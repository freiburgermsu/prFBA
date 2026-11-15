#!/usr/bin/env python
"""S8 — prefetch_pgfams.py: two-pass PGFam gene-set prefetch (DESIGN §7).

A thin, resumable wrapper over ``prefetch_gene_families.py`` (the prFBA module
that owns ``fetch_pgfams`` / ``populate_cache``).  It populates the single
prFBA-wide PGFam cache (``PATHS.pgfam_cache``, ``gid -> sorted pgfam list``) for
exactly the genomes this validation study consults, so both truth and predicted
gene sets resolve out of *one* namespace (apples-to-apples recall/precision in
``score_genecap.py``).

Two passes into the same cache (DESIGN §7 "Threaded, resumable prefetch"):

  1. **truth** — the 10k source ``genome_id``s in ``PATHS.selection_json``
     (the S1 selection manifest).
  2. **predicted reps** — the *distinct selected* representative ``genome_id``s
     in ``PATHS.selection_out`` (the S6 include-self selection), if present.
     Skipped with a note when S6 has not run yet (run me again after S6).

Both passes share one resumable cache: already-cached genomes are skipped, so a
re-run only fetches the gap.  Each genome is fetched from the **hyphenated**
host (``https://www.bv-brc.org``) with ``feature_type=CDS``,
``annotation=PATRIC`` (BV-BRC's own annotation, never the imported RefSeq one),
``select(feature_id,pgfam_id)``, sorted by feature_id and **paged**
(``limit(25000,start)`` until a short page; ``prefetch_gene_families.
fetch_patric_cds_rows``, which raises rather than return a partial set), 12
threads, 3x retry/backoff per page, atomic save-every-200.

The cap flag (since 2026-11-13).  Up to 2026-11-13 this was ONE unfiltered
``limit(25000)`` request with no paging, and a genome hitting the cap was
flagged as truncated.  Now every genome is fetched in full, so nothing is
truncated; a genome is still flagged when its PATRIC CDS rows exceed one
25000-row page, because any PGFam set fetched for it before 2026-11-13 was
truncated.  The flag covers only genomes fetched in this run (cached ones are
skipped), and not those truncated before 2026-11-13 by RefSeq rows filling the
old unfiltered window (PATRIC <= 25000 < both annotations).  The truncated
entries of the shared cache are repaired, each change logged in
``bvbrc_cache/CORRECTIONS_2026-11-13.json``.

Resumable: if every required genome is already cached, the pass short-circuits.
``--force`` re-fetches every genome (truth + reps) even if cached.

Sole interpreter: ``/home/freiburger/Documents/py_venv/bin/python``.

    ~/Documents/py_venv/bin/python region_validation/scripts/prefetch_pgfams.py
    ~/Documents/py_venv/bin/python region_validation/scripts/prefetch_pgfams.py --force
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

# --- import contract module + the prFBA prefetch machinery ------------------ #
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # for _common
import _common as C  # noqa: E402  (PATHS / done — never hardcode a data path)

# prefetch_gene_families lives at the prFBA repo root; reuse its fetch helpers
# (fetch_pgfams / load_cache / save_cache / API) rather than re-implementing.
_PRFBA_ROOT = "/home/freiburger/Documents/prFBA"
sys.path.insert(0, _PRFBA_ROOT)
import prefetch_gene_families as PGF  # noqa: E402

# Pinned, design-mandated fetch parameters (DESIGN §7 "PGFam fetch plan").
HOST = PGF.API                 # hyphenated host: https://www.bv-brc.org/api/genome_feature/
FAMILY = "pgfam_id"
LIMIT = 25000                  # rows per page; a genome over one page is flagged (fetched in full)
WORKERS = 12                   # BV-BRC Solr is the limiter, not cores (DESIGN §7)
RETRIES = 3
UA = "prFBA/1.0 (academic; metabolic-model benchmarking)"   # same string as prefetch_gene_families.py
SAVE_EVERY = 200
PROGRESS_EVERY = 200

_GID_RE = re.compile(r"^\d+\.\d+$")


# --------------------------------------------------------------------------- #
# Genome-id extraction from the two input files
# --------------------------------------------------------------------------- #
def _valid_gids(values):
    """Filter an iterable of candidate ids to canonical BV-BRC genome_ids."""
    out = set()
    for v in values:
        s = str(v).strip()
        if _GID_RE.match(s):
            out.add(s)
    return out


def source_genome_ids(path):
    """Distinct source genome_ids from the S1 selection manifest (``selection_json``).

    Tolerant of the manifest shape S1 may emit (DESIGN §2.6 lists per-genome
    records ``{genome_id, taxon_id, ...}``): accepts a list of records, a dict of
    records, or a dict keyed by genome_id.  Falls back to harvesting every
    ``genome_id`` field anywhere in the structure.
    """
    with open(path) as fh:
        data = json.load(fh)

    gids = set()

    def _harvest(obj, key=None):
        if isinstance(obj, dict):
            if "genome_id" in obj:
                gids.update(_valid_gids([obj["genome_id"]]))
            for k, v in obj.items():
                if k == "_meta":
                    continue
                # a dict keyed BY genome_id (key itself is the id)
                gids.update(_valid_gids([k]))
                _harvest(v, k)
        elif isinstance(obj, list):
            for item in obj:
                _harvest(item)

    _harvest(data)
    return gids


def rep_genome_ids(path):
    """Distinct *selected representative* genome_ids from the S6 selection output.

    The selection JSON (``select_references.py`` output, mirrored by S6) is keyed
    by amplicon (plus a ``_meta`` entry); each amplicon record carries a
    ``selected`` list whose items each have a ``genome_id``.  We union those.
    Falls back to harvesting ``genome_id`` from any ``selected`` list found.
    """
    with open(path) as fh:
        data = json.load(fh)

    gids = set()
    records = data.values() if isinstance(data, dict) else data
    for rec in records:
        if not isinstance(rec, dict):
            continue
        sel = rec.get("selected")
        if isinstance(sel, list):
            gids.update(_valid_gids(r.get("genome_id") for r in sel
                                    if isinstance(r, dict) and r.get("genome_id")))
    return gids


# --------------------------------------------------------------------------- #
# Fetch with the explicit row count, so a genome over one 25000-row page is flagged
# --------------------------------------------------------------------------- #
def fetch_pgfams_capaware(gid, *, family=FAMILY, limit=LIMIT, timeout=60, retries=RETRIES):
    """Fetch a genome's CDS PGFam set, returning ``(sorted_pgfams, n_rows, capped)``.

    All of the genome's BV-BRC PATRIC CDS rows (``eq(annotation,PATRIC)``), in
    pages of ``limit`` rows sorted by feature_id: the prFBA pager
    ``prefetch_gene_families.fetch_patric_cds_rows`` (same host and User-Agent
    as ``fetch_pgfams``), which raises -- nothing is cached -- unless the rows
    provably equal the server's Content-Range total.  ``n_rows`` is the number of
    PATRIC CDS rows.  ``capped`` is ``n_rows > limit``: the set is complete, but
    the genome exceeds one page, so a set fetched for it before 2026-11-13 by the
    old single unfiltered ``limit(25000)`` request was truncated.  (Until
    2026-11-13 ``capped`` was ``n_rows >= limit`` of that one request, and the
    set itself was truncated.)
    """
    rows = PGF.fetch_patric_cds_rows(gid, [family], page_size=limit, timeout=timeout,
                                     retries=retries)
    pgfams = sorted({r[family] for r in rows if r.get(family)})
    return pgfams, len(rows), (len(rows) > limit)


# --------------------------------------------------------------------------- #
# One resumable, threaded pass over an explicit genome-id set
# --------------------------------------------------------------------------- #
def run_pass(label, gids, cache, cache_path, *, force=False):
    """Fetch every uncached gid into ``cache`` (threaded, save-every-200, atomic).

    Returns ``(n_fetched, capped_gids)``.  ``cache`` is mutated in place and
    flushed to ``cache_path`` incrementally + at the end.  Resumable: skips ids
    already in ``cache`` unless ``force``.
    """
    gids = sorted(gids)
    todo = gids if force else [g for g in gids if g not in cache]
    cached = len(gids) - len(todo)
    print(f"[{label}] {len(gids)} genome_ids | cached={cached} to_fetch={len(todo)} "
          f"workers={WORKERS}", flush=True)
    if not todo:
        print(f"[{label}] nothing to fetch (resumed).", flush=True)
        return 0, []

    capped = []
    t0 = time.perf_counter()
    done = empty = errors = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(fetch_pgfams_capaware, g): g for g in todo}
        for fut in as_completed(futs):
            g = futs[fut]
            try:
                pgfams, n_rows, is_capped = fut.result()
                cache[g] = pgfams
                empty += not pgfams
                if is_capped:
                    capped.append(g)
            except Exception:
                # Do NOT persist a failure as an empty result: an empty PGFam set is
                # indistinguishable from "this genome has no genes" downstream. Leave the
                # key ABSENT so the next run retries it.
                pass
                errors += 1
            done += 1
            if done % SAVE_EVERY == 0:
                PGF.save_cache(cache, cache_path)
                el = time.perf_counter() - t0
                print(f"[{label}] {done}/{len(todo)} ({el:.0f}s, {done/max(el,1):.1f}/s) "
                      f"empty={empty} err={errors} capped={len(capped)}", flush=True)
    PGF.save_cache(cache, cache_path)
    el = time.perf_counter() - t0
    print(f"[{label}] fetched {len(todo)} in {el:.0f}s -> {cache_path} "
          f"(empty={empty} err={errors} capped={len(capped)})", flush=True)
    return len(todo), capped


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="re-fetch every genome (truth + reps) even if cached")
    ap.add_argument("--source", default=C.PATHS.selection_json,
                    help="S1 selection manifest (10k source genome_ids)")
    ap.add_argument("--selection", default=C.PATHS.selection_out,
                    help="S6 selection output (distinct selected rep genome_ids)")
    ap.add_argument("--cache", default=C.PATHS.pgfam_cache,
                    help="prFBA-wide PGFam cache (gid -> sorted pgfam list)")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.cache), exist_ok=True)
    cache = PGF.load_cache(args.cache)
    print(f"[cache] {args.cache}: {len(cache)} genomes already cached", flush=True)

    # --- Pass 1: truth (10k source genomes) -------------------------------- #
    if not C.done(args.source):
        print(f"[truth] ERROR: source manifest missing: {args.source}\n"
              f"        run S1 (select_10k.py) first.", file=sys.stderr)
        sys.exit(2)
    src = source_genome_ids(args.source)
    print(f"[truth] parsed {len(src)} distinct source genome_ids from {args.source}",
          flush=True)
    _, capped_src = run_pass("truth", src, cache, args.cache, force=args.force)

    # --- Pass 2: predicted reps (distinct selected reps), if S6 ran -------- #
    reps = set()
    capped_rep = []
    if C.done(args.selection):
        reps = rep_genome_ids(args.selection)
        # truth genomes are already cached; run_pass skips them resumably.
        print(f"[reps] parsed {len(reps)} distinct selected rep genome_ids "
              f"from {args.selection}", flush=True)
        _, capped_rep = run_pass("reps", reps, cache, args.cache, force=args.force)
    else:
        print(f"[reps] selection output not present yet ({args.selection}); "
              f"skipping pass 2 — re-run after S6 (select_excludeself.py).", flush=True)

    # --- Coverage + cap report --------------------------------------------- #
    wanted = sorted(src | reps)
    sizes = [len(cache[g]) for g in wanted if g in cache]
    nz = [s for s in sizes if s]
    median = sorted(nz)[len(nz) // 2] if nz else 0
    capped_all = sorted(set(capped_src) | set(capped_rep))
    print(f"[coverage] {len(nz)}/{len(sizes)} requested genomes have a CDS PGFam set; "
          f"{len(sizes) - len(nz)} empty (16S-only / no assembly / fetch error -> "
          f"score_genecap treats as truth_empty/pgfam_missing); median {median} PGFams/genome",
          flush=True)
    if capped_all:
        print(f"[CAP] {len(capped_all)} genome(s) fetched this run have > {LIMIT} PATRIC CDS "
              f"rows (fetched in full by paging; a PGFam set fetched for them before "
              f"2026-11-13 by the unpaged query was TRUNCATED -- see "
              f"bvbrc_cache/CORRECTIONS_2026-11-13.json): "
              f"{', '.join(capped_all[:20])}{' ...' if len(capped_all) > 20 else ''}",
              flush=True)
    else:
        print(f"[CAP] no genome fetched this run has > {LIMIT} PATRIC CDS rows.", flush=True)


if __name__ == "__main__":
    main()
