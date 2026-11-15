#!/usr/bin/env python
"""
prefetch_gene_families.py — populate the BV-BRC PGFam gene-set cache for exactly the
genomes that ASV selection consults, so `select_references.py --gene-provider bvbrc`
computes EXACT gene-set novelty (strategies 1 & 3) instead of the taxonomy estimator.

Only the genome_ids that survive the gates to the marginal-gain stage (stages 0-6 are
gene-independent, so this candidate set is fixed) are fetched — ~4.5k of the ~21.7k
top-20 genome_ids (a ~79% saving). Each genome's CDS PGFam (global protein family) set
is fetched once and cached to a JSON the provider reads. The run is threaded, resumable
(skips cached ids), and writes incrementally.

    python prefetch_gene_families.py \
        --hits /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_top20_alignment_hits.json \
        --cache bvbrc_cache/genome_gene_families.json
    # then:  select_references.py ... --gene-provider bvbrc --gene-cache genome_gene_families.json

SELECTION RULE (--select-mode): the stage-7 candidate set depends on the reducer. The
default (exact-tie) enumerates the candidates of select_references' shipped default, the
exact-tie union. `--select-mode legacy` enumerates them under the LEGACY ordered reducer
(knobs {"select_all": False}) instead -- the set a `select_references.py --select-mode
legacy` run will consult. Without it, a legacy selection lazy-fetches every miss and
rewrites the whole shared cache once per miss.

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import statistics as st
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import select_references as S

API = "https://www.bv-brc.org/api/genome_feature/"
# stages a hit reaches only if it entered stage 7 (the marginal-gain candidate pool)
STAGE7 = {"marginal_gain", "contamination_budget", "tier_cap", "complete", "no_gene_set"}


def knobs_for_select_mode(select_mode):
    """'exact-tie' / 'exact_tie' / None -> None (the shipped default); 'legacy' -> the
    legacy ordered reducer's knobs. Mirrors select_references.main's --select-mode."""
    if select_mode in (None, "exact-tie", "exact_tie"):
        return None
    if select_mode == "legacy":
        return {"select_all": False}
    raise ValueError(f"unsupported select_mode {select_mode!r} (exact-tie | legacy)")


def candidate_genome_ids(hits, knobs=None) -> set:
    """Distinct genome_ids that reach the marginal-gain stage across all ASVs (the only
    genomes whose gene sets affect any decision)."""
    audit, _ = S.build_audit(hits, knobs=knobs)
    cand = set()
    for r in audit.values():
        for d in r["decisions"]:
            if d["disposition"] == "selected" or d["stage"] in STAGE7:
                cand.add(str(d["genome_id"]))
    return cand


ANNOTATION = "PATRIC"   # BV-BRC's own (RASTtk) annotation; never the imported RefSeq one
# BV-BRC 403s the bare urllib User-Agent (verified 2026-10-09); every request sends this one.
USER_AGENT = "prFBA/1.0 (academic; metabolic-model benchmarking)"
# Rows per request: limit(PAGE_SIZE,start). Up to 2026-11-13 every CDS query was ONE request with
# limit(25000) and no paging, which silently truncated any genome over 25,000 rows (see
# fetch_patric_cds_rows). A module constant so tests can shrink it.
PAGE_SIZE = 25000


class IncompleteFetch(RuntimeError):
    """A paged BV-BRC fetch whose rows could not be shown to be complete. Raised instead of
    returning a partial row set; the caller must not cache anything for that genome."""


def _content_range_total(value):
    """The total of a BV-BRC Content-Range header ('items 0-24999/39609' -> 39609; 'items
    0-0/0' -> 0), or None when the header is absent or carries no numeric total."""
    if not value or "/" not in value:
        return None
    tail = value.rsplit("/", 1)[1].strip()
    return int(tail) if tail.isdigit() else None


def _get(url, timeout=60, retries=3):
    """One BV-BRC GET with USER_AGENT: (parsed JSON body, Content-Range total or None).
    Transient errors (network, 5xx) are retried with backoff; after `retries` failed attempts
    the last error propagates. The transport of fetch_patric_cds_rows (tests replace it)."""
    last = None
    for attempt in range(retries):
        try:
            rq = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(rq, timeout=timeout) as resp:
                body = json.load(resp)
                headers = getattr(resp, "headers", None)
                total = _content_range_total(headers.get("Content-Range") if headers else None)
            return body, total
        except Exception as exc:  # transient network / 5xx
            last = exc
            if isinstance(exc, urllib.error.HTTPError):
                exc.close()  # an HTTPError holds the open response; leaking it exhausts fds
            time.sleep(1.5 * (attempt + 1))
    raise last


def patric_cds_query(gid, fields):
    """The genome_feature query (without paging) for a genome's BV-BRC PATRIC CDS rows,
    selecting feature_id (needed by the overlap guard) plus `fields`, sorted by feature_id so
    that offset pages are stable (an unsorted capped query returned a different subset on
    different calls; 2026-11-13)."""
    sel = ",".join(["feature_id"] + [f for f in fields if f != "feature_id"])
    return (f"eq(genome_id,{gid})&eq(feature_type,CDS)&eq(annotation,{ANNOTATION})"
            f"&select({sel})&sort(+feature_id)")


def fetch_patric_cds_rows(gid, fields, *, page_size=None, timeout=60, retries=3) -> list:
    """Every BV-BRC PATRIC-annotation CDS row of a genome (eq(feature_type,CDS) and
    eq(annotation,PATRIC); each row has feature_id plus `fields`), fetched in pages of
    limit(page_size,start) (default PAGE_SIZE = 25,000) sorted by feature_id, until a short page
    or the server's Content-Range total is reached.

    Never returns a partial set. It raises when:
      - a page fails after `retries` attempts (the transport error propagates);
      - a response has no Content-Range total, or the total changes between pages
        (IncompleteFetch);
      - the rows fetched differ in number from that total, e.g. because the server returned a
        short page early (IncompleteFetch);
      - a feature_id is missing or repeated, i.e. pages overlap (IncompleteFetch).
    A genome with no PATRIC CDS (total 0) returns [].

    Why: up to 2026-11-13 fast_build_genome and fetch_pgfams sent ONE unpaged request with
    limit(25000), and until that day without the annotation filter, so RefSeq rows also counted
    toward it. Every genome over that window was silently truncated: 13 genome-cache files, 167
    entries of the shared PGFam cache, and copies of them. The 2026-11-13 repair replaced those 13
    files (and 2 more holding RefSeq rows with another genome's patric_id: 13 + 2) and 1,006
    shared-cache entries (the 167 truncated ones plus 839 empty lists that were stored fetch
    failures), and 22 / 35 entries of the EmilyKin / codiffusion d100 tie-set caches; every change
    is logged, with the old value, in bvbrc_cache/CORRECTIONS_2026-11-13.json."""
    page_size = page_size or PAGE_SIZE
    base = patric_cds_query(gid, fields)
    rows, start, total = [], 0, None
    while True:
        url = f"{API}?{base}&limit({page_size},{start})&http_accept=application/json"
        page, n = _get(url, timeout=timeout, retries=retries)
        where = f"genome {gid}, rows {start}+ of {base!r}"
        if not isinstance(page, list):
            raise IncompleteFetch(f"{where}: response is not a list of rows")
        if n is None:
            raise IncompleteFetch(f"{where}: no Content-Range total, so completeness cannot be "
                                  f"checked")
        if total is None:
            total = n
        elif n != total:
            raise IncompleteFetch(f"{where}: Content-Range total changed between pages "
                                  f"({total} -> {n})")
        if len(page) > page_size:
            raise IncompleteFetch(f"{where}: page of {len(page)} rows > limit {page_size}")
        rows.extend(page)
        start += len(page)
        if len(page) < page_size or start >= total:
            break
    if len(rows) != total:
        raise IncompleteFetch(f"genome {gid}: fetched {len(rows)} rows but the server reports "
                              f"{total} ({base!r}); nothing returned")
    ids = [r.get("feature_id") if isinstance(r, dict) else None for r in rows]
    if None in ids or len(set(ids)) != len(ids):
        raise IncompleteFetch(f"genome {gid}: {len(ids) - len(set(ids))} repeated or missing "
                              f"feature_id over {len(ids)} rows (pages overlap); nothing returned")
    return rows


GENOME_SEQUENCE_API = "https://www.bv-brc.org/api/genome_sequence/"


def genome_sequence_query(gid, fields):
    """The genome_sequence query (without paging) for a genome's contigs, selecting sequence_id
    (the unique key the overlap guard needs) plus `fields`, sorted by sequence_id so offset pages
    are stable."""
    sel = ",".join(["sequence_id"] + [f for f in fields if f != "sequence_id"])
    return f"eq(genome_id,{gid})&select({sel})&sort(+sequence_id)"


def fetch_genome_sequence_rows(gid, fields=("accession", "sequence"), *, page_size=None,
                               timeout=180, retries=4) -> list:
    """Every BV-BRC genome_sequence row (contig) of a genome, each with sequence_id plus `fields`,
    paged and checked exactly as fetch_patric_cds_rows pages CDS rows: limit(page_size,start)
    sorted by sequence_id until a short page or the Content-Range total; raises IncompleteFetch
    rather than return a partial set (no total, a total that changes, a row count that differs
    from it, or a repeated / missing sequence_id). A genome with no contigs returns [].

    Why (2026-11-15, re-score decision D6): the MAG-arm contig fetchers sent ONE request with
    limit(50000), and the server caps a response at 25,000 rows, so genomes with more contigs were
    written silently truncated: random-panel 1733.9197 (31,169 contigs on BV-BRC) and 1733.9535
    (36,001), and the sludge-DB 1733.9522 (25,454) and 1733.9680 (36,062) all hold exactly 25,000.
    Rows come back sorted by sequence_id (the old unpaged order was the server's)."""
    page_size = page_size or PAGE_SIZE
    base = genome_sequence_query(gid, fields)
    rows, start, total = [], 0, None
    while True:
        url = f"{GENOME_SEQUENCE_API}?{base}&limit({page_size},{start})&http_accept=application/json"
        page, n = _get(url, timeout=timeout, retries=retries)
        where = f"genome {gid}, contigs {start}+ of {base!r}"
        if not isinstance(page, list):
            raise IncompleteFetch(f"{where}: response is not a list of rows")
        if n is None:
            raise IncompleteFetch(f"{where}: no Content-Range total, so completeness cannot be "
                                  f"checked")
        if total is None:
            total = n
        elif n != total:
            raise IncompleteFetch(f"{where}: Content-Range total changed between pages "
                                  f"({total} -> {n})")
        if len(page) > page_size:
            raise IncompleteFetch(f"{where}: page of {len(page)} rows > limit {page_size}")
        rows.extend(page)
        start += len(page)
        if len(page) < page_size or start >= total:
            break
    if len(rows) != total:
        raise IncompleteFetch(f"genome {gid}: fetched {len(rows)} contigs but the server reports "
                              f"{total} ({base!r}); nothing returned")
    ids = [r.get("sequence_id") if isinstance(r, dict) else None for r in rows]
    if None in ids or len(set(ids)) != len(ids):
        raise IncompleteFetch(f"genome {gid}: {len(ids) - len(set(ids))} repeated or missing "
                              f"sequence_id over {len(ids)} contigs (pages overlap); nothing returned")
    return rows


def fetch_pgfams(gid, family="pgfam_id", timeout=60, retries=3) -> list:
    """Return the sorted set of PGFam ids over a genome's CDS features (BV-BRC API), BV-BRC
    PATRIC annotation only (eq(annotation,PATRIC), since 2026-11-13), over ALL of the genome's
    PATRIC CDS rows (fetch_patric_cds_rows: paged, raises rather than return a partial set).

    For a genome whose CDS rows fit the old 25,000-row window, the annotation filter does not
    change the values: RefSeq CDSs carry no pgfam_id (checked 2026-11-13 on 200 random genomes of
    the shared cache, none over the window: 199 equal the cached list; the 200th was cached as []
    and BV-BRC now returns 1,799 PGFams for it with and without the filter; no RefSeq row carried
    a pgfam_id). Above the window the cached values were wrong: up to 2026-11-13 this was one
    unfiltered request with limit(25000) and no paging, so a genome whose CDS rows over BOTH
    annotations exceeded 25,000 got a truncated set (RefSeq rows filled part of the window, or
    the PATRIC rows alone exceeded it). In the shared cache that was 167 genomes (105 by RefSeq
    rows, 62 over 25,000 PATRIC CDSs), missing a median 296 PGFams (6.3%), at most 5,417 (813.301)
    and up to 46.9% (165433.68); e.g. 1773.17156: 6,521 cached vs 7,771 complete. Selection (gene
    counts, novelty) and the validations read these sets, so for those genomes they read
    truncated sets. The 2026-11-13 repair replaced 1,006 shared-cache entries atomically with this
    function's values: the 167 truncated ones and 839 empty lists that were stored fetch failures
    (the prefetchers of that time persisted any failed fetch as []; BV-BRC has PGFams for all 839;
    the other 963 [] entries are BV-BRC's value and were kept), plus 22 / 35 entries of the
    EmilyKin / codiffusion d100 tie-set caches. Every change is logged in
    bvbrc_cache/CORRECTIONS_2026-11-13.json, which keeps each replaced list under
    backup_old_values. (RefSeq rows carry no pgfam_id, so no PGFam set holds a RefSeq value, also
    for the genomes whose RefSeq rows carry another genome's patric_id.)"""
    rows = fetch_patric_cds_rows(gid, [family], timeout=timeout, retries=retries)
    return sorted({r[family] for r in rows if r.get(family)})


def _open(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def load_cache(path) -> dict:
    if os.path.exists(path):
        with _open(path, "rt") as f:
            return json.load(f)
    return {}


def save_cache(cache, path):
    tmp = path + ".tmp"
    with _open(tmp, "wt") as f:
        json.dump(cache, f)
    os.replace(tmp, path)


def populate_cache(hits, cache_path, *, workers=12, family="pgfam_id", limit=0,
                   verbose=True, progress_every=2000, select_mode=None):
    """Ensure the PGFam cache covers every stage-7 candidate genome in ``hits``; fetch the
    missing ones (threaded, resumable, incremental). Returns ``(cache_dict, candidates)``.
    Reusable by the build orchestrator so the whole pipeline stays exact + cached.

    ``select_mode`` picks the reducer whose candidates are enumerated (None/'exact-tie' =
    the shipped default, 'legacy' = the ordered reducer). ``progress_every`` is also the
    save cadence: every save rewrites the whole (multi-GB) cache, so it was raised from
    200 to 2000 -- a crash loses at most that many fetches, which the next run retries."""
    cand = sorted(candidate_genome_ids(hits, knobs=knobs_for_select_mode(select_mode)))
    cache = load_cache(cache_path)
    todo = [g for g in cand if g not in cache]
    if limit:
        todo = todo[:limit]
    if verbose:
        print(f"[prefetch] {len(cand)} candidate genome_ids | cached={len(cand) - len(todo)} "
              f"to_fetch={len(todo)} workers={workers}")
    if not todo:
        return cache, cand
    t0 = time.perf_counter()
    done = empty = errors = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fetch_pgfams, g, family): g for g in todo}
        for fut in as_completed(futs):
            g = futs[fut]
            try:
                cache[g] = fut.result()
                empty += not cache[g]
            except Exception:
                # Never persist a failure: [] is indistinguishable downstream from a
                # genuinely gene-less genome, and this cache is shared across panels.
                # Leaving the key absent lets the next run retry it.
                errors += 1
            done += 1
            if verbose and done % progress_every == 0:
                save_cache(cache, cache_path)
                el = time.perf_counter() - t0
                print(f"[prefetch] {done}/{len(todo)} ({el:.0f}s, {done/max(el,1):.1f}/s) "
                      f"empty={empty} err={errors}", flush=True)
    save_cache(cache, cache_path)
    if verbose:
        print(f"[prefetch] fetched {len(todo)} in {time.perf_counter()-t0:.0f}s "
              f"-> {cache_path} (empty={empty} err={errors})")
    return cache, cand


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--cache", default="bvbrc_cache/genome_gene_families.json",
                    help="output cache (.json or .json.gz); read by select_references bvbrc provider")
    ap.add_argument("--family", default="pgfam_id", choices=["pgfam_id", "plfam_id"])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=0, help="fetch only the first N uncached candidates (sample)")
    ap.add_argument("--select-mode", default="exact-tie", choices=["exact-tie", "exact_tie", "legacy"],
                    help="reducer whose stage-7 candidates are prefetched (default exact-tie = "
                         "select_references' default; legacy = --select-mode legacy)")
    ap.add_argument("--save-every", type=int, default=2000,
                    help="rewrite the cache every N fetches (default 2000; each save rewrites the "
                         "whole file)")
    args = ap.parse_args()

    hits = json.load(open(args.hits))
    print(f"[prefetch] select_mode={args.select_mode}", flush=True)
    cache, cand = populate_cache(hits, args.cache, workers=args.workers,
                                 family=args.family, limit=args.limit,
                                 progress_every=args.save_every, select_mode=args.select_mode)
    sizes = [len(cache[g]) for g in cand if g in cache]
    nz = [s for s in sizes if s]
    print(f"  coverage: {len(nz)}/{len(sizes)} genomes have CDS gene sets; "
          f"{len(sizes)-len(nz)} empty (16S-only / no assembly -> contribute 0 genes; the "
          f"legacy reducer drops them, the default union keeps them as taxonomic reps flagged "
          f"gene_set_empty); median {int(st.median(nz)) if nz else 0} PGFams/genome")


if __name__ == "__main__":
    main()
