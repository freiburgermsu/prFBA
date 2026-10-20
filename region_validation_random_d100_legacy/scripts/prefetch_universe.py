#!/usr/bin/env python
"""
prefetch_universe.py — prefetch PGFam sets for exactly the genomes the gene-capture
scorer will consult, BEFORE score_genecap runs, so it never lazy-fetches each miss +
rewrites the whole cache per miss (an O(N^2) death-spiral; the pilot lesson).

BOUNDED (not the full top-20 universe): the set needed is
  (a) source genomes (truth T per amplicon),
  (b) candidate_genome_ids(hits)            -> the include-self selection's stage-7 candidates,
  (c) candidate_genome_ids(filtered hits)   -> the exclude-self (tau=0.987) re-selection's
      stage-7 candidates (drop each amplicon's self / same-taxon / >=tau hits first).
Final reps in either selection are a subset of their stage-7 candidates, so this covers
every provider() call score_genecap makes — at a fraction of the full top-20 universe
(measured ~3.3x smaller on the pilot, and it saturates far below the universe at scale).

Threaded, resumable (skips cached), atomic incremental writes (prefetch_gene_families).

--select-mode {exact-tie,legacy} (default exact-tie): which reducer's stage-7 candidates
(b) and (c) enumerate. Pass legacy when selection / score_genecap run with the legacy
ordered reducer (select_references --select-mode legacy, score_genecap --legacy); the
default enumerates the exact-tie union only, and every legacy-only candidate would then be
lazy-fetched by score_genecap with a whole-cache rewrite per miss.
Saves every SAVE_EVERY=2000 fetches (each save rewrites the multi-GB shared cache).

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import prefetch_gene_families as PF  # noqa: E402
import select_references as S  # noqa: E402  (candidate_genome_ids lives in PF; S kept for parity)
import _common  # noqa: E402

TAU = 0.987  # gene-capture exclude-self threshold (matches score_genecap self_mode)
SAVE_EVERY = 2000  # was 200; every save rewrites the whole shared cache


def _filtered_hits(hits):
    """Drop each amplicon's self / same-taxon / >=TAU-identity hits, re-rank, recompute
    best_identity (never None) — the exact record the exclude-self re-selection scores."""
    out = {}
    for key, rec in hits.items():
        src = key.split("__", 1)[0]
        try:
            src_tax = int(src.split(".")[0])
        except Exception:
            src_tax = None
        kept = [
            h for h in rec.get("top20", [])
            if str(h.get("genome_id")) != src
            and h.get("taxon_id") != src_tax
            and (h.get("identity") or 0.0) < TAU
        ]
        if not kept:
            continue
        fr = dict(rec)
        fr["top20"] = [{**h, "rank": i + 1} for i, h in enumerate(kept)]
        fr["best_identity"] = kept[0].get("identity") or 0.0
        out[key] = fr
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--select-mode", default="exact-tie", choices=["exact-tie", "exact_tie", "legacy"],
                    help="reducer whose stage-7 candidates are enumerated (default exact-tie)")
    args = ap.parse_args()
    knobs = PF.knobs_for_select_mode(args.select_mode)
    print(f"[universe] select_mode={args.select_mode} knobs={knobs}", flush=True)
    hits = json.load(open(_common.PATHS.hits_json))
    t = time.perf_counter()
    need = set(PF.candidate_genome_ids(hits, knobs=knobs))                  # include-self candidates
    n_incl = len(need)
    need |= set(PF.candidate_genome_ids(_filtered_hits(hits), knobs=knobs))  # exclude-self candidates
    need |= {k.split("__", 1)[0] for k in hits}                     # source genomes (truth)
    need = sorted(g for g in need if g)
    print(f"[universe] bounded rep set = {len(need)} genomes "
          f"(include-self candidates {n_incl} + exclude-self candidates + sources; "
          f"computed in {time.perf_counter()-t:.0f}s)", flush=True)

    cache_path = _common.PATHS.pgfam_cache
    cache = PF.load_cache(cache_path)
    todo = [g for g in need if g not in cache]
    print(f"[universe] cached={len(need) - len(todo)} to_fetch={len(todo)} workers=16", flush=True)
    if not todo:
        print("[universe] nothing to fetch (resumed).")
        return
    t0 = time.perf_counter()
    done = empty = errors = 0
    with ThreadPoolExecutor(max_workers=16) as ex:
        futs = {ex.submit(PF.fetch_pgfams, g): g for g in todo}
        for fut in as_completed(futs):
            g = futs[fut]
            try:
                cache[g] = fut.result()
                empty += not cache[g]
            except Exception:
                # Never persist a failure: [] reads downstream as a genuinely
                # gene-less genome, and this cache is shared across panels.
                pass
                errors += 1
            done += 1
            if done % SAVE_EVERY == 0 or done == len(todo):
                PF.save_cache(cache, cache_path)  # atomic (tmp + replace)
                el = time.perf_counter() - t0
                print(f"[universe] {done}/{len(todo)} ({el:.0f}s, {done/max(el,1):.1f}/s) "
                      f"empty={empty} err={errors}", flush=True)
    PF.save_cache(cache, cache_path)
    print(f"[universe] fetched {len(todo)} in {time.perf_counter()-t0:.0f}s -> {cache_path} "
          f"(empty={empty} err={errors})")


if __name__ == "__main__":
    main()
