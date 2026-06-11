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

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import statistics as st
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import select_references as S

API = "https://www.bv-brc.org/api/genome_feature/"
# stages a hit reaches only if it entered stage 7 (the marginal-gain candidate pool)
STAGE7 = {"marginal_gain", "contamination_budget", "tier_cap", "complete", "no_gene_set"}


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


def fetch_pgfams(gid, family="pgfam_id", timeout=60, retries=3) -> list:
    """Return the sorted set of PGFam ids over a genome's CDS features (BV-BRC API)."""
    q = (f"eq(genome_id,{gid})&eq(feature_type,CDS)&select({family})"
         f"&limit(25000)&http_accept=application/json")
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(API + "?" + q, timeout=timeout) as resp:
                rows = json.load(resp)
            return sorted({r[family] for r in rows if r.get(family)})
        except Exception as exc:  # transient network / 5xx
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise last


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
                   verbose=True, progress_every=200):
    """Ensure the PGFam cache covers every stage-7 candidate genome in ``hits``; fetch the
    missing ones (threaded, resumable, incremental). Returns ``(cache_dict, candidates)``.
    Reusable by the build orchestrator so the whole pipeline stays exact + cached."""
    cand = sorted(candidate_genome_ids(hits))
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
                cache[g] = []
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
    args = ap.parse_args()

    hits = json.load(open(args.hits))
    cache, cand = populate_cache(hits, args.cache, workers=args.workers,
                                 family=args.family, limit=args.limit)
    sizes = [len(cache[g]) for g in cand if g in cache]
    nz = [s for s in sizes if s]
    print(f"  coverage: {len(nz)}/{len(sizes)} genomes have CDS gene sets; "
          f"{len(sizes)-len(nz)} empty (16S-only / no assembly -> auto-excluded by selection); "
          f"median {int(st.median(nz)) if nz else 0} PGFams/genome")


if __name__ == "__main__":
    main()
