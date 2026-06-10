#!/usr/bin/env python
"""
build_synthetic_genomes_parallel.py — maximally-parallel driver for the per-ASV synthetic
KBase genomes, reproducing the data of the prior codiffusion_bioreactor/genome_objects run
(functions + KBase aliases + RAST ontology + per-function probability, NO sequences) but
keyed per ASV.

It reuses the tested single-genome + per-ASV merge code in build_synthetic_genomes.py and
parallelizes the two expensive, independent axes on this 64-core box:

  1. EXACT selection  — build the reference-selection audit with the BV-BRC PGFam gene
     provider (exact gene-set novelty), NOT the estimator. The estimated audit on disk is
     refused by design; this regenerates the authoritative gene_data="bvbrc" audit.
  2. WARM the genome cache (threaded, network-bound) — fetch every SELECTED genome_id's
     fast source dict once (2 BV-BRC calls each), cached to bvbrc_cache/kbase_genome_cache/{gid}.json.
  3. MERGE per ASV (process pool, CPU-bound) — union the selected genomes' genes by FUNCTION
     via LocalGenomeConverter.create_synthetic_genome, tag core/accessory provenance, and
     write {asv}.json + genome_objects/{asv}.json + taxonomies/{asv}.json + provenance/{asv}.json.

Caches are warmed BEFORE the process pool so workers are read-only on the cache (no fetch
race). Abstained ASVs (n_selected == 0) get a manifest entry, no genome.

    python build_synthetic_genomes_parallel.py \
        --hits    /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_top20_alignment_hits.json \
        --out-dir /home/freiburger/Documents/EmilyKin/synthetic_genomes \
        --procs 56 --io-workers 32          # drop nothing for the full 3,950-ASV run
        # --limit 5                          # smoke test on the first 5 ASVs

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

import build_synthetic_genomes as BSG
import select_references as S

# ---- worker globals (set per process in _winit) ----
_CACHE_DIR = None


def _winit(out_dir, cache_dir):
    """Per-worker: silence the headless matplotlib taxonomy viz, contain the merger's
    hardcoded relative side-effect dirs under out_dir, and stash the genome cache dir."""
    global _CACHE_DIR
    _CACHE_DIR = cache_dir
    BSG._patch_taxonomy_viz()
    os.chdir(out_dir)
    for d in ("genome_objects", "taxonomies", "provenance", "ASVset_taxonomies"):
        os.makedirs(d, exist_ok=True)


def _wbuild(task):
    """Build one ASV's synthetic genome (runs in a worker, cwd == out_dir)."""
    asv, rec = task
    try:
        return asv, BSG.build_for_asv(asv, rec, genome_cache_dir=_CACHE_DIR, builder="fast")
    except Exception as exc:  # never let one ASV kill the pool
        return asv, {"status": "failed", "tier": rec.get("tier"), "error": str(exc)[:300]}


def exact_select(hits, gene_cache):
    """Authoritative selection with exact gene-set novelty (BV-BRC PGFam provider)."""
    provider = S.bvbrc_gene_provider(gene_cache)  # loads the warm 4.5k-genome PGFam cache
    print(f"[select] exact gene-set novelty over {len(hits)} ASVs "
          f"(gene cache: {gene_cache}) ...", flush=True)
    t0 = time.perf_counter()
    audit, summary = S.build_audit(hits, gene_provider=provider)
    print(f"[select] done in {time.perf_counter()-t0:.0f}s | "
          f"abstain={summary['abstained']} mean_handful={summary['mean_handful']} "
          f"tiers={summary['tiers']}", flush=True)
    return audit, summary


def warm_cache(gids, cache_dir, io_workers):
    """Fetch every selected genome's fast source dict once (threaded, resumable)."""
    os.makedirs(cache_dir, exist_ok=True)
    todo = [g for g in gids if not os.path.exists(os.path.join(cache_dir, f"{g}.json"))]
    print(f"[warm] {len(gids)} selected genomes | cached={len(gids)-len(todo)} "
          f"to_fetch={len(todo)} threads={io_workers}", flush=True)
    if not todo:
        return 0, 0
    t0 = time.perf_counter()
    done = err = 0
    with ThreadPoolExecutor(max_workers=io_workers) as ex:
        futs = {ex.submit(BSG.build_genome, g, cache_dir, "fast"): g for g in todo}
        for fut in as_completed(futs):
            try:
                fut.result()
            except Exception:
                err += 1
            done += 1
            if done % 200 == 0 or done == len(todo):
                el = time.perf_counter() - t0
                print(f"[warm] {done}/{len(todo)} ({el:.0f}s, {done/max(el,1):.1f}/s) err={err}",
                      flush=True)
    print(f"[warm] fetched {len(todo)} in {time.perf_counter()-t0:.0f}s (err={err})", flush=True)
    return len(todo), err


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True, help="asv_top20_alignment_hits.json (SW hits)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--gene-cache", default="bvbrc_cache/genome_gene_families.json",
                    help="warm PGFam cache from prefetch_gene_families.py (exact novelty)")
    ap.add_argument("--genome-cache-dir", default="bvbrc_cache/kbase_genome_cache")
    ap.add_argument("--selection-out", default=None,
                    help="where to write the exact audit (default: alongside --hits)")
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 8) - 4),
                    help="merge worker processes")
    ap.add_argument("--io-workers", type=int, default=32, help="cache-warm threads")
    ap.add_argument("--limit", type=int, default=0, help="first N ASVs (smoke test)")
    args = ap.parse_args()

    hits = json.load(open(os.path.abspath(args.hits)))
    if args.limit:
        hits = {a: hits[a] for a in list(hits)[:args.limit]}
    out_dir = os.path.abspath(args.out_dir)
    cache_dir = os.path.abspath(args.genome_cache_dir)
    gene_cache = os.path.abspath(args.gene_cache)
    selection_out = os.path.abspath(args.selection_out) if args.selection_out else (
        os.path.join(os.path.dirname(os.path.abspath(args.hits)), "asv_reference_selection.json"))

    # ---- 1. exact selection ----
    audit, summary = exact_select(hits, gene_cache)
    if not args.limit:  # only persist the authoritative full audit; preserve any estimated one
        if os.path.exists(selection_out):
            try:
                old = json.load(open(selection_out))
                if (old.get("_meta", {}) or {}).get("gene_data") == "estimated":
                    bak = selection_out.replace(".json", ".estimated.json")
                    shutil.move(selection_out, bak)
                    print(f"[select] backed up prior estimated audit -> {bak}", flush=True)
            except Exception:
                pass
        with open(selection_out, "w") as f:
            json.dump({"_meta": {"gene_data": "bvbrc",
                                 "generated_from": os.path.abspath(args.hits),
                                 "summary": summary}, **audit}, f, indent=1)
        print(f"[select] wrote exact audit -> {selection_out}", flush=True)

    # ---- 2. warm the genome cache for every selected genome (threaded) ----
    gids = sorted({str(s["genome_id"]) for r in audit.values() for s in r["selected"]})
    warm_cache(gids, cache_dir, args.io_workers)

    # ---- 3. parallel per-ASV merge ----
    os.makedirs(out_dir, exist_ok=True)
    for d in ("genome_objects", "taxonomies", "provenance", "ASVset_taxonomies"):
        os.makedirs(os.path.join(out_dir, d), exist_ok=True)
    tasks = [(a, audit[a]) for a in audit if audit[a]["n_selected"] > 0]
    abstains = {a: {"status": "abstain", "tier": audit[a]["tier"], "flags": audit[a]["flags"],
                    "midas_taxonomy": audit[a].get("midas_taxonomy")}
                for a in audit if audit[a]["n_selected"] == 0}
    print(f"[build] {len(tasks)} ASVs to build | {len(abstains)} abstain | "
          f"procs={args.procs} | out={out_dir}", flush=True)

    manifest, built, failed = {}, 0, 0
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.procs, initializer=_winit,
                             initargs=(out_dir, cache_dir)) as ex:
        for i, (asv, entry) in enumerate(ex.map(_wbuild, tasks, chunksize=4), 1):
            manifest[asv] = entry
            built += entry["status"] == "built"
            failed += entry["status"] == "failed"
            if i % 100 == 0 or i == len(tasks):
                print(f"[build] {i}/{len(tasks)} built={built} failed={failed} "
                      f"({time.perf_counter()-t0:.0f}s, {i/max(time.perf_counter()-t0,1):.1f}/s)",
                      flush=True)

    manifest.update(abstains)
    man_path = os.path.join(out_dir, "synthetic_genome_manifest.json")
    with open(man_path, "w") as f:
        json.dump({"_meta": {"selection": selection_out if not args.limit else "(limit run)",
                             "gene_data": "bvbrc", "builder": "fast",
                             "n_asvs": len(audit), "built": built,
                             "abstain": len(abstains), "failed": failed}, **manifest}, f, indent=1)
    dt = time.perf_counter() - t0
    print(f"[build] DONE: built={built} abstain={len(abstains)} failed={failed} in {dt:.0f}s "
          f"-> {man_path}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
