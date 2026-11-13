#!/usr/bin/env python
"""
build_synthetic_genomes_parallel.py — maximally-parallel driver for the per-ASV synthetic
KBase genomes, reproducing the shape of the prior codiffusion_bioreactor/genome_objects run
(functions + KBase aliases + RAST ontology + per-function probability, NO sequences) but
keyed per ASV.

It reuses the tested single-genome + per-ASV merge code in build_synthetic_genomes.py and
parallelizes the two expensive, independent axes on this 64-core box:

  1. EXACT selection  — build the reference-selection audit with the BV-BRC PGFam gene
     provider (exact gene-set novelty), NOT the estimator. The estimated audit on disk is
     refused by design; this regenerates the authoritative gene_data="bvbrc" audit.
  2. WARM the genome cache (threaded, network-bound) — fetch every SELECTED genome_id's
     fast source dict once (2 BV-BRC calls each; BV-BRC PATRIC annotation only), cached to
     bvbrc_cache/kbase_genome_cache/{gid}.json (new files only; existing ones are filtered on read).
  3. MERGE per ASV (process pool, CPU-bound) — union the selected genomes' genes by FUNCTION
     via LocalGenomeConverter.create_synthetic_genome, tag core/accessory provenance, and
     write {asv}.json + genome_objects/{asv}.json + taxonomies/{asv}.json + provenance/{asv}.json.

Caches are warmed BEFORE the process pool so workers are read-only on the cache (no fetch
race). Abstained ASVs (n_selected == 0) get a manifest entry, no genome.

--selection AUDIT builds from a prebuilt exact audit instead of --hits (skips step 1; an
'estimated' audit is refused unless --allow-estimated, as in build_synthetic_genomes.py).

BV-BRC PATRIC annotation only (since 2026-11-13; no option to include RefSeq): every source
genome dict, cached or fetched by the warm step, passes through build_synthetic_genomes.
bvbrc_only (see that module's docstring); the manifest _meta records 'annotation'. Builds made
up to git 3ce12d0 also merged the imported RefSeq CDSs and are reproducible only from git at or
before 3ce12d0.

--pfunc e3 (the DEFAULT) replaces k/n by E3 for ASVs with >= 2 built source genomes (see
build_synthetic_genomes.py); --pfunc kn keeps k/n everywhere (genome, provenance and taxonomy
files byte-identical to the same code's build without E3; the manifest _meta gains a minimal
'pfunc' key). In every
mode the manifest _meta also gains 'selection_stamp' (sha256/size/mtime_ns of the audit; None
for a --hits --limit run, which writes no audit) and 'source_digests', the record of
OUT/pfunc/source_digests.json ({genome_id: source_digest} of every built source genome, computed
in the workers on the dicts they loaded): what reweight_build_e3.py checks for drift. An e3 run
first checks the E3 clade catalogue (exists and matches the PGFam cache; else it stops before
selection with the rebuild command). Before the merge pool, the parent collects the genomes
E3 needs (neighbours of every n_selected >= min_n ASV from the audit + clade members of every
rep, since the prior is anchored on the first rep that builds), streams their PGFam
sets out of the shared cache once, fetches the missing ones into OUT/pfunc/pgfam_fetched.json,
saves the slim store to OUT/pfunc/pgfam_store.json, and hands the store to the workers through
the pool initializer of a fork-context pool (inherited copy-on-write, read-only; nothing is
re-loaded per worker).

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
import multiprocessing
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

import build_synthetic_genomes as BSG
import select_references as S

# ---- worker globals (set per process in _winit) ----
_CACHE_DIR = None
_PFUNC = None


def _winit(out_dir, cache_dir, pfunc=None):
    """Per-worker: silence the headless matplotlib taxonomy viz, contain the merger's
    hardcoded relative side-effect dirs under out_dir, and stash the genome cache dir and the
    (read-only) pfunc context. With the fork context used in e3 mode, `pfunc` is inherited,
    not pickled."""
    global _CACHE_DIR, _PFUNC
    _CACHE_DIR = cache_dir
    _PFUNC = pfunc
    BSG._patch_taxonomy_viz()
    os.chdir(out_dir)
    for d in ("genome_objects", "taxonomies", "provenance", "ASVset_taxonomies"):
        os.makedirs(d, exist_ok=True)


def _wbuild(task):
    """Build one ASV's synthetic genome (runs in a worker, cwd == out_dir). Returns
    (asv, entry, {genome_id: source_digest} of its built source genomes)."""
    asv, rec = task
    dg = {}
    try:
        return asv, BSG.build_for_asv(asv, rec, genome_cache_dir=_CACHE_DIR, builder="fast",
                                      pfunc=_PFUNC, digests=dg), dg
    except Exception as exc:  # never let one ASV kill the pool
        return asv, {"status": "failed", "tier": rec.get("tier"), "error": str(exc)[:300]}, {}


def exact_select(hits, gene_cache, gene_tiebreak=True):
    """Authoritative selection with exact gene-set novelty (BV-BRC PGFam provider)."""
    provider = S.bvbrc_gene_provider(gene_cache)  # loads the warm 4.5k-genome PGFam cache
    print(f"[select] exact gene-set novelty over {len(hits)} ASVs "
          f"(gene cache: {gene_cache}) | P3 gene-tiebreak={'on' if gene_tiebreak else 'OFF'} ...", flush=True)
    t0 = time.perf_counter()
    # P3: gene-richest-conspecific species-dedup tie-break (fetch-free, from the loaded cache)
    gcf = S.gene_count_from_provider(provider) if gene_tiebreak else None
    audit, summary = S.build_audit(hits, gene_provider=provider, gene_count_fn=gcf)
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
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--hits", help="asv_top20_alignment_hits.json (SW hits) -> exact selection + build")
    src.add_argument("--selection", help="a prebuilt exact selection audit to build from (skips selection)")
    ap.add_argument("--allow-estimated", action="store_true",
                    help="(with --selection) permit an audit whose _meta.gene_data is 'estimated'")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--gene-cache", default="bvbrc_cache/genome_gene_families.json",
                    help="warm PGFam cache from prefetch_gene_families.py (exact novelty)")
    ap.add_argument("--genome-cache-dir", default="bvbrc_cache/kbase_genome_cache")
    ap.add_argument("--selection-out", default=None,
                    help="where to write the exact audit (default: alongside --hits)")
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 8) - 4),
                    help="merge worker processes")
    ap.add_argument("--io-workers", type=int, default=32, help="cache-warm threads")
    ap.add_argument("--no-gene-tiebreak", action="store_true",
                    help="disable the P3 gene-richest-conspecific species-dedup tie-break "
                         "(reproduces the pre-P3 'old' selection for comparison)")
    ap.add_argument("--limit", type=int, default=0, help="first N ASVs (smoke test)")
    BSG.add_pfunc_args(ap)
    args = ap.parse_args()

    out_dir = os.path.abspath(args.out_dir)
    cache_dir = os.path.abspath(args.genome_cache_dir)
    gene_cache = os.path.abspath(args.gene_cache)
    # E3 (the default) needs a current clade catalogue: fail now, not after selection/warm
    preflight = BSG.pfunc_preflight(args.pfunc, args.pfunc_catalogue, args.pfunc_pgfam_cache,
                                    allow_stale=args.pfunc_allow_stale_catalogue)
    if args.selection:
        # ---- 1'. a prebuilt audit (no selection run) ----
        selection_out = os.path.abspath(args.selection)
        sel_stamp = BSG.file_stamp(selection_out)      # the audit as read
        audit = BSG._load(selection_out)
        sel_meta = audit.pop("_meta", {}) or {}
        gene_data = sel_meta.get("gene_data", "unknown")
        if gene_data == "estimated" and not args.allow_estimated:
            raise SystemExit("[build] refusing to build from an ESTIMATED selection audit "
                             "(pass --hits for the exact pipeline, or --allow-estimated)")
        if args.limit:
            audit = {a: audit[a] for a in list(audit)[:args.limit]}
        print(f"[select] prebuilt audit {selection_out} (gene_data={gene_data}, "
              f"{len(audit)} ASVs)", flush=True)
    else:
        hits = json.load(open(os.path.abspath(args.hits)))
        if args.limit:
            hits = {a: hits[a] for a in list(hits)[:args.limit]}
        selection_out = os.path.abspath(args.selection_out) if args.selection_out else (
            os.path.join(os.path.dirname(os.path.abspath(args.hits)), "asv_reference_selection.json"))
        gene_data = "bvbrc"
        sel_stamp = None      # set once the audit is written (a --limit run writes none)

    # ---- 1. exact selection ----
    if not args.selection:
        audit, summary = exact_select(hits, gene_cache, gene_tiebreak=not args.no_gene_tiebreak)
    if not args.selection and not args.limit:  # only persist the authoritative full audit; preserve any estimated one
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
        sel_stamp = BSG.file_stamp(selection_out)      # the audit as written
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
    # ---- 3a. --pfunc e3: one PGFam store for the whole run, built before the pool forks ----
    pfunc, pfunc_meta = BSG.prepare_pfunc(
        dict(tasks), out_dir, mode=args.pfunc, min_n=args.pfunc_min_n, dedup=args.pfunc_dedup,
        pgfam_cache=os.path.abspath(args.pfunc_pgfam_cache),
        catalogue=os.path.abspath(args.pfunc_catalogue), fetch=not args.pfunc_no_fetch,
        fetch_workers=args.pfunc_fetch_workers,
        allow_stale=args.pfunc_allow_stale_catalogue, preflight=preflight,
        code_files=BSG.PFUNC_PROVENANCE_FILES + [os.path.abspath(__file__)])
    print(f"[build] {len(tasks)} ASVs to build | {len(abstains)} abstain | "
          f"procs={args.procs} | pfunc={args.pfunc} | out={out_dir}", flush=True)

    manifest, built, failed = {}, 0, 0
    run_digests, dig_conflicts = {}, set()
    t0 = time.perf_counter()
    pool_kw = {}
    if pfunc.e3:
        # fork: the workers inherit the (read-only) PGFam store and catalogue instead of
        # unpickling a copy each
        pool_kw["mp_context"] = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(max_workers=args.procs, initializer=_winit,
                             initargs=(out_dir, cache_dir, pfunc if pfunc.e3 else None),
                             **pool_kw) as ex:
        for i, (asv, entry, dg) in enumerate(ex.map(_wbuild, tasks, chunksize=4), 1):
            manifest[asv] = entry
            built += entry["status"] == "built"
            failed += entry["status"] == "failed"
            if entry["status"] == "built":
                BSG.merge_source_digests(run_digests, dg, dig_conflicts)
            if i % 100 == 0 or i == len(tasks):
                print(f"[build] {i}/{len(tasks)} built={built} failed={failed} "
                      f"({time.perf_counter()-t0:.0f}s, {i/max(time.perf_counter()-t0,1):.1f}/s)",
                      flush=True)

    manifest.update(abstains)
    if pfunc.e3:
        pfunc_meta.update(BSG.pfunc_counts(manifest.values()))
    dig_meta = BSG.write_source_digests(out_dir, run_digests, dig_conflicts)
    man_path = os.path.join(out_dir, "synthetic_genome_manifest.json")
    with open(man_path, "w") as f:
        json.dump({"_meta": {"selection": (selection_out if (args.selection or not args.limit)
                                           else "(limit run)"),
                             "selection_stamp": sel_stamp,
                             "gene_data": gene_data, "annotation": BSG.ANNOTATION_STAMP,
                             "builder": "fast",
                             "n_asvs": len(audit), "built": built,
                             "abstain": len(abstains), "failed": failed,
                             "source_digests": dig_meta, "pfunc": pfunc_meta},
                   **manifest}, f, indent=1)
    dt = time.perf_counter() - t0
    print(f"[build] DONE: built={built} abstain={len(abstains)} failed={failed} in {dt:.0f}s "
          f"-> {man_path}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
