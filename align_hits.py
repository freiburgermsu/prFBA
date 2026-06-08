#!/usr/bin/env python
"""
align_hits.py — hardware-aware dispatcher that picks HOW to hit experimental
sequences (16S ASVs) against the BV-BRC reference database, then runs that method.

The alignment scheme (local Smith-Waterman, match +2 / mismatch -3 / gap_open -5 /
gap_extend -2) is identical across every method, so results are directly comparable;
only the *engine* and the *prefilter depth* change with the available compute.

Decision tree (see `select_method`)
-----------------------------------
1. CUDA available .......... GPU Smith-Waterman via the standalone `gpusw` package
                             (../gpuSW): exhaustive all-reference scoring on the GPU,
                             then a Biopython rescoring of the per-ASV top candidates
                             for % identity / aligned coords.
2. CPU < 32 cores .......... edlib (infix / "HW") edit-distance prefilter, then
                             Biopython local SW on the shortlist. The shortlist depth
                             scales with cores (and disk, for the cache):
                               * 16-32 cores AND >= 10 GB free -> top 5000, cached to disk
                               * 16-32 cores AND  < 10 GB free -> top 2000 (no cache)
                               *  8-16 cores                   -> top 2000
                               *  < 8 cores                    -> top 1000
3. CPU >= 32 cores (no GPU)  full exhaustive Biopython local SW against every reference
                             (no prefilter) — slow but prefilter-free; only chosen when
                             there are plenty of cores and no GPU.

All three paths funnel into `edlib_biopython_hits.enrich_and_write`, so every method
emits the same per-ASV hits JSON + summary CSV, plus a `method_selection.json`
recording the hardware probe and the chosen method.

    python align_hits.py                 # auto-detect and run on all ASVs
    python align_hits.py --explain       # just print the probe + decision, run nothing
    python align_hits.py --limit 5       # smoke test on the first 5 ASVs
    python align_hits.py --force-method biopython_full

Interpreter: ~/Documents/py_venv/bin/python   (per project convention)
"""
from __future__ import annotations
import argparse, json, os, shutil, time
from dataclasses import asdict, dataclass
from multiprocessing import Pool

import numpy as np

import edlib_biopython_hits as P   # loaders, scoring constants, stage cores, enrichment

GPUSW_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "gpuSW")


# ============================================================================ hardware
@dataclass
class Hardware:
    cuda_available: bool
    n_cores: int
    free_disk_gb: float


@dataclass
class Decision:
    method: str           # "gpu" | "edlib_biopython" | "biopython_full"
    k1: int | None        # edlib shortlist depth passed to Biopython (edlib_biopython only)
    save_edlib: bool      # cache the edlib shortlist to disk before rescoring
    reason: str


def detect_hardware(path: str = ".") -> Hardware:
    """Probe CUDA (via gpusw), logical CPU cores, and free disk at `path`."""
    cuda = False
    try:
        import gpusw
        cuda = bool(gpusw.gpu_available())
    except Exception:
        cuda = False
    n_cores = os.cpu_count() or 1
    free_gb = shutil.disk_usage(path).free / 1e9
    return Hardware(cuda, n_cores, free_gb)


def _edlib_tier(n_cores: int, free_disk_gb: float, min_disk_gb: float):
    """Shortlist depth + on-disk-cache flag + reason for the CPU(<32)/edlib path."""
    if n_cores >= 16:   # 16-32 cores
        if free_disk_gb >= min_disk_gb:
            return 5000, True, (
                f"{n_cores} cores (16-32) and {free_disk_gb:.0f} GB free "
                f"(>= {min_disk_gb:.0f}) -> edlib prefilter, cache + pass top 5000 "
                f"hits/ASV to Biopython")
        return 2000, False, (
            f"{n_cores} cores (16-32) but only {free_disk_gb:.0f} GB free "
            f"(< {min_disk_gb:.0f}) -> edlib prefilter, pass top 2000 hits/ASV to "
            f"Biopython (shortlist not cached)")
    if n_cores >= 8:    # 8-16 cores
        return 2000, False, (
            f"{n_cores} cores (8-16) -> edlib prefilter, pass top 2000 hits/ASV to "
            f"Biopython")
    return 1000, False, (  # < 8 cores
        f"{n_cores} cores (<8) -> edlib prefilter, pass top 1000 hits/ASV to Biopython")


def select_method(cuda_available: bool, n_cores: int, free_disk_gb: float,
                  *, min_disk_gb: float = 10.0) -> Decision:
    """Pure mapping from a hardware probe to the alignment method + parameters."""
    if cuda_available:
        return Decision("gpu", None, False,
                        "CUDA device available -> GPU Smith-Waterman via the gpusw package")
    if n_cores < 32:
        k1, save, why = _edlib_tier(n_cores, free_disk_gb, min_disk_gb)
        return Decision("edlib_biopython", k1, save, why)
    return Decision("biopython_full", None, False,
                    f"{n_cores} cores (>= 32) and no CUDA -> full exhaustive Biopython "
                    f"local SW against every reference (no prefilter)")


def decide(hw: Hardware, *, force_method: str = "auto",
           min_disk_gb: float = 10.0) -> Decision:
    """Apply `select_method`, or honour an explicit `--force-method` override."""
    if force_method == "auto":
        return select_method(hw.cuda_available, hw.n_cores, hw.free_disk_gb,
                             min_disk_gb=min_disk_gb)
    if force_method == "gpu":
        return Decision("gpu", None, False, "forced: gpu")
    if force_method == "biopython_full":
        return Decision("biopython_full", None, False, "forced: biopython_full")
    if force_method == "edlib_biopython":
        k1, save, why = _edlib_tier(hw.n_cores, hw.free_disk_gb, min_disk_gb)
        return Decision("edlib_biopython", k1, save, "forced: " + why)
    raise ValueError(f"unknown force_method {force_method!r}")


# ============================================================================ backends
def _run_pool(worker, tasks, workers, chunksize, n_asv, label):
    """Run `worker` over `tasks` in a fork Pool (P.init_worker seeds each ALIGNER).

    Each worker yields (asv_idx, records, *rest); returns {asv_idx: records} and any
    extra per-task payload as {asv_idx: rest}. Prints progress + a parallel speedup.
    """
    results, extra = {}, {}
    cpu_edlib = cpu_bio = 0.0
    t0 = time.perf_counter()
    with Pool(workers, initializer=P.init_worker) as pool:
        done = 0
        for asv_idx, recs, *rest in pool.imap_unordered(worker, tasks, chunksize=chunksize):
            results[asv_idx] = recs
            # the two timing floats are always the last two elements of `rest`
            te, tb = rest[-2], rest[-1]
            cpu_edlib += te; cpu_bio += tb
            if len(rest) > 2:
                extra[asv_idx] = rest[0]
            done += 1
            if done % 100 == 0 or done == n_asv:
                el = time.perf_counter() - t0
                eta = el / done * (n_asv - done)
                print(f"[{label}] {done}/{n_asv}  elapsed={el:6.0f}s  eta={eta:6.0f}s  "
                      f"(cpu edlib={cpu_edlib:,.0f}s bio={cpu_bio:,.0f}s)", flush=True)
    return results, extra, cpu_edlib, cpu_bio, time.perf_counter() - t0


def _save_edlib_cache(outdir, asv_ids, extra, k1, min_ref_len):
    """Persist the per-ASV edlib shortlist (ref_idx + distance, padded to k1) to .npz."""
    n = len(asv_ids)
    idx = np.full((n, k1), -1, np.int32)
    dist = np.full((n, k1), -1, np.int32)
    for i in range(n):
        cand = extra.get(i, [])          # list of (dist, ref_idx) ascending
        for j, (d, ri) in enumerate(cand[:k1]):
            dist[i, j] = d; idx[i, j] = ri
    path = os.path.join(outdir, "edlib_shortlist.npz")
    np.savez_compressed(path, ref_idx=idx, edlib_distance=dist,
                        asv_ids=np.array(asv_ids), k1=np.int32(k1),
                        min_ref_len=np.int32(min_ref_len))
    print(f"[cache] wrote {path} ({os.path.getsize(path)/1e6:.0f} MB) — "
          f"ref_idx indexes load_references(min_ref_len={min_ref_len}) order", flush=True)


def run_edlib_biopython(outdir, asv_ids, asv_seqs, md5s, *, k1, k2, workers,
                        save_edlib, min_ref_len, chunksize=2):
    n_asv, n_ref = len(asv_ids), len(md5s)
    worker = P.process_asv_savecand if save_edlib else P.process_asv
    tasks = [(i, asv_seqs[i], k1, k2) for i in range(n_asv)]
    print(f"[edlib_biopython] {n_asv} ASVs x {n_ref:,} refs | k1={k1} k2={k2} | "
          f"workers={workers} | cache={'on' if save_edlib else 'off'}", flush=True)
    results, extra, ce, cb, wall = _run_pool(worker, tasks, workers, chunksize, n_asv,
                                             "edlib_biopython")
    if save_edlib:
        _save_edlib_cache(outdir, asv_ids, extra, k1, min_ref_len)
    P.enrich_and_write(outdir, asv_ids, asv_seqs, md5s, results,
                       k2=k2, n_cand_per_asv=min(k1, n_ref),
                       search_label="edlib_HW_prefilter+biopython_local_SW")
    return {"cpu_seconds_edlib": round(ce, 1), "cpu_seconds_biopython": round(cb, 1),
            "wall_seconds_scan": round(wall, 1), "k1": k1, "k2": k2}


def run_biopython_full(outdir, asv_ids, asv_seqs, md5s, *, k2, workers, chunksize=1):
    n_asv, n_ref = len(asv_ids), len(md5s)
    tasks = [(i, asv_seqs[i], 0, k2) for i in range(n_asv)]
    print(f"[biopython_full] EXHAUSTIVE {n_asv} ASVs x {n_ref:,} refs | k2={k2} | "
          f"workers={workers}  (no prefilter — this is the slow path)", flush=True)
    results, _extra, ce, cb, wall = _run_pool(P.process_asv_biopython_full, tasks,
                                              workers, chunksize, n_asv, "biopython_full")
    P.enrich_and_write(outdir, asv_ids, asv_seqs, md5s, results, k2=k2,
                       search_label="biopython_full_exhaustive_local_SW")
    return {"cpu_seconds_biopython": round(cb, 1), "wall_seconds_scan": round(wall, 1),
            "k2": k2}


# --- GPU enrichment worker globals (set before the Pool; inherited via fork) ---
_G_ASV = None      # asv_seqs
_G_CAND = None     # (n_asv, cand) int32 reference indices from gpusw top_k
_G_K2 = None


def _gpu_enrich(qi):
    """Biopython-rescore one ASV's GPU candidate shortlist -> top-k2 record dicts."""
    cand = [(None, int(j)) for j in _G_CAND[qi] if j >= 0]
    return qi, P._rescore(_G_ASV[qi], cand, _G_K2)


def run_gpu(outdir, asv_ids, asv_seqs, md5s, *, k2, workers, cand, query_batch,
            chunksize=4):
    """GPU exhaustive SW (gpusw package) -> per-ASV top-`cand`, then Biopython rescoring."""
    import gpusw
    global _G_ASV, _G_CAND, _G_K2
    n_asv, n_ref = len(asv_ids), len(md5s)
    cand = max(cand, k2)
    print(f"[gpu] gpusw {gpusw.__version__} from {os.path.dirname(gpusw.__file__)} | "
          f"{n_asv} ASVs x {n_ref:,} refs | top-{cand} GPU candidates -> top-{k2}",
          flush=True)

    # decode the shared reference blob into strings for the gpusw funnel (one-time)
    refs = [P.REF_BLOB[P.REF_OFF[i]:P.REF_OFF[i + 1]].decode("ascii") for i in range(n_ref)]
    t0 = time.perf_counter()
    al = gpusw.Aligner(gpusw.schemes.DNA).index(refs)           # bit-exact prFBA scheme
    res = al.top_k(asv_seqs, k=cand, query_batch=query_batch)   # AlignResult w/ topk_idx
    gpu_wall = time.perf_counter() - t0
    gcups = al.gcups
    print(f"[gpu] exhaustive top-{cand} done in {gpu_wall:.0f}s "
          f"({gcups:.0f} GCUPS last launch)", flush=True)

    # free the device before forking the (CPU-only) enrichment pool
    cand_idx = np.ascontiguousarray(res.topk_idx, dtype=np.int32)
    del res, al, refs
    try:
        import cupy as cp
        cp.get_default_memory_pool().free_all_blocks()
    except Exception:
        pass

    _G_ASV, _G_CAND, _G_K2 = asv_seqs, cand_idx, k2
    print(f"[gpu] enrich {n_asv} ASVs (Biopython identity/coords) across {workers} workers ...",
          flush=True)
    results = {}
    te = time.perf_counter()
    with Pool(workers, initializer=P.init_worker) as pool:
        done = 0
        for qi, recs in pool.imap_unordered(_gpu_enrich, range(n_asv), chunksize=chunksize):
            results[qi] = recs; done += 1
            if done % 500 == 0 or done == n_asv:
                print(f"[gpu] enriched {done}/{n_asv} ({time.perf_counter()-te:.0f}s)", flush=True)
    enrich_wall = time.perf_counter() - te

    P.enrich_and_write(outdir, asv_ids, asv_seqs, md5s, results, k2=k2,
                       n_cand_per_asv=cand,
                       search_label="gpu_exhaustive_local_SW_gpusw+biopython_rescore")
    return {"gpu_exhaustive_wall_s": round(gpu_wall, 1), "gpu_gcups": round(gcups, 1),
            "enrich_wall_s": round(enrich_wall, 1), "gpu_candidates": cand, "k2": k2}


# ================================================================================ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=os.path.join(P.EMILYKIN, "bvbrc_alignment_hits"))
    ap.add_argument("--force-method", default="auto",
                    choices=["auto", "gpu", "edlib_biopython", "biopython_full"],
                    help="override the hardware-based selection")
    ap.add_argument("--k2", type=int, default=20, help="final hits kept per ASV")
    ap.add_argument("--min-ref-len", type=int, default=200)
    ap.add_argument("--min-disk-gb", type=float, default=10.0,
                    help="free-disk threshold gating the cached top-5000 edlib tier")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 1) - 4))
    ap.add_argument("--limit", type=int, default=0, help="run on first N ASVs (0=all)")
    ap.add_argument("--gpu-cand", type=int, default=40,
                    help="GPU candidates pulled per ASV before Biopython rescoring")
    ap.add_argument("--query-batch", type=int, default=256,
                    help="GPU: ASVs scored per top_k launch (device-memory bound)")
    ap.add_argument("--explain", action="store_true",
                    help="probe hardware, print the decision, and exit without running")
    args = ap.parse_args()

    hw = detect_hardware(args.outdir if os.path.isdir(args.outdir) else ".")
    decision = decide(hw, force_method=args.force_method, min_disk_gb=args.min_disk_gb)
    print("=" * 78)
    print(f"[hardware] cuda={hw.cuda_available}  cores={hw.n_cores}  "
          f"free_disk={hw.free_disk_gb:.0f} GB")
    print(f"[method]   {decision.method}"
          + (f"  (k1={decision.k1}, cache={decision.save_edlib})"
             if decision.method == "edlib_biopython" else ""))
    print(f"[why]      {decision.reason}")
    print("=" * 78, flush=True)

    if decision.method == "gpu" and not hw.cuda_available:
        raise SystemExit("[abort] method 'gpu' requested but no CUDA device is available")

    os.makedirs(args.outdir, exist_ok=True)
    sel = {"hardware": asdict(hw), "decision": asdict(decision),
           "params": {"k2": args.k2, "min_ref_len": args.min_ref_len,
                      "min_disk_gb": args.min_disk_gb, "workers": args.workers,
                      "limit": args.limit, "gpu_cand": args.gpu_cand}}
    with open(os.path.join(args.outdir, "method_selection.json"), "w") as fh:
        json.dump(sel, fh, indent=2)

    if args.explain:
        print("[explain] wrote method_selection.json; no alignment run.", flush=True)
        return

    # ---- load references (into the shared P globals) + ASVs ----
    wall0 = time.perf_counter()
    print(f"[load] references from {os.path.basename(P.DB_MD5_SEQ)} ...", flush=True)
    md5s, P.REF_BLOB, P.REF_OFF, P.REF_LEN = P.load_references(args.min_ref_len)
    n_ref = len(md5s)
    print(f"[load] {n_ref:,} unique refs (>= {args.min_ref_len} bp), "
          f"blob {len(P.REF_BLOB)/1e6:.0f} MB", flush=True)
    asv_ids, asv_seqs = P.load_asvs(P.ASVS_FASTA)
    if args.limit:
        asv_ids, asv_seqs = asv_ids[:args.limit], asv_seqs[:args.limit]
    print(f"[load] {len(asv_ids):,} ASVs", flush=True)

    if decision.method == "gpu":
        run_stats = run_gpu(args.outdir, asv_ids, asv_seqs, md5s, k2=args.k2,
                            workers=args.workers, cand=args.gpu_cand,
                            query_batch=args.query_batch)
    elif decision.method == "edlib_biopython":
        run_stats = run_edlib_biopython(args.outdir, asv_ids, asv_seqs, md5s,
                                        k1=decision.k1, k2=args.k2, workers=args.workers,
                                        save_edlib=decision.save_edlib,
                                        min_ref_len=args.min_ref_len)
    else:  # biopython_full
        run_stats = run_biopython_full(args.outdir, asv_ids, asv_seqs, md5s,
                                       k2=args.k2, workers=args.workers)

    sel["run_stats"] = run_stats
    sel["wall_seconds_total"] = round(time.perf_counter() - wall0, 1)
    sel["n_asvs"] = len(asv_ids)
    sel["n_unique_refs"] = n_ref
    with open(os.path.join(args.outdir, "method_selection.json"), "w") as fh:
        json.dump(sel, fh, indent=2)
    print(f"[done] method={decision.method}  total wall={sel['wall_seconds_total']:.0f}s  "
          f"-> {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
