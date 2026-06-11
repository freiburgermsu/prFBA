#!/usr/bin/env python
"""
align_hits.py — alignment-scoring strategy for hitting 16S ASVs against the BV-BRC 16S
reference set, to pick the representative reference genome(s) per ASV.

STRATEGY (hardware-independent prefilter; hardware-accelerated scoring):
  Stage 1 — edlib (infix / "HW") edit-distance prefilter over EVERY unique reference,
            keeping the top-500 lowest-distance candidates per ASV. ALWAYS, on every
            machine. This is cheap and discards the ~99.9% of obviously-irrelevant
            references before any expensive Smith-Waterman scoring.
  Stage 2 — Smith-Waterman SCORING of just those top-500 candidates (match +2 / mismatch
            -3 / gap_open -5 / gap_extend -2, local), on:
              * GPU  if CUDA (CuPy SW kernel from gpu_align.py) or Metal is present, else
              * CPU  (Biopython PairwiseAligner).
  Stage 3 — Biopython re-alignment of the surviving top candidates for exact % identity,
            aligned length and reference coordinates; emit the top-k2 (default 20).

Why top-500 and not exhaustive: a deep prefilter captures the best representative(s) for
each ASV — including the large equal-Smith-Waterman-score tie clusters that decide which
near-identical reference genome is named the representative (empirically up to ~200 genomes
tie at the top score) — WITHOUT paying to score the whole 459k-reference database. Keeping
500 leaves wide margin above any observed tie cluster.

    python align_hits.py --explain          # print the hardware probe + chosen backend, run nothing
    python align_hits.py                     # auto-detect backend and run on all ASVs
    python align_hits.py --fasta X.fasta --taxonomy-csv X.csv --outdir DIR
    python align_hits.py --force-backend cpu # override the scoring backend (cuda|metal|cpu)
    python align_hits.py --prefilter-k 500   # edlib shortlist depth (default 500)

Outputs (in --outdir): asv_top20_alignment_hits.json + asv_alignment_summary.csv (same
schema as before) + method_selection.json (probe, decision, run stats).

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations
import argparse, json, os, shutil, time
from dataclasses import asdict, dataclass
from multiprocessing import Pool

import numpy as np

import edlib_biopython_hits as P   # loaders, scoring constants, stage cores, enrichment


# ============================================================================ hardware
@dataclass
class Hardware:
    cuda: bool
    metal: bool
    n_cores: int
    free_disk_gb: float


@dataclass
class Decision:
    prefilter_k: int      # edlib shortlist depth passed to the scorer (always set)
    backend: str          # "cuda" | "metal" | "cpu" — the SW-scoring engine for the shortlist
    reason: str


def detect_hardware(path: str = ".") -> Hardware:
    """Probe CUDA (CuPy), Apple Metal/MPS, logical cores, and free disk."""
    cuda = False
    try:
        import cupy as cp
        cuda = cp.cuda.runtime.getDeviceCount() > 0
    except Exception:
        cuda = False
    metal = False
    try:
        import platform
        if platform.system() == "Darwin":
            try:
                import torch
                metal = bool(torch.backends.mps.is_available())
            except Exception:
                metal = False
    except Exception:
        metal = False
    return Hardware(cuda, metal, os.cpu_count() or 1, shutil.disk_usage(path).free / 1e9)


def choose_backend(hw: Hardware, force: str = "auto") -> Decision:
    """Pure mapping from hardware to the shortlist-scoring backend. Prefilter is ALWAYS
    edlib top-`prefilter_k`; only the Stage-2 scoring engine depends on hardware."""
    if force != "auto":
        return Decision(0, force, f"forced backend: {force}")
    if hw.cuda:
        return Decision(0, "cuda", "CUDA device present -> GPU Smith-Waterman scoring of the top-500 shortlist")
    if hw.metal:
        return Decision(0, "metal", "Apple Metal/MPS present -> GPU scoring of the top-500 shortlist")
    return Decision(0, "cpu",
                    f"no GPU (CUDA/Metal) -> Biopython CPU Smith-Waterman scoring of the top-500 shortlist "
                    f"({hw.n_cores} cores)")


# ============================================================================ workers
def _edlib_only_worker(task):
    """Stage 1: edlib HW edit-distance prefilter for one ASV -> top-k1 (dist, ref_idx)."""
    i, seq, k1 = task
    cand, _t = P._scan_asv(seq, k1)        # uses fork-inherited P.REF_BLOB/REF_OFF/REF_LEN
    return i, cand


def _rescore_worker(task):
    """Stage 3: Biopython rescore a candidate shortlist -> top-k2 record dicts (with coords)."""
    i, seq, cand, k2 = task
    return i, P._rescore(seq, cand, k2)     # ALIGNER set by P.init_worker


# ============================================================================ CPU scoring path
def run_cpu(outdir, asv_ids, asv_seqs, md5s, *, k1, k2, workers, chunksize=2):
    """edlib top-k1 prefilter -> Biopython SW scoring of the shortlist -> top-k2 (one fused pass)."""
    n_asv, n_ref = len(asv_ids), len(md5s)
    tasks = [(i, asv_seqs[i], k1, k2) for i in range(n_asv)]
    print(f"[cpu] edlib top-{k1} prefilter + Biopython SW scoring | {n_asv} ASVs x {n_ref:,} refs | "
          f"workers={workers}", flush=True)
    results, ce, cb = {}, 0.0, 0.0
    t0 = time.perf_counter()
    with Pool(workers, initializer=P.init_worker) as pool:
        done = 0
        for i, recs, te, tb in pool.imap_unordered(P.process_asv, tasks, chunksize=chunksize):
            results[i] = recs; ce += te; cb += tb; done += 1
            if done % 200 == 0 or done == n_asv:
                el = time.perf_counter() - t0
                print(f"[cpu] {done}/{n_asv} elapsed={el:.0f}s eta={el/done*(n_asv-done):.0f}s", flush=True)
    P.enrich_and_write(outdir, asv_ids, asv_seqs, md5s, results, k2=k2, n_cand_per_asv=min(k1, n_ref),
                       search_label=f"edlib_top{k1}_prefilter+biopython_SW")
    return {"backend": "cpu", "prefilter_k": k1, "k2": k2,
            "cpu_seconds_edlib": round(ce, 1), "cpu_seconds_biopython": round(cb, 1),
            "wall_seconds": round(time.perf_counter() - t0, 1)}


# ============================================================================ GPU scoring path
def run_gpu(outdir, asv_ids, asv_seqs, md5s, *, k1, k2, workers, gpu_cand=40,
            pairs_batch=4_000_000, chunksize=2):
    """edlib top-k1 prefilter (CPU, parallel) -> GPU Smith-Waterman scoring of ONLY those
    k1 candidates per ASV (CuPy kernel) -> Biopython rescore the GPU-top candidates for
    identity/coords -> top-k2. Scores the shortlist, never the whole reference set."""
    import cupy as cp
    import gpu_align as G                  # CuPy SW kernel (CUDA); imported lazily
    n_asv, n_ref = len(asv_ids), len(md5s)

    # ---- Stage 1: edlib prefilter -> per-ASV top-k1 (edit_dist, ref_idx) ----
    print(f"[gpu] stage1 edlib top-{k1} prefilter | {n_asv} ASVs x {n_ref:,} refs | workers={workers}",
          flush=True)
    cands = {}
    t0 = time.perf_counter()
    with Pool(workers) as pool:            # edlib-only: P globals inherited via fork, no ALIGNER needed
        done = 0
        for i, cand in pool.imap_unordered(_edlib_only_worker,
                                           [(i, asv_seqs[i], k1) for i in range(n_asv)], chunksize=chunksize):
            cands[i] = cand; done += 1
            if done % 200 == 0 or done == n_asv:
                el = time.perf_counter() - t0
                print(f"[gpu] prefilter {done}/{n_asv} elapsed={el:.0f}s eta={el/done*(n_asv-done):.0f}s", flush=True)
    edlib_wall = time.perf_counter() - t0

    # ---- Stage 2: GPU-score the flattened (ASV, candidate-ref) pairs ----
    ref_codes, ref_off = G.encode_blob(P.REF_BLOB, P.REF_OFF)
    gpu = G.GpuSW(asv_seqs, ref_codes, ref_off)
    qi, rj, spans = [], [], []
    for i in range(n_asv):
        s = len(qi)
        for _dist, ref_idx in cands[i]:
            qi.append(i); rj.append(ref_idx)
        spans.append((s, len(qi)))
    qi = np.asarray(qi, np.int32); rj = np.asarray(rj, np.int32)
    print(f"[gpu] stage2 GPU SW scoring of {len(qi):,} (ASV,candidate) pairs "
          f"(~{len(qi)/max(n_asv,1):.0f}/ASV) — NOT the {n_ref:,}-ref exhaustive product", flush=True)
    t1 = time.perf_counter()
    scores = np.empty(len(qi), np.int32)
    for b in range(0, len(qi), pairs_batch):
        sc = gpu.score_pairs(qi[b:b + pairs_batch], rj[b:b + pairs_batch])
        scores[b:b + sc.shape[0]] = cp.asnumpy(sc)
    gpu_wall = time.perf_counter() - t1
    del gpu
    try:
        cp.get_default_memory_pool().free_all_blocks()
    except Exception:
        pass

    # ---- Stage 3: per-ASV GPU-top-`gpu_cand` -> Biopython rescore for identity/coords ----
    keep = max(gpu_cand, k2)
    rescore_tasks = []
    for i in range(n_asv):
        s, e = spans[i]
        if s == e:
            rescore_tasks.append((i, asv_seqs[i], [], k2)); continue
        order = np.argsort(-scores[s:e])[:keep]      # best GPU scores first
        c = cands[i]
        rescore_tasks.append((i, asv_seqs[i], [(c[j][0], c[j][1]) for j in order], k2))
    print(f"[gpu] stage3 Biopython rescore top-{keep} per ASV (identity/coords) | workers={workers}", flush=True)
    results = {}
    t2 = time.perf_counter()
    with Pool(workers, initializer=P.init_worker) as pool:
        done = 0
        for i, recs in pool.imap_unordered(_rescore_worker, rescore_tasks, chunksize=4):
            results[i] = recs; done += 1
            if done % 500 == 0 or done == n_asv:
                print(f"[gpu] rescore {done}/{n_asv} ({time.perf_counter()-t2:.0f}s)", flush=True)
    P.enrich_and_write(outdir, asv_ids, asv_seqs, md5s, results, k2=k2, n_cand_per_asv=min(k1, n_ref),
                       search_label=f"edlib_top{k1}_prefilter+gpu_SW_score+biopython_coords")
    return {"backend": "cuda", "prefilter_k": k1, "k2": k2, "gpu_cand": keep,
            "wall_seconds_edlib_prefilter": round(edlib_wall, 1),
            "wall_seconds_gpu_score": round(gpu_wall, 1),
            "wall_seconds_biopython_rescore": round(time.perf_counter() - t2, 1),
            "n_scored_pairs": int(len(qi))}


# ================================================================================ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=os.path.join(P.EMILYKIN, "bvbrc_alignment_hits"))
    ap.add_argument("--fasta", default=P.ASVS_FASTA, help="ASV FASTA to align (default = EmilyKin asvs.fasta)")
    ap.add_argument("--taxonomy-csv", default=None,
                    help="per-ASV MiDAS lineage + rel_ab CSV (cols seq,Kingdom..Species,rel_ab)")
    ap.add_argument("--prefilter-k", type=int, default=500,
                    help="edlib edit-distance shortlist depth kept per ASV before SW scoring (default 500)")
    ap.add_argument("--k2", type=int, default=20, help="final hits kept per ASV")
    ap.add_argument("--min-ref-len", type=int, default=200)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 1) - 4))
    ap.add_argument("--gpu-cand", type=int, default=100,
                    help="GPU path: top GPU-scored candidates Biopython-rescored per ASV for final "
                         "identity/coords (margin above k2 to absorb equal-score tie clusters)")
    ap.add_argument("--force-backend", default="auto", choices=["auto", "cuda", "metal", "cpu"],
                    help="override the SW-scoring backend (prefilter is always edlib top-k)")
    ap.add_argument("--limit", type=int, default=0, help="run on first N ASVs (0=all)")
    ap.add_argument("--explain", action="store_true", help="probe hardware, print the decision, exit")
    args = ap.parse_args()

    hw = detect_hardware(args.outdir if os.path.isdir(args.outdir) else ".")
    decision = choose_backend(hw, args.force_backend)
    decision.prefilter_k = args.prefilter_k
    # Metal currently has no SW kernel here (the GPU kernel is CuPy/CUDA-only) -> score on CPU.
    effective = decision.backend
    note = ""
    if decision.backend == "metal":
        note = ("Metal/MPS detected, but the GPU SW kernel is CuPy/CUDA-only and no Metal kernel is "
                "implemented yet -> scoring the top-500 shortlist on CPU (Biopython). Detection + "
                "strategy are in place for a future Metal backend.")
        effective = "cpu"
    if decision.backend == "cuda" and not hw.cuda:
        raise SystemExit("[abort] backend 'cuda' requested but no CUDA device is available")

    print("=" * 80)
    print(f"[hardware] cuda={hw.cuda} metal={hw.metal} cores={hw.n_cores} free_disk={hw.free_disk_gb:.0f}GB")
    print(f"[strategy] ALWAYS edlib top-{decision.prefilter_k} prefilter -> "
          f"{effective.upper()} Smith-Waterman scoring of the shortlist -> top-{args.k2}")
    print(f"[why]      {decision.reason}")
    if note:
        print(f"[note]     {note}")
    print("=" * 80, flush=True)

    if args.taxonomy_csv:
        P.TAXONOMY_CSV = args.taxonomy_csv
    os.makedirs(args.outdir, exist_ok=True)
    sel = {"hardware": asdict(hw), "decision": asdict(decision), "effective_backend": effective,
           "note": note, "params": {"prefilter_k": args.prefilter_k, "k2": args.k2,
                                     "min_ref_len": args.min_ref_len, "workers": args.workers,
                                     "gpu_cand": args.gpu_cand, "limit": args.limit, "fasta": args.fasta}}
    with open(os.path.join(args.outdir, "method_selection.json"), "w") as fh:
        json.dump(sel, fh, indent=2)
    if args.explain:
        print("[explain] wrote method_selection.json; no alignment run.", flush=True)
        return

    wall0 = time.perf_counter()
    print(f"[load] references from {os.path.basename(P.DB_MD5_SEQ)} ...", flush=True)
    md5s, P.REF_BLOB, P.REF_OFF, P.REF_LEN = P.load_references(args.min_ref_len)
    print(f"[load] {len(md5s):,} unique refs (>= {args.min_ref_len} bp)", flush=True)
    asv_ids, asv_seqs = P.load_asvs(args.fasta)
    if args.limit:
        asv_ids, asv_seqs = asv_ids[:args.limit], asv_seqs[:args.limit]
    print(f"[load] {len(asv_ids):,} ASVs", flush=True)

    if effective == "cuda":
        run_stats = run_gpu(args.outdir, asv_ids, asv_seqs, md5s, k1=args.prefilter_k, k2=args.k2,
                            workers=args.workers, gpu_cand=args.gpu_cand)
    else:
        run_stats = run_cpu(args.outdir, asv_ids, asv_seqs, md5s, k1=args.prefilter_k, k2=args.k2,
                            workers=args.workers)

    sel["run_stats"] = run_stats
    sel["wall_seconds_total"] = round(time.perf_counter() - wall0, 1)
    sel["n_asvs"] = len(asv_ids); sel["n_unique_refs"] = len(md5s)
    with open(os.path.join(args.outdir, "method_selection.json"), "w") as fh:
        json.dump(sel, fh, indent=2)
    print(f"[done] backend={effective} prefilter=top-{args.prefilter_k} "
          f"total wall={sel['wall_seconds_total']:.0f}s -> {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
