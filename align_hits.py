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
  Stage 3 — Biopython re-alignment of the best --gpu-cand GPU-scored candidates (default
            200) for exact % identity, aligned length and reference coordinates; emit the
            top-k2 (default 100).

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
    python align_hits.py --ref-db REF.json --ref-id-db REF_ID.json   # restricted reference DB

RESTRICTED REFERENCE DB (--ref-db / --ref-id-db)
  By default the reference set is the whole BV-BRC 16S database, hard-coded in
  edlib_biopython_hits as DB_MD5_SEQ / DB_MD5_ID. `--ref-db` (md5 -> sequence JSON) and
  `--ref-id-db` (md5 -> BV-BRC header JSON) replace those two, so a caller can align
  against an arbitrary SUBSET of the references. This is what makes cross-validation
  possible: in a k-fold benchmark the held-out genomes and their near-identical
  neighbours must be absent from the reference *at alignment time*. Deleting their hits
  afterwards is not equivalent — the record has already been truncated to the top `--k2`,
  so excising near neighbours empties it rather than demoting them (the measured
  k2-truncation artifact behind 94.6% of self-exclusion abstentions). The two flags
  rebind `P.DB_MD5_SEQ` / `P.DB_MD5_ID` (both are read at call time) and are additionally
  passed explicitly to `P.load_references` / `P.enrich_and_write`, so the override cannot
  depend on import order. Output schema is unchanged; `method_selection.json["params"]`
  records `ref_db` / `ref_id_db` (null when unused) so every run says what it aligned to.
  Pair them with `--taxonomy-csv` pointing at a header-only MiDAS CSV when no external
  taxonomic prior may enter the records.

GPU PRUNING AND THE GUARD: --gpu-cand >= --k2 (exit 2 otherwise)
  The CUDA path Biopython-rescores only the best `max(--gpu-cand, --k2)` GPU-scored
  candidates, picked under the rescore's own total order (score desc, edlib distance asc,
  reference index asc). The GPU score is bit-exact to Biopython's, so the pick contains the
  top-k2 of the whole shortlist under that order, and the record is the one a rescore of
  all --prefilter-k candidates gives, at any --gpu-cand >= --k2: the GPU prunes without
  changing the record, and a record still depends only on its sequence. (An unstable
  argsort pick, used before 2026-10-27, cut large equal-score tie clusters arbitrarily:
  it would have changed 2,180 of the 60,461 depth-100 records at --gpu-cand 200; see
  region_validation_d100/GPU_CAND_200.md.) The margin above --k2 guards against any
  GPU/Biopython score disagreement; with `gpu_cand < k2` there would be none, so the
  mistake aborts.

SHORTLIST CACHE (--save-shortlist / --shortlist-in DIR; CUDA path only; both off by default)
  The edlib prefilter is 70-90% of wall time and depends only on (ASV sequences, reference
  set, --prefilter-k, --min-ref-len), never on --k2 or --gpu-cand. `--save-shortlist`
  persists it so a later run at another depth costs only the GPU + Biopython stages:
    edlib_top{k1}_shortlist.npz  idx/dist int32[n_asv, k1] (rows padded with -1; `n`
                                 int32[n_asv] holds each row's true length), written
                                 atomically right after stage 1;
    gpu_scores.npy               int32 flattened GPU SW scores, row-major over the
                                 shortlist (row i occupies n[0..i-1].sum() onwards),
                                 written right after stage 2;
    shortlist_meta.json          k1, min_ref_len, ref_db, fasta, n_asv, n_ref and sha256
                                 digests of the ASV (id, seq) list and the reference md5
                                 order -- the reload refuses (exit 2) on any content
                                 mismatch (paths are recorded, not compared).
  `--shortlist-in DIR` loads them and skips stage 1 (and stage 2 when gpu_scores.npy is
  present and matches); rows are restored in the original (dist, ref_idx) ascending order,
  so the GPU argsort and the final (-score, dist, idx) sort see identical input and the
  output is byte-identical to a fresh run. DIR may be --outdir itself (resume after a
  crash in stage 2/3). Without either flag nothing is written beyond the usual outputs.

Outputs (in --outdir): asv_top20_alignment_hits.json + asv_alignment_summary.csv (same
schema as before) + method_selection.json (probe, decision, run stats).

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, sys, time
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
def run_cpu(outdir, asv_ids, asv_seqs, md5s, *, k1, k2, workers, chunksize=2,
            id_db_path=None, midas_csv=None):
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
                       search_label=f"edlib_top{k1}_prefilter+biopython_SW",
                       id_db_path=id_db_path, midas_csv=midas_csv)
    return {"backend": "cpu", "prefilter_k": k1, "k2": k2,
            "cpu_seconds_edlib": round(ce, 1), "cpu_seconds_biopython": round(cb, 1),
            "wall_seconds": round(time.perf_counter() - t0, 1)}


# ============================================================================ shortlist cache
SHORTLIST_NPZ = "edlib_top{k1}_shortlist.npz"
GPU_SCORES_NPY = "gpu_scores.npy"
SHORTLIST_META = "shortlist_meta.json"
SHORTLIST_CHECK_KEYS = ("k1", "min_ref_len", "limit", "n_asv", "n_ref",
                        "asv_sha256", "ref_md5_order_sha256")


def _digest_lines(items) -> str:
    h = hashlib.sha256()
    for x in items:
        h.update(x.encode("ascii", "ignore")); h.update(b"\n")
    return h.hexdigest()


def shortlist_identity(asv_ids, asv_seqs, md5s, *, k1, min_ref_len, ref_db, fasta):
    """Everything the edlib shortlist depends on; a reload must match it exactly."""
    return {"k1": int(k1), "min_ref_len": int(min_ref_len),
            "ref_db": (os.path.abspath(ref_db) if ref_db else None),
            "fasta": (os.path.abspath(fasta) if fasta else None),
            "n_asv": len(asv_ids), "n_ref": len(md5s),
            "asv_sha256": _digest_lines(a + "\t" + q for a, q in zip(asv_ids, asv_seqs)),
            "ref_md5_order_sha256": _digest_lines(md5s)}


def _atomic_write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=2)
    os.replace(tmp, path)


def save_shortlist(outdir, cands, n_asv, k1, ident, extra=None):
    """Stage-1 result -> edlib_top{k1}_shortlist.npz (+ shortlist_meta.json), atomically."""
    idx = np.full((n_asv, k1), -1, np.int32)
    dist = np.full((n_asv, k1), -1, np.int32)
    n = np.zeros(n_asv, np.int32)
    for i in range(n_asv):
        c = cands[i]
        n[i] = len(c)
        if c:
            a = np.asarray(c, dtype=np.int64)          # [(dist, ref_idx), ...] ascending
            dist[i, :len(c)] = a[:, 0]; idx[i, :len(c)] = a[:, 1]
    path = os.path.join(outdir, SHORTLIST_NPZ.format(k1=k1))
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:                        # file handle: savez must not append .npz
        np.savez(fh, idx=idx, dist=dist, n=n)
    os.replace(tmp, path)
    _atomic_write_json(os.path.join(outdir, SHORTLIST_META),
                       {**ident, **(extra or {}), "npz": os.path.basename(path),
                        "gpu_scores": None, "n_pairs": int(n.sum())})
    print(f"[gpu] saved edlib shortlist -> {path} ({int(n.sum()):,} pairs)", flush=True)


def save_gpu_scores(outdir, scores):
    path = os.path.join(outdir, GPU_SCORES_NPY)
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        np.save(fh, np.asarray(scores, np.int32))
    os.replace(tmp, path)
    mp = os.path.join(outdir, SHORTLIST_META)
    if os.path.exists(mp):
        meta = json.load(open(mp))
        meta["gpu_scores"] = GPU_SCORES_NPY
        _atomic_write_json(mp, meta)
    print(f"[gpu] saved GPU scores -> {path} ({len(scores):,} pairs)", flush=True)


def load_shortlist(indir, ident):
    """Reload a saved shortlist; abort (exit 2) unless it was built from identical inputs.
    Returns (cands {i: [(dist, ref_idx), ...]}, scores int32[] or None, meta dict)."""
    mp = os.path.join(indir, SHORTLIST_META)
    if not os.path.exists(mp):
        print(f"[abort] --shortlist-in {indir}: no {SHORTLIST_META}", file=sys.stderr, flush=True)
        raise SystemExit(2)
    meta = json.load(open(mp))
    # content checks only: paths (fasta, ref_db) are recorded but may legitimately differ
    # (e.g. a byte-identical FASTA copied into another tree); the digests catch real drift.
    bad = [k for k in SHORTLIST_CHECK_KEYS if meta.get(k) != ident.get(k)]
    if bad:
        print(f"[abort] --shortlist-in {indir} was built from different inputs: "
              + "; ".join(f"{k}: saved={meta.get(k)!r} now={ident[k]!r}" for k in bad),
              file=sys.stderr, flush=True)
        raise SystemExit(2)
    z = np.load(os.path.join(indir, meta["npz"]))
    idx, dist, n = z["idx"], z["dist"], z["n"]
    if idx.shape != (ident["n_asv"], ident["k1"]) or n.shape != (ident["n_asv"],):
        print(f"[abort] --shortlist-in {indir}: array shape {idx.shape} does not match "
              f"({ident['n_asv']}, {ident['k1']})", file=sys.stderr, flush=True)
        raise SystemExit(2)
    cands = {}
    for i in range(ident["n_asv"]):
        m = int(n[i])                                   # python ints, as _scan_asv emits
        cands[i] = list(zip(dist[i, :m].tolist(), idx[i, :m].tolist()))
    scores = None
    sp = os.path.join(indir, GPU_SCORES_NPY)
    if os.path.exists(sp):
        s = np.load(sp)
        if s.dtype == np.int32 and s.shape == (int(n.sum()),):
            scores = s
        else:
            print(f"[gpu] {sp} has shape {s.shape} dtype {s.dtype}, expected "
                  f"({int(n.sum())},) int32 -> recomputing stage 2", flush=True)
    print(f"[gpu] loaded edlib shortlist from {indir} ({int(n.sum()):,} pairs; "
          f"GPU scores {'reused' if scores is not None else 'absent -> stage 2 reruns'})",
          flush=True)
    return cands, scores, meta


# ============================================================================ GPU scoring path
def run_gpu(outdir, asv_ids, asv_seqs, md5s, *, k1, k2, workers, gpu_cand=40,
            pairs_batch=4_000_000, chunksize=2, id_db_path=None, midas_csv=None,
            save_shortlist_to=None, shortlist_in=None, shortlist_ident=None):
    """edlib top-k1 prefilter (CPU, parallel) -> GPU Smith-Waterman scoring of ONLY those
    k1 candidates per ASV (CuPy kernel) -> Biopython rescore the GPU-top candidates for
    identity/coords -> top-k2. Scores the shortlist, never the whole reference set.

    `save_shortlist_to` (a directory) persists the stage-1 shortlist and stage-2 scores;
    `shortlist_in` (a directory) reloads them and skips the stages they cover. Both need
    `shortlist_ident` (see shortlist_identity). With neither, behaviour is unchanged."""
    import cupy as cp
    import gpu_align as G                  # CuPy SW kernel (CUDA); imported lazily
    n_asv, n_ref = len(asv_ids), len(md5s)

    # ---- Stage 1: edlib prefilter -> per-ASV top-k1 (edit_dist, ref_idx) ----
    loaded_scores, loaded_meta = None, None
    t0 = time.perf_counter()
    if shortlist_in:
        cands, loaded_scores, loaded_meta = load_shortlist(shortlist_in, shortlist_ident)
    else:
        print(f"[gpu] stage1 edlib top-{k1} prefilter | {n_asv} ASVs x {n_ref:,} refs | workers={workers}",
              flush=True)
        cands = {}
        with Pool(workers) as pool:            # edlib-only: P globals inherited via fork, no ALIGNER needed
            done = 0
            for i, cand in pool.imap_unordered(_edlib_only_worker,
                                               [(i, asv_seqs[i], k1) for i in range(n_asv)], chunksize=chunksize):
                cands[i] = cand; done += 1
                if done % 200 == 0 or done == n_asv:
                    el = time.perf_counter() - t0
                    print(f"[gpu] prefilter {done}/{n_asv} elapsed={el:.0f}s eta={el/done*(n_asv-done):.0f}s", flush=True)
    edlib_wall = time.perf_counter() - t0
    same_dir = bool(shortlist_in and save_shortlist_to
                    and os.path.realpath(shortlist_in) == os.path.realpath(save_shortlist_to))
    if save_shortlist_to and not same_dir:      # fresh stage 1, or a reload into a new outdir
        save_shortlist(save_shortlist_to, cands, n_asv, k1, shortlist_ident,
                       extra=({"edlib_wall_seconds": round(edlib_wall, 1), "edlib_workers": workers}
                              if not shortlist_in else
                              {k: loaded_meta.get(k) for k in ("edlib_wall_seconds", "edlib_workers")}))

    # ---- Stage 2: GPU-score the flattened (ASV, candidate-ref) pairs ----
    qi, rj, spans = [], [], []
    for i in range(n_asv):
        s = len(qi)
        for _dist, ref_idx in cands[i]:
            qi.append(i); rj.append(ref_idx)
        spans.append((s, len(qi)))
    qi = np.asarray(qi, np.int32); rj = np.asarray(rj, np.int32)
    if loaded_scores is not None:
        t1 = time.perf_counter()
        scores = loaded_scores
        print(f"[gpu] stage2 skipped: reusing {len(scores):,} saved GPU scores", flush=True)
    else:
        ref_codes, ref_off = G.encode_blob(P.REF_BLOB, P.REF_OFF)
        gpu = G.GpuSW(asv_seqs, ref_codes, ref_off)
        print(f"[gpu] stage2 GPU SW scoring of {len(qi):,} (ASV,candidate) pairs "
              f"(~{len(qi)/max(n_asv,1):.0f}/ASV) — NOT the {n_ref:,}-ref exhaustive product", flush=True)
        t1 = time.perf_counter()               # scoring only, as before (setup excluded)
        scores = np.empty(len(qi), np.int32)
        for b in range(0, len(qi), pairs_batch):
            sc = gpu.score_pairs(qi[b:b + pairs_batch], rj[b:b + pairs_batch])
            scores[b:b + sc.shape[0]] = cp.asnumpy(sc)
        del gpu
        try:
            cp.get_default_memory_pool().free_all_blocks()
        except Exception:
            pass
        if save_shortlist_to:
            save_gpu_scores(save_shortlist_to, scores)
    gpu_wall = time.perf_counter() - t1

    # ---- Stage 3: per-ASV GPU-top-`gpu_cand` -> Biopython rescore for identity/coords ----
    keep = max(gpu_cand, k2)
    rescore_tasks = []
    for i in range(n_asv):
        s, e = spans[i]
        if s == e:
            rescore_tasks.append((i, asv_seqs[i], [], k2)); continue
        c = cands[i]
        # best first under the rescore's own order (score desc, edlib dist asc, ref idx asc),
        # so the pick holds the shortlist's top-k2 whatever the size of an equal-score tie
        dist = np.fromiter((x[0] for x in c), np.int64, len(c))
        ref = np.fromiter((x[1] for x in c), np.int64, len(c))
        order = np.lexsort((ref, dist, -scores[s:e].astype(np.int64)))[:keep]
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
                       search_label=f"edlib_top{k1}_prefilter+gpu_SW_score+biopython_coords",
                       id_db_path=id_db_path, midas_csv=midas_csv)
    stats = {"backend": "cuda", "prefilter_k": k1, "k2": k2, "gpu_cand": keep,
             "wall_seconds_edlib_prefilter": round(edlib_wall, 1),
             "wall_seconds_gpu_score": round(gpu_wall, 1),
             "wall_seconds_biopython_rescore": round(time.perf_counter() - t2, 1),
             "n_scored_pairs": int(len(qi))}
    if shortlist_in or save_shortlist_to:
        stats["shortlist_in"] = shortlist_in
        stats["stage1_skipped"] = bool(shortlist_in)
        stats["stage2_skipped"] = loaded_scores is not None
        stats["shortlist_saved_to"] = save_shortlist_to
        if loaded_meta is not None:      # timing of the run that built the reused shortlist
            stats["shortlist_edlib_wall_seconds"] = loaded_meta.get("edlib_wall_seconds")
            stats["shortlist_edlib_workers"] = loaded_meta.get("edlib_workers")
    return stats


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
    ap.add_argument("--k2", type=int, default=100,
                    help="final hits kept per ASV (default 100, the standard record depth since the "
                         "depth-100 rerun; region_validation_d100/RUNBOOK_d100.md D1)")
    ap.add_argument("--min-ref-len", type=int, default=200)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 1) - 4))
    ap.add_argument("--gpu-cand", type=int, default=200,
                    help="GPU path: best GPU-scored candidates Biopython-rescored per ASV for final "
                         "identity/coords (default 200: the GPU prunes the 500-candidate shortlist "
                         "to 200; picked under the rescore's own tie-break, so the record equals a "
                         "rescore of the whole shortlist; see GPU PRUNING AND THE GUARD)")
    ap.add_argument("--ref-db", default=None,
                    help="md5->seq JSON replacing the default BV-BRC reference DB")
    ap.add_argument("--ref-id-db", default=None,
                    help="md5->BV-BRC header JSON matching --ref-db")
    ap.add_argument("--save-shortlist", action="store_true",
                    help="CUDA path: persist the edlib top-k1 shortlist (after stage 1) and the GPU "
                         "scores (after stage 2) in --outdir, so a later run at another --k2/--gpu-cand "
                         "can skip them with --shortlist-in (see SHORTLIST CACHE)")
    ap.add_argument("--shortlist-in", default=None, metavar="DIR",
                    help="CUDA path: reuse the shortlist (and GPU scores if present) saved in DIR by "
                         "--save-shortlist; aborts with exit 2 unless DIR was built from the same "
                         "FASTA, reference set, --prefilter-k and --min-ref-len")
    ap.add_argument("--force-backend", default="auto", choices=["auto", "cuda", "metal", "cpu"],
                    help="override the SW-scoring backend (prefilter is always edlib top-k)")
    ap.add_argument("--limit", type=int, default=0, help="run on first N ASVs (0=all)")
    ap.add_argument("--explain", action="store_true", help="probe hardware, print the decision, exit")
    args = ap.parse_args()
    # Restricted reference DB: both are read at call time from the module namespace, so
    # rebinding is sufficient; they are ALSO passed explicitly below (see module docstring).
    if args.ref_db:
        P.DB_MD5_SEQ = args.ref_db
    if args.ref_id_db:
        P.DB_MD5_ID = args.ref_id_db
    # Tie-margin guard (see module docstring): the GPU path rescores max(gpu_cand, k2)
    # candidates, so gpu_cand < k2 silently narrows the top-score tie cluster.
    if args.gpu_cand < args.k2:
        print(f"[abort] --gpu-cand ({args.gpu_cand}) < --k2 ({args.k2}): the GPU path "
              f"Biopython-rescores only max(gpu_cand, k2) candidates, so the top-score tie "
              f"cluster would be cut by the GPU argsort. Raise --gpu-cand to >= --k2.",
              file=sys.stderr, flush=True)
        raise SystemExit(2)          # exit code 2 is the documented contract for this misuse

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
    if (args.save_shortlist or args.shortlist_in) and effective != "cuda":
        print(f"[abort] --save-shortlist / --shortlist-in are implemented for the CUDA path only "
              f"(effective backend here: {effective})", file=sys.stderr, flush=True)
        raise SystemExit(2)
    if args.shortlist_in and not os.path.isdir(args.shortlist_in):
        print(f"[abort] --shortlist-in {args.shortlist_in}: not a directory", file=sys.stderr, flush=True)
        raise SystemExit(2)

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
                                     "gpu_cand": args.gpu_cand, "limit": args.limit, "fasta": args.fasta,
                                     "ref_db": args.ref_db, "ref_id_db": args.ref_id_db,
                                     "save_shortlist": args.save_shortlist,
                                     "shortlist_in": args.shortlist_in}}
    with open(os.path.join(args.outdir, "method_selection.json"), "w") as fh:
        json.dump(sel, fh, indent=2)
    if args.explain:
        print("[explain] wrote method_selection.json; no alignment run.", flush=True)
        return

    wall0 = time.perf_counter()
    print(f"[load] references from {os.path.basename(P.DB_MD5_SEQ)} ...", flush=True)
    md5s, P.REF_BLOB, P.REF_OFF, P.REF_LEN = P.load_references(args.min_ref_len, db_path=args.ref_db)
    print(f"[load] {len(md5s):,} unique refs (>= {args.min_ref_len} bp)", flush=True)
    asv_ids, asv_seqs = P.load_asvs(args.fasta)
    if args.limit:
        asv_ids, asv_seqs = asv_ids[:args.limit], asv_seqs[:args.limit]
    print(f"[load] {len(asv_ids):,} ASVs", flush=True)

    if effective == "cuda":
        ident = None
        if args.save_shortlist or args.shortlist_in:
            ident = shortlist_identity(asv_ids, asv_seqs, md5s, k1=args.prefilter_k,
                                       min_ref_len=args.min_ref_len, ref_db=args.ref_db,
                                       fasta=args.fasta)
            ident["limit"] = args.limit
        run_stats = run_gpu(args.outdir, asv_ids, asv_seqs, md5s, k1=args.prefilter_k, k2=args.k2,
                            workers=args.workers, gpu_cand=args.gpu_cand,
                            id_db_path=args.ref_id_db, midas_csv=args.taxonomy_csv,
                            save_shortlist_to=(args.outdir if args.save_shortlist else None),
                            shortlist_in=args.shortlist_in, shortlist_ident=ident)
    else:
        run_stats = run_cpu(args.outdir, asv_ids, asv_seqs, md5s, k1=args.prefilter_k, k2=args.k2,
                            workers=args.workers,
                            id_db_path=args.ref_id_db, midas_csv=args.taxonomy_csv)

    sel["run_stats"] = run_stats
    sel["wall_seconds_total"] = round(time.perf_counter() - wall0, 1)
    sel["n_asvs"] = len(asv_ids); sel["n_unique_refs"] = len(md5s)
    with open(os.path.join(args.outdir, "method_selection.json"), "w") as fh:
        json.dump(sel, fh, indent=2)
    print(f"[done] backend={effective} prefilter=top-{args.prefilter_k} "
          f"total wall={sel['wall_seconds_total']:.0f}s -> {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
