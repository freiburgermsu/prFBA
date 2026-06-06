#!/usr/bin/env python
"""
edlib_vs_biopython.py — benchmark edlib (edit distance) vs the Biopython %identity ground truth for
1000 ASVs against 97,624 region-matched BV-BRC V4-V5 inserts. PARALLELIZED across CPU cores.

Ground truth = Biopython local %identity. Brute-forcing it over the whole DB for 1000 ASVs is ~100 h,
so the Biopython ground-truth ranking per ASV is computed over a broad, partly edlib-INDEPENDENT pool =
(edlib top-50) ∪ (k-mer-8 top-50) ∪ (100 random refs); the k-mer + random parts detect edlib retrieval
misses. edlib runs over the FULL DB for every ASV.

Parallelism: one ASV per task across mp.Pool workers; the big reference arrays are loaded once in the
parent and inherited by forked workers via copy-on-write (no pickling). Reports BOTH wall time (with N
cores) and intrinsic single-core per-pair time so the edlib-vs-Biopython speedup is core-count-independent.

Outputs: edlib_mappings.json, biopython_groundtruth.json, edlib_vs_biopython_stats.json
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"                                   # avoid BLAS oversubscription under mp
import json, time, hashlib, multiprocessing as mp
from pathlib import Path
import numpy as np, pandas as pd, edlib
from Bio.Align import PairwiseAligner
from sklearn.feature_extraction.text import TfidfVectorizer
from scipy.stats import spearmanr

REPO = Path("/home/freiburger/Documents/codiffusion_bioreactor"); PRFBA = Path("/home/freiburger/Documents/prFBA")
STORE = PRFBA / "v4v5_store"
SAMPLE = 1000; EDLIB_TOPK = 50; KMER_TOPK = 50; RAND = 100; REPORT_TOPK = 20
NCORES = max(1, min(60, mp.cpu_count() - 4))

_aln = PairwiseAligner(); _aln.mode = "local"; _aln.match_score = 1; _aln.mismatch_score = 0
_aln.open_gap_score = -1; _aln.extend_gap_score = -0.5
def pid(a, b):
    if not b: return -1.0
    c = _aln.align(a, b)[0].counts(); d = c.identities + c.mismatches + c.gaps
    return 100.0 * c.identities / d if d else -1.0

# ---- globals populated in main(), inherited by forked workers (copy-on-write) ----
ASV_SEQ = INSERTS = INSERTS_B = IDX = VEC = RK = VALID_ROWS = None
N = 0

def _work(a):
    s = ASV_SEQ[a]; sb = s.encode()
    # edlib over the FULL DB (timed)
    t0 = time.perf_counter()
    ed = np.fromiter((edlib.align(sb, ib, mode="HW", task="distance")["editDistance"] if ib else 10**6
                      for ib in INSERTS_B), dtype=np.int32, count=N)
    e_order = np.argpartition(ed, EDLIB_TOPK)[:EDLIB_TOPK]
    e_order = e_order[np.argsort(ed[e_order], kind="stable")]
    t_edlib = time.perf_counter() - t0
    # candidate pool: edlib ∪ k-mer ∪ random (per-ASV reproducible rng)
    km = np.asarray((VEC.transform([s]) @ RK.T).todense()).ravel()
    ktop = np.argpartition(-km, KMER_TOPK)[:KMER_TOPK]
    r = np.random.default_rng(int.from_bytes(hashlib.md5(a.encode()).digest()[:8], "little"))
    rand = r.choice(VALID_ROWS, size=RAND, replace=False)
    pool = np.unique(np.concatenate([e_order, ktop, rand]))
    # Biopython %identity over the pool (timed)
    t0 = time.perf_counter()
    pp = {int(rr): pid(s, INSERTS[int(rr)]) for rr in pool if INSERTS[int(rr)]}
    t_bio = time.perf_counter() - t0
    gt_order = sorted(pp, key=lambda rr: -pp[rr])
    # comparison flags (edlib edit-distance ranking vs Biopython %identity ranking)
    e5, e20 = set(e_order[:5].tolist()), set(e_order[:REPORT_TOPK].tolist())
    gt5 = set(gt_order[:5])
    common = [r for r in pp]
    sp = spearmanr(-ed[common], [pp[r] for r in common]).correlation if len(common) > 5 else np.nan
    res = {
        "asv": a, "asv_len": len(s),
        "top1": int(e_order[0] == gt_order[0]),
        "top5_overlap": len(gt5 & e5) / 5.0,
        "all_gt5_in_e5": int(gt5 <= e5), "all_gt5_in_e20": int(gt5 <= e20),
        "gt1_in_e20": int(gt_order[0] in e20),
        "edlib_miss": int(gt_order[0] not in set(e_order.tolist())),
        "spearman": sp, "t_edlib": t_edlib, "t_bio": t_bio, "bio_pairs": len(pp),
        "edlib_top": [(int(r), int(d)) for r, d in zip(e_order[:REPORT_TOPK], ed[e_order[:REPORT_TOPK]])],
        "gt_top": [(int(r), float(pp[r])) for r in gt_order[:REPORT_TOPK]],
    }
    return res


def _s(v): return None if (v is None or (isinstance(v, float) and pd.isna(v))) else str(v)
def hitrec(row, score_name, score):
    m = IDX.iloc[int(row)]
    return {"row_id": int(row), "feature_id": _s(m["feature_id"]), "organism": _s(m["organism"]),
            "genome_id": _s(m["genome_id"]),
            "taxon_id": (None if pd.isna(m["taxon_id"]) else int(m["taxon_id"])),
            score_name: round(float(score), 3)}


def main():
    global ASV_SEQ, INSERTS, INSERTS_B, IDX, VEC, RK, VALID_ROWS, N
    def read_fasta(p):
        o={};nm=None;cur=[]
        for ln in Path(p).read_text().splitlines():
            if ln.startswith('>'):
                if nm: o[nm]=''.join(cur)
                nm=ln[1:].strip();cur=[]
            elif ln.strip():cur.append(ln.strip().upper())
        if nm:o[nm]=''.join(cur)
        return o
    ASV_SEQ = read_fasta(REPO / "model_inputs" / "dna-sequences.fasta")
    IDX = pd.read_parquet(STORE / "index.parquet"); N = len(IDX)
    refs = json.load(open(PRFBA / "v4v5_refs.json"))
    f2i = {(h.split("|")[1].strip() if "|" in h else h): s for h, s in refs.items()}; del refs
    INSERTS = [f2i.get(f) for f in IDX["feature_id"].values]
    INSERTS_B = [s.encode() if s else b"" for s in INSERTS]
    VALID_ROWS = np.where(np.array([bool(s) for s in INSERTS]))[0]
    print("building k-mer-8 TF-IDF index ...", flush=True)
    VEC = TfidfVectorizer(analyzer="char", ngram_range=(8, 8), lowercase=False, norm="l2", dtype=np.float32)
    RK = VEC.fit_transform([s if s else "N" for s in INSERTS])

    rng = np.random.default_rng(0)
    sample = list(rng.choice([a for a in ASV_SEQ if ASV_SEQ[a]], size=min(SAMPLE, len(ASV_SEQ)), replace=False))
    print(f"{len(sample)} ASVs vs {N:,} inserts on {NCORES} cores ...", flush=True)

    t = time.time()
    with mp.Pool(NCORES) as pool:
        results = pool.map(_work, sample, chunksize=4)
    wall = time.time() - t
    R = {r["asv"]: r for r in results}; n = len(results)

    json.dump({a: {"asv_len": R[a]["asv_len"], "top20": [hitrec(r, "edit_distance", d) for r, d in R[a]["edlib_top"]]}
               for a in sample}, open("edlib_mappings.json", "w"))
    json.dump({a: {"asv_len": R[a]["asv_len"], "top20": [hitrec(r, "pct_identity", p) for r, p in R[a]["gt_top"]]}
               for a in sample}, open("biopython_groundtruth.json", "w"))

    sumf = lambda key: sum(R[a][key] for a in sample)
    t_edlib_cpu = sumf("t_edlib"); t_bio_cpu = sumf("t_bio"); bio_pairs = sumf("bio_pairs")
    sp = np.array([R[a]["spearman"] for a in sample]); sp = sp[~np.isnan(sp)]
    bio_per_pair_ms = 1000 * t_bio_cpu / max(bio_pairs, 1)
    edlib_per_pair_us = 1e6 * t_edlib_cpu / (n * N)
    edlib_per_query_s = t_edlib_cpu / n                       # single-core, full DB
    bio_full_db_per_query_s = bio_per_pair_ms / 1000 * N
    stats = {
        "n_asv": n, "n_refs": int(N), "cores_used": NCORES,
        "pool": "edlib_top50 ∪ kmer8_top50 ∪ 100 random",
        "accuracy_capture": {
            "top1_agreement_pct": round(100 * sumf("top1") / n, 1),
            "mean_top5_overlap": round(np.mean([R[a]["top5_overlap"] for a in sample]), 3),
            "all_gt_top5_in_edlib_top5_pct": round(100 * sumf("all_gt5_in_e5") / n, 1),
            "all_gt_top5_in_edlib_top20_pct": round(100 * sumf("all_gt5_in_e20") / n, 1),
            "gt_top1_in_edlib_top20_pct": round(100 * sumf("gt1_in_e20") / n, 1),
            "edlib_retrieval_miss_pct(GT#1 not in edlib top50)": round(100 * sumf("edlib_miss") / n, 2),
            "mean_spearman_negEdit_vs_pctId": round(float(np.mean(sp)), 3),
        },
        "timing": {
            "wall_s_parallel": round(wall, 1), "cores_used": NCORES,
            "edlib_cpu_s_total": round(t_edlib_cpu, 1), "edlib_per_pair_us": round(edlib_per_pair_us, 4),
            "edlib_per_query_s_fullDB_1core": round(edlib_per_query_s, 3),
            "biopython_cpu_s_pool": round(t_bio_cpu, 1), "biopython_pairs": int(bio_pairs),
            "biopython_per_pair_ms": round(bio_per_pair_ms, 3),
            "biopython_per_query_s_fullDB_1core_extrap": round(bio_full_db_per_query_s, 1),
            "biopython_fullDB_total_h_1core_extrap": round(bio_full_db_per_query_s * n / 3600, 1),
            "edlib_speedup_vs_biopython_per_pair": round(bio_per_pair_ms * 1000 / edlib_per_pair_us, 1),
        },
    }
    json.dump(stats, open("edlib_vs_biopython_stats.json", "w"), indent=2)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
