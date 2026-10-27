#!/usr/bin/env python
"""Time stage 3 (Biopython rescore) per ASV with all 500 shortlisted candidates vs the
canonical GPU-top-200, on a random sample (seed 7) of 150 depth-100 benchmark amplicons,
single process. The difference is 300 extra ALIGNER.score calls per ASV; the remainder is
200 score calls plus the 100 ALIGNER.align tracebacks, so the two give the per-call costs
reported in GPU_CAND_200.md.

    python region_validation_d100/gpu_cand_timing.py
"""
import json, sys, time
import numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import edlib_biopython_hits as P  # noqa: E402
D = os.path.join(HERE, "data", "hits", "raw")
meta = json.load(open(D + "/method_selection.json"))["params"]
md5s, P.REF_BLOB, P.REF_OFF, P.REF_LEN = P.load_references(meta["min_ref_len"])
asv_ids, asv_seqs = P.load_asvs(meta["fasta"])
P.init_worker()
z = np.load(D + "/edlib_top500_shortlist.npz")
idx, dist, n = z["idx"], z["dist"], z["n"]
sc = np.load(D + "/gpu_scores.npy")
off = np.concatenate([[0], np.cumsum(n)])
rng = np.random.default_rng(7)
sample = rng.choice(len(asv_ids), 150, replace=False)
tt = {}
same = 0
for keep in (500, 200):
    t0 = time.perf_counter()
    outs = []
    for i in sample:
        k = int(n[i]); s = sc[off[i]:off[i] + k]; d = dist[i, :k]; r = idx[i, :k]
        order = np.lexsort((r, d, -s))[:keep]
        outs.append(P._rescore(asv_seqs[i], [(int(d[j]), int(r[j])) for j in order], 100))
    tt[keep] = (time.perf_counter() - t0, outs)
same = sum(a == b for a, b in zip(tt[500][1], tt[200][1]))
print(f"stage-3 time for {len(sample)} ASVs: 500 cand {tt[500][0]:.1f}s, 200 cand {tt[200][0]:.1f}s "
      f"(ratio {tt[200][0] / tt[500][0]:.2f}); identical records {same}/{len(sample)}")
