#!/usr/bin/env python
"""Per-call cost of the Biopython stage-3 work, measured on an otherwise idle machine.

On a random sample (seed 7) of 150 depth-100 benchmark amplicons, in one process, with the
candidates in the rescore's order (GPU score desc, edlib distance asc, reference index asc):
  - ALIGNER.score on candidate ranks 1-100, 101-200 and 201-500 (the re-score), and
  - ALIGNER.align plus identity/coordinate extraction on ranks 1-20 and 1-100 (the
    traceback a record needs).
The stage-3 cost of each setting is then a sum of these parts:
  depth-20 run (re-score 100, align 20)      = s[1-100] + a[1-20]
  depth-100 run (re-score 500, align 100)    = s[1-100] + s[101-200] + s[201-500] + a[1-100]
  --gpu-cand 200 re-score (align 100)        = s[1-100] + s[101-200] + a[1-100]
  traceback only (align 100, the standard)   = a[1-100]

    python region_validation_d100/gpu_cand_timing.py
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import edlib_biopython_hits as P  # noqa: E402

D = os.path.join(HERE, "data", "hits", "raw")
meta = json.load(open(os.path.join(D, "method_selection.json")))["params"]
md5s, P.REF_BLOB, P.REF_OFF, P.REF_LEN = P.load_references(meta["min_ref_len"])
asv_ids, asv_seqs = P.load_asvs(meta["fasta"])
P.init_worker()
z = np.load(os.path.join(D, "edlib_top500_shortlist.npz"))
idx, dist, n = z["idx"], z["dist"], z["n"]
sc = np.load(os.path.join(D, "gpu_scores.npy"))
off = np.concatenate([[0], np.cumsum(n)])
sample = np.random.default_rng(7).choice(len(asv_ids), 150, replace=False)
ref = lambda j: P.REF_BLOB[P.REF_OFF[j]:P.REF_OFF[j + 1]].decode("ascii")  # noqa: E731

T = {"s1_100": 0.0, "s101_200": 0.0, "s201_500": 0.0, "a1_20": 0.0, "a1_100": 0.0}
N = dict.fromkeys(T, 0)
for i in sample:
    k = int(n[i])
    s, d, r = sc[off[i]:off[i] + k].astype(np.int64), dist[i, :k], idx[i, :k]
    order = [int(r[j]) for j in np.lexsort((r, d, -s))]
    q = asv_seqs[i]
    for key, band in (("s1_100", order[:100]), ("s101_200", order[100:200]),
                      ("s201_500", order[200:500])):
        refs = [ref(j) for j in band]
        t0 = time.perf_counter()
        for rs in refs:
            P.ALIGNER.score(q, rs)
        T[key] += time.perf_counter() - t0
        N[key] += len(refs)
    for key, m in (("a1_20", 20), ("a1_100", 100)):
        refs = [ref(j) for j in order[:m]]
        t0 = time.perf_counter()
        for rs in refs:
            P._identity_and_coords(P.ALIGNER.align(q, rs)[0])
        T[key] += time.perf_counter() - t0
        N[key] += len(refs)
per = {k: 1000 * T[k] / N[k] for k in T}           # ms per call
asv = {k: T[k] / len(sample) for k in T}           # s per amplicon
cost = {"depth-20 run (re-score 100, align 20)": asv["s1_100"] + asv["a1_20"],
        "depth-100 run (re-score 500, align 100)": asv["s1_100"] + asv["s101_200"]
        + asv["s201_500"] + asv["a1_100"],
        "--gpu-cand 200 re-score (align 100)": asv["s1_100"] + asv["s101_200"] + asv["a1_100"],
        "traceback only (align 100)": asv["a1_100"]}
print(json.dumps({"n_amplicons": len(sample), "ms_per_call": {k: round(v, 2) for k, v in per.items()},
                  "s_per_amplicon": {k: round(v, 3) for k, v in cost.items()}}, indent=1))
