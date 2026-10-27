#!/usr/bin/env python
"""Can the GPU prune the depth-100 shortlist (--gpu-cand 200) without changing a record?

Reads the depth-100 run's saved stage-1 shortlist and stage-2 GPU scores
(region_validation_d100/data/hits/raw: edlib_top500_shortlist.npz, gpu_scores.npy) and its
stored records (asv_top20_alignment_hits.json, rescored over all 500 candidates):

  1. exactness: is the GPU score equal to the stored Biopython align_score on every hit?
  2. canonical order: is every stored record exactly the top-k2 of its shortlist under the
     rescore's total order (score desc, edlib distance asc, reference index asc), computed
     from the GPU scores?
  3. the old pick: how many records would a rescore of np.argsort(-gpu_score)[:gpu_cand]
     (align_hits.py before 2026-10-27) have changed, for several gpu_cand?

If 1 and 2 hold for every record, rescoring the canonical GPU-top-gpu_cand (align_hits.py
since 2026-10-27) reproduces every record for any gpu_cand >= k2: the pick then contains
the shortlist's top-k2 under the very order the rescore applies. GPU_CAND_200.md reports
the run of 2026-10-27 and the end-to-end rerun that confirms it.

    python region_validation_d100/gpu_cand_check.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import edlib_biopython_hits as P  # noqa: E402

D = os.path.join(HERE, "data", "hits", "raw")
params = json.load(open(os.path.join(D, "method_selection.json")))["params"]
K2 = params["k2"]
md5s, *_ = P.load_references(params["min_ref_len"])
asv_ids, _seqs = P.load_asvs(params["fasta"])
ref_of = {m: i for i, m in enumerate(md5s)}
z = np.load(os.path.join(D, "edlib_top500_shortlist.npz"))
idx, dist, n = z["idx"], z["dist"], z["n"]
sc = np.load(os.path.join(D, "gpu_scores.npy"))
off = np.concatenate([[0], np.cumsum(n)])
rec = json.load(open(os.path.join(D, "asv_top20_alignment_hits.json")))

PICKS = (100, 150, 200, 250, 300)
bad_score = bad_order = n_hits = 0
changed = {g: [0, 0] for g in PICKS}          # records changed, hits replaced
for i, aid in enumerate(asv_ids):
    hits = rec[aid]["top20"]
    k = int(n[i])
    s, d, r = sc[off[i]:off[i] + k], dist[i, :k], idx[i, :k]
    canon = np.lexsort((r, d, -s.astype(np.int64)))
    got = [ref_of[h["md5"]] for h in hits]
    n_hits += len(hits)
    bad_order += got != [int(r[j]) for j in canon[:len(hits)]]
    pos = {int(r[j]): j for j in range(k)}
    bad_score += sum(float(s[pos[g]]) != float(h["align_score"]) for h, g in zip(hits, got))
    top = canon[:K2]
    for g in PICKS:
        if k > g:
            old = set(np.argsort(-s)[:g].tolist())
            miss = sum(1 for j in top if j not in old)
            if miss:
                changed[g][0] += 1
                changed[g][1] += miss
print(f"{len(asv_ids):,} records, {n_hits:,} hits")
print(f"1. GPU score != stored Biopython score: {bad_score} hits")
print(f"2. stored record != canonical top-{K2} from GPU scores: {bad_order} records")
for g in PICKS:
    print(f"3. old argsort pick, gpu_cand={g}: {changed[g][0]:,} records changed "
          f"({changed[g][1]:,} hits replaced)")
