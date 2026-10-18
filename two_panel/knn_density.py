#!/usr/bin/env python
"""Local density of each clusterable 16S sequence, in the FULL 1024-d embedding.

Density must not be read off the UMAP projection: UMAP preserves neighbourhoods,
not densities, so a visually dense island there can be sparse in the space k-means
actually optimized. The sparsity measure here is the cosine distance to the k-th
nearest neighbour among the 210,779 clusterable sequences -- the same population,
the same metric (unit-norm dot product) that select_kmeans.py clustered on.

k = 21 is the principled choice: mean cluster occupancy is 210,779/10,000 ~= 21,
so the k-th neighbour sits at roughly the scale of one cluster. k = 50 is computed
too and the two rankings are compared, so the decile assignment is not an artifact
of one k.
"""
import json, os, sys
import numpy as np, pandas as pd, torch

R = os.path.expanduser("~/Documents/prFBA")
sys.path.insert(0, R + "/kmeans_panel")
import select_kmeans as SK

OUT = R + "/two_panel"
KS = (21, 50)
def log(m): print("[knn] %s" % m, flush=True)

df = pd.read_parquet(R + "/index.parquet")
sub = df[(df.seq_len >= 1400) & (~df.atypical_len)].reset_index(drop=True)
tids = sorted({int(x) for x in sub.taxon_id.dropna().unique()})
ok = {t for t in tids if SK.lineage(t) is not None}
clust = (sub.taxon_id.notna() & sub.taxon_id.astype("Int64").isin(ok)).to_numpy()
want = json.load(open(R + "/kmeans_panel/selection_kmeans.json"))["meta"]["n_clusterable_sequences"]
assert int(clust.sum()) == int(want), (clust.sum(), want)
log("clusterable %d (matches the clustering input)" % clust.sum())

rows = sub.row_id.to_numpy()[clust]
emb = np.load(R + "/embeddings.f16.npy", mmap_mode="r")
X = torch.from_numpy(np.ascontiguousarray(emb[rows])).cuda().float()
X = torch.nn.functional.normalize(X, dim=1).half()
n = X.shape[0]
log("on GPU %s (%.1f GB)" % (tuple(X.shape), X.element_size() * X.nelement() / 1e9))

kmax = max(KS) + 1                       # +1: the point is its own nearest neighbour
out = {k: np.empty(n, dtype=np.float32) for k in KS}
CH = 2048
for s in range(0, n, CH):
    sim = (X[s:s + CH] @ X.T).float()    # cosine similarity, unit-norm rows
    top = sim.topk(kmax, dim=1, largest=True).values
    for k in KS:
        out[k][s:s + CH] = (1.0 - top[:, k]).clamp_(min=0).cpu().numpy()
    if (s // CH) % 20 == 0:
        log("  %d/%d" % (s, n))
del X; torch.cuda.empty_cache()

from scipy.stats import spearmanr
rho = spearmanr(out[21], out[50]).statistic
log("Spearman(k=21, k=50) sparsity rankings = %.4f" % rho)
d21 = pd.qcut(out[21], 10, labels=False, duplicates="drop")
d50 = pd.qcut(out[50], 10, labels=False, duplicates="drop")
agree = float((d21 == d50).mean())
log("same decile under both k: %.3f" % agree)

np.savez_compressed(OUT + "/knn_density.npz", clusterable=clust,
                    knn21=out[21], knn50=out[50], decile=d21.astype(np.int8))
json.dump(dict(n=int(n), ks=list(KS), spearman_21_50=float(rho),
               same_decile_frac=agree), open(OUT + "/knn_density.json", "w"), indent=1)
log("wrote knn_density.npz")
