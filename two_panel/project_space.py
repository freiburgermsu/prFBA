#!/usr/bin/env python
"""2D projection of the BV-BRC 16S embedding space, for mapping panel sampling.

PCA (GPU, 1024 -> 50) then UMAP on the 50-d scores: the standard recipe, and
necessary here because the raw space is globally concentrated (random sequence
pairs average cosine 0.90), so a purely linear projection renders one blob.

Projects ALL 210,779 clusterable (full-length, admissible) sequences, so every
k-means pick and every original-panel pick that is clusterable appears in the
same embedding rather than a subsample of it.
"""
import json, os, time
import numpy as np, pandas as pd, torch

R = os.path.expanduser("~/Documents/prFBA")
OUT = R + "/two_panel"
SEED = 1729


def log(m): print("[proj] %s" % m, flush=True)


df = pd.read_parquet(R + "/index.parquet")
mask = (df.seq_len >= 1400) & (~df.atypical_len)
sub = df[mask].reset_index(drop=True)
km = {g["genome_id"] for g in json.load(open(R + "/kmeans_panel/selection_kmeans.json"))["genomes"]}
og = {g["genome_id"] for g in json.load(open(R + "/region_validation/data/selection_10k.json"))["genomes"]}
log("clusterable %d sequences / %d genomes" % (len(sub), sub.genome_id.nunique()))

emb = np.load(R + "/embeddings.f16.npy", mmap_mode="r")
X = torch.from_numpy(np.ascontiguousarray(emb[sub.row_id.to_numpy()])).cuda().float()
X = torch.nn.functional.normalize(X, dim=1)
log("on GPU %s" % (tuple(X.shape),))

# PCA via SVD on the centred matrix
mu = X.mean(0, keepdim=True)
U, S, V = torch.pca_lowrank(X - mu, q=50, niter=4)
Z = ((X - mu) @ V[:, :50]).cpu().numpy().astype(np.float32)
tot = ((X - mu) ** 2).sum().item()          # total variance BEFORE truncation
ev = (S[:50] ** 2).cpu().numpy() / tot
log("PCA 1024->50: first 10 components %.1f%% of TOTAL variance, all 50 %.1f%%"
    % (100 * ev[:10].sum(), 100 * ev.sum()))
del X; torch.cuda.empty_cache()

import umap
t0 = time.time()
# NOTE: random_state is deliberately NOT set. umap-learn forces n_jobs=1 when a seed
# is given, which on 218,571 points means hours single-threaded. Reproducibility is
# preserved by ARCHIVING the coordinates (space_xy.npy) rather than by re-deriving
# them from a seed; the figure is drawn from that saved array.
red = umap.UMAP(n_neighbors=30, min_dist=0.15, metric="euclidean", verbose=True)
E = red.fit_transform(Z)
log("UMAP done in %.0fs" % (time.time() - t0))

np.save(OUT + "/space_xy.npy", E.astype(np.float32))
sub[["row_id", "genome_id", "seq_len"]].assign(
    in_kmeans=sub.genome_id.isin(km).values,
    in_original=sub.genome_id.isin(og).values,
).to_parquet(OUT + "/space_meta.parquet")
log("wrote space_xy.npy (%s) and space_meta.parquet" % (E.shape,))
log("  k-means sequences: %d ; original-panel sequences: %d"
    % (int(sub.genome_id.isin(km).sum()), int(sub.genome_id.isin(og).sum())))
