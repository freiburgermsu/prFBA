#!/usr/bin/env python
"""Assemble the data for the panel-sampling map.

The background is the set that was actually eligible for clustering: full-length
(>= 1400 bp, non-atypical) 16S sequences belonging to an admissible organism
(domain in {Archaea, Bacteria} with a named phylum). project_space.py projected a
SUPERSET of this -- its mask omitted the lineage gate -- so the gate is re-applied
here and the recovered count is asserted against the number select_kmeans.py
recorded, which is what makes this the same population the panel was drawn from.

One point per picked genome, located at its centroid-nearest sequence
(selection_kmeans.json's best_16s_md5 IS that sequence: select_kmeans.py writes
`best_16s_md5=str(r.md5)` for the row it walked to from the centroid). The
original panel is placed by its own best_16s_md5 for a like-for-like comparison.
"""
import json, os, sys
import numpy as np, pandas as pd

R = os.path.expanduser("~/Documents/prFBA")
sys.path.insert(0, R + "/kmeans_panel")
import select_kmeans as SK

OUT = R + "/two_panel"
def log(m): print("[figdata] %s" % m, flush=True)

xy = np.load(OUT + "/space_xy.npy")
meta = pd.read_parquet(OUT + "/space_meta.parquet")
assert len(xy) == len(meta)
log("projected %d sequences" % len(xy))

df = pd.read_parquet(R + "/index.parquet")
sub = df[(df.seq_len >= 1400) & (~df.atypical_len)].reset_index(drop=True)
assert len(sub) == len(meta), (len(sub), len(meta))
assert (sub.row_id.to_numpy() == meta.row_id.to_numpy()).all(), "projection row order drifted"

tids = sorted({int(x) for x in sub.taxon_id.dropna().unique()})
log("resolving %d distinct taxa through the same lineage gate ..." % len(tids))
ok = {t for t in tids if SK.lineage(t) is not None}
clust = (sub.taxon_id.notna() & sub.taxon_id.astype("Int64").isin(ok)).to_numpy()

want = json.load(open(R + "/kmeans_panel/selection_kmeans.json"))["meta"]["n_clusterable_sequences"]
log("clusterable recovered %d / recorded %d" % (clust.sum(), want))
assert int(clust.sum()) == int(want), "failed to reproduce the clustering input"

# --- picks, one point per genome, at its centroid-nearest sequence -------------
pos = {m: i for i, m in enumerate(sub.md5.astype(str).to_numpy())}

# The two selections do not share an md5 convention: index.parquet (and so
# selection_kmeans.json, which copies from it) keys on md5(sequence.upper()),
# while selection_10k.json carries the md5 of the raw lower-case sequence from
# 16S_md5_seq.json. Re-hash through the sequence to bridge them -- joining on the
# md5 string alone silently locates 0/10,000 of the original panel.
import hashlib
_raw = json.load(open("/home/freiburger/Documents/codiffusion_bioreactor/"
                      "model_inputs/16S_md5_seq.json"))
RAW2UP = {k: hashlib.md5(v.upper().encode()).hexdigest() for k, v in _raw.items()}
del _raw

def locate(selfile, label):
    gs = json.load(open(selfile))["genomes"]
    idx, miss = [], 0
    for g in gs:
        m = str(g["best_16s_md5"])
        i = pos.get(m)
        if i is None:
            i = pos.get(RAW2UP.get(m, ""))
        if i is None: miss += 1
        else: idx.append(i)
    idx = np.array(idx, dtype=np.int64)
    log("%s: %d genomes -> %d located (%d md5 not in the full-length set), %d on the clusterable set"
        % (label, len(gs), len(idx), miss, int(clust[idx].sum())))
    return idx

km = locate(R + "/kmeans_panel/selection_kmeans.json", "k-means")
og = locate(R + "/region_validation/data/selection_10k.json", "original")
# The k-means panel is on the clustering input by construction, so a miss there is
# a bug. The original panel is NOT: it was drawn by taxonomy-stratified sampling
# over all 16S, so some of its organisms' best sequence is short/atypical or their
# lineage fails the gate. Those picks have no position in this embedding and are
# reported as off-frame rather than dropped silently -- the asymmetry is a result.
assert len(km) == 10000 and clust[km].all(), "a k-means pick is off the clustering input"
og_on = og[clust[og]]
off_frame = 10000 - len(og_on)
log("original off-frame: %d of 10,000 (%d best-16S not full-length/typical, %d fail the lineage gate)"
    % (off_frame, 10000 - len(og), len(og) - len(og_on)))

np.savez_compressed(OUT + "/figdata_space.npz",
                    xy=xy.astype(np.float32), clusterable=clust, km=km, og=og_on)
json.dump(dict(n_projected=int(len(xy)), n_clusterable=int(clust.sum()),
               n_km_on_frame=int(len(km)),
               n_og_on_frame=int(len(og_on)), n_og_off_frame=int(off_frame),
               n_og_not_full_length=int(10000 - len(og)),
               n_og_lineage_gate=int(len(og) - len(og_on))),
          open(OUT + "/figdata_space.json", "w"), indent=1)
log("wrote figdata_space.npz")
