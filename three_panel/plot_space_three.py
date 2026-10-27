#!/usr/bin/env python
"""Where each of the three panels samples the 16S embedding space.

The two-panel figure (two_panel/plot_space.py -> F_panel_space) extended with the
uniform-random panel:

A, B, C: the 210,779 clusterable sequences in UMAP projection, shaded by LOCAL DENSITY
         measured in the full 1024-d embedding (not in the projection), with each
         panel's picks overlaid, in the order the panels were built -- taxonomy-
         stratified original, k-means, uniform random.
D:       share of each panel's on-frame picks per density decile, against the share of
         all clusterable sequences in that decile.

The population share is drawn per decile, not as a flat 10% line: the 21st-neighbour
cosine distance takes only 166 distinct values over 210,779 sequences, so the decile
bins are tied and unequal (7.7-12.6% of sequences each).  The two-panel figure's flat
10% line was an approximation.

Each pick is placed at its organism's best_16s_md5 (the k-means panel's centroid-nearest
sequence; the original and random panels' representative 16S), bridging the raw and
upper-case md5 conventions as figdata_space.py does.  Picks whose sequence is not in the
clusterable set are off frame and reported, not dropped silently.

    python three_panel/plot_space_three.py            # F_panel_space_three (manuscript)
    python three_panel/plot_space_three.py --slide    # F_panel_space_three_slide.png (deck)

--slide draws the same data at slide size for presentations/make_deck.py: placed 1:1 on a
13.33-in slide, larger type, no suptitle or sub-captions (the slide title and speaker notes
carry them), and the taxonomy-stratified panel named "stratified" as the deck names it.
"""
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.path.expanduser("~/Documents/prFBA")
TP = R + "/two_panel"
OUT = R + "/three_panel"
BLUE, ORANGE, GREEN = "#0072B2", "#B07600", "#009E73"
SLIDE = "--slide" in sys.argv[1:]
INK, MUTED, GRID = "#1a1a1a", "#555555", "#d8d8d8"

sp = np.load(TP + "/figdata_space.npz")
kd = np.load(TP + "/knn_density.npz")
stats = json.load(open(TP + "/figdata_space.json"))
assert (sp["clusterable"] == kd["clusterable"]).all(), "density and projection disagree on the frame"

clust = sp["clusterable"]
xy = sp["xy"][clust]
knn = kd["knn21"]
dec = kd["decile"].astype(int)            # 0 = smallest knn = DENSEST
n = len(xy)
pos = np.full(len(clust), -1, dtype=np.int64)
pos[np.flatnonzero(clust)] = np.arange(n)
km, og = pos[sp["km"]], pos[sp["og"]]
assert (km >= 0).all() and (og >= 0).all()

# --- the random panel, located exactly as figdata_space.locate() does ----------------
df = pd.read_parquet(R + "/index.parquet")
sub = df[(df.seq_len >= 1400) & (~df.atypical_len)].reset_index(drop=True)
assert len(sub) == len(clust)
md5_to_row = {m: i for i, m in enumerate(sub.md5.astype(str).to_numpy())}
_raw = json.load(open("/home/freiburger/Documents/codiffusion_bioreactor/"
                      "model_inputs/16S_md5_seq.json"))
RAW2UP = {k: hashlib.md5(v.upper().encode()).hexdigest() for k, v in _raw.items()}
del _raw
rnd_all, rnd_miss = [], 0
for g in json.load(open(R + "/random_panel/selection_random.json"))["genomes"]:
    m = str(g["best_16s_md5"])
    i = md5_to_row.get(m)
    if i is None:
        i = md5_to_row.get(RAW2UP.get(m, ""))
    if i is None:
        rnd_miss += 1
    else:
        rnd_all.append(i)
rnd_all = np.array(rnd_all, dtype=np.int64)
rnd = pos[rnd_all[clust[rnd_all]]]
n_rnd_off = 10000 - len(rnd)
print("random panel: %d of 10,000 on frame (%d best-16S not full-length/typical, %d fail the "
      "lineage gate)" % (len(rnd), rnd_miss, len(rnd_all) - len(rnd)))

order = np.argsort(np.argsort(knn))
dens_pct = 1.0 - order / (n - 1)          # dense = dark

FS = (dict(base=12, title=13, label=12, legend=11, fig=(12.13, 4.2), top=0.9, bottom=0.16,
           left=0.0, right=0.995, wspace=0.1) if SLIDE else
      dict(base=9, title=9.5, label=8.5, legend=7.4, fig=(16.4, 4.9), top=0.82, bottom=0.11,
           left=0.022, right=0.99, wspace=0.13))
plt.rcParams.update({"font.size": FS["base"], "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "figure.facecolor": "white", "axes.facecolor": "white",
                     "svg.fonttype": "none"})
fig = plt.figure(figsize=FS["fig"])
gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.25], wspace=FS["wspace"],
                      left=FS["left"], right=FS["right"], top=FS["top"], bottom=FS["bottom"])
axes = [fig.add_subplot(gs[0, i]) for i in range(4)]


def mapanel(ax, idx, color, title, sub_):
    ax.hexbin(xy[:, 0], xy[:, 1], C=dens_pct, reduce_C_function=np.mean,
              gridsize=150, cmap="Greys", vmin=-0.15, vmax=1.0,
              linewidths=0, rasterized=True)
    ax.scatter(xy[idx, 0], xy[idx, 1], color=color, s=1.1, linewidths=0,
               alpha=0.75, rasterized=True)
    ax.set_title(title, fontsize=FS["title"], color=INK, loc="left", pad=6 if SLIDE else 30)
    if not SLIDE:
        ax.text(0.0, 1.006, sub_, transform=ax.transAxes, fontsize=7.0, color=MUTED,
                ha="left", va="bottom", linespacing=1.35)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    if SLIDE:
        ax.set_anchor("N")        # titles level with D's, not floated by the equal aspect
    for s_ in ax.spines.values():
        s_.set_visible(False)


# chronological order: the original panel, then k-means, then random
mapanel(axes[1], km, BLUE, "B  k-means" if SLIDE else "B  k-means panel (embedding-uniform)",
        "10,000 picks, one per cluster\nall on frame")
mapanel(axes[0], og, ORANGE, "A  taxonomy-stratified" if SLIDE else "A  original panel (taxonomy-stratified)",
        "{:,} of 10,000 picks on frame\n{:,} off frame (short/atypical 16S or lineage gate)"
        .format(stats["n_og_on_frame"], stats["n_og_off_frame"]))
mapanel(axes[2], rnd, GREEN, "C  uniform random" if SLIDE else "C  random panel (uniform over organisms)",
        "10,000 picks, one per organism\nall on frame" if n_rnd_off == 0 else
        "{:,} of 10,000 picks on frame\n{:,} off frame (short/atypical 16S or lineage gate)"
        .format(len(rnd), n_rnd_off))

# --- D: share of picks per density decile ------------------------------------------
ax = axes[3]
lab = 10 - dec                            # 1 = sparsest ... 10 = densest
xs = np.arange(1, 11)
share = lambda idx: np.array([(lab[idx] == d).sum() for d in xs], float) / len(idx)  # noqa: E731
kshare, oshare, rshare = share(km), share(og), share(rnd)
pshare = np.array([(lab == d).mean() for d in xs])
w = 0.27
ax.bar(xs - w, 100 * oshare, w, color=ORANGE, linewidth=0,
       label="stratified" if SLIDE else "original")
ax.bar(xs, 100 * kshare, w, color=BLUE, linewidth=0, label="k-means")
ax.bar(xs + w, 100 * rshare, w, color=GREEN, linewidth=0, label="random")
ax.hlines(100 * pshare, xs - 1.6 * w, xs + 1.6 * w, color=INK, lw=1.4, zorder=3,
          label="all sequences" if SLIDE else "all clusterable sequences")
ax.set_title("D  sampling vs local density", fontsize=FS["title"], color=INK, loc="left",
             pad=6 if SLIDE else 30)
ax.set_xlabel("density decile  (1 = sparsest)" if SLIDE else
              "density decile of the embedding space  (1 = sparsest, 10 = densest)",
              fontsize=FS["label"])
ax.set_ylabel("% of on-frame picks" if SLIDE else "share of panel's on-frame picks (%)",
              fontsize=FS["label"])
ax.set_xticks(xs); ax.set_xlim(0.4, 10.6)
ax.legend(frameon=False, fontsize=FS["legend"], loc="upper right")
ax.yaxis.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)

if not SLIDE:
    fig.suptitle("Panel sampling across the BV-BRC 16S embedding space",
                 fontsize=10.5, y=0.975, x=0.022, ha="left", color=INK)

# the shares must reproduce the deciles the three-panel comparison already reports
ref = json.load(open(OUT + "/space_deciles_three.json"))
for key, arr in (("km", kshare), ("og", oshare), ("rnd", rshare), ("population", pshare)):
    d = float(np.abs(np.array(ref[key]) - arr).max())
    assert d < 1e-3, "%s decile shares drift from space_deciles_three.json by %.4f" % (key, d)
json.dump(dict(decile=list(map(int, xs)), km_share=kshare.tolist(), og_share=oshare.tolist(),
               rnd_share=rshare.tolist(), population_share=pshare.tolist(),
               n_km=int(len(km)), n_og=int(len(og)), n_rnd=int(len(rnd)),
               n_rnd_off_frame=int(n_rnd_off), n_rnd_not_full_length=int(rnd_miss),
               n_rnd_lineage_gate=int(len(rnd_all) - len(rnd))),
          open(OUT + "/space_deciles_three_plot.json", "w"), indent=1)
for ext in (("png",) if SLIDE else ("pdf", "png")):
    fig.savefig(OUT + "/F_panel_space_three" + ("_slide." if SLIDE else ".") + ext, dpi=220,
                facecolor="white")
for nm, a in (("k-means", kshare), ("original", oshare), ("random", rshare), ("population", pshare)):
    print("%-10s %s  two sparsest %.1f%%  densest %.1f%%" % (nm, np.round(100 * a, 1),
                                                             100 * a[:2].sum(), 100 * a[9]))
print("wrote F_panel_space_three_slide.png" if SLIDE else "wrote F_panel_space_three.pdf/.png")
