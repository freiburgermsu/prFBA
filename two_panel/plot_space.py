#!/usr/bin/env python
"""Where each panel samples the 16S embedding space.

A, B: the 210,779 clusterable sequences in UMAP projection, shaded by LOCAL DENSITY
      measured in the full 1024-d embedding (not in the projection -- UMAP does not
      preserve density), with each panel's picks overlaid.
C:    share of each panel's picks per density decile. Shares, not counts, because
      the original panel places only 9,097 of its 10,000 organisms in this space at
      all; raw counts would read as a uniform ~9% deficit everywhere.

Density is per SEQUENCE, not per genome: a genome with eight 16S operons
contributes eight points. That is the population k-means clustered over, so it is
the right denominator here, but it is not an organism count.
"""
import json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.path.expanduser("~/Documents/prFBA")
OUT = R + "/two_panel"
BLUE, ORANGE = "#0072B2", "#B07600"
INK, MUTED, GRID = "#1a1a1a", "#555555", "#d8d8d8"

sp = np.load(OUT + "/figdata_space.npz")
kd = np.load(OUT + "/knn_density.npz")
stats = json.load(open(OUT + "/figdata_space.json"))
kstats = json.load(open(OUT + "/knn_density.json"))
assert (sp["clusterable"] == kd["clusterable"]).all(), "density and projection disagree on the frame"

clust = sp["clusterable"]
xy = sp["xy"][clust]                      # projection, clusterable rows only
knn = kd["knn21"]                         # cosine distance to the 21st neighbour
dec = kd["decile"].astype(int)            # 0 = smallest knn = DENSEST
n = len(xy)

# map a pick's row in the full-length set -> its row in the clusterable subset
pos = np.full(len(clust), -1, dtype=np.int64)
pos[np.flatnonzero(clust)] = np.arange(n)
km, og = pos[sp["km"]], pos[sp["og"]]
assert (km >= 0).all() and (og >= 0).all()

# dense = dark: shade by density percentile (the complement of sparsity rank)
order = np.argsort(np.argsort(knn))
dens_pct = 1.0 - order / (n - 1)

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "figure.facecolor": "white", "axes.facecolor": "white",
                     "svg.fonttype": "none"})
fig = plt.figure(figsize=(12.4, 5.0))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.12], wspace=0.16,
                      left=0.035, right=0.985, top=0.84, bottom=0.09)
axes = [fig.add_subplot(gs[0, i]) for i in range(3)]


def mapanel(ax, idx, color, title, sub):
    # Hexes shaded by the MEAN 1024-d density percentile of the sequences in them:
    # a readable field that still reports embedding density, not projected density.
    ax.hexbin(xy[:, 0], xy[:, 1], C=dens_pct, reduce_C_function=np.mean,
              gridsize=150, cmap="Greys", vmin=-0.15, vmax=1.0,
              linewidths=0, rasterized=True)
    ax.scatter(xy[idx, 0], xy[idx, 1], color=color, s=1.1, linewidths=0,
               alpha=0.75, rasterized=True)
    ax.set_title(title, fontsize=9.5, color=INK, loc="left", pad=30)
    ax.text(0.0, 1.006, sub, transform=ax.transAxes, fontsize=7.0, color=MUTED,
            ha="left", va="bottom", linespacing=1.35)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    for sp_ in ax.spines.values():
        sp_.set_visible(False)



mapanel(axes[0], km, BLUE, "A  k-means panel (embedding-uniform)",
        "10,000 picks, one per cluster\nall on frame")
mapanel(axes[1], og, ORANGE, "B  original panel (taxonomy-stratified)",
        "%s of 10,000 picks on frame\n%s off frame (short/atypical 16S or lineage gate)"
        % ("{:,}".format(stats["n_og_on_frame"]), "{:,}".format(stats["n_og_off_frame"])))

# --- C: share of picks per density decile ------------------------------------
ax = axes[2]
lab = 10 - dec                            # 1 = sparsest ... 10 = densest
xs = np.arange(1, 11)
kshare = np.array([(lab[km] == d).sum() for d in xs], float) / len(km)
oshare = np.array([(lab[og] == d).sum() for d in xs], float) / len(og)
w = 0.38
ax.bar(xs - w / 2, 100 * kshare, w, color=BLUE, linewidth=0, label="k-means")
ax.bar(xs + w / 2, 100 * oshare, w, color=ORANGE, linewidth=0, label="original")
ax.axhline(10, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=0)
ax.text(10.45, 10, "population\nshare", fontsize=6.8, color=MUTED, va="center", ha="left")
ax.set_title("C  sampling vs local density", fontsize=9.5, color=INK, loc="left", pad=30)
ax.set_xlabel("density decile of the embedding space  (1 = sparsest, 10 = densest)",
              fontsize=8.5)
ax.set_ylabel("share of panel's on-frame picks (%)", fontsize=8.5)
ax.set_xticks(xs); ax.set_xlim(0.4, 11.6)
ax.legend(frameon=False, fontsize=7.6, loc="upper right")
ax.yaxis.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.suptitle("Panel sampling across the BV-BRC 16S embedding space",
             fontsize=10.5, y=0.975, x=0.035, ha="left", color=INK)
# The method note that used to sit here is dropped on purpose: at \linewidth it
# renders near 3.5pt, and the LaTeX caption carries the same text legibly.

json.dump(dict(decile=list(map(int, xs)), km_share=kshare.tolist(),
               og_share=oshare.tolist(), n_km=int(len(km)), n_og=int(len(og))),
          open(OUT + "/space_deciles.json", "w"), indent=1)
for ext in ("pdf", "png"):
    fig.savefig(OUT + "/F_panel_space." + ext, dpi=220, facecolor="white")
print("decile shares (%%): k-means", np.round(100 * kshare, 2))
print("decile shares (%%): original", np.round(100 * oshare, 2))
print("wrote F_panel_space.pdf/.png")
