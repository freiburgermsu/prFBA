#!/usr/bin/env python
"""Three-panel mechanism figure for the MAG-matching validation.

Each panel has ONE y-axis and at most three series. Palette
(#0072B2 / #B07600 / #009E73) was validated computationally: all-pair OKLab dE
>= 18.4 under normal vision and >= 10.0 under simulated protan/deutan/tritan
(Machado 2009, severity 1.0), contrast 3.4-5.2:1 on white. Marker shape carries
identity as well as color, so the series are distinguishable without hue.
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MV = os.path.expanduser("~/Documents/prFBA/mag_validation")
BLUE, ORANGE, GREEN = "#0072B2", "#B07600", "#009E73"
INK, MUTED, GRID = "#1a1a1a", "#555555", "#d8d8d8"
COMP = [100, 90, 70, 50, 30]
X = [1.0, 0.9, 0.7, 0.5, 0.3]

panel = json.load(open(MV + "/out/mag_summary_full.json"))["summary"]
gtdb = json.load(open(MV + "/out_gtdb/gtdb_summary.json"))["summary"]
pair = json.load(open(MV + "/out/paired_comparison.json"))

afq = [gtdb["c%g_x0_include" % c]["median_af_query"] for c in X]
afr = [gtdb["c%g_x0_include" % c]["median_af_ref"] for c in X]
res_gi = [gtdb["c%g_x0_include" % c]["resolve_rate"] for c in X]
res_ge = [gtdb["c%g_x0_exclude" % c]["resolve_rate"] for c in X]
res_pe = [panel["c%g_x0_exclude" % c]["coverage_answered"] for c in X]
tier = {1.0: "c1_x0", 0.9: "c0.9_x0", 0.7: "c0.7_x0", 0.5: "c0.5_x0", 0.3: "c0.3_x0"}
mag_f1 = [pair["exclude|%s|FullLength16S|f1" % tier[c]]["conditional"]["mag"] for c in X]
fl_f1 = [pair["exclude|%s|FullLength16S|f1" % tier[c]]["conditional"]["amp"] for c in X]
v4_f1 = [pair["exclude|%s|V4|f1" % tier[c]]["conditional"]["amp"] for c in X]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "figure.facecolor": "white", "axes.facecolor": "white",
                     "svg.fonttype": "none"})
fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.9))


def style(ax, title, ylab, ylim):
    ax.set_title(title, fontsize=9.5, color=INK, loc="left", pad=8)
    ax.set_xlabel("MAG completeness (%)", fontsize=8.5)
    ax.set_ylabel(ylab, fontsize=8.5)
    ax.set_xticks(COMP); ax.set_xlim(24, 106); ax.set_ylim(*ylim)
    ax.grid(True, color=GRID, lw=0.6, axis="y"); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.8)


def series(ax, y, color, marker, label):
    ax.plot(COMP, y, color=color, lw=1.8, marker=marker, ms=6.5, mfc="white",
            mew=1.6, mec=color, label=label, zorder=3, clip_on=False)


# A - the degradation is real, and the gate's input cannot see it
ax = axes[0]
series(ax, afr, BLUE, "o", "AF$_{ref}$ (fraction of the REFERENCE aligned)")
series(ax, afq, ORANGE, "s", "AF$_{query}$ (fraction of the MAG aligned)")
style(ax, "A  Only AF$_{ref}$ tracks completeness", "median aligned fraction (%)", (0, 126))
ax.text(29, 107, "AF$_{query}$ - what the selector gates on.\nFlat, so completeness is invisible to it.",
        fontsize=7.2, color=ORANGE, va="bottom")
ax.text(63, 34, "AF$_{ref}$ $\\approx$ 0.99 $\\times$ completeness\n- a free completeness estimator\nthe pipeline discards",
        fontsize=7.2, color=BLUE, va="bottom")
ax.legend(frameon=False, fontsize=7.0, loc="lower left", bbox_to_anchor=(-0.02, -0.01))

# B - resolve rate is flat; only the reference universe moves it
ax = axes[1]
series(ax, res_gi, GREEN, "^", "GTDB r226, self-include")
series(ax, res_ge, BLUE, "o", "GTDB r226, self-exclude")
series(ax, res_pe, ORANGE, "s", "10k panel, self-exclude")
style(ax, "B  Resolve rate is flat in completeness", "fraction resolving to a reference", (0, 1.12))
ax.text(31, 0.74, "the reference universe sets the level,\nnot MAG quality: 0.56 vs 0.31",
        fontsize=7.2, color=MUTED, va="center")
ax.legend(frameon=False, fontsize=7.0, loc="lower left", bbox_to_anchor=(-0.02, -0.01))

# C - the MAG advantage over amplicons holds at every completeness
ax = axes[2]
series(ax, mag_f1, BLUE, "o", "MAG (whole-genome ANI)")
series(ax, fl_f1, ORANGE, "s", "full-length 16S amplicon")
series(ax, v4_f1, GREEN, "^", "V4 amplicon")
style(ax, "C  MAG beats amplicons at every completeness", "conditional gene-capture F1", (0.55, 0.93))
ax.text(31, 0.90, "+0.064 vs full-length 16S, +0.133 vs V4\npaired within organism, $p<10^{-140}$",
        fontsize=7.2, color=BLUE, va="top")
ax.legend(frameon=False, fontsize=7.0, loc="lower left", bbox_to_anchor=(-0.02, -0.01))

fig.suptitle("MAG reference selection is insensitive to completeness from 100% down to 30%",
             fontsize=10.5, y=1.0, x=0.008, ha="left", color=INK)
fig.text(0.008, -0.04,
         "9,997 BV-BRC genomes x 5 completeness levels = 49,985 synthetic MAGs shown (panel B's GTDB series from a 1,000-genome subsample). Self-exclude removes the source genome, its taxon\n"
         "and every hit at >=95 ANI. Panel C is paired within organism against the 16S arm's own cells on the same genomes, n=1,510-1,573; it covers the 31% of organisms for which ANI\n"
         "finds a qualifying reference, which are also +0.06 easier for the 16S arm. Contaminated tiers are excluded throughout: the synthesizer's dose control fails at low completeness.",
         fontsize=6.8, color=MUTED, ha="left", va="top")
fig.tight_layout(rect=[0, 0.02, 1, 0.94])
for ext in ("pdf", "png"):
    fig.savefig(MV + "/out/F_mag_completeness." + ext, dpi=200, bbox_inches="tight",
                facecolor="white")
print("wrote out/F_mag_completeness.pdf and .png")
