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
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out")
BLUE, ORANGE, GREEN = "#0072B2", "#B07600", "#009E73"
INK, MUTED, GRID = "#1a1a1a", "#555555", "#d8d8d8"
COMP = [100, 90, 70, 50, 30]
X = [1.0, 0.9, 0.7, 0.5, 0.3]

panel = json.load(open(MV + "/out/mag_summary_full.json"))["summary"]
gtdb = json.load(open(MV + "/out_gtdb/gtdb_summary_2000.json"))["summary"]   # SCRATCH PATCH: the 2,000-organism subsample the text quotes
pair = json.load(open(MV + "/out/paired_comparison_d100.json"))   # SCRATCH PATCH: depth-100 amplicon cells, legacy reducer

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
_c1fl = pair["exclude|c1_x0|FullLength16S|f1"]["conditional"]; _c1v4 = pair["exclude|c1_x0|V4|f1"]["conditional"]   # SCRATCH PATCH
_dfl = round(_c1fl["mag"] - _c1fl["amp"], 3); _dv4 = round(_c1v4["mag"] - _c1v4["amp"], 3)
ax.text(31, 0.90, "%+.3f vs full-length 16S, %+.3f vs V4 (paired means)\npaired within organism, $p<10^{-164}$" % (_dfl, _dv4),
        fontsize=7.2, color=BLUE, va="top")
# SCRATCH PATCH: at depth 100 the V4 series (~0.615) runs under a lower-left legend;
# place the legend in the empty band between the MAG (~0.81) and full-length (~0.69) series.
ax.legend(frameon=False, fontsize=7.0, loc="center left", bbox_to_anchor=(-0.02, 0.53))

_ns = [pair["exclude|%s|%s|f1" % (tier[c], r)]["conditional"]["n"] for c in X for r in ("FullLength16S", "V4")]   # SCRATCH PATCH
_nrange = "%s-%s" % (format(min(_ns), ","), format(max(_ns), ","))
fig.suptitle("MAG reference selection is insensitive to completeness from 100% down to 30%",
             fontsize=10.5, y=1.0, x=0.008, ha="left", color=INK)
fig.text(0.008, -0.04,
         "9,997 BV-BRC genomes x 5 completeness levels = 49,985 synthetic MAGs shown (panel A and panel B's GTDB series from the 2,000-genome GTDB subsample). Self-exclude removes the source genome, its taxon\n"
         "and every hit at >=95 ANI. Panel C is paired within organism against the 16S arm's own cells (depth-100 records, exact-tie union) on the same genomes, n=%s; it covers the 31%% of organisms for which ANI\n"
         "finds a qualifying reference, which are also +0.05 easier for the 16S arm. Contaminated tiers are excluded throughout: the synthesizer's dose control fails at low completeness." % _nrange,
         fontsize=6.8, color=MUTED, ha="left", va="top")
fig.tight_layout(rect=[0, 0.02, 1, 0.94])
for ext in ("pdf", "png"):
    fig.savefig("%s/F_mag_completeness." % OUT_DIR + ext, dpi=200, bbox_inches="tight",
                facecolor="white")
print("wrote out/F_mag_completeness.pdf and .png")
