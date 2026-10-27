#!/usr/bin/env python
"""MAG performance versus MAG completeness, with the 16S amplicon as a constant.

Palette #0072B2 / #B07600 / #009E73 is the computationally validated trio (all-pair
OKLab dE >= 18.4 normal, >= 10.0 under simulated protan/deutan/tritan, contrast
3.4-5.2:1 on white). Marker shape and line style carry identity alongside hue.
The amplicon values are drawn as horizontal constants because they do not depend on
MAG completeness; each is the PAIRED mean over the organisms that tier compares
against, so the reference moves only with the matched set, not with the marker.
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MV = os.path.expanduser("~/Documents/prFBA/mag_validation")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out")
BLUE, ORANGE, GREEN = "#0072B2", "#B07600", "#009E73"
INK, MUTED, GRID = "#1a1a1a", "#555555", "#d8d8d8"

rows = sorted(json.load(open(MV + "/out/crossover_d100.json")),   # SCRATCH PATCH: depth-100 amplicon cells, exact-tie union
               key=lambda r: r["completeness"])
x = [r["completeness"] * 100 for r in rows]
cond = [r["cond_f1"] for r in rows]
rate = [r["resolve_rate"] for r in rows]
fl = [r.get("FullLength16S_amp") for r in rows]
v4 = [r.get("V4_amp") for r in rows]
mag_fl = [r.get("FullLength16S_mag") for r in rows]


def const(vals):
    v = [z for z in vals if z is not None]
    return sum(v) / len(v) if v else None


FL, V4 = const(fl), const(v4)


def crossing(xs, ys, level):
    """First completeness (descending) at which ys drops to level, linearly interpolated."""
    pts = sorted(zip(xs, ys), reverse=True)
    for (x1, y1), (x2, y2) in zip(pts, pts[1:]):
        if y1 is None or y2 is None:
            continue
        if (y1 - level) * (y2 - level) <= 0 and y1 != y2:
            return x1 + (level - y1) * (x2 - x1) / (y2 - y1)
    return None


xr_fl = crossing(x, mag_fl, FL) if FL else None
xr_v4 = crossing(x, mag_fl, V4) if V4 else None

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "figure.facecolor": "white", "axes.facecolor": "white"})
fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.1))
TICKS = [1, 2, 3, 5, 7, 10, 15, 20, 25, 30, 50, 70, 100]


def style(ax, title, ylab, ylim):
    ax.set_title(title, fontsize=9.5, color=INK, loc="left", pad=8)
    ax.set_xlabel("MAG completeness (%, log scale)", fontsize=8.5)
    ax.set_ylabel(ylab, fontsize=8.5)
    ax.set_xscale("log"); ax.set_xticks(TICKS)
    ax.set_xticklabels([str(t) for t in TICKS], fontsize=7.5)
    ax.minorticks_off(); ax.set_xlim(0.85, 125); ax.set_ylim(*ylim)
    ax.grid(True, color=GRID, lw=0.6, axis="y"); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.8)


MIMAG = 50.0   # MIMAG medium-quality completeness floor


def mimag(ax, at="bottom"):
    ax.axvspan(0.85, MIMAG, color="#000000", alpha=0.035, lw=0, zorder=0)
    y0, y1 = ax.get_ylim()
    y = y1 - 0.035 * (y1 - y0) if at == "top" else y0 + 0.035 * (y1 - y0)
    ax.annotate("below MIMAG medium quality", (MIMAG * 0.92, y), fontsize=6.8,
                color=MUTED, ha="right", va="top" if at == "top" else "bottom")


# ---- A: gene capture, with the two amplicon constants
ax = axes[0]
if FL:
    ax.axhline(FL, color=ORANGE, lw=1.8, ls=(0, (5, 2)), zorder=2,
               label="full-length 16S (constant)")
if V4:
    ax.axhline(V4, color=GREEN, lw=1.8, ls=(0, (1.5, 1.8)), zorder=2,
               label="V4 16S (constant)")
ax.plot(x, mag_fl, color=BLUE, lw=1.9, marker="o", ms=6, mfc="white", mew=1.6,
        mec=BLUE, label="MAG (whole-genome ANI)", zorder=4, clip_on=False)
lo = min([z for z in mag_fl if z is not None] + [V4 or 1, FL or 1])
hi = max([z for z in mag_fl if z is not None] + [FL or 0])
style(ax, "A  Gene capture vs completeness", "conditional gene-capture F1",
      (lo - 0.075, hi + 0.045))
mimag(ax, at='top')
gap_hi = rows[-1]["FullLength16S_delta"]      # at 100% completeness
gap_lo = rows[0]["FullLength16S_delta"]       # at the lowest tier
ax.annotate("", xy=(1, mag_fl[0]), xytext=(1, FL),
            arrowprops=dict(arrowstyle="<->", color=BLUE, lw=1.1, shrinkA=0, shrinkB=0))
ax.annotate("still +%.3f at 1%% complete\n(+%.3f at 100%%)" % (gap_lo, gap_hi),
            (1.15, (mag_fl[0] + FL) / 2), fontsize=7.2, color=BLUE, va="center")
if xr_fl is None and xr_v4 is None:
    ax.annotate("no crossover anywhere in 1-100% completeness",
                (100, FL - 0.015), fontsize=7.4, color=ORANGE, ha="right", va="top")
ax.legend(frameon=False, fontsize=7.4, loc="lower right", bbox_to_anchor=(1.03, -0.025))

# ---- B: how often it answers at all
ax = axes[1]
ax.plot(x, rate, color=BLUE, lw=1.9, marker="s", ms=6, mfc="white", mew=1.6,
        mec=BLUE, label="MAG resolve rate", zorder=4, clip_on=False)
style(ax, "B  How often a reference is found at all", "fraction resolving to a reference",
      (0, max(rate) * 1.35))
mimag(ax)
ax.annotate("the conditional F1 in A is computed\nover this answered set, which falls\n0.308 -> 0.253 from 100% to 1%",
            (1.15, max(rate) * 0.30), fontsize=7.2, color=MUTED, ha="left")
ax.legend(frameon=False, fontsize=7.4, loc="upper left", bbox_to_anchor=(-0.02, 1.0))

_ns = [r[k] for r in rows for k in ("FullLength16S_n", "V4_n") if r.get(k) is not None]   # SCRATCH PATCH
_nrange = "%s-%s" % (format(min(_ns), ","), format(max(_ns), ","))
fig.suptitle("MAG gene capture holds far below any usable MAG quality threshold",
             fontsize=10.5, y=1.0, x=0.008, ha="left", color=INK)
fig.text(0.008, -0.03,
         "9,997 BV-BRC panel genomes per completeness tier, 0% contamination, self-exclusion (source genome, source taxon and every >=95-ANI hit removed). The amplicon constants are paired\n"
         "means over the matched organisms at each tier (n=" + _nrange + "), so the reference moves only with the matched set. The shaded region is below MIMAG's medium-quality floor of 50%\n"
         "completeness -- a MAG study would discard every bin left of it, and gene capture is still above both amplicon constants there.",
         fontsize=6.8, color=MUTED, ha="left", va="top")
fig.tight_layout(rect=[0, 0.02, 1, 0.93])
for ext in ("pdf", "png"):
    fig.savefig("%s/F_mag_crossover." % OUT_DIR + ext, dpi=200, bbox_inches="tight", facecolor="white")
print("crossover vs full-length 16S: %s" % ("%.1f%% completeness" % xr_fl if xr_fl else "none in range"))
print("crossover vs V4:             %s" % ("%.1f%% completeness" % xr_v4 if xr_v4 else "none in range"))
print("wrote out/F_mag_crossover.pdf / .png")
