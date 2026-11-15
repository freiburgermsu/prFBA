#!/usr/bin/env python
"""Render the manuscript's F13 (fig:picrust2): prFBA vs PICRUSt2 at depth 100, as shipped.

    cd manuscript_figures/p2 && ~/Documents/py_venv/bin/python render_comparison_d100.py
    -> ../out/F13_prfba_vs_picrust2.{pdf,png}; copy the PDF to
       prFBA_paper/images/prfba/F13_prfba_vs_picrust2.pdf

A copy of region_validation_d100/picrust2_comparison/scripts/render_comparison.py (which reads the
`prfba_exclude` key: the shipped exact-tie union, source excluded) that COMPUTES the two numbers
that script hard-codes, so a re-scored comparison_summary.json renders a consistent figure:
  * the suptitle n: EC-scoreable amplicons, by_region.ec.all.all.prfba_exclude.n (asserted equal
    to PICRUSt2's n);
  * panel C's median-dF1 range: min and max over regions of paired_tests_ec_f1
    "<region>|prfba_exclude_vs_picrust2|conditional".median_diff, at two decimals (the text quotes
    the same values at three, L521), and "p ~ 0" only if every Holm p is 0 (else the largest Holm p).
The original scratch producer of F13 is lost; this script reproduces the published PDF pixel for
pixel from the committed summary (re-score W4 control: pdftoppm -r 150, 0 differing pixels).

Three panels, one measure per axis (never two scales on one axis):
  A  EC precision vs PICRUSt2's NSTI distance bins -- where PICRUSt2 degrades
  B  prFBA coverage over the same bins -- what it does instead of degrading
  C  per-region EC F1 when each method answers (conditional)

Colors are the first two slots of the documented categorical palette, used
unchanged (blue = prFBA, orange = PICRUSt2).
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PRFBA, PICRUST2 = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "0.88"
NSTI_ORDER = ["<0.01", "0.01-0.05", "0.05-0.15", "0.15-0.5", ">=0.5"]
IDENT_ORDER = ["<0.865", "0.865-0.945", "0.945-0.987", "0.987-0.995", ">=0.995"]
DPI = 300
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
SUMMARY = os.path.join(REPO, "region_validation_d100", "picrust2_comparison", "data",
                       "comparison_summary.json")
METHOD = "prfba_exclude"
N_GENOMES = 10000        # the panel size (fixed by design; not stored in the summary)
STEM = "F13_prfba_vs_picrust2"


def value(block, scope, metric):
    v = (block or {}).get(scope, {}).get(metric)
    return v["value"] if v else np.nan


def labels(s, method=METHOD):
    """The numbers the figure prints in text: scoreable n, region count, the median-dF1 range and
    the Wilcoxon p wording, all from the summary."""
    allb = s["by_region"]["ec"]["all"]["all"]
    n = allb[method]["n"]
    if allb["picrust2"]["n"] != n:
        raise SystemExit(f"scoreable n differs: {method} {n} vs picrust2 {allb['picrust2']['n']}")
    regions = [r for r in s["by_region"]["ec"] if r != "all"]
    pt = s["paired_tests_ec_f1"]
    keys = [f"{r}|{method}_vs_picrust2|conditional" for r in regions]
    missing = [k for k in keys if k not in pt]
    if missing:
        raise SystemExit(f"paired tests missing: {missing}")
    med = [pt[k]["median_diff"] for k in keys]
    ph = [pt[k].get("p_holm", pt[k]["p_value"]) for k in keys]
    lo, hi = min(med), max(med)
    sign = lambda v: f"{v:+.2f}"   # noqa: E731
    p_text = "p ≈ 0" if max(ph) == 0 else f"Holm p ≤ {max(ph):.1g}"
    return {"n": n, "n_text": f"{n:,}", "n_regions": len(regions),
            "median_dF1_min": lo, "median_dF1_max": hi,
            "dF1_text": f"{sign(lo)} to {sign(hi)}", "p_text": p_text,
            "all_regions_positive": all(m > 0 for m in med)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", default=SUMMARY)
    ap.add_argument("--outdir", default=os.path.join(HERE, "..", "out"))
    ap.add_argument("--print-labels", action="store_true", help="print the computed text and exit")
    args = ap.parse_args()
    s = json.load(open(args.summary))
    lab = labels(s)
    print(f"[render_comparison_d100] {args.summary}: {json.dumps(lab, ensure_ascii=False)}")
    if args.print_labels:
        return
    os.makedirs(args.outdir, exist_ok=True)
    nsti = s["by_stratum"]["ec"]["nsti_bin"]
    regions = [r for r in s["by_region"]["ec"] if r != "all"]
    regions.sort(key=lambda r: -value(s["by_region"]["ec"][r]["all"][METHOD],
                                     "conditional", "f1"))

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2,
                         "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": INK2, "ytick.color": INK2})
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.1), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.0, 1.0, 1.35]})

    # --- A: precision vs novelty ------------------------------------------- #
    ax = axes[0]
    x = np.arange(len(NSTI_ORDER))
    for method, color, label in ((("picrust2"), PICRUST2, "PICRUSt2"),
                                 (METHOD, PRFBA, "prFBA (exclude-self)")):
        y = [value(nsti[b].get(method), "conditional", "precision") for b in NSTI_ORDER]
        ax.plot(x, y, marker="o", lw=1.8, ms=5, color=color, label=label)
        ax.annotate(f"{y[-1]:.2f}", (x[-1], y[-1]), textcoords="offset points",
                    xytext=(6, -2), color=color, fontsize=8, fontweight="bold")
        ax.annotate(f"{y[0]:.2f}", (x[0], y[0]), textcoords="offset points",
                    xytext=(-2, 7), color=color, fontsize=8)
    ax.set_xticks(x, NSTI_ORDER, rotation=20, ha="right")
    ax.set_ylim(0, 1)
    ax.set_xlabel("PICRUSt2 NSTI (phylogenetic distance to reference)")
    ax.set_ylabel("EC precision, when the method answers")
    ax.set_title("A  Precision collapses with novelty\n(PICRUSt2 answers regardless)",
                 loc="left", fontsize=10)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower left", fontsize=8)

    # --- B: prFBA coverage on its own novelty axis ------------------------- #
    # Identity to the nearest genuinely different organism -- the axis prFBA's
    # family floor is defined on. (Coverage is lowest in the >=0.995 bin because
    # self-exclusion strips exactly those near-identical references.)
    ax = axes[1]
    ident = s["by_stratum"]["ec"]["ident_bin"]
    ib = [b for b in IDENT_ORDER if b in ident]
    xb = np.arange(len(ib))
    cov = [ident[b][METHOD]["coverage"] for b in ib]
    prec = [value(ident[b]["picrust2"], "conditional", "precision") for b in ib]
    ax.bar(xb, cov, 0.62, color=PRFBA)
    for xi, c in zip(xb, cov):
        ax.annotate(f"{c:.2f}", (xi, c), textcoords="offset points", xytext=(0, 3),
                    ha="center", fontsize=8, color=INK2)
    ax.axhline(1.0, color=PICRUST2, lw=1.8, ls="--")
    ax.annotate(f"PICRUSt2 always answers (1.00); in the floor bin\nits precision is "
                f"{prec[0]:.2f}, where prFBA answers nothing",
                (0.97, 0.955), xycoords="axes fraction", ha="right", va="top",
                fontsize=8, color=PICRUST2)
    ax.set_xticks(xb, ib, rotation=20, ha="right")
    ax.set_ylim(0, 1.18)
    ax.set_xlabel("identity to nearest different organism (prFBA's novelty axis)")
    ax.set_ylabel("prFBA coverage (fraction answered)")
    ax.set_title("B  prFBA abstains instead\n(below its 0.865 family floor it answers "
                 "nothing)", loc="left", fontsize=10)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)

    # --- C: per-region conditional F1 -------------------------------------- #
    ax = axes[2]
    y = np.arange(len(regions))
    h = 0.38
    f_prf = [value(s["by_region"]["ec"][r]["all"][METHOD], "conditional", "f1")
             for r in regions]
    f_p2 = [value(s["by_region"]["ec"][r]["all"]["picrust2"], "conditional", "f1")
            for r in regions]
    ax.barh(y - h / 2, f_prf, h, color=PRFBA, label="prFBA (exclude-self)")
    ax.barh(y + h / 2, f_p2, h, color=PICRUST2, label="PICRUSt2")
    for yi, (a, b) in enumerate(zip(f_prf, f_p2)):
        ax.annotate(f"{a:.3f}", (a, yi - h / 2), textcoords="offset points",
                    xytext=(3, -3), fontsize=7.5, color=INK2)
        ax.annotate(f"{b:.3f}", (b, yi + h / 2), textcoords="offset points",
                    xytext=(3, -3), fontsize=7.5, color=INK2)
    ax.set_yticks(y, regions)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.2)   # headroom beyond the longest bar for value labels + legend
    ax.set_xlabel("EC F1 when the method answers (conditional)")
    ax.set_title("C  prFBA is more accurate in every region when it answers\n"
                 f"(paired Wilcoxon {lab['p_text']}; median ΔF1 {lab['dF1_text']})",
                 loc="left", fontsize=10)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    # below the axis: every row carries two bars, so no in-panel spot is free
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12),
              ncol=2, fontsize=8)

    fig.suptitle(f"prFBA vs PICRUSt2 — EC-space functional prediction on {lab['n_text']} "
                 f"benchmark amplicons ({N_GENOMES:,} genomes × {lab['n_regions']} 16S regions)",
                 fontsize=11, x=0.005, ha="left")
    png = os.path.join(args.outdir, STEM + ".png")
    pdf = os.path.join(args.outdir, STEM + ".pdf")
    fig.savefig(png, dpi=DPI, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    print(f"[render_comparison_d100] -> {png}, {pdf}")


if __name__ == "__main__":
    main()
