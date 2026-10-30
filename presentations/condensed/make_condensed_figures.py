#!/usr/bin/env python
"""Figures for the condensed prFBA deck (presentations/make_deck_condensed.py).

Every value is read from the pipeline's own files, through make_deck.load() where the full
deck already reads it, so the two decks cannot disagree. Each figure is drawn at slide size
(placed 1:1 on the 13.33-in slide) in the deck's palette and a Calibri-metric face.

    python presentations/condensed/make_condensed_figures.py [names...]

  workflow   F_workflow.png    one real amplicon through align, select, merge and score
  arms       F_arms.png        the four validation arms: query, reference, held out, truth
  sankey     F_sankey.png      67,991 amplicons: region -> genus outcome, source present/removed
  picrust2   F_picrust2.png    how PICRUSt2 was run, as shipped and in the controlled test
  heatmap    F_heatmap.png     region x measure, source present vs removed
  space      F_space.png       where in 16S space prediction fails once the source is removed

worked_example.json comes from extract_worked_example.py (one pass over the 3.5 GB record).
"""
import csv
import hashlib
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.patches import FancyBboxPatch, PathPatch, Rectangle  # noqa: E402
from matplotlib.path import Path                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "presentations"))
import make_deck as MD                                # noqa: E402

INK, DEEP, TEAL, SEAFOAM, MINT = "#10292E", "#05343D", "#028090", "#00A896", "#02C39A"
AMBER, AMBER_TXT, MUTED, LIGHT = "#E8871A", "#C96E0A", "#5B7179", "#EFF5F5"
SAND, MIST, GRID = "#F6D3A8", "#C9D6D9", "#D8E2E4"
PANEL = {"original": "#B07600", "kmeans": "#0072B2", "random": "#009E73"}
REG = MD.REG_LABEL
W_IN = 12.13                                          # slide width less the two margins

from matplotlib import font_manager  # noqa: E402
for _f in ("Regular", "Bold", "Italic", "BoldItalic"):
    _p = "/usr/share/fonts/truetype/crosextra/Carlito-%s.ttf" % _f
    if os.path.exists(_p):
        font_manager.fontManager.addfont(_p)
plt.rcParams.update({
    "font.family": "Carlito", "font.size": 13, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "figure.facecolor": "white", "axes.facecolor": "white", "axes.titleweight": "bold",
    "axes.titlesize": 14, "axes.titlelocation": "left", "axes.titlecolor": DEEP,
    "legend.frameon": False})


def save(fig, name):
    p = os.path.join(HERE, name)
    fig.savefig(p, dpi=220, facecolor="white")
    plt.close(fig)
    print("[figs] wrote", p)


def clean(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)


def box(ax, x, y, w, h, fc, *, ec="none", r=0.025, lw=1.0, z=1, transform=None):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=%g" % r, fc=fc,
                       ec=ec, lw=lw, zorder=z, transform=transform or ax.transAxes,
                       mutation_aspect=1)
    ax.add_patch(p)
    return p


def band(ax, x0, x1, y0a, y0b, y1a, y1b, color, alpha=0.55):
    """A Sankey ribbon from [y0a, y0b] at x0 to [y1a, y1b] at x1 (data coordinates)."""
    xm = (x0 + x1) / 2
    verts = [(x0, y0a), (xm, y0a), (xm, y1a), (x1, y1a), (x1, y1b), (xm, y1b), (xm, y0b),
             (x0, y0b), (x0, y0a)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4,
             Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), fc=color, ec="none", alpha=alpha, zorder=1))


def fmt(n):
    return f"{n:,}"


# ------------------------------------------------------------------ data, shared
_D = None


def deck_data():
    global _D
    if _D is None:
        MD.configure_depth(100, "union", None)
        _D = MD.load()
    return _D


def genus_outcomes():
    """Genus outcome of every amplicon x self-mode, classified exactly as the full deck's
    denominator slide (make_deck.s09_denominators) does."""
    dd = os.path.join(MD.AMP_RV, "data")
    truth = json.load(open(os.path.join(dd, "benchmark_truth.json")))
    named = lambda v: bool(v) and not str(v).startswith("unclassified")  # noqa: E731
    out = {}
    with open(os.path.join(dd, "asv_region_concordance.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["predictor"] != "anchor" or r["self_mode"] not in ("include", "exclude_0.987"):
                continue
            if not named(truth[r["amplicon"]]["truth_lineage"].get("Genus")):
                k = "unscoreable"
            elif r["abstained"] == "1":
                k = "abstained"
            else:
                k = {"1": "correct", "0": "wrong"}.get(r["correct_Genus"], "unnamed")
            out[(r["amplicon"], r["self_mode"])] = (r["region"], r["src_genome_id"], k)
    return out


OUTCOME = [("correct", "correct genus", TEAL), ("wrong", "wrong genus", AMBER),
           ("unnamed", "no named genus", SAND), ("abstained", "abstained", MIST),
           ("unscoreable", "truth has no genus", "#E4E9EA")]


# ------------------------------------------------------------------ 1. workflow
def short_org(name, genus):
    """'Halalkalicoccus sp. GCM10025704 VKM B-1617' -> 'H. sp. GCM10025704'."""
    t = name.replace(genus, genus[0] + ".").split()
    return " ".join(t[:3])


def fig_workflow():
    w = json.load(open(os.path.join(HERE, "worked_example.json")))
    D = deck_data()
    hits, src = w["hits"], w["source"]
    exc = w["selection"]["exclude_0.987"]
    reps = exc["reps"]
    tau = 0.987
    srcgen, srcphy = w["truth_lineage"]["Genus"], w["truth_lineage"]["Phylum"]
    fig = plt.figure(figsize=(W_IN, 5.05))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.2, 1], width_ratios=[1.12, 1.18, 0.74],
                          left=0.062, right=0.99, top=0.985, bottom=0.115, wspace=0.5,
                          hspace=0.16)
    raw = D["bvbrc_raw"] / 1e6
    STEPS = [("1  Align   2  Select",
              f"{raw:.2f}M records → 459,301 unique → 500 → 100 → ties"),
             ("3  Merge", "function union · p = k / n"),
             ("4  Score", "against the source genome")]
    # A: the 100-hit record, as the selector sees it with the source removed
    ax = fig.add_subplot(gs[1, 0])
    rk = np.array([h["rank"] for h in hits])
    idn = np.array([h["identity"] for h in hits])
    removed = np.array([(h["genome_id"] == src) or (h["identity"] >= tau) for h in hits])
    chosen = np.array([h["genome_id"] in reps for h in hits])
    def colour(h):
        g = h["genus"] or ""
        if g == srcgen or (h["organism"] or "").startswith(srcgen):
            return TEAL
        return SEAFOAM if h.get("phylum") == srcphy else AMBER
    col = np.array([colour(h) for h in hits])
    ax.axhspan(0.945, 0.987, color=LIGHT, zorder=0)
    for y, lab in ((0.987, "species  0.987"), (0.945, "genus  0.945")):
        ax.axhline(y, color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=1)
        ax.text(101, y + 0.0012, lab, fontsize=10, color=MUTED, ha="right", va="bottom")
    keep = ~removed & ~chosen
    ax.scatter(rk[keep], idn[keep], s=20, c=col[keep], lw=0, zorder=3)
    ax.scatter(rk[removed], idn[removed], s=46, marker="x", c=MUTED, lw=1.6, zorder=3)
    ax.scatter(rk[chosen], idn[chosen], s=26, c=TEAL, lw=0, zorder=4)
    ax.scatter([rk[chosen].mean()], [idn[chosen][0]], s=900, facecolors="none",
               edgecolors=DEEP, lw=1.6, zorder=4)
    ax.annotate(f"source + {int(removed.sum()) - 1} hits ≥ τ: removed",
                (rk[removed][-1] + 1, idn[removed][-1]), (14, 0.9965), fontsize=10.5,
                color=MUTED, va="center",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    ax.annotate(f"exact tie at {exc['best_identity']:.4f}: {len(reps)} genomes kept",
                (rk[chosen].max() + 2, idn[chosen][0]), (16, 0.975), fontsize=10.5,
                color=DEEP, weight="bold", va="center",
                arrowprops=dict(arrowstyle="-", color=DEEP, lw=0.9))
    for lab, c, y in ((srcgen, TEAL, 0.968), ("other " + srcphy, SEAFOAM, 0.9625),
                      ("other phyla", AMBER, 0.957)):
        ax.scatter([58], [y], s=30, c=c, lw=0)
        ax.text(60.5, y, lab, fontsize=10.5, va="center", color=INK)
    ax.set_xlim(0, 102)
    ax.set_ylim(0.918, 1.004)
    ax.set_xlabel("rank in the 100-hit record")
    ax.set_ylabel("identity to the amplicon")
    ax.set_title(f"A  a {REG.get(w['amplicon'].split('__')[1], 'V4')} amplicon of "
                 f"{srcgen}, source removed", fontsize=13)
    clean(ax)
    # B: the function matrix; each p group gets a readable width (true counts labelled)
    ax = fig.add_subplot(gs[1, 1])
    sets = [set(w["ec"][g]) for g in reps]
    truth = set(w["ec"][src])
    n = len(reps)
    union = sorted(set().union(*sets),
                   key=lambda e: (-sum(e in s for s in sets), e not in truth, e))
    k = np.array([sum(e in s for s in sets) for e in union])
    groups = [j for j in range(n, 0, -1) if (k == j).any()]
    share = {groups[0]: 0.46}
    rest = [j for j in groups[1:]]
    for j in rest:
        share[j] = 0.54 / len(rest)
    edges, x0 = [], 0.0
    for j in groups:
        m = int((k == j).sum())
        edges.extend(list(x0 + share[j] * np.arange(m) / m))
        x0 += share[j]
    edges.append(1.0)
    edges = np.array(edges)
    rows = [[e in s for e in union] for s in sets]
    Z = np.array(rows + [[False] * len(union)] + [[e in truth for e in union]], float)
    Z[n] = np.nan
    Z[n + 1] *= 2
    cmap = matplotlib.colors.ListedColormap(["#EDF2F2", TEAL, INK])
    ax.pcolormesh(edges, np.arange(n + 3), Z, cmap=cmap, vmin=0, vmax=2, shading="flat",
                  rasterized=True)
    ax.set_ylim(n + 2, -1.05)
    ax.set_xlim(0, 1)
    x0 = 0.0
    for j in groups:
        m = int((k == j).sum())
        ax.text(x0 + share[j] / 2, -0.62, "p = 1" if j == n else f"p = {j}/{n}", ha="center",
                fontsize=11.5, color=DEEP, weight="bold")
        ax.text(x0 + share[j] / 2, -0.18, f"{m} ECs", ha="center", fontsize=10.5,
                color=MUTED)
        x0 += share[j]
        if x0 < 1:
            ax.axvline(x0, color="white", lw=3)
    org = {h["genome_id"]: h["organism"] for h in hits}
    ax.set_yticks(list(np.arange(n) + 0.5) + [n + 1.5])
    ax.set_yticklabels([short_org(org[g], srcgen) for g in reps] + ["truth: the source"],
                       fontsize=10.5)
    ax.get_yticklabels()[-1].set_weight("bold")
    ax.set_xticks([])
    ax.tick_params(left=False)
    ax.set_title(f"B  synthetic genome: {len(union)} ECs the {n} genomes carry", fontsize=13,
                 pad=4)
    clean(ax, keep=())
    # C: how often an EC at each p is truly in the source
    ax = fig.add_subplot(gs[1, 2])
    present = [float(np.mean([e in truth for e, kk in zip(union, k) if kk == j]))
               for j in groups]
    yy = np.arange(len(groups))
    ax.barh(yy, [1] * len(groups), color="#EDF2F2", height=0.6)
    ax.barh(yy, present, color=TEAL, height=0.6)
    for i, (j, p) in enumerate(zip(groups, present)):
        ax.text(p - 0.025, i, f"{p:.2f}", ha="right", va="center", color="white",
                fontsize=12.5, weight="bold")
    ax.set_yticks(yy)
    ax.set_yticklabels(["p = 1" if j == n else f"p = {j}/{n}" for j in groups], fontsize=11.5)
    ax.tick_params(left=False)
    ax.set_xlim(0, 1)
    ax.set_ylim(len(groups) - 0.5, -0.5)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xlabel("share in the source genome")
    tp = len(set(union) & truth)
    prec, rec = tp / len(union), tp / len(truth)
    f1 = 2 * prec * rec / (prec + rec)
    ax.set_title(f"C  EC F1 {f1:.2f} vs truth", fontsize=13)
    clean(ax, keep=("bottom",))
    hdr = fig.add_subplot(gs[0, :])
    hdr.set_axis_off()
    pos = [a.get_position() for a in fig.axes if a is not hdr][:3]
    y0, y1 = hdr.get_position().y0, hdr.get_position().y1
    lefts = [0.008, pos[1].x0 - 0.115, pos[2].x0 - 0.06]
    rights = [pos[1].x0 - 0.128, pos[2].x0 - 0.073, 0.995]
    for (h, sub), xl, xr in zip(STEPS, lefts, rights):
        fig.patches.append(FancyBboxPatch((xl, y0), xr - xl, y1 - y0,
                                          boxstyle="round,pad=0,rounding_size=0.008",
                                          fc=LIGHT, ec="none", transform=fig.transFigure))
        fig.text(xl + 0.008, y0 + 0.66 * (y1 - y0), h, fontsize=14.5, weight="bold",
                 color=TEAL, va="center")
        fig.text(xl + 0.008, y0 + 0.27 * (y1 - y0), sub, fontsize=11, color=MUTED,
                 va="center")
    save(fig, "F_workflow.png")
    return dict(amplicon=w["amplicon"], f1=f1, precision=prec, recall=rec, n_union=len(union),
                n_truth=len(truth), n_removed=int(removed.sum()),
                n_by_p={f"{j}/{n}": int((k == j).sum()) for j in groups},
                present_by_p={f"{j}/{n}": p for j, p in zip(groups, present)},
                best_identity=exc["best_identity"], reps=reps,
                hits_other_phyla=int(sum(c == AMBER for c in col)))


# ------------------------------------------------------------------ 2. validation arms
def fig_arms():
    D = deck_data()
    tax, meta = D["tax"]["regions"], D["tax"]["_meta"]
    regs = list(tax)
    g = D["kg"]
    mg = g["meta"]
    fp = json.load(open(os.path.join(MD.SR, "data", "fold_plan.json")))
    drafts = D["ret"]["clean"]
    n_drafts = sum(r["n"] for r in drafts)
    p2 = D["cmp"]["_meta"]["picrust2"]
    n_amp_reg = {}
    with open(os.path.join(MD.AMP_RV, "data", "asv_region_concordance.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["self_mode"] == "include" and r["predictor"] == "anchor":
                n_amp_reg[r["region"]] = n_amp_reg.get(r["region"], 0) + 1
    assert sum(n_amp_reg.values()) == meta["n_total"]
    dec = json.load(open(os.path.join(REPO, "three_panel", "space_deciles_three.json")))
    arms = [
        ("In-silico amplicons", TEAL,
         [("query", f"{fmt(meta['n_total'])} amplicons, {len(regs)} regions"),
          ("reference", "all 459,301 BV-BRC 16S"),
          ("held out", "source + every hit ≥ 0.987"),
          ("truth", "the source's lineage, PGFams")]),
        ("vs PICRUSt2, as shipped", AMBER,
         [("query", "the same amplicons"),
          ("reference", f"each its own: {fmt(p2['n_reference_genomes'])} vs ~315k"),
          ("held out", "prFBA only (PICRUSt2 cannot)"),
          ("truth", "BV-BRC ECs: prFBA's vocabulary")]),
        ("Controlled, shared reference", DEEP,
         [("query", f"{fmt(mg['n_eligible'])} panel genomes' 16S"),
          ("reference", f"one for both: {fmt(g['headline']['reference_density_genomes'])} (median)"),
          ("held out", f"whole genera, {mg['n_folds']} folds"),
          ("truth", "one answer key + ATP, gapfill")]),
        ("Synthetic MAGs", SEAFOAM,
         [("query", f"{fmt(n_drafts)} drafts, {len(drafts)} levels"),
          ("reference", "BV-BRC + GTDB, by skani ANI"),
          ("held out", "source + every hit ≥ 95 ANI"),
          ("truth", "the source genome")]),
        ("Three panels", PANEL["kmeans"],
         [("query", "3 × 10,000 organisms"),
          ("reference", "the same database"),
          ("held out", "source + every hit ≥ 0.987"),
          ("truth", "is coverage a panel artefact?")]),
    ]
    fig = plt.figure(figsize=(W_IN, 5.0))
    n = len(arms)
    gap, x0 = 0.012, 0.0
    cw = (1 - gap * (n - 1)) / n
    for i, (name, c, rows) in enumerate(arms):
        x = x0 + i * (cw + gap)
        fig.patches.append(FancyBboxPatch((x, 0.0), cw, 0.995, boxstyle="round,pad=0,rounding_size=0.01",
                                          fc=LIGHT, ec="none", transform=fig.transFigure, zorder=-2))
        fig.patches.append(FancyBboxPatch((x, 0.89), cw, 0.105, boxstyle="round,pad=0,rounding_size=0.01",
                                          fc=c, ec="none", transform=fig.transFigure, zorder=-1))
        fig.text(x + cw / 2, 0.942, name, ha="center", va="center", color="white",
                 fontsize=13.5, weight="bold")
        inset = 0.062 if i in (0, 1, 2, 4) else 0.015
        ax = fig.add_axes([x + inset, 0.5, cw - inset - 0.012, 0.32])
        ax.set_zorder(2)
        if i == 0:
            vals = [n_amp_reg[r] for r in regs]
            ax.bar(range(len(regs)), vals, color=TEAL, width=0.7)
            ax.set_xticks(range(len(regs)))
            ax.set_xticklabels([REG[r].replace("Full-length", "FL").replace(" 944R", "*")
                                for r in regs], fontsize=8.5, rotation=90)
            ax.set_ylabel("amplicons", fontsize=10)
            ax.tick_params(axis="y", labelsize=9)
        elif i == 1:
            ax.barh([1, 0], [p2["n_reference_genomes"], 315000], color=[AMBER, TEAL], height=0.6)
            ax.set_yticks([1, 0])
            ax.set_yticklabels(["PICRUSt2", "prFBA"], fontsize=10.5)
            ax.set_xlabel("reference genomes", fontsize=10)
            ax.set_xticks([0, 150000, 300000])
            ax.set_xticklabels(["0", "150k", "300k"], fontsize=9)
            ax.set_ylim(-0.6, 1.6)
        elif i == 2:
            sizes = {}
            for gg in fp["genome_group"].values():
                sizes[gg] = sizes.get(gg, 0) + 1
            for k, f in enumerate(fp["folds"]):
                ss = sorted((sizes[t] for t in f["test_groups"]), reverse=True)
                left = 0
                for j, v in enumerate(ss):
                    ax.barh(k, v, left=left, color=DEEP if j % 2 == 0 else "#2F6B75",
                            height=0.7, lw=0)
                    left += v
            ax.set_yticks(range(len(fp["folds"])))
            ax.set_yticklabels([f"fold {k + 1}" for k in range(len(fp["folds"]))], fontsize=9.5)
            ax.invert_yaxis()
            ax.set_xlabel("held-out genomes; stripes = genera", fontsize=9.5)
            ax.tick_params(axis="x", labelsize=9)
        elif i == 3:
            rng_ = np.random.default_rng(7)
            for k, cpl in enumerate([1.0, 0.7, 0.4, 0.1]):
                pieces = rng_.integers(2, 6, 18)
                keep = rng_.random(18) < cpl if cpl < 1 else np.ones(18, bool)
                left = 0
                for L, kk in zip(pieces, keep):
                    ax.barh(k, L - 0.25, left=left, color=SEAFOAM if kk else "#DDE6E7",
                            height=0.55, lw=0)
                    left += L
                ax.text(-1, k, f"{cpl:.0%}", ha="right", va="center", fontsize=10)
            ax.set_xlim(-9, left + 1)
            ax.invert_yaxis()
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel("fragmented, then thinned (schematic)", fontsize=9.5)
            clean(ax, keep=())
        else:
            xs = np.arange(1, 11)
            for key, colr, lab in (("og", PANEL["original"], "stratified"),
                                   ("km", PANEL["kmeans"], "k-means"),
                                   ("rnd", PANEL["random"], "random")):
                ax.plot(xs, 100 * np.array(dec[key]), color=colr, lw=2, label=lab)
            ax.plot(xs, 100 * np.array(dec["population"]), color=INK, lw=1.2, ls=":",
                    label="database")
            ax.set_xlabel("density decile (1 = sparsest)", fontsize=9.5)
            ax.set_ylabel("% of picks", fontsize=10)
            ax.tick_params(labelsize=9)
            ax.legend(fontsize=8.5, loc="upper right", handlelength=1.2)
        ax.set_facecolor(LIGHT)
        if i != 3:
            clean(ax)
        for j, (lab, txt) in enumerate(rows):
            yy = 0.33 - j * 0.083
            fig.text(x + 0.012, yy, lab.upper(), fontsize=9, color=MUTED, weight="bold",
                     va="center")
            fig.text(x + 0.012, yy - 0.033, txt, fontsize=11, color=INK, va="center")
    save(fig, "F_arms.png")
    return dict(n_drafts=n_drafts, n_amp_reg=n_amp_reg)


# ------------------------------------------------------------------ 3. sankey
def fig_sankey():
    D = deck_data()
    regs = list(D["tax"]["regions"])
    out = genus_outcomes()
    cnt = {m: {} for m in ("include", "exclude_0.987")}
    for (amp, mode), (reg, _src, k) in out.items():
        cnt[mode][(reg, k)] = cnt[mode].get((reg, k), 0) + 1
    n_tot = sum(cnt["include"].values())
    assert n_tot == D["tax"]["_meta"]["n_total"]
    # the denominator slide's pooled split must come out identically
    ex = {k: sum(v for (r, kk), v in cnt["exclude_0.987"].items() if kk == k) for k, _, _ in OUTCOME}
    scoreable = n_tot - ex["unscoreable"]
    gsum = sum(D["tax"]["regions"][r]["anchor"]["exclude_0.987"]["all"]["all"]["per_rank"]["Genus"]["n"]
               for r in regs)
    assert scoreable == gsum, (scoreable, gsum)
    fig, ax = plt.subplots(figsize=(W_IN, 5.0))
    fig.subplots_adjust(left=0.0, right=1.0, top=0.93, bottom=0.0)
    ax.set_xlim(-0.02, 1.02)
    gapr, gapo = 0.012, 0.03
    H = 1.0 - gapr * (len(regs) - 1)
    scale = H / n_tot
    # region nodes in the middle
    ry = {}
    y = 1.0
    for r in regs:
        h = sum(v for (rr, k), v in cnt["include"].items() if rr == r) * scale
        ry[r] = [y - h, y]
        y -= h + gapr
    def onodes(mode):
        H2 = 1.0 - gapo * (len(OUTCOME) - 1)
        s2 = H2 / n_tot
        nodes, y = {}, 1.0
        for k, _, _ in OUTCOME:
            h = sum(v for (r, kk), v in cnt[mode].items() if kk == k) * s2
            nodes[k] = [y - h, y]
            y -= h + gapo
        return nodes, s2
    xr0, xr1 = 0.462, 0.538
    for mode, side in (("include", -1), ("exclude_0.987", 1)):
        nodes, s2 = onodes(mode)
        xo0, xo1 = (0.11, 0.135) if side < 0 else (0.865, 0.89)
        rcur = {r: ry[r][1] for r in regs}
        ocur = {k: nodes[k][1] for k in nodes}
        # each region's slice of the region node is split half to each side: draw at
        # full height on both sides (both modes cover every amplicon once)
        for r in regs:
            rc = ry[r][1]
            for k, _, colr in OUTCOME:
                v = cnt[mode].get((r, k), 0)
                if not v:
                    continue
                hr, ho = v * scale, v * s2
                if side < 0:
                    band(ax, xo1, xr0, ocur[k] - ho, ocur[k], rc - hr, rc, colr, 0.6)
                else:
                    band(ax, xr1, xo0, rc - hr, rc, ocur[k] - ho, ocur[k], colr, 0.6)
                rc -= hr
                ocur[k] -= ho
        tot = {k: sum(v for (r, kk), v in cnt[mode].items() if kk == k) for k in nodes}
        ly, last = {}, None
        for k, _, _ in OUTCOME:                     # label centres, at least 0.085 apart
            y0, y1 = nodes[k]
            if y1 - y0 <= 0:
                continue
            c = (y0 + y1) / 2
            if last is not None and last - c < 0.085:
                c = last - 0.085
            ly[k] = last = c
        for k, lab, colr in OUTCOME:
            y0, y1 = nodes[k]
            if y1 - y0 <= 0:
                continue
            ax.add_patch(Rectangle((xo0, y0), xo1 - xo0, y1 - y0, fc=colr, ec="none", zorder=3))
            sc = tot[k] / (n_tot - tot["unscoreable"]) if k != "unscoreable" else tot[k] / n_tot
            txt = f"{lab}\n{fmt(tot[k])} · {sc:.0%}" + ("" if k != "unscoreable" else " of all")
            if side < 0:
                ax.text(xo0 - 0.008, ly[k], txt, ha="right", va="center", fontsize=11)
            else:
                ax.text(xo1 + 0.008, ly[k], txt, ha="left", va="center", fontsize=11)
    for r in regs:
        y0, y1 = ry[r]
        ax.add_patch(Rectangle((xr0, y0), xr1 - xr0, y1 - y0, fc=DEEP, ec="none", zorder=3))
        ax.text(0.5, (y0 + y1) / 2, REG[r], ha="center", va="center", color="white",
                fontsize=10.5, weight="bold", zorder=4)
    ax.text(0.123, 1.035, "SOURCE PRESENT", ha="center", fontsize=12.5, weight="bold", color=TEAL)
    ax.text(0.5, 1.035, f"{fmt(n_tot)} AMPLICONS", ha="center", fontsize=12.5, weight="bold",
            color=DEEP)
    ax.text(0.877, 1.035, "SOURCE REMOVED (τ 0.987)", ha="center", fontsize=12.5,
            weight="bold", color=AMBER_TXT)
    ax.set_ylim(-0.01, 1.07)
    ax.set_axis_off()
    save(fig, "F_sankey.png")
    inc = {k: sum(v for (r, kk), v in cnt["include"].items() if kk == k) for k, _, _ in OUTCOME}
    # abstentions whose whole record was stripped by self-exclusion (n_removed_self >= depth)
    emptied = 0
    with open(os.path.join(MD.AMP_RV, "data", "asv_region_concordance.csv")) as fh:
        for r in csv.DictReader(fh):
            if (r["predictor"] == "anchor" and r["self_mode"] == "exclude_0.987"
                    and out.get((r["amplicon"], "exclude_0.987"), ("", "", ""))[2] == "abstained"
                    and int(r["n_removed_self"] or 0) >= MD.DEPTH):
                emptied += 1
    return dict(n=n_tot, scoreable=scoreable, include=inc, exclude=ex, abstain_emptied=emptied)


# ------------------------------------------------------------------ 4. PICRUSt2
def fig_picrust2():
    D = deck_data()
    p2 = D["cmp"]["_meta"]["picrust2"]
    m = D["cmp"]["_meta"]
    g = D["kg"]
    fz = g["meta"]["frozen"]["picrust2"]
    ref = g["headline"]["reference_density_genomes"]
    fig, ax = plt.subplots(figsize=(W_IN, 5.0))
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    def node(x, y, w, h, head, sub, fc, tc="white", sc=None, hs=13):
        box(ax, x, y, w, h, fc, r=0.012, transform=ax.transData)
        ax.text(x + w / 2, y + h * (0.64 if sub else 0.5), head, ha="center", va="center",
                fontsize=hs, weight="bold", color=tc)
        if sub:
            ax.text(x + w / 2, y + h * 0.3, sub, ha="center", va="center", fontsize=10.5,
                    color=sc or tc, linespacing=1.15)

    def arrow(x0, y0, x1, y1, c=MUTED):
        ax.annotate("", (x1, y1), (x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=1.4, mutation_scale=14))

    for lane, y, lab, sub in (("A", 0.545, "A  As shipped", "each system, its own reference"),
                              ("B", 0.045, "B  Controlled", "one reference, rebuilt per fold")):
        box(ax, 0.0, y - 0.02, 1.0, 0.44, LIGHT, r=0.012, transform=ax.transData)
        ax.text(0.012, y + 0.385, lab, fontsize=14, weight="bold", color=DEEP, va="center")
        ax.text(0.012, y + 0.335, sub, fontsize=11, color=MUTED, va="center")
    # lane A
    yA = 0.56
    node(0.012, yA + 0.08, 0.12, 0.17, "Amplicon", f"{fmt(m['n_amplicons'])}\nin-silico", INK, hs=12.5)
    node(0.17, yA + 0.165, 0.43, 0.12, "PICRUSt2 " + p2["version"],
         f"HMMER align → EPA-ng place → castor {p2['hsp'].upper()} over "
         f"{fmt(p2['n_reference_genomes'])} genomes", AMBER, hs=12.5)
    node(0.17, yA + 0.015, 0.43, 0.12, "prFBA",
         "align to 459,301 BV-BRC 16S → exact ties → PGFam union", TEAL, hs=12.5)
    node(0.64, yA + 0.08, 0.15, 0.17, "EC map", f"PGFam → product → EC\n{fmt(m['ec_universe_shared'])} shared ECs",
         "white", tc=DEEP, sc=MUTED, hs=12.5)
    node(0.83, yA + 0.08, 0.16, 0.17, "BV-BRC truth", "prFBA's own\nvocabulary", "white",
         tc=AMBER_TXT, sc=MUTED, hs=12.5)
    arrow(0.132, yA + 0.2, 0.17, yA + 0.225)
    arrow(0.132, yA + 0.13, 0.17, yA + 0.075)
    arrow(0.6, yA + 0.225, 0.64, yA + 0.19)
    arrow(0.6, yA + 0.075, 0.64, yA + 0.13)
    arrow(0.79, yA + 0.165, 0.83, yA + 0.165)
    # lane B
    yB = 0.06
    node(0.012, yB + 0.08, 0.12, 0.17, f"{fmt(g['meta']['n_eligible'])}", "panel genomes,\ngenus-blocked",
         INK, hs=12.5)
    node(0.165, yB + 0.08, 0.155, 0.17, f"{g['meta']['n_folds']} folds",
         "train genera →\nreference", DEEP, hs=12.5)
    node(0.355, yB + 0.08, 0.15, 0.17, f"{fmt(ref)} genomes", "one species-\ndereplicated set", DEEP, hs=12.5)
    node(0.54, yB + 0.165, 0.26, 0.12, "PICRUSt2, rebuilt",
         f"FastTree {' '.join(g['meta']['frozen']['fasttree_args'])} · HSP {fz['hsp_method']}",
         AMBER, hs=12.5)
    node(0.54, yB + 0.015, 0.26, 0.12, "prFBA, same genomes", "nested CV picks the rule", TEAL,
         hs=12.5)
    node(0.83, yB + 0.08, 0.16, 0.17, "One answer key", "held-out genera\n+ ATP, gapfill", "white",
         tc=DEEP, sc=MUTED, hs=12.5)
    for a, b in ((0.132, 0.165), (0.32, 0.355)):
        arrow(a, yB + 0.165, b, yB + 0.165)
    arrow(0.505, yB + 0.19, 0.54, yB + 0.225)
    arrow(0.505, yB + 0.14, 0.54, yB + 0.075)
    arrow(0.8, yB + 0.225, 0.83, yB + 0.19)
    arrow(0.8, yB + 0.075, 0.83, yB + 0.13)
    save(fig, "F_picrust2.png")
    return dict(ref=ref, fasttree=g["meta"]["frozen"]["fasttree_args"], hsp=fz["hsp_method"])


# ------------------------------------------------------------------ 5. heatmap
def fig_heatmap():
    D = deck_data()
    tax = D["tax"]["regions"]
    cap = D["cap"]["_meta"]
    regs = list(tax)
    pr = lambda r, m, rank: tax[r]["anchor"][m]["all"]["all"]["per_rank"][rank]  # noqa: E731
    cols = [("Genus", "present", lambda r: pr(r, "include", "Genus")["correct"]),
            ("Genus", "removed\n÷ answered", lambda r: pr(r, "exclude_0.987", "Genus")["correct"]),
            ("Genus", "removed\n÷ all", lambda r: pr(r, "exclude_0.987", "Genus")["correct_population"]),
            ("Family", "present", lambda r: pr(r, "include", "Family")["correct"]),
            ("Family", "removed\n÷ answered", lambda r: pr(r, "exclude_0.987", "Family")["correct"]),
            ("Family", "removed\n÷ all", lambda r: pr(r, "exclude_0.987", "Family")["correct_population"]),
            ("Gene capture", "present", lambda r: cap["headline_macro_recall_pct"]["include"][r] / 100),
            ("Gene capture", "removed\n÷ answered", lambda r: cap["headline_macro_recall_pct"]["exclude_0.987"][r] / 100),
            ("Gene capture", "removed\n÷ all", lambda r: cap["headline_macro_recall_pct_population"]["exclude_0.987"][r] / 100)]
    Z = np.array([[f(r) for _, _, f in cols] for r in regs])
    fig, ax = plt.subplots(figsize=(W_IN, 4.95))
    fig.subplots_adjust(left=0.1, right=0.995, top=0.8, bottom=0.02)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("t", ["#FFFFFF", "#BFE3E0", TEAL, DEEP])
    xe = []
    x = 0
    for j in range(len(cols)):
        xe.append(x)
        x += 1 + (0.25 if j in (2, 5) else 0)
    for i, r in enumerate(regs):
        for j, xj in enumerate(xe):
            v = Z[i, j]
            ax.add_patch(Rectangle((xj, i), 0.96, 0.92, fc=cmap(v), ec="none"))
            ax.text(xj + 0.48, i + 0.46, f"{v:.2f}", ha="center", va="center", fontsize=12,
                    color="white" if v > 0.55 else INK, weight="bold" if j % 3 == 2 else None)
    ax.set_xlim(0, x)
    ax.set_ylim(len(regs), 0)
    ax.set_yticks(np.arange(len(regs)) + 0.46)
    ax.set_yticklabels([REG[r] for r in regs], fontsize=12)
    ax.set_xticks([xj + 0.48 for xj in xe])
    ax.set_xticklabels([c[1] for c in cols], fontsize=10.5)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for k, nm in enumerate(("Genus accuracy", "Family accuracy", "Gene capture (macro recall)")):
        ax.text(xe[3 * k] + 1.5, -1.25, nm, ha="center", fontsize=13.5, weight="bold",
                color=DEEP)
    clean(ax, keep=())
    save(fig, "F_heatmap.png")
    return dict(cols=[c[0] + " " + c[1].replace("\n", " ") for c in cols],
                min=Z.min(0).tolist(), max=Z.max(0).tolist())


# ------------------------------------------------------------------ 6. 16S space
def fig_space():
    import pandas as pd
    D = deck_data()
    TP = os.path.join(REPO, "two_panel")
    sp = np.load(TP + "/figdata_space.npz")
    kd = np.load(TP + "/knn_density.npz")
    clust = sp["clusterable"]
    xy = sp["xy"][clust]
    dec = kd["decile"].astype(int)
    n = len(xy)
    pos = np.full(len(clust), -1, dtype=np.int64)
    pos[np.flatnonzero(clust)] = np.arange(n)
    df = pd.read_parquet(REPO + "/index.parquet")
    sub = df[(df.seq_len >= 1400) & (~df.atypical_len)].reset_index(drop=True)
    assert len(sub) == len(clust)
    md5_to_row = {m: i for i, m in enumerate(sub.md5.astype(str).to_numpy())}
    raw = json.load(open("/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_seq.json"))
    raw2up = {k: hashlib.md5(v.upper().encode()).hexdigest() for k, v in raw.items()}
    del raw
    where = {}
    for gg in json.load(open(os.path.join(MD.AMP_RV, "data", "selection_10k.json")))["genomes"]:
        m = str(gg["best_16s_md5"])
        i = md5_to_row.get(m, md5_to_row.get(raw2up.get(m, ""), None))
        if i is not None and clust[i]:
            where[gg["genome_id"]] = pos[i]
    out = genus_outcomes()
    order = np.argsort(np.argsort(kd["knn21"]))
    dens_pct = 1.0 - order / (n - 1)
    fig = plt.figure(figsize=(W_IN, 4.95))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], left=0.0, right=0.99, top=0.9,
                          bottom=0.13, wspace=0.08)
    colmap = {k: c for k, _, c in OUTCOME}
    stats = {}
    for j, reg in enumerate(("FullLength16S", "V4")):
        ax = fig.add_subplot(gs[0, j])
        ax.hexbin(xy[:, 0], xy[:, 1], C=dens_pct, reduce_C_function=np.mean, gridsize=150,
                  cmap="Greys", vmin=-0.4, vmax=1.6, linewidths=0, rasterized=True)
        pts = [(where[src], k) for (amp, mode), (r, src, k) in out.items()
               if mode == "exclude_0.987" and r == reg and src in where and k != "unscoreable"]
        perm = np.random.default_rng(1729).permutation(len(pts))  # no outcome drawn on top
        idx = np.array([pts[i][0] for i in perm], dtype=int)
        ax.scatter(xy[idx, 0], xy[idx, 1], s=2.4, c=[colmap[pts[i][1]] for i in perm], lw=0,
                   rasterized=True, alpha=0.95)
        ax.set_title(f"{'ABC'[j]}  {REG[reg]}: genus call, source removed", fontsize=13)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal"); ax.set_anchor("N")
        clean(ax, keep=())
        stats[reg] = len(pts)
    # C: by density decile
    ax = fig.add_subplot(gs[0, 2])
    lab = 10 - dec
    xs = np.arange(1, 11)
    for reg, ls in (("FullLength16S", "-"), ("V4", "--")):
        per = {d: {"correct": 0, "n": 0, "abst": 0} for d in xs}
        for (amp, mode), (r, src, k) in out.items():
            if mode != "exclude_0.987" or r != reg or src not in where or k == "unscoreable":
                continue
            d = lab[where[src]]
            per[d]["n"] += 1
            per[d]["correct"] += k == "correct"
            per[d]["abst"] += k == "abstained"
        acc = [per[d]["correct"] / max(per[d]["n"], 1) for d in xs]
        ab = [per[d]["abst"] / max(per[d]["n"], 1) for d in xs]
        ax.plot(xs, acc, color=TEAL, lw=2.4, ls=ls, marker="o", ms=4,
                label=f"correct genus, {REG[reg]}")
        ax.plot(xs, ab, color=MUTED, lw=2.0, ls=ls, marker="o", ms=4,
                label=f"abstained, {REG[reg]}")
        stats[reg + "_acc_by_decile"] = acc
        stats[reg + "_abst_by_decile"] = ab
        stats[reg + "_n_by_decile"] = [per[d]["n"] for d in xs]
    ax.set_xticks(xs)
    ax.set_xlabel("local density decile of the source (1 = sparsest)")
    ax.set_ylabel("share of scoreable amplicons")
    ax.set_ylim(0, 1)
    ax.yaxis.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
    ax.legend(fontsize=10, loc="upper left", ncol=2, handlelength=2.2)
    ax.set_title("C  densest space: records empty, prFBA abstains", fontsize=13)
    clean(ax)
    # outcome legend under A and B
    hs = [plt.Line2D([], [], marker="o", ls="", color=colmap[k], ms=7) for k in
          ("correct", "wrong", "unnamed", "abstained")]
    fig.legend(hs, ["correct genus", "wrong genus", "no named genus", "abstained"], ncol=4,
               loc="lower left", bbox_to_anchor=(0.03, 0.0), fontsize=11.5)
    save(fig, "F_space.png")
    stats["n_located"] = len(where)
    return stats


FIGS = {"workflow": fig_workflow, "arms": fig_arms, "sankey": fig_sankey,
        "picrust2": fig_picrust2, "heatmap": fig_heatmap, "space": fig_space}

if __name__ == "__main__":
    want = sys.argv[1:] or list(FIGS)
    stats = {}
    sp = os.path.join(HERE, "figure_stats.json")
    if os.path.exists(sp):
        stats = json.load(open(sp))
    for nm in want:
        stats[nm] = FIGS[nm]()
    json.dump(stats, open(sp, "w"), indent=1, sort_keys=True, default=float)
