#!/usr/bin/env python
"""strategy_tradeoff.py — genus-level gene-capture vs F1/MCC tradeoff across the
selection strategies, for the most prominent 16S regions.

The question: how does the choice of selection strategy (keep only the best-hit
genome -> keep all best-hitting genomes -> band unions -> all reliable hits ->
legacy reducer) trade **% of the source organism's genes captured** (recall)
against the **precision-balanced quality** of the predicted gene set (F1 and MCC)?

This is only meaningful in the SELF-EXCLUSION (novel-organism) scenario: with the
source organism in the database every strategy captures ~100% (best-hit = the
source), so the strategies only separate once the source + its species-level
neighbours are removed and the organism must be represented by its **genus-level
relatives**. We therefore (a) drop the source organism and every >=tau (0.987)
near-identical hit, (b) keep the GENUS-tier amplicons (best remaining relative in
[0.945, 0.987)), and (c) score each strategy's unioned PGFam set against the
removed source's PGFam set.

Metrics per (genome, region, strategy) cell, macro-averaged over the genus-tier
amplicons of each region:
  recall    = i/t                      (= "% genes captured")
  precision = i/p
  F1        = 2i/(t+p)
  MCC       = (TP*TN - FP*FN)/sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
              over the global PGFam universe U (TP=i, FP=p-i, FN=t-i, TN=U-(t+p-i)).

Parallelism: the PGFam cache (~112k genomes) is loaded + converted ONCE in the
parent and shared copy-on-write with a fork ProcessPoolExecutor; the
(strategy x region) cells run concurrently. Each cell is pure (re-select +
set ops), so wall-clock is ~one cache load plus the slowest cell.

Outputs (dedicated folder region_validation/figures_strategy_tradeoff/):
  strategy_tradeoff.csv      per (region, strategy) recall/precision/F1/MCC + n
  strategy_tradeoff.png/.pdf per-region recall-vs-F1 and recall-vs-MCC frontier

Interpreter: /home/freiburger/Documents/py_venv/bin/python
"""
from __future__ import annotations

import csv
import json
import math
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import _common as C  # noqa: E402
from _common import PATHS  # noqa: E402
import select_references as SR  # noqa: E402
import score_genecap as GC  # noqa: E402

# Most prominent survey regions (+ full-length control), short->long.
REGIONS = ["V4", "V3-V4", "V4-V5", "FullLength16S"]

TAU = 0.987          # self-exclusion threshold (species floor); source + >=tau dropped
GENUS_LO, GENUS_HI = 0.945, 0.987   # genus band for the genus-tier stratum

# Strategies, ordered by increasing union breadth (the tradeoff axis); "legacy" is
# the off-axis reference (the ordered reducer pipeline).
STRATEGIES = [
    ("best hit only",     dict(select_all=True, select_all_window=0.0,  select_all_species_dedup=True, select_all_cap=1)),
    ("all best-hitting",  dict(select_all=True, select_all_window=0.0,  select_all_species_dedup=True)),
    ("band union 0.005",  dict(select_all=True, select_all_window=0.005, select_all_species_dedup=True)),
    ("band union 0.02",   dict(select_all=True, select_all_window=0.02,  select_all_species_dedup=True)),
    ("all reliable",      dict(select_all=True, select_all_window=1.0,   select_all_species_dedup=True)),
    ("legacy reducer",    dict(select_all=False)),
]
BREADTH_ORDER = ["best hit only", "all best-hitting", "band union 0.005",
                 "band union 0.02", "all reliable"]   # connected frontier; legacy is standalone

# ----- globals shared via fork -----
CACHE: dict = {}
GENUS: dict = {}     # region -> list of (key, src_gid, excl_record)
U = 1


def _gp(gid):
    """gene_provider over the in-memory cache (no network); empty set if absent."""
    return CACHE.get(str(gid)) or frozenset()


def score_cell(strategy_name, knobs, region):
    """Macro-mean recall/precision/F1/MCC over the region's genus-tier amplicons."""
    rec_s = prec_s = f1_s = mcc_s = 0.0
    n_prec = n = 0
    for key, src_gid, excl_rec in GENUS[region]:
        sel = SR.select_representatives(excl_rec, knobs=knobs, gene_provider=_gp)
        P = set()
        for s in (sel.get("selected") or []):
            P |= _gp(s["genome_id"])
        T = _gp(src_gid)
        t, p = len(T), len(P)
        if t == 0:
            continue
        n += 1
        i = len(T & P)
        rec_s += i / t
        f1_s += (2 * i / (t + p)) if (t + p) else 0.0
        if p:
            prec_s += i / p
            n_prec += 1
        TP, FP, FN = i, p - i, t - i
        TN = U - (t + p - i)
        d2 = (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
        mcc_s += ((TP * TN - FP * FN) / math.sqrt(d2)) if d2 > 0 else 0.0
    if n == 0:
        return strategy_name, region, None
    return strategy_name, region, dict(
        recall=rec_s / n, precision=(prec_s / n_prec if n_prec else float("nan")),
        f1=f1_s / n, mcc=mcc_s / n, n=n)


def build_globals():
    global CACHE, GENUS, U
    t0 = time.perf_counter()
    print("[tradeoff] loading PGFam cache ...", flush=True)
    raw = json.load(open(PATHS.pgfam_cache))
    CACHE = {g: frozenset(v) for g, v in raw.items()}
    uni = set()
    for v in CACHE.values():
        uni |= v
    U = max(len(uni), 1)
    del raw, uni
    print(f"[tradeoff] cache={len(CACHE)} genomes, universe U={U} ({time.perf_counter()-t0:.0f}s)", flush=True)

    print("[tradeoff] loading hits + truth, building genus-tier exclude sets ...", flush=True)
    hits = json.load(open(PATHS.hits_json))
    truth = json.load(open(PATHS.truth_json))
    GENUS = {r: [] for r in REGIONS}
    for key, t in truth.items():
        region = key.rpartition("__")[2]
        if region not in GENUS:
            continue
        src_gid = str(t["src_genome_id"])
        if src_gid not in CACHE or not CACHE[src_gid]:
            continue                      # need the source's true gene set
        src_tid = GC._src_taxon_id(src_gid)
        hr = hits.get(key)
        if not hr:
            continue
        excl_rec, _ = GC.exclude_self_record(hr, src_gid, src_tid, TAU)
        b = excl_rec.get("best_identity") or 0.0
        if GENUS_LO <= b < GENUS_HI:      # genus-tier relative remains
            GENUS[region].append((key, src_gid, excl_rec))
    for r in REGIONS:
        print(f"[tradeoff]   {r}: {len(GENUS[r])} genus-tier amplicons", flush=True)


def main():
    outdir = os.path.join(C.BASE, "figures_strategy_tradeoff")
    if "--replot" in sys.argv:
        render(load_results_csv(outdir), outdir)
        print(f"[tradeoff] re-rendered figures from CSV -> {outdir}")
        return
    build_globals()
    cells = [(name, knobs, r) for (name, knobs) in STRATEGIES for r in REGIONS]
    print(f"[tradeoff] scoring {len(cells)} (strategy x region) cells in parallel ...", flush=True)
    results = {}
    ctx = mp.get_context("fork")           # inherit CACHE/GENUS/U copy-on-write
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=min(len(cells), (os.cpu_count() or 8)), mp_context=ctx) as ex:
        futs = [ex.submit(score_cell, n, k, r) for (n, k, r) in cells]
        for fut in as_completed(futs):
            name, region, m = fut.result()
            if m:
                results[(name, region)] = m
    print(f"[tradeoff] done in {time.perf_counter()-t0:.0f}s", flush=True)

    outdir = os.path.join(C.BASE, "figures_strategy_tradeoff")
    os.makedirs(outdir, exist_ok=True)
    # ---- CSV ----
    with open(os.path.join(outdir, "strategy_tradeoff.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["region", "strategy", "n_amplicons", "pct_genes_captured", "precision", "f1", "mcc"])
        for (name, _) in STRATEGIES:
            for r in REGIONS:
                m = results.get((name, r))
                if m:
                    w.writerow([r, name, m["n"], round(100*m["recall"], 3),
                                round(m["precision"], 4), round(m["f1"], 4), round(m["mcc"], 4)])
    render(results, outdir)
    print(f"[tradeoff] wrote {outdir}/strategy_tradeoff.{{csv,png,pdf}}")


SHORT = {"best hit only": "best-hit", "all best-hitting": "all-best",
         "band union 0.005": "band .005", "band union 0.02": "band .02",
         "all reliable": "all-reliable", "legacy reducer": "legacy"}


def load_results_csv(outdir):
    """Reload the per-(region,strategy) metrics from the CSV (for fast --replot)."""
    res = {}
    with open(os.path.join(outdir, "strategy_tradeoff.csv")) as fh:
        for row in csv.DictReader(fh):
            res[(row["strategy"], row["region"])] = dict(
                recall=float(row["pct_genes_captured"]) / 100.0,
                precision=float(row["precision"]), f1=float(row["f1"]),
                mcc=float(row["mcc"]), n=int(row["n_amplicons"]))
    return res


def render(results, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    f1c, mccc, prc = "#1f77b4", "#d62728", "#7f7f7f"
    # shared x-range so the frontier shapes are directly comparable across regions
    allx = [100 * results[(n, r)]["recall"] for (n, _knobs) in STRATEGIES for r in REGIONS
            if (n, r) in results]
    xlo, xhi = min(allx) - 2, max(allx) + 4

    ncol = 2
    nrow = (len(REGIONS) + 1) // 2
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.8 * ncol, 5.8 * nrow), squeeze=False)
    for ax, region in zip(axes.flat, REGIONS):
        present = [n for n in BREADTH_ORDER if (n, region) in results]
        xs = [100 * results[(n, region)]["recall"] for n in present]
        f1s = [results[(n, region)]["f1"] for n in present]
        mcs = [results[(n, region)]["mcc"] for n in present]
        prs = [results[(n, region)]["precision"] for n in present]
        ax.plot(xs, prs, "--^", color=prc, alpha=0.6, zorder=2)
        ax.plot(xs, f1s, "-o", color=f1c, zorder=3)
        ax.plot(xs, mcs, "-s", color=mccc, zorder=3)
        # shade the high-quality corner (peak F1/MCC, moderate recall): best-hit..all-best
        corner = [100 * results[(n, region)]["recall"] for n in present[:2]]
        if len(corner) == 2:
            ax.axvspan(corner[0] - 1.2, corner[1] + 1.2, color="#2ca02c", alpha=0.06, zorder=0)
        # label each breadth strategy at its F1 point; last point right-aligned to avoid clipping
        for j, n in enumerate(present):
            last = (j == len(present) - 1)
            ax.annotate(SHORT[n], (xs[j], f1s[j]), fontsize=7.5,
                        ha="right" if last else "center",
                        xytext=(-4 if last else 0, -15 if j % 2 == 0 else -27),
                        textcoords="offset points", color="#333",
                        arrowprops=dict(arrowstyle="-", color="#bbb", lw=0.5))
        if ("legacy reducer", region) in results:
            m = results[("legacy reducer", region)]
            ax.scatter([100*m["recall"]], [m["f1"]], marker="*", s=170, color=f1c, edgecolor="k", zorder=5)
            ax.scatter([100*m["recall"]], [m["mcc"]], marker="*", s=170, color=mccc, edgecolor="k", zorder=5)
            ax.annotate("legacy", (100*m["recall"], m["f1"]), fontsize=7.5, ha="center",
                        xytext=(0, 9), textcoords="offset points", color="k", fontweight="bold")
        ax.set_xlim(xlo, xhi)
        ax.set_title(f"{region}   (n={results[(present[0], region)]['n']} genus-tier amplicons)", fontsize=11)
        ax.set_xlabel("% genes captured  (recall, self-exclusion)")
        ax.set_ylabel("score")
        ax.grid(alpha=0.3)
    for ax in axes.flat[len(REGIONS):]:
        ax.axis("off")
    handles = [Line2D([], [], color=prc, ls="--", marker="^", alpha=0.7, label="precision"),
               Line2D([], [], color=f1c, marker="o", label="F1"),
               Line2D([], [], color=mccc, marker="s", label="MCC"),
               Line2D([], [], color="w", marker="*", markerfacecolor="0.6", markeredgecolor="k",
                      markersize=12, label="legacy reducer (★ on F1 & MCC)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9, frameon=True,
               bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Genus-level gene capture vs F1 / MCC across selection strategies\n"
                 "(self-exclusion / novel-organism; breadth best-hit → all-reliable; green band = peak-quality corner)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0.045, 1, 0.95))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"strategy_tradeoff.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
