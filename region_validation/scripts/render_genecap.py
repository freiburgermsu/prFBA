#!/usr/bin/env python
"""S12 (gene-capture slice) — render the gene-capture + accuracy-coupling figures.

Renders DESIGN §8 figures **F4, F4b, F5/F5b, F8, F9, F12**:

  F4   gene-capture recall / precision / F1 violins, **tier-faceted**       (3xT panels)
  F4b  "% gene capture per region" bar + 95% CI                             (the literal answer)
  F5   amplicon median length vs genus taxonomic concordance               (scatter + fit + rho)
  F5b  amplicon median length vs gene-capture mean F1                       (scatter + fit + rho)
  F8   recall-precision scatter, point per amplicon, hue=region,           (PR operating points)
       F1 iso-contours, per-region centroid diamond
  F9   accuracy (best non-self identity, continuous) vs gene-capture F1,    (decoupling test)
       per-region regression + Spearman rho (+ CI)
  F12  capstone: genus-concordance x gene-capture-F1 bubble, size=self      (one-figure decision aid)
       recovery, colour/label=region+length

This script is **pure rendering**.  It reads ONLY the flat, pre-computed CSV
tables under ``region_validation/data/`` (produced by the scoring stages
S9 ``score_taxacc.py`` / S10 ``score_genecap.py`` / S11 ``build_region_meta.py``)
and writes ``region_validation/figures/<stem>.png`` (300 dpi) + ``<stem>.pdf``.
It never touches the alignment / selection / PGFam pipeline.

DUAL SELF-MODE (this renderer's core contract).  The gene-capture cell CSV
carries a ``self_mode`` column.  Per the task, the **honest exclude-source series
(``exclude_0.987`` — drop hits/reps at the τ=0.987 species threshold / source
organism) is the PRIMARY series**, drawn solid / fully saturated, and the
diagnostic ceiling ``include`` series is the **secondary overlay** (lighter,
hatched, or a thin marker behind the primary).  We resolve the two series by a
tolerant alias map so the figure renders whether the producer wrote
``exclude_0.987``/``exclude_source``/``near_identical`` for the primary and
``include``/``include_source`` for the overlay.  If only one series exists, it is
drawn as the primary and the overlay is silently omitted.

Resumability: each figure is skipped if BOTH its ``.png`` and ``.pdf`` already
exist; pass ``--force`` to re-render, or ``--only F4,F9`` to render a subset.
A figure whose input CSV is missing/empty/under-columned is skipped with a clear
message (a partial pipeline still renders what it can).

Consistent ordering + palette come from ``_common`` (never hard-coded):
  * ``REGION_ORDER``   — regions short -> long by amplicon length (one axis order).
  * ``REGION_PALETTE`` — one stable colour per region, used in every figure.
  * ``REGION_AMP_BP``  — published amplicon length per region (F5/F12 fallback).

Sole interpreter: ``/home/freiburger/Documents/py_venv/bin/python``.

    ~/Documents/py_venv/bin/python scripts/render_genecap.py [--force] [--only F4,F4b,...]

================================================================================
INPUT CSV CONTRACT  (the columns this renderer consumes; producers must honor)
================================================================================
PATHS.genecap_cells_csv  (F4, F4b, F5b, F8, F9, F12; falls back to
                          PATHS.gene_capture_long_csv if absent)
    region, domain, tier, self_mode (a.k.a. mode), recall, precision, f1,
    jaccard, best_nonself_identity (a.k.a. best_identity / best_cosine), outcome
      * self_mode primary  in {exclude_0.987, exclude_source, near_identical, exclude}
      * self_mode overlay  in {include, include_source}
      * recall is the "% gene capture" (stored as a fraction OR a percent; this
        renderer auto-detects and normalises to a fraction in [0,1]).
      * outcome in {ok, abstained, pred_empty, truth_empty, pgfam_missing, ...};
        truth_empty / pgfam_missing / no_amplicon cells are dropped from
        distributions (poisoned / undefined).

PATHS.asv_region_concordance_csv  (F5, F12; falls back to concordance_long_csv)
    region, rank, predictor, self_level, domain, concordance, n
      * F5/F12 use the Genus-rank concordance (anchor/best predictor).

PATHS.self_recovery_csv  (F12 bubble size, optional)
    region [, domain], one of {self_recovery, rate_any20, rate_top1,
    self_anchor_rate, recovered_frac}.

PATHS.region_meta_csv  (F5/F5b/F12 x-axis, optional)
    region [, domain], amplicon_len_median  (else panel REGION_AMP_BP fallback).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless; must precede pyplot import
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# --------------------------------------------------------------------------- #
# REPO-hop boilerplate: make ``import _common`` work whether this file is run
# from scripts/, from the repo root, or by absolute path.
# --------------------------------------------------------------------------- #
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import _common as C  # noqa: E402

# Pull shared constants (REGION_PALETTE values are matplotlib-ready RGBA tuples).
REGION_ORDER = list(C.REGION_ORDER)
REGION_PALETTE = dict(C.REGION_PALETTE)
REGION_AMP_BP = dict(C.REGION_AMP_BP)
DOMAIN_ORDER = ["Bacteria", "Archaea"]
TIER_ORDER = ["species", "genus", "family"]
DEFAULT_COLOR = (0.6, 0.6, 0.6, 1.0)
DPI = 300

PATHS = C.PATHS
FIGS = Path(C.FIGS)

# Map each figure id to its output stem (figures dir is FIGS).
FIG_STEMS = {
    "F4": "F4_genecapture_violins",
    "F4b": "F4b_percent_genecapture_by_region",
    "F5": "F5_amplicon_length_vs_accuracy",
    "F5b": "F5b_amplicon_length_vs_genecapture",
    "F8": "F8_recall_precision_scatter",
    "F9": "F9_accuracy_vs_genecapture_scatter",
    "F12": "F12_capstone_bubble",
}

# --------------------------------------------------------------------------- #
# Self-mode resolution: which raw label is the PRIMARY (honest exclude-source)
# series and which is the SECONDARY (include-self ceiling) overlay.  Order is
# preference order; the first present wins.
# --------------------------------------------------------------------------- #
PRIMARY_MODE_ALIASES = ["exclude_0.987", "exclude_source", "near_identical",
                        "exclude", "noself"]
OVERLAY_MODE_ALIASES = ["include", "include_source", "include_self", "withself"]

PRIMARY_LABEL = "exclude-self (τ=0.987, honest)"
OVERLAY_LABEL = "include-self (ceiling)"

# Outcomes that are undefined / poisoned and must be dropped from distributions.
DROP_OUTCOMES = {"truth_empty", "pgfam_missing", "no_amplicon"}


# --------------------------------------------------------------------------- #
# Small generic helpers
# --------------------------------------------------------------------------- #
def color_for(region: str):
    """Stable colour for a region (fallback grey for unexpected names)."""
    return REGION_PALETTE.get(region, DEFAULT_COLOR)


def order_regions(present) -> list:
    """Regions present in the data, in canonical short->long order; extras appended."""
    present = list(dict.fromkeys(present))
    head = [r for r in REGION_ORDER if r in present]
    tail = [r for r in present if r not in REGION_ORDER]
    return head + tail


def _pick(df: pd.DataFrame, *names, required=False, default=None):
    """Return the first column in ``names`` that exists in ``df`` (case-insensitive).

    The scoring stages accept a few alias names for each field; fail loudly only
    when a *required* field is wholly absent.
    """
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n in df.columns:
            return n
        if n.lower() in lower:
            return lower[n.lower()]
    if required:
        raise KeyError(f"none of {names!r} found in columns {list(df.columns)!r}")
    return default


def _read_csv(path):
    """Read a CSV if it exists and is non-empty, else return None (figure skipped)."""
    p = str(path)
    if not C.done(p):
        return None
    try:
        df = pd.read_csv(p)
    except Exception as exc:  # pragma: no cover - defensive
        print(f"  ! failed to read {p}: {exc}")
        return None
    return df if len(df) else None


def _mean_ci(values, conf: float = 0.95):
    """Mean + a 95% CI for a vector of per-genome proportions.

    Recall/precision/F1 cells are continuous in [0,1] (a *macro* mean over
    genomes, not a single binomial), so the honest interval is a t-based CI on
    the mean; falls back to the value itself when n<=1.  Returns (m, lo, hi, n).
    """
    v = np.asarray([x for x in values if x is not None and np.isfinite(x)], dtype=float)
    n = v.size
    if n == 0:
        return np.nan, np.nan, np.nan, 0
    m = float(v.mean())
    if n == 1:
        return m, m, m, 1
    from scipy import stats

    sem = v.std(ddof=1) / np.sqrt(n)
    half = float(stats.t.ppf(0.5 + conf / 2.0, n - 1) * sem)
    return m, max(0.0, m - half), min(1.0, m + half), n


def _spearman_ci(x, y, conf: float = 0.95, n_boot: int = 2000, seed: int = 1729):
    """Spearman rho with a bootstrap CI (deterministic seed).  Returns (rho, p, lo, hi)."""
    from scipy import stats

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = x.size
    if n < 3:
        return float("nan"), float("nan"), float("nan"), float("nan")
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.all(x[idx] == x[idx][0]) or np.all(y[idx] == y[idx][0]):
            boots[b] = np.nan
            continue
        boots[b] = stats.spearmanr(x[idx], y[idx]).correlation
    boots = boots[np.isfinite(boots)]
    if boots.size:
        lo, hi = np.percentile(boots, [100 * (1 - conf) / 2, 100 * (1 + conf) / 2])
    else:
        lo = hi = float("nan")
    return float(rho), float(p), float(lo), float(hi)


def _finish(fig, fig_id: str):
    """Tight-layout + write <stem>.png (300 dpi) and <stem>.pdf, then close."""
    stem = FIG_STEMS[fig_id]
    fig.tight_layout()
    png = FIGS / f"{stem}.png"
    pdf = FIGS / f"{stem}.pdf"
    fig.savefig(png, dpi=DPI, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png}")
    print(f"  wrote {pdf}")


def _skip_blank(fig_id: str, msg: str):
    """Write a labelled placeholder so a partial pipeline still produces the file."""
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes,
            fontsize=11, wrap=True)
    ax.set_axis_off()
    ax.set_title(f"{fig_id}  (no data)")
    _finish(fig, fig_id)


def _outputs_exist(fig_id: str) -> bool:
    stem = FIG_STEMS[fig_id]
    return C.done(FIGS / f"{stem}.png") and C.done(FIGS / f"{stem}.pdf")


# --------------------------------------------------------------------------- #
# Data loaders (gene-capture cells, accuracy, region meta, self-recovery)
# --------------------------------------------------------------------------- #
def load_genecap_cells():
    """Per-(genome,region,self_mode) gene-capture cells.

    Primary source ``gene_capture_cells.csv`` (DESIGN §7); falls back to the
    tidy ``gene_capture_long.csv`` (pivoted back to wide on ``metric``).
    Normalises to canonical columns:
        region, domain, tier, self_mode, recall, precision, f1, jaccard,
        best_identity, outcome
    Recall/precision/F1/Jaccard are auto-normalised to fractions in [0,1].
    """
    df = _read_csv(PATHS.genecap_cells_csv)
    source = "cells"
    if df is None:
        df = _read_csv(PATHS.gene_capture_long_csv)
        source = "long"
    if df is None:
        return None

    # The tidy long form is (cell x metric x value); pivot it back to wide so the
    # rest of this renderer sees one row per cell with recall/precision/f1/jaccard.
    if source == "long" and _pick(df, "metric") and _pick(df, "value"):
        idcols = [c for c in df.columns if c.lower() not in ("metric", "value")]
        mcol, vcol = _pick(df, "metric"), _pick(df, "value")
        df["value"] = pd.to_numeric(df[vcol], errors="coerce")
        df = (df.pivot_table(index=idcols, columns=mcol, values="value",
                             aggfunc="mean").reset_index())
        df.columns.name = None

    out = pd.DataFrame()
    out["region"] = df[_pick(df, "region", required=True)].astype(str)
    c = _pick(df, "domain")
    out["domain"] = df[c].astype(str) if c else "Bacteria"
    c = _pick(df, "tier", "best_tier", "confidence_tier")
    out["tier"] = df[c].astype(str).str.lower() if c else "all"
    c = _pick(df, "self_mode", "mode", "exclude_mode")
    out["self_mode"] = df[c].astype(str) if c else "include_source"
    for canon, *aliases in (
        ("recall", "recall", "recall_frac", "pct_gene_capture", "gene_capture"),
        ("precision", "precision", "prec"),
        ("f1", "f1", "F1", "f1_score"),
        ("jaccard", "jaccard", "jacc"),
        ("best_identity", "best_nonself_identity", "best_identity",
         "nonself_identity", "best_cosine", "identity"),
    ):
        col = _pick(df, *aliases)
        out[canon] = pd.to_numeric(df[col], errors="coerce") if col else np.nan
    c = _pick(df, "outcome", "edge_outcome", "cell_outcome")
    out["outcome"] = df[c].astype(str).str.lower() if c else ""

    # Recall etc. may be stored as a percentage; rescale to a fraction.
    for m in ("recall", "precision", "f1", "jaccard"):
        col = out[m].dropna()
        if len(col) and col.max() > 1.5:
            out[m] = out[m] / 100.0
    return out


def resolve_modes(cells: pd.DataFrame):
    """Return ``(primary_mode, overlay_mode_or_None)`` raw labels present in data.

    PRIMARY = honest exclude-source series; OVERLAY = include-self ceiling.
    Either may be None if not present in the CSV.
    """
    modes = set(cells["self_mode"].dropna().astype(str).unique())
    primary = next((m for m in PRIMARY_MODE_ALIASES if m in modes), None)
    overlay = next((m for m in OVERLAY_MODE_ALIASES if m in modes), None)
    if primary is None and overlay is None:
        primary = sorted(modes)[0] if modes else None
    if primary is None and overlay is not None:
        # Only the include series exists — promote it to primary (best effort).
        primary, overlay = overlay, None
    return primary, overlay


def _clean(cells: pd.DataFrame, mode) -> pd.DataFrame:
    """Rows for one self_mode with poisoned/undefined outcomes dropped."""
    if mode is None:
        return cells.iloc[0:0]
    sub = cells[cells["self_mode"].astype(str) == str(mode)].copy()
    if "outcome" in sub:
        sub = sub[~sub["outcome"].isin(DROP_OUTCOMES)]
    return sub


def load_accuracy():
    """Per-region taxonomic concordance for the Genus rank (F5, F12).

    Source ``asv_region_concordance.csv`` (region x rank x predictor x self-level);
    falls back to ``concordance_long.csv``.  Normalises to:
        region, rank, predictor, self_level, domain, concordance, n
    """
    df = _read_csv(PATHS.asv_region_concordance_csv)
    if df is None:
        df = _read_csv(PATHS.concordance_long_csv)
    if df is None:
        return None
    # concordance_long carries per-phylum breakdown rows; drop them (headline only).
    if "phylum" in df.columns:
        df = df[df["phylum"].isna() | (df["phylum"].astype(str).str.strip() == "")]
    out = pd.DataFrame()
    out["region"] = df[_pick(df, "region", required=True)].astype(str)
    c = _pick(df, "rank")
    out["rank"] = df[c].astype(str) if c else ""
    c = _pick(df, "predictor")
    out["predictor"] = df[c].astype(str) if c else "anchor"
    c = _pick(df, "self_level", "self_mode", "exclusion", "level")
    out["self_level"] = df[c].astype(str) if c else ""
    c = _pick(df, "domain")
    out["domain"] = df[c].astype(str) if c else "Bacteria"
    col = _pick(df, "concordance", "correct_frac", "accuracy", "value", "correct")
    out["concordance"] = pd.to_numeric(df[col], errors="coerce") if col else np.nan
    c = _pick(df, "n_evaluable", "n", "n_not_none", "denom")
    out["n"] = pd.to_numeric(df[c], errors="coerce") if c else np.nan
    return out


def genus_concordance_by_region(acc) -> dict:
    """{region: genus concordance}; prefers include-self (this study's headline)
    anchor/best predictor, domain='all' rows where present."""
    if acc is None:
        return {}
    df = acc.copy()
    df = df[df["rank"].str.lower() == "genus"]
    if df.empty:
        return {}
    if "domain" in df.columns and (df["domain"].str.lower() == "all").any():
        df = df[df["domain"].str.lower() == "all"]
    # Prefer the include-self level (USER DECISION 4) then any honest level.
    levels = set(df["self_level"].str.lower().unique())
    for lv in ("include_self", "include_source", "include", "near_identical",
               "exclude_source", "exclude_0.987"):
        if lv in levels:
            df = df[df["self_level"].str.lower() == lv]
            break
    preds = set(df["predictor"].str.lower().unique())
    for pr in ("best", "anchor", "any20"):
        if pr in preds:
            df = df[df["predictor"].str.lower() == pr]
            break
    return df.groupby("region")["concordance"].mean().to_dict()


def load_region_meta():
    """Region metadata: median amplicon length per region (DESIGN §11)."""
    df = _read_csv(PATHS.region_meta_csv)
    if df is None:
        return None
    out = pd.DataFrame()
    out["region"] = df[_pick(df, "region", required=True)].astype(str)
    col = _pick(df, "amplicon_len_median", "median_amplicon_len", "amp_len_median",
                "amplicon_median_bp", "median_len")
    if col:
        out["amplicon_len_median"] = pd.to_numeric(df[col], errors="coerce")
    else:
        out["amplicon_len_median"] = out["region"].map(REGION_AMP_BP)
    c = _pick(df, "domain")
    out["domain"] = df[c].astype(str) if c else "Bacteria"
    return out


def region_lengths() -> dict:
    """{region: median amplicon bp} from region_meta if available, else panel amp_bp."""
    meta = load_region_meta()
    if meta is not None:
        lens = (meta.dropna(subset=["amplicon_len_median"])
                .groupby("region")["amplicon_len_median"].median().to_dict())
        if lens:
            # Backfill any region missing a measured median with the panel value.
            for r, v in REGION_AMP_BP.items():
                lens.setdefault(r, v)
            return lens
    return dict(REGION_AMP_BP)


def load_self_recovery() -> dict:
    """{region: self-recovery fraction} for the F12 bubble size (optional)."""
    df = _read_csv(PATHS.self_recovery_csv)
    if df is None:
        return {}
    rcol = _pick(df, "region", required=True)
    if "domain" in df.columns and (df["domain"].astype(str).str.lower() == "all").any():
        df = df[df["domain"].astype(str).str.lower() == "all"]
    vcol = _pick(df, "self_recovery", "recovered_frac", "rate_any20",
                 "self_anytop20_rate", "self_anchor_rate", "rate_top1",
                 "any_top20_rate", "fraction", "recovery")
    if not vcol:
        return {}
    return (df.groupby(rcol)[vcol]
            .apply(lambda s: pd.to_numeric(s, errors="coerce").mean()).to_dict())


# --------------------------------------------------------------------------- #
# Violin drawing (primary solid + overlay outline)
# --------------------------------------------------------------------------- #
def _violin_dual(ax, regions, prim_by_region, over_by_region, ylabel, title):
    """Per-region violins: PRIMARY filled (+ quartile strip) with the OVERLAY
    drawn as a translucent, narrower outline behind it (the include-self ceiling)."""
    xs = list(range(1, len(regions) + 1))
    drew_overlay = False

    # Overlay first (behind): narrower, no fill highlight, hatched outline.
    if over_by_region:
        over = [np.asarray(over_by_region.get(r, []), dtype=float) for r in regions]
        over = [s[np.isfinite(s)] for s in over]
        ov = [(x, s) for x, s in zip(xs, over) if s.size]
        if ov:
            drew_overlay = True
            parts = ax.violinplot([s for _, s in ov], positions=[x for x, _ in ov],
                                  widths=0.92, showextrema=False)
            for body in parts["bodies"]:
                body.set_facecolor("none")
                body.set_edgecolor("0.45")
                body.set_linewidth(1.0)
                body.set_linestyle((0, (3, 2)))
                body.set_alpha(0.9)

    # Primary on top: region-coloured fill + quartile/median strip.
    prim = [np.asarray(prim_by_region.get(r, []), dtype=float) for r in regions]
    prim = [s[np.isfinite(s)] for s in prim]
    pr = [(x, s, r) for x, s, r in zip(xs, prim, regions) if s.size]
    if not pr and not drew_overlay:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        ax.set_xticks(xs)
        ax.set_xticklabels(regions, rotation=40, ha="right", fontsize=8)
        return
    if pr:
        parts = ax.violinplot([s for _, s, _ in pr], positions=[x for x, _, _ in pr],
                              widths=0.78, showextrema=False)
        for body, (_, _, r) in zip(parts["bodies"], pr):
            body.set_facecolor(color_for(r))
            body.set_edgecolor("0.2")
            body.set_alpha(0.8)
        for x, s, _ in pr:
            q1, med, q3 = np.percentile(s, [25, 50, 75])
            ax.vlines(x, q1, q3, color="0.12", lw=4, zorder=3)
            ax.scatter([x], [med], color="white", edgecolor="0.1", s=22, zorder=4)

    ax.set_xticks(xs)
    ax.set_xticklabels(regions, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(axis="y", ls=":", alpha=0.4)


# --------------------------------------------------------------------------- #
# F4 — gene-capture recall/precision/F1 violins, tier-faceted
# --------------------------------------------------------------------------- #
def fig_F4(cells: pd.DataFrame):
    """Gene-capture recall/precision/F1 violins, tier-faceted; PRIMARY exclude-self
    filled, include-self ceiling overlaid as a dashed outline."""
    primary, overlay = resolve_modes(cells)
    prim = _clean(cells, primary)
    over = _clean(cells, overlay)
    if prim.empty and over.empty:
        return _skip_blank("F4", "gene-capture cells present but no scorable rows")

    metrics = [("recall", "recall (% gene capture)"),
               ("precision", "precision"),
               ("f1", "F1")]
    present_tiers = set(prim["tier"].unique()) | set(over["tier"].unique())
    tiers = [t for t in TIER_ORDER if t in present_tiers]
    if not tiers:
        tiers = sorted(present_tiers) or ["all"]
    regions = order_regions(pd.concat([prim["region"], over["region"]]).unique())

    nrow, ncol = len(metrics), len(tiers)
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(max(4.5, 0.85 * len(regions)) * ncol,
                                      3.0 * nrow),
                             squeeze=False, sharex=True)
    for i, (mkey, mlabel) in enumerate(metrics):
        for j, tier in enumerate(tiers):
            ax = axes[i][j]
            ptier = prim if tier == "all" else prim[prim["tier"] == tier]
            otier = over if tier == "all" else over[over["tier"] == tier]
            prim_by = {r: ptier[ptier["region"] == r][mkey].dropna().values
                       for r in regions}
            over_by = {r: otier[otier["region"] == r][mkey].dropna().values
                       for r in regions} if not over.empty else None
            title = f"tier={tier}" if i == 0 else ""
            _violin_dual(ax, regions, prim_by, over_by,
                         ylabel=mlabel if j == 0 else "", title=title)
            if i != nrow - 1:
                ax.set_xticklabels([])
    handles = [Patch(facecolor="0.5", edgecolor="0.2", alpha=0.8,
                     label=f"PRIMARY · {PRIMARY_LABEL}")]
    if overlay is not None and not over.empty:
        handles.append(Line2D([], [], color="0.45", lw=1.2, ls=(0, (3, 2)),
                              label=f"overlay · {OVERLAY_LABEL}"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("F4 · Gene-capture distributions by region "
                 f"(primary={primary}; tier-faceted)", y=1.0, fontsize=12)
    _finish(fig, "F4")


# --------------------------------------------------------------------------- #
# F4b — % gene capture per region (bar + Wilson/mean CI)
# --------------------------------------------------------------------------- #
def fig_F4b(cells: pd.DataFrame):
    """% gene capture per region: PRIMARY exclude-self bar + 95% CI, with the
    include-self ceiling shown as an open marker (the literal headline answer)."""
    primary, overlay = resolve_modes(cells)
    prim = _clean(cells, primary)
    over = _clean(cells, overlay)
    if prim.empty and over.empty:
        return _skip_blank("F4b", "gene-capture cells present but no scorable rows")
    regions = order_regions(pd.concat([prim["region"], over["region"]]).unique())

    fig, ax = plt.subplots(figsize=(max(6.0, 0.95 * len(regions) + 2.5), 4.8))
    xs = np.arange(len(regions))
    means, los, his, ns = [], [], [], []
    over_means = []
    for r in regions:
        m, lo, hi, n = _mean_ci(prim[prim["region"] == r]["recall"].dropna().values)
        means.append(100.0 * m if np.isfinite(m) else 0.0)
        los.append(100.0 * (m - lo) if np.isfinite(m) else 0.0)
        his.append(100.0 * (hi - m) if np.isfinite(m) else 0.0)
        ns.append(n)
        om, *_ = _mean_ci(over[over["region"] == r]["recall"].dropna().values) \
            if not over.empty else (np.nan,)
        over_means.append(100.0 * om if np.isfinite(om) else np.nan)

    ax.bar(xs, means, color=[color_for(r) for r in regions],
           edgecolor="0.2", width=0.74, zorder=2,
           label=f"PRIMARY · {PRIMARY_LABEL}")
    ax.errorbar(xs, means, yerr=[los, his], fmt="none", ecolor="0.12",
                elinewidth=1.4, capsize=4, zorder=3)
    if overlay is not None and any(np.isfinite(over_means)):
        ax.scatter(xs, over_means, marker="_", s=420, color="0.1",
                   linewidth=2.2, zorder=4, label=f"{OVERLAY_LABEL}")
    for x, m, n in zip(xs, means, ns):
        ax.text(x, m + 1.5, f"n={n}", ha="center", va="bottom",
                fontsize=7, color="0.25")
    ax.set_xticks(xs)
    ax.set_xticklabels(regions, rotation=40, ha="right", fontsize=9)
    ax.set_ylabel("mean % gene capture (recall)")
    top = max([v for v in (means + over_means) if np.isfinite(v)] + [1.0])
    ax.set_ylim(0, min(105, top + 12))
    ax.set_title(f"F4b · % gene capture per region (primary={primary}, mean ± 95% CI)")
    ax.grid(axis="y", ls=":", alpha=0.4)
    ax.legend(fontsize=8, loc="best", framealpha=0.9)
    _finish(fig, "F4b")


# --------------------------------------------------------------------------- #
# F5 / F5b — amplicon length vs accuracy / gene-capture (scatter + fit)
# --------------------------------------------------------------------------- #
def _scatter_len_vs(ax, lengths: dict, yvals: dict, ylabel: str, title: str):
    """Per-region point at (median amplicon bp, y); fit a line; annotate corr."""
    regions = order_regions([r for r in yvals
                             if r in lengths and np.isfinite(yvals[r])
                             and np.isfinite(lengths[r])])
    if not regions:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return False
    x = np.array([lengths[r] for r in regions], dtype=float)
    y = np.array([yvals[r] for r in regions], dtype=float)
    for r, xi, yi in zip(regions, x, y):
        ax.scatter([xi], [yi], color=color_for(r), s=130, edgecolor="0.2", zorder=3)
        ax.annotate(r, (xi, yi), xytext=(5, 4), textcoords="offset points", fontsize=8)
    if len(x) >= 2:
        b, a = np.polyfit(x, y, 1)
        xx = np.linspace(x.min(), x.max(), 50)
        ax.plot(xx, a + b * xx, color="0.3", ls="--", lw=1.4, zorder=2)
        from scipy import stats

        pr, pp = stats.pearsonr(x, y)
        rho, rp = stats.spearmanr(x, y)
        ax.text(0.03, 0.97,
                f"Pearson r={pr:.2f} (p={pp:.2g})\nSpearman ρ={rho:.2f} (p={rp:.2g})",
                transform=ax.transAxes, va="top", ha="left", fontsize=8,
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    ax.set_xlabel("median amplicon length (bp)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(ls=":", alpha=0.4)
    return True


def fig_F5(cells: pd.DataFrame, acc):
    """F5: amplicon median length vs genus taxonomic concordance (per region)."""
    lengths = region_lengths()
    yvals = genus_concordance_by_region(acc)
    if not yvals:
        return _skip_blank("F5", "no genus-concordance table (asv_region_concordance.csv) "
                                 "available yet")
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    _scatter_len_vs(ax, lengths, yvals, "genus taxonomic concordance",
                    "F5 · Amplicon length vs taxonomic accuracy")
    _finish(fig, "F5")


def fig_F5b(cells: pd.DataFrame):
    """F5b: amplicon median length vs gene-capture mean F1 (per region, PRIMARY mode)."""
    lengths = region_lengths()
    primary, _ = resolve_modes(cells)
    prim = _clean(cells, primary)
    if prim.empty:
        return _skip_blank("F5b", "gene-capture cells present but no scorable rows")
    yvals = prim.groupby("region")["f1"].mean().to_dict()
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    _scatter_len_vs(ax, lengths, yvals, f"mean gene-capture F1 (primary={primary})",
                    "F5b · Amplicon length vs gene capture")
    _finish(fig, "F5b")


# --------------------------------------------------------------------------- #
# F8 — recall-precision scatter (point per amplicon, hue=region)
# --------------------------------------------------------------------------- #
def fig_F8(cells: pd.DataFrame):
    """F8: recall-precision scatter, point per amplicon, hue=region, F1 iso-contours,
    per-region centroid diamond.  PRIMARY exclude-self; include-self centroid as
    an open diamond for the leakage gap."""
    primary, overlay = resolve_modes(cells)
    prim = _clean(cells, primary).dropna(subset=["recall", "precision"])
    over = _clean(cells, overlay).dropna(subset=["recall", "precision"])
    if prim.empty and over.empty:
        return _skip_blank("F8", "no (recall, precision) pairs available")
    regions = order_regions(pd.concat([prim["region"], over["region"]]).unique())

    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    # F1 iso-contours.
    gr = np.linspace(0.001, 1.0, 200)
    R, P = np.meshgrid(gr, gr)
    F1 = 2 * R * P / (R + P)
    cs = ax.contour(P, R, F1, levels=[0.2, 0.4, 0.6, 0.8], colors="0.7",
                    linewidths=0.8, linestyles=":")
    ax.clabel(cs, inline=True, fontsize=7, fmt="F1=%.1f")

    for r in regions:
        rs = prim[prim["region"] == r]
        ax.scatter(rs["precision"], rs["recall"], s=10, alpha=0.28,
                   color=color_for(r), edgecolor="none", zorder=2)
    handles = []
    for r in regions:
        rs = prim[prim["region"] == r]
        if not rs.empty:
            cx, cy = rs["precision"].mean(), rs["recall"].mean()
            ax.scatter([cx], [cy], marker="D", s=140, color=color_for(r),
                       edgecolor="black", linewidth=1.1, zorder=5)
        os_ = over[over["region"] == r]
        if not os_.empty:
            ox, oy = os_["precision"].mean(), os_["recall"].mean()
            ax.scatter([ox], [oy], marker="D", s=140, facecolor="none",
                       edgecolor=color_for(r), linewidth=2.0, zorder=4)
        handles.append(Line2D([], [], marker="D", linestyle="none",
                              markerfacecolor=color_for(r), markeredgecolor="black",
                              markersize=9, label=r))
    extra = [Line2D([], [], marker="D", linestyle="none", markerfacecolor="0.4",
                    markeredgecolor="black", markersize=9,
                    label=f"◆ {PRIMARY_LABEL}")]
    if overlay is not None and not over.empty:
        extra.append(Line2D([], [], marker="D", linestyle="none", markerfacecolor="none",
                            markeredgecolor="0.4", markersize=9,
                            label=f"◇ {OVERLAY_LABEL}"))
    ax.set_xlabel("precision  (i / |P|)")
    ax.set_ylabel("recall  (i / |T| = % gene capture)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"F8 · Recall–precision per amplicon (primary={primary}); "
                 "◆ region centroid")
    leg1 = ax.legend(handles=handles, fontsize=8, title="region", loc="lower left",
                     framealpha=0.9, ncol=2)
    ax.add_artist(leg1)
    ax.legend(handles=extra, fontsize=8, loc="upper right", framealpha=0.9)
    ax.grid(ls=":", alpha=0.4)
    _finish(fig, "F8")


# --------------------------------------------------------------------------- #
# F9 — accuracy (best non-self identity) vs gene-capture F1 + Spearman
# --------------------------------------------------------------------------- #
def fig_F9(cells: pd.DataFrame):
    """F9: best non-self identity (continuous taxonomic resolvability) vs gene-capture
    F1, per-region regression + Spearman rho (with bootstrap CI for the overall rho)."""
    primary, _ = resolve_modes(cells)
    sub = _clean(cells, primary).dropna(subset=["best_identity", "f1"])
    if sub.empty:
        return _skip_blank("F9", "no (best_identity, F1) pairs available "
                                 "(needs best_nonself_identity column)")
    regions = order_regions(sub["region"].unique())

    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    from scipy import stats

    rho_all, p_all, lo_all, hi_all = _spearman_ci(sub["best_identity"].values,
                                                  sub["f1"].values)

    handles = []
    for r in regions:
        rs = sub[sub["region"] == r]
        if rs.empty:
            continue
        col = color_for(r)
        ax.scatter(rs["best_identity"], rs["f1"], s=10, alpha=0.28, color=col,
                   edgecolor="none", zorder=2)
        if len(rs) >= 3:
            x = rs["best_identity"].values
            y = rs["f1"].values
            b, a = np.polyfit(x, y, 1)
            xx = np.linspace(x.min(), x.max(), 30)
            ax.plot(xx, a + b * xx, color=col, lw=1.8, zorder=4)
            rho, _ = stats.spearmanr(x, y)
            lab = f"{r} (ρ={rho:.2f})"
        else:
            lab = r
        handles.append(Line2D([], [], color=col, lw=2.2, label=lab))
    ax.set_xlabel("best non-self reference identity  (taxonomic resolvability)")
    ax.set_ylabel(f"gene-capture F1 (primary={primary})")
    ax.set_ylim(-0.02, 1.02)
    ci_txt = (f" [95% CI {lo_all:.2f}, {hi_all:.2f}]"
              if np.isfinite(lo_all) and np.isfinite(hi_all) else "")
    ax.set_title("F9 · Accuracy vs gene capture — are resolvable amplicons also "
                 f"clean gene-capture?\noverall Spearman ρ={rho_all:.2f}"
                 f"{ci_txt} (p={p_all:.2g})")
    ax.legend(handles=handles, fontsize=7.5, title="region (per-region ρ)",
              loc="best", framealpha=0.9, ncol=2)
    ax.grid(ls=":", alpha=0.4)
    _finish(fig, "F9")


# --------------------------------------------------------------------------- #
# F12 — capstone bubble (genus concordance x gene-capture F1, size=self-recovery)
# --------------------------------------------------------------------------- #
def fig_F12(cells: pd.DataFrame, acc):
    """F12 capstone: genus concordance x gene-capture F1 bubble, size=self-recovery,
    colour/label=region+length; the one-figure decision aid (Pareto-dominant region)."""
    primary, _ = resolve_modes(cells)
    prim = _clean(cells, primary)
    f1_by_region = prim.groupby("region")["f1"].mean().to_dict() if not prim.empty else {}
    genus_by_region = genus_concordance_by_region(acc)
    selfrec = load_self_recovery()
    lengths = region_lengths()

    regions = order_regions([r for r in f1_by_region
                             if r in genus_by_region
                             and np.isfinite(f1_by_region[r])
                             and np.isfinite(genus_by_region[r])])
    if not regions:
        return _skip_blank("F12", "needs BOTH gene-capture F1 (cells) and genus "
                                  "concordance (asv_region_concordance.csv)")

    fig, ax = plt.subplots(figsize=(8.6, 7.2))
    sizes = []
    for r in regions:
        sr = selfrec.get(r)
        sizes.append(200 + 1400 * sr if sr is not None and np.isfinite(sr) else 600)
    xs = [genus_by_region[r] for r in regions]
    ys = [f1_by_region[r] for r in regions]
    for r, x, y, s in zip(regions, xs, ys, sizes):
        ax.scatter([x], [y], s=s, color=color_for(r), edgecolor="0.15",
                   linewidth=1.2, alpha=0.85, zorder=3)
        lab = f"{r}\n({int(lengths[r])} bp)" if r in lengths and np.isfinite(lengths[r]) else r
        ax.annotate(lab, (x, y), xytext=(7, 6), textcoords="offset points",
                    fontsize=8, zorder=5)

    ax.axhline(float(np.nanmean(ys)), color="0.8", ls=":", lw=1)
    ax.axvline(float(np.nanmean(xs)), color="0.8", ls=":", lw=1)
    ax.set_xlabel("genus taxonomic concordance")
    ax.set_ylabel(f"gene-capture F1 (primary={primary})")
    ax.set_title("F12 · Capstone: accuracy × gene capture "
                 "(bubble size = self-recovery; upper-right = Pareto-dominant)")
    ax.grid(ls=":", alpha=0.4)
    if any(selfrec.get(r) is not None for r in regions):
        for sval, slab in ((0.5, "0.5"), (0.9, "0.9")):
            ax.scatter([], [], s=200 + 1400 * sval, color="0.6", edgecolor="0.2",
                       alpha=0.6, label=f"self-recovery {slab}")
        ax.legend(scatterpoints=1, fontsize=8, labelspacing=1.6,
                  borderpad=1.0, loc="lower right", framealpha=0.9)
    _finish(fig, "F12")


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
# Loaders each figure needs (so a figure is skipped cleanly when its inputs are
# missing rather than crashing the whole render).
_NEEDS = {
    "F4": ("cells",),
    "F4b": ("cells",),
    "F5": ("cells", "acc"),
    "F5b": ("cells",),
    "F8": ("cells",),
    "F9": ("cells",),
    "F12": ("cells", "acc"),
}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Render gene-capture figures F4/F4b/F5/F5b/F8/F9/F12 "
                    "(primary=exclude-self τ=0.987, include-self overlay) from "
                    "the flat CSVs in _common.PATHS.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="re-render even if a figure's .png/.pdf already exist")
    ap.add_argument("--only", default="",
                    help="comma list of figure ids (e.g. F4,F9); default all")
    args = ap.parse_args(argv)

    C.ensure_dirs()
    requested = [f.strip() for f in args.only.split(",") if f.strip()] or list(FIG_STEMS)
    unknown = [f for f in requested if f not in FIG_STEMS]
    if unknown:
        ap.error(f"unknown figure id(s): {unknown}; valid: {list(FIG_STEMS)}")

    print(f"render_genecap → {FIGS}")
    print(f"  REGION_ORDER (short→long): {REGION_ORDER}")
    print(f"  figures: {requested}{'  | --force' if args.force else ''}")

    # Lazy-load shared tables once.
    loaded: dict = {}

    def get(name):
        if name not in loaded:
            loaded[name] = {"cells": load_genecap_cells, "acc": load_accuracy}[name]()
        return loaded[name]

    rendered, skipped, missing = [], [], []
    for fid in requested:
        if not args.force and _outputs_exist(fid):
            print(f"  [{fid}] skip (exists) — {FIG_STEMS[fid]}")
            skipped.append(fid)
            continue
        needs = _NEEDS[fid]
        data = {n: get(n) for n in needs}
        if data.get("cells") is None and "cells" in needs:
            missing.append(fid)
            print(f"  [{fid}] SKIP — gene-capture cells not available yet "
                  f"({PATHS.genecap_cells_csv})")
            continue
        if fid in ("F5", "F12") and data.get("acc") is None:
            # These need accuracy too; render a labelled placeholder.
            print(f"  [{fid}] rendering placeholder (accuracy table absent)")
        print(f"[{fid}] rendering → {FIG_STEMS[fid]}")
        try:
            if fid == "F4":
                fig_F4(data["cells"])
            elif fid == "F4b":
                fig_F4b(data["cells"])
            elif fid == "F5":
                fig_F5(data["cells"], data.get("acc"))
            elif fid == "F5b":
                fig_F5b(data["cells"])
            elif fid == "F8":
                fig_F8(data["cells"])
            elif fid == "F9":
                fig_F9(data["cells"])
            elif fid == "F12":
                fig_F12(data["cells"], data.get("acc"))
            rendered.append(fid)
        except Exception as exc:
            missing.append(fid)
            print(f"  ! {fid} failed: {exc}")
            import traceback

            traceback.print_exc()

    print("\n=== render_genecap summary ===")
    print(f"rendered: {rendered or '-'}")
    print(f"skipped (already present, use --force): {skipped or '-'}")
    print(f"not rendered (missing inputs / error): {missing or '-'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
