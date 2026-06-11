#!/usr/bin/env python
"""S12 (accuracy slice) — render the taxonomic-accuracy figure suite.

Renders DESIGN.md §8 figures **F1, F2, F3, F6, F7, F10, F11, F13, F14**:

  F1  Accuracy heatmap region x rank        (1x3: anchor / consensus / coverage)
  F2  Per-rank concordance by region        (grouped lines + markers + Wilson CI)
  F3  Deepest-correct-rank distribution      (100%-stacked horizontal bar)
  F6  Self-recovery rate by region           (bar + Wilson CI, top-1 vs any-top-20)
  F7  Genus-failure confusion                (truth x pred genus heatmap, LogNorm)
  F10 Per-phylum accuracy                     (phylum x region heatmap, n<10 masked)
  F11 Self-inclusion vs self-exclusion delta  (per region x rank; the leakage gap)
  F13 Amplification-rate / selection-bias     (grouped bar, region x phylum, Bact/Arch)
  F14 Tie/ambiguity cluster size              (ECDF per region, per-amplicon)

This is a **standalone renderer** (matplotlib ``Agg``).  It reads ONLY the flat,
pre-computed CSVs addressed by ``_common.PATHS`` (written by the scoring stages
S9 ``score_taxacc.py`` and S11 ``build_region_meta.py``) and writes
``FIGS/F*.png`` @300 dpi + ``F*.pdf``.  It performs **no** scoring / alignment /
selection itself.

Self-mode handling (DESIGN USER DECISIONS banner item 4 — report BOTH scenarios).
The realistic, include-self ("source organism in the DB") pipeline is the PRIMARY
series; a self-exclusion ("exclude_0.987" / near-identical) series, when present in
the source CSVs, is drawn as a secondary overlay, and **F11 = include − exclude**
delta per region.  The current S9 scorer (``score_taxacc.py``) emits an
**include-self only** table (no ``self_mode`` column), so the exclude series is
absent: F11 then short-circuits with a clear message and the other figures render
the include-self series alone.  The renderer auto-detects whichever self-mode
column the producers actually emit (``self_mode`` / ``self_level`` / ``mode``) so
it works unchanged the day an exclude-self branch is added upstream.

Predictors.  ``score_taxacc.py`` scores two predictors, ``anchor`` (selected[0],
the realistic pipeline call) and ``consensus`` (per-rank majority).  ``anchor`` is
the PRIMARY series everywhere; ``consensus`` is the F1 secondary panel.

Shared constants come from ``_common`` (never hard-coded):
  * ``REGION_ORDER``   — regions short -> long by amplicon length (one axis order).
  * ``REGION_PALETTE`` — one stable color per region, used in every figure.
  * ``DEPTH_RAMP`` / ``DEPTH_LABELS`` — the deepest-correct stacked-bar ramp.
  * ``SCORED``         — Phylum..Species (the scored ranks).
  * ``wilson_ci``      — Wilson 95% CI (fallback when a CSV omits a CI column).

Resumable: each figure is skipped if BOTH its ``.png`` and ``.pdf`` already
exist (pass ``--force`` to re-render, or ``--only F1,F3`` to render a subset).
A figure whose input CSV is missing / empty / under-columned is skipped with a
clear message (so a partial pipeline still renders what it can).

Sole interpreter: ``/home/freiburger/Documents/py_venv/bin/python``.

    $PY scripts/render_accuracy.py [--force] [--only F1,F2,...]

================================================================================
INPUT CSV CONTRACT  (the columns this renderer consumes; producers honor these —
verified against scripts/score_taxacc.py (S9) and scripts/build_region_meta.py (S11))
================================================================================
PATHS.concordance_long_csv   (F1, F2, F11)   <- score_taxacc.py
    region, domain, tier, predictor, rank,
    n, n_called, n_correct,
    correct, correct_ci_lo, correct_ci_hi,
    coverage, coverage_ci_lo, coverage_ci_hi
      * predictor in {anchor, consensus}
      * rank in SCORED (Phylum..Species)
      * domain in {Bacteria, Archaea, all};  tier in {species,genus,family,below,all}
      * ``correct`` = n_correct / n_called (the per-rank concordance)
      * (optional) a self-mode column {self_mode|self_level|mode}; absent today.
      * (optional) a ``phylum`` column for the F10 per-phylum breakdown; absent today.

PATHS.asv_region_concordance_csv   (F3, F14)   <- score_taxacc.py
    amplicon, region, domain, tier, predictor,
    src_genome_id, src_taxon_id, best_identity, n_selected,
    self_in_top20, self_selected, self_is_anchor, tie_cluster_size,
    correct_<rank> (one per SCORED rank; "" = no-call),
    deepest_correct_rank, first_wrong_rank, misassign_kind,
    overconfident, abstained
      * one row per (amplicon, predictor); F3 derives the deepest-correct
        fractions from ``deepest_correct_rank`` ("" -> "none"); F14 reads the
        per-amplicon ``tie_cluster_size``.

PATHS.self_recovery_csv   (F6)   <- score_taxacc.py
    region, domain, tier, n,
    n_self_in_top20, n_self_selected, n_self_anchor,
    self_top20_rate, self_top20_ci_lo, self_top20_ci_hi,
    self_anchor_rate, self_anchor_ci_lo, self_anchor_ci_hi,
    n_selected_given_top20, selected_given_top20_rate,
    selected_given_top20_ci_lo, selected_given_top20_ci_hi,
    n_anchor_given_top20, anchor_given_top20_rate
      * F6: recovered-given-available = selected_given_top20_rate (any-top-20 ->
        selected) vs anchor_given_top20_rate (top-1 -> anchor), each over
        n_self_in_top20 with a Wilson CI.

PATHS.confusion_genus_csv   (F7)   <- score_taxacc.py
    region, true_genus, pred_genus, count
      * one row per observed (true_genus -> pred_genus) anchor misassignment.

PATHS.amplification_by_phylum_csv   (F13)   <- build_region_meta.py
    region, domain, phylum, n_amplified, n_source, extract_rate
      * extract_rate = n_amplified / n_source in [0,1]; domain in {Bacteria,Archaea,...}.

PATHS.region_meta_csv   (axis annotation only, optional)   <- build_region_meta.py
    region, ..., amp_len_median, extract_rate, <version stamps...>
"""
from __future__ import annotations

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")  # headless; must precede pyplot import
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# --- REPO-hop: make `import _common` work no matter the cwd ----------------- #
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import _common as C  # noqa: E402
from _common import (  # noqa: E402
    PATHS,
    SCORED,
    REGION_ORDER,
    REGION_PALETTE,
    DEPTH_LABELS,
    DEPTH_RAMP,
)

DPI = 300

# Primary predictor = the realistic include-self pipeline call (selected[0]);
# consensus is the secondary overlay panel.
PRIMARY_PREDICTOR = "anchor"
SECONDARY_PREDICTOR = "consensus"

# Aliases for the (optional) self-mode column the producers may emit, and the
# value we treat as the realistic-include / honest-exclude series.  None of
# these are present in the current include-self-only tables; the renderer simply
# falls back to "one series".
SELFMODE_COLS = ("self_mode", "self_level", "mode", "exclusion", "level")
INCLUDE_TOKENS = ("include", "include_self", "withself", "with_self")
EXCLUDE_TOKENS = (
    "exclude_0.987", "exclude_0_987", "near_identical", "exclude_source",
    "exclude", "noself", "no_self",
)

# Short, presentation-friendly rank labels.
RANK_SHORT = {"Phylum": "P", "Class": "C", "Order": "O",
              "Family": "F", "Genus": "G", "Species": "S"}

# Archaeal-phylum markers (for the F10 Bacteria/Archaea split + F13 ordering).
_ARCH_MARKERS = ("archaeo", "halobacterota", "methanobacteriota", "nanoarchaeota",
                 "thermoplasmatota", "crenarchaeota", "euryarchaeota",
                 "thaumarchaeota", "nanoarchaeota", "asgardarchaeota")


# --------------------------------------------------------------------------- #
# Small shared helpers
# --------------------------------------------------------------------------- #
def _pick(df: pd.DataFrame, *names, required=False, default=None):
    """First column in ``names`` present in ``df`` (case-insensitive), else default.

    The scoring producers use stable names but a couple of fields have a few
    historical aliases; this keeps the renderer robust to either spelling and
    fails loudly only when a *required* column is wholly absent."""
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n in df.columns:
            return n
        if n.lower() in lower:
            return lower[n.lower()]
    if required:
        raise KeyError(f"none of {names!r} in columns {list(df.columns)!r}")
    return default


def _regions_present(values) -> list:
    """REGION_ORDER intersected with the regions actually present (order kept);
    any unexpected region names are appended after, deterministically sorted."""
    have = list(dict.fromkeys(str(v) for v in values))
    head = [r for r in REGION_ORDER if r in have]
    tail = sorted(r for r in have if r not in REGION_ORDER)
    return head + tail


def _selfmode_col(df: pd.DataFrame):
    """Return the name of the self-mode column if the CSV carries one, else None."""
    return _pick(df, *SELFMODE_COLS)


def _filter_selfmode(df: pd.DataFrame, want: str):
    """Keep rows of one self-mode (``include`` or ``exclude``).

    No-op when the CSV has no self-mode column (the current include-self-only
    schema): every row already *is* the include series."""
    col = _selfmode_col(df)
    if col is None:
        return df if want == "include" else df.iloc[0:0]
    tokens = INCLUDE_TOKENS if want == "include" else EXCLUDE_TOKENS
    vals = df[col].astype(str).str.lower()
    mask = pd.Series(False, index=df.index)
    for t in tokens:
        mask |= vals == t.lower()
    if not mask.any() and want == "include":
        # No recognized include token but a self-mode col exists: assume single mode.
        return df
    return df[mask]


def _rollup(df: pd.DataFrame, col: str, prefer: str):
    """Filter ``df[col]`` to the rollup value ``prefer`` (e.g. domain/tier == 'all')
    when present; otherwise leave ``df`` untouched (CSV may not carry that column)."""
    if col in df.columns and (df[col].astype(str) == prefer).any():
        return df[df[col].astype(str) == prefer]
    return df


def _have_cols(df: pd.DataFrame, cols, fig: str) -> bool:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        print(f"  [{fig}] SKIP — input missing columns {missing} "
              f"(have {list(df.columns)})")
        return False
    return True


def _load(path: str, fig: str, quiet: bool = False):
    """Read a CSV checkpoint; return None (with a message) if absent/empty.

    ``quiet`` suppresses the not-found message (used when the caller will try a
    fallback table and report its own combined skip)."""
    if not C.done(path):
        if not quiet:
            print(f"  [{fig}] SKIP — input not found / empty: {path}")
        return None
    try:
        df = pd.read_csv(path)
    except Exception as exc:  # pragma: no cover - defensive
        print(f"  [{fig}] SKIP — could not read {path}: {exc}")
        return None
    if df.empty:
        print(f"  [{fig}] SKIP — input has no rows: {path}")
        return None
    return df


def _save(fig, stem: str) -> None:
    """Write FIGS/<stem>.png @DPI and .pdf, then close."""
    png = os.path.join(C.FIGS, f"{stem}.png")
    pdf = os.path.join(C.FIGS, f"{stem}.pdf")
    fig.savefig(png, dpi=DPI, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {png}")
    print(f"  wrote {pdf}")


def _outputs_exist(stem: str) -> bool:
    return (C.done(os.path.join(C.FIGS, f"{stem}.png"))
            and C.done(os.path.join(C.FIGS, f"{stem}.pdf")))


def _color(region: str):
    return REGION_PALETTE.get(region, "0.6")


def _short_phylum(p) -> str:
    """Shorten long phylum/family labels (Proteobacteria convention from donor)."""
    if not isinstance(p, str):
        return str(p)
    return (p.replace("proteobacteria", "-proteo.")
             .replace("Proteobacteria", "Proteo."))


def _is_archaea_phylum(p) -> bool:
    return isinstance(p, str) and any(m in p.lower() for m in _ARCH_MARKERS)


def _resolve_predictors(present):
    """(primary, secondary) predictor names from those actually in the data.

    Prefers the real S9 predictors (anchor, then consensus); falls back to the
    fixture/legacy names (best -> any20 -> abw).  Returns (primary, secondary|None);
    secondary is None when only one predictor exists."""
    present = list(dict.fromkeys(str(p) for p in present))
    order = [PRIMARY_PREDICTOR, SECONDARY_PREDICTOR, "best", "any20", "abw"]
    ranked = [p for p in order if p in present] + \
             [p for p in present if p not in order]
    if not ranked:
        return None, None
    primary = ranked[0]
    secondary = ranked[1] if len(ranked) > 1 else None
    return primary, secondary


def _dedup_predictor(df: pd.DataFrame) -> pd.DataFrame:
    """Keep one predictor's rows so per-amplicon tables become one-row-per-amplicon.

    Picks the resolved primary predictor (anchor, else the first available); a
    no-op when there is no ``predictor`` column."""
    if "predictor" not in df.columns:
        return df
    primary, _sec = _resolve_predictors(df["predictor"])
    if primary is not None and (df["predictor"] == primary).any():
        return df[df["predictor"] == primary]
    return df


def _concordance_col(df: pd.DataFrame) -> str:
    """Name of the per-rank concordance column (``correct`` in the S9 schema,
    with a couple of tolerant aliases)."""
    return _pick(df, "correct", "concordance", "correct_frac", "accuracy",
                 required=True)


def _ci_cols(df: pd.DataFrame):
    """(lo, hi) column names for the concordance CI, or (None, None)."""
    lo = _pick(df, "correct_ci_lo", "ci_lo", "concordance_ci_lo")
    hi = _pick(df, "correct_ci_hi", "ci_hi", "concordance_ci_hi")
    return (lo, hi) if (lo and hi) else (None, None)


def _headline_concordance(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize the long concordance table to the headline grid used by F1/F2/F11.

    Returns a frame with canonical columns ``region, predictor, rank,
    concordance, ci_lo, ci_hi``.  When the source carries a ``phylum`` column
    (the F10 per-phylum breakdown shares the table) the per-phylum rows are
    dropped so the (region, rank, predictor) grid is unique — otherwise the
    headline pivots/reindexes would see duplicate (region, rank) labels."""
    if "phylum" in df.columns:
        ph = df["phylum"].astype("string").fillna("").str.strip()
        df = df[ph == ""]
    out = pd.DataFrame()
    out["region"] = df[_pick(df, "region", required=True)].astype(str)
    out["predictor"] = df[_pick(df, "predictor", default="anchor")].astype(str) \
        if _pick(df, "predictor") else "anchor"
    out["rank"] = df[_pick(df, "rank", required=True)].astype(str)
    out["concordance"] = pd.to_numeric(df[_concordance_col(df)], errors="coerce")
    lo, hi = _ci_cols(df)
    out["ci_lo"] = pd.to_numeric(df[lo], errors="coerce") if lo else np.nan
    out["ci_hi"] = pd.to_numeric(df[hi], errors="coerce") if hi else np.nan
    cov = _pick(df, "coverage")
    out["coverage"] = pd.to_numeric(df[cov], errors="coerce") if cov else np.nan
    # carry any self-mode column through unchanged so callers can split on it.
    smc = _selfmode_col(df)
    if smc is not None:
        out[smc] = df[smc].values
    return out


def _is_long_grid(df: pd.DataFrame) -> bool:
    """True iff ``df`` is the long (region, rank, predictor, scalar-concordance)
    grid rather than the per-amplicon ``asv_region_concordance`` table (whose
    concordance lives in per-rank ``correct_<rank>`` columns, not a scalar)."""
    if "rank" not in df.columns:
        return False
    return _pick(df, "correct", "concordance", "correct_frac", "accuracy") is not None


def load_concordance_grid(fig: str):
    """Load the long per-rank concordance grid for F1/F2/F11.

    The grid is produced by S9 in ``concordance_long.csv`` (real schema: a
    scalar ``correct`` per (region, domain, tier, predictor, rank)).  The
    fixture/legacy harness instead ships it inside ``asv_region_concordance.csv``
    as (region, rank, predictor, self_level, domain, concordance).  Try the
    canonical file first, then fall back to the asv table only when *that* file
    is shaped like the long grid.  Returns ``(raw_df, source_path)`` or
    ``(None, None)``."""
    raw = _load(PATHS.concordance_long_csv, fig, quiet=True)
    if raw is not None and _is_long_grid(raw):
        return raw, PATHS.concordance_long_csv
    alt = _load(PATHS.asv_region_concordance_csv, fig, quiet=True)
    if alt is not None and _is_long_grid(alt):
        return alt, PATHS.asv_region_concordance_csv
    if raw is not None:  # present but not a long grid -> let caller report
        return raw, PATHS.concordance_long_csv
    if alt is not None:
        return alt, PATHS.asv_region_concordance_csv
    print(f"  [{fig}] SKIP — no concordance grid found "
          f"({os.path.basename(PATHS.concordance_long_csv)} / "
          f"{os.path.basename(PATHS.asv_region_concordance_csv)})")
    return None, None


# =========================================================================== #
# F1 — Accuracy heatmap region x rank  (anchor / consensus / coverage)
# =========================================================================== #
def render_F1(force: bool) -> None:
    stem = "F1_accuracy_heatmap"
    if not force and _outputs_exist(stem):
        print(f"  [F1] skip (exists) — {stem}")
        return
    raw, _src = load_concordance_grid("F1")
    if raw is None:
        return
    if not _have_cols(raw, ["region", "rank"], "F1"):
        return
    raw = _rollup(_rollup(_filter_selfmode(raw, "include"), "domain", "all"),
                  "tier", "all")
    df = _headline_concordance(raw)

    regions = _regions_present(df["region"])
    if not regions:
        print("  [F1] SKIP — no regions present")
        return
    # x = rank in genus->phylum reading order (deep first), the F1 convention.
    ranks = [r for r in reversed(SCORED) if r in set(df["rank"])]
    if not ranks:
        print("  [F1] SKIP — no SCORED ranks present")
        return

    primary, secondary = _resolve_predictors(df["predictor"])
    if primary is None:
        print("  [F1] SKIP — no predictor rows")
        return
    panels = [(primary, "concordance", f"{primary} concordance", (0.0, 1.0),
               "viridis")]
    if secondary is not None:
        panels.append((secondary, "concordance", f"{secondary} concordance",
                       (0.0, 1.0), "viridis"))
    # Final panel: primary coverage (the (correct, coverage) pair, DESIGN §6) —
    # only when a real coverage column exists (the fixture grid omits it).
    if df["coverage"].notna().any():
        panels.append((primary, "coverage", f"{primary} coverage", (0.0, 1.0),
                       "magma"))

    fig, axes = plt.subplots(
        1, len(panels), figsize=(4.6 * len(panels), 0.55 * len(regions) + 2.4),
        squeeze=False, constrained_layout=True)
    axes = axes[0]
    last_im = {}
    for ax, (pred, value, lab, (vmin, vmax), cmap) in zip(axes, panels):
        sub = df[df["predictor"] == pred]
        grid = (sub.pivot_table(index="region", columns="rank", values=value,
                                aggfunc="mean")
                .reindex(index=regions, columns=ranks))
        M = grid.values.astype(float)
        im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=vmin, vmax=vmax,
                       aspect="auto")
        last_im[cmap] = im
        ax.set_xticks(range(len(ranks)))
        ax.set_xticklabels([RANK_SHORT.get(r, r) for r in ranks])
        ax.set_yticks(range(len(regions)))
        ax.set_yticklabels(regions if ax is axes[0] else [])
        ax.set_title(lab, fontsize=10)
        ax.set_xlabel("rank")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = M[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                            fontsize=7, color="white" if v < 0.55 else "black")
    axes[0].set_ylabel("region (short → long amplicon)")
    fig.suptitle("F1 · Taxonomic concordance & coverage by region × rank "
                 "(include-self)", fontsize=12)
    # one colorbar per distinct colormap used
    for cmap, im in last_im.items():
        cb = fig.colorbar(im, ax=list(axes), fraction=0.02, pad=0.01)
        cb.set_label("coverage" if cmap == "magma" else "concordance")
    _save(fig, stem)


# =========================================================================== #
# F2 — Per-rank concordance by region (grouped lines + markers + Wilson CI)
# =========================================================================== #
def render_F2(force: bool) -> None:
    stem = "F2_perrank_concordance"
    if not force and _outputs_exist(stem):
        print(f"  [F2] skip (exists) — {stem}")
        return
    raw, _src = load_concordance_grid("F2")
    if raw is None:
        return
    if not _have_cols(raw, ["region", "rank"], "F2"):
        return
    raw = _rollup(_rollup(_filter_selfmode(raw, "include"), "domain", "all"),
                  "tier", "all")
    df = _headline_concordance(raw)
    primary, _sec = _resolve_predictors(df["predictor"])
    if primary is not None:
        df = df[df["predictor"] == primary]

    regions = _regions_present(df["region"])
    ranks = [r for r in SCORED if r in set(df["rank"])]  # phylum->species l-to-r
    if not regions or not ranks:
        print("  [F2] SKIP — no regions/ranks present")
        return

    fig, ax = plt.subplots(figsize=(8.6, 5.2), constrained_layout=True)
    x = np.arange(len(ranks))
    for region in regions:
        sub = df[df["region"] == region].set_index("rank").reindex(ranks)
        y = sub["concordance"].values.astype(float)
        color = _color(region)
        ax.plot(x, y, marker="o", lw=1.8, ms=5, color=color, label=region)
        lo = sub["ci_lo"].values.astype(float)
        hi = sub["ci_hi"].values.astype(float)
        if np.isfinite(lo).any() and np.isfinite(hi).any():
            ax.fill_between(x, lo, hi, color=color, alpha=0.13, lw=0)
    ax.set_xticks(x)
    ax.set_xticklabels(ranks, rotation=30, ha="right")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel(f"{primary or 'predictor'} concordance (correct / called)")
    ax.set_xlabel("taxonomic rank")
    ax.grid(True, axis="y", ls=":", alpha=0.5)
    ax.set_title(f"F2 · Per-rank concordance by region "
                 f"({primary or 'predictor'}, include-self)")
    ax.legend(title="region (short→long)", fontsize=8, ncol=2,
              loc="lower left", framealpha=0.9)
    _save(fig, stem)


# =========================================================================== #
# F3 — Deepest-correct-rank distribution (100%-stacked horizontal bar)
# =========================================================================== #
def render_F3(force: bool) -> None:
    stem = "F3_deepest_correct_stack"
    if not force and _outputs_exist(stem):
        print(f"  [F3] skip (exists) — {stem}")
        return
    df = _load(PATHS.asv_region_concordance_csv, "F3")
    if df is None:
        return
    if not _have_cols(df, ["region", "deepest_correct_rank"], "F3"):
        return
    df = _filter_selfmode(df, "include")
    df = _rollup(df, "domain", "all")  # no-op (per-amplicon rows carry a real domain)
    df = _dedup_predictor(df)  # one row per amplicon

    df = df.copy()
    dc = df["deepest_correct_rank"].astype("string").fillna("").str.strip()
    df["deepest"] = dc.where(dc != "", "none")

    regions = _regions_present(df["region"])
    if not regions:
        print("  [F3] SKIP — no regions present")
        return
    seg_labels = [s for s in DEPTH_LABELS if s in set(df["deepest"])]
    if not seg_labels:
        print("  [F3] SKIP — no deepest-correct segments present")
        return

    counts = (df.groupby(["region", "deepest"]).size()
              .unstack(fill_value=0)
              .reindex(index=regions, columns=seg_labels, fill_value=0))
    rowsum = counts.sum(axis=1).replace(0, np.nan)
    grid = counts.div(rowsum, axis=0).fillna(0.0)

    fig, ax = plt.subplots(figsize=(8.8, 0.5 * len(regions) + 1.9),
                           constrained_layout=True)
    y = np.arange(len(regions))
    left = np.zeros(len(regions))
    for seg in seg_labels:
        w = grid[seg].values.astype(float)
        ax.barh(y, w, left=left, color=DEPTH_RAMP.get(seg, "0.6"),
                edgecolor="white", lw=0.5, label=seg)
        left += w
    ax.set_yticks(y)
    ax.set_yticklabels(regions)
    ax.invert_yaxis()  # short amplicon at top
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of scored amplicons")
    ax.set_ylabel("region (short → long amplicon)")
    ax.set_title("F3 · Deepest correct rank by region (anchor, include-self)")
    ax.legend(title="deepest correct rank", ncol=len(seg_labels), fontsize=8,
              loc="upper center", bbox_to_anchor=(0.5, -0.12))
    _save(fig, stem)


# =========================================================================== #
# F6 — Self-recovery rate by region (bar + Wilson CI, top-1 vs any-top-20)
# =========================================================================== #
def render_F6(force: bool) -> None:
    stem = "F6_self_recovery"
    if not force and _outputs_exist(stem):
        print(f"  [F6] skip (exists) — {stem}")
        return
    df = _load(PATHS.self_recovery_csv, "F6")
    if df is None:
        return
    ncol = _pick(df, "n_self_in_top20", "self_in_top20")
    if ncol is None or not _have_cols(df, ["region"], "F6"):
        print("  [F6] SKIP — no n_self_in_top20 column")
        return
    df = _rollup(_rollup(df, "domain", "all"), "tier", "all")

    regions = _regions_present(df["region"])
    if not regions:
        print("  [F6] SKIP — no regions present")
        return
    df = df.drop_duplicates(subset=["region"]).set_index("region")

    # recovered-given-available (DESIGN §6): of the amplicons whose source was in
    # top20, how often was it (a) the anchor top-1, (b) selected at all (any-20)?
    rate_top1 = _pick(df, "anchor_given_top20_rate", "rate_top1", "self_anchor_rate")
    rate_any = _pick(df, "selected_given_top20_rate", "rate_any20",
                     "self_selected_rate")
    num_top1 = _pick(df, "n_anchor_given_top20", "n_self_anchor")
    num_any = _pick(df, "n_selected_given_top20", "n_self_selected")

    def _bar_vals(rate_col, num_col):
        vals, elo, ehi = [], [], []
        for region in regions:
            n = float(df.loc[region, ncol])
            v = float(df.loc[region, rate_col]) if rate_col and np.isfinite(
                df.loc[region, rate_col]) else np.nan
            if not np.isfinite(v):
                vals.append(np.nan); elo.append(0.0); ehi.append(0.0); continue
            if num_col and num_col in df.columns and np.isfinite(df.loc[region, num_col]):
                k = int(round(float(df.loc[region, num_col])))
            else:
                k = int(round(v * n))
            lo, hi = C.wilson_ci(k, int(round(n))) if n > 0 else (v, v)
            vals.append(v); elo.append(max(v - lo, 0.0)); ehi.append(max(hi - v, 0.0))
        return vals, elo, ehi

    fig, ax = plt.subplots(figsize=(8.8, 5.2), constrained_layout=True)
    x = np.arange(len(regions))
    w = 0.38
    colors = [_color(r) for r in regions]
    series = [
        (rate_top1, num_top1, "top-1 (anchor | in top20)", "", 0.95, -0.5),
        (rate_any, num_any, "any-top-20 (selected | in top20)", "//", 0.62, 0.5),
    ]
    for rate_col, num_col, lab, hatch, alpha, off in series:
        if rate_col is None:
            continue
        vals, elo, ehi = _bar_vals(rate_col, num_col)
        ax.bar(x + off * w, vals, w, yerr=[elo, ehi], color=colors,
               edgecolor="black", lw=0.5, hatch=hatch, capsize=2.5,
               alpha=alpha, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=30, ha="right")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("fraction of self-in-top20 recovered")
    ax.set_title("F6 · Self-recovery rate by region (include-self round-trip)")
    ax.grid(True, axis="y", ls=":", alpha=0.5)
    legend_handles = [
        Patch(facecolor="0.6", edgecolor="black", label="top-1 (anchor)"),
        Patch(facecolor="0.6", edgecolor="black", hatch="//", alpha=0.62,
              label="any-top-20 (selected)"),
    ]
    ax.legend(handles=legend_handles, fontsize=9, loc="lower right")
    _save(fig, stem)


# =========================================================================== #
# F7 — Genus-failure confusion (truth x pred genus heatmap, LogNorm) per region
# =========================================================================== #
def render_F7(force: bool) -> None:
    stem_base = "F7_genus_confusion"
    df = _load(PATHS.confusion_genus_csv, "F7")
    if df is None:
        return
    # Producer emits (true_genus, pred_genus); accept the family alias too.
    tcol = _pick(df, "true_genus", "true_family")
    pcol = _pick(df, "pred_genus", "pred_family")
    ccol = _pick(df, "count", "n")
    if not (tcol and pcol and ccol and "region" in df.columns):
        print(f"  [F7] SKIP — confusion CSV missing truth/pred/count/region "
              f"(have {list(df.columns)})")
        return
    df = _filter_selfmode(df, "include")
    df = _rollup(df, "domain", "all")

    regions = _regions_present(df["region"])
    if not regions:
        print("  [F7] SKIP — no regions present")
        return
    label = "genus" if tcol.endswith("genus") else "family"

    wrote_any = False
    for region in regions:
        stem = f"{stem_base}_{region.replace('/', '-')}"
        if not force and _outputs_exist(stem):
            print(f"  [F7] skip (exists) — {stem}")
            wrote_any = True
            continue
        sub = df[df["region"] == region]
        top_true = (sub.groupby(tcol)[ccol].sum()
                    .sort_values(ascending=False).head(30).index.tolist())
        top_pred = (sub.groupby(pcol)[ccol].sum()
                    .sort_values(ascending=False).head(30).index.tolist())
        if not top_true or not top_pred:
            print(f"  [F7] {region}: no confusion cells, skip")
            continue
        grid = (sub.pivot_table(index=tcol, columns=pcol, values=ccol,
                                aggfunc="sum")
                .reindex(index=top_true, columns=top_pred).fillna(0.0))
        M = grid.values.astype(float)
        vmax = M.max() if M.max() > 0 else 1.0
        fig, ax = plt.subplots(
            figsize=(max(6, 0.32 * len(top_pred) + 3),
                     max(5, 0.32 * len(top_true) + 2.5)),
            constrained_layout=True)
        cmap = plt.get_cmap("magma").copy()
        cmap.set_bad("white")
        im = ax.imshow(np.ma.masked_less(M, 1.0), cmap=cmap,
                       norm=LogNorm(vmin=1, vmax=max(vmax, 1.0)), aspect="auto")
        ax.set_xticks(range(len(top_pred)))
        ax.set_xticklabels([_short_phylum(p) for p in top_pred],
                           rotation=80, ha="right", fontsize=6)
        ax.set_yticks(range(len(top_true)))
        ax.set_yticklabels([_short_phylum(p) for p in top_true], fontsize=6)
        ax.set_xlabel(f"predicted {label}")
        ax.set_ylabel(f"true {label}")
        ax.set_title(f"F7 · {label.capitalize()}-failure confusion — {region} "
                     f"(anchor, include-self)", fontsize=10)
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label("misassignment count (log)")
        _save(fig, stem)
        wrote_any = True
    if not wrote_any:
        print("  [F7] nothing rendered")


# =========================================================================== #
# F10 — Per-phylum accuracy (phylum x region heatmap, n<10 masked)
# =========================================================================== #
def render_F10(force: bool) -> None:
    stem = "F10_perphylum_accuracy"
    if not force and _outputs_exist(stem):
        print(f"  [F10] skip (exists) — {stem}")
        return
    # F10 needs the per-phylum Genus breakdown, which may live in either
    # concordance_long.csv (canonical) or the long-grid asv table (fixture).
    df = None
    for path in (PATHS.concordance_long_csv, PATHS.asv_region_concordance_csv):
        cand = _load(path, "F10") if C.done(path) else None
        if cand is not None and "phylum" in cand.columns and _is_long_grid(cand):
            df = cand
            break
    if df is None:
        print("  [F10] SKIP — no long concordance table carries a 'phylum' column "
              "(the per-phylum Genus breakdown is not produced by score_taxacc.py; "
              "add (region, phylum, rank, predictor) Genus rows upstream to enable "
              "F10)")
        return
    df = _filter_selfmode(df, "include")
    df = df[df["rank"].astype(str) == "Genus"]
    df = _dedup_predictor(df)
    df = df[df["phylum"].notna() & (df["phylum"].astype(str).str.strip() != "")]
    if df.empty:
        print("  [F10] SKIP — no Genus-rank per-phylum rows")
        return
    conc_col = _concordance_col(df)
    nev_col = _pick(df, "n_called", "n_evaluable", "n")

    regions = _regions_present(df["region"])
    conc = (df.pivot_table(index="phylum", columns="region", values=conc_col,
                           aggfunc="mean").reindex(columns=regions))
    nev = (df.pivot_table(index="phylum", columns="region",
                          values=nev_col, aggfunc="sum").reindex(columns=regions)
           if nev_col else None)
    M = conc.values.astype(float)
    if nev is not None:
        Nm = nev.values.astype(float)
        M = np.where(np.isfinite(Nm) & (Nm >= 10), M, np.nan)  # n<10 mask (F10)

    phyla = list(conc.index)
    mean_acc = np.where(np.all(np.isnan(M), axis=1), -1.0, np.nanmean(M, axis=1))
    arch = sorted([(p, a) for p, a in zip(phyla, mean_acc) if _is_archaea_phylum(p)],
                  key=lambda t: (-t[1], t[0]))
    bact = sorted([(p, a) for p, a in zip(phyla, mean_acc) if not _is_archaea_phylum(p)],
                  key=lambda t: (-t[1], t[0]))
    order = [p for p, _ in arch] + [p for p, _ in bact]
    idx_map = {p: i for i, p in enumerate(phyla)}
    M = M[[idx_map[p] for p in order], :]
    keep = ~np.all(np.isnan(M), axis=1)
    order = [p for p, k in zip(order, keep) if k]
    M = M[keep, :]
    if M.size == 0:
        print("  [F10] SKIP — every (phylum,region) cell has n<10")
        return

    fig, ax = plt.subplots(
        figsize=(0.55 * len(regions) + 3, 0.28 * len(order) + 2.5),
        constrained_layout=True)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("lightgrey")
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=1,
                   aspect="auto")
    ax.set_xticks(range(len(regions)))
    ax.set_xticklabels(regions, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([_short_phylum(p) for p in order], fontsize=6)
    n_arch = sum(1 for p in order if _is_archaea_phylum(p))
    if 0 < n_arch < len(order):
        ax.axhline(n_arch - 0.5, color="red", lw=1.2, ls="--", alpha=0.7)
    ax.set_xlabel("region (short → long amplicon)")
    ax.set_ylabel("phylum  (Archaea above red line)")
    ax.set_title("F10 · Per-phylum genus concordance (n≥10; anchor, include-self)",
                 fontsize=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("genus concordance")
    _save(fig, stem)


# =========================================================================== #
# F11 — Self-inclusion vs self-exclusion concordance delta per region × rank
# =========================================================================== #
def render_F11(force: bool) -> None:
    stem = "F11_self_inclusion_exclusion_delta"
    if not force and _outputs_exist(stem):
        print(f"  [F11] skip (exists) — {stem}")
        return
    raw, src = load_concordance_grid("F11")
    if raw is None:
        return
    if not _have_cols(raw, ["region", "rank"], "F11"):
        return
    col = _selfmode_col(raw)
    if col is None:
        print(f"  [F11] SKIP — {os.path.basename(src)} carries no self-mode column "
              f"({SELFMODE_COLS}); the current pipeline is include-self only "
              "(DESIGN USER DECISION 4 — pipeline unaltered). No leakage delta to "
              "render until an exclude-self branch (e.g. self_mode='exclude_0.987') "
              "is produced by the scorer.")
        return

    base = _rollup(_rollup(raw, "domain", "all"), "tier", "all")
    inc = _headline_concordance(_filter_selfmode(base, "include"))
    exc = _headline_concordance(_filter_selfmode(base, "exclude"))
    if inc.empty or exc.empty:
        print("  [F11] SKIP — need BOTH include and exclude self-mode rows "
              f"(include n={len(inc)}, exclude n={len(exc)})")
        return
    primary, _sec = _resolve_predictors(set(inc["predictor"]) & set(exc["predictor"])
                                        or inc["predictor"])
    if primary is not None:
        inc = inc[inc["predictor"] == primary]
        exc = exc[exc["predictor"] == primary]

    merged = inc.merge(exc, on=["region", "rank"], suffixes=("_inc", "_exc"))
    merged["delta"] = merged["concordance_inc"] - merged["concordance_exc"]
    regions = _regions_present(merged["region"])
    ranks = [r for r in SCORED if r in set(merged["rank"])]
    if not regions or not ranks:
        print("  [F11] SKIP — no regions/ranks after include∩exclude merge")
        return

    grid = (merged.pivot_table(index="region", columns="rank", values="delta",
                               aggfunc="mean")
            .reindex(index=regions, columns=ranks))
    M = grid.values.astype(float)
    vmax = np.nanmax(np.abs(M)) if np.isfinite(M).any() else 1.0
    vmax = max(vmax, 1e-3)

    fig, ax = plt.subplots(
        figsize=(0.9 * len(ranks) + 3, 0.55 * len(regions) + 2.4),
        constrained_layout=True)
    im = ax.imshow(np.ma.masked_invalid(M), cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="auto")
    ax.set_xticks(range(len(ranks)))
    ax.set_xticklabels(ranks, rotation=30, ha="right")
    ax.set_yticks(range(len(regions)))
    ax.set_yticklabels(regions)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7,
                        color="black" if abs(v) < 0.6 * vmax else "white")
    ax.set_xlabel("taxonomic rank")
    ax.set_ylabel("region (short → long amplicon)")
    ax.set_title(f"F11 · Include−exclude self concordance Δ (leakage / "
                 f"generalization gap; {primary or 'predictor'})", fontsize=11)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Δ concordance (include − exclude)")
    _save(fig, stem)


# =========================================================================== #
# F13 — Amplification-rate / selection-bias (region x phylum, Bact/Arch)
# =========================================================================== #
def render_F13(force: bool) -> None:
    stem = "F13_amplification_by_phylum"
    if not force and _outputs_exist(stem):
        print(f"  [F13] skip (exists) — {stem}")
        return
    df = _load(PATHS.amplification_by_phylum_csv, "F13")
    if df is None:
        return
    if not _have_cols(df, ["region", "phylum", "extract_rate"], "F13"):
        return
    df = df.copy()
    df["extract_rate"] = pd.to_numeric(df["extract_rate"], errors="coerce")

    regions = _regions_present(df["region"])
    if not regions:
        print("  [F13] SKIP — no regions present")
        return
    # Domain facets: Bacteria then Archaea (DESIGN F13 Bact-vs-Arch).
    if "domain" in df.columns:
        present = set(df["domain"].astype(str))
        domains = [d for d in ("Bacteria", "Archaea") if d in present]
        domains += sorted(d for d in present if d not in ("Bacteria", "Archaea"))
    else:
        domains = ["all"]

    weight_col = _pick(df, "n_source", "n_attempted")
    if weight_col:
        top_phyla = (df.groupby("phylum")[weight_col].sum()
                     .sort_values(ascending=False).head(15).index.tolist())
    else:
        top_phyla = df["phylum"].value_counts().head(15).index.tolist()
    top_phyla = [p for p in top_phyla if isinstance(p, str) and p]
    if not top_phyla:
        print("  [F13] SKIP — no phyla present")
        return

    fig, axes = plt.subplots(
        len(domains), 1,
        figsize=(max(9, 0.7 * len(top_phyla) + 3), 3.6 * len(domains) + 0.6),
        squeeze=False, constrained_layout=True)
    axes = axes[:, 0]
    nreg = len(regions)
    bw = 0.8 / max(nreg, 1)
    for ax, dom in zip(axes, domains):
        sub = df if dom == "all" else df[df["domain"].astype(str) == dom]
        grid = (sub.pivot_table(index="phylum", columns="region",
                                values="extract_rate", aggfunc="mean")
                .reindex(index=top_phyla, columns=regions))
        x = np.arange(len(top_phyla))
        for j, region in enumerate(regions):
            vals = grid[region].values.astype(float) if region in grid.columns \
                else np.full(len(top_phyla), np.nan)
            ax.bar(x + (j - (nreg - 1) / 2) * bw, np.nan_to_num(vals, nan=0.0),
                   bw, color=_color(region), edgecolor="none",
                   label=region if ax is axes[0] else None)
        ax.set_xticks(x)
        ax.set_xticklabels([_short_phylum(p) for p in top_phyla],
                           rotation=70, ha="right", fontsize=7)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("extract rate")
        ax.set_title(f"{dom}", fontsize=10, loc="left")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
    axes[0].legend(title="region (short→long)", fontsize=7, ncol=2,
                   loc="upper right", framealpha=0.9)
    fig.suptitle("F13 · Amplification rate by phylum × region "
                 "(survivorship-bias guard)", fontsize=12)
    _save(fig, stem)


# =========================================================================== #
# F14 — Tie/ambiguity cluster size (per-amplicon ECDF per region)
# =========================================================================== #
def render_F14(force: bool) -> None:
    stem = "F14_ambiguity_ecdf"
    if not force and _outputs_exist(stem):
        print(f"  [F14] skip (exists) — {stem}")
        return
    # Per-amplicon tie-cluster sizes live in asv_region_concordance.csv (the
    # ambiguity.csv producer only emits per-(region,domain,tier) summaries, which
    # cannot draw an ECDF).  Fall back to ambiguity.csv if the per-amplicon table
    # is absent (then plot the available summary points instead of an ECDF).
    df = _load(PATHS.asv_region_concordance_csv, "F14")
    source = "per_amplicon"
    if df is None or "tie_cluster_size" not in df.columns:
        df = _load(PATHS.ambiguity_csv, "F14")
        source = "summary"
        if df is None:
            return
        if not _have_cols(df, ["region"], "F14"):
            return

    if source == "per_amplicon":
        df = _filter_selfmode(df, "include")
        df = _dedup_predictor(df)  # one row per amplicon
        df = df.copy()
        df["tie_cluster_size"] = pd.to_numeric(df["tie_cluster_size"],
                                               errors="coerce")
        regions = _regions_present(df["region"])
        if not regions:
            print("  [F14] SKIP — no regions present")
            return
        fig, ax = plt.subplots(figsize=(8.6, 5.4), constrained_layout=True)
        drew = False
        for region in regions:
            vals = (df.loc[df["region"] == region, "tie_cluster_size"]
                    .dropna().astype(float).values)
            vals = vals[vals >= 1]
            if vals.size == 0:
                continue
            xs = np.sort(vals)
            ys = np.arange(1, xs.size + 1) / xs.size
            ax.step(np.concatenate([[xs[0]], xs]),
                    np.concatenate([[0.0], ys]), where="post",
                    color=_color(region), lw=1.8, label=region)
            drew = True
        if not drew:
            print("  [F14] SKIP — no tie-cluster sizes to plot")
            plt.close(fig)
            return
        ax.set_xscale("log")
        ax.set_xlabel("tie-cluster size  (#refs within ε of best identity)")
        ax.set_ylabel("ECDF (fraction of amplicons ≤ x)")
        ax.set_ylim(0, 1.02)
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.set_title("F14 · Tie/ambiguity cluster size by region (include-self)")
        ax.legend(title="region (short→long)", fontsize=8, ncol=2,
                  loc="lower right", framealpha=0.9)
        _save(fig, stem)
        return

    # --- summary fallback: bar of median (+ p90 whisker) tie size per region ---
    med_col = _pick(df, "tie_median", "median")
    p90_col = _pick(df, "tie_p90", "p90")
    if med_col is None:
        print("  [F14] SKIP — ambiguity summary lacks a tie_median column")
        return
    df = _rollup(_rollup(df, "domain", "all"), "tier", "all")
    regions = _regions_present(df["region"])
    df = df.drop_duplicates(subset=["region"]).set_index("region")
    fig, ax = plt.subplots(figsize=(8.6, 5.0), constrained_layout=True)
    x = np.arange(len(regions))
    med = [float(df.loc[r, med_col]) if r in df.index else np.nan for r in regions]
    ax.bar(x, med, color=[_color(r) for r in regions], edgecolor="black", lw=0.5)
    if p90_col:
        p90 = [float(df.loc[r, p90_col]) if r in df.index else np.nan for r in regions]
        err = [max(p - m, 0.0) if np.isfinite(p) and np.isfinite(m) else 0.0
               for p, m in zip(p90, med)]
        ax.errorbar(x, med, yerr=[np.zeros(len(err)), err], fmt="none",
                    ecolor="0.2", capsize=3, elinewidth=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=30, ha="right")
    ax.set_ylabel("tie-cluster size (median; whisker→p90)")
    ax.set_title("F14 · Tie/ambiguity cluster size by region "
                 "(summary fallback; include-self)")
    ax.grid(True, axis="y", ls=":", alpha=0.4)
    _save(fig, stem)


# =========================================================================== #
# Dispatch
# =========================================================================== #
FIGURES = {
    "F1": render_F1, "F2": render_F2, "F3": render_F3, "F6": render_F6,
    "F7": render_F7, "F10": render_F10, "F11": render_F11,
    "F13": render_F13, "F14": render_F14,
}


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Render accuracy figures F1/F2/F3/F6/F7/F10/F11/F13/F14 from "
                    "the flat CSVs in _common.PATHS.")
    ap.add_argument("--force", action="store_true",
                    help="re-render even if a figure's .png/.pdf already exist")
    ap.add_argument("--only", default="",
                    help="comma-separated subset, e.g. 'F1,F3,F14' (default: all)")
    args = ap.parse_args()

    C.ensure_dirs()
    if args.only.strip():
        want = [f.strip().upper() for f in args.only.split(",") if f.strip()]
        unknown = [f for f in want if f not in FIGURES]
        if unknown:
            ap.error(f"unknown figure(s) {unknown}; choose from {list(FIGURES)}")
    else:
        want = list(FIGURES)

    print(f"render_accuracy → {C.FIGS}")
    print(f"  REGION_ORDER (short→long): {REGION_ORDER}")
    print(f"  figures: {want}  | primary predictor: {PRIMARY_PREDICTOR}"
          f"{'  | --force' if args.force else ''}")
    for fig in want:
        print(f"[{fig}]")
        try:
            FIGURES[fig](args.force)
        except Exception as exc:
            import traceback
            print(f"  [{fig}] ERROR — {exc}")
            traceback.print_exc()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
