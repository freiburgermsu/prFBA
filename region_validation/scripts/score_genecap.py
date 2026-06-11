#!/usr/bin/env python
"""S10 -- gene-capture scoring (DESIGN.md Section 7), reporting BOTH self_modes.

DESIGN decision #4: the pipeline RUN is never altered (the source organism stays
in the reference DB); from that single run we report, per ``(genome, region)``
cell, TWO numbers under two ``self_mode`` values:

  * ``include``        -- realistic upper bound.  P is the union of PGFam sets of
                          the AS-IS selected reps (``PATHS.selection_out``); the
                          source genome, if among the reps, is kept in P.
  * ``exclude_0.987``  -- honest generalization (novel-organism).  P is rebuilt
                          from a SELF-EXCLUSION RE-SELECTION: drop top-hits whose
                          ``genome_id == src_genome_id`` OR ``taxon_id ==
                          src_taxon_id`` OR ``identity >= tau`` (tau = 0.987, the
                          species threshold), re-rank 1..n, recompute
                          ``best_identity`` (0.0 if nothing survives -> clean
                          abstain, never None), and re-run the PURE
                          ``select_references.select_representatives`` on the
                          filtered record.  P is the union of those reps' PGFams.

Every cell / long-form row carries a ``self_mode`` column, so the include->exclude
drop is the leakage / generalization gap (figure F11 territory).

For each (genome, region, self_mode) cell:

    Truth  T = PGF(src_genome_id)                       -- the source's PGFam set
    Pred   P = union of PGF(r) over r in selected reps  -- per self_mode (above)

with ``t=|T|``, ``p=|P|``, ``i=|T & P|``:

    recall    = i/t   (= "% gene capture", reported as 100*i/t)
    precision = i/p
    F1        = 2i/(t+p)
    Jaccard   = i/(t+p-i)

PGFam sets for BOTH truth and predicted reps come from the *one* namespace
``select_references.bvbrc_gene_provider`` over ``PATHS.pgfam_cache`` so the two
sides are apples-to-apples.  This stage is a *pure post* step: it never hits the
network.  A rep genome_id absent from the cache is classified ``pgfam_missing``
and poisons the whole cell (NOT silently fetched, NOT silently dropped from P --
which would inflate precision, verifier FINDING 6).

Edge outcomes (DESIGN Section 7, every cell classified, never silently dropped):
  * ``truth_empty``    -- T = empty (16S-only / no usable CDS source) -> all NaN,
                          excluded from the metric means.
  * ``abstained``      -- n_selected == 0 -> recall/F1/Jaccard = 0, precision NaN.
  * ``pred_empty``     -- every selected rep has an empty PGFam set -> recall 0.
  * ``pgfam_missing``  -- the source's set OR any selected rep's set is missing
                          from the cache -> poison the cell (flagged, excluded
                          from every denominator).
  * ``ok``             -- a scorable cell.

The exclude-self re-selection needs the per-amplicon top hits, so when both
``PATHS.hits_json`` and a usable selection record exist, ``exclude_0.987`` is
computed; otherwise that mode is recorded as ``no_hits`` (poisoned, excluded)
and ``include`` is still fully reported.

Stratification: by selection ``tier`` (species/genus/family/below) and by
best-non-self-hit identity bins (the recall-vs-identity axis = distance to the
nearest *genuinely different* reference; from ``PATHS.hits_json`` when present).

Resumability: if all three outputs exist and are non-empty, short-circuit unless
``--force``.

I/O contract
------------
Reads:
  * PATHS.selection_out  (selection.json: _meta + ampliconKey -> selection record)
  * PATHS.pgfam_cache    (genome_gene_families.json: genome_id -> [pgfam_id,...])
  * PATHS.hits_json      (asv_top20_alignment_hits.json: ampliconKey -> top20;
                          REQUIRED for the exclude_0.987 re-selection + non-self bins)
Writes:
  * PATHS.genecap_cells_csv      (one row per (genome,region,self_mode) cell)
  * PATHS.gene_capture_long_csv  (tidy long form for the renderers)
  * PATHS.genecap_summary_json   (macro + micro recall/precision/F1/Jaccard per
                                  region x self_mode, with Wilson CIs + per-outcome
                                  denominators)

Interpreter: /home/freiburger/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys

# Make `import _common` work regardless of cwd (agent cwd resets between calls).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common as C  # noqa: E402
from _common import PATHS, wilson_ci  # noqa: E402

# select_references provides the single gene-provider namespace shared by truth
# and predicted reps (apples-to-apples) AND the PURE select_representatives used
# for the exclude-self re-selection.  prFBA root must be importable.
sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import select_references as SR  # noqa: E402


# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
GENOME_ID_RE = re.compile(r"^\d+\.\d+$")
TIER_ORDER = ["species", "genus", "family", "below"]

# Self-exclusion threshold for the honest "exclude" mode (species threshold).
TAU = 0.987

# The two self_modes reported on every cell (DESIGN decision #4).
MODE_INCLUDE = "include"
MODE_EXCLUDE = "exclude_0.987"
SELF_MODES = [MODE_INCLUDE, MODE_EXCLUDE]

# Best-non-self-hit identity bins (distance to nearest genuinely-different ref).
# Anchored on the Yarza identity floors used elsewhere in the panel.
IDENT_BIN_EDGES = [0.0, 0.865, 0.945, 0.987, 0.995, 1.0001]
IDENT_BIN_LABELS = ["<0.865", "0.865-0.945", "0.945-0.987", "0.987-0.995", ">=0.995"]
IDENT_BIN_UNKNOWN = "unknown"  # hits unavailable / no non-self hit

# Edge outcomes (precedence order in score_cell()).
OUTCOME_TRUTH_EMPTY = "truth_empty"
OUTCOME_PGFAM_MISSING = "pgfam_missing"
OUTCOME_ABSTAINED = "abstained"
OUTCOME_PRED_EMPTY = "pred_empty"
OUTCOME_NO_HITS = "no_hits"   # exclude mode could not be re-selected (no hit record)
OUTCOME_OK = "ok"
OUTCOMES = [OUTCOME_OK, OUTCOME_ABSTAINED, OUTCOME_PRED_EMPTY,
            OUTCOME_TRUTH_EMPTY, OUTCOME_PGFAM_MISSING, OUTCOME_NO_HITS]

METRICS = ["recall", "precision", "f1", "jaccard"]


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #
def load_selection(path):
    """Load selection.json; return {ampliconKey: record} stripped of _meta."""
    with open(path) as fh:
        obj = json.load(fh)
    return {k: v for k, v in obj.items() if k != "_meta"}


def load_cache_keys(path):
    """Set of genome_ids present in the PGFam cache.

    Used to classify cache-absent genomes as ``pgfam_missing`` *without* a live
    fetch (this stage is pure post).  Supports a ``.json.gz`` sibling the same
    way ``bvbrc_gene_provider`` does.
    """
    p = path if os.path.exists(path) else (path + ".gz")
    if not os.path.exists(p):
        return set()
    if p.endswith(".gz"):
        import gzip
        with gzip.open(p, "rt") as fh:
            return set(json.load(fh).keys())
    with open(p) as fh:
        return set(json.load(fh).keys())


def load_hits(path):
    """Load the per-amplicon top-hit records (ampliconKey -> record), or {}."""
    if not C.done(path):
        return {}
    with open(path) as fh:
        hits = json.load(fh)
    return {k: v for k, v in hits.items() if k != "_meta" and isinstance(v, dict)}


def best_nonself_identity(hit_record, src_gid, src_tid):
    """Best identity among non-self top hits of one record, or None.

    'Non-self' = a hit whose ``genome_id`` differs from the source genome_id AND
    whose ``taxon_id`` differs from the source taxon_id (organism-level).
    """
    if not hit_record:
        return None
    best = None
    for h in hit_record.get("top20", []) or []:
        hid = str(h.get("genome_id")) if h.get("genome_id") is not None else None
        htid = h.get("taxon_id")
        if hid is not None and hid == src_gid:
            continue
        if src_tid is not None and htid is not None and int(htid) == src_tid:
            continue
        ident = h.get("identity")
        if ident is None:
            continue
        if best is None or ident > best:
            best = ident
    return best


# --------------------------------------------------------------------------- #
# Source-genome helpers
# --------------------------------------------------------------------------- #
def split_key(key):
    """Split ampliconKey '{genome_id}__{region}' -> (src_genome_id, region)."""
    src_gid, _, region = key.partition("__")
    return src_gid, region


def _src_taxon_id(src_gid):
    """int(taxon_id) from a numeric BV-BRC genome_id, else None (MAG/non-numeric)."""
    if not GENOME_ID_RE.match(str(src_gid)):
        return None
    try:
        return int(str(src_gid).split(".")[0])
    except (ValueError, IndexError):
        return None


_DOMAIN_CACHE: dict = {}


def domain_for(taxon_id):
    """Resolve a taxon_id to 'Bacteria'/'Archaea' (else None).

    ``_common.lineage_for`` folds domain into ``Kingdom`` as the modern-NCBI
    realm name (e.g. 'Pseudomonadati'), which is *not* the Bacteria/Archaea split
    we stratify on.  We therefore read the literal ``domain`` rank from the same
    taxdump.  Cached.
    """
    if taxon_id is None:
        return None
    if taxon_id in _DOMAIN_CACHE:
        return _DOMAIN_CACHE[taxon_id]
    dom = None
    try:
        import taxopy
        rd = taxopy.Taxon(int(taxon_id), C._get_taxdb()).rank_name_dictionary
        dom = rd.get("domain") or rd.get("superkingdom")
    except Exception:
        dom = None
    _DOMAIN_CACHE[taxon_id] = dom
    return dom


def ident_bin(best_nonself):
    """Map best-non-self identity to its bin label (or 'unknown')."""
    if best_nonself is None:
        return IDENT_BIN_UNKNOWN
    for lab, lo, hi in zip(IDENT_BIN_LABELS, IDENT_BIN_EDGES[:-1], IDENT_BIN_EDGES[1:]):
        if lo <= best_nonself < hi:
            return lab
    return IDENT_BIN_LABELS[-1]


# --------------------------------------------------------------------------- #
# Exclude-self re-selection (DESIGN Section 5.4 -> Section 7 exclude_source)
# --------------------------------------------------------------------------- #
def exclude_self_record(hit_record, src_gid, src_tid, tau=TAU):
    """Build a self-EXCLUDED copy of one amplicon's hit record for re-selection.

    Drops every top-hit that is the source organism or a near-identical neighbour
    (``genome_id == src_gid`` OR ``taxon_id == src_tid`` OR ``identity >= tau``),
    then re-ranks the survivors 1..n and recomputes ``best_identity`` (0.0 if
    nothing survives -> a clean 'below' abstain in select_representatives, never
    None which would crash ``_tier``).  Returns (record, n_removed).

    The pipeline RUN is never altered: this is a pure post-hoc transform of the
    SAME hits, exactly the in-memory ``exclude_self -> select_representatives``
    path mandated by DESIGN Section 5.4 (no JSON-rewriter, no re-alignment).
    """
    top = hit_record.get("top20", []) or []
    kept = []
    for h in top:
        hid = str(h.get("genome_id")) if h.get("genome_id") is not None else None
        htid = h.get("taxon_id")
        ident = h.get("identity")
        if hid is not None and hid == src_gid:
            continue
        if src_tid is not None and htid is not None and int(htid) == src_tid:
            continue
        if ident is not None and ident >= tau:
            continue
        kept.append(h)
    n_removed = len(top) - len(kept)

    # Re-rank 1..n on the survivors (rank is the only positional field
    # select_representatives reads off each hit besides identity/score).
    reranked = []
    for i, h in enumerate(sorted(kept, key=lambda x: (-(x.get("identity") or 0.0),
                                                      -(x.get("align_score") or 0.0))), 1):
        hh = dict(h)
        hh["rank"] = i
        reranked.append(hh)

    best_id = reranked[0]["identity"] if reranked else 0.0

    rec = dict(hit_record)
    rec["top20"] = reranked
    rec["best_identity"] = best_id
    # asv_len is required by select_representatives; carry it through unchanged.
    return rec, n_removed


# --------------------------------------------------------------------------- #
# Per-cell scoring
# --------------------------------------------------------------------------- #
def selected_rep_gids(sel_record):
    """Ordered, de-duplicated list of selected rep genome_ids from a selection record."""
    seen, out = set(), []
    for s in sel_record.get("selected", []) or []:
        gid = s.get("genome_id")
        if gid is None:
            continue
        gid = str(gid)
        if gid not in seen:
            seen.add(gid)
            out.append(gid)
    return out


def _isnan(x):
    try:
        return math.isnan(x)
    except (TypeError, ValueError):
        return False


def score_one(key, self_mode, sel_record, src_gid, src_tid, region, domain,
              provider, cache_keys, *, available):
    """Score one (genome, region, self_mode) cell against the source's PGFam set.

    ``sel_record`` is the selection record to use for THIS mode:
      * include       -> the as-is PATHS.selection_out record.
      * exclude_0.987 -> the re-selection over the self-excluded hit record.
    ``available`` is False only for exclude mode when no hit record existed to
    re-select from -> outcome ``no_hits`` (poisoned, excluded from denominators).
    """
    nan = float("nan")
    tier = (sel_record or {}).get("tier")
    n_selected = int((sel_record or {}).get("n_selected", 0) or 0)
    rep_gids = selected_rep_gids(sel_record or {})

    # --- Truth set T -----------------------------------------------------
    # A non-numeric (MAG) source has no BV-BRC genome to fetch -> truth missing.
    src_in_cache = (src_tid is not None) and (src_gid in cache_keys)
    T = provider(src_gid) if src_in_cache else set()
    truth_present = src_in_cache  # cache membership == fetched truth available

    # --- Predicted set P (union over selected reps) ----------------------
    reps_missing = [g for g in rep_gids if g not in cache_keys]
    P = set()
    for g in rep_gids:
        if g in cache_keys:
            P |= provider(g)

    source_self_only = (len(rep_gids) == 1 and rep_gids[0] == src_gid)
    n_reps_geneless = sum(1 for g in rep_gids if g in cache_keys and not provider(g))

    # --- Classify outcome (precedence order) -----------------------------
    if not available:
        # exclude mode with no hit record to re-select from: cannot honestly
        # build P -> poison this mode's cell (NOT counted as an abstain).
        outcome = OUTCOME_NO_HITS
    elif not truth_present:
        # Source PGFam set never fetched (MAG/non-numeric source, or S8 never
        # prefetched it): truth undefined -> recall uncomputable -> poison.
        outcome = OUTCOME_PGFAM_MISSING
    elif len(T) == 0:
        outcome = OUTCOME_TRUTH_EMPTY
    elif reps_missing:
        outcome = OUTCOME_PGFAM_MISSING
    elif n_selected == 0:
        outcome = OUTCOME_ABSTAINED
    elif len(P) == 0:
        outcome = OUTCOME_PRED_EMPTY
    else:
        outcome = OUTCOME_OK

    t, p = len(T), len(P)
    i = len(T & P) if (t and p) else 0

    # --- Metric values per outcome (NaN where undefined) -----------------
    if outcome == OUTCOME_OK:
        recall = i / t
        precision = i / p
        f1 = (2 * i / (t + p)) if (t + p) else nan
        jaccard = (i / (t + p - i)) if (t + p - i) else nan
    elif outcome in (OUTCOME_ABSTAINED, OUTCOME_PRED_EMPTY):
        # nothing predicted (or all reps gene-less) -> 0 capture, precision undefined
        recall, precision, f1, jaccard = 0.0, nan, 0.0, 0.0
    else:
        # truth_empty / pgfam_missing / no_hits: poisoned, all NaN (excluded from means)
        recall = precision = f1 = jaccard = nan

    return dict(
        ampliconKey=key,
        src_genome_id=src_gid,
        src_taxon_id=src_tid if src_tid is not None else "",
        region=region,
        domain=domain or "",
        tier=tier or "",
        self_mode=self_mode,
        n_selected=n_selected,
        n_reps=len(rep_gids),
        n_reps_in_cache=sum(1 for g in rep_gids if g in cache_keys),
        n_reps_missing=len(reps_missing),
        n_reps_geneless=n_reps_geneless,
        source_self_only=source_self_only,
        best_identity=(sel_record or {}).get("best_identity"),
        truth_size=t,
        pred_size=p,
        intersection=i,
        recall=recall,
        recall_pct=(100.0 * recall) if not _isnan(recall) else nan,
        precision=precision,
        f1=f1,
        jaccard=jaccard,
        outcome=outcome,
    )


def score_amplicon(key, inc_record, hit_record, provider, cache_keys, tau=TAU):
    """Score both self_modes for one amplicon; returns a list of 2 cell dicts.

    ``inc_record`` -- the as-is include-self selection record (PATHS.selection_out).
    ``hit_record`` -- the per-amplicon top-hit record (PATHS.hits_json) or None.
    """
    src_gid, region = split_key(key)
    src_tid = _src_taxon_id(src_gid)
    domain = domain_for(src_tid)

    # --- include mode: score the as-is selection ------------------------
    inc_cell = score_one(
        key, MODE_INCLUDE, inc_record, src_gid, src_tid, region, domain,
        provider, cache_keys, available=True)
    inc_cell["n_removed_self"] = 0

    # --- exclude_0.987 mode: re-select on the self-excluded hit record --
    if hit_record is not None and "top20" in hit_record:
        excl_rec, n_removed = exclude_self_record(hit_record, src_gid, src_tid, tau)
        excl_sel = SR.select_representatives(
            excl_rec, gene_provider=provider)
        excl_cell = score_one(
            key, MODE_EXCLUDE, excl_sel, src_gid, src_tid, region, domain,
            provider, cache_keys, available=True)
        excl_cell["n_removed_self"] = n_removed
    else:
        # No hits to re-select from -> exclude mode is poisoned (no_hits).
        excl_cell = score_one(
            key, MODE_EXCLUDE, None, src_gid, src_tid, region, domain,
            provider, cache_keys, available=False)
        excl_cell["n_removed_self"] = 0

    return [inc_cell, excl_cell]


# --------------------------------------------------------------------------- #
# Aggregation
# --------------------------------------------------------------------------- #
def _macro(values):
    """Macro (genome-equal) mean + normal-approx CI of finite values in [0,1].

    The values are bounded proportions, so a normal-approx mean CI (NOT a Wilson
    binomial CI) is the right interval for a mean of proportions.  Returns
    (mean, lo, hi, n).
    """
    vals = [v for v in values if v is not None and not _isnan(v)]
    n = len(vals)
    if n == 0:
        return (None, None, None, 0)
    mean = sum(vals) / n
    if n == 1:
        return (mean, mean, mean, 1)
    var = sum((v - mean) ** 2 for v in vals) / (n - 1)
    se = math.sqrt(var / n)
    half = 1.96 * se
    return (mean, max(0.0, mean - half), min(1.0, mean + half), n)


def _micro(cells):
    """Pooled micro metrics over OK cells, with Wilson CIs on the count ratios.

    Micro recall is a ratio of summed counts, so a Wilson interval on (Sum i,
    Sum t) is appropriate; precision pools on (Sum i, Sum p); F1 = 2*Si/(St+Sp);
    Jaccard = Si/(St+Sp-Si).
    """
    si = sum(c["intersection"] for c in cells)
    st = sum(c["truth_size"] for c in cells)
    sp = sum(c["pred_size"] for c in cells)
    out = {}
    if st > 0:
        lo, hi = wilson_ci(si, st)
        out["recall"] = dict(value=si / st, lo=lo, hi=hi, num=si, den=st)
    else:
        out["recall"] = dict(value=None, lo=None, hi=None, num=si, den=st)
    if sp > 0:
        lo, hi = wilson_ci(si, sp)
        out["precision"] = dict(value=si / sp, lo=lo, hi=hi, num=si, den=sp)
    else:
        out["precision"] = dict(value=None, lo=None, hi=None, num=si, den=sp)
    denom_f1 = st + sp
    out["f1"] = dict(value=(2 * si / denom_f1) if denom_f1 else None,
                     num=2 * si, den=denom_f1)
    denom_j = st + sp - si
    out["jaccard"] = dict(value=(si / denom_j) if denom_j else None,
                          num=si, den=denom_j)
    return out


def _macro_block(cells):
    """Macro means+CIs for all four metrics over the given cells (OK only)."""
    block = {}
    for m in METRICS:
        mean, lo, hi, n = _macro([c[m] for c in cells])
        block[m] = dict(value=mean, lo=lo, hi=hi, n=n)
    return block


def _outcome_counts(cells):
    return {o: sum(1 for c in cells if c["outcome"] == o) for o in OUTCOMES}


def _mode_block(cells, best_nonself):
    """Summary block for one self_mode's cells of one region (macro+micro+strata)."""
    ok = [c for c in cells if c["outcome"] == OUTCOME_OK]
    block = {
        "n_cells": len(cells),
        "outcome_counts": _outcome_counts(cells),
        "n_scorable_ok": len(ok),
        "macro": _macro_block(ok),
        "micro": _micro(ok),
        "by_domain": {},
        "by_tier": {},
        "by_identity_bin": {},
    }
    for dom in sorted({c["domain"] for c in cells if c["domain"]}):
        dok = [c for c in ok if c["domain"] == dom]
        drc = [c for c in cells if c["domain"] == dom]
        block["by_domain"][dom] = dict(
            n_cells=len(drc), n_scorable_ok=len(dok),
            outcome_counts=_outcome_counts(drc),
            macro=_macro_block(dok), micro=_micro(dok))
    for tier in TIER_ORDER:
        tok = [c for c in ok if c["tier"] == tier]
        trc = [c for c in cells if c["tier"] == tier]
        if not trc:
            continue
        block["by_tier"][tier] = dict(
            n_cells=len(trc), n_scorable_ok=len(tok),
            outcome_counts=_outcome_counts(trc),
            macro=_macro_block(tok), micro=_micro(tok))
    for lab in IDENT_BIN_LABELS + [IDENT_BIN_UNKNOWN]:
        bok = [c for c in ok if ident_bin(best_nonself.get(c["ampliconKey"])) == lab]
        brc = [c for c in cells if ident_bin(best_nonself.get(c["ampliconKey"])) == lab]
        if not brc:
            continue
        block["by_identity_bin"][lab] = dict(
            n_cells=len(brc), n_scorable_ok=len(bok),
            outcome_counts=_outcome_counts(brc),
            macro=_macro_block(bok), micro=_micro(bok))
    return block


def summarize(cells, best_nonself):
    """Per-region x self_mode summary (macro + micro + strata + denoms), + leakage gap."""
    regions = sorted({c["region"] for c in cells})
    summary = {"_meta": {
        "self_modes": SELF_MODES,
        "tau": TAU,
        "n_cells_total": len(cells),
        "outcome_counts_total": {
            m: _outcome_counts([c for c in cells if c["self_mode"] == m])
            for m in SELF_MODES},
        "metric_defs": {
            "recall": "i/t (% gene capture = 100*recall)",
            "precision": "i/p", "f1": "2i/(t+p)", "jaccard": "i/(t+p-i)",
        },
        "ident_bins": IDENT_BIN_LABELS,
        "tiers": TIER_ORDER,
    }, "regions": {}}

    for region in regions:
        rcells = [c for c in cells if c["region"] == region]
        rblock = {}
        for mode in SELF_MODES:
            mcells = [c for c in rcells if c["self_mode"] == mode]
            rblock[mode] = _mode_block(mcells, best_nonself)
        # leakage gap: include macro recall minus exclude macro recall.
        inc_r = rblock[MODE_INCLUDE]["macro"]["recall"]["value"]
        exc_r = rblock[MODE_EXCLUDE]["macro"]["recall"]["value"]
        rblock["leakage_gap_macro_recall"] = (
            None if (inc_r is None or exc_r is None) else round(inc_r - exc_r, 6))
        summary["regions"][region] = rblock

    # headline: macro recall % per region for each self_mode (the literal
    # "% gene capture by region" answer; exclude_0.987 = the honest number).
    headline = {}
    for mode in SELF_MODES:
        headline[mode] = {}
        for region in regions:
            v = summary["regions"][region][mode]["macro"]["recall"]["value"]
            headline[mode][region] = (None if v is None else round(100.0 * v, 3))
    summary["_meta"]["headline_macro_recall_pct"] = headline
    return summary


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #
CELL_FIELDS = [
    "ampliconKey", "src_genome_id", "src_taxon_id", "region", "domain", "tier",
    "self_mode", "n_selected", "n_reps", "n_reps_in_cache", "n_reps_missing",
    "n_reps_geneless", "n_removed_self", "source_self_only", "best_identity",
    "ident_bin", "best_nonself_identity", "truth_size", "pred_size",
    "intersection", "recall", "recall_pct", "precision", "f1", "jaccard",
    "outcome",
]


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, float) and _isnan(v):
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, float):
        return f"{v:.6g}"
    return v


def write_cells_csv(path, cells, best_nonself):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CELL_FIELDS, extrasaction="ignore")
        w.writeheader()
        for c in cells:
            bn = best_nonself.get(c["ampliconKey"])
            row = dict(c)
            row["best_nonself_identity"] = bn
            row["ident_bin"] = ident_bin(bn)
            w.writerow({k: _fmt(row.get(k)) for k in CELL_FIELDS})


LONG_FIELDS = ["ampliconKey", "src_genome_id", "region", "domain", "tier",
               "ident_bin", "self_mode", "outcome", "truth_size", "pred_size",
               "intersection", "metric", "value"]


def write_long_csv(path, cells, best_nonself):
    """Tidy long form: one row per (cell, metric) -- the renderers read this."""
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LONG_FIELDS)
        w.writeheader()
        for c in cells:
            base = dict(
                ampliconKey=c["ampliconKey"], src_genome_id=c["src_genome_id"],
                region=c["region"], domain=c["domain"], tier=c["tier"],
                ident_bin=ident_bin(best_nonself.get(c["ampliconKey"])),
                self_mode=c["self_mode"], outcome=c["outcome"],
                truth_size=c["truth_size"], pred_size=c["pred_size"],
                intersection=c["intersection"])
            for m in METRICS:
                row = dict(base, metric=m, value=_fmt(c[m]))
                w.writerow(row)


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def run(force=False):
    outs = (PATHS.genecap_cells_csv, PATHS.gene_capture_long_csv,
            PATHS.genecap_summary_json)
    if not force and all(C.done(p) for p in outs):
        return ("[score_genecap] all outputs present; skipping (use --force). "
                f"-> {', '.join(outs)}")

    if not C.done(PATHS.selection_out):
        raise SystemExit(f"[score_genecap] missing selection: {PATHS.selection_out}")

    selection = load_selection(PATHS.selection_out)
    cache_keys = load_cache_keys(PATHS.pgfam_cache)
    # One namespace for truth + predicted reps (apples-to-apples).  We never let
    # it fetch live here: every genome we query is pre-checked against cache_keys.
    provider = SR.bvbrc_gene_provider(PATHS.pgfam_cache)
    hits = load_hits(PATHS.hits_json)

    if not cache_keys:
        print(f"[score_genecap] WARNING: empty/absent PGFam cache "
              f"({PATHS.pgfam_cache}); every cell will be pgfam_missing. "
              f"Run S8 prefetch first.", file=sys.stderr)
    if not hits:
        print(f"[score_genecap] WARNING: absent/empty hits ({PATHS.hits_json}); "
              f"the {MODE_EXCLUDE} mode will be 'no_hits' for every amplicon. "
              f"Run S5 alignment first.", file=sys.stderr)

    # best-non-self identity per amplicon (the recall-vs-distance axis + bins).
    best_nonself = {}
    cells = []
    for key, inc_record in sorted(selection.items()):
        src_gid, _ = split_key(key)
        src_tid = _src_taxon_id(src_gid)
        hr = hits.get(key)
        best_nonself[key] = best_nonself_identity(hr, src_gid, src_tid)
        cells.extend(score_amplicon(key, inc_record, hr, provider, cache_keys))

    write_cells_csv(PATHS.genecap_cells_csv, cells, best_nonself)
    write_long_csv(PATHS.gene_capture_long_csv, cells, best_nonself)
    summary = summarize(cells, best_nonself)
    with open(PATHS.genecap_summary_json, "w") as fh:
        json.dump(summary, fh, indent=1)

    n_amplicons = len(selection)
    oc = summary["_meta"]["outcome_counts_total"]
    head = summary["_meta"]["headline_macro_recall_pct"]
    return (
        f"[score_genecap] {n_amplicons} amplicons -> {len(cells)} cells "
        f"({len(SELF_MODES)} self_modes) across {len(summary['regions'])} regions.\n"
        f"[score_genecap] include outcomes: {json.dumps(oc[MODE_INCLUDE])}\n"
        f"[score_genecap] exclude_0.987 outcomes: {json.dumps(oc[MODE_EXCLUDE])}\n"
        f"[score_genecap] headline macro recall %% (gene capture) per region:\n"
        f"  include       = {json.dumps(head[MODE_INCLUDE])}\n"
        f"  exclude_0.987 = {json.dumps(head[MODE_EXCLUDE])}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="recompute even if outputs exist")
    args = ap.parse_args()
    print(run(force=args.force))


if __name__ == "__main__":
    main()
