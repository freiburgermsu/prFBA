#!/usr/bin/env python
"""score_taxacc.py — S9 taxonomic-accuracy scoring (DESIGN.md §6).

Reports BOTH self scenarios for every amplicon (USER DECISION #4: the pipeline
RUN is never altered; exclude-self is post-hoc analysis on the SAME hits):

  * ``include``            — selection/scoring on the real, unfiltered pipeline
                             selection (``PATHS.selection_out``); the source
                             organism is in the reference DB — the realistic,
                             primary number.
  * ``exclude_0.987``      — generalization: from ``PATHS.hits_json`` drop every
  * ``exclude_0.995``        hit whose ``genome_id == src`` OR ``taxon_id ==
                             src_taxon`` OR ``identity >= tau`` (tau in
                             {0.987 species threshold, 0.995}), set
                             ``best_identity = kept[0].identity if kept else 0.0``
                             (0.0 -> ``below`` tier -> clean abstain; never None),
                             then re-run the PURE ``select_references.
                             select_representatives`` on the filtered record and
                             score it — performance as if the source organism
                             were NOT in the reference set.

Every per-record row AND every region aggregate carries a ``self_mode`` column
(``include`` | ``exclude_0.987`` | ``exclude_0.995``), plus ``domain`` and
``tier`` (DESIGN §6).

For every amplicon (keyed ``{genome_id}__{region}`` from the in-silico PCR) and
every self_mode, two predictors are scored against the source organism's true
lineage:

  * **anchor**    = ``selected[0]`` (``role=="anchor"``); ``n_selected==0`` -> ABSTAIN
  * **consensus** = per-rank majority over ``selected[]``; a tie -> AMBIG (no-call)

BOTH predictors are RE-RESOLVED through ``_common.lineage_for(taxon_id)`` (taxopy
+ prFBA taxdump), the *same* function that resolves the truth lineage, so the two
sides compare under one rank vocabulary (the "one comparability rule").  Names are
normalized with ``_common._norm`` before ``==``.

Per amplicon, per self_mode, per SCORED rank (Phylum..Species) the three-state
cell is:

    rank_correct(truth_v, pred_v):
        if pred_v is None or AMBIG:   return None   # no-call (abstention)
        if not truth_v:               return None   # truth undefined here
        return _norm(pred_v) == _norm(truth_v)      # True / False (misassign)

``None`` != ``False``: ``None`` is a no-call (excluded from %-correct numerator/
denominator, counted in coverage); ``False`` is an active misassignment.  The
**deepest correct rank** is the deepest SCORED rank that is ``True`` AND has a
contiguous True prefix (a deep match over a wrong Family does not count).

ANTI-INFLATION RULE (DESIGN §6 / verifier FINDING 7): every self_mode scores the
SAME amplicon population, so an exclude-self abstention counts as a no-call
against the full include-self denominator (coverage's ``n`` never shrinks toward
easy cases) — the include->exclude drop therefore reflects reality, not
denominator drift.

Aggregates are computed per ``region x self_mode x domain x tier`` (plus a
domain="all" / tier="all" rollup) for each predictor:
  * (correct, coverage) pair per rank, each with a Wilson 95% CI,
  * self-recovery (self genome in top20 / self genome selected) — include only,
  * abstention rate split by reason,
  * misassignment kinds + first-wrong-rank histogram + mean off_by_ranks,
  * deepest-correct-rank histogram,
  * tie-cluster size (ambiguity) summary.

Resumable: if every output exists and is non-empty, the script short-circuits
unless ``--force`` is given.

Sole interpreter: ``/home/freiburger/Documents/py_venv/bin/python``.

I/O CONTRACT
  reads : PATHS.hits_json        (asv_top20_alignment_hits.json — per amplicon top20)
          PATHS.selection_out    (select_references audit — per amplicon tier/selected[])
          PATHS.truth_json        (benchmark_truth.json — per amplicon src/domain/truth_lineage)
          + select_references.py  (in-memory exclude-self re-selection)
  writes: PATHS.taxacc_jsonl                  (one JSON record per scored amplicon, all self_modes nested)
          PATHS.concordance_long_csv          (region x self_mode x domain x tier x predictor x rank)
          PATHS.asv_region_concordance_csv    (one row per amplicon x self_mode x predictor)
          PATHS.self_recovery_csv             (region x self_mode x domain x tier self-recovery)
          PATHS.confusion_genus_csv           (region x self_mode genus true->pred confusion counts)
          PATHS.ambiguity_csv                 (region x self_mode x domain x tier tie-cluster sizes)
          PATHS.taxacc_summary_json           (nested region->predictor->self_mode->domain->tier metrics)
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict

# Import the shared contract module by absolute path so the script runs from any cwd.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# select_references lives in prFBA (the pipeline's selection module); add it so the
# exclude-self branch can call the *pure* select_representatives in-memory (DESIGN §5.4).
sys.path.insert(0, "/home/freiburger/Documents/prFBA")

import _common as C  # noqa: E402
from _common import PATHS, SCORED, _norm, done, lineage_for, wilson_ci  # noqa: E402
import select_references as SR  # noqa: E402

AMBIG = "__AMBIG__"  # sentinel: a per-rank consensus tie (treated as a no-call)

# Tie-cluster epsilon: references within this identity of the best are "tied".
TIE_EPS = 0.002

# Self-exclusion thresholds (DESIGN §5.4 / decision #4).  ``include`` keeps the
# real pipeline selection; each exclude_<tau> re-runs select_representatives on
# the source-organism + near-identical-neighbour-filtered hits.
INCLUDE_MODE = "include"
EXCLUDE_TAUS = (0.987, 0.995)
SELF_MODES = (INCLUDE_MODE,) + tuple(f"exclude_{tau:g}" for tau in EXCLUDE_TAUS)

PREDICTORS = ("anchor", "consensus")

# Selection regime for the exclude-self RE-SELECTION.  The committed validation
# (and the manuscript tables) were produced by the LEGACY ordered reducer, before
# the exact-tie union became select_references.DEFAULT_KNOBS.  To reproduce those
# numbers (so the only change is the added correct_population metric) run with
# --legacy, which pins the re-selection to select_all=False.  Default = whatever
# select_references currently defaults to (the exact-tie union).
LEGACY_SELECT = False


def _select_knobs():
    """Knobs override for the exclude-self re-selection (None = pipeline default)."""
    return {"select_all": False} if LEGACY_SELECT else None


# --------------------------------------------------------------------------- #
# Loading + record assembly
# --------------------------------------------------------------------------- #
def _strip_meta(d: dict) -> dict:
    """Drop the ``_meta`` audit header that select_references prepends."""
    return {k: v for k, v in d.items() if k != "_meta"}


def _src_from_key(key: str):
    """Recover (src_genome_id, region) from a ``{genome_id}__{region}`` amplicon key."""
    gid, _, region = key.rpartition("__")
    return gid, region


def _src_taxon(genome_id):
    """``int(genome_id.split('.')[0])`` for numeric BV-BRC ids, else None (MAG/unscoreable)."""
    try:
        return int(str(genome_id).split(".")[0])
    except (ValueError, AttributeError):
        return None


# --------------------------------------------------------------------------- #
# Exclude-self: filter hits + re-run the pure select_representatives
# --------------------------------------------------------------------------- #
def filter_self_hits(hit_record, src_genome_id, src_taxon_id, tau):
    """Drop source-organism / near-identical hits from one amplicon's top20.

    A hit is dropped iff (DESIGN decision #4):
        hit.genome_id == src_genome_id  OR
        hit.taxon_id  == src_taxon_id   OR
        hit.identity  >= tau            (the near-identical species/near floor)

    Surviving hits are re-ranked 1..n (``rank`` rewritten) so the downstream
    selector sees a clean top-N.  Returns (kept_hits, removed_hits).
    """
    kept, removed = [], []
    for h in hit_record.get("top20") or []:
        gid = str(h.get("genome_id"))
        tid = h.get("taxon_id")
        ident = h.get("identity") or 0.0
        if (
            gid == str(src_genome_id)
            or (src_taxon_id is not None and tid == src_taxon_id)
            or ident >= tau
        ):
            removed.append(h)
        else:
            kept.append(h)
    # re-rank survivors 1..n so select_representatives' anchor logic is clean
    out = []
    for i, h in enumerate(kept, start=1):
        hh = dict(h)
        hh["rank"] = i
        out.append(hh)
    return out, removed


def exclude_self_selection(hit_record, src_genome_id, src_taxon_id, tau):
    """Build the exclude-self selection record by re-running the PURE selector.

    Drops self + near-identical hits (``filter_self_hits``), sets
    ``best_identity = kept[0].identity if kept else 0.0`` (0.0 -> ``below`` tier
    -> clean abstain; NEVER None, which would crash ``_tier``), then calls the
    unmodified ``select_references.select_representatives`` (no gene_provider:
    the estimator path — re-selection is a pure function of the record).

    Returns (selection_record, n_removed).
    """
    kept, removed = filter_self_hits(hit_record, src_genome_id, src_taxon_id, tau)
    asv_record = {
        "asv_len": hit_record.get("asv_len"),
        "best_identity": (kept[0]["identity"] if kept else 0.0),
        "midas_taxonomy": hit_record.get("midas_taxonomy"),
        "rel_ab": hit_record.get("rel_ab"),
        "top20": kept,
    }
    sel = SR.select_representatives(asv_record, knobs=_select_knobs())
    return sel, len(removed)


# --------------------------------------------------------------------------- #
# Predictors
# --------------------------------------------------------------------------- #
def anchor_lineage(sel_record):
    """Anchor predictor: re-resolve ``selected[0]``'s taxon_id via lineage_for.

    Returns (lineage_dict_or_None, anchor_genome_id, anchor_taxon_id).
    ``None`` lineage encodes ABSTAIN (no selected genome).
    """
    selected = sel_record.get("selected") or []
    if not selected:
        return None, None, None
    a = selected[0]  # role == "anchor" by construction (stage-7 first inclusion)
    tid = a.get("taxon_id")
    return lineage_for(tid), a.get("genome_id"), tid


def consensus_lineage(sel_record):
    """Consensus predictor: per-rank majority over all ``selected[]`` re-resolved lineages.

    Each selected genome's taxon_id is re-resolved through ``lineage_for`` (same
    function as truth), then per SCORED rank we take the normalized-name majority.
    A tie (no strict winner) -> ``AMBIG`` (a no-call).  Empty selection -> all-None.

    Returns (lineage_dict {rank: name|AMBIG|None}, strength {rank: top/n}).
    """
    selected = sel_record.get("selected") or []
    lineages = [lineage_for(g.get("taxon_id")) for g in selected]
    cons, strength = {}, {}
    for rank in SCORED:
        # Count by normalized name but remember a representative display name.
        norm_counts = Counter()
        display = {}
        n_present = 0
        for lin in lineages:
            v = lin.get(rank)
            if not v:
                continue  # rank undefined for this rep — not a vote
            n_present += 1
            nv = _norm(v)
            norm_counts[nv] += 1
            display.setdefault(nv, v)
        if not norm_counts:
            cons[rank] = None
            strength[rank] = None
            continue
        ranked = norm_counts.most_common()
        top_norm, top_n = ranked[0]
        # strict majority requires a unique top count
        if len(ranked) > 1 and ranked[1][1] == top_n:
            cons[rank] = AMBIG
        else:
            cons[rank] = display[top_norm]
        strength[rank] = round(top_n / n_present, 4) if n_present else None
    return cons, strength


# --------------------------------------------------------------------------- #
# Three-state per-rank cell + deepest contiguous-correct rank
# --------------------------------------------------------------------------- #
def rank_correct(truth_v, pred_v):
    """Three-state cell (DESIGN §6).

    None  -> no-call (pred absent / AMBIG, OR truth undefined at this rank)
    False -> active misassignment
    True  -> correct
    """
    if pred_v is None or pred_v == AMBIG:
        return None
    if not truth_v:
        return None
    return _norm(pred_v) == _norm(truth_v)


def score_predictor(truth_lin, pred_lin):
    """Score one predictor's lineage vs truth across SCORED ranks.

    Returns dict with:
      cells          {rank: True/False/None}
      deepest_correct_rank   deepest SCORED rank with a contiguous True prefix (or None)
      first_wrong_rank       shallowest SCORED rank that is False (or None)
      shared_to              deepest SCORED rank that is True (non-contiguous ok) (or None)
      off_by_ranks           depth of truth - depth of deepest_correct (when misassigned)
      overconfident          pred named a rank deeper than truth defines (a no-call, flagged)
    """
    cells = {}
    if pred_lin is None:  # ABSTAIN: predictor made no call at any rank
        for rank in SCORED:
            cells[rank] = None
    else:
        for rank in SCORED:
            cells[rank] = rank_correct(truth_lin.get(rank), pred_lin.get(rank))

    # deepest contiguous-correct rank: walk shallow->deep, stop at first non-True
    deepest_correct = None
    for rank in SCORED:
        if cells[rank] is True:
            deepest_correct = rank
        else:
            break

    first_wrong = next((r for r in SCORED if cells[r] is False), None)
    shared_to = None
    for rank in SCORED:
        if cells[rank] is True:
            shared_to = rank

    # how deep does truth itself go (deepest non-empty truth rank)
    truth_depth = None
    for rank in SCORED:
        if truth_lin.get(rank):
            truth_depth = rank
    off_by = None
    if first_wrong is not None and truth_depth is not None:
        dc_idx = SCORED.index(deepest_correct) if deepest_correct else -1
        off_by = SCORED.index(truth_depth) - dc_idx

    # Overconfident is a SEPARATE signal from a misassignment: the three-state
    # rule makes "pred named a rank where truth is undefined" a no-call (None,
    # never True — DESIGN §6), not a False.  We flag it so the predictor never
    # gets silent credit for over-deep names: any rank deeper than truth's
    # deepest defined rank where the predictor still named something.
    overconfident = False
    if pred_lin is not None and truth_depth is not None:
        td_idx = SCORED.index(truth_depth)
        for rank in SCORED[td_idx + 1:]:
            pv = pred_lin.get(rank)
            if pv and pv != AMBIG:
                overconfident = True
                break

    return dict(
        cells=cells,
        deepest_correct_rank=deepest_correct,
        first_wrong_rank=first_wrong,
        shared_to=shared_to,
        off_by_ranks=off_by,
        overconfident=overconfident,
    )


def misassign_kind(score, self_in_top20):
    """Classify an active misassignment (a False cell at ``first_wrong_rank``).

    self_displaced self genome was in top20 but anchor went to a wrong neighbour
                   (the selector-quality KPI — takes precedence)
    sibling        wrong but shares the parent rank (diverges one level down)
    cross_lineage  diverges at the shallowest scored rank (no shared parent)

    Returns None when there is no active misassignment.  "overconfident" is NOT
    returned here — it is a distinct no-call signal carried on ``score`` and is
    counted separately (a predicted rank deeper than truth defines never scores
    True and is never an active False, so it cannot be a ``first_wrong_rank``).
    """
    fw = score["first_wrong_rank"]
    if fw is None:
        return None
    if self_in_top20:
        # the source was available as a hit yet the predictor picked a wrong neighbour
        return "self_displaced"
    fw_idx = SCORED.index(fw)
    # sibling: agrees through the rank immediately above first-wrong
    if fw_idx == 0:
        return "cross_lineage"
    parent = SCORED[fw_idx - 1]
    if score["cells"].get(parent) is True:
        return "sibling"
    return "cross_lineage"


# --------------------------------------------------------------------------- #
# Tie-cluster size (ambiguity) + self-recovery (top20)
# --------------------------------------------------------------------------- #
def tie_cluster_size(hit_record, eps=TIE_EPS):
    """#refs in top20 within ``eps`` identity of the best hit (ambiguity proxy)."""
    top = hit_record.get("top20") or []
    if not top:
        return 0, None
    best = max((h.get("identity") or 0.0) for h in top)
    n = sum(1 for h in top if (h.get("identity") or 0.0) >= best - eps)
    return n, round(best, 4)


def self_in_top20(hit_record, src_genome_id, src_taxon_id):
    """True iff the source organism appears in the top20 (by genome_id or taxon_id)."""
    for h in hit_record.get("top20") or []:
        if str(h.get("genome_id")) == str(src_genome_id):
            return True
        if src_taxon_id is not None and h.get("taxon_id") == src_taxon_id:
            return True
    return False


def self_selected(sel_record, src_genome_id, src_taxon_id):
    """True iff the source organism is among the selected reps."""
    for g in sel_record.get("selected") or []:
        if str(g.get("genome_id")) == str(src_genome_id):
            return True
        if src_taxon_id is not None and g.get("taxon_id") == src_taxon_id:
            return True
    return False


def self_is_anchor(sel_record, src_genome_id, src_taxon_id):
    """True iff the source organism is the anchor (selected[0])."""
    selected = sel_record.get("selected") or []
    if not selected:
        return False
    a = selected[0]
    if str(a.get("genome_id")) == str(src_genome_id):
        return True
    return src_taxon_id is not None and a.get("taxon_id") == src_taxon_id


def abstention_reason(sel_record):
    """Map a select_references audit record's flags to an abstention reason, or None."""
    if (sel_record.get("n_selected") or 0) > 0:
        return None
    flags = sel_record.get("flags") or []
    if any(f.startswith("abstain_below_family") for f in flags):
        return "below_family"
    if any(f.startswith("abstain_reliability") for f in flags):
        return "reliability"
    if sel_record.get("tier") == "below":
        return "below_family"
    return "other"


# --------------------------------------------------------------------------- #
# Aggregation accumulators
# --------------------------------------------------------------------------- #
class CovCorr:
    """Accumulate (correct, coverage) over a stream of three-state cells.

    ``n`` is the full amplicon count (coverage denominator); it is identical
    across self_modes because every self_mode scores the SAME amplicon
    population — the anti-inflation rule (DESIGN §6): an exclude-self abstention
    is a no-call against the full include-self denominator, never a shrunk one.
    """

    __slots__ = ("n", "called", "correct")

    def __init__(self):
        self.n = 0        # all amplicons (denominator for coverage)
        self.called = 0   # cells that are not None (denominator for correct)
        self.correct = 0  # cells that are True

    def add(self, cell):
        self.n += 1
        if cell is not None:
            self.called += 1
            if cell is True:
                self.correct += 1

    def metrics(self):
        corr = self.correct / self.called if self.called else None
        cov = self.called / self.n if self.n else None
        # correct_population: n_correct / n over the FULL amplicon population, i.e.
        # a no-call / abstention counts as NOT-correct (== correct * coverage).
        # This is the honest whole-population accuracy; `correct` alone is the
        # accuracy CONDITIONAL on the pipeline making a call, and under
        # self-exclusion the called-denominator shrinks toward the easy cases, so
        # `correct` overstates novel-organism generalization -- read it with
        # coverage or use correct_population.
        corr_pop = self.correct / self.n if self.n else None
        clo, chi = wilson_ci(self.correct, self.called)
        vlo, vhi = wilson_ci(self.called, self.n)
        plo, phi = wilson_ci(self.correct, self.n)
        return dict(
            n=self.n, n_called=self.called, n_correct=self.correct,
            correct=corr, correct_ci=[clo, chi],
            coverage=cov, coverage_ci=[vlo, vhi],
            correct_population=corr_pop, correct_population_ci=[plo, phi],
        )


def _nest():
    """region -> self_mode -> domain -> tier -> predictor -> rank -> CovCorr."""
    return defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(
        lambda: defaultdict(lambda: defaultdict(CovCorr))))))


def _nest_counter():
    """region -> self_mode -> domain -> tier -> predictor -> Counter."""
    return defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(
        lambda: defaultdict(Counter)))))


def _nest_list():
    """region -> self_mode -> domain -> tier -> predictor -> list."""
    return defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(
        lambda: defaultdict(list)))))


# --------------------------------------------------------------------------- #
# Main scoring pass
# --------------------------------------------------------------------------- #
def run(force: bool = False):
    outputs = [
        PATHS.taxacc_jsonl, PATHS.concordance_long_csv, PATHS.asv_region_concordance_csv,
        PATHS.self_recovery_csv, PATHS.confusion_genus_csv, PATHS.ambiguity_csv,
        PATHS.taxacc_summary_json,
    ]
    if not force and all(done(p) for p in outputs):
        print(f"[score_taxacc] all {len(outputs)} outputs present — skipping (use --force)")
        return

    for p in (PATHS.hits_json, PATHS.selection_out, PATHS.truth_json):
        if not done(p):
            raise SystemExit(f"[score_taxacc] required input missing/empty: {p}")

    print("[score_taxacc] loading hits / selection / truth ...")
    hits = json.load(open(PATHS.hits_json))
    selection = _strip_meta(json.load(open(PATHS.selection_out)))
    truth = json.load(open(PATHS.truth_json))

    # --- accumulators (all keyed by self_mode as the 2nd axis) ---
    cov = _nest()                          # [region][sm][dom][tier][pred][rank] = CovCorr
    # deepest-correct histograms: deep[region][sm][dom][tier][pred] = Counter(rank|"none")
    deep = _nest_counter()
    kinds = _nest_counter()                # misassign kind + overconfident counts
    firstwrong = _nest_counter()           # first-wrong-rank histogram
    offby = _nest_list()                   # off_by_ranks values
    # self-recovery counters: self_ct[region][sm][dom][tier] = Counter(keys)
    self_ct = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(Counter))))
    # abstention: abst[region][sm][dom][tier] = Counter(reason|"called")
    abst = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(Counter))))
    # ambiguity (tie-cluster sizes): amb[region][sm][dom][tier] = list[int]
    amb = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    # genus confusion: conf[region][sm][(true_genus, pred_genus)] -> count (anchor predictor)
    conf = defaultdict(lambda: defaultdict(Counter))

    jsonl_fh = open(PATHS.taxacc_jsonl, "w")
    asv_rows = []  # asv_region_concordance.csv (one row per amplicon x self_mode x predictor)

    n_total = n_scored = n_unscoreable = n_no_selection = 0
    # exclude-self bookkeeping (DESIGN §5.4 top-20 depth guard)
    n_removed_total = {sm: 0 for sm in SELF_MODES if sm != INCLUDE_MODE}
    n_empty_after = {sm: 0 for sm in SELF_MODES if sm != INCLUDE_MODE}

    for key, tinfo in truth.items():
        n_total += 1
        truth_lin = tinfo.get("truth_lineage") or {}
        src_gid = tinfo.get("src_genome_id")
        src_tid = tinfo.get("src_taxon_id")
        region = tinfo.get("region")
        domain = tinfo.get("domain") or "Unknown"
        if region is None or src_gid is None:
            g2, r2 = _src_from_key(key)
            region = region or r2
            src_gid = src_gid or g2
        if src_tid is None:
            src_tid = _src_taxon(src_gid)

        # Unscoreable: truth lineage is all-None (non-numeric / MAG source) -> log + skip.
        if not any(truth_lin.get(r) for r in SCORED):
            n_unscoreable += 1
            continue

        hit_record = hits.get(key, {})
        inc_record = selection.get(key)
        if inc_record is None:
            # No selection record for this amplicon (did not amplify / no hits).
            n_no_selection += 1
            inc_record = {"selected": [], "n_selected": 0, "tier": "below", "flags": []}

        n_scored += 1

        # --- self / top20 facts (computed on the unfiltered hits — include-self) ---
        s_top20 = self_in_top20(hit_record, src_gid, src_tid)
        s_sel = self_selected(inc_record, src_gid, src_tid)
        s_anchor = self_is_anchor(inc_record, src_gid, src_tid)
        tcs, best_id = tie_cluster_size(hit_record)

        # --- build per self_mode selection records ---
        # include = the real pipeline selection; exclude_<tau> = re-run the pure selector
        sel_by_mode = {INCLUDE_MODE: inc_record}
        removed_by_mode = {INCLUDE_MODE: 0}
        for tau in EXCLUDE_TAUS:
            sm = f"exclude_{tau:g}"
            sel_x, nrem = exclude_self_selection(hit_record, src_gid, src_tid, tau)
            sel_by_mode[sm] = sel_x
            removed_by_mode[sm] = nrem
            n_removed_total[sm] += nrem
            if (sel_x.get("n_selected") or 0) == 0:
                n_empty_after[sm] += 1

        # --- per-record JSONL: all self_modes nested ---
        rec = dict(
            amplicon=key, src_genome_id=src_gid, src_taxon_id=src_tid,
            region=region, domain=domain,
            self_in_top20=s_top20, self_selected=s_sel, self_is_anchor=s_anchor,
            tie_cluster_size=tcs, best_hit_identity=best_id,
            truth_lineage={r: truth_lin.get(r) for r in SCORED},
            self_modes={},
        )

        for sm in SELF_MODES:
            sel_record = sel_by_mode[sm]
            tier = sel_record.get("tier") or "below"
            reason = abstention_reason(sel_record)

            a_lin, a_gid, a_tid = anchor_lineage(sel_record)
            c_lin, c_strength = consensus_lineage(sel_record)
            pred_lins = {"anchor": a_lin, "consensus": c_lin}
            pred_scores = {
                "anchor": score_predictor(truth_lin, a_lin),
                "consensus": score_predictor(truth_lin, c_lin),
            }

            sm_block = dict(
                tier=tier,
                best_identity=sel_record.get("best_identity"),
                n_selected=sel_record.get("n_selected", len(sel_record.get("selected") or [])),
                n_removed_self=removed_by_mode[sm],
                abstention_reason=reason,
                anchor=dict(
                    genome_id=a_gid, taxon_id=a_tid,
                    lineage=(None if a_lin is None else {r: a_lin.get(r) for r in SCORED}),
                ),
                consensus=dict(
                    lineage={r: c_lin.get(r) for r in SCORED},
                    strength=c_strength,
                ),
                scores={},
            )

            for pred in PREDICTORS:
                sc = pred_scores[pred]
                kind = misassign_kind(sc, s_top20) if pred_lins[pred] is not None else None
                sm_block["scores"][pred] = dict(
                    cells={r: sc["cells"][r] for r in SCORED},
                    deepest_correct_rank=sc["deepest_correct_rank"],
                    first_wrong_rank=sc["first_wrong_rank"],
                    shared_to=sc["shared_to"],
                    off_by_ranks=sc["off_by_ranks"],
                    misassign_kind=kind,
                    overconfident=sc["overconfident"],
                    abstained=(pred_lins[pred] is None),
                )

                # --- aggregate into region x self_mode x domain x tier (and rollups) ---
                for dkey in (domain, "all"):
                    for tkey in (tier, "all"):
                        for rank in SCORED:
                            cov[region][sm][dkey][tkey][pred][rank].add(sc["cells"][rank])
                        deep[region][sm][dkey][tkey][pred][sc["deepest_correct_rank"] or "none"] += 1
                        if kind is not None:
                            kinds[region][sm][dkey][tkey][pred][kind] += 1
                        if sc["overconfident"]:
                            kinds[region][sm][dkey][tkey][pred]["overconfident"] += 1
                        if sc["first_wrong_rank"] is not None:
                            firstwrong[region][sm][dkey][tkey][pred][sc["first_wrong_rank"]] += 1
                        if sc["off_by_ranks"] is not None:
                            offby[region][sm][dkey][tkey][pred].append(sc["off_by_ranks"])

                # asv_region_concordance.csv row (carries self_mode + domain + tier)
                asv_rows.append(dict(
                    amplicon=key, region=region, self_mode=sm, domain=domain, tier=tier,
                    predictor=pred, src_genome_id=src_gid, src_taxon_id=src_tid,
                    best_identity=sel_record.get("best_identity"),
                    n_selected=sm_block["n_selected"], n_removed_self=removed_by_mode[sm],
                    self_in_top20=int(s_top20), self_selected=int(s_sel),
                    self_is_anchor=int(s_anchor), tie_cluster_size=tcs,
                    **{f"correct_{r}": (
                        "" if sc["cells"][r] is None else int(sc["cells"][r])) for r in SCORED},
                    deepest_correct_rank=sc["deepest_correct_rank"] or "",
                    first_wrong_rank=sc["first_wrong_rank"] or "",
                    misassign_kind=kind or "",
                    overconfident=int(sc["overconfident"]),
                    abstained=int(pred_lins[pred] is None),
                ))

            # --- self-recovery + abstention + ambiguity strata (predictor-independent) ---
            for dkey in (domain, "all"):
                for tkey in (tier, "all"):
                    sc_ct = self_ct[region][sm][dkey][tkey]
                    sc_ct["n"] += 1
                    # self_selected / self_anchor only meaningful for include (exclude
                    # drops the source by construction); recorded for completeness.
                    sm_self_sel = self_selected(sel_record, src_gid, src_tid)
                    sm_self_anchor = self_is_anchor(sel_record, src_gid, src_tid)
                    sc_ct["self_in_top20"] += int(s_top20)
                    sc_ct["self_selected"] += int(sm_self_sel)
                    sc_ct["self_anchor"] += int(sm_self_anchor)
                    sc_ct["selected_given_top20"] += int(sm_self_sel and s_top20)
                    sc_ct["anchor_given_top20"] += int(sm_self_anchor and s_top20)
                    abst[region][sm][dkey][tkey][reason or "called"] += 1
                    amb[region][sm][dkey][tkey].append(tcs)

            # --- genus confusion (anchor predictor; only where genus is misassigned) ---
            a_genus = None if a_lin is None else a_lin.get("Genus")
            t_genus = truth_lin.get("Genus")
            if t_genus and a_genus and a_genus != AMBIG:
                if _norm(a_genus) != _norm(t_genus):
                    conf[region][sm][(t_genus, a_genus)] += 1

            # attach this self_mode's full block to the per-record JSONL
            rec["self_modes"][sm] = sm_block

        jsonl_fh.write(json.dumps(rec) + "\n")

    jsonl_fh.close()

    print(f"[score_taxacc] amplicons: {n_total} total, {n_scored} scored, "
          f"{n_unscoreable} unscoreable (all-None truth), "
          f"{n_no_selection} with no selection record")
    for sm in n_removed_total:
        empty = n_empty_after[sm]
        rate = (empty / n_scored) if n_scored else 0.0
        print(f"[score_taxacc] {sm}: removed {n_removed_total[sm]} self/near hits, "
              f"{empty}/{n_scored} amplicons empty after filter ({rate*100:.1f}%)")

    # ----------------------------------------------------------------------- #
    # concordance_long.csv  (region x self_mode x domain x tier x predictor x rank)
    # ----------------------------------------------------------------------- #
    with open(PATHS.concordance_long_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "region", "self_mode", "domain", "tier", "predictor", "rank",
            "n", "n_called", "n_correct",
            "correct", "correct_ci_lo", "correct_ci_hi",
            "coverage", "coverage_ci_lo", "coverage_ci_hi",
            "correct_population", "correct_population_ci_lo", "correct_population_ci_hi",
        ])
        for region in sorted(cov):
            for sm in SELF_MODES:
                if sm not in cov[region]:
                    continue
                for dkey in sorted(cov[region][sm]):
                    for tkey in sorted(cov[region][sm][dkey]):
                        for pred in PREDICTORS:
                            for rank in SCORED:
                                m = cov[region][sm][dkey][tkey][pred][rank].metrics()
                                w.writerow([
                                    region, sm, dkey, tkey, pred, rank,
                                    m["n"], m["n_called"], m["n_correct"],
                                    _f(m["correct"]), _f(m["correct_ci"][0]), _f(m["correct_ci"][1]),
                                    _f(m["coverage"]), _f(m["coverage_ci"][0]), _f(m["coverage_ci"][1]),
                                    _f(m["correct_population"]),
                                    _f(m["correct_population_ci"][0]), _f(m["correct_population_ci"][1]),
                                ])

    # ----------------------------------------------------------------------- #
    # asv_region_concordance.csv  (one row per amplicon x self_mode x predictor)
    # ----------------------------------------------------------------------- #
    if asv_rows:
        fields = list(asv_rows[0].keys())
        with open(PATHS.asv_region_concordance_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(asv_rows)
    else:
        open(PATHS.asv_region_concordance_csv, "w").write("")

    # ----------------------------------------------------------------------- #
    # self_recovery.csv  (region x self_mode x domain x tier)
    # ----------------------------------------------------------------------- #
    with open(PATHS.self_recovery_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "region", "self_mode", "domain", "tier", "n",
            "n_self_in_top20", "n_self_selected", "n_self_anchor",
            "self_top20_rate", "self_top20_ci_lo", "self_top20_ci_hi",
            "self_anchor_rate", "self_anchor_ci_lo", "self_anchor_ci_hi",
            "n_selected_given_top20", "selected_given_top20_rate",
            "selected_given_top20_ci_lo", "selected_given_top20_ci_hi",
            "n_anchor_given_top20", "anchor_given_top20_rate",
        ])
        for region in sorted(self_ct):
            for sm in SELF_MODES:
                if sm not in self_ct[region]:
                    continue
                for dkey in sorted(self_ct[region][sm]):
                    for tkey in sorted(self_ct[region][sm][dkey]):
                        c = self_ct[region][sm][dkey][tkey]
                        n = c["n"]
                        top20 = c["self_in_top20"]
                        t20_lo, t20_hi = wilson_ci(top20, n)
                        a_lo, a_hi = wilson_ci(c["self_anchor"], n)
                        sg_lo, sg_hi = wilson_ci(c["selected_given_top20"], top20)
                        w.writerow([
                            region, sm, dkey, tkey, n,
                            top20, c["self_selected"], c["self_anchor"],
                            _f(_rate(top20, n)), _f(t20_lo), _f(t20_hi),
                            _f(_rate(c["self_anchor"], n)), _f(a_lo), _f(a_hi),
                            c["selected_given_top20"], _f(_rate(c["selected_given_top20"], top20)),
                            _f(sg_lo), _f(sg_hi),
                            c["anchor_given_top20"], _f(_rate(c["anchor_given_top20"], top20)),
                        ])

    # ----------------------------------------------------------------------- #
    # confusion_genus.csv  (region x self_mode true->pred genus counts)
    # ----------------------------------------------------------------------- #
    with open(PATHS.confusion_genus_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["region", "self_mode", "true_genus", "pred_genus", "count"])
        for region in sorted(conf):
            for sm in SELF_MODES:
                if sm not in conf[region]:
                    continue
                for (tg, pg), n in sorted(conf[region][sm].items(), key=lambda kv: -kv[1]):
                    w.writerow([region, sm, tg, pg, n])

    # ----------------------------------------------------------------------- #
    # ambiguity.csv  (region x self_mode x domain x tier tie-cluster size distribution)
    # ----------------------------------------------------------------------- #
    with open(PATHS.ambiguity_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "region", "self_mode", "domain", "tier", "n",
            "tie_mean", "tie_median", "tie_p90", "tie_max", "frac_tie_gt1",
        ])
        for region in sorted(amb):
            for sm in SELF_MODES:
                if sm not in amb[region]:
                    continue
                for dkey in sorted(amb[region][sm]):
                    for tkey in sorted(amb[region][sm][dkey]):
                        vals = amb[region][sm][dkey][tkey]
                        w.writerow([
                            region, sm, dkey, tkey, len(vals),
                            _f(_mean(vals)), _f(_median(vals)),
                            _f(_pct(vals, 90)), (max(vals) if vals else ""),
                            _f(_rate(sum(1 for v in vals if v > 1), len(vals))),
                        ])

    # ----------------------------------------------------------------------- #
    # taxacc_region_summary.json  (region->predictor->self_mode->domain->tier metrics)
    # ----------------------------------------------------------------------- #
    summary = {
        "_meta": dict(
            scenario="include_and_exclude_self",
            self_modes=list(SELF_MODES), exclude_taus=list(EXCLUDE_TAUS),
            headline_self_mode=f"exclude_{EXCLUDE_TAUS[0]:g}",
            n_total=n_total, n_scored=n_scored,
            n_unscoreable=n_unscoreable, n_no_selection=n_no_selection,
            n_removed_self=n_removed_total, n_empty_after_exclude=n_empty_after,
            scored_ranks=SCORED, predictors=list(PREDICTORS), tie_eps=TIE_EPS,
            select_mode=("legacy" if LEGACY_SELECT else "default_union"),
            metric_defs=dict(
                coverage="n_called / n (fraction of amplicons that got a call); "
                         "n is the full scored-amplicon population and is identical "
                         "across self_modes at domain=all/tier=all, so an exclude-self "
                         "abstention is a no-call, never a removed sample",
                correct="n_correct / n_called -- accuracy CONDITIONAL on a call. Under "
                        "self-exclusion the called-denominator shrinks toward easy cases, "
                        "so this OVERSTATES novel-organism generalization; read it with "
                        "coverage, or use correct_population",
                correct_population="n_correct / n -- accuracy over the FULL population "
                                   "(a no-call/abstention counts as not-correct); == "
                                   "correct * coverage; the honest generalization number",
            ),
            headline_metric="correct_population",
            inputs=dict(hits=PATHS.hits_json, selection=PATHS.selection_out,
                        truth=PATHS.truth_json),
        ),
        "regions": {},
    }
    for region in sorted(cov):
        rblock = summary["regions"].setdefault(region, {})
        for pred in PREDICTORS:
            pblock = rblock.setdefault(pred, {})
            for sm in SELF_MODES:
                if sm not in cov[region]:
                    continue
                smblock = pblock.setdefault(sm, {})
                for dkey in sorted(cov[region][sm]):
                    dblock = smblock.setdefault(dkey, {})
                    for tkey in sorted(cov[region][sm][dkey]):
                        ranks = {rank: cov[region][sm][dkey][tkey][pred][rank].metrics()
                                 for rank in SCORED}
                        obv = offby[region][sm][dkey][tkey][pred]
                        dblock[tkey] = dict(
                            per_rank=ranks,
                            deepest_correct_hist=dict(deep[region][sm][dkey][tkey][pred]),
                            misassign_kinds=dict(kinds[region][sm][dkey][tkey][pred]),
                            first_wrong_hist=dict(firstwrong[region][sm][dkey][tkey][pred]),
                            mean_off_by_ranks=_f(_mean(obv)),
                            n_misassigned=len(obv),
                        )
        # predictor-independent strata attached once per region, keyed by self_mode
        rblock["self_recovery"] = _nested_self(self_ct[region])
        rblock["abstention"] = _nested_abst(abst[region])
        rblock["ambiguity"] = _nested_amb(amb[region])

    json.dump(summary, open(PATHS.taxacc_summary_json, "w"), indent=1)

    print(f"[score_taxacc] wrote:\n  " + "\n  ".join(outputs))


# --------------------------------------------------------------------------- #
# small numeric / formatting helpers
# --------------------------------------------------------------------------- #
def _f(x):
    return "" if x is None else round(float(x), 6)


def _rate(k, n):
    return (k / n) if n else None


def _mean(xs):
    return (sum(xs) / len(xs)) if xs else None


def _median(xs):
    if not xs:
        return None
    s = sorted(xs)
    m = len(s) // 2
    return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2


def _pct(xs, p):
    if not xs:
        return None
    s = sorted(xs)
    idx = min(len(s) - 1, int(round((p / 100.0) * (len(s) - 1))))
    return s[idx]


def _nested_self(node):
    """self_mode -> domain -> tier -> self-recovery block."""
    out = {}
    for sm in SELF_MODES:
        if sm not in node:
            continue
        out[sm] = {}
        for dkey in sorted(node[sm]):
            out[sm][dkey] = {}
            for tkey in sorted(node[sm][dkey]):
                c = node[sm][dkey][tkey]
                n, top20 = c["n"], c["self_in_top20"]
                out[sm][dkey][tkey] = dict(
                    n=n, self_in_top20=top20, self_selected=c["self_selected"],
                    self_anchor=c["self_anchor"],
                    self_top20_rate=_f(_rate(top20, n)),
                    self_selected_rate=_f(_rate(c["self_selected"], n)),
                    self_anchor_rate=_f(_rate(c["self_anchor"], n)),
                    selected_given_top20=_f(_rate(c["selected_given_top20"], top20)),
                    anchor_given_top20=_f(_rate(c["anchor_given_top20"], top20)),
                )
    return out


def _nested_abst(node):
    out = {}
    for sm in SELF_MODES:
        if sm not in node:
            continue
        out[sm] = {}
        for dkey in sorted(node[sm]):
            out[sm][dkey] = {}
            for tkey in sorted(node[sm][dkey]):
                c = node[sm][dkey][tkey]
                n = sum(c.values())
                n_abst = sum(v for k, v in c.items() if k != "called")
                out[sm][dkey][tkey] = dict(
                    n=n, n_abstained=n_abst,
                    abstention_rate=_f(_rate(n_abst, n)),
                    by_reason={k: v for k, v in c.items() if k != "called"},
                )
    return out


def _nested_amb(node):
    out = {}
    for sm in SELF_MODES:
        if sm not in node:
            continue
        out[sm] = {}
        for dkey in sorted(node[sm]):
            out[sm][dkey] = {}
            for tkey in sorted(node[sm][dkey]):
                vals = node[sm][dkey][tkey]
                out[sm][dkey][tkey] = dict(
                    n=len(vals), tie_mean=_f(_mean(vals)), tie_median=_f(_median(vals)),
                    tie_p90=_f(_pct(vals, 90)), tie_max=(max(vals) if vals else None),
                    frac_tie_gt1=_f(_rate(sum(1 for v in vals if v > 1), len(vals))),
                )
    return out


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="recompute even if all outputs already exist")
    ap.add_argument("--legacy", action="store_true",
                    help="pin the exclude-self re-selection to the LEGACY ordered "
                         "reducer (select_all=False) to reproduce the committed "
                         "validation / manuscript tables")
    args = ap.parse_args()
    global LEGACY_SELECT
    LEGACY_SELECT = args.legacy
    C.ensure_dirs()
    run(force=args.force)


if __name__ == "__main__":
    main()
