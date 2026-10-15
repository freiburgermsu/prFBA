#!/usr/bin/env python
"""Score the shared-reference k-fold benchmark: every arm, one truth, one vocabulary.

WHY THIS COMPONENT EXISTS
-------------------------
The predecessor head-to-head (``picrust2_comparison/``) is confounded two ways: the
answer key was BV-BRC's annotation, which is also the universe prFBA's predictions are
drawn from (scoring the SAME predictions against each system's own annotation flips the
ranking, +0.050 vs -0.204 F1), and the two reference databases differ by 4x in size.
The shared-reference design removes both: one annotation (the BV-BRC EC projection)
is both methods' gene content AND the answer key, and the reference genome set is
identical between methods and rebuilt de novo per fold.  This script is the only place
the resulting predictions become numbers, so it is also the only place the remaining
confounds can be made visible rather than silently averaged away.  It therefore:

  * scores EVERY arm against the held-out genome's own EC set in ONE fixed vocabulary
    (``panel_ec_universe.json``), in EC space and in ModelSEED-reaction space (via
    ``function_maps.ec_rxn``, applied identically to truth and to every arm);
  * reports the population mean (an abstention is a zero) as primary and the
    conditional mean only on the subset BOTH methods answered, with its denominator
    printed, because conditional-only means compare two differently-selected
    populations (D7: letting PICRUSt2 abstain on its own worst 24.9% by NSTI shrinks
    prFBA's advantage from +0.070 to +0.017);
  * gives each method a risk-coverage curve swept on ITS OWN knob (prFBA: rank-1
    identity on the forced-answer arm; PICRUSt2: NSTI), never on the score being
    measured;
  * carries five baselines -- constant/most-prevalent, nearest-surviving-neighbour,
    family donor, clade-max oracle, and random-abstention -- so an absolute F1 is
    interpretable against a floor and a copy-a-relative ceiling (skill score);
  * clusters all inference by holdout group (measured ICC 0.66 genome / 0.48 genus /
    0.35 family: 67k amplicons carry ~3.8-13.4k independent observations), and never
    treats the k fold means as replicates;
  * bins BOTH methods by the neutral per-amplicon novelty covariates
    ``rank1_identity`` / ``rank1_rank`` -- properties of the shared reference, not of a
    method -- because NSTI is tree-specific (0.42x branch-length scale under our
    recipe) and ``best_identity`` is prFBA's own confidence;
  * excludes ``--min_align``-dropped amplicons and empty-truth cells from BOTH methods'
    denominators identically (C7, C10), so no shared data limitation is charged to one
    method;
  * writes every control verdict as a mechanical gate with an explicit falsifying
    number into ``controls.json``, merging with the verdicts other invocations
    (recipe control, include-self, plumbing) recorded.

GUARANTEES
----------
Stage 1 writes ``cells.csv.gz`` (one row per amplicon x arm x space, the expensive
pass).  Stage 2 derives the whole summary by streaming that file twice, so re-running
with an existing cell table (no ``--force``) reproduces the summary without rescoring,
and every number in ``summary_<rank>.json`` is auditable from the cell table.  Output
is deterministic: all iteration is over sorted keys and every random choice is seeded
from ``frozen.seed``.

CELL-TABLE CONVENTIONS (relied on by the promoted audit helpers and by phenotype_fold)
  * metric columns are BLANK, never 0, for a row that made no prediction
    (``abstained`` / ``pred_empty``) and for an unscoreable row; the population view is
    "blank f1 counts as 0.0", the conditional view is "blank f1 is excluded".
  * ``n_truth`` is always filled (it is 0 exactly for ``unscoreable_empty_truth``).
  * ``ident_bin`` bins the ARM's own best identity (blank for arms without one);
    ``novelty_bin`` bins the neutral ``rank1_identity``.  Both use
    ``score_genecap.ident_bin`` edges.
"""
import argparse
import collections
import csv
import glob
import gzip
import json
import math
import os
import random
import re
import statistics as st
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np
import orjson

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.abspath(os.path.join(HERE, ".."))             # region_validation/shared_reference
RV = os.path.abspath(os.path.join(SR, ".."))               # region_validation
REPO = os.path.abspath(os.path.join(RV, ".."))             # prFBA
P2C = os.path.join(RV, "picrust2_comparison")
for p in (os.path.join(RV, "scripts"), os.path.join(P2C, "scripts"), HERE, REPO):
    if p not in sys.path:
        sys.path.insert(0, p)

import _common as C                                        # noqa: E402
import compare as CMP                                      # noqa: E402
import score_genecap as SG                                 # noqa: E402
from annotation_ceiling import summarize                    # noqa: E402
from compare import Agg, load_picrust2, metrics            # noqa: E402

# SPEC Section 13: the baselines, the ICC/n_eff computation and the risk-coverage sweep
# are the promoted audit scripts, imported -- not re-derived -- so the harness and the
# audit report one body of code.
try:
    import design_audit_baselines as BASE
    import design_audit_icc as ICC
    import design_audit_riskcoverage as RC
except ImportError as exc:                                  # pragma: no cover
    raise SystemExit(
        "[score_folds] cannot import the promoted design-audit modules (%s). They must "
        "live beside this script as design_audit_{baselines,icc,riskcoverage}.py "
        "(SPEC Section 13)." % exc)

# --------------------------------------------------------------------------- #
# Section 1 -- frozen configuration, recorded verbatim into every output JSON
# --------------------------------------------------------------------------- #
FROZEN = {
    "tau_amplicon": 0.987,
    "tau_full16s": 0.987,
    "k_folds": 10,
    "seed": 1729,
    "prfba": {"k2": 200, "gpu_cand": 400, "prefilter_k": 1000, "min_ref_len": 200,
              "force_backend": "cuda",
              "select_modes": {
                  "legacy": {"select_all": False},
                  "exact_tie": {},
                  "graduated": {"select_all_window": 0.005,
                                "select_all_unconditional_window": 0.0},
                  "exact_tie_guards": {"select_all_coherence": True,
                                       "select_all_contam": True,
                                       "select_all_cap": 6}},
              "forced_answer_knobs": {"t_family": 0.0, "family_floor": 0.0},
              "t_species": 0.987, "t_genus": 0.945, "t_family": 0.865,
              "family_floor": 0.865, "decay_tau": 0.01},
    "picrust2": {"hsp_method": "mp", "edge_exponent": 0.5, "min_align": 0.80,
                 "placement_tool": "epa-ng", "max_nsti_reported": [None, 2.0],
                 "calculate_nsti": True},
    "reference_qc": {"len_min": 1200, "len_max": 1800, "occupancy_min": 0.90,
                     "ec_min": 1},
    "fasttree_args": ["-nt", "-gtr", "-gamma"],
    "density_fractions": [0.25, 0.5, 1.0],
    "sesoi_delta_f1": 0.02,
    "label_prefix": "g_",
}
SELECT_MODES = ("legacy", "exact_tie", "graduated", "exact_tie_guards")
PRIMARY_DENSITY = 1.0
COVERAGE_GRID = (0.5, 0.75, 0.9, 1.0)
BASELINE_ARMS = ("base_constant", "base_nn1", "base_family_donor", "base_clade_max",
                 "base_random_abstain")
RANK_ORDER = ("species", "genus", "family", "order", "class", "phylum", "domain")
# the oracle donor rank under each blocked rank (D13): the blocked rank has no survivor
BLOCKED_TO_DONOR = {"genus": "family", "family": "order", "order": "class",
                    "md5only": "genus"}

CELL_FIELDS = ["amplicon", "md5", "region", "domain", "src_genome_id", "fold",
               "group_id", "density", "arm", "space", "outcome", "n_truth", "n_pred",
               "tp", "recall", "precision", "f1", "jaccard", "mcc", "n_reps", "tier",
               "abstained", "nsti", "nsti_bin", "picrust2_domain", "routing",
               "rank1_identity", "rank1_rank", "ident_bin", "novelty_bin",
               "n_truth_accessory", "f1_accessory", "n_truth_impossible",
               "template_visible_frac", "min_align_dropped", "unscoreable_reason",
               # --- appended after the SPEC 8.1 contract (never reordered, so the
               # compare.py renderer keeps working) ------------------------------- #
               # why an arm predicted nothing, verbatim from the runner (or one of
               # MISSING_RECORD_REASONS when the run did not carry this md5 at all)
               "abstain_reason",
               # C10's second ceiling: truth ECs carried by no surviving genome of
               # the cell's OWN domain reference -- what forced-routed PICRUSt2 can
               # physically reach, as opposed to prFBA's bac+arc union
               "n_truth_impossible_own_domain"]

# A prediction table that does not carry an md5 the fold plan says was queried is a
# RUN DEFECT (truncated / resumed-stale / plan-mismatched), never a method
# abstention.  These sentinel reasons are written to the cells table's
# ``abstain_reason`` column and counted by C14; they must stay disjoint from the
# runner-emitted vocabulary (run_prfba_fold writes abstain_below_family /
# abstain_reliability).
MISSING_RECORD_REASONS = ("no_record", "no_prediction", "no_forced_record")

# The smallest cluster bootstrap that may stand behind a published CI.  --bootstrap 0
# produced a quotable artefact whose primary sentence read "SESOI NOT exceeded" with
# CI [None, None] at EXIT 0 -- i.e. "no interval was computed" was published in the
# same words as "the effect is below the smallest effect size of interest" -- and
# --bootstrap 5 produced a CI [-0.2895, -0.2224] from five resamples with no warning.
MIN_BOOTSTRAP = 200

# C12's gate number: a pipeline whose sanity gate excludes more than this fraction of
# ALL models is miscalibrated, not protective.  A gated TRUTH model is zero-tolerance
# (it removes that genome from every arm's denominator).
C12_MAX_GATED_FRACTION = 0.10

# Every file whose BYTES decide a number in the artefact, beyond the run directories:
# the ANSWER KEY itself (--functions, the truth EC set behind every F1) and the shared
# vocabulary (--maps, --ec-universe, --fanout, --templates, --annotation-ceiling).  A
# truth file with 2,822 ECs dropped from 20 held-out genomes re-published the shipped
# headline verbatim through the no---force reuse path, and with --force moved picrust2
# 0.6955 -> 0.6872 with nothing in the artefact naming the key that produced it.
ANSWER_KEY_INPUTS = ("functions", "maps", "ec_universe", "fanout",
                     "annotation_ceiling", "templates")

# The scorer's OWN code identity, fingerprinted from a FIXED list (never from
# sys.modules, which differs between the pre-scoring reuse check and the post-scoring
# sidecar because load_template_visible imports phenotype_arm lazily).  A modified
# scorer (picrust2 f1 * 0.5 in build_row) re-published 0.6955 / delta -0.2610 through
# the no---force reuse path while its own scoring gives 0.3478 / +0.0867 -- the primary
# conclusion's SIGN.  The sidecar recorded git_commit but nobody compared it, and
# rev-parse HEAD on a dirty tree cannot see an uncommitted edit anyway.
SCORER_CODE_FILES = ("score_folds.py", "design_audit_baselines.py",
                     "design_audit_icc.py", "design_audit_riskcoverage.py",
                     "score_genecap.py", "compare.py", "annotation_ceiling.py",
                     "_common.py", "phenotype_arm.py")

# Controls every mode must EMIT, and absence voids the run.  An absent control was
# indistinguishable from a passing one: a silently unemitted C14 read as green.  Only
# controls whose emission depends on NOTHING but the mode are listed -- a gate that
# cried wolf on a legitimately absent control would be worse than the hole, so the
# data-conditional ones are REPORTED instead (CONDITIONAL_CONTROLS below).
EXPECTED_CONTROLS = {
    "fold": ("C3", "C4", "C6", "C7", "C10", "C14"),
    "density": ("C3", "C4", "C6", "C7", "C10", "C14"),
    "include_self": ("C2", "C3", "C4", "C6", "C7", "C10", "C14"),
    "plumbing": ("C3", "C9"),
    "recipe_control": ("C1",),
}
CONDITIONAL_CONTROLS = {
    "C1": "--mode recipe_control",
    "C2": "--mode include_self",
    "C5": "a pruned stock reference (recipe arm A) under --refs",
    "C8": "at least one prFBA run (--prfba)",
    "C9": "a plumbing reference run under --picrust2",
    "C11": "phenotype_summary.json present (and bindable)",
    "C12": "phenotype_summary.json present (and bindable)",
    "C13": "a reference manifest carrying tip counts",
    "C15": "a reference manifest carrying the trait-column screen",
}
ALL_CONTROLS = tuple("C%d" % i for i in range(1, 16))

TEMPLATE_DIR = "/home/freiburger/Documents/ModelSEEDTemplates/templates/v7.0/"
TEMPLATE_FILES = {"gramneg": "GramNegModelTemplateV7.json",
                  "grampos": "GramPosModelTemplateV7.json",
                  "archaea": "ArchaeaTemplateV6.json"}
EC_IN_ROLE = re.compile(r"\(EC:?\s*([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)\)", re.I)


# --------------------------------------------------------------------------- #
# small utilities
# --------------------------------------------------------------------------- #
def jload(path):
    with open(path, "rb") as fh:
        return orjson.loads(fh.read())


def jdump(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(orjson.dumps(obj, option=orjson.OPT_INDENT_2 |
                              orjson.OPT_SERIALIZE_NUMPY | orjson.OPT_NON_STR_KEYS))
    os.replace(tmp, path)


def now_iso():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def git_commit():
    try:
        return subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"],
                              capture_output=True, text=True, timeout=20).stdout.strip()
    except Exception:
        return ""


def dkey(d):
    """Canonical density key: 1.0 -> '1.0', 0.5 -> '0.5', 0.25 -> '0.25'."""
    return str(float(d))


def norm_ec(e):
    """'EC:1.1.1.1' -> '1.1.1.1'; partial ECs (containing '-') are dropped."""
    e = str(e).strip()
    if e.upper().startswith("EC:"):
        e = e[3:].strip()
    return None if ("-" in e or not e) else e


def rnd(v, k=4):
    return None if v is None else round(float(v), k)


def fmt(v, k=4):
    """CSV cell: blank for None, rounded otherwise."""
    return "" if v is None else (round(float(v), k) if isinstance(v, float) else v)


def stable_hash(*parts):
    """Process-independent hash. ``hash()`` is salted by PYTHONHASHSEED, so using it
    for a seed would make every run's donor and clade-max draws different."""
    import hashlib
    h = hashlib.blake2b("\x1f".join(str(p) for p in parts).encode(), digest_size=8)
    return int.from_bytes(h.digest(), "big")


def seeded_rng(*parts):
    return random.Random((FROZEN["seed"] * 1000003 + stable_hash(*parts)) & 0x7FFFFFFF)


def gzip_text_writer(path):
    """gzip text writer with mtime=0, so an unchanged run reproduces the exact bytes."""
    import io
    raw = open(path, "wb")
    gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return io.TextIOWrapper(gz, encoding="utf-8", newline=""), (gz, raw)


def mcc_from_counts(tp, fp, fn, n_universe):
    tn = n_universe - tp - fp - fn
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return ((tp * tn - fp * fn) / den) if den > 0 else None


def macro(vals, n_total, population, extra_zeros=0):
    """{'value','n','sd'} over ``vals``, padded with zeros in the population view."""
    vals = list(vals)
    if population:
        vals = vals + [0.0] * extra_zeros
    if not vals:
        return None
    return {"value": round(st.mean(vals), 4), "n": len(vals),
            "sd": round(st.pstdev(vals), 4) if len(vals) > 1 else 0.0}


def lsq_slope(xs, ys):
    if len(xs) < 2 or len(set(xs)) < 2:
        return None
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    return float(np.polyfit(x, y, 1)[0])


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values, order preserved; None passes through."""
    idx = [i for i, p in enumerate(pvals) if p is not None]
    out = [None] * len(pvals)
    if not idx:
        return out
    try:
        from statsmodels.stats.multitest import multipletests
        q = multipletests([pvals[i] for i in idx], method="fdr_bh")[1]
    except Exception:
        m = len(idx)
        order = sorted(idx, key=lambda i: pvals[i])
        q, running = [0.0] * m, 1.0
        tmp = {}
        for rank, i in enumerate(reversed(order), start=1):
            running = min(running, pvals[i] * m / (m - rank + 1))
            tmp[i] = running
        q = [tmp[i] for i in idx]
    for i, v in zip(idx, q):
        out[i] = float(v)
    return out


# --------------------------------------------------------------------------- #
# ArmAgg -- compare.Agg plus population precision, MCC and template visibility
# --------------------------------------------------------------------------- #
class ArmAgg(Agg):
    """``compare.Agg`` extended with the two quantities SPEC Section 8.5 adds.

    ``Agg`` pads only recall/f1/jaccard with zeros for the population view; the spec's
    population block also carries precision and MCC, so those are tracked here and
    padded identically (an abstention has precision 0 and MCC 0 -- it predicts nothing
    and therefore correlates with nothing).
    """

    __slots__ = ("mcc_vals", "tv_vals")

    def __init__(self):
        super().__init__()
        self.mcc_vals = []
        self.tv_vals = []

    def add_cell(self, truth_n, pred_n, tp, recall, precision, f1, jaccard, mcc, tv):
        self.add(truth_n, pred_n, tp, recall, precision, f1, jaccard)
        if pred_n and mcc is not None:
            self.mcc_vals.append(mcc)
        if tv is not None:
            self.tv_vals.append(tv)

    def _pop(self, key):
        vals = self.mcc_vals if key == "mcc" else self.vals[key]
        return macro(vals, self.n, True, extra_zeros=self.n - self.called)

    def _cond(self, key):
        vals = self.mcc_vals if key == "mcc" else self.vals[key]
        return macro(vals, self.called, False)

    def population_block(self):
        return {k: self._pop(k) for k in ("recall", "precision", "f1", "jaccard", "mcc")}

    def conditional_block(self):
        """Conditional means with BOTH denominators printed: ``n_matched`` is the size
        of the matched subset this block was restricted to, ``n_answered`` how many of
        those cells this arm actually answered (they differ for an arm that abstains
        inside the matched set, e.g. a donor baseline with no surviving donor)."""
        b = {k: self._cond(k) for k in ("recall", "precision", "f1", "jaccard", "mcc")}
        b["n_matched"] = self.n
        b["n_answered"] = self.called
        return b

    def micro_block(self):
        r = self.tp / (self.tp + self.fn) if (self.tp + self.fn) else None
        p = self.tp / (self.tp + self.fp) if (self.tp + self.fp) else None
        f1 = (2 * r * p / (r + p)) if (r and p) else (0.0 if r is not None else None)
        return {"recall": rnd(r), "precision": rnd(p), "f1": rnd(f1)}

    def coverage(self):
        return round(self.called / self.n, 4) if self.n else None

    def tv_mean(self):
        return round(st.mean(self.tv_vals), 4) if self.tv_vals else None


# --------------------------------------------------------------------------- #
# input loading
# --------------------------------------------------------------------------- #
def load_panel(path):
    d = jload(path)
    return d.get("genomes", {}), d.get("meta", {})


def load_ec_universe(path):
    """The FIXED scoring vocabulary. ``ec_index`` is the builder's column order, so
    every arm, every fold and every trait table share one integer id per EC."""
    d = jload(path)
    ec_list = list(d["ec"])
    if d.get("ec_index"):
        ec_index = {e: int(i) for e, i in d["ec_index"].items()}
    else:
        ec_index = {e: i for i, e in enumerate(ec_list)}
    prevalence = d.get("prevalence", {})
    return ec_list, ec_index, prevalence, d.get("meta", {})


def load_fanout(path):
    """{cell key -> md5} from the md5 -> [cells] fanout."""
    fan = jload(path)
    out = {}
    for md5, cells in fan.items():
        for c in cells:
            out[c] = md5
    return out


def load_functions(path, keep_ids, ec_index):
    """{genome_id -> frozenset(ec int)} restricted to ``keep_ids`` and the universe."""
    raw = jload(path)
    out = {}
    for gid in keep_ids:
        rec = raw.get(gid)
        if not rec:
            continue
        ecs = set()
        for e in rec.get("ec", ()):
            e = norm_ec(e)
            if e is not None and e in ec_index:
                ecs.add(ec_index[e])
        out[gid] = frozenset(ecs)
    del raw
    return out


def load_ec_rxn(path, ec_index):
    maps = jload(path)
    rxn_index = {}
    ec_rxn = {}
    for e, rxns in sorted(maps.get("ec_rxn", {}).items()):
        e = norm_ec(e)
        if e is None or e not in ec_index:
            continue
        ids = set()
        for r in rxns:
            ids.add(rxn_index.setdefault(r, len(rxn_index)))
        ec_rxn[ec_index[e]] = frozenset(ids)
    return ec_rxn, rxn_index


def project(ec_set, ec_rxn):
    out = set()
    for e in ec_set:
        out |= ec_rxn.get(e, frozenset())
    return frozenset(out)


def load_template_visible(template_dir, ec_index):
    """{template key -> frozenset(ec int)} named by at least one template role.

    Read straight from the template JSONs with the same regex ``phenotype_arm`` uses,
    so no modelseedpy import (and no 20 s template build) is needed just to log D16's
    template-visible fraction.
    """
    out = {}
    try:
        import phenotype_arm as PA           # noqa: F401  (reuse its regex if available)
        rx = PA.EC_RE
    except Exception:
        rx = EC_IN_ROLE
    for key, fname in sorted(TEMPLATE_FILES.items()):
        path = os.path.join(template_dir, fname)
        if not os.path.exists(path):
            continue
        try:
            d = jload(path)
        except Exception:
            continue
        ecs = set()
        for r in d.get("roles", ()):
            for e in rx.findall(r.get("name", "")):
                e = norm_ec(e)
                if e is not None and e in ec_index:
                    ecs.add(ec_index[e])
        out[key] = frozenset(ecs)
    return out


def dom3(x):
    """'Archaea' / 'arc' / 'ARCHAEA' -> 'arc'; 'Bacteria' -> 'bac'; anything that is
    not a single domain ('both', 'all', None, '') -> '' (meaning "covers everything").
    """
    s = str(x or "").strip().lower()
    if s in ("", "both", "all", "none", "any", "union"):
        return ""
    return s[:3]


def _meta_of(path):
    try:
        return jload(path).get("meta", {})
    except Exception:
        return {}


def discover_picrust2(root, ec_index, extend=False):
    """[{dir, name, fold, domain, density, routing, pred, meta}] over */predictions.json.

    Falls back to ``compare.load_picrust2`` for a directory that holds only the two
    tsv.gz tables (the raw hsp.py products) and no predictions.json.

    With ``extend=False`` (the fold arm) an EC outside ``panel_ec_universe`` is dropped:
    the vocabulary is FIXED and identical for every arm, so a prediction outside it is
    not expressible.  With ``extend=True`` (the recipe control, whose vocabulary is the
    stock trait table rather than our panel) unseen ECs are added to the index instead,
    because dropping them would silently inflate precision.
    """
    runs = []
    if not root or not os.path.isdir(root):
        return runs
    for d in sorted(glob.glob(os.path.join(root, "*"))):
        if not os.path.isdir(d):
            continue
        name = os.path.basename(d)
        pj = os.path.join(d, "predictions.json")
        if os.path.exists(pj):
            blob = jload(pj)
            meta = blob.get("meta", {})
            pred = {}
            for md5, rec in blob.get("pred", {}).items():
                ecs = set()
                for e in rec.get("ec", ()):
                    e = norm_ec(e)
                    if e is None:
                        continue
                    if extend:
                        ecs.add(ec_index.setdefault(e, len(ec_index)))
                    elif e in ec_index:
                        ecs.add(ec_index[e])
                pred[str(md5)] = {"ec": frozenset(ecs), "nsti": rec.get("nsti"),
                                  "best_domain": rec.get("best_domain"),
                                  "closest": rec.get("closest_reference_genome")}
        elif os.path.exists(os.path.join(d, "combined_EC_predicted.tsv.gz")):
            local = dict(ec_index)
            raw, p2meta = load_picrust2(d, local)
            inv = {v: k for k, v in local.items()}
            pred = {}
            for md5, cols in raw.items():
                ecs = {c for c in cols if inv.get(c) in ec_index}
                nsti, dom, closest = p2meta.get(md5, (None, None, None))
                pred[str(md5)] = {"ec": frozenset(ecs), "nsti": nsti,
                                  "best_domain": dom, "closest": closest}
            meta = {"ref_name": name, "routing": "forced"}
        else:
            continue
        # run_picrust2_fold writes meta.ref_name as a STRING under --routing forced
        # and as a LIST of both domain references under --routing nsti.  Normalise to
        # a list here: every downstream consumer (manifest lookup, --mode selection,
        # the plumbing control) must work for both shapes, and dict.get(list) raises.
        rn = meta.get("ref_name") or name
        ref_names = list(rn) if isinstance(rn, (list, tuple)) else [rn]
        ref_names = [str(x) for x in ref_names if x]
        runs.append({"dir": d, "name": name, "meta": meta, "pred": pred,
                     "fold": meta.get("fold"), "domain": meta.get("domain"),
                     "density": float(meta.get("density") or PRIMARY_DENSITY),
                     "routing": meta.get("routing") or "forced",
                     "ref_name": rn, "ref_names": ref_names or [name]})
    return runs


def discover_prfba(root, ec_index):
    """[{dir, name, fold, density, forced, modes:{mode:{md5:rec}}, meta}]."""
    runs = []
    if not root or not os.path.isdir(root):
        return runs
    for d in sorted(glob.glob(os.path.join(root, "*"))):
        pj = os.path.join(d, "predictions.json")
        if not os.path.isdir(d) or not os.path.exists(pj):
            continue
        blob = jload(pj)
        meta = blob.get("meta", {})
        modes = {}
        for mode, per_md5 in blob.get("pred", {}).items():
            cur = {}
            for md5, rec in per_md5.items():
                ecs = set()
                for e in rec.get("ec", ()):
                    e = norm_ec(e)
                    if e is not None and e in ec_index:
                        ecs.add(ec_index[e])
                cur[str(md5)] = {"ec": frozenset(ecs), "reps": rec.get("reps") or [],
                                 "n_reps": rec.get("n_reps",
                                                    len(rec.get("reps") or [])),
                                 "tier": rec.get("tier"),
                                 "abstained": bool(rec.get("abstained")),
                                 "abstain_reason": rec.get("abstain_reason") or "",
                                 "best_identity": rec.get("best_identity"),
                                 "rank1_genome": rec.get("rank1_genome"),
                                 "rank1_identity": rec.get("rank1_identity"),
                                 "rank1_rank": rec.get("rank1_rank")}
            modes[mode] = cur
        ref_names = meta.get("ref_name") or []
        if isinstance(ref_names, str):
            ref_names = [ref_names]
        runs.append({"dir": d, "name": os.path.basename(d), "meta": meta,
                     "modes": modes, "fold": meta.get("fold"),
                     "density": float(meta.get("density") or PRIMARY_DENSITY),
                     "forced": bool(meta.get("forced_answer")),
                     "ref_names": ref_names})
    return runs


def load_ref_manifests(refs_dir):
    """{dir name -> manifest dict} for every refs/*/manifest.json."""
    out = {}
    if not refs_dir or not os.path.isdir(refs_dir):
        return out
    for p in sorted(glob.glob(os.path.join(refs_dir, "*", "manifest.json"))):
        try:
            out[os.path.basename(os.path.dirname(p))] = jload(p)
        except Exception:
            continue
    return out


def reinstate_full_prevalence_ec(p2_runs, manifests, ec_index):
    """Credit PICRUSt2 with the ECs its trait table could not physically carry.

    build_reference.py's §5.4 screen must drop a column whose prevalence is 0 OR n
    among the surviving reference tips, because castor's ``mp`` dies on either with
    ``invalid 'ncol' value`` (verified directly on a 39-tip archaeal reference: an
    all-zero column crashes, an all-ONES column crashes identically, and a table of
    only mixed columns succeeds).

    SPEC §5.4/§8.1 then say every dropped EC is "scored as predicted-0 for PICRUSt2".
    That is right for ``reason == 'zero'`` -- castor would have reconstructed absence
    anyway -- but wrong for ``reason == 'full'``: a trait carried by *every* reference
    tip reconstructs as present at every node, so castor would have predicted 1 for
    every query.  Scoring it 0 is a guaranteed false negative on an EC the held-out
    genome very likely has, i.e. a systematic recall penalty against PICRUSt2 that
    grows as the reference thins.  Measured here: 7 such ECs on the 39-tip archaeal
    reference vs 0 on the 143-tip bacterial one -- biting hardest in exactly the
    domain where the two methods are closest.

    D12's stated intent is to drop "only the crashing columns" without shrinking the
    predictable vocabulary, so the honest repair is on the scoring side: re-add the
    ``'full'`` ECs to the prediction set.  Returns per-run counts for the summary.

    The reinstated set is resolved PER DOMAIN, never per run: the 'full' set is a
    property of one domain's trait table and the two domains are very unequal (7/10/25
    archaeal vs 0/1/2 bacterial columns in the smoke references).  A forced-routing
    run has exactly one reference, so its set applies to every record; an nsti-routing
    run carries BOTH domain references (``meta.ref_name`` is a list), and each record
    must be credited with the set of the reference that actually produced it
    (``best_domain``) -- crediting a bacterium-routed query with the archaeal
    reference's 'full' ECs would manufacture predictions out of the other domain.
    """
    added = {}
    for r in p2_runs:
        names = r.get("ref_names") or [r.get("name")]
        per_dom, dom_names = {}, {}
        for name in names:
            man = manifests.get(name) or {}
            if not man and len(names) == 1:
                man = manifests.get(r.get("name")) or {}
            reasons = man.get("dropped_reason") or {}
            full = set()
            for ec, why in reasons.items():
                if why != "full":
                    continue
                e = norm_ec(ec)
                if e is not None and e in ec_index:
                    full.add(ec_index[e])
            dom = dom3(man.get("domain")) or dom3(r.get("domain"))
            # two references of the same domain (or an unlabelled one) can only be
            # merged, but that never happens for a routed run: one per domain.
            per_dom[dom] = per_dom.get(dom, set()) | full
            dom_names[dom] = name
        n_unresolved = 0
        if len(per_dom) <= 1:
            full = next(iter(per_dom.values()), set())
            if full:
                for rec in r["pred"].values():
                    rec["ec"] = frozenset(rec["ec"] | full)
            added[r["name"]] = {dom_names.get(next(iter(per_dom), ""), r["name"]):
                                len(full)} if full else {}
        else:
            for rec in r["pred"].values():
                full = per_dom.get(dom3(rec.get("best_domain")))
                if full is None:
                    n_unresolved += 1
                    continue
                if full:
                    rec["ec"] = frozenset(rec["ec"] | full)
            added[r["name"]] = {dom_names[d]: len(s)
                                for d, s in sorted(per_dom.items()) if s}
        r["full_prevalence_ec_reinstated"] = added[r["name"]]
        r["n_full_prevalence_ec_reinstated"] = sum(added[r["name"]].values())
        r["n_full_prevalence_domain_unresolved"] = n_unresolved
        if n_unresolved:
            print("[score_folds] WARNING: %s has %d prediction(s) whose best_domain "
                  "does not match any of its references %s; their prevalence-n ECs "
                  "were NOT reinstated (conservative: a recall penalty against "
                  "PICRUSt2, never a credit)" % (r["name"], n_unresolved,
                                                 sorted(per_dom)))
    return added


def run_kind(ref_names):
    """Which reference family a run used -> the --mode it belongs to.

    ``ref_names`` is always the normalised LIST (a --routing nsti run names both
    domain references); ``str(ref_name)`` on the raw field would stringify the list
    to "['f00_bac_d100', ...]", match no prefix, and silently classify an nsti
    control run as a fold run.
    """
    names = [str(n) for n in (ref_names or ()) if n]
    kinds = set()
    for n in names:
        if n.startswith("incself"):
            kinds.add("include_self")
        elif n.startswith("plumbing"):
            kinds.add("plumbing")
        elif n.startswith("recipe"):
            kinds.add("recipe_control")
        else:
            kinds.add("fold")
    if not kinds:
        return "fold"
    if len(kinds) > 1:
        raise SystemExit("[score_folds] a run mixes reference families %s (%s): the "
                         "two domain references of one run must belong to the same "
                         "--mode" % (sorted(kinds), names))
    return kinds.pop()


def select_runs(mode, p2_runs, pf_runs):
    """Keep only the runs a given --mode is about (the reference family they used)."""
    want = "fold" if mode == "density" else mode
    p2 = [r for r in p2_runs if run_kind(r.get("ref_names")) == want]
    pf = [r for r in pf_runs if run_kind(r.get("ref_names")) == want]
    return p2, pf


# --------------------------------------------------------------------------- #
# Stage 0 -- the run INVENTORY gate (C14 legs 6-7)
# --------------------------------------------------------------------------- #
# Every silent-population defect this harness has shipped has the same shape: the
# scorer scores whatever run directories it happens to FIND, so every way a run can
# be missing, duplicated, stale or under-specified shrinks the scored population
# while the headline stays valid at EXIT 0.  Patching leg by leg did not converge
# (a killed (fold, domain) PICRUSt2 run with its refs/ still present; a stale retry
# directory merged into the good run -- that one moved the primary delta from
# -0.2610 to -0.2039; a never-launched archaeal arm, -10.4% of the population; a
# fold launched with a shorter --select-modes; a whole lost density tier).
#
# So the check is INVERTED here: the expected inventory is DERIVED from the fold
# plan and the frozen config (folds x densities x domains x arm/mode) plus this
# invocation's own declared scope, never from what is on disk, and the discovered
# run directories must equal it EXACTLY.  The derivation is recorded beside the
# verdict so a reader can see which axis came from where.
def _ikey(arm, fold, dk, domain="", extra=""):
    return "%s|fold=%s|d=%s|domain=%s|%s" % (arm, fold, dk, domain or "-", extra)


def file_fp(path):
    """{size, sha256} of a file the scorer is about to trust, or None.

    A reviewer lost half a day to a concurrent session rewriting a predictions.json
    mid-review (it moved picrust2 0.6955 -> 0.6526 with every control green), so
    every input the numbers rest on is fingerprinted INTO the artefact.
    """
    import hashlib
    try:
        h, n = hashlib.sha256(), 0
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
                n += len(chunk)
        return {"size": n, "sha256": h.hexdigest()}
    except OSError:
        return None


def answer_key_paths(args, ec_path=None):
    """The resolved paths of the answer key and the shared vocabulary, by tag.

    Resolved here rather than at each use site so the artefact can NAME the truth
    file it scored against (``inventory.input_paths``); ``--ec-universe`` defaults to
    a file beside ``--panel``, which no reader could otherwise reconstruct.
    """
    ecp = ec_path or getattr(args, "ec_universe", None) or os.path.join(
        os.path.dirname(os.path.abspath(getattr(args, "panel", "") or ".")),
        "panel_ec_universe.json")
    out = {"functions": getattr(args, "functions", None),
           "maps": getattr(args, "maps", None),
           "ec_universe": ecp,
           "fanout": getattr(args, "fanout", None),
           "annotation_ceiling": getattr(args, "annotation_ceiling", None),
           "templates": (None if getattr(args, "no_templates", False)
                         else getattr(args, "templates", None))}
    return {k: (os.path.abspath(v) if v else None) for k, v in sorted(out.items())}


def _answer_key_fingerprints(args, ec_path=None):
    """{tag -> {size, sha256}} for every answer-key / vocabulary input.

    ``--templates`` is a DIRECTORY, so the exact files ``load_template_visible``
    reads (``TEMPLATE_FILES``, never a glob -- a glob would make the fingerprint
    depend on unrelated files landing in the template directory) are fingerprinted
    individually under ``templates/<role>``.
    """
    paths = answer_key_paths(args, ec_path)
    out = {}
    for tag in ANSWER_KEY_INPUTS:
        p = paths.get(tag)
        if not p:
            continue
        if tag == "templates":
            for role, fname in sorted(TEMPLATE_FILES.items()):
                fp = file_fp(os.path.join(p, fname))
                if fp:
                    out["templates/%s" % role] = fp
            continue
        fp = file_fp(p)
        if fp:
            out[tag] = fp
    return out


def collect_fingerprints(args, ec_path=None):
    """Size + sha256 of EVERY input these numbers rest on, by section.

    One collector for the Stage 0 inventory and for the post-cells re-verification:
    they were two copies of the same loop, and a file added to one but not the other
    would be fingerprinted at the start and never re-checked.
    """
    fps = {"predictions": {}, "manifests": {}, "inputs": {},
           "answer_key": _answer_key_fingerprints(args, ec_path)}
    for root, tag in ((getattr(args, "picrust2", None), "picrust2"),
                      (getattr(args, "prfba", None), "prfba")):
        if not root or not os.path.isdir(root):
            continue
        for p in sorted(glob.glob(os.path.join(root, "*", "predictions.json"))):
            fps["predictions"]["%s/%s" % (tag, os.path.basename(os.path.dirname(p)))] \
                = file_fp(p)
    refs = getattr(args, "refs", None)
    if refs and os.path.isdir(refs):
        for p in sorted(glob.glob(os.path.join(refs, "*", "manifest.json"))):
            fps["manifests"][os.path.basename(os.path.dirname(p))] = file_fp(p)
    for tag in ("fold_plan", "panel", "truth", "neighbors"):
        p = getattr(args, tag, None)
        if p and os.path.exists(p):
            fps["inputs"][tag] = file_fp(p)
    return fps


def fp_digest(obj):
    import hashlib
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def scorer_code_fp():
    """(digest, {name -> {path, size, sha256}}) of the code that computes the numbers.

    Fixed file list, fixed search order, paths recorded relative to the repo root so
    the digest is identical from any checkout location.  This is what binds a reused
    ``--cells`` table to the scorer that wrote it; ``git_commit`` cannot (it is
    ``rev-parse HEAD`` on a tree that may be dirty).
    """
    files = {}
    for name in SCORER_CODE_FILES:
        for d in (HERE, os.path.join(RV, "scripts"), os.path.join(P2C, "scripts"),
                  REPO):
            p = os.path.join(d, name)
            if os.path.exists(p):
                fp = file_fp(p) or {}
                files[name] = {"path": os.path.relpath(p, REPO),
                               "size": fp.get("size"), "sha256": fp.get("sha256")}
                break
    return fp_digest(files), files


def _parse_dens_list(s):
    out = []
    for tok in str(s or "").replace(" ", "").split(","):
        if tok:
            out.append(dkey(float(tok)))
    return sorted(set(out), key=float)


_P2_NAME = re.compile(r"^f(\d+)_(bac|arc)_d(\d+)$")
_PF_NAME = re.compile(r"^f(\d+)_d(\d+)(_forced)?$")


def plan_domains_by_fold(plan, panel, truth_meta):
    """{fold -> {'bac','arc'}} from the PLAN's own test cells.

    The domains that must have a forced-routing PICRUSt2 run are the domains the
    plan's test cells belong to -- not the domains whose reference happens to be on
    disk.  Deriving them from refs/ would miss the case where the archaeal reference
    build died too: no ref, no run, every archaeal cell silently
    ``unscoreable_method_not_run`` for every arm, equal denominators, EXIT 0.
    """
    out = {}
    for f in plan.get("folds", []):
        if f.get("fold") is None:
            continue
        doms = set()
        for cell in (f.get("test_cells") or ()):
            t = truth_meta.get(cell) or {}
            d = t.get("domain") or (panel.get(t.get("src_genome_id") or
                                              cell.rsplit("__", 1)[0]) or {}
                                    ).get("domain") or ""
            d3 = dom3(d)
            if d3:
                doms.add(d3)
        out[int(f["fold"])] = doms
    return out


def run_inventory(args, plan, panel, truth_meta, manifests, p2_runs, pf_runs,
                  frozen, ec_path=None):
    """DERIVE the expected run set, ASSERT the discovered set equals it, FINGERPRINT
    every input.  Returns the ``accounting.completeness.inventory`` block.

    Expected, at full scale (SPEC 11.2): 60 PICRUSt2 forced (10 folds x 2 domains x
    3 densities) + 10 nsti + 30 prFBA fold + 10 prFBA forced-answer, each EXACTLY
    once.  Derivation, axis by axis:

    * folds     -- ``fold_plan.folds`` (every planned fold, always).
    * densities -- ``--expect-densities`` when declared; else the frozen
      ``density_fractions`` (``--mode density``, where the sweep IS the claim); else
      the frozen fractions intersected with the density tiers that have a reference
      BUILD (``--mode fold``: a tier nobody built was not part of this invocation's
      scope, and the frozen-vs-built gap is reported and printed).
    * domains   -- ``--expect-domains`` when declared; else the domains of the
      plan's own test cells, per fold (``plan_domains_by_fold``).
    * arms      -- one forced PICRUSt2 run per (fold, density, domain); one prFBA
      fold run per (fold, density); one prFBA forced-answer run per fold at the
      primary density (``--expect-forced-answer``); nsti runs only when
      ``--expect-nsti`` declares them (default: none expected, and then a
      discovered nsti run is an EXTRA, because an inventory inferred from disk is
      the defect, not the gate).

    The gate is enforced in ``--mode fold`` / ``--mode density``.  The control
    reference families (include_self, plumbing) are single-fold/single-density by
    construction, so there the expected set is taken from the discovered runs and
    the block is REPORTED with ``ok: null`` -- except for duplicates, which are a
    hard failure in every mode: nothing legitimately ships two directories for one
    (fold, density, domain, arm) key, and the merge is silent (PICRUSt2 per md5 by
    min NSTI, prFBA by dict.update last-wins).
    """
    primary_dk = dkey(PRIMARY_DENSITY)
    gated = args.mode in ("fold", "density")
    plan_folds = sorted(int(f["fold"]) for f in plan.get("folds", [])
                        if f.get("fold") is not None)
    frozen_dens = _parse_dens_list(",".join(
        str(x) for x in (frozen.get("density_fractions") or [PRIMARY_DENSITY])))
    # density tiers with a FOLD-family reference build (refs/ is the builder's
    # record, not a run's self-report)
    built_dens, built_grid = set(), collections.defaultdict(set)
    for name, man in sorted(manifests.items()):
        if (man.get("mode") or "fold") != "fold" or man.get("fold") is None:
            continue
        m = _P2_NAME.match(name)
        d = man.get("density")
        dk = dkey(float(d)) if d is not None else (
            dkey(int(m.group(3)) / 100.0) if m else None)
        if dk is None:
            continue
        built_dens.add(dk)
        built_grid[(int(man["fold"]), dk)].add(dom3(man.get("domain")))
    built_dens = sorted(built_dens, key=float)

    declared_dens = _parse_dens_list(args.expect_densities) if args.expect_densities \
        else []
    dens_source = "--expect-densities"
    if declared_dens:
        exp_dens = declared_dens
    elif args.mode == "density":
        exp_dens, dens_source = list(frozen_dens), "frozen.density_fractions"
    elif args.mode == "fold":
        exp_dens = [d for d in frozen_dens if d in set(built_dens)] or [primary_dk]
        dens_source = "frozen.density_fractions INTERSECT reference builds in --refs"
    else:
        exp_dens = sorted({dkey(r["density"]) for r in p2_runs} |
                          {dkey(r["density"]) for r in pf_runs}, key=float) or \
            [primary_dk]
        dens_source = "discovered runs (control reference family; report-only)"

    declared_doms = sorted({dom3(x) for x in
                            str(args.expect_domains or "").replace(" ", "").split(",")
                            if x}) if args.expect_domains else []
    doms_by_fold = plan_domains_by_fold(plan, panel, truth_meta)
    dom_source = ("--expect-domains" if declared_doms else
                  "domains of the plan's own test cells, per fold")

    # ---- the expected set ------------------------------------------------- #
    declared_folds = sorted({int(x) for x in
                             str(args.expect_folds or "").replace(" ", "").split(",")
                             if x}) if args.expect_folds else []
    unknown_folds = [f for f in declared_folds if f not in set(plan_folds)]
    if unknown_folds:
        raise SystemExit("[score_folds] --expect-folds names fold(s) %s that the "
                         "fold plan does not contain (%s)"
                         % (unknown_folds, plan_folds))
    expected = set()
    if gated:
        folds_scope = declared_folds or list(plan_folds)
        for f in folds_scope:
            doms = sorted(declared_doms or doms_by_fold.get(f, set()))
            for dk in exp_dens:
                for dom in doms:
                    expected.add(("picrust2", f, dk, dom, "forced"))
                expected.add(("prfba", f, dk, "", "fold"))
            if args.expect_nsti == "primary":
                expected.add(("picrust2", f, primary_dk, "", "nsti"))
            elif args.expect_nsti == "all":
                for dk in exp_dens:
                    expected.add(("picrust2", f, dk, "", "nsti"))
            if args.expect_forced_answer == "primary":
                expected.add(("prfba", f, primary_dk, "", "forced_answer"))
            elif args.expect_forced_answer == "all":
                for dk in exp_dens:
                    expected.add(("prfba", f, dk, "", "forced_answer"))
    else:
        folds_scope = sorted({int(r["fold"]) for r in list(p2_runs) + list(pf_runs)
                              if r.get("fold") is not None})

    # ---- the discovered set ----------------------------------------------- #
    discovered = collections.defaultdict(list)
    scope_mismatch = []
    for r in p2_runs:
        meta = r.get("meta") or {}
        routing = str(r.get("routing") or "forced")
        d3 = dom3(meta.get("domain"))
        key = ("picrust2", r.get("fold"), dkey(r["density"]),
               d3 if routing == "forced" else "", routing)
        discovered[key].append(r["dir"])
        if gated:
            if routing == "forced" and not d3:
                scope_mismatch.append({
                    "run": r["name"], "dir": r["dir"],
                    "defect": "forced_run_without_a_domain",
                    "declared": {"domain": meta.get("domain"), "routing": routing},
                    "planned": "every forced-routing PICRUSt2 run names its own "
                               "domain (bac|arc); an unlabelled run is scored as "
                               "covering both, which hides a dead domain"})
            if (meta.get("mode") or "fold") != "fold":
                scope_mismatch.append({
                    "run": r["name"], "dir": r["dir"],
                    "defect": "declared_mode_is_not_fold",
                    "declared": {"mode": meta.get("mode")}, "planned": "fold"})
            m = _P2_NAME.match(r["name"])
            if m and (int(m.group(1)) != (r.get("fold") if r.get("fold") is not None
                                          else -1) or
                      m.group(2) != d3 or
                      dkey(int(m.group(3)) / 100.0) != dkey(r["density"])):
                scope_mismatch.append({
                    "run": r["name"], "dir": r["dir"],
                    "defect": "declared_scope_contradicts_directory_name",
                    "declared": {"fold": r.get("fold"), "domain": d3,
                                 "density": dkey(r["density"])},
                    "planned": {"fold": int(m.group(1)), "domain": m.group(2),
                                "density": dkey(int(m.group(3)) / 100.0)}})
    sm_by_run = {}
    for r in pf_runs:
        meta = r.get("meta") or {}
        kind = "forced_answer" if r.get("forced") else "fold"
        key = ("prfba", r.get("fold"), dkey(r["density"]), "", kind)
        discovered[key].append(r["dir"])
        sm_by_run[r["name"]] = sorted({str(m) for m in (meta.get("select_modes") or ())})
        if gated:
            if not sm_by_run[r["name"]]:
                scope_mismatch.append({
                    "run": r["name"], "dir": r["dir"],
                    "defect": "select_modes_absent",
                    "declared": {"select_modes": meta.get("select_modes")},
                    "planned": "meta.select_modes names the run's own scope; without "
                               "it C14 leg 1 has nothing to compare the shipped mode "
                               "blocks against and a vanished mode block passes "
                               "vacuously (it cost prfba_nestedcv_forced 1922 -> "
                               "1258 cells)"})
            m = _PF_NAME.match(r["name"])
            if m and (int(m.group(1)) != (r.get("fold") if r.get("fold") is not None
                                          else -1) or
                      dkey(int(m.group(2)) / 100.0) != dkey(r["density"]) or
                      bool(m.group(3)) != bool(r.get("forced"))):
                scope_mismatch.append({
                    "run": r["name"], "dir": r["dir"],
                    "defect": "declared_scope_contradicts_directory_name",
                    "declared": {"fold": r.get("fold"),
                                 "density": dkey(r["density"]),
                                 "forced_answer": bool(r.get("forced"))},
                    "planned": {"fold": int(m.group(1)),
                                "density": dkey(int(m.group(2)) / 100.0),
                                "forced_answer": bool(m.group(3))}})
    # prFBA select-mode scope must be the SAME on every fold: one fold launched with
    # a shorter --select-modes is internally self-consistent and plan-complete, and
    # it silently drops that arm's whole fold (zero rows, so even the sentinel
    # counters are blind).
    distinct_sm = sorted({tuple(v) for v in sm_by_run.values() if v})
    sm_consistent = len(distinct_sm) <= 1
    if gated and not sm_consistent:
        scope_mismatch.append({
            "run": "(prfba, all folds)", "dir": args.prfba,
            "defect": "select_modes_differ_across_runs",
            "declared": {k: v for k, v in sorted(sm_by_run.items())},
            "planned": "every prFBA run of a campaign declares the SAME select-mode "
                       "set; a fold launched with fewer modes drops that arm's whole "
                       "fold from the population"})

    duplicates = [{"key": _ikey(*k), "dirs": sorted(v)}
                  for k, v in sorted(discovered.items(), key=lambda kv: _ikey(*kv[0]))
                  if len(v) > 1]
    missing = sorted(_ikey(*k) for k in (expected - set(discovered))) if gated else []
    unexpected = ([{"key": _ikey(*k), "dirs": sorted(v)}
                   for k, v in sorted(discovered.items(),
                                      key=lambda kv: _ikey(*kv[0]))
                   if k not in expected] if gated else [])

    # ---- fingerprints (provenance, independent of the gate) ---------------- #
    # The ANSWER KEY is in here now (fps["answer_key"]): --functions is the truth EC
    # set behind every F1, and it was fingerprinted nowhere, re-verified nowhere and
    # named nowhere, so a tampered truth file both passed the no---force reuse gate
    # with the old headline and, with --force, moved every number with no record of
    # which key produced it.
    fps = collect_fingerprints(args, ec_path)
    scope_args = {"mode": args.mode, "rank_blocked": args.rank_blocked,
                  "space": args.space, "limit_cells": args.limit_cells,
                  "primary_select_mode": args.primary_select_mode,
                  "no_templates": bool(args.no_templates),
                  "clade_max_cap": args.clade_max_cap}
    # NOTE: the answer key enters the digest through fps["answer_key"] -- its
    # CONTENT (size + sha256), never its path.  Hashing the paths made every
    # fingerprint directory-dependent, so a tree copied to review it refused to
    # reuse its own cells table; the paths are recorded, unhashed, in input_paths.
    inventory_fp = fp_digest({"files": fps, "args": scope_args,
                              "runs": sorted(_ikey(*k) for k in discovered)})
    code_fp, code_files = scorer_code_fp()

    reasons = []
    if duplicates:
        reasons.append("duplicate_run_directories")
    if missing:
        reasons.append("run_inventory_incomplete")
    if unexpected:
        reasons.append("unexpected_run_directory")
    if scope_mismatch:
        reasons.append("run_scope_mismatch")
    ok = (not (missing or unexpected or scope_mismatch or duplicates)) if gated \
        else (False if duplicates else None)
    counts = collections.Counter(k[0] if k[0] == "picrust2" else "prfba"
                                 for k in discovered)
    return {
        "gated": gated,
        "ok": ok,
        "reasons": reasons,
        "derivation": {
            "spec": "SPEC 11.2: 60 PICRUSt2 forced (10 folds x 2 domains x 3 "
                    "densities) + 10 nsti + 30 prFBA fold + 10 prFBA "
                    "forced-answer, each exactly once",
            "folds_source": (("--expect-folds (subset of fold_plan.folds)"
                              if declared_folds else "fold_plan.folds") if gated
                             else "discovered runs (control reference family)"),
            "folds_in_plan": [str(f) for f in plan_folds],
            "folds_expected": [str(f) for f in folds_scope],
            "densities_source": dens_source,
            "densities_expected": exp_dens,
            "densities_in_frozen_plan": frozen_dens,
            "densities_with_a_reference_build": built_dens,
            "densities_in_frozen_plan_not_built": [d for d in frozen_dens
                                                   if d not in set(built_dens)],
            "domains_source": dom_source,
            "domains_expected_by_fold": {str(f): sorted(declared_doms or
                                                        doms_by_fold.get(f, set()))
                                         for f in folds_scope},
            "nsti_runs_expected": args.expect_nsti,
            "forced_answer_runs_expected": args.expect_forced_answer,
            "primary_density": primary_dk,
        },
        "n_expected": len(expected),
        "n_discovered": len(discovered),
        "n_run_directories": sum(len(v) for v in discovered.values()),
        "counts_discovered": {"picrust2": counts.get("picrust2", 0),
                              "prfba": counts.get("prfba", 0)},
        "expected": sorted(_ikey(*k) for k in expected),
        "discovered": {_ikey(*k): sorted(v)
                       for k, v in sorted(discovered.items(),
                                          key=lambda kv: _ikey(*kv[0]))},
        "missing": missing,
        "duplicates": duplicates,
        "unexpected": unexpected,
        "scope_mismatch": scope_mismatch,
        "select_modes_by_run": {k: v for k, v in sorted(sm_by_run.items())},
        "select_modes_consistent_across_runs": sm_consistent,
        "select_modes_in_frozen_plan": sorted(
            (frozen.get("prfba") or {}).get("select_modes") or {}),
        "fingerprints": fps,
        "inventory_fp": inventory_fp,
        # WHICH ANSWER KEY, by name and by bytes, and WHICH CODE scored it
        "input_paths": answer_key_paths(args, ec_path),
        "answer_key_fp": fp_digest(fps.get("answer_key") or {}),
        "scorer_code_fp": code_fp,
        "scorer_code_files": code_files,
    }


def _recipe_fingerprints(args):
    """{predictions|manifests|inputs -> {name: {size, sha256}}} for the recipe arm."""
    out = {"predictions": {}, "manifests": {}, "inputs": {}}
    if args.picrust2 and os.path.isdir(args.picrust2):
        for p in sorted(glob.glob(os.path.join(args.picrust2, "*",
                                               "predictions.json"))):
            out["predictions"][os.path.basename(os.path.dirname(p))] = file_fp(p)
    if args.refs and os.path.isdir(args.refs):
        for p in sorted(glob.glob(os.path.join(args.refs, "*", "manifest.json"))):
            out["manifests"][os.path.basename(os.path.dirname(p))] = file_fp(p)
    if args.truth and os.path.exists(args.truth):
        out["inputs"]["truth"] = file_fp(args.truth)
    out["answer_key"] = _answer_key_fingerprints(args)
    out["answer_key_paths"] = answer_key_paths(args)
    out["answer_key_fp"] = fp_digest(out["answer_key"])
    out["scorer_code_fp"], out["scorer_code_files"] = scorer_code_fp()
    return out


def reverify_fingerprints(args, inv, ec_path=None):
    """Re-read every fingerprinted input AFTER the cells pass and compare.

    ``file_fp`` hashes a second read, so on its own it cannot see a rewrite that
    happened between ``discover_*`` and the hash, nor one during the (long) scoring
    pass.  A concurrent session rewriting ONE predictions.json mid-review moved
    picrust2 0.6955 -> 0.6526 with every control green, so the comparison is made
    again at the end and a difference voids the run (C14 leg 6) instead of being
    merely recorded.
    """
    changed = []
    now = collect_fingerprints(args, ec_path)
    before = inv.get("fingerprints") or {}
    for kind in sorted(set(before) | set(now)):
        b, a = before.get(kind) or {}, now.get(kind) or {}
        for name in sorted(set(b) | set(a)):
            if b.get(name) != a.get(name):
                changed.append({"file": "%s:%s" % (kind, name),
                                "before": b.get(name), "after": a.get(name)})
    inv["fingerprints_reverified"] = True
    inv["inputs_changed_during_scoring"] = changed
    if changed:
        inv["ok"] = False
        if "inputs_changed_during_scoring" not in (inv.get("reasons") or []):
            inv.setdefault("reasons", []).append("inputs_changed_during_scoring")
        print("=" * 78)
        print("[score_folds] *** INPUT BYTES CHANGED WHILE SCORING -- %d "
              "fingerprinted input(s) differ from the bytes this run started with. "
              "The numbers are not reproducible from the tree as it now stands and "
              "are stamped invalid. ***" % len(changed))
        for c in changed[:10]:
            print("  %-48s %s -> %s" % (c["file"], c["before"], c["after"]))
        print("=" * 78)
    return changed


def inventory_banner(inv, summary_path=None):
    """Print the inventory verdict the moment it is known -- an operator tailing a
    12-16 h unattended run must not have to wait for the scoring to finish."""
    if inv.get("ok") is not False:
        return
    print("=" * 78)
    print("[score_folds] *** RUN INVENTORY MISMATCH -- the discovered run "
          "directories are NOT the set the fold plan and the frozen config call "
          "for. The numbers below are over a different population than the plan's "
          "and are stamped invalid. ***")
    print("  reasons                = %s" % (inv.get("reasons"),))
    print("  expected / discovered  = %s / %s run key(s)"
          % (inv.get("n_expected"), inv.get("n_discovered")))
    for k in (inv.get("missing") or []):
        print("  MISSING   %s" % k)
    for d in (inv.get("duplicates") or []):
        print("  DUPLICATE %s" % d["key"])
        for p in d["dirs"]:
            print("            %s" % p)
    for d in (inv.get("unexpected") or []):
        print("  UNEXPECTED %s %s" % (d["key"], d["dirs"]))
    for s in (inv.get("scope_mismatch") or []):
        print("  SCOPE     %s %s declared=%s" % (s["defect"], s["run"],
                                                 s.get("declared")))
    der = inv.get("derivation") or {}
    print("  derived from: folds %s (%s), densities %s (%s), domains %s (%s), "
          "nsti=%s, forced_answer=%s"
          % (der.get("folds_expected"), der.get("folds_source"),
             der.get("densities_expected"), der.get("densities_source"),
             der.get("domains_expected_by_fold"), der.get("domains_source"),
             der.get("nsti_runs_expected"), der.get("forced_answer_runs_expected")))
    print("  Fix the inventory (rerun the dead job, delete the stale retry "
          "directory, or declare the scope with --expect-densities / "
          "--expect-domains / --expect-nsti / --expect-forced-answer), then rescore "
          "with --force.")
    print("=" * 78)


def cells_fp_path(cells_path):
    return cells_path + ".fp"


def cells_reuse_verdict(args, inv):
    """May the existing --cells be reused without --force?

    Only when the sidecar fingerprint matches BOTH the cells bytes and the inventory
    that was just inspected.  Reusing on existence alone published stale numbers
    (n=1722, delta -0.2443) from an artefact that openly contradicted its own
    completeness block, with every control passing and EXIT 0.
    """
    fpp = cells_fp_path(args.cells)
    cur = file_fp(args.cells)
    rec = None
    if os.path.exists(fpp):
        try:
            rec = jload(fpp)
        except Exception:
            rec = None
    if rec is None:
        return False, ("no fingerprint sidecar %s: the cells table cannot be shown "
                       "to have been written from the runs now on disk" % fpp)
    if (rec.get("cells") or {}).get("sha256") != (cur or {}).get("sha256"):
        return False, ("cells sha256 %s != the sidecar's %s: %s was rewritten after "
                       "it was scored" % ((cur or {}).get("sha256"),
                                          (rec.get("cells") or {}).get("sha256"),
                                          args.cells))
    if rec.get("inventory_fp") != inv.get("inventory_fp"):
        return False, ("inventory fingerprint %s != the sidecar's %s: the runs, "
                       "manifests, the ANSWER KEY (%s vs %s) or the scoring scope "
                       "changed since the cells table was written"
                       % (inv.get("inventory_fp"), rec.get("inventory_fp"),
                          inv.get("answer_key_fp"), rec.get("answer_key_fp")))
    # ...and to the CODE, which no fingerprint of the inputs can cover: a scorer
    # edited to halve picrust2's f1 re-published 0.6955 / delta -0.2610 from a cells
    # table its own build_row would now score at 0.3478 / +0.0867.
    if not rec.get("scorer_code_fp"):
        return False, ("the fingerprint sidecar %s records no scorer_code_fp (it "
                       "predates the code-identity contract), so the cells table "
                       "cannot be shown to have been written by THIS scorer" % fpp)
    if rec.get("scorer_code_fp") != inv.get("scorer_code_fp"):
        return False, ("scorer code fingerprint %s != the sidecar's %s: "
                       "score_folds.py or one of the promoted modules it computes "
                       "with changed since the cells table was written"
                       % (inv.get("scorer_code_fp"), rec.get("scorer_code_fp")))
    return True, ("fingerprints match (inventory_fp %s, answer_key_fp %s, "
                  "scorer_code_fp %s)"
                  % (inv.get("inventory_fp"), (inv.get("answer_key_fp") or "")[:12],
                     (inv.get("scorer_code_fp") or "")[:12]))


# --------------------------------------------------------------------------- #
# Stage 1 -- score every cell of every fold
# --------------------------------------------------------------------------- #
def reference_set(fold_rec, density, manifests, mode="fold"):
    """Surviving reference genomes for (fold, density): the refs manifests are truth.

    At density < 1 the surviving set is the seeded phylum-stratified subsample recorded
    by the builder, NOT ``fold_plan.reference_genomes``; the union over the domain
    references is exactly what prFBA's restricted DB contains (C4).
    """
    f = int(fold_rec["fold"])
    # a control reference has its own name, and must not borrow the leave-out fold's
    # manifests: the include-self reference contains the test genomes by construction,
    # so prevalence, base_constant and C10 would otherwise be computed from the wrong
    # surviving set.
    pat = {"include_self": "incself", "plumbing": "plumbing",
           "recipe_control": "recipe"}.get(mode, "f%02d_" % f)
    want = "_d%03d" % int(round(density * 100))
    ids, names = set(), []
    by_dom = collections.defaultdict(set)
    for name, man in sorted(manifests.items()):
        if name.startswith(pat) and name.endswith(want) and man.get("reference_genomes"):
            ids |= set(man["reference_genomes"])
            names.append(name)
            by_dom[dom3(man.get("domain"))] |= set(man["reference_genomes"])
    if not ids and pat != "f%02d_" % f:          # control refs may be density-less
        for name, man in sorted(manifests.items()):
            if name.startswith(pat) and man.get("reference_genomes"):
                ids |= set(man["reference_genomes"])
                names.append(name)
                by_dom[dom3(man.get("domain"))] |= set(man["reference_genomes"])
    if not ids:
        ids = set(fold_rec.get("reference_genomes") or ())
        by_dom = {}
    # a domain reference the builder never labelled cannot be told apart from the
    # union, so C10's own-domain ceiling falls back to the union for it (reported)
    return ids, names, {d: v for d, v in by_dom.items() if d}


class FoldContext:
    """Everything derived from one (fold, density) reference set, computed once."""

    def __init__(self, fold_rec, density, ref_ids, ref_names, g_ec, g_rxn,
                 panel, rank_blocked, ref_ids_by_domain=None):
        self.fold = int(fold_rec["fold"])
        self.density = float(density)
        self.ref_ids = sorted(g for g in ref_ids if g in g_ec)
        self.ref_names = ref_names
        self.n_ref = len(self.ref_ids)
        self.sets = {"ec": g_ec, "rxn": g_rxn}
        self.prev = {}
        self.present = {}
        self.accessory = {}
        self.constant = {}
        for space, store in self.sets.items():
            cnt = collections.Counter()
            for gid in self.ref_ids:
                cnt.update(store.get(gid, ()))
            n = max(1, self.n_ref)
            self.prev[space] = {k: v / n for k, v in cnt.items()}
            self.present[space] = frozenset(cnt)
            self.accessory[space] = frozenset(
                BASE.accessory_vocab(self.prev[space], 0.05, 0.95))
            sizes = [len(store.get(g, ())) for g in self.ref_ids] or [0]
            k = int(st.median(sizes))
            self.constant[space] = frozenset(BASE.constant_baseline(cnt, k))
            self.constant_k = k
        # C10, second ceiling: what a FORCED-routed PICRUSt2 query can physically
        # reach is its own domain's tree and trait table, not the bac+arc union that
        # prFBA's restricted ref_db really is.  present_by_domain['arc']['ec'] is the
        # archaeal reference's vocabulary; a domain with no labelled manifest falls
        # back to the union (and is listed in own_domain_fallback).
        self.present_by_domain = {}
        self.own_domain_fallback = []
        for dom, ids in sorted((ref_ids_by_domain or {}).items()):
            sel = [g for g in ids if g in g_ec]
            if not sel:
                continue
            self.present_by_domain[dom] = {
                space: frozenset(e for g in sel for e in store.get(g, ()))
                for space, store in self.sets.items()}
            self.present_by_domain[dom]["n_ref"] = len(sel)
        self.by_rank = {}
        for rank in RANK_ORDER:
            idx = collections.defaultdict(list)
            for gid in self.ref_ids:
                v = (panel.get(gid) or {}).get(rank)
                if v:
                    idx[v].append(gid)
            self.by_rank[rank] = idx
        self.donor_rank = BLOCKED_TO_DONOR.get(rank_blocked, "family")
        self._donor_cache = {}
        self._clademax_cache = {}

    def family_donor(self, src, panel):
        if src in self._donor_cache:
            return self._donor_cache[src]
        rng = seeded_rng("donor", self.fold, self.density, src)
        val = (panel.get(src) or {}).get(self.donor_rank)
        donor, _ = BASE.donor_baseline(val, self.by_rank[self.donor_rank],
                                       {g: self.sets["ec"].get(g, frozenset())
                                        for g in self.by_rank[self.donor_rank].get(val, ())},
                                       {src}, rng)
        self._donor_cache[src] = donor
        return donor

    def clade_max(self, src, panel, truth_ec, cap):
        """Max-F1 surviving genome in the FINEST rank that still has a survivor."""
        if src in self._clademax_cache:
            return self._clademax_cache[src]
        cands = []
        for rank in RANK_ORDER:
            if rank == "species":
                continue
            val = (panel.get(src) or {}).get(rank)
            if not val:
                continue
            pool = [g for g in self.by_rank[rank].get(val, ()) if g != src]
            if pool:
                cands = pool
                break
        if len(cands) > cap:
            rng = seeded_rng("clademax", self.fold, self.density, src)
            cands = sorted(rng.sample(cands, cap))
        best, _ = BASE.clade_max_baseline(set(truth_ec), cands, self.sets["ec"])
        self._clademax_cache[src] = best
        return best


def score_cells(args, ctx_inputs, writer):
    """Write every cell row; return the light bookkeeping the controls need."""
    (plan, panel, g_ec, g_rxn, ec_rxn, ec_index, truth_meta, md5_of_cell, neighbors,
     manifests, p2_runs, pf_runs, tv_sets, universe) = ctx_inputs

    folds = plan.get("folds", [])
    genome_group = plan.get("genome_group", {})
    rank_blocked = args.rank_blocked

    # --- index the runs by (fold, density) ---------------------------------- #
    p2_by = collections.defaultdict(list)
    for r in p2_runs:
        p2_by[(r["fold"], dkey(r["density"]), r["routing"])].append(r)
    pf_by = collections.defaultdict(list)
    for r in pf_runs:
        pf_by[(r["fold"], dkey(r["density"]), r["forced"])].append(r)

    densities = sorted({dkey(r["density"]) for r in p2_runs} |
                       {dkey(r["density"]) for r in pf_runs}) or [dkey(PRIMARY_DENSITY)]
    bookkeeping = {"per_fold": {}, "min_align": collections.Counter(),
                   "min_align_by_region": collections.Counter(),
                   "n_rows": 0, "densities": densities,
                   "select_modes": set(), "arms": set(),
                   "folds_scored": set(), "folds_unscored": set(),
                   # (fold, density, arm) -> md5s an arm's run never carried (C14)
                   "missing_record": {}, "missing_record_n": collections.Counter(),
                   "mixed_domain_md5": collections.Counter()}

    spaces = [s for s in args.space.split(",") if s]
    for fold_rec in folds:
        fold = int(fold_rec["fold"])
        cells = sorted(fold_rec.get("test_cells") or ())
        if args.limit_cells:
            cells = cells[:args.limit_cells]
        for dk in densities:
            density = float(dk)
            ref_ids, ref_names, ref_by_dom = reference_set(fold_rec, density,
                                                           manifests, args.mode)
            if not ref_ids:
                continue
            ctx = FoldContext(fold_rec, density, ref_ids, ref_names, g_ec, g_rxn,
                              panel, rank_blocked, ref_ids_by_domain=ref_by_dom)
            # --- PICRUSt2: merge the domain references of this (fold, density) --- #
            p2_pred, p2_dropped, p2_meta_runs = {}, set(), []
            # Which domains does a PICRUSt2 run actually exist for?  Under forced
            # routing a cell is only ever a query of its own truth-domain reference,
            # so a fold whose archaeal reference was never run must NOT have its
            # archaeal cells scored as PICRUSt2 zeros -- that would charge PICRUSt2
            # for a run the operator did not launch.  Those rows are omitted and the
            # omission is reported in meta.notes.picrust2_domains_not_run.
            # Under FORCED routing the merge is keyed on (md5, domain), not md5: an
            # amplicon md5 carried by both a bacterial and an archaeal genome is a
            # query of BOTH domain references, and the min-NSTI merge would hand such
            # a cell the other domain's prediction -- silently making the primary
            # forced arm the nsti-routed arm for those cells, which D9 makes a
            # SECONDARY arm precisely because the two trees' branch lengths are not
            # commensurable.  The min-NSTI merge is kept for routing == "nsti", which
            # is what that arm is defined to be.
            p2_cov, p2_mixed = {}, {}
            for routing in ("forced", "nsti"):
                merged = {}
                doms, domain_unknown = set(), False
                keyed_by_domain = (routing == "forced")
                for r in p2_by.get((fold, dk, routing), ()):
                    p2_meta_runs.append(r)
                    d3 = dom3(r.get("domain"))
                    if d3:
                        doms.add(d3)
                    else:
                        # "both" (an nsti run) or an unlabelled run covers everything
                        domain_unknown = True
                    if keyed_by_domain and not d3:
                        keyed_by_domain = False
                    for md5, rec in r["pred"].items():
                        key = (str(md5), d3) if (routing == "forced" and d3) else str(md5)
                        cur = merged.get(key)
                        if cur is None or ((rec.get("nsti") is not None) and
                                           (cur.get("nsti") is None or
                                            rec["nsti"] < cur["nsti"])):
                            merged[key] = rec
                    for md5 in (r["meta"].get("min_align_dropped_md5") or ()):
                        p2_dropped.add(str(md5))
                if routing == "forced" and not keyed_by_domain:
                    # an unlabelled forced run cannot be domain-keyed; fall back to
                    # the md5 merge and say so rather than dropping its predictions
                    merged = {(k[0] if isinstance(k, tuple) else k): v
                              for k, v in merged.items()}
                if merged:
                    p2_pred[routing] = merged
                    # an unrecorded domain means "assume it covers everything", so an
                    # older run without the field is never silently omitted
                    p2_cov[routing] = None if (domain_unknown or not doms) else doms
                    if routing == "forced":
                        md5s = [k[0] for k in merged if isinstance(k, tuple)]
                        p2_mixed[routing] = len(md5s) - len(set(md5s))
                        bookkeeping["mixed_domain_md5"][(fold, dk)] = \
                            p2_mixed[routing]
            # --- prFBA: shipped + forced-answer arms ---------------------------- #
            pf_pred = {}
            for forced in (False, True):
                for r in pf_by.get((fold, dk, forced), ()):
                    for mode, per in r["modes"].items():
                        # the forced runner already names its keys "<mode>_forced"
                        # (SPEC 6.2); never append the suffix twice.
                        arm = "prfba_%s" % mode
                        if forced and not mode.endswith("_forced"):
                            arm += "_forced"
                        pf_pred.setdefault(arm, {}).update(per)
                        bookkeeping["select_modes"].add(
                            mode[:-len("_forced")] if mode.endswith("_forced") else mode)
            shipped_modes = sorted(m for m in bookkeeping["select_modes"]
                                   if ("prfba_%s" % m) in pf_pred)
            primary_mode = (args.primary_select_mode if args.primary_select_mode in
                            shipped_modes else (shipped_modes[0] if shipped_modes
                                                else None))
            # Coverage of the shipped primary arm -> base_random_abstain's budget.
            # Denominator is the SCOREABLE cells only (empty-truth and min_align
            # drops are excluded from every arm's denominator), so the masked arm
            # lands at the same coverage the shipped arm is credited with.
            # Which domains has a PICRUSt2 forced run?  A cell of a domain whose run
            # was never launched is unscoreable for EVERY arm (see below), so it is
            # out of base_random_abstain's budget too.
            def _method_not_run(cell, _cov=p2_cov, _panel=panel, _tm=truth_meta):
                cov = _cov.get("forced")
                if cov is None:
                    return ""
                t = _tm.get(cell) or {}
                d = t.get("domain") or (_panel.get(
                    t.get("src_genome_id") or cell.rsplit("__", 1)[0]) or {}
                    ).get("domain") or ""
                # EXACTLY the dom_key rule used below, so a cell with no resolvable
                # domain is treated as bacterial here too (never as a third domain
                # that no run covers, which would make it unscoreable everywhere)
                k = "arc" if str(d).lower().startswith("arch") else "bac"
                return "" if k in cov else k

            shipped_cov, rand_keep = None, None
            if primary_mode:
                recs = pf_pred.get("prfba_%s" % primary_mode, {})
                elig = [c for c in cells
                        if g_ec.get((truth_meta.get(c) or {}).get("src_genome_id") or
                                    c.rsplit("__", 1)[0])
                        and md5_of_cell.get(c) not in p2_dropped
                        and not _method_not_run(c)]
                answered = sum(1 for c in elig
                               if (recs.get(md5_of_cell.get(c)) or {}).get("ec")
                               and not (recs.get(md5_of_cell.get(c)) or {}).get(
                                   "abstained", True))
                shipped_cov = (answered / len(elig)) if elig else 0.0
                # exactly, not in expectation: the promoted audit helper picks which
                # cells an uninformative abstention policy would answer, so
                # base_random_abstain lands on the shipped arm's coverage to the cell.
                rand_keep = RC.random_abstain_mask(
                    elig, shipped_cov, FROZEN["seed"] + fold * 10 + int(density * 100))

            # A fold with no method run at this density is NOT scored at all: the
            # baselines need no run and would otherwise be aggregated over folds the
            # methods never saw, which silently puts the skill score's floor and
            # ceiling on a different population from its numerator.
            if not p2_pred and not pf_pred:
                bookkeeping["folds_unscored"].add((fold, dk))
                continue
            bookkeeping["folds_scored"].add((fold, dk))
            arms = []
            for arm in sorted(pf_pred):
                arms.append(arm)
            if "forced" in p2_pred:
                arms.append("picrust2")
            if "nsti" in p2_pred:
                arms.append("picrust2_nsti_routed")
            arms += [a for a in BASELINE_ARMS
                     if a != "base_random_abstain" or
                     ("prfba_%s_forced" % primary_mode) in pf_pred]
            bookkeeping["arms"].update(arms)

            fb = bookkeeping["per_fold"].setdefault((fold, dk), collections.Counter())
            fb["n_reference_genomes"] = ctx.n_ref
            for cell in cells:
                t = truth_meta.get(cell)
                if t is None:
                    continue
                src = t.get("src_genome_id") or cell.rsplit("__", 1)[0]
                region = t.get("region") or cell.rsplit("__", 1)[1]
                domain = t.get("domain") or (panel.get(src) or {}).get("domain") or ""
                md5 = md5_of_cell.get(cell, "")
                group_id = genome_group.get(src, "")
                dropped = md5 in p2_dropped
                truth = {"ec": frozenset(g_ec.get(src, frozenset())),
                         "rxn": frozenset(g_rxn.get(src, frozenset()))}
                # neutral novelty covariate: prFBA's rank-1 hit against the SHARED,
                # fold-restricted DB; the DB is identical for both methods, so this is
                # a property of the data (SPEC 8.2).  Fall back to neighbors.json.
                rank1_gid = rank1_id = rank1_rank = None
                if primary_mode:
                    r1 = (pf_pred.get("prfba_%s" % primary_mode, {}).get(md5) or
                          pf_pred.get("prfba_%s_forced" % primary_mode, {}).get(md5))
                    if r1:
                        rank1_gid = r1.get("rank1_genome")
                        rank1_id = r1.get("rank1_identity")
                        rank1_rank = r1.get("rank1_rank")
                if rank1_gid is None and neighbors:
                    for gid, ident in (neighbors.get("rank1", {}).get(md5) or ()):
                        if gid in ctx.ref_ids and gid != src:
                            rank1_gid, rank1_id = gid, ident
                            break
                if rank1_rank is None and rank1_gid:
                    rank1_rank = first_shared_rank(panel.get(src), panel.get(rank1_gid))
                novelty_bin = SG.ident_bin(rank1_id)

                fb["n_cells"] += 1
                if dropped:
                    bookkeeping["min_align"][dk] += 1
                    bookkeeping["min_align_by_region"][region] += 1
                if not truth["ec"]:
                    fb["n_truth_empty"] += 1

                # --- the arm predictions for this cell -------------------------- #
                dom_key = "arc" if str(domain).lower().startswith("arch") else "bac"
                # A domain whose PICRUSt2 run was never launched used to drop that
                # domain's cells FOR PICRUSt2 ONLY, which left every reported F1,
                # skill score and risk-coverage curve standing on a different
                # population from prFBA's -- the printed F1s then did not subtract to
                # the printed delta.  Such a cell is now unscoreable for EVERY arm
                # (outcome unscoreable_method_not_run), so all denominators stay
                # equal and the omission is visible in cells.csv.gz.  This is NOT the
                # same as a launched run missing an md5, which is a run defect (C14).
                not_run = _method_not_run(cell)
                cell_arms = list(arms)
                preds = {}
                if not_run:
                    fb["n_method_not_run"] += 1
                else:
                    for arm in cell_arms:
                        preds[arm] = arm_prediction(
                            arm, md5, src, truth["ec"], ctx, panel, pf_pred, p2_pred,
                            rank1_gid, shipped_cov, cell, args, rand_keep,
                            dom_key=dom_key)
                        reason = preds[arm].get("reason") or ""
                        # A cell whose md5 PICRUSt2 legitimately dropped at
                        # --min_align 0.80 has no record BY DESIGN (SPEC 7-C7: the
                        # md5 is excluded from BOTH methods' denominators, and leg 2
                        # of C14 subtracts the drops from the expected set).  Such a
                        # row still carries the informative abstain_reason sentinel,
                        # but it is NOT a missing-record run defect: counting it
                        # fired C14 leg 3 on a valid run and voided the headline.
                        if reason in MISSING_RECORD_REASONS and not dropped:
                            fb["n_missing_record"] += 1
                            miss = bookkeeping["missing_record"].setdefault(
                                (fold, dk, arm), [])
                            if len(miss) < 10:
                                miss.append(md5)
                            bookkeeping["missing_record_n"][(fold, dk, arm)] += 1
                    if primary_mode and preds.get("prfba_%s" % primary_mode, {}).get(
                            "answered") and preds.get("picrust2", {}).get("answered"):
                        fb["n_matched"] += 1
                    if primary_mode and not preds.get("prfba_%s" % primary_mode,
                                                      {}).get("answered"):
                        fb["n_prfba_abstained"] += 1
                    if preds.get("picrust2", {}).get("answered"):
                        fb["n_picrust2_emitted"] += 1

                for space in spaces:
                    t_set = truth[space]
                    acc_vocab = ctx.accessory[space]
                    present = ctx.present[space]
                    n_imposs = len(t_set - present)
                    # C10's second ceiling: against the cell's OWN domain reference,
                    # which is all a forced-routed PICRUSt2 query can reach
                    own = (ctx.present_by_domain.get(dom_key) or {}).get(space)
                    n_imposs_own = len(t_set - own) if own is not None else n_imposs
                    if space == "ec":
                        fb["impossible_sum"] += n_imposs
                        fb["truth_sum"] += len(t_set)
                    uni = universe[space]
                    for arm in cell_arms:
                        info = preds.get(arm)
                        p_set = frozenset() if info is None else \
                            info["ec" if space == "ec" else "rxn"]
                        row = build_row(cell, md5, region, domain, src, fold, group_id,
                                        dk, arm, space, t_set, p_set, info, dropped,
                                        acc_vocab, n_imposs, uni, rank1_id, rank1_rank,
                                        novelty_bin, tv_sets, panel.get(src),
                                        n_imposs_own=n_imposs_own,
                                        not_run=not_run)
                        writer.writerow(row)
                        bookkeeping["n_rows"] += 1
    bookkeeping["select_modes"] = sorted(bookkeeping["select_modes"])
    bookkeeping["arms"] = sorted(bookkeeping["arms"])
    for k in ("folds_scored", "folds_unscored"):
        bookkeeping[k] = sorted("fold_%s_d%s" % fd for fd in bookkeeping[k])
    if bookkeeping["folds_unscored"]:
        print("[score_folds] %d (fold, density) cell(s) skipped: no method run "
              "(first few: %s)" % (len(bookkeeping["folds_unscored"]),
                                   bookkeeping["folds_unscored"][:5]))
    if bookkeeping["missing_record_n"]:
        print("[score_folds] *** INCOMPLETE PREDICTION TABLE(S): an md5 the fold "
              "plan queried is absent from a LAUNCHED run, which would be scored as "
              "a method abstention. C14 fails. ***")
        for (fold, dk, arm), n in sorted(bookkeeping["missing_record_n"].items()):
            print("    fold %s d%s arm %-28s %d cell(s) with no record (e.g. %s)"
                  % (fold, dk, arm, n,
                     bookkeeping["missing_record"][(fold, dk, arm)][:3]))
    if any(bookkeeping["mixed_domain_md5"].values()):
        print("[score_folds] mixed-domain md5s (queried by both domain references, "
              "routed by truth domain under --routing forced): %s"
              % {"fold_%s_d%s" % k: v
                 for k, v in sorted(bookkeeping["mixed_domain_md5"].items()) if v})
    return bookkeeping


def first_shared_rank(a, b):
    """Finest rank at which two panel records agree -- the method-neutral novelty label."""
    if not a or not b:
        return None
    for rank in RANK_ORDER:
        va, vb = a.get(rank), b.get(rank)
        if va and vb and va == vb:
            return rank
    return "none"


def arm_prediction(arm, md5, src, truth_ec, ctx, panel, pf_pred, p2_pred, rank1_gid,
                   shipped_cov, cell, args, rand_keep=None, dom_key=None):
    """-> {ec, rxn, answered, abstained, n_reps, tier, nsti, p2_domain, routing,
    best_identity, reason}."""
    base = {"ec": frozenset(), "rxn": frozenset(), "answered": False,
            "abstained": False, "n_reps": None, "tier": "", "nsti": None,
            "p2_domain": "", "routing": "", "best_identity": None, "reason": ""}
    ec_rxn = CTX_EC_RXN[0]

    def finish(ec_set, **kw):
        out = dict(base)
        out.update(kw)
        out["ec"] = frozenset(ec_set)
        out["rxn"] = project(out["ec"], ec_rxn)
        out["answered"] = bool(out["ec"]) and not out["abstained"]
        return out

    if arm.startswith("prfba_"):
        rec = pf_pred.get(arm, {}).get(md5)
        if rec is None:
            return finish(frozenset(), abstained=True, reason="no_record")
        return finish(rec["ec"], abstained=bool(rec["abstained"]) or not rec["ec"],
                      n_reps=rec["n_reps"], tier=rec["tier"] or "",
                      best_identity=rec["best_identity"],
                      reason=rec["abstain_reason"] or "")
    if arm in ("picrust2", "picrust2_nsti_routed"):
        routing = "forced" if arm == "picrust2" else "nsti"
        table = p2_pred.get(routing, {})
        # forced routing is keyed on (md5, truth domain) so a mixed-domain md5 gets
        # ITS OWN domain reference's prediction, not the lower-NSTI one (D9); an
        # unlabelled legacy run is keyed on md5 alone and falls through
        rec = table.get((md5, dom_key)) if (routing == "forced" and dom_key) else None
        if rec is None:
            rec = table.get(md5)
        if rec is None:
            return finish(frozenset(), reason="no_prediction", routing=routing)
        return finish(rec["ec"], nsti=rec.get("nsti"),
                      p2_domain=rec.get("best_domain") or "", routing=routing)
    if arm == "base_constant":
        return finish(ctx.constant["ec"])
    if arm == "base_nn1":
        if not rank1_gid:
            return finish(frozenset(), reason="no_surviving_neighbour")
        return finish(ctx.sets["ec"].get(rank1_gid, frozenset()), n_reps=1)
    if arm == "base_family_donor":
        donor = ctx.family_donor(src, panel)
        if not donor:
            return finish(frozenset(), reason="no_%s_donor" % ctx.donor_rank)
        return finish(ctx.sets["ec"].get(donor, frozenset()), n_reps=1)
    if arm == "base_clade_max":
        best = ctx.clade_max(src, panel, truth_ec, args.clade_max_cap)
        if not best:
            return finish(frozenset(), reason="no_surviving_clade")
        return finish(ctx.sets["ec"].get(best, frozenset()), n_reps=1)
    if arm == "base_random_abstain":
        forced_arm = "prfba_%s_forced" % args.primary_select_mode
        table = pf_pred.get(forced_arm)
        if table is None:
            forced = sorted(a for a in pf_pred if a.endswith("_forced"))
            table = pf_pred.get(forced[0], {}) if forced else {}
        rec = table.get(md5)
        if rec is None:
            return finish(frozenset(), abstained=True, reason="no_forced_record")
        if rand_keep is not None and cell not in rand_keep:
            return finish(frozenset(), abstained=True, reason="random_mask")
        return finish(rec["ec"], n_reps=rec["n_reps"], tier=rec["tier"] or "")
    return finish(frozenset(), reason="unknown_arm")


# EC -> ModelSEED reaction map, set once in main(); a module-level handle keeps
# arm_prediction a plain function of the cell rather than a closure over the loader.
CTX_EC_RXN = [{}]


def build_row(cell, md5, region, domain, src, fold, group_id, dk, arm, space, t_set,
              p_set, info, dropped, acc_vocab, n_imposs, universe, rank1_id, rank1_rank,
              novelty_bin, tv_sets, panel_rec, n_imposs_own=None, not_run=""):
    info = info or {"ec": frozenset(), "rxn": frozenset(), "answered": False,
                    "abstained": False, "n_reps": None, "tier": "", "nsti": None,
                    "p2_domain": "", "routing": "", "best_identity": None,
                    "reason": ""}
    tv = None
    if tv_sets and info["ec"]:
        key = (panel_rec or {}).get("template")
        if key in tv_sets:
            tv = len(info["ec"] & tv_sets[key]) / len(info["ec"])
    row = {k: "" for k in CELL_FIELDS}
    row.update({"amplicon": cell, "md5": md5, "region": region, "domain": domain,
                "src_genome_id": src, "fold": fold, "group_id": group_id,
                "density": dk, "arm": arm, "space": space,
                "n_truth": len(t_set), "n_pred": len(p_set),
                "n_reps": "" if info["n_reps"] is None else info["n_reps"],
                "tier": info["tier"] or "",
                "abstained": int(bool(info["abstained"])),
                "nsti": fmt(info["nsti"], 6),
                "nsti_bin": CMP.bin_of(info["nsti"], CMP.NSTI_BIN_EDGES,
                                       CMP.NSTI_BIN_LABELS) if info["nsti"] is not None
                else "",
                "picrust2_domain": info["p2_domain"], "routing": info["routing"],
                "rank1_identity": fmt(rank1_id, 6), "rank1_rank": rank1_rank or "",
                "ident_bin": (SG.ident_bin(info["best_identity"])
                              if info["best_identity"] is not None else ""),
                "novelty_bin": novelty_bin,
                "n_truth_impossible": n_imposs,
                "n_truth_impossible_own_domain": (n_imposs if n_imposs_own is None
                                                  else n_imposs_own),
                "template_visible_frac": fmt(tv),
                "min_align_dropped": int(bool(dropped)),
                "unscoreable_reason": "",
                # verbatim runner reason, or a MISSING_RECORD_REASONS sentinel; C14
                # counts the sentinels, and a reader can tell a tier abstention from
                # a truncated prediction table without rerunning anything
                "abstain_reason": info.get("reason") or ""})
    if dropped:
        row["outcome"] = "unscoreable_min_align"
        row["unscoreable_reason"] = "min_align_0.80"
        return row
    if not_run:
        row["outcome"] = "unscoreable_method_not_run"
        row["unscoreable_reason"] = "picrust2_domain_%s_not_run" % not_run
        return row
    if not t_set:
        row["outcome"] = "unscoreable_empty_truth"
        row["unscoreable_reason"] = "truth_empty_in_universe"
        return row
    t_acc = t_set & acc_vocab
    row["n_truth_accessory"] = len(t_acc)
    if not p_set:
        row["outcome"] = "abstained" if info["abstained"] else "pred_empty"
        return row
    row["outcome"] = "scored"
    tp, r, p, f1, j = metrics(t_set, p_set)
    row.update({"tp": tp, "recall": fmt(r), "precision": fmt(p), "f1": fmt(f1),
                "jaccard": fmt(j)})
    row["mcc"] = fmt(mcc_from_counts(tp, len(p_set) - tp, len(t_set) - tp,
                                     len(universe)))
    if t_acc:
        row["f1_accessory"] = fmt(BASE.f1_set(t_acc, p_set & acc_vocab))
    return row


# --------------------------------------------------------------------------- #
# Stage 2 -- aggregate the cell table
# --------------------------------------------------------------------------- #
def iter_cells(path):
    with gzip.open(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh):
            yield row


def f1_pop(row):
    """Population F1: a blank (abstained / empty prediction) counts as 0.0."""
    return float(row["f1"]) if row["f1"] else 0.0


def scoreable(row):
    return row["outcome"] in ("scored", "abstained", "pred_empty")


NESTED_ARM = "prfba_nestedcv"


def split_prfba_arm(arm):
    """'prfba_legacy_forced' -> ('legacy', True); a non-prFBA arm -> (None, False)."""
    if not arm.startswith("prfba_"):
        return None, False
    rest = arm[len("prfba_"):]
    if rest.endswith("_forced"):
        return rest[:-len("_forced")], True
    return rest, False


def cells_header(path):
    with gzip.open(path, "rt", newline="") as fh:
        return next(csv.reader(fh), [])


def aggregate(args, cells_path, plan, panel, manifests, p2_runs, pf_runs, run_meta,
              run_comp=None, inventory=None):
    primary_dk = dkey(PRIMARY_DENSITY)
    # A cells table written before the C14 columns existed cannot be audited for
    # completeness -- and quietly treating its missing abstain_reason as "" is exactly
    # the laundering C14 exists to stop.  Refuse it instead.
    head = cells_header(cells_path)
    need = [c for c in ("abstain_reason", "n_truth_impossible_own_domain")
            if c not in head]
    if need:
        raise SystemExit(
            "[score_folds] %s predates the %s column(s): rescore with --force. A "
            "cells table without abstain_reason cannot distinguish a tier abstention "
            "from a truncated prediction run (C14)." % (cells_path, ", ".join(need)))
    # ---- pass A: matched cells, knobs, paired vectors, fold accounting ------ #
    answered_fa = collections.defaultdict(set)     # (dk, fold, arm) -> cells answered
    seen_arms, seen_spaces, seen_dens = set(), set(), set()
    per_fold_by_arm = collections.defaultdict(collections.Counter)
    cell_meta = {}
    fold_mode_f1 = collections.defaultdict(lambda: collections.defaultdict(list))
    p2_domains = collections.defaultdict(set)      # (fold, dk) -> domains with a p2 row
    all_domains = collections.defaultdict(set)
    # C14 leg 3 and the per-space unscoreable audit, over EVERY space and arm
    missing_by_arm = collections.defaultdict(collections.Counter)   # arm -> reason/fold
    unscoreable_by_space = collections.defaultdict(collections.Counter)
    # C10: per (arm, fold, domain) truth / impossible sums, union and own-domain
    ceil_by_fd = collections.defaultdict(collections.Counter)
    nsti_by_cell = {}                              # cell -> PICRUSt2's NSTI (D7 knob)
    seen_fold_dk = set()                           # C14 leg 5: the scored grid
    for row in iter_cells(cells_path):
        arm, space, dk = row["arm"], row["space"], row["density"]
        seen_arms.add(arm)
        seen_spaces.add(space)
        seen_dens.add(dk)
        seen_fold_dk.add((str(row["fold"]), dk))
        if row["outcome"].startswith("unscoreable"):
            unscoreable_by_space[(space, dk)][
                (arm, row["outcome"])] += 1
        # C14 leg 3 counts rows scored against a MISSING record.  A min_align drop
        # is not one: that md5 is excluded from every arm's denominator by design
        # (SPEC 7-C7) and leg 2 subtracts it from the expected set, so counting it
        # here would make leg 3 contradict leg 2 and void a valid run.  The
        # exemption is deliberately narrow -- ONLY unscoreable_min_align; a missing
        # record on an empty-truth cell is still a run defect.
        if (row["abstain_reason"] in MISSING_RECORD_REASONS
                and row["outcome"] != "unscoreable_min_align"):
            missing_by_arm[arm][(space, dk, row["fold"],
                                 row["abstain_reason"])] += 1
        if space != "ec":
            continue
        if dk == primary_dk and arm == "picrust2" and row["nsti"]:
            nsti_by_cell[row["amplicon"]] = float(row["nsti"])
        if scoreable(row):
            c2 = ceil_by_fd[(arm, row["fold"], row["domain"])]
            c2["truth_sum"] += int(row["n_truth"] or 0)
            c2["impossible_sum"] += int(row["n_truth_impossible"] or 0)
            c2["impossible_own_sum"] += int(
                row["n_truth_impossible_own_domain"] or 0)
            c2["n"] += 1
        if scoreable(row) and row["outcome"] == "scored":
            answered_fa[(dk, row["fold"], arm)].add(row["amplicon"])
        all_domains[(row["fold"], dk)].add(row["domain"])
        # a domain whose PICRUSt2 run was never launched is now unscoreable for every
        # arm (equal denominators); the omission is disclosed from those rows
        if row["outcome"] == "unscoreable_method_not_run":
            p2_domains[(row["fold"], dk)].add(row["domain"])
        if dk == primary_dk:
            cell_meta[row["amplicon"]] = (row["group_id"], row["src_genome_id"],
                                          row["region"], row["domain"],
                                          row["novelty_bin"], row["fold"])
            mode, forced = split_prfba_arm(arm)
            if mode and not forced and scoreable(row):
                fold_mode_f1[row["fold"]][mode].append(f1_pop(row))
            fk = row["fold"]
            c = per_fold_by_arm[(arm, fk)]
            c["n_cells"] += 1
            if row["outcome"] == "unscoreable_empty_truth":
                c["n_truth_empty"] += 1
            if row["outcome"] == "unscoreable_min_align":
                c["n_min_align"] += 1
            if row["outcome"] == "unscoreable_method_not_run":
                c["n_method_not_run"] += 1
            if (row["abstain_reason"] in MISSING_RECORD_REASONS
                    and row["outcome"] != "unscoreable_min_align"):
                c["n_no_record"] += 1          # same min_align exemption as leg 3
            if row["outcome"] == "scored":
                c["n_emitted"] += 1
            if scoreable(row):
                c["n_scoreable"] += 1
                c["impossible_sum"] += int(row["n_truth_impossible"] or 0)
                c["impossible_own_sum"] += int(
                    row["n_truth_impossible_own_domain"] or 0)
                c["truth_sum"] += int(row["n_truth"] or 0)

    # ---- nested CV (D8): the per-fold winner is chosen on the OTHER folds --- #
    nested, winners = nested_cv(fold_mode_f1)
    # The PRIMARY prFBA arm is the nested-CV composite: fold i is scored with the
    # select mode that won on folds != i, so the headline is not the shipped default
    # fitted on these same labels.  The shipped arm survives as a labelled arm.
    composite_cells = {(f, "prfba_%s" % m) for f, m in winners.items() if m}

    # accounting is anchored on whichever arm actually covers the cells (a prFBA-only
    # run has no picrust2 rows, and a single-domain PICRUSt2 run covers only part)
    anchor = next((a for a in ("picrust2", NESTED_ARM,
                               "prfba_%s" % args.primary_select_mode)
                   if any(k[0] == a for k in per_fold_by_arm)),
                  (sorted(per_fold_by_arm)[0][0] if per_fold_by_arm else None))
    per_fold = collections.defaultdict(collections.Counter)
    # iterate a snapshot and never index the defaultdict: looking up a (picrust2, fold)
    # key that does not exist would insert it mid-iteration, which is exactly what
    # happens on a prFBA-only run or a fold PICRUSt2 never covered.
    for (arm, fk), c in sorted(per_fold_by_arm.items()):
        if arm == anchor:
            per_fold[fk].update(c)
            per_fold[fk]["n_picrust2_emitted"] = per_fold_by_arm.get(
                ("picrust2", fk), {}).get("n_emitted", 0)

    answered = collections.defaultdict(set)
    for (dk, fk, arm), cells in answered_fa.items():
        answered[(dk, arm)] |= cells
        mode, forced = split_prfba_arm(arm)
        if mode and (fk, "prfba_%s" % mode) in composite_cells:
            answered[(dk, NESTED_ARM + ("_forced" if forced else ""))] |= cells
    primary_prfba = pick_primary_arm(seen_arms | {NESTED_ARM} if composite_cells
                                     else seen_arms, args.primary_select_mode)
    matched = (answered.get((primary_dk, primary_prfba), set()) &
               answered.get((primary_dk, "picrust2"), set())) if primary_prfba else set()

    # ---- pass B: blocks, knob tables, paired vectors ------------------------ #
    blocks = collections.defaultdict(ArmAgg)        # (arm, space, kind, value) -> agg
    matched_blocks = collections.defaultdict(ArmAgg)
    knob = collections.defaultdict(dict)            # (arm, space) -> cell -> (knob, f1)
    knob_cross = collections.defaultdict(dict)      # same, on the OTHER method's knob
    paired = collections.defaultdict(dict)          # (arm, space) -> cell -> f1
    abstain_counts = collections.Counter()
    for row in iter_cells(cells_path):
        if not scoreable(row):
            continue
        arm, space, dk = row["arm"], row["space"], row["density"]
        t_n, p_n = int(row["n_truth"] or 0), int(row["n_pred"] or 0)
        tp = int(row["tp"]) if row["tp"] != "" else 0
        r = float(row["recall"]) if row["recall"] else None
        p = float(row["precision"]) if row["precision"] else None
        f1 = float(row["f1"]) if row["f1"] else None
        j = float(row["jaccard"]) if row["jaccard"] else None
        mcc = float(row["mcc"]) if row["mcc"] else None
        tv = float(row["template_visible_frac"]) if row["template_visible_frac"] else None
        args_tuple = (t_n, p_n, tp, r, p, f1, j, mcc, tv)
        mode, forced = split_prfba_arm(arm)
        arm_names = [arm]
        if mode and (row["fold"], "prfba_%s" % mode) in composite_cells:
            arm_names.append(NESTED_ARM + ("_forced" if forced else ""))
        for name in arm_names:
            if row["abstained"] == "1":
                abstain_counts[(dk, name)] += 1
            if dk == primary_dk:
                for kind, value in (("all", "all"), ("region", row["region"]),
                                    ("novelty_bin", row["novelty_bin"]),
                                    ("domain", row["domain"]), ("fold", row["fold"])):
                    blocks[(name, space, kind, value)].add_cell(*args_tuple)
                if row["amplicon"] in matched:
                    matched_blocks[(name, space)].add_cell(*args_tuple)
                paired[(name, space)][row["amplicon"]] = f1_pop(row)
                k = k_cross = None
                if name.startswith("prfba_"):
                    k = float(row["rank1_identity"]) if row["rank1_identity"] else None
                    # D7 sweeps each method on its OWN knob; the companion curve
                    # (each method on the other's knob) is what separates "better
                    # abstention policy" from "better confidence signal"
                    k_cross = nsti_by_cell.get(row["amplicon"])
                elif name.startswith("picrust2"):
                    k = float(row["nsti"]) if row["nsti"] else None
                    k_cross = (float(row["rank1_identity"])
                               if row["rank1_identity"] else None)
                knob[(name, space)][row["amplicon"]] = (k, f1_pop(row))
                knob_cross[(name, space)][row["amplicon"]] = (k_cross, f1_pop(row))
            blocks[(name, space, "density", dk)].add_cell(*args_tuple)
    arms = sorted(seen_arms | {a for (a, _s, _k, _v) in blocks})
    omitted = {"fold_%s_d%s" % (f, d): sorted(p2_domains[(f, d)])
               for (f, d) in sorted(all_domains) if p2_domains[(f, d)]}
    # ---- C14 leg 5: did EVERY (fold, density) the plan lists get scored? --- #
    # A fold absent from BOTH methods' runs is not a laundered abstention -- it is
    # simply gone: both denominators shrink by the same cells, so the leg-4
    # denominator tripwire stays silent and legs 1-3 only inspect LAUNCHED runs.
    # The headline would then be quoted over a silently reduced population (at 10
    # folds, one OOM-killed fold = 90% of the planned cells).  The planned density
    # set is the set the LAUNCHED runs cover (the frozen fractions are not all run
    # at smoke scale), so the grid is plan folds x launched densities.
    plan_folds = sorted({str(int(f["fold"])) for f in plan.get("folds", [])
                         if f.get("fold") is not None})
    # --expect-folds declares THIS invocation's fold scope (SPEC 11.1's smoke recipe
    # deliberately launches one fold of three), so leg 5 is gated on the declared
    # subset -- the same scope leg 6 derives its inventory from.  Without the flag
    # the scope is every planned fold, i.e. unchanged.
    declared_folds = sorted({str(int(x)) for x in
                             str(args.expect_folds or "").replace(" ", "").split(",")
                             if x}, key=int) if args.expect_folds else []
    if declared_folds:
        plan_folds = [f for f in plan_folds if f in set(declared_folds)]
    run_dens = sorted({dkey(r["density"]) for r in p2_runs} |
                      {dkey(r["density"]) for r in pf_runs}) or [primary_dk]
    # include_self / plumbing / recipe_control reference families are single-fold by
    # construction (the incself run covers fold 0 only), so the grid is gated for
    # the fold-family modes only and merely reported for the others.
    cov_gated = args.mode in ("fold", "density")
    scored_folds = sorted({f for f, _d in seen_fold_dk})
    grid_missing = sorted({(f, d) for f in plan_folds for d in run_dens} -
                          seen_fold_dk)
    plan_coverage = {
        "gated": cov_gated,
        "folds_in_plan": plan_folds,
        "folds_declared": declared_folds,
        "folds_scored": scored_folds,
        "folds_in_plan_not_scored": sorted(set(plan_folds) - set(scored_folds)),
        "densities_in_launched_runs": run_dens,
        "densities_scored": sorted({d for _f, d in seen_fold_dk}),
        "grid_missing": ["fold_%s_d%s" % p for p in grid_missing],
        "n_grid_missing": len(grid_missing),
        "ok": (not grid_missing) if cov_gated else None,
    }
    completeness = {
        # Stage 0: the derived-vs-discovered run inventory (C14 leg 6).  It
        # SUBSUMES the leg-by-leg checks below -- a dead, duplicated, stale or
        # under-specified run directory is caught by set equality against the plan
        # rather than by whichever symptom it happened to produce -- while leaving
        # every one of them reporting.
        "inventory": inventory or {},
        "plan_coverage": plan_coverage,
        "missing_record_rows_by_arm": {
            arm: {"n_rows": sum(c.values()),
                  "detail": {"space=%s|density=%s|fold=%s|%s" % k: v
                             for k, v in sorted(c.items())}}
            for arm, c in sorted(missing_by_arm.items())},
        "n_missing_record_rows": sum(sum(c.values())
                                     for c in missing_by_arm.values()),
        "unscoreable_by_space": {
            "space=%s|density=%s" % k: {"%s|%s" % kk: vv
                                        for kk, vv in sorted(c.items())}
            for k, c in sorted(unscoreable_by_space.items())},
        "unscoreable_uniform_across_arms": _unscoreable_uniform(unscoreable_by_space),
        # C14 legs 1-2, computed from the runs' own metas and the fold plan
        "runs": (run_comp or {}).get("per_run", {}),
        "n_runs_incomplete": (run_comp or {}).get("n_runs_incomplete", 0),
        "runs_complete": (run_comp or {}).get("ok"),
    }
    return finalize(args, blocks, matched_blocks, knob, paired, per_fold, cell_meta,
                    nested, winners, arms, sorted(seen_spaces),
                    sorted(seen_dens), primary_prfba, matched, plan, panel, manifests,
                    p2_runs, pf_runs, run_meta, abstain_counts, anchor, omitted,
                    knob_cross=knob_cross, ceil_by_fd=ceil_by_fd,
                    per_fold_by_arm=per_fold_by_arm, completeness=completeness)


def _unscoreable_uniform(unscoreable_by_space):
    """Every arm must see the SAME unscoreable cells: min_align drops, empty truth and
    an unlaunched domain are properties of the cell, never of a method.  A False here
    means some arm's denominator differs from the others' (R2F1's defect class)."""
    for _key, c in unscoreable_by_space.items():
        per_arm = collections.defaultdict(dict)
        for (arm, outcome), n in c.items():
            per_arm[arm][outcome] = n
        vals = list(per_arm.values())
        if any(v != vals[0] for v in vals[1:]):
            return False
    return True


def nested_cv(fold_mode_f1):
    """D8: per fold, the select mode with the best mean F1 on the OTHER folds.

    -> (summary block, {fold: winning mode}).  The shipped default was itself chosen
    by maximising gene-capture F1 on these very amplicons, so scoring the headline
    with it would report a number fitted on every fold's labels; the winner is
    therefore re-selected per fold using training folds only.
    """
    fold_means = {f: {m: (st.mean(v) if v else None) for m, v in d.items()}
                  for f, d in fold_mode_f1.items()}
    all_modes = sorted({m for d in fold_means.values() for m in d})
    nested, winners = {}, {}
    for f in sorted(fold_means):
        inner = {}
        for m in all_modes:
            vals = [fold_means[g][m] for g in fold_means
                    if g != f and fold_means[g].get(m) is not None]
            inner[m] = rnd(st.mean(vals)) if vals else None
        # with a single fold there is no training set: fall back to this fold's own
        # ranking and say so, rather than silently reporting a fitted winner
        if all(v is None for v in inner.values()):
            inner = {m: rnd(fold_means[f].get(m)) for m in all_modes}
            nested.setdefault("_meta", {})["single_fold_fallback"] = True
        best = max((m for m in all_modes if inner.get(m) is not None),
                   key=lambda m: inner[m], default=None)
        winners[f] = best
        nested["fold_%s" % f] = {"inner_winner": best, "inner_f1": inner}
    if winners:
        modal = collections.Counter(winners.values()).most_common(1)[0][0]
        nested["winner_instability"] = rnd(
            sum(1 for v in winners.values() if v != modal) / len(winners))
    else:
        nested["winner_instability"] = None
    return nested, winners


def pick_primary_arm(arms, preferred):
    """The prFBA arm the headline and the matched subset are defined against.

    The nested-CV composite wins when it exists (SPEC D8 makes it the primary and
    demotes the shipped select mode to a labelled sensitivity arm).
    """
    cands = sorted(a for a in arms if a.startswith("prfba_") and
                   not a.endswith("_forced"))
    if not cands:
        return None
    for want in (NESTED_ARM, "prfba_%s" % preferred):
        if want in cands:
            return want
    return cands[0]


def cluster_bootstrap(diffs, groups, n_boot, seed=FROZEN["seed"]):
    """Cluster bootstrap over holdout groups -> (delta, ci95, p, n_groups, ci_note).

    ``ci_note`` is None when the interval is a real 2000+-resample interval,
    ``bootstrap_disabled`` / ``fewer_than_2_holdout_groups`` when there is NO
    interval, and ``bootstrap_below_<MIN_BOOTSTRAP>`` when one was computed from too
    few resamples to stand behind a published claim.

    Groups -- not amplicons -- are resampled with replacement, because the per-amplicon
    difference is clustered at ICC ~0.35-0.66 and amplicon-level resampling would give
    a CI several times too narrow.
    """
    if not diffs:
        return None, [None, None], None, 0, "no_paired_cells"
    by = collections.defaultdict(list)
    for key, v in diffs.items():
        by[groups.get(key, key)].append(v)
    keys = sorted(by)
    sums = np.array([sum(by[k]) for k in keys], dtype=float)
    cnts = np.array([len(by[k]) for k in keys], dtype=float)
    delta = float(sums.sum() / cnts.sum())
    # WHY there is no interval travels with the absent interval: "no CI was
    # computed" and "the effect is below the SESOI" were published in the same words.
    if len(keys) < 2 or n_boot <= 0:
        return delta, [None, None], None, len(keys), (
            "bootstrap_disabled" if n_boot <= 0 else "fewer_than_2_holdout_groups")
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(keys), size=(n_boot, len(keys)))
    boot = sums[idx].sum(axis=1) / cnts[idx].sum(axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    p = 2.0 * min(float((boot <= 0).mean()), float((boot >= 0).mean()))
    p = max(min(p, 1.0), 1.0 / n_boot)
    return delta, [float(lo), float(hi)], p, len(keys), (
        None if n_boot >= MIN_BOOTSTRAP else "bootstrap_below_%d" % MIN_BOOTSTRAP)


def glmm_p(diffs, groups, covariate=None, cap=20000):
    """MixedLM on the paired difference with a holdout-group random intercept.

    -> ``(p_value, model_label)``.  The novelty covariate is dropped when it carries
    no variance in the subset (which makes the design matrix singular); the model
    actually fitted is reported rather than silently returning no p-value at all.
    """
    try:
        import pandas as pd
        from statsmodels.regression.mixed_linear_model import MixedLM
    except Exception:
        return None, None
    keys = sorted(diffs)
    if len(keys) < 30:
        return None, "n<30"
    if covariate is not None:
        vals = [covariate.get(k) for k in keys if covariate.get(k) is not None]
        if len(vals) < len(keys) or len(set(vals)) < 3 or st.pstdev(vals) < 1e-9:
            covariate = None
    if len(keys) > cap:
        rng = random.Random(FROZEN["seed"])
        keys = sorted(rng.sample(keys, cap))
    rows = {"y": [diffs[k] for k in keys],
            "g": [groups.get(k, k) for k in keys]}
    if covariate:
        rows["x"] = [covariate.get(k) if covariate.get(k) is not None else 0.0
                     for k in keys]
    df = pd.DataFrame(rows)
    import warnings
    attempts = ([("intercept+rank1_identity", True), ("intercept", False)]
                if covariate else [("intercept", False)])
    for label, use_cov in attempts:
        try:
            exog = np.ones((len(df), 2 if use_cov else 1))
            if use_cov:
                exog[:, 1] = df["x"].to_numpy(dtype=float)
            with warnings.catch_warnings():
                # a near-zero group variance is informative, not an error: it is
                # reported through icc_group, and the fixed effect is still wanted.
                warnings.simplefilter("ignore")
                res = MixedLM(df["y"].to_numpy(dtype=float), exog,
                              groups=df["g"]).fit()
            return float(res.pvalues[0]), label
        except Exception:
            continue
    return None, "failed"


def rc_block(arm, space, knob, arms_present, cross=False):
    """Risk-coverage for one arm, swept on its own knob (D7).

    A prFBA arm is swept on the FORCED-answer arm when one exists, because a
    risk-coverage curve needs a full-coverage run to sweep -- the shipped arm has
    already spent its abstention budget.

    With ``cross=True`` the arm is swept on the OTHER method's knob (prFBA on NSTI,
    PICRUSt2 on rank1_identity), which is the companion curve that separates "this
    method's abstention policy is better" from "this method's confidence signal is
    better".  D7 keeps the own-knob curves as the co-primary endpoint; the cross
    curves are reported beside them.
    """
    source = (arm, space)
    if arm.startswith("prfba_") and not arm.endswith("_forced"):
        forced = ("%s_forced" % arm, space)
        if forced in knob:
            source = forced
    table = knob.get(source) or knob.get((arm, space)) or {}
    if not table:
        return {"aurc": None, "at": {str(c): None for c in COVERAGE_GRID},
                "knob": None, "source_arm": None}
    cells = sorted(table)
    knobs = [table[c][0] for c in cells]
    scores = [table[c][1] for c in cells]
    own_is_identity = not arm.startswith("picrust2")
    # rank1_identity: higher = more confident; NSTI: lower = more confident
    on_identity = own_is_identity if not cross else (not own_is_identity)
    b = RC.risk_coverage(scores, knobs, higher_is_better=on_identity,
                         grid=COVERAGE_GRID, seed=FROZEN["seed"])
    return {"aurc": rnd(b["aurc"]), "at": {str(c): rnd(b["at"][c]) for c in COVERAGE_GRID},
            "knob": ("rank1_identity" if on_identity else "nsti"),
            "n_knob_missing": sum(1 for k in knobs if k is None),
            "source_arm": source[0]}


def finalize(args, blocks, matched_blocks, knob, paired, per_fold, cell_meta,
             nested, winners, arms, spaces, densities, primary_prfba, matched, plan,
             panel, manifests, p2_runs, pf_runs, run_meta, abstain_counts, anchor,
             omitted_picrust2_domains, knob_cross=None, ceil_by_fd=None,
             per_fold_by_arm=None, completeness=None):
    groups = {c: m[0] for c, m in cell_meta.items()}
    genomes = {c: m[1] for c, m in cell_meta.items()}
    rank_blocked = args.rank_blocked
    primary_dk = dkey(PRIMARY_DENSITY)

    def pop_f1(arm, space, kind="all", value="all"):
        b = blocks.get((arm, space, kind, value))
        if not b:
            return None
        blk = b.population_block()["f1"]
        return blk["value"] if blk else None

    # ---- tips per density, for the density slope --------------------------- #
    # A reference at one density is split across the two domain trees, so the tip
    # count of a (fold, density) reference is the SUM over its domain manifests; the
    # density axis is then the mean over folds.  Control references (recipe, include
    # self, plumbing) are excluded -- they are not points on this curve.
    tips_fd = collections.defaultdict(int)
    refs_fd = {}
    mode_prefix = {"include_self": "incself", "plumbing": "plumbing",
                   "recipe_control": "recipe"}.get(args.mode)
    for name, man in sorted(manifests.items()):
        if mode_prefix:
            if not name.startswith(mode_prefix):
                continue
        elif (man.get("mode") or "fold") != "fold" or man.get("fold") is None:
            continue
        key = (int(man.get("fold") or 0),
               dkey(man.get("density") or PRIMARY_DENSITY))
        tips_fd[key] += int(man.get("n_tips") or 0)
        refs_fd.setdefault(key, set()).update(man.get("reference_genomes") or ())
    tips_mean = {}
    for d in densities:
        vals = [v for (f, dd), v in tips_fd.items() if dd == d and v]
        tips_mean[d] = st.mean(vals) if vals else None
    primary_refs = {f: len(v) for (f, dd), v in refs_fd.items() if dd == primary_dk}
    ref_counts = sorted(primary_refs.values())
    if not ref_counts:
        ref_counts = [len(f.get("reference_genomes") or ())
                      for f in plan.get("folds", [])] or [0]

    # ---- arms block -------------------------------------------------------- #
    out_arms = {}
    for arm in arms:
        for space in spaces:
            b = blocks.get((arm, space, "all", "all"))
            if b is None:
                continue
            base_const = pop_f1("base_constant", space)
            base_ceil = pop_f1("base_clade_max", space)
            f1v = pop_f1(arm, space)
            skill = None
            if None not in (f1v, base_const, base_ceil) and base_ceil != base_const:
                skill = (f1v - base_const) / (base_ceil - base_const)
            rc = rc_block(arm, space, knob, arms)
            rc_x = rc_block(arm, space, knob_cross or {}, arms, cross=True)
            by_density = {}
            for d in densities:
                bd = blocks.get((arm, space, "density", d))
                if bd is not None:
                    blk = bd.population_block()["f1"]
                    by_density[d] = {"f1": blk["value"] if blk else None,
                                     "n": bd.n, "n_tips": tips_mean.get(d)}
            xs, ys = [], []
            for d, v in sorted(by_density.items()):
                if v["f1"] is not None and v.get("n_tips"):
                    xs.append(math.log10(max(v["n_tips"], 1)))
                    ys.append(v["f1"])
            mb = matched_blocks.get((arm, space))
            out_arms.setdefault(arm, {})[space] = {
                "population": b.population_block(),
                "conditional_matched": (mb.conditional_block() if mb is not None else
                                        {"n_matched": 0}),
                "micro": b.micro_block(),
                "coverage": b.coverage(),
                "skill_score": rnd(skill),
                "aurc": rc["aurc"],
                "f1_at_coverage": rc["at"],
                "by_region": strat_block(blocks, arm, space, "region"),
                "by_novelty_bin": strat_block(blocks, arm, space, "novelty_bin"),
                "by_domain": strat_block(blocks, arm, space, "domain"),
                "by_density": by_density,
                "density_slope_f1_per_log10_tips": rnd(lsq_slope(xs, ys)),
                "template_visible_frac_mean": b.tv_mean(),
                "n_abstained": abstain_counts.get((primary_dk, arm), 0),
                "risk_coverage_knob": rc["knob"],
                "risk_coverage_source_arm": rc["source_arm"],
                # companion curve on the other method's knob (reported, not the D7
                # endpoint): prFBA on NSTI, PICRUSt2 on rank1_identity
                "aurc_cross_knob": rc_x["aurc"],
                "f1_at_coverage_cross_knob": rc_x["at"],
                "risk_coverage_knob_cross": rc_x["knob"],
                "n_cross_knob_missing": rc_x.get("n_knob_missing"),
            }

    # ---- contrasts --------------------------------------------------------- #
    pairs = []
    for arm in arms:
        if arm.startswith("picrust2"):
            continue
        pairs.append((arm, "picrust2"))
    if primary_prfba:
        for b in BASELINE_ARMS:
            if b in arms:
                pairs.append((primary_prfba, b))
    if "base_nn1" in arms:
        pairs.append(("picrust2", "base_nn1"))
    contrasts = {}
    secondary_p = []
    for a, b in pairs:
        if a == b:
            continue
        for space in spaces:
            da, db = paired.get((a, space), {}), paired.get((b, space), {})
            common = sorted(set(da) & set(db))
            if not common:
                continue
            diffs = {c: da[c] - db[c] for c in common}
            delta, ci, p, n_groups, ci_note = cluster_bootstrap(
                diffs, groups, args.bootstrap)
            icc_g = ICC.icc(diffs, lambda c: genomes.get(c))
            icc_grp = ICC.icc(diffs, lambda c: groups.get(c))
            cov = {c: (knob.get((a, space), {}).get(c) or (None,))[0] for c in common}
            rc_a = out_arms.get(a, {}).get(space, {}).get("f1_at_coverage", {})
            rc_b = out_arms.get(b, {}).get(space, {}).get("f1_at_coverage", {})
            key = "%s__vs__%s" % (a, b)
            p_glmm, glmm_model = glmm_p(diffs, groups, cov)
            blk = {
                "delta_f1_population": rnd(delta),
                "ci95": [rnd(ci[0]), rnd(ci[1])],
                "bootstrap_n_groups": n_groups,
                "icc_genome": rnd(icc_g["icc"]),
                "icc_group": rnd(icc_grp["icc"]),
                "n_eff": int(round(icc_grp["n_eff"])) if icc_grp["n_eff"] else 0,
                "p_cluster_bootstrap": rnd(p, 6),
                "p_glmm": rnd(p_glmm, 6),
                "p_glmm_model": glmm_model,
                "q_bh": None,
                # THREE outcomes, never two: exceeded / not exceeded / INDETERMINATE.
                # exceeds_sesoi alone reported "no interval was computed" and "the
                # effect is below the smallest effect size of interest" with the
                # same False, and the headline sentence said "SESOI not exceeded"
                # for both -- so --bootstrap 0 inverted the campaign's primary claim
                # at EXIT 0.  exceeds_sesoi keeps its meaning (a decided, exceeded
                # claim) for consumers that already read it.
                "exceeds_sesoi": bool(
                    delta is not None and abs(delta) >= FROZEN["sesoi_delta_f1"] and
                    ci[0] is not None and (ci[0] > 0 or ci[1] < 0)),
                "ci_available": ci[0] is not None,
                "ci_unavailable_reason": (ci_note if ci[0] is None else None),
                "ci_underpowered": bool(ci[0] is not None and ci_note),
                "ci_underpowered_reason": (ci_note if ci[0] is not None else None),
                "min_bootstrap": MIN_BOOTSTRAP,
                "sesoi_verdict": (
                    "indeterminate_no_ci" if (delta is None or ci[0] is None) else
                    ("exceeded" if (abs(delta) >= FROZEN["sesoi_delta_f1"] and
                                    (ci[0] > 0 or ci[1] < 0)) else "not_exceeded")),
                "delta_f1_at_matched_coverage": {
                    str(c): rnd((rc_a.get(str(c)) - rc_b.get(str(c)))
                                if (rc_a.get(str(c)) is not None and
                                    rc_b.get(str(c)) is not None) else None)
                    for c in (0.5, 0.75, 0.9)},
                "n_pairs": len(common),
                "by_region": {}, "by_domain": {},
            }
            for kind, dest in (("region", "by_region"), ("domain", "by_domain")):
                sub = collections.defaultdict(dict)
                for c in common:
                    v = cell_meta.get(c)
                    if not v:
                        continue
                    label = v[2] if kind == "region" else v[3]
                    sub[label][c] = diffs[c]
                for label in sorted(sub):
                    d2, _ci2, p2, _ng, _note2 = cluster_bootstrap(
                        sub[label], groups, min(args.bootstrap, 500))
                    blk[dest][label] = {"delta_f1_population": rnd(d2), "q_bh": None,
                                        "p_cluster_bootstrap": rnd(p2, 6),
                                        "n_pairs": len(sub[label])}
                    secondary_p.append((key, space, dest, label, p2))
            contrasts.setdefault(key, {})[space] = blk
            is_primary = (a == primary_prfba and b == "picrust2" and space == "ec")
            if not is_primary:
                secondary_p.append((key, space, None, None, p))
    qs = bh_fdr([t[4] for t in secondary_p])
    for (key, space, dest, label, _p), q in zip(secondary_p, qs):
        blk = contrasts.get(key, {}).get(space)
        if blk is None:
            continue
        if dest is None:
            blk["q_bh"] = rnd(q, 6)
        else:
            blk[dest][label]["q_bh"] = rnd(q, 6)

    # ---- accounting -------------------------------------------------------- #
    ceil_by_fd = ceil_by_fd or {}
    per_fold_by_arm = per_fold_by_arm or {}
    acc_rows = []
    for f in sorted(per_fold, key=lambda x: (len(x), x)):
        c = per_fold[f]
        matched_fold = sum(1 for cell, m in cell_meta.items()
                           if m[5] == f and cell in matched)
        truth_sum = c.get("truth_sum", 0)
        # C10 per domain: the union ceiling is prFBA-reachable, the own-domain
        # ceiling is what a forced-routed PICRUSt2 query can reach (its domain tree
        # and trait table alone).  The two differ by ~3-5 recall points in the thin
        # archaeal stratum, which SPEC 8.4 requires be reported separately.
        by_dom = {}
        for (a, fk, dom), cc in sorted(ceil_by_fd.items()):
            if a != anchor or fk != f or not dom:
                continue
            ts = cc.get("truth_sum", 0)
            by_dom[dom] = {
                "n_cells_scoreable": cc.get("n", 0),
                "truth_sum": ts,
                "structural_recall_ceiling_union": rnd(
                    1 - cc.get("impossible_sum", 0) / ts) if ts else None,
                "structural_recall_ceiling_own_domain": rnd(
                    1 - cc.get("impossible_own_sum", 0) / ts) if ts else None}
        acc_rows.append({
            "fold": int(f) if str(f).isdigit() else f,
            "n_cells": c.get("n_cells", 0),
            "n_prfba_abstained": fold_abstained(blocks, primary_prfba, f),
            "n_picrust2_emitted": c.get("n_picrust2_emitted", 0),
            "n_truth_empty": c.get("n_truth_empty", 0),
            "n_matched": matched_fold,
            "n_impossible_ec_mean": rnd(c.get("impossible_sum", 0) /
                                        max(1, c.get("n_scoreable", 0))),
            "structural_recall_ceiling": rnd(
                1 - c.get("impossible_sum", 0) / truth_sum) if truth_sum else None,
            "structural_recall_ceiling_picrust2_own_tree": rnd(
                1 - c.get("impossible_own_sum", 0) / truth_sum)
            if truth_sum else None,
            "structural_recall_ceiling_by_domain": by_dom,
            "n_min_align_dropped": c.get("n_min_align", 0),
            "n_method_not_run": c.get("n_method_not_run", 0),
            # C14: an md5 the fold plan queried that a LAUNCHED run never carried.
            # Non-zero means some arm's cells were scored against a truncated or
            # stale prediction table, which reads out as lower coverage.
            "n_no_record": {a: int(cc.get("n_no_record", 0))
                            for (a, fk), cc in sorted(per_fold_by_arm.items())
                            if fk == f and cc.get("n_no_record")},
            "n_reference_genomes": primary_refs.get(
                int(f) if str(f).lstrip("-").isdigit() else f,
                next((len(x.get("reference_genomes") or ())
                      for x in plan.get("folds", []) if str(x["fold"]) == str(f)), 0)),
        })
    # the trait-column screen that build_reference applied, per reference: the 'full'
    # (prevalence-n) columns are reinstated as PICRUSt2 predictions, the 'zero' ones
    # are genuinely unpredictable, and 'singleton'/'crash' (the SPEC 5.4 aggressive
    # retry) are NOT repairable -- a fold that used the retry carries a disclosed,
    # unquantified PICRUSt2 recall penalty, which C15 reports.
    trait_screen = {}
    for name, man in sorted(manifests.items()):
        reasons = man.get("dropped_reason") or {}
        if not reasons and not man.get("n_dropped"):
            continue
        per_reason = collections.Counter(reasons.values())
        trait_screen[name] = {
            "fold": man.get("fold"), "domain": man.get("domain"),
            "n_tips": man.get("n_tips"), "n_dropped": man.get("n_dropped"),
            "by_reason": {k: per_reason[k] for k in sorted(per_reason)},
            "retry_used": bool(man.get("retry_used")),
            "n_full_prevalence_ec_reinstated": per_reason.get("full", 0)}
    reinstated_by_run = {r["name"]: r.get("full_prevalence_ec_reinstated")
                         for r in p2_runs
                         if r.get("n_full_prevalence_ec_reinstated")}

    # ---- annotation ceiling + headline ------------------------------------- #
    ceiling = {"recall": None, "precision": None, "source": args.annotation_ceiling}
    if os.path.exists(args.annotation_ceiling):
        try:
            cj = jload(args.annotation_ceiling)["ceiling"]["shared"]
            ceiling.update({"recall": cj["recall"]["mean"],
                            "precision": cj["precision"]["mean"]})
        except Exception:
            pass
    med_ref = int(st.median(ref_counts)) if ref_counts else 0
    cov = (out_arms.get(primary_prfba, {}).get("ec", {}).get("coverage")
           if primary_prfba else None)
    prim_key = "%s__vs__picrust2" % primary_prfba if primary_prfba else None
    prim = (contrasts.get(prim_key, {}).get("ec") or {}) if prim_key else {}
    # DENOMINATOR TRIPWIRE: the two population F1s in the sentence must subtract to
    # the paired delta beside them.  They only do so when every arm was scored over
    # the SAME cells, so this one assertion catches the whole class of
    # unequal-denominator defects (a domain's run absent, an arm skipped on some
    # cells, a stale cells table) generically.  5e-4 is the rounding slack of three
    # 4-decimal values.
    f1_prim = pop_f1(primary_prfba, "ec") if primary_prfba else None
    f1_p2 = pop_f1("picrust2", "ec")
    delta_pop = prim.get("delta_f1_population")
    denom_gap = None
    if None not in (f1_prim, f1_p2, delta_pop):
        denom_gap = (f1_prim - f1_p2) - delta_pop
    denominators_ok = denom_gap is None or abs(denom_gap) <= 5e-4
    sentence = (
        "leave-%s-out, %d species-dereplicated reference genomes (median per fold at "
        "density %s), coverage %s: %s population EC-F1 %s vs picrust2 %s, "
        "delta %s (95%% CI %s, cluster bootstrap over %s holdout groups, n_eff %s); "
        "skill score vs constant/clade-max floor-ceiling %s; SESOI %s %s."
        % (rank_blocked, med_ref, primary_dk,
           ("%.3f" % cov) if cov is not None else "n/a",
           primary_prfba or "prfba(absent)",
           pop_f1(primary_prfba, "ec") if primary_prfba else "n/a",
           pop_f1("picrust2", "ec"),
           prim.get("delta_f1_population"), prim.get("ci95"),
           prim.get("bootstrap_n_groups"), prim.get("n_eff"),
           (out_arms.get(primary_prfba, {}).get("ec", {}).get("skill_score")
            if primary_prfba else None),
           FROZEN["sesoi_delta_f1"],
           {"exceeded": "exceeded", "not_exceeded": "not exceeded"}.get(
               prim.get("sesoi_verdict"),
               "INDETERMINATE -- no CI was computed (%s), so this run cannot say "
               "whether the SESOI is exceeded or not"
               % (prim.get("ci_unavailable_reason") or "no primary contrast"))))
    # An incomplete prediction table does NOT make the denominators differ (the
    # missing md5's cells are scored 0.0 for that arm, so n stays equal) -- it moves
    # the arm's COVERAGE instead.  C14 fails on it and score_folds exits non-zero,
    # but the sentence itself is stamped too, because summary_<rank>.json is what the
    # write-up reads and a quoted sentence must not survive a void run.
    comp = completeness or {}
    n_sentinel_rows = int(comp.get("n_missing_record_rows") or 0)
    runs_complete = comp.get("runs_complete")
    plan_cov = comp.get("plan_coverage") or {}
    plan_cov_bad = plan_cov.get("ok") is False
    # ---- C14 leg 6: the derived-vs-discovered run inventory --------------- #
    inv = comp.get("inventory") or {}
    inv_bad = inv.get("ok") is False
    inv_gated = bool(inv.get("gated"))
    # ---- C14 leg 7: EVERY arm on the SAME denominator --------------------- #
    # The leg-4 tripwire is pairwise (primary prFBA vs picrust2), so a NON-primary
    # arm -- including prfba_nestedcv_forced, itself a headline number and SPEC D7's
    # risk-coverage/AURC source -- could sit on a different population at EXIT 0
    # (observed: n=1258 against picrust2's 1922).  Every arm is scored over the same
    # planned cells by construction, so one set-size assertion covers the class.
    n_by_arm = {}
    n_by_space = collections.defaultdict(set)
    for _arm, _spaces in sorted(out_arms.items()):
        for _space, _blk in sorted(_spaces.items()):
            _n = ((_blk.get("population") or {}).get("f1") or {}).get("n")
            if _n:
                n_by_arm["%s|%s" % (_arm, _space)] = int(_n)
                n_by_space[_space].add(int(_n))
    arm_denoms_unequal = any(len(v) > 1 for v in n_by_space.values())
    comp["arm_denominators"] = {
        "gated": inv_gated,
        "n_by_arm": n_by_arm,
        "distinct_n_by_space": {k: sorted(v) for k, v in sorted(n_by_space.items())},
        "ok": (not arm_denoms_unequal) if inv_gated else None}
    arm_denoms_bad = inv_gated and arm_denoms_unequal
    # A PICRUSt2 (fold, domain) that was never LAUNCHED keeps every denominator
    # equal by design (SPEC D22), which is why it slipped past legs 1-5: it costs
    # 10.4% of the population at EXIT 0 behind a fully formed headline that never
    # mentions Archaea.  In the gated modes the omission now voids the run unless
    # the operator declared the narrower scope with --expect-domains.
    dom_not_run = omitted_picrust2_domains or {}
    dom_not_run_bad = bool(dom_not_run) and inv_gated and not args.expect_domains
    # ---- the primary claim must be DECIDABLE ------------------------------ #
    # A run whose primary contrast has no confidence interval cannot state the
    # campaign's claim either way; publishing "SESOI not exceeded" for it inverted
    # the conclusion at EXIT 0.  An interval from fewer than MIN_BOOTSTRAP resamples
    # is quotable but not supportable, so it is flagged the same way.
    claim_indeterminate = bool(prim) and prim.get("sesoi_verdict") == \
        "indeterminate_no_ci"
    claim_underpowered = bool(prim.get("ci_underpowered"))
    claim_bad = inv_gated and (claim_indeterminate or claim_underpowered)
    run_incomplete = (bool(n_sentinel_rows) or runs_complete is False or plan_cov_bad
                      or inv_bad or arm_denoms_bad or dom_not_run_bad)
    if run_incomplete:
        sentence = (
            "RUN_INCOMPLETE -- these numbers are not valid: %d run(s) shipped a "
            "prediction table that does not cover every md5 the fold plan queried, "
            "%d cell row(s) were scored against a missing record (abstain_reason "
            "in %s), %d planned (fold, density) cell(s) were never scored at all "
            "because NEITHER method has a run for them%s, the run INVENTORY is %s "
            "(missing %s, duplicate %s, unexpected %s, scope %s, input bytes "
            "changed while scoring %s), the per-arm "
            "denominators are %s%s, and PICRUSt2 (fold, domain) run(s) never "
            "launched: %s -- the scored population is then not the planned one. See "
            "control C14 and accounting.completeness; rerun the affected run(s) and "
            "rescore with --force. The provisional sentence was: %s"
            % (comp.get("n_runs_incomplete", 0), n_sentinel_rows,
               list(MISSING_RECORD_REASONS), plan_cov.get("n_grid_missing") or 0,
               (" %s" % (plan_cov.get("grid_missing") or [])[:5]
                if plan_cov_bad else ""),
               "INCOMPLETE" if inv_bad else "complete",
               (inv.get("missing") or [])[:5],
               [d.get("dirs") for d in (inv.get("duplicates") or [])][:5],
               [d.get("key") for d in (inv.get("unexpected") or [])][:5],
               [d.get("defect") for d in (inv.get("scope_mismatch") or [])][:5],
               [d.get("file") for d in
                (inv.get("inputs_changed_during_scoring") or [])][:5],
               "UNEQUAL" if arm_denoms_bad else "equal",
               (" %s" % comp["arm_denominators"]["distinct_n_by_space"]
                if arm_denoms_bad else ""),
               dom_not_run if dom_not_run_bad else {},
               sentence))
        print("[score_folds] *** %s" % sentence)
    if not denominators_ok:
        n_prim = ((blocks.get((primary_prfba, "ec", "all", "all")) or
                   ArmAgg()).n if primary_prfba else 0)
        n_p2 = (blocks.get(("picrust2", "ec", "all", "all")) or ArmAgg()).n
        sentence = (
            "DENOMINATORS_DIFFER -- headline withheld: %s population EC-F1 %s is "
            "over %d cells, picrust2 %s over %d, so they do not subtract to the "
            "paired delta %s (gap %s). Fix the missing arm/run coverage (see "
            "meta.notes.picrust2_domains_not_run and accounting.completeness) and "
            "rescore; the suppressed sentence was: %s"
            % (primary_prfba, f1_prim, n_prim, f1_p2, n_p2, delta_pop,
               rnd(denom_gap, 6), sentence))
        print("[score_folds] *** %s" % sentence)

    # Per-ARM validity (operator trap): the arms table is written in full even when
    # the run is declared void, so a consumer reading arms.<arm>.<space>.population.f1
    # without checking headline.valid would quote a number the scorer itself
    # invalidated.  Mark it where that number lives, and mirror the flags in meta.
    run_valid = bool(denominators_ok) and not run_incomplete and not claim_bad
    invalid_reasons = []
    # Only a GATED mode's claim voids the run (the control reference families carry
    # no primary claim); every mode still reports the verdict per contrast.
    if claim_bad and claim_indeterminate:
        invalid_reasons.append("primary_claim_indeterminate_%s"
                               % (prim.get("ci_unavailable_reason") or "no_ci"))
    if claim_bad and claim_underpowered:
        invalid_reasons.append("primary_claim_%s"
                               % (prim.get("ci_underpowered_reason") or
                                  "ci_underpowered"))
    if not denominators_ok:
        invalid_reasons.append("denominators_differ")
    if n_sentinel_rows:
        invalid_reasons.append("cells_scored_without_a_record")
    if runs_complete is False:
        invalid_reasons.append("launched_run_incomplete")
    if plan_cov_bad:
        invalid_reasons.append("planned_fold_density_never_scored")
    for _r in (inv.get("reasons") or []):   # includes inputs_changed_during_scoring
        if _r not in invalid_reasons:
            invalid_reasons.append(_r)
    if arm_denoms_bad:
        invalid_reasons.append("arm_denominators_unequal")
    if dom_not_run_bad:
        invalid_reasons.append("picrust2_domain_never_launched")
    # (stamped into every arm block, meta and headline by stamp_invalid below, once
    # the summary dict exists -- one writer, so the three copies cannot disagree)

    n_unscoreable_empty = sum(c.get("n_truth_empty", 0) for c in per_fold.values())
    n_unscoreable_min = sum(c.get("n_min_align", 0) for c in per_fold.values())
    n_unscoreable_notrun = sum(c.get("n_method_not_run", 0)
                               for c in per_fold.values())
    summary = {
        "meta": {
            # mirror of headline.valid, so a reader of meta alone cannot miss it
            "run_valid": run_valid,
            "run_invalid_reasons": list(invalid_reasons),
            "frozen": run_meta["frozen"], "rank_blocked": rank_blocked,
            "n_folds": len(plan.get("folds", [])),
            "n_eligible": len(panel),
            "n_reference_genomes_median": med_ref,
            "n_cells_scored": sum(c.get("n_scoreable", 0) for c in per_fold.values()),
            "n_unscoreable_empty_truth": n_unscoreable_empty,
            "n_unscoreable_min_align": n_unscoreable_min,
            # a cell of a domain whose PICRUSt2 run was never launched: unscoreable
            # for EVERY arm, so all denominators stay equal
            "n_unscoreable_method_not_run": n_unscoreable_notrun,
            # the two scalars above are the EC-space anchor-arm counts; the rxn space
            # excludes different cells (an EC set that projects to no ModelSEED
            # reaction), which a reader cannot otherwise reconstruct
            "n_unscoreable_by_space": (completeness or {}).get(
                "unscoreable_by_space", {}),
            "unscoreable_uniform_across_arms": (completeness or {}).get(
                "unscoreable_uniform_across_arms"),
            "generated": now_iso(), "git_commit": run_meta["git_commit"],
            "notes": {
                "mode": args.mode, "script": "score_folds.py",
                "primary_arm": primary_prfba,
                "primary_arm_definition":
                    ("nested-CV composite: fold i scored with the select mode that won "
                     "on folds != i (SPEC D8); per-fold winners in nested_cv"
                     if primary_prfba == NESTED_ARM else
                     "single select mode (no nested-CV composite available)"),
                "nested_cv_winners": {str(k): v for k, v in sorted(winners.items())},
                "accounting_anchor_arm": anchor,
                "folds_scored": sorted({str(m[5]) for m in cell_meta.values()}),
                "n_folds_in_plan": len(plan.get("folds", [])),
                "picrust2_domains_not_run":
                    omitted_picrust2_domains or {},
                "primary_density": primary_dk,
                "densities": densities, "spaces": spaces,
                "matched_subset": "cells where %s and picrust2 both predicted"
                                  % (primary_prfba or "the primary prFBA arm"),
                "matched_n": len(matched),
                "population_view": "a blank f1 cell counts as 0.0; empty-truth and "
                                   "min_align-dropped cells are excluded from every "
                                   "denominator for every arm",
                "contrast_density": primary_dk,
                "bootstrap": args.bootstrap,
                "min_bootstrap": MIN_BOOTSTRAP,
                "primary_claim": {
                    "sesoi": FROZEN["sesoi_delta_f1"],
                    "sesoi_verdict": prim.get("sesoi_verdict"),
                    "ci_available": prim.get("ci_available"),
                    "ci_unavailable_reason": prim.get("ci_unavailable_reason"),
                    "ci_underpowered": prim.get("ci_underpowered"),
                    "gate": "a run whose primary contrast has no CI cannot state "
                            "the claim either way; 'indeterminate_no_ci' is NOT "
                            "'the effect is below the SESOI'"},
                "cells": os.path.abspath(args.cells),
                "n_picrust2_runs": len(p2_runs), "n_prfba_runs": len(pf_runs),
                # provenance: which run keys, and which BYTES, produced these
                # numbers (accounting.completeness.inventory.fingerprints has the
                # per-file size + sha256 of every predictions.json/manifest.json)
                "run_inventory_fp": inv.get("inventory_fp"),
                # the answer key is NAMED here, not merely hashed: every F1 in this
                # file is against the EC sets in meta.notes.answer_key_paths
                "answer_key_fp": inv.get("answer_key_fp"),
                "answer_key_paths": inv.get("input_paths") or {},
                "scorer_code_fp": inv.get("scorer_code_fp"),
                "run_inventory_gated": inv.get("gated"),
                "n_run_directories": inv.get("n_run_directories"),
                "denominators_consistent": denominators_ok,
                "denominator_check": {
                    "population_f1_primary": f1_prim,
                    "population_f1_picrust2": f1_p2,
                    "delta_f1_population": delta_pop,
                    "gap": rnd(denom_gap, 6), "tolerance": 5e-4},
            },
        },
        "accounting": {"per_fold": acc_rows,
                       "completeness": completeness or {},
                       "trait_column_screen": trait_screen,
                       "full_prevalence_ec_reinstated_by_run": reinstated_by_run},
        "arms": out_arms,
        "contrasts": contrasts,
        "nested_cv": nested,
        "annotation_ceiling": ceiling,
        "headline": {"sentence": sentence, "rank_blocked": rank_blocked,
                     "reference_density_genomes": med_ref,
                     "coverage": rnd(cov),
                     "denominators_consistent": denominators_ok,
                     "prediction_runs_complete": not run_incomplete,
                     "plan_coverage_complete": not plan_cov_bad,
                     "run_inventory_complete": not inv_bad,
                     "arm_denominators_equal": not arm_denoms_bad,
                     "run_inventory_fp": inv.get("inventory_fp"),
                     # WHICH answer key and WHICH code produced this sentence
                     "answer_key_fp": inv.get("answer_key_fp"),
                     "scorer_code_fp": inv.get("scorer_code_fp"),
                     "sesoi_verdict": prim.get("sesoi_verdict"),
                     "primary_claim_decidable": not claim_indeterminate,
                     "valid": run_valid,
                     # filled by stamp_invalid in main() once every hard control has
                     # been evaluated; a hard control failure makes this list
                     # non-empty and valid false
                     "hard_controls_failed": []},
    }
    stamp_invalid(summary, invalid_reasons)
    return summary, {"primary_arm": primary_prfba, "matched": matched,
                     "per_fold": per_fold, "arms": arms, "blocks": blocks,
                     "paired": paired, "groups": groups}


def hard_control_failures(fresh):
    """The HARD controls that failed in THIS invocation's own verdicts.

    ONE rule, used by the exit status and by the validity stamp so they cannot
    disagree: a control is HARD unless its own ``gate`` text begins with "reported"
    (C7/C10/C13/C15; a future reported-only control MUST start its gate text with
    that word).
    """
    return sorted(cid for cid, v in (fresh or {}).items()
                  if v.get("status") == "fail" and not
                  str(v.get("gate") or "").strip().lower().startswith("reported"))


def stamp_invalid(summary, reasons):
    """Write ONE validity verdict into every place a consumer reads it.

    ``headline.valid`` / ``meta.run_valid`` / ``arms.<arm>.<space>.run_valid`` used to
    reflect only C14 and the inventory, so a C4/C5/C6/C9/C11 failure exited 3 while
    leaving a fully formed ``"valid": true`` artefact whose published skill score had
    flipped sign (-2.0449 -> +1.5523).  Any consumer reading summary_*.json rather
    than the exit code was unprotected, so validity is stamped HERE, from one reason
    list, after every hard control has been evaluated.
    """
    meta = summary.setdefault("meta", {})
    head = summary.setdefault("headline", {})
    out = list(meta.get("run_invalid_reasons") or [])
    for r in reasons:
        if r not in out:
            out.append(r)
    valid = not out
    meta["run_valid"], meta["run_invalid_reasons"] = valid, out
    head["valid"], head["run_invalid_reasons"] = valid, out
    for _arm, _spaces in (summary.get("arms") or {}).items():
        for _space, _blk in (_spaces or {}).items():
            if isinstance(_blk, dict):
                _blk["run_valid"] = valid
                _blk["run_invalid_reasons"] = list(out)
    sent = head.get("sentence") or ""
    if not valid and not sent.startswith(("RUN_INCOMPLETE", "DENOMINATORS_DIFFER",
                                          "INVALID")):
        head["sentence"] = ("INVALID (%s) -- these numbers are not valid and must "
                            "not be quoted; see controls and "
                            "accounting.completeness. The provisional sentence was: "
                            "%s" % (", ".join(out), sent))
    return valid, out


def fold_abstained(blocks, arm, fold):
    if not arm:
        return 0
    b = blocks.get((arm, "ec", "fold", str(fold)))
    return (b.n - b.called) if b is not None else 0


def strat_block(blocks, arm, space, kind):
    out = {}
    for (a, s, k, v), b in blocks.items():
        if a == arm and s == space and k == kind:
            out[v] = {"population": b.population_block(), "n": b.n,
                      "coverage": b.coverage(),
                      "micro": b.micro_block()}
    return {k: out[k] for k in sorted(out)}


# --------------------------------------------------------------------------- #
# controls
# --------------------------------------------------------------------------- #
def phenotype_binding(args, plan, frozen, inv, cells_fp_rec):
    """Bind ``phenotype_summary.json`` to THIS run, or say loudly that it is not.

    The phenotype summary was auto-discovered by EXISTENCE and trusted: a DIFFERENT
    campaign's 12-genome summary dropped in beside ``--cells`` scored C11/C12 pass,
    and a truncated one made C11/C12 vanish silently (the whole block was wrapped in
    ``except Exception: pass``, so an unreadable phenotype arm was indistinguishable
    from a tree where the arm never ran).

    ``phenotype_fold.py`` now stamps ``meta.binding`` -- the sha256 of the cells
    table it scored, of the fold plan, the panel and the ANSWER KEY it used, the run
    directories it read, and its own genome list.  Every one of those is re-checked
    here against this invocation's inputs.  A summary that predates the contract
    carries no ``meta.binding``: it is reported as ``bound: false`` (and fails under
    ``--require-phenotype-binding``, which SPEC 11.2's campaign invocation passes)
    rather than being silently accepted as this run's phenotype arm.
    """
    path = args.phenotype_summary
    out = {"path": (os.path.abspath(path) if path else None),
           "discovered": ("--phenotype-summary" if getattr(args, "_pheno_explicit",
                                                            False)
                          else "auto-discovered beside --cells"),
           "present": bool(path and os.path.exists(path)),
           "bound": False, "defects": [], "readable": None,
           "require_binding": bool(getattr(args, "require_phenotype_binding", False)),
           "fingerprints": {}, "declared": {}}
    if not out["present"]:
        return out
    out["fingerprints"]["phenotype_summary"] = file_fp(path)
    try:
        pj = jload(path)
        if not isinstance(pj, dict):
            raise ValueError("not a JSON object")
        meta = pj.get("meta")
        if not isinstance(meta, dict):
            raise ValueError("no meta block")
    except Exception as exc:
        out["readable"] = False
        out["defects"].append({"defect": "phenotype_summary_unreadable",
                               "detail": "%s: %s" % (type(exc).__name__, exc)})
        return out
    out["readable"] = True
    out["declared"] = {
        "n_genomes": meta.get("n_genomes"), "n_models": meta.get("n_models"),
        "select_mode": meta.get("select_mode"), "prfba_arm": meta.get("prfba_arm"),
        "density": meta.get("density"), "generated": meta.get("generated"),
        "git_commit": meta.get("git_commit"),
        "n_gated_out": meta.get("n_gated_out"),
        "identity_control": meta.get("identity_control")}
    # ---- checks that need no binding block (a foreign campaign's config) --- #
    if meta.get("frozen") and frozen and meta["frozen"] != frozen:
        out["defects"].append({
            "defect": "frozen_config_differs_from_this_fold_plan",
            "detail": "the phenotype arm was run under a different frozen config"})
    # Every run directory THIS TREE contains, per root -- not the mode-selected
    # subset: the phenotype arm legitimately reads the include-self run too, so
    # comparing against --mode fold's selection would flag a correct summary.  A root
    # this invocation was not given at all (--mode plumbing takes no --prfba) can
    # neither confirm nor deny its runs, so those are recorded as UNCHECKED rather
    # than reported as foreign.
    run_dirs = {}
    for kind, root in (("picrust2", args.picrust2), ("prfba", args.prfba)):
        if root and os.path.isdir(root):
            run_dirs[kind] = {os.path.basename(os.path.dirname(q)) for q in
                              glob.glob(os.path.join(root, "*", "predictions.json"))}
    declared, unchecked, foreign = set(), set(), []
    for key, kind in (("picrust2_runs", "picrust2"), ("prfba_runs", "prfba")):
        for r in (meta.get(key) or ()):
            if not isinstance(r, dict):
                continue
            d = os.path.basename(os.path.dirname(str(r.get("path") or "")))
            if not d:
                continue
            declared.add(d)
            if kind not in run_dirs:
                unchecked.add("%s/%s" % (kind, d))
            elif d not in run_dirs[kind]:
                foreign.append("%s/%s" % (kind, d))
    out["declared"]["run_dirs"] = sorted(declared)
    out["declared"]["run_dirs_unchecked_root_not_given"] = sorted(unchecked)
    if foreign:
        out["defects"].append({
            "defect": "phenotype_read_run_directories_this_scoring_does_not_cover",
            "detail": sorted(foreign)})
    # ---- the binding block itself ----------------------------------------- #
    b = meta.get("binding")
    if not isinstance(b, dict) or not b.get("contract"):
        out["bound"] = False
        out["unbound_reason"] = (
            "phenotype_summary.json carries no meta.binding block: it predates the "
            "binding contract, so nothing ties it to this cells table, this fold "
            "plan or this answer key. Re-run phenotype_fold.py to bind it (or pass "
            "--require-phenotype-binding to make this a C11 failure).")
        if out["require_binding"]:
            out["defects"].append({"defect": "phenotype_summary_unbound",
                                   "detail": out["unbound_reason"]})
        return out
    out["binding_contract"] = b.get("contract")
    out["declared"]["binding"] = {k: v for k, v in sorted(b.items())
                                  if k not in ("genomes",)}
    # taken from the Stage 0 fingerprints rather than re-hashed: one read per input
    # per run (--functions is 1.7 GB), and the comparison is then against exactly the
    # bytes THIS scoring used, which reverify_fingerprints already re-checks
    _ifp = (inv.get("fingerprints") or {}).get("inputs") or {}
    _akfp = (inv.get("fingerprints") or {}).get("answer_key") or {}
    ours = {"cells": cells_fp_rec, "fold_plan": _ifp.get("fold_plan"),
            "panel": _ifp.get("panel"), "functions": _akfp.get("functions")}
    for tag, mine in sorted(ours.items()):
        theirs = b.get(tag)
        if not isinstance(theirs, dict) or not theirs.get("sha256"):
            out["defects"].append({"defect": "phenotype_binding_incomplete",
                                   "detail": "meta.binding.%s is missing" % tag})
        elif mine and theirs.get("sha256") != mine.get("sha256"):
            out["defects"].append({
                "defect": "phenotype_bound_to_a_different_%s" % tag,
                "detail": {"phenotype_arm": theirs.get("sha256"),
                           "this_run": mine.get("sha256")}})
    # the phenotype arm's own outputs must still be the bytes it summarised
    for tag, key in (("per_genome_jsonl", "per_genome_jsonl"),
                     ("identity_control_json", "identity_control_json")):
        theirs = b.get(key)
        if not isinstance(theirs, dict):
            continue
        # The recorded path is absolute and the whole tree legitimately gets copied
        # (every review of this harness runs on a copy of data_smoke_e2e), so the
        # file BESIDE THIS SUMMARY is checked FIRST: it is the phenotype arm this
        # tree actually ships.  Preferring the recorded path instead let an edited
        # local per_genome.jsonl pass by fingerprinting the ORIGINAL tree's copy.
        rec_path = theirs.get("path") or ""
        local = os.path.join(os.path.dirname(os.path.abspath(path)),
                             os.path.basename(str(b.get("outdir") or "")),
                             os.path.basename(rec_path))
        cand = local if os.path.exists(local) else rec_path
        out["fingerprints"]["%s_checked" % tag] = cand
        now = file_fp(cand)
        out["fingerprints"][tag] = now
        if now is None:
            out["defects"].append({"defect": "phenotype_output_missing",
                                   "detail": theirs.get("path")})
        elif now.get("sha256") != theirs.get("sha256"):
            out["defects"].append({"defect": "phenotype_output_rewritten",
                                   "detail": {"file": theirs.get("path"),
                                              "summarised": theirs.get("sha256"),
                                              "now": now.get("sha256")}})
    # the genomes it scored must be held-out genomes of THIS plan
    test_g = set()
    for f in plan.get("folds", []):
        test_g |= set(f.get("test_genomes") or ())
    declared_g = [g for g in (b.get("genomes") or ())]
    out["declared"]["n_bound_genomes"] = len(declared_g)
    outside = sorted(set(declared_g) - test_g)
    if outside:
        out["defects"].append({
            "defect": "phenotype_scored_genomes_outside_this_plans_holdout",
            "detail": outside[:10]})
    out["bound"] = not out["defects"]
    return out


def phenotype_controls(pb):
    """C11 (identity) and C12 (phenotype validity gate), from a BOUND summary.

    C12 used to be written as ``control(True, ...)``: it could not fail, in any tree,
    for any input.  Its SPEC row asked for the gated-out count to be REPORTED while
    7.1 listed it HARD, so the gate is now stated: a gated-out TRUTH model removes
    that genome from every arm's denominator (the silent-population class this
    harness exists to stop), and a pipeline gating more than
    ``C12_MAX_GATED_FRACTION`` of all models is miscalibrated rather than protective.
    """
    out = {}
    defects = pb.get("defects") or []
    decl = pb.get("declared") or {}
    ident = decl.get("identity_control") or {}
    gated = decl.get("n_gated_out")
    n_models = decl.get("n_models") or 0
    binding_measured = {"phenotype_summary": pb.get("path"),
                        "bound_to_this_run": pb.get("bound"),
                        "binding_defects": defects,
                        "unbound_reason": pb.get("unbound_reason"),
                        "fingerprints": pb.get("fingerprints")}
    if defects:
        # An unreadable / foreign / unbound-and-required phenotype summary is a
        # FAILURE of the phenotype controls, never their absence: C11/C12 vanishing
        # from controls.json read exactly like a tree where the arm never ran.
        out["C11"] = control(False, dict(binding_measured, identity_control=ident),
                             "phenotype_summary.json is bound to THIS cells table, "
                             "fold plan and answer key, and its identity control "
                             "passes (Jaccard 1.000)",
                             "an unbound or unreadable phenotype summary lets "
                             "another campaign's numbers (or none at all) stand in "
                             "for this run's phenotype arm")
        out["C12"] = control(False, dict(binding_measured, n_gated_out=gated),
                             "phenotype_summary.json is bound to THIS run, no TRUTH "
                             "model was gated out, and fewer than %.0f%% of all "
                             "models were gated" % (100 * C12_MAX_GATED_FRACTION),
                             "a phenotype validity gate that cannot be tied to this "
                             "run cannot falsify anything")
        return out
    if ident:
        out["C11"] = control(bool(ident.get("pass")),
                             dict(ident, binding=binding_measured),
                             "two builds from the same EC set give identical "
                             "reaction sets and gapfill solution (Jaccard 1.000), "
                             "from a phenotype summary bound to this run",
                             "degenerate optima are polluting every phenotype "
                             "number")
    if isinstance(gated, dict):
        truth_gated = int(gated.get("truth") or 0)
        total_gated = sum(int(v or 0) for v in gated.values())
        frac = (total_gated / n_models) if n_models else 0.0
        out["C12"] = control(
            truth_gated == 0 and frac <= C12_MAX_GATED_FRACTION,
            {"n_gated_out": gated, "n_truth_models_gated": truth_gated,
             "n_models_gated": total_gated, "n_models": n_models,
             "gated_fraction": rnd(frac), "max_gated_fraction":
                 C12_MAX_GATED_FRACTION, "binding": binding_measured},
            "biomass and ATPM on water+protons both == 0 and yield below the "
            "calibrated ceiling, with ZERO truth models gated out and at most "
            "%.0f%% of all models gated" % (100 * C12_MAX_GATED_FRACTION),
            "the raw draft yields 1000 ATP from water via rxn15249_c0 + "
            "rxn20583_c0 cycling and loopless_solution does NOT remove it; a GATED "
            "TRUTH model silently removes that genome from every arm's denominator")
    elif gated is not None:
        out["C12"] = control(False, {"n_gated_out": gated,
                                     "binding": binding_measured},
                             "meta.n_gated_out is a per-arm mapping",
                             "a malformed gate report cannot falsify anything")
    return out


def control(status, measured, gate, falsifies):
    return {"status": "pass" if status else "fail", "measured": measured,
            "gate": gate, "falsifies": falsifies}


def _base_mode(name):
    """'graduated_forced' -> 'graduated'; the forced runner suffixes its table keys
    (SPEC 6.2) while meta.select_modes records the bare mode names."""
    name = str(name)
    return name[:-len("_forced")] if name.endswith("_forced") else name


def run_completeness(plan, panel, truth_meta, md5_of_cell, p2_runs, pf_runs):
    """C14 legs 1 and 2: is every launched run's prediction table COMPLETE?

    Leg 1 (self-consistency) -- a run's own meta counts must match the table it
    shipped: ``len(pred) == n_placed == n_query - n_min_align_dropped`` for PICRUSt2,
    and one record per query per select mode for prFBA, over EXACTLY the mode set
    ``meta.select_modes`` names (a missing or empty mode block is a fail, not a
    vacuous pass).
    Leg 2 (against the fold plan) -- every md5 the plan says was queried, minus that
    run's own min_align drops, must appear in the table.  A forced-routing PICRUSt2
    run is only responsible for its own domain's md5s; an nsti run (domain "both")
    for all of them.

    Without this, a truncated / resumed-stale / plan-mismatched run is laundered into
    method abstentions: the md5 is simply absent, so every one of its cells scores
    0.0 and the arm reads out at lower coverage.  Measured cost of not checking:
    deleting 60 of 561 md5s from one fold's bacterial run moved the primary delta
    from -0.2610 to -0.2347 with every control still passing.
    """
    dom_of_md5 = collections.defaultdict(set)
    md5s_of_fold = {}
    for f in plan.get("folds", []):
        fk = int(f["fold"])
        md5s = set()
        for cell in (f.get("test_cells") or ()):
            md5 = md5_of_cell.get(cell)
            if not md5:
                continue
            md5s.add(md5)
            t = truth_meta.get(cell) or {}
            d = t.get("domain") or (panel.get(t.get("src_genome_id") or
                                              cell.rsplit("__", 1)[0]) or {}
                                    ).get("domain") or ""
            dom_of_md5[md5].add(dom3(d))
        md5s_of_fold[fk] = md5s | set(f.get("test_md5s") or ())
    out, ok = {}, True
    for r in list(p2_runs) + list(pf_runs):
        meta = r.get("meta") or {}
        is_p2 = "pred" in r
        if is_p2:
            tables = {"pred": set(k[0] if isinstance(k, tuple) else k
                                  for k in r["pred"])}
        else:
            tables = {m: set(v) for m, v in (r.get("modes") or {}).items()}
        # The EXPECTED prFBA table set comes from the run's own meta.select_modes,
        # never from the tables it shipped: deriving it from the output made a
        # MISSING mode block vacuously complete (both legs loop over `tables`, so a
        # deleted key leaves `missing` empty), and a vanished select mode silently
        # shrinks that arm -- prfba_nestedcv_forced fell 1922 -> 1258 cells against
        # PICRUSt2's 1922 with C14 passing.  Compared on BASE mode names because the
        # forced runner names its keys "<mode>_forced" (SPEC 6.2) and score_cells
        # appends the suffix itself when the runner did not.
        want_modes = {_base_mode(m) for m in (meta.get("select_modes") or ())}
        got_modes = {_base_mode(m) for m in tables}
        n_query = meta.get("n_query")
        dropped = {str(m) for m in (meta.get("min_align_dropped_md5") or ())}
        rec = {"mode": meta.get("mode"), "fold": r.get("fold"),
               "domain": meta.get("domain"), "routing": r.get("routing"),
               "n_query": n_query, "n_placed": meta.get("n_placed"),
               "n_min_align_dropped": meta.get("n_min_align_dropped"),
               "n_rows": {k: len(v) for k, v in sorted(tables.items())},
               "self_consistent": True, "plan_complete": None, "n_missing": 0,
               "missing_examples": []}
        if not is_p2:
            rec["select_modes"] = sorted(want_modes)
            rec["tables_present"] = sorted(got_modes)
            rec["tables_missing"] = sorted(want_modes - got_modes)
            rec["tables_unexpected"] = sorted(got_modes - want_modes)
            # An ABSENT meta.select_modes used to short-circuit this guard (empty
            # want_modes), which reopened exactly the blind spot it was added to
            # close: a resumed/older/hand-assembled run block could drop a mode
            # and pass vacuously.  No declared scope == not self-consistent.
            if not tables or not want_modes or got_modes != want_modes:
                rec["self_consistent"] = False
        # ---- leg 1 ---------------------------------------------------------- #
        for name, tbl in sorted(tables.items()):
            if is_p2:
                if meta.get("n_placed") is not None and len(tbl) != meta["n_placed"]:
                    rec["self_consistent"] = False
                if n_query is not None and len(tbl) != (
                        int(n_query) - int(meta.get("n_min_align_dropped") or 0)):
                    rec["self_consistent"] = False
            elif n_query is not None and len(tbl) != int(n_query):
                rec["self_consistent"] = False
        # ---- leg 2 ---------------------------------------------------------- #
        fk = r.get("fold")
        if fk is not None and (meta.get("mode") or "fold") != "recipe_control":
            want = set(md5s_of_fold.get(int(fk), ()))
            dkey_dom = dom3(meta.get("domain")) if is_p2 else ""
            if dkey_dom:
                want = {m for m in want if dkey_dom in dom_of_md5.get(m, {dkey_dom})}
            want -= dropped
            missing = set()
            for tbl in tables.values():
                missing |= (want - tbl)
            rec["n_expected_from_plan"] = len(want)
            rec["plan_complete"] = not missing
            rec["n_missing"] = len(missing)
            rec["missing_examples"] = sorted(missing)[:5]
        if not rec["self_consistent"] or rec["plan_complete"] is False:
            ok = False
        out[r["name"]] = rec
    return {"per_run": out, "ok": ok,
            "n_runs_incomplete": sum(1 for v in out.values()
                                     if not v["self_consistent"] or
                                     v["plan_complete"] is False)}


def control_provenance(args, inv, cells_fp_rec, mode):
    """What a control verdict must carry to be attributable to a run.

    Control entries carried no timestamp, no cells sha and no inventory_fp, and
    ``merge_controls`` preserves whatever was already in the file -- so five
    pre-seeded "pass" entries made a fresh run read 15/15 pass with four of them
    never evaluated by that invocation at all.
    """
    return {"mode": mode, "evaluated_at": now_iso(),
            "cells": os.path.abspath(args.cells) if args.cells else None,
            "cells_sha256": (cells_fp_rec or {}).get("sha256"),
            "inventory_fp": inv.get("inventory_fp"),
            "answer_key_fp": inv.get("answer_key_fp"),
            "scorer_code_fp": inv.get("scorer_code_fp"),
            "git_commit": git_commit()}


def controls_provenance_report(merged, fresh, prov, mode):
    """Which entries in controls.json THIS invocation actually evaluated.

    Reported, never a gate by itself: C1 legitimately comes from
    ``--mode recipe_control``, C2 from include_self and C9 from plumbing, so a
    foreign-mode entry is normal.  What was NOT distinguishable before is a stale or
    pre-seeded entry from an unknown run, which is now named.
    """
    rep = {"this_run": prov, "fresh": sorted(fresh),
           "from_other_modes": [], "without_provenance": [],
           "stale_same_mode": [], "missing_expected": [],
           "expected_this_mode": sorted(EXPECTED_CONTROLS.get(mode, ())),
           "all_control_ids": list(ALL_CONTROLS)}
    for cid in sorted(merged):
        if cid in fresh:
            continue
        pv = (merged[cid] or {}).get("provenance")
        if not isinstance(pv, dict):
            rep["without_provenance"].append(cid)
            continue
        # Classify by MODE first.  inventory_fp can never match across modes (the
        # mode and the mode-selected run set are inside it), so comparing
        # fingerprints first labelled every legitimate cross-mode control -- C1 from
        # recipe_control, C2 from include_self, C9 from plumbing -- "different
        # inputs/code", which at full scale is a cry-wolf in every campaign log.
        # What the ANSWER KEY and the CODE agreeing does tell us is reported beside
        # each entry instead.
        rec = {"id": cid, "mode": pv.get("mode"),
               "evaluated_at": pv.get("evaluated_at"),
               "cells_sha256": pv.get("cells_sha256"),
               "inventory_fp": pv.get("inventory_fp"),
               "answer_key_matches":
                   pv.get("answer_key_fp") == prov.get("answer_key_fp"),
               "scorer_code_matches":
                   pv.get("scorer_code_fp") == prov.get("scorer_code_fp")}
        if pv.get("mode") != mode:
            rep["from_other_modes"].append(rec)
        else:
            # this mode recorded it once and did NOT re-evaluate it now
            rep["stale_same_mode"].append(rec)
    rep["missing_expected"] = sorted(set(EXPECTED_CONTROLS.get(mode, ())) -
                                     set(fresh))
    rep["conditional_not_emitted"] = {
        cid: why for cid, why in sorted(CONDITIONAL_CONTROLS.items())
        if cid not in fresh}
    rep["n_absent_entirely"] = sorted(set(ALL_CONTROLS) - set(merged))
    return rep


def merge_controls(path, new, mode=None, provenance=None):
    """Read-modify-write: recipe_control, include_self, plumbing and fold invocations
    each contribute a different subset of C1..C13 to the SAME file.

    The merge is FAILURE-PRESERVING.  C3/C4/C6/C7/C8 are evaluated in every mode, so a
    plain dict update would let the include-self invocation's pass silently overwrite
    the fold arm's failure on the same gate -- i.e. the run would look clean because a
    control reference happened to be wired correctly.  A fail therefore survives a
    later pass, with the passing mode's measurement kept beside it.

    The protection is CROSS-mode only.  Re-scoring the SAME mode after fixing what the
    control caught (the operator path C14's banner prints: rerun the run, rescore with
    --force) must be able to clear that mode's own fail, or the stale entry would make
    score_folds exit 2 forever until controls.json were deleted by hand.
    """
    cur = {}
    if os.path.exists(path):
        try:
            cur = jload(path)
        except Exception:
            cur = {}
    for cid, entry in sorted(new.items()):
        if mode:
            entry["measured"] = dict(entry["measured"], mode=mode)
        # every FRESH verdict is stamped with the run that produced it, so a stale or
        # foreign "pass" is no longer indistinguishable from a just-evaluated one
        if provenance:
            entry["provenance"] = dict(provenance)
        old = cur.get(cid)
        same_mode = ((old or {}).get("measured") or {}).get("mode") == mode
        if (old or {}).get("status") == "fail" and entry["status"] == "pass" \
                and not same_mode:
            old.setdefault("measured", {}).setdefault("passed_in_other_modes", {})[
                mode or "unknown"] = entry["measured"]
            continue
        cur[cid] = entry
    jdump(path, {k: cur[k] for k in sorted(cur)})
    return cur


def fold_controls(args, plan, panel, manifests, p2_runs, pf_runs, neighbors, summary,
                  extras, p2_all=None):
    out = {}
    # ---- C3 tip-label round trip ------------------------------------------ #
    rt = {r["name"]: r["meta"].get("unknown_rows_match_queries")
          for r in p2_runs if "unknown_rows_match_queries" in r["meta"]}
    if rt:
        out["C3"] = control(all(bool(v) for v in rt.values()),
                            {"per_run": rt,
                             "n_unknown_rows": {r["name"]: r["meta"].get("n_unknown_rows")
                                                for r in p2_runs},
                             "n_query": {r["name"]: r["meta"].get("n_query")
                                         for r in p2_runs}},
                            "n_unknown_output_rows == n_query_sequences for every run",
                            "castor_nsti.R read numeric BV-BRC ids as numeric and "
                            "deleted reference genomes from the NSTI known set "
                            "(3.61% of panel ids end in a trailing zero); the g_ "
                            "label prefix is what prevents it")
    # ---- C4 artefact set equality ----------------------------------------- #
    eq = {n: bool(m.get("set_equality_ok")) for n, m in manifests.items()}
    prfba_eq = {}
    for r in pf_runs:
        want = set()
        for name in r["ref_names"]:
            want |= set((manifests.get(name) or {}).get("reference_genomes") or ())
        n_md5 = r["meta"].get("n_reference_md5")
        prfba_eq[r["name"]] = {"n_reference_md5": n_md5, "n_manifest_union": len(want),
                               "ok": (not want) or n_md5 == len(want)}
    if eq or prfba_eq:
        out["C4"] = control(all(eq.values()) and all(v["ok"] for v in prfba_eq.values()),
                            {"manifests": eq, "prfba": prfba_eq},
                            "set(fna ids) == set(tre tips) == set(trait rows) == "
                            "set(marker rows) == set(prFBA ref_db genomes)",
                            "one differing id silently turns reference genomes into "
                            "study sequences")
    # ---- C6 per-fold leakage ---------------------------------------------- #
    pre = plan.get("leakage_preflight") or {}
    recomputed = {}
    if neighbors:
        for f in plan.get("folds", []):
            ref = set(f.get("reference_genomes") or ())
            test = set(f.get("test_genomes") or ())
            md5_share = amp_ge = full_ge = tax = 0
            for g in sorted(test):
                for pair in (neighbors.get("md5_share", {}).get(g) or ()):
                    if pair and pair[0] in ref:
                        md5_share += 1
                for gid, ident in (neighbors.get("full16s", {}).get(g) or ()):
                    if gid in ref and ident >= FROZEN["tau_full16s"]:
                        full_ge += 1
                for _region, lst in sorted((neighbors.get("amplicon", {}).get(g)
                                            or {}).items()):
                    for gid, ident in lst:
                        if gid in ref and ident >= FROZEN["tau_amplicon"]:
                            amp_ge += 1
                tid = (panel.get(g) or {}).get("taxon_id")
                if tid:
                    tax += sum(1 for r2 in ref
                               if (panel.get(r2) or {}).get("taxon_id") == tid)
            recomputed["fold_%s" % f["fold"]] = {
                "md5_shared_with_reference": md5_share,
                "amplicon_identity_ge_tau": amp_ge,
                "full16s_identity_ge_tau": full_ge,
                "taxon_id_shared": tax}
    canary = {r["name"]: (r["meta"].get("assertions") or {}).get(
        "exclude_self_n_removed_total") for r in pf_runs}
    counters = [v for d in list(pre.values()) + list(recomputed.values())
                for k, v in d.items() if k != "genus_shared"]
    if pre or recomputed or canary:
        ok = all(int(v or 0) == 0 for v in counters) and \
            all(int(v or 0) == 0 for v in canary.values() if v is not None)
        out["C6"] = control(ok, {"plan_preflight": pre, "recomputed": recomputed,
                                 "exclude_self_canary": canary},
                            "md5_shared == amplicon_identity_ge_tau == "
                            "full16s_identity_ge_tau == taxon_id_shared == "
                            "exclude_self_n_removed_total == 0",
                            "the group graph is wrong and both methods got a free "
                            "near-exact neighbour")
    # ---- C7 placement-loss accounting (reported) --------------------------- #
    dropped = {r["name"]: r["meta"].get("n_min_align_dropped", 0) for r in p2_runs}
    out["C7"] = control(True, {"per_run": dropped,
                               "total": sum(int(v or 0) for v in dropped.values()),
                               "excluded_from_both_denominators": True,
                               "n_unscoreable_min_align":
                                   summary["meta"]["n_unscoreable_min_align"]},
                        "reported, not gated",
                        "scoring min_align-dropped amplicons as PICRUSt2 abstentions "
                        "would charge PICRUSt2 for a shared data limitation")
    # ---- C8 prFBA leak-adjacent hooks -------------------------------------- #
    asserts = {r["name"]: (r["meta"].get("assertions") or {}) for r in pf_runs}
    if asserts:
        ok = all(int(a.get("midas_nonempty_records", 0) or 0) == 0 and
                 int(a.get("reps_outside_reference", 0) or 0) == 0
                 for a in asserts.values())
        out["C8"] = control(ok, {"per_run": asserts},
                            "midas_nonempty_records == 0 and "
                            "reps_outside_reference == 0",
                            "a non-empty MiDAS lineage hands prFBA an external "
                            "taxonomic prior for the held-out organism")
    # ---- C10 structural recall ceiling (reported) -------------------------- #
    # TWO ceilings, because the two methods do not reach the same reference.  prFBA's
    # restricted ref_db IS the bac+arc union, so the union figure is its ceiling.
    # Under the primary --routing forced arm (D9) a PICRUSt2 query is placed on its
    # OWN domain's tree with that domain's trait table alone, so its ceiling is the
    # own-domain figure -- 3-5 recall points lower in the thin archaeal stratum,
    # which is exactly the subpopulation SPEC 8.4 insists be reported separately.
    # Quoting the union number as "the shared ceiling" overstates what PICRUSt2 can
    # physically predict there.
    ceil_rows = {str(r["fold"]): r["structural_recall_ceiling"]
                 for r in summary["accounting"]["per_fold"]}
    own_rows = {str(r["fold"]): r.get("structural_recall_ceiling_picrust2_own_tree")
                for r in summary["accounting"]["per_fold"]}
    by_dom = {str(r["fold"]): r.get("structural_recall_ceiling_by_domain") or {}
              for r in summary["accounting"]["per_fold"]}
    gaps = []
    for f, doms in sorted(by_dom.items()):
        for dom, blk in sorted(doms.items()):
            u, o = (blk.get("structural_recall_ceiling_union"),
                    blk.get("structural_recall_ceiling_own_domain"))
            if u is not None and o is not None:
                gaps.append({"fold": f, "domain": dom, "union": u, "own_domain": o,
                             "gap": rnd(u - o)})
    vals = [v for v in ceil_rows.values() if v is not None]
    own_vals = [v for v in own_rows.values() if v is not None]
    out["C10"] = control(True, {"per_fold": ceil_rows,
                                "per_fold_picrust2_own_tree": own_rows,
                                "mean": rnd(st.mean(vals)) if vals else None,
                                "mean_picrust2_own_tree":
                                    rnd(st.mean(own_vals)) if own_vals else None,
                                "distribution": summarize(vals),
                                "by_domain": by_dom,
                                "union_minus_own_domain": gaps,
                                "routing_of_primary_picrust2_arm": "forced",
                                "n_impossible_ec_mean": {
                                    str(r["fold"]): r["n_impossible_ec_mean"]
                                    for r in summary["accounting"]["per_fold"]}},
                         "reported per fold, per domain, against BOTH reference "
                         "vocabularies",
                         "ECs carried by zero surviving reference genomes cannot be "
                         "predicted: against the bac+arc UNION for prFBA (whose "
                         "ref_db is the union) and against the cell's OWN domain "
                         "reference for forced-routed PICRUSt2. The union figure "
                         "alone would overstate PICRUSt2's archaeal ceiling")
    # ---- C13 tree-quality descriptives (reported) -------------------------- #
    trees = {n: {"n_tips": m.get("n_tips"),
                 "genus_monophyly_frac": m.get("genus_monophyly_frac"),
                 "family_monophyly_frac": m.get("family_monophyly_frac"),
                 "mean_pairwise_tip_distance": m.get("mean_pairwise_tip_distance"),
                 "msa_width": m.get("msa_width"),
                 "msa_mean_gap_frac": m.get("msa_mean_gap_frac")}
             for n, m in sorted(manifests.items())}
    nsti = {}
    for arm in ("picrust2", "picrust2_nsti_routed"):
        b = summary["arms"].get(arm, {}).get("ec")
        if b:
            nsti[arm] = {"coverage": b["coverage"], "aurc": b["aurc"]}
    nsti_dist = {}
    for r in p2_runs:
        vals = [rec.get("nsti") for rec in r["pred"].values()
                if rec.get("nsti") is not None]
        if vals:
            nsti_dist[r["name"]] = summarize(vals)
    if trees:
        out["C13"] = control(True, {"references": trees, "picrust2_arms": nsti,
                                    "nsti_distribution": nsti_dist,
                                    "stock_reference_nsti_median_as_shipped": 0.0105,
                                    "recipe_nsti_ratio_ours_over_stock": 0.377},
                             "reported, not a gate",
                             "a poorly built tree would handicap PICRUSt2; genus / "
                             "family monophyly and the NSTI distribution are the "
                             "quantitative evidence it does not")
    # ---- C2 include-self ceiling ------------------------------------------- #
    if args.mode == "include_self":
        src_rec = extras.get("include_self_source_recovery")
        p2f1 = summary["arms"].get("picrust2", {}).get("ec", {}) \
                              .get("population", {}).get("f1")
        # an ABSENT measurement is not a pass: a failed or never-launched prFBA
        # include-self run would otherwise report the hard gate as satisfied
        hard_ok = src_rec is not None and abs(src_rec - 1.0) < 1e-9
        out["C2"] = control(hard_ok,
                            {"prfba_source_recovery": rnd(src_rec),
                             "prfba_source_recovery_measured": src_rec is not None,
                             "prfba_n_twin_free": extras.get("include_self_n"),
                             # SPEC 7-C2: twin cells (md5 shared by >1 panel genome)
                             # are reported against the analytic source-union-twin
                             # EC ceiling, never gated -- an exact-amplicon twin may
                             # legitimately be selected instead of the source
                             "prfba_twin_cells": extras.get("include_self_twin"),
                             "picrust2_mean_ec_f1": (p2f1 or {}).get("value"),
                             "picrust2_soft_threshold": 0.85,
                             "picrust2_soft_flag":
                                 (p2f1 or {}).get("value") is not None and
                                 (p2f1 or {}).get("value", 1.0) < 0.85,
                             "reference_measurements": {
                                 "custom_240tip_panel_dialect": 0.8155,
                                 "stock_8600tip_own_dialect": 0.9343}},
                            "HARD: prFBA source-recovery is MEASURED and == 1.00; "
                            "SOFT (flag only): PICRUSt2 mean EC F1 < 0.85",
                            "prFBA source-recovery < 1.00 means the restricted DB, "
                            "the synthetic md5->header join or the select call is "
                            "wrong; an absent measurement means the include-self "
                            "run did not happen, which is equally a failure; the "
                            "PICRUSt2 value is the density- and dialect-conditional "
                            "ceiling, never a setup failure")
    # ---- C9 plumbing ------------------------------------------------------- #
    # Evaluated from ANY plumbing run discovered under --picrust2, whatever --mode
    # this invocation is: the spec's smoke sequence scores the plumbing reference in
    # the same `--mode fold` call as the folds, so gating C9 on --mode plumbing would
    # mean the control never fires in the run the orchestrator actually makes.
    plumb = [r for r in (p2_all if p2_all is not None else p2_runs)
             if run_kind(r.get("ref_names") or [r["name"]]) == "plumbing"]
    if plumb or args.mode == "plumbing":
        src = plumb or p2_runs
        n_q = sum(int(r["meta"].get("n_query") or 0) for r in src)
        n_p = sum(len(r["pred"]) for r in src)
        # the gate names only the method it measures: these are PICRUSt2 runs.  No
        # prFBA plumbing run exists (none is specified), and prFBA's end-to-end
        # wiring is covered substantively by C2's include-self run at coverage 1.000.
        out["C9"] = control(n_q > 0 and n_p >= n_q,
                            {"n_query": n_q, "n_predicted": n_p,
                             "methods_measured": ["picrust2"],
                             "prfba_wiring_covered_by": "C2 (include_self)",
                             "runs": sorted(r["name"] for r in src)},
                            "every smoke query has a PICRUSt2 prediction and the "
                            "trait/md5 joins and domain routing resolve",
                            "wiring only -- a pass licenses nothing about tree "
                            "quality (D6)")
    # ---- C5 branch-length preservation (recipe arm A ONLY) ----------------- #
    # Arm A is the stock RAxML tree PRUNED with preserve_branch_length -- the only
    # arm where pruning happens and therefore the only one this gate can falsify.
    # Arm B is rebuilt de novo from the retained sequences, so its before == after by
    # construction; recording B (it sorts last) made the gate tautological.  Every
    # candidate is accumulated and a fail can never be overwritten by a later pass.
    c5_measured, c5_ok, c5_tested = {}, True, 0
    for name, man in sorted(manifests.items()):
        before = man.get("mean_pairwise_tip_distance_before_prune")
        after = man.get("mean_pairwise_tip_distance_after_prune")
        if not (before and after):
            continue
        ratio = after / before
        pruned = (str(man.get("recipe_arm") or "").upper() == "A")
        rec = {"before": before, "after": after, "ratio": ratio,
               "recipe_arm": man.get("recipe_arm"),
               "gated": pruned,
               "note": ("pruned stock tree: the real test" if pruned else
                        "built de novo, never pruned -- before == after by "
                        "construction, so it cannot falsify the gate and is not "
                        "gated on")}
        if pruned:
            c5_tested += 1
            if abs(ratio - 1.0) >= 1e-6:
                c5_ok = False
        c5_measured[name] = rec
    if c5_measured:
        out["C5"] = control(c5_ok and c5_tested > 0,
                            {"per_reference": c5_measured,
                             "n_pruned_references_tested": c5_tested,
                             "gated_arm": "A"},
                            "abs(pruned/full - 1) < 1e-6 over a 200-tip sample, on "
                            "the PRUNED arm (recipe_arm A) and on at least one such "
                            "reference",
                            "ete3's default pruning yields 0.537x and would "
                            "fabricate the whole NSTI/novelty relationship; a C5 "
                            "recorded from the unpruned arm B cannot fail at all")
    # ---- C11/C12 phenotype controls, bound to THIS run --------------------- #
    # Computed in main() (it needs the cells fingerprint and the discovered runs) and
    # handed in through extras; a present-but-unbindable summary FAILS here instead
    # of disappearing.
    out.update(phenotype_controls(extras.get("phenotype_binding") or {}))
    # ---- C14 prediction-run completeness (HARD) ---------------------------- #
    # Three independent legs, any one of which fails the gate:
    #  (1) every launched run's table is self-consistent with its own meta counts;
    #  (2) every md5 the fold plan queried (minus that run's own min_align drops,
    #      domain-filtered for a forced-routing PICRUSt2 run) is in the table;
    #  (3) zero cells were scored against a missing record -- no cells.csv.gz row
    #      carries a MISSING_RECORD_REASONS sentinel in abstain_reason;
    #  plus (4) the printed population F1s subtract to the printed paired delta,
    #      i.e. every arm was scored over the same denominator.
    comp = extras.get("run_completeness") or {}
    cells_comp = (summary.get("accounting", {}) or {}).get("completeness") or {}
    notes = (summary.get("meta", {}) or {}).get("notes", {}) or {}
    n_sentinel = int(cells_comp.get("n_missing_record_rows") or 0)
    denom_ok = notes.get("denominators_consistent", True)
    uniform = cells_comp.get("unscoreable_uniform_across_arms", True)
    plan_cov = cells_comp.get("plan_coverage") or {}
    inv = cells_comp.get("inventory") or extras.get("inventory") or {}
    arm_den = cells_comp.get("arm_denominators") or {}
    if comp or cells_comp:
        ok = (comp.get("ok", True) and n_sentinel == 0 and bool(denom_ok) and
              uniform is not False and plan_cov.get("ok") is not False and
              inv.get("ok") is not False and arm_den.get("ok") is not False and
              not (inv.get("gated") and notes.get("picrust2_domains_not_run")
                   and not args.expect_domains))
        out["C14"] = control(
            ok,
            {"runs": comp.get("per_run", {}),
             "plan_coverage": plan_cov,
             # leg 6: the derived-vs-discovered run inventory.  The fingerprints
             # live in the summary (accounting.completeness.inventory), not here.
             "inventory": {k: v for k, v in sorted(inv.items())
                           if k != "fingerprints"},
             # leg 7: every arm over the same denominator
             "arm_denominators": arm_den,
             "n_runs_incomplete": comp.get("n_runs_incomplete", 0),
             "n_cells_scored_without_a_record": n_sentinel,
             "missing_record_rows_by_arm":
                 cells_comp.get("missing_record_rows_by_arm", {}),
             "unscoreable_uniform_across_arms": uniform,
             "denominators_consistent": denom_ok,
             "denominator_check": notes.get("denominator_check", {}),
             "picrust2_domains_not_run": notes.get("picrust2_domains_not_run", {})},
            "for every launched run: len(pred) == n_placed == n_query - "
            "n_min_align_dropped (PICRUSt2) / one record per query for exactly the "
            "meta.select_modes mode set (prFBA); set(fold_plan.test_md5s) - "
            "min_align_dropped is covered; zero cells scored with abstain_reason in "
            "%s and outcome != unscoreable_min_align; set(folds scored) x "
            "set(launched densities) == set(fold_plan folds) x set(launched "
            "densities); population F1(primary) - F1(picrust2) == "
            "delta_f1_population; the DISCOVERED run directories equal the "
            "inventory DERIVED from the fold plan, the frozen config and this "
            "invocation's declared scope EXACTLY (no missing key, no duplicate or "
            "unexpected directory, every run's own declared scope -- prFBA "
            "select_modes equal across folds, PICRUSt2 domain/routing -- matching "
            "the plan); every arm's population.f1.n equal within a space; and no "
            "PICRUSt2 (fold, domain) left unlaunched unless --expect-domains "
            "declared it" % (list(MISSING_RECORD_REASONS),),
            "a truncated, resumed-stale or plan-mismatched prediction run is "
            "otherwise laundered into method abstentions: the md5 is simply absent, "
            "every cell of it scores 0.0, and the arm reads out at lower COVERAGE -- "
            "the quantity the headline rests on. Measured: deleting 60 of 561 md5s "
            "from one fold's bacterial run moved the primary delta -0.2610 -> "
            "-0.2347 with every other control still passing. A domain whose run was "
            "never LAUNCHED is treated differently PER CELL -- those cells are "
            "unscoreable_method_not_run for every arm, so the denominators stay "
            "equal -- but in the gated modes it is now a leg-6/7 FAILURE rather "
            "than a note, because equal denominators over a 10.4%-smaller "
            "population behind a headline that never mentions Archaea is precisely "
            "the laundering this control exists to stop (pass --expect-domains to "
            "declare a deliberately narrower scope, SPEC D22). An md5 dropped "
            "by --min_align 0.80 is excluded from both denominators (SPEC 7-C7), so "
            "its sentinel row is not counted by leg 3. A planned fold absent from "
            "BOTH methods is caught by leg 5: both denominators shrink equally, so "
            "legs 1-4 are silent while the headline covers 2 of 3 folds. Legs 6-7 "
            "invert the question, because legs 1-5 all ask it of the runs that "
            "happen to be ON DISK: a killed (fold, domain) job whose refs/ survives "
            "(-30% of the population at EXIT 0), a stale retry directory merged in "
            "silently (primary delta -0.2610 -> -0.0345 for PICRUSt2, -0.2039 for "
            "prFBA), a never-launched archaeal arm (-10.4%, delta -0.2443 against a "
            "frozen SESOI of 0.02), a fold launched with a shorter --select-modes "
            "(that arm at n=1258 vs 1922) and a whole lost density tier are all one "
            "defect -- the scored inventory is not the planned one -- and set "
            "equality against the DERIVED plan catches them by construction rather "
            "than by whichever symptom each happens to leave behind")
    # ---- C15 trait-column screen (reported) -------------------------------- #
    screen = (summary.get("accounting", {}) or {}).get("trait_column_screen") or {}
    if screen:
        retry = sorted(n for n, v in screen.items() if v.get("retry_used"))
        not_repairable = {n: {k: v["by_reason"][k] for k in ("singleton", "crash")
                              if k in v.get("by_reason", {})}
                          for n, v in sorted(screen.items())
                          if set(v.get("by_reason", {})) & {"singleton", "crash"}}
        out["C15"] = control(True,
                             {"per_reference": screen,
                              "references_using_aggressive_retry": retry,
                              "columns_dropped_but_not_repairable": not_repairable,
                              "reinstated_by_run": (summary.get("accounting", {})
                                                    or {}).get(
                                  "full_prevalence_ec_reinstated_by_run", {})},
                             "reported, not gated",
                             "a prevalence-n ('full') column IS reinstated as a "
                             "PICRUSt2 prediction (castor would reconstruct it "
                             "present everywhere); a prevalence-0 ('zero') column is "
                             "genuinely unpredictable; but a column dropped by the "
                             "SPEC 5.4 aggressive retry ('singleton'/'crash') is not "
                             "reconstructable in either direction, so a fold that "
                             "used the retry carries an undisclosed PICRUSt2 recall "
                             "penalty unless it is named here")
    return out


# --------------------------------------------------------------------------- #
# recipe control (C1)
# --------------------------------------------------------------------------- #
def run_recipe_control(args):
    """Score the recipe arms against the stock truth and settle C1.

    Arm A = stock RAxML tree pruned with preserve_branch_length + stock .model on
    the stock MSA; arm B = our hmmbuild/FastTree/.model recipe on that SAME stock
    MSA; arm C = the same retained sequences degapped and run through the full
    §5.1 recipe (hmmalign + RF match columns + hmmbuild + FastTree + derived
    model), i.e. what every fold reference actually is.  Traits, truth and queries
    are identical across arms, so the paired differences decompose the build
    recipe: C-A is the whole dialect (GATED), B-A the tree+model half with the MSA
    held fixed, C-B the alignment half with the tree recipe held fixed, and per
    query (C-A) == (C-B) + (B-A) exactly.  The arm labels extend the cell table's
    ``arm`` vocabulary, which is why this table is written to its own file
    (``cells_recipe.csv.gz``).
    """
    fps = _recipe_fingerprints(args)
    truth_raw = jload(args.truth)
    ec_index = {}
    truth = {}
    for label, ecs in sorted(truth_raw.items()):
        s = set()
        for e in ecs:
            e = norm_ec(e)
            if e is None:
                continue
            s.add(ec_index.setdefault(e, len(ec_index)))
        truth[str(label)] = frozenset(s)
    # the recipe control scores in the STOCK trait vocabulary, so a predicted EC the
    # stock truth happens never to carry still counts as a false positive (extend=True)
    runs = discover_picrust2(args.picrust2, ec_index, extend=True)
    runs = [r for r in runs if "recipe" in r["name"].lower()]
    if not runs:
        raise SystemExit("[score_folds] --mode recipe_control found no recipe_* run "
                         "under %s" % args.picrust2)
    universe = frozenset(ec_index.values())
    os.makedirs(os.path.dirname(os.path.abspath(args.cells)) or ".", exist_ok=True)
    per_arm = {}
    fh, handles = gzip_text_writer(args.cells)
    try:
        w = csv.DictWriter(fh, fieldnames=CELL_FIELDS)
        w.writeheader()
        for r in runs:
            arm = "picrust2_" + re.sub(r"[^A-Za-z0-9_]", "_", r["name"])
            scores = {}
            for label in sorted(truth):
                t = truth[label]
                rec = r["pred"].get(label)
                p = rec["ec"] if rec else frozenset()
                row = {k: "" for k in CELL_FIELDS}
                row.update({"amplicon": label, "md5": label, "region": "V4",
                            "domain": "Bacteria", "src_genome_id": label, "fold": "",
                            "group_id": label, "density": "", "arm": arm, "space": "ec",
                            "n_truth": len(t), "n_pred": len(p),
                            "nsti": fmt(rec.get("nsti") if rec else None, 6),
                            "routing": "forced", "abstained": 0,
                            "min_align_dropped": 0, "novelty_bin": "",
                            "n_truth_impossible": "", "unscoreable_reason": ""})
                if not t:
                    row["outcome"] = "unscoreable_empty_truth"
                    row["unscoreable_reason"] = "truth_empty_in_universe"
                elif not p:
                    row["outcome"] = "pred_empty"
                else:
                    row["outcome"] = "scored"
                    tp, rc_, pr_, f1, j = metrics(t, p)
                    row.update({"tp": tp, "recall": fmt(rc_), "precision": fmt(pr_),
                                "f1": fmt(f1), "jaccard": fmt(j),
                                "mcc": fmt(mcc_from_counts(tp, len(p) - tp,
                                                           len(t) - tp, len(universe)))})
                    scores[label] = f1 or 0.0
                w.writerow(row)
            per_arm[arm] = scores
    finally:
        fh.close()
        for h in handles:
            h.close()
    arms = sorted(per_arm)
    by_letter = {}
    for letter in ("a", "b", "c"):
        hit = next((a for a in arms if a.lower().endswith("_" + letter)), None)
        if hit is not None:
            by_letter[letter.upper()] = hit
    arm_a = by_letter.get("A", arms[0])
    arm_b = by_letter.get("B", arms[-1])
    arm_c = by_letter.get("C")
    # An arm C reference that was built but never run would silently drop C1 back
    # to the tree-only contrast it used to be -- the exact defect this covers.
    if arm_c is None:
        c_manifest = os.path.join(
            os.path.dirname(os.path.abspath(args.truth)), "recipe_C",
            "manifest.json")
        if os.path.exists(c_manifest):
            raise SystemExit(
                "[score_folds] %s exists but no recipe_C run was found under %s: "
                "C1 would silently fall back to the tree+model-only contrast. Run "
                "run_picrust2_fold.py --mode recipe_control on refs/recipe_C "
                "first." % (c_manifest, args.picrust2))

    def paired(a_arm, b_arm):
        """mean / Wilcoxon of F1(b) - F1(a) over the queries both arms scored."""
        keys = sorted(set(per_arm[a_arm]) & set(per_arm[b_arm]))
        d = [per_arm[b_arm][k] - per_arm[a_arm][k] for k in keys]
        mean_d = st.mean(d) if d else None
        p_w = None
        if len(d) >= 10 and any(x != 0 for x in d):
            try:
                from scipy.stats import wilcoxon
                p_w = float(wilcoxon([per_arm[b_arm][k] for k in keys],
                                     [per_arm[a_arm][k] for k in keys],
                                     zero_method="zsplit")[1])
            except Exception:
                p_w = None
        return keys, mean_d, p_w

    contrasts = {}
    pairs = [(arm_a, arm_b, "tree_model_only_msa_held_fixed")]
    if arm_c:
        pairs = [(arm_a, arm_c, "full_build_recipe_alignment_tree_model"),
                 (arm_a, arm_b, "tree_model_only_msa_held_fixed"),
                 (arm_b, arm_c, "alignment_only_tree_recipe_held_fixed")]
    for a_arm, b_arm, covers in pairs:
        keys, mean_d, p_w = paired(a_arm, b_arm)
        contrasts["%s__vs__%s" % (b_arm, a_arm)] = {"ec": {
            "mean_paired_delta_f1": rnd(mean_d), "wilcoxon_p": rnd(p_w, 6),
            "n_pairs": len(keys), "covers": covers,
            "arm_b_mean_f1": rnd(st.mean([per_arm[b_arm][k] for k in keys]))
            if keys else None,
            "arm_a_mean_f1": rnd(st.mean([per_arm[a_arm][k] for k in keys]))
            if keys else None}}
    # the gated contrast: the whole recipe when arm C exists, else tree+model only
    gate_a, gate_b = (arm_a, arm_c) if arm_c else (arm_a, arm_b)
    common, mean_d, p_w = paired(gate_a, gate_b)
    coverage = ("alignment+tree+model" if arm_c else "tree_model_only")
    summary = {
        "meta": {"frozen": FROZEN, "mode": "recipe_control", "generated": now_iso(),
                 "git_commit": git_commit(), "n_queries": len(truth),
                 "n_paired": len(common), "arms": arms,
                 "arm_a": arm_a, "arm_b": arm_b, "arm_c": arm_c,
                 "gated_contrast": "%s__vs__%s" % (gate_b, gate_a),
                 "gated_contrast_covers": coverage,
                 "ec_universe": len(universe),
                 "note": "arm labels extend the cells.csv.gz arm vocabulary; this "
                         "table is written separately for that reason",
                 "cells": os.path.abspath(args.cells)},
        "arms": {a: {"ec": {"mean_f1": rnd(st.mean(list(v.values()))) if v else None,
                            "n": len(v)}} for a, v in sorted(per_arm.items())},
        "contrasts": contrasts,
        # provenance, same as the fold path: size + sha256 of every prediction
        # table and manifest these numbers were computed from, so a later reader
        # can prove which bytes produced them (a concurrent session rewriting a
        # predictions.json mid-review moved picrust2 0.6955 -> 0.6526 with every
        # control green, and nothing in the artefact could show it)
        "accounting": {"fingerprints": fps},
    }
    jdump(args.summary, summary)
    gate_letter = "C" if arm_c else "B"
    measured = {"gated_contrast": "%s__vs__%s" % (gate_b, gate_a),
                "gated_contrast_covers": coverage,
                "mean_paired_delta_f1_%s_minus_A" % gate_letter: rnd(mean_d),
                "wilcoxon_p": rnd(p_w, 6), "n_pairs": len(common),
                "arm_a": arm_a, "arm_b": arm_b, "arm_c": arm_c,
                "decomposition": {k: v["ec"] for k, v in contrasts.items()},
                "not_covered": [
                    "ssu-align's upstream curation of the raw stock sequences "
                    "(only the aligned bac_ref.fna is distributed, so arm C "
                    "re-aligns the degapped stock residues)",
                    "count-vs-binary traits: all recipe arms use the stock EC "
                    "counts, while the fold arm uses binary traits"],
                "audit_reference_n40": {"paired_dF1_B_minus_A": -0.0016,
                                        "wilcoxon_p": 0.6368}}
    ctl = {"C1": control(mean_d is not None and abs(mean_d) < 0.01, measured,
                         "abs(mean paired dF1(%s-A)) < 0.01" % gate_letter,
                         "dF1 <= -0.01 means our full hmmalign/RF/hmmbuild/FastTree "
                         "recipe handicaps PICRUSt2 relative to its own stock "
                         "ssu-align MSA + RAxML tree + stock model, and every "
                         "PICRUSt2 number must carry that band as an uncertainty; "
                         "the B-A and C-B rows of the decomposition attribute it to "
                         "the tree+model or the alignment half")}
    merge_controls(args.controls, ctl, mode="recipe_control",
                   provenance={"mode": "recipe_control",
                               "evaluated_at": now_iso(),
                               "cells": os.path.abspath(args.cells),
                               "cells_sha256": (file_fp(args.cells) or {}).get(
                                   "sha256"),
                               "inventory_fp": None,
                               "answer_key_fp": fps.get("answer_key_fp"),
                               "scorer_code_fp": fps.get("scorer_code_fp"),
                               "git_commit": git_commit()})
    print("[score_folds] C1 recipe calibration [%s]: mean paired dF1(%s-A) = %s "
          "(wilcoxon p=%s, n=%d) -> %s" % (coverage, gate_letter, rnd(mean_d),
                                           rnd(p_w, 4), len(common),
                                           ctl["C1"]["status"]))
    for k, v in contrasts.items():
        print("[score_folds]   %-44s dF1=%-10s p=%-8s n=%-4s (%s)"
              % (k, rnd(v["ec"]["mean_paired_delta_f1"]),
                 rnd(v["ec"]["wilcoxon_p"], 4), v["ec"]["n_pairs"],
                 v["ec"]["covers"]))
    print("[score_folds] -> %s, %s, %s" % (args.cells, args.summary, args.controls))
    # C1 is HARD, so it must move the exit status too (same rule as the fold modes,
    # SPEC 7.1): a recipe calibration that failed used to exit 0 and be visible only
    # inside controls.json.
    if ctl["C1"]["status"] == "fail":
        print("=" * 78)
        print("[score_folds] *** HARD CONTROL C1 FAILED: mean paired dF1(%s-A) = %s "
              "is outside the +-0.01 band -- our build recipe handicaps PICRUSt2 "
              "relative to its own stock build. EXIT 3. ***"
              % (gate_letter, rnd(mean_d)))
        print("=" * 78)
        return 3
    return 0


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    D = os.path.join(SR, "data")
    ap.add_argument("--mode", required=True,
                    choices=["fold", "recipe_control", "include_self", "density",
                             "plumbing"])
    ap.add_argument("--fold-plan", default=os.path.join(D, "fold_plan.json"))
    ap.add_argument("--panel", default=os.path.join(D, "panel_eligible.json"))
    ap.add_argument("--picrust2", required=True)
    ap.add_argument("--prfba", default=None)
    ap.add_argument("--truth", default=C.PATHS.truth_json)
    ap.add_argument("--functions", default=os.path.join(P2C, "data",
                                                        "genome_function_sets.json"))
    ap.add_argument("--maps", default=os.path.join(P2C, "data", "function_maps.json"))
    ap.add_argument("--neighbors", default=os.path.join(D, "neighbors.json"))
    ap.add_argument("--refs", required=True)
    ap.add_argument("--cells", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--controls", required=True)
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--rank-blocked", default="genus")
    ap.add_argument("--space", default="ec,rxn")
    # --- additions (optional; the spec's CLI works verbatim without them) --- #
    ap.add_argument("--ec-universe", default=None,
                    help="panel_ec_universe.json (default: beside --panel)")
    ap.add_argument("--fanout", default=C.PATHS.amplicon_fanout)
    ap.add_argument("--annotation-ceiling",
                    default=os.path.join(P2C, "data", "annotation_ceiling.json"))
    ap.add_argument("--phenotype-summary", default=None,
                    help="phenotype_summary.json for C11/C12 (default: beside --cells)")
    ap.add_argument("--templates", default=TEMPLATE_DIR,
                    help="ModelSEED template dir for the template-visible fraction")
    ap.add_argument("--no-templates", action="store_true")
    ap.add_argument("--primary-select-mode", default="exact_tie")
    ap.add_argument("--clade-max-cap", type=int, default=3000)
    ap.add_argument("--limit-cells", type=int, default=0,
                    help="smoke only: score at most N cells per fold")
    ap.add_argument("--force", action="store_true",
                    help="rescore even when --cells already exists")
    ap.add_argument("--require-phenotype-binding", action="store_true",
                    help="fail C11/C12 when phenotype_summary.json carries no "
                         "meta.binding block tying it to this cells table, fold "
                         "plan and answer key (SPEC 11.2 passes this; a summary "
                         "written before the binding contract has none)")
    ap.add_argument("--allow-underpowered-bootstrap", action="store_true",
                    help="run with --bootstrap below %d (or 0) in a gated mode; the "
                         "primary claim is then stamped indeterminate/underpowered "
                         "and the run is NOT valid" % MIN_BOOTSTRAP)
    # --- declared scope of THIS invocation (the inventory gate, C14 leg 6) --- #
    # The expected run set is derived from the fold plan and the frozen config; the
    # axes below are the ones that legitimately depend on an operator choice, so
    # they are DECLARED here rather than inferred from the directories on disk.
    ap.add_argument("--expect-folds", default=None,
                    help="comma list of fold indices this invocation covers "
                         "(default: every fold in the fold plan; declare a subset "
                         "for a deliberately partial smoke run)")
    ap.add_argument("--expect-densities", default=None,
                    help="comma list of density fractions this invocation covers "
                         "(default: frozen.density_fractions in --mode density; "
                         "frozen.density_fractions INTERSECT the reference builds "
                         "in --refs in --mode fold)")
    ap.add_argument("--expect-domains", default=None,
                    help="comma list of domains every fold's PICRUSt2 arm covers "
                         "(default: the domains of the plan's own test cells; pass "
                         "e.g. 'bac' for a deliberately bacteria-only run, SPEC D22)")
    ap.add_argument("--expect-nsti", default="none",
                    choices=["none", "primary", "all"],
                    help="how many --routing nsti PICRUSt2 runs the campaign "
                         "launched: none, one per fold at the primary density "
                         "(SPEC 11.2's 10 runs), or one per fold x density")
    ap.add_argument("--expect-forced-answer", default="primary",
                    choices=["none", "primary", "all"],
                    help="how many prFBA --forced-answer runs the campaign "
                         "launched: none, one per fold at the primary density "
                         "(SPEC 11.2's 10 runs), or one per fold x density")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    t0 = time.time()
    if args.mode == "recipe_control":
        return run_recipe_control(args)
    # A gated mode publishes the primary claim, and the claim is the CI: refuse a
    # bootstrap too small to support one rather than printing "SESOI not exceeded"
    # over a missing interval.
    if args.mode in ("fold", "density") and args.bootstrap < MIN_BOOTSTRAP \
            and not args.allow_underpowered_bootstrap:
        raise SystemExit(
            "[score_folds] REFUSING --bootstrap %d in --mode %s: the primary claim "
            "is a cluster-bootstrap CI over the holdout groups, and fewer than %d "
            "resamples cannot support one (--bootstrap 0 computes NO interval at "
            "all, and the headline then reads 'SESOI not exceeded' -- which is NOT "
            "what 'no interval was computed' means). Use --bootstrap %d or more "
            "(SPEC 11.2 uses 2000), or pass --allow-underpowered-bootstrap to run "
            "anyway with the claim stamped indeterminate and the run stamped "
            "invalid." % (args.bootstrap, args.mode, MIN_BOOTSTRAP, MIN_BOOTSTRAP))
    if not args.prfba:
        print("[score_folds] no --prfba given: scoring PICRUSt2 and the baselines only")

    args._pheno_explicit = bool(args.phenotype_summary)
    if args._pheno_explicit and not os.path.exists(args.phenotype_summary):
        # Silently treating it as "the phenotype arm did not run" is the same
        # silent-loss shape as the vanishing C11/C12: the operator NAMED a file.
        raise SystemExit("[score_folds] --phenotype-summary %s does not exist. An "
                         "explicitly named phenotype summary that is absent is an "
                         "error, not an arm that did not run (C11/C12 would simply "
                         "vanish from %s, which reads exactly like a clean run)."
                         % (args.phenotype_summary, args.controls))
    if not args.phenotype_summary:
        cand = os.path.join(os.path.dirname(os.path.abspath(args.cells)),
                            "phenotype_summary.json")
        args.phenotype_summary = cand if os.path.exists(cand) else None

    plan = jload(args.fold_plan)
    panel, panel_meta = load_panel(args.panel)
    frozen = (plan.get("meta", {}) or {}).get("frozen") or FROZEN
    ec_path = args.ec_universe or os.path.join(os.path.dirname(
        os.path.abspath(args.panel)), "panel_ec_universe.json")
    ec_list, ec_index, _prev, _ecmeta = load_ec_universe(ec_path)
    ec_rxn, rxn_index = load_ec_rxn(args.maps, ec_index)
    CTX_EC_RXN[0] = ec_rxn
    universe = {"ec": frozenset(ec_index.values()),
                "rxn": frozenset(r for e in ec_index.values()
                                 for r in ec_rxn.get(e, ()))}
    print("[score_folds] universe: %d EC, %d reactions; panel %d genomes"
          % (len(universe["ec"]), len(universe["rxn"]), len(panel)))

    truth_meta = jload(args.truth)
    md5_of_cell = load_fanout(args.fanout)
    neighbors = jload(args.neighbors) if os.path.exists(args.neighbors) else {}
    manifests = load_ref_manifests(args.refs)
    p2_runs = discover_picrust2(args.picrust2, ec_index)
    pf_runs = discover_prfba(args.prfba, ec_index)
    reinstated = reinstate_full_prevalence_ec(p2_runs, manifests, ec_index)
    if any(reinstated.values()):
        print("[score_folds] reinstated prevalence-n ('full') trait columns as "
              "PICRUSt2 predictions: %s"
              % {k: v for k, v in sorted(reinstated.items()) if v})
    p2_all = list(p2_runs)
    p2_runs, pf_runs = select_runs(args.mode, p2_runs, pf_runs)
    print("[score_folds] %d PICRUSt2 run(s), %d prFBA run(s), %d reference manifest(s)"
          % (len(p2_runs), len(pf_runs), len(manifests)))
    if not p2_runs and not pf_runs:
        raise SystemExit("[score_folds] no predictions found for --mode %s" % args.mode)

    # ---- Stage 0: the run INVENTORY gate, BEFORE any scoring -------------- #
    # Derived from the fold plan + the frozen config + this invocation's declared
    # scope; asserted against the discovered directories; every input fingerprinted.
    inventory = run_inventory(args, plan, panel, truth_meta, manifests, p2_runs,
                              pf_runs, frozen, ec_path=ec_path)
    print("[score_folds] inventory: %d expected / %d discovered run key(s) in %d "
          "directory(ies) [%s]; fp %s"
          % (inventory["n_expected"], inventory["n_discovered"],
             inventory["n_run_directories"],
             "GATED" if inventory["gated"] else "reported only",
             inventory["inventory_fp"][:12]))
    der = inventory["derivation"]
    print("[score_folds] inventory derived from: folds %s (%s); densities %s (%s); "
          "domains %s; nsti=%s; forced_answer=%s"
          % (der["folds_expected"], der["folds_source"], der["densities_expected"],
             der["densities_source"], der["domains_source"],
             der["nsti_runs_expected"], der["forced_answer_runs_expected"]))
    if der["densities_in_frozen_plan_not_built"]:
        print("[score_folds] NOTE: %s of the frozen density_fractions %s have NO "
              "reference build under --refs, so the density sweep this invocation "
              "can see is %s. Pass --expect-densities to gate the full frozen set "
              "(--mode density gates it by default)."
              % (der["densities_in_frozen_plan_not_built"],
                 der["densities_in_frozen_plan"], der["densities_expected"]))
    # The ANSWER KEY, by name and by bytes: a reader of the artefact (and an
    # operator tailing the log) must be able to see WHICH truth file produced these
    # numbers.  A tampered genome_function_sets.json republished the shipped headline
    # verbatim through the reuse path and moved every number with --force.
    _akp = inventory.get("input_paths") or {}
    _akf = (inventory.get("fingerprints") or {}).get("answer_key") or {}
    print("[score_folds] answer key fp %s / scorer code fp %s"
          % ((inventory.get("answer_key_fp") or "")[:12],
             (inventory.get("scorer_code_fp") or "")[:12]))
    for _tag in sorted(_akp):
        _fp = _akf.get(_tag) or {}
        if _tag == "templates":
            _fp = {"size": sum((_akf.get("templates/%s" % r) or {}).get("size") or 0
                               for r in TEMPLATE_FILES),
                   "sha256": "(per role)"}
        print("[score_folds]   %-18s %s%s" % (_tag, _akp[_tag] or "(none)",
                                              (" [%s %s]" % (_fp.get("size"),
                                                             str(_fp.get("sha256"))[:12])
                                               if _fp else " [MISSING]")))
    _missing_ak = [t for t in ANSWER_KEY_INPUTS
                   if _akp.get(t) and not (_akf.get(t) or
                                           any(k.startswith(t + "/") for k in _akf))]
    if _missing_ak:
        print("[score_folds] NOTE: answer-key input(s) %s are declared but absent "
              "from disk, so nothing fingerprints them" % _missing_ak)
    inventory_banner(inventory)

    needed = set(panel)
    for f in plan.get("folds", []):
        needed |= set(f.get("reference_genomes") or ())
        needed |= set(f.get("test_genomes") or ())
    g_ec = load_functions(args.functions, sorted(needed), ec_index)
    g_rxn = {g: project(v, ec_rxn) for g, v in g_ec.items()}
    tv_sets = {} if args.no_templates else load_template_visible(args.templates,
                                                                 ec_index)
    if tv_sets:
        print("[score_folds] template-visible EC sets: %s"
              % {k: len(v) for k, v in sorted(tv_sets.items())})

    os.makedirs(os.path.dirname(os.path.abspath(args.cells)) or ".", exist_ok=True)
    bookkeeping = None
    # Reuse of an existing cells table without --force is now conditional on the
    # fingerprints matching: a stale table from a DIFFERENT inventory published
    # n=1722 / delta -0.2443 as valid while the completeness block it was shipped
    # with listed every run as present.  No sidecar, different cells bytes or a
    # different inventory -> recompute; if the destination is not writable, refuse.
    reuse_ok, reuse_why = (False, "--force")
    if not args.force and os.path.exists(args.cells):
        reuse_ok, reuse_why = cells_reuse_verdict(args, inventory)
        if not reuse_ok:
            cells_dir = os.path.dirname(os.path.abspath(args.cells)) or "."
            if not os.access(cells_dir, os.W_OK):
                raise SystemExit(
                    "[score_folds] REFUSING to reuse %s: %s -- and %s is not "
                    "writable, so it cannot be recomputed here. Rescore with "
                    "--force into a writable directory (or make this one writable); "
                    "scoring a cells table that cannot be tied to the runs on disk "
                    "is how a stale population (n=1722 vs the planned 1922) was "
                    "published as valid." % (args.cells, reuse_why, cells_dir))
            print("[score_folds] NOT reusing %s: %s -- recomputing"
                  % (args.cells, reuse_why))
    if args.force or not os.path.exists(args.cells) or not reuse_ok:
        tmp = args.cells + ".tmp"
        fh, handles = gzip_text_writer(tmp)
        try:
            w = csv.DictWriter(fh, fieldnames=CELL_FIELDS)
            w.writeheader()
            bookkeeping = score_cells(
                args, (plan, panel, g_ec, g_rxn, ec_rxn, ec_index, truth_meta,
                       md5_of_cell, neighbors, manifests, p2_runs, pf_runs, tv_sets,
                       universe), w)
        finally:
            fh.close()
            for h in handles:
                h.close()
        os.replace(tmp, args.cells)
        print("[score_folds] wrote %d cell rows -> %s"
              % (bookkeeping["n_rows"], args.cells))
        fp_rec = {"cells": file_fp(args.cells),
                  "inventory_fp": inventory["inventory_fp"],
                  # the answer key and the CODE this table was scored with: the two
                  # identities a later no---force rescore must match before it may
                  # republish these numbers
                  "answer_key_fp": inventory.get("answer_key_fp"),
                  "answer_key_paths": inventory.get("input_paths"),
                  "scorer_code_fp": inventory.get("scorer_code_fp"),
                  "scorer_code_files": inventory.get("scorer_code_files"),
                  "mode": args.mode, "generated": now_iso(),
                  "script": "score_folds.py", "git_commit": git_commit(),
                  "fingerprints": inventory["fingerprints"]}
        try:
            jdump(cells_fp_path(args.cells), fp_rec)
        except OSError as e:
            print("[score_folds] WARNING: could not write the fingerprint sidecar "
                  "%s (%s): a later no---force rescore will refuse to reuse this "
                  "cells table rather than trust it"
                  % (cells_fp_path(args.cells), e))
    else:
        print("[score_folds] reusing existing %s: %s (use --force to rescore)"
              % (args.cells, reuse_why))

    # Every fingerprinted input is re-read now: the bytes the numbers rest on must
    # be the bytes this invocation started with (C14 leg 6).
    reverify_fingerprints(args, inventory, ec_path=ec_path)

    run_meta = {"frozen": frozen, "git_commit": git_commit()}
    # C14 legs 1-2 first, so the summary itself can be stamped: are the prediction
    # tables we are about to score actually complete?
    run_comp = run_completeness(plan, panel, truth_meta, md5_of_cell, p2_runs,
                                pf_runs)
    summary, extras = aggregate(args, args.cells, plan, panel, manifests, p2_runs,
                                pf_runs, run_meta, run_comp=run_comp,
                                inventory=inventory)
    # include-self source recovery (C2): did prFBA select the source genome itself?
    if args.mode == "include_self":
        rec = include_self_recovery(args, plan, pf_runs, md5_of_cell, g_ec)
        extras.update(rec)
    extras["run_completeness"] = run_comp
    extras["inventory"] = inventory
    # The phenotype arm is BOUND to this run (or declared unbound) before its
    # controls are evaluated: it was auto-discovered by existence and trusted, so
    # another campaign's summary scored C11/C12 pass and a truncated one made both
    # controls vanish inside an `except Exception: pass`.
    cells_fp_now = file_fp(args.cells)
    pheno_bind = phenotype_binding(args, plan, frozen, inventory, cells_fp_now)
    extras["phenotype_binding"] = pheno_bind
    summary.setdefault("accounting", {})["phenotype_binding"] = pheno_bind
    if pheno_bind.get("present"):
        print("[score_folds] phenotype arm: %s (%s) -- bound=%s%s"
              % (pheno_bind["path"], pheno_bind["discovered"],
                 pheno_bind.get("bound"),
                 ("; defects %s" % [d["defect"] for d in pheno_bind["defects"]])
                 if pheno_bind.get("defects") else ""))
        if not pheno_bind.get("bound") and not pheno_bind.get("defects"):
            print("[score_folds] WARNING: %s" % pheno_bind.get("unbound_reason"))
    jdump(args.summary, summary)                    # artefacts first, always
    ctl = fold_controls(args, plan, panel, manifests, p2_runs, pf_runs, neighbors,
                        summary, extras, p2_all=p2_all)
    # this invocation's own verdict, before the failure-preserving merge folds in
    # another mode's older fail: the exit status must describe THIS run
    fresh = {k: dict(v) for k, v in ctl.items()}
    hard_failed_pre = hard_control_failures(fresh)
    prov = control_provenance(args, inventory, cells_fp_now, args.mode)
    ctl = merge_controls(args.controls, ctl, mode=args.mode, provenance=prov)
    cprov = controls_provenance_report(ctl, set(fresh), prov, args.mode)
    summary["accounting"]["controls_provenance"] = cprov
    summary["headline"]["hard_controls_failed"] = list(hard_failed_pre)
    summary["headline"]["controls_evaluated_this_run"] = sorted(fresh)
    # VALIDITY reflects every hard control now, not just C14 + the inventory: a
    # C4/C5/C6/C9/C11 failure exited 3 while publishing "valid": true beside a skill
    # score that had flipped sign.
    _vreasons = ["hard_control_failed:%s" % c for c in hard_failed_pre]
    if cprov["missing_expected"]:
        _vreasons.append("controls_not_emitted:%s"
                         % ",".join(cprov["missing_expected"]))
    stamp_invalid(summary, _vreasons)
    jdump(args.summary, summary)

    print("[score_folds] %s" % summary["headline"]["sentence"])
    for arm in sorted(summary["arms"]):
        b = summary["arms"][arm].get("ec")
        if not b:
            continue
        pop = (b["population"]["f1"] or {}).get("value")
        cond = (b["conditional_matched"].get("f1") or {}).get("value")
        print("  %-28s n=%-6s coverage=%-6s popF1=%-7s condF1(matched n=%s)=%s"
              % (arm, (b["population"]["f1"] or {}).get("n"), b["coverage"], pop,
                 b["conditional_matched"].get("n_matched"), cond))
    _stale = {r["id"]: r for r in cprov["stale_same_mode"]}
    _other = {r["id"]: r for r in cprov["from_other_modes"]}
    for cid in sorted(ctl):
        tag = ""
        if cid not in fresh:
            if cid in _stale:
                tag = ("   [NOT re-evaluated by this --mode %s run; recorded %s]"
                       % (args.mode, _stale[cid].get("evaluated_at")))
            elif cid in _other:
                tag = ("   [from --mode %s @ %s%s]"
                       % (_other[cid].get("mode"), _other[cid].get("evaluated_at"),
                          "" if (_other[cid].get("scorer_code_matches") and
                                 _other[cid].get("answer_key_matches"))
                          else " -- ANSWER KEY OR CODE DIFFERS"))
            else:
                tag = "   [NOT evaluated by this run and carries NO provenance]"
        print("  control %-4s %s%s" % (cid, ctl[cid]["status"], tag))
    _suspect = ([r["id"] for r in cprov["from_other_modes"]
                 if not (r.get("scorer_code_matches") and
                         r.get("answer_key_matches"))] +
                [r["id"] for r in cprov["stale_same_mode"]])
    if cprov["without_provenance"] or _suspect:
        print("[score_folds] NOTE: control verdict(s) in %s that this invocation "
              "cannot vouch for: no provenance at all %s; recorded against a "
              "different answer key or a different scorer, or by this same mode and "
              "not re-evaluated now, %s. A pre-seeded 'pass' used to be "
              "indistinguishable from a fresh one."
              % (args.controls, cprov["without_provenance"], _suspect))
    if cprov["missing_expected"]:
        print("=" * 78)
        print("[score_folds] *** CONTROL(S) %s WERE NOT EMITTED by --mode %s, which "
              "always evaluates them. An absent control reads exactly like a passing "
              "one in %s, so the run is stamped invalid. ***"
              % (cprov["missing_expected"], args.mode, args.controls))
        print("=" * 78)
    print("[score_folds] -> %s, %s, %s (%.1fs)"
          % (args.cells, args.summary, args.controls, time.time() - t0))
    # The exit status reflects EVERY hard control of THIS invocation (a control is
    # hard unless its own gate text says "reported"): a C5/C9/C12 fail used to leave
    # EXIT=0, so an orchestrator that trusted the return code never saw it.
    # C14 keeps its own code (2) and its own banner, because its failure means THESE
    # NUMBERS ARE NOT THE NUMBERS -- an incomplete prediction table silently lowers
    # an arm's coverage, and coverage is what the headline rests on.  Any other hard
    # control failing exits 3.  Every artefact is written first: the operator needs
    # them to debug.
    hard_failed = hard_control_failures(fresh)
    if (fresh.get("C14") or {}).get("status") == "fail":
        m = fresh["C14"]["measured"]
        bad = {n: {k: v for k, v in r.items()
                   if k in ("mode", "fold", "domain", "n_query", "n_placed",
                            "n_rows", "n_expected_from_plan", "n_missing",
                            "missing_examples", "self_consistent", "plan_complete")}
               for n, r in (m.get("runs") or {}).items()
               if not r.get("self_consistent") or r.get("plan_complete") is False}
        print("=" * 78)
        print("[score_folds] *** C14 FAIL -- THE SCORED POPULATION IS NOT THE "
              "PLANNED ONE. The scores in %s are NOT valid: an md5 missing from a "
              "LAUNCHED run is scored as a method abstention (which lowers that "
              "arm's coverage and moves the headline), and a planned (fold, "
              "density) absent from BOTH methods shrinks the population itself. ***"
              % args.summary)
        print("  runs_incomplete        = %s" % (m.get("n_runs_incomplete"),))
        print("  cells_without_a_record = %s %s"
              % (m.get("n_cells_scored_without_a_record"),
                 m.get("missing_record_rows_by_arm") or ""))
        print("  denominators_consistent= %s" % (m.get("denominators_consistent"),))
        for n, r in sorted(bad.items()):
            print("  %-24s %s" % (n, r))
        print("  plan_coverage          = %s" % (m.get("plan_coverage") or {},))
        _inv = m.get("inventory") or {}
        print("  inventory_ok           = %s reasons=%s (expected %s / discovered "
              "%s key(s))" % (_inv.get("ok"), _inv.get("reasons"),
                              _inv.get("n_expected"), _inv.get("n_discovered")))
        for _k in (_inv.get("missing") or []):
            print("    MISSING   %s" % _k)
        for _d in (_inv.get("duplicates") or []):
            print("    DUPLICATE %s -> %s" % (_d.get("key"), _d.get("dirs")))
        for _d in (_inv.get("unexpected") or []):
            print("    UNEXPECTED %s -> %s" % (_d.get("key"), _d.get("dirs")))
        for _d in (_inv.get("scope_mismatch") or []):
            print("    SCOPE     %s %s" % (_d.get("defect"), _d.get("run")))
        print("  arm_denominators       = %s"
              % ((m.get("arm_denominators") or {}).get("distinct_n_by_space"),))
        print("  domains_not_run        = %s"
              % (m.get("picrust2_domains_not_run") or {},))
        print("  inventory_fp           = %s" % (_inv.get("inventory_fp"),))
        print("  Re-run the affected run(s), then rescore with --force.")
        if [c for c in hard_failed if c != "C14"]:
            print("  OTHER HARD CONTROL(S) ALSO FAILED: %s"
                  % ([c for c in hard_failed if c != "C14"],))
        print("=" * 78)
        return 2
    if hard_failed:
        print("=" * 78)
        print("[score_folds] *** HARD CONTROL(S) FAILED: %s -- see %s. The exit "
              "status is 3 (not 2: C14 passed, so the numbers are over the planned "
              "population, but a hard setup/leakage control did not hold). ***"
              % (hard_failed, args.controls))
        for cid in hard_failed:
            print("  %-4s gate: %s" % (cid, (fresh[cid].get("gate") or "")[:140]))
        print("=" * 78)
        return 3
    # C14 passed and every hard control held, but the artefact may still be stamped
    # invalid -- an indeterminate primary claim (no CI) or a control this mode must
    # emit missing.  An invalid artefact must not leave EXIT 0 behind it.
    if summary["headline"].get("valid") is False:
        print("=" * 78)
        print("[score_folds] *** THE ARTEFACT IS STAMPED INVALID: %s. EXIT 3. ***"
              % (summary["meta"].get("run_invalid_reasons"),))
        print("  headline.sentence = %s" % summary["headline"]["sentence"])
        print("=" * 78)
        return 3
    return 0


def include_self_recovery(args, plan, pf_runs, md5_of_cell, g_ec=None):
    """C2's hard gate: on twin-free cells, prFBA must select the source genome.

    A cell is twin-free when its amplicon md5 maps to exactly one panel genome; twin
    cells are reported separately against the analytic ``source union twin`` ceiling
    and are never gated (an exact-amplicon twin can legitimately be selected instead).
    SPEC 7-C2 requires that twin report, so it is computed here rather than letting
    the twin cells vanish: for each twin cell, prFBA's EC set against the union of
    the EC sets of every genome that md5 maps to (the analytic ceiling) and against
    the source genome alone.
    """
    fan = jload(args.fanout)
    twins = {md5: {c.rsplit("__", 1)[0] for c in cells} for md5, cells in fan.items()}
    g_ec = g_ec or {}
    hits = total = twin_cells = 0
    twin_f1_union, twin_f1_source, twin_recovery = [], [], 0
    for r in pf_runs:
        per = r["modes"].get(args.primary_select_mode) or next(
            iter(r["modes"].values()), {})
        for f in plan.get("folds", []):
            if r["fold"] not in (None, f["fold"]):
                continue
            for cell in sorted(f.get("test_cells") or ()):
                md5 = md5_of_cell.get(cell)
                rec = per.get(md5)
                if rec is None:
                    continue
                src = cell.rsplit("__", 1)[0]
                sibs = twins.get(md5, {src})
                if len(sibs) > 1:
                    twin_cells += 1
                    pred = set(rec.get("ec") or ())
                    ceil_set = set()
                    for g in sorted(sibs):
                        ceil_set |= set(g_ec.get(g, ()))
                    own = set(g_ec.get(src, ()))
                    if pred and ceil_set:
                        twin_f1_union.append(BASE.f1_set(ceil_set, pred))
                    if pred and own:
                        twin_f1_source.append(BASE.f1_set(own, pred))
                    if set(rec.get("reps") or ()) & sibs:
                        twin_recovery += 1
                    continue
                total += 1
                if src in (rec.get("reps") or ()):
                    hits += 1
    return {"include_self_source_recovery": (hits / total) if total else None,
            "include_self_n": total, "include_self_twin_cells": twin_cells,
            "include_self_twin": {
                "n_twin_cells": twin_cells,
                "n_selected_a_twin_or_the_source": twin_recovery,
                "mean_f1_vs_source_union_twin_ceiling":
                    rnd(st.mean(twin_f1_union)) if twin_f1_union else None,
                "mean_f1_vs_source_only":
                    rnd(st.mean(twin_f1_source)) if twin_f1_source else None,
                "gate": "reported, never gated (SPEC 7-C2)"}}


if __name__ == "__main__":
    sys.exit(main())
