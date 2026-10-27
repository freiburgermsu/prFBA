#!/usr/bin/env python
"""Phenotype arm of the shared-reference k-fold benchmark: does an EC-level difference
between prFBA and PICRUSt2 propagate to metabolic phenotype?

WHY THIS COMPONENT EXISTS
-------------------------
Every EC-level endpoint in this benchmark is an agreement score against one annotation
(BV-BRC's EC projection).  Two EC sets differing by 20% may yield identical flux
predictions or wildly different ones, and no F1 can tell which.  Phenotype is
annotation-neutral: it is computed by a fixed reconstruction procedure that knows
nothing about which method produced the input.  This component builds one ModelSEED
model per (held-out genome, arm) by an IDENTICAL procedure whose ONLY varying input is
the EC set, then compares each arm's model to the TRUTH model of the same genome.

WHAT IS HELD FIXED (so the only live variable is gene content)
  * EC -> model input: a template role is granted iff EVERY complete EC number in its
    name is present in the input EC set (phenotype_arm.genome_from_ec); each granted
    role becomes one synthetic MSFeature carrying the role name in ontology_terms["RAST"]
    (D10 -- the curated role->complex->reaction route, NOT the uncurated ec_rxn fan-out
    used for set scoring in score_folds.py; the two routes are never interchanged).
  * `annotate_with_rast=False` everywhere -- the modelseedpy default True would
    re-annotate the input and erase the very difference this arm measures.
  * Template is a pure function of the HELD-OUT genome's TRUE lineage, read verbatim
    from panel_eligible.json:genomes.<id>.template (never auto_select_template, which
    would make the template depend on predicted content).
  * No gapfilling for the primary readout (D11): un-gapfilled drafts give bio1 = 0.0
    with every exchange open, complete-medium-gapfilled models grow on 0/167 Biolog
    minimal media, and minimal-medium gapfilling produces models that fail the
    thermodynamic gate.  The stated growth panel returns an all-zero vector for EVERY
    arm, so it cannot discriminate.

READOUTS
  atp          PRIMARY.  54-medium core ATP-yield vector (phenotype_arm.core_atp_vector,
               MSATPCorrection on the core-template-restricted model, no gapfilling).
               Binary Hamming at 1e-6, signed set difference (media lost vs gained),
               Spearman of continuous yields, and ATP on the `empty` medium as the
               validity assertion.
  structure    draft reaction-set Jaccard, gene-associated reaction Jaccard, counts.
  growth       SECONDARY.  complete-medium-gapfilled bio1 ratio and log-ratio.
  essentiality SECONDARY.  single-reaction deletion at 1% of wild-type, evaluated ONLY
               where the truth model grows; shared-reaction disagreements and
               |essential in truth, absent from prediction| reported separately, never
               merged.

PIPELINES (A4 fatal flaw 3 -- the pipeline is itself an input-dependent actor)
  as_shipped   each model gets its own MSATPCorrection media selection and tests.
  controlled   the TRUTH model's `selected_media` and `build_tests()` output are imposed
               on every arm, so the only difference is gene content.
  The gap between them is the pipeline-reaction effect and must not be attributed to the
  methods.  This requires inlining MSBuilder.build_and_gapfill's body (msbuilder.py
  ~l.977-1050) rather than calling it -- which is also how the per-media ATP vector it
  computes and discards (`media_eval`) is captured.
  NOTE: the PRIMARY (atp) readout has no input-dependent pipeline decision at all --
  core_atp_vector evaluates all 54 media with no selection step -- so the ATP vector is
  pipeline-independent by construction and is computed once per arm and reported under
  both pipelines.  meta.atp_pipeline_independent records this, and pipeline_effect on
  atp_hamming is therefore identically 0.0 unless a gapfilled readout is requested.

WHAT THE SMOKE RUN ESTABLISHED ABOUT READOUT QUALITY (governs how any full run reads)
  Measured on the 8-genome as-shipped panel -- where EVERY Wilcoxon p on the ATP
  contrast was >= 0.20 -- plus the 12-genome e2e panel, whose primary p is 0.19 once the
  D19 prediction-source leak is fixed (the pre-fix p of 0.002 was self-leakage: the
  include_self ceiling control had overwritten every fold-0 prFBA prediction, handing
  prFBA 7 of 12 held-out genomes' own annotations).  n is tiny and nothing below is a
  result -- it is the readout's dynamic range, which is what decides what a full run can
  detect at all.
  * The 167-medium Biolog growth panel is DEAD on draft-quality reconstructions: 3
    genomes x 3 arms x 167 Carbon-* minimal media = 1,503 growth calls produced 0
    positives, for TRUTH as well as for both methods, while all 9 models grow on the
    complete medium and pass the sanity gate.  Arm-vs-truth Hamming over that panel is
    identically 0, i.e. it declares every method perfect.  Run it only as a one-paragraph
    negative result, never as a discriminating endpoint (D11).
  * The 54-medium core ATP panel -- the PRIMARY readout -- is narrow: only 19 of 54 media
    are ever positive across 112 models, and 9 of those 19 are the co-linear Pyr.* block,
    so the vector carries roughly 5 effective dimensions.
  * 2 of 8 truth genomes are phenotypically silent (all-zero truth ATP vector) -- 3 and
    4 of 12 on the two e2e panels -- and are therefore ATP-unscoreable; expect a
    comparable ~25% fraction at scale and read meta.n_unscoreable_atp before quoting any
    ATP number.
  * Structure is where the two methods actually differ and it does NOT favour prFBA:
    draft reaction-set Jaccard vs truth prFBA 0.600 vs PICRUSt2 0.645 (0.634 vs 0.674
    post-gapfill), with BOTH predicted networks ~26% larger than truth (truth 636
    metabolic reactions, prFBA 819, PICRUSt2 801) because both over-predict ECs.
  * Consequence: at this scale the phenotype arm cannot rank the methods.  It can report
    over-prediction, the pipeline effect, and the calibration ladders.

GATES (hard, per the component contract)
  atp.gate_ok            = vector["empty"] < 1e-6.  The loop-controlled core check:
                           the raw draft otherwise yields 1000 ATP from water via
                           rxn15249_c0 + rxn20583_c0 cycling cpd19001/cpd19006, and
                           loopless_solution does NOT remove it.  Biomass is not
                           evaluated here (the core model has no biomass objective).
  gapfilled.sanity_gate_ok = phenotype_arm.sanity_gate(): biomass AND ATPM on
                           water+protons alone must both be 0.
  A row failing its gate gets status "gated", is excluded from every `arms` mean, and is
  counted in meta.n_gated_out.  A gated TRUTH model leaves nothing to compare against,
  but the two gates guard different readouts.  A truth that fails the ATP empty-medium
  gate (a DRAFT-model property) makes every other arm of that genome status
  "unscoreable".  A truth that fails only the gapfill sanity gate (a property of ONE
  pipeline's GAPFILLED model) voids only that pipeline's gapfilled comparisons: the other
  arms keep status "ok", their growth/essentiality/gapfill fields are set to None, and
  vs_truth.gapfill_unscoreable = "truth_gapfill_gated".  Their draft-level ATP and
  structure comparisons stand.  Abstentions and
  missing predictions produce NO model and are counted in meta.n_unscoreable -- never
  scored as a zero phenotype.
  identity_control (C11) builds the truth EC set twice inside the same worker batch; the
  two reaction sets must have Jaccard 1.000, else degenerate optima are polluting every
  phenotype number.
  ATP PRIMARY DENOMINATOR: a genome whose TRUTH ATP vector is identically zero is
  phenotypically silent through this route, so atp_hamming against it is not a
  measurement -- every arm, ladder rung and size-matched null scores 0, and the only way
  to score above 0 is to predict a capability.  Such genomes are UNSCOREABLE for ATP
  (meta.n_unscoreable_atp, their own denominator) exactly as score_folds.py treats an
  empty truth cell; the contract-named primary (arms[*].atp_hamming,
  contrasts[*].delta_atp_hamming/wilcoxon_p/n_pairs, every ladder rung and the EC-dF1
  Spearman) is the informative subset, and the all-genome value is kept beside it as a
  labelled sensitivity (atp_all_genomes_sensitivity, delta_atp_hamming_all_genomes).
  Structure, growth and essentiality keep the full denominator -- they are well defined
  on a silent genome.
  PREDICTION SOURCE: load_prfba accepts only runs whose reference EXCLUDES the test
  genomes (mode fold / as_shipped).  The include_self and plumbing runs are C2 ceiling
  controls whose selected rep is the source genome itself; one of them silently
  overwrote every fold-0 prFBA prediction before this filter existed.

CALIBRATION LADDERS (so an EC-level dF1 becomes an expected phenotype distance)
  dropout 5/10/20/40%, false-positive addition 10/20/30/50%, donor substitution
  (nearest surviving relative / family / order), and a SIZE-MATCHED NULL per arm: a
  random EC set of the same cardinality as that arm's prediction drawn from the
  template-visible universe.  The null is what prevents a phenotype "win" that merely
  reflects predicting a larger set.  Measured on the 8-genome smoke panel: per unit of
  perturbation false negatives move the phenotype ~2x more than random false positives
  (dropout 1.125/2.50/2.75/4.00 at 5/10/20/40% vs false-positive 1.75/1.75/1.875/2.29 at
  10/20/30/50%), while the actual predictions are false-positive dominated.  Every rung
  shares the primary endpoint's denominator (see GATES), or the calibration and the arms
  would be averaged over different genome sets.

Everything expensive is per-genome and resumable: completed genomes are skipped unless
--force.  Output rows are rewritten in a deterministic (genome, arm, pipeline) order at
the end of every run, so the artefact does not depend on worker scheduling.
"""
import argparse
import csv
import glob
import gzip
import json
import math
import os
import random
import re
import subprocess
import sys
import time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
P2_SCRIPTS = os.path.join(REPO, "region_validation", "picrust2_comparison", "scripts")
if P2_SCRIPTS not in sys.path:
    sys.path.insert(0, P2_SCRIPTS)

MSD_DEFAULT = "/home/freiburger/Documents/ModelSEEDDatabase"
ATP_EMPTY_MEDIUM = "empty"
ATP_THRESH = 1e-6
GROWTH_THRESH = 1e-5
LADDER_DROPOUT = (0.05, 0.10, 0.20, 0.40)
LADDER_FP = (0.10, 0.20, 0.30, 0.50)
LADDER_DONORS = ("nn1", "family", "order")
ALL_LADDERS = ("dropout", "false_positive", "donor", "size_matched_null")
EC_OK = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

_PA = None


def pa():
    """Lazy import: keeps cobra/modelseedpy out of the parent process (spawn workers
    re-import anyway) and out of the aggregation-only code path."""
    global _PA
    if _PA is None:
        import phenotype_arm as _m
        _PA = _m
    return _PA


# --------------------------------------------------------------------------- helpers
def git_commit():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, check=True,
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        return "unknown"


def now_iso():
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def file_fp(path):
    """{size, sha256} of a file this arm's numbers rest on, or None."""
    import hashlib
    try:
        h, n = hashlib.sha256(), 0
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
                n += len(chunk)
        return {"size": n, "sha256": h.hexdigest(), "path": os.path.abspath(path)}
    except OSError:
        return None


def binding_block(args, rows, jsonl, ident_path, ec_path):
    """What ties THIS phenotype summary to the run that produced it (score_folds C11).

    ``phenotype_summary.json`` was auto-discovered by EXISTENCE beside ``--cells``
    and trusted: a DIFFERENT campaign's 12-genome summary dropped in there scored
    C11/C12 pass, and a truncated one made both controls vanish.  Nothing in the file
    said which cells table, which fold plan or which ANSWER KEY it was computed from,
    so nothing could be checked.  Everything needed to check it is recorded here, and
    ``score_folds.phenotype_binding`` re-verifies every entry.
    """
    code = {}
    for name, d in (("phenotype_fold.py", HERE), ("phenotype_arm.py", P2_SCRIPTS)):
        fp = file_fp(os.path.join(d, name))
        if fp:
            code[name] = {"sha256": fp["sha256"], "size": fp["size"]}
    run_dirs = sorted({os.path.basename(os.path.dirname(os.path.abspath(q)))
                       for root in (args.picrust2, args.prfba)
                       if root and os.path.isdir(root)
                       for q in glob.glob(os.path.join(root, "*",
                                                       "predictions.json"))})
    return {"contract": 1,
            "cells": file_fp(args.cells),
            "fold_plan": file_fp(args.fold_plan),
            "panel": file_fp(args.panel),
            "functions": file_fp(args.functions),
            "ec_universe": file_fp(ec_path) if ec_path else None,
            "outdir": os.path.abspath(args.outdir),
            "per_genome_jsonl": file_fp(jsonl),
            "identity_control_json": file_fp(ident_path),
            "run_directories_on_disk": run_dirs,
            "genomes": sorted({r["src_genome_id"] for r in rows}),
            "density": args.density, "select_mode": args.select_mode,
            "n_genomes_requested": args.n_genomes,
            "generated": now_iso(), "git_commit": git_commit(),
            "code": code}


def template_provenance(template_dir):
    """Which template bytes this run read: the directory, the sha256 of each template
    file the harness can load, and the git commit of the checkout when there is one."""
    m = pa()
    out = {"dir": template_dir, "files": {}, "git_commit": None}
    for key, name in sorted(m.TEMPLATE_FILES.items()):
        fp = file_fp(os.path.join(template_dir, name))
        out["files"][name] = fp["sha256"] if fp else None
    try:
        out["git_commit"] = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=template_dir, check=True,
            capture_output=True, text=True).stdout.strip()
    except Exception:
        pass
    return out


def library_provenance():
    """Where the libraries that build and solve the models were imported from. For one
    imported from a git checkout (the editable ModelSEEDpy fork, cobrakbase) this is the
    checkout, its commit and whether tracked files differ from that commit: an editable
    install follows whatever its working tree holds, so its version string says nothing."""
    import importlib
    import sysconfig
    site = os.path.abspath(sysconfig.get_paths()["purelib"])
    out = {}
    for name in ("modelseedpy", "cobrakbase", "cobra", "optlang", "gurobipy"):
        try:
            mod = importlib.import_module(name)
        except Exception:
            out[name] = None
            continue
        path = os.path.dirname(os.path.abspath(mod.__file__))
        rec = {"path": path, "version": str(getattr(mod, "__version__", None)),
               "git_commit": None, "git_dirty": None}
        if not path.startswith(site):
            try:
                git = lambda *a: subprocess.run(["git", *a], cwd=path, check=True,
                                                capture_output=True, text=True).stdout.strip()
                rec["git_commit"] = git("rev-parse", "HEAD")
                rec["git_dirty"] = bool(git("status", "--porcelain", "--untracked-files=no"))
            except Exception:
                pass
        out[name] = rec
    return out


def jload(path):
    import orjson
    with open(path, "rb") as fh:
        return orjson.loads(fh.read())


def norm_ec(values):
    """One normalisation for every arm: strip an EC: prefix, drop partial ECs (with '-')
    and anything that is not a 4-field EC number."""
    out = set()
    for v in values or ():
        s = str(v).strip()
        if s.upper().startswith("EC:"):
            s = s[3:].strip()
        if "-" in s or not EC_OK.match(s):
            continue
        out.add(s)
    return out


def mean(xs):
    xs = [x for x in xs if x is not None]
    return (sum(xs) / len(xs)) if xs else None


def median(xs):
    xs = sorted(x for x in xs if x is not None)
    if not xs:
        return None
    n = len(xs)
    return float(xs[n // 2]) if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


def stat_block(xs):
    xs = [x for x in xs if x is not None]
    return {"mean": mean(xs), "median": median(xs), "n": len(xs)}


def bh(pvals):
    """Benjamini-Hochberg q-values, None-safe (None stays None)."""
    idx = [i for i, p in enumerate(pvals) if p is not None and p == p]
    out = [None] * len(pvals)
    if not idx:
        return out
    try:
        from statsmodels.stats.multitest import multipletests
        q = multipletests([pvals[i] for i in idx], method="fdr_bh")[1]
        for j, i in enumerate(idx):
            out[i] = float(q[j])
    except Exception:
        m = len(idx)
        order = sorted(idx, key=lambda i: pvals[i])
        prev = 1.0
        for rank, i in enumerate(reversed(order), start=1):
            prev = min(prev, pvals[i] * m / (m - rank + 1))
            out[i] = float(prev)
    return out


F1_CONSISTENCY_TOL = 0.02       # score_folds rounds and may clip a different universe


def set_f1(pred, truth):
    """Symmetric set F1, the same quantity score_folds.py writes into cells.csv."""
    if not pred and not truth:
        return None
    tp = len(pred & truth)
    return 2.0 * tp / (len(pred) + len(truth))


def rng_for(seed, *parts):
    """Per-genome deterministic RNG: resume order and worker scheduling can never change
    a ladder sample."""
    return random.Random("{}:{}".format(seed, ":".join(str(p) for p in parts)))


# ----------------------------------------------------------------------- input loading
def load_cells(path, density, space="ec"):
    """-> {(genome_id, region): {"md5":.., "fold":.., "group_id":.., "domain":..,
                                 "arms": {arm: (outcome, f1, n_truth)}}}"""
    cells = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh):
            if (row.get("space") or "ec") != space:
                continue
            if density is not None and row.get("density"):
                try:
                    if abs(float(row["density"]) - density) > 1e-9:
                        continue
                except ValueError:
                    pass
            gid, region = row.get("src_genome_id"), row.get("region")
            if not gid:
                continue
            key = (gid, region)
            rec = cells.setdefault(key, {"md5": row.get("md5"), "region": region,
                                         "amplicon": row.get("amplicon"),
                                         "fold": row.get("fold"),
                                         "group_id": row.get("group_id"),
                                         "domain": row.get("domain"), "arms": {}})
            f1 = row.get("f1")
            try:
                f1 = float(f1) if f1 not in (None, "", "None") else None
            except ValueError:
                f1 = None
            try:
                n_truth = int(float(row.get("n_truth") or 0))
            except ValueError:
                n_truth = 0
            rec["arms"][row.get("arm")] = (row.get("outcome"), f1, n_truth)
    return cells


# run_prfba_fold.py --mode include_self / plumbing builds a reference that CONTAINS the
# test genomes: it is the C2 ceiling control (SPEC 7-C2), its selected reps ARE the
# source genome, and so its EC set IS the truth set.  Such a record must never reach a
# method arm.  load_picrust2 has always filtered on `mode`; this loader did not, and
# because `prfba/incself_d100/` sorts AFTER `prfba/f00_d100/` and covers the same md5s
# (measured: 604 of 604 on the e2e smoke panel), the ceiling control silently overwrote
# every fold-0 prediction -- giving prFBA the held-out genome's own annotation on one
# fold in k.  Allow-list the modes that are genuinely held out; `as_shipped` is
# adapt_as_shipped.py's fixture dialect (no folds, self-exclusion applied upstream).
PRFBA_MODES_OK = (None, "fold", "as_shipped")


def load_prfba(root, select_mode, density, forced=False):
    """-> {md5: {"ec": set, "reps": [...], "rank1_genome": str, "abstained": bool,
                 "tier": str}} for the matching density / forced-answer flag.
    Self-inclusive control runs are skipped, and two accepted runs may never claim the
    same md5 (fold prediction blocks are disjoint by construction, so a collision means
    the roots have been mixed -- the failure that used to pass silently)."""
    out, metas, src_of = {}, [], {}
    for path in sorted(glob.glob(os.path.join(root, "*", "predictions.json"))):
        d = jload(path)
        m = d.get("meta", {})
        if m.get("mode") not in PRFBA_MODES_OK:
            print("[pheno] skip {} (mode={}: its reference contains the test genomes, "
                  "so it is a ceiling control, not a prediction)".format(
                      path, m.get("mode")), flush=True)
            continue
        if density is not None and m.get("density") is not None:
            if abs(float(m["density"]) - density) > 1e-9:
                continue
        if bool(m.get("forced_answer", False)) != bool(forced):
            continue
        key = select_mode + ("_forced" if forced else "")
        pred = d.get("pred", {})
        block = pred.get(key) or pred.get(select_mode) or {}
        for md5, rec in block.items():
            if md5 in out and rec != out[md5]:
                raise SystemExit(
                    "[pheno] two accepted prFBA runs disagree on md5 {}: {} vs {} -- "
                    "fold prediction blocks must be disjoint".format(
                        md5, src_of[md5], path))
            out[md5], src_of[md5] = rec, path
        metas.append({"path": path, "fold": m.get("fold"), "density": m.get("density"),
                      "mode": m.get("mode"),
                      "prediction_source": m.get("prediction_source"),
                      "n_reference_md5": m.get("n_reference_md5")})
    return out, metas


def load_picrust2(root, density, routing="forced"):
    """-> {(domain, md5): rec}; domain is the reference domain the run used."""
    out, metas = {}, []
    for path in sorted(glob.glob(os.path.join(root, "*", "predictions.json"))):
        d = jload(path)
        m = d.get("meta", {})
        if m.get("mode") not in (None, "fold"):
            continue
        if density is not None and m.get("density") is not None:
            if abs(float(m["density"]) - density) > 1e-9:
                continue
        if routing and m.get("routing") and m["routing"] != routing:
            continue
        dom = m.get("domain")
        for md5, rec in d.get("pred", {}).items():
            out[(dom, md5)] = rec
        metas.append({"path": path, "fold": m.get("fold"), "domain": dom,
                      "density": m.get("density")})
    return out, metas


def dom_key(domain):
    return "arc" if (domain or "").lower().startswith("arch") else "bac"


def ec_universe(path, panel, functions):
    if path and os.path.exists(path):
        d = jload(path)
        if isinstance(d, dict) and d.get("ec"):
            return norm_ec(d["ec"])
    uni = set()
    for gid in panel.get("genomes", {}):
        fs = functions.get(gid)
        if fs:
            uni |= norm_ec(fs.get("ec"))
    return uni


# ------------------------------------------------------------------ arm EC-set assembly
def fold_reference_genomes(fold_plan):
    out = {}
    for f in fold_plan.get("folds", []):
        out[str(f["fold"])] = f.get("reference_genomes", [])
    return out


def constant_arm(ref_genomes, functions, universe):
    """Top-K most prevalent ECs among the fold's TRAINING genomes, K = median EC-set
    size over those genomes.  The floor every method must beat."""
    counts = defaultdict(int)
    sizes = []
    n = 0
    for gid in ref_genomes:
        fs = functions.get(gid)
        if not fs:
            continue
        ecs = norm_ec(fs.get("ec")) & universe
        if not ecs:
            continue
        n += 1
        sizes.append(len(ecs))
        for e in ecs:
            counts[e] += 1
    if not n:
        return set()
    k = int(median(sizes) or 0)
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return {e for e, _ in ranked[:k]}


def donor_pick(rng, ref_genomes, panel_g, gid, rank, exclude_rank=None):
    """A random surviving genome sharing `rank` with the test genome (D13).  Under
    leave-genus-out no same-genus survivor exists by construction, so the oracle is
    defined per blocked rank: family, then order, then class."""
    me = panel_g.get(gid, {})
    want = me.get(rank)
    if not want:
        return None
    pool = []
    for other in ref_genomes:
        o = panel_g.get(other)
        if not o or other == gid:
            continue
        if o.get(rank) != want:
            continue
        if exclude_rank and o.get(exclude_rank) and o.get(exclude_rank) == me.get(exclude_rank):
            continue
        pool.append(other)
    return rng.choice(sorted(pool)) if pool else None


# ----------------------------------------------------------- worker: model construction
_WORKER = {"ready": False, "integrate_reversals": False}


def use_templates(template_dir):
    """Read every ModelSEED template from `template_dir` (a checkout of
    ModelSEEDTemplates' templates/v7.0/) instead of phenotype_arm.MST. SPEC amendment
    A-C12 (2026-10-27) runs both phenotype passes on commit 21a8084, the five-reaction
    fix; with no --templates the harness reads phenotype_arm.MST exactly as before."""
    m = pa()
    if template_dir:
        m.MST = os.path.abspath(template_dir).rstrip("/") + "/"
        m._tcache.clear()
    return m.MST


def worker_init(msd_path, quiet=True, template_dir=None, integrate_reversals=False):
    import logging
    if quiet:
        for name in ("modelseedpy", "cobra", "optlang"):
            logging.getLogger(name).setLevel(logging.ERROR)
    m = pa()
    use_templates(template_dir)
    _WORKER["integrate_reversals"] = bool(integrate_reversals)
    m.preload_biochemistry(msd_path)
    try:                              # one LP thread per worker; workers are the parallelism
        import gurobipy
        gurobipy.setParam("Threads", 1)
        gurobipy.setParam("OutputFlag", 0)
    except Exception:
        pass
    m.load_template("core")
    _WORKER["ready"] = True


def safe_id(*parts):
    return re.sub(r"[^A-Za-z0-9]+", "_", "_".join(str(p) for p in parts))


def visible_ecs(role_ecs):
    out = set()
    for _name, ecs in role_ecs.values():
        out |= set(ecs)
    return out


def atp_block(m, ec, tmpl_key, model_id):
    vec = m.core_atp_vector(sorted(ec), tmpl_key, model_id, no_gapfilling=True)
    empty = float(vec.get(ATP_EMPTY_MEDIUM, 0.0))
    return {"vector": {k: round(float(v), 6) for k, v in sorted(vec.items())},
            "n_positive": sum(1 for v in vec.values() if v > ATP_THRESH),
            "empty_medium_atp": round(empty, 6),
            "gate_ok": bool(empty < ATP_THRESH)}


def draft_block(m, ec, tmpl_key, model_id):
    mdl, _b, n_feat = m.build_draft(sorted(ec), tmpl_key, model_id)
    core, gpr = m.rxn_sets(mdl)
    return mdl, n_feat, core, gpr, {
        "n_rxn": len(mdl.reactions), "n_metabolic_rxn": len(core),
        "n_gene_rxn": len(gpr), "n_genes": len(mdl.genes)}


def gapfill_block(m, ec, tmpl_key, model_id, gf_media, gf_cpds, pipeline, truth_ctx,
                  need_ess):
    """Inlined MSBuilder.build_and_gapfill (msbuilder.py ~977-1050) so that (a) the
    per-media ATP vector it discards is captured and (b) the `controlled` pipeline can
    impose the TRUTH model's selected media and tests on every arm."""
    from modelseedpy import MSGapfill
    from modelseedpy.core.msbuilder import MSBuilder
    from modelseedpy.core.msatpcorrection import MSATPCorrection, load_default_medias
    from modelseedpy.core.msmodel import get_reaction_constraints_from_direction

    tmpl, role_ecs = m.load_template(tmpl_key)
    core_t = m.load_template("core")[0]
    g = m.genome_from_ec(sorted(ec), role_ecs, model_id)
    b = MSBuilder(g, tmpl, ontology_term="RAST", template_core=core_t)
    # keyword-explicit: annotate_with_rast MUST be False (the modelseedpy default True
    # would re-annotate the input and erase the difference this arm measures)
    base = b.build_base_model(model_id, "0", allow_all_non_grp_reactions=False,
                              annotate_with_rast=False)
    MSBuilder.add_atpm(base)
    atpc = MSATPCorrection(base, core_t, load_default_medias(), compartment="c0",
                           atp_hydrolysis_id="ATPM_c0", load_default_medias=False)
    media_eval = atpc.evaluate_growth_media()
    controlled = (pipeline == "controlled" and truth_ctx is not None)
    if controlled:
        want = set(truth_ctx["selected_media_ids"])
        atpc.selected_media = [md for md in atpc.media_gapfill_stats if md.id in want]
    else:
        atpc.determine_growth_media()
    atpc.apply_growth_media_gapfilling()
    atpc.expand_model_to_genome_scale()
    tests = truth_ctx["tests"] if controlled else atpc.build_tests()

    gf = MSGapfill(base, default_gapfill_templates=[tmpl], test_conditions=tests,
                   default_target="bio1")
    infeasible = ""
    try:
        res = gf.run_gapfilling(gf_media)
    except Exception as exc:                       # noqa: BLE001 - must never crash a fold
        res, infeasible = None, "gapfill_exception:{}".format(type(exc).__name__)
    gap_sol = {}
    if res and res.get("new"):
        for rid, d in res["new"].items():
            if rid[:-1] in tmpl.reactions:
                gap_sol[rid[:-1]] = get_reaction_constraints_from_direction(d)
    elif res is None and not infeasible:
        infeasible = ("gapfill_infeasible_under_truth_tests" if controlled
                      else "gapfill_infeasible")
    proposed_rev = dict((res or {}).get("reversed") or {})
    mdl = base.copy()
    MSBuilder.integrate_gapfill_solution(tmpl, mdl, gap_sol)
    # The gapfiller also opens existing reactions in the other direction (`reversed`).
    # MSBuilder.build_and_gapfill in the ModelSEEDpy fork this harness inlines drops
    # them (upstream MSGapfill.integrate_gapfill_solution applies them), which left
    # models unable to grow (EGC_C12_REPORT.md, root cause 3). SPEC amendment A-C12
    # applies them, exactly as egc_rebuild.py's validated `reversed_only` variant does.
    applied_rev = []
    if _WORKER.get("integrate_reversals"):
        for rid, d in sorted(proposed_rev.items()):
            if rid in mdl.reactions:
                r = mdl.reactions.get_by_id(rid)
                if d in (">", "="):
                    r.upper_bound = max(r.upper_bound, 1000.0)
                if d in ("<", "="):
                    r.lower_bound = min(r.lower_bound, -1000.0)
                applied_rev.append(rid)

    gate = m.sanity_gate(mdl)
    yc = m.yield_check(mdl, gf_cpds)
    ess, wt = (None, 0.0)
    if need_ess:
        ess, wt = m.essentiality(mdl, gf_cpds)
    if not need_ess or wt <= 0.0:
        gv = m.growth_vector(mdl, {"complete": gf_cpds}, thresh=GROWTH_THRESH)
        wt = max(wt, float(gv.get("complete", 0.0)))
    gcore, ggpr = m.rxn_sets(mdl)
    block = {"n_new_rxn": len(gap_sol), "new_rxn": sorted(gap_sol),
             "reversed_proposed": {k: proposed_rev[k] for k in sorted(proposed_rev)},
             "reversed_applied": applied_rev,
             "growth": round(float(wt), 6), "sanity_gate_ok": bool(gate["pass"]),
             "yield_ratio": (None if yc["ratio_carbon"] in (float("inf"), float("-inf"))
                             else round(float(yc["ratio_carbon"]), 6)),
             "n_essential": (sum(1 for v in (ess or {}).values() if v) if ess else 0),
             "essential": sorted(k for k, v in (ess or {}).items() if v)}
    ctx = {"selected_media_ids": [md.id for md in (atpc.selected_media or [])],
            "tests": tests, "media_eval": media_eval}
    return mdl, block, ess, gcore, ggpr, ctx, infeasible


def vs_truth_block(m, arm, truth):
    """Every comparison is arm-vs-truth for the SAME genome, same pipeline."""
    out = {"atp_hamming": None, "atp_lost": [], "atp_gained": [], "atp_spearman": None,
           "atp_n_truth_positive": None, "atp_informative": None,
           "growth_ratio": None, "growth_log_ratio": None, "growth_agree": None,
           "ess_shared_disagree": None, "ess_n_shared": None,
           "ess_shared_disagree_frac": None, "ess_truth_only_absent": None,
           "rxn_jaccard_draft": None, "rxn_jaccard_gapfilled": None,
           "gene_rxn_jaccard": None, "gapfill_burden_delta": None}
    tv, av = (truth.get("atp") or {}).get("vector"), (arm.get("atp") or {}).get("vector")
    if tv and av:
        d, _n = m.hamming(tv, av, thresh=ATP_THRESH)
        out["atp_hamming"] = int(d)
        out["atp_lost"] = sorted(k for k in tv
                                 if tv[k] > ATP_THRESH and av.get(k, 0.0) <= ATP_THRESH)
        out["atp_gained"] = sorted(k for k in av
                                   if av[k] > ATP_THRESH and tv.get(k, 0.0) <= ATP_THRESH)
        keys = sorted(set(tv) | set(av))
        try:
            import warnings
            from scipy.stats import spearmanr
            with warnings.catch_warnings():     # an all-zero (phenotypically silent)
                warnings.simplefilter("ignore")  # vector is expected, not an error
                rho = spearmanr([tv.get(k, 0.0) for k in keys],
                                [av.get(k, 0.0) for k in keys]).statistic
            out["atp_spearman"] = None if rho != rho else round(float(rho), 6)
        except Exception:
            out["atp_spearman"] = None
        # A genome whose TRUTH ATP vector is all-zero is phenotypically silent through
        # the EC->role->core-template route: EVERY arm (and every ladder rung, and the
        # size-matched null) scores Hamming 0 against it, so including it in the primary
        # mean pulls every arm toward 0 by the same amount and shrinks the contrast.
        # Flag it per row so the summary can report the informative-subset mean beside
        # the all-genome mean instead of silently averaging dead cells in.
        out["atp_n_truth_positive"] = sum(1 for v in tv.values() if v > ATP_THRESH)
        out["atp_informative"] = bool(out["atp_n_truth_positive"] > 0)
    if arm.get("_draft_core") is not None and truth.get("_draft_core") is not None:
        out["rxn_jaccard_draft"] = round(m.jaccard(truth["_draft_core"],
                                                   arm["_draft_core"]), 6)
        out["gene_rxn_jaccard"] = round(m.jaccard(truth["_draft_gpr"],
                                                  arm["_draft_gpr"]), 6)
    tg, ag = truth.get("gapfilled"), arm.get("gapfilled")
    if tg and ag:
        out["gapfill_burden_delta"] = int(ag["n_new_rxn"] - tg["n_new_rxn"])
        if arm.get("_gf_core") is not None and truth.get("_gf_core") is not None:
            out["rxn_jaccard_gapfilled"] = round(m.jaccard(truth["_gf_core"],
                                                            arm["_gf_core"]), 6)
        t_g, a_g = float(tg["growth"]), float(ag["growth"])
        out["growth_ratio"] = round(a_g / t_g, 6) if t_g > 0 else None
        out["growth_log_ratio"] = round(math.log10((a_g + 1e-9) / (t_g + 1e-9)), 6)
        out["growth_agree"] = bool((a_g > GROWTH_THRESH) == (t_g > GROWTH_THRESH))
        te, ae = truth.get("_ess"), arm.get("_ess")
        if te and ae and t_g > GROWTH_THRESH:   # essentiality only where truth grows
            shared = set(te) & set(ae)
            out["ess_shared_disagree"] = sum(1 for k in shared if te[k] != ae[k])
            # The disagreement COUNT is not an agreement rate: models differ in size
            # (measured 182-287 essential reactions over 8 genomes), so 45 disagreements
            # means something different in a 600-reaction model than in a 1,000-reaction
            # one.  Carry the denominator so the two are never confused downstream.
            out["ess_n_shared"] = len(shared)
            out["ess_shared_disagree_frac"] = (
                round(out["ess_shared_disagree"] / len(shared), 6) if shared else None)
            out["ess_truth_only_absent"] = sum(1 for k, v in te.items()
                                               if v and k not in ae)
    return out


def ladder_sets(job, truth_ec, universe, visible):
    """Deterministic calibration ladders.  Seeded per genome AND per ladder step."""
    out = {}
    truth_sorted = sorted(truth_ec)
    pool = sorted((universe & visible) - truth_ec) or sorted(universe - truth_ec)
    if "dropout" in job["ladders"]:
        # ONE shuffle per genome, prefixes taken: the rungs are NESTED
        # (40% dropout removes a superset of what 20% removes), so a non-monotone
        # phenotype distance is a real cliff edge and not a resampling artefact.
        shuffled = list(truth_sorted)
        rng_for(job["seed"], job["genome_id"], "dropout").shuffle(shuffled)
        for f in LADDER_DROPOUT:
            keep = max(0, len(shuffled) - int(round(f * len(shuffled))))
            out["ladder_dropout_{:.2f}".format(f)] = set(shuffled[:keep])
    if "false_positive" in job["ladders"]:
        fp_pool = list(pool)                 # nested for the same reason
        rng_for(job["seed"], job["genome_id"], "fp").shuffle(fp_pool)
        for f in LADDER_FP:
            k = min(len(fp_pool), int(round(f * len(truth_sorted))))
            out["ladder_false_positive_{:.2f}".format(f)] = truth_ec | set(fp_pool[:k])
    if "donor" in job["ladders"]:
        for tag in LADDER_DONORS:
            ecs = job["donors"].get(tag)
            if ecs:
                out["ladder_donor_" + tag] = norm_ec(ecs) & universe
    if "size_matched_null" in job["ladders"]:
        vis_pool = sorted(universe & visible) or sorted(universe)
        for arm, ecs in job["arms"].items():
            if arm == "truth":
                continue
            r = rng_for(job["seed"], job["genome_id"], "null", arm)
            k = min(len(vis_pool), len(ecs))
            out["ladder_size_matched_null_" + arm] = set(r.sample(vis_pool, k))
    return out


def run_genome(job):
    """One genome, every arm, every pipeline, every ladder -- all in ONE process so the
    `controlled` pipeline can reuse the truth model's media selection and tests, and so
    the identity control is genuinely 'in the same batch'."""
    t_all = time.time()
    m = pa()
    tmpl_key = job["template"]
    tmpl, role_ecs = m.load_template(tmpl_key)
    visible = visible_ecs(role_ecs)
    universe = set(job["ec_universe"])
    gf_media, gf_cpds = (m.complete_media()
                         if ({"growth", "essentiality"} & set(job["readouts"]))
                         else (None, None))
    want_atp = "atp" in job["readouts"]
    want_struct = "structure" in job["readouts"]
    want_gf = bool({"growth", "essentiality"} & set(job["readouts"]))
    need_ess = "essentiality" in job["readouts"]

    truth_ec = norm_ec(job["arms"]["truth"])
    ec_by_arm = {a: norm_ec(v) for a, v in job["arms"].items()}
    ec_by_arm.update(ladder_sets(job, truth_ec, universe, visible))
    # truth is always first: every vs_truth and the controlled pipeline depend on it.
    arm_order = ["truth"] + [a for a in job["arms"] if a != "truth"] + \
                [a for a in ec_by_arm if a.startswith("ladder_")]

    rows, built, identity = [], {}, None
    truth_ctx = {}
    for arm in arm_order:
        ec = ec_by_arm[arm]
        res = {"atp": None, "gapfilled": None, "_draft_core": None, "_draft_gpr": None,
               "_gf_core": None, "_ess": None, "draft": None, "n_roles": 0,
               "n_features": 0, "infeasible": {}, "_gf": {}}
        t0 = time.time()
        try:
            g = m.genome_from_ec(sorted(ec), role_ecs, "probe")
            res["n_roles"] = res["n_features"] = len(g.features)
            res["template_visible_frac"] = (len(ec & visible) / len(ec)) if ec else None
            if want_atp:
                res["atp"] = atp_block(m, ec, tmpl_key, safe_id(job["genome_id"], arm, "atp"))
            if want_struct:
                mdl, n_feat, core, gpr, counts = draft_block(
                    m, ec, tmpl_key, safe_id(job["genome_id"], arm, "draft"))
                res["draft"], res["_draft_core"], res["_draft_gpr"] = counts, core, gpr
                if job.get("dump_models"):
                    import cobra.io
                    os.makedirs(job["dump_models"], exist_ok=True)
                    cobra.io.save_json_model(mdl, os.path.join(
                        job["dump_models"], "{}__{}__draft.json".format(job["genome_id"], arm)))
            for pipe in job["pipelines"]:
                if not want_gf:
                    continue
                # The calibration ladders (9.3) are an ATP-Hamming calibration BY DESIGN
                # -- summarize() reads only vs_truth.atp_hamming for every rung, and
                # nothing reads a ladder row's `gapfilled` block.  Gapfilling them is
                # therefore pure cost, and it is the dominant cost: measured 1,640 s per
                # genome for 14 arms x 2 pipelines gapfilled vs 415 s for the 3 real arms
                # x 2 pipelines (~91 core-hours vs ~46 at the orchestrator's 200
                # genomes).  Skip it, so `--readouts atp,growth,essentiality,structure`
                # with the default --ladders is affordable.
                if arm.startswith("ladder_"):
                    continue
                ctx = truth_ctx.get(pipe) if arm != "truth" else None
                # The truth arm DEFINES the pipeline, so its own gapfill is identical under
                # both variants: build it once and alias the second pipeline to it.
                if arm == "truth" and truth_ctx:
                    prior = next(iter(truth_ctx))
                    res["_gf"][pipe] = res["_gf"][prior]
                    res["infeasible"][pipe] = res["infeasible"].get(prior, "")
                    truth_ctx[pipe] = truth_ctx[prior]
                    continue
                mdl, block, ess, gcore, ggpr, newctx, infeas = gapfill_block(
                    m, ec, tmpl_key, safe_id(job["genome_id"], arm, pipe), gf_media,
                    gf_cpds, pipe, ctx, need_ess)
                res["_gf"][pipe] = {"block": block, "ess": ess,
                                    "core": gcore, "gpr": ggpr}
                res["infeasible"][pipe] = infeas
                if arm == "truth":
                    truth_ctx[pipe] = newctx
        except Exception as exc:                   # noqa: BLE001
            # Isolate the arm, not the genome: one exotic arm (a ladder, a
            # baseline) must not discard the truth/prFBA/PICRUSt2 rows.  A
            # TRUTH failure has nothing to compare against, so it cascades.
            if arm == "truth":
                raise
            built[arm] = {"_failed": "build_exception:{}:{}".format(
                type(exc).__name__, exc)[:300]}
            continue
        res["secs"] = round(time.time() - t0, 2)
        built[arm] = res

        if arm == "truth" and job.get("identity_control"):
            ident = {"jaccard": None, "pass": None, "atp_identical": None}
            _mdl2, _nf2, core2, _gpr2, _c2 = draft_block(
                m, ec, tmpl_key, safe_id(job["genome_id"], "truth", "identity")) \
                if want_struct else (None, None, None, None, None)
            if core2 is not None:
                ident["jaccard"] = round(m.jaccard(res["_draft_core"], core2), 6)
            if want_atp:
                v2 = atp_block(m, ec, tmpl_key, safe_id(job["genome_id"], "idatp"))
                ident["atp_identical"] = bool(v2["vector"] == res["atp"]["vector"])
            if want_gf and res["_gf"]:
                pipe0 = job["pipelines"][0]
                _m2, b2, _e2, _c, _g, _ctx, _i = gapfill_block(
                    m, ec, tmpl_key, safe_id(job["genome_id"], "idgf"), gf_media,
                    gf_cpds, "as_shipped", None, False)
                ident["gapfill_identical"] = bool(
                    b2["new_rxn"] == res["_gf"][pipe0]["block"]["new_rxn"])
            ident["pass"] = bool((ident["jaccard"] in (None, 1.0))
                                 and ident.get("atp_identical") in (None, True)
                                 and ident.get("gapfill_identical", True))
            identity = ident

    # ------------------------------------------------------------------ emit the rows
    # The two truth gates guard different readouts (see the module docstring).  Before
    # 2026-10-18 only the ATP gate was propagated, so a truth model that failed its
    # gapfill sanity gate stayed the reference for every other arm's gapfilled readouts.
    truth_atp_ok = True
    if want_atp and built["truth"]["atp"] is not None:
        truth_atp_ok = built["truth"]["atp"]["gate_ok"]

    def truth_gf_ok(pipe):
        tg = built["truth"]["_gf"].get(pipe)
        return not (tg and not tg["block"]["sanity_gate_ok"])

    for arm in arm_order:
        res = built[arm]
        if res.get("_failed"):                 # arm-level failure: fully-keyed rows,
            for pipe in job["pipelines"]:      # never a zero phenotype
                rows.append(empty_row(job, arm, pipe, "unscoreable", res["_failed"],
                                      n_ec=len(ec_by_arm[arm])))
            continue
        for pipe in job["pipelines"]:
            gfp = res["_gf"].get(pipe)
            arm_view = {"atp": res["atp"],
                        "gapfilled": gfp["block"] if gfp else None,
                        "_draft_core": res["_draft_core"], "_draft_gpr": res["_draft_gpr"],
                        "_gf_core": gfp["core"] if gfp else None,
                        "_ess": gfp["ess"] if gfp else None}
            tgfp = built["truth"]["_gf"].get(pipe)
            truth_view = {"atp": built["truth"]["atp"],
                          "gapfilled": tgfp["block"] if tgfp else None,
                          "_draft_core": built["truth"]["_draft_core"],
                          "_draft_gpr": built["truth"]["_draft_gpr"],
                          "_gf_core": tgfp["core"] if tgfp else None,
                          "_ess": tgfp["ess"] if tgfp else None}
            vt = vs_truth_block(m, arm_view, truth_view)
            if arm == "truth":              # self-comparison: makes C11 trivially checkable
                vt["atp_hamming"] = 0 if vt["atp_hamming"] is not None else None
                vt["atp_lost"], vt["atp_gained"] = [], []
            gate_ok = True if res["atp"] is None else res["atp"]["gate_ok"]
            if gfp and not gfp["block"]["sanity_gate_ok"]:
                gate_ok = False
            status, reason = "ok", ""
            if not gate_ok:
                status, reason = "gated", ("atp_empty_medium_nonzero"
                                           if (res["atp"] and not res["atp"]["gate_ok"])
                                           else "sanity_gate_failed")
            elif arm != "truth" and not truth_atp_ok:
                status, reason = "unscoreable", "truth_model_gated"
            if res["infeasible"].get(pipe):
                status, reason = "unscoreable", res["infeasible"][pipe]
            if arm != "truth" and status == "ok" and not truth_gf_ok(pipe):
                for k in GAPFILL_VS_TRUTH_KEYS:
                    vt[k] = None
                vt["gapfill_unscoreable"] = "truth_gapfill_gated"
            rows.append({
                "src_genome_id": job["genome_id"], "fold": job["fold"],
                "group_id": job.get("group_id"), "region": job.get("region"),
                "md5": job.get("md5"), "arm": arm, "pipeline": pipe,
                "template": job["template"], "domain": job["domain"],
                "n_ec_input": len(ec_by_arm[arm]),
                "template_visible_frac": (None if res["template_visible_frac"] is None
                                          else round(res["template_visible_frac"], 6)),
                "n_roles": res["n_roles"], "n_features": res["n_features"],
                "draft": res["draft"] or {"n_rxn": None, "n_metabolic_rxn": None,
                                          "n_gene_rxn": None, "n_genes": None},
                "atp": res["atp"] or {"vector": {}, "n_positive": None,
                                      "empty_medium_atp": None, "gate_ok": None},
                "gapfilled": gfp["block"] if gfp else {
                    "n_new_rxn": None, "new_rxn": [], "growth": None,
                    "sanity_gate_ok": None, "yield_ratio": None,
                    "n_essential": None, "essential": []},
                "vs_truth": vt,
                "ec_f1": job.get("ec_f1", {}).get(arm),
                "readouts": sorted(job["readouts"]),
                "ladders": sorted(job["ladders"]),
                "wall_seconds": res["secs"], "status": status,
                "status_reason": reason})
    return {"genome_id": job["genome_id"], "rows": rows, "identity_control": identity,
            "wall_seconds": round(time.time() - t_all, 2)}


def empty_row(job, arm, pipe, status, reason, n_ec=None):
    """A fully-keyed row for a (genome, arm, pipeline) that produced no model, so the
    artefact schema never depends on whether a build succeeded."""
    return {"src_genome_id": job["genome_id"], "fold": job["fold"],
            "group_id": job.get("group_id"), "region": job.get("region"),
            "md5": job.get("md5"), "arm": arm, "pipeline": pipe,
            "template": job["template"], "domain": job["domain"],
            "n_ec_input": n_ec, "template_visible_frac": None,
            "n_roles": None, "n_features": None,
            "draft": {"n_rxn": None, "n_metabolic_rxn": None, "n_gene_rxn": None,
                      "n_genes": None},
            "atp": {"vector": {}, "n_positive": None, "empty_medium_atp": None,
                    "gate_ok": None},
            "gapfilled": {"n_new_rxn": None, "new_rxn": [], "growth": None,
                          "sanity_gate_ok": None, "yield_ratio": None,
                          "n_essential": None, "essential": []},
            "vs_truth": {k: v for k, v in
                         (("atp_hamming", None), ("atp_lost", []), ("atp_gained", []),
                          ("atp_spearman", None), ("growth_ratio", None),
                          ("growth_log_ratio", None), ("growth_agree", None),
                          ("ess_shared_disagree", None), ("ess_truth_only_absent", None),
                          ("rxn_jaccard_draft", None), ("rxn_jaccard_gapfilled", None),
                          ("gene_rxn_jaccard", None), ("gapfill_burden_delta", None))},
            "ec_f1": job.get("ec_f1", {}).get(arm),
            "readouts": sorted(job["readouts"]), "ladders": sorted(job["ladders"]),
            "wall_seconds": 0.0, "status": status, "status_reason": reason}


def run_genome_safe(job):
    """One bad genome must not kill a 1,000-genome batch: the failure is recorded as
    unscoreable rows (never as a zero phenotype) and the run continues."""
    try:
        return run_genome(job)
    except Exception as exc:                       # noqa: BLE001
        import traceback
        traceback.print_exc()
        reason = "build_exception:{}:{}".format(type(exc).__name__, exc)[:300]
        rows = [empty_row(job, arm, pipe, "unscoreable", reason)
                for arm in job["arms"] for pipe in job["pipelines"]]
        return {"genome_id": job["genome_id"], "rows": rows, "identity_control": None,
                "wall_seconds": 0.0, "error": reason}


# --------------------------------------------------------------------------- selection
def quintile(rank, n, nbins=5):
    if n <= 1:
        return 0
    nbins = max(1, min(nbins, n))
    return min(nbins - 1, int(nbins * rank / n))


def choose_genomes(cands, n_target, seed, experimental):
    """Stratified by domain x template x Delta(prFBA, PICRUSt2) EC-F1 quintile (A1
    amendment G): the arm can only test whether EC disagreement propagates if the sample
    spans the disagreement range.  Genomes with experimental data are taken first."""
    if n_target >= len(cands):
        return sorted(cands, key=lambda c: c["genome_id"])
    ordered = sorted(cands, key=lambda c: ((c["delta_f1"] if c["delta_f1"] is not None
                                            else 0.0), c["genome_id"]))
    n = len(ordered)
    strata = defaultdict(list)
    for i, c in enumerate(ordered):
        c["quintile"] = quintile(i, n)
        strata[(c["domain"], c["template"], c["quintile"])].append(c)
    picked, pool = [], []
    exp = set(experimental or ())
    for c in ordered:
        (picked if c["genome_id"] in exp else pool).append(c)
    picked = picked[:n_target]
    remaining = n_target - len(picked)
    if remaining > 0:
        by_s = defaultdict(list)
        for c in pool:
            by_s[(c["domain"], c["template"], c["quintile"])].append(c)
        keys = sorted(by_s)
        quota, rema = {}, []
        for k in keys:
            exact = remaining * len(by_s[k]) / len(pool)
            quota[k] = int(exact)
            rema.append((exact - int(exact), k))
        for _frac, k in sorted(rema, reverse=True):
            if sum(quota.values()) >= remaining:
                break
            quota[k] += 1
        for k in keys:
            r = random.Random("{}:{}".format(seed, k))
            members = sorted(by_s[k], key=lambda c: c["genome_id"])
            r.shuffle(members)
            picked.extend(members[:quota[k]])
        if len(picked) < n_target:          # top-up when strata were short
            have = {c["genome_id"] for c in picked}
            extra = [c for c in ordered if c["genome_id"] not in have]
            r = random.Random("{}:topup".format(seed))
            r.shuffle(extra)
            picked.extend(extra[:n_target - len(picked)])
    return sorted(picked[:n_target], key=lambda c: c["genome_id"])


# -------------------------------------------------------------------------- aggregation
METRIC_KEYS = (("atp_hamming", lambda r: r["vs_truth"]["atp_hamming"]),
               ("rxn_jaccard_draft", lambda r: r["vs_truth"]["rxn_jaccard_draft"]),
               ("growth_log_ratio", lambda r: r["vs_truth"]["growth_log_ratio"]),
               ("ess_shared_disagree", lambda r: r["vs_truth"]["ess_shared_disagree"]),
               ("ess_shared_disagree_frac",
                lambda r: r["vs_truth"].get("ess_shared_disagree_frac")),
               ("gapfill_burden", lambda r: None if r["vs_truth"].get("gapfill_unscoreable")
                else r["gapfilled"].get("n_new_rxn")))

# vs_truth fields that compare GAPFILLED models.  Voided when the truth model of that
# pipeline fails its gapfill sanity gate, so every gapfilled readout shares one denominator.
GAPFILL_VS_TRUTH_KEYS = ("growth_ratio", "growth_log_ratio", "growth_agree",
                         "ess_shared_disagree", "ess_n_shared", "ess_shared_disagree_frac",
                         "ess_truth_only_absent", "rxn_jaccard_gapfilled",
                         "gapfill_burden_delta")


def dedupe_rows(rows):
    """One row per (genome, arm, pipeline).

    per_genome.jsonl accumulates rows from every invocation that wrote into it -- the ATP
    pass (atp,structure) and the gapfilled pass (atp,growth,essentiality,structure) share
    one file -- so a genome in both passes has two rows per arm and pipeline.  Averaging
    over the raw list counted those genomes twice in every per-arm statistic and ladder
    rung, while the contrasts' {genome: row} dicts kept whichever row came last in file
    order.  Keep the row with the most readouts (a superset of the other's; the ATP
    vector is identical across passes, the identity control guarantees it), ties to the
    later row, so every consumer sees the same single row.
    """
    best = {}
    for i, r in enumerate(rows):
        k = (r["src_genome_id"], r["arm"], r["pipeline"])
        rank = (len(r.get("readouts") or []), i)
        if k not in best or rank > best[k][0]:
            best[k] = (rank, r)
    return [v[1] for v in best.values()]

# The ATP metrics.  Their PRIMARY denominator is the genomes whose TRUTH vector has at
# least one ATP-positive medium (see ATP_PRIMARY_DENOM below); the all-genome value is
# kept beside them as a labelled sensitivity, never as the contract number.
INFORMATIVE_KEYS = (("atp_hamming", lambda r: r["vs_truth"]["atp_hamming"]),
                    ("atp_spearman", lambda r: r["vs_truth"].get("atp_spearman")))

# Metric names in METRIC_KEYS whose denominator is the informative subset.  The
# structural / growth / essentiality metrics are NOT degenerate on a phenotypically
# silent genome (its reaction network and biomass flux are perfectly well defined), so
# they keep the full denominator.
ATP_METRICS = ("atp_hamming",)
ATP_PRIMARY_DENOM = "informative_truth_only (>=1 ATP-positive medium in the truth vector)"


def is_informative(r):
    """True iff this row's truth ATP vector was not all-zero.  Falls back to recomputing
    from the row's own vector for rows written before `atp_informative` existed (an arm
    row carries its own vector, so only the truth row can answer this; the summary
    therefore resolves it per genome, see summarize())."""
    v = r["vs_truth"].get("atp_informative")
    return None if v is None else bool(v)


def summarize(rows, arms, pipelines, frozen, meta_extra, ladders_run):
    ok = dedupe_rows([r for r in rows if r["status"] == "ok"])
    by = defaultdict(list)
    for r in ok:
        by[(r["arm"], r["pipeline"])].append(r)

    # Informativeness is a property of the GENOME (its truth ATP vector), so resolve it
    # ONCE from the truth rows -- which also makes it available for rows written before
    # vs_truth.atp_informative existed -- and use that single predicate everywhere a
    # per-genome atp_hamming is read (arms, by_domain, contrasts, the EC-dF1 Spearman
    # and every ladder rung), so the primary endpoint and its calibration share one
    # denominator.  A genome whose truth ATP vector is identically zero is
    # phenotypically silent through the EC->role->core route: EVERY arm, ladder rung and
    # size-matched null scores Hamming 0 against it, so atp_hamming there is not a
    # measurement of the method but a penalty for predicting any capability at all.
    # Such genomes are unscoreable for ATP (meta.n_unscoreable_atp), exactly as the EC
    # scorer treats an empty truth cell.
    truth_npos = {r["src_genome_id"]: (r["atp"] or {}).get("n_positive")
                  for r in rows if r["arm"] == "truth"}
    informative = {g: (n is not None and n > 0) for g, n in truth_npos.items()}

    def inf_rows(rs):
        """Rows on the ATP primary denominator (unknown genomes are kept, as before)."""
        return [r for r in rs if informative.get(r["src_genome_id"], True)]

    arms_out = {}
    for arm in arms:
        arms_out[arm] = {}
        for pipe in pipelines:
            rs = by.get((arm, pipe), [])
            inf = inf_rows(rs)
            # ATP metrics on the informative denominator (PRIMARY), everything else on
            # the full one.
            block = {name: stat_block([fn(r) for r in (inf if name in ATP_METRICS
                                                       else rs)])
                     for name, fn in METRIC_KEYS}
            block["atp_primary_denominator"] = ATP_PRIMARY_DENOM
            # alias of the primary (same rows), kept so the swap is self-documenting
            block["atp_informative_only"] = {
                name: stat_block([fn(r) for r in inf]) for name, fn in INFORMATIVE_KEYS}
            block["atp_informative_only"]["n_genomes"] = len(inf)
            # labelled SENSITIVITY: the old all-genome value, silent truth included
            block["atp_all_genomes_sensitivity"] = {
                name: stat_block([fn(r) for r in rs]) for name, fn in INFORMATIVE_KEYS}
            block["atp_all_genomes_sensitivity"]["n_genomes"] = len(rs)
            doms = defaultdict(list)
            for r in rs:
                doms[r["domain"] or "unknown"].append(r)
            block["by_domain"] = {
                d: {name: stat_block([fn(r) for r in (inf_rows(drs)
                                                      if name in ATP_METRICS else drs)])
                    for name, fn in METRIC_KEYS}
                for d, drs in sorted(doms.items())}
            arms_out[arm][pipe] = block

    # ----- contrasts: paired per organism, Wilcoxon signed-rank, McNemar per medium
    pairs = [(a, b) for i, a in enumerate(arms) for b in arms[i + 1:]
             if a != "truth" and b != "truth"]
    contrasts, pvals, ptags = {}, [], []
    for a, b in pairs:
        key = "{}__vs__{}".format(a, b)
        contrasts[key] = {}
        for pipe in pipelines:
            ra = {r["src_genome_id"]: r for r in by.get((a, pipe), [])}
            rb = {r["src_genome_id"]: r for r in by.get((b, pipe), [])}
            rt = {r["src_genome_id"]: r for r in by.get(("truth", pipe), [])}
            shared = sorted(set(ra) & set(rb))
            # PRIMARY: paired per organism over the ATP-informative genomes only.
            inf_shared = [g for g in shared if informative.get(g, True)]

            def paired_diffs(gs):
                d = [(ra[g]["vs_truth"]["atp_hamming"],
                      rb[g]["vs_truth"]["atp_hamming"]) for g in gs]
                return [x - y for x, y in d if x is not None and y is not None]

            def wilcox(dd):
                if len(dd) >= 2 and any(x != 0 for x in dd):
                    try:
                        from scipy.stats import wilcoxon
                        return float(wilcoxon(dd).pvalue)
                    except Exception:
                        return None
                return None

            diffs = paired_diffs(inf_shared)        # contract primary
            diffs_all = paired_diffs(shared)        # labelled sensitivity
            p, p_all = wilcox(diffs), wilcox(diffs_all)

            # EC dF1 vs phenotype dHamming: does annotation disagreement propagate?
            # On a silent genome the phenotype difference is identically 0 whatever the
            # EC difference is, so this shares the primary's denominator; the
            # all-genome value is reported beside it.
            def pheno_spearman(gs):
                xs = [(ra[g]["ec_f1"], rb[g]["ec_f1"], ra[g]["vs_truth"]["atp_hamming"],
                       rb[g]["vs_truth"]["atp_hamming"]) for g in gs]
                xs = [(f1a - f1b, ha - hb) for f1a, f1b, ha, hb in xs
                      if None not in (f1a, f1b, ha, hb)]
                if len(xs) < 3:
                    return None
                try:
                    from scipy.stats import spearmanr
                    rr = spearmanr([x for x, _ in xs], [y for _, y in xs]).statistic
                    return None if rr != rr else round(float(rr), 6)
                except Exception:
                    return None

            rho, rho_all = pheno_spearman(inf_shared), pheno_spearman(shared)
            # McNemar per ATP medium: b = arm A agrees with truth where B does not.
            # Deliberately over ALL shared genomes: the unit here is the MEDIUM, and on
            # a silent genome an arm predicting ATP where truth has none is a genuine
            # medium-level false positive, not a degenerate cell.
            media = sorted({k for g in shared for k in rt.get(g, {}).get("atp", {})
                            .get("vector", {})})
            mc, mp = {}, []
            for med in media:
                b_ct = c_ct = both = neither = 0
                for g in shared:
                    tv = rt.get(g, {}).get("atp", {}).get("vector", {})
                    if med not in tv:
                        continue
                    t_pos = tv[med] > ATP_THRESH
                    a_ok = (ra[g]["atp"]["vector"].get(med, 0.0) > ATP_THRESH) == t_pos
                    b_ok = (rb[g]["atp"]["vector"].get(med, 0.0) > ATP_THRESH) == t_pos
                    if a_ok and b_ok:
                        both += 1
                    elif a_ok and not b_ok:
                        b_ct += 1
                    elif b_ok and not a_ok:
                        c_ct += 1
                    else:
                        neither += 1
                pm = None
                if b_ct + c_ct > 0:
                    try:
                        from statsmodels.stats.contingency_tables import mcnemar
                        pm = float(mcnemar([[both, b_ct], [c_ct, neither]],
                                           exact=True).pvalue)
                    except Exception:
                        pm = None
                mc[med] = {"p": pm, "q_bh": None, "b": b_ct, "c": c_ct}
                mp.append(pm)
            for med, q in zip(media, bh(mp)):
                mc[med]["q_bh"] = q
            d_inf = round(mean(diffs), 6) if diffs else None
            contrasts[key][pipe] = {
                # ---- contract primary (ATP-informative genomes only)
                "delta_atp_hamming": d_inf,
                "wilcoxon_p": p, "q_bh": None, "n_pairs": len(diffs),
                "atp_denominator": ATP_PRIMARY_DENOM,
                # aliases of the primary, kept so the swap is self-documenting
                "delta_atp_hamming_informative": d_inf,
                "n_pairs_informative": len(diffs),
                "wilcoxon_p_informative": p,
                # ---- labelled sensitivity: every genome, silent truth included
                "delta_atp_hamming_all_genomes": (round(mean(diffs_all), 6)
                                                  if diffs_all else None),
                "wilcoxon_p_all_genomes": p_all,
                "n_pairs_all_genomes": len(diffs_all),
                "mcnemar_per_medium": mc,
                "ec_f1_delta_vs_phenotype_spearman": rho,
                "ec_f1_delta_vs_phenotype_spearman_all_genomes": rho_all}
            pvals.append(p)          # BH over the PRIMARY p-values
            ptags.append((key, pipe))
    for (key, pipe), q in zip(ptags, bh(pvals)):
        contrasts[key][pipe]["q_bh"] = q

    # ----- ladders (always full-skeleton, null where the ladder was not run)
    def lad(prefix, tags):
        """Rungs are the calibration for the primary endpoint, so they must share its
        denominator -- a silent genome scores 0 on every rung by construction, which
        would compress the whole ladder toward 0 and make the comparison with the arms
        apples-to-oranges."""
        out = {}
        for tag in tags:
            rs = [r for r in ok if r["arm"] == "{}{}".format(prefix, tag)]
            vals = [v for v in (r["vs_truth"]["atp_hamming"] for r in inf_rows(rs))
                    if v is not None]
            vals_all = [v for v in (r["vs_truth"]["atp_hamming"] for r in rs)
                        if v is not None]
            out[str(tag)] = {"atp_hamming_mean": (round(mean(vals), 6) if vals else None),
                             "n": len(vals),
                             "atp_denominator": ATP_PRIMARY_DENOM,
                             "atp_hamming_mean_all_genomes": (round(mean(vals_all), 6)
                                                              if vals_all else None),
                             "n_all_genomes": len(vals_all)}
        return out

    ladders = {
        "dropout": lad("ladder_dropout_", ["{:.2f}".format(f) for f in LADDER_DROPOUT]),
        "false_positive": lad("ladder_false_positive_",
                              ["{:.2f}".format(f) for f in LADDER_FP]),
        "donor": lad("ladder_donor_", list(LADDER_DONORS)),
        "size_matched_null": lad("ladder_size_matched_null_",
                                 [a for a in arms if a != "truth"])}
    # The nn1 rung (and the base_nn1 arm) copy prfba.rank1_genome.  In a fold run that
    # is the alignment rank-1 hit against the fold-restricted shared reference DB, i.e.
    # a genuine nearest surviving relative; in the adapt_as_shipped.py fixture it is
    # sorted(reps)[0], i.e. prFBA's first selected rep, which is NOT an independent
    # neighbour.  Carry the provenance wherever the number is read.
    if "nn1" in ladders.get("donor", {}):
        ladders["donor"]["nn1"]["donor_source"] = meta_extra.get("donor_nn1_source")

    # The pipeline-reaction effect (D11/9.2-6): how much of an arm's apparent distance
    # from truth comes from MSATPCorrection's input-dependent media selection rather
    # than from gene content.  atp_hamming is identically 0 here because the core ATP
    # vector has no selection step (meta.atp_pipeline_independent); the gapfilled
    # metrics are where this effect is real, so they are reported too.
    pipeline_effect = {}
    for arm in arms:
        block = {}
        for name, fn in METRIC_KEYS:
            v = {}
            for pipe in pipelines:
                vals = [fn(r) for r in by.get((arm, pipe), [])]
                v[pipe] = mean([x for x in vals if x is not None])
            x, y = v.get("as_shipped"), v.get("controlled")
            block["{}_as_shipped_minus_controlled".format(name)] = (
                round(x - y, 6) if (x is not None and y is not None) else None)
        pipeline_effect[arm] = block

    gated = defaultdict(int)
    unscoreable = defaultdict(int)
    for r in rows:
        if r["status"] == "gated":
            gated[r["arm"]] += 1
        elif r["status"] == "unscoreable":
            unscoreable[r["arm"]] += 1
    # Genomes whose TRUTH ATP vector is all-zero are phenotypically silent through this
    # route (no ATP-positive medium at all): every arm, rung and null scores Hamming 0
    # against them, so they carry no method signal and are UNSCOREABLE for ATP -- their
    # own denominator here, not a diluted primary mean.
    silent = sorted(g for g, n in truth_npos.items() if n == 0)
    n_inf = sum(1 for v in informative.values() if v)
    undefined = sorted(g for g, n in truth_npos.items() if n is None)
    meta = {"frozen": frozen,
            "ladders_run": sorted(ladders_run),
            "atp_primary_denominator": ATP_PRIMARY_DENOM,
            # n_atp_informative + truth_atp_allzero + truth_atp_undefined == the number
            # of genomes carrying a truth row, by construction of one predicate.
            "n_unscoreable_atp": {
                "truth_atp_allzero": len(silent),
                "truth_atp_undefined": len(undefined),
                "n_atp_scoreable_genomes": n_inf,
                "n_truth_rows": len(truth_npos),
                "genomes": silent,
                "denominator": ATP_PRIMARY_DENOM},
            "n_truth_atp_allzero": len(silent),
            "truth_atp_allzero_genomes": silent,
            "n_atp_informative": n_inf,
            "n_genomes": len({r["src_genome_id"] for r in rows}),
            "n_models": len(rows),
            "n_gated_out": dict(sorted(gated.items())),
            "n_unscoreable": dict(sorted(unscoreable.items())),
            "generated": now_iso(), "git_commit": git_commit()}
    meta.update(meta_extra)
    return {"meta": meta, "arms": arms_out, "contrasts": contrasts,
            "ladders": ladders, "pipeline_effect": pipeline_effect}


# --------------------------------------------------------------------------------- main
def build_jobs(args, panel, fold_plan, functions, universe, cells, prf, p2):
    panel_g = panel.get("genomes", {})
    ref_by_fold = fold_reference_genomes(fold_plan)
    genome_fold = {}
    for f in fold_plan.get("folds", []):
        for gid in f.get("test_genomes", []):
            genome_fold[gid] = str(f["fold"])
    prfba_arm = args.prfba_arm
    const_cache = {}

    # ---- candidates: one amplicon per genome, both methods scored, truth big enough
    cands = {}
    for (gid, region), rec in cells.items():
        if gid not in panel_g:
            continue
        a_p = rec["arms"].get(prfba_arm)
        a_2 = rec["arms"].get("picrust2")
        if not a_p or not a_2 or a_p[0] != "scored" or a_2[0] != "scored":
            continue
        truth_ec = norm_ec((functions.get(gid) or {}).get("ec")) & universe
        if len(truth_ec) < args.min_truth_ec:
            continue
        md5 = rec["md5"]
        if md5 not in prf:
            continue
        dom = rec["domain"] or panel_g[gid].get("domain")
        if (dom_key(dom), md5) not in p2:
            continue
        fold_val = rec["fold"] if rec["fold"] not in (None, "") else genome_fold.get(gid)
        try:
            fold_val = int(float(fold_val))
        except (TypeError, ValueError):
            pass
        cand = {"genome_id": gid, "region": region, "md5": md5,
                "amplicon": rec["amplicon"],
                "fold": fold_val,
                "group_id": rec["group_id"], "domain": dom,
                "template": panel_g[gid].get("template")
                            or pa_template_fallback(panel_g[gid]),
                "delta_f1": (None if (a_p[1] is None or a_2[1] is None)
                             else a_p[1] - a_2[1]),
                "f1": {prfba_arm: a_p[1], "picrust2": a_2[1]}}
        prev = cands.get(gid)
        if prev is None or region_rank(region, args.region_pref) < \
                region_rank(prev["region"], args.region_pref):
            cands[gid] = cand
    chosen = choose_genomes(list(cands.values()), args.n_genomes, args.seed,
                            load_experimental(args.experimental_list))

    jobs = []
    for c in chosen:
        gid, md5 = c["genome_id"], c["md5"]
        fold_key = str(c["fold"])               # fold_plan keys are strings
        fold = c["fold"]                        # ...but the artefact carries the int
        truth_ec = norm_ec((functions.get(gid) or {}).get("ec")) & universe
        refs = ref_by_fold.get(fold_key, [])
        arms, ec_f1 = {"truth": sorted(truth_ec)}, {}
        prec = prf.get(md5, {})
        p2rec = p2.get((dom_key(c["domain"]), md5), {})
        if prfba_arm in args.arms:
            arms[prfba_arm] = sorted(norm_ec(prec.get("ec")) & universe)
            ec_f1[prfba_arm] = c["f1"].get(prfba_arm)
        if "picrust2" in args.arms:
            arms["picrust2"] = sorted(norm_ec(p2rec.get("ec")) & universe)
            ec_f1["picrust2"] = c["f1"].get("picrust2")
        # The cell's F1 was scored by score_folds.py from the same prediction run; if
        # the EC set recomputed here disagrees, the artefact feeding the model is NOT
        # the run that produced the cell (this is the guard that catches a ceiling- or
        # density-control run leaking into a method arm -- it fires at 1.0 Jaccard
        # against truth while the cell still records the honest fold F1).
        for arm_name, recorded in sorted(ec_f1.items()):
            got = set_f1(set(arms.get(arm_name) or ()), truth_ec)
            if recorded is None or got is None:
                continue
            if abs(got - float(recorded)) > F1_CONSISTENCY_TOL:
                print("[pheno] WARNING {} {}: EC F1 recomputed from the prediction "
                      "artefact ({:.4f}) != cells.csv ({:.4f}) -- the prediction source "
                      "may not be the run that produced this cell".format(
                          gid, arm_name, got, float(recorded)), flush=True)
        nn1 = prec.get("rank1_genome")
        rng = rng_for(args.seed, gid, "donor")
        fam = donor_pick(rng, refs, panel_g, gid, "family") or \
            donor_pick(rng, refs, panel_g, gid, "order")
        order_donor = donor_pick(rng_for(args.seed, gid, "donor_order"), refs, panel_g,
                                 gid, "order", exclude_rank="family")
        donors = {}
        for tag, src in (("nn1", nn1), ("family", fam), ("order", order_donor)):
            if src and src in functions:
                donors[tag] = sorted(norm_ec(functions[src]["ec"]) & universe)
        if "base_nn1" in args.arms and "nn1" in donors:
            arms["base_nn1"] = donors["nn1"]
        if "base_family_donor" in args.arms and "family" in donors:
            arms["base_family_donor"] = donors["family"]
        if "base_constant" in args.arms:
            if fold_key not in const_cache:
                const_cache[fold_key] = sorted(constant_arm(refs, functions, universe))
            arms["base_constant"] = const_cache[fold_key]
        jobs.append({"genome_id": gid, "fold": fold, "group_id": c["group_id"],
                     "region": c["region"], "md5": md5, "domain": c["domain"],
                     "template": c["template"], "arms": arms, "donors": donors,
                     "ec_universe": sorted(universe), "ladders": args.ladders,
                     "pipelines": args.pipelines, "readouts": args.readouts,
                     "ec_f1": ec_f1, "seed": args.seed,
                     "dump_models": args.dump_models,
                     "identity_control": False})
    if args.ladder_genomes >= 0:
        for i, j in enumerate(jobs):
            if i >= args.ladder_genomes:
                j["ladders"] = []
    return jobs


def job_covered(job, rows, args):
    """True iff `rows` already answer this invocation for this genome."""
    if not rows:
        return False
    want_r, want_l = set(job["readouts"]), set(job["ladders"])
    by = {}
    for r in rows:
        if str(r.get("status_reason", "")).startswith("build_exception"):
            return False                      # crashed last time -> retry the genome
        by[(r.get("arm"), r.get("pipeline"))] = r
    for arm in job["arms"]:
        for pipe in job["pipelines"]:
            r = by.get((arm, pipe))
            if r is None or not want_r <= set(r.get("readouts") or []):
                return False
    for fam in want_l:
        if not any(str(a).startswith("ladder_" + fam) for a, _p in by):
            return False
    return True


def pa_template_fallback(rec):
    """panel_eligible.json is the single source of `template` (spec 9.2); this fallback
    exists only so a hand-made fixture without the field still runs."""
    return pa().template_for(rec)


def region_rank(region, pref):
    return (0, "") if region == pref else (1, region or "")


def load_experimental(path):
    if not path or not os.path.exists(path):
        return set()
    with open(path) as fh:
        return {l.strip() for l in fh if l.strip() and not l.startswith("#")}


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", required=True)
    ap.add_argument("--fold-plan", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--prfba", required=True, help="root dir globbed for */predictions.json")
    ap.add_argument("--picrust2", required=True, help="root dir globbed for */predictions.json")
    ap.add_argument("--functions", default=os.path.join(
        REPO, "region_validation", "picrust2_comparison", "data",
        "genome_function_sets.json"),
        help="BV-BRC EC projection; default is the spec'd immutable input, so the "
             "SPEC section 11.1 smoke command (which omits the flag) runs verbatim")
    ap.add_argument("--outdir", required=True, help="per_genome.jsonl lives here")
    ap.add_argument("--n-genomes", type=int, default=200)
    ap.add_argument("--select-mode", default="exact_tie")
    ap.add_argument("--arms", default="truth,prfba,picrust2,base_nn1,base_family_donor,base_constant")
    ap.add_argument("--pipeline", default="as_shipped,controlled")
    ap.add_argument("--readouts", default="atp,growth,essentiality,structure")
    ap.add_argument("--experimental-list", default=None,
                    help="one genome_id per line; over-sampled into the subset")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--seed", type=int, default=1729)
    ap.add_argument("--min-truth-ec", type=int, default=200)
    # --- additive options (defaults reproduce the spec'd behaviour) ---
    ap.add_argument("--density", type=float, default=1.0,
                    help="which density arm of cells/predictions to read")
    ap.add_argument("--region-pref", default="V4",
                    help="preferred amplicon region per genome; else lexicographic first")
    ap.add_argument("--ec-universe", default=None,
                    help="panel_ec_universe.json (default: beside --panel)")
    ap.add_argument("--ladders", default=",".join(ALL_LADDERS),
                    help="dropout,false_positive,donor,size_matched_null or none")
    ap.add_argument("--ladder-genomes", type=int, default=-1,
                    help="run ladders on the first N sampled genomes (-1 = all)")
    ap.add_argument("--summary", default=None,
                    help="default: <dirname of --outdir>/phenotype_summary.json")
    ap.add_argument("--dump-models", default=None, help="dir for optional model JSON dumps")
    ap.add_argument("--msd", default=MSD_DEFAULT)
    ap.add_argument("--force", action="store_true", help="ignore and overwrite existing rows")
    ap.add_argument("--templates", default=None, metavar="DIR",
                    help="read the ModelSEED templates from DIR (a ModelSEEDTemplates "
                         "templates/v7.0/ checkout) instead of phenotype_arm.MST; SPEC "
                         "amendment A-C12 uses commit 21a8084")
    ap.add_argument("--integrate-reversals", action="store_true",
                    help="apply the gapfiller's reversed reactions as well as its new ones "
                         "(SPEC amendment A-C12; off reproduces the pre-registered runs)")
    ap.add_argument("--only-genomes", default=None, metavar="ID,ID",
                    help="after the usual stratified selection of --n-genomes, build only "
                         "these genome ids (targeted reruns and checks; the selection and "
                         "every per-genome seed are unchanged)")
    ap.add_argument("--expect-modelseedpy", default=None, metavar="COMMIT",
                    help="stop unless modelseedpy is imported from a git checkout at COMMIT "
                         "(a prefix will do) with no change to its tracked files; SPEC "
                         "amendment A-C12 pins d884a50, the pre-registered run's commit, "
                         "by putting a detached clone first on PYTHONPATH")
    ap.add_argument("--maxtasksperchild", type=int, default=None, metavar="N",
                    help="replace each worker after N genomes, which bounds the memory a "
                         "worker accumulates over a long pass (rows do not depend on which "
                         "worker builds a genome; default keeps workers for the whole pass)")
    args = ap.parse_args(argv)

    args.arms_raw = [a for a in args.arms.split(",") if a]
    args.prfba_arm = "prfba_{}".format(args.select_mode)
    args.arms = ["truth"] + [args.prfba_arm if a == "prfba" else a
                             for a in args.arms_raw if a != "truth"]
    seen, uniq = set(), []
    for a in args.arms:
        if a not in seen:
            seen.add(a)
            uniq.append(a)
    args.arms = uniq
    args.pipelines = [p for p in args.pipeline.split(",") if p]
    args.readouts = [r for r in args.readouts.split(",") if r]
    args.ladders = ([] if args.ladders in ("", "none")
                    else [l for l in args.ladders.split(",") if l in ALL_LADDERS])
    args.summary = args.summary or os.path.join(
        os.path.dirname(os.path.abspath(args.outdir.rstrip("/"))),
        "phenotype_summary.json")
    return args


def main(argv=None):
    args = parse_args(argv)
    libs = library_provenance()
    msp = libs.get("modelseedpy") or {}
    print("[pheno] modelseedpy {} at {} (tracked files changed: {})".format(
        msp.get("path"), msp.get("git_commit"), msp.get("git_dirty")), flush=True)
    if args.expect_modelseedpy and not (
            (msp.get("git_commit") or "").startswith(args.expect_modelseedpy)
            and msp.get("git_dirty") is False):
        raise SystemExit("[pheno] modelseedpy is not a clean checkout of {} (see above); "
                         "put the pinned clone first on PYTHONPATH".format(
                             args.expect_modelseedpy))
    os.makedirs(args.outdir, exist_ok=True)
    jsonl = os.path.join(args.outdir, "per_genome.jsonl")

    panel = jload(args.panel)
    fold_plan = jload(args.fold_plan)
    functions = jload(args.functions)
    _ec_path = args.ec_universe or os.path.join(
        os.path.dirname(os.path.abspath(args.panel)), "panel_ec_universe.json")
    universe = ec_universe(_ec_path, panel, functions)
    cells = load_cells(args.cells, args.density)
    prf, prf_meta = load_prfba(args.prfba, args.select_mode, args.density)
    p2, p2_meta = load_picrust2(args.picrust2, args.density)
    print("[pheno] universe={} cells={} prfba_md5={} picrust2_md5={}".format(
        len(universe), len(cells), len(prf), len(p2)), flush=True)

    # Provenance of the nn1 donor (the base_nn1 arm and the donor ladder's nn1 rung):
    # a fold run's rank1_genome is the alignment rank-1 hit against the fold-restricted
    # shared reference DB (run_prfba_fold.py: top20[0]), which is the intended nearest
    # surviving relative; adapt_as_shipped.py instead sets it to sorted(reps)[0], i.e.
    # prFBA's own first selected rep, which is not an independent neighbour.  Label it
    # rather than silently calling both "nn1".
    pf_src = {(mm.get("prediction_source") or mm.get("mode") or "fold")
              for mm in prf_meta}
    donor_nn1_source = (
        "prfba_reps_first: adapt_as_shipped.py sets rank1_genome = sorted(reps)[0], so "
        "the nn1 rung is 'copy prFBA's first rep', NOT an independent nearest neighbour"
        if "as_shipped" in pf_src else
        "alignment_rank1_hit: run_prfba_fold.py top20[0] against the fold-restricted "
        "shared reference DB (the intended nearest surviving relative)")

    template_dir = use_templates(args.templates)
    jobs = build_jobs(args, panel, fold_plan, functions, universe, cells, prf, p2)
    if args.only_genomes:
        keep = {g for g in args.only_genomes.split(",") if g}
        jobs = [j for j in jobs if j["genome_id"] in keep]
    if not jobs:
        print("[pheno] no eligible genomes -- nothing to do", flush=True)

    done = {}
    if os.path.exists(jsonl) and not args.force:
        with open(jsonl) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                done.setdefault(r["src_genome_id"], []).append(r)
    elif args.force and os.path.exists(jsonl):
        os.remove(jsonl)
    # A genome counts as done only if its stored rows COVER what this invocation asks
    # for: the orchestrator runs this arm twice into the same per_genome.jsonl (ATP over
    # 1,000 genomes, then atp+growth+essentiality+structure over 200), and a
    # presence-only check would silently skip the second, heavier pass.
    todo = []
    for j in jobs:
        if not job_covered(j, done.get(j["genome_id"], []), args):
            done.pop(j["genome_id"], None)     # stale rows are replaced, not duplicated
            todo.append(j)
    if todo:
        todo[0]["identity_control"] = True          # C11, inside a real worker batch
    print("[pheno] {} genomes selected, {} already done, {} to build".format(
        len(jobs), len(jobs) - len(todo), len(todo)), flush=True)

    rows = [r for rs in done.values() for r in rs]
    ident_path = os.path.join(args.outdir, "identity_control.json")
    identity = (jload(ident_path) if (os.path.exists(ident_path) and not args.force)
                else None)

    if todo:
        if not os.path.isdir(args.msd):
            raise SystemExit("[pheno] --msd is not a directory: {} (ModelSEEDDatabase "
                             "is required by MSATPCorrection)".format(args.msd))
        import multiprocessing as mp
        ctx = mp.get_context("spawn")               # Gurobi envs do not survive fork
        nw = max(1, min(args.workers, len(todo)))
        fh = open(jsonl, "a")
        try:
            winit = (args.msd, True, args.templates, args.integrate_reversals)
            if nw == 1:
                worker_init(*winit)
                results = (run_genome_safe(j) for j in todo)
            else:
                pool = ctx.Pool(nw, initializer=worker_init, initargs=winit,
                                maxtasksperchild=args.maxtasksperchild)
                results = pool.imap_unordered(run_genome_safe, todo)
            for res in results:
                if res.get("identity_control"):
                    identity = res["identity_control"]
                    with open(ident_path, "w") as ifh:
                        json.dump(identity, ifh, indent=1, sort_keys=True)
                for r in res["rows"]:
                    fh.write(json.dumps(r, sort_keys=True) + "\n")
                fh.flush()
                rows.extend(res["rows"])
                print("[pheno] {} done in {}s ({} rows)".format(
                    res["genome_id"], res["wall_seconds"], len(res["rows"])), flush=True)
        finally:
            fh.close()
            if nw > 1:
                pool.close()
                pool.join()

    # deterministic artefact regardless of worker scheduling
    rows.sort(key=lambda r: (r["src_genome_id"], r["arm"], r["pipeline"]))
    with open(jsonl, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r, sort_keys=True) + "\n")

    frozen = (panel.get("meta", {}).get("frozen")
              or fold_plan.get("meta", {}).get("frozen") or {})
    meta_extra = {
        "identity_control": identity or {"jaccard": None, "pass": None},
        "atp_pipeline_independent": True,
        "select_mode": args.select_mode, "prfba_arm": args.prfba_arm,
        "density": args.density, "region_pref": args.region_pref,
        "readouts": args.readouts, "pipelines": args.pipelines,
        "arms_requested": args.arms, "ladders_requested": args.ladders,
        "min_truth_ec": args.min_truth_ec, "n_ec_universe": len(universe),
        "n_candidates_eligible": len(jobs),
        "prfba_runs": prf_meta, "picrust2_runs": p2_meta,
        "donor_nn1_source": donor_nn1_source,
        "templates": template_provenance(template_dir),
        "libraries": libs,
        "integrate_reversals": bool(args.integrate_reversals),
        "only_genomes": (sorted({g for g in args.only_genomes.split(",") if g})
                         if args.only_genomes else None),
        # the only thing that can tie this file to a run (score_folds C11/C12)
        "binding": binding_block(args, rows, jsonl, ident_path, _ec_path),
        "script": "phenotype_fold.py"}
    summary = summarize(rows, args.arms, args.pipelines, frozen, meta_extra,
                        args.ladders)
    os.makedirs(os.path.dirname(os.path.abspath(args.summary)), exist_ok=True)
    with open(args.summary, "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
    print("[pheno] wrote {} ({} rows) and {}".format(jsonl, len(rows), args.summary),
          flush=True)
    ic = summary["meta"]["identity_control"]
    if ic.get("pass") is False:
        print("[pheno] WARNING C11 identity control FAILED: {}".format(ic), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
