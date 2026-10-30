#!/usr/bin/env python
"""In-silico test of prFBA's weighted objective (WEIGHTED_OBJECTIVE_TEST.md, 2026-10-29).

QUESTION
  prFBA minimizes  sum_r max(min_prob, 1 - p_r^prob_exp) |v_r|  +  sum_ex ex_weight |v_ex|
  with biomass >= 95% of its maximum (prFBA.tex eq:objective; MSProbability.prFBA).
  Where the selection holds n >= 2 genomes, p is fractional.  Does the gradation move
  flux off reactions the source organism lacks, and does it shrink the optimal face
  (degeneracy)?  The paper states this has not been measured.

ONE AMPLICON, ONE WORKER (run_amplicon)
  1. Gapfilled model of the selection's UNION EC set, built exactly like the phenotype
     harness's gapfilled pass: phenotype_fold.gapfill_block (as_shipped pipeline, complete
     medium, gapfiller reversals applied), templates of the truth genome's lineage
     (panel_eligible.json:template), ModelSEEDpy pinned by PYTHONPATH.
  2. Per-role probability, the EC-route analogue of the production merge
     (bvbrc_to_kbase_genome.create_synthetic_genome: a function's counter is incremented
     once per selected genome carrying it, then divided by the number of genomes):
        p(role) = #{members g : every EC in the role's name is in EC(g)} / n
     The harness's features ARE template roles (genome_from_ec), so a model gene id is a
     role id; that is checked per amplicon.
  3. Per-reaction probability by the production rule (build_prfba_community_models.py
     _build_one; load_genomes_to_modelseedpy.py):  p(R) = max over R's genes of p(gene),
     0.0 when R has no genes, stored as str in rxn.notes['probability'].  Applied AFTER
     gapfilling, so gapfilled reactions get 0.0 -- the value the rule gives every gene-less
     reaction, and the paper's 'for example a weight of 0'.
  4. Truth support: the truth genome's own base model (same template, same
     build_base_model arguments, no gapfilling); a reaction is truth-supported iff it is
     gene-associated there.
  5. Arms (all through prfba_released_copy, each on a fresh copy of the gapfilled model):
        real        production p                              (prob_exp 1, ex_weight 100)
        two_tier    gene-backed p = 1, gene-less 0            what production gives at n = 1
        uniform     every rxn p = 1 (plain parsimony + exchange term)
        shuffled    production p permuted among the gene-backed reactions (seeded)
        oracle      gene-backed p = 1 if truth-supported else 0 (positive control)
        real_exp2   production p, prob_exp 2                  (sensitivity)
        real_ex1 / two_tier_ex1   ex_weight 1                 (exploratory)
  6. Primary, vertex-robust: on the numerically optimal face (objective <= opt + 1e-6, ONE
     absolute slack for every arm) the min and max of the total truth-absent flux, each
     divided by an arm-independent per-amplicon constant D (the two-tier arm's gene-backed
     flux at its returned optimum): taf_face_min / taf_face_mid / taf_face_max.
     Secondary, 'as returned': the point metrics at the solver's vertex, with their spread
     over Gurobi Method 0/1/2 re-solves of the same LP.  Degeneracy (four arms): FVA of the
     gene-backed reactions on the optimal face, confirmed on a 1e-7 face; descriptive
     near-optimal faces (int1e-3, obj1e-3).  WEIGHTED_OBJECTIVE_TEST.md 'Review' says why.

THE OBJECTIVE CODE
  prfba_released_copy() is MSProbability.prFBA (ModelSEEDpy c2e31a0, msprobability.py
  l.155-201) with only the edits needed for the single-model path to run.  The released
  function cannot run on one model: `len(model_s_)` raises TypeError on a cobra Model and
  MSModelUtil(list) raises AttributeError; with that unwrapped, the second
  ElementUptakePkg (l.174) raises PackageError.  Every bound, constraint and coefficient
  is verbatim.  Each row also re-solves the real arm through the RELEASED SOURCE itself
  (exec'd from the imported module with the two textual edits, asserted unique) and
  records whether the objective agrees (released_check).

Interpreter: /home/freiburger/Documents/py_venv/bin/python (the only interpreter).
PYTHONPATH must put the pinned ModelSEEDpy clone first; --expect-modelseedpy enforces it.
"""
import argparse
import contextlib
import gzip
import hashlib
import inspect
import io
import json
import math
import os
import random
import subprocess
import sys
import textwrap
import time
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)
D = os.path.join(SR, "data")
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import phenotype_fold as PF          # noqa: E402  (no cobra at import time)

DTAG = {1.0: "d100", 0.5: "d050", 0.25: "d025"}
STRATA = ("n2", "n3", "n4+", "n1")          # job interleave order; n1 = degenerate control
MULTI = ("n2", "n3", "n4+")
ARMS = ("two_tier", "real", "uniform", "shuffled", "oracle", "real_exp2",
        "real_ex1", "two_tier_ex1")       # two_tier first: it sets the common slack
ACTIVE = 1e-6                                 # |v| above this = active
SESOI_TAF = 0.01                              # declared in WEIGHTED_OBJECTIVE_TEST.md
SESOI_DEGEN_REL = 0.10
OPT_SLACK = 1e-6                              # the optimal face: one absolute slack, all arms
OPT7_SLACK = 1e-7                             # confirmation face of the degeneracy count
CONFIRM_REL = 0.05                            # flag: opt vs opt7 count differ by more than this
LP_METHODS = (0, 1, 2)                        # Gurobi primal simplex, dual simplex, barrier
EX100_ARMS = ("two_tier", "real", "uniform", "shuffled", "oracle", "real_exp2")
TEMPLATES_DEFAULT = "/home/freiburger/Documents/ModelSEEDTemplates-c12run/templates/v7.0/"


def stratum_of(n):
    return "n1" if n == 1 else "n2" if n == 2 else "n3" if n == 3 else "n4+"


def now_iso():
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit(path):
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=path, check=True,
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        return None


# ============================================================ the objective (released)
def prfba_released_copy(model, environment=None, min_prob=0.01, prob_exp=1, ex_weight=100):
    """MSProbability.prFBA, ModelSEEDpy c2e31a0 modelseedpy/core/msprobability.py l.155-201,
    single-model path.  Differences from the released text -- each either required for the
    function to run on one model at all, or output-only:
      (1) MSModelUtil(model) for `MSModelUtil(model_s_ if len(model_s_) == 1 else ...)`:
          a cobra Model has no len() (TypeError) and MSModelUtil([model]) fails on
          `self.id = model.id` (AttributeError).
      (2) the second `ElementUptakePkg(mdlUtil.model).build_package({"C": 100})` (l.174)
          is dropped: it raises PackageError("Package with name ElementUptakePkg already
          in model!").  The first build (l.168) is kept verbatim.
      (3) `print([cons.name ...])` (l.194) and the printLP branch are dropped (stdout/file
          only); the mscommunity imports and the MSCommunity rows are community-only.
    Every bound, constraint and objective coefficient is the released code's.  Like the
    released code it maximizes WHATEVER objective the model carries to set the 95% floor;
    the caller sets bio1 (a gapfilled model as built carries ATPM_c0).
    Returns (solution, MSModelUtil, coef, maxBioSol); the model is modified in place."""
    from modelseedpy.core.msmodelutl import MSModelUtil
    from modelseedpy.fbapkg.elementuptakepkg import ElementUptakePkg
    from optlang.symbolics import Zero

    mdlUtil = MSModelUtil(model)                                            # edit (1)
    mdlUtil.add_medium(environment)
    # constrain carbon consumption and community composition
    elepkg = ElementUptakePkg(mdlUtil.model)  ;  elepkg.build_package({"C": 100})
    # (edit 2: the duplicate ElementUptakePkg build is not repeated)
    # constrain the model to 95% of the optimum growth
    maxBioSol = mdlUtil.model.slim_optimize()
    mdlUtil.add_minimal_objective_cons(maxBioSol*.95)
    # weight internal reactions based on their probabilities
    coef = {}
    for rxn in mdlUtil.model.reactions:
        if "rxn" == rxn.id[0:3]:
            coef.update({rxn.forward_variable: max(min_prob, (1-float(rxn.notes["probability"])**prob_exp))})
            coef.update({rxn.reverse_variable: max(min_prob, (1-float(rxn.notes["probability"])**prob_exp))})
        elif "EX_" == rxn.id[0:3]:
            coef.update({rxn.forward_variable: ex_weight})
            coef.update({rxn.reverse_variable: ex_weight})
    mdlUtil.add_objective(Zero, "min", coef)
    # (edit 3: no constraint-name print, no LP file)
    return mdlUtil.model.optimize(), mdlUtil, coef, maxBioSol


_RELEASED = {}
RELEASED_EDITS = (
    ("mdlUtil = MSModelUtil(model_s_ if len(model_s_) == 1 else",
     "model_s_ = model_s_ if isinstance(model_s_, (list, tuple)) else [model_s_]\n"
     "    mdlUtil = MSModelUtil(model_s_[0] if len(model_s_) == 1 else"),
    ("    # constrain carbon consumption\n"
     "    elepkg = ElementUptakePkg(mdlUtil.model)  ;  elepkg.build_package({\"C\": 100})\n",
     ""),
)


def released_prfba():
    """The imported module's own prFBA source with the two edits above (each asserted to
    occur exactly once), exec'd in the module's namespace.  Used per row as the check that
    prfba_released_copy() is the released objective."""
    if "fn" not in _RELEASED:
        import modelseedpy.core.msprobability as MP
        src = textwrap.dedent(inspect.getsource(MP.MSProbability.prFBA)).replace(
            "@staticmethod\n", "")
        _RELEASED["sha256"] = hashlib.sha256(src.encode()).hexdigest()
        for old, new in RELEASED_EDITS:
            if src.count(old) != 1:
                raise RuntimeError("released prFBA source changed: edit target found {} "
                                   "times: {!r}".format(src.count(old), old[:60]))
            src = src.replace(old, new)
        ns = dict(vars(MP))
        exec(compile(src, "<released prFBA + edits 1,2>", "exec"), ns)
        _RELEASED["fn"] = ns["prFBA"]
    return _RELEASED["fn"]


def purge_caches():
    """MSModelUtil.mdlutls and MSPackageManager.pkgmgrs are class-level dicts keyed by the
    model object: every copy would otherwise stay alive for the life of the worker."""
    from modelseedpy.core.msmodelutl import MSModelUtil
    from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
    MSModelUtil.mdlutls.clear()
    MSPackageManager.pkgmgrs.clear()


# ============================================================ population
def load_population(args):
    """Scored EC cells of prfba_<select_mode> at --density -> one candidate per
    (source genome, stratum), the selection read from the fold's predictions.json."""
    import numpy as np  # noqa: F401
    import orjson
    import pandas as pd
    t0 = time.time()
    panel = PF.jload(args.panel)["genomes"]
    uni = PF.norm_ec(PF.jload(args.ec_universe)["ec"])
    pec = PF.jload(args.panel_ec)
    arm = "prfba_" + args.select_mode
    usecols = ["amplicon", "md5", "region", "src_genome_id", "fold", "group_id", "density",
               "arm", "space", "outcome", "n_truth", "f1", "n_reps", "novelty_bin",
               "rank1_identity"]
    parts = []
    for ch in pd.read_csv(args.cells, dtype=str, usecols=usecols, keep_default_na=False,
                          chunksize=500_000):
        ch = ch[(ch.arm == arm) & (ch.space == "ec") & (ch.outcome == "scored")]
        ch = ch[ch.density.astype(float).sub(args.density).abs() < 1e-9]
        parts.append(ch)
    cells = pd.concat(parts, ignore_index=True)
    A = Counter()
    A["scored_cells"] = len(cells)
    md5_srcs = cells.groupby("md5").src_genome_id.nunique()
    twin = set(md5_srcs[md5_srcs > 1].index)
    A["twin_md5_scored"] = len(twin)
    by_fold = {f: g for f, g in cells.groupby("fold")}
    cands = []
    for fold in sorted(by_fold, key=int):
        sub = by_fold[fold]
        path = os.path.join(args.prfba, "f{:02d}_{}".format(int(fold), DTAG[args.density]),
                            "predictions.json")
        with open(path, "rb") as fh:
            P = orjson.loads(fh.read())
        m = P["meta"]
        if m.get("mode") not in PF.PRFBA_MODES_OK or bool(m.get("forced_answer")) or \
                abs(float(m.get("density", args.density)) - args.density) > 1e-9:
            raise SystemExit("[wobj] {} is not a held-out fold run at density {} "
                             "(mode={}, forced={})".format(path, args.density, m.get("mode"),
                                                           m.get("forced_answer")))
        block = P["pred"][args.select_mode]
        del P
        for row in sub.itertuples(index=False):
            rec = block.get(row.md5)
            if rec is None:
                A["md5_missing_in_predictions"] += 1
                continue
            reps = list(rec["reps"])
            n = len(reps)
            if n == 0 or n != int(row.n_reps):
                A["n_reps_mismatch_or_empty"] += 1
                continue
            src = row.src_genome_id
            if src not in panel or src not in pec or any(r not in pec for r in reps):
                A["genome_outside_panel_ec"] += 1
                continue
            if src in reps:
                A["leak_source_in_reps"] += 1
                continue
            truth = PF.norm_ec(pec[src]["ec"]) & uni
            if len(truth) < args.min_truth_ec:
                A["truth_below_min_ec"] += 1
                continue
            member = {r: sorted(PF.norm_ec(pec[r]["ec"]) & uni) for r in reps}
            union = set().union(*map(set, member.values()))
            pred_ec = PF.norm_ec(rec.get("ec")) & uni
            f1 = PF.set_f1(union, truth)
            cands.append({
                "key": "{}__{}".format(src, row.region), "src_genome_id": src,
                "region": row.region, "md5": row.md5, "fold": int(fold),
                "group_id": row.group_id, "template": panel[src].get("template")
                or PF.pa_template_fallback(panel[src]),
                "genus": panel[src].get("genus"), "n_reps": n, "stratum": stratum_of(n),
                "reps": reps, "twin_md5": row.md5 in twin,
                "novelty_bin": row.novelty_bin, "rank1_identity": row.rank1_identity,
                "tier": rec.get("tier"),
                "checks": {"union_eq_pred_ec": union == pred_ec,
                           "f1_cells": float(row.f1), "f1_recomputed": round(f1, 6),
                           "f1_agrees": abs(f1 - float(row.f1)) <= 1e-3,
                           "n_truth_cells": int(row.n_truth), "n_truth": len(truth),
                           "n_truth_agrees": int(row.n_truth) == len(truth)},
                "_member": member, "_truth": sorted(truth)})
        del block
        print("[wobj] fold {}: {} candidate cells so far ({:.0f}s)".format(
            fold, len(cands), time.time() - t0), flush=True)
    return cands, A


def sample_population(cands, args):
    """One amplicon per source genome and per md5 across the whole sample; strata drawn
    rarest first (n4+, n3, n2, n1); within a stratum a seeded shuffle of the sorted source
    genomes, each source contributing its --region-pref cell (else lexicographic first);
    twin md5s excluded (an identical model would be scored against two truths)."""
    sizes = args.n_per_stratum
    by = defaultdict(lambda: defaultdict(list))
    elig = Counter()
    for c in cands:
        if c["twin_md5"]:
            continue
        by[c["stratum"]][c["src_genome_id"]].append(c)
    for s in by:
        elig[s] = len(by[s])                 # distinct source genomes in the stratum
    used_src, used_md5, chosen = set(), set(), {s: [] for s in STRATA}
    for s in ("n4+", "n3", "n2", "n1"):
        srcs = sorted(by.get(s, {}))
        random.Random("{}:wobj:{}".format(args.seed, s)).shuffle(srcs)
        for src in srcs:
            if len(chosen[s]) >= sizes.get(s, 0):
                break
            if src in used_src:
                continue
            opts = sorted(by[s][src], key=lambda c: PF.region_rank(c["region"],
                                                                 args.region_pref))
            pick = next((c for c in opts if c["md5"] not in used_md5), None)
            if pick is None:
                continue
            used_src.add(src)
            used_md5.add(pick["md5"])
            chosen[s].append(pick)
    jobs = []
    for i in range(max(len(v) for v in chosen.values())):
        for s in STRATA:
            if i < len(chosen[s]):
                jobs.append(chosen[s][i])
    return jobs, dict(elig)


def build_or_load_population(args):
    path = os.path.join(args.outdir, "population.json")
    inputs = {k: {"path": os.path.abspath(p), "sha256": sha256(p)} for k, p in (
        ("cells", args.cells), ("panel", args.panel), ("panel_ec", args.panel_ec),
        ("ec_universe", args.ec_universe), ("fold_plan", args.fold_plan))}
    want = {"select_mode": args.select_mode, "density": args.density,
            "n_per_stratum": args.n_per_stratum, "seed": args.seed,
            "min_truth_ec": args.min_truth_ec, "region_pref": args.region_pref,
            "inputs": inputs}
    if os.path.exists(path):
        pop = PF.jload(path)
        if pop["meta"]["spec"] != json.loads(json.dumps(want)):
            raise SystemExit("[wobj] {} was built from different inputs or settings; use a "
                             "new --outdir".format(path))
        print("[wobj] population reloaded: {} jobs".format(len(pop["jobs"])), flush=True)
        return pop
    cands, A = load_population(args)
    jobs, elig = sample_population(cands, args)
    allc = Counter(c["stratum"] for c in cands)
    meta = {"spec": want, "generated": now_iso(), "git_commit": git_commit(REPO),
            "candidate_cells_by_stratum": dict(allc),
            "eligible_source_genomes_by_stratum_excl_twin": elig,
            "sampled_by_stratum": dict(Counter(j["stratum"] for j in jobs)),
            "assertions": dict(A),
            "candidate_checks": {
                "union_ne_pred_ec": sum(1 for c in cands if not c["checks"]["union_eq_pred_ec"]),
                "f1_disagrees": sum(1 for c in cands if not c["checks"]["f1_agrees"]),
                "n_truth_disagrees": sum(1 for c in cands if not c["checks"]["n_truth_agrees"])},
            "n_candidates": len(cands)}
    out = []
    for j in jobs:
        j = dict(j)
        j["member_ec"] = j.pop("_member")
        j["truth_ec"] = j.pop("_truth")
        out.append(j)
    pop = {"meta": meta, "jobs": out}
    os.makedirs(args.outdir, exist_ok=True)
    with open(path + ".tmp", "w") as fh:
        json.dump(pop, fh, sort_keys=True)
    os.replace(path + ".tmp", path)
    print("[wobj] population: {} candidates, sampled {} ({})".format(
        len(cands), len(out), meta["sampled_by_stratum"]), flush=True)
    return pop


# ============================================================ worker
_W = {}


def worker_init(msd, template_dir, integrate_reversals, quiet=True):
    PF.worker_init(msd, quiet, template_dir, integrate_reversals)
    m = PF.pa()
    _W["gf_media"], _W["gf_cpds"] = m.complete_media()
    # the environment passed to prFBA: the gapfilling medium's compounds at its uptake
    # bound (100), as a {cpd: uptake} dict (MSModelUtil.add_medium maps it to EX_*_e0).
    # An MSMedia object would instead be converted with FBAHelper.convert_kbase_media's
    # default uniform_uptake=1000.
    _W["env"] = {c: 100.0 for c in _W["gf_cpds"]}
    released_prfba()


def role_probabilities(gene_ids, role_ecs, member_ec):
    """p(role) = share of members whose own EC set grants the role (all ECs of its name)."""
    n = len(member_ec)
    sets = [set(v) for v in member_ec.values()]
    p, unknown = {}, []
    for g in gene_ids:
        rec = role_ecs.get(g)
        if rec is None:
            unknown.append(g)
            continue
        ecs = rec[1]
        p[g] = sum(1 for s in sets if ecs <= s) / n
    return p, unknown


def arm_probs(arm, rxns, prod, gene_backed, supported, seed, key):
    """{rxn_id: probability} for the rxn-prefixed reactions under one arm."""
    if arm in ("real", "real_exp2", "real_ex1"):
        return dict(prod)
    if arm in ("two_tier", "two_tier_ex1"):
        return {r: (1.0 if r in gene_backed else 0.0) for r in rxns}
    if arm == "uniform":
        return {r: 1.0 for r in rxns}
    if arm == "oracle":
        return {r: ((1.0 if r in supported else 0.0) if r in gene_backed else 0.0)
                for r in rxns}
    if arm == "shuffled":
        gb = sorted(gene_backed)
        vals = [prod[r] for r in gb]
        PF.rng_for(seed, key, "shuffle").shuffle(vals)
        out = dict(prod)
        out.update(zip(gb, vals))
        return out
    raise ValueError(arm)


# Which faces each arm is measured on (WEIGHTED_OBJECTIVE_TEST.md section 8 and 'Review').
#   opt      objective <= opt + OPT_SLACK (1e-6 absolute, the same for every arm): the
#            numerically optimal face, i.e. the alternative optima.  Every arm: min/max of
#            the truth-absent flux (2 LPs; the vertex-robust primary).  FVA_ARMS also: FVA of
#            the scan-variable gene-backed reactions (the degeneracy primary, n_var_1e-3).
#   opt7     objective <= opt + 1e-7: confirmation of the opt count (FVA_ARMS).  The review
#            found the count plateaus between slacks 1e-8 and 1e-6 and inflates 3-7x above.
#   int1e-3  ONE absolute slack per amplicon, 0.1% of the two-tier arm's internal term:
#            near-optimal, descriptive (sum of ranges).  Its truth-absent interval is stored
#            but not read: at slack s the LP can add futile forward+reverse flux on the
#            cheapest truth-absent reaction (coefficient c_min), i.e. up to s / c_min of
#            'absent flux' (futile_offset_bound); on the smoke that was most of the width.
#   obj1e-3  0.1% of the arm's whole objective (real only): the diagnostic that a tolerance
#            on the objective is swamped by the exchange term.
FVA_ARMS = ("two_tier", "real", "uniform", "shuffled")
FACES = {"real": ("opt", "opt7", "int1e-3", "obj1e-3"), "two_tier": ("opt", "opt7", "int1e-3"),
         "uniform": ("opt", "opt7", "int1e-3"), "shuffled": ("opt", "opt7", "int1e-3"),
         "oracle": ("opt",), "real_exp2": ("opt",), "real_ex1": ("opt",),
         "two_tier_ex1": ("opt",)}

ARM_PARAMS = {"real": (1, 100), "two_tier": (1, 100), "uniform": (1, 100),
              "shuffled": (1, 100), "oracle": (1, 100), "real_exp2": (2, 100),
              "real_ex1": (1, 1), "two_tier_ex1": (1, 1)}


def set_probs(model, probs):
    for r in model.reactions:
        if r.id[:3] == "rxn":
            r.notes["probability"] = str(probs[r.id])


def mag(model):
    """|v| per reaction as forward + reverse primal (the quantity the objective prices)."""
    return {r.id: max(0.0, r.forward_variable.primal) + max(0.0, r.reverse_variable.primal)
            for r in model.reactions}


def feasible_scan(M, targets):
    """Min/max net flux of each target over the prFBA constraint set alone (medium,
    carbon cap, biomass floor; no objective bound).  A reaction whose range there is 0 has
    range 0 on every optimal face of every arm (all arms share the constraint set), so the
    per-arm FVA skips it.  Returns ({id: (min, max)}, n_failed)."""
    from optlang.symbolics import Zero
    nan = float("nan")
    lo, hi, fail = {}, {}, 0
    with M:
        M.objective = M.problem.Objective(Zero, direction="min")
        M.solver.update()
        for direction, store in (("min", lo), ("max", hi)):
            M.solver.objective.direction = direction
            for r in targets:
                M.solver.objective.set_linear_coefficients(
                    {r.forward_variable: 1, r.reverse_variable: -1})
                v = M.slim_optimize(error_value=nan)
                fail += int(v != v)
                store[r.id] = v
                M.solver.objective.set_linear_coefficients(
                    {r.forward_variable: 0, r.reverse_variable: 0})
    return {k: (lo[k], hi[k]) for k in lo}, fail


def _solve_retry(M, nan):
    """slim_optimize; on failure, one cold retry (gurobipy Model.reset).  The tight faces
    (slack 1e-6, 1e-7) occasionally fail one LP from a warm start (first amended smoke: 4 of
    about 2 x 10^4 face LPs).  Returns (value, status of the first failure or None)."""
    v = M.slim_optimize(error_value=nan)
    if v == v:
        return v, None
    st = str(M.solver.status)
    try:
        M.solver.problem.reset()
    except Exception:                                 # noqa: BLE001
        pass
    return M.slim_optimize(error_value=nan), st


def face_block(M, opt, slack, targets, intervals, fva=True):
    """On {objective <= opt + slack}: FVA of net flux over `targets` (when `fva`), and the
    min and max of the total flux (fwd + rev) through each reaction set in `intervals`
    ({name: [reactions]}, stored as <name>_flux_min / <name>_flux_max).  Restores the model
    on exit."""
    from optlang.symbolics import Zero
    nan = float("nan")
    out = {"slack": slack}
    with M:
        c = M.problem.Constraint(M.solver.objective.expression, ub=opt + slack,
                                 name="wobj_face")
        M.add_cons_vars([c])
        M.objective = M.problem.Objective(Zero, direction="min")
        M.solver.update()
        lo, hi, fail, retried, first_status = {}, {}, 0, 0, Counter()
        for direction, store in (("min", lo), ("max", hi)):
            M.solver.objective.direction = direction
            for r in (targets if fva else ()):
                M.solver.objective.set_linear_coefficients(
                    {r.forward_variable: 1, r.reverse_variable: -1})
                v, st = _solve_retry(M, nan)
                if st is not None:
                    retried += 1
                    first_status[st] += 1
                if v != v:
                    fail += 1
                store[r.id] = v
                M.solver.objective.set_linear_coefficients(
                    {r.forward_variable: 0, r.reverse_variable: 0})
        iv, ifail = {}, 0
        for name, rset in intervals.items():
            coefs = {}
            for r in rset:
                coefs[r.forward_variable] = 1
                coefs[r.reverse_variable] = 1
            if not coefs:
                iv[name] = {"min": 0.0, "max": 0.0}
                continue
            M.solver.objective.set_linear_coefficients(coefs)
            d = {}
            for direction in ("min", "max"):
                M.solver.objective.direction = direction
                d[direction], st = _solve_retry(M, nan)
                if st is not None:
                    retried += 1
                    first_status[st] += 1
                ifail += int(d[direction] != d[direction])
            M.solver.objective.set_linear_coefficients({k: 0 for k in coefs})
            iv[name] = d
    rng = {k: hi[k] - lo[k] for k in lo if lo[k] == lo[k] and hi[k] == hi[k]}
    out.update({"n_targets": len(targets) if fva else 0, "n_lp_failed": fail,
                "n_interval_lp_failed": ifail, "n_lp_retried": retried,
                "lp_first_fail_status": dict(first_status),
                "n_var_1e-6": sum(1 for x in rng.values() if x > 1e-6) if fva else None,
                "n_var_1e-3": sum(1 for x in rng.values() if x > 1e-3) if fva else None,
                "sum_range": round(sum(rng.values()), 6) if fva else None,
                "_range": rng})
    for name, d in iv.items():
        out["{}_flux_min".format(name)] = d["min"]
        out["{}_flux_max".format(name)] = d["max"]
    return out


def vertex_resolves(M, gb, absn, D):
    """The arm's LP re-solved from scratch (gurobipy Model.reset) with Method 0 (primal
    simplex), 1 (dual simplex) and 2 (barrier with crossover): the same optimum, possibly at
    another optimal vertex.  Returns the point measures at each.  Leaves Method -1
    (automatic) set and the LP re-solved before returning."""
    nan = float("nan")
    try:
        import gurobipy
        gp = M.solver.problem
        if not isinstance(gp, gurobipy.Model):
            return {"skipped": "solver is not gurobi"}
    except Exception as exc:                          # noqa: BLE001
        return {"skipped": "{}: {}".format(type(exc).__name__, exc)}
    out = {}
    try:
        for meth in LP_METHODS:
            gp.reset()
            gp.setParam("Method", meth)
            obj = M.slim_optimize(error_value=nan)
            if obj != obj:
                out[str(meth)] = {"status": str(M.solver.status)}
                continue
            v = mag(M)
            fgb = sum(v[i] for i in gb)
            fab = sum(v[i] for i in absn)
            out[str(meth)] = {"objective": obj, "flux_gene_backed": fgb, "flux_absent": fab,
                              "taf_gb": (fab / fgb) if fgb > 0 else None,
                              "taf_point_D": (fab / D) if D else None}
    finally:
        gp.reset()
        gp.setParam("Method", -1)
        M.slim_optimize(error_value=nan)
    return out


def _spread(xs):
    xs = [x for x in xs if x is not None and x == x]
    return (max(xs) - min(xs)) if len(xs) > 1 else None


def arm_metrics(M, sol, coef, maxbio, cls, prod, probs, args, face_tags, row, arm):
    """Everything measured on one arm's optimum."""
    v = mag(M)
    rx = [r for r in M.reactions if r.id[:3] == "rxn"]
    tot = lambda ids: sum(v[i] for i in ids)
    gb, sup, absn, unann = cls["gene_backed"], cls["supported"], cls["absent"], cls["unann"]
    gapf = cls["gapfilled"]
    unc = {r for r in gb if prod[r] < 1.0}           # fractional production p
    flux_rxn = tot(r.id for r in rx)
    flux_gb = tot(gb)
    flux_abs = tot(absn)
    flux_sup = tot(i for i in (r.id for r in rx) if i in sup)
    flux_unann = tot(unann)
    if arm == "two_tier":
        # D, the per-amplicon normaliser of the primary: it does not depend on the arm
        # (two_tier runs first; see run_amplicon)
        cls["ref_flux_gb"] = flux_gb
        row["model"]["D_flux_gb_two_tier"] = flux_gb
    D = cls.get("ref_flux_gb")
    pv = {var.name: var.primal for var in M.variables}
    internal = exch = 0.0
    split = Counter()
    for var, c in coef.items():
        x = max(0.0, pv[var.name]) * c
        rid = var.name.split("_reverse_")[0] if "_reverse_" in var.name else var.name
        if rid.startswith("EX_"):
            exch += x
        else:
            internal += x
            split["gene_backed_p1" if (rid in gb and probs.get(rid, 0) >= 1.0) else
                  "gene_backed_plt1" if rid in gb else
                  "gapfilled" if rid in gapf else "geneless_draft"] += x
    act = {r.id for r in rx if v[r.id] > ACTIVE}
    act_gb = act & gb
    act_abs = act & absn
    res = {
        "status": sol.status, "objective": sol.objective_value,
        "max_biomass": maxbio, "biomass": sol.fluxes.get("bio1"),
        "internal_term": internal, "exchange_term": exch,
        "internal_share": (internal / (internal + exch)) if (internal + exch) > 0 else None,
        "internal_term_split": dict(split),
        "flux_rxn": flux_rxn, "flux_gene_backed": flux_gb, "flux_absent": flux_abs,
        "flux_supported": flux_sup, "flux_unannotated": flux_unann,
        "flux_gapfilled": tot(gapf),
        "taf_gb": (flux_abs / flux_gb) if flux_gb > 0 else None,
        "taf_point_D": (flux_abs / D) if D else None,
        "taf_all": (flux_abs / flux_rxn) if flux_rxn > 0 else None,
        # gene-backed on both sides (precision_flux = 1 - taf_gb): gene-less truth-supported
        # reactions (gapfilled, spontaneous) are on neither side
        "precision_flux": ((flux_gb - flux_abs) / flux_gb) if flux_gb > 0 else None,
        "precision_count": (len(act_gb - act_abs) / len(act_gb)) if act_gb else None,
        "n_active": len(act), "n_active_gb": len(act_gb), "n_active_absent": len(act_abs),
        "n_active_gapfilled": len(act & gapf),
        "flux_plt1": tot(unc), "flux_plt1_absent": tot(unc & absn),
        "flux_p1_absent": tot(absn - unc),
        "plt1_share_gb": (tot(unc) / flux_gb) if flux_gb > 0 else None,
    }
    res["_v"] = {k: round(x, 9) for k, x in v.items() if x > ACTIVE}
    # secondary, 'as returned': the same LP at other optimal vertices.  After the point
    # metrics (C1 compares the default vertex), before the scan and the face LPs.
    t1 = time.time()
    vr = vertex_resolves(M, gb, absn, D)
    meth = {k: x for k, x in vr.items() if isinstance(x, dict) and "objective" in x}
    res["vertex_resolves"] = vr
    res["n_vertex_resolves"] = len(meth)
    res["taf_gb_vertex_spread"] = _spread([res["taf_gb"]] + [x["taf_gb"] for x in meth.values()])
    res["taf_point_D_vertex_spread"] = _spread(
        [res["taf_point_D"]] + [x["taf_point_D"] for x in meth.values()])
    res["flux_gb_vertex_spread"] = _spread(
        [flux_gb] + [x["flux_gene_backed"] for x in meth.values()])
    res["objective_vertex_rel_diff_max"] = (max(
        abs(x["objective"] - sol.objective_value) for x in meth.values()) / max(
        1.0, abs(sol.objective_value))) if meth else None
    row["timing"]["vertex_" + arm] = round(time.time() - t1, 2)
    if arm == "two_tier":
        cls["ref_internal"] = internal
        row["model"]["slack_int1e-3_abs"] = max(1e-6, args.fva_tol * abs(internal))
    if face_tags and arm in FVA_ARMS and "scan_variable" not in cls:
        t1 = time.time()
        rng0, nfail = feasible_scan(M, [r for r in rx if r.id in gb])
        cls["scan_variable"] = {k for k, (a, b) in rng0.items()
                                if a == a and b == b and b - a > 1e-9}
        row["model"]["n_gene_backed_scan_variable"] = len(cls["scan_variable"])
        row["model"]["n_gene_backed_scan_fixed"] = len(gb) - len(cls["scan_variable"])
        row["model"]["scan_lp_failed"] = nfail
        row["model"]["scan_arm"] = arm
        row["timing"]["scan"] = round(time.time() - t1, 2)
    if face_tags:
        fva = arm in FVA_ARMS
        targets = [r for r in rx if r.id in cls["scan_variable"]] if fva else []
        absent_rx = [r for r in rx if r.id in absn]
        # the cheapest truth-absent coefficient: at slack s, futile fwd+rev flux on that
        # reaction can add up to s / c_min of 'absent flux' to a face maximum
        cmin = min((coef.get(r.forward_variable) for r in absent_rx
                    if coef.get(r.forward_variable) is not None), default=None)
        allsl = {"opt": OPT_SLACK, "opt7": OPT7_SLACK,
                 "int1e-3": max(1e-6, args.fva_tol * abs(cls.get("ref_internal") or internal)),
                 "obj1e-3": args.fva_tol * abs(sol.objective_value)}
        faces = {}
        for tag in face_tags:
            sl = allsl[tag]
            iv = {"absent": absent_rx}
            if arm == "two_tier" and tag == "opt":
                iv["gene_backed"] = [r for r in rx if r.id in gb]     # D's face interval
            fb = face_block(M, sol.objective_value, sl, targets, iv, fva=fva)
            rng = fb.pop("_range")
            if fva:
                fb["n_var_plt1_1e-3"] = sum(1 for k, x in rng.items() if k in unc and x > 1e-3)
                fb["n_var_absent_1e-3"] = sum(1 for k, x in rng.items()
                                              if k in absn and x > 1e-3)
            fb["futile_offset_bound"] = (sl / cmin) if cmin else None
            faces[tag] = fb
        res["faces"] = faces
        o = faces.get("opt") or {}
        amin, amax = o.get("absent_flux_min"), o.get("absent_flux_max")
        good = bool(D) and amin is not None and amax is not None and amin == amin and \
            amax == amax
        res["taf_face_min"] = (amin / D) if good else None
        res["taf_face_max"] = (amax / D) if good else None
        res["taf_face_mid"] = ((amin + amax) / 2.0 / D) if good else None
        res["taf_face_width"] = ((amax - amin) / D) if good else None
        if arm == "two_tier" and D and o.get("gene_backed_flux_max") is not None:
            row["model"]["D_face_rel_width"] = (
                o["gene_backed_flux_max"] - o["gene_backed_flux_min"]) / D
        if fva and "opt7" in faces and "opt" in faces:
            a, b = faces["opt"]["n_var_1e-3"], faces["opt7"]["n_var_1e-3"]
            rel = abs(a - b) / max(1, a)
            faces["opt7"]["n_var_1e-3_rel_diff_vs_opt"] = rel
            res["opt7_count_flag"] = bool(rel > CONFIRM_REL)
    return res


def run_amplicon(job):
    t_all = time.time()
    m = PF.pa()
    from modelseedpy.core.msbuilder import MSBuilder
    key = job["key"]
    row = {k: job[k] for k in ("key", "src_genome_id", "region", "md5", "fold",
                               "group_id", "template", "genus", "n_reps", "stratum",
                               "reps", "novelty_bin", "rank1_identity", "tier", "checks")}
    row.update({"status": "ok", "reason": "", "arms": {}, "timing": {}})
    try:
        tmpl, role_ecs = m.load_template(job["template"])
        core_t = m.load_template("core")[0]
        union = set().union(*map(set, job["member_ec"].values()))
        # ---- 1. gapfilled model of the union (phenotype harness, as_shipped pipeline)
        t0 = time.time()
        mdl, block, _ess, _gc, _gg, _ctx, infeas = PF.gapfill_block(
            m, union, job["template"], PF.safe_id(key, "wobj"), _W["gf_media"],
            _W["gf_cpds"], "as_shipped", None, False)
        row["timing"]["gapfill"] = round(time.time() - t0, 2)
        row["gapfill"] = {k: block[k] for k in ("n_new_rxn", "growth", "sanity_gate_ok",
                                                "yield_ratio")}
        row["gapfill"]["n_reversed_applied"] = len(block["reversed_applied"])
        row["gapfill"]["infeasible"] = infeas
        if infeas:
            row["status"], row["reason"] = "unscoreable", infeas
        elif not block["sanity_gate_ok"]:
            row["status"], row["reason"] = "gated", "sanity_gate_failed"
        elif block["growth"] <= PF.GROWTH_THRESH:
            row["status"], row["reason"] = "unscoreable", "no_growth_on_complete"
        if row["status"] != "ok":
            return finish(row, t_all)
        gapf = {r.id for r in mdl.reactions if r.id[:-1] in set(block["new_rxn"])}
        # ---- 2-3. role and reaction probabilities (production rule)
        genes = sorted(g.id for g in mdl.genes)
        p_role, unknown = role_probabilities(genes, role_ecs, job["member_ec"])
        rx_ids = [r.id for r in mdl.reactions if r.id[:3] == "rxn"]
        gene_backed = {r.id for r in mdl.reactions if r.id[:3] == "rxn" and r.genes}
        prod = {}
        for r in mdl.reactions:
            if r.id[:3] != "rxn":
                continue
            prod[r.id] = (max((p_role.get(g.id, 0.0) for g in r.genes), default=0.0)
                          if r.genes else 0.0)
        union_only = sorted(g for g, p in p_role.items() if p == 0.0)
        # ---- 4. truth support
        t0 = time.time()
        tg = m.genome_from_ec(sorted(job["truth_ec"]), role_ecs, PF.safe_id(key, "truth"))
        tb = MSBuilder(tg, tmpl, ontology_term="RAST", template_core=core_t)
        tbase = tb.build_base_model(PF.safe_id(key, "truth"), "0",
                                    allow_all_non_grp_reactions=False,
                                    annotate_with_rast=False)
        truth_gpr = {r.id for r in tbase.reactions if r.gene_reaction_rule}
        row["timing"]["truth"] = round(time.time() - t0, 2)
        supported = {r for r in rx_ids if r in truth_gpr}
        absent = gene_backed - supported
        unann = {r for r in rx_ids if r not in gene_backed and r not in supported}
        cls = {"gene_backed": gene_backed, "supported": supported, "absent": absent,
               "unann": unann, "gapfilled": gapf}
        pv = Counter(round(prod[r], 6) for r in gene_backed)
        unc = {r for r in gene_backed if prod[r] < 1.0}
        row["model"] = {
            "n_union_ec": len(union), "n_truth_ec": len(job["truth_ec"]),
            "n_rxn": len(mdl.reactions), "n_rxn_prefixed": len(rx_ids),
            "n_gene_backed": len(gene_backed), "n_gapfilled": len(gapf),
            "n_geneless_draft": len(set(rx_ids) - gene_backed - gapf),
            "n_genes": len(genes), "genes_not_template_roles": unknown,
            "n_union_only_roles": len(union_only),
            "n_gene_backed_plt1": len(unc),
            "p_values_gene_backed": {str(k): v for k, v in sorted(pv.items())},
            "n_truth_gpr": len(truth_gpr), "n_supported": len(supported),
            "n_absent": len(absent),
            "static_precision_gb": (len(gene_backed & supported) / len(gene_backed))
            if gene_backed else None,
            "static_precision_plt1": (len(unc & supported) / len(unc)) if unc else None,
            "static_precision_p1": (len((gene_backed - unc) & supported) /
                                    len(gene_backed - unc)) if (gene_backed - unc) else None,
            "n_gapfilled_supported": len(gapf & supported),
        }
        # ---- 5-6. arms
        t0 = time.time()
        v_ref = {}
        for arm in sorted(job["arms"], key=ARMS.index):      # two_tier first
            prob_exp, ex_weight = ARM_PARAMS[arm]
            probs = arm_probs(arm, rx_ids, prod, gene_backed, supported, job["seed"], key)
            M = mdl.copy()
            set_probs(M, probs)
            M.objective = "bio1"            # see prfba_released_copy: the floor's objective
            sol, _u, coef, maxbio = prfba_released_copy(
                M, _W["env"], min_prob=job["min_prob"], prob_exp=prob_exp,
                ex_weight=ex_weight)
            if sol.status != "optimal":
                row["arms"][arm] = {"status": sol.status}
                purge_caches()
                continue
            if arm == "two_tier":
                cls["ref_internal"] = None      # set from this arm's point metrics below
            res = arm_metrics(M, sol, coef, maxbio, cls, prod, probs, job["fva_args"],
                              FACES.get(arm, ()), row, arm)
            v_ref[arm] = res.pop("_v")
            if arm == "real" and job.get("check_released"):
                M2 = mdl.copy()
                set_probs(M2, probs)
                M2.objective = "bio1"
                with contextlib.redirect_stdout(io.StringIO()):
                    s2 = released_prfba()(M2, environment=_W["env"],
                                          min_prob=job["min_prob"], prob_exp=1,
                                          ex_weight=100)
                rel = abs(s2.objective_value - sol.objective_value) / max(
                    1.0, abs(sol.objective_value))
                row["released_check"] = {"objective_copy": sol.objective_value,
                                         "objective_released": s2.objective_value,
                                         "rel_diff": rel, "pass": bool(rel <= 1e-6),
                                         "released_source_sha256": _RELEASED.get("sha256")}
                del M2
            row["arms"][arm] = res
            del M
            purge_caches()
        # flux-vector differences between arms (L1 over all reactions)
        def l1(a, b):
            ks = set(v_ref[a]) | set(v_ref[b])
            return sum(abs(v_ref[a].get(k, 0.0) - v_ref[b].get(k, 0.0)) for k in ks)
        row["l1"] = {"{}-{}".format(a, b): round(l1(a, b), 9)
                     for a, b in (("real", "two_tier"), ("real", "shuffled"),
                                  ("real", "uniform"), ("oracle", "two_tier"),
                                  ("real_exp2", "two_tier"), ("real_ex1", "two_tier_ex1"))
                     if a in v_ref and b in v_ref}
        ex = [row["arms"][a]["exchange_term"] for a in EX100_ARMS
              if "exchange_term" in (row["arms"].get(a) or {})]
        row["exchange_term_rel_spread_ex100"] = (
            (max(ex) - min(ex)) / max(abs(max(ex)), 1e-12)) if len(ex) > 1 else None
        row["identical_real_two_tier"] = (
            ("real" in v_ref and "two_tier" in v_ref and v_ref["real"] == v_ref["two_tier"]
             and row["arms"]["real"]["objective"] == row["arms"]["two_tier"]["objective"]))
        row["timing"]["arms"] = round(time.time() - t0, 2)
    except Exception as exc:                          # noqa: BLE001
        row["status"], row["reason"] = "error", "build_exception:{}:{}".format(
            type(exc).__name__, str(exc)[:300])
    return finish(row, t_all)


def finish(row, t_all):
    purge_caches()
    row["wall_seconds"] = round(time.time() - t_all, 2)
    try:
        import resource
        row["worker_maxrss_gb"] = round(
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0, 3)
    except Exception:
        pass
    return row


# ============================================================ statistics
def cluster_boot_mean(d, clusters, B, seed):
    """Mean of d with a percentile cluster bootstrap (score_folds' draw:
    default_rng(seed).integers(0, K, size=(B, K)) over the sorted cluster keys)."""
    import numpy as np
    keys = sorted(set(clusters))
    idx = {k: i for i, k in enumerate(keys)}
    K = len(keys)
    S = np.zeros(K)
    N = np.zeros(K)
    for x, c in zip(d, clusters):
        S[idx[c]] += x
        N[idx[c]] += 1
    draw = np.random.default_rng(seed).integers(0, K, size=(B, K))
    boots = S[draw].sum(1) / N[draw].sum(1)
    return float(S.sum() / N.sum()), boots


def wilcoxon(d):
    """Signed-rank, zero_method='pratt' (exact zeros are expected: identical optima).
    Effect size: matched-pairs rank-biserial over the non-zero differences."""
    import numpy as np
    d = np.asarray(d, float)
    nz = d[d != 0]
    out = {"n": int(len(d)), "n_nonzero": int(len(nz)), "p": None, "rank_biserial": None}
    if len(nz) == 0:
        return out
    from scipy.stats import rankdata, wilcoxon as W
    r = rankdata(np.abs(nz))
    out["rank_biserial"] = float((r[nz > 0].sum() - r[nz < 0].sum()) / r.sum())
    try:
        out["p"] = float(W(d, zero_method="pratt").pvalue)
    except ValueError:
        pass
    return out


CONTRASTS = (("real", "two_tier"), ("real", "shuffled"), ("real", "uniform"),
             ("oracle", "two_tier"), ("real_exp2", "two_tier"),
             ("real_ex1", "two_tier_ex1"), ("shuffled", "two_tier"))


def metric_value(row, arm, metric):
    a = row["arms"].get(arm) or {}
    if metric.startswith("faces."):
        _f, tag, name = metric.split(".", 2)
        return ((a.get("faces") or {}).get(tag) or {}).get(name)
    return a.get(metric)


METRICS = (
    # primary, vertex-robust: total truth-absent flux on the optimal face (1e-6 absolute
    # slack, every arm) divided by D, the two-tier arm's point gene-backed flux
    "taf_face_mid", "taf_face_min", "taf_face_max", "taf_face_width",
    # secondary, as returned (the solver's default vertex) and its spread over re-solves
    "taf_gb", "taf_point_D", "taf_gb_vertex_spread", "taf_all", "precision_count",
    "n_active_gb", "n_active_absent", "flux_absent", "plt1_share_gb", "internal_share",
    "exchange_term", "biomass",
    # degeneracy (FVA_ARMS): the count on the optimal face decides H2, the 1e-7 face
    # confirms it; near-optimal faces are descriptive (on any face with a positive slack
    # every reaction that can move moves, so their count is the same for every arm)
    "faces.opt.n_var_1e-3", "faces.opt7.n_var_1e-3", "faces.opt.sum_range",
    "faces.int1e-3.sum_range", "faces.int1e-3.n_var_1e-3", "faces.obj1e-3.sum_range",
    "faces.obj1e-3.n_var_1e-3")


def _clean(xs):
    return [float(x) for x in xs if x is not None and x == x]


def smean(xs):
    xs = _clean(xs)
    return (sum(xs) / len(xs)) if xs else None


def smedian(xs):
    return PF.median(_clean(xs))


def smin(xs):
    xs = _clean(xs)
    return min(xs) if xs else None


def smax(xs):
    xs = _clean(xs)
    return max(xs) if xs else None


def summarize(rows, pop_meta, args):
    import numpy as np
    ok = [r for r in rows if r["status"] == "ok" and r.get("arms")]
    weights = {}
    elig = pop_meta.get("eligible_source_genomes_by_stratum_excl_twin", {})
    tw = sum(elig.get(s, 0) for s in MULTI)
    for s in MULTI:
        weights[s] = (elig.get(s, 0) / tw) if tw else None
    out = {"meta": {"generated": now_iso(), "n_rows": len(rows),
                    "status": {"{}|{}".format(k[0], k[1]): v for k, v in sorted(
                        Counter((r["status"], r["reason"].split(":")[0]) for r in rows).items())},
                    "n_ok_by_stratum": dict(Counter(r["stratum"] for r in ok)),
                    "multi_stratum_weights": weights, "sesoi_taf": SESOI_TAF,
                    "sesoi_degeneracy_rel": SESOI_DEGEN_REL, "bootstrap": args.bootstrap,
                    "seed": args.seed},
           "arms": {}, "contrasts": {}, "checks": {}, "diagnostics": {}}
    strata = {s: [r for r in ok if r["stratum"] == s] for s in STRATA}
    strata["multi"] = [r for r in ok if r["stratum"] in MULTI]
    # ---- per-arm descriptive
    for s, rs in strata.items():
        for arm in ARMS:
            for met in METRICS:
                xs = [metric_value(r, arm, met) for r in rs]
                xs = [x for x in xs if x is not None and x == x]
                if xs:
                    out["arms"].setdefault(s, {}).setdefault(arm, {})[met] = {
                        "mean": float(np.mean(xs)), "median": float(np.median(xs)),
                        "n": len(xs)}
    # ---- paired contrasts
    for a, b in CONTRASTS:
        ck = "{}-{}".format(a, b)
        for met in METRICS:
            per = {}
            for s in MULTI + ("n1",):
                pairs = [(metric_value(r, a, met), metric_value(r, b, met), r)
                         for r in strata[s]]
                pairs = [(x, y, r) for x, y, r in pairs if x is not None and y is not None
                         and x == x and y == y]
                if not pairs:
                    continue
                d = [x - y for x, y, _ in pairs]
                # one seed per stratum (the pooled CI combines strata, so their draws must
                # be independent); within a stratum every contrast and metric shares the
                # same draws, so contrasts are paired
                sseed = args.seed + 1000 * STRATA.index(s)
                est, bg = cluster_boot_mean(d, [r["group_id"] for *_, r in pairs],
                                            args.bootstrap, sseed)
                _e2, bs = cluster_boot_mean(d, [r["src_genome_id"] for *_, r in pairs],
                                            args.bootstrap, sseed)
                w = wilcoxon(d)
                base = float(np.mean([y for _, y, _ in pairs]))
                per[s] = {"n": len(d), "mean_delta": est, "median_delta": float(np.median(d)),
                          "ci_group": [float(np.percentile(bg, 2.5)), float(np.percentile(bg, 97.5))],
                          "ci_src": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                          "mean_comparator": base,
                          "frac_lt0": sum(1 for x in d if x < 0) / len(d),
                          "frac_eq0": sum(1 for x in d if x == 0) / len(d),
                          "frac_gt0": sum(1 for x in d if x > 0) / len(d),
                          "wilcoxon_p_pratt": w["p"], "rank_biserial": w["rank_biserial"],
                          "n_nonzero": w["n_nonzero"], "_boot_group": bg}
            # pooled multi: stratum-weighted by eligible source genomes, stratified draws
            if all(s in per for s in MULTI) and all(weights[s] is not None for s in MULTI):
                est = sum(weights[s] * per[s]["mean_delta"] for s in MULTI)
                boots = sum(weights[s] * per[s]["_boot_group"] for s in MULTI)
                dall = [x for s in MULTI for x in
                        [metric_value(r, a, met) - metric_value(r, b, met) for r in strata[s]
                         if metric_value(r, a, met) is not None and
                         metric_value(r, b, met) is not None]]
                w = wilcoxon(dall)
                per["multi"] = {"n": len(dall), "mean_delta_weighted": est,
                                "mean_comparator_weighted": sum(
                                    weights[s] * per[s]["mean_comparator"] for s in MULTI),
                                "ci_group": [float(np.percentile(boots, 2.5)),
                                             float(np.percentile(boots, 97.5))],
                                "mean_delta_unweighted": float(np.mean(dall)),
                                "wilcoxon_p_pratt_unweighted": w["p"],
                                "rank_biserial_unweighted": w["rank_biserial"],
                                "n_nonzero": w["n_nonzero"]}
            for s in per:
                per[s].pop("_boot_group", None)
            if per:
                out["contrasts"].setdefault(ck, {})[met] = per
    # ---- BH for the per-stratum reporting only (the verdict is one pooled test per
    # hypothesis): taf_face_mid x {real-two_tier, real-shuffled} x {n2, n3, n4+}, 6 tests
    fam = []
    for ck in ("real-two_tier", "real-shuffled"):
        for s in MULTI:
            c = ((out["contrasts"].get(ck) or {}).get("taf_face_mid") or {}).get(s)
            if c:
                fam.append((ck, s, c.get("wilcoxon_p_pratt")))
    qs = PF.bh([p for *_, p in fam])
    for (ck, s, _p), q in zip(fam, qs):
        out["contrasts"][ck]["taf_face_mid"][s]["bh_q_per_stratum_family"] = q
    # ---- checks
    n1 = strata["n1"]
    out["checks"]["n1_real_identical_two_tier"] = {
        "n": len(n1), "identical": sum(1 for r in n1 if r.get("identical_real_two_tier"))}
    rc = [r["released_check"] for r in ok if r.get("released_check")]
    out["checks"]["released_check"] = {"n": len(rc), "pass": sum(1 for x in rc if x["pass"]),
                                       "max_rel_diff": max((x["rel_diff"] for x in rc),
                                                           default=None)}
    out["checks"]["genes_not_template_roles"] = sum(
        len(r["model"]["genes_not_template_roles"]) for r in ok)
    out["checks"]["face_lp_failures"] = sum(
        f.get("n_lp_failed", 0) + f.get("n_interval_lp_failed", 0) for r in ok
        for a in r["arms"].values() for f in (a.get("faces") or {}).values())
    out["checks"]["face_lp_retried"] = sum(
        f.get("n_lp_retried", 0) for r in ok
        for a in r["arms"].values() for f in (a.get("faces") or {}).values())
    # n1: the same LP in both arms, so the face values agree up to LP tolerance (the two arms
    # reach their face LPs from different warm starts; the first amended smoke: <= 2e-11)
    d1 = []
    for r in n1:
        a, b = r["arms"].get("real") or {}, r["arms"].get("two_tier") or {}
        if a.get("taf_face_mid") is not None and b.get("taf_face_mid") is not None:
            d1.append(max(abs(a[k] - b[k]) for k in ("taf_face_min", "taf_face_max")))
    out["checks"]["n1_taf_face_agree_real_two_tier"] = {
        "n": len(n1), "n_compared": len(d1), "tol": 1e-9,
        "agree": sum(1 for x in d1 if x <= 1e-9), "max_abs_diff": max(d1) if d1 else None}
    flags, vtx = {}, {}
    for arm in ARMS:
        fl = [(r["arms"].get(arm) or {}).get("opt7_count_flag") for r in ok]
        fl = [x for x in fl if x is not None]
        if fl:
            flags[arm] = {"n": len(fl), "flagged_gt_5pct": sum(1 for x in fl if x),
                          "keys_flagged": [r["key"] for r in ok
                                           if (r["arms"].get(arm) or {}).get("opt7_count_flag")]}
        aa = [r["arms"].get(arm) or {} for r in ok]
        sp = [a.get("taf_gb_vertex_spread") for a in aa if a.get("taf_gb_vertex_spread")
              is not None]
        if aa:
            vtx[arm] = {"n_rows_resolved": sum(1 for a in aa if a.get("n_vertex_resolves")),
                        "n_rows_skipped": sum(1 for a in aa if "skipped" in
                                              (a.get("vertex_resolves") or {})),
                        "n_taf_gb_spread_gt_1e-6": sum(1 for x in sp if x > 1e-6),
                        "max_taf_gb_spread": max(sp) if sp else None,
                        "max_taf_point_D_spread": smax(a.get("taf_point_D_vertex_spread")
                                                       for a in aa),
                        "max_objective_rel_diff": smax(a.get("objective_vertex_rel_diff_max")
                                                       for a in aa)}
    out["checks"]["opt7_count_flags"] = flags
    out["checks"]["vertex_resolves"] = vtx
    exs = [r.get("exchange_term_rel_spread_ex100") for r in ok
           if r.get("exchange_term_rel_spread_ex100") is not None]
    out["checks"]["exchange_term_rel_spread_ex100"] = {
        "n": len(exs), "max": max(exs) if exs else None,
        "n_le_1e-6": sum(1 for x in exs if x <= 1e-6)}
    dw = [r["model"].get("D_face_rel_width") for r in ok
          if r["model"].get("D_face_rel_width") is not None]
    out["checks"]["D_face_rel_width_max"] = max(dw) if dw else None
    out["checks"]["candidate_checks_failed_in_sample"] = {
        k: sum(1 for r in rows if not r["checks"][k]) for k in
        ("union_eq_pred_ec", "f1_agrees", "n_truth_agrees")}
    # ---- diagnostics: addressable set, flux differences, oracle headroom, composition
    for s, rs in strata.items():
        if not rs:
            continue
        tt = [r["arms"].get("two_tier") or {} for r in rs]
        absn = [x.get("flux_absent") for x in tt if x.get("flux_absent") is not None]
        addr = [x["flux_plt1_absent"] / x["flux_absent"] for x in tt
                if x.get("flux_absent") and x["flux_absent"] > 0]
        dg = {"n": len(rs),
              "n_gene_backed_plt1_median": smedian(r["model"]["n_gene_backed_plt1"] for r in rs),
              "n_gene_backed_median": smedian(r["model"]["n_gene_backed"] for r in rs),
              "addressable_share_of_absent_flux_two_tier_mean": smean(addr),
              "frac_any_flux_difference_real_vs_two_tier": smean(
                  1.0 if r.get("l1", {}).get("real-two_tier", 0) > 1e-6 else 0.0 for r in rs),
              "internal_share_real_median": smedian(
                  (r["arms"].get("real") or {}).get("internal_share") for r in rs),
              "internal_share_real_min": smin(
                  (r["arms"].get("real") or {}).get("internal_share") for r in rs),
              "internal_share_real_max": smax(
                  (r["arms"].get("real") or {}).get("internal_share") for r in rs),
              "static_precision_plt1_mean": smean(r["model"]["static_precision_plt1"] for r in rs),
              "static_precision_p1_mean": smean(r["model"]["static_precision_p1"] for r in rs),
              "flux_absent_two_tier_mean": smean(absn)}
        for met, name in (("taf_face_mid", "oracle_headroom_captured"),
                          ("taf_gb", "oracle_headroom_captured_point")):
            c_or = ((out["contrasts"].get("oracle-two_tier") or {}).get(met) or {}).get(s)
            c_re = ((out["contrasts"].get("real-two_tier") or {}).get(met) or {}).get(s)
            if c_or and c_re:
                o = c_or.get("mean_delta", c_or.get("mean_delta_weighted"))
                rr = c_re.get("mean_delta", c_re.get("mean_delta_weighted"))
                dg[name] = (rr / o) if o not in (None, 0) else None
        out["diagnostics"][s] = dg
    out["verdict"] = verdict(out)
    return out


def verdict(out):
    """The decision rules of WEIGHTED_OBJECTIVE_TEST.md section 10 (as amended in its
    'Review' section), applied mechanically to the pooled multi-genome contrasts."""
    def ci(ck, met, s="multi"):
        c = ((out["contrasts"].get(ck) or {}).get(met) or {}).get(s) or {}
        m = c.get("mean_delta", c.get("mean_delta_weighted"))
        return m, c.get("ci_group")
    v = {}
    chk = out["checks"]
    # ---- C1 and C2 gate every other number
    defects, notes = [], []
    c1 = chk.get("n1_real_identical_two_tier") or {}
    if not c1.get("n"):
        notes.append("C1 not evaluable: no n1 rows")
    elif c1.get("identical", 0) < c1["n"]:
        defects.append("C1: {} of {} n1 rows differ between real and two_tier".format(
            c1["n"] - c1.get("identical", 0), c1["n"]))
    c2 = chk.get("released_check") or {}
    if not c2.get("n"):
        defects.append("C2: released_check was not run")
    elif c2.get("pass", 0) < c2["n"]:
        defects.append("C2: {} of {} released_check rows fail".format(
            c2["n"] - c2.get("pass", 0), c2["n"]))
    v["notes"] = notes
    if defects:
        v["harness_defect"] = defects
        for k in ("positive_control_oracle_expressible", "H1", "H1s_specific_to_p", "H2"):
            v[k] = "not read: harness defect"
        return v
    E = ("taf_face_mid", "taf_face_min", "taf_face_max")

    def three(ck):
        return {e: ci(ck, e) for e in E}

    def ub_below(d, e):
        return bool(d[e][1] and d[e][1][1] < 0)

    def lb_above(d, e):
        return bool(d[e][1] and d[e][1][0] > 0)

    def inside(d, e):
        return bool(d[e][1] and -SESOI_TAF < d[e][1][0] and d[e][1][1] < SESOI_TAF)

    orc, rt, rs = three("oracle-two_tier"), three("real-two_tier"), three("real-shuffled")
    v["numbers"] = {ck: {e: {"mean": d[e][0], "ci_group": d[e][1]} for e in E}
                    for ck, d in (("oracle-two_tier", orc), ("real-two_tier", rt),
                                  ("real-shuffled", rs))}
    mp, cp = ci("real-two_tier", "taf_gb")
    v["numbers"]["real-two_tier"]["taf_gb_point_as_returned"] = {"mean": mp, "ci_group": cp}
    # ---- C3: can ANY reweighting of gene-backed reactions move truth-absent flux?
    m_o = orc["taf_face_mid"][0]
    expressible = bool(ub_below(orc, "taf_face_mid") and m_o is not None and
                       m_o <= -SESOI_TAF and ub_below(orc, "taf_face_min") and
                       ub_below(orc, "taf_face_max"))
    v["positive_control_oracle_expressible"] = expressible
    # ---- H1 (decided on the face midpoint; the sign must hold at both face endpoints)
    m, c = rt["taf_face_mid"]
    ends_below = ub_below(rt, "taf_face_min") and ub_below(rt, "taf_face_max")
    ends_above = lb_above(rt, "taf_face_min") and lb_above(rt, "taf_face_max")
    ends_inside = inside(rt, "taf_face_min") and inside(rt, "taf_face_max")
    spec = ub_below(rs, "taf_face_mid")
    v["H1s_specific_to_p"] = spec
    structural = ("" if expressible else "; the oracle also fails, so the objective cannot "
                  "express a reweighting of gene-backed reactions (structural null)")
    if c is None:
        h1 = "not evaluable"
    elif c[1] < 0 and m <= -SESOI_TAF and ends_below and spec:
        h1 = ("works: real lowers truth-absent flux on the optimal face vs two-tier by >= "
              "SESOI, at both face endpoints, and beats shuffled")
    elif c[0] > 0 and m >= SESOI_TAF and ends_above:
        h1 = ("harmful: real raises truth-absent flux on the optimal face vs two-tier by >= "
              "SESOI, at both face endpoints")
    elif -SESOI_TAF < c[0] and c[1] < SESOI_TAF and ends_inside:
        h1 = (("does nothing (equivalence): the 95% CIs of the face midpoint and both "
               "endpoints lie inside +-SESOI" if c[0] <= 0 <= c[1] else
               "equivalent within SESOI, nonzero: the CIs lie inside +-SESOI but the "
               "midpoint CI excludes 0") + structural)
    elif c[1] < 0 and m <= -SESOI_TAF:
        why = []
        if not ends_below:
            why.append("not at both face endpoints (vertex-dependent)")
        if not spec:
            why.append("not specific to p (real - shuffled midpoint CI does not lie below 0)")
        h1 = "lowers truth-absent flux by >= SESOI at the face midpoint, but " + \
            "; ".join(why)
    elif c[0] > 0 and m >= SESOI_TAF:
        h1 = ("raises truth-absent flux by >= SESOI at the face midpoint, but not at both "
              "face endpoints (vertex-dependent)")
    elif c[1] < 0 or c[0] > 0:
        h1 = ("detectable but below SESOI: the midpoint CI excludes 0, |mean| < SESOI, and "
              "the CIs do not lie inside +-SESOI")
    else:
        h1 = "inconclusive: the midpoint CI crosses 0 and a SESOI bound" + structural
    v["H1"] = h1

    # ---- H2 (count on the optimal face; the 1e-7 face confirms)
    def h2(met):
        md, cd = ci("real-two_tier", met)
        cm = ((out["contrasts"].get("real-two_tier") or {}).get(met) or {}).get("multi") or {}
        base = cm.get("mean_comparator_weighted")
        _mds, cds = ci("real-shuffled", met)
        sp = ("; specific to p (real also below shuffled)" if (cds and cds[1] < 0) else
              "; NOT specific to p: shuffled costs do as well or better" if cds else "")
        if cd is None or not base:
            return "not evaluable"
        if cd[1] < 0 and -md / base >= SESOI_DEGEN_REL:
            return ("reduces degeneracy: >= 10% fewer variable gene-backed reactions on the "
                    "optimal face" + sp)
        if cd[1] < 0:
            return "reduces degeneracy, below the 10% SESOI" + sp
        if cd[0] > 0:
            return "increases degeneracy" + sp
        return "no detectable change in degeneracy" + sp
    v["H2"] = h2("faces.opt.n_var_1e-3")
    v["H2_at_1e-7_face"] = h2("faces.opt7.n_var_1e-3")
    v["H2_confirmed_at_1e-7_face"] = v["H2"] == v["H2_at_1e-7_face"]
    return v


# ============================================================ driver
def library_provenance():
    libs = PF.library_provenance()
    try:
        import mscommunity
        p = os.path.dirname(os.path.abspath(mscommunity.__file__))
        try:
            dirty = bool(subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
                                        cwd=p, check=True, capture_output=True,
                                        text=True).stdout.strip())
        except Exception:
            dirty = None
        # imported by the released prFBA (and so by released_check) but never called on
        # the single-model path; recorded because it is another session's moving checkout
        libs["mscommunity"] = {"path": p, "git_commit": git_commit(p), "git_dirty": dirty}
    except Exception:
        libs["mscommunity"] = None
    return libs


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", default=os.path.join(D, "cells.csv.gz"))
    ap.add_argument("--fold-plan", default=os.path.join(D, "fold_plan.json"))
    ap.add_argument("--panel", default=os.path.join(D, "panel_eligible.json"))
    ap.add_argument("--panel-ec", default=os.path.join(D, "_cache_panel_ec.json"))
    ap.add_argument("--ec-universe", default=os.path.join(D, "panel_ec_universe.json"))
    ap.add_argument("--prfba", default=os.path.join(D, "prfba"))
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--select-mode", default="exact_tie")
    ap.add_argument("--density", type=float, default=1.0)
    ap.add_argument("--n-per-stratum", default="n2=100,n3=100,n4+=100,n1=60",
                    help="sample size per stratum (one amplicon per source genome)")
    ap.add_argument("--seed", type=int, default=1729)
    ap.add_argument("--min-truth-ec", type=int, default=200)
    ap.add_argument("--region-pref", default="V4")
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--min-prob", type=float, default=0.01)
    ap.add_argument("--fva-tol", type=float, default=1e-3,
                    help="slack of the descriptive near-optimal faces: a fraction of the "
                         "two-tier arm's internal term (int1e-3) and of the real arm's whole "
                         "objective (obj1e-3).  The optimal face uses OPT_SLACK (1e-6)")
    ap.add_argument("--no-check-released", action="store_true")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None,
                    help="first N jobs of the stratum-interleaved population (smoke)")
    ap.add_argument("--only", default=None, metavar="KEY,KEY")
    ap.add_argument("--templates", default=TEMPLATES_DEFAULT)
    ap.add_argument("--no-integrate-reversals", action="store_true")
    ap.add_argument("--msd", default=PF.MSD_DEFAULT)
    ap.add_argument("--expect-modelseedpy", default=None, metavar="COMMIT")
    ap.add_argument("--maxtasksperchild", type=int, default=10)
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--summarize-only", action="store_true")
    ap.add_argument("--population-only", action="store_true")
    args = ap.parse_args(argv)
    args.n_per_stratum = {k: int(v) for k, v in
                          (kv.split("=") for kv in args.n_per_stratum.split(",") if kv)}
    args.arms = [a for a in args.arms.split(",") if a]
    bad = [a for a in args.arms if a not in ARMS]
    if bad or "real" not in args.arms or "two_tier" not in args.arms:
        raise SystemExit("[wobj] --arms must include real and two_tier; unknown: {}".format(bad))
    return args


def main(argv=None):
    args = parse_args(argv)
    libs = library_provenance()
    msp = libs.get("modelseedpy") or {}
    print("[wobj] modelseedpy {} at {} (tracked files changed: {})".format(
        msp.get("path"), msp.get("git_commit"), msp.get("git_dirty")), flush=True)
    if args.expect_modelseedpy and not (
            (msp.get("git_commit") or "").startswith(args.expect_modelseedpy)
            and msp.get("git_dirty") is False):
        raise SystemExit("[wobj] modelseedpy is not a clean checkout of {}; put the pinned "
                         "clone first on PYTHONPATH".format(args.expect_modelseedpy))
    os.makedirs(args.outdir, exist_ok=True)
    pop = build_or_load_population(args)
    if args.population_only:
        return 0
    jsonl = os.path.join(args.outdir, "per_amplicon.jsonl")
    jobs = pop["jobs"]
    if args.only:
        keep = set(args.only.split(","))
        jobs = [j for j in jobs if j["key"] in keep]
    if args.limit:
        jobs = jobs[:args.limit]
    done = {}
    if os.path.exists(jsonl):
        with open(jsonl) as fh:
            for line in fh:
                if line.strip():
                    r = json.loads(line)
                    if r["status"] != "error":
                        done[r["key"]] = r
    fva_args = argparse.Namespace(fva_tol=args.fva_tol)
    todo = []
    for j in jobs:
        if j["key"] in done:
            continue
        jj = dict(j)
        jj.update({"arms": args.arms, "seed": args.seed, "min_prob": args.min_prob,
                   "fva_args": fva_args, "check_released": not args.no_check_released})
        todo.append(jj)
    print("[wobj] {} jobs selected, {} done, {} to run".format(len(jobs), len(jobs) - len(todo),
                                                              len(todo)), flush=True)
    if todo and not args.summarize_only:
        import multiprocessing as mp
        ctx = mp.get_context("spawn")             # Gurobi environments do not survive fork
        nw = max(1, min(args.workers, len(todo)))
        winit = (args.msd, args.templates, not args.no_integrate_reversals)
        with open(jsonl, "a") as fh:
            if nw == 1:
                worker_init(*winit)
                results = map(run_amplicon, todo)
                pool = None
            else:
                pool = ctx.Pool(nw, initializer=worker_init, initargs=winit,
                                maxtasksperchild=args.maxtasksperchild)
                results = pool.imap_unordered(run_amplicon, todo)
            try:
                for r in results:
                    fh.write(json.dumps(r, sort_keys=True) + "\n")
                    fh.flush()
                    done[r["key"]] = r
                    print("[wobj] {} {} {} {}s".format(r["key"], r["stratum"], r["status"],
                                                       r["wall_seconds"]), flush=True)
            finally:
                if pool is not None:
                    pool.close()
                    pool.join()
    rows = [done[j["key"]] for j in jobs if j["key"] in done]
    released_prfba()                 # records the released source's sha256 in this process
    shas = sorted({r["released_check"]["released_source_sha256"] for r in rows
                   if r.get("released_check")})
    summary = summarize(rows, pop["meta"], args)
    summary["checks"]["released_source_sha256_in_rows"] = shas
    summary["meta"].update({
        "libraries": libs, "templates": PF.template_provenance(PF.use_templates(args.templates)),
        "integrate_reversals": not args.no_integrate_reversals,
        "released_prfba_source_sha256": _RELEASED.get("sha256"),
        "script": {"path": os.path.abspath(__file__), "sha256": sha256(os.path.abspath(__file__))},
        "phenotype_fold_sha256": sha256(os.path.join(HERE, "phenotype_fold.py")),
        "phenotype_arm_sha256": sha256(os.path.join(PF.P2_SCRIPTS, "phenotype_arm.py")),
        "git_commit": git_commit(REPO), "population": pop["meta"],
        "args": {k: (v if not isinstance(v, argparse.Namespace) else vars(v))
                 for k, v in vars(args).items()},
        "per_amplicon_jsonl": os.path.abspath(jsonl)})
    import resource
    summary["meta"]["parent_maxrss_gb"] = round(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0, 3)
    # the largest single reaped child (a pool worker), not their sum
    summary["meta"]["largest_child_maxrss_gb"] = round(
        resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1048576.0, 3)
    with open(os.path.join(args.outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True, default=str)
    print("[wobj] wrote {} rows -> {}".format(len(rows), os.path.join(args.outdir,
                                                                     "summary.json")), flush=True)
    print("[wobj] verdict: {}".format(json.dumps(summary["verdict"])), flush=True)
    print("[wobj] peak RSS: parent {} GB, largest pool worker {} GB, largest row {} GB".format(
        summary["meta"]["parent_maxrss_gb"], summary["meta"]["largest_child_maxrss_gb"],
        max((r.get("worker_maxrss_gb") or 0 for r in rows), default=None)), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
