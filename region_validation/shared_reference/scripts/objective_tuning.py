#!/usr/bin/env python
"""Stage O of PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md: tuning prFBA's probability-weighted
FBA objective (declared 2026-11-02, post hoc; the design is binding).

SAMPLE AND MODEL (identical to the weighted-objective test, WEIGHTED_OBJECTIVE_TEST.md)
  The 360 amplicons of data/weighted_objective/population.json.  Per amplicon, one worker:
  the gapfilled model of the selection's UNION EC set is rebuilt exactly as
  weighted_objective_test.run_amplicon builds it (phenotype_fold.gapfill_block, as_shipped
  pipeline, complete medium, gapfiller reversals applied; ModelSEEDpy pinned at c2e31a0 by
  PYTHONPATH and --expect-modelseedpy; templates ModelSEEDTemplates-c12run 21a8084).  The
  rebuilt model is checked against the prior run's row (data/weighted_objective/
  per_amplicon.jsonl): gapfill.n_new_rxn, gapfill.growth (1e-6), model.n_gene_backed,
  model.n_gapfilled (and n_rxn, n_supported, n_absent, n_gapfilled_supported); any mismatch
  is recorded and flagged (row.prior_match), never silently accepted.  Truth support and the
  gene-backed / gapfilled / gene-less classification are run_amplicon's.

PROBABILITY SOURCES (role level; reaction p = max over the reaction's genes; gene-less 0)
  kn            the production k/n (weighted_objective_test.role_probabilities): the prior
                run's `real` arm.
  est           Stage E's estimator from the weights file (data/pfunc_estimators/
                stage_o_weights.json): for a role whose EC set is S,
                  q = (sum_j w_j 1[S <= EC(g_j)] + alpha pi(S)) / (sum_j w_j + alpha),
                pi(S) = share of the amplicon's clade members whose EC set contains S
                (the alpha*pi term is 0 when the clade is null; alpha stays in the
                denominator: the literal formula, counted in row.est_info); EC(g) from the
                fold's data/prfba/f{fold:02d}_d100/_cache_ec_subset.json.  Never a product
                of per-EC values.  meta.isotonic, when not null, is applied to q (linear
                interpolation x -> y, clipped at the ends, as sklearn's IsotonicRegression
                .predict; isotonic.interp = "step" selects a right-continuous step map).
  est_shuffled  est permuted among the model's gene-backed reactions by the prior run's
                own helper (weighted_objective_test.arm_probs 'shuffled', rng_for(seed, key,
                'shuffle')): one permutation per amplicon, the prior shuffled arm's, reused
                in all 18 (C, G, X) cells.
  two_tier      gene-backed p = 1 (what production gives at n = 1).
  oracle        gene-backed p = 1 if truth-supported else 0 (the headroom).

COST SHAPES AND THE OTHER FACTORS (on rxn-prefixed reactions; exchange reactions X)
  C  lin = max(0.01, 1 - p) (released, prob_exp 1);  log = max(0.01, ln max(p, 0.01) / ln 0.01).
     Both give cost 1 at p = 0 and 0.01 at p = 1, so gene-less draft reactions cost 1 in
     every arm.
  G  gapfilled reactions take p = G in EVERY arm (P, two_tier and oracle alike): G0 = 0
     (released), Gmeas = the pooled truth-support rate of gapfilled reactions over the TUNE
     half of the prior run (sum n_gapfilled_supported / sum n_gapfilled; computed in the
     parent, recorded in meta), G0.2 = 0.2.  Gapfilled = gapfill_block's new_rxn (gene-less;
     asserted disjoint from the gene-backed set on every row).
  X  exchange weight 100 (released), 10, 1.
  Arms: 5 P x 2 C x 3 G x 3 X = 90, named "P/C/G/X" (shipped = kn/lin/G0/X100), each solved
  on a FRESH copy of the gapfilled model through prfba_tuned_copy() (below), which builds
  exactly prfba_released_copy's constraint set (medium at uptake 100, the ElementUptakePkg
  carbon cap, bio1 >= 0.95 x max) and takes the coefficient map from cost_map().  A fresh
  copy costs about 0.15 s; a cold re-solve on a reused copy was not bitwise identical to a
  fresh solve in the pre-build probe, so no reuse.

CHECKS ON EVERY ROW
  * two_tier/lin/G0/X100 runs first: D = its gene-backed flux at its returned optimum (as
    before), fixed for every arm of the amplicon; compared with the prior row's D.
  * kn/lin/G0/X100 is solved through prfba_tuned_copy like every arm; its coefficient map is
    then compared by variable name with weighted_objective_test.prfba_released_copy's on a
    second fresh copy (shipped_check), and the released source itself is re-solved on a
    third copy (released_check, as before; it must equal the prior row's
    released_check.objective_released).
  * Its objective, returned vertex and face values are compared with the prior `real` arm;
    two_tier/lin/G0/X100 with the prior `two_tier` arm (row.prior_match).
  * Determinism pairs (identical LPs, so bitwise-identical results expected): two_tier and
    oracle at G0 under lin vs log, every X; kn vs two_tier on n1 rows in all 18 cells.

MEASURES PER ARM
  status, objective, biomass and max biomass; the truth-absent gene-backed flux on the
  optimal face (objective <= opt + 1e-6, weighted_objective_test.face_block without FVA):
  taf_face_min / max / mid = F / D; the gapfilled-flux interval on the same face; as
  returned: taf_gb, taf_point_D; guard metrics: total exchange flux sum |v_ex| (from the
  primal, not exchange_term / X), gene-backed flux, flux precision (truth-supported share of
  gene-backed flux), gapfilled flux, internal / exchange terms and the internal-term split.
  --fva-arms ARM,...: also FVA of the scan-variable gene-backed reactions on the optimal face
  and the 1e-7 face (as before), AFTER the primary face block so the primary numbers do not
  depend on the flag (a second pass, finalists only, on the confirm half).
  --repro-check: also re-runs the prior run's `two_tier` and `real` arms through
  weighted_objective_test's own code (prfba_released_copy + arm_metrics: vertex re-solves,
  feasible scan, opt-face FVA), which should reproduce the prior row bitwise (row.repro).

SPLIT
  A holdout group is in the tune half iff int(sha256(("1729|" + group_id).encode())
  .hexdigest(), 16) % 2 == 0, else confirm; each row records its half.

OUTPUT
  <outdir>/per_amplicon.jsonl (one row per amplicon, all arms; resumable: finished keys are
  skipped), <outdir>/meta.json.  Selection and confirmation: scripts/objective_tuning_
  results.py.

Interpreter: /home/freiburger/Documents/py_venv/bin/python (the only interpreter).
PYTHONPATH must put the pinned ModelSEEDpy clone first; --expect-modelseedpy enforces it.
"""
import argparse
import bisect
import contextlib
import gc
import hashlib
import io
import json
import math
import os
import sys
import time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)
D = os.path.join(SR, "data")
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import weighted_objective_test as W   # noqa: E402  (no cobra at import time)

PF = W.PF
PRIOR_DIR = os.path.join(D, "weighted_objective")
WEIGHTS_REAL = os.path.join(D, "pfunc_estimators", "stage_o_weights.json")
OUT_REAL = os.path.join(D, "objective_tuning")

P_SOURCES = ("kn", "est", "est_shuffled", "two_tier", "oracle")
P_CANDIDATES = ("kn", "est")                  # selectable; the other three are controls
C_SHAPES = ("lin", "log")
G_LEVELS = ("G0", "Gmeas", "G0.2")
X_LEVELS = ("X100", "X10", "X1")
X_VALUE = {"X100": 100, "X10": 10, "X1": 1}   # ints, as the released default ex_weight=100
MIN_PROB = 0.01
SHIPPED = "kn/lin/G0/X100"
D_ARM = "two_tier/lin/G0/X100"
SPLIT_SALT = "1729|"


def arm_name(P, C, G, X):
    return "/".join((P, C, G, X))


def arm_factors(name):
    P, C, G, X = name.split("/")
    return {"P": P, "C": C, "G": G, "X": X}


def all_arms():
    """D_ARM first (it sets D), SHIPPED second, then the rest in factor order."""
    rest = [arm_name(P, C, G, X) for P in P_SOURCES for C in C_SHAPES for G in G_LEVELS
            for X in X_LEVELS]
    return [D_ARM, SHIPPED] + [a for a in rest if a not in (D_ARM, SHIPPED)]


def half_of(group_id):
    return "tune" if int(hashlib.sha256((SPLIT_SALT + str(group_id)).encode()).hexdigest(),
                         16) % 2 == 0 else "confirm"


def now_iso():
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def sha256(path):
    return W.sha256(path)


def maxrss_gb():
    import resource
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0, 3)


# ============================================================ coefficients (NEW; the
# released function and weighted_objective_test.prfba_released_copy are untouched)
def cost_lin(p, min_prob=MIN_PROB):
    """The released coefficient at prob_exp 1: max(min_prob, (1 - float(p) ** 1))."""
    return max(min_prob, (1 - float(p) ** 1))


def cost_log(p, min_prob=MIN_PROB):
    """max(min_prob, ln max(p, min_prob) / ln min_prob): 1 at p <= min_prob, min_prob at p = 1."""
    return max(min_prob, math.log(max(float(p), min_prob)) / math.log(min_prob))


COST = {"lin": cost_lin, "log": cost_log}


def cost_map(probs, shape, min_prob=MIN_PROB):
    """{rxn_id: objective coefficient} for the rxn-prefixed reactions of one arm."""
    f = COST[shape]
    return {r: f(p, min_prob) for r, p in probs.items()}


def prfba_tuned_copy(model, environment, rxn_cost, ex_weight):
    """weighted_objective_test.prfba_released_copy (itself MSProbability.prFBA, ModelSEEDpy
    c2e31a0 msprobability.py l.155-201, with the edits needed for one model) with ONE
    difference: the coefficient of a rxn-prefixed reaction is rxn_cost[rxn.id] (from
    cost_map) instead of max(min_prob, 1 - float(rxn.notes['probability']) ** prob_exp).
    Every bound and constraint, their order, and the exchange coefficient (ex_weight on both
    directions) are prfba_released_copy's.  The caller sets bio1 as the objective first.
    Returns (solution, MSModelUtil, coef, maxBioSol); the model is modified in place."""
    from modelseedpy.core.msmodelutl import MSModelUtil
    from modelseedpy.fbapkg.elementuptakepkg import ElementUptakePkg
    from optlang.symbolics import Zero

    mdlUtil = MSModelUtil(model)
    mdlUtil.add_medium(environment)
    elepkg = ElementUptakePkg(mdlUtil.model)  ;  elepkg.build_package({"C": 100})
    maxBioSol = mdlUtil.model.slim_optimize()
    mdlUtil.add_minimal_objective_cons(maxBioSol*.95)
    coef = {}
    for rxn in mdlUtil.model.reactions:
        if "rxn" == rxn.id[0:3]:
            c = rxn_cost[rxn.id]
            coef.update({rxn.forward_variable: c})
            coef.update({rxn.reverse_variable: c})
        elif "EX_" == rxn.id[0:3]:
            coef.update({rxn.forward_variable: ex_weight})
            coef.update({rxn.reverse_variable: ex_weight})
    mdlUtil.add_objective(Zero, "min", coef)
    return mdlUtil.model.optimize(), mdlUtil, coef, maxBioSol


# ============================================================ the estimator (role level)
def ec_mask(ecs, index):
    """Bitmask of an EC set over the universe index; None if an EC is outside it."""
    m = 0
    for e in ecs:
        i = index.get(e)
        if i is None:
            return None
        m |= 1 << i
    return m


def apply_isotonic(q, iso):
    if not iso:
        return q
    xs, ys = iso["x"], iso["y"]
    if iso.get("interp", "linear") == "step":
        i = bisect.bisect_right(xs, q) - 1
        return float(ys[max(0, min(i, len(ys) - 1))])
    if q <= xs[0]:
        return float(ys[0])
    if q >= xs[-1]:
        return float(ys[-1])
    i = bisect.bisect_right(xs, q)
    x0, x1, y0, y1 = xs[i - 1], xs[i], ys[i - 1], ys[i]
    return float(y0 if x1 == x0 else y0 + (y1 - y0) * (q - x0) / (x1 - x0))


def est_role_probabilities(gene_ids, role_ecs, est, index, iso):
    """q(role) for the model's genes (= template roles), from one amplicon's weights record
    (neighbour and clade-member EC sets as universe bitmasks).  Genes that are not template
    roles get no entry (reaction p then takes 0 from them, as kn)."""
    neigh = [(est["masks"][g], float(w)) for g, w in est["neighbours"]]
    wsum = sum(w for _, w in neigh)
    alpha = float(est["alpha"])
    members = [est["masks"][g] for g in (est["clade"] or {}).get("members", [])]
    den = wsum + alpha
    q, info = {}, Counter()
    for g in gene_ids:
        rec = role_ecs.get(g)
        if rec is None:
            info["genes_not_template_roles"] += 1
            continue
        m = ec_mask(rec[1], index)
        if m is None:
            info["roles_with_ec_outside_universe"] += 1
            q[g] = apply_isotonic(0.0, iso)
            continue
        num = sum(w for e, w in neigh if e & m == m)
        if members and alpha:
            num += alpha * (sum(1 for e in members if e & m == m) / len(members))
        q[g] = apply_isotonic(num / den, iso)
    return q, dict(info)


# ============================================================ per-arm measurement
def point_metrics(M, sol, coef, maxbio, cls, probs, D):
    """Everything read at one arm's returned optimum (the solver's vertex)."""
    v = W.mag(M)
    gb, absn, gapf = cls["gene_backed"], cls["absent"], cls["gapfilled"]
    rx = cls["rx_ids"]
    ex = [r.id for r in M.reactions if r.id[:3] == "EX_"]
    tot = lambda ids: sum(v[i] for i in ids)
    flux_gb, flux_abs = tot(gb), tot(absn)
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
    act = {i for i in rx if v[i] > W.ACTIVE}
    bio = sol.fluxes.get("bio1")
    return {
        "status": sol.status, "objective": sol.objective_value,
        "max_biomass": maxbio, "biomass": bio,
        "biomass_ratio": (bio / maxbio) if maxbio else None,
        "flux_exchange": tot(ex), "flux_gene_backed": flux_gb, "flux_absent": flux_abs,
        "flux_gapfilled": tot(gapf), "flux_rxn": tot(rx),
        "precision_flux": ((flux_gb - flux_abs) / flux_gb) if flux_gb > 0 else None,
        "taf_gb": (flux_abs / flux_gb) if flux_gb > 0 else None,
        "taf_point_D": (flux_abs / D) if D else None,
        "internal_term": internal, "exchange_term": exch,
        "internal_share": (internal / (internal + exch)) if (internal + exch) > 0 else None,
        "internal_term_split": dict(split),
        "n_active_gb": len(act & gb), "n_active_absent": len(act & absn),
        "n_active_gapfilled": len(act & gapf),
    }


def copy_sets(M, cls):
    """The reaction sets as reactions OF THIS COPY, in model order (a reaction object of
    another copy carries another solver's variables)."""
    rx = [r for r in M.reactions if r.id[:3] == "rxn"]
    return {"rx": rx, "absent": [r for r in rx if r.id in cls["absent"]],
            "gapfilled": [r for r in rx if r.id in cls["gapfilled"]],
            "gene_backed": [r for r in rx if r.id in cls["gene_backed"]]}


def primary_face(M, sol, coef, sets, D, with_gene_backed):
    """The optimal face (objective <= opt + OPT_SLACK, one absolute slack for every arm):
    min / max of the total truth-absent gene-backed flux (the primary) and of the total
    gapfilled flux; for the D arm also of the gene-backed flux (D's face robustness)."""
    absent_rx = sets["absent"]
    iv = {"absent": absent_rx, "gapfilled": sets["gapfilled"]}
    if with_gene_backed:
        iv["gene_backed"] = sets["gene_backed"]
    fb = W.face_block(M, sol.objective_value, W.OPT_SLACK, [], iv, fva=False)
    fb.pop("_range", None)
    amin, amax = fb.get("absent_flux_min"), fb.get("absent_flux_max")
    good = bool(D) and amin is not None and amax is not None and amin == amin and amax == amax
    cmin = min((coef.get(r.forward_variable) for r in absent_rx
                if coef.get(r.forward_variable) is not None), default=None)
    out = {"absent_flux_min": amin, "absent_flux_max": amax,
           "taf_face_min": (amin / D) if good else None,
           "taf_face_max": (amax / D) if good else None,
           "taf_face_mid": ((amin + amax) / 2.0 / D) if good else None,
           "taf_face_width": ((amax - amin) / D) if good else None,
           "gapfilled_flux_face_min": fb.get("gapfilled_flux_min"),
           "gapfilled_flux_face_max": fb.get("gapfilled_flux_max"),
           "face_lp_failed": fb["n_interval_lp_failed"], "face_lp_retried": fb["n_lp_retried"],
           "futile_offset_bound": (W.OPT_SLACK / cmin) if cmin else None}
    if with_gene_backed:
        out["gene_backed_flux_face_min"] = fb.get("gene_backed_flux_min")
        out["gene_backed_flux_face_max"] = fb.get("gene_backed_flux_max")
    return out


def fva_faces(M, sol, cls, sets, probs, row):
    """FVA degeneracy as in the weighted-objective test (section 8): the feasible scan once
    per amplicon, then FVA of the scan-variable gene-backed reactions on the optimal face
    (1e-6) and the 1e-7 confirmation face."""
    rx = sets["rx"]
    gb = cls["gene_backed"]
    if "scan_variable" not in cls:
        t1 = time.time()
        rng0, nfail = W.feasible_scan(M, [r for r in rx if r.id in gb])
        cls["scan_variable"] = {k for k, (a, b) in rng0.items()
                                if a == a and b == b and b - a > 1e-9}
        row["model"]["n_gene_backed_scan_variable"] = len(cls["scan_variable"])
        row["model"]["scan_lp_failed"] = nfail
        row["timing"]["scan"] = round(time.time() - t1, 2)
    targets = [r for r in rx if r.id in cls["scan_variable"]]
    unc_kn = {r for r in gb if cls["kn"][r] < 1.0}
    unc_arm = {r for r in gb if probs[r] < 1.0}
    faces = {}
    for tag, sl in (("opt", W.OPT_SLACK), ("opt7", W.OPT7_SLACK)):
        fb = W.face_block(M, sol.objective_value, sl, targets, {"absent": sets["absent"]},
                          fva=True)
        rng = fb.pop("_range")
        fb["n_var_plt1_kn_1e-3"] = sum(1 for k, x in rng.items() if k in unc_kn and x > 1e-3)
        fb["n_var_plt1_arm_1e-3"] = sum(1 for k, x in rng.items() if k in unc_arm and x > 1e-3)
        fb["n_var_absent_1e-3"] = sum(1 for k, x in rng.items()
                                      if k in cls["absent"] and x > 1e-3)
        faces[tag] = fb
    a, b = faces["opt"]["n_var_1e-3"], faces["opt7"]["n_var_1e-3"]
    rel = abs(a - b) / max(1, a)
    faces["opt7"]["n_var_1e-3_rel_diff_vs_opt"] = rel
    return faces, bool(rel > W.CONFIRM_REL)


# ============================================================ worker
_O = {}


def release(*models):
    """Free each finished copy's Gurobi model.  Every cobra Model.copy() leaves its optlang
    solver objects alive after the copy is collected (an upstream optlang/gurobipy reference
    the harness cannot break: ~6-7 MiB per copy of a 1,500-reaction model, ~10,700 objects);
    disposing the gurobipy model returns the native part (~1.5 MiB).  Called only after every
    measurement of the copy, so no number depends on it.  The rest is why the launcher runs
    --maxtasksperchild 1."""
    for M in models:
        try:
            M.solver.problem.dispose()
        except Exception:                             # noqa: BLE001
            pass


def worker_init(msd, template_dir, integrate_reversals, universe):
    W.worker_init(msd, template_dir, integrate_reversals)
    _O["index"] = {e: i for i, e in enumerate(universe)}


def _same(a, b, tol=0.0):
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= tol


def prior_match(row, prior):
    """The rebuilt model against the prior run's row: the declared keys must match."""
    if not prior:
        return {"prior_row": False, "ok": False, "mismatches": ["no prior row"]}
    pg, pm = prior.get("gapfill") or {}, prior.get("model") or {}
    g, m = row["gapfill"], row["model"]
    mm = []
    if g["n_new_rxn"] != pg.get("n_new_rxn"):
        mm.append("gapfill.n_new_rxn {} != {}".format(g["n_new_rxn"], pg.get("n_new_rxn")))
    if pg.get("growth") is None or abs(g["growth"] - pg["growth"]) > 1e-6:
        mm.append("gapfill.growth {} vs {}".format(g["growth"], pg.get("growth")))
    for k in ("n_gene_backed", "n_gapfilled"):
        if m.get(k) != pm.get(k):
            mm.append("model.{} {} != {}".format(k, m.get(k), pm.get(k)))
    extra = [k for k in ("n_rxn", "n_supported", "n_absent", "n_gapfilled_supported",
                         "n_genes", "n_gene_backed_plt1") if m.get(k) != pm.get(k)]
    return {"prior_row": True, "ok": not mm, "mismatches": mm,
            "extra_mismatches": {k: [m.get(k), pm.get(k)] for k in extra}}


PRIOR_ARM_KEYS = ("objective", "flux_gene_backed", "flux_absent", "taf_gb", "taf_face_min",
                  "taf_face_max", "taf_face_mid")


def compare_prior_arm(mine, prior_arm):
    out = {}
    for k in PRIOR_ARM_KEYS:
        a, b = mine.get(k), (prior_arm or {}).get(k)
        out[k] = {"bitwise": (a == b) if (a is not None and b is not None) else None,
                  "abs_diff": abs(a - b) if (a is not None and b is not None) else None}
    return out


def run_amplicon(job):
    t_all = time.time()
    m = PF.pa()
    from modelseedpy.core.msbuilder import MSBuilder
    key = job["key"]
    row = {k: job[k] for k in ("key", "src_genome_id", "region", "md5", "fold", "group_id",
                               "template", "genus", "n_reps", "stratum", "reps",
                               "novelty_bin", "rank1_identity", "tier")}
    row.update({"half": half_of(job["group_id"]), "status": "ok", "reason": "", "arms": {},
                "timing": {"arm_seconds": {}}})
    try:
        tmpl, role_ecs = m.load_template(job["template"])
        core_t = m.load_template("core")[0]
        union = set().union(*map(set, job["member_ec"].values()))
        # ---- 1. gapfilled model of the union (weighted_objective_test.run_amplicon, verbatim)
        t0 = time.time()
        mdl, block, _ess, _gc, _gg, _ctx, infeas = PF.gapfill_block(
            m, union, job["template"], PF.safe_id(key, "wobj"), W._W["gf_media"],
            W._W["gf_cpds"], "as_shipped", None, False)
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
        row["maxrss_gb_after_gapfill"] = maxrss_gb()
        gapf = {r.id for r in mdl.reactions if r.id[:-1] in set(block["new_rxn"])}
        # ---- 2. probabilities: kn (production rule) and est (Stage E), role level
        genes = sorted(g.id for g in mdl.genes)
        p_role, unknown = W.role_probabilities(genes, role_ecs, job["member_ec"])
        rx_ids = [r.id for r in mdl.reactions if r.id[:3] == "rxn"]
        gene_backed = {r.id for r in mdl.reactions if r.id[:3] == "rxn" and r.genes}
        rxn_genes = {r.id: [g.id for g in r.genes] for r in mdl.reactions if r.id[:3] == "rxn"}
        prod = {r: (max((p_role.get(g, 0.0) for g in gs), default=0.0) if gs else 0.0)
                for r, gs in rxn_genes.items()}
        union_only = sorted(g for g, p in p_role.items() if p == 0.0)
        q_role, est_counts = est_role_probabilities(genes, role_ecs, job["est"], _O["index"],
                                                    job["isotonic"])
        est = {r: (max((q_role.get(g, 0.0) for g in gs), default=0.0) if gs else 0.0)
               for r, gs in rxn_genes.items()}
        # ---- 3. truth support (run_amplicon, verbatim)
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
        unc = {r for r in gene_backed if prod[r] < 1.0}
        pv = Counter(round(prod[r], 6) for r in gene_backed)
        row["model"] = {
            "n_union_ec": len(union), "n_truth_ec": len(job["truth_ec"]),
            "n_rxn": len(mdl.reactions), "n_rxn_prefixed": len(rx_ids),
            "n_gene_backed": len(gene_backed), "n_gapfilled": len(gapf),
            "n_gapfilled_and_gene_backed": len(gapf & gene_backed),
            "n_geneless_draft": len(set(rx_ids) - gene_backed - gapf),
            "n_genes": len(genes), "genes_not_template_roles": unknown,
            "n_union_only_roles": len(union_only),
            "n_gene_backed_plt1": len(unc),
            "p_values_gene_backed": {str(k): v for k, v in sorted(pv.items())},
            "n_truth_gpr": len(truth_gpr), "n_supported": len(supported),
            "n_absent": len(absent), "n_gapfilled_supported": len(gapf & supported),
            "static_precision_gb": (len(gene_backed & supported) / len(gene_backed))
            if gene_backed else None,
        }
        row["prior_match"] = prior_match(row, job["prior"])
        if gapf & gene_backed:
            # the G factor is defined on gene-less gapfilled reactions; the design never
            # met an overlap (0 on the probe), so stop rather than choose silently
            raise RuntimeError("{} gapfilled reactions are gene-backed".format(
                len(gapf & gene_backed)))
        cls = {"gene_backed": gene_backed, "supported": supported, "absent": absent,
               "unann": unann, "gapfilled": gapf, "rx_ids": rx_ids, "kn": prod}
        row["est_info"] = est_descriptives(job, est, prod, gene_backed, supported, est_counts,
                                           q_role)
        # ---- 4. per-arm probabilities (gapfilled reactions take G in every arm)
        base = {"kn": prod, "est": est,
                "est_shuffled": W.arm_probs("shuffled", rx_ids, est, gene_backed, supported,
                                            job["seed"], key),
                "two_tier": W.arm_probs("two_tier", rx_ids, prod, gene_backed, supported,
                                        job["seed"], key),
                "oracle": W.arm_probs("oracle", rx_ids, prod, gene_backed, supported,
                                      job["seed"], key)}
        gval = {"G0": 0.0, "Gmeas": float(job["g_meas"]), "G0.2": 0.2}
        # ---- 5. arms.  The worker heap holds ~2.5 million long-lived objects (biochemistry,
        # templates, the gapfilled model): a full collection costs ~1.2 s, so they are
        # frozen once (gc.freeze) and each arm's copy is collected right after it (~0.02 s)
        # instead of piling up in the oldest generation.  finish() unfreezes.
        gc.collect()
        gc.freeze()
        t0 = time.time()
        Dv = None
        for arm in job["arms"]:
            ta = time.time()
            f = arm_factors(arm)
            probs = dict(base[f["P"]])
            for r in gapf:
                probs[r] = gval[f["G"]]
            M = mdl.copy()
            M.objective = "bio1"            # see prfba_released_copy: the floor's objective
            sol, _u, coef, maxbio = prfba_tuned_copy(M, W._W["env"], cost_map(probs, f["C"]),
                                                     X_VALUE[f["X"]])
            res = {"factors": f, "g_value": gval[f["G"]], "x_value": X_VALUE[f["X"]]}
            if sol.status != "optimal":
                res["status"] = sol.status
                row["arms"][arm] = res
                release(M)
                del M, sol, _u, coef
                W.purge_caches()
                continue
            if arm == D_ARM:
                v0 = W.mag(M)
                Dv = sum(v0[i] for i in gene_backed)
                row["model"]["D_flux_gb_two_tier"] = Dv
            sets = copy_sets(M, cls)
            res.update(point_metrics(M, sol, coef, maxbio, cls, probs, Dv))
            res.update(primary_face(M, sol, coef, sets, Dv, with_gene_backed=(arm == D_ARM)))
            if arm == D_ARM and Dv and res.get("gene_backed_flux_face_max") is not None:
                row["model"]["D_face_rel_width"] = (res["gene_backed_flux_face_max"] -
                                                    res["gene_backed_flux_face_min"]) / Dv
            if arm in job["fva_arms"]:
                t1 = time.time()
                res["faces"], res["opt7_count_flag"] = fva_faces(M, sol, cls, sets, probs, row)
                row["timing"].setdefault("fva", {})[arm] = round(time.time() - t1, 2)
            if arm == SHIPPED:
                shipped_checks(row, job, mdl, probs, prod, coef, sol)
            row["arms"][arm] = res
            release(M)
            del M, sol, _u, coef, sets
            W.purge_caches()
            gc.collect()
            row["timing"]["arm_seconds"][arm] = round(time.time() - ta, 3)
        row["timing"]["arms"] = round(time.time() - t0, 2)
        row["maxrss_gb_after_arms"] = maxrss_gb()
        # ---- 6. against the prior run, and the determinism pairs
        pa_ = (job["prior"] or {}).get("arms") or {}
        pm = row["prior_match"]
        if SHIPPED in row["arms"]:
            pm["shipped_vs_prior_real"] = compare_prior_arm(row["arms"][SHIPPED],
                                                            pa_.get("real"))
        if D_ARM in row["arms"]:
            pm["two_tier_vs_prior_two_tier"] = compare_prior_arm(row["arms"][D_ARM],
                                                                 pa_.get("two_tier"))
        pD = ((job["prior"] or {}).get("model") or {}).get("D_flux_gb_two_tier")
        pm["D"] = {"mine": Dv, "prior": pD, "bitwise": (Dv == pD) if pD is not None else None,
                   "abs_diff": abs(Dv - pD) if (Dv is not None and pD is not None) else None}
        row["determinism"] = determinism(row)
        # ---- 7. optional: the prior arms through the prior harness's own code
        if job.get("repro_check"):
            t1 = time.time()
            row["repro"] = repro_check(job, mdl, cls, prod, rx_ids, gene_backed, supported)
            row["timing"]["repro"] = round(time.time() - t1, 2)
    except Exception as exc:                          # noqa: BLE001
        import traceback
        row["status"], row["reason"] = "error", "build_exception:{}:{}".format(
            type(exc).__name__, str(exc)[:300])
        row["traceback"] = traceback.format_exc()[-2000:]
    return finish(row, t_all)


def shipped_checks(row, job, mdl, probs, prod, coef, sol):
    """kn/lin/G0/X100: (a) its probabilities are the production rule's (gapfilled 0);
    (b) its coefficient map equals prfba_released_copy's, by variable name, on a second
    fresh copy (that copy is the prior run's `real` arm verbatim); (c) the released source
    itself re-solves the same objective on a third copy (released_check)."""
    env = W._W["env"]
    out = {"probs_equal_production_rule": probs == prod}
    M2 = mdl.copy()
    W.set_probs(M2, probs)
    M2.objective = "bio1"
    s2, _u2, coef2, _mb2 = W.prfba_released_copy(M2, env, min_prob=MIN_PROB, prob_exp=1,
                                                 ex_weight=100)
    a = {var.name: c for var, c in coef.items()}
    b = {var.name: c for var, c in coef2.items()}
    out.update({"coef_map_identical": a == b, "n_coef": len(a),
                "n_coef_differ": sum(1 for k in set(a) | set(b) if a.get(k) != b.get(k)),
                "objective_released_copy": s2.objective_value,
                "objective_bitwise_equal": s2.objective_value == sol.objective_value})
    row["shipped_check"] = out
    release(M2)
    del M2, s2, _u2, coef2
    W.purge_caches()
    M3 = mdl.copy()
    W.set_probs(M3, probs)
    M3.objective = "bio1"
    with contextlib.redirect_stdout(io.StringIO()):
        s3 = W.released_prfba()(M3, environment=env, min_prob=MIN_PROB, prob_exp=1,
                                ex_weight=100)
    rel = abs(s3.objective_value - sol.objective_value) / max(1.0, abs(sol.objective_value))
    prior_rel = ((job["prior"] or {}).get("released_check") or {}).get("objective_released")
    row["released_check"] = {
        "objective_copy": sol.objective_value, "objective_released": s3.objective_value,
        "rel_diff": rel, "pass": bool(rel <= 1e-6),
        "released_source_sha256": W._RELEASED.get("sha256"),
        "prior_objective_released": prior_rel,
        "equals_prior_released_bitwise": (s3.objective_value == prior_rel)
        if prior_rel is not None else None,
        "abs_diff_prior_released": abs(s3.objective_value - prior_rel)
        if prior_rel is not None else None}
    release(M3)
    del M3, s3
    W.purge_caches()


DET_KEYS = ("objective", "flux_gene_backed", "flux_absent", "flux_exchange", "taf_face_min",
            "taf_face_max", "gapfilled_flux_face_min", "gapfilled_flux_face_max")


def determinism(row):
    """Pairs of arms that are the same LP (identical coefficient maps): two_tier and oracle
    at G0 under lin vs log (both shapes give 0.01 at p = 1 and 1 at p = 0), every X; and kn
    vs two_tier in all 18 cells on n1 rows (k/n = 1 for every gene-backed reaction)."""
    pairs = []
    for X in X_LEVELS:
        for P in ("two_tier", "oracle"):
            pairs.append((arm_name(P, "lin", "G0", X), arm_name(P, "log", "G0", X)))
    if row["stratum"] == "n1":
        for C in C_SHAPES:
            for G in G_LEVELS:
                for X in X_LEVELS:
                    pairs.append((arm_name("kn", C, G, X), arm_name("two_tier", C, G, X)))
    n = ident = 0
    worst, diffs = 0.0, []
    for a, b in pairs:
        A, B = row["arms"].get(a), row["arms"].get(b)
        if not A or not B or A.get("status") != "optimal" or B.get("status") != "optimal":
            continue
        n += 1
        same = all(A.get(k) == B.get(k) for k in DET_KEYS)
        ident += int(same)
        if not same:
            d = max((abs(A[k] - B[k]) for k in DET_KEYS
                     if A.get(k) is not None and B.get(k) is not None), default=float("nan"))
            worst = max(worst, d) if d == d else worst
            diffs.append([a, b, d])
    return {"n_pairs": n, "n_identical": ident, "max_abs_diff": worst, "differing": diffs}


REPRO_KEYS = PRIOR_ARM_KEYS + ("taf_point_D", "precision_flux", "internal_term",
                               "exchange_term")


def repro_check(job, mdl, cls, prod, rx_ids, gene_backed, supported):
    """The prior run's `two_tier` and `real` arms re-run through weighted_objective_test's
    own code path (arm_probs, set_probs, prfba_released_copy, arm_metrics with the optimal
    face only: vertex re-solves, the feasible scan on two_tier, FVA + truth-absent interval
    on the 1e-6 face).  The LP sequence is the prior run's, so its numbers should equal the
    prior row's bitwise."""
    cls_p = {k: cls[k] for k in ("gene_backed", "supported", "absent", "unann", "gapfilled")}
    row_p = {"model": {}, "timing": {}}
    fargs = argparse.Namespace(fva_tol=1e-3)
    prior_arms = (job["prior"] or {}).get("arms") or {}
    out = {"arms": {}}
    for arm in ("two_tier", "real"):
        probs = W.arm_probs(arm, rx_ids, prod, gene_backed, supported, job["seed"], job["key"])
        M = mdl.copy()
        W.set_probs(M, probs)
        M.objective = "bio1"
        sol, _u, coef, maxbio = W.prfba_released_copy(M, W._W["env"], min_prob=MIN_PROB,
                                                      prob_exp=1, ex_weight=100)
        if sol.status != "optimal":
            out["arms"][arm] = {"status": sol.status}
            release(M)
            continue
        if arm == "two_tier":
            cls_p["ref_internal"] = None
        res = W.arm_metrics(M, sol, coef, maxbio, cls_p, prod, probs, fargs, ("opt",),
                            row_p, arm)
        res.pop("_v", None)
        mine = {k: res.get(k) for k in REPRO_KEYS}
        mine["opt_n_var_1e-3"] = ((res.get("faces") or {}).get("opt") or {}).get("n_var_1e-3")
        mine["opt_sum_range"] = ((res.get("faces") or {}).get("opt") or {}).get("sum_range")
        pr = prior_arms.get(arm) or {}
        pr_flat = dict(pr)
        pr_flat["opt_n_var_1e-3"] = ((pr.get("faces") or {}).get("opt") or {}).get("n_var_1e-3")
        pr_flat["opt_sum_range"] = ((pr.get("faces") or {}).get("opt") or {}).get("sum_range")
        cmp_ = {}
        for k, a in mine.items():
            b = pr_flat.get(k)
            cmp_[k] = {"mine": a, "prior": b,
                       "bitwise": (a == b) if (a is not None and b is not None) else None}
        out["arms"][arm] = {"values": cmp_,
                            "all_bitwise": all(c["bitwise"] for c in cmp_.values()
                                               if c["bitwise"] is not None)}
        release(M)
        del M, sol, _u, coef
        W.purge_caches()
        gc.collect()
    pm = (job["prior"] or {}).get("model") or {}
    out["D"] = {"mine": row_p["model"].get("D_flux_gb_two_tier"),
                "prior": pm.get("D_flux_gb_two_tier")}
    out["n_gene_backed_scan_variable"] = {"mine": row_p["model"].get(
        "n_gene_backed_scan_variable"), "prior": pm.get("n_gene_backed_scan_variable")}
    out["all_bitwise"] = all(a.get("all_bitwise") for a in out["arms"].values()) and \
        out["D"]["mine"] == out["D"]["prior"]
    return out


def est_descriptives(job, est, prod, gene_backed, supported, counts, q_role):
    """What est looks like on this model, beside kn: histograms, static separation of
    truth-supported from truth-absent gene-backed reactions (AUROC, reaction-level Brier)."""
    gb = sorted(gene_backed)
    y = [1 if r in supported else 0 for r in gb]

    def auroc(p):
        pos = [p[r] for r, t in zip(gb, y) if t]
        neg = [p[r] for r, t in zip(gb, y) if not t]
        if not pos or not neg:
            return None
        allv = sorted((x, i) for i, x in enumerate(pos + neg))
        ranks = [0.0] * len(allv)
        i = 0
        while i < len(allv):
            j = i
            while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
                j += 1
            for k in range(i, j + 1):
                ranks[allv[k][1]] = (i + j) / 2.0 + 1
            i = j + 1
        rp = sum(ranks[:len(pos)])
        return (rp - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))

    def brier(p):
        return (sum((p[r] - t) ** 2 for r, t in zip(gb, y)) / len(gb)) if gb else None

    hist = Counter(round(est[r], 2) for r in gb)
    e = job["est"]
    return {"n_neighbours": len(e["neighbours"]),
            "weight_sum": sum(float(w) for _, w in e["neighbours"]),
            "alpha": e["alpha"], "clade_rank": (e["clade"] or {}).get("rank"),
            "n_clade_members": len((e["clade"] or {}).get("members", [])),
            "alpha_with_null_clade": bool(e["clade"] is None and float(e["alpha"]) > 0),
            "isotonic": bool(job["isotonic"]), "counts": counts,
            "n_roles_q": len(q_role),
            "n_gene_backed_est_lt1": sum(1 for r in gb if est[r] < 1.0),
            "n_gene_backed_est_lt0.5": sum(1 for r in gb if est[r] < 0.5),
            "n_gene_backed_est_0": sum(1 for r in gb if est[r] == 0.0),
            "mean_est_gene_backed": (sum(est[r] for r in gb) / len(gb)) if gb else None,
            "mean_kn_gene_backed": (sum(prod[r] for r in gb) / len(gb)) if gb else None,
            "est_hist_gene_backed": {str(k): v for k, v in sorted(hist.items())},
            "auroc_supported_vs_absent": {"kn": auroc(prod), "est": auroc(est)},
            "brier_vs_truth_support": {"kn": brier(prod), "est": brier(est)}}


def finish(row, t_all):
    gc.unfreeze()
    W.purge_caches()
    row["wall_seconds"] = round(time.time() - t_all, 2)
    try:
        row["worker_maxrss_gb"] = maxrss_gb()
    except Exception:                                 # noqa: BLE001
        pass
    return row


# ============================================================ parent: inputs
def load_prior(path):
    rows = {}
    with open(path) as fh:
        for line in fh:
            if line.strip():
                r = json.loads(line)
                rows[r["key"]] = r
    return rows


def slim_prior(r):
    """The parts of a prior row the worker compares against (keeps the job small)."""
    if r is None:
        return None
    arms = {}
    for a in ("real", "two_tier"):
        x = (r.get("arms") or {}).get(a) or {}
        arms[a] = {k: x.get(k) for k in REPRO_KEYS}
        arms[a]["faces"] = {"opt": {k: ((x.get("faces") or {}).get("opt") or {}).get(k)
                                    for k in ("n_var_1e-3", "sum_range")}}
    return {"gapfill": r.get("gapfill"), "model": r.get("model"), "arms": arms,
            "released_check": r.get("released_check"), "status": r.get("status")}


def g_meas_from_prior(prior):
    rows = [r for r in prior.values() if r.get("status") == "ok" and
            half_of(r["group_id"]) == "tune"]
    num = sum(r["model"]["n_gapfilled_supported"] for r in rows)
    den = sum(r["model"]["n_gapfilled"] for r in rows)
    return num / den, {"n_rows": len(rows), "sum_n_gapfilled_supported": num,
                       "sum_n_gapfilled": den, "half": "tune",
                       "rule": "sum n_gapfilled_supported / sum n_gapfilled over the prior "
                               "run's tune-half rows (status ok)"}


def load_weights(path, pop_jobs, allow_stub):
    """Read and validate the Stage E weights file against the population.  Returns
    (weights, report).  Raises SystemExit on any violation."""
    with open(path) as fh:
        Wt = json.load(fh)
    meta = Wt.get("meta") or {}
    if meta.get("stub") and not allow_stub:
        raise SystemExit("[otune] {} is a STUB (meta.stub); refusing without --allow-stub"
                         .format(path))
    amps = Wt.get("amplicons") or {}
    bad = []
    for j in pop_jobs:
        a = amps.get(j["key"])
        if a is None:
            bad.append("{}: missing".format(j["key"]))
            continue
        if a.get("md5") != j["md5"] or int(a.get("fold", -1)) != int(j["fold"]):
            bad.append("{}: md5/fold differ".format(j["key"]))
        if a.get("reps") is not None and list(a["reps"]) != list(j["reps"]):
            bad.append("{}: reps differ".format(j["key"]))
        if a.get("n_reps") is not None and int(a["n_reps"]) != int(j["n_reps"]):
            bad.append("{}: n_reps differ".format(j["key"]))
        nb = a.get("neighbours") or []
        ws = [float(w) for _, w in nb]
        if any(not (w == w) or w < 0 or math.isinf(w) for w in ws):
            bad.append("{}: a weight is negative or not finite".format(j["key"]))
        alpha = float(a.get("alpha", 0.0))
        if alpha < 0 or sum(ws) + alpha <= 0:
            bad.append("{}: sum of weights + alpha <= 0".format(j["key"]))
        cl = a.get("clade")
        gl = [g for g, _ in nb] + (list(cl.get("members") or []) if cl else [])
        if j["src_genome_id"] in gl:
            bad.append("{}: the source genome is a neighbour or clade member".format(j["key"]))
    iso = meta.get("isotonic")
    if iso:
        xs, ys = iso.get("x"), iso.get("y")
        if not xs or not ys or len(xs) != len(ys) or any(b < a for a, b in zip(xs, xs[1:])) \
                or any(b < a for a, b in zip(ys, ys[1:])):
            bad.append("meta.isotonic is not a non-decreasing map with equal-length x and y")
        if iso.get("interp", "linear") not in ("linear", "step"):
            bad.append("meta.isotonic.interp must be linear or step")
    if bad:
        raise SystemExit("[otune] weights file {} fails validation ({} problems): {}".format(
            path, len(bad), "; ".join(bad[:20])))
    # counted, not refusing: the exact-tie reps (the genomes production selects) should be
    # neighbours at w = 1 (the top-identity tie set of the stage-0 ranking)
    cnt = Counter()
    for j in pop_jobs:
        a = amps[j["key"]]
        w_of = {g: float(w) for g, w in a.get("neighbours") or []}
        cnt["reps_not_among_neighbours"] += sum(1 for g in j["reps"] if g not in w_of)
        cnt["reps_with_weight_ne_1"] += sum(1 for g in j["reps"] if g in w_of and w_of[g] != 1.0)
        cnt["rows_with_a_rep_missing_or_ne_1"] += int(any(w_of.get(g) != 1.0 for g in j["reps"]))
        if a.get("clade") is None:
            cnt["null_clade_rows"] += 1
            cnt["null_clade_rows_without_neighbours"] += int(not w_of)
        cnt["rows_without_neighbours"] += int(not w_of)
    rep = {"counts": dict(cnt),
           "path": os.path.abspath(path), "sha256": sha256(path), "stub": bool(meta.get("stub")),
           "meta": {k: v for k, v in meta.items() if k not in ("assertions",)},
           "n_amplicons_in_file": len(amps),
           "n_alpha_with_null_clade": sum(1 for j in pop_jobs if (amps[j["key"]].get("clade")
                                                                   is None and
                                                                   float(amps[j["key"]].get(
                                                                       "alpha", 0)) > 0))}
    return Wt, rep


def attach_est(jobs, Wt, universe):
    """Per job: the neighbour and clade-member EC sets as universe bitmasks, from that
    fold's _cache_ec_subset.json (loaded one fold at a time)."""
    import orjson
    index = {e: i for i, e in enumerate(universe)}
    amps = Wt["amplicons"]
    by_fold = {}
    for j in jobs:
        by_fold.setdefault(int(j["fold"]), []).append(j)
    A = Counter()
    for fold, js in sorted(by_fold.items()):
        with open(os.path.join(D, "prfba", "f{:02d}_d100".format(fold),
                               "_cache_ec_subset.json"), "rb") as fh:
            ecsub = orjson.loads(fh.read())
        for j in js:
            a = amps[j["key"]]
            gl = [g for g, _ in a["neighbours"]] + list((a.get("clade") or {}).get("members")
                                                       or [])
            masks = {}
            for g in gl:
                if g in masks:
                    continue
                ecs = ecsub.get(g)
                if ecs is None:
                    raise SystemExit("[otune] {}: genome {} is not in fold {}'s reference "
                                     "(_cache_ec_subset.json)".format(j["key"], g, fold))
                mk = 0
                for e in ecs:
                    i = index.get(e)
                    if i is None:
                        A["ec_outside_universe_dropped"] += 1
                        continue
                    mk |= 1 << i
                masks[g] = mk
            j["est"] = {"neighbours": [[g, float(w)] for g, w in a["neighbours"]],
                        "clade": a.get("clade"), "alpha": float(a.get("alpha", 0.0)),
                        "masks": masks}
            A["genomes_masked"] += len(masks)
        del ecsub
    return dict(A)


# ============================================================ driver
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--population", default=os.path.join(PRIOR_DIR, "population.json"))
    ap.add_argument("--prior", default=os.path.join(PRIOR_DIR, "per_amplicon.jsonl"))
    ap.add_argument("--weights", default=WEIGHTS_REAL)
    ap.add_argument("--allow-stub", action="store_true",
                    help="accept a weights file with meta.stub (smoke tests only)")
    ap.add_argument("--ec-universe", default=os.path.join(D, "panel_ec_universe.json"))
    ap.add_argument("--half", choices=("all", "tune", "confirm"), default="all")
    ap.add_argument("--arms", default="all",
                    help="comma list of P/C/G/X arm names (default all 90); %s is always "
                         "added (it sets D)" % D_ARM)
    ap.add_argument("--fva-arms", default="",
                    help="comma list of arms that also get FVA degeneracy (opt, opt7 faces)")
    ap.add_argument("--repro-check", action="store_true",
                    help="also re-run the prior two_tier and real arms through the prior "
                         "harness's own code (exact reproduction test; smoke)")
    ap.add_argument("--seed", type=int, default=1729)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--maxtasksperchild", type=int, default=1,
                    help="1: a fresh worker per amplicon (each of the ~92 model copies of a "
                         "task leaks ~6 MiB upstream; see release())")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--only", default=None, metavar="KEY,KEY")
    ap.add_argument("--templates", default=W.TEMPLATES_DEFAULT)
    ap.add_argument("--no-integrate-reversals", action="store_true")
    ap.add_argument("--msd", default=PF.MSD_DEFAULT)
    ap.add_argument("--expect-modelseedpy", default=None, metavar="COMMIT")
    ap.add_argument("--validate-only", action="store_true",
                    help="parent checks only (population, prior, weights file, fold EC caches, "
                         "g_meas); writes nothing, starts no worker (the launcher's pre-flight)")
    args = ap.parse_args(argv)
    every = all_arms()
    if args.arms == "all":
        args.arms = every
    else:
        want = {a for a in args.arms.split(",") if a} | {D_ARM}
        bad = sorted(want - set(every))
        if bad:
            raise SystemExit("[otune] unknown arms: {}".format(bad))
        args.arms = [a for a in every if a in want]
    args.fva_arms = [a for a in args.fva_arms.split(",") if a]
    bad = sorted(set(args.fva_arms) - set(args.arms))
    if bad:
        raise SystemExit("[otune] --fva-arms not among --arms: {}".format(bad))
    return args


def main(argv=None):
    args = parse_args(argv)
    libs = W.library_provenance()
    msp = libs.get("modelseedpy") or {}
    print("[otune] modelseedpy {} at {} (tracked files changed: {})".format(
        msp.get("path"), msp.get("git_commit"), msp.get("git_dirty")), flush=True)
    if args.expect_modelseedpy and not (
            (msp.get("git_commit") or "").startswith(args.expect_modelseedpy)
            and msp.get("git_dirty") is False):
        raise SystemExit("[otune] modelseedpy is not a clean checkout of {}; put the pinned "
                         "clone first on PYTHONPATH".format(args.expect_modelseedpy))
    out_abs = os.path.abspath(args.outdir)
    w_abs = os.path.abspath(args.weights)
    if out_abs == os.path.abspath(OUT_REAL) and (args.allow_stub or ".STUB." in w_abs):
        raise SystemExit("[otune] the real output directory never takes stub weights")
    if not args.validate_only:
        os.makedirs(args.outdir, exist_ok=True)
    with open(args.population) as fh:
        pop = json.load(fh)
    prior = load_prior(args.prior)
    g_meas, g_meas_info = g_meas_from_prior(prior)
    jobs = list(pop["jobs"])
    if args.only:
        keep = set(args.only.split(","))
        jobs = [j for j in jobs if j["key"] in keep]
    if args.half != "all":
        jobs = [j for j in jobs if half_of(j["group_id"]) == args.half]
    if args.limit:
        jobs = jobs[:args.limit]
    Wt, wrep = load_weights(args.weights, jobs, args.allow_stub)
    universe = sorted(PF.norm_ec(PF.jload(args.ec_universe)["ec"]))
    est_attach = attach_est(jobs, Wt, universe)
    iso = (Wt.get("meta") or {}).get("isotonic")
    if args.validate_only:
        print(json.dumps({"validate_only": "ok", "n_jobs": len(jobs),
                          "n_jobs_by_half": dict(Counter(half_of(j["group_id"]) for j in jobs)),
                          "g_meas": g_meas, "g_meas_info": g_meas_info,
                          "weights": {k: wrep[k] for k in ("path", "sha256", "stub",
                                                           "n_amplicons_in_file",
                                                           "n_alpha_with_null_clade",
                                                           "counts")},
                          "isotonic": bool(iso), "est_attach": est_attach,
                          "arms": len(args.arms), "fva_arms": args.fva_arms}, indent=1),
              flush=True)
        return 0
    jsonl = os.path.join(args.outdir, "per_amplicon.jsonl")
    meta_path = os.path.join(args.outdir, "meta.json")
    script_sha = sha256(os.path.abspath(__file__))
    spec = {"arms": args.arms, "fva_arms": args.fva_arms, "seed": args.seed,
            "min_prob": MIN_PROB, "g_meas": g_meas, "weights_sha256": wrep["sha256"],
            "population_sha256": sha256(args.population), "prior_sha256": sha256(args.prior),
            "templates": os.path.abspath(args.templates),
            "integrate_reversals": not args.no_integrate_reversals,
            "repro_check": bool(args.repro_check)}
    if os.path.exists(meta_path):
        old = PF.jload(meta_path)
        if old.get("spec") != json.loads(json.dumps(spec)):
            diff = sorted(k for k in set(spec) | set(old.get("spec") or {})
                          if json.loads(json.dumps(spec)).get(k) != (old.get("spec") or {}).get(k))
            raise SystemExit("[otune] {} was written with a different spec ({}); use a new "
                             "--outdir".format(meta_path, diff))
    done = {}
    if os.path.exists(jsonl):
        with open(jsonl) as fh:
            for line in fh:
                if line.strip():
                    r = json.loads(line)
                    if r["status"] != "error":
                        done[r["key"]] = r
    todo = []
    for j in jobs:
        if j["key"] in done:
            continue
        jj = dict(j)
        jj.update({"arms": args.arms, "fva_arms": args.fva_arms, "seed": args.seed,
                   "g_meas": g_meas, "isotonic": iso, "prior": slim_prior(prior.get(j["key"])),
                   "repro_check": bool(args.repro_check)})
        todo.append(jj)
    meta = {"spec": spec, "generated": now_iso(), "git_commit": W.git_commit(REPO),
            "design": os.path.join(SR, "PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md"),
            "design_sha256": sha256(os.path.join(SR, "PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md")),
            "g_meas": g_meas, "g_meas_info": g_meas_info, "weights": wrep,
            "est_attach": est_attach, "isotonic": iso,
            "factors": {"P": P_SOURCES, "P_candidates": P_CANDIDATES, "C": C_SHAPES,
                        "G": G_LEVELS, "X": X_LEVELS, "x_values": X_VALUE,
                        "g_values": {"G0": 0.0, "Gmeas": g_meas, "G0.2": 0.2},
                        "shipped": SHIPPED, "d_arm": D_ARM, "min_prob": MIN_PROB,
                        "cost": {"lin": "max(0.01, 1 - p)",
                                 "log": "max(0.01, ln(max(p, 0.01)) / ln(0.01))"}},
            "split": {"rule": "tune iff int(sha256(('1729|' + group_id).encode()).hexdigest(),"
                              " 16) % 2 == 0",
                      "n_jobs_by_half": dict(Counter(half_of(j["group_id"]) for j in jobs))},
            "n_jobs": len(jobs), "n_done_at_start": len(jobs) - len(todo),
            "args": {k: v for k, v in vars(args).items()},
            "libraries": libs,
            "templates": PF.template_provenance(PF.use_templates(args.templates)),
            "script": {"path": os.path.abspath(__file__), "sha256": script_sha},
            "weighted_objective_test_sha256": sha256(os.path.join(HERE,
                                                                  "weighted_objective_test.py")),
            "phenotype_fold_sha256": sha256(os.path.join(HERE, "phenotype_fold.py")),
            "phenotype_arm_sha256": sha256(os.path.join(PF.P2_SCRIPTS, "phenotype_arm.py")),
            "outputs": {"per_amplicon_jsonl": os.path.abspath(jsonl)}}
    if os.path.exists(meta_path):
        old = PF.jload(meta_path)
        meta["generated_first"] = old.get("generated_first", old.get("generated"))
        meta["script_sha256_history"] = sorted(set(old.get("script_sha256_history", [])) |
                                               {old["script"]["sha256"], script_sha})
    with open(meta_path + ".tmp", "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True, default=str)
    os.replace(meta_path + ".tmp", meta_path)
    print("[otune] {} jobs selected ({}), {} done, {} to run; g_meas {:.6f}; weights {} "
          "(stub: {}); {} arms".format(len(jobs), meta["split"]["n_jobs_by_half"],
                                      len(jobs) - len(todo), len(todo), g_meas,
                                      wrep["path"], wrep["stub"], len(args.arms)), flush=True)
    t_run = time.time()
    if todo:
        import multiprocessing as mp
        ctx = mp.get_context("spawn")             # Gurobi environments do not survive fork
        nw = max(1, min(args.workers, len(todo)))
        winit = (args.msd, args.templates, not args.no_integrate_reversals, universe)
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
                    pmok = (r.get("prior_match") or {}).get("ok")
                    rc = (r.get("released_check") or {}).get("pass")
                    print("[otune] {} {} {} {} {}s prior_match={} released={} rss={}".format(
                        r["key"], r["stratum"], r.get("half"), r["status"], r["wall_seconds"],
                        pmok, rc, r.get("worker_maxrss_gb")), flush=True)
            finally:
                if pool is not None:
                    pool.close()
                    pool.join()
    rows = [done[j["key"]] for j in jobs if j["key"] in done]
    import resource
    meta["run"] = {
        "finished": now_iso(), "seconds_this_invocation": round(time.time() - t_run, 1),
        "n_rows": len(rows), "status": dict(Counter(r["status"] for r in rows)),
        "prior_match_fail": [r["key"] for r in rows
                             if not (r.get("prior_match") or {}).get("ok", False)],
        "released_check_fail": [r["key"] for r in rows if r["status"] == "ok" and
                                not (r.get("released_check") or {}).get("pass", False)],
        "parent_maxrss_gb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss /
                                  1048576.0, 3),
        "largest_child_maxrss_gb": round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
                                         / 1048576.0, 3)}
    with open(meta_path + ".tmp", "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True, default=str)
    os.replace(meta_path + ".tmp", meta_path)
    print("[otune] wrote {} rows -> {}; status {}; prior_match failures {}; released_check "
          "failures {}".format(len(rows), jsonl, meta["run"]["status"],
                               len(meta["run"]["prior_match_fail"]),
                               len(meta["run"]["released_check_fail"])), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
