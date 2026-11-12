#!/usr/bin/env python
"""Equivalence gate for prfba_objective.run_prfba (the prFBA-side wrapper of
ModelSEEDpy's MSProbability.prFBA, fix branch fix-prfba-single-model).

QUESTION
  Does the library call with prFBA's tuned keywords (TUNED_PRFBA_KWARGS: ex_weight 1,
  min_prob 0.01, prob_exp 1, missing_probability 0.0, biomass_id "bio1") solve the same LP,
  and return the same optimum, as Stage O's harness copy objective_tuning.prfba_tuned_copy
  at the selected arm est/lin/G0/X1?  And does the released configuration (tuned=False,
  ex_weight 100) reproduce the shipped arm kn/lin/G0/X100?

MODELS (Stage O's, rebuilt)
  Per amplicon, objective_tuning.run_amplicon's steps 1-2, called directly (NOT through
  weighted_objective_test.worker_init / released_prfba(), which exec-edit the released
  prFBA source and raise on the fixed library): phenotype_fold.worker_init with Stage O's
  biochemistry, templates (ModelSEEDTemplates-c12run 21a8084) and gapfiller reversals; the
  complete medium at uptake 100 as the environment; phenotype_fold.gapfill_block of the
  selection's union EC set ("as_shipped"); the role-level kn (weighted_objective_test.
  role_probabilities) and E3 (objective_tuning.est_role_probabilities, from Stage O's
  weights file stage_o_weights_E3.json, attached with objective_tuning.load_weights /
  attach_est).  ModelSEEDpy is the fix worktree (asserted clean at FIX_COMMIT); it differs
  from Stage O's pinned c2e31a0 only in msprobability.py and a new test file (recorded from
  git diff --stat), so model building is Stage O's code.  The rebuilt model is compared with
  Stage O's row (n_rxn, n_gene_backed, n_gapfilled, gapfill n_new_rxn and growth) and the
  harness objective with the row's recorded objective.

ARMS, each on a fresh copy of the gapfilled model
  tuned     (a) harness: probs = est with gapfilled p = 0; M.objective = "bio1";
                objective_tuning.prfba_tuned_copy(M, env, cost_map(probs, "lin"), 1).
            (b) wrapper: prfba_objective.set_reaction_probabilities(M, q_role) (gene-level
                E3; the production rule: max over genes, gene-less and so gapfilled 0.0),
                then run_prfba(M, env) with the tuned keywords (biomass_id sets bio1).
            (c) as (b) with the notes removed from the gapfilled reactions, so
                missing_probability=0.0 prices them (in (b) every reaction is annotated and
                missing_probability is a no-op).
  released  (a) harness at kn/lin/G0/X100: prfba_tuned_copy(M, env, cost_map(kn, "lin"), 100).
            (b) wrapper: set_reaction_probabilities(M, p_role); M.objective = "bio1" (the
                released keywords carry no biomass_id; the harness sets bio1 the same way);
                run_prfba(M, env, tuned=False, ex_weight=100).
  Compared: the wrapper rule's reaction probabilities against the harness's (exact); the
  objective coefficient maps by variable name; the 95% floor (bound and variables); the
  objective value (relative difference, and bitwise); net fluxes (max |diff| by reaction)
  and variable primals (max |diff| by name).
  PASS per comparison: rule probabilities identical, coefficient maps identical, floor
  identical, objective relative difference <= 1e-9.  Flux differences are reported, not
  gated (the optimal face is degenerate), though an identical LP is expected to return the
  identical vertex.

OUTPUT  data/objective_tuning/wrapper_equivalence.json

    PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-prfba \\
        ~/Documents/py_venv/bin/python scripts/objective_wrapper_equivalence.py

Interpreter: /home/freiburger/Documents/py_venv/bin/python (the only interpreter).
"""
import argparse
import gc
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)
D = os.path.join(SR, "data")
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
for p in (HERE, REPO):
    if p not in sys.path:
        sys.path.insert(0, p)
import objective_tuning as OT          # noqa: E402  (no cobra at import time)
import weighted_objective_test as W    # noqa: E402
import phenotype_fold as PF            # noqa: E402

FIX_COMMIT = "09856f5"
STAGE_O_COMMIT = "c2e31a0"
TUNED_ARM = "est/lin/G0/X1"
RELEASED_ARM = "kn/lin/G0/X100"
REL_TOL = 1e-9
# two per stratum; the three rows that rebuild a different LP are excluded
DEFAULT_KEYS = ("2053306.99__V4", "418003.3__V4",               # n1
                "298656.3__V7-V9", "500466.3__V1-V2",            # n2
                "2287.19__V7-V9", "3125783.3__V4-V5_944R",       # n3
                "3133138.3__V4", "2837493.4__V4-V5")             # n4+ (n = 5, n = 4)
BAD_KEYS = {"1674805.4__V3-V4", "1798454.3__V4", "1017270.21__V4-V5"}


def now_iso():
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def git(path, *a):
    return subprocess.run(["git", *a], cwd=path, check=True, capture_output=True,
                          text=True).stdout.strip()


def coefficients(M):
    coefs = M.solver.objective.get_linear_coefficients(M.variables)
    return {v.name: float(c) for v, c in coefs.items() if c != 0}


def floor_of(M):
    """The 95% floor: its bounds, and whether its variables are bio1's (fwd and rev)."""
    c = M.constraints["min_value"]
    names = sorted(v.name for v in c.variables)
    bio = M.reactions.get_by_id("bio1")
    return {"lb": c.lb, "ub": c.ub, "variables": names,
            "on_bio1": set(names) == {bio.forward_variable.name, bio.reverse_variable.name}}


def solution_record(M, sol):
    return {"status": sol.status, "objective": sol.objective_value,
            "bio1": float(sol.fluxes["bio1"]),
            "fluxes": {k: float(v) for k, v in sol.fluxes.items()},
            "primals": {v.name: v.primal for v in M.variables},
            "coef": coefficients(M), "floor": floor_of(M)}


def compare(a, b):
    """Two solution_records: objective, fluxes, variable primals, coefficients, floor."""
    oa, ob = a["objective"], b["objective"]
    rel = abs(oa - ob) / max(1.0, abs(oa), abs(ob))
    fl = max(abs(a["fluxes"][k] - b["fluxes"][k]) for k in a["fluxes"])
    same_keys_f = set(a["fluxes"]) == set(b["fluxes"])
    pr = max(abs(a["primals"][k] - b["primals"].get(k, float("nan"))) for k in a["primals"])
    same_keys_p = set(a["primals"]) == set(b["primals"])
    ca, cb = a["coef"], b["coef"]
    n_differ = sum(1 for k in set(ca) | set(cb) if ca.get(k) != cb.get(k))
    out = {"status": [a["status"], b["status"]],
           "objective": [oa, ob], "objective_rel_diff": rel,
           "objective_bitwise_equal": oa == ob,
           "bio1": [a["bio1"], b["bio1"]],
           "fluxes_max_abs_diff": fl, "fluxes_same_reactions": same_keys_f,
           "fluxes_bitwise_equal": same_keys_f and fl == 0.0,
           "primals_max_abs_diff": pr, "primals_same_variables": same_keys_p,
           "coef_n": [len(ca), len(cb)], "coef_n_differ": n_differ,
           "coef_maps_identical": n_differ == 0,
           "floor": [a["floor"], b["floor"]],
           "floor_identical": a["floor"] == b["floor"]}
    out["floor_on_bio1"] = [a["floor"]["on_bio1"], b["floor"]["on_bio1"]]
    out["pass"] = bool(a["status"] == b["status"] == "optimal" and out["coef_maps_identical"]
                       and out["floor_identical"] and all(out["floor_on_bio1"])
                       and rel <= REL_TOL)
    return out


def check_library():
    libs = W.library_provenance()
    msp = libs.get("modelseedpy") or {}
    path = msp.get("path") or ""
    ok = (path.startswith("/home/freiburger/Documents/ModelSEEDpy-prfba/")
          and (msp.get("git_commit") or "").startswith(FIX_COMMIT)
          and msp.get("git_dirty") is False)
    if not ok:
        raise SystemExit("[wequiv] modelseedpy must be the clean fix worktree at {} (got {} "
                         "at {}, dirty {}); run with PYTHONPATH=/home/freiburger/Documents/"
                         "ModelSEEDpy-prfba".format(FIX_COMMIT, path, msp.get("git_commit"),
                                                    msp.get("git_dirty")))
    root = os.path.dirname(path)
    diff = git(root, "diff", "--stat", "{}..HEAD".format(STAGE_O_COMMIT))
    names = git(root, "diff", "--name-only", "{}..HEAD".format(STAGE_O_COMMIT)).split()
    return libs, {"fix_vs_stage_o_diff_stat": diff, "files_changed": names,
                  "only_msprobability_and_its_test": sorted(names) == sorted(
                      ["modelseedpy/core/msprobability.py",
                       "tests/core/test_msprobability_prfba.py"])}


def build(job, m, env_media, index):
    """objective_tuning.run_amplicon steps 1-2 (verbatim calls); returns the gapfilled model
    and the probability maps."""
    tmpl, role_ecs = m.load_template(job["template"])
    union = set().union(*map(set, job["member_ec"].values()))
    t0 = time.time()
    mdl, block, _ess, _gc, _gg, _ctx, infeas = PF.gapfill_block(
        m, union, job["template"], PF.safe_id(job["key"], "wobj"), env_media[0],
        env_media[1], "as_shipped", None, False)
    t_gf = round(time.time() - t0, 2)
    if infeas or not block["sanity_gate_ok"] or block["growth"] <= PF.GROWTH_THRESH:
        raise RuntimeError("{}: model not scoreable ({}, gate {}, growth {})".format(
            job["key"], infeas, block["sanity_gate_ok"], block["growth"]))
    gapf = {r.id for r in mdl.reactions if r.id[:-1] in set(block["new_rxn"])}
    genes = sorted(g.id for g in mdl.genes)
    p_role, _unknown = W.role_probabilities(genes, role_ecs, job["member_ec"])
    rxn_genes = {r.id: [g.id for g in r.genes] for r in mdl.reactions if r.id[:3] == "rxn"}
    prod = {r: (max((p_role.get(g, 0.0) for g in gs), default=0.0) if gs else 0.0)
            for r, gs in rxn_genes.items()}
    q_role, est_counts = OT.est_role_probabilities(genes, role_ecs, job["est"], index,
                                                   job["isotonic"])
    est = {r: (max((q_role.get(g, 0.0) for g in gs), default=0.0) if gs else 0.0)
           for r, gs in rxn_genes.items()}
    gene_backed = {r.id for r in mdl.reactions if r.id[:3] == "rxn" and r.genes}
    if gapf & gene_backed:
        raise RuntimeError("{} gapfilled reactions are gene-backed".format(len(gapf &
                                                                                gene_backed)))
    info = {"gapfill_seconds": t_gf, "n_new_rxn": block["n_new_rxn"],
            "growth": block["growth"], "n_rxn": len(mdl.reactions),
            "n_rxn_prefixed": len(rxn_genes), "n_gene_backed": len(gene_backed),
            "n_gapfilled": len(gapf), "n_genes": len(genes),
            "n_template_role_genes_kn": len(p_role), "n_template_role_genes_est": len(q_role),
            "est_counts": est_counts,
            "mean_est_gene_backed": (sum(est[r] for r in gene_backed) / len(gene_backed))
            if gene_backed else None,
            "mean_kn_gene_backed": (sum(prod[r] for r in gene_backed) / len(gene_backed))
            if gene_backed else None,
            "objective_on_entry": str(mdl.objective.expression)[:80]}
    return mdl, gapf, p_role, q_role, prod, est, info


def rule_check(rp, harness_probs):
    """The wrapper rule's probability against the harness's, on every rxn reaction."""
    keys = [r for r in rp if r[:3] == "rxn"]
    bad = [r for r in keys if rp[r] != harness_probs[r]]
    return {"n_rxn": len(keys), "n_differ": len(bad), "first_differ": bad[:5],
            "same_reaction_set": set(keys) == set(harness_probs),
            "identical": not bad and set(keys) == set(harness_probs)}


def solve_harness(mdl, env, probs, x):
    M = mdl.copy()
    M.objective = "bio1"
    sol, _u, coef, maxbio = OT.prfba_tuned_copy(M, env, OT.cost_map(probs, "lin"), x)
    rec = solution_record(M, sol)
    rec["maxBioSol"] = maxbio
    rec["harness_coef_equals_model_objective"] = (
        {v.name: float(c) for v, c in coef.items() if c != 0} == rec["coef"])
    OT.release(M)
    del M, sol, _u, coef
    W.purge_caches()
    gc.collect()
    return rec


def solve_wrapper(mdl, env, gene_probs, harness_probs, gapf, tuned, strip_gapfilled=False,
                  set_bio1=False, **overrides):
    import prfba_objective as PO
    M = mdl.copy()
    rp = PO.set_reaction_probabilities(M, gene_probs)
    rule = rule_check(rp, harness_probs)
    n_stripped = 0
    if strip_gapfilled:
        for r in gapf:
            M.reactions.get_by_id(r).notes.pop("probability", None)
            n_stripped += 1
    if set_bio1:
        M.objective = "bio1"
    sol = PO.run_prfba(M, env, tuned=tuned, **overrides)
    rec = solution_record(M, sol)
    rec["rule"] = rule
    rec["n_gapfilled_notes_removed"] = n_stripped
    rec["kwargs"] = PO.prfba_kwargs(tuned, **overrides)
    OT.release(M)
    del M, sol
    W.purge_caches()
    gc.collect()
    return rec


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--keys", default=",".join(DEFAULT_KEYS))
    ap.add_argument("--out", default=os.path.join(D, "objective_tuning",
                                                  "wrapper_equivalence.json"))
    ap.add_argument("--weights", default=os.path.join(D, "pfunc_estimators",
                                                      "stage_o_weights_E3.json"))
    args = ap.parse_args(argv)
    t_all = time.time()
    libs, libdiff = check_library()
    print("[wequiv] modelseedpy {} at {}".format(libs["modelseedpy"]["path"],
                                                 libs["modelseedpy"]["git_commit"]), flush=True)
    meta_o = PF.jload(os.path.join(D, "objective_tuning", "meta.json"))
    so_args = meta_o["args"]
    keys = [k for k in args.keys.split(",") if k]
    bad = sorted(set(keys) & BAD_KEYS)
    if bad:
        raise SystemExit("[wequiv] {} rebuild a different LP in Stage O; choose others".format(
            bad))
    pop = PF.jload(so_args["population"])
    jobs = [j for j in pop["jobs"] if j["key"] in set(keys)]
    if len(jobs) != len(keys):
        raise SystemExit("[wequiv] keys not in the population: {}".format(
            sorted(set(keys) - {j["key"] for j in jobs})))
    jobs.sort(key=lambda j: keys.index(j["key"]))
    if os.path.abspath(args.weights) != os.path.abspath(so_args["weights"]):
        raise SystemExit("[wequiv] Stage O ran with weights {}".format(so_args["weights"]))
    Wt, wrep = OT.load_weights(args.weights, jobs, False)
    iso = (Wt.get("meta") or {}).get("isotonic")
    if iso is not None:
        raise SystemExit("[wequiv] the E3 weights carry an isotonic map; Stage O's did not")
    if wrep["sha256"] != meta_o["weights"]["sha256"]:
        raise SystemExit("[wequiv] weights sha256 differs from Stage O's")
    universe = sorted(PF.norm_ec(PF.jload(so_args["ec_universe"])["ec"]))
    est_attach = OT.attach_est(jobs, Wt, universe)
    index = {e: i for i, e in enumerate(universe)}
    for j in jobs:
        j["isotonic"] = iso
    rows = {}
    with open(os.path.join(D, "objective_tuning", "per_amplicon.jsonl")) as fh:
        for line in fh:
            if line.strip():
                r = json.loads(line)
                if r["key"] in set(keys):
                    rows[r["key"]] = r
    # ---- worker setup: weighted_objective_test.worker_init WITHOUT released_prfba()
    t0 = time.time()
    PF.worker_init(so_args["msd"], True, so_args["templates"],
                   not so_args["no_integrate_reversals"])
    m = PF.pa()
    gf_media, gf_cpds = m.complete_media()
    env = {c: 100.0 for c in gf_cpds}
    t_init = round(time.time() - t0, 1)
    print("[wequiv] worker setup {}s; {} amplicons".format(t_init, len(jobs)), flush=True)
    import cobra
    solver = cobra.Configuration().solver.__name__
    results = []
    for job in jobs:
        ta = time.time()
        key = job["key"]
        row = rows[key]
        res = {"key": key, "stratum": job["stratum"], "n_reps": job["n_reps"],
               "half": OT.half_of(job["group_id"]), "template": job["template"]}
        try:
            mdl, gapf, p_role, q_role, prod, est, info = build(job, m, (gf_media, gf_cpds),
                                                                index)
            res["model"] = info
            pm, pg = row["model"], row["gapfill"]
            res["rebuild_vs_stage_o"] = {
                "n_rxn": [info["n_rxn"], pm["n_rxn"]],
                "n_gene_backed": [info["n_gene_backed"], pm["n_gene_backed"]],
                "n_gapfilled": [info["n_gapfilled"], pm["n_gapfilled"]],
                "n_new_rxn": [info["n_new_rxn"], pg["n_new_rxn"]],
                "growth": [info["growth"], pg["growth"]]}
            res["rebuild_matches_stage_o"] = (
                info["n_rxn"] == pm["n_rxn"] and info["n_gene_backed"] == pm["n_gene_backed"]
                and info["n_gapfilled"] == pm["n_gapfilled"]
                and info["n_new_rxn"] == pg["n_new_rxn"]
                and abs(info["growth"] - pg["growth"]) <= 1e-6)
            probs_est = dict(est)
            probs_kn = dict(prod)
            for r in gapf:
                probs_est[r] = 0.0
                probs_kn[r] = 0.0
            # ---- tuned
            ha = solve_harness(mdl, env, probs_est, OT.X_VALUE["X1"])
            wb = solve_wrapper(mdl, env, q_role, probs_est, gapf, tuned=True)
            wc = solve_wrapper(mdl, env, q_role, probs_est, gapf, tuned=True,
                               strip_gapfilled=True)
            so_t = row["arms"][TUNED_ARM]["objective"]
            res["tuned"] = {
                "arm": TUNED_ARM,
                "harness_vs_wrapper": compare(ha, wb),
                "harness_vs_wrapper_missing_probability": compare(ha, wc),
                "wrapper_rule_vs_harness_probs": wb["rule"],
                "n_gapfilled_notes_removed": wc["n_gapfilled_notes_removed"],
                "wrapper_kwargs": wb["kwargs"],
                "harness_coef_equals_model_objective": ha["harness_coef_equals_model_objective"],
                "stage_o_objective": so_t,
                "harness_vs_stage_o_rel_diff": abs(ha["objective"] - so_t) /
                max(1.0, abs(so_t)),
                "harness_vs_stage_o_bitwise": ha["objective"] == so_t,
                "wrapper_vs_stage_o_rel_diff": abs(wb["objective"] - so_t) /
                max(1.0, abs(so_t))}
            # ---- released configuration
            hr = solve_harness(mdl, env, probs_kn, OT.X_VALUE["X100"])
            wr = solve_wrapper(mdl, env, p_role, probs_kn, gapf, tuned=False, set_bio1=True,
                               ex_weight=100)
            so_r = row["arms"][RELEASED_ARM]["objective"]
            so_rel = (row.get("released_check") or {}).get("objective_released")
            res["released"] = {
                "arm": RELEASED_ARM,
                "harness_vs_wrapper": compare(hr, wr),
                "wrapper_rule_vs_harness_probs": wr["rule"],
                "wrapper_kwargs": wr["kwargs"],
                "objective_set_by_caller": "bio1",
                "stage_o_objective": so_r,
                "stage_o_released_check_objective": so_rel,
                "harness_vs_stage_o_rel_diff": abs(hr["objective"] - so_r) /
                max(1.0, abs(so_r)),
                "harness_vs_stage_o_bitwise": hr["objective"] == so_r,
                "wrapper_vs_stage_o_released_check_bitwise": (wr["objective"] == so_rel)
                if so_rel is not None else None}
            res["status"] = "ok"
            res["pass"] = bool(
                res["tuned"]["harness_vs_wrapper"]["pass"]
                and res["tuned"]["harness_vs_wrapper_missing_probability"]["pass"]
                and res["tuned"]["wrapper_rule_vs_harness_probs"]["identical"]
                and res["released"]["harness_vs_wrapper"]["pass"]
                and res["released"]["wrapper_rule_vs_harness_probs"]["identical"])
            del mdl
        except Exception as exc:                      # noqa: BLE001
            import traceback
            res["status"] = "error"
            res["pass"] = False
            res["error"] = "{}: {}".format(type(exc).__name__, str(exc)[:500])
            res["traceback"] = traceback.format_exc()[-3000:]
        res["seconds"] = round(time.time() - ta, 1)
        W.purge_caches()
        gc.collect()
        t = res.get("tuned", {}).get("harness_vs_wrapper", {})
        r_ = res.get("released", {}).get("harness_vs_wrapper", {})
        print("[wequiv] {} {} {} pass={} tuned rel={} flux={} | released rel={} flux={} | "
              "rebuild={} {}s".format(key, res["stratum"], res["status"], res["pass"],
                                      t.get("objective_rel_diff"), t.get("fluxes_max_abs_diff"),
                                      r_.get("objective_rel_diff"),
                                      r_.get("fluxes_max_abs_diff"),
                                      res.get("rebuild_matches_stage_o"), res["seconds"]),
              flush=True)
        results.append(res)

    def agg(path):
        vals = []
        for r in results:
            x = r
            for k in path:
                x = (x or {}).get(k)
            vals.append(x)
        return vals

    ok = [r for r in results if r["status"] == "ok"]
    summary = {
        "n_amplicons": len(results), "n_ok": len(ok),
        "strata": {s: sum(1 for r in results if r["stratum"] == s)
                   for s in ("n1", "n2", "n3", "n4+")},
        "all_pass": bool(results) and all(r["pass"] for r in results),
        "rebuild_matches_stage_o_all": all(r.get("rebuild_matches_stage_o") for r in results),
        "tuned": {
            "max_objective_rel_diff": max(agg(("tuned", "harness_vs_wrapper",
                                               "objective_rel_diff")), default=None),
            "objective_bitwise_all": all(agg(("tuned", "harness_vs_wrapper",
                                              "objective_bitwise_equal"))),
            "max_fluxes_abs_diff": max(agg(("tuned", "harness_vs_wrapper",
                                            "fluxes_max_abs_diff")), default=None),
            "fluxes_bitwise_all": all(agg(("tuned", "harness_vs_wrapper",
                                           "fluxes_bitwise_equal"))),
            "coef_maps_identical_all": all(agg(("tuned", "harness_vs_wrapper",
                                                "coef_maps_identical"))),
            "missing_probability_variant_pass_all": all(
                agg(("tuned", "harness_vs_wrapper_missing_probability", "pass"))),
            "rule_identical_all": all(agg(("tuned", "wrapper_rule_vs_harness_probs",
                                           "identical"))),
            "max_harness_vs_stage_o_rel_diff": max(agg(("tuned",
                                                        "harness_vs_stage_o_rel_diff")),
                                                   default=None),
            "harness_vs_stage_o_bitwise_all": all(agg(("tuned",
                                                       "harness_vs_stage_o_bitwise")))},
        "released": {
            "max_objective_rel_diff": max(agg(("released", "harness_vs_wrapper",
                                               "objective_rel_diff")), default=None),
            "objective_bitwise_all": all(agg(("released", "harness_vs_wrapper",
                                              "objective_bitwise_equal"))),
            "max_fluxes_abs_diff": max(agg(("released", "harness_vs_wrapper",
                                            "fluxes_max_abs_diff")), default=None),
            "fluxes_bitwise_all": all(agg(("released", "harness_vs_wrapper",
                                           "fluxes_bitwise_equal"))),
            "coef_maps_identical_all": all(agg(("released", "harness_vs_wrapper",
                                                "coef_maps_identical"))),
            "rule_identical_all": all(agg(("released", "wrapper_rule_vs_harness_probs",
                                           "identical"))),
            "max_harness_vs_stage_o_rel_diff": max(agg(("released",
                                                        "harness_vs_stage_o_rel_diff")),
                                                   default=None),
            "harness_vs_stage_o_bitwise_all": all(agg(("released",
                                                       "harness_vs_stage_o_bitwise"))),
            "wrapper_vs_stage_o_released_check_bitwise_all": all(
                agg(("released", "wrapper_vs_stage_o_released_check_bitwise")))}}
    import prfba_objective as PO
    out = {
        "generated": now_iso(),
        "script": {"path": os.path.abspath(__file__), "sha256": W.sha256(os.path.abspath(
            __file__))},
        "wrapper": {"path": os.path.abspath(PO.__file__),
                    "sha256": W.sha256(os.path.abspath(PO.__file__)),
                    "TUNED_PRFBA_KWARGS": PO.TUNED_PRFBA_KWARGS,
                    "RELEASED_PRFBA_KWARGS": PO.RELEASED_PRFBA_KWARGS},
        "prfba_git_commit": W.git_commit(REPO),
        "libraries": libs, "modelseedpy_fix_vs_stage_o": libdiff, "solver": solver,
        "stage_o": {"meta": os.path.join(D, "objective_tuning", "meta.json"),
                    "modelseedpy": meta_o["libraries"]["modelseedpy"],
                    "templates": so_args["templates"], "msd": so_args["msd"],
                    "integrate_reversals": not so_args["no_integrate_reversals"],
                    "weights": wrep["path"], "weights_sha256": wrep["sha256"],
                    "est_attach": est_attach,
                    "objective_tuning_sha256": W.sha256(os.path.join(HERE,
                                                                     "objective_tuning.py")),
                    "objective_tuning_sha256_in_stage_o_meta": meta_o["script"]["sha256"]},
        "environment": "complete medium (phenotype_arm.complete_media), every compound at "
                       "uptake 100, as a {{cpd: 100.0}} dict; {} compounds".format(len(env)),
        "pass_rule": "per comparison: both optimal, objective coefficient maps identical by "
                     "variable name, the 95% floor identical, objective relative difference "
                     "<= {}; per amplicon also the wrapper rule's reaction probabilities "
                     "identical to the harness's (tuned and released) and the "
                     "missing_probability variant passing; fluxes reported, not gated"
                     .format(REL_TOL),
        "worker_setup_seconds": t_init, "wall_seconds": round(time.time() - t_all, 1),
        "summary": summary, "amplicons": results}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out + ".tmp", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True, default=str)
    os.replace(args.out + ".tmp", args.out)
    print("[wequiv] summary " + json.dumps(summary, sort_keys=True, default=str), flush=True)
    print("[wequiv] wrote {}".format(args.out), flush=True)
    print("[wequiv] DONE all_pass={}".format(summary["all_pass"]), flush=True)
    return 0 if summary["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
