#!/usr/bin/env python
"""Which reactions let a gapfilled phenotype-arm model make ATP (or biomass) from water
and protons alone?  Runs on the models egc_rebuild.py writes.

For every model it
  1. reproduces the sanity gate (phenotype_arm.sanity_gate semantics: every exchange
     closed for uptake and open to 100 for excretion, water and H+ at +/-100) and reads
     the maximum ATPM_c0 and bio1 flux;
  2. replays the gate on each reconstruction stage (draft, after media evaluation, after
     the ATP-correction gapfill, after genome-scale expansion, final), by restricting
     the final model to that stage's reactions and bounds, so the leak is attributed to
     the stage that first opens it;
  3. enumerates minimum-cardinality reaction supports of the leak (MILP, binary
     indicators, big-M = 1000, no-good cuts), i.e. the smallest flux-carrying networks
     that make ATP from water, in increasing size;
  4. finds the minimum set of NON-DRAFT reactions (added by ATP correction, expansion or
     MSGapfill, or whose bounds were widened after the draft) whose removal -- or bound
     reversion -- abolishes the leak: an iterate-to-convergence minimum hitting set over
     the supports, verified by LP, which is exact for the restricted cut;
  5. annotates every implicated reaction with its stage, bounds, equation with compound
     names, the cobra mass/charge balance, and the ModelSEED database record.

    python egc_analyze.py --models data/egc_analysis/models --out data/egc_analysis
"""
import argparse
import glob
import json
import os
import re
import sys
import time

MSD = "/home/freiburger/Documents/ModelSEEDDatabase/Biochemistry"
THRESH = 1e-6
U = 1000.0
STAGES = ["draft", "evaluated", "atp_gapfill", "expanded", "final"]


def load_db():
    import csv
    rx, cp = {}, {}
    for f in sorted(glob.glob(os.path.join(MSD, "reaction_*.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rx[r["id"]] = {k: r.get(k) for k in ("name", "equation", "definition",
                                                      "reversibility", "status",
                                                      "is_transport", "is_obsolete",
                                                      "ec_numbers", "linked_reaction")}
    for f in sorted(glob.glob(os.path.join(MSD, "compound_*.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cp[r["id"]] = {k: r.get(k) for k in ("name", "formula", "charge")}
    return rx, cp


def empty_medium(model):
    for r in model.reactions:
        if r.id.startswith("EX_"):
            r.bounds = (0.0, 100.0)
    for c in ("cpd00001_e0", "cpd00067_e0"):
        ex = "EX_" + c
        if ex in model.reactions:
            model.reactions.get_by_id(ex).bounds = (-100.0, 100.0)


def maxflux(model, rid):
    with model:
        model.objective = rid
        v = model.slim_optimize(error_value=0.0)
    return 0.0 if (v is None or v != v) else float(v)


def provenance(stages, rid, final_bounds):
    first = next((s for s in STAGES if rid in stages[s]), None)
    if first is None:
        return "gate_exchange"
    if first == "draft":
        d = tuple(stages["draft"][rid])
        widened = (final_bounds[0] < d[0] - 1e-9) or (final_bounds[1] > d[1] + 1e-9)
        changed = tuple(final_bounds) != d
        return "draft_widened" if widened else ("draft_narrowed" if changed else "draft")
    return {"evaluated": "evaluation", "atp_gapfill": "atp_correction_gapfill",
            "expanded": "expansion", "final": "msgapfill"}[first]


def restrict_to_stage(model, stages, stage):
    snap = stages[stage]
    for r in model.reactions:
        if r.id.startswith("EX_"):
            continue
        if r.id in snap:
            r.bounds = tuple(snap[r.id])
        else:
            r.bounds = (0.0, 0.0)


def min_support(model, obj, delta, cuts=(), timeout=300):
    """Minimum number of non-exchange reactions carrying flux with obj >= delta."""
    prob = model.problem
    with model:
        model.reactions.get_by_id(obj).lower_bound = delta
        ys, cons = {}, []
        for r in model.reactions:
            if r.id.startswith("EX_") or (r.lower_bound == 0 and r.upper_bound == 0):
                continue
            y = prob.Variable("y_" + r.id, type="binary")
            ys[r.id] = y
            cons.append(prob.Constraint(r.flux_expression - U * y, ub=0))
            cons.append(prob.Constraint(r.flux_expression + U * y, lb=0))
        model.add_cons_vars(list(ys.values()) + cons)
        for S in cuts:
            vs = [ys[x] for x in S if x in ys]
            if vs:
                model.add_cons_vars(prob.Constraint(sum(vs), ub=len(vs) - 1))
        model.objective = prob.Objective(sum(ys.values()), direction="min")
        model.solver.configuration.timeout = timeout
        st = model.solver.optimize()
        if st != "optimal":
            return None, st
        sup = sorted(k for k, y in ys.items() if y.primal > 0.5)
        flux = {k: round(model.reactions.get_by_id(k).flux, 6) for k in sup}
        return {"reactions": sup, "flux": flux}, st


def leak(model, obj):
    return maxflux(model, obj) > THRESH


def min_cut(model, obj, candidates, prov, stages, max_iter=60):
    """Minimum set of candidate reactions whose removal (non-draft) or reversion to draft
    bounds (draft_widened) abolishes the leak.  Exact for this restricted cut: supports
    are generated in the current residual model until it no longer leaks."""
    import itertools  # noqa: F401
    supports, cut, it = [], [], 0
    while it < max_iter:
        it += 1
        with model:
            apply_cut(model, cut, prov, stages)
            if not leak(model, obj):
                return {"cut": cut, "n_supports": len(supports), "iterations": it,
                        "status": "exact"}
            mx = maxflux(model, obj)
            sup, st = min_support(model, obj, min(1.0, mx / 2))
        if sup is None:
            return {"cut": cut, "status": "support_failed:" + st, "iterations": it}
        cand = [r for r in sup["reactions"] if r in candidates]
        if not cand:
            return {"cut": None, "status": "leak_without_candidates",
                    "draft_only_support": sup["reactions"], "iterations": it}
        supports.append(cand)
        cut = hitting_set(supports)
    return {"cut": cut, "status": "max_iter", "iterations": it}


def apply_cut(model, cut, prov, stages):
    for rid in cut:
        r = model.reactions.get_by_id(rid)
        if prov.get(rid, "").startswith("draft"):
            r.bounds = tuple(stages["draft"][rid])
        else:
            r.bounds = (0.0, 0.0)


def hitting_set(sets):
    import gurobipy as gp
    items = sorted({x for s in sets for x in s})
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = {i: m.addVar(vtype=gp.GRB.BINARY) for i in items}
    for s in sets:
        m.addConstr(gp.quicksum(x[i] for i in s) >= 1)
    m.setObjective(gp.quicksum(x.values()), gp.GRB.MINIMIZE)
    m.optimize()
    return sorted(i for i in items if x[i].X > 0.5)


def describe(model, rid, prov, stages, db, cp, info):
    r = model.reactions.get_by_id(rid)
    base = re.sub(r"_[a-z]\d+$", "", rid)
    mb = {k: round(v, 6) for k, v in r.check_mass_balance().items()}
    eq = " + ".join("{} {}".format(-c, met.name or met.id) for met, c in r.metabolites.items() if c < 0)
    eq += " => " if r.lower_bound >= 0 else " <=> "
    eq += " + ".join("{} {}".format(c, met.name or met.id) for met, c in r.metabolites.items() if c > 0)
    stage_bounds = {s: stages[s].get(rid) for s in STAGES}
    d = db.get(base, {})
    return {"id": rid, "name": r.name, "stage": prov.get(rid), "bounds": list(r.bounds),
            "stage_bounds": stage_bounds, "equation": eq,
            "mass_charge_imbalance": mb, "gene_reaction_rule": r.gene_reaction_rule,
            "msgapfill_direction": info.get("gapfill_directions", {}).get(rid),
            "db_status": d.get("status"), "db_reversibility": d.get("reversibility"),
            "db_name": d.get("name"), "db_definition": d.get("definition"),
            "db_is_transport": d.get("is_transport"), "db_obsolete": d.get("is_obsolete")}


def analyze(path):
    import cobra.io
    try:
        import gurobipy
        gurobipy.setParam("Threads", 1)
        gurobipy.setParam("OutputFlag", 0)
    except Exception:
        pass
    t0 = time.time()
    key = os.path.basename(path).replace(".model.json", "")
    meta = json.load(open(path.replace(".model.json", ".stages.json")))
    stages, info = meta["stages"], meta["info"]
    model = cobra.io.load_json_model(path)
    model.solver = "gurobi"
    db, cp = DB
    prov = {r.id: provenance(stages, r.id, r.bounds) for r in model.reactions}
    empty_medium(model)
    out = {"key": key, "genome_id": meta["genome_id"], "arm": meta["arm"],
           "pipeline": meta["pipeline"], "template": meta["template"],
           "stored_gate_ok": info["sanity_gate_ok"], "stored_empty_atp": info["empty_atp"],
           "stored_empty_biomass": info["empty_biomass"],
           "n_reactions": len(model.reactions),
           "n_by_stage": {}, "max_empty": {}, "stage_leak": {}}
    for p in prov.values():
        out["n_by_stage"][p] = out["n_by_stage"].get(p, 0) + 1
    for obj in ("ATPM_c0", "bio1"):
        out["max_empty"][obj] = round(maxflux(model, obj), 6)
    for s in STAGES:
        with model:
            restrict_to_stage(model, stages, s)
            out["stage_leak"][s] = {o: round(maxflux(model, o), 6) for o in ("ATPM_c0", "bio1")}
    candidates = {rid for rid, p in prov.items()
                  if p not in ("draft", "draft_narrowed", "gate_exchange")
                  and not rid.startswith("EX_")}
    implicated = set()
    for obj in ("ATPM_c0", "bio1"):
        mx = out["max_empty"][obj]
        if mx <= THRESH:
            continue
        sups, cuts = [], []
        for _ in range(12):                       # minimum supports, increasing size
            sup, st = min_support(model, obj, min(1.0, mx / 2), cuts=cuts)
            if sup is None:
                break
            sups.append(sup)
            cuts.append(sup["reactions"])
        mc_all = min_cut(model, obj, candidates, prov, stages)
        msg_only = {r for r in candidates if prov[r] == "msgapfill"}
        mc_msg = min_cut(model, obj, msg_only, prov, stages)
        ko = {}
        for rid in sorted({x for s in sups for x in s["reactions"]} & candidates):
            with model:
                apply_cut(model, [rid], prov, stages)
                ko[rid] = round(maxflux(model, obj), 6)
        out[obj] = {"supports": sups, "single_knockout_leak": ko,
                    "min_cut_nondraft": mc_all, "min_cut_msgapfill_only": mc_msg}
        for s in sups:
            implicated |= set(s["reactions"])
        implicated |= set(mc_all.get("cut") or []) | set(mc_msg.get("cut") or [])
    out["reactions"] = {rid: describe(model, rid, prov, stages, db, cp, info)
                        for rid in sorted(implicated)}
    out["secs"] = round(time.time() - t0, 1)
    return out


DB = None


def init():
    global DB
    DB = load_db()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--only-gated", action="store_true")
    a = ap.parse_args()
    paths = sorted(glob.glob(os.path.join(a.models, "*.model.json")))
    if a.only_gated:
        paths = [p for p in paths if not json.load(open(
            p.replace(".model.json", ".stages.json")))["info"]["sanity_gate_ok"]]
    print("[egc] analysing {} models".format(len(paths)), flush=True)
    import multiprocessing as mp
    res = []
    with mp.get_context("spawn").Pool(a.workers, initializer=init) as pool:
        for r in pool.imap_unordered(analyze, paths):
            res.append(r)
            print("[egc] {} max_empty {} stage_leak_first {} cut {}".format(
                r["key"], r["max_empty"],
                next((s for s in STAGES if r["stage_leak"][s]["ATPM_c0"] > THRESH
                      or r["stage_leak"][s]["bio1"] > THRESH), None),
                (r.get("ATPM_c0") or {}).get("min_cut_nondraft", {}).get("cut")), flush=True)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "egc_per_model.json"), "w") as fh:
        json.dump(sorted(res, key=lambda r: r["key"]), fh, indent=1)
    print("[egc] wrote", os.path.join(a.out, "egc_per_model.json"))


if __name__ == "__main__":
    sys.exit(main())
