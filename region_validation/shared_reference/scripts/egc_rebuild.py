#!/usr/bin/env python
"""Rebuild the phenotype arm's gapfilled models with stage provenance, for the
energy-generating-cycle (EGC) analysis of the C12 sanity-gate failures.

The phenotype arm (phenotype_fold.py, --readouts atp,growth,essentiality,structure)
gates a gapfilled model when it makes biomass or ATP on water + protons alone
(phenotype_arm.sanity_gate).  Its rows keep only the final MSGapfill solution's reaction
ids, not the model, and not the reactions the ATP-correction step adds before it.  This
script rebuilds every gapfilled model of the selected genomes -- truth first, exactly as
run_genome() orders them, so the `controlled` pipeline reuses the truth model's media and
tests -- through a copy of gapfill_block() that snapshots the model after each stage:

  draft        MSBuilder.build_base_model + add_atpm
  atp_gapfill  after MSATPCorrection.apply_growth_media_gapfilling()
  expanded     after MSATPCorrection.expand_model_to_genome_scale()
  final        after integrating the MSGapfill solution

It writes, per model, the final cobra model (after the gate's ensure_exchanges, i.e. the
object the gate tested), the stage snapshots (reaction id -> bounds), the MSGapfill
solution with directions, the gate's two numbers, and a reproduction check of the
stored row (new_rxn and sanity_gate_ok must match).

    python egc_rebuild.py --genomes all_gated_and_siblings --outdir data/egc_analysis/models
"""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import phenotype_fold as PF  # noqa: E402

SR = os.path.abspath(os.path.join(HERE, ".."))
D = os.path.join(SR, "data")
GF_ARGS = [
    "--cells", os.path.join(D, "cells.csv.gz"),
    "--fold-plan", os.path.join(D, "fold_plan.json"),
    "--panel", os.path.join(D, "panel_eligible.json"),
    "--prfba", os.path.join(D, "prfba"),
    "--picrust2", os.path.join(D, "picrust2"),
    "--density", "1.0", "--select-mode", "exact_tie",
    "--outdir", os.path.join(D, "phenotype_gf"),
    "--summary", os.path.join(D, "phenotype_gf", "phenotype_summary.json"),
    "--n-genomes", "200", "--readouts", "atp,growth,essentiality,structure"]


# Counterfactual gapfill variants (2026-10-20 EGC analysis): the five reactions whose
# addition by MSGapfill re-opens an ATP-from-water leak in the 44 gated models.
CULPRITS = ["rxn10966", "rxn09180", "rxn08167", "rxn00381", "rxn00084"]
VARIANTS = {
    "none": {"blacklist": [], "block_dir": {}},
    # remove the five from the gapfill universe entirely
    "blacklist5": {"blacklist": CULPRITS, "block_dir": {}},
    # block only the thermodynamically infeasible (reverse) direction of the four the
    # gapfiller added backwards; rxn09180 is added forward, so it is blacklisted
    "dirfix": {"blacklist": ["rxn09180"],
               "block_dir": {"rxn10966": "<", "rxn08167": "<", "rxn00381": "<",
                             "rxn00084": "<"}},
    # correct what is thermodynamically wrong instead of removing reactions: the
    # database-irreversible directions for the two draft-level culprits (catalase is '<' in
    # every template against database '>' and dG' -45 kcal/mol; L-aspartate:NAD+
    # oxidoreductase is '>' against database '<') and for the three reactions MSGapfill
    # added backwards, and integrate the gapfiller's `reversed` set, which
    # MSBuilder.build_and_gapfill (and phenotype_fold, which inlines it) drops
    "principled": {"blacklist": [], "block_dir": {},
                   "template_dir": {"rxn00006": ">", "rxn05119": "<", "rxn10966": ">",
                                    "rxn00381": ">", "rxn00084": ">"},
                   "integrate_reversed": True},
    # + the NADP isoform of L-aspartate dehydrogenase, which carries the residual leak
    # in the rxn08167 family under "principled" (template '>' against database '<')
    "principled2": {"blacklist": [], "block_dir": {},
                    "template_dir": {"rxn00006": ">", "rxn05119": "<", "rxn05117": "<",
                                     "rxn10966": ">", "rxn00381": ">", "rxn00084": ">"},
                    "integrate_reversed": True},
    # systematic: every template reaction the ModelSEED database marks irreversible takes
    # the database's direction (draft bounds and GapfillDirection), rxn10966 is forced
    # forward (database '?'; the coupled NADP-dependent reaction rxn00572 is '<'), and the
    # gapfiller's reversals are integrated
    "dbdir": {"blacklist": [], "block_dir": {}, "template_dir": {"rxn10966": ">"},
              "db_align": True, "integrate_reversed": True},
    # template files as given (use with --template-dir), plus the gapfiller's reversals
    "reversed_only": {"blacklist": [], "block_dir": {}, "integrate_reversed": True},
}
VARIANT = {"name": "none"}


def apply_variant(name):
    """GapfillingPkg reads the module-level base_blacklist ({id: direction to block};
    '<' zeroes the lower bound, '>' the upper, '=' both) when it builds the gapfill
    universe, so the direction blocks are patched in per worker process."""
    import modelseedpy.fbapkg.gapfillingpkg as gpkg
    v = VARIANTS[name]
    for rid, d in v["block_dir"].items():
        gpkg.base_blacklist[rid] = d
    if v.get("db_align"):
        import csv, glob as _glob
        dbrev = {}
        for f in sorted(_glob.glob("/home/freiburger/Documents/ModelSEEDDatabase/"
                                   "Biochemistry/reaction_*.tsv")):
            with open(f) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r["reversibility"] in (">", "<"):
                        dbrev[r["id"]] = r["reversibility"]
        v = dict(v)
        v["template_dir"] = {**dbrev, **v["template_dir"]}
    if v.get("template_dir"):
        m = PF.pa()
        bounds = {">": (0.0, 1000.0), "<": (-1000.0, 0.0), "=": (-1000.0, 1000.0)}
        for key in ("gramneg", "grampos", "archaea", "core"):
            tmpl = m.load_template(key)[0]
            for rid, d in v["template_dir"].items():
                tid = rid + "_c"
                if tid in tmpl.reactions:
                    r = tmpl.reactions.get_by_id(tid)
                    r.lower_bound, r.upper_bound = bounds[d]
                    r.GapfillDirection = d
    VARIANT["name"] = name
    return v


def snap(model):
    return {r.id: [float(r.lower_bound), float(r.upper_bound)] for r in model.reactions}


def gapfill_traced(m, ec, tmpl_key, model_id, gf_media, gf_cpds, pipeline, truth_ctx):
    """phenotype_fold.gapfill_block(), statement for statement, plus snapshots.
    Essentiality is skipped: it does not touch the model the gate tests."""
    from modelseedpy import MSGapfill
    from modelseedpy.core.msbuilder import MSBuilder
    from modelseedpy.core.msatpcorrection import MSATPCorrection, load_default_medias
    from modelseedpy.core.msmodel import get_reaction_constraints_from_direction

    tmpl, role_ecs = m.load_template(tmpl_key)
    core_t = m.load_template("core")[0]
    g = m.genome_from_ec(sorted(ec), role_ecs, model_id)
    b = MSBuilder(g, tmpl, ontology_term="RAST", template_core=core_t)
    base = b.build_base_model(model_id, "0", allow_all_non_grp_reactions=False,
                              annotate_with_rast=False)
    MSBuilder.add_atpm(base)
    stages = {"draft": snap(base)}
    atpc = MSATPCorrection(base, core_t, load_default_medias(), compartment="c0",
                           atp_hydrolysis_id="ATPM_c0", load_default_medias=False)
    media_eval = atpc.evaluate_growth_media()
    controlled = (pipeline == "controlled" and truth_ctx is not None)
    if controlled:
        want = set(truth_ctx["selected_media_ids"])
        atpc.selected_media = [md for md in atpc.media_gapfill_stats if md.id in want]
    else:
        atpc.determine_growth_media()
    stages["evaluated"] = snap(base)
    atpc.apply_growth_media_gapfilling()
    stages["atp_gapfill"] = snap(base)
    atpc.expand_model_to_genome_scale()
    stages["expanded"] = snap(base)
    tests = truth_ctx["tests"] if controlled else atpc.build_tests()

    gf = MSGapfill(base, default_gapfill_templates=[tmpl], test_conditions=tests,
                   default_target="bio1",
                   blacklist=list(VARIANTS[VARIANT["name"]]["blacklist"]))
    infeasible = ""
    try:
        res = gf.run_gapfilling(gf_media)
    except Exception as exc:                       # noqa: BLE001
        res, infeasible = None, "gapfill_exception:{}".format(type(exc).__name__)
    gap_sol, gap_dir, reversed_ = {}, {}, {}
    if res and res.get("new"):
        for rid, d in res["new"].items():
            if rid[:-1] in tmpl.reactions:
                gap_sol[rid[:-1]] = get_reaction_constraints_from_direction(d)
                gap_dir[rid] = d
    if res and res.get("reversed"):
        reversed_ = {rid: d for rid, d in res["reversed"].items()}
    elif res is None and not infeasible:
        infeasible = ("gapfill_infeasible_under_truth_tests" if controlled
                      else "gapfill_infeasible")
    mdl = base.copy()
    MSBuilder.integrate_gapfill_solution(tmpl, mdl, gap_sol)
    if VARIANTS[VARIANT["name"]].get("integrate_reversed"):
        for rid, d in reversed_.items():
            if rid in mdl.reactions:
                r = mdl.reactions.get_by_id(rid)
                if d in (">", "="):
                    r.upper_bound = max(r.upper_bound, 1000.0)
                if d in ("<", "="):
                    r.lower_bound = min(r.lower_bound, -1000.0)
    stages["final"] = snap(mdl)
    gate = m.sanity_gate(mdl)
    yc = m.yield_check(mdl, gf_cpds)
    growth = m.growth_vector(mdl, {"complete": gf_cpds}, thresh=PF.GROWTH_THRESH)
    ctx = {"selected_media_ids": [md.id for md in (atpc.selected_media or [])],
           "tests": tests, "media_eval": media_eval}
    info = {"n_new_rxn": len(gap_sol), "new_rxn": sorted(gap_sol),
            "gapfill_directions": gap_dir, "gapfill_reversed_not_integrated": reversed_,
            "sanity_gate_ok": bool(gate["pass"]), "empty_biomass": gate["empty_biomass"],
            "empty_atp": gate["empty_atp"],
            "yield_ratio": (None if yc["ratio_carbon"] in (float("inf"), float("-inf"))
                            else round(float(yc["ratio_carbon"]), 6)),
            "selected_media_ids": ctx["selected_media_ids"],
            "tests": [{"media": getattr(t.get("media"), "id", None),
                       "objective": t.get("objective"), "is_max_threshold":
                       t.get("is_max_threshold"), "threshold": t.get("threshold")}
                      for t in (tests or [])],
            "infeasible": infeasible, "variant": VARIANT["name"],
            "growth_complete": float(growth.get("complete", 0.0))}
    return mdl, stages, info, ctx


def run_one(job):
    """Truth first (gapfilled once, aliased to both pipelines, as run_genome does), then
    every requested (arm, pipeline)."""
    import cobra.io
    if job.get("_template_dir"):
        # read every template from another directory (e.g. a patched ModelSEEDTemplates
        # checkout); the worker's cached core template is dropped and reloaded from it
        m0 = PF.pa()
        m0.MST = job["_template_dir"].rstrip("/") + "/"
        m0._tcache.clear()
    apply_variant(job.get("_variant", "none"))
    t_all = time.time()
    m = PF.pa()
    tmpl_key = job["template"]
    gf_media, gf_cpds = m.complete_media()
    out_dir = job["_outdir"]
    os.makedirs(out_dir, exist_ok=True)
    ec_by_arm = {a: PF.norm_ec(v) for a, v in job["arms"].items()}
    results, truth_ctx = [], {}
    pipes = job["pipelines"]
    order = ["truth"] + [a for a in job["arms"] if a != "truth"]
    for arm in order:
        if arm not in job["_want_arms"]:
            continue
        for pipe in pipes:
            if arm == "truth" and truth_ctx:
                continue                      # identical build, aliased in run_genome
            ctx = truth_ctx.get(pipe) if arm != "truth" else None
            key = "{}__{}__{}".format(job["genome_id"], arm,
                                      "shared" if arm == "truth" else pipe)
            t0 = time.time()
            try:
                mdl, stages, info, newctx = gapfill_traced(
                    m, ec_by_arm[arm], tmpl_key,
                    PF.safe_id(job["genome_id"], arm, pipe), gf_media, gf_cpds, pipe, ctx)
            except Exception as exc:          # noqa: BLE001
                results.append({"key": key, "error": "{}: {}".format(
                    type(exc).__name__, exc)[:400]})
                if arm == "truth":
                    return {"genome_id": job["genome_id"], "results": results}
                continue
            if arm == "truth":
                for p in pipes:
                    truth_ctx[p] = newctx
            cobra.io.save_json_model(mdl, os.path.join(out_dir, key + ".model.json"))
            with open(os.path.join(out_dir, key + ".stages.json"), "w") as fh:
                json.dump({"key": key, "genome_id": job["genome_id"], "arm": arm,
                           "pipeline": pipe, "template": tmpl_key, "info": info,
                           "stages": stages, "secs": round(time.time() - t0, 1)}, fh)
            results.append({"key": key, "arm": arm, "pipeline": pipe,
                            "n_new_rxn": info["n_new_rxn"], "new_rxn": info["new_rxn"],
                            "sanity_gate_ok": info["sanity_gate_ok"],
                            "empty_biomass": info["empty_biomass"],
                            "empty_atp": info["empty_atp"], "infeasible": info["infeasible"],
                            "secs": round(time.time() - t0, 1)})
    return {"genome_id": job["genome_id"], "results": results,
            "wall_seconds": round(time.time() - t_all, 1)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", default="gated",
                    help="'gated' (genomes with any gated row), 'all', or a comma list")
    ap.add_argument("--arms", default="all", help="'all' or a comma list of arm names")
    ap.add_argument("--outdir", default=os.path.join(D, "egc_analysis", "models"))
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--variant", default="none", choices=sorted(VARIANTS))
    ap.add_argument("--template-dir", default=None,
                    help="read the ModelSEED templates from this directory instead")
    a = ap.parse_args()

    rows = [json.loads(l) for l in open(os.path.join(D, "phenotype_gf", "per_genome.jsonl"))]
    gated = sorted({r["src_genome_id"] for r in rows if r["status"] == "gated"})
    want = (gated if a.genomes == "gated" else
            None if a.genomes == "all" else a.genomes.split(","))

    args = PF.parse_args(GF_ARGS)
    panel, fold_plan = PF.jload(args.panel), PF.jload(args.fold_plan)
    functions = PF.jload(args.functions)
    universe = PF.ec_universe(os.path.join(os.path.dirname(args.panel),
                                           "panel_ec_universe.json"), panel, functions)
    cells = PF.load_cells(args.cells, args.density)
    prf, _ = PF.load_prfba(args.prfba, args.select_mode, args.density)
    p2, _ = PF.load_picrust2(args.picrust2, args.density)
    jobs = PF.build_jobs(args, panel, fold_plan, functions, universe, cells, prf, p2)
    del functions
    sel = [j for j in jobs if want is None or j["genome_id"] in want]
    missing = sorted(set(want or []) - {j["genome_id"] for j in sel})
    if missing:
        raise SystemExit("genomes not in the 200-genome selection: {}".format(missing))
    arms = None if a.arms == "all" else set(a.arms.split(","))
    for j in sel:
        j["_outdir"] = a.outdir
        j["_variant"] = a.variant
        j["_template_dir"] = a.template_dir
        j["_want_arms"] = set(j["arms"]) if arms is None else (arms | {"truth"})
    print("[egc] {} jobs ({} gated genomes); arms per job {}".format(
        len(sel), len(gated), sorted(sel[0]["_want_arms"]) if sel else []), flush=True)

    import multiprocessing as mp
    os.makedirs(a.outdir, exist_ok=True)
    log = open(os.path.join(a.outdir, "..", "rebuild_results_{}.jsonl".format(a.variant)), "a")
    pool = mp.get_context("spawn").Pool(min(a.workers, len(sel)),
                                        initializer=PF.worker_init, initargs=(PF.MSD_DEFAULT,))
    for res in pool.imap_unordered(run_one, sel):
        log.write(json.dumps(res, sort_keys=True) + "\n")
        log.flush()
        print("[egc] {} done in {}s: {}".format(
            res["genome_id"], res.get("wall_seconds"),
            [(r.get("key", "").split("__", 1)[-1], r.get("sanity_gate_ok"), r.get("error"))
             for r in res["results"]]), flush=True)
    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
