#!/usr/bin/env python
"""Scope of the C12 leak: how often does a RAW genome-scale draft make ATP from nothing, under
the original templates (main 2294299) and under the five-reaction fix (21a8084)?

EGC_C12_REPORT.md found that 190 of the 218 raw drafts of the 20 gated genomes leak, that
MSATPCorrection's expansion step closes those loops, and that the final MSGapfill step
re-opens one in 44 models. This scan asks how general the draft-level loops are, over every
genome and arm the phenotype passes use, with nothing gapfilled:

  for each (genome, arm) of phenotype_fold.py's own job list (same selection, seed and EC
  sets) and each --templates directory:
    build the base model exactly as gapfill_block does (genome_from_ec ->
    MSBuilder.build_base_model(allow_all_non_grp_reactions=False, annotate_with_rast=False)
    -> add_atpm), then phenotype_arm.sanity_gate: ATPM_c0 and bio1 with every exchange
    closed for uptake and water and protons free.

    python egc_scope_scan.py <phenotype_fold.py input flags> --n-genomes 200 \
        --templates DIR [--templates DIR ...] --out data/egc_analysis/scope_scan_200.jsonl
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import phenotype_fold as PF  # noqa: E402

_W = {}


def pinned_libraries(args, out):
    """Print where modelseedpy was imported from, stop unless it is a clean checkout of
    --expect-modelseedpy when that is given, and record every library in OUT.meta.json."""
    libs = PF.library_provenance()
    msp = libs.get("modelseedpy") or {}
    print("[scope] modelseedpy %s at %s (tracked files changed: %s)" % (
        msp.get("path"), msp.get("git_commit"), msp.get("git_dirty")), flush=True)
    if args.expect_modelseedpy and not (
            (msp.get("git_commit") or "").startswith(args.expect_modelseedpy)
            and msp.get("git_dirty") is False):
        raise SystemExit("[scope] modelseedpy is not a clean checkout of %s"
                         % args.expect_modelseedpy)
    with open(out + ".meta.json", "w") as fh:
        json.dump({"libraries": libs, "argv": sys.argv[1:]}, fh, indent=1, sort_keys=True)
    return libs


def init(msd, template_dirs):
    import logging
    for name in ("modelseedpy", "cobra", "optlang"):
        logging.getLogger(name).setLevel(logging.ERROR)
    m = PF.pa()
    m.preload_biochemistry(msd)
    try:
        import gurobipy
        gurobipy.setParam("Threads", 1)
        gurobipy.setParam("OutputFlag", 0)
    except Exception:
        pass
    _W["tmpl"] = {}
    for d in template_dirs:              # every template version, loaded once per worker
        PF.use_templates(d)
        _W["tmpl"][d] = {k: m.load_template(k) for k in ("gramneg", "grampos", "archaea",
                                                          "core")}


def scan(job):
    from modelseedpy.core.msbuilder import MSBuilder
    m = PF.pa()
    rows = []
    for d, T in _W["tmpl"].items():
        tmpl, role_ecs = T[job["template"]]
        core_t = T["core"][0]
        for arm, ecs in sorted(job["arms"].items()):
            ec = PF.norm_ec(ecs)
            t0 = time.time()
            row = {"src_genome_id": job["genome_id"], "arm": arm, "template": job["template"],
                   "templates_dir": d, "n_ec": len(ec)}
            try:
                g = m.genome_from_ec(sorted(ec), role_ecs, "scan")
                b = MSBuilder(g, tmpl, ontology_term="RAST", template_core=core_t)
                base = b.build_base_model("scan", "0", allow_all_non_grp_reactions=False,
                                          annotate_with_rast=False)
                MSBuilder.add_atpm(base)
                gate = m.sanity_gate(base)
                row.update({"n_rxn": len(base.reactions),
                            "empty_atp": round(gate["empty_atp"], 6),
                            "empty_biomass": round(gate["empty_biomass"], 6),
                            "leak": not gate["pass"], "status": "ok"})
            except Exception as exc:                     # noqa: BLE001
                row.update({"status": "error", "error": "%s: %s" % (type(exc).__name__, exc)})
            row["wall_seconds"] = round(time.time() - t0, 2)
            rows.append(row)
    return rows


def main():
    argv = sys.argv[1:]
    dirs, out, rest, i = [], None, [], 0
    while i < len(argv):
        if argv[i] == "--templates":
            dirs.append(os.path.abspath(argv[i + 1]).rstrip("/") + "/"); i += 2
        elif argv[i] == "--out":
            out = argv[i + 1]; i += 2
        else:
            rest.append(argv[i]); i += 1
    if not dirs or not out:
        raise SystemExit("need --templates DIR (one or more) and --out FILE")
    args = PF.parse_args(rest + ["--outdir", os.path.dirname(os.path.abspath(out))])
    pinned_libraries(args, out)
    panel, fold_plan, functions = (PF.jload(args.panel), PF.jload(args.fold_plan),
                                   PF.jload(args.functions))
    ec_path = args.ec_universe or os.path.join(os.path.dirname(os.path.abspath(args.panel)),
                                               "panel_ec_universe.json")
    universe = PF.ec_universe(ec_path, panel, functions)
    cells = PF.load_cells(args.cells, args.density)
    prf, _ = PF.load_prfba(args.prfba, args.select_mode, args.density)
    p2, _ = PF.load_picrust2(args.picrust2, args.density)
    jobs = PF.build_jobs(args, panel, fold_plan, functions, universe, cells, prf, p2)
    del functions
    jobs = [{"genome_id": j["genome_id"], "template": j["template"],
             "arms": {a: v for a, v in j["arms"].items()}} for j in jobs]
    print("[scope] %d genomes x %d template dirs" % (len(jobs), len(dirs)), flush=True)
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    n = 0
    with open(out, "w") as fh, ctx.Pool(min(args.workers, len(jobs)), initializer=init,
                                        initargs=(args.msd, dirs)) as pool:
        for rows in pool.imap_unordered(scan, jobs):
            for r in rows:
                fh.write(json.dumps(r, sort_keys=True) + "\n")
            n += 1
            if n % 25 == 0 or n == len(jobs):
                print("[scope] %d/%d genomes" % (n, len(jobs)), flush=True)
    print("[scope] wrote", out, flush=True)


if __name__ == "__main__":
    main()
