#!/usr/bin/env python
"""Does ModelSEEDpy's own MSBuilder.build_and_gapfill reproduce the harness's gapfilled
models once it applies the gapfiller's reversals?

For each (genome, arm) given, this builds the gapfilled model twice in one process, with
the same EC set, template and complete medium:
  package  phenotype_arm.build_gapfilled, which calls MSBuilder.build_and_gapfill
  harness  phenotype_fold.gapfill_block (as_shipped), the builder inlined, with the
           reversals applied as SPEC amendment A-C12 does (--integrate-reversals)
and reports the gapfill solution, the reversed reactions and their final bounds, growth on
the complete medium and the C12 sanity gate. Run it with the ModelSEEDpy to test first on
PYTHONPATH (the fix branch reproduces the harness; an unfixed commit drops the reversals).

    python validate_builder_reversals.py <phenotype_fold.py input flags> --templates DIR \\
        --n-genomes 200 --targets targets.json --out result.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import phenotype_fold as PF  # noqa: E402


def main():
    argv = sys.argv[1:]
    targets = out = None
    rest, i = [], 0
    while i < len(argv):
        if argv[i] == "--targets":
            targets = json.load(open(argv[i + 1])); i += 2
        elif argv[i] == "--out":
            out = argv[i + 1]; i += 2
        else:
            rest.append(argv[i]); i += 1
    args = PF.parse_args(rest + ["--outdir", os.path.dirname(os.path.abspath(out))])
    libs = PF.library_provenance()
    print("[revcheck] modelseedpy", libs["modelseedpy"]["path"], libs["modelseedpy"]["git_commit"],
          "dirty", libs["modelseedpy"]["git_dirty"], flush=True)
    panel, fold_plan, functions = (PF.jload(args.panel), PF.jload(args.fold_plan),
                                   PF.jload(args.functions))
    ec_path = args.ec_universe or os.path.join(os.path.dirname(os.path.abspath(args.panel)),
                                               "panel_ec_universe.json")
    universe = PF.ec_universe(ec_path, panel, functions)
    cells = PF.load_cells(args.cells, args.density)
    prf, _ = PF.load_prfba(args.prfba, args.select_mode, args.density)
    p2, _ = PF.load_picrust2(args.picrust2, args.density)
    jobs = {j["genome_id"]: j for j in
            PF.build_jobs(args, panel, fold_plan, functions, universe, cells, prf, p2)}
    del functions, cells
    PF.worker_init(args.msd, True, args.templates, True)
    m = PF.pa()
    gf_media, gf_cpds = m.complete_media()
    res = {"libraries": libs, "models": []}
    for gid, arm in targets:
        job = jobs[gid]
        ec, tmpl_key = job["arms"][arm], job["template"]
        mid = PF.safe_id(gid, arm, "revcheck")
        pk, gapsol, _atpc, _tests = m.build_gapfilled(sorted(ec), tmpl_key, mid, gf_media)
        hk, block, _e, _c, _g, _ctx, infeasible = PF.gapfill_block(
            m, ec, tmpl_key, mid, gf_media, gf_cpds, "as_shipped", None, False)
        rev = sorted(block["reversed_proposed"])
        entry = {
            "genome": gid, "arm": arm, "template": tmpl_key, "harness_infeasible": infeasible,
            "new_rxn": {"package": sorted(gapsol), "harness": block["new_rxn"]},
            "reversed_proposed": block["reversed_proposed"],
            "reversed_bounds": {
                "package": {r: list(pk.reactions.get_by_id(r).bounds) for r in rev
                            if r in pk.reactions},
                "harness": {r: list(hk.reactions.get_by_id(r).bounds) for r in rev
                            if r in hk.reactions}},
            "growth": {
                "package": round(float(m.growth_vector(
                    pk, {"complete": gf_cpds}, thresh=PF.GROWTH_THRESH).get("complete", 0.0)), 6),
                "harness": block["growth"]},
            "sanity_gate_ok": {"package": bool(m.sanity_gate(pk)["pass"]),
                               "harness": block["sanity_gate_ok"]},
        }
        entry["identical"] = (entry["new_rxn"]["package"] == entry["new_rxn"]["harness"]
                              and entry["reversed_bounds"]["package"]
                              == entry["reversed_bounds"]["harness"]
                              and abs(entry["growth"]["package"] - entry["growth"]["harness"]) < 1e-6
                              and entry["sanity_gate_ok"]["package"]
                              == entry["sanity_gate_ok"]["harness"])
        res["models"].append(entry)
        print("[revcheck]", gid, arm, "identical" if entry["identical"] else "DIFFER",
              "| growth package", entry["growth"]["package"], "harness", entry["growth"]["harness"],
              "| reversed", entry["reversed_bounds"], flush=True)
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)


if __name__ == "__main__":
    main()
