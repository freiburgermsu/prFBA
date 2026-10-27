#!/usr/bin/env python
"""Which reactions form the ATP-from-nothing loops of RAW genome-scale drafts?

egc_scope_scan.py found that most raw drafts (before MSATPCorrection or gapfilling) make
ATP with every exchange closed for uptake and water and protons free. This script takes the
drafts that scan flagged and, for each, enumerates up to --k minimum-cardinality reaction
supports of ATPM_c0 under that empty medium. It uses egc_analyze.min_support: a MILP with
binary flux indicators and no-good cuts, so the i-th support is the smallest flux-carrying
network that differs from the first i-1. The draft is rebuilt exactly as egc_scope_scan.py
builds it.

    python egc_draft_loops.py <phenotype_fold.py input flags> --n-genomes 200 \
        --scan data/egc_analysis/scope_scan_200.jsonl --arms truth --templates DIR \
        --k 5 --out data/egc_analysis/draft_loops_200.jsonl
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import phenotype_fold as PF  # noqa: E402
import egc_analyze as EA  # noqa: E402
import egc_scope_scan as SS  # noqa: E402


def loops(job):
    from modelseedpy.core.msbuilder import MSBuilder
    m = PF.pa()
    out = []
    for d, T in SS._W["tmpl"].items():
        tmpl, role_ecs = T[job["template"]]
        core_t = T["core"][0]
        for arm in job["want"]:
            ec = PF.norm_ec(job["arms"][arm])
            t0 = time.time()
            row = {"src_genome_id": job["genome_id"], "arm": arm, "template": job["template"],
                   "templates_dir": d}
            try:
                g = m.genome_from_ec(sorted(ec), role_ecs, "loop")
                b = MSBuilder(g, tmpl, ontology_term="RAST", template_core=core_t)
                base = b.build_base_model("loop", "0", allow_all_non_grp_reactions=False,
                                          annotate_with_rast=False)
                MSBuilder.add_atpm(base)
                m.ensure_exchanges(base, {"cpd00001", "cpd00067"})
                base.solver = "gurobi"
                EA.empty_medium(base)
                mx = EA.maxflux(base, "ATPM_c0")
                row["empty_atp"] = round(mx, 6)
                sups, cuts = [], []
                if mx > EA.THRESH:
                    for _ in range(job["k"]):
                        sup, st = EA.min_support(base, "ATPM_c0", min(1.0, mx / 2), cuts=cuts,
                                                 timeout=120)
                        if sup is None:
                            break
                        sups.append(sup["reactions"])
                        cuts.append(sup["reactions"])
                row.update({"supports": sups, "status": "ok"})
            except Exception as exc:                     # noqa: BLE001
                row.update({"status": "error", "error": "%s: %s" % (type(exc).__name__, exc)})
            row["secs"] = round(time.time() - t0, 1)
            out.append(row)
    return out


def main():
    argv = sys.argv[1:]
    dirs, out, scan, arms, k, rest, i = [], None, None, ["truth"], 5, [], 0
    while i < len(argv):
        a = argv[i]
        if a == "--templates":
            dirs.append(os.path.abspath(argv[i + 1]).rstrip("/") + "/"); i += 2
        elif a == "--out":
            out = argv[i + 1]; i += 2
        elif a == "--scan":
            scan = argv[i + 1]; i += 2
        elif a == "--arms":
            arms = argv[i + 1].split(","); i += 2
        elif a == "--k":
            k = int(argv[i + 1]); i += 2
        else:
            rest.append(a); i += 1
    leaking = None
    if scan:
        leaking = set()
        for line in open(scan):
            r = json.loads(line)
            if r.get("leak") and os.path.abspath(r["templates_dir"]).rstrip("/") + "/" in dirs:
                leaking.add((r["src_genome_id"], r["arm"]))
    args = PF.parse_args(rest + ["--outdir", os.path.dirname(os.path.abspath(out))])
    SS.pinned_libraries(args, out)
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
    todo = []
    for j in jobs:
        want = [a for a in arms if a in j["arms"]
                and (leaking is None or (j["genome_id"], a) in leaking)]
        if want:
            todo.append({"genome_id": j["genome_id"], "template": j["template"],
                         "arms": {a: j["arms"][a] for a in want}, "want": want, "k": k})
    print("[loops] %d genomes, %d drafts" % (len(todo), sum(len(t["want"]) for t in todo)),
          flush=True)
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    n = 0
    with open(out, "w") as fh, ctx.Pool(min(args.workers, len(todo)), initializer=SS.init,
                                        initargs=(args.msd, dirs)) as pool:
        for rows in pool.imap_unordered(loops, todo):
            for r in rows:
                fh.write(json.dumps(r, sort_keys=True) + "\n")
            n += 1
            if n % 20 == 0 or n == len(todo):
                print("[loops] %d/%d genomes" % (n, len(todo)), flush=True)
    print("[loops] wrote", out, flush=True)


if __name__ == "__main__":
    main()
