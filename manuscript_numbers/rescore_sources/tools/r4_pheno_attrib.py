#!/usr/bin/env python
"""R4.1 control: attribute every changed phenotype row (A-C13 configuration, old -> new) to a changed input.

Rebuilds phenotype_fold's own jobs (build_jobs: genome sample + every arm's input EC set + donors) twice, from the
PRE-REPAIR inputs and from the NEW inputs, then compares the old and new per_genome.jsonl row by row on
(genome, arm, pipeline), ignoring wall_seconds.

A row may differ only if its own input set changed or the genome's truth set changed (every row's vs_truth block and
the gapfill truth context read the truth model):
  truth, ladder_dropout_*, ladder_false_positive_*                      -> arms.truth
  prfba_exact_tie, ladder_size_matched_null_prfba_exact_tie             -> arms.prfba_exact_tie
  picrust2, ladder_size_matched_null_picrust2                           -> arms.picrust2
  base_nn1, ladder_donor_nn1, ladder_size_matched_null_base_nn1         -> donors.nn1
  base_family_donor, ladder_donor_family, ladder_size_matched_null_...  -> donors.family
  ladder_donor_order                                                    -> donors.order
  base_constant, ladder_size_matched_null_base_constant                 -> arms.base_constant
Any other changed row is UNATTRIBUTED (exit 1). Also checks that each rebuilt sample equals the genome list the
corresponding summary's meta.binding records (proves the old-input reconstruction).

usage: r4_pheno_attrib.py --pass atp|gf --out JSON
"""
import argparse
import gc
import hashlib
import json
import os
import sys
import tempfile
from collections import Counter, defaultdict

R = "/home/freiburger/Documents/prFBA"
SR = f"{R}/region_validation/shared_reference"
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
PRS = f"{W}/pre_repair/region_validation/shared_reference"
PRP = f"{W}/pre_repair/region_validation/picrust2_comparison"
sys.path.insert(0, f"{SR}/scripts")
import phenotype_fold as PF  # noqa: E402

PASSES = {"atp": ("1000", "atp,structure"), "gf": ("200", "atp,growth,essentiality,structure")}


def old_picrust2_root():
    """data/picrust2 as it was before R3.3: the moved-aside run dirs where R3.3 re-ran, the untouched ones otherwise."""
    d = tempfile.mkdtemp(prefix="r4_old_p2_", dir=f"{W}/control/r4")
    n_old = 0
    for run in sorted(os.listdir(f"{SR}/data/picrust2")):
        src_old = f"{PRS}/data/picrust2/{run}/predictions.json"
        src = src_old if os.path.exists(src_old) else f"{SR}/data/picrust2/{run}/predictions.json"
        if not os.path.exists(src):
            continue
        n_old += src is src_old
        os.makedirs(f"{d}/{run}")
        os.symlink(src, f"{d}/{run}/predictions.json")
    return d, n_old


def h(xs):
    return hashlib.sha256(json.dumps(xs, sort_keys=True).encode()).hexdigest()[:16]


def jobs_for(pass_, cells, functions, prfba, picrust2):
    n, ro = PASSES[pass_]
    args = PF.parse_args([
        "--cells", cells, "--fold-plan", f"{SR}/data/fold_plan.json", "--panel", f"{SR}/data/panel_eligible.json",
        "--prfba", prfba, "--picrust2", picrust2, "--functions", functions,
        "--outdir", f"{W}/control/r4/_unused_{pass_}", "--n-genomes", n, "--readouts", ro])
    panel, fold_plan = PF.jload(args.panel), PF.jload(args.fold_plan)
    fn = PF.jload(args.functions)
    universe = PF.ec_universe(os.path.join(os.path.dirname(args.panel), "panel_ec_universe.json"), panel, fn)
    cl = PF.load_cells(args.cells, args.density)
    prf, _ = PF.load_prfba(args.prfba, args.select_mode, args.density)
    p2, _ = PF.load_picrust2(args.picrust2, args.density)
    jobs = PF.build_jobs(args, panel, fold_plan, fn, universe, cl, prf, p2)
    out = {}
    for j in jobs:
        out[j["genome_id"]] = {"arms": {a: h(v) for a, v in j["arms"].items()},
                               "donors": {a: h(v) for a, v in j["donors"].items()},
                               "md5": j["md5"], "region": j["region"], "fold": j["fold"], "template": j["template"]}
    del fn, cl, prf, p2, panel, fold_plan, jobs
    gc.collect()
    return out


def dep(arm):
    if arm == "truth" or arm.startswith("ladder_dropout") or arm.startswith("ladder_false_positive"):
        return ("arms", "truth")
    for a in ("prfba_exact_tie", "picrust2", "base_constant"):
        if arm == a or arm == f"ladder_size_matched_null_{a}":
            return ("arms", a)
    for tag, a in (("nn1", "base_nn1"), ("family", "base_family_donor")):
        if arm in (a, f"ladder_donor_{tag}", f"ladder_size_matched_null_{a}"):
            return ("donors", tag)
    if arm == "ladder_donor_order":
        return ("donors", "order")
    return None


def rows(path):
    out = {}
    for line in open(path):
        if line.strip():
            r = json.loads(line)
            r.pop("wall_seconds", None)
            out[(r["src_genome_id"], r["arm"], r["pipeline"])] = r
    return out


def diff_fields(a, b, p=""):
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            out += diff_fields(a.get(k), b.get(k), f"{p}/{k}")
        return out
    return [] if a == b else [p]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass", dest="pass_", required=True, choices=sorted(PASSES))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    p = a.pass_
    old_dir, new_dir = f"{PRS}/data/phenotype_{p}_c13", f"{SR}/data/phenotype_{p}_c13"
    p2_old, n_p2_old = old_picrust2_root()
    old_j = jobs_for(p, f"{PRS}/data/cells.csv.gz", f"{PRP}/data/genome_function_sets.json", f"{PRS}/data/prfba", p2_old)
    new_j = jobs_for(p, f"{SR}/data/cells.csv.gz", f"{R}/region_validation/picrust2_comparison/data/genome_function_sets.json",
                     f"{SR}/data/prfba", f"{SR}/data/picrust2")
    old_s = json.load(open(f"{old_dir}/phenotype_summary.json"))["meta"]["binding"]["genomes"]
    new_s = json.load(open(f"{new_dir}/phenotype_summary.json"))["meta"]["binding"]["genomes"]
    res = {"pass": p, "old_picrust2_runs_from_pre_repair": n_p2_old,
           "rebuilt_old_sample_equals_old_binding": sorted(old_j) == sorted(old_s),
           "rebuilt_new_sample_equals_new_binding": sorted(new_j) == sorted(new_s),
           "n_old": len(old_s), "n_new": len(new_s), "n_both": len(set(old_s) & set(new_s)),
           "only_old": sorted(set(old_s) - set(new_s)), "only_new": sorted(set(new_s) - set(old_s))}
    both = set(old_j) & set(new_j)
    changed_inputs = {}
    for g in sorted(both):
        o, n = old_j[g], new_j[g]
        ch = {f"arms.{k}" for k in set(o["arms"]) | set(n["arms"]) if o["arms"].get(k) != n["arms"].get(k)}
        ch |= {f"donors.{k}" for k in set(o["donors"]) | set(n["donors"]) if o["donors"].get(k) != n["donors"].get(k)}
        ch |= {k for k in ("md5", "region", "fold", "template") if o[k] != n[k]}
        if ch:
            changed_inputs[g] = sorted(ch)
    res["genomes_in_both_with_changed_inputs"] = changed_inputs
    res["input_change_counts"] = dict(Counter(x for v in changed_inputs.values() for x in v))
    ro, rn = rows(f"{old_dir}/per_genome.jsonl"), rows(f"{new_dir}/per_genome.jsonl")
    reasons, unattr, fields = Counter(), [], Counter()
    n_cmp = n_same = 0
    for key in sorted(set(ro) & set(rn)):
        g, arm, pipe = key
        if g not in both:
            continue
        n_cmp += 1
        if ro[key] == rn[key]:
            n_same += 1
            continue
        chg = set(changed_inputs.get(g, []))
        d = dep(arm)
        if {"md5", "region", "fold", "template"} & chg:
            reasons["amplicon_or_template_changed"] += 1
        elif d and f"{d[0]}.{d[1]}" in chg:
            reasons[f"own_input_changed:{d[0]}.{d[1]}"] += 1
        elif "arms.truth" in chg:
            reasons["truth_changed"] += 1
        else:
            unattr.append({"key": key, "fields": diff_fields(ro[key], rn[key])[:12]})
        for f in diff_fields(ro[key], rn[key]):
            fields[f.split("/")[1]] += 1
    res["rows_only_old_sample_keys"] = len([k for k in ro if k[0] not in both])
    res["rows_only_new_sample_keys"] = len([k for k in rn if k[0] not in both])
    res["rows_compared_genomes_in_both"] = n_cmp
    res["rows_identical"] = n_same
    res["rows_changed"] = n_cmp - n_same
    res["changed_by_reason"] = dict(reasons)
    res["changed_top_level_fields"] = dict(fields)
    res["unattributed"] = len(unattr)
    res["unattributed_examples"] = unattr[:20]
    res["genomes_unchanged_inputs_all_rows_identical"] = all(
        ro[k] == rn[k] for k in set(ro) & set(rn) if k[0] in both and k[0] not in changed_inputs)
    res["pass_ok"] = (res["unattributed"] == 0 and res["rebuilt_old_sample_equals_old_binding"]
                      and res["rebuilt_new_sample_equals_new_binding"])
    json.dump(res, open(a.out, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k not in ("genomes_in_both_with_changed_inputs", "only_old",
                                                                  "only_new", "unattributed_examples")}, indent=1))
    sys.exit(0 if res["pass_ok"] else 1)


if __name__ == "__main__":
    main()
