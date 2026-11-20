#!/usr/bin/env python
"""R3.4 checks over the 45 re-selected runs (old = pre_repair copies, new = canonical):
- log: EXIT=0, '[skip] reusing', no '[align]', C6/C8 clean (fold) or C2 source recovery (include_self),
  inputs_identical=1;
- meta.assertions: reps_outside_reference == 0, hits_outside_reference == 0, forced_invariant (forced runs),
  n_reference_genomes_missing_function_set/pgfam_cache == 0; meta equal old->new minus run fields and assertions;
- subsets: old vs new differ ONLY in the changed eligible genomes (= fold_pgfam_staleness), and those all differ;
- predictions (W8): every changed cell has a changed genome among its old or new reps (cmp_prfba_pred logic).
"""
import json, os, re, sys, subprocess, collections
import orjson
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference"
B = f"{W}/pre_repair/region_validation/shared_reference"
C = f"{W}/control/r3"
ch = set(map(str, json.load(open(f"{W}/expected/changed_sets.json"))["changed_eligible"]))
stale = json.load(open(f"{W}/expected/rescore_work/fold_pgfam_staleness.json"))["per_subset"]
L = lambda p: orjson.loads(open(p, "rb").read())
runs = [l.rstrip("\n").split("\t")[:2] for l in open(f"{C}/r34_runs.tsv")]
out = {"runs": {}, "problems": []}
os.makedirs(f"{C}/r34_cmp", exist_ok=True)
for root, name in runs:
    d, b = f"{SR}/{root}/prfba/{name}", f"{B}/{root}/prfba/{name}"
    r = {}
    log = open(f"{W}/logs/r3/prfba_{root}_{name}.log").read()
    r["exit0"] = bool(re.search(r"^EXIT=0 ", log, re.M))
    r["hits_reused"] = "[skip] reusing" in log and not re.search(r"^\[align\]", log, re.M)
    r["inputs_identical"] = "inputs_identical=1" in log
    r["canary_line"] = ("[canary] C6/C8 clean" in log)
    m_new = L(f"{d}/predictions.json"); m_old = L(f"{b}/predictions.json")
    a = m_new["meta"]["assertions"]
    r["reps_outside_reference"] = a.get("reps_outside_reference")
    r["hits_outside_reference"] = (a.get("c6") or {}).get("hits_outside_reference")
    r["missing_fn"] = a.get("n_reference_genomes_missing_function_set")
    r["missing_pgfam"] = a.get("n_reference_genomes_missing_pgfam_cache")
    r["forced_invariant"] = a.get("forced_invariant")
    r["source_recovery"] = a.get("source_recovery")
    run_f = {"generated", "git_commit", "wall_seconds", "align_seconds", "assertions"}
    mo = {k: v for k, v in m_old["meta"].items() if k not in run_f}; mn = {k: v for k, v in m_new["meta"].items() if k not in run_f}
    r["meta_equal_minus_run_and_assertions"] = mo == mn
    ao, an = m_old["meta"]["assertions"], a
    r["assertion_diffs"] = {k: [ao.get(k), an.get(k)] for k in set(ao) | set(an) if ao.get(k) != an.get(k)}
    del m_new, m_old
    sub = {}
    for f in ("_cache_ec_subset.json", "_cache_pgfam_subset.json"):
        so, sn = L(f"{b}/{f}"), L(f"{d}/{f}")
        diff = {g for g in set(so) | set(sn) if so.get(g) != sn.get(g)}
        sub[f] = {"n": len(sn), "same_keys": set(so) == set(sn), "n_diff": len(diff),
                  "diff_not_changed_genome": sorted(diff - ch)[:10], "changed_in_subset": len(set(sn) & ch)}
    key = f"{root}/prfba/{name}/_cache_pgfam_subset.json"
    sub["expected_n_changed_pgfam"] = (stale.get(key) or {}).get("n_changed")
    r["subsets"] = sub
    cp = subprocess.run([sys.executable, f"{W}/tools/cmp_prfba_pred.py", f"{b}/predictions.json", f"{d}/predictions.json",
                         "--out", f"{C}/r34_cmp/{root}__{name}.json", "--changed-genomes", f"{W}/expected/changed_sets.json"],
                        capture_output=True, text=True)
    rep = json.load(open(f"{C}/r34_cmp/{root}__{name}.json"))
    r["pred_changed_cells"] = {arm: v["n_changed"] for arm, v in rep["arms"].items()}
    r["pred_unattributed"] = rep.get("n_unattributed_total")
    r["pred_ec_gained"] = sum(v["ec_gained_total"] for v in rep["arms"].values())
    r["pred_ec_lost"] = sum(v["ec_lost_total"] for v in rep["arms"].values())
    is_fold = not name.startswith("incself")
    ok = (r["exit0"] and r["hits_reused"] and r["inputs_identical"] and (r["canary_line"] or not is_fold)
          and r["reps_outside_reference"] == 0 and not r["hits_outside_reference"] and r["missing_fn"] == 0
          and r["missing_pgfam"] == 0 and r["meta_equal_minus_run_and_assertions"]
          and all(s["same_keys"] and not s["diff_not_changed_genome"] for s in (sub["_cache_ec_subset.json"], sub["_cache_pgfam_subset.json"]))
          and sub["_cache_pgfam_subset.json"]["n_diff"] == sub["expected_n_changed_pgfam"]
          and r["pred_unattributed"] == 0)
    if r["forced_invariant"] is not None:
        # a broken contract dies (exit 3); the emitted-verbatim line proves the guard ran
        ok = ok and "[forced-answer] contract:" in log and r["forced_invariant"].get("checked", 0) > 0
    r["ok"] = ok
    if not ok:
        out["problems"].append(f"{root}/{name}")
    out["runs"][f"{root}/{name}"] = r
    print(f"{root}/{name}: ok={ok} cells_changed={sum(r['pred_changed_cells'].values())} unattr={r['pred_unattributed']} "
          f"pgfam_diff={sub['_cache_pgfam_subset.json']['n_diff']}/{sub['expected_n_changed_pgfam']} "
          f"ec_diff={sub['_cache_ec_subset.json']['n_diff']} assert_diffs={sorted(r['assertion_diffs'])}", flush=True)
out["n_runs"] = len(out["runs"]); out["n_ok"] = sum(r["ok"] for r in out["runs"].values())
json.dump(out, open(f"{C}/r34_checks.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ("n_runs", "n_ok", "problems")}))
sys.exit(0 if not out["problems"] else 1)
