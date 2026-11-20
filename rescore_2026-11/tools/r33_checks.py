#!/usr/bin/env python
"""R3.3 checks: every re-run PICRUSt2 run exited 0 with C3/C3_preflight/C4 pass and C7 reported; placed.tre
(per domain workdir) byte-identical to the pre-repair copy (epa-ng determinism: placement re-ran only because
the stage stamp hashes traits_sha256); stage_place.json differs only in traits_sha256; unaffected runs are
untouched (pre-anything stamp vs final stamp); prediction changes counted per run (queries whose EC set moved).
"""
import json, os, re, sys, glob, hashlib
import orjson
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference"
B = f"{W}/pre_repair/region_validation/shared_reference"
C = f"{W}/control/r3"
L = lambda p: orjson.loads(open(p, "rb").read())
h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
runs = [l.rstrip("\n").split("\t")[:2] for l in open(f"{C}/r33_runs.tsv")]
rerun = {f"{r}/picrust2/{n}" for r, n in runs}
def sha(path):
    d = {}
    for line in open(path):
        x, p = line.rstrip("\n").split("  ", 1); d[p] = x
    return d
before = sha(f"{C}/r33_p2_before_gf.sha256")
after = sha(sys.argv[1]) if len(sys.argv) > 1 else sha(f"{C}/r33_p2_after_c2.sha256")
untouched_changed = sorted(p for p in set(before) | set(after) if before.get(p) != after.get(p)
                           and "/".join(p.split("/")[:3]) not in rerun)
out = {"runs": {}, "problems": [], "untouched_runs_changed_files": untouched_changed}
for root, name in runs:
    d, b = f"{SR}/{root}/picrust2/{name}", f"{B}/{root}/picrust2/{name}"
    r = {}
    log = open(f"{W}/logs/r3/p2_{root}_{name}.log").read()
    r["exit0"] = bool(re.search(r"^EXIT=0 ", log, re.M))
    pn, po = L(f"{d}/predictions.json"), L(f"{b}/predictions.json")
    c = pn["meta"]["controls"]
    st = {}
    for k in ("C3", "C3_preflight", "C4"):
        st[k] = sorted({v["status"] for v in (c.get(k) or {}).values()})
    st["C7"] = (c.get("C7") or {}).get("status")
    r["controls"] = st
    r["controls_ok"] = all(st[k] == ["pass"] for k in ("C3", "C3_preflight", "C4")) and st["C7"] == "reported"
    workdirs = [""] if "nsti" not in name else ["bac/", "arc/"]
    r["placed_identical"] = {}; r["stamp_diff_keys"] = {}; r["placement_reran"] = {}
    for w in workdirs:
        r["placed_identical"][w or "."] = h(f"{d}/{w}placed.tre") == h(f"{b}/{w}placed.tre")
        so, sn = json.load(open(f"{b}/{w}stage_place.json")), json.load(open(f"{d}/{w}stage_place.json"))
        r["stamp_diff_keys"][w or "."] = sorted(k for k in set(so) | set(sn) if so.get(k) != sn.get(k))
        r["placement_reran"][w or "."] = os.path.getmtime(f"{d}/{w}placed.tre") > os.path.getmtime(f"{b}/{w}placed.tre")
    qo, qn = po["pred"], pn["pred"]
    r["n_query"] = len(qn); r["same_queries"] = set(qo) == set(qn)
    r["n_pred_ec_changed"] = sum(1 for q in qn if set(qn[q]["ec"]) != set((qo.get(q) or {}).get("ec", [])))
    r["meta_n_placed"] = [po["meta"].get("n_placed"), pn["meta"].get("n_placed")]
    r["trait_cols"] = [po["meta"].get("n_trait_columns"), pn["meta"].get("n_trait_columns")]
    ok = (r["exit0"] and r["controls_ok"] and all(r["placed_identical"].values()) and r["same_queries"]
          and all(set(v) <= {"traits_sha256"} for v in r["stamp_diff_keys"].values())
          and r["meta_n_placed"][0] == r["meta_n_placed"][1])
    r["ok"] = ok
    if not ok:
        out["problems"].append(f"{root}/{name}")
    out["runs"][f"{root}/{name}"] = r
out["n_runs"] = len(out["runs"]); out["n_ok"] = sum(v["ok"] for v in out["runs"].values())
out["n_runs_placement_reran"] = sum(any(v["placement_reran"].values()) for v in out["runs"].values())
out["n_runs_any_prediction_changed"] = sum(v["n_pred_ec_changed"] > 0 for v in out["runs"].values())
out["n_queries_prediction_changed"] = sum(v["n_pred_ec_changed"] for v in out["runs"].values())
json.dump(out, open(f"{C}/r33_checks.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ("n_runs", "n_ok", "problems", "n_runs_placement_reran", "n_runs_any_prediction_changed",
                                      "n_queries_prediction_changed")} | {"untouched_changed": len(untouched_changed)}))
sys.exit(0 if not out["problems"] and not untouched_changed else 1)
