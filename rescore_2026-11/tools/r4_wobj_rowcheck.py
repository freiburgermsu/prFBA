#!/usr/bin/env python
"""R4.3 wobj attribution (W8 analogue): every amplicon whose population job is unchanged (old -> new population.json)
must give a row identical to the committed 10-29 row in every field except timing/memory; only rows of changed jobs
may differ. Writes control/r4/wobj_rowcheck.json; exit 1 on any unattributed difference."""
import json
import sys

SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference/data/weighted_objective"
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
OLD = f"{W}/pre_repair/region_validation/shared_reference/data/weighted_objective"
SKIP = ("seconds", "timing", "maxrss", "wall")


def walk(a, b, p, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if any(s in k for s in SKIP):
                continue
            if k not in a or k not in b:
                out.append(p + "/" + k)
                continue
            walk(a[k], b[k], p + "/" + k, out)
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, "%s[%d]" % (p, i), out)
    elif a != b or type(a) != type(b):
        out.append(p)


def rows(path):
    return {r["key"]: r for r in map(json.loads, filter(str.strip, open(path)))}


pn, po = json.load(open(f"{SR}/population.json")), json.load(open(f"{OLD}/population.json"))
jo = {j["key"]: j for j in po["jobs"]}
changed_jobs = sorted(j["key"] for j in pn["jobs"] if jo.get(j["key"]) != j)
rn, ro = rows(f"{SR}/per_amplicon.jsonl"), rows(f"{OLD}/per_amplicon.jsonl")
res = {"n_rows_new": len(rn), "n_rows_old": len(ro), "same_keys": sorted(rn) == sorted(ro),
       "status_ok_new": sum(1 for r in rn.values() if r.get("status") == "ok"),
       "changed_jobs": changed_jobs, "unchanged_job_rows_identical": 0, "unattributed": [], "changed_job_rows": {}}
for k in sorted(rn):
    d = []
    walk(rn[k], ro.get(k, {}), "", d)
    if k in changed_jobs:
        res["changed_job_rows"][k] = {"n_fields_differ": len(d), "first": d[:10]}
    elif d:
        res["unattributed"].append({"key": k, "n_fields_differ": len(d), "first": d[:10]})
    else:
        res["unchanged_job_rows_identical"] += 1
res["pass"] = res["same_keys"] and not res["unattributed"] and res["status_ok_new"] == len(rn) == 360
json.dump(res, open(f"{W}/control/r4/wobj_rowcheck.json", "w"), indent=1)
print(json.dumps({k: (v if k != "unattributed" else v[:5]) for k, v in res.items()}, indent=1))
sys.exit(0 if res["pass"] else 1)
