#!/usr/bin/env python
"""R4.3: check the new weighted-objective population.json (FM) against the expected new inputs and the
pre-repair population; write control/r4/wobj_population_check.json. Exit 1 if a hard check fails."""
import hashlib, json, sys
SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference"
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
new_p = f"{SR}/data/weighted_objective/population.json"
old_p = f"{W}/pre_repair/region_validation/shared_reference/data/weighted_objective/population.json"
EXP = {"cells": "d57709e36ffb2b675f2b48d1e68081ec02070a3fc7a1ab84fbbf2076e0b49ce7",
       "panel_ec": "8eb296b1a083d3e80869779913e1193cdf2eab86d5d5c5bbc74f5a71bb56b2c6",
       "panel": "6eafdc5eb464af2f1c582f35d128c6d72235e49a5c390bb4b979c6a8cb2dc310",
       "ec_universe": "4e9da0800802cdd7c78f73491203196bade25b2ddffa473adc62e166bee95ca2",
       "fold_plan": "fdc8734a7302dac89e97b2d350f5eae384e9a22b78b66b38d852a10362bdb8f3"}
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()
new, old = json.load(open(new_p)), json.load(open(old_p))
nm, om = new["meta"], old["meta"]
out = {"population_json": {"path": new_p, "sha256": sha(new_p)}, "checks": {}, "vs_pre_repair": {}}
c = out["checks"]
for k, v in EXP.items():
    c[f"input_{k}_sha256_expected"] = nm["spec"]["inputs"][k]["sha256"] == v
c["sampled_by_stratum_100_100_100_60"] = nm["sampled_by_stratum"] == {"n1": 60, "n2": 100, "n3": 100, "n4+": 100}
c["n_jobs_360"] = len(new["jobs"]) == 360
c["candidate_checks_all_zero"] = all(v == 0 for v in nm["candidate_checks"].values())
c["spec_settings_equal_old"] = {k: v for k, v in nm["spec"].items() if k != "inputs"} == \
    {k: v for k, v in om["spec"].items() if k != "inputs"}
ok = all(c.values())
ok_keys = {j["key"] for j in new["jobs"]}; ol_keys = {j["key"] for j in old["jobs"]}
oj = {j["key"]: j for j in old["jobs"]}
same_job = sum(1 for j in new["jobs"] if j["key"] in oj and j == oj[j["key"]])
out["vs_pre_repair"] = {
    "n_candidates": [om["n_candidates"], nm["n_candidates"]],
    "candidate_cells_by_stratum": [om["candidate_cells_by_stratum"], nm["candidate_cells_by_stratum"]],
    "eligible_source_genomes_by_stratum_excl_twin": [om["eligible_source_genomes_by_stratum_excl_twin"],
                                                      nm["eligible_source_genomes_by_stratum_excl_twin"]],
    "assertions": [om["assertions"], nm["assertions"]],
    "keys_in_both": len(ok_keys & ol_keys), "keys_only_new": len(ok_keys - ol_keys),
    "keys_only_old": len(ol_keys - ok_keys), "jobs_identical_dict": same_job,
    "job_fields": sorted(new["jobs"][0].keys()),
}
out["all_hard_checks_pass"] = ok
json.dump(out, open(f"{W}/control/r4/wobj_population_check.json", "w"), indent=1)
print(json.dumps(out, indent=1)[:3500])
sys.exit(0 if ok else 1)
