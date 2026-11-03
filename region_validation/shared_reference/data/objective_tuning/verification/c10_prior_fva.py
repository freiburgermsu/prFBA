"""Cross-checks against the prior run (weighted-objective test) and the finished FVA pass."""
import json
import os

import numpy as np

import vlib as L

rows = L.load_rows()
prior = L.load_prior()
pk = {r["key"]: r for r in prior}
WA, WB, _ = L.weights()
POOLS = {"A": WA, "B": WB}
out = {}

# ---- 1. shipped - two_tier on ALL 360 Stage O rows vs the prior study's pooled -0.006054
c = L.paired(rows, lambda r: L.v(r, L.SHIPPED, "taf_face_mid"),
             lambda r: L.v(r, "two_tier/lin/G0/X100", "taf_face_mid"), pools=POOLS)
out["stageO_all360_shipped_minus_tt_X100"] = c["pooled"]


# prior rows, real - two_tier, all and by half
def pv(r, arm, met):
    a = r["arms"].get(arm) or {}
    x = a.get(met)
    return None if x is None else float(x)


for tag, rs in (("all", prior), ("tune", [r for r in prior if L.half_of(r["group_id"]) == "tune"]),
                ("confirm", [r for r in prior if L.half_of(r["group_id"]) == "confirm"])):
    rs = [r for r in rs if r["status"] == "ok"]
    c = L.paired(rs, lambda r: pv(r, "real", "taf_face_mid"), lambda r: pv(r, "two_tier", "taf_face_mid"),
                 pools=POOLS)
    out[f"prior_{tag}_real_minus_two_tier"] = c["pooled"]
    if all("real_ex1" in r["arms"] for r in rs):
        c1 = L.paired(rs, lambda r: pv(r, "real_ex1", "taf_face_mid"),
                      lambda r: pv(r, "two_tier_ex1", "taf_face_mid"), pools=POOLS)
        out[f"prior_{tag}_real_ex1_minus_two_tier_ex1"] = c1["pooled"]
c = L.paired(rows, lambda r: L.v(r, "kn/lin/G0/X1", "taf_face_mid"),
             lambda r: L.v(r, "two_tier/lin/G0/X1", "taf_face_mid"), pools=POOLS)
out["stageO_all360_kn_X1_minus_tt_X1"] = c["pooled"]
out["prior_arms_available"] = sorted(prior[0]["arms"].keys())
# do the prior's ex1 arms equal Stage O's kn/two_tier X1 arms on the same rows?
same = {"real_ex1==kn/lin/G0/X1 (face mid, |d|<=1e-9)": 0, "two_tier_ex1==two_tier/lin/G0/X1": 0, "n": 0}
for r in rows:
    p = pk[r["key"]]
    if "real_ex1" not in p["arms"]:
        continue
    same["n"] += 1
    same["real_ex1==kn/lin/G0/X1 (face mid, |d|<=1e-9)"] += abs(p["arms"]["real_ex1"]["taf_face_mid"] - r["arms"]["kn/lin/G0/X1"]["taf_face_mid"]) <= 1e-9
    same["two_tier_ex1==two_tier/lin/G0/X1"] += abs(p["arms"]["two_tier_ex1"]["taf_face_mid"] - r["arms"]["two_tier/lin/G0/X1"]["taf_face_mid"]) <= 1e-9
out["prior_ex1_vs_stageO_X1"] = same

# ---- 2. g_meas from the prior run's tune-half rows
num = den = 0
numc = denc = 0
for r in prior:
    if r["status"] != "ok":
        continue
    m = r["model"] or {}
    if L.half_of(r["group_id"]) == "tune":
        num += m.get("n_gapfilled_supported", 0)
        den += m.get("n_gapfilled", 0)
    else:
        numc += m.get("n_gapfilled_supported", 0)
        denc += m.get("n_gapfilled", 0)
out["g_meas_recomputed"] = {"tune": [num, den, num / den], "confirm": [numc, denc, numc / denc],
                            "all": [num + numc, den + denc, (num + numc) / (den + denc)]}
meta = json.load(open(os.path.join(L.SR, "data/objective_tuning/meta.json")))
out["g_meas_meta"] = [meta["g_meas_info"]["sum_n_gapfilled_supported"], meta["g_meas_info"]["sum_n_gapfilled"],
                      meta["g_meas"]]
# any arm value that depends on g_meas differs between tune/confirm g? (irrelevant; G not decisive)

# ---- 3. FVA pass: same primary numbers?  and the degeneracy contrasts
fva = {}
with open(os.path.join(L.SR, "data/objective_tuning_fva/per_amplicon.jsonl")) as fh:
    for line in fh:
        r = json.loads(line)
        fva[r["key"]] = r
conf = {r["key"]: r for r in rows if r["half"] == "confirm"}
FARMS = ("est/lin/G0/X1", "kn/lin/G0/X100", "two_tier/lin/G0/X1", "est_shuffled/lin/G0/X1")
eq = {"n_rows": len(fva), "keys_equal_confirm_half": set(fva) == set(conf)}
mism = []
for k, r in fva.items():
    for a in FARMS:
        x, y = r["arms"][a], conf[k]["arms"][a]
        for met in ("objective", "taf_face_min", "taf_face_max", "taf_face_mid", "flux_absent",
                    "flux_gene_backed", "flux_exchange"):
            if x.get(met) != y.get(met):
                mism.append((k, a, met, x.get(met), y.get(met)))
eq["n_metric_mismatches_vs_main_run"] = len(mism)
eq["mismatch_examples"] = mism[:10]
eq["mismatch_keys"] = sorted({m[0] for m in mism})
fv_prior_fail = [k for k, r in fva.items() if not (r.get("prior_match") or {}).get("ok")]
eq["fva_prior_match_fail"] = fv_prior_fail
out["fva_vs_main"] = eq


def fv(r, arm, tag):
    a = r["arms"].get(arm) or {}
    f = (a.get("faces") or {}).get(tag) or {}
    x = f.get("n_var_1e-3")
    return None if x is None else float(x)


frows = [dict(fva[k], stratum=fva[k]["stratum"], group_id=fva[k]["group_id"]) for k in sorted(fva)]
deg = {}
for name, (a, b) in {"selected-shipped": ("est/lin/G0/X1", "kn/lin/G0/X100"),
                     "selected-two_tier": ("est/lin/G0/X1", "two_tier/lin/G0/X1"),
                     "selected-shuffled": ("est/lin/G0/X1", "est_shuffled/lin/G0/X1")}.items():
    for tag in ("opt", "opt7"):
        c = L.paired(frows, lambda r: fv(r, a, tag), lambda r: fv(r, b, tag), pools=POOLS)
        deg[f"{name}|{tag}"] = {"A": c["pooled"]["A"], "B": c["pooled"]["B"],
                                "per": {s: {"n": p["n"], "mean": p["mean"], "ci": p["ci"]}
                                        for s, p in c["per"].items()}}
# arm levels (pooled A) of the variable-reaction count
lev = {}
for a in FARMS:
    by = L.by_stratum(frows)
    lev[a] = sum(WA[s] * float(np.mean([fv(r, a, "opt") for r in by[s]])) for s in L.MULTI)
out["fva_degeneracy"] = deg
out["fva_levels_A_opt_n_var"] = lev
out["fva_n_gene_backed_scan_variable_mean"] = float(np.mean([r["model"].get("n_gene_backed_scan_variable", np.nan) for r in frows]))
json.dump(out, open(os.path.join(L.HERE, "out_c10.json"), "w"), indent=1, default=str)
for k, x in out.items():
    if k in ("fva_degeneracy",):
        for kk, vv in x.items():
            print("  deg", kk, "A", L.fmt(vv["A"]), "B", L.fmt(vv["B"]))
        continue
    if isinstance(x, dict) and "A" in x and isinstance(x["A"], dict):
        print(k, "A", L.fmt(x["A"]), "B", L.fmt(x["B"]))
    else:
        print(k, json.dumps(x, default=str)[:900])
