"""Phase 5: side-by-side of my numbers (p4_results.json) with the study's scores.json."""
import json
import os

import orjson

import vctx

mine = json.load(open(os.path.join(vctx.OUT, "p4_results.json")))
theirs = orjson.loads(open(os.path.join(vctx.DATA, "pfunc_estimators", "scores.json"), "rb").read())
T = theirs["composites"]["d100"]
MAP = {"universe_brier": "universe_brier", "universe_logloss": "universe_logloss_eps0.01",
       "union_brier": "union_brier", "union_wls_slope": "union_slope_wls",
       "union_auroc": "union_auroc", "macro_dF1_vs_union": "d_macro_F1_vs_union",
       "micro_dF1_vs_union": "d_micro_F1_vs_union"}
rows = []
worst = {"est": 0.0, "ci": 0.0}


def fmt(e, lo=None, hi=None, k=4):
    if e is None:
        return "n/a"
    s = f"{e:.{k}f}"
    if lo is not None:
        s += f" [{lo:.{k}f}, {hi:.{k}f}]"
    return s


out = {}
for st in ("all", "n1", "n2", "n3", "n4+", "multi"):
    M = mine["strata"][st]
    t = T[st]
    for nm in ("E0", "E1", "E2", "E3", "E4(E0)", "E4(E3)"):
        for m, tm in MAP.items():
            me = M["est"][nm].get(m)
            mci = M["ci"][nm].get(m)
            te = t[nm].get(tm, {})
            if te is None or te.get("est") is None or me is None:
                continue
            de = abs(me - te["est"])
            dci = (max(abs(mci[0] - te["lo"]), abs(mci[1] - te["hi"]))
                   if (mci and te.get("lo") is not None) else None)
            out[f"{st}|{nm}|{m}"] = {"mine": [me] + (mci or []),
                                     "theirs": [te["est"], te.get("lo"), te.get("hi")],
                                     "d_est": de, "d_ci": dci}
            if not nm.startswith("E4"):
                worst["est"] = max(worst["est"], de)
                if dci is not None and m != "union_auroc":
                    worst["ci"] = max(worst["ci"], dci)
    for nm in ("E1", "E2", "E3", "E4(E0)", "E4(E3)"):
        for m, tm in MAP.items():
            mc = M["contrasts"].get(f"{nm} - E0", {}).get(m)
            tc = t["contrasts_vs_E0"][nm].get(tm)
            if mc is None or tc is None or tc.get("est") is None:
                continue
            out[f"{st}|{nm}-E0|{m}"] = {"mine": mc, "theirs": [tc["est"], tc.get("lo"), tc.get("hi")],
                                        "d_est": abs(mc[0] - tc["est"])}
    for a, b in (("E3", "E1"), ("E4(E3)", "E3")):
        for m, tm in MAP.items():
            mc = M["contrasts"].get(f"{a} - {b}", {}).get(m)
            tc = t["contrasts_pairwise"].get(f"{a} - {b}", {}).get(tm)
            if mc is None or tc is None or tc.get("est") is None:
                continue
            out[f"{st}|{a}-{b}|{m}"] = {"mine": mc, "theirs": [tc["est"], tc.get("lo"), tc.get("hi")],
                                        "d_est": abs(mc[0] - tc["est"])}
json.dump({"worst_nonE4_abs_diff": worst, "rows": out},
          open(os.path.join(vctx.OUT, "p5_compare.json"), "w"), indent=1)
print("worst non-E4 abs diffs (est, CI bound excl. AUROC):", worst)
keys = [k for k in out if k.split("|")[0] in ("all", "n1", "multi") and
        (k.split("|")[1] in ("E0", "E3", "E3-E0", "E3-E1", "E4(E3)-E3", "E1-E0", "E4(E0)-E0"))]
for k in keys:
    v = out[k]
    mm = v["mine"]
    tt = v["theirs"]
    print(f"{k:45s} mine {fmt(*mm) if len(mm)==3 else fmt(mm[0])}   theirs {fmt(*tt) if tt[1] is not None else fmt(tt[0])}")
