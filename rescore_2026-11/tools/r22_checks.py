#!/usr/bin/env python
"""R2.2 checks on FM after install: row counts; resolve rates (ok / (ok + abstained)) old vs new per
(regime, completeness) at contamination 0 (main + lowc; tau99) against expected/resolve_impact.json (the
11-13 prediction: old, new, outcome changes); outcome transitions."""
import csv, glob, json, collections, sys
R = "/home/freiburger/Documents/prFBA"; W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
S, P = W + "/start", W + "/pre_repair"


def rows_of(paths):
    out = {}
    for p in paths:
        for r in csv.DictReader(open(p)):
            out[(r["mag"], r["regime"])] = r
    return out


def batches(d):
    return sorted(glob.glob(d + "/b*/mag_cells.csv"))


def outcome(r):
    if r["outcome"] == "ok" and r["precision"] in ("", "None"):
        return "pred_empty"
    return r["outcome"]


arms = {
    "tau95": (rows_of([S + "/mag_validation/out/mag_cells_all.csv"] + batches(P + "/mag_validation/lowc_out_batches")),
              rows_of([R + "/mag_validation/out/mag_cells_all.csv"] + batches(R + "/mag_validation/lowc_out_batches"))),
    "tau99": (rows_of([S + "/mag_validation/out/mag_cells_tau99.csv"]), rows_of([R + "/mag_validation/out/mag_cells_tau99.csv"])),
}
exp = json.load(open(W + "/expected/resolve_impact.json"))
res, ok_all = {"row_counts": {}}, True
for name, path in (("main", R + "/mag_validation/out/mag_cells_all.csv"), ("tau99", R + "/mag_validation/out/mag_cells_tau99.csv")):
    res["row_counts"][name] = sum(1 for _ in open(path)) - 1
res["row_counts"]["n500"] = sum(sum(1 for _ in open(p)) - 1 for p in batches(R + "/mag_validation/out_batches_n500"))
res["row_counts"]["lowc_batches"] = len(batches(R + "/mag_validation/lowc_out_batches"))
res["row_counts"]["km"] = sum(sum(1 for _ in open(p)) - 1 for p in batches(R + "/mag_validation_km/out_batches"))
res["row_counts"]["rnd"] = sum(sum(1 for _ in open(p)) - 1 for p in batches(R + "/mag_validation_rnd/out_batches"))
want = {"main": 199940, "tau99": 199940, "n500": 99970, "lowc_batches": 20}
res["row_counts_ok"] = all(res["row_counts"][k] == v for k, v in want.items())
ok_all &= res["row_counts_ok"]
for arm, (O, N) in arms.items():
    blk = {}
    for key, e in exp[arm].items():
        reg, c = key.split("|"); c = float(c[1:])
        ks = [k for k in O if k[1] == reg and float(O[k]["completeness"]) == c and float(O[k]["contamination"]) == 0]
        def rate(T):
            sc = [k for k in ks if T[k]["outcome"] != "truth_empty"]
            n_ok = sum(outcome(T[k]) == "ok" for k in sc)   # resolved = ok with a non-empty prediction (aggregate.py, resolve_impact)
            return [n_ok, len(sc), round(n_ok / len(sc), 4)]
        ro, rn = rate(O), rate(N)
        chg = collections.Counter(f"chg_{outcome(O[k])}->{outcome(N[k])}" for k in ks if outcome(O[k]) != outcome(N[k]))
        match = ro == e["old"] and rn == e["new"] and dict(chg) == e["changes"]
        blk[key] = {"old": ro, "new": rn, "expected_old": e["old"], "expected_new": e["new"],
                    "changes": dict(chg), "expected_changes": e["changes"], "match": match}
        ok_all &= match
    res[arm] = blk
res["all_match"] = ok_all
json.dump(res, open(W + "/control/r22_checks.json", "w"), indent=1)
print(json.dumps(res["row_counts"]), "row_counts_ok", res["row_counts_ok"])
for arm in arms:
    for k, v in res[arm].items():
        print(arm, k, v["old"], "->", v["new"], "expected", v["expected_new"], "MATCH" if v["match"] else "DIFF", "" if v["match"] else (v["changes"], v["expected_changes"]))
print("ALL MATCH" if ok_all else "NOT ALL MATCH")
