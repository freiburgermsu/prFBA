"""Claims 2 (eligibility) and 3 (selection), tune half only, own code."""
import json
import os
from itertools import product

import numpy as np

import vlib as L

rows = L.load_rows()
WA, WB, _e = L.weights()
tune = [r for r in rows if r["half"] == "tune"]
conf = [r for r in rows if r["half"] == "confirm"]
assert all(L.half_of(r["group_id"]) == "tune" for r in tune)
by = L.by_stratum(tune)
P_C, C_S, G_L, X_L = ("kn", "est"), ("lin", "log"), ("G0", "Gmeas", "G0.2"), ("X100", "X10", "X1")
CAND = ["/".join(t) for t in product(P_C, C_S, G_L, X_L)]
assert len(CAND) == 36
out = {"weights_A": WA, "weights_B": WB}


def med(rs, arm, met):
    xs = [L.v(r, arm, met) for r in rs]
    xs = [x for x in xs if x is not None]
    return float(np.median(xs)), len(xs)


def level(byS, arm, met, W, renorm=False):
    per = {}
    for s in W:
        xs = [L.v(r, arm, met) for r in byS.get(s, [])]
        xs = [x for x in xs if x is not None]
        if xs:
            per[s] = float(np.mean(xs))
    if renorm:
        tw = sum(W[s] for s in per)
        return (sum(W[s] / tw * per[s] for s in per) if per else None), per
    if len(per) < len(W):
        return None, per
    return sum(W[s] * per[s] for s in W), per


def deps(a):
    return sum(x != y for x, y in zip(a.split("/"), L.SHIPPED.split("/")))


# ---------------------------------------------------------------- eligibility
sh = {m: med(tune, L.SHIPPED, m)[0] for m in ("flux_exchange", "flux_gene_backed", "precision_flux")}
elig = {}
for c in CAND:
    rat = {m: med(tune, c, m)[0] / sh[m] for m in sh}
    nopt = sum(1 for r in tune if r["arms"][c]["status"] == "optimal")
    floor = all(abs(r["arms"][c]["biomass_ratio"] - 0.95) <= 1e-6 for r in tune)
    rowrat = {}
    for m in ("flux_exchange", "flux_gene_backed"):
        rr = np.array([L.v(r, c, m) / L.v(r, L.SHIPPED, m) for r in tune])
        rowrat[m] = float(np.median(rr))
    ok = (abs(rat["flux_exchange"] - 1) <= 0.10 and abs(rat["flux_gene_backed"] - 1) <= 0.10
          and nopt == len(tune) and floor)
    elig[c] = {"ratio_of_medians": rat, "median_row_ratio": rowrat, "n_optimal": nopt,
               "floor": floor, "eligible": ok}
E = [c for c in CAND if elig[c]["eligible"]]
out["eligibility"] = {"n_eligible": len(E), "shipped_medians": sh,
                      "max_abs_ratio_dev": {m: max(abs(elig[c]["ratio_of_medians"][m] - 1)
                                                   for c in CAND)
                                            for m in ("flux_exchange", "flux_gene_backed")},
                      "max_abs_rowratio_dev": {m: max(abs(elig[c]["median_row_ratio"][m] - 1)
                                                      for c in CAND)
                                               for m in ("flux_exchange", "flux_gene_backed")},
                      "per_config": elig}


# per-row ratio distribution for the selected arm (tune and confirm)
def ratio_dist(rs, arm, m):
    rr = np.array([L.v(r, arm, m) / L.v(r, L.SHIPPED, m) for r in rs])
    q = np.percentile(rr, [0, 1, 5, 25, 50, 75, 95, 99, 100])
    return {"quantiles_0_1_5_25_50_75_95_99_100": [float(x) for x in q],
            "frac_outside_10pct": float(np.mean(np.abs(rr - 1) > 0.10)),
            "frac_outside_5pct": float(np.mean(np.abs(rr - 1) > 0.05)),
            "frac_outside_25pct": float(np.mean(np.abs(rr - 1) > 0.25)),
            "n_outside_10pct": int(np.sum(np.abs(rr - 1) > 0.10)), "n": len(rr),
            "frac_exactly_1": float(np.mean(rr == 1.0))}


out["selected_row_ratios"] = {h: {m: ratio_dist(rs, L.SELECTED, m) for m in
                                  ("flux_exchange", "flux_gene_backed", "flux_rxn",
                                   "flux_gapfilled")}
                              for h, rs in (("tune", tune), ("confirm", conf))}

# ---------------------------------------------------------------- levels
lev = {}
allarms = sorted(tune[0]["arms"].keys())
for a in allarms:
    la, pa = level(by, a, "taf_face_mid", WA)
    lb, _ = level(by, a, "taf_face_mid", WB)
    lev[a] = {"A": la, "B": lb, "per": pa}
rankA = sorted(E, key=lambda a: lev[a]["A"])
rankB = sorted(E, key=lambda a: lev[a]["B"])
out["levels_A"] = {a: lev[a]["A"] for a in rankA}
out["levels_B_top10"] = [(a, lev[a]["B"]) for a in rankB[:10]]
out["controls_A"] = {a: lev[a]["A"] for a in allarms if a not in CAND}

# lin vs log at every (P, G, X)
lvl = {}
for P, G, X in product(P_C, G_L, X_L):
    a, b = f"{P}/lin/{G}/{X}", f"{P}/log/{G}/{X}"
    lvl[f"{P}/*/{G}/{X}"] = {"lin": lev[a]["A"], "log": lev[b]["A"],
                              "log_minus_lin": lev[b]["A"] - lev[a]["A"]}
out["log_vs_lin_tune_A"] = lvl
out["log_ever_better_than_lin"] = [k for k, x in lvl.items() if x["log_minus_lin"] < 0]

# ---------------------------------------------------------------- G reading
fam = {}
for P, C, X in product(P_C, C_S, X_L):
    arms = [f"{P}/{C}/{G}/{X}" for G in G_L]
    vary = []
    for r in tune:
        if r["stratum"] not in L.MULTI:
            continue
        mids = [(L.v(r, a, "gapfilled_flux_face_min") + L.v(r, a, "gapfilled_flux_face_max")) / 2
                for a in arms]
        Dv = r["model"]["D_flux_gb_two_tier"]
        if max(mids) - min(mids) > 1e-6 * max(1.0, Dv):
            vary.append(r)
    vb = L.by_stratum(vary)
    vals = {}
    for a in arms:
        vals[a] = level(vb, a, "taf_face_mid", WA, renorm=True)[0] if vary else None
    rec = {"n_vary": {s: len(vb.get(s, [])) for s in L.MULTI}, "vary_keys": [r["key"] for r in vary]}
    if vary:
        best = min(vals.values())
        ties = [a for a in arms if vals[a] <= best + 0.002]
        rec["A_on_vary"] = vals
        rec["max_minus_min_A_on_vary"] = max(vals.values()) - min(vals.values())
        rec["choice"] = min(ties, key=lambda a: (deps(a), vals[a]))
    else:
        rec["choice"] = f"{P}/{C}/G0/{X}"
    allA = {a: lev[a]["A"] for a in arms}
    rec["A_all_rows"] = allA
    rec["max_minus_min_A_all_rows"] = max(allA.values()) - min(allA.values())
    # rows where the TAF face mid differs across G (any stratum)
    rec["n_rows_taf_mid_differs_across_G_gt1e-9"] = sum(
        1 for r in tune if max(L.v(r, a, "taf_face_mid") for a in arms) -
        min(L.v(r, a, "taf_face_mid") for a in arms) > 1e-9)
    fam[f"{P}/{C}/*/{X}"] = rec
out["G_reading"] = fam
reps = [f["choice"] for f in fam.values()]
A = {a: lev[a]["A"] for a in reps}
best = min(A.values())
ties = [a for a in A if A[a] <= best + 0.002]
sel = min(ties, key=lambda a: (deps(a), A[a]))
out["selection"] = {"selected": sel, "pooled_A": A[sel], "pooled_B": lev[sel]["B"],
                    "shipped_A": lev[L.SHIPPED]["A"], "shipped_B": lev[L.SHIPPED]["B"],
                    "tie_set": ties, "representatives_sorted": sorted(A.items(), key=lambda t: t[1]),
                    "margin_to_next_family": sorted(A.values())[1] - best}
# naive joint selection over all eligible
AN = {a: lev[a]["A"] for a in E}
bn = min(AN.values())
tn = [a for a in AN if AN[a] <= bn + 0.002]
out["naive_selection"] = {"selected": min(tn, key=lambda a: (deps(a), AN[a])), "tie_set": tn}
# pooling-B selection (not used; sensitivity)
ABv = {a: lev[a]["B"] for a in E}
bb = min(ABv.values())
tb = [a for a in ABv if ABv[a] <= bb + 0.002]
out["B_selection_sensitivity"] = {"selected": min(tb, key=lambda a: (deps(a), ABv[a])),
                                  "tie_set": tb}

# headroom (tune), ratio of pooled-A levels
s_, tt, orc = lev[sel]["A"], lev["two_tier/lin/G0/X1"]["A"], lev["oracle/lin/G0/X1"]["A"]
tt0, orc0, sh0 = lev["two_tier/lin/G0/X100"]["A"], lev["oracle/lin/G0/X100"]["A"], lev[L.SHIPPED]["A"]
out["headroom_tune"] = {"own_cell": (s_ - tt) / (orc - tt), "shipped_frame": (s_ - tt0) / (orc0 - tt0),
                        "shipped": (sh0 - tt0) / (orc0 - tt0)}
# tune-half paired contrast selected - shipped (winner's curse check)
c = L.paired(tune, lambda r: L.v(r, sel, "taf_face_mid"), lambda r: L.v(r, L.SHIPPED, "taf_face_mid"),
             pools={"A": WA, "B": WB})
out["tune_selected_minus_shipped"] = {"A": c["pooled"]["A"], "B": c["pooled"]["B"],
                                      "per": c["per"]}

json.dump(out, open(os.path.join(L.HERE, "out_c2_c3.json"), "w"), indent=1, default=str)
print("eligible", len(E), "max|ratio-1|", out["eligibility"]["max_abs_ratio_dev"],
      "rowratio", out["eligibility"]["max_abs_rowratio_dev"])
print("selection", json.dumps(out["selection"], indent=0, default=str))
print("naive", out["naive_selection"], "B-sens", out["B_selection_sensitivity"])
print("rank A top 12:")
for a in rankA[:12]:
    print("  %-22s A=%.6f B=%.6f" % (a, lev[a]["A"], lev[a]["B"]))
print("shipped A=%.6f B=%.6f" % (lev[L.SHIPPED]["A"], lev[L.SHIPPED]["B"]))
print("log ever better than lin:", out["log_ever_better_than_lin"])
for k, rec in fam.items():
    print("  G %-18s vary=%s choice=%s maxdiff_vary=%s maxdiff_all=%.6f n_rows_taf_differs=%d" % (
        k, rec["n_vary"], rec["choice"], ("%.5f" % rec["max_minus_min_A_on_vary"]) if
        "max_minus_min_A_on_vary" in rec else "-", rec["max_minus_min_A_all_rows"],
        rec["n_rows_taf_mid_differs_across_G_gt1e-9"]))
print("headroom_tune", out["headroom_tune"])
print("tune sel-ship A", L.fmt(out["tune_selected_minus_shipped"]["A"]), "B",
      L.fmt(out["tune_selected_minus_shipped"]["B"]))
for h in ("tune", "confirm"):
    for m in ("flux_exchange", "flux_gene_backed", "flux_rxn", "flux_gapfilled"):
        d = out["selected_row_ratios"][h][m]
        print(h, m, ["%.4f" % x for x in d["quantiles_0_1_5_25_50_75_95_99_100"]],
              "out10%%=%d/%d out5%%=%.3f out25%%=%.3f eq1=%.3f" % (
                  d["n_outside_10pct"], d["n"], d["frac_outside_5pct"], d["frac_outside_25pct"],
                  d["frac_exactly_1"]))
