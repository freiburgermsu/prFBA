"""Robustness: (1) precision / taf_gb decomposition; (2) leave-one-group-out on the confirm
contrasts; (3) the dominant group grp00000 in and out; (4) other bootstrap seeds and the joint
bootstrap; (5) excluding the rows whose rebuilt model differs from the prior run."""
import json
import os

import numpy as np

import vlib as L

rows = L.load_rows()
WA, WB, _ = L.weights()
POOLS = {"A": WA, "B": WB}
conf = [r for r in rows if r["half"] == "confirm"]
tune = [r for r in rows if r["half"] == "tune"]
SEL, SH = L.SELECTED, L.SHIPPED
out = {}


def C(rs, a, b, met="taf_face_mid", seed=L.SEED, pools=POOLS):
    return L.paired(rs, lambda r: L.v(r, a, met), lambda r: L.v(r, b, met), seed_base=seed,
                    pools=pools)


# ---- (1) precision decomposition (taf_gb = 1 - precision)
dec = {}
for name, (a, b) in {"selected-shipped": (SEL, SH), "est@X100-shipped": ("est/lin/G0/X100", SH),
                     "est X1-X100": (SEL, "est/lin/G0/X100"), "kn X1-X100": ("kn/lin/G0/X1", SH),
                     "two_tier X1-X100": ("two_tier/lin/G0/X1", "two_tier/lin/G0/X100"),
                     "oracle X1-X100": ("oracle/lin/G0/X1", "oracle/lin/G0/X100")}.items():
    dec[name] = {m: {p: C(conf, a, b, m)["pooled"][p] for p in ("A", "B")}
                 for m in ("precision_flux", "taf_gb", "taf_point_D", "taf_face_mid")}
out["precision_decomposition_confirm"] = dec


# ---- (2) leave-one-group-out, pooled A (point means; weights fixed)
def pooled_mean(rs, a, b, met, W):
    by = L.by_stratum(rs)
    tot = 0.0
    for s, w in W.items():
        d = [L.v(r, a, met) - L.v(r, b, met) for r in by.get(s, [])]
        if not d:
            return None
        tot += w * float(np.mean(d))
    return tot


logo = {}
for name, (a, b) in {"selected-shipped": (SEL, SH), "selected-two_tier": (SEL, "two_tier/lin/G0/X1"),
                     "selected-shuffled": (SEL, "est_shuffled/lin/G0/X1"),
                     "est@X100-shipped": ("est/lin/G0/X100", SH),
                     "shipped-two_tier@X100": (SH, "two_tier/lin/G0/X100")}.items():
    full = pooled_mean(conf, a, b, "taf_face_mid", WA)
    groups = sorted({r["group_id"] for r in conf})
    vals = []
    for g in groups:
        rs = [r for r in conf if r["group_id"] != g]
        vals.append((pooled_mean(rs, a, b, "taf_face_mid", WA), g))
    vv = sorted(vals)
    infl = sorted(((x - full, g) for x, g in vals), key=lambda t: -abs(t[0]))
    fullB = pooled_mean(conf, a, b, "taf_face_mid", WB)
    valsB = [(pooled_mean([r for r in conf if r["group_id"] != g], a, b, "taf_face_mid", WB), g)
             for g in groups]
    logo[name] = {"full_A": full, "min_A": vv[0], "max_A": vv[-1],
                  "top5_influence_A": infl[:5],
                  "n_groups": len(groups),
                  "any_drop_above_-0.01_A": [g for x, g in vals if x > -0.01],
                  "full_B": fullB, "min_B": min(valsB), "max_B": max(valsB)}
out["logo_confirm"] = logo

# ---- (3) grp00000 in / out, with CIs
G0 = "grp00000"
wo = [r for r in conf if r["group_id"] != G0]
only = [r for r in conf if r["group_id"] == G0]
g0 = {}
for name, (a, b) in {"selected-shipped": (SEL, SH), "selected-two_tier": (SEL, "two_tier/lin/G0/X1"),
                     "selected-shuffled": (SEL, "est_shuffled/lin/G0/X1"),
                     "oracle-two_tier@X1": ("oracle/lin/G0/X1", "two_tier/lin/G0/X1"),
                     "shipped-two_tier@X100": (SH, "two_tier/lin/G0/X100"),
                     "est@X100-shipped": ("est/lin/G0/X100", SH),
                     "kn X1-X100": ("kn/lin/G0/X1", SH)}.items():
    rec = {}
    for met in ("taf_face_mid", "taf_face_min", "taf_face_max"):
        cw = C(wo, a, b, met)
        rec[met] = {"without_grp00000_A": cw["pooled"]["A"], "without_grp00000_B": cw["pooled"]["B"],
                    "without_per": {s: {"n": p["n"], "K": p["K"], "mean": p["mean"], "ci": p["ci"]}
                                    for s, p in cw["per"].items()}}
        # grp00000 alone: plain means per stratum, pooled A
        by = L.by_stratum(only)
        per = {s: float(np.mean([L.v(r, a, met) - L.v(r, b, met) for r in by[s]])) for s in L.STRATA
               if by.get(s)}
        rec[met]["grp00000_only_per_stratum"] = per
        rec[met]["grp00000_only_A"] = sum(WA[s] * per[s] for s in L.MULTI)
    g0[name] = rec
out["grp00000"] = {"n_rows": len(only), "by_stratum": {s: len(L.by_stratum(only).get(s, []))
                                                        for s in L.STRATA},
                   "contrasts": g0}

# ---- (4) bootstrap seeds and joint bootstrap for the headline
seeds = {}
for sd in (1729, 424242, 1, 2026, 99991):
    c = C(conf, SEL, SH, seed=sd)
    seeds[str(sd)] = {"A": c["pooled"]["A"]["ci"], "B": c["pooled"]["B"]["ci"],
                      "n1": c["per"]["n1"]["ci"]}
jA = L.joint_boot(conf, lambda r: L.v(r, SEL, "taf_face_mid"), lambda r: L.v(r, SH, "taf_face_mid"), WA)
jB = L.joint_boot(conf, lambda r: L.v(r, SEL, "taf_face_mid"), lambda r: L.v(r, SH, "taf_face_mid"), WB)
# bigger B for stability of percentile ends
big = C(conf, SEL, SH, seed=1729)
bigB = L.paired(conf, lambda r: L.v(r, SEL, "taf_face_mid"), lambda r: L.v(r, SH, "taf_face_mid"),
                seed_base=1729, pools=POOLS, B_=20000)
out["bootstrap_sensitivity_selected_shipped_mid"] = {
    "seeds": seeds, "joint_A": jA, "joint_B": jB,
    "B20000_seed1729": {"A": bigB["pooled"]["A"]["ci"], "B": bigB["pooled"]["B"]["ci"]}}
# how many grp00000 copies per stratum draw, and the conditional means (lumpiness)
by = L.by_stratum(conf)
lump = {}
for i, s in enumerate(L.STRATA):
    keys = sorted({r["group_id"] for r in by[s]})
    K = len(keys)
    j = keys.index(G0)
    dr = np.random.default_rng(L.SEED + 1000 * i).integers(0, K, size=(L.B, K))
    cnt = (dr == j).sum(1)
    lump[s] = {"K": K, "share_rows_grp00000": sum(1 for r in by[s] if r["group_id"] == G0) / len(by[s]),
               "draws_with_0_copies": float(np.mean(cnt == 0)), "draws_with_ge2": float(np.mean(cnt >= 2))}
out["grp00000_bootstrap_presence"] = lump

# ---- (5) excluding rows whose rebuilt model differs from the prior run
EXC2 = {"1674805.4__V3-V4", "1798454.3__V4"}
EXC3 = EXC2 | {"1017270.21__V4-V5"}
EXC4 = EXC3 | {"2053306.99__V4"}
exc = {}
for tag, ex in (("excl_2_prior_match_false", EXC2), ("excl_3_D_not_bitwise", EXC3),
                ("excl_4_face_diff_gt_1e-9", EXC4)):
    t2 = [r for r in tune if r["key"] not in ex]
    c2 = [r for r in conf if r["key"] not in ex]
    # selection on the reduced tune half: lowest pooled A among the 36 (G0 rule kept simple:
    # the representative of each family is its G0 arm unless another G level is lower by >0.002)
    by = L.by_stratum(t2)
    lev = {}
    for a in t2[0]["arms"]:
        lev[a] = sum(WA[s] * float(np.mean([L.v(r, a, "taf_face_mid") for r in by[s]])) for s in L.MULTI)
    cands = [a for a in lev if a.split("/")[0] in ("kn", "est")]
    best = min(cands, key=lambda a: lev[a])
    rec = {"n_tune": len(t2), "n_confirm": len(c2), "lowest_A_candidate_tune": best,
           "lowest_A": lev[best], "selected_A": lev[SEL], "shipped_A": lev[SH],
           "est/lin/G0/X10_A": lev["est/lin/G0/X10"]}
    for name, (a, b) in {"selected-shipped": (SEL, SH), "selected-two_tier": (SEL, "two_tier/lin/G0/X1"),
                         "selected-shuffled": (SEL, "est_shuffled/lin/G0/X1"),
                         "shipped-two_tier@X100": (SH, "two_tier/lin/G0/X100")}.items():
        cc = C(c2, a, b)
        rec[name] = {"A": cc["pooled"]["A"], "B": cc["pooled"]["B"]}
    exc[tag] = rec
out["exclusions"] = exc

json.dump(out, open(os.path.join(L.HERE, "out_c8.json"), "w"), indent=1, default=str)
print("== precision decomposition (confirm, A / B)")
for k, d in dec.items():
    print("  %-20s prec A %s B %s | taf_gb A %s | taf_point_D A %s" % (
        k, L.fmt(d["precision_flux"]["A"]), L.fmt(d["precision_flux"]["B"]), L.fmt(d["taf_gb"]["A"]),
        L.fmt(d["taf_point_D"]["A"])))
print("== LOGO (pooled A face mid)")
for k, d in logo.items():
    print("  %-22s full %.4f range [%.4f (%s), %.4f (%s)] top: %s ; drops above -0.01: %s ; B full %.4f range [%.4f, %.4f]" % (
        k, d["full_A"], d["min_A"][0], d["min_A"][1], d["max_A"][0], d["max_A"][1],
        [(round(x, 4), g) for x, g in d["top5_influence_A"][:3]], d["any_drop_above_-0.01_A"],
        d["full_B"], d["min_B"][0], d["max_B"][0]))
print("== grp00000:", out["grp00000"]["n_rows"], out["grp00000"]["by_stratum"])
for k, d in g0.items():
    m = d["taf_face_mid"]
    print("  %-22s without A %s B %s | grp00000-only A %.4f per %s" % (
        k, L.fmt(m["without_grp00000_A"]), L.fmt(m["without_grp00000_B"]), m["grp00000_only_A"],
        {s: round(x, 4) for s, x in m["grp00000_only_per_stratum"].items()}))
    print("  %-22s without: min A %s ; max A %s ; per-stratum mid %s" % (
        "", L.fmt(d["taf_face_min"]["without_grp00000_A"]), L.fmt(d["taf_face_max"]["without_grp00000_A"]),
        {s: (p["n"], p["K"], round(p["mean"], 4), [round(x, 4) for x in p["ci"]]) for s, p in m["without_per"].items()}))
print("== seeds", json.dumps(out["bootstrap_sensitivity_selected_shipped_mid"], default=str))
print("== grp00000 presence in draws", lump)
print("== exclusions")
for k, d in exc.items():
    print("  ", k, {kk: (L.fmt(vv["A"]) if isinstance(vv, dict) and "A" in vv else vv) for kk, vv in d.items()})
