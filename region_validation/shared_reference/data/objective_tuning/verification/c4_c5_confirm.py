"""Claims 4 (confirmation) and 5 (decomposition), confirm half, own code; then a
comparison of every number with the author's summary.json."""
import json
import os

import numpy as np

import vlib as L

rows = L.load_rows()
WA, WB, _ = L.weights()
POOLS = {"A": WA, "B": WB}
conf = [r for r in rows if r["half"] == "confirm"]
assert all(L.half_of(r["group_id"]) == "confirm" for r in conf)
SEL, SH = L.SELECTED, L.SHIPPED
PAIRS = {"selected-shipped": (SEL, SH),
         "selected-two_tier": (SEL, "two_tier/lin/G0/X1"),
         "selected-shuffled": (SEL, "est_shuffled/lin/G0/X1"),
         "oracle-two_tier@selected_cell": ("oracle/lin/G0/X1", "two_tier/lin/G0/X1"),
         "shipped-two_tier@shipped_cell": (SH, "two_tier/lin/G0/X100"),
         "oracle-two_tier@shipped_cell": ("oracle/lin/G0/X100", "two_tier/lin/G0/X100")}
DECOMP = {"est_alone@X100 (est/lin/G0/X100 - kn/lin/G0/X100)": ("est/lin/G0/X100", SH),
          "X1_alone@kn (kn/lin/G0/X1 - kn/lin/G0/X100)": ("kn/lin/G0/X1", SH),
          "X1_given_est (est/lin/G0/X1 - est/lin/G0/X100)": (SEL, "est/lin/G0/X100"),
          "est_given_X1 (est/lin/G0/X1 - kn/lin/G0/X1)": (SEL, "kn/lin/G0/X1"),
          "log_vs_lin@est/X1 (est/log/G0/X1 - est/lin/G0/X1)": ("est/log/G0/X1", SEL),
          "Gmeas_vs_G0@est/lin/X1": ("est/lin/Gmeas/X1", SEL),
          "G0.2_vs_G0@est/lin/X1": ("est/lin/G0.2/X1", SEL),
          "two_tier X1 - X100": ("two_tier/lin/G0/X1", "two_tier/lin/G0/X100"),
          "oracle X1 - X100": ("oracle/lin/G0/X1", "oracle/lin/G0/X100"),
          "est_shuffled/lin/G0/X1 - shipped": ("est_shuffled/lin/G0/X1", SH),
          "est_shuffled/lin/G0/X1 - two_tier/lin/G0/X1": ("est_shuffled/lin/G0/X1",
                                                          "two_tier/lin/G0/X1")}
METS = ("taf_face_mid", "taf_face_min", "taf_face_max", "taf_gb", "taf_point_D",
        "flux_exchange", "flux_gene_backed", "precision_flux", "flux_gapfilled",
        "taf_face_width")


def contrast(rows_, a, b, met, seed=L.SEED):
    return L.paired(rows_, lambda r: L.v(r, a, met), lambda r: L.v(r, b, met),
                    seed_base=seed, pools=POOLS)


res = {}
for name, (a, b) in {**PAIRS, **DECOMP}.items():
    res[name] = {}
    for met in METS:
        c = contrast(conf, a, b, met)
        c2 = contrast(conf, a, b, met, seed=L.ALT_SEED)
        res[name][met] = {"per": c["per"], "pooled": c["pooled"],
                          "pooled_alt_seed": c2["pooled"],
                          "per_alt_seed": {s: c2["per"][s]["ci"] for s in c2["per"]}}
        if met in ("taf_face_mid", "taf_face_min", "taf_face_max"):
            jp = L.joint_boot(conf, lambda r: L.v(r, a, met), lambda r: L.v(r, b, met), WA)
            res[name][met]["joint_A"] = {"mean": jp[0], "ci": jp[1], "K": jp[2],
                                         "n_stratum_missing_draws": jp[3]}
            jb = L.joint_boot(conf, lambda r: L.v(r, a, met), lambda r: L.v(r, b, met), WB)
            res[name][met]["joint_B"] = {"mean": jb[0], "ci": jb[1]}


# headroom (ratio of pooled-A means)
def pm(name, pool="A", met="taf_face_mid"):
    return res[name][met]["pooled"][pool]["mean"]


hr = {pool: {"selected_own_cell": pm("selected-two_tier", pool) / pm("oracle-two_tier@selected_cell", pool),
             "shipped": pm("shipped-two_tier@shipped_cell", pool) / pm("oracle-two_tier@shipped_cell", pool),
             "selected_in_shipped_frame": (pm("selected-shipped", pool) + pm("shipped-two_tier@shipped_cell", pool))
             / pm("oracle-two_tier@shipped_cell", pool)} for pool in ("A", "B")}

# additivity check of the decomposition
add = {"est@X100 + X1@kn": pm("est_alone@X100 (est/lin/G0/X100 - kn/lin/G0/X100)") +
       pm("X1_alone@kn (kn/lin/G0/X1 - kn/lin/G0/X100)"),
       "joint (selected-shipped)": pm("selected-shipped"),
       "interaction": pm("selected-shipped") - pm("est_alone@X100 (est/lin/G0/X100 - kn/lin/G0/X100)")
       - pm("X1_alone@kn (kn/lin/G0/X1 - kn/lin/G0/X100)")}

# ------------------------------------------------ compare with the author's summary.json
S = json.load(open(os.path.join(L.SR, "data/objective_tuning/summary.json")))
SC = S["confirm"]["contrasts"]
maxd, n, worst = 0.0, 0, None
for name in PAIRS:
    for met in METS:
        if met not in SC[name]["metrics"]:
            continue
        theirs = SC[name]["metrics"][met]
        mine = res[name][met]
        for pool in ("A", "B"):
            t = theirs["pooled"][pool]
            m = mine["pooled"][pool]
            for x, y in ((m["mean"], t["mean"]), (m["ci"][0], t["ci"][0]), (m["ci"][1], t["ci"][1])):
                n += 1
                if abs(x - y) > maxd:
                    maxd, worst = abs(x - y), (name, met, pool, x, y)
        for s, t in theirs["per_stratum"].items():
            m = mine["per"][s]
            for x, y in ((m["mean"], t["mean"]), (m["ci"][0], t["ci"][0]), (m["ci"][1], t["ci"][1])):
                n += 1
                if abs(x - y) > maxd:
                    maxd, worst = abs(x - y), (name, met, s, x, y)
cmp = {"n_numbers_compared": n, "max_abs_diff": maxd, "worst": worst,
       "author_headroom_confirm": S["confirm"]["headroom_confirm"], "mine": hr}

out = {"contrasts": res, "headroom": hr, "additivity": add, "cmp_vs_summary": cmp}
json.dump(out, open(os.path.join(L.HERE, "out_c4_c5.json"), "w"), indent=1, default=str)

print("compared", n, "numbers with summary.json; max abs diff", maxd, worst)
print("headroom", json.dumps(hr, indent=0))
print("author headroom", S["confirm"]["headroom_confirm"])
for name in list(PAIRS) + list(DECOMP):
    r = res[name]
    print("%-58s mid %s | min %s | max %s" % (name, L.fmt(r["taf_face_mid"]["pooled"]["A"]),
                                             L.fmt(r["taf_face_min"]["pooled"]["A"]),
                                             L.fmt(r["taf_face_max"]["pooled"]["A"])))
    print("%-58s B mid %s ; altseed A %s ; jointA %s ; jointB %s" % (
        "", L.fmt(r["taf_face_mid"]["pooled"]["B"]),
        "[%+.4f, %+.4f]" % tuple(r["taf_face_mid"]["pooled_alt_seed"]["A"]["ci"]),
        "[%+.4f, %+.4f]" % tuple(r["taf_face_mid"]["joint_A"]["ci"]),
        "[%+.4f, %+.4f]" % tuple(r["taf_face_mid"]["joint_B"]["ci"])))
    print("%-58s per-stratum mid: %s" % ("", "; ".join(
        "%s %s (n=%d,K=%d, <0:%d =0:%d >0:%d)" % (s, L.fmt(p), p["n"], p["K"], p["n_lt0"], p["n_eq0"], p["n_gt0"])
        for s, p in r["taf_face_mid"]["per"].items())))
print("additivity", add)
g = res["selected-shipped"]
for met in ("flux_exchange", "flux_gene_backed", "precision_flux", "flux_gapfilled", "taf_gb",
            "taf_point_D", "taf_face_width"):
    print("guard sel-ship %-18s A %s  B %s" % (met, L.fmt(g[met]["pooled"]["A"]), L.fmt(g[met]["pooled"]["B"])))
