"""Check the filled doc's guard table (3 dp) and checks table against own numbers."""
import json, os, re
import vlib as L
rows = L.load_rows(); prior = {r["key"]: r for r in L.load_prior()}
o = json.load(open(os.path.join(L.HERE, "out_c2_c3.json")))
md = open(os.path.join(L.SR, "data/objective_tuning/OBJECTIVE_TUNING_RESULTS.filled.md")).read()
bad = []
n = 0
for line in md.splitlines():
    m = re.match(r"\| ((?:kn|est)/(?:lin|log)/G[^/]+/X\d+) \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \| (yes|no) \| ([\d.]+) \| ([\d.]+) \|", line)
    if not m:
        continue
    a = m.group(1); e = o["eligibility"]["per_config"][a]
    mine = ("%.3f" % e["ratio_of_medians"]["flux_exchange"], "%.3f" % e["ratio_of_medians"]["flux_gene_backed"],
            "%.3f" % e["ratio_of_medians"]["precision_flux"], "yes" if e["eligible"] else "no",
            "%.4f" % o["levels_A"][a], None)
    theirs = m.groups()[1:6]
    n += 1
    if mine[:5] != theirs:
        bad.append((a, mine[:5], theirs))
print("guard-table rows checked", n, "mismatches", bad)
# checks table
pm = [r["prior_match"] for r in rows]
def cnt(name, k):
    xs = [((x.get(name) or {}).get(k) or {}) for x in pm]
    return sum(1 for y in xs if y.get("bitwise"))
print("shipped vertex gb/absent bitwise", cnt("shipped_vs_prior_real", "flux_gene_backed"), cnt("shipped_vs_prior_real", "flux_absent"))
# own: compare with prior file directly
vb = sum(1 for r in rows if r["arms"][L.SHIPPED]["flux_gene_backed"] == prior[r["key"]]["arms"]["real"]["flux_gene_backed"])
va = sum(1 for r in rows if r["arms"][L.SHIPPED]["flux_absent"] == prior[r["key"]]["arms"]["real"]["flux_absent"])
tto = sum(1 for r in rows if r["arms"]["two_tier/lin/G0/X100"]["objective"] == prior[r["key"]]["arms"]["two_tier"]["objective"])
ttf = sum(1 for r in rows if all(abs(r["arms"]["two_tier/lin/G0/X100"][k] - prior[r["key"]]["arms"]["two_tier"][k]) <= 1e-9 for k in ("taf_face_min", "taf_face_mid", "taf_face_max")))
ttmax = max(abs(r["arms"]["two_tier/lin/G0/X100"]["taf_face_mid"] - prior[r["key"]]["arms"]["two_tier"]["taf_face_mid"]) for r in rows)
shf = sum(1 for r in rows if all(abs(r["arms"][L.SHIPPED][k] - prior[r["key"]]["arms"]["real"][k]) <= 1e-9 for k in ("taf_face_min", "taf_face_mid", "taf_face_max")))
print("own: shipped vertex gb bitwise", vb, "absent bitwise", va, "| tt objective bitwise", tto, "tt face within 1e-9", ttf, "max mid diff", ttmax, "| shipped face within 1e-9", shf)
print("retried", sum(x["face_lp_retried"] or 0 for r in rows for x in r["arms"].values()))
print("gapf&gb", sum(r["model"]["n_gapfilled_and_gene_backed"] for r in rows), "alpha>0 & null clade", sum(1 for r in rows if r["est_info"]["alpha_with_null_clade"]))
print("released pass", sum(1 for r in rows if (r["released_check"] or {}).get("pass")), "released==prior bitwise", sum(1 for r in rows if (r["released_check"] or {}).get("equals_prior_released_bitwise")),
      "coef map identical", sum(1 for r in rows if (r["shipped_check"] or {}).get("coef_map_identical")))
print("which released != prior:", [r["key"] for r in rows if not (r["released_check"] or {}).get("equals_prior_released_bitwise")])
# B ranking order check
lev = {}
WA, WB, _ = L.weights()
by = L.by_stratum([r for r in rows if r["half"] == "tune"])
import numpy as np
for a in ("est/lin/G0/X1", "est/lin/Gmeas/X1", "est/lin/G0.2/X1"):
    lev[a] = sum(WB[s] * float(np.mean([L.v(r, a, "taf_face_mid") for r in by[s]])) for s in WB)
print("B levels", {k: "%.10f" % v for k, v in lev.items()})
