"""Claim 1 (split) and claim 6 (health), plus row-internal consistency checks."""
import json
import os
from collections import Counter, defaultdict

import numpy as np

import vlib as L

rows = L.load_rows()
prior = {r["key"]: r for r in L.load_prior()}
pop = json.load(open(os.path.join(L.SR, "data/weighted_objective/population.json")))
popjobs = {j["key"]: j for j in pop["jobs"]}
out = {}

# ---------------------------------------------------------------- claim 1: split
mism_half = [r["key"] for r in rows if L.half_of(r["group_id"]) != r["half"]]
cnt = Counter(L.half_of(r["group_id"]) for r in rows)
g_half = defaultdict(set)
for r in rows:
    g_half[r["group_id"]].add(r["half"])
span = [g for g, h in g_half.items() if len(h) > 1]
gid_vs_pop = [r["key"] for r in rows if popjobs.get(r["key"], {}).get("group_id") != r["group_id"]]
keys_missing_pop = [r["key"] for r in rows if r["key"] not in popjobs]
by_half_stratum = Counter((r["half"], r["stratum"]) for r in rows)
groups_by_half = Counter(L.half_of(g) for g in g_half)
# do other identities cross halves? (genus / source genome / amplicon md5 / fold)
def cross(field):
    m = defaultdict(set)
    for r in rows:
        m[r[field]].add(r["half"])
    return sorted(k for k, h in m.items() if len(h) > 1)
out["claim1"] = {
    "rows_by_half_recomputed": dict(cnt),
    "rows_half_field_mismatch": mism_half,
    "n_groups": len(g_half), "groups_by_half": dict(groups_by_half),
    "groups_spanning_halves": span,
    "group_id_differs_from_population_json": gid_vs_pop,
    "keys_missing_from_population": keys_missing_pop,
    "rows_by_half_stratum": {f"{h}|{s}": n for (h, s), n in sorted(by_half_stratum.items())},
    "groups_by_half_stratum": {f"{h}|{s}": len({r['group_id'] for r in rows
                                                if r['half'] == h and r['stratum'] == s})
                               for h in ("tune", "confirm") for s in L.STRATA},
    "genera_in_both_halves": cross("genus"),
    "src_genomes_in_both_halves": cross("src_genome_id"),
    "md5_in_both_halves": cross("md5"),
    "folds_by_half": {h: dict(Counter(r["fold"] for r in rows if r["half"] == h))
                      for h in ("tune", "confirm")},
}
# groups that hold rows in >1 stratum (bootstrap independence across strata)
g_str = defaultdict(set)
for r in rows:
    g_str[r["group_id"]].add(r["stratum"])
out["claim1"]["groups_with_rows_in_multiple_strata"] = {
    h: sum(1 for g, s in g_str.items() if len(s) > 1 and L.half_of(g) == h)
    for h in ("tune", "confirm")}
# largest group sizes per half
gsz = Counter(r["group_id"] for r in rows)
out["claim1"]["largest_groups"] = {h: sorted(((n, g) for g, n in gsz.items()
                                              if L.half_of(g) == h), reverse=True)[:5]
                                   for h in ("tune", "confirm")}

# ---------------------------------------------------------------- claim 6: health
st = Counter(r["status"] for r in rows)
narms = Counter(len(r["arms"]) for r in rows)
solves = [(r["key"], a, x["status"]) for r in rows for a, x in r["arms"].items()]
nonopt = [s for s in solves if s[2] != "optimal"]
br = [x["biomass_ratio"] for r in rows for x in r["arms"].values()]
pm_fail = [r for r in rows if not (r.get("prior_match") or {}).get("ok")]
pm = []
for r in pm_fail:
    p = prior.get(r["key"])
    pm.append({"key": r["key"], "half": r["half"], "stratum": r["stratum"],
               "group": r["group_id"],
               "mismatches": r["prior_match"].get("mismatches"),
               "extra_mismatches": r["prior_match"].get("extra_mismatches"),
               "gapfill_now": r["gapfill"], "gapfill_prior": (p or {}).get("gapfill"),
               "D_now": r["model"].get("D_flux_gb_two_tier"),
               "D_prior": ((p or {}).get("model") or {}).get("D_flux_gb_two_tier")})
# D bitwise / shipped vs prior real, by my own comparison with the prior file
dmis, ship_obj_mis, ship_mid = [], [], []
for r in rows:
    p = prior[r["key"]]
    Dn, Dp = r["model"].get("D_flux_gb_two_tier"), (p["model"] or {}).get("D_flux_gb_two_tier")
    if Dn != Dp:
        dmis.append((r["key"], r["half"], Dn, Dp))
    a, b = r["arms"][L.SHIPPED], p["arms"].get("real") or {}
    if a["objective"] != b.get("objective"):
        ship_obj_mis.append((r["key"], r["half"], a["objective"], b.get("objective")))
    if b.get("taf_face_mid") is not None:
        ship_mid.append((abs(a["taf_face_mid"] - b["taf_face_mid"]), r["key"], r["half"]))
ship_mid.sort(reverse=True)
out["claim6"] = {
    "status": dict(st), "arms_per_row": dict(narms), "n_solves": len(solves),
    "n_nonoptimal": len(nonopt),
    "biomass_ratio_min_max": [min(br), max(br)],
    "n_off_floor_gt_1e-6": sum(1 for x in br if abs(x - 0.95) > 1e-6),
    "face_lp_failed": sum(x["face_lp_failed"] or 0 for r in rows for x in r["arms"].values()),
    "prior_match_fail": pm,
    "D_not_bitwise_prior": dmis,
    "shipped_objective_not_bitwise_prior_real": ship_obj_mis,
    "shipped_face_mid_absdiff_vs_prior_top5": ship_mid[:5],
    "n_shipped_face_mid_absdiff_gt_1e-9": sum(1 for x in ship_mid if x[0] > 1e-9),
    "determinism": {"pairs": sum((r["determinism"] or {}).get("n_pairs", 0) for r in rows),
                    "identical": sum((r["determinism"] or {}).get("n_identical", 0) for r in rows)},
}

# ---------------------------------------------------------------- row-internal consistency
bad = Counter()
maxdev = Counter()
for r in rows:
    D = r["model"]["D_flux_gb_two_tier"]
    for a, x in r["arms"].items():
        mid = (x["absent_flux_min"] + x["absent_flux_max"]) / 2 / D
        dv = abs(mid - x["taf_face_mid"])
        maxdev["mid"] = max(maxdev["mid"], dv)
        dv2 = abs(x["absent_flux_min"] / D - x["taf_face_min"])
        maxdev["min"] = max(maxdev["min"], dv2)
        # vertex inside face?
        if not (x["absent_flux_min"] - 1e-7 <= x["flux_absent"] <= x["absent_flux_max"] + 1e-7):
            bad["vertex_outside_face"] += 1
        # exchange_term / X vs flux_exchange
        dex = abs(x["exchange_term"] / x["x_value"] - x["flux_exchange"])
        maxdev["ex_rel"] = max(maxdev["ex_rel"], dex / max(1.0, x["flux_exchange"]))
        sp = x["internal_term_split"] or {}
        geneless = x["flux_rxn"] - x["flux_gene_backed"] - x["flux_gapfilled"]
        maxdev["geneless_vs_split"] = max(maxdev["geneless_vs_split"],
                                          abs(geneless - sp.get("geneless_draft", 0.0)))
        if x["precision_flux"] is not None:
            maxdev["prec"] = max(maxdev["prec"], abs(
                x["precision_flux"] - (x["flux_gene_backed"] - x["flux_absent"]) /
                x["flux_gene_backed"]))
        if x["factors"] != dict(zip("PCGX", a.split("/"))):
            bad["factor_label_mismatch"] += 1
        if x["taf_face_min"] > x["taf_face_max"] + 1e-12:
            bad["min_gt_max"] += 1
out["consistency"] = {"bad": dict(bad), "maxdev": dict(maxdev)}
# D arm: two_tier/lin/G0/X100 gene-backed flux == D
dd = max(abs(r["arms"]["two_tier/lin/G0/X100"]["flux_gene_backed"] - r["model"]["D_flux_gb_two_tier"])
         for r in rows)
out["consistency"]["D_equals_tt_lin_G0_X100_flux_gb_maxdiff"] = dd

json.dump(out, open(os.path.join(L.HERE, "out_c1_c6.json"), "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str)[:12000])
