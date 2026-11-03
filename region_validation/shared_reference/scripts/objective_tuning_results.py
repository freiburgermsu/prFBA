#!/usr/bin/env python
"""Stage O selection (tune half) and confirmation (confirm half), PFUNC_ESTIMATORS_AND_
OBJECTIVE_TUNING.md section "Stage O".  Reads <outdir>/{per_amplicon.jsonl, meta.json}
written by scripts/objective_tuning.py; no solver.

SELECTION (tune half; rows whose holdout group hashes to the tune half)
  Candidates: P in {kn, est} x C in {lin, log} x G in {G0, Gmeas, G0.2} x X in {X100, X10,
  X1} = 36 configurations (est_shuffled, two_tier and oracle are controls).  Shipped:
  kn/lin/G0/X100.
  1. Eligibility (guards): a configuration is eligible iff its median total exchange flux
     and its median gene-backed flux (as returned; medians over all tune-half rows, n1
     included) are each within 10% of the shipped arm's median: |median_c / median_shipped
     - 1| <= 0.10.  Flux precision, gapfilled flux, the internal share, the biomass floor
     (bio1 = 0.95 x max on every arm) and solver status are reported beside it.
  2. Level: the pooled TAF face midpoint (taf_face_mid = truth-absent gene-backed flux on
     the optimal face, mid of min and max, / D), pooled
       A: n >= 2 strata weighted by eligible source genomes (n2 0.6975, n3 0.1869, n4+
          0.1156: 2,474 / 663 / 410, recomputed from population.json), the primary and the
          selection criterion;
       B: all answered amplicons, n1 at 0.868 and the n >= 2 strata sharing 0.132 in A's
          proportions; reported, never used to select.
  3. The G factor is read only on rows where gapfilled flux varies across G's levels:
     within each (P, C, X) family the G level is chosen by pooled A over the n >= 2 rows
     whose gapfilled-flux face midpoint differs across the three G levels by more than
     1e-6 x max(1, D) (stratum weights renormalised over the strata that have such rows),
     ties within 0.002 broken toward fewer departures from the shipped configuration (G0).
     A family with no such row keeps G0 (G is not read).  Implemented here because the
     design's sentence fixes where G is read, not how; the joint selection over all 36
     configurations on all rows is reported as a sensitivity (`naive_selection`).
  4. Across the families' representatives: the lowest pooled A; ties within 0.002 of the
     lowest -> the fewest departures from the shipped configuration (factors P, C, G, X),
     then the lowest pooled A, then factor order.
  5. Oracle headroom captured: (sel - two_tier) / (oracle - two_tier), the two controls at
     the selected (C, G, X) cell; the shipped-cell version (two_tier and oracle at
     lin/G0/X100) beside it.

CONFIRMATION (confirm half)
  Paired by amplicon: selected - shipped; selected - its own two_tier (same C, G, X);
  selected - est_shuffled at the same cell (its own shuffled arm when P = est; when P = kn
  there is no kn-shuffled arm in the design, so est_shuffled is shown with a flag);
  oracle - two_tier at the selected cell; and, as context, shipped - two_tier and
  oracle - two_tier at the shipped cell (the weighted-objective test's H1 and C3 on this
  half).  TAF face min / mid / max (+ as returned, and the guard metrics), per stratum and
  pooled A and B; 95% percentile cluster bootstrap over holdout groups, 2,000 draws, seed
  1729 + 1000 x stratum index (n2 0, n3 1, n4+ 2, n1 3), weighted_objective_test.
  cluster_boot_mean (score_folds' draw); pooled CIs combine the strata's draws with the
  weights.  Labels against the SESOI 0.01 (the weighted-objective test's H1 rules, without
  the specificity clause).  Degeneracy (FVA counts on the 1e-6 and 1e-7 faces) when a
  second pass with --fva-arms has been run (--fva-dir), else "not run".

OUTPUTS
  <outdir>/summary.json (or --summary-out) and the filled draft of OBJECTIVE_TUNING_
  RESULTS.md (template: --md-template, default the repository's OBJECTIVE_TUNING_
  RESULTS.md; output --md-out, default <outdir>/OBJECTIVE_TUNING_RESULTS.filled.md).  The
  template itself is never overwritten.
  --plumbing-all-halves (smoke only) treats every row as both tune and confirm, so the
  whole code path runs on a handful of rows; its output is labelled PLUMBING.

Run: ~/Documents/py_venv/bin/python scripts/objective_tuning_results.py [--outdir DIR]
"""
import argparse
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)
D = os.path.join(SR, "data")
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import weighted_objective_test as W   # noqa: E402  (cluster_boot_mean; no cobra at import)
import objective_tuning as OT          # noqa: E402  (arm names and factors; no cobra)

STRATA = ("n2", "n3", "n4+", "n1")     # the bootstrap seed index order of both harnesses
MULTI = ("n2", "n3", "n4+")
DECLARED_W_A = {"n2": 0.6975, "n3": 0.1869, "n4+": 0.1156}
N1_SHARE_B = 0.868
SESOI = 0.01
TIE = 0.002
GUARD = 0.10
GAPF_TOL = 1e-6
FLOOR = 0.95
B = 2000
SEED = 1729
FACE = ("taf_face_mid", "taf_face_min", "taf_face_max")
GUARD_METRICS = ("flux_exchange", "flux_gene_backed")
REPORT_GUARDS = ("flux_exchange", "flux_gene_backed", "precision_flux", "flux_gapfilled",
                 "internal_share")
CONFIRM_METRICS = FACE + ("taf_gb", "taf_point_D") + REPORT_GUARDS
FVA_METRICS = ("faces.opt.n_var_1e-3", "faces.opt7.n_var_1e-3")
CANDIDATES = [OT.arm_name(P, C, G, X) for P in OT.P_CANDIDATES for C in OT.C_SHAPES
              for G in OT.G_LEVELS for X in OT.X_LEVELS]
ORDER = {a: i for i, a in enumerate(OT.all_arms())}


# ============================================================ helpers
def val(row, arm, met):
    a = (row.get("arms") or {}).get(arm) or {}
    if a.get("status") != "optimal":
        return None
    if met.startswith("faces."):
        _f, tag, name = met.split(".", 2)
        x = ((a.get("faces") or {}).get(tag) or {}).get(name)
    else:
        x = a.get(met)
    return None if (x is None or x != x) else float(x)


def gapf_mid(row, arm):
    lo, hi = val(row, arm, "gapfilled_flux_face_min"), val(row, arm, "gapfilled_flux_face_max")
    return None if lo is None or hi is None else (lo + hi) / 2.0


def departures(name):
    f, s = OT.arm_factors(name), OT.arm_factors(OT.SHIPPED)
    return sum(f[k] != s[k] for k in "PCGX")


def pct(b):
    return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


def by_stratum(rows):
    out = defaultdict(list)
    for r in rows:
        out[r["stratum"]].append(r)
    return out


def pooled_level(by, arm, met, weights, renorm=False):
    """sum_s w_s mean_s; None if a weighted stratum has no value (unless renorm: the
    weights are then renormalised over the strata that have values)."""
    per = {}
    for s in weights:
        xs = [val(r, arm, met) for r in by.get(s, [])]
        xs = [x for x in xs if x is not None]
        if xs:
            per[s] = {"mean": float(np.mean(xs)), "n": len(xs)}
    if renorm:
        tw = sum(weights[s] for s in per)
        if not per or tw <= 0:
            return None, per
        return sum(weights[s] / tw * per[s]["mean"] for s in per), per
    if not all(s in per for s in weights):
        return None, per
    return sum(weights[s] * per[s]["mean"] for s in weights), per


def contrast(by, a, b, met, weights_by_name):
    """Paired a - b per stratum (cluster bootstrap over holdout groups, the stratum's own
    seed) and pooled under each weighting (the strata's draws combined with the weights)."""
    per, boots = {}, {}
    for s in STRATA:
        pr = [(val(r, a, met), val(r, b, met), r) for r in by.get(s, [])]
        pr = [(x, y, r) for x, y, r in pr if x is not None and y is not None]
        if not pr:
            continue
        d = [x - y for x, y, _ in pr]
        est, bg = W.cluster_boot_mean(d, [r["group_id"] for *_, r in pr], B,
                                      SEED + 1000 * STRATA.index(s))
        per[s] = {"n": len(d), "n_groups": len({r["group_id"] for *_, r in pr}),
                  "mean": est, "ci": pct(bg), "median": float(np.median(d)),
                  "n_lt0": sum(1 for x in d if x < 0), "n_eq0": sum(1 for x in d if x == 0),
                  "n_gt0": sum(1 for x in d if x > 0),
                  "arm_mean": float(np.mean([x for x, _, _ in pr])),
                  "comparator_mean": float(np.mean([y for _, y, _ in pr]))}
        boots[s] = bg
    pooled = {}
    for name, wts in weights_by_name.items():
        miss = [s for s in wts if s not in boots]
        if miss:
            pooled[name] = {"not_evaluable": "no rows in stratum {}".format(",".join(miss))}
            continue
        pooled[name] = {"mean": sum(wts[s] * per[s]["mean"] for s in wts),
                        "ci": pct(sum(wts[s] * boots[s] for s in wts)),
                        "arm_mean": sum(wts[s] * per[s]["arm_mean"] for s in wts),
                        "comparator_mean": sum(wts[s] * per[s]["comparator_mean"] for s in wts),
                        "n": sum(per[s]["n"] for s in wts)}
    return {"per_stratum": per, "pooled": pooled}


def sesoi_label(c):
    """The weighted-objective test's H1 rules (section 10) without the specificity clause,
    on one contrast's pooled face mid / min / max."""
    try:
        m, ci = c["taf_face_mid"]["mean"], c["taf_face_mid"]["ci"]
        lo, hi = c["taf_face_min"]["ci"], c["taf_face_max"]["ci"]
    except (KeyError, TypeError):
        return "not evaluable"
    ends_below = lo[1] < 0 and hi[1] < 0
    ends_above = lo[0] > 0 and hi[0] > 0
    ends_inside = all(-SESOI < x[0] and x[1] < SESOI for x in (lo, hi))
    if ci[1] < 0 and m <= -SESOI and ends_below:
        return ("lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= "
                "-0.01, and both face endpoints' CIs below 0")
    if ci[0] > 0 and m >= SESOI and ends_above:
        return ("raises truth-absent flux by >= SESOI (harmful): midpoint CI above 0 with "
                "mean >= +0.01, and both endpoints' CIs above 0")
    if -SESOI < ci[0] and ci[1] < SESOI and ends_inside:
        return ("equivalent within SESOI: the midpoint and both endpoint CIs lie inside "
                "+-0.01" + ("" if ci[0] <= 0 <= ci[1] else ", nonzero (midpoint CI excludes 0)"))
    if ci[1] < 0 and m <= -SESOI:
        return "lowers by >= SESOI at the face midpoint, but not at both face endpoints"
    if ci[0] > 0 and m >= SESOI:
        return "raises by >= SESOI at the face midpoint, but not at both face endpoints"
    if ci[1] < 0 or ci[0] > 0:
        return ("detectable but below SESOI: the midpoint CI excludes 0, |mean| < 0.01, and "
                "the CIs do not lie inside +-0.01")
    return "inconclusive: the midpoint CI crosses 0 and a SESOI bound"


# ============================================================ health
def health(rows, meta):
    ok = [r for r in rows if r["status"] == "ok"]
    h = {"n_rows": len(rows), "status": dict(Counter(r["status"] for r in rows))}
    pm = [r.get("prior_match") or {} for r in ok]
    h["prior_model_match"] = {"n": len(pm), "ok": sum(1 for x in pm if x.get("ok")),
                              "keys_failing": [r["key"] for r, x in zip(ok, pm)
                                               if not x.get("ok")],
                              "extra_mismatch_rows": sum(1 for x in pm
                                                         if x.get("extra_mismatches"))}

    def bitwise_count(name, k):
        xs = [((x.get(name) or {}).get(k) or {}) for x in pm]
        xs = [y for y in xs if y]
        return {"n": len(xs), "bitwise": sum(1 for y in xs if y.get("bitwise")),
                "max_abs_diff": max((y["abs_diff"] for y in xs if y.get("abs_diff")
                                     is not None), default=None)}
    for name, lab in (("shipped_vs_prior_real", "shipped_vs_prior_real"),
                      ("two_tier_vs_prior_two_tier", "two_tier_vs_prior_two_tier")):
        h[lab] = {k: bitwise_count(name, k) for k in OT.PRIOR_ARM_KEYS}
        h[lab]["taf_face_within_1e-9"] = sum(
            1 for x in pm if all(((x.get(name) or {}).get(k) or {}).get("abs_diff") is not None
                                 and x[name][k]["abs_diff"] <= 1e-9 for k in FACE))
    Dm = [x.get("D") or {} for x in pm]
    h["D_vs_prior"] = {"n": sum(1 for d in Dm if d.get("prior") is not None),
                       "bitwise": sum(1 for d in Dm if d.get("bitwise")),
                       "max_abs_diff": max((d["abs_diff"] for d in Dm
                                            if d.get("abs_diff") is not None), default=None)}
    rc = [r.get("released_check") or {} for r in ok if r.get("released_check")]
    h["released_check"] = {"n": len(rc), "pass": sum(1 for x in rc if x.get("pass")),
                           "equals_prior_released_bitwise": sum(
                               1 for x in rc if x.get("equals_prior_released_bitwise")),
                           "max_rel_diff": max((x["rel_diff"] for x in rc), default=None),
                           "released_source_sha256": sorted({x.get("released_source_sha256")
                                                             for x in rc})}
    sc = [r.get("shipped_check") or {} for r in ok if r.get("shipped_check")]
    h["shipped_coefficients_equal_released_copy"] = {
        "n": len(sc), "coef_map_identical": sum(1 for x in sc if x.get("coef_map_identical")),
        "objective_bitwise": sum(1 for x in sc if x.get("objective_bitwise_equal")),
        "probs_equal_production_rule": sum(1 for x in sc if x.get("probs_equal_production_rule"))}
    dt = [r.get("determinism") or {} for r in ok]
    h["determinism_pairs"] = {"n_pairs": sum(x.get("n_pairs", 0) for x in dt),
                              "identical": sum(x.get("n_identical", 0) for x in dt),
                              "max_abs_diff": max((x.get("max_abs_diff") or 0.0 for x in dt),
                                                  default=None),
                              "differing": [d for x in dt for d in x.get("differing", [])][:20]}
    rp = [r.get("repro") for r in ok if r.get("repro")]
    h["repro_check"] = ({"n": len(rp), "all_bitwise": sum(1 for x in rp if x.get("all_bitwise")),
                         "keys_not_bitwise": [r["key"] for r in ok if r.get("repro") and
                                              not r["repro"].get("all_bitwise")]}
                        if rp else "not run")
    arms = [(r, a, x) for r in ok for a, x in (r.get("arms") or {}).items()]
    h["arm_solves"] = {"n": len(arms), "optimal": sum(1 for _, _, x in arms
                                                      if x.get("status") == "optimal")}
    br = [x["biomass_ratio"] for _, _, x in arms if x.get("biomass_ratio") is not None]
    h["biomass_ratio"] = {"min": min(br, default=None), "max": max(br, default=None),
                          "n_off_floor_gt_1e-6": sum(1 for x in br if abs(x - FLOOR) > 1e-6)}
    mb = []
    for r in ok:
        xs = [x.get("max_biomass") for x in r["arms"].values() if x.get("max_biomass")]
        if xs:
            mb.append((max(xs) - min(xs)) / max(xs))
    h["max_biomass_rel_spread_across_arms_max"] = max(mb, default=None)
    h["face_lp"] = {"failed": sum(x.get("face_lp_failed", 0) for _, _, x in arms),
                    "retried": sum(x.get("face_lp_retried", 0) for _, _, x in arms)}
    h["gapfilled_and_gene_backed"] = sum(r["model"].get("n_gapfilled_and_gene_backed", 0)
                                         for r in ok)
    h["genes_not_template_roles"] = sum(len(r["model"].get("genes_not_template_roles") or [])
                                        for r in ok)
    h["est"] = {"alpha_with_null_clade_rows": sum(1 for r in ok if (r.get("est_info") or {})
                                                  .get("alpha_with_null_clade")),
                "auroc_mean": {p: _mean([((r.get("est_info") or {}).get(
                    "auroc_supported_vs_absent") or {}).get(p) for r in ok]) for p in ("kn", "est")},
                "brier_mean": {p: _mean([((r.get("est_info") or {}).get(
                    "brier_vs_truth_support") or {}).get(p) for r in ok]) for p in ("kn", "est")}}
    w = [r["wall_seconds"] for r in ok]
    h["time"] = {"wall_mean_s": _mean(w), "wall_max_s": max(w, default=None),
                 "gapfill_mean_s": _mean([r["timing"].get("gapfill") for r in ok]),
                 "arms_mean_s": _mean([r["timing"].get("arms") for r in ok]),
                 "worker_maxrss_gb_max": max((r.get("worker_maxrss_gb") or 0 for r in ok),
                                             default=None)}
    return h


def _mean(xs):
    xs = [float(x) for x in xs if x is not None]
    return (sum(xs) / len(xs)) if xs else None


# ============================================================ selection
def select(rows, WA, WB):
    by = by_stratum(rows)
    arms_present = sorted({a for r in rows for a in (r.get("arms") or {})}, key=ORDER.get)
    out = {"n_rows": len(rows), "n_by_stratum": {s: len(by.get(s, [])) for s in STRATA},
           "n_groups_by_stratum": {s: len({r["group_id"] for r in by.get(s, [])})
                                   for s in STRATA}}
    # ---- 1. guards (medians over all tune rows, n1 included)
    med = {}
    for a in arms_present:
        med[a] = {}
        for m in REPORT_GUARDS:
            xs = [val(r, a, m) for r in rows]
            xs = [x for x in xs if x is not None]
            med[a][m] = float(np.median(xs)) if xs else None
        med[a]["n_optimal"] = sum(1 for r in rows if val(r, a, "objective") is not None)
        br = [val(r, a, "biomass_ratio") for r in rows]
        med[a]["floor_met_all"] = all(x is not None and abs(x - FLOOR) <= 1e-6 for x in br)
    ship = med.get(OT.SHIPPED) or {}
    guard = {}
    for c in CANDIDATES:
        if c not in med:
            continue
        ratios = {m: (med[c][m] / ship[m]) if (med[c].get(m) is not None and ship.get(m))
                  else None for m in REPORT_GUARDS}
        row_ratio = {}
        for m in GUARD_METRICS:
            rr = [val(r, c, m) / val(r, OT.SHIPPED, m) for r in rows
                  if val(r, c, m) is not None and val(r, OT.SHIPPED, m)]
            row_ratio[m] = float(np.median(rr)) if rr else None
        elig = (all(ratios[m] is not None and abs(ratios[m] - 1.0) <= GUARD
                    for m in GUARD_METRICS) and med[c]["n_optimal"] == len(rows)
                and med[c]["floor_met_all"])
        guard[c] = {"median": {m: med[c][m] for m in REPORT_GUARDS},
                    "ratio_of_medians_to_shipped": ratios,
                    "median_of_row_ratios_to_shipped": row_ratio,
                    "n_optimal": med[c]["n_optimal"], "floor_met_all": med[c]["floor_met_all"],
                    "eligible": bool(elig)}
    E = [c for c in CANDIDATES if guard.get(c, {}).get("eligible")]
    out["guards"] = guard
    out["eligible"] = E
    out["n_eligible"] = len(E)
    out["eligible_by_X"] = {X: sum(1 for c in E if OT.arm_factors(c)["X"] == X)
                            for X in OT.X_LEVELS}
    out["control_medians"] = {a: med[a] for a in arms_present if a not in CANDIDATES}
    # ---- 2. levels of every arm (controls too)
    levels = {}
    for a in arms_present:
        la, pa = pooled_level(by, a, "taf_face_mid", WA)
        lb, pb = pooled_level(by, a, "taf_face_mid", WB)
        levels[a] = {"A": la, "B": lb, "per_stratum": {**pa, **pb},
                     "A_min": pooled_level(by, a, "taf_face_min", WA)[0],
                     "A_max": pooled_level(by, a, "taf_face_max", WA)[0]}
    out["levels"] = levels
    # ---- 3. the G factor, read per (P, C, X) family on G-varying rows
    fams = {}
    for P in OT.P_CANDIDATES:
        for C in OT.C_SHAPES:
            for X in OT.X_LEVELS:
                fam = [OT.arm_name(P, C, G, X) for G in OT.G_LEVELS]
                fam_e = [a for a in fam if a in E]
                if not fam_e:
                    continue
                vary = []
                for r in rows:
                    if r["stratum"] not in MULTI:
                        continue
                    mids = [gapf_mid(r, a) for a in fam]
                    mids = [x for x in mids if x is not None]
                    Dv = (r.get("model") or {}).get("D_flux_gb_two_tier") or 0.0
                    if len(mids) >= 2 and max(mids) - min(mids) > GAPF_TOL * max(1.0, Dv):
                        vary.append(r)
                vby = by_stratum(vary)
                rec = {"family": "{}/{}/*/{}".format(P, C, X), "eligible_G": fam_e,
                       "n_vary_rows": len(vary),
                       "n_vary_by_stratum": {s: len(vby.get(s, [])) for s in MULTI}}
                if len(fam_e) == 1:
                    rec["choice"], rec["why"] = fam_e[0], "only one eligible G level"
                elif not vary:
                    g0 = OT.arm_name(P, C, "G0", X)
                    if g0 in fam_e:
                        rec["choice"], rec["why"] = g0, "gapfilled flux varies on no row: G not read"
                    else:
                        vals = {a: levels[a]["A"] for a in fam_e if levels[a]["A"] is not None}
                        rec["choice"] = (min(vals, key=lambda a: (vals[a], ORDER[a]))
                                         if vals else fam_e[0])
                        rec["why"] = ("gapfilled flux varies on no row and G0 is ineligible: "
                                      "the lowest all-row pooled A among eligible G levels")
                else:
                    vals = {}
                    for a in fam_e:
                        v_, per = pooled_level(vby, a, "taf_face_mid", WA, renorm=True)
                        if v_ is not None:
                            vals[a] = v_
                    rec["pooled_A_on_vary_rows"] = vals
                    if vals:
                        best = min(vals.values())
                        ties = [a for a in vals if vals[a] <= best + TIE]
                        rec["tie_set"] = ties
                        rec["choice"] = min(ties, key=lambda a: (departures(a), vals[a], ORDER[a]))
                        rec["why"] = ("lowest pooled A on the G-varying rows; ties within {} -> "
                                      "fewer departures".format(TIE))
                    else:
                        rec["choice"], rec["why"] = fam_e[0], "no evaluable G level on vary rows"
                fams[rec["family"]] = rec
    out["g_reading"] = fams
    reps = [f["choice"] for f in fams.values()]
    # ---- 4. across families (all tune rows of the n >= 2 strata, pooling A)
    A = {a: levels[a]["A"] for a in reps if levels[a]["A"] is not None}
    out["representatives"] = {a: A.get(a) for a in reps}
    if A:
        best = min(A.values())
        ties = [a for a in A if A[a] <= best + TIE]
        sel = min(ties, key=lambda a: (departures(a), A[a], ORDER[a]))
        out["selection"] = {"selected": sel, "factors": OT.arm_factors(sel),
                            "pooled_A": A[sel], "pooled_B": levels[sel]["B"],
                            "lowest_pooled_A": best,
                            "lowest_arm": min(A, key=lambda a: (A[a], ORDER[a])),
                            "tie_set": sorted(ties, key=ORDER.get),
                            "departures_from_shipped": departures(sel),
                            "shipped_pooled_A": levels.get(OT.SHIPPED, {}).get("A"),
                            "shipped_pooled_B": levels.get(OT.SHIPPED, {}).get("B")}
    else:
        out["selection"] = {"selected": None, "why": "no eligible configuration has a pooled "
                                                     "A level (a stratum without rows?)"}
    # ---- sensitivity: the joint selection over every eligible configuration, all rows
    AN = {a: levels[a]["A"] for a in E if levels[a]["A"] is not None}
    if AN:
        best = min(AN.values())
        ties = [a for a in AN if AN[a] <= best + TIE]
        naive = min(ties, key=lambda a: (departures(a), AN[a], ORDER[a]))
        out["naive_selection"] = {"selected": naive, "tie_set": sorted(ties, key=ORDER.get),
                                  "agrees": naive == out["selection"].get("selected")}
    # pooling B ranking (reported, never used)
    AB = {a: levels[a]["B"] for a in E if levels[a]["B"] is not None}
    out["ranking_B"] = sorted(AB, key=lambda a: (AB[a], ORDER[a]))[:10]
    out["ranking_A"] = sorted(AN, key=lambda a: (AN[a], ORDER[a]))[:10]
    # ---- 5. headroom on the tune half
    sel = out["selection"].get("selected")
    if sel:
        out["headroom_tune"] = headroom(levels, sel)
        out["contrast_vs_shipped_tune"] = {
            met: contrast(by, sel, OT.SHIPPED, met, {"A": WA, "B": WB})["pooled"]
            for met in FACE}
    return out


def headroom(levels, sel):
    f = OT.arm_factors(sel)
    tt = OT.arm_name("two_tier", f["C"], f["G"], f["X"])
    orc = OT.arm_name("oracle", f["C"], f["G"], f["X"])
    out = {}
    for pool in ("A", "B"):
        def g(a):
            return (levels.get(a) or {}).get(pool)
        vals = {k: g(a) for k, a in (("sel", sel), ("tt", tt), ("orc", orc),
                                     ("ship", OT.SHIPPED), ("tt0", OT.D_ARM),
                                     ("orc0", OT.arm_name("oracle", "lin", "G0", "X100")))}
        if any(v is None for v in vals.values()):
            out[pool] = "not evaluable"
            continue
        own = (vals["sel"] - vals["tt"]) / (vals["orc"] - vals["tt"]) \
            if vals["orc"] != vals["tt"] else None
        frame0 = (vals["sel"] - vals["tt0"]) / (vals["orc0"] - vals["tt0"]) \
            if vals["orc0"] != vals["tt0"] else None
        ship = (vals["ship"] - vals["tt0"]) / (vals["orc0"] - vals["tt0"]) \
            if vals["orc0"] != vals["tt0"] else None
        out[pool] = {"captured_own_cell": own, "captured_shipped_cell_frame": frame0,
                     "shipped_captured": ship, "levels": vals}
    return out


# ============================================================ confirmation
def confirm(rows, sel, WA, WB, fva_rows=None):
    by = by_stratum(rows)
    f = OT.arm_factors(sel)
    tt = OT.arm_name("two_tier", f["C"], f["G"], f["X"])
    orc = OT.arm_name("oracle", f["C"], f["G"], f["X"])
    shuf = OT.arm_name("est_shuffled", f["C"], f["G"], f["X"])
    tt0, orc0 = OT.D_ARM, OT.arm_name("oracle", "lin", "G0", "X100")
    pairs = {"selected-shipped": (sel, OT.SHIPPED), "selected-two_tier": (sel, tt),
             "selected-shuffled": (sel, shuf), "oracle-two_tier@selected_cell": (orc, tt),
             "shipped-two_tier@shipped_cell": (OT.SHIPPED, tt0),
             "oracle-two_tier@shipped_cell": (orc0, tt0)}
    wts = {"A": WA, "B": WB}
    out = {"selected": sel, "n_rows": len(rows),
           "n_by_stratum": {s: len(by.get(s, [])) for s in STRATA},
           "n_groups_by_stratum": {s: len({r["group_id"] for r in by.get(s, [])})
                                   for s in STRATA},
           "shuffled_comparator": shuf,
           "shuffled_is_own": f["P"] == "est",
           "shuffled_note": ("" if f["P"] == "est" else
                             "P = {} has no own shuffled arm in the design; est_shuffled at "
                             "the same cell is shown, and the weighted-objective test's "
                             "real - shuffled (kn, lin/G0/X100: -0.0120 [-0.0179, -0.0069] "
                             "pooled) is the only kn specificity contrast".format(f["P"])),
           "contrasts": {}}
    for name, (a, b) in pairs.items():
        out["contrasts"][name] = {"arms": [a, b], "metrics": {
            met: contrast(by, a, b, met, wts) for met in CONFIRM_METRICS}}
    labels = {}
    for name in pairs:
        for pool in ("A", "B"):
            c = {met: out["contrasts"][name]["metrics"][met]["pooled"].get(pool) or {}
                 for met in FACE}
            labels["{}|{}".format(name, pool)] = (
                sesoi_label(c) if all("ci" in c[m] for m in FACE) else "not evaluable")
    out["labels"] = labels
    # headroom on the confirm half (point estimates of pooled means)
    hr = {}
    for pool in ("A", "B"):
        s_t = out["contrasts"]["selected-two_tier"]["metrics"]["taf_face_mid"]["pooled"].get(pool) or {}
        o_t = out["contrasts"]["oracle-two_tier@selected_cell"]["metrics"]["taf_face_mid"]["pooled"].get(pool) or {}
        sh_t = out["contrasts"]["shipped-two_tier@shipped_cell"]["metrics"]["taf_face_mid"]["pooled"].get(pool) or {}
        o0 = out["contrasts"]["oracle-two_tier@shipped_cell"]["metrics"]["taf_face_mid"]["pooled"].get(pool) or {}
        hr[pool] = {"selected_own_cell": (s_t["mean"] / o_t["mean"]) if s_t.get("mean") is not None
                    and o_t.get("mean") else None,
                    "shipped": (sh_t["mean"] / o0["mean"]) if sh_t.get("mean") is not None
                    and o0.get("mean") else None}
    out["headroom_confirm"] = hr
    # degeneracy: a second pass with --fva-arms (finalists only)
    src = fva_rows if fva_rows is not None else rows
    fby = by_stratum(src)
    have = any(val(r, sel, FVA_METRICS[0]) is not None and val(r, OT.SHIPPED, FVA_METRICS[0])
               is not None for r in src)
    if have:
        out["degeneracy"] = {name: {met: contrast(fby, a, b, met, wts) for met in FVA_METRICS}
                             for name, (a, b) in pairs.items()
                             if any(val(r, a, FVA_METRICS[0]) is not None and
                                    val(r, b, FVA_METRICS[0]) is not None for r in src)}
    else:
        out["degeneracy"] = "not run (second pass: objective_tuning.py --half confirm --arms "\
                            "<finalists> --fva-arms <finalists> --outdir <new dir>; then "\
                            "--fva-dir here)"
    return out


# ============================================================ rendering
def f4(x, nd=4):
    if x is None:
        return "n/a"
    if isinstance(x, str):
        return x
    return "{:+.{}f}".format(x, nd) if nd else str(x)


def fl(x, nd=4):
    return "n/a" if x is None else ("{:.{}f}".format(x, nd) if not isinstance(x, str) else x)


def ci_str(c):
    if not c or "ci" not in c:
        return (c or {}).get("not_evaluable", "n/a") if isinstance(c, dict) else "n/a"
    return "{} [{}, {}]".format(f4(c["mean"]), f4(c["ci"][0]), f4(c["ci"][1]))


def render(template, summ):
    """Fill {{placeholders}} of the results template from the summary."""
    h, t, c = summ["checks"], summ.get("tune") or {}, summ.get("confirm") or {}
    m = summ["meta"]
    vals = {}
    vals["status_banner"] = ("**PLUMBING / SMOKE OUTPUT: NOT RESULTS.** " if m.get("plumbing")
                             or m.get("weights_stub") else "") + (
        "Weights file is a STUB (meta.stub), not Stage E's estimator. " if m.get("weights_stub")
        else "") + ("Every row was treated as both tune and confirm "
                    "(--plumbing-all-halves). " if m.get("plumbing") else "")
    vals["generated"] = m["generated"]
    vals["n_rows"] = str(h["n_rows"])
    vals["status_counts"] = json.dumps(h["status"])
    vals["halves"] = json.dumps(m["n_by_half"])
    vals["weights_path"] = str(m["weights"]["path"])
    vals["weights_sha256"] = str(m["weights"]["sha256"])[:12]
    vals["estimator"] = str((m["weights"].get("meta") or {}).get("estimator", "see weights meta"))
    vals["g_meas"] = "{:.6f}".format(m["g_meas"])
    vals["g_meas_info"] = "{} / {} over {} tune-half rows of the prior run".format(
        m["g_meas_info"]["sum_n_gapfilled_supported"], m["g_meas_info"]["sum_n_gapfilled"],
        m["g_meas_info"]["n_rows"])
    vals["weights_A"] = ", ".join("{} {:.4f}".format(s, w) for s, w in m["weights_A"].items())
    vals["weights_B"] = ", ".join("{} {:.4f}".format(s, w) for s, w in m["weights_B"].items())
    # checks table
    pmm = h["prior_model_match"]
    sv, tv = h["shipped_vs_prior_real"], h["two_tier_vs_prior_two_tier"]
    rc, sc, dp = h["released_check"], h["shipped_coefficients_equal_released_copy"], h["determinism_pairs"]
    rp = h["repro_check"]
    lines = ["| check | result |", "|---|---|",
             "| rows ok / total | {} / {} |".format(h["status"].get("ok", 0), h["n_rows"]),
             "| rebuilt model = prior row (n_new_rxn, growth 1e-6, n_gene_backed, n_gapfilled) | "
             "{} / {} |".format(pmm["ok"], pmm["n"]),
             "| D (two_tier/lin/G0/X100 gene-backed flux) bitwise = prior | {} / {} |".format(
                 h["D_vs_prior"]["bitwise"], h["D_vs_prior"]["n"]),
             "| kn/lin/G0/X100 objective bitwise = prior `real` | {} / {} |".format(
                 sv["objective"]["bitwise"], sv["objective"]["n"]),
             "| kn/lin/G0/X100 returned vertex (gene-backed, truth-absent flux) bitwise = prior | "
             "{} / {}, {} / {} |".format(sv["flux_gene_backed"]["bitwise"],
                                         sv["flux_gene_backed"]["n"],
                                         sv["flux_absent"]["bitwise"], sv["flux_absent"]["n"]),
             "| kn/lin/G0/X100 face min/mid/max within 1e-9 of prior `real` | {} / {} (max abs "
             "diff mid {}) |".format(sv["taf_face_within_1e-9"], sv["taf_face_mid"]["n"],
                                     sv["taf_face_mid"]["max_abs_diff"]),
             "| two_tier/lin/G0/X100 objective bitwise = prior `two_tier` | {} / {} |".format(
                 tv["objective"]["bitwise"], tv["objective"]["n"]),
             "| two_tier/lin/G0/X100 face within 1e-9 of prior | {} / {} (max abs diff mid {}) |"
             .format(tv["taf_face_within_1e-9"], tv["taf_face_mid"]["n"],
                     tv["taf_face_mid"]["max_abs_diff"]),
             "| coefficient map of the new function = prfba_released_copy (shipped arm) | "
             "{} / {} |".format(sc["coef_map_identical"], sc["n"]),
             "| released_check pass; = prior objective_released bitwise | {} / {}; {} / {} |"
             .format(rc["pass"], rc["n"], rc["equals_prior_released_bitwise"], rc["n"]),
             "| determinism pairs identical (same LP) | {} / {} |".format(dp["identical"],
                                                                       dp["n_pairs"]),
             "| prior arms re-run through the prior code, bitwise (repro) | {} |".format(
                 "{} / {}".format(rp["all_bitwise"], rp["n"]) if isinstance(rp, dict) else rp),
             "| arm solves optimal | {} / {} |".format(h["arm_solves"]["optimal"],
                                                       h["arm_solves"]["n"]),
             "| biomass at 0.95 x max (min, max ratio; off floor > 1e-6) | {}, {}; {} |".format(
                 h["biomass_ratio"]["min"], h["biomass_ratio"]["max"],
                 h["biomass_ratio"]["n_off_floor_gt_1e-6"]),
             "| face LP failures (after one cold retry); retried | {}; {} |".format(
                 h["face_lp"]["failed"], h["face_lp"]["retried"]),
             "| gapfilled reactions that are gene-backed | {} |".format(
                 h["gapfilled_and_gene_backed"]),
             "| rows with alpha > 0 and a null clade (literal formula) | {} |".format(
                 h["est"]["alpha_with_null_clade_rows"]),
             "| wall mean / max per amplicon; largest worker RSS | {:.1f} s / {:.1f} s; {} GiB |"
             .format(h["time"]["wall_mean_s"] or 0, h["time"]["wall_max_s"] or 0,
                     h["time"]["worker_maxrss_gb_max"])]
    vals["checks_table"] = "\n".join(lines)
    # tune
    vals["tune_n"] = json.dumps(t.get("n_by_stratum", {}))
    vals["n_eligible"] = str(t.get("n_eligible", "n/a"))
    vals["eligible_by_X"] = json.dumps(t.get("eligible_by_X", {}))
    gl = ["| configuration | exchange flux ratio | gene-backed flux ratio | precision ratio | "
          "eligible | pooled A | pooled B |", "|---|---|---|---|---|---|---|"]
    lev = t.get("levels", {})
    for cfg in CANDIDATES:
        g = (t.get("guards") or {}).get(cfg)
        if not g:
            continue
        r_ = g["ratio_of_medians_to_shipped"]
        gl.append("| {} | {} | {} | {} | {} | {} | {} |".format(
            cfg, fl(r_.get("flux_exchange"), 3), fl(r_.get("flux_gene_backed"), 3),
            fl(r_.get("precision_flux"), 3), "yes" if g["eligible"] else "no",
            fl((lev.get(cfg) or {}).get("A")), fl((lev.get(cfg) or {}).get("B"))))
    vals["guard_table"] = "\n".join(gl)
    fr = ["| family | eligible G | rows where gapfilled flux varies (n2/n3/n4+) | choice | why |",
          "|---|---|---|---|---|"]
    for fam, rec in (t.get("g_reading") or {}).items():
        fr.append("| {} | {} | {} | {} | {} |".format(
            fam, ", ".join(OT.arm_factors(a)["G"] for a in rec["eligible_G"]),
            "/".join(str(rec["n_vary_by_stratum"].get(s, 0)) for s in MULTI),
            rec.get("choice"), rec.get("why")))
    vals["g_reading_table"] = "\n".join(fr)
    sel = (t.get("selection") or {})
    vals["selected"] = str(sel.get("selected"))
    vals["selection_paragraph"] = (
        "Selected: **{}** (pooled A {} against the shipped {}; lowest pooled A {} at {}; tie "
        "set within {}: {}; {} departure(s) from the shipped configuration). The joint "
        "selection over all eligible configurations on all rows (sensitivity) picks {} ({})."
        .format(sel.get("selected"), fl(sel.get("pooled_A")), fl(sel.get("shipped_pooled_A")),
                fl(sel.get("lowest_pooled_A")), sel.get("lowest_arm"), TIE,
                ", ".join(sel.get("tie_set") or []), sel.get("departures_from_shipped"),
                (t.get("naive_selection") or {}).get("selected"),
                "agrees" if (t.get("naive_selection") or {}).get("agrees") else "DIFFERS")
        if sel.get("selected") else "No configuration selected: {}".format(sel.get("why")))
    vals["ranking_A"] = ", ".join("{} ({})".format(a, fl(lev[a]["A"])) for a in t.get("ranking_A", []))
    vals["ranking_B"] = ", ".join("{} ({})".format(a, fl(lev[a]["B"])) for a in t.get("ranking_B", []))
    ht = (t.get("headroom_tune") or {}).get("A")
    vals["headroom_tune"] = ("{} of the oracle's headroom at its own cell; {} in the shipped "
                             "cell's frame (the shipped configuration: {})".format(
                                 fl(ht["captured_own_cell"], 3), fl(ht["captured_shipped_cell_frame"], 3),
                                 fl(ht["shipped_captured"], 3))
                             if isinstance(ht, dict) else "not evaluable")
    # confirm
    cr = ["| contrast | pooling | face mid | face min | face max | as returned (taf_gb) | label |",
          "|---|---|---|---|---|---|---|"]
    for name in ("selected-shipped", "selected-two_tier", "selected-shuffled",
                 "oracle-two_tier@selected_cell", "shipped-two_tier@shipped_cell",
                 "oracle-two_tier@shipped_cell"):
        cc = (c.get("contrasts") or {}).get(name)
        if not cc:
            continue
        for pool in ("A", "B"):
            mm = cc["metrics"]
            cr.append("| {} ({} - {}) | {} | {} | {} | {} | {} | {} |".format(
                name, cc["arms"][0], cc["arms"][1], pool,
                ci_str(mm["taf_face_mid"]["pooled"].get(pool)),
                ci_str(mm["taf_face_min"]["pooled"].get(pool)),
                ci_str(mm["taf_face_max"]["pooled"].get(pool)),
                ci_str(mm["taf_gb"]["pooled"].get(pool)),
                (c.get("labels") or {}).get("{}|{}".format(name, pool), "")))
    vals["confirm_table"] = "\n".join(cr)
    vals["confirm_n"] = json.dumps(c.get("n_by_stratum", {}))
    vals["shuffled_note"] = c.get("shuffled_note") or "P = est: est_shuffled is its own shuffled arm."
    st = ["| contrast | stratum | n (groups) | face mid | face min | face max |",
          "|---|---|---|---|---|---|"]
    for name in ("selected-shipped", "selected-two_tier", "selected-shuffled"):
        cc = (c.get("contrasts") or {}).get(name)
        if not cc:
            continue
        for s in STRATA:
            p = cc["metrics"]["taf_face_mid"]["per_stratum"].get(s)
            if not p:
                continue
            st.append("| {} | {} | {} ({}) | {} | {} | {} |".format(
                name, s, p["n"], p["n_groups"], ci_str(p),
                ci_str(cc["metrics"]["taf_face_min"]["per_stratum"].get(s)),
                ci_str(cc["metrics"]["taf_face_max"]["per_stratum"].get(s))))
    vals["confirm_strata_table"] = "\n".join(st)
    gd = ["| contrast | pooling | exchange flux | gene-backed flux | flux precision | gapfilled flux |",
          "|---|---|---|---|---|---|"]
    for name in ("selected-shipped",):
        cc = (c.get("contrasts") or {}).get(name)
        if cc:
            for pool in ("A", "B"):
                gd.append("| {} | {} | {} | {} | {} | {} |".format(
                    name, pool, *[ci_str(cc["metrics"][g_]["pooled"].get(pool)) for g_ in
                                  ("flux_exchange", "flux_gene_backed", "precision_flux",
                                   "flux_gapfilled")]))
    vals["confirm_guard_table"] = "\n".join(gd)
    hc = (c.get("headroom_confirm") or {}).get("A") or {}
    vals["headroom_confirm"] = "selected (own cell) {}; shipped {}".format(
        fl(hc.get("selected_own_cell"), 3), fl(hc.get("shipped"), 3))
    dg = c.get("degeneracy")
    vals["degeneracy"] = dg if isinstance(dg, str) else "\n".join(
        "- {}: opt face {}; 1e-7 face {}".format(
            name, ci_str(d["faces.opt.n_var_1e-3"]["pooled"].get("A")),
            ci_str(d["faces.opt7.n_var_1e-3"]["pooled"].get("A"))) for name, d in dg.items())
    out = re.sub(r"<!-- TEMPLATE-BANNER-START -->.*?<!-- TEMPLATE-BANNER-END -->\n*",
                 "> Filled by `scripts/objective_tuning_results.py` from `{}` on {}. Check every "
                 "number by an independent re-derivation before quoting it.\n\n".format(
                     m["outdir"], m["generated"]), template, flags=re.S)
    for k, v in vals.items():
        out = out.replace("{{" + k + "}}", v)
    left = sorted(set(re.findall(r"\{\{([a-zA-Z0-9_]+)\}\}", out)))
    return out, left


# ============================================================ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=os.path.join(D, "objective_tuning"))
    ap.add_argument("--population", default=os.path.join(D, "weighted_objective",
                                                         "population.json"))
    ap.add_argument("--fva-dir", default=None,
                    help="outdir of a second pass with --fva-arms (confirm half, finalists)")
    ap.add_argument("--summary-out", default=None)
    ap.add_argument("--md-template", default=os.path.join(D, "objective_tuning", "OBJECTIVE_TUNING_RESULTS.template.md"))
    ap.add_argument("--md-out", default=None)
    ap.add_argument("--plumbing-all-halves", action="store_true",
                    help="SMOKE ONLY: every row is both tune and confirm")
    args = ap.parse_args(argv)
    rows = [json.loads(l) for l in open(os.path.join(args.outdir, "per_amplicon.jsonl"))
            if l.strip()]
    meta = json.load(open(os.path.join(args.outdir, "meta.json")))
    elig = json.load(open(args.population))["meta"]["eligible_source_genomes_by_stratum_excl_twin"]
    tw = sum(elig[s] for s in MULTI)
    WA = {s: elig[s] / tw for s in MULTI}
    assert all(round(WA[s], 4) == DECLARED_W_A[s] for s in MULTI), WA
    WB = {"n1": N1_SHARE_B, **{s: (1 - N1_SHARE_B) * WA[s] for s in MULTI}}
    ok = [r for r in rows if r["status"] == "ok"]
    if args.plumbing_all_halves:
        tune_rows, conf_rows = ok, ok
    else:
        tune_rows = [r for r in ok if r["half"] == "tune"]
        conf_rows = [r for r in ok if r["half"] == "confirm"]
    fva_rows = None
    if args.fva_dir:
        fr = {json.loads(l)["key"]: json.loads(l) for l in
              open(os.path.join(args.fva_dir, "per_amplicon.jsonl")) if l.strip()}
        fva_rows = [fr[r["key"]] for r in conf_rows if r["key"] in fr and
                    fr[r["key"]]["status"] == "ok"]
    summ = {"meta": {"generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                     "outdir": os.path.abspath(args.outdir),
                     "plumbing": bool(args.plumbing_all_halves),
                     "weights": meta.get("weights", {}), "weights_stub":
                     bool((meta.get("weights") or {}).get("stub")),
                     "g_meas": meta["g_meas"], "g_meas_info": meta["g_meas_info"],
                     "weights_A": WA, "weights_B": WB, "sesoi": SESOI, "tie": TIE,
                     "guard": GUARD, "gapf_tol_rel_D": GAPF_TOL, "bootstrap": B, "seed": SEED,
                     "n_by_half": dict(Counter(r["half"] for r in ok)),
                     "script": {"path": os.path.abspath(__file__),
                                "sha256": W.sha256(os.path.abspath(__file__))},
                     "harness_meta": {k: meta.get(k) for k in ("spec", "script", "git_commit",
                                                               "generated", "run")},
                     "interpretation": [
                         "eligibility = ratio of tune-half medians (all strata, n1 included) "
                         "of total exchange flux and of gene-backed flux, each within 10% of "
                         "the shipped arm's; every arm optimal and at the 0.95 floor on every row",
                         "G is read per (P, C, X) family on the n >= 2 rows where the gapfilled-"
                         "flux FACE MIDPOINT differs across G's levels by > 1e-6 x max(1, D); "
                         "stratum weights renormalised over strata with such rows; ties -> G0",
                         "selection across families on all tune rows of the n >= 2 strata, "
                         "pooling A; ties within 0.002 -> fewest departures, then lowest A",
                         "headroom captured relative to the selected cell's own two_tier and "
                         "oracle; the shipped-cell version beside it",
                         "when P = kn is selected, the shuffled comparator is est_shuffled at "
                         "the same cell, flagged (the design has no kn-shuffled arm)"]},
            "checks": health(rows, meta)}
    summ["tune"] = select(tune_rows, WA, WB) if tune_rows else "no tune rows"
    sel = (summ["tune"].get("selection") or {}).get("selected") if isinstance(
        summ["tune"], dict) else None
    summ["confirm"] = (confirm(conf_rows, sel, WA, WB, fva_rows) if (sel and conf_rows)
                       else "not evaluable (no selection or no confirm rows)")
    sp = args.summary_out or os.path.join(args.outdir, "summary.json")
    with open(sp, "w") as fh:
        json.dump(summ, fh, indent=1, sort_keys=True, default=str)
    msg = {"summary": sp}
    if os.path.exists(args.md_template):
        tpl = open(args.md_template).read()
        if isinstance(summ["confirm"], str):
            summ["confirm"] = {}
        if isinstance(summ["tune"], str):
            summ["tune"] = {}
        filled, left = render(tpl, summ)
        mo = args.md_out or os.path.join(args.outdir, "OBJECTIVE_TUNING_RESULTS.filled.md")
        if os.path.abspath(mo) == os.path.abspath(args.md_template):
            raise SystemExit("[otune-results] --md-out must not overwrite the template")
        with open(mo, "w") as fh:
            fh.write(filled)
        msg.update({"md_out": mo, "unfilled_placeholders": left})
    msg["selected"] = sel
    msg["checks"] = {k: summ["checks"][k] for k in ("status", "prior_model_match",
                                                     "released_check", "determinism_pairs",
                                                     "arm_solves", "biomass_ratio")}
    print(json.dumps(msg, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
