"""Recompute the weighted-objective test's contrasts from per_amplicon.jsonl, check them
against summary.json, re-apply the section-10 decision rules, compute the descriptive
numbers WEIGHTED_OBJECTIVE_RESULTS.md quotes, and draw its two figures.

Reads only data/weighted_objective/{per_amplicon.jsonl, population.json, summary.json}.
Writes data/weighted_objective/recomputed_contrasts.json and fig_wobj_{flux,degeneracy}.{png,pdf}.
No solver. Run: ~/Documents/py_venv/bin/python scripts/wobj_results.py
"""
import json
import os
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data", "weighted_objective")
rows = [json.loads(l) for l in open(os.path.join(D, "per_amplicon.jsonl"))]
S = json.load(open(os.path.join(D, "summary.json")))
elig = json.load(open(os.path.join(D, "population.json")))["meta"][
    "eligible_source_genomes_by_stratum_excl_twin"]

MULTI = ("n2", "n3", "n4+")
STRATA = ("n2", "n3", "n4+", "n1")          # seed index order of the harness
B, SEED, SESOI, SESOI_DEG = 2000, 1729, 0.01, 0.10
ALT_SEED = 424242                            # an independent draw, sensitivity only
W = {s: elig[s] / sum(elig[t] for t in MULTI) for s in MULTI}

ok = [r for r in rows if r["status"] == "ok" and r.get("arms")]
by = defaultdict(list)
for r in ok:
    by[r["stratum"]].append(r)


def val(r, arm, met):
    a = r["arms"].get(arm) or {}
    if met.startswith("faces."):
        _f, tag, name = met.split(".", 2)
        return ((a.get("faces") or {}).get(tag) or {}).get(name)
    return a.get(met)


def boot(d, cl, seed):
    keys = sorted(set(cl))
    ix = {k: i for i, k in enumerate(keys)}
    s, n = np.zeros(len(keys)), np.zeros(len(keys))
    for x, c in zip(d, cl):
        s[ix[c]] += x
        n[ix[c]] += 1
    dr = np.random.default_rng(seed).integers(0, len(keys), size=(B, len(keys)))
    return float(s.sum() / n.sum()), s[dr].sum(1) / n[dr].sum(1), len(keys)


def pct(b):
    return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


CON = [("real", "two_tier"), ("real", "shuffled"), ("real", "uniform"), ("oracle", "two_tier"),
       ("shuffled", "two_tier"), ("real_exp2", "two_tier"), ("real_ex1", "two_tier_ex1")]
MET = ["taf_face_mid", "taf_face_min", "taf_face_max", "taf_gb", "taf_point_D",
       "faces.opt.n_var_1e-3", "faces.opt7.n_var_1e-3"]

rec = {}
for a, b in CON:
    ck = f"{a}-{b}"
    for met in MET:
        per, bg, bg2 = {}, {}, {}
        for s in STRATA:
            pr = [(val(r, a, met), val(r, b, met), r) for r in by[s]]
            pr = [(x, y, r) for x, y, r in pr if x is not None and y is not None]
            if not pr:
                continue
            d = np.array([x - y for x, y, _ in pr])
            cl = [r["group_id"] for *_, r in pr]
            m, bg[s], K = boot(d, cl, SEED + 1000 * STRATA.index(s))
            _m, bg2[s], _ = boot(d, cl, ALT_SEED + 1000 * STRATA.index(s))
            per[s] = {"n": len(d), "K_groups": K, "mean": m, "median": float(np.median(d)),
                      "ci": pct(bg[s]), "ci_alt_seed": pct(bg2[s]),
                      "comparator_mean": float(np.mean([y for _, y, _ in pr])),
                      "arm_mean": float(np.mean([x for x, _, _ in pr])),
                      "n_lt0": int((d < 0).sum()), "n_eq0": int((d == 0).sum()),
                      "n_gt0": int((d > 0).sum())}
        if all(s in per for s in MULTI):
            per["pooled"] = {
                "mean": sum(W[s] * per[s]["mean"] for s in MULTI),
                "ci": pct(sum(W[s] * bg[s] for s in MULTI)),
                "ci_alt_seed": pct(sum(W[s] * bg2[s] for s in MULTI)),
                "comparator_mean": sum(W[s] * per[s]["comparator_mean"] for s in MULTI),
                "arm_mean": sum(W[s] * per[s]["arm_mean"] for s in MULTI)}
        rec.setdefault(ck, {})[met] = per

# ---------------------------------------------------------------- agreement with summary.json
cmp_n, cmp_bad, maxdiff = 0, [], 0.0
for ck, mets in rec.items():
    for met, per in mets.items():
        for s, v in per.items():
            sk = "multi" if s == "pooled" else s
            ref = S["contrasts"][ck][met][sk]
            refm = ref.get("mean_delta", ref.get("mean_delta_weighted"))
            for mine, theirs in ((v["mean"], refm), (v["ci"][0], ref["ci_group"][0]),
                                 (v["ci"][1], ref["ci_group"][1])):
                cmp_n += 1
                dd = abs(mine - theirs)
                maxdiff = max(maxdiff, dd)
                if dd > 1e-12:
                    cmp_bad.append((ck, met, s, mine, theirs))
alt_shift = max(abs(v["ci_alt_seed"][i] - v["ci"][i]) for mets in rec.values()
                for met, per in mets.items() if met.startswith("taf_face")
                for v in per.values() for i in (0, 1))


# ---------------------------------------------------------------- decision rules, re-applied
def P(ck, met):
    return rec[ck][met]["pooled"]


E = ("taf_face_mid", "taf_face_min", "taf_face_max")
orc = {e: P("oracle-two_tier", e) for e in E}
rt = {e: P("real-two_tier", e) for e in E}
rs = {e: P("real-shuffled", e) for e in E}
c3 = (orc["taf_face_mid"]["ci"][1] < 0 and orc["taf_face_mid"]["mean"] <= -SESOI and
      orc["taf_face_min"]["ci"][1] < 0 and orc["taf_face_max"]["ci"][1] < 0)
m, c = rt["taf_face_mid"]["mean"], rt["taf_face_mid"]["ci"]
ends_below = rt["taf_face_min"]["ci"][1] < 0 and rt["taf_face_max"]["ci"][1] < 0
ends_above = rt["taf_face_min"]["ci"][0] > 0 and rt["taf_face_max"]["ci"][0] > 0
ends_inside = all(-SESOI < rt[e]["ci"][0] and rt[e]["ci"][1] < SESOI for e in E[1:])
spec = rs["taf_face_mid"]["ci"][1] < 0
if c[1] < 0 and m <= -SESOI and ends_below and spec:
    h1 = "works"
elif c[0] > 0 and m >= SESOI and ends_above:
    h1 = "harmful"
elif -SESOI < c[0] and c[1] < SESOI and ends_inside:
    h1 = "does nothing (equivalence)" if c[0] <= 0 <= c[1] else "equivalent within SESOI, nonzero"
elif c[1] < 0 and m <= -SESOI:
    h1 = "lowers by >= SESOI at the midpoint, but ..."
elif c[0] > 0 and m >= SESOI:
    h1 = "raises by >= SESOI at the midpoint, but ..."
elif c[1] < 0 or c[0] > 0:
    h1 = "detectable but below SESOI"
else:
    h1 = "inconclusive"


def h2(met):
    p, q = P("real-two_tier", met), P("real-shuffled", met)
    rel = -p["mean"] / p["comparator_mean"]
    lab = ("reduces degeneracy (>= 10%)" if p["ci"][1] < 0 and rel >= SESOI_DEG else
           "reduces, below 10%" if p["ci"][1] < 0 else
           "increases" if p["ci"][0] > 0 else "no detectable change")
    return lab + ("; specific to p" if q["ci"][1] < 0 else "; NOT specific to p"), rel


h2_opt, rel_opt = h2("faces.opt.n_var_1e-3")
h2_opt7, rel_opt7 = h2("faces.opt7.n_var_1e-3")
verdict = {"C3_oracle_expressible": bool(c3), "H1": h1, "H1s_specific_to_p": bool(spec),
           "H2": h2_opt, "H2_rel_reduction": rel_opt, "H2_1e-7": h2_opt7,
           "H2_1e-7_rel_reduction": rel_opt7,
           "headroom_captured_face_mid": rt["taf_face_mid"]["mean"] / orc["taf_face_mid"]["mean"],
           "headroom_captured_point": P("real-two_tier", "taf_gb")["mean"] /
           P("oracle-two_tier", "taf_gb")["mean"]}
for name, mine in (("H1", h1.split(":")[0]), ("H1s_specific_to_p", spec),
                   ("positive_control_oracle_expressible", c3)):
    theirs = S["verdict"][name]
    verdict[f"agrees_summary_{name}"] = (theirs == mine if isinstance(mine, bool)
                                         else str(theirs).startswith(mine))

# ---------------------------------------------------------------- descriptive numbers
multi = [r for r in ok if r["stratum"] in MULTI]
n1 = by["n1"]
desc = {}
dm = np.array([val(r, "real", "taf_face_mid") - val(r, "two_tier", "taf_face_mid") for r in multi])
desc["concentration_real_minus_two_tier_face_mid"] = {
    "n": len(dm), "n_le_-0.01": int((dm <= -0.01).sum()), "n_ge_+0.01": int((dm >= 0.01).sum()),
    "quantiles_5_25_50_75_95": [float(np.percentile(dm, q)) for q in (5, 25, 50, 75, 95)],
    "median_by_stratum": {s: rec["real-two_tier"]["taf_face_mid"][s]["median"] for s in STRATA},
    "top30_reductions_share_of_net": float(np.sort(dm)[:30].sum() / dm.sum())}
ish = lambda r, arm: r["arms"][arm]["internal_share"]
plt1 = lambda r, arm: (r["arms"][arm]["internal_term_split"].get("gene_backed_plt1", 0.0) /
                       r["arms"][arm]["objective"])
desc["composition"] = {
    "exchange_share_real_median_all360": float(np.median([1 - ish(r, "real") for r in ok])),
    "exchange_share_real_min_max_all360": [float(min(1 - ish(r, "real") for r in ok)),
                                           float(max(1 - ish(r, "real") for r in ok))],
    "exchange_share_real_median_multi": float(np.median([1 - ish(r, "real") for r in multi])),
    "internal_share_two_tier_median_multi": float(np.median([ish(r, "two_tier") for r in multi])),
    "internal_share_uniform_median_multi": float(np.median([ish(r, "uniform") for r in multi])),
    "internal_share_real_ex1_median_multi": float(np.median([ish(r, "real_ex1") for r in multi])),
    "plt1_term_share_of_objective_real_median_multi": float(np.median([plt1(r, "real") for r in multi])),
    "plt1_term_share_of_objective_real_median_all360": float(np.median([plt1(r, "real") for r in ok])),
    "geneless_draft_share_of_internal_real_median_multi": float(np.median(
        [r["arms"]["real"]["internal_term_split"]["geneless_draft"] /
         r["arms"]["real"]["internal_term"] for r in multi])),
    "exchange_spread_ex100_le_1e-6": int(sum(r["exchange_term_rel_spread_ex100"] <= 1e-6 for r in ok)),
    "exchange_spread_ex100_max": float(max(r["exchange_term_rel_spread_ex100"] for r in ok))}
desc["addressable_share_two_tier_mean"] = {
    s: float(np.mean([r["arms"]["two_tier"]["flux_plt1_absent"] / r["arms"]["two_tier"]["flux_absent"]
                      for r in by[s] if r["arms"]["two_tier"]["flux_absent"] > 0])) for s in MULTI}
desc["static_precision_multi_mean"] = {
    "plt1": float(np.mean([r["model"]["static_precision_plt1"] for r in multi
                           if r["model"]["static_precision_plt1"] is not None])),
    "p1": float(np.mean([r["model"]["static_precision_p1"] for r in multi]))}
desc["model_size_median"] = {
    "n_gene_backed_multi": float(np.median([r["model"]["n_gene_backed"] for r in multi])),
    "n_gene_backed_plt1_multi": float(np.median([r["model"]["n_gene_backed_plt1"] for r in multi])),
    **{f"n_gene_backed_plt1_{s}": float(np.median([r["model"]["n_gene_backed_plt1"] for r in by[s]]))
       for s in MULTI}}
desc["face_width_mean"] = {s: {arm: float(np.mean([val(r, arm, "taf_face_width") for r in by[s]]))
                               for arm in ("real", "two_tier")} for s in MULTI}
desc["arm_levels_pooled"] = {}
for arm in ("two_tier", "uniform", "real", "shuffled", "oracle"):
    desc["arm_levels_pooled"][arm] = {
        met: sum(W[s] * float(np.mean([val(r, arm, met) for r in by[s]])) for s in MULTI)
        for met in ("taf_face_mid", "faces.opt.n_var_1e-3")
        if val(by["n2"][0], arm, met) is not None}
# growth
ratios = [r["arms"][a]["biomass"] / r["arms"][a]["max_biomass"] for r in ok for a in r["arms"]]
desc["growth"] = {"n_arm_solves": len(ratios), "ratio_min": float(min(ratios)),
                  "ratio_max": float(max(ratios)),
                  "n_optimal": sum(r["arms"][a]["status"] == "optimal" for r in ok for a in r["arms"]),
                  "max_biomass_rel_spread_across_arms_max": float(max(
                      (max(a["max_biomass"] for a in r["arms"].values()) -
                       min(a["max_biomass"] for a in r["arms"].values())) /
                      max(a["max_biomass"] for a in r["arms"].values()) for r in ok))}
# controls and health
rc = [r["released_check"] for r in ok]
desc["checks"] = {
    "status": dict(Counter(r["status"] for r in rows)),
    "released_check_pass": sum(x["pass"] for x in rc), "released_check_n": len(rc),
    "released_check_max_rel_diff": max(x["rel_diff"] for x in rc),
    "released_sha256": sorted({x["released_source_sha256"] for x in rc}),
    "n1_identical": sum(bool(r["identical_real_two_tier"]) for r in n1), "n1_n": len(n1),
    "n1_face_max_abs_diff": float(max(max(abs(val(r, "real", k) - val(r, "two_tier", k))
                                          for k in ("taf_face_min", "taf_face_max")) for r in n1)),
    "n1_face_n_within_1e-9": int(sum(max(abs(val(r, "real", k) - val(r, "two_tier", k))
                                         for k in ("taf_face_min", "taf_face_max")) <= 1e-9 for r in n1)),
    "face_lp_failures": sum(f.get("n_lp_failed", 0) + f.get("n_interval_lp_failed", 0)
                            for r in ok for a in r["arms"].values()
                            for f in (a.get("faces") or {}).values()),
    "face_lp_retried": sum(f.get("n_lp_retried", 0) for r in ok for a in r["arms"].values()
                           for f in (a.get("faces") or {}).values()),
    "gapfill_sanity_ok": sum(bool(r["gapfill"]["sanity_gate_ok"]) for r in ok),
    "gapfill_n_new_rxn_median": float(np.median([r["gapfill"]["n_new_rxn"] for r in ok])),
    "models_with_reversals": sum(r["gapfill"]["n_reversed_applied"] > 0 for r in ok),
    "reversals_total": sum(r["gapfill"]["n_reversed_applied"] for r in ok),
    "D_face_rel_width_max": float(max(r["model"]["D_face_rel_width"] for r in ok)),
    "wall_mean_s": float(np.mean([r["wall_seconds"] for r in ok])),
    "wall_max_s": float(max(r["wall_seconds"] for r in ok)),
    "worker_maxrss_gb_max": float(max(r["worker_maxrss_gb"] for r in ok)),
    "opt7_flags": {a: sum(bool(r["arms"][a].get("opt7_count_flag")) for r in ok)
                   for a in ("two_tier", "real", "uniform", "shuffled")},
    "taf_gb_vertex_dependent_rows": {a: int(sum(r["arms"][a]["taf_gb_vertex_spread"] > 1e-6 for r in ok))
                                     for a in ("real", "two_tier", "oracle")},
    "taf_gb_vertex_spread_max": {a: float(max(r["arms"][a]["taf_gb_vertex_spread"] for r in ok))
                                 for a in ("real", "two_tier", "oracle")},
    "n1_bitwise_returned_identical_note": "identical_real_two_tier per row"}
desc["deployment_share_n1"] = "86.8% (prFBA.tex l.210; not recomputed here)"

out = {"inputs": {"rows": len(rows), "ok": len(ok), "by_stratum": {s: len(by[s]) for s in STRATA},
                  "groups_by_stratum": {s: len({r["group_id"] for r in by[s]}) for s in STRATA},
                  "weights": W, "bootstrap": B, "seed": SEED},
       "agreement_with_summary": {"cells_compared": cmp_n, "cells_differing_gt_1e-12": len(cmp_bad),
                                  "max_abs_diff": maxdiff, "first_disagreements": cmp_bad[:10],
                                  "alt_seed_max_ci_shift_taf_face": alt_shift},
       "verdict_reapplied": verdict, "contrasts": rec, "descriptive": desc}
json.dump(out, open(os.path.join(D, "recomputed_contrasts.json"), "w"), indent=1)

# ---------------------------------------------------------------- figures
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
COL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # default categorical slots 1-4
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF})
GROUPS = ["n2", "n3", "n4+", "pooled"]
GLAB = ["n = 2", "n = 3", "n ≥ 4", "pooled n ≥ 2"]


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels(GLAB)
    ax.axvline(2.5, color=GRID, lw=0.8)


# Figure 1: flux contrasts on the optimal-face midpoint
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.axhspan(-SESOI, SESOI, color="#f0efec", zorder=0)
ax.axhline(0, color=INK2, lw=0.8)
F = [("real-two_tier", "real − two-tier"), ("real-shuffled", "real − shuffled"),
     ("oracle-two_tier", "oracle − two-tier")]
for j, (ck, lab) in enumerate(F):
    xs = np.arange(len(GROUPS)) + (j - 1) * 0.22
    ms = [rec[ck]["taf_face_mid"][g]["mean"] for g in GROUPS]
    lo = [rec[ck]["taf_face_mid"][g]["ci"][0] for g in GROUPS]
    hi = [rec[ck]["taf_face_mid"][g]["ci"][1] for g in GROUPS]
    ax.errorbar(xs, ms, yerr=[np.subtract(ms, lo), np.subtract(hi, ms)], fmt="o", ms=5,
                color=COL[j], ecolor=COL[j], elinewidth=1.6, capsize=0, label=lab)
style(ax)
ax.text(3.45, SESOI * 0.55, "±SESOI", color=INK2, fontsize=7.5, ha="right", va="center")
ax.set_ylabel("Δ truth-absent flux, face midpoint\n(fraction of two-tier gene-backed flux)")
ax.legend(frameon=False, fontsize=8, loc="lower left", ncol=3, bbox_to_anchor=(0, 1.0))
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(D, f"fig_wobj_flux.{ext}"), dpi=200)
plt.close(fig)

# Figure 2: degeneracy, variable gene-backed reactions on the optimal face
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ARMS = [("two_tier", "two-tier"), ("uniform", "uniform"), ("real", "real p"),
        ("shuffled", "shuffled p")]
met = "faces.opt.n_var_1e-3"
lev = {}
for arm, _ in ARMS:
    for g in GROUPS:
        if g == "pooled":
            lev[arm, g] = sum(W[s] * lev[arm, s] for s in MULTI)
        else:
            lev[arm, g] = float(np.mean([val(r, arm, met) for r in by[g]]))
for j, (arm, lab) in enumerate(ARMS):
    xs = np.arange(len(GROUPS)) + (j - 1.5) * 0.19
    ax.bar(xs, [lev[arm, g] for g in GROUPS], width=0.17, color=COL[j], label=lab,
           edgecolor=SURF, linewidth=1.0)
    if arm in ("real", "shuffled"):     # CI of the paired difference vs two-tier, at the arm's bar
        ck = f"{arm}-two_tier" if arm == "shuffled" else "real-two_tier"
        for x, g in zip(xs, GROUPS):
            ci = rec[ck][met][g]["ci"]
            ax.plot([x, x], [lev["two_tier", g] + ci[0], lev["two_tier", g] + ci[1]],
                    color=INK, lw=1.2)
style(ax)
ax.set_ylabel("variable gene-backed reactions\non the optimal face (range > 1e-3)")
ax.legend(frameon=False, fontsize=8, loc="lower left", ncol=4, bbox_to_anchor=(0, 1.0))
ax.text(0.99, -0.2, "Black ticks: two-tier mean + 95% CI of the paired difference",
        transform=ax.transAxes, ha="right", fontsize=7, color=INK2)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(D, f"fig_wobj_degeneracy.{ext}"), dpi=200)
plt.close(fig)

print(json.dumps({"agreement": out["agreement_with_summary"], "verdict": verdict}, indent=1))
