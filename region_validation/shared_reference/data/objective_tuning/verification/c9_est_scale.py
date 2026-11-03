"""How est (E3) and kn place gene-backed reactions on the cost scale (from est_info's
2-decimal histogram and model.p_values_gene_backed; the per-reaction values are not stored).
lin cost = max(0.01, 1 - p): the floor 0.01 is reached at p >= 0.99.  log cost = max(0.01,
ln p / ln 0.01): the floor is reached at p >= 0.01**0.01 = 0.9550."""
import json
import math
import os

import numpy as np

import vlib as L

rows = L.load_rows()
P_LOG_FLOOR = 0.01 ** 0.01


def lin(p):
    return max(0.01, 1 - p)


def logc(p):
    return max(0.01, math.log(max(p, 0.01)) / math.log(0.01))


per = []
for r in rows:
    e, m = r["est_info"], r["model"]
    n = m["n_gene_backed"]
    h = {float(k): v for k, v in e["est_hist_gene_backed"].items()}
    assert sum(h.values()) == n
    kn = {float(k): v for k, v in m["p_values_gene_backed"].items()}
    assert sum(kn.values()) == n
    exact1 = n - e["n_gene_backed_est_lt1"]
    b100 = h.get(1.0, 0)                       # p >= 0.995 -> lin cost 0.01 for sure
    b099 = h.get(0.99, 0)                      # 0.985 <= p < 0.995: cost 0.010-0.015
    rec = {"key": r["key"], "stratum": r["stratum"], "half": r["half"], "n_gb": n,
           "est_exact1_share": exact1 / n, "est_bin1.0_share": b100 / n,
           "est_ge0.985_share": (b100 + b099) / n,
           "est_lt0.5_share": e["n_gene_backed_est_lt0.5"] / n,
           "kn_p1_share": kn.get(1.0, 0) / n,
           "est_mean": e["mean_est_gene_backed"], "kn_mean": e["mean_kn_gene_backed"],
           "est_lin_cost_mean_hist": sum(v * lin(k) for k, v in h.items()) / n,
           "kn_lin_cost_mean": sum(v * lin(k) for k, v in kn.items()) / n,
           "est_log_floor_share_hist": sum(v for k, v in h.items() if k >= 0.96) / n,
           "kn_log_floor_share": sum(v for k, v in kn.items() if k >= P_LOG_FLOOR) / n,
           "est_log_cost_mean_hist": sum(v * logc(k) for k, v in h.items()) / n,
           "kn_log_cost_mean": sum(v * logc(k) for k, v in kn.items()) / n,
           "auroc_est": e["auroc_supported_vs_absent"]["est"], "auroc_kn": e["auroc_supported_vs_absent"]["kn"],
           "brier_est": e["brier_vs_truth_support"]["est"], "brier_kn": e["brier_vs_truth_support"]["kn"],
           "n_neighbours": e["n_neighbours"], "clade_rank": e["clade_rank"], "alpha": e["alpha"]}
    per.append(rec)


def summ(rs, k):
    xs = np.array([x[k] for x in rs if x[k] is not None], float)
    return {"mean": float(xs.mean()), "median": float(np.median(xs)), "min": float(xs.min()),
            "max": float(xs.max()), "n": len(xs)}


KEYS = ("est_exact1_share", "est_bin1.0_share", "est_ge0.985_share", "est_lt0.5_share",
        "kn_p1_share", "est_mean", "kn_mean", "est_lin_cost_mean_hist", "kn_lin_cost_mean",
        "est_log_floor_share_hist", "kn_log_floor_share", "est_log_cost_mean_hist",
        "kn_log_cost_mean", "auroc_est", "auroc_kn", "brier_est", "brier_kn")
out = {"all": {k: summ(per, k) for k in KEYS}}
for s in L.STRATA:
    out[s] = {k: summ([x for x in per if x["stratum"] == s], k) for k in KEYS}
# pooled reaction counts
tot = sum(x["n_gb"] for x in per)
out["pooled_counts"] = {
    "n_gene_backed_total": tot,
    "est_exact_p1": int(round(sum(x["est_exact1_share"] * x["n_gb"] for x in per))),
    "est_p_ge_0.995": int(round(sum(x["est_bin1.0_share"] * x["n_gb"] for x in per))),
    "est_p_ge_0.985": int(round(sum(x["est_ge0.985_share"] * x["n_gb"] for x in per))),
    "kn_p1": int(round(sum(x["kn_p1_share"] * x["n_gb"] for x in per))),
    "rows_with_zero_est_exact1": sum(1 for x in per if x["est_exact1_share"] == 0),
    "null_clade_rows(alpha0)": sum(1 for x in per if not x["alpha"]),
}
json.dump({"summary": out, "per_row": per}, open(os.path.join(L.HERE, "out_c9.json"), "w"), indent=1)
for grp in ("all",) + L.STRATA:
    o = out[grp]
    print("%-4s est exact p=1 %.3f (median %.3f) | est p>=0.995 %.3f | p>=0.985 %.3f | kn p=1 %.3f | "
          "mean p est %.3f kn %.3f | mean lin cost est %.4f kn %.4f | log-floor share est %.3f kn %.3f | "
          "AUROC est %.3f kn %.3f" % (
              grp, o["est_exact1_share"]["mean"], o["est_exact1_share"]["median"],
              o["est_bin1.0_share"]["mean"], o["est_ge0.985_share"]["mean"], o["kn_p1_share"]["mean"],
              o["est_mean"]["mean"], o["kn_mean"]["mean"], o["est_lin_cost_mean_hist"]["mean"],
              o["kn_lin_cost_mean"]["mean"], o["est_log_floor_share_hist"]["mean"],
              o["kn_log_floor_share"]["mean"], o["auroc_est"]["mean"], o["auroc_kn"]["mean"]))
print(out["pooled_counts"])
