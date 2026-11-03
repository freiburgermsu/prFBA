"""Phase 4: metrics, strata and cluster-bootstrap CIs (holdout groups; score_folds' draw:
default_rng(seed).integers(0, K, size=(B, K)) over the stratum's sorted group keys) for the
composites of phase 2; non-nested grid / probe point estimates from pass A."""
import collections
import json
import os
import pickle
import sys

import numpy as np

import vctx

B = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
SEED = 1729
NU = 2235
cells = vctx.load_cells("1.0")
nC = len(cells)
nreps = cells.n_reps.astype(int).to_numpy()
UQn = np.load(os.path.join(vctx.CACHE, "union_quantities.npz"))
Ns, T, POS = UQn["Ns"], UQn["T"], UQn["POS"]
F1U = 2 * POS / (Ns + T)
with open(os.path.join(vctx.CACHE, "passB.pkl"), "rb") as fh:
    PB = pickle.load(fh)
comps = PB["comps"]
ST = PB["ST"]                      # nC x ncomp x nstat
summ = json.load(open(os.path.join(vctx.OUT, "p2_summary.json")))
STATS = summ["stats"]
S = {k: j for j, k in enumerate(STATS)}
UQ, UY, UC = PB["UQ"], PB["UY"], PB["UC"]

strata = {"all": np.ones(nC, bool), "n1": nreps == 1, "n2": nreps == 2, "n3": nreps == 3,
          "n4+": nreps >= 4, "multi": nreps >= 2}


def metrics(Ssum, n_cells):
    """Ssum: dict of summed columns (scalars or arrays over replicates)."""
    o = {}
    with np.errstate(all="ignore"):
        for c in comps:
            g = lambda k: Ssum[f"{c}:{k}"]  # noqa: E731
            o[c] = {
                "universe_brier": g("SQ") / (n_cells * NU),
                "universe_logloss": g("LL") / (n_cells * NU),
                "union_brier": g("SQs") / Ssum["Ns"],
                "union_wls_slope": (g("Sqys") - g("Sqs") * Ssum["POS"] / Ssum["Ns"]) /
                                   (g("Sqqs") - g("Sqs") ** 2 / Ssum["Ns"]),
                "macro_dF1_vs_union": (g("F1") - Ssum["F1U"]) / n_cells,
                "micro_dF1_vs_union": 2 * g("TP") / (g("PRED") + Ssum["T"]) -
                                      2 * Ssum["POS"] / (Ssum["Ns"] + Ssum["T"]),
                "macro_dF1_vs_union_restricted_to_union": (g("F1u") - Ssum["F1U"]) / n_cells,
                "support_per_cell": g("NSUP") / n_cells,
                "pred_q_ge_0.5_outside_union_per_cell": (g("PRED") - g("PREDu")) / n_cells,
                "tp_outside_union_per_cell": (g("TP") - g("TPu")) / n_cells,
            }
    return o


def auroc_exact(q, y):
    order = np.argsort(q, kind="mergesort")
    qs = q[order]
    ys = y[order]
    # average ranks for ties
    bound = np.flatnonzero(np.r_[True, qs[1:] != qs[:-1], True])
    starts = bound[:-1]
    ends = bound[1:]
    avg = (starts + ends + 1) / 2.0           # 1-based average rank of each tie block
    lens = ends - starts
    ranks = np.repeat(avg, lens)
    P = ys.sum()
    N = len(ys) - P
    return float((ranks[ys].sum() - P * (P + 1) / 2) / (P * N))


def auroc_hist(Hp, Hn):
    above = np.cumsum(Hp[..., ::-1], -1)[..., ::-1] - Hp
    num = (Hn * (above + 0.5 * Hp)).sum(-1)
    with np.errstate(all="ignore"):
        return num / (Hp.sum(-1) * Hn.sum(-1))


res = {}
gid = cells.group_id.to_numpy()
for sn, mask in strata.items():
    idx = np.flatnonzero(mask)
    keys, inv = np.unique(gid[idx], return_inverse=True)
    K = len(keys)
    rng = np.random.default_rng(SEED)
    draws = rng.integers(0, K, size=(B, K))
    W = np.zeros((B, K))
    for b in range(B):
        W[b] = np.bincount(draws[b], minlength=K)
    cols = {"Ns": Ns, "POS": POS, "T": T, "F1U": F1U, "NCELL": np.ones(nC)}
    for c, nm in enumerate(comps):
        for st in STATS:
            cols[f"{nm}:{st}"] = ST[:, c, S[st]]
    names = list(cols)
    M = np.stack([np.bincount(inv, weights=cols[k][idx], minlength=K) for k in names], 1)
    tot = M.sum(0)
    rep = W @ M
    pt_sum = {k: tot[j] for j, k in enumerate(names)}
    bt_sum = {k: rep[:, j] for j, k in enumerate(names)}
    pt = metrics(pt_sum, pt_sum["NCELL"])
    bt = metrics(bt_sum, bt_sum["NCELL"])
    # AUROC on the union support (exact estimate, 1,001-bin replicates)
    cmask = mask[UC]
    yy = UY[cmask]
    lab = inv[np.searchsorted(idx, UC[cmask])]
    for c, nm in enumerate(comps):
        q = UQ[c, cmask].astype(np.float64)
        pt[nm]["union_auroc"] = auroc_exact(q, yy) if (yy.any() and (~yy).any()) else float("nan")
        bins = np.where(q >= 1.0, 1000, np.minimum((q * 1000).astype(np.int64), 999))
        key = lab * 1001 + bins
        Hp = np.bincount(key[yy], minlength=K * 1001).reshape(K, 1001).astype(float)
        Hn = np.bincount(key[~yy], minlength=K * 1001).reshape(K, 1001).astype(float)
        pt[nm]["union_auroc_binned"] = float(auroc_hist(Hp.sum(0), Hn.sum(0)))
        bt[nm]["union_auroc"] = auroc_hist(W @ Hp, W @ Hn)
    out = {"n_cells": int(len(idx)), "n_groups": int(K), "est": {}, "ci": {}, "contrasts": {}}
    metr = list(next(iter(pt.values())).keys())
    for nm in comps:
        out["est"][nm] = {m: float(pt[nm][m]) for m in metr if m in pt[nm]}
        out["ci"][nm] = {}
        for m in pt[nm]:
            if m == "union_auroc_binned":
                continue
            r = np.asarray(bt[nm][m], float)
            r = r[np.isfinite(r)]
            out["ci"][nm][m] = ([float(x) for x in np.percentile(r, (2.5, 97.5))]
                                if len(r) >= 10 else None)
    pairs = [(a, "E0") for a in comps if a != "E0"] + [
        ("E3", "E1"), ("E3", "E4(E0)"), ("E1", "E4(E0)"), ("E4(E3)", "E3"),
        ("E3", "E1|K=100|tau=0.02"), ("E3x|tau=0.05|family|alpha=2.0", "E3"),
        ("E3x|tau=inf|family|alpha=2.0", "E3"), ("E3", "PREV_foldref"),
        ("E3", "CLADE_family_only"), ("E1x|K=100|tau=inf", "E1")]
    for a, b in pairs:
        if a not in comps or b not in comps:
            continue
        d = {}
        for m in ("universe_brier", "universe_logloss", "union_brier", "union_wls_slope",
                  "union_auroc", "macro_dF1_vs_union", "micro_dF1_vs_union"):
            e = pt[a][m] - pt[b][m]
            r = np.asarray(bt[a][m], float) - np.asarray(bt[b][m], float)
            r = r[np.isfinite(r)]
            d[m] = [float(e)] + ([float(x) for x in np.percentile(r, (2.5, 97.5))]
                                 if len(r) >= 10 else [None, None])
        out["contrasts"][f"{a} - {b}"] = d
    res[sn] = out
    print(f"[analyze] {sn}: {len(idx)} cells, {K} groups", flush=True)

# --------------------------------------------------------------- non-nested grid (pass A)
A = np.load(os.path.join(vctx.CACHE, "passA.npz"))["ST"]
V = summ["variants"]
grid = {}
for sn in ("all", "n1", "multi"):
    m = strata[sn]
    n = m.sum()
    g = {}
    for j, nm in enumerate(V):
        s = A[m, j].sum(0)
        g[nm] = {"universe_brier": s[S["SQ"]] / (n * NU),
                 "universe_logloss": s[S["LL"]] / (n * NU),
                 "union_brier": s[S["SQs"]] / Ns[m].sum(),
                 "union_wls_slope": (s[S["Sqys"]] - s[S["Sqs"]] * POS[m].sum() / Ns[m].sum()) /
                                    (s[S["Sqqs"]] - s[S["Sqs"]] ** 2 / Ns[m].sum()),
                 "macro_dF1_vs_union": (s[S["F1"]] - F1U[m].sum()) / n,
                 "macro_dF1_vs_union_restricted_to_union": (s[S["F1u"]] - F1U[m].sum()) / n,
                 "support_per_cell": s[S["NSUP"]] / n}
    grid[sn] = g
# per-fold training log-loss of the probes vs the declared optimum (E3 family, E1)
fold = cells.fold.to_numpy()
LLf = np.stack([A[fold == f][:, :, S["LL"]].sum(0) for f in range(5)], 1)
NCf = np.array([(fold == f).sum() for f in range(5)], float)
train = {}
for j, nm in enumerate(V):
    train[nm] = {f: float(LLf[j, [x for x in range(5) if x != f]].sum() /
                          (NCf[[x for x in range(5) if x != f]].sum() * NU)) for f in range(5)}
with open(os.path.join(vctx.OUT, "p4_results.json"), "w") as fh:
    json.dump({"B": B, "seed": SEED, "strata": res, "grid_pooled_not_nested": grid,
               "train_ll_by_variant_and_outer_fold": train}, fh, indent=1, default=float)
print("[analyze] written p4_results.json")
