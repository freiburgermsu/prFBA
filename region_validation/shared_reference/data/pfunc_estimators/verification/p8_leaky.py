"""Counterfactual: E3 (K100, tau .02, family start, alpha 2) with the clade prior computed
over the FULL panel (held-out genomes included) instead of the fold reference. If the study
had leaked this way, its E3 would equal this, not the fold-reference E3."""
import json, os
import numpy as np
import p2_pass as P
import vctx
P.setup()
G = P.G
C, E = G["C"], G["E"]
inv = {v: k for k, v in G["clade_key"].items()}
lab = {rk: np.array([x or "" for x in C.lab[rk]], object) for rk in P.CHAIN}
pi_full = {}
for j, (f, rk, L) in inv.items():
    mem = np.flatnonzero(lab[rk] == L)
    pi_full[j] = E[mem].mean(0)
nC = len(G["cells"])
SQ = {"ref": 0.0, "full": 0.0}
SQs = {"ref": 0.0, "full": 0.0}
Ns = 0.0
src_in_full = 0
for i in range(nC):
    f = int(G["fold"][i])
    g1, i1 = G["l100"][i]
    w = np.exp(-(i1[0] - i1) / 0.02)
    num = w @ G["Ef"][f][g1].astype(np.float64)
    den = w.sum()
    y = E[G["src"][i]].astype(np.float64)
    u = G["Ef"][f][G["reps"][i]].max(0) > 0
    j = G["cl_of"]["family"][i]
    for nm in ("ref", "full"):
        if j >= 0:
            p = G["PI"][j] if nm == "ref" else pi_full[j]
            q = (num + 2.0 * p) / (den + 2.0)
        else:
            q = num / den
        q[q >= 1 - 1e-9] = 1.0
        q[q <= 1e-12] = 0.0
        SQ[nm] += ((q - y) ** 2).sum()
        SQs[nm] += ((q[u] - y[u]) ** 2).sum()
    Ns += u.sum()
    if j >= 0:
        rk = inv[j][1]
        src_in_full += int(lab[rk][G["src"][i]] == inv[j][2])
out = {"E3_universe_brier_fold_reference_prior": SQ["ref"] / (nC * C.NU),
       "E3_universe_brier_full_panel_prior_counterfactual": SQ["full"] / (nC * C.NU),
       "E3_union_brier_fold_reference_prior": SQs["ref"] / Ns,
       "E3_union_brier_full_panel_prior_counterfactual": SQs["full"] / Ns,
       "cells_whose_source_would_be_in_the_full_panel_clade": src_in_full,
       "reported_E3_universe_brier": 0.05778561}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(vctx.OUT, "p8_leaky.json"), "w"), indent=1)
