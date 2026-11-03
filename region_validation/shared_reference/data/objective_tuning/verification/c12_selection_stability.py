"""Selection stability: resample tune-half holdout groups (per stratum, the harness's seeds;
and jointly) and re-select among the 36 candidates by pooled A (lowest; the 0.002 tie rule
toward fewer departures applied).  Also tune vs confirm per-stratum selected-shipped."""
import json, os
from itertools import product
import numpy as np
import vlib as L
rows = L.load_rows(); WA, WB, _ = L.weights()
tune = [r for r in rows if r["half"] == "tune"]; conf = [r for r in rows if r["half"] == "confirm"]
CAND = ["/".join(t) for t in product(("kn", "est"), ("lin", "log"), ("G0", "Gmeas", "G0.2"), ("X100", "X10", "X1"))]
deps = lambda a: sum(x != y for x, y in zip(a.split("/"), L.SHIPPED.split("/")))
by = L.by_stratum(tune)
def draws_per_stratum(seed_base, joint=False, B=2000):
    lev = np.zeros((B, len(CAND)))
    if joint:
        groups = sorted({r["group_id"] for r in tune if r["stratum"] in L.MULTI})
        pos = {g: i for i, g in enumerate(groups)}
        dr = np.random.default_rng(seed_base).integers(0, len(groups), size=(B, len(groups)))
        cnt = np.zeros((B, len(groups)))
        for b in range(B):
            cnt[b] = np.bincount(dr[b], minlength=len(groups))
    for i, s in enumerate(L.MULTI):
        rs = by[s]
        if joint:
            S = np.zeros((len(groups), len(CAND))); N = np.zeros(len(groups))
            for r in rs:
                S[pos[r["group_id"]]] += [r["arms"][a]["taf_face_mid"] for a in CAND]; N[pos[r["group_id"]]] += 1
            num = cnt @ S; den = cnt @ N
            lev += WA[s] * num / den[:, None]
        else:
            keys = sorted({r["group_id"] for r in rs}); pos2 = {g: j for j, g in enumerate(keys)}
            S = np.zeros((len(keys), len(CAND))); N = np.zeros(len(keys))
            for r in rs:
                S[pos2[r["group_id"]]] += [r["arms"][a]["taf_face_mid"] for a in CAND]; N[pos2[r["group_id"]]] += 1
            dr = np.random.default_rng(seed_base + 1000 * i).integers(0, len(keys), size=(B, len(keys)))
            lev += WA[s] * S[dr].sum(1) / N[dr].sum(1)[:, None]
    winners = []
    for b in range(B):
        best = lev[b].min()
        ties = [CAND[j] for j in range(len(CAND)) if lev[b, j] <= best + 0.002]
        winners.append(min(ties, key=lambda a: (deps(a), lev[b, CAND.index(a)])))
    from collections import Counter
    c = Counter(winners)
    fam = Counter("/".join((w.split("/")[0], w.split("/")[1], w.split("/")[3])) for w in winners)
    return {"winner_freq": dict(c.most_common(8)), "P_est_share": sum(v for k, v in c.items() if k.startswith("est")) / B,
            "X1_share": sum(v for k, v in c.items() if k.endswith("X1")) / B, "family_freq": dict(fam.most_common(6)),
            "est_lin_X1_any_G_share": sum(v for k, v in c.items() if k.startswith("est/lin") and k.endswith("X1")) / B}
out = {"per_stratum_seed1729": draws_per_stratum(1729), "joint_seed1729": draws_per_stratum(1729, joint=True)}
# tune vs confirm, per stratum, selected - shipped
for h, rs in (("tune", tune), ("confirm", conf)):
    c = L.paired(rs, lambda r: L.v(r, L.SELECTED, "taf_face_mid"), lambda r: L.v(r, L.SHIPPED, "taf_face_mid"), pools={"A": WA, "B": WB})
    out[f"{h}_selected_minus_shipped"] = {"per": {s: [round(p["mean"], 4), [round(x, 4) for x in p["ci"]], p["n"], p["K"]] for s, p in c["per"].items()},
                                          "A": L.fmt(c["pooled"]["A"]), "B": L.fmt(c["pooled"]["B"])}
json.dump(out, open(os.path.join(L.HERE, "out_c12.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
