"""Own helpers for the Stage O verification (written from the design + the bootstrap
convention stated in the brief; no import of the author's scripts)."""
import hashlib
import json
import pickle
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference"
STRATA = ("n2", "n3", "n4+", "n1")     # seed index order
MULTI = ("n2", "n3", "n4+")
B = 2000
SEED = 1729
ALT_SEED = 424242
SHIPPED = "kn/lin/G0/X100"
SELECTED = "est/lin/G0/X1"


def load_rows():
    return pickle.load(open(os.path.join(HERE, "rows.pkl"), "rb"))


def load_prior():
    return pickle.load(open(os.path.join(HERE, "prior.pkl"), "rb"))


def half_of(gid):
    h = int(hashlib.sha256(("1729|" + str(gid)).encode()).hexdigest(), 16)
    return "tune" if h % 2 == 0 else "confirm"


def weights():
    pop = json.load(open(os.path.join(SR, "data/weighted_objective/population.json")))
    e = pop["meta"]["eligible_source_genomes_by_stratum_excl_twin"]
    tot = sum(e[s] for s in MULTI)
    WA = {s: e[s] / tot for s in MULTI}
    WB = {"n1": 0.868}
    for s in MULTI:
        WB[s] = 0.132 * WA[s]
    return WA, WB, e


def v(row, arm, met):
    a = row["arms"].get(arm)
    if a is None or a.get("status") != "optimal":
        return None
    x = a.get(met)
    if x is None:
        return None
    x = float(x)
    return None if x != x else x


def boot_stratum(d, groups, seed, B_=B):
    """Ratio-of-sums cluster bootstrap over sorted unique group keys."""
    keys = sorted(set(groups))
    pos = {k: i for i, k in enumerate(keys)}
    K = len(keys)
    S = np.zeros(K)
    N = np.zeros(K)
    for x, g in zip(d, groups):
        S[pos[g]] += x
        N[pos[g]] += 1
    dr = np.random.default_rng(seed).integers(0, K, size=(B_, K))
    return float(S.sum() / N.sum()), S[dr].sum(1) / N[dr].sum(1), K


def ci(b):
    return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


def by_stratum(rows):
    out = defaultdict(list)
    for r in rows:
        out[r["stratum"]].append(r)
    return out


def paired(rows, fa, fb=None, seed_base=SEED, pools=None, B_=B):
    """fa, fb: functions row -> value (or None).  d = fa - fb (fb None -> fa alone).
    Per-stratum mean + CI (seed_base + 1000*idx), pooled per weights dict in pools."""
    by = by_stratum(rows)
    per, boots = {}, {}
    for i, s in enumerate(STRATA):
        d, g = [], []
        for r in by.get(s, []):
            x = fa(r)
            y = fb(r) if fb is not None else 0.0
            if x is None or y is None:
                continue
            d.append(x - y)
            g.append(r["group_id"])
        if not d:
            continue
        m, bs, K = boot_stratum(d, g, seed_base + 1000 * i, B_)
        per[s] = {"n": len(d), "K": K, "mean": m, "ci": ci(bs),
                  "n_lt0": int(sum(1 for x in d if x < 0)),
                  "n_eq0": int(sum(1 for x in d if x == 0)),
                  "n_gt0": int(sum(1 for x in d if x > 0))}
        boots[s] = bs
    pooled = {}
    for name, W in (pools or {}).items():
        if not all(s in per for s in W):
            pooled[name] = None
            continue
        pooled[name] = {"mean": sum(W[s] * per[s]["mean"] for s in W),
                        "ci": ci(sum(W[s] * boots[s] for s in W))}
    return {"per": per, "pooled": pooled, "boots": boots}


def joint_boot(rows, fa, fb, W, seed=SEED, B_=B):
    """Sensitivity: resample holdout groups ONCE across all strata (a group's rows in
    different strata move together); per draw, each stratum's ratio-of-sums mean over the
    drawn groups' rows in that stratum, pooled with W.  Draws where a stratum gets no rows
    are re-weighted over present strata (counted)."""
    by = by_stratum(rows)
    groups = sorted({r["group_id"] for r in rows if r["stratum"] in W})
    pos = {k: i for i, k in enumerate(groups)}
    K = len(groups)
    S = {s: np.zeros(K) for s in W}
    N = {s: np.zeros(K) for s in W}
    for s in W:
        for r in by.get(s, []):
            x, y = fa(r), (fb(r) if fb is not None else 0.0)
            if x is None or y is None:
                continue
            S[s][pos[r["group_id"]]] += x - y
            N[s][pos[r["group_id"]]] += 1
    dr = np.random.default_rng(seed).integers(0, K, size=(B_, K))
    tot = np.zeros(B_)
    wsum = np.zeros(B_)
    nmiss = 0
    for s in W:
        num = S[s][dr].sum(1)
        den = N[s][dr].sum(1)
        okm = den > 0
        nmiss += int((~okm).sum())
        mean_s = np.where(okm, num / np.where(okm, den, 1), 0.0)
        tot += W[s] * mean_s * okm
        wsum += W[s] * okm
    boots = tot / wsum
    point = sum(W[s] * S[s].sum() / N[s].sum() for s in W)
    return float(point), ci(boots), K, nmiss


def fmt(c, nd=4):
    if c is None:
        return "n/a"
    return "{:+.{n}f} [{:+.{n}f}, {:+.{n}f}]".format(c["mean"], c["ci"][0], c["ci"][1], n=nd)
