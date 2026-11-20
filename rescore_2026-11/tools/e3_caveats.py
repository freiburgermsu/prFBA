#!/usr/bin/env python
"""Re-derive the PFUNC_E3_GATE.md caveat numbers that neither pfunc_e3_gate.json nor
review_recompute.json stores (re-score R2.1), from the gate's per_amplicon.csv.gz + the E3 catalogue.

    e3_caveats.py --csv per_amplicon.csv.gz --catalogue e3_clade_catalogue.json [--out J]

exclude_0.987, scored n >= 2 amplicons. Pooled union Brier difference E3(species) - k/n per stratum,
from per-amplicon Brier x union_size (cells); 95% percentile CI from the gate's cluster bootstrap
(calibrate_pfunc.boot_weights: default_rng(1729).integers(0, K, (2000, K)) over the stratum's sorted
cluster keys). Strata: n = 2, n = 3, n >= 4, n >= 10; amplicons with no neighbour weight lost to
missing or empty PGFam sets (w_missing_sp == 0 and n_nb_empty_sp == 0); clusters = source genome
(the gate's unit), source family, source order (catalogue genome -> taxon fallback by prefix taxid).
"""
import argparse, csv, gzip, json
import numpy as np


def boot_weights(K, B, seed):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, K, size=(B, K))
    flat = (idx + (np.arange(B, dtype=np.int64) * K)[:, None]).ravel()
    return np.bincount(flat, minlength=B * K).reshape(B, K).astype(np.float64)


def stat(rows, cluster_of, B=2000, seed=1729):
    if not rows:
        return None
    c = np.array([r["c"] for r in rows]); dk = np.array([r["d"] for r in rows])
    point = dk.sum() / c.sum()
    keys = sorted({cluster_of(r) for r in rows})
    pos = {k: i for i, k in enumerate(keys)}
    inv = np.array([pos[cluster_of(r)] for r in rows])
    K = len(keys)
    C = np.bincount(inv, weights=c, minlength=K); D = np.bincount(inv, weights=dk, minlength=K)
    W = boot_weights(K, B, seed)
    bs = (W @ D) / (W @ C)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return {"amplicons": len(rows), "clusters": K, "dbrier": round(float(point), 6),
            "ci95": [round(float(lo), 6), round(float(hi), 6)], "ci_width": round(float(hi - lo), 6)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True); ap.add_argument("--catalogue", required=True); ap.add_argument("--out")
    a = ap.parse_args()
    cat = json.load(open(a.catalogue))
    G, T = cat["genomes"], cat["taxa"]

    def ranks(g):
        v = G.get(g)
        if v is None:
            v = T.get(g.split(".")[0])
        return (None, None, None) if v is None else tuple(v)
    rows = []
    with gzip.open(a.csv, "rt", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["mode"] != "exclude_0.987" or r["scored"] != "1":
                continue
            n = int(r["n"])
            if n < 2:
                continue
            c = int(r["union_size"])
            rows.append({"src": r["src_genome_id"], "n": n, "c": c,
                         "d": (float(r["brier_e3s"]) - float(r["brier_kn"])) * c,
                         "lost0": float(r["w_missing_sp"]) == 0 and int(r["n_nb_empty_sp"]) == 0,
                         "fam": ranks(r["src_genome_id"])[1], "ord": ranks(r["src_genome_id"])[2]})
    src = lambda r: r["src"]
    fam = lambda r: r["fam"] if r["fam"] else "src:" + r["src"]
    ordr = lambda r: r["ord"] if r["ord"] else "src:" + r["src"]
    out = {"n_amplicons_n2plus": len(rows)}
    out["all_n2plus"] = stat(rows, src)
    out["n=2"] = stat([r for r in rows if r["n"] == 2], src)
    out["n=3"] = stat([r for r in rows if r["n"] == 3], src)
    out["n>=4"] = stat([r for r in rows if r["n"] >= 4], src)
    out["n>=10"] = stat([r for r in rows if r["n"] >= 10], src)
    out["no_weight_lost"] = stat([r for r in rows if r["lost0"]], src)
    out["cluster_family"] = stat(rows, fam)
    out["cluster_order"] = stat(rows, ordr)
    w0 = out["all_n2plus"]["ci_width"]
    out["ci_width_ratio_family"] = round(out["cluster_family"]["ci_width"] / w0, 3)
    out["ci_width_ratio_order"] = round(out["cluster_order"]["ci_width"] / w0, 3)
    out["n_sources_without_family"] = len({r["src"] for r in rows if not r["fam"]})
    out["n_sources_without_order"] = len({r["src"] for r in rows if not r["ord"]})
    s = json.dumps(out, indent=1)
    print(s)
    if a.out:
        open(a.out, "w").write(s + "\n")


if __name__ == "__main__":
    main()
