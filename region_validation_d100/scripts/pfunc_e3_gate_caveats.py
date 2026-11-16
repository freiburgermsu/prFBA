#!/usr/bin/env python
"""Caveat numbers of PFUNC_E3_GATE.md that neither data/pfunc_e3_gate.json nor
data/pfunc_e3_gate/review_recompute.json stores: the by-n bootstrap CIs, n >= 10, the
clustering-unit sensitivity (source family / order) and the no-lost-neighbour-weight subset.

    python region_validation_d100/scripts/pfunc_e3_gate_caveats.py
        [--csv data/pfunc_e3_gate/per_amplicon.csv.gz] [--catalogue bvbrc_cache/e3_clade_catalogue.json]
        [--out data/pfunc_e3_gate/caveats_recompute.json]

Scope: exclude_0.987, scored amplicons with n >= 2 (the gate's primary cell). Per stratum, the pooled
union Brier difference E3(species) - k/n = sum((brier_e3s - brier_kn) * union_size) / sum(union_size)
over the gate's per-amplicon rows; 95% percentile CI from the gate's own cluster bootstrap
(calibrate_pfunc.boot_weights: numpy default_rng(1729).integers(0, K, (2000, K)) over the stratum's
sorted cluster keys, each stratum drawn on its own). Strata: n = 2, n = 3, n >= 4, n >= 10; amplicons
with no neighbour weight lost to missing or empty PGFam sets (w_missing_sp == 0 and n_nb_empty_sp == 0).
Cluster units: the source genome (the gate's unit), and the source's family / order from the E3
catalogue (genome entry, else the taxon entry of its prefix taxid; a source without that rank is its
own cluster).

Control (2026-11-15, pre-repair inputs: the committed 2026-11-12 per_amplicon.csv.gz and the
pre-repair catalogue): reproduces the gate's own pooled CI exactly (-0.0530 [-0.0547, -0.0513]), the
11,557-amplicon subset exactly (-0.0488 [-0.0507, -0.0468]), the 9,096 n = 2 amplicons, every by-n
point estimate and the 1.9x family widening. The 2026-11-12 reviewer's ad hoc bootstrap gave
n = 2 [-0.1040, -0.0976], n >= 4 [-0.0254, -0.0226], family [-0.0564, -0.0499] and order 2.7x; this
script's draw gives [-0.1040, -0.0974], [-0.0253, -0.0226], [-0.0562, -0.0496] and 2.6x.
"""
import argparse
import csv
import gzip
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
ROOT = os.path.dirname(os.path.dirname(HERE))
B_BOOT, SEED = 2000, 1729


def boot_weights(K, B, seed):
    """calibrate_pfunc.boot_weights verbatim (the gate's draw)."""
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, K, size=(B, K))
    flat = (idx + (np.arange(B, dtype=np.int64) * K)[:, None]).ravel()
    return np.bincount(flat, minlength=B * K).reshape(B, K).astype(np.float64)


def stat(rows, cluster_of, B=B_BOOT, seed=SEED):
    if not rows:
        return None
    c = np.array([r["c"] for r in rows], np.float64)
    dk = np.array([r["d"] for r in rows], np.float64)
    keys = sorted({cluster_of(r) for r in rows})
    pos = {k: i for i, k in enumerate(keys)}
    inv = np.array([pos[cluster_of(r)] for r in rows])
    K = len(keys)
    C = np.bincount(inv, weights=c, minlength=K)
    D = np.bincount(inv, weights=dk, minlength=K)
    W = boot_weights(K, B, seed)
    lo, hi = np.percentile((W @ D) / (W @ C), [2.5, 97.5])
    return {"amplicons": len(rows), "clusters": K, "dbrier": round(float(dk.sum() / c.sum()), 6),
            "ci95": [round(float(lo), 6), round(float(hi), 6)], "ci_width": round(float(hi - lo), 6)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=os.path.join(DATA, "pfunc_e3_gate", "per_amplicon.csv.gz"))
    ap.add_argument("--catalogue", default=os.path.join(ROOT, "bvbrc_cache", "e3_clade_catalogue.json"))
    ap.add_argument("--out", default=os.path.join(DATA, "pfunc_e3_gate", "caveats_recompute.json"))
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
            if r["mode"] != "exclude_0.987" or r["scored"] != "1" or int(r["n"]) < 2:
                continue
            c = int(r["union_size"])
            _, fam, order = ranks(r["src_genome_id"])
            rows.append({"src": r["src_genome_id"], "n": int(r["n"]), "c": c,
                         "d": (float(r["brier_e3s"]) - float(r["brier_kn"])) * c,
                         "lost0": float(r["w_missing_sp"]) == 0 and int(r["n_nb_empty_sp"]) == 0,
                         "fam": fam, "ord": order})
    src = lambda r: r["src"]  # noqa: E731
    out = {"_meta": {"what": "PFUNC_E3_GATE.md caveat numbers not stored by the gate or the review "
                             "(exclude_0.987, scored n >= 2, E3 species - k/n pooled union Brier)",
                     "csv": os.path.abspath(a.csv), "catalogue": os.path.abspath(a.catalogue),
                     "bootstrap": {"B": B_BOOT, "seed": SEED, "draw": "calibrate_pfunc.boot_weights per stratum"}},
           "n_amplicons_n2plus": len(rows)}
    out["all_n2plus"] = stat(rows, src)
    out["n=2"] = stat([r for r in rows if r["n"] == 2], src)
    out["n=3"] = stat([r for r in rows if r["n"] == 3], src)
    out["n>=4"] = stat([r for r in rows if r["n"] >= 4], src)
    out["n>=10"] = stat([r for r in rows if r["n"] >= 10], src)
    out["no_neighbour_weight_lost"] = stat([r for r in rows if r["lost0"]], src)
    out["cluster_source_family"] = stat(rows, lambda r: r["fam"] or "src:" + r["src"])
    out["cluster_source_order"] = stat(rows, lambda r: r["ord"] or "src:" + r["src"])
    w0 = out["all_n2plus"]["ci_width"]
    out["ci_width_ratio_family"] = round(out["cluster_source_family"]["ci_width"] / w0, 3)
    out["ci_width_ratio_order"] = round(out["cluster_source_order"]["ci_width"] / w0, 3)
    out["n_sources_without_family"] = len({r["src"] for r in rows if not r["fam"]})
    out["n_sources_without_order"] = len({r["src"] for r in rows if not r["ord"]})
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: v for k, v in out.items() if k != "_meta"}, indent=1))


if __name__ == "__main__":
    main()
