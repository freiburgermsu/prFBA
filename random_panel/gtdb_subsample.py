#!/usr/bin/env python
"""Draw the phylum-stratified GTDB-arm subsample for an alternative panel.

The original panel's mag_validation/gtdb_sub_ids.txt is a 20% PROPORTIONAL allocation
by phylum (verified: Bacillota 2718->543, Pseudomonadota 1032->206, Actinomycetota
861->173, Methanobacteriota 818->163, Bacteroidota 535->107 -- all exactly 0.20, and
160 of 206 phyla are represented because phyla with fewer than 3 members round to 0).
This reproduces that rule so the k-means and random panels get subsamples built the
same way, and the three arms differ only in the panel they are drawn from.

Cross-panel rates from this arm must still be reported phylum-standardized: the three
panels have different phylum compositions by construction, so a proportional subsample
inherits that difference and only the standardized column is attributable to the
selection scheme.
"""
import argparse, collections, json
import numpy as np

SEED = 1729


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--frac", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=SEED)
    a = ap.parse_args()

    genomes = json.load(open(a.selection))["genomes"]
    by = collections.defaultdict(list)
    for e in genomes:
        by[e["phylum"]].append(str(e["genome_id"]))

    rng = np.random.default_rng(a.seed)
    picked = []
    for ph in sorted(by):                      # sorted => reproducible
        ids = sorted(by[ph])
        k = int(round(a.frac * len(ids)))       # proportional allocation, banker-free round
        if k <= 0:
            continue
        idx = rng.choice(len(ids), size=min(k, len(ids)), replace=False)
        picked.extend(ids[i] for i in sorted(idx))

    picked = sorted(set(picked))
    with open(a.out, "w") as fh:
        fh.write("".join(g + "\n" for g in picked))
    cov = sum(1 for ph in by if int(round(a.frac * len(by[ph]))) > 0)
    print("[sub] %s: %d genomes from %d phyla (%d/%d phyla represented) -> %s"
          % (a.selection.rsplit("/", 1)[-1], len(picked), len(by), cov, len(by), a.out))


if __name__ == "__main__":
    main()
