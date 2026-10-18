#!/usr/bin/env python
"""Three-panel MAG comparison on the GTDB dense-reference arm.

The GTDB arm is the cross-panel-comparable MAG measurement: every panel's synthetic MAGs are
searched against the same 732,475-accession GTDB r226 sketch, so reference depth is identical
and only the query panel differs. (The panel-as-sketch arm's resolve rate is NOT comparable:
each panel is its own reference, so a less diverse panel is a shallower reference.)

Matched on everything but the panel: same synth_mags.py and seed, same skani 0.3.2 binary and
sketch, -n 100 --min-af 0, contamination 0, completeness {100,90,70,50,30}%, the same 20%
proportional-by-phylum subsample rule, and score_gtdb.py's exclusion rule (ANI >= 95).

Ladder as in compare_three_panels.py: crude -> common support (full-length best 16S) ->
domain-std -> phylum-std. Only the phylum-standardised column is attributable to the selection
scheme.
"""
import collections, json, math, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_three_panels import PANELS, FULL_LEN, B, SEED, load_truth

ROOT = os.path.expanduser("~/Documents/prFBA")
CELLS = {
    "original": ROOT + "/mag_validation/out_gtdb/gtdb_cells_2000.csv",
    "kmeans":   ROOT + "/mag_validation_km/out_gtdb/gtdb_cells.csv",
    "random":   ROOT + "/mag_validation_rnd/out_gtdb/gtdb_cells.csv",
}
COMPL = [1.0, 0.9, 0.7, 0.5, 0.3]


def main():
    rng = np.random.default_rng(SEED)
    T = {p: load_truth(r) for p, r in PANELS.items()}
    C = {}
    for p, f in CELLS.items():
        d = pd.read_csv(f, dtype={"genome_id": str})
        d = d[d.contamination.astype(float) == 0]
        d["Phylum"] = d.genome_id.map(lambda g: (T[p].get(g) or {}).get("Phylum"))
        d["domain"] = d.genome_id.map(lambda g: (T[p].get(g) or {}).get("domain"))
        d["best_len"] = d.genome_id.map(lambda g: (T[p].get(g) or {}).get("best_len") or 0)
        C[p] = d
        print("[load] %-8s %5d cells, %4d genomes, %4d with truth lineage"
              % (p, len(d), d.genome_id.nunique(), d.Phylum.notna().groupby(d.genome_id).first().sum()), flush=True)

    def counts(p, key):
        g = C[p].drop_duplicates("genome_id")
        g = g[(g.best_len >= FULL_LEN) & g[key].notna()]
        return collections.Counter(g[key])
    phy = {p: counts(p, "Phylum") for p in CELLS}
    dom = {p: counts(p, "domain") for p in CELLS}
    common_phy = set.intersection(*(set(v) for v in phy.values()))
    common_dom = set.intersection(*(set(v) for v in dom.values()))
    wphy = {f: sum(phy[p][f] for p in CELLS) for f in common_phy}
    wdom = {d: sum(dom[p][d] for p in CELLS) for d in common_dom}
    cover = {p: round(sum(phy[p][f] for f in common_phy) / max(1, sum(phy[p].values())), 4) for p in CELLS}
    print("[std] common phyla %d; share of each panel's genomes inside them %s" % (len(common_phy), cover), flush=True)

    def std(sub, key, w):
        num = den = 0.0
        for s, wt in w.items():
            v = sub.loc[sub[key] == s, "resolves"]
            if len(v):
                num += wt * v.mean(); den += wt
        return num / den if den else float("nan")

    def boot(sub):
        groups = {f: sub.loc[sub.Phylum == f, "resolves"].to_numpy() for f in wphy if (sub.Phylum == f).any()}
        wsum = sum(wphy[f] for f in groups)
        out = np.empty(B)
        for b in range(B):
            out[b] = sum(wphy[f] * v[rng.integers(0, len(v), len(v))].mean() for f, v in groups.items()) / wsum
        return out

    res = dict(meta=dict(n_common_phyla=len(common_phy), std_coverage=cover, bootstrap=B, seed=SEED,
                         n_genomes={p: int(C[p].genome_id.nunique()) for p in CELLS}), cells=[])
    for regime in ("include", "exclude"):
        for comp in COMPL:
            cell = dict(regime=regime, completeness=comp)
            boots = {}
            for p in CELLS:
                s = C[p][(C[p].regime == regime) & (C[p].completeness.astype(float) == comp)]
                c = s[s.best_len >= FULL_LEN]
                boots[p] = boot(c)
                cell[p] = dict(n=int(len(s)), crude=float(s.resolves.mean()),
                               common=float(c.resolves.mean()),
                               domstd=float(std(c, "domain", wdom)), phystd=float(std(c, "Phylum", wphy)),
                               phystd_ci=[float(np.percentile(boots[p], 2.5)), float(np.percentile(boots[p], 97.5))],
                               median_best_ani=float(s.best_ani.median()) if s.best_ani.notna().any() else None)
            for x, y in (("kmeans", "original"), ("random", "original"), ("random", "kmeans")):
                dd = boots[x] - boots[y]
                lo, hi = np.percentile(dd, [2.5, 97.5])
                cell["delta_%s_minus_%s" % (x, y)] = dict(est=round(cell[x]["phystd"] - cell[y]["phystd"], 4),
                                                         ci=[round(float(lo), 4), round(float(hi), 4)],
                                                         excludes_zero=bool(lo > 0 or hi < 0))
            res["cells"].append(cell)
    out = ROOT + "/three_panel/three_panel_mag_gtdb.json"
    json.dump(res, open(out, "w"), indent=1)
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()
