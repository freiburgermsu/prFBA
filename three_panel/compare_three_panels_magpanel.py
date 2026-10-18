#!/usr/bin/env python
"""Three-panel MAG comparison on the PANEL-AS-SKETCH arm (taxonomy + gene capture).

In this arm each panel's synthetic MAGs are searched against that panel's own 10,000-genome
sketch at -n 500, contamination 0, tau_ANI = 95. Its RESOLVE RATE is not comparable across panels
-- the reference IS the panel, so a less diverse panel is a shallower reference (5,048 vs 10,000
distinct taxids for k-means vs original) -- and is reported per panel only. Use the GTDB arm
(compare_three_panels_mag.py) for a comparable resolve rate. What IS comparable is the
CONDITIONAL gene-capture F1: how much of the true gene content the selected reference recovers,
given that one was found.

Cell semantics are aggregate.py's, copied from two_panel/compare_mag_panels.load_cells:
  * the population denominator is outcome in {ok, abstained}; truth-unavailable cells are dropped
  * a cell is resolved iff outcome == ok AND precision is not blank (blank precision = pred_empty)
"""
import collections, glob, json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_three_panels import PANELS, FULL_LEN, B, SEED, load_truth

ROOT = os.path.expanduser("~/Documents/prFBA")
CELLS = {
    "original": ROOT + "/mag_validation/out_batches_n500/b*/mag_cells.csv",
    "kmeans":   ROOT + "/mag_validation_km/out_batches/b*/mag_cells.csv",
    "random":   ROOT + "/mag_validation_rnd/out_batches/b*/mag_cells.csv",
}


def load_cells(pat):
    d = pd.concat([pd.read_csv(f, dtype={"genome_id": str}) for f in sorted(glob.glob(pat))], ignore_index=True)
    d = d[d.contamination.astype(float) == 0].copy()
    d = d[d.outcome.isin(("ok", "abstained"))].copy()
    d["resolved"] = d.outcome.eq("ok") & d.precision.notna()
    return d


def main():
    rng = np.random.default_rng(SEED)
    T = {p: load_truth(r) for p, r in PANELS.items()}
    C = {}
    for p, pat in CELLS.items():
        d = load_cells(pat)
        d["Phylum"] = d.genome_id.map(lambda g: (T[p].get(g) or {}).get("Phylum"))
        d["best_len"] = d.genome_id.map(lambda g: (T[p].get(g) or {}).get("best_len") or 0)
        C[p] = d
        print("[load] %-8s %6d cells over %5d genomes, %d batch files" % (p, len(d), d.genome_id.nunique(), len(glob.glob(pat))), flush=True)

    phy = {}
    for p in CELLS:
        g = C[p].drop_duplicates("genome_id")
        g = g[(g.best_len >= FULL_LEN) & g.Phylum.notna()]
        phy[p] = collections.Counter(g.Phylum)
    common = set.intersection(*(set(v) for v in phy.values()))
    w = {f: sum(phy[p][f] for p in CELLS) for f in common}
    print("[std] common phyla %d" % len(common), flush=True)

    def std_f1(sub):
        num = den = 0.0
        for f, wt in w.items():
            v = sub.loc[(sub.Phylum == f) & sub.resolved, "f1"]
            if len(v):
                num += wt * v.mean(); den += wt
        return num / den if den else float("nan")

    def boot_f1(sub):
        groups = {f: sub.loc[(sub.Phylum == f) & sub.resolved, "f1"].to_numpy() for f in w}
        groups = {f: v for f, v in groups.items() if len(v)}
        wsum = sum(w[f] for f in groups)
        out = np.empty(B)
        for b in range(B):
            out[b] = sum(w[f] * v[rng.integers(0, len(v), len(v))].mean() for f, v in groups.items()) / wsum
        return out

    res = dict(meta=dict(n_common_phyla=len(common), depth=500, contamination=0, tau_ani=95.0, bootstrap=B,
                         note="resolve_rate is per-panel only (reference = the panel); cond_f1 is comparable"),
               cells=[])
    comps = sorted(set.intersection(*(set(C[p].completeness.astype(float)) for p in CELLS)), reverse=True)
    for regime in ("include", "exclude"):
        for comp in comps:
            cell = dict(regime=regime, completeness=comp)
            boots = {}
            for p in CELLS:
                s = C[p][(C[p].regime == regime) & (C[p].completeness.astype(float) == comp)]
                c = s[s.best_len >= FULL_LEN]
                r = s.resolved
                boots[p] = boot_f1(c)
                cell[p] = dict(n=int(len(s)), resolve_rate_NOT_COMPARABLE=float(r.mean()),
                               cond_f1_crude=float(s.f1[r].mean()) if r.any() else None,
                               cond_f1_common=float(c.f1[c.resolved].mean()) if c.resolved.any() else None,
                               cond_f1_phystd=float(std_f1(c)),
                               cond_f1_phystd_ci=[float(np.percentile(boots[p], 2.5)), float(np.percentile(boots[p], 97.5))])
            for x, y in (("kmeans", "original"), ("random", "original"), ("random", "kmeans")):
                dd = boots[x] - boots[y]
                lo, hi = np.percentile(dd, [2.5, 97.5])
                cell["delta_cond_f1_%s_minus_%s" % (x, y)] = dict(
                    est=round(cell[x]["cond_f1_phystd"] - cell[y]["cond_f1_phystd"], 4),
                    ci=[round(float(lo), 4), round(float(hi), 4)], excludes_zero=bool(lo > 0 or hi < 0))
            res["cells"].append(cell)
    out = ROOT + "/three_panel/three_panel_mag_panelarm.json"
    json.dump(res, open(out, "w"), indent=1)
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()
