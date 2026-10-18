#!/usr/bin/env python
"""Cross-panel MAG comparison: k-means panel vs the taxonomy-stratified original.

Matched on everything that is not the selection scheme:
  * skani retention depth -n 500 on BOTH arms (the original's main run was -n 100;
    out_batches_n500 is the matched re-run)
  * contamination 0 on both (the k-means arm was synthesized contamination-free)
  * the same exclusion rule, tau_ANI = 95

The ladder is the same as compare_panels.py's: crude -> domain-standardized ->
phylum-standardized. Only the phylum-standardized column is attributable to the
selection scheme; the crude column also carries the panels' different taxon mixes.
"""
import collections, glob, json, os, sys
import pandas as pd

sys.path.insert(0, os.path.expanduser("~/Documents/prFBA/two_panel"))
from compare_panels import load_truth, standardize

ROOT = os.path.expanduser("~/Documents/prFBA")
OUT = ROOT + "/two_panel"


def load_cells(pat):
    """Load MAG cells using aggregate.py's EXACT definitions, not an approximation.

    Two traps, both of which this function previously fell into:
      * a cell with outcome 'ok' and a BLANK precision is a pred_empty -- a
        reference was selected but predicted no gene families. aggregate.py
        excludes it from the answered set; counting it as resolved reproduces the
        very bug the score_mag.py audit fixed, and deflates conditional F1
        (0.797 instead of 0.807 on the original panel).
      * the coverage denominator is outcome in {ok, abstained}, NOT every cell.
        Cells whose truth is unavailable (a gene-less source genome) are dropped
        from both numerator and denominator. The original panel has 161 gene-less
        sources against the k-means panel's 45, so keeping them penalises the
        original roughly three times as hard -- an asymmetry larger than the
        coverage difference being measured.
    """
    d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(pat))], ignore_index=True)
    d = d[d.contamination.astype(float) == 0].copy()
    d = d[d.outcome.isin(("ok", "abstained"))].copy()          # the pop denominator
    d["resolved"] = d.outcome.eq("ok") & d.precision.notna()   # blank precision = pred_empty
    return d


def ladder(d, truth, wdom, wphy):
    """crude / domain-std / phylum-std resolve rate, plus conditional + population F1."""
    def strata(keyfn):
        agg = collections.defaultdict(lambda: [0, 0])
        for g, ok in zip(d.genome_id.astype(str), d.resolved):
            k = keyfn(g)
            if k is None:
                continue
            agg[k][1] += 1
            agg[k][0] += int(ok)
        return agg
    dom = strata(lambda g: (truth.get(g) or {}).get("domain"))
    phy = strata(lambda g: ((truth.get(g) or {}).get("truth_lineage") or {}).get("Phylum"))
    r = d.resolved
    return dict(n=int(len(d)), crude=float(r.mean()),
                domstd=standardize(dom, wdom), phystd=standardize(phy, wphy),
                condF1=float(d.f1[r].mean()) if r.any() else None,
                popF1=float(d.f1.where(r, 0).mean()),
                n_resolved=int(r.sum()))


def main():
    ot = load_truth(ROOT + "/region_validation/data/benchmark_truth.json")
    kt = load_truth(ROOT + "/region_validation_kmeans/data/benchmark_truth.json")
    ophy = collections.Counter(v["truth_lineage"].get("Phylum") for v in ot.values())
    kphy = collections.Counter(v["truth_lineage"].get("Phylum") for v in kt.values())
    common_phyla = {p for p in (set(ophy) & set(kphy)) if p}
    wphy = {p: ophy[p] + kphy[p] for p in common_phyla}
    odom = collections.Counter(v.get("domain") for v in ot.values())
    kdom = collections.Counter(v.get("domain") for v in kt.values())
    common_dom = {x for x in (set(odom) & set(kdom)) if x}
    wdom = {x: odom[x] + kdom[x] for x in common_dom}
    print("common phyla %d (orig %d, km %d); domains %s"
          % (len(common_phyla), len(ophy), len(kphy), sorted(common_dom)), flush=True)

    og = load_cells(ROOT + "/mag_validation/out_batches_n500/b*/mag_cells.csv")
    km = load_cells(ROOT + "/mag_validation_km/out_batches/b*/mag_cells.csv")
    print("cells: original %d, k-means %d" % (len(og), len(km)), flush=True)

    res = {"meta": dict(n_common_phyla=len(common_phyla), depth=500, contamination=0,
                        tau_ani=95.0, n_orig_cells=int(len(og)), n_km_cells=int(len(km)))}
    rows = []
    for regime in ("exclude", "include"):
        for comp in sorted(set(og.completeness) & set(km.completeness), reverse=True):
            o = ladder(og[(og.regime == regime) & (og.completeness == comp)], ot, wdom, wphy)
            k = ladder(km[(km.regime == regime) & (km.completeness == comp)], kt, wdom, wphy)
            rows.append(dict(regime=regime, completeness=float(comp),
                             **{"orig_" + a: b for a, b in o.items()},
                             **{"km_" + a: b for a, b in k.items()}))
            f = lambda x: "  n/a" if x is None else "%.4f" % x
            print("%-8s c=%.2f  coverage crude %s/%s  dom-std %s/%s  phy-std %s/%s  condF1 %s/%s  popF1 %s/%s"
                  % (regime, comp, f(k["crude"]), f(o["crude"]), f(k["domstd"]), f(o["domstd"]),
                     f(k["phystd"]), f(o["phystd"]), f(k["condF1"]), f(o["condF1"]),
                     f(k["popF1"]), f(o["popF1"])), flush=True)
    res["rows"] = rows
    json.dump(res, open(OUT + "/mag_panel_comparison.json", "w"), indent=1)
    print("wrote mag_panel_comparison.json  (k-means / original in each pair)")


if __name__ == "__main__":
    main()
