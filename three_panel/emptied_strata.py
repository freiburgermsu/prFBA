#!/usr/bin/env python
"""Emptied-record strata: how much of each panel's self-exclusion abstention is the
record running out (tab:emptied), at any record depth.

Promoted from the ad hoc scratchpad em2.py, which hard-coded a 20-hit record
(`'20' if n >= 20`). Here the record depth is read per panel from
<root>/data/hits/raw/method_selection.json (params.k2) and the strata follow it:

  depth 20   0-4 | 5-9 | 10-14 | 15-19 | 20 (emptied)        (em2.py's bins, reproduced)
  depth 100  0-19 | 20-49 | 50-99 | 100 (emptied)
  other D    0-19 | 20-49 | 50-(D-1) | D (emptied)

Population: FullLength16S, self_mode exclude_0.987, outcome != truth_empty, one row per
organism from gene_capture_cells.csv, binned on n_removed_self (hits removed by
self-exclusion). "Emptied" = n_removed_self >= depth, i.e. the whole record was removed;
this assumes every record holds exactly k2 hits, which holds for the full reference
database (the 500-deep shortlist always fills, so len(record) = min(k2, 500)).

Per panel and stratum: n, share of organisms, non-ok rate (outcome != ok) and strict
abstention rate (outcome == abstained). Survivors (all non-emptied strata) get a Wilson
CI; panels are compared on survivors (pairwise chi2 + Fisher, 3-way chi2), per stratum
(Fisher, original vs each), with a CMH test stratified on the survivor strata and a
directly standardised non-ok rate over the pooled survivor-stratum mix.

    python three_panel/emptied_strata.py                       # cited depth-20 trees
    python three_panel/emptied_strata.py --suffix _d100        # depth-100, exact-tie union
    python three_panel/emptied_strata.py --suffix _d100_legacy # depth-100, legacy reducer

Writes --out (default three_panel/emptied_strata{suffix}.json); `emptied_pct` per panel
is the EMPTIED constant of presentations/make_deck.py.
"""
import argparse, collections, csv, json, os, sys

from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.proportion import proportion_confint

ROOT = os.path.expanduser("~/Documents/prFBA")
PANELS = {"original": ROOT + "/region_validation",
          "kmeans": ROOT + "/region_validation_kmeans",
          "random": ROOT + "/region_validation_random"}
REGION, REGIME = "FullLength16S", "exclude_0.987"


def record_depth(root):
    p = os.path.join(root, "data", "hits", "raw", "method_selection.json")
    try:
        return int(json.load(open(p))["params"]["k2"])
    except (OSError, KeyError, ValueError) as exc:
        raise SystemExit("[emptied] cannot read the record depth from %s (%s)" % (p, exc))


def strata_for(depth):
    """[(label, lo, hi_inclusive)], last one = emptied (n >= depth)."""
    if depth == 20:
        edges = [0, 5, 10, 15]
    else:
        edges = [e for e in (0, 20, 50) if e < depth]
    out = [("%d-%d" % (lo, (edges[i + 1] if i + 1 < len(edges) else depth) - 1), lo,
            (edges[i + 1] if i + 1 < len(edges) else depth) - 1) for i, lo in enumerate(edges)]
    return out + [("%d" % depth, depth, 10**9)]


def load(root, strata):
    path = os.path.join(root, "data", "gene_capture_cells.csv")
    bins = collections.OrderedDict((lab, [0, 0, 0]) for lab, _, _ in strata)   # n, non-ok, abstained
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["region"] != REGION or r["self_mode"] != REGIME or r["outcome"] == "truth_empty":
                continue
            n = int(r["n_removed_self"])
            lab = next(l for l, lo, hi in strata if lo <= n <= hi)
            b = bins[lab]
            b[0] += 1; b[1] += r["outcome"] != "ok"; b[2] += r["outcome"] == "abstained"
    return bins


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--suffix", default="", help="append to each default panel root (_d100, _d100_legacy)")
    ap.add_argument("--roots", default=None, help="original=DIR,kmeans=DIR,random=DIR")
    ap.add_argument("--out", default=None, help="default three_panel/emptied_strata{suffix}.json")
    a = ap.parse_args()
    roots = {p: r + a.suffix for p, r in PANELS.items()}
    if a.roots:
        roots = dict(kv.split("=", 1) for kv in a.roots.split(","))
        if set(roots) != set(PANELS):
            raise SystemExit("--roots must name exactly %s" % sorted(PANELS))
        if a.out is None and not a.suffix:
            raise SystemExit("--roots without --suffix: pass --out explicitly")
    out_path = a.out or os.path.join(ROOT, "three_panel", "emptied_strata%s.json" % a.suffix)

    depth = {p: record_depth(r) for p, r in roots.items()}
    if len(set(depth.values())) != 1:
        raise SystemExit("[emptied] panels have different record depths: %s" % depth)
    D = next(iter(depth.values()))
    strata = strata_for(D)
    labs = [l for l, _, _ in strata]
    surv_labs, empt = labs[:-1], labs[-1]
    print("[emptied] record depth %d; strata %s; roots %s" % (D, labs, roots), flush=True)

    res = dict(meta=dict(depth=D, strata=labs, emptied_stratum=empt, region=REGION, regime=REGIME,
                         roots=roots, population="outcome != truth_empty, one row per organism"),
               panels={}, tests={})
    B = {}
    for p, r in roots.items():
        bins = load(r, strata); B[p] = bins
        N = sum(v[0] for v in bins.values())
        rows = {}
        print("%s  N=%d" % (p, N))
        for lab in labs:
            n, a_, s = bins[lab]
            rows[lab] = dict(n=n, share=round(n / N, 4) if N else None,
                             non_ok_rate=round(a_ / n, 4) if n else None,
                             strict_abstain_rate=round(s / n, 4) if n else None)
            print("  %-7s n=%-6d share=%.3f non-ok=%s strict=%s" % (
                lab, n, n / N if N else 0, "%.3f" % (a_ / n) if n else "n/a",
                "%.3f" % (s / n) if n else "n/a"))
        sn = sum(bins[l][0] for l in surv_labs); sa = sum(bins[l][1] for l in surv_labs)
        ss = sum(bins[l][2] for l in surv_labs)
        lo, hi = proportion_confint(sa, sn, method="wilson") if sn else (None, None)
        print("  survivors n=%d non-ok=%d (%.4f, Wilson %.4f-%.4f) strict=%d (%.4f)"
              % (sn, sa, sa / sn, lo, hi, ss, ss / sn) if sn else "  survivors n=0")
        res["panels"][p] = dict(N=N, strata=rows,
                                emptied_n=bins[empt][0],
                                emptied_pct=round(100 * bins[empt][0] / N, 1) if N else None,
                                survivors=dict(n=sn, non_ok=sa, non_ok_rate=round(sa / sn, 4) if sn else None,
                                               wilson=[lo, hi], strict=ss,
                                               strict_rate=round(ss / sn, 4) if sn else None))

    sv = {p: (sum(B[p][l][1] for l in surv_labs), sum(B[p][l][0] for l in surv_labs)) for p in B}
    pair = {}
    for x, y in (("original", "kmeans"), ("original", "random"), ("kmeans", "random")):
        t = [[sv[x][0], sv[x][1] - sv[x][0]], [sv[y][0], sv[y][1] - sv[y][0]]]
        pair["%s_vs_%s" % (x, y)] = dict(chi2_p=round(float(chi2_contingency(t)[1]), 4),
                                         fisher_p=round(float(fisher_exact(t)[1]), 4))
        print(x, y, pair["%s_vs_%s" % (x, y)])
    res["tests"]["survivors_pairwise"] = pair
    res["tests"]["survivors_3way_chi2_p"] = float(chi2_contingency([[a_, n - a_] for a_, n in sv.values()])[1])
    print("3way", res["tests"]["survivors_3way_chi2_p"])
    per = {}
    for lab in surv_labs:
        for y in ("kmeans", "random"):
            o, r = B["original"][lab], B[y][lab]
            per["%s|original_vs_%s" % (lab, y)] = round(float(
                fisher_exact([[o[1], o[0] - o[1]], [r[1], r[0] - r[1]]])[1]), 3)
            print(lab, y, per["%s|original_vs_%s" % (lab, y)])
    res["tests"]["per_stratum_fisher_p"] = per
    cmh = {}
    for y in ("kmeans", "random"):
        tabs = [[[B["original"][l][1], B["original"][l][0] - B["original"][l][1]],
                 [B[y][l][1], B[y][l][0] - B[y][l][1]]] for l in surv_labs]
        st = StratifiedTable(tabs)
        cmh["original_vs_%s" % y] = dict(odds_ratio_pooled=float(st.oddsratio_pooled),
                                         p=float(st.test_null_odds().pvalue))
        print("CMH orig vs", y, st.oddsratio_pooled, st.test_null_odds().pvalue)
    res["tests"]["cmh_survivor_strata"] = cmh
    tot = {l: sum(B[p][l][0] for p in B) for l in surv_labs}; T = sum(tot.values())
    std = {p: sum(tot[l] / T * B[p][l][1] / B[p][l][0] for l in surv_labs if B[p][l][0])
           for p in B}
    for p in B:
        print("std", p, std[p])
    res["tests"]["standardised_non_ok_rate"] = std
    json.dump(res, open(out_path, "w"), indent=1)
    print("[emptied] emptied_pct %s -> %s" % ({p: v["emptied_pct"] for p, v in res["panels"].items()},
                                              out_path))


if __name__ == "__main__":
    main()
