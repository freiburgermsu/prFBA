#!/usr/bin/env python
"""Is the k-means panel's 16S coverage deficit a retention-depth artifact?

Hypothesis: the k-means panel is NOT species-dereplicated (it samples the embedding
uniformly), so a query often has many same-species panel-mates. The 16S arm keeps
only the top 20 alignment hits per query, so when a query has >= 20 same-species
mates, self-exclusion (tau = 0.987) removes every retained hit and the query
abstains for want of record depth rather than for want of a relative. The MAG arm
keeps 500, so the same query still has survivors.

Prediction: 16S abstention rises steeply with same-species mate count and
approaches 1 at >= 20; MAG abstention does not.

Falsifiable the other way: if 16S abstention is flat in mate count, the deficit is
something else and the manuscript needs a different sentence.
"""
import collections, csv, glob, json, os, sys
import pandas as pd

ROOT = os.path.expanduser("~/Documents/prFBA")
sys.path.insert(0, ROOT + "/two_panel")
from compare_panels import load_truth, load_genecap

OUT = ROOT + "/two_panel"
BINS = [(1, 1), (2, 2), (3, 4), (5, 9), (10, 19), (20, 10**9)]


def mates(truth):
    """genome_id -> number of panel organisms sharing its truth species (incl. itself)."""
    bysp = collections.Counter(
        (v.get("truth_lineage") or {}).get("Species") for v in truth.values())
    return {g: bysp[(v.get("truth_lineage") or {}).get("Species")] for g, v in truth.items()}


def mag_abstain(pat):
    """genome_id -> abstained? at completeness 1.0, self-exclusion, contamination 0."""
    d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(pat))], ignore_index=True)
    d = d[(d.regime == "exclude") & (d.completeness.astype(float) == 1.0)
          & (d.contamination.astype(float) == 0)]
    res = d.outcome.eq("ok") & d.n_selected.gt(0)
    return dict(zip(d.genome_id.astype(str), ~res))


def report(label, truth, genecap, magab):
    m = mates(truth)
    agg = {b: [0, 0, 0, 0] for b in BINS}      # n16, abstain16, nmag, abstainmag
    for g, k in m.items():
        b = next(b for b in BINS if b[0] <= k <= b[1])
        gc = genecap.get(g)
        # TRUE abstention only: no reference survived self-exclusion. Earlier this
        # lumped in pgfam_missing, which was a cache gap (now closed, 0 cells) and
        # not a retrieval failure -- that inflated k-means abstention in every bin.
        if gc is not None and gc[0] != "truth_empty":
            agg[b][0] += 1
            agg[b][1] += int(gc[0] == "abstained")
        ma = magab.get(g)
        if ma is not None:
            agg[b][2] += 1
            agg[b][3] += int(ma)
    print("\n=== %s ===" % label)
    print("%-10s %7s %10s %7s %10s" % ("mates", "n(16S)", "abstain16S", "n(MAG)", "abstainMAG"))
    rows = []
    for b in BINS:
        n1, a1, n2, a2 = agg[b]
        lab = "%d" % b[0] if b[0] == b[1] else ("%d+" % b[0] if b[1] > 10**8 else "%d-%d" % b)
        print("%-10s %7d %10s %7d %10s"
              % (lab, n1, "%.4f" % (a1 / n1) if n1 else "n/a",
                 n2, "%.4f" % (a2 / n2) if n2 else "n/a"))
        rows.append(dict(panel=label, bin=lab, n_16s=n1, abstain_16s=(a1 / n1) if n1 else None,
                         n_mag=n2, abstain_mag=(a2 / n2) if n2 else None))
    tot = sum(1 for v in m.values() if v >= 20)
    print("organisms with >= 20 same-species panel-mates: %d of %d (%.1f%%)"
          % (tot, len(m), 100 * tot / len(m)))
    return rows, tot, len(m)


res = {}
rows = []
for label, tdir, cells in (
        ("k-means", ROOT + "/region_validation_kmeans", ROOT + "/mag_validation_km/out_batches"),
        ("original", ROOT + "/region_validation", ROOT + "/mag_validation/out_batches_n500")):
    truth = load_truth(tdir + "/data/benchmark_truth.json")
    gc = load_genecap(tdir + "/data/gene_capture_cells.csv", "FullLength16S")
    ma = mag_abstain(cells + "/b*/mag_cells.csv")
    r, tot, n = report(label, truth, gc, ma)
    rows += r
    res[label] = dict(n_ge20=tot, n=n)
json.dump(dict(meta=res, rows=rows), open(OUT + "/mechanism_depth.json", "w"), indent=1)
print("\nwrote mechanism_depth.json")
