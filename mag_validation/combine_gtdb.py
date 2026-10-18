#!/usr/bin/env python
"""Combine the GTDB dense-reference arm over the FULL phylum-stratified 2,000.

The first four batches (b000-b003, 1,000 genomes) were run earlier with the
contamination grid; b004-b007 were run contamination-free. Only contamination 0
is combined, which is the regime both halves share -- and the only one the
contamination claims survived the synthesizer audit.

Why this matters: gtdb_sub_ids.txt is not randomly ordered with respect to phylum.
The first 1,000 ids span 85 phyla and the second 1,000 span 123, so results from
the completed half alone over-weight the common phyla. Every figure quoted from
this arm should come from all 2,000.
"""
import glob, json, os, sys
import numpy as np, pandas as pd

MV = os.path.expanduser("~/Documents/prFBA/mag_validation")
sys.path.insert(0, os.path.expanduser("~/Documents/prFBA/region_validation/scripts"))
from _common import wilson_ci

old = pd.read_csv(MV + "/out_gtdb/gtdb_cells.csv")
new = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(MV + "/out_gtdb_batches/*/gtdb_cells.csv"))],
                ignore_index=True)
old = old[old.contamination.astype(float) == 0]
new = new[new.contamination.astype(float) == 0]
assert not (set(old.genome_id) & set(new.genome_id)), "halves overlap"
d = pd.concat([old, new], ignore_index=True)
print("combined: %d genomes, %d cells (old %d + new %d)"
      % (d.genome_id.nunique(), len(d), old.genome_id.nunique(), new.genome_id.nunique()))

summary = {}
for (comp, reg), g in d.groupby(["completeness", "regime"]):
    n = len(g); r = int(g.resolves.sum())
    lo, hi = wilson_ci(r, n)
    summary["c%g_x0_%s" % (comp, reg)] = dict(
        completeness=float(comp), contamination=0, regime=reg, n=n, resolves=r,
        resolve_rate=r / n, ci=[lo, hi],
        median_best_ani=float(g.best_ani.median()),
        median_af_query=float(g.best_af_query.median()),
        median_af_ref=float(g.best_af_ref.median()),
        median_n_af_qualified=float(g.n_af_qualified.median()),
        mean_n_af_qualified=float(g.n_af_qualified.mean()),
        frac_zero_af_qualified=float((g.n_af_qualified == 0).mean()),
        ani_classes={k: int(v) for k, v in g.ani_class.value_counts().items()})

out = MV + "/out_gtdb/gtdb_summary_2000.json"
json.dump(dict(n_genomes=int(d.genome_id.nunique()), contamination=0, summary=summary),
          open(out, "w"), indent=1)
d.to_csv(MV + "/out_gtdb/gtdb_cells_2000.csv", index=False)
print("\n%-10s %8s %8s   %s" % ("cell", "n", "resolve", "vs the 1,000-genome half"))
for c in sorted(summary, key=lambda k: (summary[k]["regime"], -summary[k]["completeness"])):
    s = summary[c]
    o = old[(old.completeness == s["completeness"]) & (old.regime == s["regime"])]
    print("%-10s %8d %8.4f   %.4f" % (c, s["n"], s["resolve_rate"], o.resolves.mean()))
print("\nwrote", out)
