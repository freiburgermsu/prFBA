#!/usr/bin/env python
"""R4.1: why the phenotype genome samples changed. Rebuilds phenotype_fold's candidate list (one amplicon per genome,
truth >= 200 ECs, both methods scored) from the PRE-REPAIR and the NEW inputs, by capturing the arguments
build_jobs passes to choose_genomes, and reports, per pass size (1000 ATP, 200 gapfilled):
  candidates old/new, candidates whose delta EC F1 changed, strata (domain x template x delta-F1 quintile) whose
  membership changed, and how many sampled genomes come from strata with identical membership (those picks are
  identical by construction: the per-stratum shuffle is seeded by the stratum key and sees the same list).
usage: r4_pheno_sample_diff.py OUT_JSON
"""
import gc
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, "/home/freiburger/Documents/prFBA_rescore_2026-11/tools")
import r4_pheno_attrib as A  # noqa: E402  (paths, old picrust2 root)

PF = A.PF


def cands(cells, functions, prfba, picrust2):
    cap = {}
    real = PF.choose_genomes

    def spy(c, n_target, seed, experimental):
        cap["c"] = [dict(x) for x in c]
        return real(c, n_target, seed, experimental)
    PF.choose_genomes = spy
    try:
        A.jobs_for("gf", cells, functions, prfba, picrust2)
    finally:
        PF.choose_genomes = real
    gc.collect()
    return cap["c"]


def strata(cs, n_target):
    picked = {c["genome_id"] for c in PF.choose_genomes([dict(x) for x in cs], n_target, 1729, None)}
    ordered = sorted(cs, key=lambda c: ((c["delta_f1"] if c["delta_f1"] is not None else 0.0), c["genome_id"]))
    st = defaultdict(set)
    for i, c in enumerate(ordered):
        st[(c["domain"], c["template"], PF.quintile(i, len(ordered)))].add(c["genome_id"])
    return st, picked


p2_old, _ = A.old_picrust2_root()
old = cands(f"{A.PRS}/data/cells.csv.gz", f"{A.PRP}/data/genome_function_sets.json", f"{A.PRS}/data/prfba", p2_old)
new = cands(f"{A.SR}/data/cells.csv.gz",
            f"{A.R}/region_validation/picrust2_comparison/data/genome_function_sets.json",
            f"{A.SR}/data/prfba", f"{A.SR}/data/picrust2")
oc, nc = {c["genome_id"]: c for c in old}, {c["genome_id"]: c for c in new}
both = set(oc) & set(nc)
res = {"n_candidates": [len(oc), len(nc)], "only_old": len(set(oc) - set(nc)), "only_new": len(set(nc) - set(oc)),
       "delta_f1_changed": sum(1 for g in both if oc[g]["delta_f1"] != nc[g]["delta_f1"]),
       "amplicon_changed": sum(1 for g in both if (oc[g]["region"], oc[g]["md5"]) != (nc[g]["region"], nc[g]["md5"]))}
for n in (1000, 200):
    so, po = strata(old, n)
    sn, pn = strata(new, n)
    keys = set(so) | set(sn)
    same = {k for k in keys if so.get(k) == sn.get(k)}
    res[f"n{n}"] = {"strata": len(keys), "strata_identical_membership": len(same),
                    "picked_both": len(po & pn), "picked_only_old": len(po - pn),
                    "picked_from_identical_strata_new": sum(1 for k in same for g in sn[k] if g in pn),
                    "picked_from_identical_strata_equal": all((so[k] & po) == (sn[k] & pn) for k in same)}
json.dump(res, open(sys.argv[1], "w"), indent=1)
print(json.dumps(res, indent=1))
