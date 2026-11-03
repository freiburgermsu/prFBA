"""Phase 0: population and reference-level sanity checks (no estimator math)."""
import collections
import json
import os
import time

import numpy as np

import vctx

t0 = time.time()
C = vctx.Ctx(load_truth=True)
out = collections.OrderedDict()
out["panel_genomes"] = C.NG
out["ec_universe"] = C.NU
out["genomes_without_truth_record"] = C.missing_truth
out["panel_species_labels_distinct"] = len(set(C.lab["species"]))
out["panel_species_norm_distinct"] = len(set(vctx.sp_norm(s) for s in C.lab["species"]))
out["panel_taxon_ids_distinct"] = len(set(C.lab["taxon"]))
out["panel_genomes_without_group"] = sum(1 for g in C.lab["group"] if g is None)
# genera spanning more than one holdout group (genus blocking)
gg = collections.defaultdict(set)
for gen, grp in zip(C.lab["genus"], C.lab["group"]):
    if gen:
        gg[gen].add(grp)
out["genera_spanning_more_than_one_group"] = sum(1 for v in gg.values() if len(v) > 1)
out["genomes_with_missing_genus_label"] = sum(1 for v in C.lab["genus"] if not v)
out["genomes_with_placeholder_genus_label"] = sum(1 for v in C.lab["genus"]
                                                  if v and v.startswith("unclassified_"))
# placeholder genus labels shared by several genomes across groups?
ph = collections.defaultdict(set)
for gen, grp in zip(C.lab["genus"], C.lab["group"]):
    if gen and gen.startswith("unclassified_"):
        ph[gen].add(grp)
out["placeholder_genus_labels"] = len(ph)
out["placeholder_genus_labels_spanning_groups"] = sum(1 for v in ph.values() if len(v) > 1)

cells = vctx.load_cells("1.0")
out["cells"] = len(cells)
out["sources"] = int(cells.src_genome_id.nunique())
out["groups"] = int(cells.group_id.nunique())
nr = cells.n_reps.astype(int)
out["n1"] = int((nr == 1).sum())
out["n2"] = int((nr == 2).sum())
out["n3"] = int((nr == 3).sum())
out["n4+"] = int((nr >= 4).sum())
out["cells_by_fold"] = {int(k): int(v) for k, v in cells.fold.value_counts().sort_index().items()}
out["records_by_fold"] = {int(f): int(s.md5.nunique()) for f, s in cells.groupby("fold")}
# cell group == fold_plan group of the source
out["cell_group_ne_fold_plan_group"] = int(sum(
    1 for s, g in zip(cells.src_genome_id, cells.group_id) if C.ggroup.get(s) != g))
# truth size equals cells n_truth
T = C.E.sum(1)
out["n_truth_mismatch"] = int(sum(1 for s, n in zip(cells.src_genome_id, cells.n_truth)
                                  if int(T[C.gi[s]]) != int(n)))

# fold reference vs fold plan
fp = C.fold_plan
test_of_fold = {f["fold"]: set(f["test_genomes"]) for f in fp["folds"]}
test_groups_of_fold = {f["fold"]: set(f["test_groups"]) for f in fp["folds"]}
fold_of_group = {}
for f, gs in test_groups_of_fold.items():
    for g in gs:
        fold_of_group[g] = f
out["cell_fold_ne_fold_plan_fold_of_group"] = int(sum(
    1 for g, f in zip(cells.group_id, cells.fold) if fold_of_group.get(g) != f))
refchk = {}
for f in range(5):
    ref = set(vctx.jl(os.path.join(vctx.DATA, "prfba", f"f{f:02d}_d100",
                                   "_cache_ec_subset.json")).keys())
    allg = set(C.gids)
    fs = [x for x in fp["folds"] if x["fold"] == f][0]
    refchk[f] = {
        "ref_size": len(ref), "fold_plan_reference_genomes": len(fs["reference_genomes"]), "ref_equals_fold_plan_reference": ref == set(fs["reference_genomes"]),
        "ref_intersect_test_genomes": len(ref & test_of_fold[f]),
        "ref_equals_panel_minus_test": ref == (allg - test_of_fold[f]),
        "ref_outside_panel": len(ref - allg),
        "ref_genera_in_test_genera": len({C.panel[g]["genus"] for g in ref if C.panel[g]["genus"]} &
                                         {C.panel[g]["genus"] for g in test_of_fold[f]
                                          if C.panel[g]["genus"]}),
        "ref_groups_in_test_groups": len({C.ggroup[g] for g in ref} & test_groups_of_fold[f]),
        "cells_sources_in_ref": int(sum(1 for s in cells[cells.fold == f].src_genome_id
                                        if s in ref)),
    }
out["fold_reference"] = refchk
out["seconds"] = round(time.time() - t0, 1)
print(json.dumps(out, indent=1))
with open(os.path.join(vctx.OUT, "p0_sanity.json"), "w") as fh:
    json.dump(out, fh, indent=1)
