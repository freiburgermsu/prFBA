"""Synthesize minimal component-1 outputs so build_reference.py can be exercised.

plan_folds.py (component 1) is being written in parallel, so its outputs do not
exist yet.  This builds the smallest honest stand-ins in the §4.2 formats from
already-verified repo data (kfold_audit/{panel_16s.fna,panel_16s_meta.json,
panel_ec_sets.json} + selection_10k.json).
"""
import json, os, sys, collections
sys.path.insert(0, "/home/freiburger/Documents/prFBA/region_validation/shared_reference/scripts")
from build_reference import FROZEN

RV = "/home/freiburger/Documents/prFBA/region_validation"
AUD = os.path.join(RV, "picrust2_comparison/data/kfold_audit")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RV, "shared_reference/data_selftest")
os.makedirs(OUT, exist_ok=True)

sel = {r["genome_id"]: r for r in json.load(open(os.path.join(RV, "data/selection_10k.json")))["genomes"]}
ecs = json.load(open(os.path.join(AUD, "panel_ec_sets.json")))
seqs = {}
cur = None
for l in open(os.path.join(AUD, "panel_16s.fna")):
    if l[0] == ">":
        cur = l[1:].strip(); seqs[cur] = []
    else:
        seqs[cur].append(l.strip())
seqs = {k: "".join(v).upper() for k, v in seqs.items()}

def eligible(g):
    r = sel.get(g)
    if not r or g not in seqs or g not in ecs or not ecs[g]["ec"]: return False
    if not (1200 <= len(seqs[g]) <= 1800): return False
    return set(seqs[g]) <= set("ACGTN")

bac = sorted(g for g in seqs if eligible(g) and sel[g]["domain"] == "Bacteria")
arc = sorted(g for g in seqs if eligible(g) and sel[g]["domain"] == "Archaea")
# Deterministic, phylum-spread pick: round-robin over phyla.
def spread(ids, n):
    byp = collections.OrderedDict()
    for g in ids: byp.setdefault(sel[g]["phylum"], []).append(g)
    out, keys = [], list(byp)
    i = 0
    while len(out) < n and keys:
        k = keys[i % len(keys)]
        if byp[k]: out.append(byp[k].pop(0))
        else: keys.remove(k); i -= 1
        i += 1
    return sorted(out)

# Keep a couple of same-genus pairs so monophyly has something to test.
pick_b = spread(bac, 30); pick_a = spread(arc, 8)
genera = collections.Counter(sel[g]["genus"] for g in pick_b)
for g in bac:
    if len(pick_b) >= 34: break
    if g not in pick_b and genera.get(sel[g]["genus"], 0) >= 1:
        pick_b.append(g); genera[sel[g]["genus"]] += 1
pick_b = sorted(set(pick_b))
chosen = pick_b + pick_a
print("bac", len(pick_b), "arc", len(pick_a), "genera with >=2:",
      sum(1 for k, v in collections.Counter(sel[g]["genus"] for g in chosen).items() if v >= 2))

GRAMPOS = {"Bacillota", "Firmicutes", "Actinomycetota", "Actinobacteria",
           "Chloroflexota", "Deinococcota", "Thermotogota", "Mycoplasmatota", "Tenericutes"}
def template(r):
    if r["domain"] == "Archaea": return "archaea"
    return "grampos" if (r["phylum"] or "").split("_")[0] in GRAMPOS else "gramneg"

genomes = {}
for g in chosen:
    r = sel[g]
    genomes[g] = {"label": FROZEN["label_prefix"] + g, "taxon_id": r["taxon_id"],
                  "md5": r["best_16s_md5"], "len": len(seqs[g]), "domain": r["domain"],
                  "phylum": r["phylum"], "class": r["class"], "order": r["order"],
                  "family": r["family"], "genus": r["genus"], "species": r["species"],
                  "n_full_length_16s": r.get("n_full_length_16s") or 1,
                  "n_ec": len(set(ecs[g]["ec"])), "occupancy": 0.99, "template": template(r)}
json.dump({"meta": {"frozen": FROZEN, "generated": "selftest", "git_commit": "selftest",
                    "n_panel": len(sel), "n_eligible": len(genomes),
                    "script": "sr_selftest_setup.py",
                    "drop_counts": {"no_function_set": 0, "empty_ec": 0, "length": 0,
                                    "occupancy": 0, "bad_chars": 0}},
           "genomes": genomes, "dropped": {}},
          open(os.path.join(OUT, "panel_eligible.json"), "w"), indent=1)

for dom, ids in (("bac", pick_b), ("arc", pick_a)):
    with open(os.path.join(OUT, f"panel_{dom}.fna"), "w") as fh:
        for g in ids: fh.write(f">{FROZEN['label_prefix']}{g}\n{seqs[g]}\n")

prev = collections.Counter()
for g in chosen: prev.update({e for e in ecs[g]["ec"] if "-" not in e})
order = sorted(prev)
json.dump({"meta": {"n_genomes": len(chosen), "n_ec": len(order)}, "ec": order,
           "prevalence": {e: prev[e] / len(chosen) for e in order},
           "ec_index": {e: i for i, e in enumerate(order)}},
          open(os.path.join(OUT, "panel_ec_universe.json"), "w"))
json.dump({g: {"ec": sorted({e for e in ecs[g]["ec"] if "-" not in e})} for g in chosen},
          open(os.path.join(OUT, "panel_ec_sets.json"), "w"))

# k=2 fold plan: alternate genomes between folds within each domain.
folds = []
for i in range(2):
    test = sorted([g for j, g in enumerate(pick_b) if j % 2 == i] +
                  [g for j, g in enumerate(pick_a) if j % 2 == i])
    ref = sorted(set(chosen) - set(test))
    folds.append({"fold": i, "test_groups": [f"grp{i}"], "test_genomes": test,
                  "reference_genomes": ref, "test_cells": [], "test_md5s": [],
                  "n_test_genomes": len(test), "n_reference_genomes": len(ref),
                  "n_reference_bac": sum(1 for g in ref if sel[g]["domain"] == "Bacteria"),
                  "n_reference_arc": sum(1 for g in ref if sel[g]["domain"] == "Archaea"),
                  "domain_frac": {}, "top5_phylum_frac": {}})
json.dump({"meta": {"frozen": FROZEN, "rank_blocked": "genus", "k": 2,
                    "n_eligible": len(chosen), "n_components": 2, "n_singletons": 0,
                    "largest_component": 1, "generated": "selftest", "git_commit": "selftest"},
           "groups": {"grp0": folds[0]["test_genomes"], "grp1": folds[1]["test_genomes"]},
           "genome_group": {g: ("grp0" if g in folds[0]["test_genomes"] else "grp1") for g in chosen},
           "folds": folds, "leakage_preflight": {}},
          open(os.path.join(OUT, "fold_plan.json"), "w"), indent=1)
print("WROTE", OUT, "ecu", len(order), "fold0 refs", folds[0]["n_reference_genomes"])
