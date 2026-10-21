#!/usr/bin/env python
"""Summarise the template-file validation rebuilds (egc_rebuild.py --template-dir) of the
218 gapfilled models behind the C12 gate, one row per variant.

For each variant directory under data/egc_analysis/ it counts, against the original
gapfilled models (data/phenotype_gf/per_genome.jsonl for growth; data/egc_analysis/models
for the gate):
  gated_now_ok       originally leaking models that now pass the gate
  leaking, new_leaks models failing the gate now / of those, ones that passed originally
  no_growth          models with zero growth on the complete gapfilling medium
  lost_growth        models that grew originally and now do not, split by cause:
                     gapfill infeasible, or a gapfill solution found but its reversals
                     dropped (only possible when the variant does not integrate them)
  passing_*          over the 164 models that passed the gate and grew originally: the
                     median growth ratio to the original (zeros included, and over the
                     models that still grow), and how many moved by more than 10%

    python template_validation_summary.py            # writes data/template_direction_audit/
"""
import csv
import glob
import json
import os
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.abspath(os.path.join(HERE, ".."))
EGC = os.path.join(SR, "data/egc_analysis")
OUT = os.path.join(SR, "data/template_direction_audit")
THRESH = 1e-5

# variant -> the variant it is compared with for gapfill-solution changes
PARENT = {"PR1": "baseline", "PR1+rev": "baseline", "B": "PR1", "Bcur": "PR1", "Bcur+rev": "PR1+rev",
          "BC_dir": "B", "BC_gf": "B", "BC_full": "B", "BC_full+rev": "PR1+rev", "BC_core": "B",
          "PEN": "PR1", "PEN+rev": "PR1+rev", "PENcur+rev": "PR1+rev", "mainPEN": "baseline"}
# name -> (directory, what the templates are, reversals integrated?)
VARIANTS = [
    ("baseline", "models", "original v7.0 templates (main 2294299)", False),
    ("PR1", "models_tmplfix_none", "six-reaction fix (fcfe00d)", False),
    ("PR1+rev", "models_tmplfix_reversed_only", "six-reaction fix (fcfe00d)", True),
    ("B", "models_tv_B", "fcfe00d + all 38 tier-B flips (0fa1bd8)", False),
    ("Bcur", "models_tv_Bcur", "fcfe00d + 5 curated tier-B flips", False),
    ("Bcur+rev", "models_tv_Bcur_rev", "fcfe00d + 5 curated tier-B flips", True),
    ("BC_dir", "models_tv_BC_dir", "0fa1bd8 + tier C, direction field only", False),
    ("BC_gf", "models_tv_BC_gf", "0fa1bd8 + tier C, GapfillDirection field only", False),
    ("BC_full", "models_tv_BC_full", "0fa1bd8 + tier C (048a984): full database alignment", False),
    ("BC_full+rev", "models_tv_BC_full_rev", "0fa1bd8 + tier C (048a984): full database alignment", True),
    ("BC_core", "models_tv_BC_core", "0fa1bd8 + tier C in Core only", False),
    ("PEN", "models_tv_PEN", "fcfe00d + every negative reverse_penalty negated", False),
    ("PEN+rev", "models_tv_PEN_rev", "fcfe00d + every negative reverse_penalty negated", True),
    ("PENcur+rev", "models_tv_PENcur_rev", "fcfe00d + negative reverse_penalty negated where forward is favoured", True),
    ("mainPEN", "models_cf_mainPEN", "main + every negative reverse_penalty negated (16 genomes only)", False),
]


def load(d):
    out = {}
    for f in glob.glob(os.path.join(EGC, d, "*.stages.json")):
        s = json.load(open(f))
        out[s["key"]] = s["info"]
    return out


def main():
    orig = {}
    for line in open(os.path.join(SR, "data/phenotype_gf/per_genome.jsonl")):
        r = json.loads(line)
        orig[(r["src_genome_id"], r["arm"], r["pipeline"])] = r["gapfilled"].get("growth") or 0.0

    def g0(k):
        g, a, p = k.split("__")
        return orig[(g, a, "as_shipped" if p == "shared" else p)]

    base = load("models")
    assert len(base) == 218
    loaded = {}
    summary, per_model = {}, {k: {"model": k, "orig_gate_ok": base[k]["sanity_gate_ok"],
                                   "orig_growth": round(g0(k), 6)} for k in sorted(base)}
    for name, d, desc, rev in VARIANTS:
        cf = load(d)
        if not cf:
            print("skip %s: no outputs in %s" % (name, d))
            continue
        keys = [k for k in base if k in cf]
        loaded[name] = cf

        def g(k):
            v = cf[k].get("growth_complete")
            return g0(k) if v is None else v
        lost = [k for k in keys if g0(k) > THRESH and g(k) <= THRESH]
        infeasible = [k for k in lost if cf[k].get("infeasible")]
        dropped = [k for k in lost if not cf[k].get("infeasible") and cf[k].get("gapfill_reversed_not_integrated")]
        passing = [k for k in keys if base[k]["sanity_gate_ok"] and g0(k) > 0]
        ratio = [g(k) / g0(k) for k in passing]
        growing = [x for x in ratio if x > THRESH / 10]
        summary[name] = {
            "templates": desc, "reversals_integrated": rev, "dir": d, "n": len(keys),
            "gated_now_ok": sum(1 for k in keys if not base[k]["sanity_gate_ok"] and cf[k]["sanity_gate_ok"]),
            "leaking": sum(1 for k in keys if not cf[k]["sanity_gate_ok"]),
            "new_leaks": sum(1 for k in keys if base[k]["sanity_gate_ok"] and not cf[k]["sanity_gate_ok"]),
            "new_leak_models": sorted(k for k in keys if base[k]["sanity_gate_ok"] and not cf[k]["sanity_gate_ok"]),
            "no_growth": sum(1 for k in keys if g(k) <= THRESH),
            "lost_growth": len(lost), "lost_infeasible": len(infeasible),
            "lost_dropped_reversal": len(dropped),
            "lost_growth_genomes": sorted({k.split("__")[0] for k in lost}),
            "n_passing_growing": len(passing),
            "passing_growth_median_ratio": round(st.median(ratio), 3) if ratio else None,
            "passing_growth_median_ratio_still_growing": round(st.median(growing), 3) if growing else None,
            "passing_changed_gt10pct": sum(abs(x - 1) > 0.1 for x in ratio),
        }
        par = PARENT.get(name)
        if par in loaded:
            pc = loaded[par]
            both = [k for k in keys if k in pc and cf[k]["n_new_rxn"] and pc[k]["n_new_rxn"]]
            jac = [len(set(cf[k]["new_rxn"]) & set(pc[k]["new_rxn"])) /
                   len(set(cf[k]["new_rxn"]) | set(pc[k]["new_rxn"])) for k in both]

            def gp(k):
                v = pc[k].get("growth_complete")
                return g0(k) if v is None else v
            vs = [k for k in passing if k in pc and gp(k) > THRESH]
            summary[name]["vs_parent"] = {
                "parent": par, "n_both_gapfilled": len(both),
                "identical_gapfill_solution": sum(1 for x in jac if x == 1.0),
                "median_jaccard": round(st.median(jac), 3) if jac else None,
                "passing_changed_gt10pct_vs_parent": sum(abs(g(k) / gp(k) - 1) > 0.1 for k in vs),
                "n_passing_growing_in_parent": len(vs)}
        for k in keys:
            per_model[k]["%s_gate_ok" % name] = cf[k]["sanity_gate_ok"]
            per_model[k]["%s_growth" % name] = round(g(k), 6)
        s = summary[name]
        print("%-11s n=%3d gatedOK=%2d leak=%2d new=%d nogrow=%2d lost=%2d (infeas %d, dropped-rev %d) "
              "med=%.3f med_growing=%s >10%%=%d/%d" % (
                  name, s["n"], s["gated_now_ok"], s["leaking"], s["new_leaks"], s["no_growth"],
                  s["lost_growth"], s["lost_infeasible"], s["lost_dropped_reversal"],
                  s["passing_growth_median_ratio"], s["passing_growth_median_ratio_still_growing"],
                  s["passing_changed_gt10pct"], s["n_passing_growing"]))
    json.dump(summary, open(os.path.join(OUT, "validation_summary.json"), "w"), indent=1)
    cols = sorted({c for r in per_model.values() for c in r}, key=lambda c: (c != "model", c))
    with open(os.path.join(OUT, "validation_per_model.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(per_model.values())


if __name__ == "__main__":
    main()
