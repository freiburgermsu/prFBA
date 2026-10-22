#!/usr/bin/env python
"""Summarise the template-file validation rebuilds (egc_rebuild.py --template-dir) of the
218 gapfilled models behind the C12 gate, one entry per variant.

Every variant is compared with a reference built the same way, so that a template effect
is not confused with the integration of the gapfiller's reversals:
  - variants without reversal integration ('none') against the original runs
    (data/egc_analysis/models; growth from data/phenotype_gf/per_genome.jsonl);
  - variants with reversals integrated ('reversed_only') against the original templates
    rebuilt with reversals integrated (base+rev).
The gate verdict always refers to the original runs. The reversals only loosen bounds of
the final model, and the baseline rebuild with reversals has the same 44 leaks.

Per variant:
  gated_now_ok        originally leaking models that now pass the gate
  leaking, new_leaks  models failing the gate now / of those, ones that passed originally
  no_growth           models with zero growth on the complete gapfilling medium
  lost_growth         models that grew in the reference and now do not, split by cause:
                      gapfill infeasible, or a gapfill solution found but its reversals
                      dropped
  passing_*           over models that passed the gate originally and grow in the
                      reference: the median growth ratio to the reference (zeros included,
                      and over the models that still grow), and how many moved by > 10%
  vs_parent           gapfill-solution identity and growth changes against the variant
                      it is stacked on

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

# name -> (directory, templates, reversals integrated?, parent for gapfill-solution comparison)
# ModelSEEDTemplates commits: main 2294299; five-reaction fix cf63a8f; its curated
# follow-ups 0805094 (five opposite-direction reactions) and 8064b97 (penalty sign);
# blanket penalty da0c441; reference alignment c6b14d6 (41 flips) and ea88229 (complete).
# The earlier six-reaction variant 7d104ac set NadX (rxn05119/rxn05117) '<' instead of
# rxn08167 '>'; its rows are kept for the record.
VARIANTS = [
    ("baseline", "models", "original v7.0 (main 2294299), original runs", False, None),
    ("base+rev", "models_tv_base_rev", "original v7.0 (main 2294299)", True, None),
    ("PR1v2", "models_tv_PR1v2", "five-reaction fix (cf63a8f)", False, "baseline"),
    ("PR1v2+rev", "models_tv_PR1v2_rev", "five-reaction fix (cf63a8f)", True, "base+rev"),
    ("Bcur2+rev", "models_tv_Bcur2_rev", "cf63a8f + 5 opposite-direction reactions (0805094)", True, "PR1v2+rev"),
    ("PENcur2+rev", "models_tv_PENcur2_rev", "cf63a8f + penalty sign where forward is favoured (8064b97)", True, "PR1v2+rev"),
    ("PEN2+rev", "models_tv_PEN2_rev", "cf63a8f + every negative reverse_penalty negated (da0c441)", True, "PR1v2+rev"),
    ("BC_full+rev", "models_tv_BC_full_rev", "complete database alignment (ea88229; same files as on 7d104ac)", True, "PR1v2+rev"),
    ("BC_full", "models_tv_BC_full", "complete database alignment (ea88229; same files as on 7d104ac)", False, "baseline"),
    ("mainPEN", "models_cf_mainPEN", "main + every negative reverse_penalty negated (142 models of 16 genomes)", False, "baseline"),
    # earlier six-reaction variant (7d104ac, NadX '<') and variants stacked on it
    ("PR1", "models_tmplfix_none", "six-reaction variant (7d104ac)", False, "baseline"),
    ("PR1+rev", "models_tmplfix_reversed_only", "six-reaction variant (7d104ac)", True, "base+rev"),
    ("B", "models_tv_B", "7d104ac + all 38 opposite-direction flips", False, "PR1"),
    ("Bcur", "models_tv_Bcur", "7d104ac + 5 opposite-direction reactions", False, "PR1"),
    ("Bcur+rev", "models_tv_Bcur_rev", "7d104ac + 5 opposite-direction reactions", True, "PR1+rev"),
    ("BC_dir", "models_tv_BC_dir", "7d104ac + 38 + tier C, direction field only", False, "B"),
    ("BC_gf", "models_tv_BC_gf", "7d104ac + 38 + tier C, GapfillDirection field only", False, "B"),
    ("BC_core", "models_tv_BC_core", "7d104ac + 38 + tier C in Core only", False, "B"),
    ("PEN", "models_tv_PEN", "7d104ac + every negative reverse_penalty negated", False, "PR1"),
    ("PEN+rev", "models_tv_PEN_rev", "7d104ac + every negative reverse_penalty negated", True, "PR1+rev"),
    ("PENcur+rev", "models_tv_PENcur_rev", "7d104ac + penalty sign where forward is favoured", True, "PR1+rev"),
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

    def g_orig(k):
        g, a, p = k.split("__")
        return orig[(g, a, "as_shipped" if p == "shared" else p)]

    base = load("models")
    assert len(base) == 218
    loaded = {}

    def growth(name, k):
        if name == "baseline":
            return g_orig(k)
        v = loaded[name][k].get("growth_complete")
        return g_orig(k) if v is None else v

    summary, per_model = {}, {k: {"model": k, "orig_gate_ok": base[k]["sanity_gate_ok"],
                                   "orig_growth": round(g_orig(k), 6)} for k in sorted(base)}
    loaded["baseline"] = base
    for name, d, desc, rev, parent in VARIANTS:
        cf = load(d)
        if not cf:
            print("skip %s: no outputs in %s" % (name, d))
            continue
        loaded[name] = cf
        if name == "baseline":
            keys = sorted(base)
        else:
            keys = [k for k in base if k in cf]
        ref = "base+rev" if rev and name != "base+rev" else "baseline"
        if ref not in loaded:
            print("skip %s: reference %s not loaded" % (name, ref))
            continue
        keys = [k for k in keys if k in loaded[ref]]
        gr = {k: growth(ref, k) for k in keys}
        g = {k: growth(name, k) for k in keys}
        lost = [k for k in keys if gr[k] > THRESH and g[k] <= THRESH]
        infeasible = [k for k in lost if cf[k].get("infeasible")]
        dropped = [k for k in lost if not cf[k].get("infeasible") and cf[k].get("gapfill_reversed_not_integrated")]
        passing = [k for k in keys if base[k]["sanity_gate_ok"] and gr[k] > THRESH]
        ratio = [g[k] / gr[k] for k in passing]
        growing = [x for x in ratio if x > 0]
        s = summary[name] = {
            "templates": desc, "reversals_integrated": rev, "dir": d, "growth_reference": ref,
            "n": len(keys),
            "gated_now_ok": sum(1 for k in keys if not base[k]["sanity_gate_ok"] and cf[k]["sanity_gate_ok"]),
            "gated_now_ok_via_infeasible_gapfill": sum(
                1 for k in keys if not base[k]["sanity_gate_ok"] and cf[k]["sanity_gate_ok"] and cf[k].get("infeasible")),
            "leaking": sum(1 for k in keys if not cf[k]["sanity_gate_ok"]),
            "new_leaks": sum(1 for k in keys if base[k]["sanity_gate_ok"] and not cf[k]["sanity_gate_ok"]),
            "new_leak_models": sorted(k for k in keys if base[k]["sanity_gate_ok"] and not cf[k]["sanity_gate_ok"]),
            "no_growth": sum(1 for k in keys if g[k] <= THRESH),
            "lost_growth": len(lost), "lost_infeasible": len(infeasible),
            "lost_dropped_reversal": len(dropped),
            "lost_growth_genomes": sorted({k.split("__")[0] for k in lost}),
            "n_passing_growing": len(passing),
            "passing_growth_median_ratio": round(st.median(ratio), 3) if ratio else None,
            "passing_growth_median_ratio_still_growing": round(st.median(growing), 3) if growing else None,
            "passing_changed_gt10pct": sum(abs(x - 1) > 0.1 for x in ratio),
        }
        if parent in loaded:
            pc = loaded[parent]
            both = [k for k in keys if k in pc and cf[k].get("n_new_rxn") and pc[k].get("n_new_rxn")]
            jac = [len(set(cf[k]["new_rxn"]) & set(pc[k]["new_rxn"])) /
                   len(set(cf[k]["new_rxn"]) | set(pc[k]["new_rxn"])) for k in both]
            vs = [k for k in passing if k in pc and growth(parent, k) > THRESH]
            s["vs_parent"] = {
                "parent": parent, "n_both_gapfilled": len(both),
                "identical_gapfill_solution": sum(1 for x in jac if x == 1.0),
                "median_jaccard": round(st.median(jac), 3) if jac else None,
                "passing_changed_gt10pct_vs_parent": sum(abs(g[k] / growth(parent, k) - 1) > 0.1 for k in vs),
                "n_passing_growing_in_parent": len(vs)}
        if name != "baseline":
            for k in keys:
                per_model[k]["%s_gate_ok" % name] = cf[k]["sanity_gate_ok"]
                per_model[k]["%s_growth" % name] = round(g[k], 6)
        print("%-12s n=%3d gatedOK=%2d leak=%2d new=%d nogrow=%2d lost=%2d (infeas %d, dropped-rev %d) "
              "med=%s med_growing=%s >10%%=%d/%d" % (
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
