#!/usr/bin/env python
"""Self-test for score_folds.py's run-INVENTORY gate (C14 legs 6-7).

Pure function, no fixtures, seconds to run: it builds the FULL-SCALE shape SPEC 11.2
describes (10 folds x 2 domains x 3 densities) in memory and asserts that
``run_inventory`` DERIVES exactly 60 PICRUSt2 forced + 10 nsti + 30 prFBA fold + 10
prFBA forced-answer = 110 run keys from the fold plan, the frozen config and the
declared scope -- and that each of the campaign failure modes the reviewers
demonstrated (a killed (fold, domain) job whose refs/ survives, a never-launched
archaeal arm, a lost density tier, a leftover retry directory, a fold launched with a
shorter --select-modes, a dead nsti run) comes back ok=False with the offending key
named.  Run it before the launch; it needs no data on disk.

    ~/Documents/py_venv/bin/python scripts/_selftest_score_folds_inventory.py
"""
import os, sys, types
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score_folds as SF

FOLDS = range(10)
DENS = [0.25, 0.5, 1.0]
plan = {"meta": {"frozen": {"density_fractions": DENS}},
        "folds": [{"fold": f, "test_cells": ["g%d_%d__V4" % (f, i) for i in range(4)]}
                  for f in FOLDS]}
panel = {}
truth = {}
for f in FOLDS:                       # two bacterial + two archaeal cells per fold
    for i in range(4):
        truth["g%d_%d__V4" % (f, i)] = {"domain": "Bacteria" if i < 2 else "Archaea"}
manifests = {}
for f in FOLDS:
    for dom in ("bac", "arc"):
        for d in DENS:
            manifests["f%02d_%s_d%03d" % (f, dom, round(d * 100))] = {
                "mode": "fold", "fold": f, "domain": dom, "density": d,
                "reference_genomes": ["x"]}

def p2(f, dom, d, routing="forced", name=None):
    return {"dir": "/runs/picrust2/" + (name or "f%02d_%s_d%03d"
                                        % (f, dom, round(d * 100))),
            "name": name or "f%02d_%s_d%03d" % (f, dom, round(d * 100)),
            "meta": {"mode": "fold", "fold": f, "domain": dom, "density": d,
                     "routing": routing},
            "pred": {}, "fold": f, "domain": dom, "density": d, "routing": routing,
            "ref_names": ["f%02d_%s_d%03d" % (f, dom, round(d * 100))]}

def pf(f, d, forced=False, name=None):
    nm = name or ("f%02d_d%03d%s" % (f, round(d * 100), "_forced" if forced else ""))
    return {"dir": "/runs/prfba/" + nm, "name": nm,
            "meta": {"mode": "fold", "fold": f, "density": d,
                     "forced_answer": forced,
                     "select_modes": ["legacy", "exact_tie", "graduated",
                                      "exact_tie_guards"]},
            "modes": {}, "fold": f, "density": d, "forced": forced,
            "ref_names": ["f%02d_bac_d%03d" % (f, round(d * 100))]}

p2_runs = [p2(f, dom, d) for f in FOLDS for dom in ("bac", "arc") for d in DENS]
p2_runs += [p2(f, "both", 1.0, "nsti", name="f%02d_nsti" % f) for f in FOLDS]
pf_runs = [pf(f, d) for f in FOLDS for d in DENS]
pf_runs += [pf(f, 1.0, True) for f in FOLDS]

args = types.SimpleNamespace(
    mode="fold", picrust2="/runs/picrust2", prfba="/runs/prfba", refs="/runs/refs",
    fold_plan="/runs/fold_plan.json", panel="/runs/panel.json", truth="/runs/t.json",
    neighbors="/runs/n.json", rank_blocked="genus", space="ec,rxn", limit_cells=0,
    primary_select_mode="exact_tie", no_templates=False, clade_max_cap=3000,
    expect_folds=None, expect_densities=None, expect_domains=None,
    expect_nsti="primary", expect_forced_answer="primary")

inv = SF.run_inventory(args, plan, panel, truth, manifests, p2_runs, pf_runs,
                       plan["meta"]["frozen"])
exp = inv["expected"]
n_forced = sum(1 for k in exp if k.startswith("picrust2") and k.endswith("forced"))
n_nsti = sum(1 for k in exp if k.endswith("nsti"))
n_pf = sum(1 for k in exp if k.startswith("prfba") and k.endswith("|fold"))
n_pfa = sum(1 for k in exp if k.endswith("forced_answer"))
print("expected: picrust2 forced %d, nsti %d, prfba fold %d, prfba forced-answer %d, "
      "total %d" % (n_forced, n_nsti, n_pf, n_pfa, inv["n_expected"]))
print("discovered %d keys in %d directories; ok=%s reasons=%s missing=%s "
      "unexpected=%s duplicates=%s"
      % (inv["n_discovered"], inv["n_run_directories"], inv["ok"], inv["reasons"],
         inv["missing"], inv["unexpected"], inv["duplicates"]))
assert (n_forced, n_nsti, n_pf, n_pfa) == (60, 10, 30, 10), (n_forced, n_nsti, n_pf,
                                                             n_pfa)
assert inv["ok"] is True and inv["n_expected"] == 110 == inv["n_discovered"]
print("SPEC 11.2 inventory 60+10+30+10 == 110 DERIVED and matched: OK")

# now break it the five ways the reviewers broke the campaign
def check(label, p2r, pfr, want_reasons):
    i = SF.run_inventory(args, plan, panel, truth, manifests, p2r, pfr,
                         plan["meta"]["frozen"])
    print("%-34s ok=%-5s reasons=%-60s detail=%s"
          % (label, i["ok"], i["reasons"],
             (i["missing"][:2] or [d["key"] for d in i["duplicates"]][:2] or
              [d["key"] for d in i["unexpected"]][:2] or
              [d["defect"] for d in i["scope_mismatch"]][:2])))
    assert i["ok"] is False and set(i["reasons"]) == set(want_reasons), i["reasons"]

check("one (fold,domain) job killed", [r for r in p2_runs if r["name"] !=
                                       "f07_arc_d100"], pf_runs,
      ["run_inventory_incomplete"])
check("archaeal arm never launched", [r for r in p2_runs if "arc" not in r["name"]],
      pf_runs, ["run_inventory_incomplete"])
check("whole d050 tier lost", [r for r in p2_runs if "d050" not in r["name"]],
      [r for r in pf_runs if "d050" not in r["name"]], ["run_inventory_incomplete"])
check("stale retry directory", p2_runs + [p2(0, "bac", 1.0,
                                             name="f00_bac_d100_retry")], pf_runs,
      ["duplicate_run_directories"])
short = [r for r in pf_runs if r["name"] != "f00_d100_forced"]
r0 = pf(0, 1.0, True)
r0["meta"]["select_modes"] = ["exact_tie"]
check("fold 0 short --select-modes", p2_runs, short + [r0], ["run_scope_mismatch"])
nsti_gone = [r for r in p2_runs if r["name"] != "f03_nsti"]
check("one nsti run dead", nsti_gone, pf_runs, ["run_inventory_incomplete"])
print("all five full-scale failure modes caught")

# the one residual case: the d050 tier's REFERENCES died too, so fold mode cannot
# expect runs for it -- reported in fold mode, GATED in --mode density
man2 = {k: v for k, v in manifests.items() if "d050" not in k}
p2b = [r for r in p2_runs if "d050" not in r["name"]]
pfb = [r for r in pf_runs if "d050" not in r["name"]]
i = SF.run_inventory(args, plan, panel, truth, man2, p2b, pfb,
                     plan["meta"]["frozen"])
print("refs+runs of d050 both gone, --mode fold: ok=%s expected_dens=%s "
      "frozen_not_built=%s" % (i["ok"], i["derivation"]["densities_expected"],
                               i["derivation"]["densities_in_frozen_plan_not_built"]))
assert i["ok"] is True and i["derivation"]["densities_in_frozen_plan_not_built"] == ["0.5"]
args.mode = "density"
i = SF.run_inventory(args, plan, panel, truth, man2, p2b, pfb,
                     plan["meta"]["frozen"])
print("                           --mode density: ok=%s reasons=%s missing=%d "
      "(e.g. %s)" % (i["ok"], i["reasons"], len(i["missing"]), i["missing"][:1]))
assert i["ok"] is False and len(i["missing"]) == 30
args.mode = "fold"
args.expect_densities = "0.25,0.5,1.0"
i = SF.run_inventory(args, plan, panel, truth, man2, p2b, pfb,
                     plan["meta"]["frozen"])
print("            --mode fold --expect-densities 0.25,0.5,1.0: ok=%s missing=%d"
      % (i["ok"], len(i["missing"])))
assert i["ok"] is False and len(i["missing"]) == 30
print("density axis comes from the PLAN, never from the launched runs")

print("== ALL INVENTORY CHECKS PASSED ==")
