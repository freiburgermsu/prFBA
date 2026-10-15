#!/usr/bin/env python
"""Pure-function selftest for score_folds.py's PROVENANCE contracts.

Every check here is a round-4 review finding reduced to an assertion, so the
behaviour cannot regress silently (the end-to-end proofs live in SPEC 11.1's smoke
gate and in the review log; these are the parts that need no data on disk):

  1. the answer key is fingerprinted by CONTENT, never by path (hashing the paths
     made a copied tree refuse to reuse its own cells table);
  2. the scorer's code identity is a FIXED file list (sys.modules introspection
     differs between the pre-scoring reuse check and the post-scoring sidecar,
     because load_template_visible imports phenotype_arm lazily);
  3. cluster_bootstrap says WHY there is no CI, and the SESOI verdict is
     three-valued (no CI != below the SESOI);
  4. a control is HARD unless its gate text starts with "reported", by ONE rule;
  5. validity is stamped into every place a consumer reads it, from one reason list;
  6. a control this mode must emit, absent, is reported as missing;
  7. the phenotype arm is bound to the run, and C12 can actually FAIL.
"""
import json
import os
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score_folds as SF                                     # noqa: E402

FAIL = []


def check(name, got, want):
    ok = got == want
    print("  [%s] %-62s got %r" % ("PASS" if ok else "FAIL", name, got))
    if not ok:
        FAIL.append("%s: got %r want %r" % (name, got, want))


print("1. answer-key fingerprints are CONTENT, not paths")
a = types.SimpleNamespace(functions=os.path.join(HERE, "score_folds.py"),
                          maps=None, ec_universe=None, fanout=None,
                          annotation_ceiling=None, templates=None,
                          no_templates=True, panel="/nowhere/panel.json")
b = types.SimpleNamespace(**dict(vars(a)))
b.panel = "/elsewhere/panel.json"          # a different tree, same answer key file
check("same content -> same answer_key digest",
      SF.fp_digest(SF._answer_key_fingerprints(a)) ==
      SF.fp_digest(SF._answer_key_fingerprints(b)), True)
check("the digest is over {size, sha256} only",
      sorted(SF._answer_key_fingerprints(a)["functions"]), ["sha256", "size"])
check("every declared answer-key input has a tag",
      sorted(SF.ANSWER_KEY_INPUTS),
      ["annotation_ceiling", "ec_universe", "fanout", "functions", "maps",
       "templates"])
paths = SF.answer_key_paths(a)
check("--ec-universe is resolved beside --panel (a reader cannot guess it)",
      paths["ec_universe"], "/nowhere/panel_ec_universe.json")

print("2. scorer code identity is a fixed list, stable across calls")
fp1, files1 = SF.scorer_code_fp()
fp2, files2 = SF.scorer_code_fp()
check("stable", fp1 == fp2, True)
check("score_folds.py itself is in it", "score_folds.py" in files1, True)
check("the promoted audit modules are in it",
      all(n in files1 for n in ("design_audit_baselines.py", "design_audit_icc.py",
                                "design_audit_riskcoverage.py")), True)
check("paths are repo-relative (location independent)",
      files1["score_folds.py"]["path"].startswith("region_validation/"), True)
check("a module absent from disk is simply not listed",
      all(os.path.exists(os.path.join(SF.REPO, v["path"]))
          for v in files1.values()), True)

print("3. no CI is not 'below the SESOI'")
diffs = {("c%d" % i): (-0.3 if i % 2 else -0.2) for i in range(20)}
groups = {k: ("g%d" % (i % 5)) for i, k in enumerate(sorted(diffs))}
d, ci, p, ng, note = SF.cluster_bootstrap(diffs, groups, 500)
check("a real CI carries no note", (ci[0] is not None, note), (True, None))
check("n_groups", ng, 5)
d, ci, p, ng, note = SF.cluster_bootstrap(diffs, groups, 0)
check("--bootstrap 0 -> no CI, reason named", (ci, note),
      ([None, None], "bootstrap_disabled"))
one = {"c0": -0.3}
d, ci, p, ng, note = SF.cluster_bootstrap(one, {"c0": "g0"}, 500)
check("one holdout group -> no CI, reason named", (ci, note),
      ([None, None], "fewer_than_2_holdout_groups"))
d, ci, p, ng, note = SF.cluster_bootstrap(diffs, groups, 5)
check("a CI from 5 resamples is flagged underpowered", note,
      "bootstrap_below_%d" % SF.MIN_BOOTSTRAP)
d, ci, p, ng, note = SF.cluster_bootstrap({}, {}, 500)
check("no paired cells -> reason named", note, "no_paired_cells")
check("MIN_BOOTSTRAP is the smoke recipe's own 200 (SPEC 11.1)",
      SF.MIN_BOOTSTRAP, 200)

print("4. ONE definition of a HARD control")
fresh = {
    "C14": {"status": "fail", "gate": "for every launched run: ..."},
    "C7": {"status": "fail", "gate": "reported, not gated"},
    "C10": {"status": "fail", "gate": "  REPORTED, never gated (SPEC 7-C10)"},
    "C4": {"status": "pass", "gate": "set(fna ids) == ..."},
    "C5": {"status": "fail", "gate": "abs(pruned/full - 1) < 1e-6 ..."}}
check("hard failures", SF.hard_control_failures(fresh), ["C14", "C5"])
check("no fail -> nothing", SF.hard_control_failures(
    {"C4": {"status": "pass", "gate": "x"}}), [])

print("5. validity is stamped everywhere, from one list")
summ = {"meta": {"run_valid": True, "run_invalid_reasons": []},
        "headline": {"valid": True, "sentence": "leave-genus-out, ... SESOI 0.02 "
                                                "exceeded."},
        "arms": {"picrust2": {"ec": {"population": {"f1": {"value": 0.6955}}}},
                 "prfba_nestedcv": {"ec": {}, "rxn": {}}}}
valid, reasons = SF.stamp_invalid(summ, ["hard_control_failed:C4"])
check("valid", valid, False)
check("meta", summ["meta"]["run_valid"], False)
check("headline", summ["headline"]["valid"], False)
check("every arm x space block", all(
    blk.get("run_valid") is False for sp in summ["arms"].values()
    for blk in sp.values()), True)
check("the sentence is stamped", summ["headline"]["sentence"].startswith(
    "INVALID (hard_control_failed:C4)"), True)
before = summ["headline"]["sentence"]
SF.stamp_invalid(summ, ["hard_control_failed:C4"])
check("idempotent (no double prefix, no duplicate reason)",
      (summ["headline"]["sentence"] == before,
       summ["meta"]["run_invalid_reasons"]), (True, ["hard_control_failed:C4"]))
summ2 = {"meta": {}, "headline": {"sentence": "RUN_INCOMPLETE -- ..."}, "arms": {}}
SF.stamp_invalid(summ2, ["cells_scored_without_a_record"])
check("an already-stamped sentence is not re-prefixed",
      summ2["headline"]["sentence"].startswith("RUN_INCOMPLETE"), True)
clean = {"meta": {"run_invalid_reasons": []}, "headline": {"sentence": "ok."},
         "arms": {}}
check("no reasons -> valid, sentence untouched",
      (SF.stamp_invalid(clean, [])[0], clean["headline"]["sentence"]),
      (True, "ok."))

print("6. an absent control is not a passing one")
prov = {"mode": "fold", "inventory_fp": "i", "answer_key_fp": "a",
        "scorer_code_fp": "s", "evaluated_at": "now"}
merged = {"C1": {"status": "pass", "provenance": dict(prov, mode="recipe_control")},
          "C2": {"status": "pass"},
          "C3": {"status": "pass", "provenance": dict(prov, inventory_fp="OLD")},
          "C14": {"status": "pass", "provenance": dict(prov)}}
rep = SF.controls_provenance_report(merged, {"C14"}, prov, "fold")
check("missing_expected names what fold must emit", rep["missing_expected"],
      ["C10", "C3", "C4", "C6", "C7"])
check("a foreign-mode entry is reported as such",
      [r["id"] for r in rep["from_other_modes"]], ["C1"])
check("an entry with no provenance is named", rep["without_provenance"], ["C2"])
check("an entry this SAME mode recorded but did not re-evaluate is stale",
      [r["id"] for r in rep["stale_same_mode"]], ["C3"])
check("a cross-mode entry is classified by MODE, not by inventory_fp (which can "
      "never match across modes)",
      [(r["id"], r["scorer_code_matches"], r["answer_key_matches"])
       for r in rep["from_other_modes"]], [("C1", True, True)])
rep2 = SF.controls_provenance_report(
    {"C1": {"status": "pass", "provenance": dict(prov, mode="recipe_control",
                                                 scorer_code_fp="OTHER")}},
    {"C14"}, prov, "fold")
check("a cross-mode entry from a DIFFERENT scorer is flagged",
      [(r["id"], r["scorer_code_matches"]) for r in rep2["from_other_modes"]],
      [("C1", False)])
check("C14 is in every gated mode's must-emit set",
      all("C14" in SF.EXPECTED_CONTROLS[m] for m in ("fold", "density",
                                                     "include_self")), True)
check("the data-conditional controls are reported, never required",
      sorted(set(SF.CONDITIONAL_CONTROLS) &
             set(SF.EXPECTED_CONTROLS["fold"])), [])

print("7. the phenotype arm is bound, and C12 can fail")
pb_unbound = {"path": "/t/phenotype_summary.json", "present": True, "bound": False,
              "defects": [], "declared": {
                  "identity_control": {"pass": True, "jaccard": 1.0},
                  "n_gated_out": {}, "n_models": 91}}
c = SF.phenotype_controls(pb_unbound)
check("an unbound legacy summary keeps its own verdicts",
      (c["C11"]["status"], c["C12"]["status"]), ("pass", "pass"))
c = SF.phenotype_controls(dict(pb_unbound, defects=[
    {"defect": "phenotype_summary_unreadable"}]))
check("a binding defect FAILS both, never makes them vanish",
      (c["C11"]["status"], c["C12"]["status"]), ("fail", "fail"))
c = SF.phenotype_controls(dict(pb_unbound, declared=dict(
    pb_unbound["declared"], n_gated_out={"truth": 1})))
check("one gated TRUTH model fails C12 (it shrinks every arm's denominator)",
      c["C12"]["status"], "fail")
c = SF.phenotype_controls(dict(pb_unbound, declared=dict(
    pb_unbound["declared"], n_gated_out={"prfba_exact_tie": 20})))
check("gating 22%% of 91 models fails C12", c["C12"]["status"], "fail")
c = SF.phenotype_controls(dict(pb_unbound, declared=dict(
    pb_unbound["declared"], n_gated_out={"prfba_exact_tie": 5})))
check("gating 5 of 91 passes (the gate working as designed)",
      c["C12"]["status"], "pass")
c = SF.phenotype_controls(dict(pb_unbound, declared=dict(
    pb_unbound["declared"], identity_control={"pass": False, "jaccard": 0.9})))
check("a failing identity control fails C11", c["C11"]["status"], "fail")
check("C12's threshold is recorded in the artefact",
      SF.phenotype_controls(pb_unbound)["C12"]["measured"]["max_gated_fraction"],
      SF.C12_MAX_GATED_FRACTION)
check("no phenotype summary -> no C11/C12 at all (the arm did not run)",
      SF.phenotype_controls({"present": False}), {})

print()
if FAIL:
    print("== %d PROVENANCE CHECK(S) FAILED ==" % len(FAIL))
    for f in FAIL:
        print("  " + f)
    sys.exit(1)
print("== ALL PROVENANCE CHECKS PASSED ==")
