#!/usr/bin/env python
"""Builds prFBA/manuscript_numbers/claims_spec.json (W7 extractor spec) from the plan's
claims_by_stage.json plus the per-item extractors written below.

Run: ~/Documents/py_venv/bin/python $W/tools/build_claims_spec.py
Each item id is "<inventory>:<key>". Values: see manuscript_numbers/collect.py docstring.
"""
import json
import os

W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
R = "/home/freiburger/Documents/prFBA"
CLAIMS = os.path.join(W, "expected/plan/claims_by_stage.json")
OUT = os.path.join(R, "manuscript_numbers/claims_spec.json")

R8 = ["FullLength16S", "V1-V2", "V1-V3", "V3-V4", "V4", "V4-V5", "V4-V5_944R", "V7-V9"]
V3ON = ["V3-V4", "V4", "V4-V5", "V4-V5_944R", "V7-V9"]
V1A = ["FullLength16S", "V1-V2", "V1-V3"]

# ---------------------------------------------------------------- artefacts
CS = "region_validation_d100/picrust2_comparison/data/comparison_summary.json"
CS20 = "region_validation/picrust2_comparison/data/comparison_summary.json"
CC = "region_validation_d100/picrust2_comparison/data/comparison_cells.csv.gz"
DIAL = "region_validation_d100/picrust2_comparison/data/annotation_dialect.json"
CEIL = "region_validation/picrust2_comparison/data/annotation_ceiling.json"
PROBE = "region_validation/picrust2_comparison/data/depth_probe/depth_probe_results.json"
GC = "region_validation_d100/data/gene_capture_by_region.json"
GCL = "region_validation_d100_legacy/data/gene_capture_by_region.json"
GC20 = "region_validation/data/gene_capture_by_region.json"
GCC = "region_validation_d100/data/gene_capture_cells.csv"
GCCL = "region_validation_d100_legacy/data/gene_capture_cells.csv"
TAU = "region_validation_d100/data/gene_capture_tau_compare.json"
SRD = "region_validation_d100_legacy/data/self_recovery_diagnosis.json"
E3 = "region_validation_d100/data/pfunc_e3_gate.json"
RR = "region_validation_d100/data/pfunc_e3_gate/review_recompute.json"
SR = "region_validation/shared_reference/"
SG = SR + "data/summary_genus.json"
SF = SR + "data_family/summary_family.json"
CTRL = SR + "data/controls.json"
PATP = SR + "data/phenotype_atp_c13/phenotype_summary.json"
PGF = SR + "data/phenotype_gf_c13/phenotype_summary.json"
PGF0 = SR + "data/phenotype_gf/phenotype_summary.json"
WOBJ = SR + "data/weighted_objective/summary.json"
OT = SR + "data/objective_tuning/summary.json"
OTM = SR + "data/objective_tuning/meta.json"
WEQ = SR + "data/objective_tuning/wrapper_equivalence.json"
CAL = SR + "data/calibration/scores.json"
PFE = SR + "data/pfunc_estimators/scores.json"
TP = "three_panel/three_panel_comparison_d100.json"
PA = "three_panel/three_panel_mag_panelarm.json"
EM = "three_panel/emptied_strata_d100.json"
GT = "three_panel/three_panel_mag_gtdb.json"
PC = "mag_validation/out/paired_comparison_d100.json"
CX = "mag_validation/out/crossover_d100.json"
MSF = "mag_validation/out/mag_summary_full.json"
MST = "mag_validation/out/mag_summary_tau99.json"
N500 = "mag_validation/out/summary_n500.json"
AFG = "mag_validation/out/af_gate_counterfactual.json"
MAGC = "mag_validation/out/mag_cells_all.csv"
MPC = "two_panel/mag_panel_comparison.json"
DEP = "deployment_d100/stats_d100_emilykin.json"
DEP_NEW = "deployment_d100/stats_d100_emilykin_bvbrc.json"     # written by plan step R1.K/K1
BRS = "deployment_d100/bvbrc_rebuild_summary.json"
EDL = "edlib_vs_gpu_assessment/assessment_summary.json"


# ---------------------------------------------------------------- helpers
def V(label, art, sel, printed, fmt="f3", **kw):
    d = {"label": label, "artifact": art, "selector": sel, "printed": printed, "fmt": fmt}
    d.update(kw)
    return d


def RANGE(label, art, pre, keys, post, lo, hi, fmt="f3", lo_red="min", hi_red="max", **kw):
    """Two values (lo/hi) over pre + {each: keys} + post."""
    sel = pre + [{"each": keys}] + post
    return [V(f"{label} (low)", art, sel, lo, fmt, reduce=lo_red, **kw),
            V(f"{label} (high)", art, sel, hi, fmt, reduce=hi_red, **kw)]


def X(label, expr, vars_, printed, fmt="f3", **kw):
    d = {"label": label, "expr": expr, "vars": vars_, "printed": printed, "fmt": fmt}
    d.update(kw)
    return d


def S(art, sel, **kw):
    d = {"artifact": art, "selector": sel}
    d.update(kw)
    return d


def CSV(label, art, recipe, args, printed, fmt="f3", **kw):
    d = {"label": label, "artifact": art, "recipe": recipe, "args": args, "printed": printed, "fmt": fmt}
    d.update(kw)
    return d


def D(label, key, printed, fmt="f3", **kw):
    """A W6 derived number (manuscript_numbers/derived_numbers.{pre,post}.json)."""
    d = {"label": label, "artifact": f"derived:{key}", "printed": printed, "fmt": fmt}
    d.update(kw)
    return d


# comparison_summary shortcuts
def cs(*path):
    return list(path)


ALL = ["by_region", "ec", "all", "all"]


def csm(method, *rest, art=CS):
    return ALL + [method] + list(rest)


def nsti(b, m, *rest):
    return ["by_stratum", "ec", "nsti_bin", b, m] + list(rest)


def ident(b, m, *rest):
    return ["by_stratum", "ec", "ident_bin", b, m] + list(rest)


def strat(flag, m, *rest):
    return ["by_stratum", "ec", "src_in_picrust2_ref", flag, m] + list(rest)


def pt(region):
    return ["paired_tests_ec_f1", f"{region}|prfba_exclude_vs_picrust2|conditional"]


# three-panel shortcuts
def tp(region, regime, endpoint, *rest):
    return ["cells", {"where": {"region": region, "regime": regime, "endpoint": endpoint}}] + list(rest)


def pa(regime, comp, *rest):
    return ["cells", {"where": {"regime": regime, "completeness": comp}}] + list(rest)


def cx(comp, *rest):
    if comp == 0.01:
        return [-1] + list(rest)
    return [{"where": {"completeness": comp}}] + list(rest)


def pc(region, metric="f1", block="conditional", comp="c1_x0", *rest):
    return [f"exclude|{comp}|{region}|{metric}", block] + list(rest)


# shared-reference shortcuts
def arm(a, *rest, space="ec"):
    return ["arms", a, space] + list(rest)


def con(c, *rest, space="ec"):
    return ["contrasts", c, space] + list(rest)


POP = ["population", "f1", "value"]
CM = ["conditional_matched", "f1", "value"]

ITEMS = {}


def item(iid, values=None, manual=None):
    assert iid not in ITEMS, iid
    ITEMS[iid] = {"values": values or [], "manual": manual}


# ================================================================ R4H post-hoc studies
OTA = ["confirm", "contrasts", "selected-shipped", "metrics"]
item("front:abs-objective-posthoc", [], manual="qualitative abstract sentence (no number printed at L61); the Stage O numbers "
     "are under methods:stageO_results")
WC = ["contrasts"]
item("methods:wobj_results", [
    V("p_func vs two-tier face-mid delta", WOBJ, WC + ["real-two_tier", "taf_face_mid", "multi", "mean_delta_weighted"], "0.0061", "absf4"),
    V("  CI bound (near)", WOBJ, WC + ["real-two_tier", "taf_face_mid", "multi", "ci_group", 1], "0.0025", "absf4"),
    V("  CI bound (far)", WOBJ, WC + ["real-two_tier", "taf_face_mid", "multi", "ci_group", 0], "0.0097", "absf4"),
    V("two-tier comparator", WOBJ, WC + ["real-two_tier", "taf_face_mid", "multi", "mean_comparator_weighted"], "0.130"),
    X("share of two-tier (4.7%)", "abs(a)/b", {"a": S(WOBJ, WC + ["real-two_tier", "taf_face_mid", "multi", "mean_delta_weighted"]),
                                              "b": S(WOBJ, WC + ["real-two_tier", "taf_face_mid", "multi", "mean_comparator_weighted"])}, "4.7", "pct1"),
    V("face minimum delta", WOBJ, WC + ["real-two_tier", "taf_face_min", "multi", "mean_delta_weighted"], "-0.0025", "f4"),
    V("  CI low", WOBJ, WC + ["real-two_tier", "taf_face_min", "multi", "ci_group", 0], "-0.0053", "f4"),
    V("  CI high", WOBJ, WC + ["real-two_tier", "taf_face_min", "multi", "ci_group", 1], "+0.0011", "+f4"),
    V("beat shuffled", WOBJ, WC + ["real-shuffled", "taf_face_mid", "multi", "mean_delta_weighted"], "0.0120", "absf4"),
    V("  CI (near)", WOBJ, WC + ["real-shuffled", "taf_face_mid", "multi", "ci_group", 1], "0.0069", "absf4"),
    V("  CI (far)", WOBJ, WC + ["real-shuffled", "taf_face_mid", "multi", "ci_group", 0], "0.0179", "absf4"),
    V("shuffled vs two-tier", WOBJ, WC + ["shuffled-two_tier", "taf_face_mid", "multi", "mean_delta_weighted"], "+0.0060", "+f4"),
    V("oracle vs two-tier", WOBJ, WC + ["oracle-two_tier", "taf_face_mid", "multi", "mean_delta_weighted"], "0.0249", "absf4"),
    V("  CI (near)", WOBJ, WC + ["oracle-two_tier", "taf_face_mid", "multi", "ci_group", 1], "0.0210", "absf4"),
    V("  CI (far)", WOBJ, WC + ["oracle-two_tier", "taf_face_mid", "multi", "ci_group", 0], "0.0291", "absf4"),
    V("p_func share of oracle", WOBJ, ["diagnostics", "multi", "oracle_headroom_captured"], "24", "pct0"),
    V("free gene-backed reactions two-tier", WOBJ, WC + ["real-two_tier", "faces.opt.n_var_1e-3", "multi", "mean_comparator_weighted"], "41.5", "f1"),
    V("free gene-backed reactions delta", WOBJ, WC + ["real-two_tier", "faces.opt.n_var_1e-3", "multi", "mean_delta_weighted"], "-8.9", "f1"),
    V("  CI low", WOBJ, WC + ["real-two_tier", "faces.opt.n_var_1e-3", "multi", "ci_group", 0], "-10.3", "f1"),
    V("  CI high", WOBJ, WC + ["real-two_tier", "faces.opt.n_var_1e-3", "multi", "ci_group", 1], "-7.5", "f1"),
    V("shuffled free reactions", WOBJ, WC + ["real-shuffled", "faces.opt.n_var_1e-3", "multi", "mean_comparator_weighted"], "27.9", "f1"),
    V("exchange-term n<=1e-6 spread", WOBJ, ["checks", "exchange_term_rel_spread_ex100", "n_le_1e-6"], "357", "int"),
    V("n=1 identical", WOBJ, ["checks", "n1_real_identical_two_tier", "identical"], "60", "int"),
], manual="exchange sum 98.9% of objective, weighted term median 0.11%, p_func 32.7 and parsimony 41.8 have no single stored key "
          "(arms.<arm>.<stratum>.internal_share medians; WEIGHTED_OBJECTIVE_RESULTS.md section 1)")
item("methods:msprob_released_equivalence", [
    V("released_check pass", WOBJ, ["checks", "released_check", "pass"], "360", "int"),
    V("released_check n", WOBJ, ["checks", "released_check", "n"], "360", "int"),
])
item("methods:stageO_design_and_split", [
    V("g_meas", OTM, ["g_meas"], "0.084"),
    V("tune half amplicons", OTM, ["split", "n_jobs_by_half", "tune"], "167", "int"),
    V("confirm half amplicons", OTM, ["split", "n_jobs_by_half", "confirm"], "193", "int"),
], manual="121/80 groups per half are distinct group_id counts in per_amplicon.jsonl (not stored)")
item("methods:stageO_results", [
    V("selected-shipped truth-absent flux", OT, OTA + ["taf_face_mid", "pooled", "A", "mean"], "0.0172", "absf4"),
    V("  CI (near)", OT, OTA + ["taf_face_mid", "pooled", "A", "ci", 1], "0.0118", "absf4"),
    V("  CI (far)", OT, OTA + ["taf_face_mid", "pooled", "A", "ci", 0], "0.0221", "absf4"),
    V("face minimum", OT, OTA + ["taf_face_min", "pooled", "A", "mean"], "0.0076", "absf4"),
    V("beat shuffled", OT, ["confirm", "contrasts", "selected-shuffled", "metrics", "taf_face_mid", "pooled", "A", "mean"], "0.0282", "absf4"),
    V("headroom selected", OT, ["confirm", "headroom_confirm", "A", "selected_own_cell"], "69", "pct0"),
    V("headroom shipped", OT, ["confirm", "headroom_confirm", "A", "shipped"], "23", "pct0"),
    V("free reactions selected", OT, ["confirm", "degeneracy", "selected-shipped", "faces.opt.n_var_1e-3", "pooled", "A", "arm_mean"], "13.6", "f1"),
    V("free reactions shipped", OT, ["confirm", "degeneracy", "selected-shipped", "faces.opt.n_var_1e-3", "pooled", "A", "comparator_mean"], "35.4", "f1"),
    V("free reactions shuffled", OT, ["confirm", "degeneracy", "selected-shuffled", "faces.opt.n_var_1e-3", "pooled", "A", "comparator_mean"], "2.1", "f1"),
    V("exchange flux", OT, OTA + ["flux_exchange", "pooled", "A", "mean"], "+2.5", "+f1"),
], manual="face-min CI, shuffled CI, shuffled-vs-released +0.0110, estimate-alone 0.0128, exchange-weight 0.0043, log-cost 0.0039, "
          "G effect, n=1 0.0048/0.0124, 96-of-193 group and 0.0164 sensitivity are in OBJECTIVE_TUNING_RESULTS.md / "
          "verification/MANUSCRIPT_VERIFICATION_7dc189e.md, not as stable summary keys")
item("methods:stageO_wrapper_equivalence", [
    V("amplicons passing (bitwise)", WEQ, ["summary", "n_ok"], "8", "int", tex_text="eight amplicons"),
    V("amplicons compared", WEQ, ["summary", "n_amplicons"], "8", "int", tex_text="eight amplicons"),
])
item("val16s:pfunc_calibration", [], manual="calibration/scores.json is 60 MB (git-ignored) and nests the WLS slope per unit/stratum; "
     "0.78 [0.76, 0.79] and +0.0007 are quoted from calibration/README.md (key results 1 and 5); add a selector once R4.3 regenerates it")
item("val16s:stage_e_estimator", [], manual="pfunc_estimators/scores.json rows are a list per estimator/stratum/fold; the 13 quoted "
     "numbers (0.0889 -> 0.0578, 0.209 -> 0.127, 0.543 -> 0.816, ...) are PFUNC_ESTIMATORS_RESULTS.md section 1 pooled values")
item("panels:disc.objective_test1", [
    V("lowers truth-absent flux", WOBJ, WC + ["real-two_tier", "taf_face_mid", "multi", "mean_delta_weighted"], "0.0061", "absf4"),
    V("  CI (near)", WOBJ, WC + ["real-two_tier", "taf_face_mid", "multi", "ci_group", 1], "0.0025", "absf4"),
    V("  CI (far)", WOBJ, WC + ["real-two_tier", "taf_face_mid", "multi", "ci_group", 0], "0.0097", "absf4"),
    X("free gene-backed reactions cut", "-a/b", {"a": S(WOBJ, WC + ["real-two_tier", "faces.opt.n_var_1e-3", "multi", "mean_delta_weighted"]),
                                                 "b": S(WOBJ, WC + ["real-two_tier", "faces.opt.n_var_1e-3", "multi", "mean_comparator_weighted"])}, "21", "pct0"),
    V("share of oracle", WOBJ, ["diagnostics", "multi", "oracle_headroom_captured"], "24", "pct0", tex_text="quarter"),
], manual="shuffled cut 33% and exchange 98.9% at the median: no single stored key")
item("panels:disc.objective_test2", [
    V("selected-shipped", OT, OTA + ["taf_face_mid", "pooled", "A", "mean"], "0.0172", "absf4"),
    V("  CI (near)", OT, OTA + ["taf_face_mid", "pooled", "A", "ci", 1], "0.0118", "absf4"),
    V("  CI (far)", OT, OTA + ["taf_face_mid", "pooled", "A", "ci", 0], "0.0221", "absf4"),
], manual="121 of 201 / 80 groups: per_amplicon.jsonl group counts")
item("panels:lim.reconstruction_tests", [], manual="design restatement (300/60 amplicons, 121 of 201 groups, slope 0.78, ~7,000 genomes): "
     "calibration slope from calibration/README.md; see methods:wobj_results, stageO_design_and_split")

# ================================================================ R2.3 cross-panel + MAG pairing
item("front:abs-mag-vs-16s", [
    X("MAG-16S conditional delta at 100% (low, FL)", "a-b", {"a": S(PC, pc("FullLength16S") + ["mag"]), "b": S(PC, pc("FullLength16S") + ["amp"])}, "+0.10", "+f2"),
    X("MAG-16S conditional delta at 100% (high, V4)", "a-b", {"a": S(PC, pc("V4") + ["mag"]), "b": S(PC, pc("V4") + ["amp"])}, "+0.21", "+f2"),
    V("FL delta at 1%", CX, cx(0.01, "FullLength16S_delta"), "+0.08", "+f2"),
    V("V4 delta at 1%", CX, cx(0.01, "V4_delta"), "+0.19", "+f2"),
    D("positive in every region at every level", "magvs16s_positive_everywhere", "True", "str"),
    V("resolve rate at 100% ('~30%')", CX, cx(1.0, "resolve_rate"), "0.3", "f1", note="printed as '~30%' of organisms", tex_text="30\\%"),
])
item("front:abs-threepanel-genecap-quality", [
    V("max |phystd delta| cond F1 (random V7-V9)", TP, tp("V7-V9", "exclude_0.987", "genecap_cond_f1", "delta_random_minus_original", "est"), "0.028", "absf3"),
])
item("front:abs-threepanel-genecap-coverage", [
    V("random-original genus coverage (min |.|, V1-V2: the abstract's '0.20')", TP, tp("V1-V2", "exclude_0.987", "genus_coverage", "delta_random_minus_original", "est"), "0.20", "absf2",
      note="the abstract's 'by 0.20 to 0.29' spans genus and gene-capture coverage; gene capture alone is 0.24-0.29"),
    V("random-original genecap coverage (max |.|, V4-V5_944R)", TP, tp("V4-V5_944R", "exclude_0.987", "genecap_coverage", "delta_random_minus_original", "est"), "0.29", "absf2"),
])
item("front:abs-selfincl-panel-dependence", [
    V("include genecap cond F1 random-original V4", TP, tp("V4", "include", "genecap_cond_f1", "delta_random_minus_original", "est"), "0.11", "absf2"),
    V("include genecap cond F1 random-original FL", TP, tp("FullLength16S", "include", "genecap_cond_f1", "delta_random_minus_original", "est"), "0.02", "absf2"),
])
MAGVS = {"FullLength16S": ("1,916", "0.813", "0.710", "+0.103", "+0.056", "0.304", "0.921", "0.248", "0.605"),
         "V1-V2": ("1,992", "0.814", "0.698", "+0.116", "+0.059", "0.298", "0.920", "0.242", "0.593"),
         "V1-V3": ("2,020", "0.813", "0.704", "+0.110", "+0.052", "0.298", "0.927", "0.243", "0.602"),
         "V3-V4": ("1,885", "0.815", "0.629", "+0.186", "+0.148", "0.298", "0.875", "0.242", "0.516"),
         "V4": ("2,299", "0.810", "0.598", "+0.212", "+0.171", "0.313", "0.863", "0.253", "0.483"),
         "V4-V5": ("2,334", "0.809", "0.616", "+0.192", "+0.157", "0.314", "0.868", "0.254", "0.504"),
         "V4-V5_944R": ("2,255", "0.811", "0.623", "+0.188", "+0.157", "0.318", "0.873", "0.257", "0.510"),
         "V7-V9": ("2,301", "0.808", "0.641", "+0.168", "+0.126", "0.326", "0.880", "0.264", "0.523")}


def magvs_values(with_cov=True):
    out = []
    for r, (n, mag, amp, dl, md, cm, ca, pm, pa_) in MAGVS.items():
        out += [V(f"{r} n", PC, pc(r) + ["n"], n, "int"),
                V(f"{r} MAG", PC, pc(r) + ["mag"], mag),
                V(f"{r} 16S", PC, pc(r) + ["amp"], amp),
                X(f"{r} delta", "a-b", {"a": S(PC, pc(r) + ["mag"]), "b": S(PC, pc(r) + ["amp"])}, dl, "+f3"),
                V(f"{r} median delta", PC, pc(r) + ["median_delta"], md, "+f3"),
                V(f"{r} pop MAG", PC, pc(r, block="population") + ["mag"], pm),
                V(f"{r} pop 16S", PC, pc(r, block="population") + ["amp"], pa_)]
        if with_cov:
            out += [D(f"{r} coverage MAG", f"magvs16s_{r}_cov_mag", cm),
                    D(f"{r} coverage 16S", f"magvs16s_{r}_cov_amp", ca)]
    return out


item("front:tab:magvs16s", magvs_values(), manual=None)
MAGVAL = {1.0: ("0.308", "0.807", "0.813", "0.710", "+0.103", "0.810", "0.598", "+0.212", "1,916"),
          0.5: ("0.307", "0.807", None, None, None, None, None, None, "1,927"),
          0.25: ("0.309", "0.806", None, None, None, None, None, None, "1,932"),
          0.1: ("0.295", "0.803", None, None, None, None, None, None, "1,845"),
          0.05: ("0.292", "0.797", None, None, None, None, None, None, "1,821"),
          0.02: ("0.283", "0.791", None, None, None, None, None, None, "1,791"),
          0.01: ("0.253", "0.783", "0.783", "0.706", "+0.077", "0.783", "0.593", "+0.190", "1,603")}


def magval_values():
    out = []
    for c, (rr, cf, fm, fa, fd, vm, va, vd, n) in MAGVAL.items():
        out += [V(f"{c:g} resolve rate", CX, cx(c, "resolve_rate"), rr), V(f"{c:g} cond F1", CX, cx(c, "cond_f1"), cf),
                V(f"{c:g} FL n", CX, cx(c, "FullLength16S_n"), n, "int")]
        if fm:
            out += [V(f"{c:g} FL MAG", CX, cx(c, "FullLength16S_mag"), fm), V(f"{c:g} FL 16S", CX, cx(c, "FullLength16S_amp"), fa),
                    V(f"{c:g} FL delta", CX, cx(c, "FullLength16S_delta"), fd, "+f3"),
                    V(f"{c:g} V4 MAG", CX, cx(c, "V4_mag"), vm), V(f"{c:g} V4 16S", CX, cx(c, "V4_amp"), va),
                    V(f"{c:g} V4 delta", CX, cx(c, "V4_delta"), vd, "+f3")]
    return out


item("front:tab:magval", magval_values() + [
    V("V4 pairing n (min over tab levels)", CX, [{"each": [0, 3, 5, 8, 10, 12, 13]}, "V4_n"], "1,893", "int", reduce="min",
      note="rows 1.0/0.5/0.25/0.1/0.05/0.02/0.01 by list position"),
    V("V4 pairing n (max over tab levels)", CX, [{"each": [0, 3, 5, 8, 10, 12, 13]}, "V4_n"], "2,311", "int", reduce="max"),
], manual="9,997 panel organisms is the panel size (selection_10k), not re-scored")
item("front:fig:magcomp", [], manual="figure (F14); its numbers are inventoried under tab:magval / tab:magvs16s; re-render in R5.2")
item("front:fig:magcross", [], manual="figure (F15); numbers in tab:magval; re-render in R5.2")
TP2 = {"FullLength16S": ("0.921", "0.750", "0.650", "0.657", "0.645", "0.641"),
       "V1-V3": ("0.927", "0.757", "0.664", "0.649", "0.641", "0.636"),
       "V4": ("0.863", "0.647", "0.536", "0.559", "0.541", "0.536"),
       "V7-V9": ("0.880", "0.659", "0.547", "0.594", "0.571", "0.565")}


def twopanel16s_values():
    out = []
    for r, vals in TP2.items():
        for i, p in enumerate(("original", "kmeans", "random")):
            out.append(V(f"{r} coverage {p}", TP, tp(r, "exclude_0.987", "genecap_coverage", p, "crude"), vals[i]))
        for i, p in enumerate(("original", "kmeans", "random")):
            out.append(V(f"{r} cond F1 {p}", TP, tp(r, "exclude_0.987", "genecap_cond_f1", p, "crude"), vals[3 + i]))
    return out


item("front:tab:twopanel-16S", twopanel16s_values())
TPM = {1.0: ("0.309", "0.342", "0.428", "0.807", "0.797", "0.791"),
       0.5: ("0.309", "0.343", "0.428", "0.806", "0.798", "0.792"),
       0.3: ("0.310", "0.345", "0.430", "0.804", "0.796", "0.795")}


def twopanelmag_values():
    out = []
    for c, vals in TPM.items():
        for i, p in enumerate(("original", "kmeans", "random")):
            out.append(V(f"{c:g} resolve {p}", PA, pa("exclude", c, p, "resolve_rate_NOT_COMPARABLE"), vals[i]))
        for i, p in enumerate(("original", "kmeans", "random")):
            out.append(V(f"{c:g} cond F1 {p}", PA, pa("exclude", c, p, "cond_f1_crude"), vals[3 + i]))
    return out


item("front:tab:twopanel-MAG", twopanelmag_values())


def emptied_values(cols=("n", "share", "unc")):
    out = []
    pv = {"original": ("7,215", ("84.4", "7.2", "3.8", "4.7"), ("0.037", "0.017", "0.015", "1.000")),
          "kmeans": ("8,815", ("44.0", "17.9", "15.2", "22.9"), ("0.028", "0.022", "0.032", "1.000")),
          "random": ("8,762", ("37.3", "14.8", "14.5", "33.4"), ("0.023", "0.021", "0.028", "1.000"))}
    for p, (n, shares, unc) in pv.items():
        if "n" in cols:
            out.append(V(f"{p} N", EM, ["panels", p, "N"], n, "int"))
        for s, sh, u in zip(("0-19", "20-49", "50-99", "100"), shares, unc):
            if "share" in cols:
                out.append(V(f"{p} {s} share", EM, ["panels", p, "strata", s, "share"], sh, "pct1"))
            if "unc" in cols:
                out.append(V(f"{p} {s} uncovered", EM, ["panels", p, "strata", s, "non_ok_rate"], u))
    return out


item("front:tab:emptied", emptied_values())
TP3 = {"FullLength16S": ("0.926", "0.767", "0.682", "-0.243", "-0.255", "-0.233", "-0.019", "-0.031", "-0.005"),
       "V4": ("0.855", "0.665", "0.576", "-0.279", "-0.291", "-0.268", "-0.024", "-0.031", "-0.017"),
       "V7-V9": ("0.870", "0.681", "0.594", "-0.275", "-0.287", "-0.263", "-0.028", "-0.035", "-0.020")}


def threepanel_values():
    out = []
    for r, v in TP3.items():
        for i, p in enumerate(("original", "kmeans", "random")):
            out.append(V(f"{r} phystd coverage {p}", TP, tp(r, "exclude_0.987", "genecap_coverage", p, "phystd"), v[i]))
        out += [V(f"{r} delta coverage", TP, tp(r, "exclude_0.987", "genecap_coverage", "delta_random_minus_original", "est"), v[3]),
                V(f"{r} delta coverage CI lo", TP, tp(r, "exclude_0.987", "genecap_coverage", "delta_random_minus_original", "ci", 0), v[4]),
                V(f"{r} delta coverage CI hi", TP, tp(r, "exclude_0.987", "genecap_coverage", "delta_random_minus_original", "ci", 1), v[5]),
                V(f"{r} delta cond F1", TP, tp(r, "exclude_0.987", "genecap_cond_f1", "delta_random_minus_original", "est"), v[6]),
                V(f"{r} delta cond F1 CI lo", TP, tp(r, "exclude_0.987", "genecap_cond_f1", "delta_random_minus_original", "ci", 0), v[7]),
                V(f"{r} delta cond F1 CI hi", TP, tp(r, "exclude_0.987", "genecap_cond_f1", "delta_random_minus_original", "ci", 1), v[8])]
    return out


item("front:tab:threepanel-genecap", threepanel_values())
PAR = {1.0: ("0.809", "0.788", "0.783", "-0.026", "-0.032", "-0.020"),
       0.3: ("0.806", "0.786", "0.789", "-0.017", "-0.023", "-0.011")}


def panelref_values():
    out = []
    for c, v in PAR.items():
        for i, p in enumerate(("original", "kmeans", "random")):
            out.append(V(f"{c:g} phystd cond F1 {p}", PA, pa("exclude", c, p, "cond_f1_phystd"), v[i]))
        out += [V(f"{c:g} delta random-original", PA, pa("exclude", c, "delta_cond_f1_random_minus_original", "est"), v[3]),
                V(f"{c:g} delta CI lo", PA, pa("exclude", c, "delta_cond_f1_random_minus_original", "ci", 0), v[4]),
                V(f"{c:g} delta CI hi", PA, pa("exclude", c, "delta_cond_f1_random_minus_original", "ci", 1), v[5])]
    return out


item("front:tab:threepanel-panelref", panelref_values())
item("val16s:emptied_vs_uncovered_panels", [
    V("FL emptied original", EM, ["panels", "original", "emptied_pct"], "4.7", "pp1"),
    V("FL emptied k-means", EM, ["panels", "kmeans", "emptied_pct"], "22.9", "pp1"),
    V("FL survivors uncovered original", EM, ["panels", "original", "survivors", "non_ok_rate"], "3.4", "pct1"),
    V("FL survivors uncovered k-means", EM, ["panels", "kmeans", "survivors", "non_ok_rate"], "2.8", "pct1"),
    D("pooled emptied original", "panels_emptied_pooled_original", "8.3", "pct1", note="no W6 entry at the time of writing"),
    D("pooled emptied k-means", "panels_emptied_pooled_kmeans", "29.4", "pct1", note="no W6 entry at the time of writing"),
    D("pooled uncovered survivors original", "panels_surv_uncovered_pooled_original", "3.1", "pct1"),
    D("pooled uncovered survivors k-means", "panels_surv_uncovered_pooled_kmeans", "2.4", "pct1"),
])
item("mag:abstract_mag_margin", [
    X("FL delta at 100%", "a-b", {"a": S(PC, pc("FullLength16S") + ["mag"]), "b": S(PC, pc("FullLength16S") + ["amp"])}, "+0.10", "+f2"),
    X("V4 delta at 100%", "a-b", {"a": S(PC, pc("V4") + ["mag"]), "b": S(PC, pc("V4") + ["amp"])}, "+0.21", "+f2"),
    V("FL delta at 1%", CX, cx(0.01, "FullLength16S_delta"), "+0.08", "+f2"),
    V("V4 delta at 1%", CX, cx(0.01, "V4_delta"), "+0.19", "+f2"),
    D("positive in every region at every level", "magvs16s_positive_everywhere", "True", "str"),
])
item("mag:magval_completeness_null", [
    *[V(f"resolve {c:g}", CX, cx(c, "resolve_rate"), rr) for c, rr in ((1.0, "0.308"), (0.5, "0.307"), (0.25, "0.309"), (0.1, "0.295"),
                                                                         (0.05, "0.292"), (0.02, "0.283"), (0.01, "0.253"))],
    *[V(f"cond F1 {c:g}", CX, cx(c, "cond_f1"), cf) for c, cf in ((1.0, "0.807"), (0.5, "0.807"), (0.25, "0.806"), (0.1, "0.803"),
                                                                    (0.05, "0.797"), (0.02, "0.791"), (0.01, "0.783"))],
    V("resolve 100-30% (low)", CX, [{"each": [0, 1, 2, 3, 4]}, "resolve_rate"], "0.306", reduce="min"),
    V("resolve 100-30% (high)", CX, [{"each": [0, 1, 2, 3, 4]}, "resolve_rate"], "0.308", reduce="max"),
    V("cond F1 100-30% (low)", CX, [{"each": [0, 1, 2, 3, 4]}, "cond_f1"], "0.805", reduce="min"),
    V("cond F1 100-30% (high)", CX, [{"each": [0, 1, 2, 3, 4]}, "cond_f1"], "0.807", reduce="max"),
])
item("mag:magval_paired_main", [
    X("FL lead at 100%", "a-b", {"a": S(PC, pc("FullLength16S") + ["mag"]), "b": S(PC, pc("FullLength16S") + ["amp"])}, "0.103"),
    X("V4 lead at 100%", "a-b", {"a": S(PC, pc("V4") + ["mag"]), "b": S(PC, pc("V4") + ["amp"])}, "0.212"),
    *[V(f"FL n {c:g}", CX, cx(c, "FullLength16S_n"), n, "int") for c, n in ((1.0, "1,916"), (0.5, "1,927"), (0.25, "1,932"), (0.1, "1,845"),
                                                                             (0.05, "1,821"), (0.02, "1,791"), (0.01, "1,603"))],
    V("FL MAG 1%", CX, cx(0.01, "FullLength16S_mag"), "0.783"), V("FL 16S 1%", CX, cx(0.01, "FullLength16S_amp"), "0.706"),
    V("V4 MAG 1%", CX, cx(0.01, "V4_mag"), "0.783"), V("V4 16S 1%", CX, cx(0.01, "V4_amp"), "0.593"),
    V("V4 delta 1%", CX, cx(0.01, "V4_delta"), "+0.190", "+f3"),
])
item("mag:magval_no_crossover", [
    V("FL lead at 1%", CX, cx(0.01, "FullLength16S_delta"), "0.077"),
    V("V4 lead at 1%", CX, cx(0.01, "V4_delta"), "0.190"),
    V("MAG wins", CX, cx(0.01, "FullLength16S_mag_wins"), "1,129", "int"),
    V("16S wins", CX, cx(0.01, "FullLength16S_amp_wins"), "473", "int"),
    V("p", CX, cx(0.01, "FullLength16S_p"), "1e-79", "e0"),
    X("cond F1 cost 100% -> 1%", "a-b", {"a": S(CX, cx(1.0, "cond_f1")), "b": S(CX, cx(0.01, "cond_f1"))}, "0.024"),
    V("resolve 100%", CX, cx(1.0, "resolve_rate"), "0.308"), V("resolve 1%", CX, cx(0.01, "resolve_rate"), "0.253"),
], manual="37 kb (median retained source sequence at 1%) is quoted from out/RESULTS_v2.md")


def delta_range(label, regions, metric, lo, hi):
    vars_ = {}
    for i, r in enumerate(regions):
        vars_[f"m{i}"] = S(PC, pc(r, metric) + ["mag"])
        vars_[f"a{i}"] = S(PC, pc(r, metric) + ["amp"])
    diffs = "[" + ", ".join(f"m{i}-a{i}" for i in range(len(regions))) + "]"
    return [X(f"{label} (low)", f"min({diffs})", vars_, lo), X(f"{label} (high)", f"max({diffs})", vars_, hi)]


item("mag:magval_vs_every_region", [
    *delta_range("F1 delta FL/V1-V2/V1-V3", V1A, "f1", "+0.103", "+0.116"),
    *delta_range("F1 delta V3-onward", V3ON, "f1", "+0.168", "+0.212"),
    V("MAG recall (low)", PC, [{"each": [f"exclude|c1_x0|{r}|recall" for r in R8]}, "conditional", "mag"], "0.816", reduce="min"),
    V("MAG recall (high)", PC, [{"each": [f"exclude|c1_x0|{r}|recall" for r in R8]}, "conditional", "mag"], "0.824", reduce="max"),
    V("16S recall (low)", PC, [{"each": [f"exclude|c1_x0|{r}|recall" for r in R8]}, "conditional", "amp"], "0.703", reduce="min"),
    V("16S recall (high)", PC, [{"each": [f"exclude|c1_x0|{r}|recall" for r in R8]}, "conditional", "amp"], "0.768", reduce="max"),
    *delta_range("recall delta FL/V1", V1A, "recall", "0.055", "0.083"),
    *delta_range("recall delta V3-onward", V3ON, "recall", "0.097", "0.114"),
    V("max F1 p", PC, [{"each": [f"exclude|c1_x0|{r}|f1" for r in R8]}, "conditional", "p"], "1e-161", "lt", reduce="max", skip_null=True,
      note="V4, V4-V5, V4-V5_944R store p=null (underflow)"),
    V("max recall p", PC, [{"each": [f"exclude|c1_x0|{r}|recall" for r in R8]}, "conditional", "p"], "1e-29", "lt", reduce="max"),
])
item("mag:tab_magvs16s", magvs_values())
item("mag:magval_population_view", [
    V("FL pop 16S", PC, pc("FullLength16S", block="population") + ["amp"], "0.605"),
    V("FL pop MAG", PC, pc("FullLength16S", block="population") + ["mag"], "0.248"),
    V("V4 pop 16S", PC, pc("V4", block="population") + ["amp"], "0.483"),
    V("V4 pop MAG", PC, pc("V4", block="population") + ["mag"], "0.253"),
    X("amplicon coverage (low)", "min([" + ", ".join(f"c{i}" for i in range(8)) + "])",
      {f"c{i}": {"artifact": f"derived:magvs16s_{r}_cov_amp"} for i, r in enumerate(R8)}, "86", "pct0"),
    X("amplicon coverage (high)", "max([" + ", ".join(f"c{i}" for i in range(8)) + "])",
      {f"c{i}": {"artifact": f"derived:magvs16s_{r}_cov_amp"} for i, r in enumerate(R8)}, "93", "pct0"),
])
item("mag:magval_af_gate_counterfactual", [
    D("qualifying share at gate 0.15, floor 0.83", "mag_af_shipped", "0.311"),
    D("qualifying share with floor 0", "mag_af_floor0", "0.316"),
    D("qualifying share at gate 0.02", "mag_af_gate2", "0.452"),
    V("newly admitted with prediction", AFG, ["summary", "n_newly_admitted"], "1,367", "int"),
    V("mean F1 newly admitted", AFG, ["summary", "mean_f1_newly"], "0.599"),
    V("mean F1 already answered", AFG, ["summary", "mean_f1_already"], "0.807"),
    V("resolve rate", CX, cx(1.0, "resolve_rate"), "0.308"),
], manual="1,386 newly admitted, 26.7% above 0.710, +45% coverage, pooled 0.742: af_gate_counterfactual.json summary sub-keys (check names after R2.3)")
item("mag:fig_magcomp_magcross", [], manual="figures F14/F15; numbers inventoried in tab:magval / tab:magvs16s; re-render in R5.2")
item("mag:twopanel_mag_panelarm", twopanelmag_values()[:6] + panelref_values()[:6])
item("mag:twopanel_mag_abstention_strata", [
    D("singleton abstention k-means", "panels_mag_single_abst_kmeans", "0.6895", "f4"),
    D("singleton abstention original", "panels_mag_single_abst_original", "0.6883", "f4"),
    D(">=20 members abstention k-means", "panels_mag_ge20_abst_kmeans", "0.559"),
    D("overall k-means", "panels_mag_species_level_abst_kmeans", "0.658"),
    D("overall original", "panels_mag_species_level_abst_original", "0.688"),
    D("k-means without >=20 stratum", "panels_mag_without_ge20_abst_kmeans", "0.687"),
])
item("panels:panels.condF1_crude_bounds", [
    X("random-original crude cond F1 V7-V9", "a-b", {"a": S(TP, tp("V7-V9", "exclude_0.987", "genecap_cond_f1", "random", "crude")),
                                                     "b": S(TP, tp("V7-V9", "exclude_0.987", "genecap_cond_f1", "original", "crude"))}, "0.029", "absf3"),
    X("k-means-original crude cond F1 V7-V9", "a-b", {"a": S(TP, tp("V7-V9", "exclude_0.987", "genecap_cond_f1", "kmeans", "crude")),
                                                      "b": S(TP, tp("V7-V9", "exclude_0.987", "genecap_cond_f1", "original", "crude"))}, "0.023", "absf3"),
], manual="'within 0.016 / 0.010 at every MAG completeness' and 'original ahead in all eight regions and five levels' span every cell; "
          "check by eye in number_diff (panelarm cond_f1_crude)")
item("panels:panels.coverage_disagrees", [
    V("MAG resolve k-means 100%", PA, pa("exclude", 1.0, "kmeans", "resolve_rate_NOT_COMPARABLE"), "0.342"),
    V("MAG resolve random 100%", PA, pa("exclude", 1.0, "random", "resolve_rate_NOT_COMPARABLE"), "0.428"),
    V("MAG resolve original 100%", PA, pa("exclude", 1.0, "original", "resolve_rate_NOT_COMPARABLE"), "0.309"),
    V("FL coverage k-means", TP, tp("FullLength16S", "exclude_0.987", "genecap_coverage", "kmeans", "crude"), "0.750"),
    V("FL coverage random", TP, tp("FullLength16S", "exclude_0.987", "genecap_coverage", "random", "crude"), "0.650"),
    V("FL coverage original", TP, tp("FullLength16S", "exclude_0.987", "genecap_coverage", "original", "crude"), "0.921"),
])
item("panels:panels.tab_twopanel_16S", twopanel16s_values())
item("panels:panels.tab_twopanel_MAG", twopanelmag_values())
item("panels:panels.fl_gap_standardized", [
    V("phystd k-means-original FL coverage", TP, tp("FullLength16S", "exclude_0.987", "genecap_coverage", "delta_kmeans_minus_original", "est"), "0.159", "absf3"),
    X("crude k-means-original", "a-b", {"a": S(TP, tp("FullLength16S", "exclude_0.987", "genecap_coverage", "kmeans", "crude")),
                                        "b": S(TP, tp("FullLength16S", "exclude_0.987", "genecap_coverage", "original", "crude"))}, "0.171", "absf3"),
    V("common phyla", TP, ["meta", "n_common_phyla"], "113", "int"),
])
item("panels:panels.surviving_uncovered", [
    V("FL survivors uncovered original", EM, ["panels", "original", "survivors", "non_ok_rate"], "3.4", "pct1"),
    V("FL survivors uncovered k-means", EM, ["panels", "kmeans", "survivors", "non_ok_rate"], "2.8", "pct1"),
    V("FL survivors uncovered random", EM, ["panels", "random", "survivors", "non_ok_rate"], "2.3", "pct1"),
    V("FL strata uncovered (low)", EM, ["panels", {"each": ["original", "kmeans", "random"]}, "strata"], "1.5", "pct1", reduce=None,
      note="placeholder replaced below"),
    D("pooled original", "panels_surv_uncovered_pooled_original", "3.1", "pct1"),
    D("pooled k-means", "panels_surv_uncovered_pooled_kmeans", "2.4", "pct1"),
    D("pooled random", "panels_surv_uncovered_pooled_random", "1.8", "pct1"),
    D("emptied / surviving uncovered rate (low, fold)", "panels_emptied_fold_min", "29", "f0"),
    D("emptied / surviving uncovered rate (high, fold)", "panels_emptied_fold_max", "43", "f0"),
    D("pooled spread (points)", "panels_surv_uncovered_pooled_spread_pts", "1.3", "f1"),
    X("original-random FL difference (points)", "a-b", {"a": S(EM, ["panels", "original", "survivors", "non_ok_rate"]),
                                                       "b": S(EM, ["panels", "random", "survivors", "non_ok_rate"])}, "1.1", "pct1"),
    V("CMH original vs random p", EM, ["tests", "cmh_survivor_strata", "original_vs_random", "p"], "0.004"),
    V("CMH original vs k-means p", EM, ["tests", "cmh_survivor_strata", "original_vs_kmeans", "p"], "0.15", "f2"),
], manual="'29- to 43-fold below the emptied-record rate' is a ratio over strata; '1.3 pooled' from the derived pooled rates")
# replace the strata placeholder with explicit low/high over the 9 non-emptied FL strata
_strata = [S(EM, ["panels", p, "strata", s, "non_ok_rate"]) for p in ("original", "kmeans", "random") for s in ("0-19", "20-49", "50-99")]
_sv = {f"s{i}": v for i, v in enumerate(_strata)}
_lst = "[" + ", ".join(_sv) + "]"
ITEMS["panels:panels.surviving_uncovered"]["values"][3] = X("FL strata uncovered (low)", f"min({_lst})", _sv, "1.5", "pct1")
ITEMS["panels:panels.surviving_uncovered"]["values"].insert(4, X("FL strata uncovered (high)", f"max({_lst})", _sv, "3.7", "pct1"))
item("panels:panels.emptied_arithmetic", [
    V("emptied random", EM, ["panels", "random", "emptied_pct"], "33.4", "pp1"),
    V("emptied k-means", EM, ["panels", "kmeans", "emptied_pct"], "22.9", "pp1"),
    V("emptied original", EM, ["panels", "original", "emptied_pct"], "4.7", "pp1"),
    X("1-emptied random", "1-a/100", {"a": S(EM, ["panels", "random", "emptied_pct"])}, "0.666"),
    X("1-emptied k-means", "1-a/100", {"a": S(EM, ["panels", "kmeans", "emptied_pct"])}, "0.771"),
    X("1-emptied original", "1-a/100", {"a": S(EM, ["panels", "original", "emptied_pct"])}, "0.953"),
    V("observed random", TP, tp("FullLength16S", "exclude_0.987", "genecap_coverage", "random", "crude"), "0.650"),
    V("observed k-means", TP, tp("FullLength16S", "exclude_0.987", "genecap_coverage", "kmeans", "crude"), "0.750"),
    V("observed original", TP, tp("FullLength16S", "exclude_0.987", "genecap_coverage", "original", "crude"), "0.921"),
])
item("panels:panels.surviving_pred_empty", [
    X(f"{p} non-strict uncovered survivors", "a-b", {"a": S(EM, ["panels", p, "survivors", "non_ok"]),
                                                     "b": S(EM, ["panels", p, "survivors", "strict"])}, n, "int")
    for p, n in (("original", "107"), ("kmeans", "76"), ("random", "62"))])
item("panels:panels.tab_emptied_shares", emptied_values(("n", "share")))
item("panels:panels.tab_emptied_uncovered", emptied_values(("unc",)))


_gv = {}
for i, (ep, r) in enumerate([(e, r) for e in ("genus_coverage", "genecap_coverage") for r in R8]):
    _gv[f"d{i}"] = S(TP, tp(r, "exclude_0.987", ep, "delta_random_minus_original", "est"))
    _gv[f"r{i}"] = S(TP, tp(r, "exclude_0.987", ep, "random", "crude"))
    _gv[f"o{i}"] = S(TP, tp(r, "exclude_0.987", ep, "original", "crude"))
_gl = "[" + ", ".join(f"d{i}/(r{i}-o{i})" for i in range(16)) + "]"
item("panels:panels.std_gap_share", [
    X("standardized / crude gap (low)", f"min({_gl})", _gv, "76", "pct0"),
    X("standardized / crude gap (high)", f"max({_gl})", _gv, "92", "pct0"),
], manual="'interval excludes zero' in all 16: delta_random_minus_original.excludes_zero")
item("panels:panels.cond_f1_random_gap", [
    X("phystd random-original cond F1 (min |.|, V1-V2)", "abs(a)", {"a": S(TP, tp("V1-V2", "exclude_0.987", "genecap_cond_f1", "delta_random_minus_original", "est"))}, "0.013"),
    X("phystd random-original cond F1 (max |.|, V7-V9)", "abs(a)", {"a": S(TP, tp("V7-V9", "exclude_0.987", "genecap_cond_f1", "delta_random_minus_original", "est"))}, "0.028"),
])
item("panels:panels.tab_threepanel_genecap_block", threepanel_values())
item("panels:panels.tab_threepanel_magpanel_block", panelref_values())
_tax = []
for r, lo, hi in (("FullLength16S", "-0.251", "-0.189"), ("V1-V3", "-0.238", "-0.191"), ("V3-V4", "-0.274", "-0.215"),
                  ("V4", "-0.276", "-0.238"), ("V7-V9", "-0.277", "-0.237")):
    _tax += [V(f"genus coverage {r} CI lo", TP, tp(r, "exclude_0.987", "genus_coverage", "delta_random_minus_original", "ci", 0), lo),
             V(f"genus coverage {r} CI hi", TP, tp(r, "exclude_0.987", "genus_coverage", "delta_random_minus_original", "ci", 1), hi)]
for r, lo, hi in (("FullLength16S", "+0.013", "+0.061"), ("V1-V3", "+0.039", "+0.087"), ("V3-V4", "-0.040", "+0.019"),
                  ("V4", "-0.015", "+0.046"), ("V7-V9", "-0.044", "+0.011")):
    _tax += [V(f"genus cond acc {r} CI lo", TP, tp(r, "exclude_0.987", "genus_cond_acc", "delta_random_minus_original", "ci", 0), lo, "+f3"),
             V(f"genus cond acc {r} CI hi", TP, tp(r, "exclude_0.987", "genus_cond_acc", "delta_random_minus_original", "ci", 1), hi, "+f3")]
item("panels:panels.taxonomy_bootstrap_CIs", _tax, manual="the region order of the tab:threepanel genus block is assumed FL, V1-V3, V3-V4, V4, V7-V9 "
     "(the table's rows); re-check against the tex if a value does not match")
item("panels:panels.selfincl_genecap_penalty", [
    V("V4 original phystd", TP, tp("V4", "include", "genecap_cond_f1", "original", "phystd"), "0.837"),
    V("V4 random phystd", TP, tp("V4", "include", "genecap_cond_f1", "random", "phystd"), "0.730"),
    V("V4 delta", TP, tp("V4", "include", "genecap_cond_f1", "delta_random_minus_original", "est"), "-0.107"),
    V("V4 delta CI lo", TP, tp("V4", "include", "genecap_cond_f1", "delta_random_minus_original", "ci", 0), "-0.114"),
    V("V4 delta CI hi", TP, tp("V4", "include", "genecap_cond_f1", "delta_random_minus_original", "ci", 1), "-0.100"),
    V("V7-V9 original phystd", TP, tp("V7-V9", "include", "genecap_cond_f1", "original", "phystd"), "0.906"),
    V("V7-V9 random phystd", TP, tp("V7-V9", "include", "genecap_cond_f1", "random", "phystd"), "0.805"),
    V("FL delta", TP, tp("FullLength16S", "include", "genecap_cond_f1", "delta_random_minus_original", "est"), "-0.016"),
])
item("panels:panels.what_this_licenses", [
    V("genus coverage random-original (min |.|)", TP, tp("V1-V2", "exclude_0.987", "genus_coverage", "delta_random_minus_original", "est"), "0.20", "absf2"),
    V("genecap coverage random-original (max |.|)", TP, tp("V4-V5_944R", "exclude_0.987", "genecap_coverage", "delta_random_minus_original", "est"), "0.29", "absf2"),
    V("GTDB resolve-rate shift", GT, ["cells", {"where": {"regime": "exclude", "completeness": 1.0}}, "delta_random_minus_original", "est"], "0.14", "absf2"),
    V("max |pairwise cond F1 delta| (random V7-V9)", TP, tp("V7-V9", "exclude_0.987", "genecap_cond_f1", "delta_random_minus_original", "est"), "0.028", "absf3"),
])

# ================================================================ R1.3U union POST
HL = ["_meta", "headline_macro_recall_pct"]
HR = ["_meta", "headline_macro_recall_pct_resolved", "exclude_0.987"]
HP = ["_meta", "headline_macro_recall_pct_population", "exclude_0.987"]
item("front:abs-selfincl-genecap", [
    V("include capture (min)", GC, HL + ["include", {"each": R8}], "99", "gt", reduce="min"),
    V("V4 include precision", GC, ["regions", "V4", "include", "macro", "precision", "value"], "0.79", "f2"),
])
item("front:abs-selfexcl-genecap", [
    V("resolved capture (min)", GC, HR + [{"each": R8}], "66", "pp0", reduce="min"),
    V("resolved capture (max)", GC, HR + [{"each": R8}], "70", "pp0", reduce="max"),
    V("population capture (min)", GC, HP + [{"each": R8}], "58", "pp0", reduce="min"),
    V("population capture (max)", GC, HP + [{"each": R8}], "64", "pp0", reduce="max"),
])
TIERF = {"self_mode": "exclude_0.987", "outcome": "ok", "tier": "family"}
TIERG = {"self_mode": "exclude_0.987", "outcome": "ok", "tier": "genus"}
item("front:abs-tier-genecap", [
    CSV("family-tier mean recall", GCC, "csv_mean", {"filter": TIERF, "column": "recall", "usecols": ["self_mode", "outcome", "tier", "recall"]}, "58", "pct0"),
    CSV("genus-tier mean recall", GCC, "csv_mean", {"filter": TIERG, "column": "recall", "usecols": ["self_mode", "outcome", "tier", "recall"]}, "69", "pct0"),
    V("family-tier by region (min)", GC, ["regions", {"each": R8}, "exclude_0.987", "by_tier", "family", "macro", "recall", "value"], "56", "pct0", reduce="min"),
    V("family-tier by region (max)", GC, ["regions", {"each": R8}, "exclude_0.987", "by_tier", "family", "macro", "recall", "value"], "61", "pct0", reduce="max"),
])
item("front:tab:main-genecap", [
    V("include (min)", GC, HL + ["include", {"each": R8}], "99.1", "pp1", reduce="min"),
    V("include (max)", GC, HL + ["include", {"each": R8}], "99.98", "pp2", reduce="max"),
    V("excl cond (min)", GC, HR + [{"each": R8}], "66.4", "pp1", reduce="min"),
    V("excl cond (max)", GC, HR + [{"each": R8}], "69.9", "pp1", reduce="max"),
    V("excl pop (min)", GC, HP + [{"each": R8}], "57.6", "pp1", reduce="min"),
    V("excl pop (max)", GC, HP + [{"each": R8}], "64.4", "pp1", reduce="max"),
])
item("front:fig:genecap", [
    V("pop (min)", GC, HP + [{"each": R8}], "58", "pp0", reduce="min"),
    V("pop (max)", GC, HP + [{"each": R8}], "64", "pp0", reduce="max"),
    V("resolved (min)", GC, HR + [{"each": R8}], "66", "pp0", reduce="min"),
    V("resolved (max)", GC, HR + [{"each": R8}], "70", "pp0", reduce="max"),
], manual="figure F4b; re-render in R5.2")
item("front:fig:leakage-genecap", [
    V("gap (min)", GC, ["regions", {"each": R8}, "leakage_gap_macro_recall"], "30", "pct0", reduce="min"),
    V("gap (max)", GC, ["regions", {"each": R8}, "leakage_gap_macro_recall"], "33", "pct0", reduce="max"),
])
_tm = [V(f"include {r}", GC, HL + ["include", r], p, "pp" + str(len(p.split(".")[1]))) for r, p in (("FullLength16S", "99.98"), ("V4", "99.1"))]
_tm += [V(f"excl cond {r}", GC, HR + [r], p, "pp1") for r, p in (("V1-V2", "69.9"), ("V1-V3", "69.5"), ("FullLength16S", "68.7"), ("V4-V5", "66.4"))]
_tm += [V(f"excl pop {r}", GC, HP + [r], p, "pp1") for r, p in (("FullLength16S", "63.3"), ("V4", "58.1"), ("V4-V5", "57.6"))]
_tm += [V("excl pop V3-onward (max)", GC, HP + [{"each": V3ON}], "59.2", "pp1", reduce="max"),
        V("excl pop FL/V1 (min)", GC, HP + [{"each": V1A}], "63.3", "pp1", reduce="min"),
        V("excl pop FL/V1 (max)", GC, HP + [{"each": V1A}], "64.4", "pp1", reduce="max")]
item("val16s:tabmain_genecap", _tm)
item("val16s:genecap_leakage_gap", [
    V("gap V1-V2", GC, ["regions", "V1-V2", "leakage_gap_macro_recall"], "29.6", "pct1"),
    V("gap V4-V5", GC, ["regions", "V4-V5", "leakage_gap_macro_recall"], "32.9", "pct1"),
    V("gap (min)", GC, ["regions", {"each": R8}, "leakage_gap_macro_recall"], "30", "pct0", reduce="min"),
    V("gap (max)", GC, ["regions", {"each": R8}, "leakage_gap_macro_recall"], "33", "pct0", reduce="max"),
])
item("val16s:genecap_by_tier", [
    CSV("family-tier mean recall", GCC, "csv_mean", {"filter": TIERF, "column": "recall", "usecols": ["self_mode", "outcome", "tier", "recall"]}, "58", "pct0"),
    CSV("genus-tier mean recall", GCC, "csv_mean", {"filter": TIERG, "column": "recall", "usecols": ["self_mode", "outcome", "tier", "recall"]}, "69", "pct0"),
    V("family-tier by region (min)", GC, ["regions", {"each": R8}, "exclude_0.987", "by_tier", "family", "macro", "recall", "value"], "56", "pct0", reduce="min"),
    V("family-tier by region (max)", GC, ["regions", {"each": R8}, "exclude_0.987", "by_tier", "family", "macro", "recall", "value"], "61", "pct0", reduce="max"),
])
item("mag:methods_union_size_amplicon", [
    CSV("include ok share n_selected>1", GCC, "csv_share", {"filter": {"self_mode": "include", "outcome": "ok"}, "where": {"n_selected": [">", 1]},
                                                            "usecols": ["self_mode", "outcome", "n_selected"]}, "34.0", "pct1"),
])
_pr = ["regions", {"each": V3ON}, "include", "macro", "precision", "value"]
_pr1 = ["regions", {"each": ["V1-V2", "V1-V3"]}, "include", "macro", "precision", "value"]
item("panels:disc.region_tradeoff_genecap", [
    V("include capture (min)", GC, HL + ["include", {"each": R8}], "99", "gt", reduce="min"),
    V("include precision V3-onward (min)", GC, _pr, "0.79", "f2", reduce="min"),
    V("include precision V3-onward (max)", GC, _pr, "0.88", "f2", reduce="max"),
    V("include precision V1-V2/V1-V3 (min)", GC, _pr1, "0.92", "f2", reduce="min"),
    V("include precision V1-V2/V1-V3 (max)", GC, _pr1, "0.94", "f2", reduce="max"),
    V("excl pop V3-onward (min)", GC, HP + [{"each": V3ON}], "57.6", "pp1", reduce="min"),
    V("excl pop V3-onward (max)", GC, HP + [{"each": V3ON}], "59.2", "pp1", reduce="max"),
    V("excl pop FL/V1 (min)", GC, HP + [{"each": V1A}], "63.3", "pp1", reduce="min"),
    V("excl pop FL/V1 (max)", GC, HP + [{"each": V1A}], "64.4", "pp1", reduce="max"),
])

# ================================================================ R2.5 d20 as-shipped, ceiling, probe
item("front:abs-depth-ec", [
    V("d100 cond EC F1", CS, csm("prfba_exclude", "conditional", "f1", "value"), "0.802"),
    V("d20 cond EC F1", CS20, csm("prfba_exclude", "conditional", "f1", "value"), "0.798"),
    V("probe recovered EC F1", PROBE, ["ec_space_of_recovered", "f1", "mean"], "0.783"),
    V("probe recovered precision", PROBE, ["ec_space_of_recovered", "precision", "mean"], "0.747"),
    V("d20 answered precision", CS20, csm("prfba_exclude", "conditional", "precision", "value"), "0.790"),
])
item("front:tab:picrust2", [
    *[V(f"All PICRUSt2 {k}", CS, csm("picrust2", *p), v) for k, p, v in (
        ("coverage", ["coverage"], "1.000"), ("recall", ["conditional", "recall", "value"], "0.849"),
        ("precision", ["conditional", "precision", "value"], "0.669"), ("F1", ["conditional", "f1", "value"], "0.729"),
        ("pop F1", ["population", "f1", "value"], "0.729"))],
    *[V(f"All prFBA excl {k}", CS, csm("prfba_exclude", *p), v) for k, p, v in (
        ("coverage", ["coverage"], "0.889"), ("recall", ["conditional", "recall", "value"], "0.847"),
        ("precision", ["conditional", "precision", "value"], "0.787"), ("F1", ["conditional", "f1", "value"], "0.802"),
        ("pop F1", ["population", "f1", "value"], "0.713"))],
    *[V(f"All prFBA incl {k}", CS, csm("prfba_include", *p), v) for k, p, v in (
        ("coverage", ["coverage"], "1.000"), ("recall", ["conditional", "recall", "value"], "0.999"),
        ("precision", ["conditional", "precision", "value"], "0.936"), ("F1", ["conditional", "f1", "value"], "0.959"),
        ("pop F1", ["population", "f1", "value"], "0.959"))],
    V("n scoreable", CS, csm("picrust2", "n"), "67,038", "int"),
    V("ceiling recall", CEIL, ["ceiling", "shared", "recall", "mean"], "0.904"),
    V("ceiling precision", CEIL, ["ceiling", "shared", "precision", "mean"], "0.743"),
    V("ceiling n genomes", CEIL, ["ceiling", "shared", "recall", "n"], "2,689", "int"),
], manual="in-reference / absent strata rows are under picrust:AS-20; 43% share (more precise than the reference annotation) has no stored key")
item("val16s:depth_probe_quality", [
    V("recovered EC recall", PROBE, ["ec_space_of_recovered", "recall", "mean"], "0.854"),
    V("recovered EC precision", PROBE, ["ec_space_of_recovered", "precision", "mean"], "0.747"),
    V("recovered EC F1", PROBE, ["ec_space_of_recovered", "f1", "mean"], "0.783"),
    V("d20 answered recall", CS20, csm("prfba_exclude", "conditional", "recall", "value"), "0.835"),
    V("d20 answered precision", CS20, csm("prfba_exclude", "conditional", "precision", "value"), "0.790"),
    V("d20 answered F1", CS20, csm("prfba_exclude", "conditional", "f1", "value"), "0.798"),
    V("recovered family-correct", PROBE, ["anchor_accuracy_of_recovered", "family", "correct"], "80.2", "pct1"),
    V("recovered genus-correct", PROBE, ["anchor_accuracy_of_recovered", "genus", "correct"], "60.1", "pct1"),
], manual="0.651 / 0.800 answered-set anchor accuracy: no stored key located")
item("val16s:conditional_ec_d20_vs_d100", [
    V("d100 F1", CS, csm("prfba_exclude", "conditional", "f1", "value"), "0.802"),
    V("d20 F1", CS20, csm("prfba_exclude", "conditional", "f1", "value"), "0.798"),
    V("d100 precision", CS, csm("prfba_exclude", "conditional", "precision", "value"), "0.787"),
    V("d20 precision", CS20, csm("prfba_exclude", "conditional", "precision", "value"), "0.790"),
    V("d100 recall", CS, csm("prfba_exclude", "conditional", "recall", "value"), "0.847"),
    V("d20 recall", CS20, csm("prfba_exclude", "conditional", "recall", "value"), "0.835"),
])
item("val16s:population_ec_and_projections", [
    V("d20 pop EC F1", CS20, csm("prfba_exclude", "population", "f1", "value"), "0.600"),
    V("d100 pop EC F1", CS, csm("prfba_exclude", "population", "f1", "value"), "0.713"),
    V("PICRUSt2 pop EC F1", CS, csm("picrust2", "population", "f1", "value"), "0.729"),
], manual="projections onto the 15,026 records (~0.96 / ~0.49 / ~0.75) are computed in picrust2_comparison/README.md, not stored")
item("picrust:AS-03", [
    V("shared EC universe", CS, ["_meta", "ec_universe_shared"], "2,132", "int"),
    V("shared rxn universe", CS, ["_meta", "rxn_universe_shared"], "7,042", "int"),
])
item("picrust:AS-07", [
    V("in-reference scoreable genomes", CEIL, ["ceiling", "shared", "recall", "n"], "2,689", "int"),
    V("in-reference amplicons", CS, strat("True", "picrust2", "n"), "19,686", "int"),
    V("source-absent amplicons", CS, strat("False", "picrust2", "n"), "47,352", "int"),
    X("in-reference share", "a/(a+b)", {"a": S(CS, strat("True", "picrust2", "n")), "b": S(CS, strat("False", "picrust2", "n"))}, "29", "pct0"),
], manual="2,720 of 9,746 source genomes in PICRUSt2's reference: no stored key")
item("picrust:AS-19", [
    V("ceiling recall", CEIL, ["ceiling", "shared", "recall", "mean"], "90.4", "pct1"),
    V("ceiling precision", CEIL, ["ceiling", "shared", "precision", "mean"], "74.3", "pct1"),
    V("ceiling Jaccard", CEIL, ["ceiling", "shared", "jaccard", "mean"], "0.688"),
    X("in-ref recall / ceiling recall", "a/b", {"a": S(CS, strat("True", "picrust2", "conditional", "recall", "value")),
                                                "b": S(CEIL, ["ceiling", "shared", "recall", "mean"])}, "96", "pct0"),
    V("PICRUSt2 in-ref recall", CS, strat("True", "picrust2", "conditional", "recall", "value"), "0.868"),
], manual="43% of in-reference amplicons more precise than the reference annotation: per_genome x comparison_cells, no stored key")
item("picrust:XR-02", [
    V("d20 pop EC F1", CS20, csm("prfba_exclude", "population", "f1", "value"), "0.600"),
    V("d100 pop EC F1", CS, csm("prfba_exclude", "population", "f1", "value"), "0.713"),
    V("PICRUSt2 pop EC F1", CS, csm("picrust2", "population", "f1", "value"), "0.729"),
    X("d20 population gap", "a-b", {"a": S(CS20, csm("prfba_exclude", "population", "f1", "value")),
                                    "b": S(CS20, csm("picrust2", "population", "f1", "value"))}, "-0.129"),
])
item("panels:lim.depth_probe_quality", [
    V("recovered genus-correct", PROBE, ["anchor_accuracy_of_recovered", "genus", "correct"], "60.1", "pct1"),
    V("recovered family-correct", PROBE, ["anchor_accuracy_of_recovered", "family", "correct"], "80.2", "pct1"),
    V("recovered EC F1", PROBE, ["ec_space_of_recovered", "f1", "mean"], "0.783"),
])
item("panels:lim.depth20_vs_100", [
    V("d100 cond EC F1", CS, csm("prfba_exclude", "conditional", "f1", "value"), "0.802"),
    V("d20 cond EC F1", CS20, csm("prfba_exclude", "conditional", "f1", "value"), "0.798"),
], manual="legacy genus accuracy within 0.01, macro recall +0.7..+1.3 and gene-capture F1 -0.1..-1.2 points: values under "
          "val16s:legacy_genecap_d20_vs_d100 (R1.4)")

# ================================================================ R2.4 d100 as-shipped
OC = {"usecols": ["method", "space", "outcome", "n_reps", "nsti_bin", "best_nonself_identity", "amplicon", "n_pred", "n_truth",
                  "precision", "nsti", "src_in_picrust2_ref", "domain", "region", "f1"]}


def cc(label, recipe, filt, printed, fmt="int", **kw):
    a = dict(OC)
    a["filter"] = filt
    a.update({k: kw.pop(k) for k in list(kw) if k in ("column", "where", "within", "derive")})
    return CSV(label, CC, recipe, a, printed, fmt, **kw)


EXC = {"method": "prfba_exclude", "space": "ec"}
item("front:abs-picrust2-asshipped", [
    V("cond EC F1 prFBA", CS, csm("prfba_exclude", "conditional", "f1", "value"), "0.802"),
    D("matched both-answer PICRUSt2 F1", "as10_picrust2_f1", "0.727"),
    X("no EC prediction share", "1-a", {"a": S(CS, csm("prfba_exclude", "coverage"))}, "11.1", "pct1"),
    X("abstain share", "a/b", {"a": S(GC, ["_meta", "outcome_counts_total", "exclude_0.987", "abstained"]),
                               "b": S(CS, csm("prfba_exclude", "n"))}, "9.6", "pct1"),
])
item("front:fig:picrust2", [
    V("PICRUSt2 precision NSTI<0.01", CS, nsti("<0.01", "picrust2", "conditional", "precision", "value"), "0.735"),
    V("PICRUSt2 precision NSTI>=0.5", CS, nsti(">=0.5", "picrust2", "conditional", "precision", "value"), "0.235"),
    V("prFBA precision NSTI>=0.5", CS, nsti(">=0.5", "prfba_exclude", "conditional", "precision", "value"), "0.606"),
    V("PICRUSt2 precision below family floor", CS, ident("<0.865", "picrust2", "conditional", "precision", "value"), "0.32", "f2"),
    V("family-floor bin n", CS, ident("<0.865", "prfba_exclude", "n"), "547", "int"),
    cc("prFBA EC non-answers", "csv_count", dict(EXC, outcome="abstained"), "7,422"),
    cc("outright abstentions (no reps)", "csv_count", dict(EXC, outcome="abstained", n_reps=["==", 0]), "6,407"),
    X("floor bin share of non-answers", "a/b", {"a": S(CS, ident("<0.865", "prfba_exclude", "n")),
                                                "b": dict(artifact=CC, recipe="csv_count", args=dict(OC, filter=dict(EXC, outcome="abstained")))}, "7", "pct0"),
], manual="145 further / 692 (10.8%) need prfba_selection.json surviving-hit identities (UNION_VERIFICATION.md L336)")
item("val16s:ec_coverage_bins_and_noanswer", [
    V("coverage 0.945-0.987", CS, ident("0.945-0.987", "prfba_exclude", "coverage"), "0.982"),
    V("coverage >=0.995", CS, ident(">=0.995", "prfba_exclude", "coverage"), "0.829"),
    V("answers", CS, csm("prfba_exclude", "coverage"), "89", "pct0"),
    X("no EC prediction", "1-a", {"a": S(CS, csm("prfba_exclude", "coverage"))}, "11.1", "pct1"),
    X("abstain", "a/b", {"a": S(GC, ["_meta", "outcome_counts_total", "exclude_0.987", "abstained"]), "b": S(CS, csm("prfba_exclude", "n"))}, "9.6", "pct1"),
    X("reps without gene families", "(a-c)/b", {"a": dict(artifact=CC, recipe="csv_count", args=dict(OC, filter=dict(EXC, outcome="abstained"))),
                                                "c": S(GC, ["_meta", "outcome_counts_total", "exclude_0.987", "abstained"]),
                                                "b": S(CS, csm("prfba_exclude", "n"))}, "1.5", "pct1"),
])
item("picrust:AS-04", [
    cc("truth_empty amplicons (prfba_exclude ec)", "csv_count", dict(EXC, outcome="truth_empty"), "953"),
    V("scoreable n", CS, csm("picrust2", "n"), "67,038", "int"),
    cc("truth_empty source genomes", "csv_nunique", dict(EXC, outcome="truth_empty"), "161", column="src", derive={"src": ["split", "amplicon", "__", 0]}),
])
item("picrust:AS-05", [
    V("truth_empty (stated in the disclosure)", GC, ["_meta", "outcome_counts_total", "exclude_0.987", "truth_empty"], "953", "int",
      note="after the re-score this becomes 278; L517 must be rewritten as a past-tense note (plan section 5)"),
    V("pred_empty", GC, ["_meta", "outcome_counts_total", "exclude_0.987", "pred_empty"], "1,060", "int"),
], manual="675/818/167 and 0.75-0.86/0.05 points come from CORRECTIONS_2026-11-13.json and genecap_repaired_impact/impact.json (historical record)")
item("picrust:AS-08", [
    X("PICRUSt2 below alignment threshold (scoreable)", "a-b", {"a": S(CS, csm("picrust2", "n")), "b": S(CS, csm("picrust2", "n_called"))}, "26", "int"),
    V("all amplicons without PICRUSt2 prediction", CS, ["_meta", "n_amplicons_without_picrust2_prediction"], "27", "int"),
])
item("picrust:AS-09", [
    V("PICRUSt2 coverage", CS, csm("picrust2", "coverage"), "1.000"),
    V("PICRUSt2 recall", CS, csm("picrust2", "conditional", "recall", "value"), "0.849"),
    V("PICRUSt2 precision", CS, csm("picrust2", "conditional", "precision", "value"), "0.669"),
    V("PICRUSt2 F1", CS, csm("picrust2", "conditional", "f1", "value"), "0.729"),
    V("PICRUSt2 pop F1", CS, csm("picrust2", "population", "f1", "value"), "0.729"),
    V("PICRUSt2 answered n", CS, csm("picrust2", "n_called"), "67,012", "int"),
    V("prFBA excl coverage", CS, csm("prfba_exclude", "coverage"), "88.9", "pct1"),
    V("prFBA excl recall", CS, csm("prfba_exclude", "conditional", "recall", "value"), "0.847"),
    V("prFBA excl precision", CS, csm("prfba_exclude", "conditional", "precision", "value"), "0.787"),
    V("prFBA excl cond F1", CS, csm("prfba_exclude", "conditional", "f1", "value"), "0.802"),
    V("prFBA excl pop F1", CS, csm("prfba_exclude", "population", "f1", "value"), "0.713"),
    V("prFBA incl recall", CS, csm("prfba_include", "conditional", "recall", "value"), "0.999"),
    V("prFBA incl precision", CS, csm("prfba_include", "conditional", "precision", "value"), "0.936"),
    V("prFBA incl F1", CS, csm("prfba_include", "conditional", "f1", "value"), "0.959"),
    X("population delta", "a-b", {"a": S(CS, csm("prfba_exclude", "population", "f1", "value")),
                                  "b": S(CS, csm("picrust2", "population", "f1", "value"))}, "-0.016"),
])
item("picrust:AS-11", [
    X("no EC prediction", "1-a", {"a": S(CS, csm("prfba_exclude", "coverage"))}, "11.1", "pct1"),
    X("abstain", "a/b", {"a": S(GC, ["_meta", "outcome_counts_total", "exclude_0.987", "abstained"]), "b": S(CS, csm("prfba_exclude", "n"))}, "9.6", "pct1"),
    X("reps without gene families", "(a-c)/b", {"a": dict(artifact=CC, recipe="csv_count", args=dict(OC, filter=dict(EXC, outcome="abstained"))),
                                                "c": S(GC, ["_meta", "outcome_counts_total", "exclude_0.987", "abstained"]),
                                                "b": S(CS, csm("prfba_exclude", "n"))}, "1.5", "pct1"),
])
item("picrust:AS-12", [
    V("median dF1 (low)", CS, ["paired_tests_ec_f1", {"each": [f"{r}|prfba_exclude_vs_picrust2|conditional" for r in R8]}, "median_diff"], "+0.072", "+f3", reduce="min"),
    V("median dF1 (high)", CS, ["paired_tests_ec_f1", {"each": [f"{r}|prfba_exclude_vs_picrust2|conditional" for r in R8]}, "median_diff"], "+0.095", "+f3", reduce="max"),
    V("Holm p (max)", CS, ["paired_tests_ec_f1", {"each": [f"{r}|prfba_exclude_vs_picrust2|conditional" for r in R8]}, "p_holm"], "1e-10", "lt", reduce="max",
      note="printed 'p ~ 0'", tex_text="p \\approx 0"),
    X("wins per loss (min)", "min([" + ", ".join(f"w{i}/l{i}" for i in range(8)) + "])",
      {**{f"w{i}": S(CS, pt(r) + ["prfba_wins"]) for i, r in enumerate(R8)}, **{f"l{i}": S(CS, pt(r) + ["picrust2_wins"]) for i, r in enumerate(R8)}}, "3", "f0"),
    X("wins per loss (max)", "max([" + ", ".join(f"w{i}/l{i}" for i in range(8)) + "])",
      {**{f"w{i}": S(CS, pt(r) + ["prfba_wins"]) for i, r in enumerate(R8)}, **{f"l{i}": S(CS, pt(r) + ["picrust2_wins"]) for i, r in enumerate(R8)}}, "7", "f0"),
])
item("picrust:AS-13", [
    V("conditional dF1 vs BV-BRC truth", DIAL, ["_meta", "prfba_minus_picrust2_f1", "conditional|truth=bvbrc"], "+0.046", "+f3"),
    V("conditional dF1 vs PICRUSt2 truth", DIAL, ["_meta", "prfba_minus_picrust2_f1", "conditional|truth=picrust2"], "-0.202"),
    V("amplicons in both databases", DIAL, ["_meta", "n_amplicons"], "19,686", "int"),
])
_pfr = {r: S(CS, ["by_region", "ec", r, "all", "prfba_exclude", "population", "f1", "value"]) for r in R8}
item("picrust:AS-14", [
    *[V(f"prFBA pop F1 {r}", CS, ["by_region", "ec", r, "all", "prfba_exclude", "population", "f1", "value"], p) for r, p in
      (("FullLength16S", "0.768"), ("V1-V2", "0.760"), ("V1-V3", "0.766"), ("V4", "0.668"))],
    *[V(f"PICRUSt2 pop F1 {r}", CS, ["by_region", "ec", r, "all", "picrust2", "population", "f1", "value"], p) for r, p in
      (("FullLength16S", "0.745"), ("V1-V2", "0.731"), ("V1-V3", "0.732"))],
    X("regions where prFBA pop F1 < PICRUSt2 (count)", "sum([" + ", ".join(f"1 if a{i} < b{i} else 0" for i in range(8)) + "])",
      {**{f"a{i}": S(CS, ["by_region", "ec", r, "all", "prfba_exclude", "population", "f1", "value"]) for i, r in enumerate(R8)},
       **{f"b{i}": S(CS, ["by_region", "ec", r, "all", "picrust2", "population", "f1", "value"]) for i, r in enumerate(R8)}}, "5", "int",
      tex_text="five regions"),
])
item("picrust:AS-15", [
    V("rxn cond F1 prFBA", CS, ["by_region", "rxn", "all", "all", "prfba_exclude", "conditional", "f1", "value"], "0.805"),
    V("rxn cond F1 PICRUSt2", CS, ["by_region", "rxn", "all", "all", "picrust2", "conditional", "f1", "value"], "0.721"),
], manual="divergences (0.026 archaeal, 0.014 distant novelty bins, 0.05 / 0.03 in the 9- and 8-amplicon cells) need per-stratum "
          "rxn-ec differences incl. a domain split of comparison_cells; recompute by hand from number_diff + cells")
item("picrust:AS-16", [
    V("PICRUSt2 precision <0.01", CS, nsti("<0.01", "picrust2", "conditional", "precision", "value"), "0.735"),
    V("PICRUSt2 precision >=0.5", CS, nsti(">=0.5", "picrust2", "conditional", "precision", "value"), "0.235"),
    V("PICRUSt2 recall <0.01", CS, nsti("<0.01", "picrust2", "conditional", "recall", "value"), "0.870"),
    V("PICRUSt2 recall >=0.5", CS, nsti(">=0.5", "picrust2", "conditional", "recall", "value"), "0.769"),
    V("prFBA precision <0.01", CS, nsti("<0.01", "prfba_exclude", "conditional", "precision", "value"), "0.822"),
    V("prFBA precision >=0.5", CS, nsti(">=0.5", "prfba_exclude", "conditional", "precision", "value"), "0.606"),
    V("prFBA coverage >=0.5", CS, nsti(">=0.5", "prfba_exclude", "coverage"), "86.0", "pct1"),
    V("prFBA coverage <0.01", CS, nsti("<0.01", "prfba_exclude", "coverage"), "92.1", "pct1"),
    V("prFBA answered at >=0.5", CS, nsti(">=0.5", "prfba_exclude", "n_called"), "603", "int"),
    cc("PICRUSt2 precision on the 603", "csv_mean", {"method": "picrust2", "space": "ec", "nsti_bin": ">=0.5"}, "0.230", "f3",
       column="precision", within={"filter": dict(EXC, outcome="ok", nsti_bin=">=0.5"), "on": "amplicon"}),
    cc("prFBA ECs per answered amplicon", "csv_mean", dict(EXC, outcome="ok", nsti_bin=">=0.5"), "141", "f0", column="n_pred"),
    cc("PICRUSt2 ECs on the same amplicons", "csv_mean", {"method": "picrust2", "space": "ec", "nsti_bin": ">=0.5"}, "412", "f0",
       column="n_pred", within={"filter": dict(EXC, outcome="ok", nsti_bin=">=0.5"), "on": "amplicon"}),
    cc("truth ECs on the same amplicons", "csv_mean", dict(EXC, outcome="ok", nsti_bin=">=0.5"), "123", "f0", column="n_truth"),
])
item("picrust:AS-17", [
    cc("high-NSTI answered with a >=0.995 non-self in the record", "csv_count",
       dict(EXC, outcome="ok", nsti_bin=">=0.5", best_nonself_identity=[">=", 0.995]), "248"),
    V("high-NSTI answered", CS, nsti(">=0.5", "prfba_exclude", "n_called"), "603", "int"),
], manual="'all 603 rest on references below 0.987 (median 0.956)' needs prfba_selection.json best_identity (git-ignored, 19 MB)")
item("picrust:AS-18", [
    V("prFBA coverage below floor", CS, ident("<0.865", "prfba_exclude", "coverage"), "0", "int"),
    V("bin n", CS, ident("<0.865", "prfba_exclude", "n"), "547", "int"),
    V("PICRUSt2 precision in bin", CS, ident("<0.865", "picrust2", "conditional", "precision", "value"), "0.316"),
    cc("EC non-answers", "csv_count", dict(EXC, outcome="abstained"), "7,422"),
    cc("outright abstentions", "csv_count", dict(EXC, outcome="abstained", n_reps=["==", 0]), "6,407"),
], manual="145 / 692 (10.8%) need prfba_selection.json (UNION_VERIFICATION.md L336)")
item("picrust:AS-20", [
    *[V(f"in-ref PICRUSt2 {k}", CS, strat("True", "picrust2", *p), v) for k, p, v in (
        ("coverage", ["coverage"], "1.000"), ("recall", ["conditional", "recall", "value"], "0.868"),
        ("precision", ["conditional", "precision", "value"], "0.735"), ("F1", ["conditional", "f1", "value"], "0.792"),
        ("pop F1", ["population", "f1", "value"], "0.792"))],
    *[V(f"in-ref prFBA {k}", CS, strat("True", "prfba_exclude", *p), v) for k, p, v in (
        ("coverage", ["coverage"], "0.932"), ("recall", ["conditional", "recall", "value"], "0.871"),
        ("precision", ["conditional", "precision", "value"], "0.828"), ("F1", ["conditional", "f1", "value"], "0.838"),
        ("pop F1", ["population", "f1", "value"], "0.781"))],
    *[V(f"absent PICRUSt2 {k}", CS, strat("False", "picrust2", *p), v) for k, p, v in (
        ("coverage", ["coverage"], "0.999"), ("recall", ["conditional", "recall", "value"], "0.840"),
        ("precision", ["conditional", "precision", "value"], "0.642"), ("F1", ["conditional", "f1", "value"], "0.703"),
        ("pop F1", ["population", "f1", "value"], "0.703"))],
    *[V(f"absent prFBA {k}", CS, strat("False", "prfba_exclude", *p), v) for k, p, v in (
        ("coverage", ["coverage"], "0.872"), ("recall", ["conditional", "recall", "value"], "0.837"),
        ("precision", ["conditional", "precision", "value"], "0.769"), ("F1", ["conditional", "f1", "value"], "0.786"),
        ("pop F1", ["population", "f1", "value"], "0.685"))],
    V("in-ref n", CS, strat("True", "picrust2", "n"), "19,686", "int"),
    V("absent n", CS, strat("False", "picrust2", "n"), "47,352", "int"),
    cc("median NSTI of in-reference placements", "csv_median", {"method": "picrust2", "space": "ec", "src_in_picrust2_ref": 1}, "0.010", "f3",
       column="nsti"),
])
item("picrust:AS-21", [
    V("include-self EC Jaccard", CS, csm("prfba_include", "conditional", "jaccard", "value"), "0.94", "f2"),
    V("agreement Jaccard (low)", CS, ["agreement_prfba_include_vs_picrust2_jaccard", {"each": [f"ec|{r}" for r in R8]}, "mean"], "0.56", "f2", reduce="min"),
    V("agreement Jaccard (high)", CS, ["agreement_prfba_include_vs_picrust2_jaccard", {"each": [f"ec|{r}" for r in R8]}, "mean"], "0.61", "f2", reduce="max"),
    V("PICRUSt2 Jaccard by region (low)", CS, ["by_region", "ec", {"each": R8}, "all", "picrust2", "conditional", "jaccard", "value"], "0.58", "f2", reduce="min"),
    V("PICRUSt2 Jaccard by region (high)", CS, ["by_region", "ec", {"each": R8}, "all", "picrust2", "conditional", "jaccard", "value"], "0.61", "f2", reduce="max"),
    V("ceiling Jaccard", CEIL, ["ceiling", "shared", "jaccard", "mean"], "0.688"),
])
item("picrust:AS-23", [], manual="median selection size 1 from prfba_selection.json (git-ignored, 19 MB): exclude_0.987.n_selected median; "
     "add a JSON-stream recipe if it is ever in doubt")
item("picrust:AS-24", [], manual="figure F13 (W4 render_comparison_d100.py computes its suptitle n and panel C range); numbers in AS-12/16/18")
item("picrust:XR-01", [
    V("coverage 0.945-0.987", CS, ident("0.945-0.987", "prfba_exclude", "coverage"), "0.982"),
    V("coverage >=0.995", CS, ident(">=0.995", "prfba_exclude", "coverage"), "0.829"),
])
item("panels:disc.picrust2_failure_modes", [
    X("PICRUSt2 not answering (scoreable)", "a-b", {"a": S(CS, csm("picrust2", "n")), "b": S(CS, csm("picrust2", "n_called"))}, "26", "int"),
    V("scoreable n", CS, csm("picrust2", "n"), "67,038", "int"),
    V("PICRUSt2 precision <0.01", CS, nsti("<0.01", "picrust2", "conditional", "precision", "value"), "0.735"),
    V("PICRUSt2 precision >=0.5", CS, nsti(">=0.5", "picrust2", "conditional", "precision", "value"), "0.235"),
    V("PICRUSt2 recall <0.01", CS, nsti("<0.01", "picrust2", "conditional", "recall", "value"), "0.870"),
    V("PICRUSt2 recall >=0.5", CS, nsti(">=0.5", "picrust2", "conditional", "recall", "value"), "0.769"),
])
item("panels:disc.prfba_answer_rate", [
    V("answers", CS, csm("prfba_exclude", "coverage"), "89", "pct0"),
    X("non-answer", "1-a", {"a": S(CS, csm("prfba_exclude", "coverage"))}, "11", "pct0"),
    X("abstain", "a/b", {"a": S(GC, ["_meta", "outcome_counts_total", "exclude_0.987", "abstained"]), "b": S(CS, csm("prfba_exclude", "n"))}, "9.6", "pct1"),
    X("reps without gene families", "(a-c)/b", {"a": dict(artifact=CC, recipe="csv_count", args=dict(OC, filter=dict(EXC, outcome="abstained"))),
                                                "c": S(GC, ["_meta", "outcome_counts_total", "exclude_0.987", "abstained"]),
                                                "b": S(CS, csm("prfba_exclude", "n"))}, "1.5", "pct1"),
    V("coverage below family floor", CS, ident("<0.865", "prfba_exclude", "coverage"), "0", "int"),
])

# ================================================================ D1 (unchanged under Path A)
item("front:abs-sharedref-size", [
    V("median reference genomes", SG, ["meta", "n_reference_genomes_median"], "7,175", "int"),
    V("min reference genomes per fold", SG, ["accounting", "per_fold", {"each": "*"}, "n_reference_genomes"], "6,788", "int", reduce="min"),
])
item("picrust:SR-01", [
    V("median reference genomes", SG, ["meta", "n_reference_genomes_median"], "7,175", "int"),
    V("min reference genomes per fold", SG, ["accounting", "per_fold", {"each": "*"}, "n_reference_genomes"], "6,788", "int", reduce="min"),
    V("eligible genomes", SG, ["meta", "n_eligible"], "8,872", "int"),
    V("scored cells", SG, ["meta", "n_cells_scored"], "63,871", "int"),
    V("holdout groups", SG, con("prfba_nestedcv__vs__picrust2", "bootstrap_n_groups"), "2,376", "int"),
    V("effective n", SG, con("prfba_nestedcv__vs__picrust2", "n_eff"), "6,716", "int"),
    V("ICC (group)", SG, con("prfba_nestedcv__vs__picrust2", "icc_group"), "0.33", "f2"),
])
item("picrust:FR-01", [
    V("median reference genomes (family)", SF, ["meta", "n_reference_genomes_median"], "8,468", "int"),
    V("scored cells (family)", SF, ["meta", "n_cells_scored"], "14,692", "int"),
    V("holdout groups (family)", SF, con("prfba_nestedcv__vs__picrust2", "bootstrap_n_groups"), "171", "int"),
    V("effective n (family)", SF, con("prfba_nestedcv__vs__picrust2", "n_eff"), "598", "int"),
], manual="2,059 held-out genomes from data_family/fold_plan_family.json")
item("picrust:FR-04", [], manual="1,653 of 8,872 unclassified-family placeholders: count over panel_eligible.json genomes[*].family (unchanged under Path A)")
item("panels:lim.database_sizes", [
    V("shared reference median", SG, ["meta", "n_reference_genomes_median"], "7,175", "int"),
], manual="~315,000 BV-BRC genomes / 459,301 unique 16S / 27,870 PICRUSt2 genomes are database facts (no re-score)")
item("panels:lim.holdout_power", [
    V("family holdout groups", SF, con("prfba_nestedcv__vs__picrust2", "bootstrap_n_groups"), "171", "int"),
    V("genus holdout groups", SG, con("prfba_nestedcv__vs__picrust2", "bootstrap_n_groups"), "2,376", "int"),
    V("family n_eff", SF, con("prfba_nestedcv__vs__picrust2", "n_eff"), "598", "int"),
    V("genus n_eff", SG, con("prfba_nestedcv__vs__picrust2", "n_eff"), "6,716", "int"),
    V("eligible", SG, ["meta", "n_eligible"], "8,872", "int"),
], manual="18.6% (1,653) placeholders from panel_eligible.json")

# ================================================================ R3 shared reference (genus)
ARMS_POP = [("picrust2", "0.779"), ("prfba_nestedcv", "0.687"), ("prfba_exact_tie", "0.685"),
            ("prfba_nestedcv_forced", "0.747"), ("prfba_exact_tie_forced", "0.746")]


def ci(c, i, sp="ec"):
    return con(c, "ci95", i, space=sp)


item("front:abs-controlled-genus", [
    *[V(f"{a} pop F1", SG, arm(a, *POP), p) for a, p in ARMS_POP],
    V("delta nestedcv", SG, con("prfba_nestedcv__vs__picrust2", "delta_f1_population"), "-0.092"),
    V("  CI lo", SG, ci("prfba_nestedcv__vs__picrust2", 0), "-0.107"),
    V("  CI hi", SG, ci("prfba_nestedcv__vs__picrust2", 1), "-0.081"),
    V("delta exact-tie", SG, con("prfba_exact_tie__vs__picrust2", "delta_f1_population"), "-0.094"),
    V("  CI lo", SG, ci("prfba_exact_tie__vs__picrust2", 0), "-0.109"),
    V("  CI hi", SG, ci("prfba_exact_tie__vs__picrust2", 1), "-0.083"),
    V("delta forced", SG, con("prfba_nestedcv_forced__vs__picrust2", "delta_f1_population"), "-0.031"),
    V("delta exact-tie forced", SG, con("prfba_exact_tie_forced__vs__picrust2", "delta_f1_population"), "-0.033"),
])
item("front:abs-dialect-density", [
    V("dialect conditional dF1 (BV-BRC truth)", DIAL, ["_meta", "prfba_minus_picrust2_f1", "conditional|truth=bvbrc"], "+0.046", "+f3"),
    V("dialect conditional dF1 (PICRUSt2 truth)", DIAL, ["_meta", "prfba_minus_picrust2_f1", "conditional|truth=picrust2"], "-0.202"),
    V("amplicons", DIAL, ["_meta", "n_amplicons"], "19,686", "int"),
    V("density slope prFBA", SG, arm("prfba_nestedcv", "density_slope_f1_per_log10_tips"), "+0.107", "+f3"),
    V("density slope exact-tie", SG, arm("prfba_exact_tie", "density_slope_f1_per_log10_tips"), "+0.106", "+f3"),
    V("density slope PICRUSt2", SG, arm("picrust2", "density_slope_f1_per_log10_tips"), "+0.030", "+f3"),
])
item("front:abs-region-penalty", [
    V("as-shipped FL pop F1", CS, ["by_region", "ec", "FullLength16S", "all", "prfba_exclude", "population", "f1", "value"], "0.768"),
    V("as-shipped V4 pop F1", CS, ["by_region", "ec", "V4", "all", "prfba_exclude", "population", "f1", "value"], "0.668"),
    V("controlled FL pop F1", SG, arm("prfba_nestedcv", "by_region", "FullLength16S", *POP), "0.692"),
    V("controlled V4 pop F1", SG, arm("prfba_nestedcv", "by_region", "V4", *POP), "0.710"),
    V("exact-tie FL pop F1", SG, arm("prfba_exact_tie", "by_region", "FullLength16S", *POP), "0.690"),
    V("exact-tie V4 pop F1", SG, arm("prfba_exact_tie", "by_region", "V4", *POP), "0.708"),
    V("controlled FL coverage", SG, arm("prfba_nestedcv", "by_region", "FullLength16S", "coverage"), "0.890"),
    V("controlled V4 coverage", SG, arm("prfba_nestedcv", "by_region", "V4", "coverage"), "0.933"),
])
SRTAB = [("picrust2", "1.000", "0.779", "0.796", "0.92"), ("picrust2_nsti_routed", "1.000", "0.775", "0.793", "0.90"),
         ("prfba_nestedcv", "0.893", "0.687", "0.769", "0.28"), ("prfba_nestedcv_forced", "0.9995", "0.747", "0.769", "0.70"),
         ("prfba_exact_tie", "0.893", "0.685", "0.767", "0.27"), ("prfba_exact_tie_forced", "0.9995", "0.746", "0.767", "0.69"),
         ("base_nn1", "1.000", "0.745", "0.766", "0.69"), ("base_random_abstain", "0.892", "0.666", "0.767", "0.14"),
         ("base_family_donor", "0.679", "0.516", "0.765", "-0.90"), ("base_constant", "1.000", "0.646", "0.665", "0.00"),
         ("base_clade_max", "1.000", "0.790", "0.805", "1.00")]


def srtab_values():
    out = []
    for a, cov, pop, cond, sk in SRTAB:
        out += [V(f"{a} coverage", SG, arm(a, "coverage"), cov, "f4" if len(cov) == 6 else "f3"),
                V(f"{a} pop F1", SG, arm(a, *POP), pop), V(f"{a} cond (matched) F1", SG, arm(a, *CM), cond),
                V(f"{a} skill", SG, arm(a, "skill_score"), sk, "f2")]
    out += [V("matched cells", SG, arm("prfba_nestedcv", "conditional_matched", "n_answered"), "57,021", "int"),
            V("random-abstain matched cells", SG, arm("base_random_abstain", "conditional_matched", "n_answered"), "50,946", "int"),
            V("family-donor matched cells", SG, arm("base_family_donor", "conditional_matched", "n_answered"), "41,257", "int"),
            V("scored cells", SG, ["meta", "n_cells_scored"], "63,871", "int"),
            V("reference genomes", SG, ["meta", "n_reference_genomes_median"], "7,175", "int")]
    return out


item("front:tab:sharedref", srtab_values())
item("methods:sharedref_selection_sizes", [
    D("fold 0 single-genome answered", "sr_fold0_single", "9,846", "int"),
    D("fold 0 answered", "sr_fold0_answered", "11,751", "int"),
    D("pooled single-genome share", "sr_single_share", "86.8", "pct1"),
    D("pooled answered", "sr_answered", "50,842", "int"),
    D("mean selection size", "sr_mean_size", "1.199"),
])
item("methods:sharedref_ngenomes_vs_identity", [], manual="P(n_genomes>=2) by novelty bin (0.100 -> 0.180) from SR cells.csv.gz; W6 controlled tie-set block")
item("val16s:controlled_forced_recovery", [
    V("forced coverage", SG, arm("prfba_nestedcv_forced", "coverage"), "99.95", "pct2"),
    X("recovered by forcing", "a-b", {"a": S(SG, arm("prfba_nestedcv_forced", *POP)), "b": S(SG, arm("prfba_nestedcv", *POP))}, "+0.061", "+f3"),
    V("deficit", SG, con("prfba_nestedcv__vs__picrust2", "delta_f1_population"), "0.092", "absf3"),
    V("exact-tie deficit", SG, con("prfba_exact_tie__vs__picrust2", "delta_f1_population"), "0.094", "absf3"),
    V("residual forced", SG, con("prfba_nestedcv_forced__vs__picrust2", "delta_f1_population"), "-0.031"),
    V("  CI lo", SG, ci("prfba_nestedcv_forced__vs__picrust2", 0), "-0.037"),
    V("  CI hi", SG, ci("prfba_nestedcv_forced__vs__picrust2", 1), "-0.028"),
    V("residual exact-tie forced", SG, con("prfba_exact_tie_forced__vs__picrust2", "delta_f1_population"), "-0.033"),
    V("  CI lo", SG, ci("prfba_exact_tie_forced__vs__picrust2", 0), "-0.039"),
    V("  CI hi", SG, ci("prfba_exact_tie_forced__vs__picrust2", 1), "-0.029"),
])
item("val16s:controlled_selection_sizes", [
    D("fold 0 single", "sr_fold0_single", "9,846", "int"),
    D("fold 0 answered", "sr_fold0_answered", "11,751", "int"),
    D("fold 0 two", "sr_fold0_two", "1,354", "int"),
    D("fold 0 three+", "sr_fold0_three_plus", "551", "int"),
    D("pooled single share", "sr_single_share", "86.8", "pct1"),
    D("pooled single", "sr_single_n", "44,132", "int"),
    D("pooled answered", "sr_answered", "50,842", "int"),
    D("mean size", "sr_mean_size", "1.199"),
    D("legacy single share", "sr_legacy_single_share", "61.7", "pct1"),
    D("cross-phylum tie sets", "sr_tieset_span_phyla_unique", "6.4", "pct1"),
    D("same-genus tie sets", "sr_tieset_share_genus_unique", "35.8", "pct1"),
], manual="4.8% of multi-genome selections (L1236) and the monotone P(n>=2) by identity are part of the W6 tie-set block")
item("picrust:SR-02", [
    V("nested-CV winners legacy (count)", SG, ["meta", "notes", "nested_cv_winners", {"each": "*"}], "5", "int",
      reduce="len", note="all five are 'legacy' at the tag; the count of folds is what is checked"),
    V("nestedcv pop F1", SG, arm("prfba_nestedcv", *POP), "0.687"),
    V("exact-tie pop F1", SG, arm("prfba_exact_tie", *POP), "0.685"),
    X("exact-tie minus nestedcv", "b-a", {"a": S(SG, arm("prfba_nestedcv", *POP)), "b": S(SG, arm("prfba_exact_tie", *POP))}, "0.002", "absf3"),
    X("inner margin legacy - best other (min)", "min([" + ", ".join(
        f"l{i}-max(e{i},g{i},t{i})" for i in range(5)) + "])",
      {**{f"l{i}": S(SG, ["nested_cv", f"fold_{i}", "inner_f1", "legacy"]) for i in range(5)},
       **{f"e{i}": S(SG, ["nested_cv", f"fold_{i}", "inner_f1", "exact_tie"]) for i in range(5)},
       **{f"g{i}": S(SG, ["nested_cv", f"fold_{i}", "inner_f1", "graduated"]) for i in range(5)},
       **{f"t{i}": S(SG, ["nested_cv", f"fold_{i}", "inner_f1", "exact_tie_guards"]) for i in range(5)}}, "0.001"),
], manual="'10/10 at both ranks' includes data_family nested_cv; margins 0.001-0.003 are inner-F1 gaps (upper end by hand)")
item("picrust:SR-03", [
    *[V(f"{a} pop F1", SG, arm(a, *POP), p) for a, p in ARMS_POP],
    V("prFBA coverage", SG, arm("prfba_nestedcv", "coverage"), "0.893"),
    V("forced coverage", SG, arm("prfba_nestedcv_forced", "coverage"), "0.9995", "f4"),
    *[V(f"{c} delta", SG, con(f"{c}__vs__picrust2", "delta_f1_population"), d) for c, d in (
        ("prfba_nestedcv", "-0.092"), ("prfba_exact_tie", "-0.094"), ("prfba_nestedcv_forced", "-0.031"), ("prfba_exact_tie_forced", "-0.033"))],
    *[V(f"{c} CI {i}", SG, ci(f"{c}__vs__picrust2", i), d) for c, i, d in (
        ("prfba_nestedcv", 0, "-0.107"), ("prfba_nestedcv", 1, "-0.081"), ("prfba_exact_tie", 0, "-0.109"), ("prfba_exact_tie", 1, "-0.083"),
        ("prfba_nestedcv_forced", 0, "-0.037"), ("prfba_nestedcv_forced", 1, "-0.028"),
        ("prfba_exact_tie_forced", 0, "-0.039"), ("prfba_exact_tie_forced", 1, "-0.029"))],
    V("SESOI verdict", SG, con("prfba_nestedcv__vs__picrust2", "sesoi_verdict"), "exceeded", "str"),
])
_sr4 = []
for c, lo, hi in (("prfba_nestedcv", "-0.059", "-0.182"), ("prfba_nestedcv_forced", "-0.029", "-0.034"),
                  ("prfba_exact_tie", "-0.060", "-0.184"), ("prfba_exact_tie_forced", "-0.031", "-0.036")):
    sel = con(f"{c}__vs__picrust2", "by_region", {"each": R8}, "delta_f1_population")
    _sr4 += [V(f"{c} by region (smallest deficit)", SG, sel, lo, reduce="max"), V(f"{c} by region (largest deficit)", SG, sel, hi, reduce="min")]
for c, b, a in (("prfba_nestedcv", "-0.095", "-0.052"), ("prfba_nestedcv_forced", "-0.032", "-0.019"),
                ("prfba_exact_tie", "-0.097", "-0.053"), ("prfba_exact_tie_forced", "-0.034", "-0.020")):
    _sr4 += [V(f"{c} Bacteria", SG, con(f"{c}__vs__picrust2", "by_domain", "Bacteria", "delta_f1_population"), b),
             V(f"{c} Archaea", SG, con(f"{c}__vs__picrust2", "by_domain", "Archaea", "delta_f1_population"), a)]
_sr4.append(V("max BH q over regions (nestedcv)", SG, con("prfba_nestedcv__vs__picrust2", "by_region", {"each": R8}, "q_bh"), "0.01", "le", reduce="max"))
item("picrust:SR-04", _sr4)
item("picrust:SR-05", [
    V("PICRUSt2 AURC", SG, arm("picrust2", "aurc"), "0.812"),
    V("prFBA AURC", SG, arm("prfba_nestedcv", "aurc"), "0.788"),
    V("exact-tie AURC", SG, arm("prfba_exact_tie", "aurc"), "0.787"),
    V("PICRUSt2 cross-knob AURC", SG, arm("picrust2", "aurc_cross_knob"), "0.813"),
    V("prFBA cross-knob AURC", SG, arm("prfba_nestedcv", "aurc_cross_knob"), "0.788"),
    V("exact-tie cross-knob AURC", SG, arm("prfba_exact_tie", "aurc_cross_knob"), "0.786"),
    X("family swap (prFBA aurc - cross)", "a-b", {"a": S(SF, arm("prfba_nestedcv", "aurc")), "b": S(SF, arm("prfba_nestedcv", "aurc_cross_knob"))}, "0.007",
      note="'~0.007 in PICRUSt2's favour' at leave-family-out; the exact pairing used for the text is not recorded"),
])
item("picrust:SR-06", [
    V("FL pop F1", SG, arm("prfba_nestedcv", "by_region", "FullLength16S", *POP), "0.692"),
    V("V4 pop F1", SG, arm("prfba_nestedcv", "by_region", "V4", *POP), "0.710"),
    V("exact-tie FL", SG, arm("prfba_exact_tie", "by_region", "FullLength16S", *POP), "0.690"),
    V("exact-tie V4", SG, arm("prfba_exact_tie", "by_region", "V4", *POP), "0.708"),
    V("FL coverage", SG, arm("prfba_nestedcv", "by_region", "FullLength16S", "coverage"), "0.890"),
    V("V4 coverage", SG, arm("prfba_nestedcv", "by_region", "V4", "coverage"), "0.933"),
    V("exact-tie FL coverage", SG, arm("prfba_exact_tie", "by_region", "FullLength16S", "coverage"), "0.890"),
    V("exact-tie V4 coverage", SG, arm("prfba_exact_tie", "by_region", "V4", "coverage"), "0.933"),
])
item("picrust:SR-07", [
    V("PICRUSt2 rxn pop F1", SG, arm("picrust2", *POP, space="rxn"), "0.780"),
    V("prFBA rxn pop F1", SG, arm("prfba_nestedcv", *POP, space="rxn"), "0.691"),
    V("delta", SG, con("prfba_nestedcv__vs__picrust2", "delta_f1_population", space="rxn"), "-0.089"),
    V("  CI lo", SG, ci("prfba_nestedcv__vs__picrust2", 0, "rxn"), "-0.105"),
    V("  CI hi", SG, ci("prfba_nestedcv__vs__picrust2", 1, "rxn"), "-0.078"),
    V("exact-tie delta", SG, con("prfba_exact_tie__vs__picrust2", "delta_f1_population", space="rxn"), "-0.092"),
    V("  CI lo", SG, ci("prfba_exact_tie__vs__picrust2", 0, "rxn"), "-0.107"),
    V("  CI hi", SG, ci("prfba_exact_tie__vs__picrust2", 1, "rxn"), "-0.081"),
    V("forced delta", SG, con("prfba_nestedcv_forced__vs__picrust2", "delta_f1_population", space="rxn"), "-0.028"),
    V("exact-tie forced delta", SG, con("prfba_exact_tie_forced__vs__picrust2", "delta_f1_population", space="rxn"), "-0.031"),
    V("max BH q over regions", SG, con("prfba_nestedcv__vs__picrust2", "by_region", {"each": R8}, "q_bh", space="rxn"), "0.002", reduce="max"),
    V("PICRUSt2 rxn AURC", SG, arm("picrust2", "aurc", space="rxn"), "0.815"),
    V("prFBA rxn AURC", SG, arm("prfba_nestedcv", "aurc", space="rxn"), "0.794"),
    V("exact-tie rxn AURC", SG, arm("prfba_exact_tie", "aurc", space="rxn"), "0.792"),
])
item("picrust:SR-08", [
    V("constant pop F1", SG, arm("base_constant", *POP), "0.646"),
    V("clade-max pop F1", SG, arm("base_clade_max", *POP), "0.790"),
    V("PICRUSt2 skill", SG, arm("picrust2", "skill_score"), "0.92", "f2"),
    V("forced skill", SG, arm("prfba_nestedcv_forced", "skill_score"), "0.70", "f2"),
    V("exact-tie forced skill", SG, arm("prfba_exact_tie_forced", "skill_score"), "0.69", "f2"),
    V("nn1 skill", SG, arm("base_nn1", "skill_score"), "0.69", "f2"),
    V("forced pop F1", SG, arm("prfba_nestedcv_forced", *POP), "0.747"),
    V("exact-tie forced pop F1", SG, arm("prfba_exact_tie_forced", *POP), "0.746"),
    V("nn1 pop F1", SG, arm("base_nn1", *POP), "0.745"),
    V("nestedcv vs nn1", SG, con("prfba_nestedcv__vs__base_nn1", "delta_f1_population"), "-0.058"),
    V("  CI lo", SG, ci("prfba_nestedcv__vs__base_nn1", 0), "-0.069"),
    V("  CI hi", SG, ci("prfba_nestedcv__vs__base_nn1", 1), "-0.051"),
    V("  q", SG, con("prfba_nestedcv__vs__base_nn1", "q_bh"), "0.002"),
    X("exact-tie vs nn1", "a-b", {"a": S(SG, arm("prfba_exact_tie", *POP)), "b": S(SG, arm("base_nn1", *POP))}, "-0.060",
      note="from rounded 4-dp arm values; the text's -0.060135 is a cells.csv.gz recomputation"),
    V("PICRUSt2 vs nn1", SG, con("picrust2__vs__base_nn1", "delta_f1_population"), "+0.034", "+f3"),
    V("  CI lo", SG, ci("picrust2__vs__base_nn1", 0), "0.030"),
    V("  CI hi", SG, ci("picrust2__vs__base_nn1", 1), "0.040"),
])
item("picrust:SR-09", [
    V("vs random abstain", SG, con("prfba_nestedcv__vs__base_random_abstain", "delta_f1_population"), "+0.021", "+f3"),
    V("  CI lo", SG, ci("prfba_nestedcv__vs__base_random_abstain", 0), "0.014"),
    V("  CI hi", SG, ci("prfba_nestedcv__vs__base_random_abstain", 1), "0.030"),
    V("  q", SG, con("prfba_nestedcv__vs__base_random_abstain", "q_bh"), "0.002"),
    X("declined share", "1-a", {"a": S(SG, arm("prfba_nestedcv", "coverage"))}, "10.7", "pct1"),
    V("random-abstain matched cells", SG, arm("base_random_abstain", "conditional_matched", "n_answered"), "50,946", "int"),
    X("exact-tie vs its own random mask (approx.)", "a-b", {"a": S(SG, arm("prfba_exact_tie", *POP)), "b": S(SG, arm("base_random_abstain", *POP))}, "+0.019", "+f3",
      note="the text's like-for-like value is a cells.csv.gz recomputation with the exact-tie mask; this proxy uses the nested-CV mask arm"),
    X("answering everywhere gains", "a-b", {"a": S(SG, arm("prfba_nestedcv_forced", *POP)), "b": S(SG, arm("prfba_nestedcv", *POP))}, "+0.061", "+f3"),
    V("forced pop F1 refused band", SG, arm("prfba_nestedcv_forced", "by_novelty_bin", "<0.865", *POP), "0.567"),
    V("exact-tie forced refused band", SG, arm("prfba_exact_tie_forced", "by_novelty_bin", "<0.865", *POP), "0.567"),
    V("PICRUSt2 refused band", SG, arm("picrust2", "by_novelty_bin", "<0.865", *POP), "0.634"),
], manual="legacy-masked forced 0.668 and the like-for-like +0.019 are cells.csv.gz recomputations (UNION_THROUGHOUT_CHANGES.md L2530, L2748)")
item("picrust:SR-10", [
    V("PICRUSt2 precision 0.945-0.987", SG, arm("picrust2", "by_novelty_bin", "0.945-0.987", "population", "precision", "value"), "0.862"),
    V("PICRUSt2 precision <0.865", SG, arm("picrust2", "by_novelty_bin", "<0.865", "population", "precision", "value"), "0.680"),
    V("prFBA precision 0.945-0.987", SG, arm("prfba_nestedcv", "by_novelty_bin", "0.945-0.987", "population", "precision", "value"), "0.783"),
    V("prFBA precision 0.865-0.945", SG, arm("prfba_nestedcv", "by_novelty_bin", "0.865-0.945", "population", "precision", "value"), "0.740"),
    V("exact-tie precision 0.945-0.987", SG, arm("prfba_exact_tie", "by_novelty_bin", "0.945-0.987", "population", "precision", "value"), "0.803"),
    V("exact-tie precision 0.865-0.945", SG, arm("prfba_exact_tie", "by_novelty_bin", "0.865-0.945", "population", "precision", "value"), "0.751"),
    V("prFBA coverage <0.865", SG, arm("prfba_nestedcv", "by_novelty_bin", "<0.865", "coverage"), "0", "int"),
    V("PICRUSt2 pop F1 <0.865", SG, arm("picrust2", "by_novelty_bin", "<0.865", *POP), "0.634"),
])
item("picrust:SR-11", srtab_values() + [
    X("matched conditional gap", "a-b", {"a": S(SG, arm("prfba_nestedcv", *CM)), "b": S(SG, arm("picrust2", *CM))}, "-0.027"),
    X("exact-tie matched gap", "a-b", {"a": S(SG, arm("prfba_exact_tie", *CM)), "b": S(SG, arm("picrust2", *CM))}, "-0.029"),
])
item("picrust:SR-12", [
    V("slope prFBA", SG, arm("prfba_nestedcv", "density_slope_f1_per_log10_tips"), "+0.107", "+f3"),
    V("slope exact-tie", SG, arm("prfba_exact_tie", "density_slope_f1_per_log10_tips"), "+0.106", "+f3"),
    V("slope PICRUSt2", SG, arm("picrust2", "density_slope_f1_per_log10_tips"), "+0.030", "+f3"),
    X("slope ratio", "a/b", {"a": S(SG, arm("prfba_nestedcv", "density_slope_f1_per_log10_tips")),
                             "b": S(SG, arm("picrust2", "density_slope_f1_per_log10_tips"))}, "3.5", "f1"),
    X("quartering cost prFBA", "a-b", {"a": S(SG, arm("prfba_nestedcv", "by_density", "1.0", "f1")), "b": S(SG, arm("prfba_nestedcv", "by_density", "0.25", "f1"))}, "0.064"),
    X("quartering cost exact-tie", "a-b", {"a": S(SG, arm("prfba_exact_tie", "by_density", "1.0", "f1")), "b": S(SG, arm("prfba_exact_tie", "by_density", "0.25", "f1"))}, "0.064"),
    X("quartering cost PICRUSt2", "a-b", {"a": S(SG, arm("picrust2", "by_density", "1.0", "f1")), "b": S(SG, arm("picrust2", "by_density", "0.25", "f1"))}, "0.018"),
])
item("picrust:XR-03", [
    D("pooled single share", "sr_single_share", "86.8", "pct1"),
    D("pooled answered", "sr_answered", "50,842", "int"),
    D("pooled single", "sr_single_n", "44,132", "int"),
    D("fold 0 single", "sr_fold0_single", "9,846", "int"),
    D("fold 0 answered", "sr_fold0_answered", "11,751", "int"),
    D("fold 0 two", "sr_fold0_two", "1,354", "int"),
    D("fold 0 three+", "sr_fold0_three_plus", "551", "int"),
    D("mean size", "sr_mean_size", "1.199"),
    D("legacy single share", "sr_legacy_single_share", "61.7", "pct1"),
    D("cross-phylum tie sets", "sr_tieset_span_phyla_unique", "6.4", "pct1"),
    D("same-genus tie sets", "sr_tieset_share_genus_unique", "35.8", "pct1"),
])
item("panels:disc.selection_size_controlled", [
    D("pooled single share", "sr_single_share", "86.8", "pct1"),
    D("pooled single", "sr_single_n", "44,132", "int"),
    D("pooled answered", "sr_answered", "50,842", "int"),
    D("mean size", "sr_mean_size", "1.199"),
    D("legacy single share", "sr_legacy_single_share", "61.7", "pct1"),
])
item("panels:disc.margin_reversal", [
    V("as-shipped cond F1 prFBA", CS, csm("prfba_exclude", "conditional", "f1", "value"), "0.802"),
    D("as-shipped matched PICRUSt2", "as10_picrust2_f1", "0.727"),
    D("as-shipped matched margin", "as10_delta", "+0.075", "+f3"),
    V("controlled matched prFBA", SG, arm("prfba_nestedcv", *CM), "0.769"),
    V("controlled matched PICRUSt2", SG, arm("picrust2", *CM), "0.796"),
    V("controlled matched exact-tie", SG, arm("prfba_exact_tie", *CM), "0.767"),
    X("matched reversal", "a-b", {"a": S(SG, arm("prfba_nestedcv", *CM)), "b": S(SG, arm("picrust2", *CM))}, "-0.027"),
    X("exact-tie matched reversal", "a-b", {"a": S(SG, arm("prfba_exact_tie", *CM)), "b": S(SG, arm("picrust2", *CM))}, "-0.029"),
    X("as-shipped population deficit", "a-b", {"a": S(CS, csm("prfba_exclude", "population", "f1", "value")),
                                              "b": S(CS, csm("picrust2", "population", "f1", "value"))}, "-0.016"),
    V("as-shipped prFBA pop", CS, csm("prfba_exclude", "population", "f1", "value"), "0.713"),
    V("as-shipped PICRUSt2 pop", CS, csm("picrust2", "population", "f1", "value"), "0.729"),
    X("depth-20 population deficit", "a-b", {"a": S(CS20, csm("prfba_exclude", "population", "f1", "value")),
                                            "b": S(CS20, csm("picrust2", "population", "f1", "value"))}, "-0.129"),
    V("controlled deficit", SG, con("prfba_nestedcv__vs__picrust2", "delta_f1_population"), "-0.092"),
    V("controlled exact-tie deficit", SG, con("prfba_exact_tie__vs__picrust2", "delta_f1_population"), "-0.094"),
])
item("panels:disc.vocab_breadth", [
    V("ceiling recall", CEIL, ["ceiling", "shared", "recall", "mean"], "90.4", "pct1"),
    V("slope prFBA", SG, arm("prfba_nestedcv", "density_slope_f1_per_log10_tips"), "0.107"),
    V("slope exact-tie", SG, arm("prfba_exact_tie", "density_slope_f1_per_log10_tips"), "0.106"),
    V("slope PICRUSt2", SG, arm("picrust2", "density_slope_f1_per_log10_tips"), "0.030"),
    X("quartering cost prFBA", "a-b", {"a": S(SG, arm("prfba_nestedcv", "by_density", "1.0", "f1")), "b": S(SG, arm("prfba_nestedcv", "by_density", "0.25", "f1"))}, "0.064"),
    X("quartering cost PICRUSt2", "a-b", {"a": S(SG, arm("picrust2", "by_density", "1.0", "f1")), "b": S(SG, arm("picrust2", "by_density", "0.25", "f1"))}, "0.018"),
    V("dialect +0.046", DIAL, ["_meta", "prfba_minus_picrust2_f1", "conditional|truth=bvbrc"], "+0.046", "+f3"),
    V("dialect -0.202", DIAL, ["_meta", "prfba_minus_picrust2_f1", "conditional|truth=picrust2"], "-0.202"),
])

# ================================================================ R3F shared reference (family)
item("front:abs-controlled-family-nn1", [
    V("family delta", SF, con("prfba_nestedcv__vs__picrust2", "delta_f1_population"), "-0.123"),
    V("  CI lo", SF, ci("prfba_nestedcv__vs__picrust2", 0), "-0.150"),
    V("  CI hi", SF, ci("prfba_nestedcv__vs__picrust2", 1), "-0.099"),
    V("family exact-tie delta", SF, con("prfba_exact_tie__vs__picrust2", "delta_f1_population"), "-0.126"),
    V("  CI lo", SF, ci("prfba_exact_tie__vs__picrust2", 0), "-0.152"),
    V("  CI hi", SF, ci("prfba_exact_tie__vs__picrust2", 1), "-0.102"),
    V("genus nestedcv vs nn1", SG, con("prfba_nestedcv__vs__base_nn1", "delta_f1_population"), "-0.058"),
    V("  CI lo", SG, ci("prfba_nestedcv__vs__base_nn1", 0), "-0.069"),
    V("  CI hi", SG, ci("prfba_nestedcv__vs__base_nn1", 1), "-0.051"),
    X("exact-tie vs nn1", "a-b", {"a": S(SG, arm("prfba_exact_tie", *POP)), "b": S(SG, arm("base_nn1", *POP))}, "-0.060"),
    V("forced pop F1", SG, arm("prfba_nestedcv_forced", *POP), "0.747"),
    V("nn1 pop F1", SG, arm("base_nn1", *POP), "0.745"),
    V("exact-tie forced pop F1", SG, arm("prfba_exact_tie_forced", *POP), "0.746"),
])
item("front:abs-controlled-random-mask", [
    V("genus vs random mask", SG, con("prfba_nestedcv__vs__base_random_abstain", "delta_f1_population"), "+0.021", "+f3"),
    V("family vs random mask", SF, con("prfba_nestedcv__vs__base_random_abstain", "delta_f1_population"), "+0.018", "+f3"),
    V("  CI lo", SF, ci("prfba_nestedcv__vs__base_random_abstain", 0), "-0.004"),
    V("  CI hi", SF, ci("prfba_nestedcv__vs__base_random_abstain", 1), "0.040"),
], manual="exact-tie +0.019 / +0.015 are cells.csv.gz recomputations against the exact-tie mask")
item("picrust:FR-02", [
    V("PICRUSt2 pop F1", SF, arm("picrust2", *POP), "0.749"),
    V("prFBA pop F1", SF, arm("prfba_nestedcv", *POP), "0.626"),
    V("prFBA coverage", SF, arm("prfba_nestedcv", "coverage"), "0.856"),
    V("delta", SF, con("prfba_nestedcv__vs__picrust2", "delta_f1_population"), "-0.123"),
    V("  CI lo", SF, ci("prfba_nestedcv__vs__picrust2", 0), "-0.150"),
    V("  CI hi", SF, ci("prfba_nestedcv__vs__picrust2", 1), "-0.099"),
    V("exact-tie pop F1", SF, arm("prfba_exact_tie", *POP), "0.623"),
    V("exact-tie delta", SF, con("prfba_exact_tie__vs__picrust2", "delta_f1_population"), "-0.126"),
    V("  CI lo", SF, ci("prfba_exact_tie__vs__picrust2", 0), "-0.152"),
    V("  CI hi", SF, ci("prfba_exact_tie__vs__picrust2", 1), "-0.102"),
    V("forced pop F1", SF, arm("prfba_nestedcv_forced", *POP), "0.709"),
    V("forced delta", SF, con("prfba_nestedcv_forced__vs__picrust2", "delta_f1_population"), "-0.040"),
    V("  CI lo", SF, ci("prfba_nestedcv_forced__vs__picrust2", 0), "-0.048"),
    V("  CI hi", SF, ci("prfba_nestedcv_forced__vs__picrust2", 1), "-0.031"),
    V("exact-tie forced pop F1", SF, arm("prfba_exact_tie_forced", *POP), "0.707"),
    V("exact-tie forced delta", SF, con("prfba_exact_tie_forced__vs__picrust2", "delta_f1_population"), "-0.042"),
    V("  CI lo", SF, ci("prfba_exact_tie_forced__vs__picrust2", 0), "-0.050"),
    V("  CI hi", SF, ci("prfba_exact_tie_forced__vs__picrust2", 1), "-0.033"),
])
item("picrust:FR-03", [
    V("by region (smallest deficit)", SF, con("prfba_nestedcv__vs__picrust2", "by_region", {"each": R8}, "delta_f1_population"), "-0.076", reduce="max"),
    V("by region (largest deficit)", SF, con("prfba_nestedcv__vs__picrust2", "by_region", {"each": R8}, "delta_f1_population"), "-0.245", reduce="min"),
    V("exact-tie by region (smallest)", SF, con("prfba_exact_tie__vs__picrust2", "by_region", {"each": R8}, "delta_f1_population"), "-0.077", reduce="max"),
    V("exact-tie by region (largest)", SF, con("prfba_exact_tie__vs__picrust2", "by_region", {"each": R8}, "delta_f1_population"), "-0.247", reduce="min"),
    V("max q over regions", SF, con("prfba_nestedcv__vs__picrust2", "by_region", {"each": R8}, "q_bh"), "0.002", reduce="max"),
    V("forced pop F1", SF, arm("prfba_nestedcv_forced", *POP), "0.7095", "f4"),
    V("nn1 pop F1", SF, arm("base_nn1", *POP), "0.7070", "f4"),
    V("exact-tie forced pop F1", SF, arm("prfba_exact_tie_forced", *POP), "0.7072", "f4"),
    V("vs nn1", SF, con("prfba_nestedcv__vs__base_nn1", "delta_f1_population"), "-0.081"),
    V("  CI lo", SF, ci("prfba_nestedcv__vs__base_nn1", 0), "-0.107"),
    V("  CI hi", SF, ci("prfba_nestedcv__vs__base_nn1", 1), "-0.059"),
    V("vs constant", SF, con("prfba_nestedcv__vs__base_constant", "delta_f1_population"), "-0.017"),
    V("  CI lo", SF, ci("prfba_nestedcv__vs__base_constant", 0), "-0.053"),
    V("  CI hi", SF, ci("prfba_nestedcv__vs__base_constant", 1), "0.018"),
    V("skill", SF, arm("prfba_nestedcv", "skill_score"), "-0.14", "f2"),
    V("vs random mask", SF, con("prfba_nestedcv__vs__base_random_abstain", "delta_f1_population"), "+0.018", "+f3"),
    V("vs family donor", SF, con("prfba_nestedcv__vs__base_family_donor", "delta_f1_population"), "+0.107", "+f3"),
    V("  CI lo", SF, ci("prfba_nestedcv__vs__base_family_donor", 0), "0.044"),
    V("  CI hi", SF, ci("prfba_nestedcv__vs__base_family_donor", 1), "0.179"),
    X("family swap", "a-b", {"a": S(SF, arm("prfba_nestedcv", "aurc")), "b": S(SF, arm("prfba_nestedcv", "aurc_cross_knob"))}, "0.007"),
], manual="exact-tie point values (-0.083, -0.020, -0.16, +0.015, +0.104, +0.016 like-for-like) are data_family/cells.csv.gz recomputations")

# ================================================================ R4P phenotype
PC13 = ["contrasts", "prfba_exact_tie__vs__picrust2", "controlled"]
PN13 = ["contrasts", "prfba_exact_tie__vs__base_nn1", "controlled"]
item("front:abs-phenotype-readout", [
    V("ATP-informative genomes", PATP, ["meta", "n_atp_informative"], "755", "int"),
    V("hamming delta", PATP, PC13 + ["delta_atp_hamming"], "0.13", "f2"),
    V("Wilcoxon p", PATP, PC13 + ["wilcoxon_p"], "0.49", "f2"),
], manual="747 V4, EC-F1 -0.028, '2 of 22 media' and '22 of 54 media' derive from per_genome.jsonl / mcnemar_per_medium (PH-01/02/05)")
item("picrust:PH-01", [
    V("paired genomes", PATP, PC13 + ["n_pairs"], "755", "int"),
    V("hamming delta", PATP, PC13 + ["delta_atp_hamming"], "0.13", "f2"),
    V("Wilcoxon p", PATP, PC13 + ["wilcoxon_p"], "0.49", "f2"),
    V("q", PATP, PC13 + ["q_bh"], "0.49", "f2"),
    V("subsample", PATP, ["meta", "n_genomes"], "1,000", "int"),
], manual="747 V4 and EC-F1 -0.028 from phenotype_atp_c13/per_genome.jsonl (no stored key)")
item("picrust:PH-02", [
    V("Ac.NO q", PATP, PC13 + ["mcnemar_per_medium", "Ac.NO", "q_bh"], "0.002", "f3"),
    V("Etho.O2 q", PATP, PC13 + ["mcnemar_per_medium", "Etho.O2", "q_bh"], "0.002", "f3"),
], manual="'two of the 22 testable media' is a count over mcnemar_per_medium favouring PICRUSt2 at q<=0.05")
item("picrust:PH-03", [
    V("Spearman rho", PATP, PC13 + ["ec_f1_delta_vs_phenotype_spearman"], "-0.14", "f2"),
], manual="p = 2e-4 is not stored (per_genome.jsonl)")
item("picrust:PH-04", [
    V("hamming delta vs nn1", PATP, PN13 + ["delta_atp_hamming"], "0.146", "absf3"),
    V("p", PATP, PN13 + ["wilcoxon_p"], "0.017"),
    V("q", PATP, PN13 + ["q_bh"], "0.028"),
    V("p all genomes", PATP, PN13 + ["wilcoxon_p_all_genomes"], "0.28", "f2"),
    V("matched cond F1 exact-tie", SG, arm("prfba_exact_tie", *CM), "0.767"),
    V("matched cond F1 nn1", SG, arm("base_nn1", *CM), "0.766"),
], manual="EC F1 0.782 vs 0.781 (p=0.64) on the 755 genomes from per_genome.jsonl")
item("picrust:PH-05", [], manual="22 of 54 media with ATP signal: derived from phenotype_atp_c13/per_genome.jsonl atp vectors (no stored key)")
item("panels:disc.atp_vs_picrust2", [
    V("hamming delta", PATP, PC13 + ["delta_atp_hamming"], "0.13", "f2"),
    V("Wilcoxon p", PATP, PC13 + ["wilcoxon_p"], "0.49", "f2"),
    V("q", PATP, PC13 + ["q_bh"], "0.49", "f2"),
    V("n", PATP, PC13 + ["n_pairs"], "755", "int"),
    V("Ac.NO q", PATP, PC13 + ["mcnemar_per_medium", "Ac.NO", "q_bh"], "0.002", "f3"),
    V("Etho.O2 q", PATP, PC13 + ["mcnemar_per_medium", "Etho.O2", "q_bh"], "0.002", "f3"),
])
item("panels:disc.atp_vs_nn1", [
    V("hamming delta vs nn1", PATP, PN13 + ["delta_atp_hamming"], "0.146", "absf3"),
    V("p", PATP, PN13 + ["wilcoxon_p"], "0.017"),
    V("q", PATP, PN13 + ["q_bh"], "0.028"),
    V("p all", PATP, PN13 + ["wilcoxon_p_all_genomes"], "0.28", "f2"),
    V("n", PATP, PN13 + ["n_pairs"], "755", "int"),
], manual="EC F1 0.782 vs 0.781 (p=0.64) from per_genome.jsonl")
item("panels:disc.refusal_credit", [
    V("vs random mask", SG, con("prfba_nestedcv__vs__base_random_abstain", "delta_f1_population"), "+0.021", "+f3"),
    V("Spearman rho", PATP, PC13 + ["ec_f1_delta_vs_phenotype_spearman"], "-0.14", "f2"),
], manual="+0.019 exact-tie (cells.csv.gz), p=2e-4 and -0.028 (per_genome.jsonl) not stored; 'two thirds of the deficit' is "
          "(forced - abstaining)/deficit, see picrust:SR-09")
item("panels:lim.phenotype_scope", [
    V("paired genomes", PATP, PC13 + ["n_pairs"], "755", "int"),
], manual="22 of 54 media (per_genome.jsonl); 200-genome gapfilled subsample and 5 gated truth models are design/history (C12_SCOPE.md)")

# ================================================================ R4S2 score pass 2 (bound)
item("front:abs-phenotype-controls", [], manual="pre-registered C12 history (5 of 200 gated) and the A-C12/A-C13 reproduction counts (43,140; 0 of 8,614 vs 49) "
     "are in C12_SCOPE.md and c13_compare_*.json; re-derive from rescore_compare_{atp,gf}.json after R4.1")
item("picrust:PH-06", [
    V("draft rxn Jaccard PICRUSt2", PATP, ["arms", "picrust2", "controlled", "rxn_jaccard_draft", "median"], "0.737"),
    V("draft rxn Jaccard prFBA", PATP, ["arms", "prfba_exact_tie", "controlled", "rxn_jaccard_draft", "median"], "0.707"),
    V("draft rxn Jaccard nn1", PATP, ["arms", "base_nn1", "controlled", "rxn_jaccard_draft", "median"], "0.705"),
    V("essentiality disagreement PICRUSt2", PGF, ["arms", "picrust2", "controlled", "ess_shared_disagree_frac", "median"], "0.065"),
    V("essentiality disagreement prFBA", PGF, ["arms", "prfba_exact_tie", "controlled", "ess_shared_disagree_frac", "median"], "0.068"),
    V("essentiality disagreement nn1", PGF, ["arms", "base_nn1", "controlled", "ess_shared_disagree_frac", "median"], "0.069"),
    V("growth log-ratio prFBA", PGF, ["arms", "prfba_exact_tie", "controlled", "growth_log_ratio", "median"], "0.008"),
    V("growth log-ratio nn1", PGF, ["arms", "base_nn1", "controlled", "growth_log_ratio", "median"], "0.008"),
    V("growth log-ratio PICRUSt2", PGF, ["arms", "picrust2", "controlled", "growth_log_ratio", "median"], "-0.006"),
    V("pre-registered essentiality PICRUSt2", PGF0, ["arms", "picrust2", "controlled", "ess_shared_disagree_frac", "median"], "0.066"),
    V("pre-registered essentiality prFBA", PGF0, ["arms", "prfba_exact_tie", "controlled", "ess_shared_disagree_frac", "median"], "0.068"),
    V("pre-registered growth prFBA", PGF0, ["arms", "prfba_exact_tie", "controlled", "growth_log_ratio", "median"], "0.002"),
    V("pre-registered growth nn1", PGF0, ["arms", "base_nn1", "controlled", "growth_log_ratio", "median"], "0.004"),
], manual="gapfill burden 80-84 vs truth 82.5 (gapfill_burden medians per arm and truth); the arm assignment of the pre-registered pair "
          "(0.066 vs 0.068, 0.002 vs 0.004) follows the inventory and should be re-checked against C12_SCOPE.md L325")
item("picrust:AM-02", [], manual="A-C12/A-C13 properties (43,140 rows, 0 of 8,614 vs 49, 107/114 reversals, 865 -> 95 of 1,181 drafts, 138 of 200) are "
     "C12_SCOPE.md and c12fix/c13 compare JSONs: historical (D9) except what R4.1's rescore_compare files restate")
item("picrust:AM-03", [
    V("C2 source recovery", CTRL, ["C2", "measured", "prfba_source_recovery"], "1.0", "f1", tex_text="all 52{,}979"),
    V("C2 twin-free cells", CTRL, ["C2", "measured", "prfba_n_twin_free"], "52,979", "int"),
    V("C2 PICRUSt2 mean EC F1", CTRL, ["C2", "measured", "picrust2_mean_ec_f1"], "0.899"),
], manual="'nine controls evaluated, no hard failure' = controls.json {C*}.status")

# ================================================================ R1.5 tau 0.995
TR = ["_meta", "headline_macro_recall_pct_resolved"]
TPOP = ["_meta", "headline_macro_recall_pct_population"]
item("front:tab:tau-genecap", [
    V("cond 0.987 (min)", TAU, TR + ["exclude_0.987", {"each": R8}], "66.4", "pp1", reduce="min"),
    V("cond 0.987 (max)", TAU, TR + ["exclude_0.987", {"each": R8}], "69.9", "pp1", reduce="max"),
    V("cond 0.995 (min)", TAU, TR + ["exclude_0.995", {"each": R8}], "69.1", "pp1", reduce="min"),
    V("cond 0.995 (max)", TAU, TR + ["exclude_0.995", {"each": R8}], "72.9", "pp1", reduce="max"),
    V("pop 0.987 (min)", TAU, TPOP + ["exclude_0.987", {"each": R8}], "57.6", "pp1", reduce="min"),
    V("pop 0.987 (max)", TAU, TPOP + ["exclude_0.987", {"each": R8}], "64.4", "pp1", reduce="max"),
    V("pop 0.995 (min)", TAU, TPOP + ["exclude_0.995", {"each": R8}], "62.5", "pp1", reduce="min"),
    V("pop 0.995 (max)", TAU, TPOP + ["exclude_0.995", {"each": R8}], "69.0", "pp1", reduce="max"),
])
OCT = ["_meta", "outcome_counts_total"]
item("val16s:genecap_denominators", [
    X("score zero at 0.987", "(a+p)/n", {"a": S(GC, OCT + ["exclude_0.987", "abstained"]), "p": S(GC, OCT + ["exclude_0.987", "pred_empty"]),
                                        "n": S(CS, csm("prfba_exclude", "n"))}, "11.1", "pct1"),
    X("abstained at 0.987", "a/n", {"a": S(GC, OCT + ["exclude_0.987", "abstained"]), "n": S(CS, csm("prfba_exclude", "n"))}, "9.6", "pct1"),
    X("pred_empty at 0.987", "p/n", {"p": S(GC, OCT + ["exclude_0.987", "pred_empty"]), "n": S(CS, csm("prfba_exclude", "n"))}, "1.6", "pct1"),
    X("abstained at 0.995", "a/n", {"a": S(TAU, OCT + ["exclude_0.995", "abstained"]), "n": S(CS, csm("prfba_exclude", "n"))}, "6.0", "pct1"),
    X("pred_empty at 0.995", "p/n", {"p": S(TAU, OCT + ["exclude_0.995", "pred_empty"]), "n": S(CS, csm("prfba_exclude", "n"))}, "1.4", "pct1"),
], manual="denominator = scoreable amplicons (67,038 = 67,991 - truth_empty); after R1.3 the scoreable n itself moves (953 -> 278 truth_empty)")
item("val16s:tau0995_genecap", [
    V("FL cond 0.987", TAU, TR + ["exclude_0.987", "FullLength16S"], "68.7", "pp1"),
    V("FL cond 0.995", TAU, TR + ["exclude_0.995", "FullLength16S"], "71.5", "pp1"),
    V("FL pop 0.987", TAU, TPOP + ["exclude_0.987", "FullLength16S"], "63.3", "pp1"),
    V("FL pop 0.995", TAU, TPOP + ["exclude_0.995", "FullLength16S"], "68.7", "pp1"),
    V("V4 cond 0.987", TAU, TR + ["exclude_0.987", "V4"], "67.3", "pp1"),
    V("V4 cond 0.995", TAU, TR + ["exclude_0.995", "V4"], "69.7", "pp1"),
    V("V4 pop 0.987", TAU, TPOP + ["exclude_0.987", "V4"], "58.1", "pp1"),
    V("V4 pop 0.995", TAU, TPOP + ["exclude_0.995", "V4"], "63.1", "pp1"),
    X("pop gain (min)", "min([" + ", ".join(f"b{i}-a{i}" for i in range(8)) + "])",
      {**{f"a{i}": S(TAU, TPOP + ["exclude_0.987", r]) for i, r in enumerate(R8)}, **{f"b{i}": S(TAU, TPOP + ["exclude_0.995", r]) for i, r in enumerate(R8)}}, "4.1", "pp1"),
    X("pop gain (max)", "max([" + ", ".join(f"b{i}-a{i}" for i in range(8)) + "])",
      {**{f"a{i}": S(TAU, TPOP + ["exclude_0.987", r]) for i, r in enumerate(R8)}, **{f"b{i}": S(TAU, TPOP + ["exclude_0.995", r]) for i, r in enumerate(R8)}}, "5.5", "pp1"),
    X("resolved gain (min)", "min([" + ", ".join(f"b{i}-a{i}" for i in range(8)) + "])",
      {**{f"a{i}": S(TAU, TR + ["exclude_0.987", r]) for i, r in enumerate(R8)}, **{f"b{i}": S(TAU, TR + ["exclude_0.995", r]) for i, r in enumerate(R8)}}, "2.4", "pp1"),
    X("resolved gain (max)", "max([" + ", ".join(f"b{i}-a{i}" for i in range(8)) + "])",
      {**{f"a{i}": S(TAU, TR + ["exclude_0.987", r]) for i, r in enumerate(R8)}, **{f"b{i}": S(TAU, TR + ["exclude_0.995", r]) for i, r in enumerate(R8)}}, "3.0", "pp1"),
])

# ================================================================ R1.3L legacy POST
LEG_INC_OK = {"self_mode": "include", "outcome": "ok"}
item("methods:select_legacy_mean_size", [
    CSV("legacy include ok mean n_selected", GCCL, "csv_mean", {"filter": LEG_INC_OK, "column": "n_selected",
                                                                "usecols": ["self_mode", "outcome", "n_selected"]}, "1.32", "f2"),
])
item("methods:select_v4_include_precision", [
    V("union V4 include precision", GC, ["regions", "V4", "include", "macro", "precision", "value"], "0.795"),
    V("legacy V4 include precision", GCL, ["regions", "V4", "include", "macro", "precision", "value"], "0.856"),
])
item("methods:select_exclusion_precision_union_ge_legacy", [
    X("regions where union exclude precision >= legacy (count)", "sum([" + ", ".join(f"1 if u{i} >= l{i} else 0" for i in range(8)) + "])",
      {**{f"u{i}": S(GC, ["regions", r, "exclude_0.987", "macro", "precision", "value"]) for i, r in enumerate(R8)},
       **{f"l{i}": S(GCL, ["regions", r, "exclude_0.987", "macro", "precision", "value"]) for i, r in enumerate(R8)}}, "8", "int",
      tex_text="at or above"),
    X("narrowest margin (min union - legacy)", "min([" + ", ".join(f"u{i}-l{i}" for i in range(8)) + "])",
      {**{f"u{i}": S(GC, ["regions", r, "exclude_0.987", "macro", "precision", "value"]) for i, r in enumerate(R8)},
       **{f"l{i}": S(GCL, ["regions", r, "exclude_0.987", "macro", "precision", "value"]) for i, r in enumerate(R8)}}, "0", "ge",
      note="the claim is 'at or above' (margin >= 0); the V4 margin was 0.0017"),
])
item("methods:select_legacy_self_recovery_loss", [
    V("failures", SRD, ["n_failures"], "12,602", "int"),
    V("self in record", SRD, ["per_region", "ALL", "self_in_top20"], "67,466", "int"),
    V("fail rate", SRD, ["overall_fail_rate_given_top20"], "19", "pct0"),
])
item("val16s:genecap_precision", [
    V("union V4 include precision", GC, ["regions", "V4", "include", "macro", "precision", "value"], "0.795"),
    V("legacy V4 include precision", GCL, ["regions", "V4", "include", "macro", "precision", "value"], "0.856"),
    V("union FL include precision", GC, ["regions", "FullLength16S", "include", "macro", "precision", "value"], "1.00", "f2"),
    V("union V4 include precision (2 dp)", GC, ["regions", "V4", "include", "macro", "precision", "value"], "0.79", "f2"),
    V("V3-onward (min)", GC, _pr, "0.79", "f2", reduce="min"), V("V3-onward (max)", GC, _pr, "0.88", "f2", reduce="max"),
    V("V1-anchored (min)", GC, _pr1, "0.92", "f2", reduce="min"), V("V1-anchored (max)", GC, _pr1, "0.94", "f2", reduce="max"),
])
item("val16s:legacy_self_recovery", [
    CSV("legacy include ok mean n_selected", GCCL, "csv_mean", {"filter": LEG_INC_OK, "column": "n_selected",
                                                                "usecols": ["self_mode", "outcome", "n_selected"]}, "1.32", "f2"),
    V("fail rate", SRD, ["overall_fail_rate_given_top20"], "19", "pct0"),
    V("failures", SRD, ["n_failures"], "12,602", "int"),
    V("self in record", SRD, ["per_region", "ALL", "self_in_top20"], "67,466", "int"),
])

# ================================================================ R5A W6 derived numbers
item("methods:sharedref_tieset_coherence", [
    D("cross-phylum tie sets", "sr_tieset_span_phyla_unique", "6.4", "pct1"),
    D("same-genus tie sets", "sr_tieset_share_genus_unique", "35.8", "pct1"),
])
item("methods:recon_modal_fractional_weight", [D("modal k/n weight", "recon_modal_fractional_weight", "0.5", "f1")])
item("picrust:AS-10", [
    D("matched amplicons", "as10_n_matched", "59,607", "int"),
    D("prFBA matched F1", "as10_prfba_f1", "0.802"),
    D("PICRUSt2 matched F1", "as10_picrust2_f1", "0.727"),
    D("matched margin", "as10_delta", "+0.075", "+f3"),
])
item("mag:magval_identity_bins", [
    D("amplicon answer share 0.945-0.987 (low)", "mag_band_amp_min", "98", "pct0"),
    D("amplicon answer share 0.945-0.987 (high)", "mag_band_amp_max", "98", "pct0"),
    D("ANI answer share (low)", "mag_band_mag_low", "8", "pct0"),
    D("ANI answer share (high)", "mag_band_mag_high", "20", "pct0"),
])
item("mag:magval_scope_claims", [
    D("conditional coverage", "mag_resolve_rate_pct", "31", "pct0"),
    D("MAG abstains on", "mag_abstain_rate_pct", "69", "pct0"),
    D("amplicon F1 where MAG abstains, FL", "magvs16s_amp_f1_where_mag_abstains_FL", "0.64", "f2"),
    D("amplicon F1 where MAG abstains, V4", "magvs16s_amp_f1_where_mag_abstains_V4", "0.54", "f2"),
    V("FL pop 16S", PC, pc("FullLength16S", block="population") + ["amp"], "0.605"),
    V("FL pop MAG", PC, pc("FullLength16S", block="population") + ["mag"], "0.248"),
    D("MAG-answered easier for amplicon", "magvs16s_mag_answered_easier", "0.053"),
    D("MAG advantage (min over regions x levels)", "magvs16s_delta_min_all_levels", "0.08", "f2"),
    D("MAG advantage (max over regions x levels)", "magvs16s_delta_max_all_levels", "0.21", "f2"),
])
item("mag:magval_tau99_paired", [
    X("FL delta at tau 95", "a-b", {"a": S(PC, pc("FullLength16S") + ["mag"]), "b": S(PC, pc("FullLength16S") + ["amp"])}, "+0.103", "+f3"),
    D("FL delta at tau 99 (d100 pairing)", "magvs16s_tau99_delta_FL", "+0.127", "+f3",
      note="becomes mag_validation/out/paired_tau99_d100.json after R2.3 step 3"),
    V("pop MAG tau 95", PC, pc("FullLength16S", block="population") + ["mag"], "0.248"),
    D("pop MAG tau 99", "magvs16s_tau99_pop_mag_FL", "0.301"),
    D("matched-stringency gap", "magvs16s_gap_matched", "0.354"),
    V("pop amplicon", PC, pc("FullLength16S", block="population") + ["amp"], "0.605"),
    V("FL pop F1 at 0.995", TAU, ["regions", "FullLength16S", "exclude_0.995", "f1_population", "value"], "0.655"),
    X("conservative gap", "b-a", {"a": S(PC, pc("FullLength16S", block="population") + ["mag"]),
                                  "b": S(PC, pc("FullLength16S", block="population") + ["amp"])}, "0.357"),
])
item("panels:panels.welch_tests", [
    D("FL original - k-means cond F1", "panels_welch_FL_delta", "+0.011", "+f3"),
    D("FL p", "panels_welch_FL_p", "2.3e-4", "e1"),
    D("all eight regions p < 0.01", "panels_welch_all_regions_p_lt_0.01", "True", "str"),
    D("regions with p < 1e-8", "panels_welch_regions_p_lt_1e-8", "V4,V4-V5,V7-V9", "str"),
    D("every MAG completeness p < 1e-3", "panels_welch_mag_every_completeness_p_lt_1e-3", "True", "str"),
])
item("panels:panels.mag_singleton_abstention", [
    D("k-means singleton abstention", "panels_mag_single_abst_kmeans", "0.6895", "f4"),
    D("k-means n", "panels_mag_single_n_kmeans", "3,942", "int"),
    D("original singleton abstention", "panels_mag_single_abst_original", "0.6883", "f4"),
    D("original n", "panels_mag_single_n_original", "9,580", "int"),
    D("difference", "panels_mag_single_diff", "+0.001", "+f3"),
    D("CI lo", "panels_mag_single_ci_lo", "-0.016"),
    D("CI hi", "panels_mag_single_ci_hi", "+0.018", "+f3"),
])
item("panels:panels.mag_multiplicity", [
    D(">=20 members abstention", "panels_mag_ge20_abst_kmeans", "0.559"),
    D("k-means overall", "panels_mag_species_level_abst_kmeans", "0.658"),
    D("original overall", "panels_mag_species_level_abst_original", "0.688"),
    D("k-means without >=20", "panels_mag_without_ge20_abst_kmeans", "0.687"),
])

# ================================================================ R1.K deployments (cite the new *_bvbrc stats file)
CPH = ["cross_phylum"]
_dep = lambda label, sel, p, fmt="int": V(label, DEP, sel, p, fmt, artifact_new=DEP_NEW)
DEPV = [
    _dep("resolved merges", CPH + ["all_multi", "resolved"], "1,238"),
    _dep("multi-genome merges", CPH + ["all_multi", "multi_genome_merges"], "1,410"),
    _dep("two-genome cross-phylum", CPH + ["n=2", "cross_phylum"], "145"),
    _dep("two-genome resolved", CPH + ["n=2", "resolved"], "454"),
    _dep("two-genome pct", CPH + ["n=2", "cross_phylum_pct"], "31.9", "pp1"),
    _dep("larger cross-phylum", CPH + ["n>=3", "cross_phylum"], "488"),
    _dep("larger resolved", CPH + ["n>=3", "resolved"], "784"),
    _dep("larger pct", CPH + ["n>=3", "cross_phylum_pct"], "62.2", "pp1"),
    _dep("Bacteria+Archaea", CPH + ["all_multi", "bacteria_archaea"], "20"),
]
item("methods:deployment_cross_phylum", DEPV)
item("methods:e3_production_unresolved_share", [
    V("unresolved feature share", BRS, ["deployments", "emilykin", "e3", "feature_resolution_shares", "unresolved"], "11", "pct0"),
    V("multi-genome synthetic genomes", BRS, ["deployments", "emilykin", "e3", "n_applied"], "1,410", "int"),
])
item("val16s:deployment_cross_phylum", DEPV)
item("val16s:deployment_e3_unresolved_share", [
    V("unresolved feature share", BRS, ["deployments", "emilykin", "e3", "feature_resolution_shares", "unresolved"], "11", "pct0"),
    V("ASVs", DEP, ["n_asvs"], "3,950", "int", artifact_new=DEP_NEW),
])
item("panels:disc.selection_size_deployment", [
    _dep("single-genome share", ["n_genomes_1_pct_of_answered"], "62.5", "pp1"),
    _dep("single-genome ASVs", ["n_genomes_1"], "2,345"),
    _dep("answered", ["answered"], "3,755"),
    _dep("mean selection size", ["mean_n_genomes_answered"], "2.289", "f3"),
])
item("panels:disc.phylum_span", DEPV[:8], manual="controlled 6.4% / 4.8%: see methods:sharedref_tieset_coherence (W6)")
item("panels:lim.deployment_truncation", [
    _dep("single-genome share", ["n_genomes_1_pct_of_answered"], "62.5", "pp1"),
    _dep("mean selection size", ["mean_n_genomes_answered"], "2.289", "f3"),
    V("largest selection", DEP, ["n_genomes_distribution"], "77", "int", artifact_new=DEP_NEW, reduce="maxkey"),
    _dep("ASVs", ["n_asvs"], "3,950"),
], manual="25 ASVs whose 100th hit ties the best, and the depth-20 comparator (2.034, 19) are in deployment_d100/RESULTS.md")

# ================================================================ R2.1 E3 gate
E3N = ["results", "exclude_0.987", "n2+"]
item("methods:e3_gate_brier_auroc", [
    V("k/n Brier", E3, E3N + ["estimators", "kn", "brier"], "0.155"),
    V("E3 Brier", E3, E3N + ["estimators", "e3_species", "brier"], "0.102"),
    V("delta Brier", E3, E3N + ["delta_vs_kn", "e3_species", "brier"], "-0.053"),
    V("  CI lo", E3, E3N + ["delta_vs_kn", "e3_species", "brier_ci95", 0], "-0.055"),
    V("  CI hi", E3, E3N + ["delta_vs_kn", "e3_species", "brier_ci95", 1], "-0.051"),
    V("k/n AUROC", E3, E3N + ["estimators", "kn", "auroc"], "0.874"),
    V("E3 AUROC", E3, E3N + ["estimators", "e3_species", "auroc"], "0.928"),
    V("delta AUROC", E3, E3N + ["delta_vs_kn", "e3_species", "auroc"], "+0.054", "+f3"),
    V("  CI lo", E3, E3N + ["delta_vs_kn", "e3_species", "auroc_ci95", 0], "+0.052", "+f3"),
    V("  CI hi", E3, E3N + ["delta_vs_kn", "e3_species", "auroc_ci95", 1], "+0.056", "+f3"),
    V("amplicons", E3, E3N + ["n_amplicons"], "19,158", "int"),
    V("source genomes", E3, E3N + ["n_sources"], "6,776", "int"),
])
_rec = X("recalibration share of the Brier gain", "(k-s)/(-d)", {"k": S(E3, E3N + ["estimators", "kn", "brier"]),
                                                                 "s": S(RR, ["kn_saturated_recalibration", "brier"]),
                                                                 "d": S(E3, E3N + ["delta_vs_kn", "e3_species", "brier"])}, "45", "pct0",
         known_mismatch="'about 45%' is 0.0236/0.0530 from rounded inputs; the unrounded share is 44.5%")
item("methods:e3_gate_recalibration_share", [
    _rec,
    V("E3 vs saturated k/n Brier", RR, ["kn_saturated_recalibration", "e3s_minus_knsat_brier"], "0.029", "absf3"),
    V("E3 vs saturated k/n AUROC", RR, ["kn_saturated_recalibration", "e3s_minus_knsat_auroc"], "0.048"),
])
item("val16s:e3_production_gate", ITEMS["methods:e3_gate_brier_auroc"]["values"] + [
    _rec,
    V("E3 vs saturated k/n Brier", RR, ["kn_saturated_recalibration", "e3s_minus_knsat_brier"], "0.029", "absf3"),
    V("E3 vs saturated k/n AUROC", RR, ["kn_saturated_recalibration", "e3s_minus_knsat_auroc"], "0.048"),
])

# ================================================================ R1.4 d20 trees
_l20 = {f"a{i}": S(GC20, ["regions", r, "exclude_0.987", "macro", "recall", "value"]) for i, r in enumerate(R8)}
_l100 = {f"b{i}": S(GCL, ["regions", r, "exclude_0.987", "macro", "recall", "value"]) for i, r in enumerate(R8)}
_f20 = {f"c{i}": S(GC20, ["regions", r, "exclude_0.987", "macro", "f1", "value"]) for i, r in enumerate(R8)}
_f100 = {f"d{i}": S(GCL, ["regions", r, "exclude_0.987", "macro", "f1", "value"]) for i, r in enumerate(R8)}
item("val16s:legacy_genecap_d20_vs_d100", [
    X("cond macro recall d100-d20 (min, points)", "min([" + ", ".join(f"b{i}-a{i}" for i in range(8)) + "])", {**_l20, **_l100}, "0.7", "pct1"),
    X("cond macro recall d100-d20 (max, points)", "max([" + ", ".join(f"b{i}-a{i}" for i in range(8)) + "])", {**_l20, **_l100}, "1.3", "pct1"),
    X("cond gene-capture F1 d20-d100 (min drop, points)", "min([" + ", ".join(f"c{i}-d{i}" for i in range(8)) + "])", {**_f20, **_f100}, "0.1", "pct1"),
    X("cond gene-capture F1 d20-d100 (max drop, points)", "max([" + ", ".join(f"c{i}-d{i}" for i in range(8)) + "])", {**_f20, **_f100}, "1.2", "pct1"),
])

# ================================================================ D7 backend concordance (historical)
item("val16s:backend_concordance", [], manual="historical June check (D7): wording only; edlib_vs_gpu_assessment/assessment_summary.json is not re-run")

# ================================================================ R2.2 MAG re-score
item("mag:magval_selfinclusion_ceiling", [
    V("include coverage 100%", MSF, ["summary", "c1_x0_include", "coverage_answered"], "1.000"),
    V("include gene capture 100%", MSF, ["summary", "c1_x0_include", "genecap_cond_f1", "mean"], "1.000"),
    CSV("include anchor is the source", MAGC, "csv_share", {"filter": {"regime": "include", "outcome": "ok", "contamination": 0}, "where": {"anchor": ["eqcol", "genome_id"]},
                                                           "usecols": ["regime", "outcome", "anchor", "genome_id", "contamination"]}, "99.9", "pct1"),
], manual="genus accuracy 1.000 (summary.c*_x0_include Genus block); 'every level' spans all completeness keys")
item("mag:magval_tau99_mag_capture", [
    V("resolve tau 95", MSF, ["summary", "c1_x0_exclude", "coverage_answered"], "0.308"),
    V("resolve tau 99", MST, ["summary", "c1_x0_exclude", "coverage_answered"], "0.366"),
    X("relative gain", "b/a-1", {"a": S(MSF, ["summary", "c1_x0_exclude", "coverage_answered"]),
                                 "b": S(MST, ["summary", "c1_x0_exclude", "coverage_answered"])}, "+19", "+pct0"),
    V("cond F1 tau 95", MSF, ["summary", "c1_x0_exclude", "genecap_cond_f1", "mean"], "0.807"),
    V("cond F1 tau 99", MST, ["summary", "c1_x0_exclude", "genecap_cond_f1", "mean"], "0.823"),
], manual="genus accuracy 0.910 -> 0.895 and the 1,285 retained references are in RESULTS_tau99.md")
_c5 = [f"c{c}_x0_exclude" for c in ("1", "0.9", "0.7", "0.5", "0.3")]
_n100 = {**{f"r{i}": S(MSF, ["summary", k, "coverage_answered"]) for i, k in enumerate(_c5)},
         **{f"f{i}": S(MSF, ["summary", k, "genecap_cond_f1", "mean"]) for i, k in enumerate(_c5)}}
_n500 = {**{f"r{i}": S(MPC, ["rows", {"where": {"regime": "exclude", "completeness": c}}, "orig_crude"]) for i, c in enumerate((1.0, 0.9, 0.7, 0.5, 0.3))},
         **{f"f{i}": S(MPC, ["rows", {"where": {"regime": "exclude", "completeness": c}}, "orig_condF1"]) for i, c in enumerate((1.0, 0.9, 0.7, 0.5, 0.3))},
         **{f"n{i}": S(MPC, ["rows", {"where": {"regime": "exclude", "completeness": c}}, "orig_n_resolved"]) for i, c in enumerate((1.0, 0.9, 0.7, 0.5, 0.3))}}
_mean5 = "sum([r0, r1, r2, r3, r4]) / 5"
item("mag:twopanel_mag_depth500_original", [
    X("resolve depth 100 (pooled 100-30%)", _mean5, _n100, "0.3071", "f4"),
    X("resolve depth 500 (pooled 100-30%)", _mean5, _n500, "0.3089", "f4"),
    X("cond F1 depth 100 (pooled, unweighted)", "sum([f0, f1, f2, f3, f4]) / 5", _n100, "0.8066", "f4",
      note="mag_summary_full has no per-level resolved counts; equal weights (n_resolved differs by <1%)"),
    X("cond F1 depth 500 (pooled, resolved-weighted)", "sum([f0*n0, f1*n1, f2*n2, f3*n3, f4*n4]) / sum([n0, n1, n2, n3, n4])", _n500, "0.8061", "f4"),
], manual="49,185 scoreable cells per regime and the 88/87/1 cell comparison are ad hoc (paper f8047c5); re-derive from out_batches vs out_batches_n500")
item("panels:panels.mag_depth500_rerun", ITEMS["mag:twopanel_mag_depth500_original"]["values"],
     manual=ITEMS["mag:twopanel_mag_depth500_original"]["manual"])


# ---------------------------------------------------------------- assemble
def main():
    claims = json.load(open(CLAIMS))
    out, missing = {}, []
    for stage, items in claims.items():
        if stage == "NO":
            continue
        for it in items:
            iid = f"{it['inv']}:{it['key']}"
            spec = ITEMS.get(iid)
            if spec is None:
                missing.append(iid)
                spec = {"values": [], "manual": "no extractor written"}
            vals = [v for v in spec["values"] if v.get("printed") is not None]
            out[iid] = {"stage": stage, "inv": it["inv"], "key": it["key"], "tex_lines": it["tex_lines"],
                        "claim": it["claim"], "values": vals}
            if spec.get("manual"):
                out[iid]["manual"] = spec["manual"]
    extra = sorted(set(ITEMS) - set(out))
    meta = {"what": "W7 extractor spec: per dependent manuscript claim, where each cited number lives",
            "built_by": "$W/tools/build_claims_spec.py (W=/home/freiburger/Documents/prFBA_rescore_2026-11)",
            "claims_source": "re-score plan claims_by_stage.json (178 dependent items)",
            "n_items": len(out), "n_values": sum(len(v["values"]) for v in out.values()),
            "n_items_without_values": sum(1 for v in out.values() if not v["values"]),
            "derived_keys_note": "derived:<key> values are read from derived_numbers.{pre,post}.json; collect.py falls back to "
                                 "matching the W6 entry by printed value + tex lines when the key differs"}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT + ".tmp", "w") as fh:
        json.dump({"_meta": meta, "items": out}, fh, indent=1)
    os.replace(OUT + ".tmp", OUT)
    print(json.dumps(meta, indent=1))
    if missing:
        print("items without extractors:", missing)
    if extra:
        print("extractors for unknown items:", extra)


if __name__ == "__main__":
    main()
