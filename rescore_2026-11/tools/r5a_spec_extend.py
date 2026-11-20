#!/usr/bin/env python
"""R5.1: extend manuscript_numbers/claims_spec.json (W7) so that the 'manual' parts of its items are
either (a) real values (W6 derived keys, JSON selectors or expressions), (b) hand-carried values with a
named source (`manual_values`), or (c) an explicit resolution (`manual_resolution`: unchanged by
construction / historical record / figure / not recomputed, with the reason). Also adds the curated
`meaning_change` notes (claims whose meaning, not just digits, changes in the re-score).

Idempotent: re-running from the committed spec gives the same file (values are keyed by label).
    ~/Documents/py_venv/bin/python $W/tools/r5a_spec_extend.py
"""
import json
import os

R = "/home/freiburger/Documents/prFBA"
SPEC = os.path.join(R, "manuscript_numbers/claims_spec.json")

SR = "region_validation/shared_reference/"
SG = SR + "data/summary_genus.json"
SF = SR + "data_family/summary_family.json"
PGF = SR + "data/phenotype_gf_c13/phenotype_summary.json"
CTRL = SR + "data/controls.json"
WOBJ = SR + "data/weighted_objective/summary.json"
MSF = "mag_validation/out/mag_summary_full.json"
MST = "mag_validation/out/mag_summary_tau99.json"
AFG = "mag_validation/out/af_gate_counterfactual.json"
PA = "three_panel/three_panel_mag_panelarm.json"
TP = "three_panel/three_panel_comparison_d100.json"


def D(key, label, printed, fmt, **kw):
    v = {"label": label, "printed": printed, "fmt": fmt, "artifact": "derived:" + key}
    v.update(kw)
    return v


def J(art, sel, label, printed, fmt, **kw):
    v = {"label": label, "printed": printed, "fmt": fmt, "artifact": art, "selector": sel}
    v.update(kw)
    return v


def X(expr, vars_, label, printed, fmt, **kw):
    v = {"label": label, "printed": printed, "fmt": fmt, "expr": expr, "vars": vars_}
    v.update(kw)
    return v


def pop(art, arm):
    return {"artifact": art, "selector": ["arms", arm, "ec", "population", "f1", "value"]}


def MV(label, printed, new, source, note=None):
    d = {"label": label, "printed": printed, "new": new, "source": source}
    if note:
        d["note"] = note
    return d


K_NN1_ECF1 = ("printed 0.781 has no surviving recipe; the cells-table recipe gives 0.782 for both arms on the "
              "pre-repair pass (see derived_numbers KNOWN_PRE_MISMATCH)")
K_NN1_P = ("no paired test reproduces 0.64; the recipe's Wilcoxon gives 0.45 on the pre-repair pass: replace the "
           "printed p with the recipe value or drop it")

ADD = {
    # ---------------------------------------------------------------- phenotype (R4P / R4S2)
    "front:abs-phenotype-readout": dict(values=[
        D("pheno_n_paired_v4", "of which V4", "747", "int"),
        D("pheno_ecf1_delta_vs_picrust2", "EC-F1 delta prFBA - PICRUSt2 on the paired genomes", "-0.028", "+f3"),
        D("pheno_mcnemar_n_favour_picrust2", "media favouring PICRUSt2 (McNemar q<=0.05)", "2", "int"),
        D("pheno_media_with_signal", "media with ATP signal", "22", "int"),
        D("pheno_media_in_vector", "media in the ATP vector", "54", "int")], manual=None),
    "picrust:PH-01": dict(values=[
        D("pheno_n_paired_v4", "of which V4", "747", "int"),
        D("pheno_ecf1_delta_vs_picrust2", "EC-F1 delta prFBA - PICRUSt2", "-0.028", "+f3")], manual=None),
    "picrust:PH-02": dict(values=[
        D("pheno_mcnemar_n_q005", "media with McNemar q<=0.05", "2", "int"),
        D("pheno_mcnemar_n_favour_picrust2", "of which favour PICRUSt2", "2", "int"),
        D("pheno_mcnemar_media_favour_picrust2", "media favouring PICRUSt2", "Ac.NO,Etho.O2", "str"),
        D("pheno_media_with_signal", "testable media (with ATP signal)", "22", "int"),
        D("pheno_mcnemar_q_max", "largest q of those media", "0.002", "f3")], manual=None),
    "picrust:PH-03": dict(values=[D("pheno_spearman_p", "Spearman p", "2e-4", "e0")], manual=None),
    "picrust:PH-04": dict(values=[
        D("pheno_nn1_ecf1_prfba", "EC F1 prFBA on the nn1-paired genomes", "0.782", "f3"),
        D("pheno_nn1_ecf1_nn1", "EC F1 nn1 on the nn1-paired genomes", "0.781", "f3", known_mismatch=K_NN1_ECF1),
        D("pheno_nn1_ecf1_p", "EC F1 paired test p", "0.64", "f2", known_mismatch=K_NN1_P)], manual=None),
    "panels:disc.atp_vs_nn1": dict(values=[
        D("pheno_nn1_ecf1_prfba", "EC F1 prFBA on the nn1-paired genomes", "0.782", "f3"),
        D("pheno_nn1_ecf1_nn1", "EC F1 nn1 on the nn1-paired genomes", "0.781", "f3", known_mismatch=K_NN1_ECF1),
        D("pheno_nn1_ecf1_p", "EC F1 paired test p", "0.64", "f2", known_mismatch=K_NN1_P)], manual=None),
    "picrust:PH-05": dict(values=[
        D("pheno_media_with_signal", "media with ATP signal", "22", "int"),
        D("pheno_media_in_vector", "media in the ATP vector (53 + empty control)", "54", "int")], manual=None),
    "panels:disc.refusal_credit": dict(values=[
        D("sr_genus_exact_tie_vs_mask", "exact-tie vs its own masked run", "+0.019", "+f3"),
        D("pheno_spearman_p", "Spearman p", "2e-4", "e0"),
        D("pheno_ecf1_delta_vs_picrust2", "EC-F1 delta prFBA - PICRUSt2", "-0.028", "+f3"),
        X("(f-a)/(p-a)", {"f": pop(SG, "prfba_nestedcv_forced"), "a": pop(SG, "prfba_nestedcv"),
                           "p": pop(SG, "picrust2")},
          "share of the deficit recovered by answering everywhere (printed 'two thirds')", None, "f2")],
        manual=None),
    "panels:lim.phenotype_scope": dict(values=[
        D("pheno_media_with_signal", "media with ATP signal", "22", "int"),
        D("pheno_media_in_vector", "media in the ATP vector", "54", "int"),
        D("pheno_gf_n_genomes", "gapfilled subsample", "200", "int")],
        manual_resolution={"status": "historical", "parts": "5 gated truth models",
                           "why": "the pre-registered C12 failure (5 of 200 gated) is a fixed history (D9), not re-run"}),
    "front:abs-phenotype-controls": dict(values=[
        D("pheno_atp_rows", "ATP-pass rows (A-C13 configuration, re-run)", "43,140", "int"),
        D("pheno_atp_rows_empty_medium_atp", "ATP-pass rows with empty-medium ATP", "0", "int"),
        D("pheno_gf_rows", "gapfilled-pass rows", "8,614", "int"),
        D("pheno_gf_gated", "gapfilled rows gated (C12 leg c)", "0", "int")],
        manual_resolution={"status": "historical", "parts": "5 of 200 gated (pre-registered); 49 gated before A-C12",
                           "why": "D9: the pre-registered and A-C12 passes are historical records; only the A-C13 "
                                  "configuration was re-run (R4.1), on a re-drawn sample"}),
    "picrust:AM-02": dict(values=[
        D("pheno_atp_rows", "ATP-pass rows (A-C13 configuration, re-run)", "43,140", "int"),
        D("pheno_atp_rows_empty_medium_atp", "ATP-pass rows with empty-medium ATP", "0", "int"),
        D("pheno_gf_rows", "gapfilled-pass rows", "8,614", "int"),
        D("pheno_gf_gated", "gapfilled rows gated", "0", "int"),
        D("pheno_gf_models_reversal", "models applying a reversal (A-C13)", "114", "int")],
        manual_resolution={"status": "historical",
                           "parts": "49 gated before; 136 models in 50 genomes; 107 reversals under A-C12; 166 of 166; "
                                    "865 -> 95 of 1,181 drafts; 138 of 200 truth drafts; 218 rebuilt models of 20 "
                                    "genomes; 72 of 72; medians move <= 0.004; burden <= 1 (1.5 as shipped)",
                           "why": "properties of the original A-C12 / A-C13 runs and side studies on the 10-29 sample "
                                  "(D9, not re-run). The re-run A-C13 pass is on a re-drawn sample (442 of 1,000 ATP and "
                                  "62 of 200 gapfilled genomes kept), so only the counts above are current; the "
                                  "'reproduces the pre-registered pass exactly (43,140 of 43,140 rows)' sentence is a "
                                  "statement about the 10-29 sample and cannot be restated for the new one"}),
    "picrust:AM-03": dict(values=[
        J(CTRL, [c, "status"], "%s status" % c, "pass", "str") for c in
        ("C3", "C4", "C6", "C7", "C8", "C10", "C13", "C14", "C15")], manual=None),
    "picrust:PH-06": dict(values=[
        J(PGF, ["arms", "truth", "controlled", "gapfill_burden", "median"], "gapfill burden truth (median)", "82.5",
          "f1"),
        J(PGF, ["arms", {"each": ["picrust2", "prfba_exact_tie", "base_nn1", "base_family_donor"]}, "controlled",
                "gapfill_burden", "median"], "gapfill burden method arms (min median)", "80", "f0", reduce="min"),
        J(PGF, ["arms", {"each": ["picrust2", "prfba_exact_tie", "base_nn1", "base_family_donor"]}, "controlled",
                "gapfill_burden", "median"], "gapfill burden method arms (max median)", "84", "f0", reduce="max")],
        manual_resolution={"status": "historical", "parts": "pre-registered pair 0.066 vs 0.068, 0.002 vs 0.004",
                           "why": "pre-registered gapfilled pass (phenotype_gf/, D9); unchanged file"}),
    # ---------------------------------------------------------------- depth 500 (R2.2)
    "mag:twopanel_mag_depth500_original": dict(values=[
        D("mag_d500_n_scoreable", "scoreable cells per regime", "49,185", "int"),
        D("mag_d500_n_outcome_changed", "cells changing outcome", "88", "int"),
        D("mag_d500_n_newly_resolved", "newly resolved", "87", "int"),
        D("mag_d500_n_geneless", "gene-less reference", "1", "int"),
        D("mag_d500_n_selection_switches", "selection switches", "0", "int")], manual=None),
    "panels:panels.mag_depth500_rerun": dict(values=[
        D("mag_d500_n_scoreable", "scoreable cells per regime", "49,185", "int"),
        D("mag_d500_n_outcome_changed", "cells changing outcome", "88", "int"),
        D("mag_d500_n_newly_resolved", "newly resolved", "87", "int"),
        D("mag_d500_n_geneless", "gene-less reference", "1", "int"),
        D("mag_d500_n_selection_switches", "selection switches", "0", "int")], manual=None),
    # ---------------------------------------------------------------- shared reference (R3 / R3F)
    "picrust:SR-09": dict(values=[
        D("sr_genus_legacy_forced_masked", "legacy forced, random-masked", "0.668", "f3"),
        D("sr_genus_legacy_like_for_like", "like-for-like credit", "+0.019", "+f3")], manual=None),
    "front:abs-controlled-random-mask": dict(values=[
        D("sr_genus_exact_tie_vs_mask", "exact-tie vs own mask (genus)", "+0.019", "+f3"),
        D("sr_family_exact_tie_vs_mask", "exact-tie vs own mask (family)", "+0.015", "+f3")], manual=None),
    "picrust:FR-03": dict(values=[
        D("sr_family_exact_tie_vs_mask", "exact-tie vs own mask", "+0.015", "+f3"),
        D("sr_family_legacy_like_for_like", "legacy like-for-like", "+0.016", "+f3"),
        X("a-b", {"a": pop(SF, "prfba_exact_tie"), "b": pop(SF, "base_nn1")}, "exact-tie vs nn1", "-0.083", "+f3"),
        X("a-b", {"a": pop(SF, "prfba_exact_tie"), "b": pop(SF, "base_constant")}, "exact-tie vs constant",
          "-0.020", "+f3"),
        J(SF, ["arms", "prfba_exact_tie", "ec", "skill_score"], "exact-tie skill", "-0.16", "f2"),
        X("a-b", {"a": pop(SF, "prfba_exact_tie"), "b": pop(SF, "base_family_donor")}, "exact-tie vs family donor",
          "+0.104", "+f3")], manual=None),
    "methods:sharedref_ngenomes_vs_identity": dict(values=[
        D("sr_p_multi_monotone", "P(n>=2) rises monotonically", "yes", "str", field="formatted"),
        D("sr_p_multi_low_bin", "P(n>=2) lowest novelty bin (not printed)", None, "f3"),
        D("sr_p_multi_high_bin", "P(n>=2) highest novelty bin (not printed)", None, "f3")], manual=None),
    "val16s:controlled_selection_sizes": dict(values=[
        D("sr_tieset_span_phyla_records", "multi-genome selections spanning phyla", "4.8", "pct1"),
        D("sr_p_multi_monotone", "P(n>=2) rises monotonically", "yes", "str", field="formatted")], manual=None),
    "picrust:SR-02": dict(values=[
        J(SF, ["meta", "notes", "nested_cv_winners", {"each": "*"}], "family nested-CV winners legacy (count)", "5",
          "int", reduce="len")],
        manual_resolution={"status": "not recomputed", "parts": "inner margins 0.001-0.003",
                           "why": "inner-CV F1 gaps (upper end by hand); the winners are unchanged (legacy in 5 of 5 "
                                  "folds at genus, R3.5)"}),
    # ---------------------------------------------------------------- as shipped (R2.4 / R2.5)
    "picrust:AS-17": dict(values=[
        D("as17_all_below_0987", "all answers below 0.987 identity", "yes", "str", field="formatted"),
        D("as17_median_identity", "median reference identity", "0.956", "f3")], manual=None),
    "picrust:AS-18": dict(values=[
        D("as18_n_further_below_floor", "further amplicons with surviving hits all below the floor", "145", "int"),
        D("as18_n_floor_explains", "abstentions the floor explains", "692", "int"),
        D("as18_floor_share_of_abstentions", "share of outright abstentions", "10.8", "pct1")], manual=None),
    "front:fig:picrust2": dict(values=[
        D("as18_n_further_below_floor", "further amplicons with surviving hits all below the floor", "145", "int"),
        D("as18_n_floor_explains", "abstentions the floor explains", "692", "int"),
        D("as18_floor_share_of_abstentions", "share of outright abstentions", "10.8", "pct1")], manual=None),
    "picrust:AS-19": dict(values=[
        D("as19_share_more_precise_than_reference", "in-reference amplicons more precise than the reference annotation",
          "43", "pct0")], manual=None),
    "front:tab:picrust2": dict(values=[
        D("as19_share_more_precise_than_reference", "in-reference amplicons more precise than the reference annotation",
          "43", "pct0")], manual=None),
    "picrust:AS-23": dict(values=[D("as23_median_selection_size", "median self-excluded selection size", "1", "int",
                                    tex_text="a median size of one reference")],
                         manual=None),
    "picrust:AS-07": dict(manual_resolution={
        "status": "unchanged", "parts": "2,720 of 9,746 source genomes in PICRUSt2's reference",
        "why": "the source-genome set and PICRUSt2's reference are fixed: bvbrc_source_genomes.json byte-equal after "
               "R2.4 and annotation_ceiling.json _meta.n_source_assembly_in_picrust2_ref = 2,720 of 9,746 in both"}),
    "picrust:AS-05": dict(manual_resolution={
        "status": "historical",
        "parts": "675 of 953, 818 of 1,060, 167 genomes, 0.75-0.86 / 0.05 points",
        "why": "the 11-13 repair measurement (CORRECTIONS_2026-11-13.json, genecap_repaired_impact/impact.json). After the "
               "re-score the reported state is the repaired one (truth_empty 278, pred_empty 243), so L517 must become a "
               "past-tense correction note (plan section 5)"}),
    "picrust:AS-15": dict(manual_resolution={
        "status": "not recomputed",
        "parts": "0.026 archaeal, 0.014 distant novelty bins, 0.05 / 0.03 in the 9- and 8-amplicon cells",
        "why": "per-stratum rxn-minus-EC differences with a domain split of comparison_cells; no recipe written in this "
               "run. The overall rxn-EC tracking (0.805 / 0.721) is extracted above; the manuscript agent should "
               "re-derive the four divergences or soften them to 'small subsets'"}),
    "val16s:depth_probe_quality": dict(manual_resolution={
        "status": "not recomputed", "parts": "0.651 / 0.800 answered-set anchor accuracy",
        "why": "no stored key or surviving recipe located (not in depth_probe_results.json or any summary); the d20 "
               "exclude-mode taxonomy is unchanged by the re-score (R1.4), so the values are expected unchanged, "
               "but this is not verified"}),
    "val16s:population_ec_and_projections": dict(manual_resolution={
        "status": "not recomputed", "parts": "projections ~0.96 / ~0.49 / ~0.75",
        "why": "hand projections in picrust2_comparison/README.md from the probe's per-region depth-500 rates; the probe "
               "recovery counts are unchanged (283/320) but the EC quality moved (0.783 -> 0.829), so the ~0.75 EC "
               "projection should be re-derived by the manuscript agent or dropped"}),
    "mag:magval_af_gate_counterfactual": dict(values=[
        X("a/b", {"a": {"artifact": AFG, "selector": ["summary", "newly_above_full_length_16S_level",
                                                      {"each": "*"}], "reduce": "min"},
                  "b": {"artifact": AFG, "selector": ["summary", "n_newly_admitted"]}},
          "newly admitted above the FL 16S level (share; min over the three levels)", None, "pct1"),
        X("a/b", {"a": {"artifact": AFG, "selector": ["summary", "n_newly_admitted"]},
                  "b": {"artifact": AFG, "selector": ["summary", "n_already_answered"]}},
          "coverage gain", "+45", "+pct0"),
        X("(na*fa+nn*fn)/(na+nn)", {"na": {"artifact": AFG, "selector": ["summary", "n_already_answered"]},
                                    "fa": {"artifact": AFG, "selector": ["summary", "mean_f1_already"]},
                                    "nn": {"artifact": AFG, "selector": ["summary", "n_newly_admitted"]},
                                    "fn": {"artifact": AFG, "selector": ["summary", "mean_f1_newly"]}},
          "pooled mean F1", "0.742", "f3")],
        manual_resolution={"status": "not recomputed", "parts": "1,386 newly admitted incl. 19 gene-less; 26.7% above 0.710",
                           "why": "the summary stores the non-empty count (n_newly_admitted) and per-level counts; the "
                                  "share above the union FL level is 369/1,391 = 26.5% new (Run 2), printed 26.7%; the "
                                  "1,386 total (incl. gene-less picks) needs the per-MAG list"}),
    "mag:magval_selfinclusion_ceiling": dict(values=[
        J(MSF, ["summary", {"each": ["c1_x0_include", "c0.9_x0_include", "c0.7_x0_include", "c0.5_x0_include",
                                     "c0.3_x0_include"]}, "Genus", "correct"],
          "include genus accuracy (min over completeness)", "1.000", "f3", reduce="min", skip_null=True)], manual=None),
    "mag:magval_tau99_mag_capture": dict(values=[
        J(MSF, ["summary", "c1_x0_exclude", "Genus", "correct"], "genus accuracy tau 95", "0.910", "f3"),
        J(MST, ["summary", "c1_x0_exclude", "Genus", "correct"], "genus accuracy tau 99", "0.895", "f3")],
        manual_resolution={"status": "unchanged", "parts": "1,285 retained references",
                           "why": "an ANI/skani property of the tau-99 search (no gene data)"}),
    "panels:panels.std_gap_share": dict(values=[]),          # filled below from the item's own d0..d15 vars
    "panels:panels.condF1_crude_bounds": dict(values=[
        X("max([abs(r1-o1), abs(r2-o2), abs(r3-o3), abs(r4-o4), abs(r5-o5)])",
          {**{"r%d" % i: {"artifact": PA, "selector": ["cells", {"where": {"regime": "exclude", "completeness": c}},
                                                         "random", "cond_f1_crude"]} for i, c in
              enumerate((1.0, 0.9, 0.7, 0.5, 0.3), 1)},
           **{"o%d" % i: {"artifact": PA, "selector": ["cells", {"where": {"regime": "exclude", "completeness": c}},
                                                         "original", "cond_f1_crude"]} for i, c in
              enumerate((1.0, 0.9, 0.7, 0.5, 0.3), 1)}},
          "random vs original MAG crude cond F1, max |gap| over completeness (printed 'within 0.016')", "0.016", "le"),
        X("max([abs(r1-o1), abs(r2-o2), abs(r3-o3), abs(r4-o4), abs(r5-o5)])",
          {**{"r%d" % i: {"artifact": PA, "selector": ["cells", {"where": {"regime": "exclude", "completeness": c}},
                                                         "kmeans", "cond_f1_crude"]} for i, c in
              enumerate((1.0, 0.9, 0.7, 0.5, 0.3), 1)},
           **{"o%d" % i: {"artifact": PA, "selector": ["cells", {"where": {"regime": "exclude", "completeness": c}},
                                                         "original", "cond_f1_crude"]} for i, c in
              enumerate((1.0, 0.9, 0.7, 0.5, 0.3), 1)}},
          "k-means vs original MAG crude cond F1, max |gap| over completeness (printed 'within 0.010')", "0.010", "le"),
        X("min([o1-max(r1,k1), o2-max(r2,k2), o3-max(r3,k3), o4-max(r4,k4), o5-max(r5,k5)])",
          {**{"%s%d" % (p[0], i): {"artifact": PA, "selector": ["cells", {"where": {"regime": "exclude",
                                                                                       "completeness": c}},
                                                                   p, "cond_f1_crude"]}
              for p in ("original", "random", "kmeans") for i, c in enumerate((1.0, 0.9, 0.7, 0.5, 0.3), 1)}},
          "original ahead at every MAG completeness (min margin over both panels)", "0", "gt")],
        manual_resolution={"status": "partly", "parts": "'original ahead in all eight regions'",
                           "why": "amplicon side: see tab:threepanel genecap_cond_f1 crude rows (R2.3 items)"}),
    "methods:wobj_results": dict(values=[
        X("1-a", {"a": {"artifact": WOBJ, "selector": ["arms", "multi", "real", "internal_share", "median"]}},
          "exchange term share of the objective (median)", "98.9", "pct1"),
        X("c+d", {"c": {"artifact": WOBJ, "selector": ["contrasts", "real-two_tier", "faces.opt.n_var_1e-3", "multi",
                                                       "mean_comparator_weighted"]},
                  "d": {"artifact": WOBJ, "selector": ["contrasts", "real-two_tier", "faces.opt.n_var_1e-3", "multi",
                                                       "mean_delta_weighted"]}},
          "free gene-backed reactions, p_func (real)", "32.7", "f1"),
        J(WOBJ, ["contrasts", "real-uniform", "faces.opt.n_var_1e-3", "multi", "mean_comparator_weighted"],
          "free gene-backed reactions, parsimony (uniform)", "41.8", "f1")],
        manual_resolution={"status": "not recomputed", "parts": "weighted term median 0.11%",
                           "why": "a per-row share computed in WEIGHTED_OBJECTIVE_RESULTS.md section 1 (internal-term "
                                  "decomposition); the wobj rows moved only in 2 of 360 jobs (R4.3), so it is expected "
                                  "unchanged at print"}),
    "panels:disc.objective_test1": dict(values=[
        X("-d/c", {"c": {"artifact": WOBJ, "selector": ["contrasts", "shuffled-two_tier", "faces.opt.n_var_1e-3", "multi",
                                                        "mean_comparator_weighted"]},
                   "d": {"artifact": WOBJ, "selector": ["contrasts", "shuffled-two_tier", "faces.opt.n_var_1e-3", "multi",
                                                        "mean_delta_weighted"]}},
          "shuffled weights cut free reactions by (not printed: 'by more')", None, "pct0"),
        X("(-ds/cs)-(-dr/cr)",
          {"cs": {"artifact": WOBJ, "selector": ["contrasts", "shuffled-two_tier", "faces.opt.n_var_1e-3", "multi",
                                                 "mean_comparator_weighted"]},
           "ds": {"artifact": WOBJ, "selector": ["contrasts", "shuffled-two_tier", "faces.opt.n_var_1e-3", "multi",
                                                 "mean_delta_weighted"]},
           "cr": {"artifact": WOBJ, "selector": ["contrasts", "real-two_tier", "faces.opt.n_var_1e-3", "multi",
                                                 "mean_comparator_weighted"]},
           "dr": {"artifact": WOBJ, "selector": ["contrasts", "real-two_tier", "faces.opt.n_var_1e-3", "multi",
                                                 "mean_delta_weighted"]}},
          "shuffled weights cut by more than p_func (cut difference)", "0", "gt"),
        X("-d/c", {"c": {"artifact": WOBJ, "selector": ["contrasts", "real-two_tier", "faces.opt.n_var_1e-3", "multi",
                                                        "mean_comparator_weighted"]},
                   "d": {"artifact": WOBJ, "selector": ["contrasts", "real-two_tier", "faces.opt.n_var_1e-3", "multi",
                                                        "mean_delta_weighted"]}},
          "p_func weights cut free reactions by", "21", "pct0"),
        X("1-a", {"a": {"artifact": WOBJ, "selector": ["arms", "multi", "real", "internal_share", "median"]}},
          "exchange term share of the objective (median)", "98.9", "pct1")], manual=None),
    # ---------------------------------------------------------------- Stage O / Stage E / calibration (R4H)
    "methods:stageO_results": dict(manual_values=[
        MV("shuffled vs released", "+0.0110 [0.0031, 0.0222]", "+0.0111 [0.0031, 0.0223]",
           "OBJECTIVE_TUNING_RESULTS.md (re-numbered by the R4.3 TF agent from objective_tuning/summary.json); "
           "RUN4_TF_REPORT 'What Run 5 needs' 1"),
        MV("face-min CI, estimate-alone 0.0128, exchange-weight 0.0043, log-cost 0.0039, G effect, n=1 0.0048 / 0.0124, "
           "96 of 193, 0.0164 [0.0100, 0.0236], 69% / 23%", "as printed", "unchanged at print",
           "OBJECTIVE_TUNING_RESULTS.md; RUN4_TF_REPORT 'What Run 5 needs' 1 (checked by hand against the new outputs)")],
        manual=None),
    "val16s:pfunc_calibration": dict(manual_values=[
        MV("calibration slope", "0.78 [0.76, 0.79]", "0.78 [0.76, 0.79]",
           "data/calibration/README.md (regenerated by calibrate_pfunc.py, R4.3)"),
        MV("thresholding at 0.5 changes macro F1 by", "+0.0007", "+0.0007", "data/calibration/README.md key result 5"),
        MV("modal fractional weight", "0.5", "see derived recon_modal_fractional_weight", "W6")], manual=None),
    "val16s:stage_e_estimator": dict(manual_values=[
        MV("isotonic recalibration lowers log-loss by", "0.0411", "0.0412", "PFUNC_ESTIMATORS_RESULTS.md section 1 "
           "(re-numbered by the R4.3 TF agent)"),
        MV("within-union macro dF1", "+0.0085 [+0.0072, +0.0096]", "+0.0085 [+0.0072, +0.0097]",
           "PFUNC_ESTIMATORS_RESULTS.md; make_deck.py PFE_WITHIN_UNION_F1"),
        MV("0.0889 -> 0.0578 (0.0311 [0.0297, 0.0322]); 0.209 -> 0.127; 0.543 -> 0.816; 0.805; 0.0408 [0.0384, 0.0424]; "
           "0.739; 0.067 [0.062, 0.072]; +0.032; +0.0249 [+0.0228, +0.0265]; 13%; 57,021 / 49,048 cells",
           "as printed", "unchanged at print", "PFUNC_ESTIMATORS_RESULTS.md; RUN4_TF_REPORT 'What Run 5 needs' 1")],
        manual=None),
    # ---------------------------------------------------------------- explicit resolutions
    "front:abs-objective-posthoc": dict(manual_resolution={
        "status": "no number", "why": "qualitative sentence; its numbers are under methods:stageO_results"}),
    "methods:stageO_design_and_split": dict(manual_resolution={
        "status": "unchanged", "parts": "121 / 80 groups",
        "why": "design: the wobj population has the same 54,243 candidates and 360 keys (R4.3 population check), so the "
               "hash split is unchanged"}),
    "panels:disc.objective_test2": dict(manual_resolution={
        "status": "unchanged", "parts": "121 of 201 / 80 groups", "why": "design (same 360 keys)"}),
    "panels:lim.reconstruction_tests": dict(manual_resolution={
        "status": "unchanged", "parts": "300/60 amplicons, 121 of 201 groups, slope 0.78, ~7,000 genomes",
        "why": "design restatement; the calibration slope is unchanged at print (val16s:pfunc_calibration)"}),
    "front:tab:magval": dict(manual_resolution={
        "status": "unchanged", "parts": "9,997 panel organisms", "why": "panel size (selection_10k), not re-scored"}),
    "front:fig:magcomp": dict(manual_resolution={"status": "figure", "why": "F14, re-rendered in R5.2"}),
    "front:fig:magcross": dict(manual_resolution={"status": "figure", "why": "F15, re-rendered in R5.2"}),
    "mag:fig_magcomp_magcross": dict(manual_resolution={"status": "figure", "why": "F14/F15, re-rendered in R5.2"}),
    "front:fig:genecap": dict(manual_resolution={"status": "figure", "why": "F4b, re-rendered in R5.2"}),
    "picrust:AS-24": dict(manual_resolution={"status": "figure", "why": "F13, re-rendered in R5.2 (suptitle n computed)"}),
    "mag:magval_no_crossover": dict(manual_resolution={
        "status": "unchanged", "parts": "37 kb",
        "why": "median retained source sequence at 1% completeness: a property of the synthetic MAGs, which were not "
               "re-synthesised"}),
    "panels:panels.taxonomy_bootstrap_CIs": dict(manual_resolution={
        "status": "note", "why": "region-order assumption of the tab:threepanel genus block; taxonomy point estimates "
                                 "are unchanged, CIs move by bootstrap-stream noise (R2.3)"}),
    "panels:panels.surviving_uncovered": dict(manual=None),
    "panels:lim.depth20_vs_100": dict(manual_resolution={
        "status": "see item", "parts": "legacy genus within 0.01, recall +0.7..+1.3, F1 -0.1..-1.2 points",
        "why": "values under val16s:legacy_genecap_d20_vs_d100 (R1.4)"}),
    "picrust:FR-01": dict(manual_resolution={
        "status": "unchanged (Path A)", "parts": "2,059 held-out genomes", "why": "fold plan kept (D1 Path A)"}),
    "picrust:FR-04": dict(manual_resolution={
        "status": "unchanged (Path A)", "parts": "1,653 of 8,872 (18.6%)",
        "why": "panel_eligible.json unchanged (sha256 and inode checked in R3.1)"}),
    "panels:lim.database_sizes": dict(manual_resolution={
        "status": "unchanged", "why": "database facts (~315,000 genomes, 459,301 16S, 27,870 PICRUSt2 genomes)"}),
    "panels:lim.holdout_power": dict(manual_resolution={
        "status": "unchanged (Path A)", "parts": "18.6% (1,653)", "why": "panel_eligible.json unchanged"}),
    "panels:disc.phylum_span": dict(manual=None),
    "panels:lim.deployment_truncation": dict(manual_resolution={
        "status": "unchanged", "parts": "25 ASVs; 2.034 / 19 depth-20 comparator",
        "why": "deployment builds were not re-run (D4); deployment_d100/RESULTS.md"}),
    "val16s:backend_concordance": dict(manual_resolution={
        "status": "historical (D7)", "why": "June check; wording only ('on the cache of the time')"}),
    "val16s:genecap_denominators": dict(manual=None, manual_resolution={
        "status": "note", "why": "the denominator itself moves: scoreable amplicons 67,038 -> 67,713 "
                                 "(truth_empty 953 -> 278)"}),
}

# claims whose meaning (not just digits) changes; the manuscript agent must rewrite these sentences
MEANING = {
    "picrust:PH-04": "prFBA-vs-nn1 ATP advantage is no longer significant: -0.146 media (p=0.017, q=0.028) -> -0.030 "
                     "(p=0.41, q=0.51); all-genome p 0.28 -> 0.46 with the sign flipped (+0.06). The 'phenotype favors "
                     "prFBA against nn1' sentence (L635-637) and L1186 must be rewritten; the printed EC-F1 p=0.64 has no "
                     "recipe (see derived pheno_nn1_ecf1_p)",
    "panels:disc.atp_vs_nn1": "as picrust:PH-04 (L1186): no longer significant",
    "picrust:PH-02": "McNemar: PICRUSt2 ahead on FOUR media at q<=0.05, not two (Ac.NO q 0.029, Ac.O2 0.040, Etho.O2 "
                     "0.046, For.O2 0.046); 'two of the 22 testable media (Ac.NO and Etho.O2, q=0.002)' at L631, L639, "
                     "L1185 and the abstract (L79) must be rewritten",
    "front:abs-phenotype-readout": "abstract: 755 -> 759 genomes (747 -> 750 V4), 0.13 -> 0.07 media (p 0.49 -> 0.79), "
                                   "'two of 22 media' -> four; the phenotype sample was re-drawn (442 of 1,000 kept)",
    "picrust:PH-01": "sample re-drawn by the stratified draw (442 of 1,000 ATP genomes kept); disclose why",
    "picrust:PH-03": "rho -0.14 (p 2e-4) -> -0.19 (p 1.5e-7): the per-genome link is stronger",
    "front:abs-phenotype-controls": "the A-C13 counts are now those of a re-run on a re-drawn sample: 43,140 -> 43,174 "
                                    "rows, 0 of 8,614 -> 0 of 8,660 gated, 114 -> 107 reversals; 'reproduces the "
                                    "pre-registered pass exactly (43,140 of 43,140)' is a statement about the old sample",
    "picrust:AM-02": "as front:abs-phenotype-controls (L667-679): the reproduction sentences describe the 10-29 sample",
    "picrust:AS-05": "L517 cache-repair disclosure must become a past-tense correction note (953 / 1,060 no longer "
                     "describe the reported state)",
    "val16s:genecap_denominators": "L397/L419 denominators move: score-zero 11.1% -> 9.9%, pred_empty 1.6% -> 0.4% "
                                   "(0.987) and 1.4% -> 0.4% (0.995); scoreable n 67,038 -> 67,713",
    "panels:panels.std_gap_share": "check 'interval excludes zero in all 16' against the new three-panel CIs "
                                   "(21 excludes_zero verdicts flipped in R2.3: control/r2/r23/three_panel_ci_flips.json)",
    "panels:panels.taxonomy_bootstrap_CIs": "21 excludes_zero verdicts flipped in R2.3 (5 in _d100, 16 in _d100_legacy), "
                                            "all borderline; none of the 21 cells is cited by any extractor (checked "
                                            "against every three_panel selector), but blanket sentences over cells "
                                            "(L997-999, L1094-1137) should be re-read against "
                                            "control/r2/r23/three_panel_ci_flips.json",
    "val16s:emptied_vs_uncovered_panels": "L742: original FL emptied-record rate 4.7% -> 4.6% (N 7,215 -> 7,268), pooled "
                                          "8.3% -> 8.2%; 'the share of full-length records that survive exclusion yet go "
                                          "uncovered stays near 3% on both (3.4% and 2.8%; 3.1% and 2.4% pooled)' becomes "
                                          "near 2% (2.4% and 2.3%; 1.8% and 1.8% pooled)",
    "panels:panels.surviving_uncovered": "L1057-1060: the surviving-uncovered rates fall well below 3% because the repair "
                                         "turned most gene-less (pred_empty) survivors into answers: FL 3.4/2.8/2.3% -> "
                                         "2.4/2.3/1.8%, pooled 3.1/2.4/1.8% -> 1.8/1.8/1.3%, FL strata 1.5-3.7% -> "
                                         "0.0-3.0% (the original 50-99 stratum is now 0 of 274), '29- to 43-fold below' "
                                         "-> 42- to 55-fold, spreads 1.1 (1.3 pooled) -> 0.6 (0.5) points; the CMH test "
                                         "original vs random is no longer significant (p 0.004 -> 0.117; vs k-means 0.15 -> "
                                         "0.99). Rewrite any 'original panel significantly higher' sentence",
    "panels:panels.surviving_pred_empty": "L1064: non-strict uncovered survivors 107/76/62 -> 36/42/31: the original panel "
                                          "no longer has the most (k-means 42 > original 36)",
    "panels:panels.tab_emptied_uncovered": "tab:emptied: every surviving-stratum uncovered rate falls (original 0-19 0.037 "
                                           "-> 0.026; original 50-99 0.015 -> 0.000)",
    "mag:twopanel_mag_abstention_strata": "L1085-1087: without the >=20-member stratum k-means abstains at 0.686 against the "
                                          "original's 0.682 (was 0.687 vs 0.688): the stratum still accounts for the "
                                          "overall gap (0.657 vs 0.682), but 'the same as the original' becomes 'slightly "
                                          "above'",
    "picrust:PH-06": "L641-646: growth log-ratio prFBA 0.008 -> 0.020 vs nn1 0.008 -> 0.011 (PICRUSt2 -0.006 -> -0.003); "
                     "'growth log-ratio 0.008 against 0.008' no longer reads as identical; essentiality 0.065/0.068/0.069 "
                     "-> 0.059/0.063/0.063 keeps the ordering; burden 80-84 vs 82.5 -> 78-82 vs 79.5",
    "panels:panels.condF1_crude_bounds": "L1018-1019: 'random within 0.016 at every MAG completeness' no longer holds "
                                         "(max gap 0.0167 -> print 0.017)",
    "val16s:depth_probe_quality": "L728-731: the probe's recovered amplicons no longer score BELOW the depth-20 answered "
                                  "set in EC space: precision 0.747 -> 0.791 (vs 0.790), F1 0.783 -> 0.829 (vs 0.798), "
                                  "recall 0.854 -> 0.904 (vs 0.835); the sentence 'score below the depth-20 answered set ... "
                                  "which predicted conditional figures roughly unchanged or slightly lower' must be rewritten "
                                  "(genus/family anchor accuracy 60.1% / 80.2% unchanged; the 0.651 / 0.800 comparators "
                                  "were not recomputed)",
    "front:abs-depth-ec": "L67: probe-recovered answers are now at or above the answered set in EC space (F1 0.829 vs "
                          "0.798, precision 0.791 vs 0.790), not 'comparable' with lower precision",
    "panels:lim.depth_probe_quality": "L1204 'at or below the answered set': in EC space the recovered calls are now above "
                                      "it (EC F1 0.829 vs 0.798); the genus/family anchor figures are unchanged",
    "panels:disc.margin_reversal": "magnitude (not direction): the as-shipped population deficit shrinks from -0.016 to "
                                   "-0.006 (prFBA 0.713 -> 0.723 vs PICRUSt2 0.729) and the depth-20 one from -0.129 to "
                                   "-0.120; the per-region ranking is unchanged (prFBA ahead for FL, V1-V2, V1-V3 only); "
                                   "'widens the population deficit from -0.016' (L1187, L1216) still reads correctly with "
                                   "-0.006",
    "picrust:AS-07": "",
}


def main():
    spec = json.load(open(SPEC))
    items = spec["items"]
    for iid, add in ADD.items():
        it = items[iid]
        have = {v.get("label") for v in it.get("values", [])}
        for v in add.get("values", []):
            if v["label"] not in have:
                it.setdefault("values", []).append(v)
                have.add(v["label"])
        if "manual" in add:
            if add["manual"] is None:
                if it.get("manual"):
                    it["manual_resolved_r5a"] = it["manual"]
                it.pop("manual", None)
            else:
                it["manual"] = add["manual"]
        for k in ("manual_values", "manual_resolution"):
            if k in add:
                it[k] = add[k]
    # std_gap_share: 'interval excludes zero' in all 16 cells, from the item's own d0..d15 selectors
    it = items["panels:panels.std_gap_share"]
    base = next(v for v in it["values"] if v.get("expr", "").startswith("min("))
    ev = {}
    for k, v in base["vars"].items():
        if k.startswith("d"):
            sel = list(v["selector"])
            assert sel[-1] == "est", sel
            ev["e" + k[1:]] = {"artifact": v["artifact"], "selector": sel[:-1] + ["excludes_zero"]}
    lab = "standardized intervals excluding zero (of 16)"
    if lab not in {v.get("label") for v in it["values"]}:
        it["values"].append({"label": lab, "printed": "16", "fmt": "int",
                             "expr": "sum([%s])" % ", ".join("e%d" % i for i in range(len(ev))), "vars": ev})
    if it.get("manual"):
        it["manual_resolved_r5a"] = it.pop("manual")
    for iid, note in MEANING.items():
        if note:
            items[iid]["meaning_change"] = note
    spec["_meta"]["amended_by_r5a"] = ("$W/tools/r5a_spec_extend.py (2026-11-16): manual parts turned into W6 derived "
                                       "values, selectors and expressions, hand-carried manual_values with a source, or an "
                                       "explicit manual_resolution; curated meaning_change notes")
    spec["_meta"]["n_values"] = sum(len(v.get("values", [])) for v in items.values())
    spec["_meta"]["n_items_without_values"] = sum(not v.get("values") for v in items.values())
    with open(SPEC + ".tmp", "w") as fh:
        json.dump(spec, fh, indent=1)
    os.replace(SPEC + ".tmp", SPEC)
    print("items", len(items), "values", spec["_meta"]["n_values"],
          "still manual", sum(bool(v.get("manual")) for v in items.values()),
          "with manual_resolution", sum(bool(v.get("manual_resolution")) for v in items.values()),
          "with manual_values", sum(bool(v.get("manual_values")) for v in items.values()),
          "meaning_change", sum(bool(v.get("meaning_change")) for v in items.values()))


if __name__ == "__main__":
    main()
