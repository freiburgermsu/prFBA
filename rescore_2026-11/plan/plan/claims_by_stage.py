"""Join the six manuscript inventories to the rescore stages by source-artefact path.
Read-only; writes claims_by_stage.{json,md} next to this script."""
import json, os, re, collections
D = os.path.dirname(os.path.abspath(__file__))
INV = os.path.dirname(D)
FILES = ["inventory_front.json", "inventory_methods.json", "inventory_val16s.json",
         "inventory_picrust.json", "inventory_mag.json", "inventory_panels.json"]
# (stage, regex over source_artifact + producer); the LAST stage in DAG order that matches wins,
# because a claim is final only after its latest input is regenerated.
ORDER = ["R1.3U", "R1.3L", "R1.4", "R1.5", "R1.K", "R2.1", "R2.2", "R2.3", "R2.4", "R2.5",
         "R3", "R3F", "R3C2", "R4P", "R4S2", "R4H", "R5A"]
RULES = [
    ("R1.3U", r"region_validation_d100/data/gene_capture_(cells|by_region)|region_validation_d100/data/genecap_repaired_impact|headline_macro_recall|gene_capture_cells\.csv"),
    ("R1.3L", r"region_validation_d100_legacy|_d100\{,_legacy\}|legacy/data|d100\{,_legacy\}"),
    ("R1.4", r"region_validation/data/(gene_capture|taxacc|selection)"),
    ("R1.5", r"tau_compare|tau0995"),
    ("R1.K", r"deployment_d100|bvbrc_rebuild_summary|stats_d100"),
    ("R2.1", r"pfunc_e3_gate"),
    ("R2.2", r"mag_validation/out|mag_cells|mag_summary|summary_n500|af_gate|out_batches|crossover|paired_comparison|paired_tau99|RESULTS_v2|RESULTS_tau99"),
    ("R2.3", r"three_panel_comparison|emptied_strata|three_panel_mag_panelarm|mag_panel_comparison|panel_comparison|mechanism_depth|UNION_THROUGHOUT_CHANGES.md item 1\.72"),
    ("R2.4", r"region_validation_d100/picrust2_comparison"),
    ("R2.5", r"region_validation/picrust2_comparison|annotation_ceiling|depth_probe|picrust2_comparison/README"),
    ("R3", r"shared_reference/data/(summary_genus|cells|controls|prfba|fold_plan|panel_eligible)|shared_reference/data\b|shared_reference/data \(|summary_genus"),
    ("R3F", r"data_family|summary_family"),
    ("R3C2", r"data_c2|summary_incself"),
    ("R4P", r"phenotype|C12_SCOPE|c12fix|c13_|EGC_C12"),
    ("R4S2", r"phenotype_gf_c13/phenotype_summary"),
    ("R4H", r"calibration|pfunc_estimators|PFUNC_ESTIMATORS|weighted_objective|WEIGHTED_OBJECTIVE|objective_tuning|OBJECTIVE_TUNING|wrapper_equivalence"),
]
# per-key overrides where the path rule is wrong or the claim needs new code (R5A = ad hoc
# derivation script to be written in run 5; R0 = no stage, depends on a decision)
OVERRIDE = {
    "sharedref_tieset_coherence": "R5A", "recon_modal_fractional_weight": "R5A",
    "panels.welch_tests": "R5A", "panels.mag_singleton_abstention": "R5A",
    "panels.mag_multiplicity": "R5A", "magval_scope_claims": "R5A", "magval_tau99_paired": "R5A",
    "magval_identity_bins": "R5A", "AS-10": "R5A", "backend_concordance": "D7",
    "SR-01": "D1", "FR-01": "D1", "FR-04": "D1",
    "abs-phenotype-readout": "R4P", "abs-phenotype-controls": "R4S2",
    "abs-objective-posthoc": "R4H", "disc.objective_test1": "R4H", "disc.objective_test2": "R4H",
    "lim.reconstruction_tests": "R4H", "stage_e_estimator": "R4H", "pfunc_calibration": "R4H",
    "abs-sharedref-size": "D1",
    "abs-dialect-density": "R3", "abs-region-penalty": "R3",
    "e3_gate_brier_auroc": "R2.1", "e3_gate_recalibration_share": "R2.1", "e3_production_gate": "R2.1",
    "fig:picrust2": "R2.4", "AS-24": "R2.4", "tab:picrust2": "R2.5",
    "abs-depth-ec": "R2.5", "conditional_ec_d20_vs_d100": "R2.5", "population_ec_and_projections": "R2.5",
    "lim.depth20_vs_100": "R2.5", "lim.depth_probe_quality": "R2.5", "depth_probe_quality": "R2.5",
    "legacy_genecap_d20_vs_d100": "R1.4", "controlled_forced_recovery": "R3",
    "abs-mag-vs-16s": "R2.3", "abstract_mag_margin": "R2.3", "tab:magvs16s": "R2.3", "tab_magvs16s": "R2.3",
    "tab:magval": "R2.3", "fig:magcomp": "R2.3", "fig:magcross": "R2.3", "fig_magcomp_magcross": "R2.3",
    "magval_completeness_null": "R2.3", "magval_paired_main": "R2.3", "magval_no_crossover": "R2.3",
    "magval_vs_every_region": "R2.3", "magval_population_view": "R2.3", "magval_af_gate_counterfactual": "R2.3",
    "methods_union_size_amplicon": "R1.3U", "magval_selfinclusion_ceiling": "R2.2",
    "magval_tau99_mag_capture": "R2.2", "twopanel_mag_depth500_original": "R2.2", "panels.mag_depth500_rerun": "R2.2",
    "emptied_vs_uncovered_panels": "R2.3", "disc.selection_size_controlled": "R3",
    "controlled_selection_sizes": "R3", "sharedref_selection_sizes": "R3", "sharedref_ngenomes_vs_identity": "R3",
    "XR-03": "R3", "disc.atp_vs_picrust2": "R4P", "disc.atp_vs_nn1": "R4P", "lim.phenotype_scope": "R4P",
    "msprob_released_equivalence": "R4H", "stageO_design_and_split": "R4H", "stageO_results": "R4H",
    "stageO_wrapper_equivalence": "R4H", "wobj_results": "R4H",
    "deployment_cross_phylum": "R1.K", "e3_production_unresolved_share": "R1.K", "deployment_e3_unresolved_share": "R1.K",
    "disc.selection_size_deployment": "R1.K", "disc.phylum_span": "R1.K", "lim.deployment_truncation": "R1.K",
    "AS-05": "R2.4",
    "SR-02": "R3", "SR-03": "R3", "SR-04": "R3", "SR-05": "R3", "SR-06": "R3", "SR-07": "R3",
    "SR-08": "R3", "SR-09": "R3", "SR-10": "R3", "SR-11": "R3", "SR-12": "R3",
    "AM-02": "R4S2", "AM-03": "R4S2", "AS-19": "R2.5", "lim.database_sizes": "D1",
    "lim.holdout_power": "D1", "disc.vocab_breadth": "R3", "disc.refusal_credit": "R4P",
}
def stage_of(it):
    if it["key"] in OVERRIDE:
        return OVERRIDE[it["key"]]
    if it["depends_on_bvbrc_data"] == "no":
        return "NO"
    text = it["source_artifact"] + " || " + it["producer"]
    hit = None
    for st, rx in RULES:
        if re.search(rx, text):
            if hit is None or ORDER.index(st) > ORDER.index(hit):
                hit = st
    return hit or "UNMATCHED"
out = collections.defaultdict(list)
for f in FILES:
    d = json.load(open(os.path.join(INV, f)))
    items = d["items"] if isinstance(d, dict) else d
    for it in items:
        st = stage_of(it)
        if it["depends_on_bvbrc_data"] == "no" and it["key"] not in OVERRIDE:
            st = "NO"
        out[st].append({"inv": f.replace("inventory_", "").replace(".json", ""), "key": it["key"],
                        "tex_lines": it["tex_lines"], "dep": it["depends_on_bvbrc_data"],
                        "claim": it["claim"]})
json.dump(out, open(os.path.join(D, "claims_by_stage.json"), "w"), indent=1)
with open(os.path.join(D, "claims_by_stage.md"), "w") as fh:
    for st in sorted(out, key=lambda s: (ORDER + ["D1", "D7", "NO", "UNMATCHED"]).index(s)
                     if s in ORDER + ["D1", "D7", "NO", "UNMATCHED"] else 99):
        fh.write(f"\n### {st} ({len(out[st])} items)\n\n")
        for it in out[st]:
            fh.write(f"- `{it['inv']}:{it['key']}` L{it['tex_lines']}: {it['claim'][:160]}\n")
for st in out:
    print(st, len(out[st]))
