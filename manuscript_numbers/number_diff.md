# Manuscript number diff (W7)

Generated 2026-11-17T10:41:42-0600 by `manuscript_numbers/collect.py`; old = tag `rescore-2026-11-start` (`cc8014f`; git-ignored files from `$W/pre_repair` or, failing that, the live tree), new = live tree at `5dcfc0f`.

| count | value |
|---|---|
| n_items | 178 |
| n_items_spec | 142 |
| n_items_partial | 28 |
| n_items_manual | 8 |
| n_values | 1580 |
| n_old_resolved | 1580 |
| n_new_resolved | 1580 |
| n_pending | 0 |
| n_errors | 0 |
| n_changed | 1183 |
| n_changed_at_print_precision | 716 |
| n_old_matches_printed | 1556 |
| n_old_matches_printed_at_rounding_boundary | 48 |
| n_old_mismatch_printed | 9 |
| n_printed_found_in_tex | 1546 |
| n_printed_not_found_in_tex | 0 |
| n_items_changed_at_print | 117 |
| n_items_meaning_change_curated | 23 |
| n_items_meaning_flag_auto | 7 |
| n_items_meaning_change_any | 25 |
| n_items_with_manual_values | 3 |
| n_items_manual_unresolved | 0 |
| manual_resolution_status | {'no number': 1, 'not recomputed': 6, 'unchanged': 9, 'figure': 5, 'partly': 1, 'note': 2, 'see item': 1, 'historical': 5, 'unchanged (Path A)': 3, 'historical (D7)': 1} |
| n_new_from_precise_source | 263 |
| n_new_precise_disagrees_with_stored | 0 |
| n_new_print_indeterminate | 3 |
| n_new_fmt_in_tex | 1544 |
| n_new_fmt_not_in_tex | 3 |
| n_new_fmt_not_in_tex_with_exception | 3 |
| n_new_fmt_not_in_tex_unexplained | 0 |
| mismatch_reasons | {'rounding': 5, 'differs': 4} |

## Claims whose meaning changes (rewrite, not just refresh)

Curated notes (`meaning_change` in claims_spec.json) and automatic flags (sign flip, a p or q crossing 0.05 / 0.01, a categorical or yes/no value changing, a CI's 'excludes zero' flipping, a printed bound no longer holding). Automatic flags are candidates: read the row before rewriting.

| item | tex lines | curated note | automatic flags (value: flag) |
|---|---|---|---|
| val16s:emptied_vs_uncovered_panels | 742 | L742: original FL emptied-record rate 4.7% -> 4.6% (N 7,215 -> 7,268), pooled 8.3% -> 8.2%; 'the share of full-length records that survive exclusion yet go uncovered stays near 3% on both (3.4% and 2.8%; 3.1% and 2.4% pooled)' becomes near 2% (2.4% and 2.3%; 1.8% and 1.8% pooled) |  |
| mag:twopanel_mag_abstention_strata | 1085-1087 | L1085-1087: without the >=20-member stratum k-means abstains at 0.686 against the original's 0.682 (was 0.687 vs 0.688): the stratum still accounts for the overall gap (0.657 vs 0.682), but 'the same as the original' becomes 'slightly above' |  |
| panels:panels.condF1_crude_bounds | 1018-1019 | L1018-1019: 'random within 0.016 at every MAG completeness' no longer holds (max gap 0.0167 -> print 0.017) |  |
| panels:panels.surviving_uncovered | 1057-1060 | L1057-1060: the surviving-uncovered rates fall well below 3% because the repair turned most gene-less (pred_empty) survivors into answers: FL 3.4/2.8/2.3% -> 2.4/2.3/1.8%, pooled 3.1/2.4/1.8% -> 1.8/1.8/1.3%, FL strata 1.5-3.7% -> 0.0-3.0% (the original 50-99 stratum is now 0 of 274), '29- to 43-fold below' -> 42- to 55-fold, spreads 1.1 (1.3 pooled) -> 0.6 (0.5) points; the CMH test original vs random is no longer significant (p 0.004 -> 0.117; vs k-means 0.15 -> 0.99). Rewrite any 'original panel significantly higher' sentence | CMH original vs random p: significance crosses 0.05 (0.00407 -> 0.117); unstratified Fisher p original vs random: significance crosses 0.01 (0.0003 -> 0.0229) |
| panels:panels.surviving_pred_empty | 1064 | L1064: non-strict uncovered survivors 107/76/62 -> 36/42/31: the original panel no longer has the most (k-means 42 > original 36) |  |
| panels:panels.tab_emptied_uncovered | 1066-1081 | tab:emptied: every surviving-stratum uncovered rate falls (original 0-19 0.037 -> 0.026; original 50-99 0.015 -> 0.000) |  |
| panels:panels.std_gap_share | 1095-1096 | check 'interval excludes zero in all 16' against the new three-panel CIs (21 excludes_zero verdicts flipped in R2.3: control/r2/r23/three_panel_ci_flips.json) |  |
| panels:panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | 21 excludes_zero verdicts flipped in R2.3 (5 in _d100, 16 in _d100_legacy), all borderline; none of the 21 cells is cited by any extractor (checked against every three_panel selector), but blanket sentences over cells (L997-999, L1094-1137) should be re-read against control/r2/r23/three_panel_ci_flips.json |  |
| front:abs-depth-ec | 67 | L67: probe-recovered answers are now at or above the answered set in EC space (F1 0.829 vs 0.798, precision 0.791 vs 0.790), not 'comparable' with lower precision |  |
| val16s:depth_probe_quality | 67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204 | L728-731: the probe's recovered amplicons no longer score BELOW the depth-20 answered set in EC space: precision 0.747 -> 0.791 (vs 0.790), F1 0.783 -> 0.829 (vs 0.798), recall 0.854 -> 0.904 (vs 0.835); the sentence 'score below the depth-20 answered set ... which predicted conditional figures roughly unchanged or slightly lower' must be rewritten (genus/family anchor accuracy 60.1% / 80.2% unchanged; the 0.651 / 0.800 comparators were not recomputed) |  |
| panels:lim.depth_probe_quality | 1204 | L1204 'at or below the answered set': in EC space the recovered calls are now above it (EC F1 0.829 vs 0.798); the genus/family anchor figures are unchanged |  |
| picrust:AS-05 | 517 | L517 cache-repair disclosure must become a past-tense correction note (953 / 1,060 no longer describe the reported state) |  |
| picrust:AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') |  | legacy reducer minus PICRUSt2, as-shipped population EC F1 (depth 100): sign flip |
| panels:disc.margin_reversal | 1184, 1187-1188, 1215-1216 | magnitude (not direction): the as-shipped population deficit shrinks from -0.016 to -0.006 (prFBA 0.713 -> 0.723 vs PICRUSt2 0.729) and the depth-20 one from -0.129 to -0.120; the per-region ranking is unchanged (prFBA ahead for FL, V1-V2, V1-V3 only); 'widens the population deficit from -0.016' (L1187, L1216) still reads correctly with -0.006 | legacy reducer minus PICRUSt2, as-shipped population EC F1 (depth 100): sign flip |
| front:abs-phenotype-readout | 79 | abstract: 755 -> 759 genomes (747 -> 750 V4), 0.13 -> 0.07 media (p 0.49 -> 0.79), 'two of 22 media' -> four; the phenotype sample was re-drawn (442 of 1,000 kept) |  |
| picrust:PH-01 | 79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755) | sample re-drawn by the stratified draw (442 of 1,000 ATP genomes kept); disclose why |  |
| picrust:PH-02 | 79 ('two of 22 media'); 631 (Ac.NO and Etho.O2, q=0.002); 639; 1185 | McNemar: PICRUSt2 ahead on FOUR media at q<=0.05, not two (Ac.NO q 0.029, Ac.O2 0.040, Etho.O2 0.046, For.O2 0.046); 'two of the 22 testable media (Ac.NO and Etho.O2, q=0.002)' at L631, L639, L1185 and the abstract (L79) must be rewritten | Ac.NO q: significance crosses 0.01 (0.00243 -> 0.0289); Etho.O2 q: significance crosses 0.01 (0.00243 -> 0.0461); media favouring PICRUSt2: categorical change Ac.NO,Etho.O2 -> Ac.NO,Ac.O2,Etho.O2,For.O2; largest q of those media: significance crosses 0.01 (0.00243 -> 0.0461); smallest q of the media with q <= 0.05: significance crosses 0.01 (0.00243 -> 0.0289) |
| picrust:PH-03 | 634 (correlation -0.14); 1190 (Spearman rho = -0.14, p = 2x10^-4) | rho -0.14 (p 2e-4) -> -0.19 (p 1.5e-7): the per-genome link is stronger |  |
| picrust:PH-04 | 635-637; 1186 | prFBA-vs-nn1 ATP advantage is no longer significant: -0.146 media (p=0.017, q=0.028) -> -0.030 (p=0.41, q=0.51); all-genome p 0.28 -> 0.46 with the sign flipped (+0.06). The 'phenotype favors prFBA against nn1' sentence (L635-637) and L1186 must be rewritten; the printed EC-F1 p=0.64 has no recipe (see derived pheno_nn1_ecf1_p) | p: significance crosses 0.05 (0.0168 -> 0.407); q: significance crosses 0.05 (0.0285 -> 0.509) |
| panels:disc.atp_vs_picrust2 | 1185 |  | Ac.NO q: significance crosses 0.01 (0.00243 -> 0.0289); Etho.O2 q: significance crosses 0.01 (0.00243 -> 0.0461); smallest q of the media with q <= 0.05: significance crosses 0.01 (0.00243 -> 0.0289) |
| panels:disc.atp_vs_nn1 | 1186 | as picrust:PH-04 (L1186): no longer significant | p: significance crosses 0.05 (0.0168 -> 0.407); q: significance crosses 0.05 (0.0285 -> 0.509) |
| front:abs-phenotype-controls | 80 | the A-C13 counts are now those of a re-run on a re-drawn sample: 43,140 -> 43,174 rows, 0 of 8,614 -> 0 of 8,660 gated, 114 -> 107 reversals; 'reproduces the pre-registered pass exactly (43,140 of 43,140)' is a statement about the old sample |  |
| picrust:PH-06 | 641-646; 1222 (n=200) | L641-646: growth log-ratio prFBA 0.008 -> 0.020 vs nn1 0.008 -> 0.011 (PICRUSt2 -0.006 -> -0.003); 'growth log-ratio 0.008 against 0.008' no longer reads as identical; essentiality 0.065/0.068/0.069 -> 0.059/0.063/0.063 keeps the ordering; burden 80-84 vs 82.5 -> 78-82 vs 79.5 |  |
| picrust:AM-02 | 80 ('43,140 of 43,140 rows ... 0 of 8,614 rows, against 49 ... valid under the amendments'); 667-673; 676-679 (seven-reaction draft correction: 865 -> 95 of 1,181 drafts; medians move <= 0.004; burden <= 1 reaction, 1.5 as shipped; 138 of 200 truth drafts make ATP from nothing; 43,140 ATP-pass rows with no empty-medium ATP); 1193 (ATP pass reproduces every vector) | as front:abs-phenotype-controls (L667-679): the reproduction sentences describe the 10-29 sample |  |
| val16s:genecap_denominators | 397 (gene-capture part: 11.1%, 9.6%, 1.6%), 419 (9.6% -> 6.0%; 1.6% -> 1.4%) | L397/L419 denominators move: score-zero 11.1% -> 9.9%, pred_empty 1.6% -> 0.4% (0.987) and 1.4% -> 0.4% (0.995); scoreable n 67,038 -> 67,713 |  |

## New values against the current manuscript

Each new value, formatted as the manuscript prints it, looked for at the claim's lines mapped from `0214a14` onto the current tex (±2 lines). 1544 found; 3 not found, 3 of them with a recorded reason (`tex_exception`). Values from a finer source than the stored one: 263; still on a rounding boundary (printed by the half-up convention): 3.

| item | tex lines | value | new | where | reason |
|---|---|---|---|---|---|
| panels.surviving_uncovered | 1057-1060 | FL strata uncovered (low) | 0.0 | absent | printed in words: 'down to none of the original panel's 274 records in the 50--99 stratum' |
| AS-12 | 521; 576 (fig C) | wins per loss (max) | 7 | elsewhere in the file | printed in words: 'between three and seven amplicons' |
| abs-phenotype-controls | 80 | ATP-pass rows (A-C13 configuration, re-run) | 43,174 | elsewhere in the file | the abstract keeps the pre-registered pass's 43,140 on purpose (the row-for-row reproduction on the original sample); the re-run's 43,174 is printed in sec:sharedref |

Values still on a rounding boundary (no finer source; printed half-up):

| item | value | new | printed as |
|---|---|---|---|
| tab:threepanel-genecap | V7-V9 delta cond F1 CI lo | -0.0345 | -0.035 |
| panels.tab_threepanel_genecap_block | V7-V9 delta cond F1 CI lo | -0.0345 | -0.035 |
| panels.taxonomy_bootstrap_CIs | genus cond acc V3-V4 CI lo | -0.0385 | -0.039 |

## R4H

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| wobj_results | 307-311 | beat shuffled | 0.0120 | 0.0120 | 0.0121 | **CHANGED** |
| wobj_results | 307-311 | shuffled free reactions | 27.9 | 27.9 | 27.8 | **CHANGED** |
| stageO_results | 319-324 | beat shuffled | 0.0282 | 0.0282 | 0.0283 | **CHANGED** |
| disc.objective_test1 | 1172, 1225 | shuffled weights cut by more than p_func (cut difference) | 0 | 0.11601208467416319 | 0.11672154599548298 | **CHANGED** |
| wobj_results | 307-311 | p_func vs two-tier face-mid delta | 0.0061 | 0.0061 | 0.0061 | changed below print precision |
| wobj_results | 307-311 |   CI bound (near) | 0.0025 | 0.0025 | 0.0025 | changed below print precision |
| wobj_results | 307-311 |   CI bound (far) | 0.0097 | 0.0097 | 0.0097 | changed below print precision |
| wobj_results | 307-311 | two-tier comparator | 0.130 | 0.130 | 0.130 | changed below print precision |
| wobj_results | 307-311 | share of two-tier (4.7%) | 4.7 | 4.7 | 4.7 | changed below print precision |
| wobj_results | 307-311 | face minimum delta | -0.0025 | -0.0025 | -0.0025 | changed below print precision |
| wobj_results | 307-311 |   CI low | -0.0053 | -0.0053 | -0.0053 | changed below print precision |
| wobj_results | 307-311 |   CI high | +0.0011 | +0.0011 | +0.0011 | changed below print precision |
| wobj_results | 307-311 |   CI (near) | 0.0069 | 0.0069 | 0.0069 | changed below print precision |
| wobj_results | 307-311 |   CI (far) | 0.0179 | 0.0179 | 0.0179 | changed below print precision |
| wobj_results | 307-311 | shuffled vs two-tier | +0.0060 | +0.0060 | +0.0060 | changed below print precision |
| wobj_results | 307-311 | oracle vs two-tier | 0.0249 | 0.0249 | 0.0249 | changed below print precision |
| wobj_results | 307-311 |   CI (near) | 0.0210 | 0.0210 | 0.0210 | changed below print precision |
| wobj_results | 307-311 |   CI (far) | 0.0291 | 0.0291 | 0.0291 | changed below print precision |
| wobj_results | 307-311 | p_func share of oracle | 24 | 24 | 24 | changed below print precision |
| wobj_results | 307-311 | free gene-backed reactions two-tier | 41.5 | 41.5 | 41.5 | changed below print precision |
| wobj_results | 307-311 | free gene-backed reactions delta | -8.9 | -8.9 | -8.9 | changed below print precision |
| wobj_results | 307-311 |   CI low | -10.3 | -10.3 | -10.3 | changed below print precision |
| wobj_results | 307-311 |   CI high | -7.5 | -7.5 | -7.5 | changed below print precision |
| wobj_results | 307-311 | free gene-backed reactions, p_func (real) | 32.7 | 32.7 | 32.7 | changed below print precision |
| wobj_results | 307-311 | free gene-backed reactions, parsimony (uniform) | 41.8 | 41.8 | 41.8 | changed below print precision |
| stageO_results | 319-324 | selected-shipped truth-absent flux | 0.0172 | 0.0172 | 0.0172 | changed below print precision |
| stageO_results | 319-324 |   CI (near) | 0.0118 | 0.0118 | 0.0118 | changed below print precision |
| stageO_results | 319-324 |   CI (far) | 0.0221 | 0.0221 | 0.0221 | changed below print precision |
| stageO_results | 319-324 | face minimum | 0.0076 | 0.0076 | 0.0076 | changed below print precision |
| stageO_results | 319-324 | headroom selected | 69 | 69 | 69 | changed below print precision |
| stageO_results | 319-324 | headroom shipped | 23 | 23 | 23 | changed below print precision |
| stageO_results | 319-324 | free reactions shuffled | 2.1 | 2.1 | 2.1 | changed below print precision |
| stageO_results | 319-324 | shuffled vs released | +0.0110 [0.0031, 0.0222] | +0.0110 [0.0031, 0.0222] | +0.0111 [0.0031, 0.0223] | manual value; source: OBJECTIVE_TUNING_RESULTS.md (re-numbered by the R4.3 TF agent from objective_tuning/summary.json); RUN4_TF_REPORT 'What Run 5 needs' 1 |
| stageO_results | 319-324 | face-min CI, estimate-alone 0.0128, exchange-weight 0.0043, log-cost 0.0039, G effect, n=1 0.0048 / 0.0124, 96 of 193, 0.0164 [0.0100, 0.0236], 69% / 23% | as printed | as printed | unchanged at print | manual value; source: OBJECTIVE_TUNING_RESULTS.md; RUN4_TF_REPORT 'What Run 5 needs' 1 (checked by hand against the new outputs) |
| pfunc_calibration | 225 (modal fractional weight 0.5), 353, 1172, 1225 (slope 0.78; informative only for >=3 refs) | calibration slope | 0.78 [0.76, 0.79] | 0.78 [0.76, 0.79] | 0.78 [0.76, 0.79] | manual value; source: data/calibration/README.md (regenerated by calibrate_pfunc.py, R4.3) |
| pfunc_calibration | 225 (modal fractional weight 0.5), 353, 1172, 1225 (slope 0.78; informative only for >=3 refs) | thresholding at 0.5 changes macro F1 by | +0.0007 | +0.0007 | +0.0007 | manual value; source: data/calibration/README.md key result 5 |
| pfunc_calibration | 225 (modal fractional weight 0.5), 353, 1172, 1225 (slope 0.78; informative only for >=3 refs) | modal fractional weight | 0.5 | 0.5 | see derived recon_modal_fractional_weight | manual value; source: W6 |
| stage_e_estimator | 354-360, 1172 (calibration), 1225 ('tuned on a benchmark reference of about 7,000 genomes') | isotonic recalibration lowers log-loss by | 0.0411 | 0.0411 | 0.0412 | manual value; source: PFUNC_ESTIMATORS_RESULTS.md section 1 (re-numbered by the R4.3 TF agent) |
| stage_e_estimator | 354-360, 1172 (calibration), 1225 ('tuned on a benchmark reference of about 7,000 genomes') | within-union macro dF1 | +0.0085 [+0.0072, +0.0096] | +0.0085 [+0.0072, +0.0096] | +0.0085 [+0.0072, +0.0097] | manual value; source: PFUNC_ESTIMATORS_RESULTS.md; make_deck.py PFE_WITHIN_UNION_F1 |
| stage_e_estimator | 354-360, 1172 (calibration), 1225 ('tuned on a benchmark reference of about 7,000 genomes') | 0.0889 -> 0.0578 (0.0311 [0.0297, 0.0322]); 0.209 -> 0.127; 0.543 -> 0.816; 0.805; 0.0408 [0.0384, 0.0424]; 0.739; 0.067 [0.062, 0.072]; +0.032; +0.0249 [+0.0228, +0.0265]; 13%; 57,021 / 49,048 cells | as printed | as printed | unchanged at print | manual value; source: PFUNC_ESTIMATORS_RESULTS.md; RUN4_TF_REPORT 'What Run 5 needs' 1 |
| disc.objective_test1 | 1172, 1225 | lowers truth-absent flux | 0.0061 | 0.0061 | 0.0061 | changed below print precision |
| disc.objective_test1 | 1172, 1225 |   CI (near) | 0.0025 | 0.0025 | 0.0025 | changed below print precision |
| disc.objective_test1 | 1172, 1225 |   CI (far) | 0.0097 | 0.0097 | 0.0097 | changed below print precision |
| disc.objective_test1 | 1172, 1225 | free gene-backed reactions cut | 21 | 21 | 21 | changed below print precision |
| disc.objective_test1 | 1172, 1225 | share of oracle | 24 | 24 | 24 | changed below print precision |
| disc.objective_test1 | 1172, 1225 | shuffled weights cut free reactions by (not printed: 'by more') | None | 33 | 33 | changed below print precision |
| disc.objective_test1 | 1172, 1225 | p_func weights cut free reactions by | 21 | 21 | 21 | changed below print precision |
| disc.objective_test2 | 1172, 1225 | selected-shipped | 0.0172 | 0.0172 | 0.0172 | changed below print precision |
| disc.objective_test2 | 1172, 1225 |   CI (near) | 0.0118 | 0.0118 | 0.0118 | changed below print precision |
| disc.objective_test2 | 1172, 1225 |   CI (far) | 0.0221 | 0.0221 | 0.0221 | changed below print precision |
| wobj_results | 307-311 | exchange-term n<=1e-6 spread | 357 | 357 | 357 |  |
| wobj_results | 307-311 | n=1 identical | 60 | 60 | 60 |  |
| wobj_results | 307-311 | exchange term share of the objective (median) | 98.9 | 98.9 | 98.9 |  |
| msprob_released_equivalence | 285, 313-314 | released_check pass | 360 | 360 | 360 |  |
| msprob_released_equivalence | 285, 313-314 | released_check n | 360 | 360 | 360 |  |
| stageO_design_and_split | 317-318 | g_meas | 0.084 | 0.084 | 0.084 |  |
| stageO_design_and_split | 317-318 | tune half amplicons | 167 | 167 | 167 |  |
| stageO_design_and_split | 317-318 | confirm half amplicons | 193 | 193 | 193 |  |
| stageO_results | 319-324 | free reactions selected | 13.6 | 13.6 | 13.6 |  |
| stageO_results | 319-324 | free reactions shipped | 35.4 | 35.4 | 35.4 |  |
| stageO_results | 319-324 | exchange flux | +2.5 | +2.5 | +2.5 |  |
| stageO_wrapper_equivalence | 325 | amplicons passing (bitwise) | 8 | 8 | 8 |  |
| stageO_wrapper_equivalence | 325 | amplicons compared | 8 | 8 | 8 |  |
| disc.objective_test1 | 1172, 1225 | exchange term share of the objective (median) | 98.9 | 98.9 | 98.9 |  |
| abs-objective-posthoc | 61 | (manual) | | | | **no number**: qualitative sentence; its numbers are under methods:stageO_results |
| wobj_results | 307-311 | (manual) | | | | **not recomputed** (weighted term median 0.11%): a per-row share computed in WEIGHTED_OBJECTIVE_RESULTS.md section 1 (internal-term decomposition); the wobj rows moved only in 2 of 360 jobs (R4.3), so it is expected unchanged at print |
| stageO_design_and_split | 317-318 | (manual) | | | | **unchanged** (121 / 80 groups): design: the wobj population has the same 54,243 candidates and 360 keys (R4.3 population check), so the hash split is unchanged |
| disc.objective_test2 | 1172, 1225 | (manual) | | | | **unchanged** (121 of 201 / 80 groups): design (same 360 keys) |
| lim.reconstruction_tests | 1225 | (manual) | | | | **unchanged** (300/60 amplicons, 121 of 201 groups, slope 0.78, ~7,000 genomes): design restatement; the calibration slope is unchanged at print (val16s:pfunc_calibration) |

## R2.3

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| tab:magvs16s | 842-865 | FullLength16S n | 1,916 | 1,916 | 1,972 | **CHANGED** |
| tab:magvs16s | 842-865 | FullLength16S MAG | 0.813 | 0.813 | 0.814 | **CHANGED** |
| tab:magvs16s | 842-865 | FullLength16S 16S | 0.710 | 0.710 | 0.712 | **CHANGED** |
| tab:magvs16s | 842-865 | FullLength16S delta | +0.103 | +0.103 | +0.102 | **CHANGED** |
| tab:magvs16s | 842-865 | FullLength16S pop MAG | 0.248 | 0.248 | 0.250 | **CHANGED** |
| tab:magvs16s | 842-865 | FullLength16S pop 16S | 0.605 | 0.605 | 0.612 | **CHANGED** |
| tab:magvs16s | 842-865 | FullLength16S coverage MAG | 0.304 | 0.304 | 0.307 | **CHANGED** |
| tab:magvs16s | 842-865 | FullLength16S coverage 16S | 0.921 | 0.921 | 0.931 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V2 n | 1,992 | 1,992 | 2,045 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V2 MAG | 0.814 | 0.814 | 0.815 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V2 median delta | +0.059 | +0.059 | +0.060 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V2 pop MAG | 0.242 | 0.242 | 0.245 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V2 pop 16S | 0.593 | 0.593 | 0.599 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V2 coverage MAG | 0.298 | 0.298 | 0.301 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V2 coverage 16S | 0.920 | 0.920 | 0.930 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V3 n | 2,020 | 2,020 | 2,072 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V3 MAG | 0.813 | 0.813 | 0.814 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V3 pop MAG | 0.243 | 0.243 | 0.245 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V3 pop 16S | 0.602 | 0.602 | 0.609 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V3 coverage MAG | 0.298 | 0.298 | 0.301 | **CHANGED** |
| tab:magvs16s | 842-865 | V1-V3 coverage 16S | 0.927 | 0.927 | 0.938 | **CHANGED** |
| tab:magvs16s | 842-865 | V3-V4 n | 1,885 | 1,885 | 1,951 | **CHANGED** |
| tab:magvs16s | 842-865 | V3-V4 MAG | 0.815 | 0.815 | 0.816 | **CHANGED** |
| tab:magvs16s | 842-865 | V3-V4 16S | 0.629 | 0.629 | 0.632 | **CHANGED** |
| tab:magvs16s | 842-865 | V3-V4 delta | +0.186 | +0.186 | +0.184 | **CHANGED** |
| tab:magvs16s | 842-865 | V3-V4 median delta | +0.148 | +0.148 | +0.145 | **CHANGED** |
| tab:magvs16s | 842-865 | V3-V4 pop MAG | 0.242 | 0.242 | 0.245 | **CHANGED** |
| tab:magvs16s | 842-865 | V3-V4 pop 16S | 0.516 | 0.516 | 0.522 | **CHANGED** |
| tab:magvs16s | 842-865 | V3-V4 coverage MAG | 0.298 | 0.298 | 0.300 | **CHANGED** |
| tab:magvs16s | 842-865 | V3-V4 coverage 16S | 0.875 | 0.875 | 0.885 | **CHANGED** |
| tab:magvs16s | 842-865 | V4 n | 2,299 | 2,299 | 2,441 | **CHANGED** |
| tab:magvs16s | 842-865 | V4 MAG | 0.810 | 0.810 | 0.809 | **CHANGED** |
| tab:magvs16s | 842-865 | V4 16S | 0.598 | 0.598 | 0.600 | **CHANGED** |
| tab:magvs16s | 842-865 | V4 delta | +0.212 | +0.212 | +0.209 | **CHANGED** |
| tab:magvs16s | 842-865 | V4 median delta | +0.171 | +0.171 | +0.170 | **CHANGED** |
| tab:magvs16s | 842-865 | V4 pop MAG | 0.253 | 0.253 | 0.258 | **CHANGED** |
| tab:magvs16s | 842-865 | V4 pop 16S | 0.483 | 0.483 | 0.493 | **CHANGED** |
| tab:magvs16s | 842-865 | V4 coverage MAG | 0.313 | 0.313 | 0.319 | **CHANGED** |
| tab:magvs16s | 842-865 | V4 coverage 16S | 0.863 | 0.863 | 0.881 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5 n | 2,334 | 2,334 | 2,464 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5 MAG | 0.809 | 0.809 | 0.808 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5 delta | +0.192 | +0.193 | +0.192 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5 median delta | +0.157 | +0.157 | +0.158 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5 pop MAG | 0.254 | 0.254 | 0.259 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5 pop 16S | 0.504 | 0.504 | 0.512 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5 coverage MAG | 0.314 | 0.314 | 0.320 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5 coverage 16S | 0.868 | 0.868 | 0.883 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5_944R n | 2,255 | 2,255 | 2,371 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5_944R 16S | 0.623 | 0.623 | 0.624 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5_944R delta | +0.188 | +0.188 | +0.187 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5_944R pop MAG | 0.257 | 0.257 | 0.262 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5_944R pop 16S | 0.510 | 0.510 | 0.517 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5_944R coverage MAG | 0.318 | 0.318 | 0.323 | **CHANGED** |
| tab:magvs16s | 842-865 | V4-V5_944R coverage 16S | 0.873 | 0.873 | 0.886 | **CHANGED** |
| tab:magvs16s | 842-865 | V7-V9 n | 2,301 | 2,301 | 2,430 | **CHANGED** |
| tab:magvs16s | 842-865 | V7-V9 MAG | 0.808 | 0.808 | 0.807 | **CHANGED** |
| tab:magvs16s | 842-865 | V7-V9 delta | +0.168 | +0.168 | +0.167 | **CHANGED** |
| tab:magvs16s | 842-865 | V7-V9 median delta | +0.126 | +0.126 | +0.125 | **CHANGED** |
| tab:magvs16s | 842-865 | V7-V9 pop MAG | 0.264 | 0.264 | 0.269 | **CHANGED** |
| tab:magvs16s | 842-865 | V7-V9 pop 16S | 0.523 | 0.523 | 0.533 | **CHANGED** |
| tab:magvs16s | 842-865 | V7-V9 coverage MAG | 0.326 | 0.326 | 0.332 | **CHANGED** |
| tab:magvs16s | 842-865 | V7-V9 coverage 16S | 0.880 | 0.880 | 0.896 | **CHANGED** |
| tab:magval | 896-919 | 1 resolve rate | 0.308 | 0.308 | 0.313 | **CHANGED** |
| tab:magval | 896-919 | 1 FL n | 1,916 | 1,916 | 1,972 | **CHANGED** |
| tab:magval | 896-919 | 1 FL MAG | 0.813 | 0.813 | 0.814 | **CHANGED** |
| tab:magval | 896-919 | 1 FL 16S | 0.710 | 0.710 | 0.712 | **CHANGED** |
| tab:magval | 896-919 | 1 FL delta | +0.103 | +0.103 | +0.102 | **CHANGED** |
| tab:magval | 896-919 | 1 V4 MAG | 0.810 | 0.810 | 0.809 | **CHANGED** |
| tab:magval | 896-919 | 1 V4 16S | 0.598 | 0.598 | 0.600 | **CHANGED** |
| tab:magval | 896-919 | 1 V4 delta | +0.212 | +0.212 | +0.209 | **CHANGED** |
| tab:magval | 896-919 | 0.5 resolve rate | 0.307 | 0.307 | 0.312 | **CHANGED** |
| tab:magval | 896-919 | 0.5 cond F1 | 0.807 | 0.807 | 0.806 | **CHANGED** |
| tab:magval | 896-919 | 0.5 FL n | 1,927 | 1,927 | 1,980 | **CHANGED** |
| tab:magval | 896-919 | 0.25 resolve rate | 0.309 | 0.309 | 0.315 | **CHANGED** |
| tab:magval | 896-919 | 0.25 FL n | 1,932 | 1,932 | 1,987 | **CHANGED** |
| tab:magval | 896-919 | 0.1 resolve rate | 0.295 | 0.295 | 0.301 | **CHANGED** |
| tab:magval | 896-919 | 0.1 cond F1 | 0.803 | 0.803 | 0.802 | **CHANGED** |
| tab:magval | 896-919 | 0.1 FL n | 1,845 | 1,845 | 1,903 | **CHANGED** |
| tab:magval | 896-919 | 0.05 resolve rate | 0.292 | 0.292 | 0.297 | **CHANGED** |
| tab:magval | 896-919 | 0.05 cond F1 | 0.797 | 0.797 | 0.796 | **CHANGED** |
| tab:magval | 896-919 | 0.05 FL n | 1,821 | 1,821 | 1,878 | **CHANGED** |
| tab:magval | 896-919 | 0.02 resolve rate | 0.283 | 0.283 | 0.289 | **CHANGED** |
| tab:magval | 896-919 | 0.02 FL n | 1,791 | 1,791 | 1,846 | **CHANGED** |
| tab:magval | 896-919 | 0.01 resolve rate | 0.253 | 0.253 | 0.258 | **CHANGED** |
| tab:magval | 896-919 | 0.01 FL n | 1,603 | 1,603 | 1,645 | **CHANGED** |
| tab:magval | 896-919 | 0.01 FL MAG | 0.783 | 0.783 | 0.785 | **CHANGED** |
| tab:magval | 896-919 | 0.01 FL 16S | 0.706 | 0.706 | 0.708 | **CHANGED** |
| tab:magval | 896-919 | 0.01 V4 16S | 0.593 | 0.593 | 0.596 | **CHANGED** |
| tab:magval | 896-919 | 0.01 V4 delta | +0.190 | +0.190 | +0.187 | **CHANGED** |
| tab:magval | 896-919 | V4 pairing n (min over tab levels) | 1,893 | 1,893 | 2,006 | **CHANGED** |
| tab:magval | 896-919 | V4 pairing n (max over tab levels) | 2,311 | 2,311 | 2,448 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | FullLength16S coverage original | 0.921 | 0.921 | 0.931 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | FullLength16S coverage kmeans | 0.750 | 0.750 | 0.754 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | FullLength16S coverage random | 0.650 | 0.650 | 0.655 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | FullLength16S cond F1 random | 0.641 | 0.641 | 0.642 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V1-V3 coverage original | 0.927 | 0.927 | 0.938 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V1-V3 coverage kmeans | 0.757 | 0.757 | 0.761 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V1-V3 coverage random | 0.664 | 0.664 | 0.667 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V4 coverage original | 0.863 | 0.863 | 0.881 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V4 coverage kmeans | 0.647 | 0.647 | 0.652 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V4 coverage random | 0.536 | 0.536 | 0.540 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V4 cond F1 kmeans | 0.541 | 0.541 | 0.540 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V4 cond F1 random | 0.536 | 0.536 | 0.535 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V7-V9 coverage original | 0.880 | 0.880 | 0.896 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V7-V9 coverage kmeans | 0.659 | 0.659 | 0.664 | **CHANGED** |
| tab:twopanel-16S | 1025-1050 | V7-V9 coverage random | 0.547 | 0.547 | 0.552 | **CHANGED** |
| tab:twopanel-MAG | 1025-1050 | 1 resolve original | 0.309 | 0.309 | 0.314 | **CHANGED** |
| tab:twopanel-MAG | 1025-1050 | 1 cond F1 original | 0.807 | 0.807 | 0.806 | **CHANGED** |
| tab:twopanel-MAG | 1025-1050 | 1 cond F1 random | 0.791 | 0.791 | 0.790 | **CHANGED** |
| tab:twopanel-MAG | 1025-1050 | 0.5 resolve original | 0.309 | 0.309 | 0.314 | **CHANGED** |
| tab:twopanel-MAG | 1025-1050 | 0.5 cond F1 random | 0.792 | 0.792 | 0.791 | **CHANGED** |
| tab:twopanel-MAG | 1025-1050 | 0.3 resolve original | 0.310 | 0.310 | 0.316 | **CHANGED** |
| tab:emptied | 1065-1080 | original N | 7,215 | 7,215 | 7,268 | **CHANGED** |
| tab:emptied | 1065-1080 | original 0-19 share | 84.4 | 84.4 | 84.3 | **CHANGED** |
| tab:emptied | 1065-1080 | original 0-19 uncovered | 0.037 | 0.037 | 0.026 | **CHANGED** |
| tab:emptied | 1065-1080 | original 20-49 share | 7.2 | 7.2 | 7.3 | **CHANGED** |
| tab:emptied | 1065-1080 | original 20-49 uncovered | 0.017 | 0.017 | 0.011 | **CHANGED** |
| tab:emptied | 1065-1080 | original 50-99 uncovered | 0.015 | 0.015 | 0.000 | **CHANGED** |
| tab:emptied | 1065-1080 | original 100 share | 4.7 | 4.7 | 4.6 | **CHANGED** |
| tab:emptied | 1065-1080 | kmeans N | 8,815 | 8,815 | 8,817 | **CHANGED** |
| tab:emptied | 1065-1080 | kmeans 0-19 uncovered | 0.028 | 0.028 | 0.022 | **CHANGED** |
| tab:emptied | 1065-1080 | kmeans 20-49 uncovered | 0.022 | 0.022 | 0.018 | **CHANGED** |
| tab:emptied | 1065-1080 | kmeans 50-99 uncovered | 0.032 | 0.032 | 0.030 | **CHANGED** |
| tab:emptied | 1065-1080 | random N | 8,762 | 8,762 | 8,780 | **CHANGED** |
| tab:emptied | 1065-1080 | random 0-19 share | 37.3 | 37.3 | 37.4 | **CHANGED** |
| tab:emptied | 1065-1080 | random 0-19 uncovered | 0.023 | 0.023 | 0.017 | **CHANGED** |
| tab:emptied | 1065-1080 | random 20-49 uncovered | 0.021 | 0.021 | 0.019 | **CHANGED** |
| tab:emptied | 1065-1080 | random 50-99 uncovered | 0.028 | 0.028 | 0.020 | **CHANGED** |
| tab:emptied | 1065-1080 | random 100 share | 33.4 | 33.4 | 33.3 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | FullLength16S phystd coverage original | 0.926 | 0.926 | 0.934 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | FullLength16S phystd coverage kmeans | 0.767 | 0.767 | 0.771 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | FullLength16S phystd coverage random | 0.682 | 0.682 | 0.687 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | FullLength16S delta coverage | -0.243 | -0.243 | -0.247 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | FullLength16S delta coverage CI lo | -0.255 | -0.255 | -0.257 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | FullLength16S delta coverage CI hi | -0.233 | -0.233 | -0.236 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | FullLength16S delta cond F1 | -0.019 | -0.019 | -0.018 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | FullLength16S delta cond F1 CI hi | -0.005 | -0.005 | -0.003 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V4 phystd coverage original | 0.855 | 0.855 | 0.868 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V4 phystd coverage kmeans | 0.665 | 0.665 | 0.671 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V4 phystd coverage random | 0.576 | 0.576 | 0.583 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V4 delta coverage | -0.279 | -0.279 | -0.285 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V4 delta coverage CI lo | -0.291 | -0.291 | -0.296 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V4 delta coverage CI hi | -0.268 | -0.268 | -0.274 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V4 delta cond F1 CI lo | -0.031 | -0.031 | -0.030 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V4 delta cond F1 CI hi | -0.017 | -0.017 | -0.016 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V7-V9 phystd coverage original | 0.870 | 0.870 | 0.884 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V7-V9 phystd coverage kmeans | 0.681 | 0.681 | 0.689 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V7-V9 phystd coverage random | 0.594 | 0.594 | 0.600 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V7-V9 delta coverage | -0.275 | -0.275 | -0.284 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V7-V9 delta coverage CI lo | -0.287 | -0.287 | -0.295 | **CHANGED** |
| tab:threepanel-genecap | 1101-1137 | V7-V9 delta coverage CI hi | -0.263 | -0.263 | -0.272 | **CHANGED** |
| tab:threepanel-panelref | 1101-1137 | 1 phystd cond F1 kmeans | 0.788 | 0.788 | 0.787 | **CHANGED** |
| tab:threepanel-panelref | 1101-1137 | 1 delta random-original | -0.026 | -0.026 | -0.027 | **CHANGED** |
| tab:threepanel-panelref | 1101-1137 | 1 delta CI lo | -0.032 | -0.032 | -0.033 | **CHANGED** |
| tab:threepanel-panelref | 1101-1137 | 0.3 phystd cond F1 kmeans | 0.786 | 0.786 | 0.785 | **CHANGED** |
| emptied_vs_uncovered_panels | 742 | FL emptied original | 4.7 | 4.7 | 4.6 | **CHANGED** |
| emptied_vs_uncovered_panels | 742 | FL survivors uncovered original | 3.4 | 3.4 | 2.4 | **CHANGED** |
| emptied_vs_uncovered_panels | 742 | FL survivors uncovered k-means | 2.8 | 2.8 | 2.3 | **CHANGED** |
| emptied_vs_uncovered_panels | 742 | pooled emptied original | 8.3 | 8.3 | 8.2 | **CHANGED** |
| emptied_vs_uncovered_panels | 742 | pooled uncovered survivors original | 3.1 | 3.1 | 1.8 | **CHANGED** |
| emptied_vs_uncovered_panels | 742 | pooled uncovered survivors k-means | 2.4 | 2.4 | 1.8 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | resolve 1 | 0.308 | 0.308 | 0.313 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | resolve 0.5 | 0.307 | 0.307 | 0.312 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | resolve 0.25 | 0.309 | 0.309 | 0.315 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | resolve 0.1 | 0.295 | 0.295 | 0.301 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | resolve 0.05 | 0.292 | 0.292 | 0.297 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | resolve 0.02 | 0.283 | 0.283 | 0.289 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | resolve 0.01 | 0.253 | 0.253 | 0.258 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | cond F1 0.5 | 0.807 | 0.807 | 0.806 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | cond F1 0.1 | 0.803 | 0.803 | 0.802 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | cond F1 0.05 | 0.797 | 0.797 | 0.796 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | resolve 100-30% (low) | 0.306 | 0.306 | 0.311 | **CHANGED** |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | resolve 100-30% (high) | 0.308 | 0.308 | 0.314 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | FL lead at 100% | 0.103 | 0.103 | 0.102 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | V4 lead at 100% | 0.212 | 0.212 | 0.209 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | FL n 1 | 1,916 | 1,916 | 1,972 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | FL n 0.5 | 1,927 | 1,927 | 1,980 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | FL n 0.25 | 1,932 | 1,932 | 1,987 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | FL n 0.1 | 1,845 | 1,845 | 1,903 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | FL n 0.05 | 1,821 | 1,821 | 1,878 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | FL n 0.02 | 1,791 | 1,791 | 1,846 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | FL n 0.01 | 1,603 | 1,603 | 1,645 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | FL MAG 1% | 0.783 | 0.783 | 0.785 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | FL 16S 1% | 0.706 | 0.706 | 0.708 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | V4 16S 1% | 0.593 | 0.593 | 0.596 | **CHANGED** |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | V4 delta 1% | +0.190 | +0.190 | +0.187 | **CHANGED** |
| magval_no_crossover | 826-830 | V4 lead at 1% | 0.190 | 0.190 | 0.187 | **CHANGED** |
| magval_no_crossover | 826-830 | MAG wins | 1,129 | 1,129 | 1,161 | **CHANGED** |
| magval_no_crossover | 826-830 | 16S wins | 473 | 473 | 483 | **CHANGED** |
| magval_no_crossover | 826-830 | p | 1e-79 | 1e-79 | 3e-83 | **CHANGED** |
| magval_no_crossover | 826-830 | resolve 100% | 0.308 | 0.308 | 0.313 | **CHANGED** |
| magval_no_crossover | 826-830 | resolve 1% | 0.253 | 0.253 | 0.258 | **CHANGED** |
| magval_vs_every_region | 832-835 | F1 delta FL/V1-V2/V1-V3 (low) | +0.103 | 0.103 | 0.102 | **CHANGED** |
| magval_vs_every_region | 832-835 | F1 delta V3-onward (low) | +0.168 | 0.168 | 0.167 | **CHANGED** |
| magval_vs_every_region | 832-835 | F1 delta V3-onward (high) | +0.212 | 0.212 | 0.209 | **CHANGED** |
| magval_vs_every_region | 832-835 | MAG recall (high) | 0.824 | 0.824 | 0.825 | **CHANGED** |
| magval_vs_every_region | 832-835 | 16S recall (low) | 0.703 | 0.703 | 0.704 | **CHANGED** |
| magval_vs_every_region | 832-835 | recall delta FL/V1 (low) | 0.055 | 0.055 | 0.056 | **CHANGED** |
| magval_vs_every_region | 832-835 | recall delta FL/V1 (high) | 0.083 | 0.083 | 0.082 | **CHANGED** |
| magval_vs_every_region | 832-835 | recall delta V3-onward (low) | 0.097 | 0.097 | 0.095 | **CHANGED** |
| magval_vs_every_region | 832-835 | recall delta V3-onward (high) | 0.114 | 0.114 | 0.113 | **CHANGED** |
| magval_vs_every_region | 832-835 | max F1 p | 1e-161 | 6.283494966921468e-162 | 9.359008488785446e-167 | **CHANGED** |
| magval_vs_every_region | 832-835 | max recall p | 1e-29 | 5.348998555164015e-30 | 1.1445113493896564e-30 | **CHANGED** |
| tab_magvs16s | 842-868 | FullLength16S n | 1,916 | 1,916 | 1,972 | **CHANGED** |
| tab_magvs16s | 842-868 | FullLength16S MAG | 0.813 | 0.813 | 0.814 | **CHANGED** |
| tab_magvs16s | 842-868 | FullLength16S 16S | 0.710 | 0.710 | 0.712 | **CHANGED** |
| tab_magvs16s | 842-868 | FullLength16S delta | +0.103 | +0.103 | +0.102 | **CHANGED** |
| tab_magvs16s | 842-868 | FullLength16S pop MAG | 0.248 | 0.248 | 0.250 | **CHANGED** |
| tab_magvs16s | 842-868 | FullLength16S pop 16S | 0.605 | 0.605 | 0.612 | **CHANGED** |
| tab_magvs16s | 842-868 | FullLength16S coverage MAG | 0.304 | 0.304 | 0.307 | **CHANGED** |
| tab_magvs16s | 842-868 | FullLength16S coverage 16S | 0.921 | 0.921 | 0.931 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V2 n | 1,992 | 1,992 | 2,045 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V2 MAG | 0.814 | 0.814 | 0.815 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V2 median delta | +0.059 | +0.059 | +0.060 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V2 pop MAG | 0.242 | 0.242 | 0.245 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V2 pop 16S | 0.593 | 0.593 | 0.599 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V2 coverage MAG | 0.298 | 0.298 | 0.301 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V2 coverage 16S | 0.920 | 0.920 | 0.930 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V3 n | 2,020 | 2,020 | 2,072 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V3 MAG | 0.813 | 0.813 | 0.814 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V3 pop MAG | 0.243 | 0.243 | 0.245 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V3 pop 16S | 0.602 | 0.602 | 0.609 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V3 coverage MAG | 0.298 | 0.298 | 0.301 | **CHANGED** |
| tab_magvs16s | 842-868 | V1-V3 coverage 16S | 0.927 | 0.927 | 0.938 | **CHANGED** |
| tab_magvs16s | 842-868 | V3-V4 n | 1,885 | 1,885 | 1,951 | **CHANGED** |
| tab_magvs16s | 842-868 | V3-V4 MAG | 0.815 | 0.815 | 0.816 | **CHANGED** |
| tab_magvs16s | 842-868 | V3-V4 16S | 0.629 | 0.629 | 0.632 | **CHANGED** |
| tab_magvs16s | 842-868 | V3-V4 delta | +0.186 | +0.186 | +0.184 | **CHANGED** |
| tab_magvs16s | 842-868 | V3-V4 median delta | +0.148 | +0.148 | +0.145 | **CHANGED** |
| tab_magvs16s | 842-868 | V3-V4 pop MAG | 0.242 | 0.242 | 0.245 | **CHANGED** |
| tab_magvs16s | 842-868 | V3-V4 pop 16S | 0.516 | 0.516 | 0.522 | **CHANGED** |
| tab_magvs16s | 842-868 | V3-V4 coverage MAG | 0.298 | 0.298 | 0.300 | **CHANGED** |
| tab_magvs16s | 842-868 | V3-V4 coverage 16S | 0.875 | 0.875 | 0.885 | **CHANGED** |
| tab_magvs16s | 842-868 | V4 n | 2,299 | 2,299 | 2,441 | **CHANGED** |
| tab_magvs16s | 842-868 | V4 MAG | 0.810 | 0.810 | 0.809 | **CHANGED** |
| tab_magvs16s | 842-868 | V4 16S | 0.598 | 0.598 | 0.600 | **CHANGED** |
| tab_magvs16s | 842-868 | V4 delta | +0.212 | +0.212 | +0.209 | **CHANGED** |
| tab_magvs16s | 842-868 | V4 median delta | +0.171 | +0.171 | +0.170 | **CHANGED** |
| tab_magvs16s | 842-868 | V4 pop MAG | 0.253 | 0.253 | 0.258 | **CHANGED** |
| tab_magvs16s | 842-868 | V4 pop 16S | 0.483 | 0.483 | 0.493 | **CHANGED** |
| tab_magvs16s | 842-868 | V4 coverage MAG | 0.313 | 0.313 | 0.319 | **CHANGED** |
| tab_magvs16s | 842-868 | V4 coverage 16S | 0.863 | 0.863 | 0.881 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5 n | 2,334 | 2,334 | 2,464 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5 MAG | 0.809 | 0.809 | 0.808 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5 delta | +0.192 | +0.193 | +0.192 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5 median delta | +0.157 | +0.157 | +0.158 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5 pop MAG | 0.254 | 0.254 | 0.259 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5 pop 16S | 0.504 | 0.504 | 0.512 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5 coverage MAG | 0.314 | 0.314 | 0.320 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5 coverage 16S | 0.868 | 0.868 | 0.883 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5_944R n | 2,255 | 2,255 | 2,371 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5_944R 16S | 0.623 | 0.623 | 0.624 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5_944R delta | +0.188 | +0.188 | +0.187 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5_944R pop MAG | 0.257 | 0.257 | 0.262 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5_944R pop 16S | 0.510 | 0.510 | 0.517 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5_944R coverage MAG | 0.318 | 0.318 | 0.323 | **CHANGED** |
| tab_magvs16s | 842-868 | V4-V5_944R coverage 16S | 0.873 | 0.873 | 0.886 | **CHANGED** |
| tab_magvs16s | 842-868 | V7-V9 n | 2,301 | 2,301 | 2,430 | **CHANGED** |
| tab_magvs16s | 842-868 | V7-V9 MAG | 0.808 | 0.808 | 0.807 | **CHANGED** |
| tab_magvs16s | 842-868 | V7-V9 delta | +0.168 | +0.168 | +0.167 | **CHANGED** |
| tab_magvs16s | 842-868 | V7-V9 median delta | +0.126 | +0.126 | +0.125 | **CHANGED** |
| tab_magvs16s | 842-868 | V7-V9 pop MAG | 0.264 | 0.264 | 0.269 | **CHANGED** |
| tab_magvs16s | 842-868 | V7-V9 pop 16S | 0.523 | 0.523 | 0.533 | **CHANGED** |
| tab_magvs16s | 842-868 | V7-V9 coverage MAG | 0.326 | 0.326 | 0.332 | **CHANGED** |
| tab_magvs16s | 842-868 | V7-V9 coverage 16S | 0.880 | 0.880 | 0.896 | **CHANGED** |
| magval_population_view | 837-840 | FL pop 16S | 0.605 | 0.605 | 0.612 | **CHANGED** |
| magval_population_view | 837-840 | FL pop MAG | 0.248 | 0.248 | 0.250 | **CHANGED** |
| magval_population_view | 837-840 | V4 pop 16S | 0.483 | 0.483 | 0.493 | **CHANGED** |
| magval_population_view | 837-840 | V4 pop MAG | 0.253 | 0.253 | 0.258 | **CHANGED** |
| magval_population_view | 837-840 | amplicon coverage (low) | 86 | 86 | 88 | **CHANGED** |
| magval_population_view | 837-840 | amplicon coverage (high) | 93 | 93 | 94 | **CHANGED** |
| magval_af_gate_counterfactual | 878-880 | qualifying share at gate 0.15, floor 0.83 | 0.311 | 0.311 | 0.314 | **CHANGED** |
| magval_af_gate_counterfactual | 878-880 | qualifying share with floor 0 | 0.316 | 0.316 | 0.319 | **CHANGED** |
| magval_af_gate_counterfactual | 878-880 | qualifying share at gate 0.02 | 0.452 | 0.452 | 0.454 | **CHANGED** |
| magval_af_gate_counterfactual | 878-880 | newly admitted with prediction | 1,367 | 1,367 | 1,391 | **CHANGED** |
| magval_af_gate_counterfactual | 878-880 | mean F1 newly admitted | 0.599 | 0.599 | 0.600 | **CHANGED** |
| magval_af_gate_counterfactual | 878-880 | resolve rate | 0.308 | 0.308 | 0.313 | **CHANGED** |
| magval_af_gate_counterfactual | 878-880 | newly admitted above the FL 16S level (share; min over the three levels) | None | 26.7 | 26.5 | **CHANGED** |
| magval_af_gate_counterfactual | 878-880 | pooled mean F1 | 0.742 | 0.742 | 0.743 | **CHANGED** |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 resolve original | 0.309 | 0.309 | 0.314 | **CHANGED** |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 cond F1 original | 0.807 | 0.807 | 0.806 | **CHANGED** |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 cond F1 random | 0.791 | 0.791 | 0.790 | **CHANGED** |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 phystd cond F1 kmeans | 0.788 | 0.788 | 0.787 | **CHANGED** |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 delta random-original | -0.026 | -0.026 | -0.027 | **CHANGED** |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 delta CI lo | -0.032 | -0.032 | -0.033 | **CHANGED** |
| twopanel_mag_abstention_strata | 1085-1087 | singleton abstention k-means | 0.6895 | 0.6895 | 0.6889 | **CHANGED** |
| twopanel_mag_abstention_strata | 1085-1087 | singleton abstention original | 0.6883 | 0.6883 | 0.6823 | **CHANGED** |
| twopanel_mag_abstention_strata | 1085-1087 | overall k-means | 0.658 | 0.658 | 0.657 | **CHANGED** |
| twopanel_mag_abstention_strata | 1085-1087 | overall original | 0.688 | 0.688 | 0.682 | **CHANGED** |
| twopanel_mag_abstention_strata | 1085-1087 | k-means without >=20 stratum | 0.687 | 0.687 | 0.686 | **CHANGED** |
| panels.condF1_crude_bounds | 1018-1019 | random-original crude cond F1 V7-V9 | 0.029 | 0.029 | 0.030 | **CHANGED** |
| panels.condF1_crude_bounds | 1018-1019 | random vs original MAG crude cond F1, max |gap| over completeness (printed 'within 0.016') | 0.016 | 0.016 | 0.017 | **CHANGED** |
| panels.condF1_crude_bounds | 1018-1019 | k-means vs original MAG crude cond F1, max |gap| over completeness (printed 'within 0.010') | 0.010 | 0.00947194126350226 | 0.009427861404042726 | **CHANGED** |
| panels.condF1_crude_bounds | 1018-1019 | original ahead at every MAG completeness (min margin over both panels) | 0 | 0.007980254886422511 | 0.007897000489042849 | **CHANGED** |
| panels.coverage_disagrees | 1021-1023 | MAG resolve original 100% | 0.309 | 0.309 | 0.314 | **CHANGED** |
| panels.coverage_disagrees | 1021-1023 | FL coverage k-means | 0.750 | 0.750 | 0.754 | **CHANGED** |
| panels.coverage_disagrees | 1021-1023 | FL coverage random | 0.650 | 0.650 | 0.655 | **CHANGED** |
| panels.coverage_disagrees | 1021-1023 | FL coverage original | 0.921 | 0.921 | 0.931 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | FullLength16S coverage original | 0.921 | 0.921 | 0.931 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | FullLength16S coverage kmeans | 0.750 | 0.750 | 0.754 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | FullLength16S coverage random | 0.650 | 0.650 | 0.655 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | FullLength16S cond F1 random | 0.641 | 0.641 | 0.642 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V1-V3 coverage original | 0.927 | 0.927 | 0.938 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V1-V3 coverage kmeans | 0.757 | 0.757 | 0.761 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V1-V3 coverage random | 0.664 | 0.664 | 0.667 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V4 coverage original | 0.863 | 0.863 | 0.881 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V4 coverage kmeans | 0.647 | 0.647 | 0.652 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V4 coverage random | 0.536 | 0.536 | 0.540 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V4 cond F1 kmeans | 0.541 | 0.541 | 0.540 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V4 cond F1 random | 0.536 | 0.536 | 0.535 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V7-V9 coverage original | 0.880 | 0.880 | 0.896 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V7-V9 coverage kmeans | 0.659 | 0.659 | 0.664 | **CHANGED** |
| panels.tab_twopanel_16S | 1025-1043 | V7-V9 coverage random | 0.547 | 0.547 | 0.552 | **CHANGED** |
| panels.tab_twopanel_MAG | 1044-1048 | 1 resolve original | 0.309 | 0.309 | 0.314 | **CHANGED** |
| panels.tab_twopanel_MAG | 1044-1048 | 1 cond F1 original | 0.807 | 0.807 | 0.806 | **CHANGED** |
| panels.tab_twopanel_MAG | 1044-1048 | 1 cond F1 random | 0.791 | 0.791 | 0.790 | **CHANGED** |
| panels.tab_twopanel_MAG | 1044-1048 | 0.5 resolve original | 0.309 | 0.309 | 0.314 | **CHANGED** |
| panels.tab_twopanel_MAG | 1044-1048 | 0.5 cond F1 random | 0.792 | 0.792 | 0.791 | **CHANGED** |
| panels.tab_twopanel_MAG | 1044-1048 | 0.3 resolve original | 0.310 | 0.310 | 0.316 | **CHANGED** |
| panels.fl_gap_standardized | 1055-1056 | phystd k-means-original FL coverage | 0.159 | 0.159 | 0.163 | **CHANGED** |
| panels.fl_gap_standardized | 1055-1056 | crude k-means-original | 0.171 | 0.171 | 0.177 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | FL survivors uncovered original | 3.4 | 3.4 | 2.4 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | FL survivors uncovered k-means | 2.8 | 2.8 | 2.3 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | FL survivors uncovered random | 2.3 | 2.4 | 1.8 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | FL strata uncovered (low) | 1.5 | 1.5 | 0.0 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | FL strata uncovered (high) | 3.7 | 3.7 | 3.0 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | pooled original | 3.1 | 3.1 | 1.8 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | pooled k-means | 2.4 | 2.4 | 1.8 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | pooled random | 1.8 | 1.8 | 1.3 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | emptied / surviving uncovered rate (low, fold) | 29 | 29 | 42 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | emptied / surviving uncovered rate (high, fold) | 43 | 43 | 55 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | pooled spread (points) | 1.3 | 1.3 | 0.5 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | original-random FL difference (points) | 1.1 | 1.1 | 0.6 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | CMH original vs random p | 0.004 | 0.004 | 0.117 | **CHANGED**; MEANING: significance crosses 0.05 (0.00407 -> 0.117) |
| panels.surviving_uncovered | 1057-1060 | CMH original vs k-means p | 0.15 | 0.15 | 0.99 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | stratum-standardized uncovered rate original | None | 3.0 | 2.0 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | stratum-standardized uncovered rate kmeans | None | 2.8 | 2.3 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | stratum-standardized uncovered rate random | None | 2.3 | 1.8 | **CHANGED** |
| panels.surviving_uncovered | 1057-1060 | unstratified Fisher p original vs random | None | 0.000 | 0.023 | **CHANGED**; MEANING: significance crosses 0.01 (0.0003 -> 0.0229) |
| panels.emptied_arithmetic | 1061-1063 | emptied random | 33.4 | 33.4 | 33.3 | **CHANGED** |
| panels.emptied_arithmetic | 1061-1063 | emptied original | 4.7 | 4.7 | 4.6 | **CHANGED** |
| panels.emptied_arithmetic | 1061-1063 | 1-emptied random | 0.666 | 0.666 | 0.667 | **CHANGED** |
| panels.emptied_arithmetic | 1061-1063 | 1-emptied original | 0.953 | 0.953 | 0.954 | **CHANGED** |
| panels.emptied_arithmetic | 1061-1063 | observed random | 0.650 | 0.650 | 0.655 | **CHANGED** |
| panels.emptied_arithmetic | 1061-1063 | observed k-means | 0.750 | 0.750 | 0.754 | **CHANGED** |
| panels.emptied_arithmetic | 1061-1063 | observed original | 0.921 | 0.921 | 0.931 | **CHANGED** |
| panels.surviving_pred_empty | 1064 | original non-strict uncovered survivors | 107 | 107 | 36 | **CHANGED** |
| panels.surviving_pred_empty | 1064 | kmeans non-strict uncovered survivors | 76 | 76 | 42 | **CHANGED** |
| panels.surviving_pred_empty | 1064 | random non-strict uncovered survivors | 62 | 62 | 31 | **CHANGED** |
| panels.tab_emptied_shares | 1066-1081 | original N | 7,215 | 7,215 | 7,268 | **CHANGED** |
| panels.tab_emptied_shares | 1066-1081 | original 0-19 share | 84.4 | 84.4 | 84.3 | **CHANGED** |
| panels.tab_emptied_shares | 1066-1081 | original 20-49 share | 7.2 | 7.2 | 7.3 | **CHANGED** |
| panels.tab_emptied_shares | 1066-1081 | original 100 share | 4.7 | 4.7 | 4.6 | **CHANGED** |
| panels.tab_emptied_shares | 1066-1081 | kmeans N | 8,815 | 8,815 | 8,817 | **CHANGED** |
| panels.tab_emptied_shares | 1066-1081 | random N | 8,762 | 8,762 | 8,780 | **CHANGED** |
| panels.tab_emptied_shares | 1066-1081 | random 0-19 share | 37.3 | 37.3 | 37.4 | **CHANGED** |
| panels.tab_emptied_shares | 1066-1081 | random 100 share | 33.4 | 33.4 | 33.3 | **CHANGED** |
| panels.tab_emptied_uncovered | 1066-1081 | original 0-19 uncovered | 0.037 | 0.037 | 0.026 | **CHANGED** |
| panels.tab_emptied_uncovered | 1066-1081 | original 20-49 uncovered | 0.017 | 0.017 | 0.011 | **CHANGED** |
| panels.tab_emptied_uncovered | 1066-1081 | original 50-99 uncovered | 0.015 | 0.015 | 0.000 | **CHANGED** |
| panels.tab_emptied_uncovered | 1066-1081 | kmeans 0-19 uncovered | 0.028 | 0.028 | 0.022 | **CHANGED** |
| panels.tab_emptied_uncovered | 1066-1081 | kmeans 20-49 uncovered | 0.022 | 0.022 | 0.018 | **CHANGED** |
| panels.tab_emptied_uncovered | 1066-1081 | kmeans 50-99 uncovered | 0.032 | 0.032 | 0.030 | **CHANGED** |
| panels.tab_emptied_uncovered | 1066-1081 | random 0-19 uncovered | 0.023 | 0.023 | 0.017 | **CHANGED** |
| panels.tab_emptied_uncovered | 1066-1081 | random 20-49 uncovered | 0.021 | 0.021 | 0.019 | **CHANGED** |
| panels.tab_emptied_uncovered | 1066-1081 | random 50-99 uncovered | 0.028 | 0.028 | 0.020 | **CHANGED** |
| panels.std_gap_share | 1095-1096 | standardized / crude gap (high) | 92 | 92 | 91 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | FullLength16S phystd coverage original | 0.926 | 0.926 | 0.934 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | FullLength16S phystd coverage kmeans | 0.767 | 0.767 | 0.771 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | FullLength16S phystd coverage random | 0.682 | 0.682 | 0.687 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | FullLength16S delta coverage | -0.243 | -0.243 | -0.247 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | FullLength16S delta coverage CI lo | -0.255 | -0.255 | -0.257 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | FullLength16S delta coverage CI hi | -0.233 | -0.233 | -0.236 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | FullLength16S delta cond F1 | -0.019 | -0.019 | -0.018 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | FullLength16S delta cond F1 CI hi | -0.005 | -0.005 | -0.003 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V4 phystd coverage original | 0.855 | 0.855 | 0.868 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V4 phystd coverage kmeans | 0.665 | 0.665 | 0.671 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V4 phystd coverage random | 0.576 | 0.576 | 0.583 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V4 delta coverage | -0.279 | -0.279 | -0.285 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V4 delta coverage CI lo | -0.291 | -0.291 | -0.296 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V4 delta coverage CI hi | -0.268 | -0.268 | -0.274 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V4 delta cond F1 CI lo | -0.031 | -0.031 | -0.030 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V4 delta cond F1 CI hi | -0.017 | -0.017 | -0.016 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V7-V9 phystd coverage original | 0.870 | 0.870 | 0.884 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V7-V9 phystd coverage kmeans | 0.681 | 0.681 | 0.689 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V7-V9 phystd coverage random | 0.594 | 0.594 | 0.600 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V7-V9 delta coverage | -0.275 | -0.275 | -0.284 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V7-V9 delta coverage CI lo | -0.287 | -0.287 | -0.295 | **CHANGED** |
| panels.tab_threepanel_genecap_block | 1120-1124 | V7-V9 delta coverage CI hi | -0.263 | -0.263 | -0.272 | **CHANGED** |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 1 phystd cond F1 kmeans | 0.788 | 0.788 | 0.787 | **CHANGED** |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 1 delta random-original | -0.026 | -0.026 | -0.027 | **CHANGED** |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 1 delta CI lo | -0.032 | -0.032 | -0.033 | **CHANGED** |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 0.3 phystd cond F1 kmeans | 0.786 | 0.786 | 0.785 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus coverage FullLength16S CI lo | -0.251 | -0.251 | -0.249 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus coverage FullLength16S CI hi | -0.189 | -0.189 | -0.192 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus coverage V1-V3 CI lo | -0.238 | -0.238 | -0.236 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus coverage V1-V3 CI hi | -0.191 | -0.191 | -0.194 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus coverage V3-V4 CI lo | -0.274 | -0.274 | -0.275 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus coverage V4 CI hi | -0.238 | -0.238 | -0.239 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus cond acc FullLength16S CI lo | +0.013 | +0.013 | +0.014 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus cond acc FullLength16S CI hi | +0.061 | +0.061 | +0.062 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus cond acc V1-V3 CI lo | +0.039 | +0.039 | +0.038 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus cond acc V1-V3 CI hi | +0.087 | +0.087 | +0.086 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus cond acc V3-V4 CI lo | -0.040 | -0.040 | -0.039 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus cond acc V3-V4 CI hi | +0.019 | +0.019 | +0.016 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus cond acc V4 CI lo | -0.015 | -0.015 | -0.012 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus cond acc V4 CI hi | +0.046 | +0.046 | +0.045 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus cond acc V7-V9 CI lo | -0.044 | -0.044 | -0.043 | **CHANGED** |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus cond acc V7-V9 CI hi | +0.011 | +0.011 | +0.010 | **CHANGED** |
| panels.selfincl_genecap_penalty | 1149-1150 | V4 original phystd | 0.837 | 0.837 | 0.836 | **CHANGED** |
| panels.selfincl_genecap_penalty | 1149-1150 | V4 delta | -0.107 | -0.107 | -0.106 | **CHANGED** |
| panels.selfincl_genecap_penalty | 1149-1150 | V4 delta CI lo | -0.114 | -0.114 | -0.113 | **CHANGED** |
| panels.selfincl_genecap_penalty | 1149-1150 | V4 delta CI hi | -0.100 | -0.100 | -0.099 | **CHANGED** |
| abs-mag-vs-16s | 63 | MAG-16S conditional delta at 100% (low, FL) | +0.10 | +0.10 | +0.10 | changed below print precision |
| abs-mag-vs-16s | 63 | MAG-16S conditional delta at 100% (high, V4) | +0.21 | +0.21 | +0.21 | changed below print precision |
| abs-mag-vs-16s | 63 | FL delta at 1% | +0.08 | +0.08 | +0.08 | changed below print precision |
| abs-mag-vs-16s | 63 | V4 delta at 1% | +0.19 | +0.19 | +0.19 | changed below print precision |
| abs-mag-vs-16s | 63 | resolve rate at 100% ('~30%') | 0.3 | 0.3 | 0.3 | changed below print precision |
| abs-threepanel-genecap-quality | 69 | max |phystd delta| cond F1 (random V7-V9) | 0.028 | 0.028 | 0.028 | changed below print precision |
| abs-threepanel-genecap-coverage | 70 | random-original genecap coverage (max |.|, V4-V5_944R) | 0.29 | 0.29 | 0.29 | changed below print precision |
| abs-selfincl-panel-dependence | 71 | include genecap cond F1 random-original V4 | 0.11 | 0.11 | 0.11 | changed below print precision |
| abs-selfincl-panel-dependence | 71 | include genecap cond F1 random-original FL | 0.02 | 0.02 | 0.02 | changed below print precision |
| tab:magvs16s | 842-865 | FullLength16S median delta | +0.056 | +0.056 | +0.056 | changed below print precision |
| tab:magvs16s | 842-865 | V1-V2 16S | 0.698 | 0.698 | 0.698 | changed below print precision |
| tab:magvs16s | 842-865 | V1-V2 delta | +0.116 | +0.116 | +0.116 | changed below print precision |
| tab:magvs16s | 842-865 | V1-V3 16S | 0.704 | 0.704 | 0.704 | changed below print precision |
| tab:magvs16s | 842-865 | V1-V3 delta | +0.110 | +0.110 | +0.110 | changed below print precision |
| tab:magvs16s | 842-865 | V1-V3 median delta | +0.052 | +0.052 | +0.052 | changed below print precision |
| tab:magvs16s | 842-865 | V4-V5 16S | 0.616 | 0.616 | 0.616 | changed below print precision |
| tab:magvs16s | 842-865 | V4-V5_944R MAG | 0.811 | 0.811 | 0.811 | changed below print precision |
| tab:magvs16s | 842-865 | V4-V5_944R median delta | +0.157 | +0.157 | +0.157 | changed below print precision |
| tab:magvs16s | 842-865 | V7-V9 16S | 0.641 | 0.641 | 0.641 | changed below print precision |
| tab:magval | 896-919 | 1 cond F1 | 0.807 | 0.807 | 0.807 | changed below print precision |
| tab:magval | 896-919 | 0.25 cond F1 | 0.806 | 0.806 | 0.806 | changed below print precision |
| tab:magval | 896-919 | 0.02 cond F1 | 0.791 | 0.791 | 0.791 | changed below print precision |
| tab:magval | 896-919 | 0.01 cond F1 | 0.783 | 0.783 | 0.783 | changed below print precision |
| tab:magval | 896-919 | 0.01 FL delta | +0.077 | +0.077 | +0.077 | changed below print precision |
| tab:magval | 896-919 | 0.01 V4 MAG | 0.783 | 0.783 | 0.783 | changed below print precision |
| tab:twopanel-16S | 1025-1050 | FullLength16S cond F1 original | 0.657 | 0.657 | 0.657 | changed below print precision |
| tab:twopanel-16S | 1025-1050 | FullLength16S cond F1 kmeans | 0.645 | 0.645 | 0.645 | changed below print precision |
| tab:twopanel-16S | 1025-1050 | V1-V3 cond F1 original | 0.649 | 0.649 | 0.649 | changed below print precision |
| tab:twopanel-16S | 1025-1050 | V1-V3 cond F1 kmeans | 0.641 | 0.641 | 0.641 | changed below print precision |
| tab:twopanel-16S | 1025-1050 | V1-V3 cond F1 random | 0.636 | 0.636 | 0.636 | changed below print precision |
| tab:twopanel-16S | 1025-1050 | V4 cond F1 original | 0.559 | 0.559 | 0.559 | changed below print precision |
| tab:twopanel-16S | 1025-1050 | V7-V9 cond F1 original | 0.594 | 0.594 | 0.594 | changed below print precision |
| tab:twopanel-16S | 1025-1050 | V7-V9 cond F1 kmeans | 0.571 | 0.571 | 0.571 | changed below print precision |
| tab:twopanel-16S | 1025-1050 | V7-V9 cond F1 random | 0.565 | 0.565 | 0.565 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 1 resolve kmeans | 0.342 | 0.342 | 0.342 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 1 resolve random | 0.428 | 0.428 | 0.428 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 1 cond F1 kmeans | 0.797 | 0.797 | 0.797 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 0.5 resolve kmeans | 0.343 | 0.343 | 0.343 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 0.5 resolve random | 0.428 | 0.428 | 0.428 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 0.5 cond F1 original | 0.806 | 0.806 | 0.806 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 0.5 cond F1 kmeans | 0.798 | 0.798 | 0.798 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 0.3 resolve kmeans | 0.345 | 0.345 | 0.345 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 0.3 resolve random | 0.430 | 0.430 | 0.430 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 0.3 cond F1 original | 0.804 | 0.804 | 0.804 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 0.3 cond F1 kmeans | 0.796 | 0.796 | 0.796 | changed below print precision |
| tab:twopanel-MAG | 1025-1050 | 0.3 cond F1 random | 0.795 | 0.795 | 0.795 | changed below print precision |
| tab:emptied | 1065-1080 | original 50-99 share | 3.8 | 3.8 | 3.8 | changed below print precision |
| tab:emptied | 1065-1080 | kmeans 0-19 share | 44.0 | 44.0 | 44.0 | changed below print precision |
| tab:emptied | 1065-1080 | kmeans 20-49 share | 17.9 | 17.9 | 17.9 | changed below print precision |
| tab:emptied | 1065-1080 | kmeans 50-99 share | 15.2 | 15.2 | 15.2 | changed below print precision |
| tab:emptied | 1065-1080 | kmeans 100 share | 22.9 | 22.9 | 22.9 | changed below print precision |
| tab:emptied | 1065-1080 | random 20-49 share | 14.8 | 14.8 | 14.8 | changed below print precision |
| tab:emptied | 1065-1080 | random 50-99 share | 14.5 | 14.5 | 14.5 | changed below print precision |
| tab:threepanel-genecap | 1101-1137 | FullLength16S delta cond F1 CI lo | -0.031 | -0.031 | -0.031 | changed below print precision |
| tab:threepanel-genecap | 1101-1137 | V4 delta cond F1 | -0.024 | -0.024 | -0.024 | changed below print precision |
| tab:threepanel-genecap | 1101-1137 | V7-V9 delta cond F1 | -0.028 | -0.028 | -0.028 | changed below print precision |
| tab:threepanel-genecap | 1101-1137 | V7-V9 delta cond F1 CI lo | -0.035 | -0.035 | -0.035 | changed below print precision |
| tab:threepanel-genecap | 1101-1137 | V7-V9 delta cond F1 CI hi | -0.020 | -0.020 | -0.020 | changed below print precision |
| tab:threepanel-panelref | 1101-1137 | 1 phystd cond F1 original | 0.809 | 0.809 | 0.809 | changed below print precision |
| tab:threepanel-panelref | 1101-1137 | 1 phystd cond F1 random | 0.783 | 0.783 | 0.783 | changed below print precision |
| tab:threepanel-panelref | 1101-1137 | 1 delta CI hi | -0.020 | -0.020 | -0.020 | changed below print precision |
| tab:threepanel-panelref | 1101-1137 | 0.3 phystd cond F1 original | 0.806 | 0.806 | 0.806 | changed below print precision |
| tab:threepanel-panelref | 1101-1137 | 0.3 phystd cond F1 random | 0.789 | 0.789 | 0.789 | changed below print precision |
| tab:threepanel-panelref | 1101-1137 | 0.3 delta random-original | -0.017 | -0.017 | -0.017 | changed below print precision |
| tab:threepanel-panelref | 1101-1137 | 0.3 delta CI lo | -0.023 | -0.023 | -0.023 | changed below print precision |
| tab:threepanel-panelref | 1101-1137 | 0.3 delta CI hi | -0.011 | -0.011 | -0.011 | changed below print precision |
| emptied_vs_uncovered_panels | 742 | pooled emptied k-means | 29.4 | 29.4 | 29.4 | changed below print precision |
| abstract_mag_margin | 63 | FL delta at 100% | +0.10 | +0.10 | +0.10 | changed below print precision |
| abstract_mag_margin | 63 | V4 delta at 100% | +0.21 | +0.21 | +0.21 | changed below print precision |
| abstract_mag_margin | 63 | FL delta at 1% | +0.08 | +0.08 | +0.08 | changed below print precision |
| abstract_mag_margin | 63 | V4 delta at 1% | +0.19 | +0.19 | +0.19 | changed below print precision |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | cond F1 1 | 0.807 | 0.807 | 0.807 | changed below print precision |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | cond F1 0.25 | 0.806 | 0.806 | 0.806 | changed below print precision |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | cond F1 0.02 | 0.791 | 0.791 | 0.791 | changed below print precision |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | cond F1 0.01 | 0.783 | 0.783 | 0.783 | changed below print precision |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | cond F1 100-30% (low) | 0.805 | 0.805 | 0.805 | changed below print precision |
| magval_completeness_null | 823; tab:magval 896-919 resolve-rate and Cond. F1 columns | cond F1 100-30% (high) | 0.807 | 0.807 | 0.807 | changed below print precision |
| magval_paired_main | 824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311) | V4 MAG 1% | 0.783 | 0.783 | 0.783 | changed below print precision |
| magval_no_crossover | 826-830 | FL lead at 1% | 0.077 | 0.077 | 0.077 | changed below print precision |
| magval_no_crossover | 826-830 | cond F1 cost 100% -> 1% | 0.024 | 0.024 | 0.024 | changed below print precision |
| magval_vs_every_region | 832-835 | F1 delta FL/V1-V2/V1-V3 (high) | +0.116 | 0.116 | 0.116 | changed below print precision |
| magval_vs_every_region | 832-835 | MAG recall (low) | 0.816 | 0.816 | 0.816 | changed below print precision |
| magval_vs_every_region | 832-835 | 16S recall (high) | 0.768 | 0.768 | 0.768 | changed below print precision |
| tab_magvs16s | 842-868 | FullLength16S median delta | +0.056 | +0.056 | +0.056 | changed below print precision |
| tab_magvs16s | 842-868 | V1-V2 16S | 0.698 | 0.698 | 0.698 | changed below print precision |
| tab_magvs16s | 842-868 | V1-V2 delta | +0.116 | +0.116 | +0.116 | changed below print precision |
| tab_magvs16s | 842-868 | V1-V3 16S | 0.704 | 0.704 | 0.704 | changed below print precision |
| tab_magvs16s | 842-868 | V1-V3 delta | +0.110 | +0.110 | +0.110 | changed below print precision |
| tab_magvs16s | 842-868 | V1-V3 median delta | +0.052 | +0.052 | +0.052 | changed below print precision |
| tab_magvs16s | 842-868 | V4-V5 16S | 0.616 | 0.616 | 0.616 | changed below print precision |
| tab_magvs16s | 842-868 | V4-V5_944R MAG | 0.811 | 0.811 | 0.811 | changed below print precision |
| tab_magvs16s | 842-868 | V4-V5_944R median delta | +0.157 | +0.157 | +0.157 | changed below print precision |
| tab_magvs16s | 842-868 | V7-V9 16S | 0.641 | 0.641 | 0.641 | changed below print precision |
| magval_af_gate_counterfactual | 878-880 | mean F1 already answered | 0.807 | 0.807 | 0.807 | changed below print precision |
| magval_af_gate_counterfactual | 878-880 | coverage gain | +45 | +45 | +45 | changed below print precision |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 resolve kmeans | 0.342 | 0.342 | 0.342 | changed below print precision |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 resolve random | 0.428 | 0.428 | 0.428 | changed below print precision |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 cond F1 kmeans | 0.797 | 0.797 | 0.797 | changed below print precision |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 phystd cond F1 original | 0.809 | 0.809 | 0.809 | changed below print precision |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 phystd cond F1 random | 0.783 | 0.783 | 0.783 | changed below print precision |
| twopanel_mag_panelarm | 1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134 | 1 delta CI hi | -0.020 | -0.020 | -0.020 | changed below print precision |
| panels.condF1_crude_bounds | 1018-1019 | k-means-original crude cond F1 V7-V9 | 0.023 | 0.023 | 0.023 | changed below print precision |
| panels.coverage_disagrees | 1021-1023 | MAG resolve k-means 100% | 0.342 | 0.342 | 0.342 | changed below print precision |
| panels.coverage_disagrees | 1021-1023 | MAG resolve random 100% | 0.428 | 0.428 | 0.428 | changed below print precision |
| panels.tab_twopanel_16S | 1025-1043 | FullLength16S cond F1 original | 0.657 | 0.657 | 0.657 | changed below print precision |
| panels.tab_twopanel_16S | 1025-1043 | FullLength16S cond F1 kmeans | 0.645 | 0.645 | 0.645 | changed below print precision |
| panels.tab_twopanel_16S | 1025-1043 | V1-V3 cond F1 original | 0.649 | 0.649 | 0.649 | changed below print precision |
| panels.tab_twopanel_16S | 1025-1043 | V1-V3 cond F1 kmeans | 0.641 | 0.641 | 0.641 | changed below print precision |
| panels.tab_twopanel_16S | 1025-1043 | V1-V3 cond F1 random | 0.636 | 0.636 | 0.636 | changed below print precision |
| panels.tab_twopanel_16S | 1025-1043 | V4 cond F1 original | 0.559 | 0.559 | 0.559 | changed below print precision |
| panels.tab_twopanel_16S | 1025-1043 | V7-V9 cond F1 original | 0.594 | 0.594 | 0.594 | changed below print precision |
| panels.tab_twopanel_16S | 1025-1043 | V7-V9 cond F1 kmeans | 0.571 | 0.571 | 0.571 | changed below print precision |
| panels.tab_twopanel_16S | 1025-1043 | V7-V9 cond F1 random | 0.565 | 0.565 | 0.565 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 1 resolve kmeans | 0.342 | 0.342 | 0.342 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 1 resolve random | 0.428 | 0.428 | 0.428 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 1 cond F1 kmeans | 0.797 | 0.797 | 0.797 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 0.5 resolve kmeans | 0.343 | 0.343 | 0.343 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 0.5 resolve random | 0.428 | 0.428 | 0.428 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 0.5 cond F1 original | 0.806 | 0.806 | 0.806 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 0.5 cond F1 kmeans | 0.798 | 0.798 | 0.798 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 0.3 resolve kmeans | 0.345 | 0.345 | 0.345 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 0.3 resolve random | 0.430 | 0.430 | 0.430 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 0.3 cond F1 original | 0.804 | 0.804 | 0.804 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 0.3 cond F1 kmeans | 0.796 | 0.796 | 0.796 | changed below print precision |
| panels.tab_twopanel_MAG | 1044-1048 | 0.3 cond F1 random | 0.795 | 0.795 | 0.795 | changed below print precision |
| panels.tab_emptied_shares | 1066-1081 | original 50-99 share | 3.8 | 3.8 | 3.8 | changed below print precision |
| panels.tab_emptied_shares | 1066-1081 | kmeans 0-19 share | 44.0 | 44.0 | 44.0 | changed below print precision |
| panels.tab_emptied_shares | 1066-1081 | kmeans 20-49 share | 17.9 | 17.9 | 17.9 | changed below print precision |
| panels.tab_emptied_shares | 1066-1081 | kmeans 50-99 share | 15.2 | 15.2 | 15.2 | changed below print precision |
| panels.tab_emptied_shares | 1066-1081 | kmeans 100 share | 22.9 | 22.9 | 22.9 | changed below print precision |
| panels.tab_emptied_shares | 1066-1081 | random 20-49 share | 14.8 | 14.8 | 14.8 | changed below print precision |
| panels.tab_emptied_shares | 1066-1081 | random 50-99 share | 14.5 | 14.5 | 14.5 | changed below print precision |
| panels.cond_f1_random_gap | 1097-1098 | phystd random-original cond F1 (max |.|, V7-V9) | 0.028 | 0.028 | 0.028 | changed below print precision |
| panels.tab_threepanel_genecap_block | 1120-1124 | FullLength16S delta cond F1 CI lo | -0.031 | -0.031 | -0.031 | changed below print precision |
| panels.tab_threepanel_genecap_block | 1120-1124 | V4 delta cond F1 | -0.024 | -0.024 | -0.024 | changed below print precision |
| panels.tab_threepanel_genecap_block | 1120-1124 | V7-V9 delta cond F1 | -0.028 | -0.028 | -0.028 | changed below print precision |
| panels.tab_threepanel_genecap_block | 1120-1124 | V7-V9 delta cond F1 CI lo | -0.035 | -0.035 | -0.035 | changed below print precision |
| panels.tab_threepanel_genecap_block | 1120-1124 | V7-V9 delta cond F1 CI hi | -0.020 | -0.020 | -0.020 | changed below print precision |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 1 phystd cond F1 original | 0.809 | 0.809 | 0.809 | changed below print precision |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 1 phystd cond F1 random | 0.783 | 0.783 | 0.783 | changed below print precision |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 1 delta CI hi | -0.020 | -0.020 | -0.020 | changed below print precision |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 0.3 phystd cond F1 original | 0.806 | 0.806 | 0.806 | changed below print precision |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 0.3 phystd cond F1 random | 0.789 | 0.789 | 0.789 | changed below print precision |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 0.3 delta random-original | -0.017 | -0.017 | -0.017 | changed below print precision |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 0.3 delta CI lo | -0.023 | -0.023 | -0.023 | changed below print precision |
| panels.tab_threepanel_magpanel_block | 1130-1133 | 0.3 delta CI hi | -0.011 | -0.011 | -0.011 | changed below print precision |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus coverage V3-V4 CI hi | -0.215 | -0.215 | -0.215 | changed below print precision |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus coverage V4 CI lo | -0.276 | -0.276 | -0.276 | changed below print precision |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus coverage V7-V9 CI lo | -0.277 | -0.277 | -0.277 | changed below print precision |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | genus coverage V7-V9 CI hi | -0.237 | -0.237 | -0.237 | changed below print precision |
| panels.selfincl_genecap_penalty | 1149-1150 | V4 random phystd | 0.730 | 0.730 | 0.730 | changed below print precision |
| panels.selfincl_genecap_penalty | 1149-1150 | V7-V9 original phystd | 0.906 | 0.906 | 0.906 | changed below print precision |
| panels.selfincl_genecap_penalty | 1149-1150 | V7-V9 random phystd | 0.805 | 0.805 | 0.805 | changed below print precision |
| panels.selfincl_genecap_penalty | 1149-1150 | FL delta | -0.016 | -0.016 | -0.016 | changed below print precision |
| panels.what_this_licenses | 1160-1163 | genecap coverage random-original (max |.|) | 0.29 | 0.29 | 0.29 | changed below print precision |
| panels.what_this_licenses | 1160-1163 | max |pairwise cond F1 delta| (random V7-V9) | 0.028 | 0.028 | 0.028 | changed below print precision |
| abs-mag-vs-16s | 63 | positive in every region at every level | True | True | True |  |
| abs-threepanel-genecap-coverage | 70 | random-original genus coverage (min |.|, V1-V2: the abstract's '0.20') | 0.20 | 0.20 | 0.20 |  |
| tab:emptied | 1065-1080 | original 100 uncovered | 1.000 | 1.000 | 1.000 |  |
| tab:emptied | 1065-1080 | kmeans 100 uncovered | 1.000 | 1.000 | 1.000 |  |
| tab:emptied | 1065-1080 | random 100 uncovered | 1.000 | 1.000 | 1.000 |  |
| emptied_vs_uncovered_panels | 742 | FL emptied k-means | 22.9 | 22.9 | 22.9 |  |
| abstract_mag_margin | 63 | positive in every region at every level | True | True | True |  |
| twopanel_mag_abstention_strata | 1085-1087 | >=20 members abstention k-means | 0.559 | 0.559 | 0.559 |  |
| panels.fl_gap_standardized | 1055-1056 | common phyla | 113 | 113 | 113 |  |
| panels.emptied_arithmetic | 1061-1063 | emptied k-means | 22.9 | 22.9 | 22.9 |  |
| panels.emptied_arithmetic | 1061-1063 | 1-emptied k-means | 0.771 | 0.771 | 0.771 |  |
| panels.tab_emptied_uncovered | 1066-1081 | original 100 uncovered | 1.000 | 1.000 | 1.000 |  |
| panels.tab_emptied_uncovered | 1066-1081 | kmeans 100 uncovered | 1.000 | 1.000 | 1.000 |  |
| panels.tab_emptied_uncovered | 1066-1081 | random 100 uncovered | 1.000 | 1.000 | 1.000 |  |
| panels.std_gap_share | 1095-1096 | standardized / crude gap (low) | 76 | 76 | 76 |  |
| panels.std_gap_share | 1095-1096 | standardized intervals excluding zero (of 16) | 16 | 16 | 16 |  |
| panels.cond_f1_random_gap | 1097-1098 | phystd random-original cond F1 (min |.|, V1-V2) | 0.013 | 0.013 | 0.013 |  |
| panels.what_this_licenses | 1160-1163 | genus coverage random-original (min |.|) | 0.20 | 0.20 | 0.20 |  |
| panels.what_this_licenses | 1160-1163 | GTDB resolve-rate shift | 0.14 | 0.14 | 0.14 |  |
| tab:magval | 896-919 | (manual) | | | | **unchanged** (9,997 panel organisms): panel size (selection_10k), not re-scored |
| fig:magcomp | 921-929 | (manual) | | | | **figure**: F14, re-rendered in R5.2 |
| fig:magcross | 931-938 | (manual) | | | | **figure**: F15, re-rendered in R5.2 |
| magval_no_crossover | 826-830 | (manual) | | | | **unchanged** (37 kb): median retained source sequence at 1% completeness: a property of the synthetic MAGs, which were not re-synthesised |
| magval_af_gate_counterfactual | 878-880 | (manual) | | | | **not recomputed** (1,386 newly admitted incl. 19 gene-less; 26.7% above 0.710): the summary stores the non-empty count (n_newly_admitted) and per-level counts; the share above the union FL level is 369/1,391 = 26.5% new (Run 2), printed 26.7%; the 1,386 total (incl. gene-less picks) needs the per-MAG list |
| fig_magcomp_magcross | 921-938 | (manual) | | | | **figure**: F14/F15, re-rendered in R5.2 |
| panels.condF1_crude_bounds | 1018-1019 | (manual) | | | | **partly** ('original ahead in all eight regions'): amplicon side: see tab:threepanel genecap_cond_f1 crude rows (R2.3 items) |
| panels.taxonomy_bootstrap_CIs | 1094, 1095, 1097, 1115-1119 | (manual) | | | | **note**: region-order assumption of the tab:threepanel genus block; taxonomy point estimates are unchanged, CIs move by bootstrap-stream noise (R2.3) |

## R1.3U

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| abs-selfincl-genecap | 64 | include capture (min) | 99 | 99.148 | 99.157 | **CHANGED** |
| abs-selfexcl-genecap | 65 | population capture (min) | 58 | 58 | 59 | **CHANGED** |
| abs-selfexcl-genecap | 65 | population capture (max) | 64 | 64 | 65 | **CHANGED** |
| tab:main-genecap | 420-442 | include (min) | 99.1 | 99.1 | 99.2 | **CHANGED** |
| tab:main-genecap | 420-442 | excl pop (min) | 57.6 | 57.6 | 58.6 | **CHANGED** |
| tab:main-genecap | 420-442 | excl pop (max) | 64.4 | 64.4 | 65.2 | **CHANGED** |
| fig:genecap | 477-483 | pop (min) | 58 | 58 | 59 | **CHANGED** |
| fig:genecap | 477-483 | pop (max) | 64 | 64 | 65 | **CHANGED** |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | include V4 | 99.1 | 99.1 | 99.2 | **CHANGED** |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | excl pop FullLength16S | 63.3 | 63.3 | 64.0 | **CHANGED** |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | excl pop V4 | 58.1 | 58.1 | 59.3 | **CHANGED** |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | excl pop V4-V5 | 57.6 | 57.6 | 58.6 | **CHANGED** |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | excl pop V3-onward (max) | 59.2 | 59.2 | 60.3 | **CHANGED** |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | excl pop FL/V1 (min) | 63.3 | 63.3 | 64.0 | **CHANGED** |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | excl pop FL/V1 (max) | 64.4 | 64.4 | 65.2 | **CHANGED** |
| genecap_leakage_gap | 404 (30-33 points; 29.6 V1-V2 to 32.9 V4-V5), 489 (fig:leakage caption, capture gap) | gap V4-V5 | 32.9 | 32.9 | 33.0 | **CHANGED** |
| methods_union_size_amplicon | 178 | include ok share n_selected>1 | 34.0 | 34.0 | 34.1 | **CHANGED** |
| disc.region_tradeoff_genecap | 1175 | include capture (min) | 99 | 99.148 | 99.157 | **CHANGED** |
| disc.region_tradeoff_genecap | 1175 | excl pop V3-onward (min) | 57.6 | 57.6 | 58.6 | **CHANGED** |
| disc.region_tradeoff_genecap | 1175 | excl pop V3-onward (max) | 59.2 | 59.2 | 60.3 | **CHANGED** |
| disc.region_tradeoff_genecap | 1175 | excl pop FL/V1 (min) | 63.3 | 63.3 | 64.0 | **CHANGED** |
| disc.region_tradeoff_genecap | 1175 | excl pop FL/V1 (max) | 64.4 | 64.4 | 65.2 | **CHANGED** |
| abs-selfincl-genecap | 64 | V4 include precision | 0.79 | 0.79 | 0.79 | changed below print precision |
| abs-selfexcl-genecap | 65 | resolved capture (min) | 66 | 66 | 66 | changed below print precision |
| abs-selfexcl-genecap | 65 | resolved capture (max) | 70 | 70 | 70 | changed below print precision |
| abs-tier-genecap | 73 | family-tier mean recall | 58 | 58 | 58 | changed below print precision |
| abs-tier-genecap | 73 | genus-tier mean recall | 69 | 69 | 69 | changed below print precision |
| abs-tier-genecap | 73 | family-tier by region (min) | 56 | 56 | 56 | changed below print precision |
| abs-tier-genecap | 73 | family-tier by region (max) | 61 | 61 | 61 | changed below print precision |
| tab:main-genecap | 420-442 | excl cond (min) | 66.4 | 66.4 | 66.4 | changed below print precision |
| tab:main-genecap | 420-442 | excl cond (max) | 69.9 | 69.9 | 69.9 | changed below print precision |
| fig:genecap | 477-483 | resolved (min) | 66 | 66 | 66 | changed below print precision |
| fig:genecap | 477-483 | resolved (max) | 70 | 70 | 70 | changed below print precision |
| fig:leakage-genecap | 485-491 | gap (min) | 30 | 30 | 30 | changed below print precision |
| fig:leakage-genecap | 485-491 | gap (max) | 33 | 33 | 33 | changed below print precision |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | excl cond V1-V2 | 69.9 | 69.9 | 69.9 | changed below print precision |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | excl cond V1-V3 | 69.5 | 69.5 | 69.5 | changed below print precision |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | excl cond FullLength16S | 68.7 | 68.7 | 68.7 | changed below print precision |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | excl cond V4-V5 | 66.4 | 66.4 | 66.4 | changed below print precision |
| genecap_leakage_gap | 404 (30-33 points; 29.6 V1-V2 to 32.9 V4-V5), 489 (fig:leakage caption, capture gap) | gap V1-V2 | 29.6 | 29.6 | 29.6 | changed below print precision |
| genecap_leakage_gap | 404 (30-33 points; 29.6 V1-V2 to 32.9 V4-V5), 489 (fig:leakage caption, capture gap) | gap (min) | 30 | 30 | 30 | changed below print precision |
| genecap_leakage_gap | 404 (30-33 points; 29.6 V1-V2 to 32.9 V4-V5), 489 (fig:leakage caption, capture gap) | gap (max) | 33 | 33 | 33 | changed below print precision |
| genecap_by_tier | 73, 409 | family-tier mean recall | 58 | 58 | 58 | changed below print precision |
| genecap_by_tier | 73, 409 | genus-tier mean recall | 69 | 69 | 69 | changed below print precision |
| genecap_by_tier | 73, 409 | family-tier by region (min) | 56 | 56 | 56 | changed below print precision |
| genecap_by_tier | 73, 409 | family-tier by region (max) | 61 | 61 | 61 | changed below print precision |
| disc.region_tradeoff_genecap | 1175 | include precision V3-onward (min) | 0.79 | 0.79 | 0.79 | changed below print precision |
| disc.region_tradeoff_genecap | 1175 | include precision V3-onward (max) | 0.88 | 0.88 | 0.88 | changed below print precision |
| disc.region_tradeoff_genecap | 1175 | include precision V1-V2/V1-V3 (min) | 0.92 | 0.92 | 0.92 | changed below print precision |
| disc.region_tradeoff_genecap | 1175 | include precision V1-V2/V1-V3 (max) | 0.94 | 0.94 | 0.94 | changed below print precision |
| tab:main-genecap | 420-442 | include (max) | 99.98 | 99.98 | 99.98 |  |
| tabmain_genecap | 64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%) | include FullLength16S | 99.98 | 99.98 | 99.98 |  |
| fig:genecap | 477-483 | (manual) | | | | **figure**: F4b, re-rendered in R5.2 |

## R2.5

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| abs-depth-ec | 67 | d20 cond EC F1 | 0.798 | 0.799 | 0.798 | **CHANGED** |
| abs-depth-ec | 67 | probe recovered EC F1 | 0.783 | 0.783 | 0.829 | **CHANGED** |
| abs-depth-ec | 67 | probe recovered precision | 0.747 | 0.747 | 0.791 | **CHANGED** |
| tab:picrust2 | 542-568 | All prFBA excl coverage | 0.889 | 0.889 | 0.902 | **CHANGED** |
| tab:picrust2 | 542-568 | All prFBA excl pop F1 | 0.713 | 0.713 | 0.723 | **CHANGED** |
| tab:picrust2 | 542-568 | All prFBA incl precision | 0.936 | 0.936 | 0.935 | **CHANGED** |
| tab:picrust2 | 542-568 | n scoreable | 67,038 | 67,038 | 67,713 | **CHANGED** |
| tab:picrust2 | 542-568 | ceiling n genomes | 2,689 | 2,689 | 2,719 | **CHANGED** |
| depth_probe_quality | 67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204 | recovered EC recall | 0.854 | 0.854 | 0.904 | **CHANGED** |
| depth_probe_quality | 67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204 | recovered EC precision | 0.747 | 0.747 | 0.791 | **CHANGED** |
| depth_probe_quality | 67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204 | recovered EC F1 | 0.783 | 0.783 | 0.829 | **CHANGED** |
| depth_probe_quality | 67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204 | d20 answered recall | 0.835 | 0.835 | 0.834 | **CHANGED** |
| depth_probe_quality | 67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204 | d20 answered F1 | 0.798 | 0.799 | 0.798 | **CHANGED** |
| conditional_ec_d20_vs_d100 | 67 (0.802 vs 0.798), 732 (F1 0.802 vs 0.798; precision 0.787 vs 0.790; recall 0.847 vs 0.835), 1204 | d20 F1 | 0.798 | 0.799 | 0.798 | **CHANGED** |
| conditional_ec_d20_vs_d100 | 67 (0.802 vs 0.798), 732 (F1 0.802 vs 0.798; precision 0.787 vs 0.790; recall 0.847 vs 0.835), 1204 | d20 recall | 0.835 | 0.835 | 0.834 | **CHANGED** |
| population_ec_and_projections | 733 (0.600 -> 0.713 vs PICRUSt2 0.729), 734-735 (projections: coverage ~0.96, population genus ~0.49, population EC F1 ~0.75) | d20 pop EC F1 | 0.600 | 0.600 | 0.609 | **CHANGED** |
| population_ec_and_projections | 733 (0.600 -> 0.713 vs PICRUSt2 0.729), 734-735 (projections: coverage ~0.96, population genus ~0.49, population EC F1 ~0.75) | d100 pop EC F1 | 0.713 | 0.713 | 0.723 | **CHANGED** |
| AS-07 | 518 (2,689 scoreable; 19,686 amplicons); 524 (29% of amplicons); 523; 536; 548; 560 (n=19,686); 563 (n=47,352); 81; 1217 | in-reference scoreable genomes | 2,689 | 2,689 | 2,719 | **CHANGED** |
| AS-07 | 518 (2,689 scoreable; 19,686 amplicons); 524 (29% of amplicons); 523; 536; 548; 560 (n=19,686); 563 (n=47,352); 81; 1217 | in-reference amplicons | 19,686 | 19,686 | 19,843 | **CHANGED** |
| AS-07 | 518 (2,689 scoreable; 19,686 amplicons); 524 (29% of amplicons); 523; 536; 548; 560 (n=19,686); 563 (n=47,352); 81; 1217 | source-absent amplicons | 47,352 | 47,352 | 47,870 | **CHANGED** |
| AS-19 | 536; 548; 565 (tab:picrust2 ceiling row); 1189 ('90.4%') | ceiling Jaccard | 0.688 | 0.688 | 0.687 | **CHANGED** |
| XR-02 | 733 ('population EC F1 rises from 0.600 to 0.713, still below PICRUSt2's 0.729'); 1187 ('-0.129 at the earlier record depth of 20') | d20 pop EC F1 | 0.600 | 0.600 | 0.609 | **CHANGED** |
| XR-02 | 733 ('population EC F1 rises from 0.600 to 0.713, still below PICRUSt2's 0.729'); 1187 ('-0.129 at the earlier record depth of 20') | d100 pop EC F1 | 0.713 | 0.713 | 0.723 | **CHANGED** |
| XR-02 | 733 ('population EC F1 rises from 0.600 to 0.713, still below PICRUSt2's 0.729'); 1187 ('-0.129 at the earlier record depth of 20') | d20 population gap | -0.129 | -0.129 | -0.120 | **CHANGED** |
| lim.depth_probe_quality | 1204 | recovered EC F1 | 0.783 | 0.783 | 0.829 | **CHANGED** |
| lim.depth20_vs_100 | 1204 | d20 cond EC F1 | 0.798 | 0.799 | 0.798 | **CHANGED** |
| abs-depth-ec | 67 | d100 cond EC F1 | 0.802 | 0.802 | 0.802 | changed below print precision |
| abs-depth-ec | 67 | d20 answered precision | 0.790 | 0.790 | 0.790 | changed below print precision |
| tab:picrust2 | 542-568 | All PICRUSt2 recall | 0.849 | 0.849 | 0.849 | changed below print precision |
| tab:picrust2 | 542-568 | All PICRUSt2 precision | 0.669 | 0.669 | 0.669 | changed below print precision |
| tab:picrust2 | 542-568 | All PICRUSt2 F1 | 0.729 | 0.729 | 0.729 | changed below print precision |
| tab:picrust2 | 542-568 | All PICRUSt2 pop F1 | 0.729 | 0.729 | 0.729 | changed below print precision |
| tab:picrust2 | 542-568 | All prFBA excl recall | 0.847 | 0.847 | 0.847 | changed below print precision |
| tab:picrust2 | 542-568 | All prFBA excl precision | 0.787 | 0.787 | 0.787 | changed below print precision |
| tab:picrust2 | 542-568 | All prFBA excl F1 | 0.802 | 0.802 | 0.802 | changed below print precision |
| tab:picrust2 | 542-568 | All prFBA incl F1 | 0.959 | 0.959 | 0.959 | changed below print precision |
| tab:picrust2 | 542-568 | All prFBA incl pop F1 | 0.959 | 0.959 | 0.959 | changed below print precision |
| tab:picrust2 | 542-568 | ceiling precision | 0.743 | 0.743 | 0.743 | changed below print precision |
| tab:picrust2 | 542-568 | in-reference amplicons more precise than the reference annotation | 43 | 43 | 43 | changed below print precision |
| depth_probe_quality | 67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204 | d20 answered precision | 0.790 | 0.790 | 0.790 | changed below print precision |
| conditional_ec_d20_vs_d100 | 67 (0.802 vs 0.798), 732 (F1 0.802 vs 0.798; precision 0.787 vs 0.790; recall 0.847 vs 0.835), 1204 | d100 F1 | 0.802 | 0.802 | 0.802 | changed below print precision |
| conditional_ec_d20_vs_d100 | 67 (0.802 vs 0.798), 732 (F1 0.802 vs 0.798; precision 0.787 vs 0.790; recall 0.847 vs 0.835), 1204 | d100 precision | 0.787 | 0.787 | 0.787 | changed below print precision |
| conditional_ec_d20_vs_d100 | 67 (0.802 vs 0.798), 732 (F1 0.802 vs 0.798; precision 0.787 vs 0.790; recall 0.847 vs 0.835), 1204 | d20 precision | 0.790 | 0.790 | 0.790 | changed below print precision |
| conditional_ec_d20_vs_d100 | 67 (0.802 vs 0.798), 732 (F1 0.802 vs 0.798; precision 0.787 vs 0.790; recall 0.847 vs 0.835), 1204 | d100 recall | 0.847 | 0.847 | 0.847 | changed below print precision |
| population_ec_and_projections | 733 (0.600 -> 0.713 vs PICRUSt2 0.729), 734-735 (projections: coverage ~0.96, population genus ~0.49, population EC F1 ~0.75) | PICRUSt2 pop EC F1 | 0.729 | 0.729 | 0.729 | changed below print precision |
| AS-07 | 518 (2,689 scoreable; 19,686 amplicons); 524 (29% of amplicons); 523; 536; 548; 560 (n=19,686); 563 (n=47,352); 81; 1217 | in-reference share | 29 | 29 | 29 | changed below print precision |
| AS-19 | 536; 548; 565 (tab:picrust2 ceiling row); 1189 ('90.4%') | ceiling precision | 74.3 | 74.3 | 74.3 | changed below print precision |
| AS-19 | 536; 548; 565 (tab:picrust2 ceiling row); 1189 ('90.4%') | in-ref recall / ceiling recall | 96 | 96 | 96 | changed below print precision |
| AS-19 | 536; 548; 565 (tab:picrust2 ceiling row); 1189 ('90.4%') | PICRUSt2 in-ref recall | 0.868 | 0.868 | 0.868 | changed below print precision |
| AS-19 | 536; 548; 565 (tab:picrust2 ceiling row); 1189 ('90.4%') | in-reference amplicons more precise than the reference annotation | 43 | 43 | 43 | changed below print precision |
| XR-02 | 733 ('population EC F1 rises from 0.600 to 0.713, still below PICRUSt2's 0.729'); 1187 ('-0.129 at the earlier record depth of 20') | PICRUSt2 pop EC F1 | 0.729 | 0.729 | 0.729 | changed below print precision |
| lim.depth20_vs_100 | 1204 | d100 cond EC F1 | 0.802 | 0.802 | 0.802 | changed below print precision |
| tab:picrust2 | 542-568 | All PICRUSt2 coverage | 1.000 | 1.000 | 1.000 |  |
| tab:picrust2 | 542-568 | All prFBA incl coverage | 1.000 | 1.000 | 1.000 |  |
| tab:picrust2 | 542-568 | All prFBA incl recall | 0.999 | 0.999 | 0.999 |  |
| tab:picrust2 | 542-568 | ceiling recall | 0.904 | 0.904 | 0.904 |  |
| depth_probe_quality | 67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204 | recovered family-correct | 80.2 | 80.2 | 80.2 |  |
| depth_probe_quality | 67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204 | recovered genus-correct | 60.1 | 60.1 | 60.1 |  |
| AS-03 | 515 | shared EC universe | 2,132 | 2,132 | 2,132 |  |
| AS-03 | 515 | shared rxn universe | 7,042 | 7,042 | 7,042 |  |
| AS-19 | 536; 548; 565 (tab:picrust2 ceiling row); 1189 ('90.4%') | ceiling recall | 90.4 | 90.4 | 90.4 |  |
| lim.depth_probe_quality | 1204 | recovered genus-correct | 60.1 | 60.1 | 60.1 |  |
| lim.depth_probe_quality | 1204 | recovered family-correct | 80.2 | 80.2 | 80.2 |  |
| depth_probe_quality | 67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204 | (manual) | | | | **not recomputed** (0.651 / 0.800 answered-set anchor accuracy): no stored key or surviving recipe located (not in depth_probe_results.json or any summary); the d20 exclude-mode taxonomy is unchanged by the re-score (R1.4), so the values are expected unchanged, but this is not verified |
| population_ec_and_projections | 733 (0.600 -> 0.713 vs PICRUSt2 0.729), 734-735 (projections: coverage ~0.96, population genus ~0.49, population EC F1 ~0.75) | (manual) | | | | **not recomputed** (projections ~0.96 / ~0.49 / ~0.75): hand projections in picrust2_comparison/README.md from the probe's per-region depth-500 rates; the probe recovery counts are unchanged (283/320) but the EC quality moved (0.783 -> 0.829), so the ~0.75 EC projection should be re-derived by the manuscript agent or dropped |
| AS-07 | 518 (2,689 scoreable; 19,686 amplicons); 524 (29% of amplicons); 523; 536; 548; 560 (n=19,686); 563 (n=47,352); 81; 1217 | (manual) | | | | **unchanged** (2,720 of 9,746 source genomes in PICRUSt2's reference): the source-genome set and PICRUSt2's reference are fixed: bvbrc_source_genomes.json byte-equal after R2.4 and annotation_ceiling.json _meta.n_source_assembly_in_picrust2_ref = 2,720 of 9,746 in both |
| lim.depth20_vs_100 | 1204 | (manual) | | | | **see item** (legacy genus within 0.01, recall +0.7..+1.3, F1 -0.1..-1.2 points): values under val16s:legacy_genecap_d20_vs_d100 (R1.4) |

## R2.4

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| abs-picrust2-asshipped | 74 | no EC prediction share | 11.1 | 11.1 | 9.8 | **CHANGED** |
| abs-picrust2-asshipped | 74 | abstain share | 9.6 | 9.6 | 9.5 | **CHANGED** |
| fig:picrust2 | 570-579 | PICRUSt2 precision NSTI>=0.5 | 0.235 | 0.235 | 0.234 | **CHANGED** |
| fig:picrust2 | 570-579 | prFBA precision NSTI>=0.5 | 0.606 | 0.606 | 0.605 | **CHANGED** |
| fig:picrust2 | 570-579 | family-floor bin n | 547 | 547 | 548 | **CHANGED** |
| fig:picrust2 | 570-579 | prFBA EC non-answers | 7,422 | 7,422 | 6,660 | **CHANGED** |
| fig:picrust2 | 570-579 | outright abstentions (no reps) | 6,407 | 6,407 | 6,423 | **CHANGED** |
| fig:picrust2 | 570-579 | floor bin share of non-answers | 7 | 7 | 8 | **CHANGED** |
| fig:picrust2 | 570-579 | further amplicons with surviving hits all below the floor | 145 | 145 | 147 | **CHANGED** |
| fig:picrust2 | 570-579 | abstentions the floor explains | 692 | 692 | 695 | **CHANGED** |
| ec_coverage_bins_and_noanswer | 74 (11.1%; 9.6%), 715 (EC coverage 0.982 and 0.829), 1179 (answers 89%), 1181 (11%; 9.6%; 1.5%) | coverage 0.945-0.987 | 0.982 | 0.982 | 0.993 | **CHANGED** |
| ec_coverage_bins_and_noanswer | 74 (11.1%; 9.6%), 715 (EC coverage 0.982 and 0.829), 1179 (answers 89%), 1181 (11%; 9.6%; 1.5%) | coverage >=0.995 | 0.829 | 0.829 | 0.841 | **CHANGED** |
| ec_coverage_bins_and_noanswer | 74 (11.1%; 9.6%), 715 (EC coverage 0.982 and 0.829), 1179 (answers 89%), 1181 (11%; 9.6%; 1.5%) | answers | 89 | 89 | 90 | **CHANGED** |
| ec_coverage_bins_and_noanswer | 74 (11.1%; 9.6%), 715 (EC coverage 0.982 and 0.829), 1179 (answers 89%), 1181 (11%; 9.6%; 1.5%) | no EC prediction | 11.1 | 11.1 | 9.8 | **CHANGED** |
| ec_coverage_bins_and_noanswer | 74 (11.1%; 9.6%), 715 (EC coverage 0.982 and 0.829), 1179 (answers 89%), 1181 (11%; 9.6%; 1.5%) | abstain | 9.6 | 9.6 | 9.5 | **CHANGED** |
| ec_coverage_bins_and_noanswer | 74 (11.1%; 9.6%), 715 (EC coverage 0.982 and 0.829), 1179 (answers 89%), 1181 (11%; 9.6%; 1.5%) | reps without gene families | 1.5 | 1.5 | 0.4 | **CHANGED** |
| AS-04 | 516; 544 (n=67,038); 1179 (67,038) | truth_empty amplicons (prfba_exclude ec) | 953 | 953 | 278 | **CHANGED** |
| AS-04 | 516; 544 (n=67,038); 1179 (67,038) | scoreable n | 67,038 | 67,038 | 67,713 | **CHANGED** |
| AS-04 | 516; 544 (n=67,038); 1179 (67,038) | truth_empty source genomes | 161 | 161 | 38 | **CHANGED** |
| AS-05 | 517 | truth_empty (stated in the disclosure) | 953 | 953 | 278 | **CHANGED** |
| AS-05 | 517 | pred_empty | 1,060 | 1,060 | 243 | **CHANGED** |
| AS-08 | 521; 1179 | PICRUSt2 below alignment threshold (scoreable) | 26 | 26 | 27 | **CHANGED** |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | PICRUSt2 answered n | 67,012 | 67,012 | 67,686 | **CHANGED** |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | prFBA excl coverage | 88.9 | 88.9 | 90.2 | **CHANGED** |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | prFBA excl pop F1 | 0.713 | 0.713 | 0.723 | **CHANGED** |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | prFBA incl precision | 0.936 | 0.936 | 0.935 | **CHANGED** |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | population delta | -0.016 | -0.016 | -0.006 | **CHANGED** |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | legacy reducer population EC F1, as shipped (depth 100) | None | 0.726 | 0.734 | **CHANGED** |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | legacy reducer population EC F1, as shipped (depth 20) | None | 0.608 | 0.616 | **CHANGED** |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | legacy reducer minus PICRUSt2, as-shipped population EC F1 (depth 100) | None | -0.002 | +0.005 | **CHANGED**; MEANING: sign flip |
| AS-11 | 74; 1181 | no EC prediction | 11.1 | 11.1 | 9.8 | **CHANGED** |
| AS-11 | 74; 1181 | abstain | 9.6 | 9.6 | 9.5 | **CHANGED** |
| AS-11 | 74; 1181 | reps without gene families | 1.5 | 1.5 | 0.4 | **CHANGED** |
| AS-12 | 521; 576 (fig C) | median dF1 (low) | +0.072 | +0.072 | +0.073 | **CHANGED** |
| AS-13 | 81 ('+0.046 to -0.202'); 523; 1217 | conditional dF1 vs BV-BRC truth | +0.046 | +0.046 | +0.047 | **CHANGED** |
| AS-13 | 81 ('+0.046 to -0.202'); 523; 1217 | conditional dF1 vs PICRUSt2 truth | -0.202 | -0.202 | -0.203 | **CHANGED** |
| AS-13 | 81 ('+0.046 to -0.202'); 523; 1217 | amplicons in both databases | 19,686 | 19,686 | 19,843 | **CHANGED** |
| AS-14 | 82 ('0.768 for full-length 16S against 0.668 for V4'); 526 | prFBA pop F1 FullLength16S | 0.768 | 0.768 | 0.777 | **CHANGED** |
| AS-14 | 82 ('0.768 for full-length 16S against 0.668 for V4'); 526 | prFBA pop F1 V1-V2 | 0.760 | 0.760 | 0.768 | **CHANGED** |
| AS-14 | 82 ('0.768 for full-length 16S against 0.668 for V4'); 526 | prFBA pop F1 V1-V3 | 0.766 | 0.766 | 0.774 | **CHANGED** |
| AS-14 | 82 ('0.768 for full-length 16S against 0.668 for V4'); 526 | prFBA pop F1 V4 | 0.668 | 0.668 | 0.681 | **CHANGED** |
| AS-14 | 82 ('0.768 for full-length 16S against 0.668 for V4'); 526 | PICRUSt2 pop F1 V1-V2 | 0.731 | 0.731 | 0.732 | **CHANGED** |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | PICRUSt2 precision >=0.5 | 0.235 | 0.235 | 0.234 | **CHANGED** |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | PICRUSt2 recall >=0.5 | 0.769 | 0.769 | 0.768 | **CHANGED** |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | prFBA precision >=0.5 | 0.606 | 0.606 | 0.605 | **CHANGED** |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | prFBA coverage >=0.5 | 86.0 | 86.0 | 89.2 | **CHANGED** |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | prFBA coverage <0.01 | 92.1 | 92.1 | 93.5 | **CHANGED** |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | prFBA answered at >=0.5 | 603 | 603 | 629 | **CHANGED** |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | PICRUSt2 precision on the 603 | 0.230 | 0.230 | 0.231 | **CHANGED** |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | PICRUSt2 ECs on the same amplicons | 412 | 412 | 410 | **CHANGED** |
| AS-17 | 531 | high-NSTI answered with a >=0.995 non-self in the record | 248 | 248 | 253 | **CHANGED** |
| AS-17 | 531 | high-NSTI answered | 603 | 603 | 629 | **CHANGED** |
| AS-18 | 532; 575 (fig B: 0.32; 7% = 547 of 7,422; 145 further; 692 = 10.8% of 6,407); 1179 ('none at all below its family floor') | bin n | 547 | 547 | 548 | **CHANGED** |
| AS-18 | 532; 575 (fig B: 0.32; 7% = 547 of 7,422; 145 further; 692 = 10.8% of 6,407); 1179 ('none at all below its family floor') | EC non-answers | 7,422 | 7,422 | 6,660 | **CHANGED** |
| AS-18 | 532; 575 (fig B: 0.32; 7% = 547 of 7,422; 145 further; 692 = 10.8% of 6,407); 1179 ('none at all below its family floor') | outright abstentions | 6,407 | 6,407 | 6,423 | **CHANGED** |
| AS-18 | 532; 575 (fig B: 0.32; 7% = 547 of 7,422; 145 further; 692 = 10.8% of 6,407); 1179 ('none at all below its family floor') | further amplicons with surviving hits all below the floor | 145 | 145 | 147 | **CHANGED** |
| AS-18 | 532; 575 (fig B: 0.32; 7% = 547 of 7,422; 145 further; 692 = 10.8% of 6,407); 1179 ('none at all below its family floor') | abstentions the floor explains | 692 | 692 | 695 | **CHANGED** |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref prFBA coverage | 0.932 | 0.932 | 0.943 | **CHANGED** |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref prFBA pop F1 | 0.781 | 0.781 | 0.791 | **CHANGED** |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent prFBA coverage | 0.872 | 0.872 | 0.884 | **CHANGED** |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent prFBA recall | 0.837 | 0.837 | 0.836 | **CHANGED** |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent prFBA F1 | 0.786 | 0.786 | 0.785 | **CHANGED** |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent prFBA pop F1 | 0.685 | 0.685 | 0.695 | **CHANGED** |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref n | 19,686 | 19,686 | 19,843 | **CHANGED** |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent n | 47,352 | 47,352 | 47,870 | **CHANGED** |
| AS-21 | 539 | include-self EC Jaccard | 0.94 | 0.94 | 0.93 | **CHANGED** |
| AS-21 | 539 | ceiling Jaccard | 0.688 | 0.688 | 0.687 | **CHANGED** |
| XR-01 | 715 ('the EC coverage of fig:picrust2B ... is 0.982 and 0.829') | coverage 0.945-0.987 | 0.982 | 0.982 | 0.993 | **CHANGED** |
| XR-01 | 715 ('the EC coverage of fig:picrust2B ... is 0.982 and 0.829') | coverage >=0.995 | 0.829 | 0.829 | 0.841 | **CHANGED** |
| disc.picrust2_failure_modes | 1179 | PICRUSt2 not answering (scoreable) | 26 | 26 | 27 | **CHANGED** |
| disc.picrust2_failure_modes | 1179 | scoreable n | 67,038 | 67,038 | 67,713 | **CHANGED** |
| disc.picrust2_failure_modes | 1179 | PICRUSt2 precision >=0.5 | 0.235 | 0.235 | 0.234 | **CHANGED** |
| disc.picrust2_failure_modes | 1179 | PICRUSt2 recall >=0.5 | 0.769 | 0.769 | 0.768 | **CHANGED** |
| disc.prfba_answer_rate | 1179-1181 | answers | 89 | 89 | 90 | **CHANGED** |
| disc.prfba_answer_rate | 1179-1181 | non-answer | 11 | 11 | 10 | **CHANGED** |
| disc.prfba_answer_rate | 1179-1181 | abstain | 9.6 | 9.6 | 9.5 | **CHANGED** |
| disc.prfba_answer_rate | 1179-1181 | reps without gene families | 1.5 | 1.5 | 0.4 | **CHANGED** |
| abs-picrust2-asshipped | 74 | cond EC F1 prFBA | 0.802 | 0.802 | 0.802 | changed below print precision |
| abs-picrust2-asshipped | 74 | matched both-answer PICRUSt2 F1 | 0.727 | 0.727 | 0.727 | changed below print precision |
| fig:picrust2 | 570-579 | PICRUSt2 precision NSTI<0.01 | 0.735 | 0.735 | 0.735 | changed below print precision |
| fig:picrust2 | 570-579 | PICRUSt2 precision below family floor | 0.32 | 0.32 | 0.32 | changed below print precision |
| fig:picrust2 | 570-579 | share of outright abstentions | 10.8 | 10.8 | 10.8 | changed below print precision |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | PICRUSt2 recall | 0.849 | 0.849 | 0.849 | changed below print precision |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | PICRUSt2 precision | 0.669 | 0.669 | 0.669 | changed below print precision |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | PICRUSt2 F1 | 0.729 | 0.729 | 0.729 | changed below print precision |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | PICRUSt2 pop F1 | 0.729 | 0.729 | 0.729 | changed below print precision |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | prFBA excl recall | 0.847 | 0.847 | 0.847 | changed below print precision |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | prFBA excl precision | 0.787 | 0.787 | 0.787 | changed below print precision |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | prFBA excl cond F1 | 0.802 | 0.802 | 0.802 | changed below print precision |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | prFBA incl F1 | 0.959 | 0.959 | 0.959 | changed below print precision |
| AS-12 | 521; 576 (fig C) | median dF1 (high) | +0.095 | +0.095 | +0.095 | changed below print precision |
| AS-12 | 521; 576 (fig C) | wins per loss (min) | 3 | 3 | 3 | changed below print precision |
| AS-12 | 521; 576 (fig C) | wins per loss (max) | 7 | 7 | 7 | changed below print precision |
| AS-14 | 82 ('0.768 for full-length 16S against 0.668 for V4'); 526 | PICRUSt2 pop F1 FullLength16S | 0.745 | 0.745 | 0.745 | changed below print precision |
| AS-14 | 82 ('0.768 for full-length 16S against 0.668 for V4'); 526 | PICRUSt2 pop F1 V1-V3 | 0.732 | 0.732 | 0.732 | changed below print precision |
| AS-15 | 527 | rxn cond F1 prFBA | 0.805 | 0.805 | 0.805 | changed below print precision |
| AS-15 | 527 | rxn cond F1 PICRUSt2 | 0.721 | 0.722 | 0.722 | changed below print precision |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | PICRUSt2 precision <0.01 | 0.735 | 0.735 | 0.735 | changed below print precision |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | PICRUSt2 recall <0.01 | 0.870 | 0.871 | 0.871 | changed below print precision |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | prFBA precision <0.01 | 0.822 | 0.822 | 0.822 | changed below print precision |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | prFBA ECs per answered amplicon | 141 | 141 | 141 | changed below print precision |
| AS-16 | 530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179 | truth ECs on the same amplicons | 123 | 123 | 123 | changed below print precision |
| AS-17 | 531 | median reference identity | 0.956 | 0.956 | 0.956 | changed below print precision |
| AS-18 | 532; 575 (fig B: 0.32; 7% = 547 of 7,422; 145 further; 692 = 10.8% of 6,407); 1179 ('none at all below its family floor') | PICRUSt2 precision in bin | 0.316 | 0.316 | 0.316 | changed below print precision |
| AS-18 | 532; 575 (fig B: 0.32; 7% = 547 of 7,422; 145 further; 692 = 10.8% of 6,407); 1179 ('none at all below its family floor') | share of outright abstentions | 10.8 | 10.8 | 10.8 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref PICRUSt2 coverage | 1.000 | 1.000 | 1.000 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref PICRUSt2 recall | 0.868 | 0.868 | 0.868 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref PICRUSt2 precision | 0.735 | 0.735 | 0.735 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref PICRUSt2 F1 | 0.792 | 0.792 | 0.792 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref prFBA recall | 0.871 | 0.871 | 0.871 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref prFBA precision | 0.828 | 0.828 | 0.828 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref prFBA F1 | 0.838 | 0.838 | 0.838 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent PICRUSt2 coverage | 0.999 | 0.999 | 0.999 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent PICRUSt2 recall | 0.840 | 0.841 | 0.841 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent PICRUSt2 precision | 0.642 | 0.642 | 0.642 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent PICRUSt2 F1 | 0.703 | 0.703 | 0.703 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent PICRUSt2 pop F1 | 0.703 | 0.703 | 0.703 | changed below print precision |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | absent prFBA precision | 0.769 | 0.769 | 0.769 | changed below print precision |
| AS-21 | 539 | agreement Jaccard (low) | 0.56 | 0.56 | 0.56 | changed below print precision |
| AS-21 | 539 | agreement Jaccard (high) | 0.61 | 0.61 | 0.61 | changed below print precision |
| AS-21 | 539 | PICRUSt2 Jaccard by region (low) | 0.58 | 0.58 | 0.58 | changed below print precision |
| AS-21 | 539 | PICRUSt2 Jaccard by region (high) | 0.61 | 0.61 | 0.61 | changed below print precision |
| disc.picrust2_failure_modes | 1179 | PICRUSt2 precision <0.01 | 0.735 | 0.735 | 0.735 | changed below print precision |
| disc.picrust2_failure_modes | 1179 | PICRUSt2 recall <0.01 | 0.870 | 0.871 | 0.871 | changed below print precision |
| AS-08 | 521; 1179 | all amplicons without PICRUSt2 prediction | 27 | 27 | 27 |  |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | PICRUSt2 coverage | 1.000 | 1.000 | 1.000 |  |
| AS-09 | 521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729') | prFBA incl recall | 0.999 | 0.999 | 0.999 |  |
| AS-12 | 521; 576 (fig C) | Holm p (max) | 1e-10 | 0.0 | 0.0 |  |
| AS-14 | 82 ('0.768 for full-length 16S against 0.668 for V4'); 526 | regions where prFBA pop F1 < PICRUSt2 (count) | 5 | 5 | 5 |  |
| AS-17 | 531 | all answers below 0.987 identity | yes | yes | yes |  |
| AS-18 | 532; 575 (fig B: 0.32; 7% = 547 of 7,422; 145 further; 692 = 10.8% of 6,407); 1179 ('none at all below its family floor') | prFBA coverage below floor | 0 | 0 | 0 |  |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | in-ref PICRUSt2 pop F1 | 0.792 | 0.792 | 0.792 |  |
| AS-20 | 537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows) | median NSTI of in-reference placements | 0.010 | 0.010 | 0.010 |  |
| AS-23 | 540 (second clause) | median self-excluded selection size | 1 | 1 | 1 |  |
| disc.prfba_answer_rate | 1179-1181 | coverage below family floor | 0 | 0 | 0 |  |
| AS-05 | 517 | (manual) | | | | **historical** (675 of 953, 818 of 1,060, 167 genomes, 0.75-0.86 / 0.05 points): the 11-13 repair measurement (CORRECTIONS_2026-11-13.json, genecap_repaired_impact/impact.json). After the re-score the reported state is the repaired one (truth_empty 278, pred_empty 243), so L517 must become a past-tense correction note (plan section 5) |
| AS-15 | 527 | (manual) | | | | **not recomputed** (0.026 archaeal, 0.014 distant novelty bins, 0.05 / 0.03 in the 9- and 8-amplicon cells): per-stratum rxn-minus-EC differences with a domain split of comparison_cells; no recipe written in this run. The overall rxn-EC tracking (0.805 / 0.721) is extracted above; the manuscript agent should re-derive the four divergences or soften them to 'small subsets' |
| AS-24 | 572-577 (fig:picrust2, images/prfba/F13_prfba_vs_picrust2.pdf) | (manual) | | | | **figure**: F13, re-rendered in R5.2 (suptitle n computed) |

## D1

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| SR-01 | 75 (6,788-7,175; median 7,175); 584; 587 (8,872 eligible; 63,871 cells); 588 (2,376 groups, n_eff 6,716, ICC 0.33); 683 (7,175; n=63,871); 1215; 1219 (7,175, 'forty-fourth'); 1220 (2,376; 6,716) | effective n | 6,716 | 6,716 | 6,715 | **CHANGED** |
| lim.holdout_power | 1220-1221 | genus n_eff | 6,716 | 6,716 | 6,715 | **CHANGED** |
| SR-01 | 75 (6,788-7,175; median 7,175); 584; 587 (8,872 eligible; 63,871 cells); 588 (2,376 groups, n_eff 6,716, ICC 0.33); 683 (7,175; n=63,871); 1215; 1219 (7,175, 'forty-fourth'); 1220 (2,376; 6,716) | ICC (group) | 0.33 | 0.33 | 0.33 | changed below print precision |
| abs-sharedref-size | 75 | median reference genomes | 7,175 | 7,175 | 7,175 |  |
| abs-sharedref-size | 75 | min reference genomes per fold | 6,788 | 6,788 | 6,788 |  |
| SR-01 | 75 (6,788-7,175; median 7,175); 584; 587 (8,872 eligible; 63,871 cells); 588 (2,376 groups, n_eff 6,716, ICC 0.33); 683 (7,175; n=63,871); 1215; 1219 (7,175, 'forty-fourth'); 1220 (2,376; 6,716) | median reference genomes | 7,175 | 7,175 | 7,175 |  |
| SR-01 | 75 (6,788-7,175; median 7,175); 584; 587 (8,872 eligible; 63,871 cells); 588 (2,376 groups, n_eff 6,716, ICC 0.33); 683 (7,175; n=63,871); 1215; 1219 (7,175, 'forty-fourth'); 1220 (2,376; 6,716) | min reference genomes per fold | 6,788 | 6,788 | 6,788 |  |
| SR-01 | 75 (6,788-7,175; median 7,175); 584; 587 (8,872 eligible; 63,871 cells); 588 (2,376 groups, n_eff 6,716, ICC 0.33); 683 (7,175; n=63,871); 1215; 1219 (7,175, 'forty-fourth'); 1220 (2,376; 6,716) | eligible genomes | 8,872 | 8,872 | 8,872 |  |
| SR-01 | 75 (6,788-7,175; median 7,175); 584; 587 (8,872 eligible; 63,871 cells); 588 (2,376 groups, n_eff 6,716, ICC 0.33); 683 (7,175; n=63,871); 1215; 1219 (7,175, 'forty-fourth'); 1220 (2,376; 6,716) | scored cells | 63,871 | 63,871 | 63,871 |  |
| SR-01 | 75 (6,788-7,175; median 7,175); 584; 587 (8,872 eligible; 63,871 cells); 588 (2,376 groups, n_eff 6,716, ICC 0.33); 683 (7,175; n=63,871); 1215; 1219 (7,175, 'forty-fourth'); 1220 (2,376; 6,716) | holdout groups | 2,376 | 2,376 | 2,376 |  |
| FR-01 | 624 (2,059 held-out genomes; 14,692 cells; median 8,468-genome reference); 627 (171 groups vs 2,376; n_eff 598 vs 6,716); 1220 | median reference genomes (family) | 8,468 | 8,468 | 8,468 |  |
| FR-01 | 624 (2,059 held-out genomes; 14,692 cells; median 8,468-genome reference); 627 (171 groups vs 2,376; n_eff 598 vs 6,716); 1220 | scored cells (family) | 14,692 | 14,692 | 14,692 |  |
| FR-01 | 624 (2,059 held-out genomes; 14,692 cells; median 8,468-genome reference); 627 (171 groups vs 2,376; n_eff 598 vs 6,716); 1220 | holdout groups (family) | 171 | 171 | 171 |  |
| FR-01 | 624 (2,059 held-out genomes; 14,692 cells; median 8,468-genome reference); 627 (171 groups vs 2,376; n_eff 598 vs 6,716); 1220 | effective n (family) | 598 | 598 | 598 |  |
| lim.database_sizes | 1214-1215, 1219 | shared reference median | 7,175 | 7,175 | 7,175 |  |
| lim.holdout_power | 1220-1221 | family holdout groups | 171 | 171 | 171 |  |
| lim.holdout_power | 1220-1221 | genus holdout groups | 2,376 | 2,376 | 2,376 |  |
| lim.holdout_power | 1220-1221 | family n_eff | 598 | 598 | 598 |  |
| lim.holdout_power | 1220-1221 | eligible | 8,872 | 8,872 | 8,872 |  |
| FR-01 | 624 (2,059 held-out genomes; 14,692 cells; median 8,468-genome reference); 627 (171 groups vs 2,376; n_eff 598 vs 6,716); 1220 | (manual) | | | | **unchanged (Path A)** (2,059 held-out genomes): fold plan kept (D1 Path A) |
| FR-04 | 627 (18.6%, 1,653 of 8,872); 1221 | (manual) | | | | **unchanged (Path A)** (1,653 of 8,872 (18.6%)): panel_eligible.json unchanged (sha256 and inode checked in R3.1) |
| lim.database_sizes | 1214-1215, 1219 | (manual) | | | | **unchanged**: database facts (~315,000 genomes, 459,301 16S, 27,870 PICRUSt2 genomes) |
| lim.holdout_power | 1220-1221 | (manual) | | | | **unchanged (Path A)** (18.6% (1,653)): panel_eligible.json unchanged |

## R3

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| abs-dialect-density | 81 | dialect conditional dF1 (BV-BRC truth) | +0.046 | +0.046 | +0.047 | **CHANGED** |
| abs-dialect-density | 81 | dialect conditional dF1 (PICRUSt2 truth) | -0.202 | -0.202 | -0.203 | **CHANGED** |
| abs-dialect-density | 81 | amplicons | 19,686 | 19,686 | 19,843 | **CHANGED** |
| abs-region-penalty | 82 | as-shipped FL pop F1 | 0.768 | 0.768 | 0.777 | **CHANGED** |
| abs-region-penalty | 82 | as-shipped V4 pop F1 | 0.668 | 0.668 | 0.681 | **CHANGED** |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | as-shipped population deficit | -0.016 | -0.016 | -0.006 | **CHANGED** |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | as-shipped prFBA pop | 0.713 | 0.713 | 0.723 | **CHANGED** |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | depth-20 population deficit | -0.129 | -0.129 | -0.120 | **CHANGED** |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | legacy reducer minus PICRUSt2, as-shipped population EC F1 (depth 100) | None | -0.002 | +0.005 | **CHANGED**; MEANING: sign flip |
| disc.vocab_breadth | 1189, 1217 | dialect +0.046 | +0.046 | +0.046 | +0.047 | **CHANGED** |
| disc.vocab_breadth | 1189, 1217 | dialect -0.202 | -0.202 | -0.202 | -0.203 | **CHANGED** |
| abs-controlled-genus | 76 | prfba_exact_tie pop F1 | 0.685 | 0.685 | 0.685 | changed below print precision |
| abs-controlled-genus | 76 | prfba_exact_tie_forced pop F1 | 0.746 | 0.746 | 0.746 | changed below print precision |
| abs-controlled-genus | 76 |   CI lo | -0.107 | -0.107 | -0.107 | changed below print precision |
| abs-controlled-genus | 76 |   CI hi | -0.081 | -0.081 | -0.081 | changed below print precision |
| abs-controlled-genus | 76 |   CI lo | -0.109 | -0.109 | -0.109 | changed below print precision |
| abs-controlled-genus | 76 |   CI hi | -0.083 | -0.083 | -0.083 | changed below print precision |
| abs-dialect-density | 81 | density slope prFBA | +0.107 | +0.107 | +0.107 | changed below print precision |
| tab:sharedref | 681-707 | picrust2 skill | 0.92 | 0.92 | 0.92 | changed below print precision |
| tab:sharedref | 681-707 | picrust2_nsti_routed cond (matched) F1 | 0.793 | 0.793 | 0.793 | changed below print precision |
| tab:sharedref | 681-707 | picrust2_nsti_routed skill | 0.90 | 0.90 | 0.90 | changed below print precision |
| tab:sharedref | 681-707 | prfba_nestedcv skill | 0.28 | 0.28 | 0.28 | changed below print precision |
| tab:sharedref | 681-707 | prfba_nestedcv_forced skill | 0.70 | 0.70 | 0.70 | changed below print precision |
| tab:sharedref | 681-707 | prfba_exact_tie pop F1 | 0.685 | 0.685 | 0.685 | changed below print precision |
| tab:sharedref | 681-707 | prfba_exact_tie_forced pop F1 | 0.746 | 0.746 | 0.746 | changed below print precision |
| tab:sharedref | 681-707 | base_nn1 skill | 0.69 | 0.69 | 0.69 | changed below print precision |
| tab:sharedref | 681-707 | base_random_abstain skill | 0.14 | 0.14 | 0.14 | changed below print precision |
| tab:sharedref | 681-707 | base_family_donor skill | -0.90 | -0.90 | -0.90 | changed below print precision |
| tab:sharedref | 681-707 | base_constant pop F1 | 0.646 | 0.646 | 0.646 | changed below print precision |
| tab:sharedref | 681-707 | base_constant cond (matched) F1 | 0.665 | 0.665 | 0.665 | changed below print precision |
| tab:sharedref | 681-707 | base_clade_max pop F1 | 0.790 | 0.790 | 0.790 | changed below print precision |
| controlled_forced_recovery | 736-737 |   CI lo | -0.037 | -0.037 | -0.037 | changed below print precision |
| controlled_forced_recovery | 736-737 |   CI hi | -0.028 | -0.028 | -0.028 | changed below print precision |
| controlled_forced_recovery | 736-737 |   CI lo | -0.039 | -0.039 | -0.039 | changed below print precision |
| controlled_forced_recovery | 736-737 |   CI hi | -0.029 | -0.029 | -0.029 | changed below print precision |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | legacy single share | 61.7 | 61.7 | 61.7 | changed below print precision |
| SR-02 | 77 ('picked legacy in every fold'); 340 (10/10, 0.001-0.003); 588-590 (legacy in all five folds at both ranks by 0.001-0.003 F1; exact-tie 0.002 lower at identical coverage) | exact-tie pop F1 | 0.685 | 0.685 | 0.685 | changed below print precision |
| SR-02 | 77 ('picked legacy in every fold'); 340 (10/10, 0.001-0.003); 588-590 (legacy in all five folds at both ranks by 0.001-0.003 F1; exact-tie 0.002 lower at identical coverage) | exact-tie minus nestedcv | 0.002 | 0.002 | 0.002 | changed below print precision |
| SR-02 | 77 ('picked legacy in every fold'); 340 (10/10, 0.001-0.003); 588-590 (legacy in all five folds at both ranks by 0.001-0.003 F1; exact-tie 0.002 lower at identical coverage) | inner margin legacy - best other (min) | 0.001 | 0.001 | 0.001 | changed below print precision |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_exact_tie pop F1 | 0.685 | 0.685 | 0.685 | changed below print precision |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_exact_tie_forced pop F1 | 0.746 | 0.746 | 0.746 | changed below print precision |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_nestedcv CI 0 | -0.107 | -0.107 | -0.107 | changed below print precision |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_nestedcv CI 1 | -0.081 | -0.081 | -0.081 | changed below print precision |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_exact_tie CI 0 | -0.109 | -0.109 | -0.109 | changed below print precision |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_exact_tie CI 1 | -0.083 | -0.083 | -0.083 | changed below print precision |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_nestedcv_forced CI 0 | -0.037 | -0.037 | -0.037 | changed below print precision |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_nestedcv_forced CI 1 | -0.028 | -0.028 | -0.028 | changed below print precision |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_exact_tie_forced CI 0 | -0.039 | -0.039 | -0.039 | changed below print precision |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_exact_tie_forced CI 1 | -0.029 | -0.029 | -0.029 | changed below print precision |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_nestedcv Archaea | -0.052 | -0.052 | -0.052 | changed below print precision |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_nestedcv_forced Archaea | -0.019 | -0.019 | -0.019 | changed below print precision |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_exact_tie_forced Archaea | -0.020 | -0.020 | -0.020 | changed below print precision |
| SR-05 | 598 (AURC 0.812 vs 0.788, exact-tie 0.787); 599 (cross-knob 0.813 and 0.788, exact-tie 0.786; family swap ~0.007) | PICRUSt2 AURC | 0.812 | 0.812 | 0.812 | changed below print precision |
| SR-05 | 598 (AURC 0.812 vs 0.788, exact-tie 0.787); 599 (cross-knob 0.813 and 0.788, exact-tie 0.786; family swap ~0.007) | prFBA cross-knob AURC | 0.788 | 0.788 | 0.788 | changed below print precision |
| SR-07 | 601-602 | PICRUSt2 rxn pop F1 | 0.780 | 0.780 | 0.780 | changed below print precision |
| SR-07 | 601-602 | exact-tie delta | -0.092 | -0.092 | -0.092 | changed below print precision |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | constant pop F1 | 0.646 | 0.646 | 0.646 | changed below print precision |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | clade-max pop F1 | 0.790 | 0.790 | 0.790 | changed below print precision |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | PICRUSt2 skill | 0.92 | 0.92 | 0.92 | changed below print precision |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | forced skill | 0.70 | 0.70 | 0.70 | changed below print precision |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | nn1 skill | 0.69 | 0.69 | 0.69 | changed below print precision |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | exact-tie forced pop F1 | 0.746 | 0.746 | 0.746 | changed below print precision |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 |   CI lo | -0.069 | -0.069 | -0.069 | changed below print precision |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 |   CI hi | -0.051 | -0.051 | -0.051 | changed below print precision |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | exact-tie vs nn1 | -0.060 | -0.060 | -0.060 | changed below print precision |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 |   CI lo | 0.030 | 0.029 | 0.029 | changed below print precision; old≠printed (rounding; known: already a rounding error before the re-score: the stored 4-dp bound 0.0295 was printed 0.030, but the unrounded bound (score_folds' bootstrap re-run, W6 sr_ci_*) is 0.02949, so 0.029; corrected in the tex 2026-11-17) |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 |   CI hi | 0.040 | 0.040 | 0.040 | changed below print precision |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') |   CI lo | 0.014 | 0.014 | 0.014 | changed below print precision |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') |   CI hi | 0.030 | 0.030 | 0.030 | changed below print precision |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') | exact-tie vs its own random mask (approx.) | +0.019 | +0.019 | +0.019 | changed below print precision |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') | legacy forced, random-masked | 0.668 | 0.668 | 0.668 | changed below print precision |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') | like-for-like credit | +0.019 | +0.019 | +0.019 | changed below print precision |
| SR-10 | 620-621 | PICRUSt2 precision 0.945-0.987 | 0.862 | 0.862 | 0.862 | changed below print precision |
| SR-10 | 620-621 | PICRUSt2 precision <0.865 | 0.680 | 0.680 | 0.680 | changed below print precision |
| SR-10 | 620-621 | prFBA precision 0.945-0.987 | 0.783 | 0.783 | 0.783 | changed below print precision |
| SR-10 | 620-621 | exact-tie precision 0.945-0.987 | 0.803 | 0.804 | 0.804 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | picrust2 skill | 0.92 | 0.92 | 0.92 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | picrust2_nsti_routed cond (matched) F1 | 0.793 | 0.793 | 0.793 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | picrust2_nsti_routed skill | 0.90 | 0.90 | 0.90 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_nestedcv skill | 0.28 | 0.28 | 0.28 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_nestedcv_forced skill | 0.70 | 0.70 | 0.70 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_exact_tie pop F1 | 0.685 | 0.685 | 0.685 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_exact_tie_forced pop F1 | 0.746 | 0.746 | 0.746 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_nn1 skill | 0.69 | 0.69 | 0.69 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_random_abstain skill | 0.14 | 0.14 | 0.14 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_family_donor skill | -0.90 | -0.90 | -0.90 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_constant pop F1 | 0.646 | 0.646 | 0.646 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_constant cond (matched) F1 | 0.665 | 0.665 | 0.665 | changed below print precision |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_clade_max pop F1 | 0.790 | 0.790 | 0.790 | changed below print precision |
| SR-12 | 81 ('+0.107, or +0.106 ... against +0.030 ... fourfold sweep'); 523 ('25, 50 and 100% subsamples ... 3.5 times as steeply (0.107, or 0.106 ... versus 0.030)'); 1189; 1217 ('quartering the shared reference costs prFBA 0.064 ... against PICRUSt2's 0.018') | slope prFBA | +0.107 | +0.107 | +0.107 | changed below print precision |
| SR-12 | 81 ('+0.107, or +0.106 ... against +0.030 ... fourfold sweep'); 523 ('25, 50 and 100% subsamples ... 3.5 times as steeply (0.107, or 0.106 ... versus 0.030)'); 1189; 1217 ('quartering the shared reference costs prFBA 0.064 ... against PICRUSt2's 0.018') | slope ratio | 3.5 | 3.5 | 3.5 | changed below print precision |
| SR-12 | 81 ('+0.107, or +0.106 ... against +0.030 ... fourfold sweep'); 523 ('25, 50 and 100% subsamples ... 3.5 times as steeply (0.107, or 0.106 ... versus 0.030)'); 1189; 1217 ('quartering the shared reference costs prFBA 0.064 ... against PICRUSt2's 0.018') | quartering cost prFBA | 0.064 | 0.064 | 0.064 | changed below print precision |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | legacy single share | 61.7 | 61.7 | 61.7 | changed below print precision |
| disc.selection_size_controlled | 1170, 1229 | legacy single share | 61.7 | 61.7 | 61.7 | changed below print precision |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | as-shipped cond F1 prFBA | 0.802 | 0.802 | 0.802 | changed below print precision |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | as-shipped matched PICRUSt2 | 0.727 | 0.727 | 0.727 | changed below print precision |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | as-shipped matched margin | +0.075 | +0.075 | +0.075 | changed below print precision |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | as-shipped PICRUSt2 pop | 0.729 | 0.729 | 0.729 | changed below print precision |
| disc.vocab_breadth | 1189, 1217 | slope prFBA | 0.107 | 0.107 | 0.107 | changed below print precision |
| disc.vocab_breadth | 1189, 1217 | quartering cost prFBA | 0.064 | 0.064 | 0.064 | changed below print precision |
| abs-controlled-genus | 76 | picrust2 pop F1 | 0.779 | 0.779 | 0.779 |  |
| abs-controlled-genus | 76 | prfba_nestedcv pop F1 | 0.687 | 0.687 | 0.687 |  |
| abs-controlled-genus | 76 | prfba_nestedcv_forced pop F1 | 0.747 | 0.747 | 0.747 |  |
| abs-controlled-genus | 76 | delta nestedcv | -0.092 | -0.092 | -0.092 |  |
| abs-controlled-genus | 76 | delta exact-tie | -0.094 | -0.094 | -0.094 |  |
| abs-controlled-genus | 76 | delta forced | -0.031 | -0.031 | -0.031 |  |
| abs-controlled-genus | 76 | delta exact-tie forced | -0.033 | -0.033 | -0.033 |  |
| abs-dialect-density | 81 | density slope exact-tie | +0.106 | +0.106 | +0.106 |  |
| abs-dialect-density | 81 | density slope PICRUSt2 | +0.030 | +0.030 | +0.030 |  |
| abs-region-penalty | 82 | controlled FL pop F1 | 0.692 | 0.692 | 0.692 |  |
| abs-region-penalty | 82 | controlled V4 pop F1 | 0.710 | 0.710 | 0.710 |  |
| abs-region-penalty | 82 | exact-tie FL pop F1 | 0.690 | 0.690 | 0.690 |  |
| abs-region-penalty | 82 | exact-tie V4 pop F1 | 0.708 | 0.708 | 0.708 |  |
| abs-region-penalty | 82 | controlled FL coverage | 0.890 | 0.890 | 0.890 |  |
| abs-region-penalty | 82 | controlled V4 coverage | 0.933 | 0.933 | 0.933 |  |
| tab:sharedref | 681-707 | picrust2 coverage | 1.000 | 1.000 | 1.000 |  |
| tab:sharedref | 681-707 | picrust2 pop F1 | 0.779 | 0.779 | 0.779 |  |
| tab:sharedref | 681-707 | picrust2 cond (matched) F1 | 0.796 | 0.796 | 0.796 |  |
| tab:sharedref | 681-707 | picrust2_nsti_routed coverage | 1.000 | 1.000 | 1.000 |  |
| tab:sharedref | 681-707 | picrust2_nsti_routed pop F1 | 0.775 | 0.775 | 0.775 |  |
| tab:sharedref | 681-707 | prfba_nestedcv coverage | 0.893 | 0.893 | 0.893 |  |
| tab:sharedref | 681-707 | prfba_nestedcv pop F1 | 0.687 | 0.687 | 0.687 |  |
| tab:sharedref | 681-707 | prfba_nestedcv cond (matched) F1 | 0.769 | 0.769 | 0.769 |  |
| tab:sharedref | 681-707 | prfba_nestedcv_forced coverage | 0.9995 | 0.9995 | 0.9995 |  |
| tab:sharedref | 681-707 | prfba_nestedcv_forced pop F1 | 0.747 | 0.747 | 0.747 |  |
| tab:sharedref | 681-707 | prfba_nestedcv_forced cond (matched) F1 | 0.769 | 0.769 | 0.769 |  |
| tab:sharedref | 681-707 | prfba_exact_tie coverage | 0.893 | 0.893 | 0.893 |  |
| tab:sharedref | 681-707 | prfba_exact_tie cond (matched) F1 | 0.767 | 0.767 | 0.767 |  |
| tab:sharedref | 681-707 | prfba_exact_tie skill | 0.27 | 0.27 | 0.27 |  |
| tab:sharedref | 681-707 | prfba_exact_tie_forced coverage | 0.9995 | 0.9995 | 0.9995 |  |
| tab:sharedref | 681-707 | prfba_exact_tie_forced cond (matched) F1 | 0.767 | 0.767 | 0.767 |  |
| tab:sharedref | 681-707 | prfba_exact_tie_forced skill | 0.69 | 0.69 | 0.69 |  |
| tab:sharedref | 681-707 | base_nn1 coverage | 1.000 | 1.000 | 1.000 |  |
| tab:sharedref | 681-707 | base_nn1 pop F1 | 0.745 | 0.745 | 0.745 |  |
| tab:sharedref | 681-707 | base_nn1 cond (matched) F1 | 0.766 | 0.766 | 0.766 |  |
| tab:sharedref | 681-707 | base_random_abstain coverage | 0.892 | 0.892 | 0.892 |  |
| tab:sharedref | 681-707 | base_random_abstain pop F1 | 0.666 | 0.666 | 0.666 |  |
| tab:sharedref | 681-707 | base_random_abstain cond (matched) F1 | 0.767 | 0.767 | 0.767 |  |
| tab:sharedref | 681-707 | base_family_donor coverage | 0.679 | 0.679 | 0.679 |  |
| tab:sharedref | 681-707 | base_family_donor pop F1 | 0.516 | 0.516 | 0.516 |  |
| tab:sharedref | 681-707 | base_family_donor cond (matched) F1 | 0.765 | 0.765 | 0.765 |  |
| tab:sharedref | 681-707 | base_constant coverage | 1.000 | 1.000 | 1.000 |  |
| tab:sharedref | 681-707 | base_constant skill | 0.00 | 0.00 | 0.00 |  |
| tab:sharedref | 681-707 | base_clade_max coverage | 1.000 | 1.000 | 1.000 |  |
| tab:sharedref | 681-707 | base_clade_max cond (matched) F1 | 0.805 | 0.805 | 0.805 |  |
| tab:sharedref | 681-707 | base_clade_max skill | 1.00 | 1.00 | 1.00 |  |
| tab:sharedref | 681-707 | matched cells | 57,021 | 57,021 | 57,021 |  |
| tab:sharedref | 681-707 | random-abstain matched cells | 50,946 | 50,946 | 50,946 |  |
| tab:sharedref | 681-707 | family-donor matched cells | 41,257 | 41,257 | 41,257 |  |
| tab:sharedref | 681-707 | scored cells | 63,871 | 63,871 | 63,871 |  |
| tab:sharedref | 681-707 | reference genomes | 7,175 | 7,175 | 7,175 |  |
| sharedref_selection_sizes | 219; 303 ("86.8% of answered amplicons on the controlled benchmark") | fold 0 single-genome answered | 9,846 | 9,846 | 9,846 |  |
| sharedref_selection_sizes | 219; 303 ("86.8% of answered amplicons on the controlled benchmark") | fold 0 answered | 11,751 | 11,751 | 11,751 |  |
| sharedref_selection_sizes | 219; 303 ("86.8% of answered amplicons on the controlled benchmark") | pooled single-genome share | 86.8 | 86.8 | 86.8 |  |
| sharedref_selection_sizes | 219; 303 ("86.8% of answered amplicons on the controlled benchmark") | pooled answered | 50,842 | 50,842 | 50,842 |  |
| sharedref_selection_sizes | 219; 303 ("86.8% of answered amplicons on the controlled benchmark") | mean selection size | 1.199 | 1.199 | 1.199 |  |
| sharedref_ngenomes_vs_identity | 228 | P(n>=2) rises monotonically | yes | yes | yes |  |
| sharedref_ngenomes_vs_identity | 228 | P(n>=2) lowest novelty bin (not printed) | None | 0.100 | 0.100 |  |
| sharedref_ngenomes_vs_identity | 228 | P(n>=2) highest novelty bin (not printed) | None | 0.180 | 0.180 |  |
| controlled_forced_recovery | 736-737 | forced coverage | 99.95 | 99.95 | 99.95 |  |
| controlled_forced_recovery | 736-737 | recovered by forcing | +0.061 | +0.061 | +0.061 |  |
| controlled_forced_recovery | 736-737 | deficit | 0.092 | 0.092 | 0.092 |  |
| controlled_forced_recovery | 736-737 | exact-tie deficit | 0.094 | 0.094 | 0.094 |  |
| controlled_forced_recovery | 736-737 | residual forced | -0.031 | -0.031 | -0.031 |  |
| controlled_forced_recovery | 736-737 | residual exact-tie forced | -0.033 | -0.033 | -0.033 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | fold 0 single | 9,846 | 9,846 | 9,846 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | fold 0 answered | 11,751 | 11,751 | 11,751 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | fold 0 two | 1,354 | 1,354 | 1,354 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | fold 0 three+ | 551 | 551 | 551 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | pooled single share | 86.8 | 86.8 | 86.8 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | pooled single | 44,132 | 44,132 | 44,132 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | pooled answered | 50,842 | 50,842 | 50,842 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | mean size | 1.199 | 1.199 | 1.199 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | cross-phylum tie sets | 6.4 | 6.4 | 6.4 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | same-genus tie sets | 35.8 | 35.8 | 35.8 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | multi-genome selections spanning phyla | 4.8 | 4.8 | 4.8 |  |
| controlled_selection_sizes | 219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets) | P(n>=2) rises monotonically | yes | yes | yes |  |
| SR-02 | 77 ('picked legacy in every fold'); 340 (10/10, 0.001-0.003); 588-590 (legacy in all five folds at both ranks by 0.001-0.003 F1; exact-tie 0.002 lower at identical coverage) | nested-CV winners legacy (count) | 5 | 5 | 5 |  |
| SR-02 | 77 ('picked legacy in every fold'); 340 (10/10, 0.001-0.003); 588-590 (legacy in all five folds at both ranks by 0.001-0.003 F1; exact-tie 0.002 lower at identical coverage) | nestedcv pop F1 | 0.687 | 0.687 | 0.687 |  |
| SR-02 | 77 ('picked legacy in every fold'); 340 (10/10, 0.001-0.003); 588-590 (legacy in all five folds at both ranks by 0.001-0.003 F1; exact-tie 0.002 lower at identical coverage) | family nested-CV winners legacy (count) | 5 | 5 | 5 |  |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | picrust2 pop F1 | 0.779 | 0.779 | 0.779 |  |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_nestedcv pop F1 | 0.687 | 0.687 | 0.687 |  |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_nestedcv_forced pop F1 | 0.747 | 0.747 | 0.747 |  |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prFBA coverage | 0.893 | 0.893 | 0.893 |  |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | forced coverage | 0.9995 | 0.9995 | 0.9995 |  |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_nestedcv delta | -0.092 | -0.092 | -0.092 |  |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_exact_tie delta | -0.094 | -0.094 | -0.094 |  |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_nestedcv_forced delta | -0.031 | -0.031 | -0.031 |  |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | prfba_exact_tie_forced delta | -0.033 | -0.033 | -0.033 |  |
| SR-03 | 76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094') | SESOI verdict | exceeded | exceeded | exceeded |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_nestedcv by region (smallest deficit) | -0.059 | -0.059 | -0.059 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_nestedcv by region (largest deficit) | -0.182 | -0.182 | -0.182 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_nestedcv_forced by region (smallest deficit) | -0.029 | -0.029 | -0.029 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_nestedcv_forced by region (largest deficit) | -0.034 | -0.034 | -0.034 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_exact_tie by region (smallest deficit) | -0.060 | -0.060 | -0.060 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_exact_tie by region (largest deficit) | -0.184 | -0.184 | -0.184 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_exact_tie_forced by region (smallest deficit) | -0.031 | -0.031 | -0.031 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_exact_tie_forced by region (largest deficit) | -0.036 | -0.036 | -0.036 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_nestedcv Bacteria | -0.095 | -0.095 | -0.095 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_nestedcv_forced Bacteria | -0.032 | -0.032 | -0.032 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_exact_tie Bacteria | -0.097 | -0.097 | -0.097 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_exact_tie Archaea | -0.053 | -0.053 | -0.053 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | prfba_exact_tie_forced Bacteria | -0.034 | -0.034 | -0.034 |  |
| SR-04 | 76 ('the deficit holds in all eight regions'); 598 | max BH q over regions (nestedcv) | 0.01 | 0.002022 | 0.002022 |  |
| SR-05 | 598 (AURC 0.812 vs 0.788, exact-tie 0.787); 599 (cross-knob 0.813 and 0.788, exact-tie 0.786; family swap ~0.007) | prFBA AURC | 0.788 | 0.788 | 0.788 |  |
| SR-05 | 598 (AURC 0.812 vs 0.788, exact-tie 0.787); 599 (cross-knob 0.813 and 0.788, exact-tie 0.786; family swap ~0.007) | exact-tie AURC | 0.787 | 0.787 | 0.787 |  |
| SR-05 | 598 (AURC 0.812 vs 0.788, exact-tie 0.787); 599 (cross-knob 0.813 and 0.788, exact-tie 0.786; family swap ~0.007) | PICRUSt2 cross-knob AURC | 0.813 | 0.813 | 0.813 |  |
| SR-05 | 598 (AURC 0.812 vs 0.788, exact-tie 0.787); 599 (cross-knob 0.813 and 0.788, exact-tie 0.786; family swap ~0.007) | exact-tie cross-knob AURC | 0.786 | 0.786 | 0.786 |  |
| SR-05 | 598 (AURC 0.812 vs 0.788, exact-tie 0.787); 599 (cross-knob 0.813 and 0.788, exact-tie 0.786; family swap ~0.007) | family swap (prFBA aurc - cross) | 0.007 | 0.007 | 0.007 |  |
| SR-06 | 82 ('0.692 vs. 0.710; exact-tie 0.690 vs. 0.708 ... coverage 0.890 against 0.933'); 600 | FL pop F1 | 0.692 | 0.692 | 0.692 |  |
| SR-06 | 82 ('0.692 vs. 0.710; exact-tie 0.690 vs. 0.708 ... coverage 0.890 against 0.933'); 600 | V4 pop F1 | 0.710 | 0.710 | 0.710 |  |
| SR-06 | 82 ('0.692 vs. 0.710; exact-tie 0.690 vs. 0.708 ... coverage 0.890 against 0.933'); 600 | exact-tie FL | 0.690 | 0.690 | 0.690 |  |
| SR-06 | 82 ('0.692 vs. 0.710; exact-tie 0.690 vs. 0.708 ... coverage 0.890 against 0.933'); 600 | exact-tie V4 | 0.708 | 0.708 | 0.708 |  |
| SR-06 | 82 ('0.692 vs. 0.710; exact-tie 0.690 vs. 0.708 ... coverage 0.890 against 0.933'); 600 | FL coverage | 0.890 | 0.890 | 0.890 |  |
| SR-06 | 82 ('0.692 vs. 0.710; exact-tie 0.690 vs. 0.708 ... coverage 0.890 against 0.933'); 600 | V4 coverage | 0.933 | 0.933 | 0.933 |  |
| SR-06 | 82 ('0.692 vs. 0.710; exact-tie 0.690 vs. 0.708 ... coverage 0.890 against 0.933'); 600 | exact-tie FL coverage | 0.890 | 0.890 | 0.890 |  |
| SR-06 | 82 ('0.692 vs. 0.710; exact-tie 0.690 vs. 0.708 ... coverage 0.890 against 0.933'); 600 | exact-tie V4 coverage | 0.933 | 0.933 | 0.933 |  |
| SR-07 | 601-602 | prFBA rxn pop F1 | 0.691 | 0.691 | 0.691 |  |
| SR-07 | 601-602 | delta | -0.089 | -0.089 | -0.089 |  |
| SR-07 | 601-602 |   CI lo | -0.105 | -0.105 | -0.105 |  |
| SR-07 | 601-602 |   CI hi | -0.078 | -0.078 | -0.078 |  |
| SR-07 | 601-602 |   CI lo | -0.107 | -0.107 | -0.107 |  |
| SR-07 | 601-602 |   CI hi | -0.081 | -0.081 | -0.081 |  |
| SR-07 | 601-602 | forced delta | -0.028 | -0.028 | -0.028 |  |
| SR-07 | 601-602 | exact-tie forced delta | -0.031 | -0.031 | -0.031 |  |
| SR-07 | 601-602 | max BH q over regions | 0.002 | 0.002 | 0.002 |  |
| SR-07 | 601-602 | PICRUSt2 rxn AURC | 0.815 | 0.815 | 0.815 |  |
| SR-07 | 601-602 | prFBA rxn AURC | 0.794 | 0.794 | 0.794 |  |
| SR-07 | 601-602 | exact-tie rxn AURC | 0.792 | 0.792 | 0.792 |  |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | exact-tie forced skill | 0.69 | 0.69 | 0.69 |  |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | forced pop F1 | 0.747 | 0.747 | 0.747 |  |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | nn1 pop F1 | 0.745 | 0.745 | 0.745 |  |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | nestedcv vs nn1 | -0.058 | -0.058 | -0.058 |  |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 |   q | 0.002 | 0.002 | 0.002 |  |
| SR-08 | 77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184 | PICRUSt2 vs nn1 | +0.034 | +0.034 | +0.034 |  |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') | vs random abstain | +0.021 | +0.021 | +0.021 |  |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') |   q | 0.002 | 0.002 | 0.002 |  |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') | declined share | 10.7 | 10.7 | 10.7 |  |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') | random-abstain matched cells | 50,946 | 50,946 | 50,946 |  |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') | answering everywhere gains | +0.061 | +0.061 | +0.061 |  |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') | forced pop F1 refused band | 0.567 | 0.567 | 0.567 |  |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') | exact-tie forced refused band | 0.567 | 0.567 | 0.567 |  |
| SR-09 | 78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit') | PICRUSt2 refused band | 0.634 | 0.634 | 0.634 |  |
| SR-10 | 620-621 | prFBA precision 0.865-0.945 | 0.740 | 0.740 | 0.740 |  |
| SR-10 | 620-621 | exact-tie precision 0.865-0.945 | 0.751 | 0.751 | 0.751 |  |
| SR-10 | 620-621 | prFBA coverage <0.865 | 0 | 0 | 0 |  |
| SR-10 | 620-621 | PICRUSt2 pop F1 <0.865 | 0.634 | 0.634 | 0.634 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | picrust2 coverage | 1.000 | 1.000 | 1.000 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | picrust2 pop F1 | 0.779 | 0.779 | 0.779 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | picrust2 cond (matched) F1 | 0.796 | 0.796 | 0.796 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | picrust2_nsti_routed coverage | 1.000 | 1.000 | 1.000 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | picrust2_nsti_routed pop F1 | 0.775 | 0.775 | 0.775 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_nestedcv coverage | 0.893 | 0.893 | 0.893 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_nestedcv pop F1 | 0.687 | 0.687 | 0.687 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_nestedcv cond (matched) F1 | 0.769 | 0.769 | 0.769 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_nestedcv_forced coverage | 0.9995 | 0.9995 | 0.9995 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_nestedcv_forced pop F1 | 0.747 | 0.747 | 0.747 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_nestedcv_forced cond (matched) F1 | 0.769 | 0.769 | 0.769 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_exact_tie coverage | 0.893 | 0.893 | 0.893 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_exact_tie cond (matched) F1 | 0.767 | 0.767 | 0.767 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_exact_tie skill | 0.27 | 0.27 | 0.27 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_exact_tie_forced coverage | 0.9995 | 0.9995 | 0.9995 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_exact_tie_forced cond (matched) F1 | 0.767 | 0.767 | 0.767 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | prfba_exact_tie_forced skill | 0.69 | 0.69 | 0.69 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_nn1 coverage | 1.000 | 1.000 | 1.000 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_nn1 pop F1 | 0.745 | 0.745 | 0.745 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_nn1 cond (matched) F1 | 0.766 | 0.766 | 0.766 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_random_abstain coverage | 0.892 | 0.892 | 0.892 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_random_abstain pop F1 | 0.666 | 0.666 | 0.666 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_random_abstain cond (matched) F1 | 0.767 | 0.767 | 0.767 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_family_donor coverage | 0.679 | 0.679 | 0.679 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_family_donor pop F1 | 0.516 | 0.516 | 0.516 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_family_donor cond (matched) F1 | 0.765 | 0.765 | 0.765 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_constant coverage | 1.000 | 1.000 | 1.000 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_constant skill | 0.00 | 0.00 | 0.00 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_clade_max coverage | 1.000 | 1.000 | 1.000 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_clade_max cond (matched) F1 | 0.805 | 0.805 | 0.805 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | base_clade_max skill | 1.00 | 1.00 | 1.00 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | matched cells | 57,021 | 57,021 | 57,021 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | random-abstain matched cells | 50,946 | 50,946 | 50,946 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | family-donor matched cells | 41,257 | 41,257 | 41,257 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | scored cells | 63,871 | 63,871 | 63,871 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | reference genomes | 7,175 | 7,175 | 7,175 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | matched conditional gap | -0.027 | -0.027 | -0.027 |  |
| SR-11 | 683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796) | exact-tie matched gap | -0.029 | -0.029 | -0.029 |  |
| SR-12 | 81 ('+0.107, or +0.106 ... against +0.030 ... fourfold sweep'); 523 ('25, 50 and 100% subsamples ... 3.5 times as steeply (0.107, or 0.106 ... versus 0.030)'); 1189; 1217 ('quartering the shared reference costs prFBA 0.064 ... against PICRUSt2's 0.018') | slope exact-tie | +0.106 | +0.106 | +0.106 |  |
| SR-12 | 81 ('+0.107, or +0.106 ... against +0.030 ... fourfold sweep'); 523 ('25, 50 and 100% subsamples ... 3.5 times as steeply (0.107, or 0.106 ... versus 0.030)'); 1189; 1217 ('quartering the shared reference costs prFBA 0.064 ... against PICRUSt2's 0.018') | slope PICRUSt2 | +0.030 | +0.030 | +0.030 |  |
| SR-12 | 81 ('+0.107, or +0.106 ... against +0.030 ... fourfold sweep'); 523 ('25, 50 and 100% subsamples ... 3.5 times as steeply (0.107, or 0.106 ... versus 0.030)'); 1189; 1217 ('quartering the shared reference costs prFBA 0.064 ... against PICRUSt2's 0.018') | quartering cost exact-tie | 0.064 | 0.064 | 0.064 |  |
| SR-12 | 81 ('+0.107, or +0.106 ... against +0.030 ... fourfold sweep'); 523 ('25, 50 and 100% subsamples ... 3.5 times as steeply (0.107, or 0.106 ... versus 0.030)'); 1189; 1217 ('quartering the shared reference costs prFBA 0.064 ... against PICRUSt2's 0.018') | quartering cost PICRUSt2 | 0.018 | 0.018 | 0.018 |  |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | pooled single share | 86.8 | 86.8 | 86.8 |  |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | pooled answered | 50,842 | 50,842 | 50,842 |  |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | pooled single | 44,132 | 44,132 | 44,132 |  |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | fold 0 single | 9,846 | 9,846 | 9,846 |  |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | fold 0 answered | 11,751 | 11,751 | 11,751 |  |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | fold 0 two | 1,354 | 1,354 | 1,354 |  |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | fold 0 three+ | 551 | 551 | 551 |  |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | mean size | 1.199 | 1.199 | 1.199 |  |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | cross-phylum tie sets | 6.4 | 6.4 | 6.4 |  |
| XR-03 | 219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229 | same-genus tie sets | 35.8 | 35.8 | 35.8 |  |
| disc.selection_size_controlled | 1170, 1229 | pooled single share | 86.8 | 86.8 | 86.8 |  |
| disc.selection_size_controlled | 1170, 1229 | pooled single | 44,132 | 44,132 | 44,132 |  |
| disc.selection_size_controlled | 1170, 1229 | pooled answered | 50,842 | 50,842 | 50,842 |  |
| disc.selection_size_controlled | 1170, 1229 | mean size | 1.199 | 1.199 | 1.199 |  |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | controlled matched prFBA | 0.769 | 0.769 | 0.769 |  |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | controlled matched PICRUSt2 | 0.796 | 0.796 | 0.796 |  |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | controlled matched exact-tie | 0.767 | 0.767 | 0.767 |  |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | matched reversal | -0.027 | -0.027 | -0.027 |  |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | exact-tie matched reversal | -0.029 | -0.029 | -0.029 |  |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | controlled deficit | -0.092 | -0.092 | -0.092 |  |
| disc.margin_reversal | 1184, 1187-1188, 1215-1216 | controlled exact-tie deficit | -0.094 | -0.094 | -0.094 |  |
| disc.vocab_breadth | 1189, 1217 | ceiling recall | 90.4 | 90.4 | 90.4 |  |
| disc.vocab_breadth | 1189, 1217 | slope exact-tie | 0.106 | 0.106 | 0.106 |  |
| disc.vocab_breadth | 1189, 1217 | slope PICRUSt2 | 0.030 | 0.030 | 0.030 |  |
| disc.vocab_breadth | 1189, 1217 | quartering cost PICRUSt2 | 0.018 | 0.018 | 0.018 |  |
| SR-02 | 77 ('picked legacy in every fold'); 340 (10/10, 0.001-0.003); 588-590 (legacy in all five folds at both ranks by 0.001-0.003 F1; exact-tie 0.002 lower at identical coverage) | (manual) | | | | **not recomputed** (inner margins 0.001-0.003): inner-CV F1 gaps (upper end by hand); the winners are unchanged (legacy in 5 of 5 folds at genus, R3.5) |

## R3F

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| abs-controlled-family-nn1 | 77 | family delta | -0.123 | -0.124 | -0.123 | **CHANGED** |
| abs-controlled-random-mask | 78 |   CI lo | -0.004 | -0.004 | -0.005 | **CHANGED** |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 | PICRUSt2 pop F1 | 0.749 | 0.749 | 0.750 | **CHANGED** |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 | delta | -0.123 | -0.124 | -0.123 | **CHANGED** |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | exact-tie by region (largest) | -0.247 | -0.247 | -0.248 | **CHANGED** |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | forced pop F1 | 0.7095 | 0.7095 | 0.7096 | **CHANGED** |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | nn1 pop F1 | 0.7070 | 0.7070 | 0.7072 | **CHANGED** |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | exact-tie forced pop F1 | 0.7072 | 0.7072 | 0.7074 | **CHANGED** |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 |   CI lo | -0.053 | -0.053 | -0.052 | **CHANGED** |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | exact-tie vs nn1 | -0.083 | -0.083 | -0.084 | **CHANGED** |
| abs-controlled-family-nn1 | 77 |   CI lo | -0.150 | -0.150 | -0.150 | changed below print precision |
| abs-controlled-family-nn1 | 77 |   CI hi | -0.099 | -0.099 | -0.099 | changed below print precision |
| abs-controlled-family-nn1 | 77 |   CI lo | -0.152 | -0.152 | -0.152 | changed below print precision |
| abs-controlled-family-nn1 | 77 |   CI hi | -0.102 | -0.102 | -0.102 | changed below print precision |
| abs-controlled-family-nn1 | 77 |   CI lo | -0.069 | -0.069 | -0.069 | changed below print precision |
| abs-controlled-family-nn1 | 77 |   CI hi | -0.051 | -0.051 | -0.051 | changed below print precision |
| abs-controlled-family-nn1 | 77 | exact-tie vs nn1 | -0.060 | -0.060 | -0.060 | changed below print precision |
| abs-controlled-family-nn1 | 77 | exact-tie forced pop F1 | 0.746 | 0.746 | 0.746 | changed below print precision |
| abs-controlled-random-mask | 78 |   CI hi | 0.040 | 0.040 | 0.040 | changed below print precision |
| abs-controlled-random-mask | 78 | exact-tie vs own mask (genus) | +0.019 | +0.019 | +0.019 | changed below print precision |
| abs-controlled-random-mask | 78 | exact-tie vs own mask (family) | +0.015 | +0.015 | +0.015 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 | prFBA pop F1 | 0.626 | 0.626 | 0.626 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 |   CI lo | -0.150 | -0.150 | -0.150 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 |   CI hi | -0.099 | -0.099 | -0.099 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 | exact-tie pop F1 | 0.623 | 0.624 | 0.624 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 |   CI lo | -0.152 | -0.152 | -0.152 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 |   CI hi | -0.102 | -0.102 | -0.102 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 | forced pop F1 | 0.709 | 0.710 | 0.710 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 | forced delta | -0.040 | -0.040 | -0.040 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 |   CI lo | -0.048 | -0.048 | -0.048 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 |   CI hi | -0.031 | -0.031 | -0.031 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 | exact-tie forced pop F1 | 0.707 | 0.707 | 0.707 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 | exact-tie forced delta | -0.042 | -0.042 | -0.042 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 |   CI lo | -0.050 | -0.050 | -0.050 | changed below print precision |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 |   CI hi | -0.033 | -0.033 | -0.033 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | by region (smallest deficit) | -0.076 | -0.076 | -0.076 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | by region (largest deficit) | -0.245 | -0.245 | -0.245 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | vs nn1 | -0.081 | -0.081 | -0.081 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 |   CI lo | -0.107 | -0.107 | -0.107 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 |   CI hi | -0.059 | -0.059 | -0.059 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 |   CI hi | 0.018 | 0.018 | 0.018 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | skill | -0.14 | -0.14 | -0.14 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | vs family donor | +0.107 | +0.107 | +0.107 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 |   CI lo | 0.044 | 0.044 | 0.044 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 |   CI hi | 0.179 | 0.179 | 0.179 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | exact-tie vs own mask | +0.015 | +0.015 | +0.015 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | legacy like-for-like | +0.016 | +0.016 | +0.016 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | exact-tie vs constant | -0.020 | -0.020 | -0.020 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | exact-tie skill | -0.16 | -0.16 | -0.16 | changed below print precision |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | exact-tie vs family donor | +0.104 | +0.104 | +0.104 | changed below print precision |
| abs-controlled-family-nn1 | 77 | family exact-tie delta | -0.126 | -0.126 | -0.126 |  |
| abs-controlled-family-nn1 | 77 | genus nestedcv vs nn1 | -0.058 | -0.058 | -0.058 |  |
| abs-controlled-family-nn1 | 77 | forced pop F1 | 0.747 | 0.747 | 0.747 |  |
| abs-controlled-family-nn1 | 77 | nn1 pop F1 | 0.745 | 0.745 | 0.745 |  |
| abs-controlled-random-mask | 78 | genus vs random mask | +0.021 | +0.021 | +0.021 |  |
| abs-controlled-random-mask | 78 | family vs random mask | +0.018 | +0.018 | +0.018 |  |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 | prFBA coverage | 0.856 | 0.856 | 0.856 |  |
| FR-02 | 77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624 | exact-tie delta | -0.126 | -0.126 | -0.126 |  |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | exact-tie by region (smallest) | -0.077 | -0.077 | -0.077 |  |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | max q over regions | 0.002 | 0.002 | 0.002 |  |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | vs constant | -0.017 | -0.017 | -0.017 |  |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | vs random mask | +0.018 | +0.018 | +0.018 |  |
| FR-03 | 78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192 | family swap | 0.007 | 0.007 | 0.007 |  |

## R4P

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| abs-phenotype-readout | 79 | ATP-informative genomes | 755 | 755 | 759 | **CHANGED** |
| abs-phenotype-readout | 79 | hamming delta | 0.13 | 0.13 | 0.07 | **CHANGED** |
| abs-phenotype-readout | 79 | Wilcoxon p | 0.49 | 0.49 | 0.79 | **CHANGED** |
| abs-phenotype-readout | 79 | of which V4 | 747 | 747 | 750 | **CHANGED** |
| abs-phenotype-readout | 79 | EC-F1 delta prFBA - PICRUSt2 on the paired genomes | -0.028 | -0.028 | -0.031 | **CHANGED** |
| abs-phenotype-readout | 79 | media favouring PICRUSt2 (McNemar q<=0.05) | 2 | 2 | 4 | **CHANGED** |
| PH-01 | 79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755) | paired genomes | 755 | 755 | 759 | **CHANGED** |
| PH-01 | 79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755) | hamming delta | 0.13 | 0.13 | 0.07 | **CHANGED** |
| PH-01 | 79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755) | Wilcoxon p | 0.49 | 0.49 | 0.79 | **CHANGED** |
| PH-01 | 79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755) | q | 0.49 | 0.49 | 0.79 | **CHANGED** |
| PH-01 | 79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755) | of which V4 | 747 | 747 | 750 | **CHANGED** |
| PH-01 | 79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755) | EC-F1 delta prFBA - PICRUSt2 | -0.028 | -0.028 | -0.031 | **CHANGED** |
| PH-01 | 79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755) | all-genome hamming delta (n 1,000) | None | 0.22 | 0.28 | **CHANGED** |
| PH-01 | 79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755) | all-genome Wilcoxon p | None | 0.110 | 0.057 | **CHANGED** |
| PH-02 | 79 ('two of 22 media'); 631 (Ac.NO and Etho.O2, q=0.002); 639; 1185 | Ac.NO q | 0.002 | 0.002 | 0.029 | **CHANGED**; MEANING: significance crosses 0.01 (0.00243 -> 0.0289) |
| PH-02 | 79 ('two of 22 media'); 631 (Ac.NO and Etho.O2, q=0.002); 639; 1185 | Etho.O2 q | 0.002 | 0.002 | 0.046 | **CHANGED**; MEANING: significance crosses 0.01 (0.00243 -> 0.0461) |
| PH-02 | 79 ('two of 22 media'); 631 (Ac.NO and Etho.O2, q=0.002); 639; 1185 | media with McNemar q<=0.05 | 2 | 2 | 4 | **CHANGED** |
| PH-02 | 79 ('two of 22 media'); 631 (Ac.NO and Etho.O2, q=0.002); 639; 1185 | of which favour PICRUSt2 | 2 | 2 | 4 | **CHANGED** |
| PH-02 | 79 ('two of 22 media'); 631 (Ac.NO and Etho.O2, q=0.002); 639; 1185 | media favouring PICRUSt2 | Ac.NO,Etho.O2 | Ac.NO,Etho.O2 | Ac.NO,Ac.O2,Etho.O2,For.O2 | **CHANGED**; MEANING: categorical change Ac.NO,Etho.O2 -> Ac.NO,Ac.O2,Etho.O2,For.O2 |
| PH-02 | 79 ('two of 22 media'); 631 (Ac.NO and Etho.O2, q=0.002); 639; 1185 | largest q of those media | 0.002 | 0.002 | 0.046 | **CHANGED**; MEANING: significance crosses 0.01 (0.00243 -> 0.0461) |
| PH-02 | 79 ('two of 22 media'); 631 (Ac.NO and Etho.O2, q=0.002); 639; 1185 | smallest q of the media with q <= 0.05 | 0.002 | 0.002 | 0.029 | **CHANGED**; MEANING: significance crosses 0.01 (0.00243 -> 0.0289) |
| PH-03 | 634 (correlation -0.14); 1190 (Spearman rho = -0.14, p = 2x10^-4) | Spearman rho | -0.14 | -0.14 | -0.19 | **CHANGED** |
| PH-03 | 634 (correlation -0.14); 1190 (Spearman rho = -0.14, p = 2x10^-4) | Spearman p | 2e-4 | 1.7e-4 | 1.5e-7 | **CHANGED**; old≠printed (differs; known: the earlier text printed one significant digit (2x10^-4 for 1.65e-4); the re-scored text prints two (1.5x10^-7)) |
| PH-04 | 635-637; 1186 | hamming delta vs nn1 | 0.146 | 0.146 | 0.030 | **CHANGED** |
| PH-04 | 635-637; 1186 | p | 0.017 | 0.02 | 0.41 | **CHANGED**; MEANING: significance crosses 0.05 (0.0168 -> 0.407) |
| PH-04 | 635-637; 1186 | q | 0.028 | 0.03 | 0.51 | **CHANGED**; MEANING: significance crosses 0.05 (0.0285 -> 0.509) |
| PH-04 | 635-637; 1186 | p all genomes | 0.28 | 0.28 | 0.46 | **CHANGED** |
| PH-04 | 635-637; 1186 | EC F1 nn1 on the nn1-paired genomes | 0.781 | 0.782 | 0.780 | **CHANGED**; old≠printed (rounding; known: printed 0.781 has no surviving recipe; the cells-table recipe gives 0.782 for both arms on the pre-repair pass (see derived_numbers KNOWN_PRE_MISMATCH)) |
| PH-04 | 635-637; 1186 | EC F1 paired test p | 0.64 | 0.45 | 0.84 | **CHANGED**; old≠printed (differs; known: no paired test reproduces 0.64; the recipe's Wilcoxon gives 0.45 on the pre-repair pass: replace the printed p with the recipe value or drop it) |
| disc.atp_vs_picrust2 | 1185 | hamming delta | 0.13 | 0.13 | 0.07 | **CHANGED** |
| disc.atp_vs_picrust2 | 1185 | Wilcoxon p | 0.49 | 0.49 | 0.79 | **CHANGED** |
| disc.atp_vs_picrust2 | 1185 | q | 0.49 | 0.49 | 0.79 | **CHANGED** |
| disc.atp_vs_picrust2 | 1185 | n | 755 | 755 | 759 | **CHANGED** |
| disc.atp_vs_picrust2 | 1185 | Ac.NO q | 0.002 | 0.002 | 0.029 | **CHANGED**; MEANING: significance crosses 0.01 (0.00243 -> 0.0289) |
| disc.atp_vs_picrust2 | 1185 | Etho.O2 q | 0.002 | 0.002 | 0.046 | **CHANGED**; MEANING: significance crosses 0.01 (0.00243 -> 0.0461) |
| disc.atp_vs_picrust2 | 1185 | smallest q of the media with q <= 0.05 | 0.002 | 0.002 | 0.029 | **CHANGED**; MEANING: significance crosses 0.01 (0.00243 -> 0.0289) |
| disc.atp_vs_nn1 | 1186 | hamming delta vs nn1 | 0.146 | 0.146 | 0.030 | **CHANGED** |
| disc.atp_vs_nn1 | 1186 | p | 0.017 | 0.02 | 0.41 | **CHANGED**; MEANING: significance crosses 0.05 (0.0168 -> 0.407) |
| disc.atp_vs_nn1 | 1186 | q | 0.028 | 0.03 | 0.51 | **CHANGED**; MEANING: significance crosses 0.05 (0.0285 -> 0.509) |
| disc.atp_vs_nn1 | 1186 | p all | 0.28 | 0.28 | 0.46 | **CHANGED** |
| disc.atp_vs_nn1 | 1186 | n | 755 | 755 | 759 | **CHANGED** |
| disc.atp_vs_nn1 | 1186 | EC F1 nn1 on the nn1-paired genomes | 0.781 | 0.782 | 0.780 | **CHANGED**; old≠printed (rounding; known: printed 0.781 has no surviving recipe; the cells-table recipe gives 0.782 for both arms on the pre-repair pass (see derived_numbers KNOWN_PRE_MISMATCH)) |
| disc.atp_vs_nn1 | 1186 | EC F1 paired test p | 0.64 | 0.45 | 0.84 | **CHANGED**; old≠printed (differs; known: no paired test reproduces 0.64; the recipe's Wilcoxon gives 0.45 on the pre-repair pass: replace the printed p with the recipe value or drop it) |
| disc.refusal_credit | 1190-1193 | Spearman rho | -0.14 | -0.14 | -0.19 | **CHANGED** |
| disc.refusal_credit | 1190-1193 | Spearman p | 2e-4 | 1.7e-4 | 1.5e-7 | **CHANGED**; old≠printed (differs; known: the earlier text printed one significant digit (2x10^-4 for 1.65e-4); the re-scored text prints two (1.5x10^-7)) |
| disc.refusal_credit | 1190-1193 | EC-F1 delta prFBA - PICRUSt2 | -0.028 | -0.028 | -0.031 | **CHANGED** |
| lim.phenotype_scope | 1222 | paired genomes | 755 | 755 | 759 | **CHANGED** |
| PH-04 | 635-637; 1186 | EC F1 prFBA on the nn1-paired genomes | 0.782 | 0.782 | 0.782 | changed below print precision |
| disc.atp_vs_nn1 | 1186 | EC F1 prFBA on the nn1-paired genomes | 0.782 | 0.782 | 0.782 | changed below print precision |
| disc.refusal_credit | 1190-1193 | exact-tie vs its own masked run | +0.019 | +0.019 | +0.019 | changed below print precision |
| abs-phenotype-readout | 79 | media with ATP signal | 22 | 22 | 22 |  |
| abs-phenotype-readout | 79 | media in the ATP vector | 54 | 54 | 54 |  |
| PH-01 | 79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755) | subsample | 1,000 | 1,000 | 1,000 |  |
| PH-02 | 79 ('two of 22 media'); 631 (Ac.NO and Etho.O2, q=0.002); 639; 1185 | testable media (with ATP signal) | 22 | 22 | 22 |  |
| PH-04 | 635-637; 1186 | matched cond F1 exact-tie | 0.767 | 0.767 | 0.767 |  |
| PH-04 | 635-637; 1186 | matched cond F1 nn1 | 0.766 | 0.766 | 0.766 |  |
| PH-05 | 79 ('only 22 of 54 media carry ATP signal'); 638 (53 media + empty control; 22 carry signal); 1222 | media with ATP signal | 22 | 22 | 22 |  |
| PH-05 | 79 ('only 22 of 54 media carry ATP signal'); 638 (53 media + empty control; 22 carry signal); 1222 | media in the ATP vector (53 + empty control) | 54 | 54 | 54 |  |
| disc.refusal_credit | 1190-1193 | vs random mask | +0.021 | +0.021 | +0.021 |  |
| disc.refusal_credit | 1190-1193 | share of the deficit recovered by answering everywhere (printed 'two thirds') | None | 0.66 | 0.66 |  |
| lim.phenotype_scope | 1222 | media with ATP signal | 22 | 22 | 22 |  |
| lim.phenotype_scope | 1222 | media in the ATP vector | 54 | 54 | 54 |  |
| lim.phenotype_scope | 1222 | gapfilled subsample | 200 | 200 | 200 |  |
| lim.phenotype_scope | 1222 | (manual) | | | | **historical** (5 gated truth models): the pre-registered C12 failure (5 of 200 gated) is a fixed history (D9), not re-run |

## R4S2

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| abs-phenotype-controls | 80 | ATP-pass rows (A-C13 configuration, re-run) | 43,140 | 43,140 | 43,174 | **CHANGED** |
| abs-phenotype-controls | 80 | gapfilled-pass rows | 8,614 | 8,614 | 8,660 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | draft rxn Jaccard PICRUSt2 | 0.737 | 0.737 | 0.735 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | draft rxn Jaccard prFBA | 0.707 | 0.707 | 0.704 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | draft rxn Jaccard nn1 | 0.705 | 0.705 | 0.700 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | essentiality disagreement PICRUSt2 | 0.065 | 0.065 | 0.059 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | essentiality disagreement prFBA | 0.068 | 0.068 | 0.063 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | essentiality disagreement nn1 | 0.069 | 0.069 | 0.063 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | growth log-ratio prFBA | 0.008 | 0.008 | 0.020 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | growth log-ratio nn1 | 0.008 | 0.008 | 0.011 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | growth log-ratio PICRUSt2 | -0.006 | -0.006 | -0.003 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | gapfill burden truth (median) | 82.5 | 82.5 | 79.5 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | gapfill burden method arms (min median) | 80 | 80 | 78 | **CHANGED** |
| PH-06 | 641-646; 1222 (n=200) | gapfill burden method arms (max median) | 84 | 84 | 82 | **CHANGED** |
| AM-02 | 80 ('43,140 of 43,140 rows ... 0 of 8,614 rows, against 49 ... valid under the amendments'); 667-673; 676-679 (seven-reaction draft correction: 865 -> 95 of 1,181 drafts; medians move <= 0.004; burden <= 1 reaction, 1.5 as shipped; 138 of 200 truth drafts make ATP from nothing; 43,140 ATP-pass rows with no empty-medium ATP); 1193 (ATP pass reproduces every vector) | ATP-pass rows (A-C13 configuration, re-run) | 43,140 | 43,140 | 43,174 | **CHANGED** |
| AM-02 | 80 ('43,140 of 43,140 rows ... 0 of 8,614 rows, against 49 ... valid under the amendments'); 667-673; 676-679 (seven-reaction draft correction: 865 -> 95 of 1,181 drafts; medians move <= 0.004; burden <= 1 reaction, 1.5 as shipped; 138 of 200 truth drafts make ATP from nothing; 43,140 ATP-pass rows with no empty-medium ATP); 1193 (ATP pass reproduces every vector) | gapfilled-pass rows | 8,614 | 8,614 | 8,660 | **CHANGED** |
| AM-02 | 80 ('43,140 of 43,140 rows ... 0 of 8,614 rows, against 49 ... valid under the amendments'); 667-673; 676-679 (seven-reaction draft correction: 865 -> 95 of 1,181 drafts; medians move <= 0.004; burden <= 1 reaction, 1.5 as shipped; 138 of 200 truth drafts make ATP from nothing; 43,140 ATP-pass rows with no empty-medium ATP); 1193 (ATP pass reproduces every vector) | models applying a reversal (A-C13) | 114 | 114 | 107 | **CHANGED** |
| abs-phenotype-controls | 80 | ATP-pass rows with empty-medium ATP | 0 | 0 | 0 |  |
| abs-phenotype-controls | 80 | gapfilled rows gated (C12 leg c) | 0 | 0 | 0 |  |
| PH-06 | 641-646; 1222 (n=200) | pre-registered essentiality PICRUSt2 | 0.066 | 0.066 | 0.066 |  |
| PH-06 | 641-646; 1222 (n=200) | pre-registered essentiality prFBA | 0.068 | 0.068 | 0.068 |  |
| PH-06 | 641-646; 1222 (n=200) | pre-registered growth prFBA | 0.002 | 0.002 | 0.002 |  |
| PH-06 | 641-646; 1222 (n=200) | pre-registered growth nn1 | 0.004 | 0.004 | 0.004 |  |
| AM-02 | 80 ('43,140 of 43,140 rows ... 0 of 8,614 rows, against 49 ... valid under the amendments'); 667-673; 676-679 (seven-reaction draft correction: 865 -> 95 of 1,181 drafts; medians move <= 0.004; burden <= 1 reaction, 1.5 as shipped; 138 of 200 truth drafts make ATP from nothing; 43,140 ATP-pass rows with no empty-medium ATP); 1193 (ATP pass reproduces every vector) | ATP-pass rows with empty-medium ATP | 0 | 0 | 0 |  |
| AM-02 | 80 ('43,140 of 43,140 rows ... 0 of 8,614 rows, against 49 ... valid under the amendments'); 667-673; 676-679 (seven-reaction draft correction: 865 -> 95 of 1,181 drafts; medians move <= 0.004; burden <= 1 reaction, 1.5 as shipped; 138 of 200 truth drafts make ATP from nothing; 43,140 ATP-pass rows with no empty-medium ATP); 1193 (ATP pass reproduces every vector) | gapfilled rows gated | 0 | 0 | 0 |  |
| AM-03 | 648-652 | C2 source recovery | 1.0 | 1.0 | 1.0 |  |
| AM-03 | 648-652 | C2 twin-free cells | 52,979 | 52,979 | 52,979 |  |
| AM-03 | 648-652 | C2 PICRUSt2 mean EC F1 | 0.899 | 0.899 | 0.899 |  |
| AM-03 | 648-652 | C3 status | pass | pass | pass |  |
| AM-03 | 648-652 | C4 status | pass | pass | pass |  |
| AM-03 | 648-652 | C6 status | pass | pass | pass |  |
| AM-03 | 648-652 | C7 status | pass | pass | pass |  |
| AM-03 | 648-652 | C8 status | pass | pass | pass |  |
| AM-03 | 648-652 | C10 status | pass | pass | pass |  |
| AM-03 | 648-652 | C13 status | pass | pass | pass |  |
| AM-03 | 648-652 | C14 status | pass | pass | pass |  |
| AM-03 | 648-652 | C15 status | pass | pass | pass |  |
| abs-phenotype-controls | 80 | (manual) | | | | **historical** (5 of 200 gated (pre-registered); 49 gated before A-C12): D9: the pre-registered and A-C12 passes are historical records; only the A-C13 configuration was re-run (R4.1), on a re-drawn sample |
| PH-06 | 641-646; 1222 (n=200) | (manual) | | | | **historical** (pre-registered pair 0.066 vs 0.068, 0.002 vs 0.004): pre-registered gapfilled pass (phenotype_gf/, D9); unchanged file |
| AM-02 | 80 ('43,140 of 43,140 rows ... 0 of 8,614 rows, against 49 ... valid under the amendments'); 667-673; 676-679 (seven-reaction draft correction: 865 -> 95 of 1,181 drafts; medians move <= 0.004; burden <= 1 reaction, 1.5 as shipped; 138 of 200 truth drafts make ATP from nothing; 43,140 ATP-pass rows with no empty-medium ATP); 1193 (ATP pass reproduces every vector) | (manual) | | | | **historical** (49 gated before; 136 models in 50 genomes; 107 reversals under A-C12; 166 of 166; 865 -> 95 of 1,181 drafts; 138 of 200 truth drafts; 218 rebuilt models of 20 genomes; 72 of 72; medians move <= 0.004; burden <= 1 (1.5 as shipped)): properties of the original A-C12 / A-C13 runs and side studies on the 10-29 sample (D9, not re-run). The re-run A-C13 pass is on a re-drawn sample (442 of 1,000 ATP and 62 of 200 gapfilled genomes kept), so only the counts above are current; the 'reproduces the pre-registered pass exactly (43,140 of 43,140 rows)' sentence is a statement about the 10-29 sample and cannot be restated for the new one |

## R1.5

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| tab:tau-genecap | 444-467 | cond 0.995 (min) | 69.1 | 69.1 | 69.0 | **CHANGED** |
| tab:tau-genecap | 444-467 | pop 0.987 (min) | 57.6 | 57.6 | 58.6 | **CHANGED** |
| tab:tau-genecap | 444-467 | pop 0.987 (max) | 64.4 | 64.4 | 65.2 | **CHANGED** |
| tab:tau-genecap | 444-467 | pop 0.995 (min) | 62.5 | 62.5 | 63.4 | **CHANGED** |
| tab:tau-genecap | 444-467 | pop 0.995 (max) | 69.0 | 69.0 | 69.7 | **CHANGED** |
| genecap_denominators | 397 (gene-capture part: 11.1%, 9.6%, 1.6%), 419 (9.6% -> 6.0%; 1.6% -> 1.4%) | score zero at 0.987 | 11.1 | 11.1 | 9.8 | **CHANGED** |
| genecap_denominators | 397 (gene-capture part: 11.1%, 9.6%, 1.6%), 419 (9.6% -> 6.0%; 1.6% -> 1.4%) | abstained at 0.987 | 9.6 | 9.6 | 9.5 | **CHANGED** |
| genecap_denominators | 397 (gene-capture part: 11.1%, 9.6%, 1.6%), 419 (9.6% -> 6.0%; 1.6% -> 1.4%) | pred_empty at 0.987 | 1.6 | 1.6 | 0.4 | **CHANGED** |
| genecap_denominators | 397 (gene-capture part: 11.1%, 9.6%, 1.6%), 419 (9.6% -> 6.0%; 1.6% -> 1.4%) | pred_empty at 0.995 | 1.4 | 1.4 | 0.4 | **CHANGED** |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | FL pop 0.987 | 63.3 | 63.3 | 64.0 | **CHANGED** |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | FL pop 0.995 | 68.7 | 68.7 | 69.2 | **CHANGED** |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | V4 pop 0.987 | 58.1 | 58.1 | 59.3 | **CHANGED** |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | V4 pop 0.995 | 63.1 | 63.1 | 64.0 | **CHANGED** |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | pop gain (max) | 5.5 | 5.5 | 5.4 | **CHANGED** |
| tab:tau-genecap | 444-467 | cond 0.987 (min) | 66.4 | 66.4 | 66.4 | changed below print precision |
| tab:tau-genecap | 444-467 | cond 0.987 (max) | 69.9 | 69.9 | 69.9 | changed below print precision |
| tab:tau-genecap | 444-467 | cond 0.995 (max) | 72.9 | 72.9 | 72.9 | changed below print precision |
| genecap_denominators | 397 (gene-capture part: 11.1%, 9.6%, 1.6%), 419 (9.6% -> 6.0%; 1.6% -> 1.4%) | abstained at 0.995 | 6.0 | 6.0 | 6.0 | changed below print precision |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | FL cond 0.987 | 68.7 | 68.7 | 68.7 | changed below print precision |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | FL cond 0.995 | 71.5 | 71.5 | 71.5 | changed below print precision |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | V4 cond 0.987 | 67.3 | 67.3 | 67.3 | changed below print precision |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | V4 cond 0.995 | 69.7 | 69.7 | 69.7 | changed below print precision |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | pop gain (min) | 4.1 | 4.1 | 4.1 | changed below print precision |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | resolved gain (min) | 2.4 | 2.4 | 2.4 | changed below print precision |
| tau0995_genecap | 394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns) | resolved gain (max) | 3.0 | 3.0 | 3.0 | changed below print precision |
| genecap_denominators | 397 (gene-capture part: 11.1%, 9.6%, 1.6%), 419 (9.6% -> 6.0%; 1.6% -> 1.4%) | (manual) | | | | **note**: the denominator itself moves: scoreable amplicons 67,038 -> 67,713 (truth_empty 953 -> 278) |

## R1.3L

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| select_v4_include_precision | 187 ("V4 0.795 against 0.856") | union V4 include precision | 0.795 | 0.795 | 0.794 | **CHANGED** |
| select_exclusion_precision_union_ge_legacy | 187 ("keeps gene-capture precision under self-exclusion at or above the ordered reducer's") | narrowest margin (min union - legacy) | 0 | 0.0016801410332314415 | 0.002484344039949238 | **CHANGED** |
| select_legacy_self_recovery_loss | 191 | failures | 12,602 | 12,602 | 12,127 | **CHANGED** |
| select_legacy_self_recovery_loss | 191 | fail rate | 19 | 19 | 18 | **CHANGED** |
| genecap_precision | 64, 187 (V4 0.795 vs 0.856), 401 (0.79 V4 vs 1.00 FL), 1175 (0.79-0.88 V3-onward; 0.92-0.94 V1-V2/V1-V3) | union V4 include precision | 0.795 | 0.795 | 0.794 | **CHANGED** |
| legacy_self_recovery | 187 (legacy mean 1.32), 191 (~19%, 12,602 of 67,466) | fail rate | 19 | 19 | 18 | **CHANGED** |
| legacy_self_recovery | 187 (legacy mean 1.32), 191 (~19%, 12,602 of 67,466) | failures | 12,602 | 12,602 | 12,127 | **CHANGED** |
| select_legacy_mean_size | 187 ("against 1.32 for the ordered reducer") | legacy include ok mean n_selected | 1.32 | 1.32 | 1.32 | changed below print precision |
| select_v4_include_precision | 187 ("V4 0.795 against 0.856") | legacy V4 include precision | 0.856 | 0.856 | 0.856 | changed below print precision |
| genecap_precision | 64, 187 (V4 0.795 vs 0.856), 401 (0.79 V4 vs 1.00 FL), 1175 (0.79-0.88 V3-onward; 0.92-0.94 V1-V2/V1-V3) | legacy V4 include precision | 0.856 | 0.856 | 0.856 | changed below print precision |
| genecap_precision | 64, 187 (V4 0.795 vs 0.856), 401 (0.79 V4 vs 1.00 FL), 1175 (0.79-0.88 V3-onward; 0.92-0.94 V1-V2/V1-V3) | union FL include precision | 1.00 | 1.00 | 1.00 | changed below print precision |
| genecap_precision | 64, 187 (V4 0.795 vs 0.856), 401 (0.79 V4 vs 1.00 FL), 1175 (0.79-0.88 V3-onward; 0.92-0.94 V1-V2/V1-V3) | union V4 include precision (2 dp) | 0.79 | 0.79 | 0.79 | changed below print precision |
| genecap_precision | 64, 187 (V4 0.795 vs 0.856), 401 (0.79 V4 vs 1.00 FL), 1175 (0.79-0.88 V3-onward; 0.92-0.94 V1-V2/V1-V3) | V3-onward (min) | 0.79 | 0.79 | 0.79 | changed below print precision |
| genecap_precision | 64, 187 (V4 0.795 vs 0.856), 401 (0.79 V4 vs 1.00 FL), 1175 (0.79-0.88 V3-onward; 0.92-0.94 V1-V2/V1-V3) | V3-onward (max) | 0.88 | 0.88 | 0.88 | changed below print precision |
| genecap_precision | 64, 187 (V4 0.795 vs 0.856), 401 (0.79 V4 vs 1.00 FL), 1175 (0.79-0.88 V3-onward; 0.92-0.94 V1-V2/V1-V3) | V1-anchored (min) | 0.92 | 0.92 | 0.92 | changed below print precision |
| genecap_precision | 64, 187 (V4 0.795 vs 0.856), 401 (0.79 V4 vs 1.00 FL), 1175 (0.79-0.88 V3-onward; 0.92-0.94 V1-V2/V1-V3) | V1-anchored (max) | 0.94 | 0.94 | 0.94 | changed below print precision |
| legacy_self_recovery | 187 (legacy mean 1.32), 191 (~19%, 12,602 of 67,466) | legacy include ok mean n_selected | 1.32 | 1.32 | 1.32 | changed below print precision |
| select_exclusion_precision_union_ge_legacy | 187 ("keeps gene-capture precision under self-exclusion at or above the ordered reducer's") | regions where union exclude precision >= legacy (count) | 8 | 8 | 8 |  |
| select_legacy_self_recovery_loss | 191 | self in record | 67,466 | 67,466 | 67,466 |  |
| legacy_self_recovery | 187 (legacy mean 1.32), 191 (~19%, 12,602 of 67,466) | self in record | 67,466 | 67,466 | 67,466 |  |

## R5A

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| AS-10 | 74; 521; 1187 ('+0.075 on matched amplicons (0.802 against 0.727)'); 1216 ('matched 0.802 against 0.727 as shipped') | matched amplicons | 59,607 | 59,607 | 61,044 | **CHANGED** |
| magval_identity_bins | 874 | amplicon answer share 0.945-0.987 (low) | 98 | 98 | 99 | **CHANGED** |
| magval_identity_bins | 874 | amplicon answer share 0.945-0.987 (high) | 98 | 98 | 99 | **CHANGED** |
| magval_scope_claims | 882-885 | amplicon F1 where MAG abstains, FL | 0.64 | 0.64 | 0.63 | **CHANGED** |
| magval_scope_claims | 882-885 | FL pop 16S | 0.605 | 0.605 | 0.612 | **CHANGED** |
| magval_scope_claims | 882-885 | FL pop MAG | 0.248 | 0.248 | 0.250 | **CHANGED** |
| magval_scope_claims | 882-885 | MAG-answered easier for amplicon | 0.053 | 0.053 | 0.055 | **CHANGED** |
| magval_tau99_paired | 892-895 | FL delta at tau 95 | +0.103 | +0.103 | +0.102 | **CHANGED** |
| magval_tau99_paired | 892-895 | pop MAG tau 95 | 0.248 | 0.248 | 0.250 | **CHANGED** |
| magval_tau99_paired | 892-895 | pop MAG tau 99 | 0.301 | 0.301 | 0.304 | **CHANGED** |
| magval_tau99_paired | 892-895 | matched-stringency gap | 0.354 | 0.354 | 0.357 | **CHANGED** |
| magval_tau99_paired | 892-895 | pop amplicon | 0.605 | 0.605 | 0.612 | **CHANGED** |
| magval_tau99_paired | 892-895 | FL pop F1 at 0.995 | 0.655 | 0.655 | 0.661 | **CHANGED** |
| magval_tau99_paired | 892-895 | conservative gap | 0.357 | 0.357 | 0.361 | **CHANGED** |
| panels.welch_tests | 1020 | FL original - k-means cond F1 | +0.011 | +0.011 | +0.012 | **CHANGED** |
| panels.welch_tests | 1020 | FL p | 2.3e-4 | 2.3e-4 | 1.5e-4 | **CHANGED** |
| panels.mag_singleton_abstention | 1083 | k-means singleton abstention | 0.6895 | 0.6895 | 0.6889 | **CHANGED** |
| panels.mag_singleton_abstention | 1083 | k-means n | 3,942 | 3,942 | 3,951 | **CHANGED** |
| panels.mag_singleton_abstention | 1083 | original singleton abstention | 0.6883 | 0.6883 | 0.6823 | **CHANGED** |
| panels.mag_singleton_abstention | 1083 | original n | 9,580 | 9,580 | 9,702 | **CHANGED** |
| panels.mag_singleton_abstention | 1083 | difference | +0.001 | +0.001 | +0.007 | **CHANGED** |
| panels.mag_singleton_abstention | 1083 | CI lo | -0.016 | -0.016 | -0.011 | **CHANGED** |
| panels.mag_singleton_abstention | 1083 | CI hi | +0.018 | +0.018 | +0.024 | **CHANGED** |
| panels.mag_multiplicity | 1084-1085 | k-means overall | 0.658 | 0.658 | 0.657 | **CHANGED** |
| panels.mag_multiplicity | 1084-1085 | original overall | 0.688 | 0.688 | 0.682 | **CHANGED** |
| panels.mag_multiplicity | 1084-1085 | k-means without >=20 | 0.687 | 0.687 | 0.686 | **CHANGED** |
| AS-10 | 74; 521; 1187 ('+0.075 on matched amplicons (0.802 against 0.727)'); 1216 ('matched 0.802 against 0.727 as shipped') | prFBA matched F1 | 0.802 | 0.802 | 0.802 | changed below print precision |
| AS-10 | 74; 521; 1187 ('+0.075 on matched amplicons (0.802 against 0.727)'); 1216 ('matched 0.802 against 0.727 as shipped') | PICRUSt2 matched F1 | 0.727 | 0.727 | 0.727 | changed below print precision |
| AS-10 | 74; 521; 1187 ('+0.075 on matched amplicons (0.802 against 0.727)'); 1216 ('matched 0.802 against 0.727 as shipped') | matched margin | +0.075 | +0.075 | +0.075 | changed below print precision |
| magval_identity_bins | 874 | ANI answer share (low) | 8 | 8 | 8 | changed below print precision |
| magval_identity_bins | 874 | ANI answer share (high) | 20 | 20 | 20 | changed below print precision |
| magval_scope_claims | 882-885 | conditional coverage | 31 | 31 | 31 | changed below print precision |
| magval_scope_claims | 882-885 | MAG abstains on | 69 | 69 | 69 | changed below print precision |
| magval_scope_claims | 882-885 | amplicon F1 where MAG abstains, V4 | 0.54 | 0.54 | 0.54 | changed below print precision |
| magval_scope_claims | 882-885 | MAG advantage (min over regions x levels) | 0.08 | 0.08 | 0.08 | changed below print precision |
| magval_scope_claims | 882-885 | MAG advantage (max over regions x levels) | 0.21 | 0.21 | 0.21 | changed below print precision |
| magval_tau99_paired | 892-895 | FL delta at tau 99 (d100 pairing) | +0.127 | +0.127 | +0.127 | changed below print precision |
| sharedref_tieset_coherence | 232 (first clause) | cross-phylum tie sets | 6.4 | 6.4 | 6.4 |  |
| sharedref_tieset_coherence | 232 (first clause) | same-genus tie sets | 35.8 | 35.8 | 35.8 |  |
| recon_modal_fractional_weight | 225 | modal k/n weight | 0.5 | 0.5 | 0.5 |  |
| panels.welch_tests | 1020 | all eight regions p < 0.01 | True | True | True |  |
| panels.welch_tests | 1020 | regions with p < 1e-8 | V4,V4-V5,V7-V9 | V4,V4-V5,V7-V9 | V4,V4-V5,V7-V9 |  |
| panels.welch_tests | 1020 | every MAG completeness p < 1e-3 | True | True | True |  |
| panels.mag_multiplicity | 1084-1085 | >=20 members abstention | 0.559 | 0.559 | 0.559 |  |

## R1.K

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| deployment_cross_phylum | 232 (deployment clause) | resolved merges | 1,238 | 1,238 | 1,238 |  |
| deployment_cross_phylum | 232 (deployment clause) | multi-genome merges | 1,410 | 1,410 | 1,410 |  |
| deployment_cross_phylum | 232 (deployment clause) | two-genome cross-phylum | 145 | 145 | 145 |  |
| deployment_cross_phylum | 232 (deployment clause) | two-genome resolved | 454 | 454 | 454 |  |
| deployment_cross_phylum | 232 (deployment clause) | two-genome pct | 31.9 | 31.9 | 31.9 |  |
| deployment_cross_phylum | 232 (deployment clause) | larger cross-phylum | 488 | 488 | 488 |  |
| deployment_cross_phylum | 232 (deployment clause) | larger resolved | 784 | 784 | 784 |  |
| deployment_cross_phylum | 232 (deployment clause) | larger pct | 62.2 | 62.2 | 62.2 |  |
| deployment_cross_phylum | 232 (deployment clause) | Bacteria+Archaea | 20 | 20 | 20 |  |
| e3_production_unresolved_share | 253 ("11% of the features") | unresolved feature share | 11 | 11 | 11 |  |
| e3_production_unresolved_share | 253 ("11% of the features") | multi-genome synthetic genomes | 1,410 | 1,410 | 1,410 |  |
| deployment_cross_phylum | 232 (1,238 of 1,410; 31.9% = 145/454; 62.2% = 488/784; 20 Bacteria+Archaea), 1171, 1236 | resolved merges | 1,238 | 1,238 | 1,238 |  |
| deployment_cross_phylum | 232 (1,238 of 1,410; 31.9% = 145/454; 62.2% = 488/784; 20 Bacteria+Archaea), 1171, 1236 | multi-genome merges | 1,410 | 1,410 | 1,410 |  |
| deployment_cross_phylum | 232 (1,238 of 1,410; 31.9% = 145/454; 62.2% = 488/784; 20 Bacteria+Archaea), 1171, 1236 | two-genome cross-phylum | 145 | 145 | 145 |  |
| deployment_cross_phylum | 232 (1,238 of 1,410; 31.9% = 145/454; 62.2% = 488/784; 20 Bacteria+Archaea), 1171, 1236 | two-genome resolved | 454 | 454 | 454 |  |
| deployment_cross_phylum | 232 (1,238 of 1,410; 31.9% = 145/454; 62.2% = 488/784; 20 Bacteria+Archaea), 1171, 1236 | two-genome pct | 31.9 | 31.9 | 31.9 |  |
| deployment_cross_phylum | 232 (1,238 of 1,410; 31.9% = 145/454; 62.2% = 488/784; 20 Bacteria+Archaea), 1171, 1236 | larger cross-phylum | 488 | 488 | 488 |  |
| deployment_cross_phylum | 232 (1,238 of 1,410; 31.9% = 145/454; 62.2% = 488/784; 20 Bacteria+Archaea), 1171, 1236 | larger resolved | 784 | 784 | 784 |  |
| deployment_cross_phylum | 232 (1,238 of 1,410; 31.9% = 145/454; 62.2% = 488/784; 20 Bacteria+Archaea), 1171, 1236 | larger pct | 62.2 | 62.2 | 62.2 |  |
| deployment_cross_phylum | 232 (1,238 of 1,410; 31.9% = 145/454; 62.2% = 488/784; 20 Bacteria+Archaea), 1171, 1236 | Bacteria+Archaea | 20 | 20 | 20 |  |
| deployment_e3_unresolved_share | 253 | unresolved feature share | 11 | 11 | 11 |  |
| deployment_e3_unresolved_share | 253 | ASVs | 3,950 | 3,950 | 3,950 |  |
| disc.selection_size_deployment | 1170, 1229 | single-genome share | 62.5 | 62.5 | 62.5 |  |
| disc.selection_size_deployment | 1170, 1229 | single-genome ASVs | 2,345 | 2,345 | 2,345 |  |
| disc.selection_size_deployment | 1170, 1229 | answered | 3,755 | 3,755 | 3,755 |  |
| disc.selection_size_deployment | 1170, 1229 | mean selection size | 2.289 | 2.289 | 2.289 |  |
| disc.phylum_span | 1171, 1236 | resolved merges | 1,238 | 1,238 | 1,238 |  |
| disc.phylum_span | 1171, 1236 | multi-genome merges | 1,410 | 1,410 | 1,410 |  |
| disc.phylum_span | 1171, 1236 | two-genome cross-phylum | 145 | 145 | 145 |  |
| disc.phylum_span | 1171, 1236 | two-genome resolved | 454 | 454 | 454 |  |
| disc.phylum_span | 1171, 1236 | two-genome pct | 31.9 | 31.9 | 31.9 |  |
| disc.phylum_span | 1171, 1236 | larger cross-phylum | 488 | 488 | 488 |  |
| disc.phylum_span | 1171, 1236 | larger resolved | 784 | 784 | 784 |  |
| disc.phylum_span | 1171, 1236 | larger pct | 62.2 | 62.2 | 62.2 |  |
| lim.deployment_truncation | 1207-1208 | single-genome share | 62.5 | 62.5 | 62.5 |  |
| lim.deployment_truncation | 1207-1208 | mean selection size | 2.289 | 2.289 | 2.289 |  |
| lim.deployment_truncation | 1207-1208 | largest selection | 77 | 77 | 77 |  |
| lim.deployment_truncation | 1207-1208 | ASVs | 3,950 | 3,950 | 3,950 |  |
| lim.deployment_truncation | 1207-1208 | (manual) | | | | **unchanged** (25 ASVs; 2.034 / 19 depth-20 comparator): deployment builds were not re-run (D4); deployment_d100/RESULTS.md |

## R2.1

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| e3_gate_brier_auroc | 254 | k/n Brier | 0.155 | 0.155 | 0.154 | **CHANGED** |
| e3_gate_brier_auroc | 254 | delta Brier | -0.053 | -0.053 | -0.052 | **CHANGED** |
| e3_gate_brier_auroc | 254 |   CI lo | -0.055 | -0.055 | -0.054 | **CHANGED** |
| e3_gate_brier_auroc | 254 |   CI hi | -0.051 | -0.051 | -0.050 | **CHANGED** |
| e3_gate_brier_auroc | 254 | k/n AUROC | 0.874 | 0.874 | 0.878 | **CHANGED** |
| e3_gate_brier_auroc | 254 | delta AUROC | +0.054 | +0.054 | +0.050 | **CHANGED** |
| e3_gate_brier_auroc | 254 |   CI lo | +0.052 | +0.052 | +0.048 | **CHANGED** |
| e3_gate_brier_auroc | 254 |   CI hi | +0.056 | +0.056 | +0.052 | **CHANGED** |
| e3_gate_brier_auroc | 254 | amplicons | 19,158 | 19,158 | 19,542 | **CHANGED** |
| e3_gate_brier_auroc | 254 | source genomes | 6,776 | 6,776 | 6,891 | **CHANGED** |
| e3_gate_recalibration_share | 255 | recalibration share of the Brier gain | 45 | 44 | 48 | **CHANGED**; old≠printed (rounding; known: 'about 45%' is 0.0236/0.0530 from rounded inputs; the unrounded share is 44.5%) |
| e3_gate_recalibration_share | 255 | E3 vs saturated k/n Brier | 0.029 | 0.029 | 0.027 | **CHANGED** |
| e3_gate_recalibration_share | 255 | E3 vs saturated k/n AUROC | 0.048 | 0.048 | 0.044 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | k/n Brier | 0.155 | 0.155 | 0.154 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | delta Brier | -0.053 | -0.053 | -0.052 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') |   CI lo | -0.055 | -0.055 | -0.054 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') |   CI hi | -0.051 | -0.051 | -0.050 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | k/n AUROC | 0.874 | 0.874 | 0.878 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | delta AUROC | +0.054 | +0.054 | +0.050 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') |   CI lo | +0.052 | +0.052 | +0.048 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') |   CI hi | +0.056 | +0.056 | +0.052 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | amplicons | 19,158 | 19,158 | 19,542 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | source genomes | 6,776 | 6,776 | 6,891 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | recalibration share of the Brier gain | 45 | 44 | 48 | **CHANGED**; old≠printed (rounding; known: 'about 45%' is 0.0236/0.0530 from rounded inputs; the unrounded share is 44.5%) |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | E3 vs saturated k/n Brier | 0.029 | 0.029 | 0.027 | **CHANGED** |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | E3 vs saturated k/n AUROC | 0.048 | 0.048 | 0.044 | **CHANGED** |
| e3_gate_brier_auroc | 254 | E3 Brier | 0.102 | 0.102 | 0.102 | changed below print precision |
| e3_gate_brier_auroc | 254 | E3 AUROC | 0.928 | 0.928 | 0.928 | changed below print precision |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | E3 Brier | 0.102 | 0.102 | 0.102 | changed below print precision |
| e3_production_gate | 254-255, 361, 1225 ('checked on the 10,000-organism panel') | E3 AUROC | 0.928 | 0.928 | 0.928 | changed below print precision |

## R1.4

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| legacy_genecap_d20_vs_d100 | 732, 1204 (conditional gene capture +0.7 to +1.3 points; conditional gene-capture F1 -0.1 to -1.2 points) | cond gene-capture F1 d20-d100 (max drop, points) | 1.2 | 1.2 | 1.1 | **CHANGED** |
| legacy_genecap_d20_vs_d100 | 732, 1204 (conditional gene capture +0.7 to +1.3 points; conditional gene-capture F1 -0.1 to -1.2 points) | cond macro recall d100-d20 (min, points) | 0.7 | 0.7 | 0.7 | changed below print precision |
| legacy_genecap_d20_vs_d100 | 732, 1204 (conditional gene capture +0.7 to +1.3 points; conditional gene-capture F1 -0.1 to -1.2 points) | cond macro recall d100-d20 (max, points) | 1.3 | 1.3 | 1.3 | changed below print precision |
| legacy_genecap_d20_vs_d100 | 732, 1204 (conditional gene capture +0.7 to +1.3 points; conditional gene-capture F1 -0.1 to -1.2 points) | cond gene-capture F1 d20-d100 (min drop, points) | 0.1 | 0.1 | 0.1 | changed below print precision |

## D7

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| backend_concordance | 384-386 | (manual) | | | | **historical (D7)**: June check; wording only ('on the cache of the time') |

## R2.2

| item | tex lines | value | printed | old | new | flag |
|---|---|---|---|---|---|---|
| magval_tau99_mag_capture | 889-891 | resolve tau 95 | 0.308 | 0.308 | 0.313 | **CHANGED** |
| magval_tau99_mag_capture | 889-891 | resolve tau 99 | 0.366 | 0.366 | 0.372 | **CHANGED** |
| magval_tau99_mag_capture | 889-891 | cond F1 tau 99 | 0.823 | 0.823 | 0.822 | **CHANGED** |
| twopanel_mag_depth500_original | 943-945; 1009-1013 | resolve depth 100 (pooled 100-30%) | 0.3071 | 0.3071 | 0.3127 | **CHANGED** |
| twopanel_mag_depth500_original | 943-945; 1009-1013 | resolve depth 500 (pooled 100-30%) | 0.3089 | 0.3089 | 0.3144 | **CHANGED** |
| twopanel_mag_depth500_original | 943-945; 1009-1013 | cond F1 depth 100 (pooled, unweighted) | 0.8066 | 0.8066 | 0.8063 | **CHANGED** |
| twopanel_mag_depth500_original | 943-945; 1009-1013 | cond F1 depth 500 (pooled, resolved-weighted) | 0.8061 | 0.8061 | 0.8059 | **CHANGED** |
| twopanel_mag_depth500_original | 943-945; 1009-1013 | scoreable cells per regime | 49,185 | 49,185 | 49,795 | **CHANGED** |
| panels.mag_depth500_rerun | 1011-1013 | resolve depth 100 (pooled 100-30%) | 0.3071 | 0.3071 | 0.3127 | **CHANGED** |
| panels.mag_depth500_rerun | 1011-1013 | resolve depth 500 (pooled 100-30%) | 0.3089 | 0.3089 | 0.3144 | **CHANGED** |
| panels.mag_depth500_rerun | 1011-1013 | cond F1 depth 100 (pooled, unweighted) | 0.8066 | 0.8066 | 0.8063 | **CHANGED** |
| panels.mag_depth500_rerun | 1011-1013 | cond F1 depth 500 (pooled, resolved-weighted) | 0.8061 | 0.8061 | 0.8059 | **CHANGED** |
| panels.mag_depth500_rerun | 1011-1013 | scoreable cells per regime | 49,185 | 49,185 | 49,795 | **CHANGED** |
| magval_selfinclusion_ceiling | 822 | include anchor is the source | 99.9 | 99.9 | 99.9 | changed below print precision |
| magval_tau99_mag_capture | 889-891 | relative gain | +19 | +19 | +19 | changed below print precision |
| magval_tau99_mag_capture | 889-891 | cond F1 tau 95 | 0.807 | 0.807 | 0.807 | changed below print precision |
| magval_selfinclusion_ceiling | 822 | include coverage 100% | 1.000 | 1.000 | 1.000 |  |
| magval_selfinclusion_ceiling | 822 | include gene capture 100% | 1.000 | 1.000 | 1.000 |  |
| magval_selfinclusion_ceiling | 822 | include genus accuracy (min over completeness) | 1.000 | 1.000 | 1.000 |  |
| magval_tau99_mag_capture | 889-891 | genus accuracy tau 95 | 0.910 | 0.910 | 0.910 |  |
| magval_tau99_mag_capture | 889-891 | genus accuracy tau 99 | 0.895 | 0.895 | 0.895 |  |
| twopanel_mag_depth500_original | 943-945; 1009-1013 | cells changing outcome | 88 | 88 | 88 |  |
| twopanel_mag_depth500_original | 943-945; 1009-1013 | newly resolved | 87 | 87 | 87 |  |
| twopanel_mag_depth500_original | 943-945; 1009-1013 | gene-less reference | 1 | 1 | 1 |  |
| twopanel_mag_depth500_original | 943-945; 1009-1013 | selection switches | 0 | 0 | 0 |  |
| panels.mag_depth500_rerun | 1011-1013 | cells changing outcome | 88 | 88 | 88 |  |
| panels.mag_depth500_rerun | 1011-1013 | newly resolved | 87 | 87 | 87 |  |
| panels.mag_depth500_rerun | 1011-1013 | gene-less reference | 1 | 1 | 1 |  |
| panels.mag_depth500_rerun | 1011-1013 | selection switches | 0 | 0 | 0 |  |
| magval_tau99_mag_capture | 889-891 | (manual) | | | | **unchanged** (1,285 retained references): an ANI/skani property of the tau-99 search (no gene data) |
