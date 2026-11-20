
### R1.3U (11 items)

- `front:abs-selfincl-genecap` L64: Source in database: gene capture over 99%; larger tie sets lower gene-capture precision on short amplicons (0.79 on V4).
- `front:abs-selfexcl-genecap` L65: Self-exclusion (0.987): gene capture 66-70% over amplicons resolved to a reference that carries gene families; 58-64% over the whole novel-organism population.
- `front:abs-tier-genecap` L73: Answers resolved only at family tier recover 58% of the true gene complement (56-61% by region) vs 69% for genus-tier answers.
- `front:tab:main-genecap` L420-442: tab:main gene capture incl. (99.1-99.98%) and excl. cond/pop (66.4-69.9% / 57.6-64.4%).
- `front:fig:genecap` L477-483: F4b_percent_genecapture_by_region.pdf: include-self capture >99% (ticks); self-exclusion population columns 58-64%; resolved 66-70%.
- `front:fig:leakage-genecap` L485-491: F11 caption: gene-capture gap (not shown) 30-33 points.
- `val16s:tabmain_genecap` L64-65, 401-403, 432-439 (gene capture incl. and excl. cond./pop. columns), 480-481 (fig:genecap caption + F4b image), 1175 (57.6-59.2%; 63.3-64.4%): Self-inclusion capture >99% (FL 99.98%, V4 99.1%); self-exclusion conditional 66-70% (V1-V2 69.9, V1-V3 69.5 vs FL 68.7; V4-V5 66.4 lowest), population 58-64% (
- `val16s:genecap_leakage_gap` L404 (30-33 points; 29.6 V1-V2 to 32.9 V4-V5), 489 (fig:leakage caption, capture gap): Gene-capture leakage gap (include minus conditional exclude) 30-33 points, 29.6 (V1-V2) to 32.9 (V4-V5)
- `val16s:genecap_by_tier` L73, 409: Family-tier answers recover 58% of the true gene complement under self-exclusion (56-61% by region) vs 69% for genus-tier answers
- `mag:methods_union_size_amplicon` L178: The amplicon arm's exact-tie union holds more than one genome for 34.0% of its self-inclusion selections on the same panel.
- `panels:disc.region_tradeoff_genecap` L1175: Self-inclusion gene capture over 99% with macro precision 0.79-0.88 (V3 onward) and 0.92-0.94 (V1-V2, V1-V3); self-excluded population gene capture 57.6-59.2% (

### R1.3L (6 items)

- `methods:select_legacy_mean_size` L187 ("against 1.32 for the ordered reducer"): The ordered (legacy) reducer selects a mean 1.32 genomes on the same panel.
- `methods:select_v4_include_precision` L187 ("V4 0.795 against 0.856"): Short-amplicon self-inclusion gene-capture precision falls under the union: V4 0.795 (union) against 0.856 (legacy).
- `methods:select_exclusion_precision_union_ge_legacy` L187 ("keeps gene-capture precision under self-exclusion at or above the ordered reducer's"): Union conditional precision under self-exclusion is >= the ordered reducer's in all 8 regions.
- `methods:select_legacy_self_recovery_loss` L191: The ordered reducer discarded the correct source in ~19% of self-recovery round-trips (12,602 of the 67,466 self-inclusion amplicons whose source is in the 100-
- `val16s:genecap_precision` L64, 187 (V4 0.795 vs 0.856), 401 (0.79 V4 vs 1.00 FL), 1175 (0.79-0.88 V3-onward; 0.92-0.94 V1-V2/V1-V3): Self-inclusion gene-capture macro precision: V4 0.795 under the union (0.856 under the ordered reducer), FL 1.00; 0.79-0.88 for V3-onward regions, 0.92-0.94 for
- `val16s:legacy_self_recovery` L187 (legacy mean 1.32), 191 (~19%, 12,602 of 67,466): The ordered (legacy) reducer discards the source in ~19% of self-recovery round-trips (12,602 of the 67,466 self-inclusion amplicons whose source is in the 100-

### R1.4 (1 items)

- `val16s:legacy_genecap_d20_vs_d100` L732, 1204 (conditional gene capture +0.7 to +1.3 points; conditional gene-capture F1 -0.1 to -1.2 points): Under the legacy reducer conditional gene capture (macro recall) rises 0.7-1.3 points from depth 20 to 100 while conditional gene-capture F1 falls 0.1-1.2 point

### R1.5 (3 items)

- `front:tab:tau-genecap` L444-467: tab:tau capture cond (66.4-69.9% / 69.1-72.9%) and pop (57.6-64.4% / 62.5-69.0%) at 0.987 / 0.995.
- `val16s:genecap_denominators` L397 (gene-capture part: 11.1%, 9.6%, 1.6%), 419 (9.6% -> 6.0%; 1.6% -> 1.4%): 11.1% of scoreable amplicons score zero under self-exclusion: selection keeps no genome (9.6%) or only references without gene families (1.6%); at 0.995 these f
- `val16s:tau0995_genecap` L394, 419 (+4.1-5.5 pop. -> 62.5-69.0%, V4 58.1->63.1; +2.4-3.0 resolved -> 69.1-72.9%), 456-463 (tab:tau capture cond./pop. columns): Gene capture at 0.987 / 0.995 per region (e.g. FL 68.7/71.5 cond., 63.3/68.7 pop.; V4 67.3/69.7, 58.1/63.1)

### R1.K (7 items)

- `methods:deployment_cross_phylum` L232 (deployment clause): Deployment: among 1,238 of 1,410 multi-genome merges with >=2 lineage-resolved members, 31.9% of two-genome (145/454) and 62.2% of larger (488/784) merges span 
- `methods:e3_production_unresolved_share` L253 ("11% of the features"): A merged function without a PGFam keeps Eq. frequency: 11% of the features of the 1,410 multi-genome synthetic genomes of the 3,950-ASV deployment.
- `val16s:deployment_cross_phylum` L232 (1,238 of 1,410; 31.9% = 145/454; 62.2% = 488/784; 20 Bacteria+Archaea), 1171, 1236: Among 1,238 of 1,410 multi-genome merges that resolve to an NCBI lineage, 31.9% of two-genome and 62.2% of larger merges span phyla; 20 unite Bacteria with Arch
- `val16s:deployment_e3_unresolved_share` L253: A merged function without a PGFam keeps p_func: 11% of the features of the 1,410 multi-genome synthetic genomes of the 3,950-ASV deployment
- `panels:disc.selection_size_deployment` L1170, 1229: 62.5% of deployment ASVs get a single genome (2,345 of 3,755 answered; mean 2.289).
- `panels:disc.phylum_span` L1171, 1236: In deployment 31.9% of two-genome (145/454) and 62.2% of larger (488/784) merges span phyla, among 1,238 of 1,410 multi-genome merges resolving to an NCBI linea
- `panels:lim.deployment_truncation` L1207-1208: Deployment: in 25 of 3,950 ASVs the 100th retained hit ties the best, so those selections (7 to 77 genomes) are lower bounds; vs the same hits cut to 20, depth 

### R2.1 (3 items)

- `methods:e3_gate_brier_auroc` L254: Released E3 beat Eq. frequency on the 10k panel at self-exclusion 0.987, n>=2: Brier 0.155 -> 0.102 (-0.053 [-0.055, -0.051]); AUROC 0.874 -> 0.928 (+0.054 [+0.
- `methods:e3_gate_recalibration_share` L255: About 45% of the Brier gain is recalibration; against a saturated in-sample recalibration of Eq. frequency the estimate still lowers Brier by 0.029 and raises A
- `val16s:e3_production_gate` L254-255, 361, 1225 ('checked on the 10,000-organism panel'): On the 10k panel, exclude_0.987, n>=2 (19,158 amplicons, 6,776 sources): union PGFam Brier 0.155 -> 0.102 (-0.053 [-0.055, -0.051]); AUROC 0.874 -> 0.928 (+0.05

### R2.2 (4 items)

- `mag:magval_selfinclusion_ceiling` L822: Under self-inclusion the arm is at ceiling at every level: coverage 1.000, genus accuracy 1.000, gene capture 1.000; selected reference is the source genome in 
- `mag:magval_tau99_mag_capture` L889-891: tau=99 lifts the resolve rate 0.308 -> 0.366 (+19% relative) and conditional gene-capture F1 0.807 -> 0.823; completeness stays flat across the five main-grid t
- `mag:twopanel_mag_depth500_original` L943-945; 1009-1013: [sec:twopanel, cross-reference] Original panel MAG arm re-run at skani depth 500: over 49,185 scoreable cells per regime, self-exclusion resolve 0.3071 -> 0.308
- `panels:panels.mag_depth500_rerun` L1011-1013: Original MAG panel re-run at depth 500: over 49,185 scoreable cells per regime the self-exclusion resolve rate moves 0.3071 -> 0.3089 and conditional F1 0.8066 

### R2.3 (42 items)

- `front:abs-mag-vs-16s` L63: ANI front end on the same panel captures more true gene complement than every amplicon region wherever it answers: +0.10 to +0.21 F1 at full completeness, +0.08
- `front:abs-threepanel-genecap-quality` L69: Phylum-standardized conditional gene-capture F1 differs from the stratified panel by at most 0.028.
- `front:abs-threepanel-genecap-coverage` L70: Self-exclusion gene-capture coverage on the random panel falls by 0.24 to 0.29 in every region (upper end of 'by 0.20 to 0.29').
- `front:abs-selfincl-panel-dependence` L71: Self-inclusion quality falls on the random panel by up to 0.11 in short regions and 0.02 on full-length 16S.
- `front:tab:magvs16s` L842-865: tab:magvs16s: per region n (1,885-2,334), conditional MAG/16S F1 (0.808-0.815 / 0.598-0.710), Delta +0.103..+0.212 (median +0.052..+0.171), coverage MAG/16S (0.
- `front:tab:magval` L896-919: tab:magval: 9,997 panel organisms; resolve rate 0.308..0.253 and cond F1 0.807..0.783 at 100%..1% completeness; MAG/16S vs FL (0.813/0.710 +0.103 .. 0.783/0.706
- `front:fig:magcomp` L921-929: F14_mag_completeness.pdf: (A) aligned fractions vs completeness (2,000-organism GTDB subsample), (B) resolve rate flat in completeness, (C) MAG arm leads both a
- `front:fig:magcross` L931-938: F15_mag_crossover.pdf: no crossover down to 1% completeness; conditional F1 vs the two amplicon constants; resolve rate degrades.
- `front:tab:twopanel-16S` L1025-1050: tab:twopanel 16S block: crude coverage (FL 0.921/0.750/0.650; V1-V3 0.927/0.757/0.664; V4 0.863/0.647/0.536; V7-V9 0.880/0.659/0.547) and conditional gene-captu
- `front:tab:twopanel-MAG` L1025-1050: tab:twopanel MAG block (depth 500, panel as reference): resolve 0.309/0.342/0.428 (100%), 0.309/0.343/0.428 (50%), 0.310/0.345/0.430 (30%); cond F1 0.807/0.797/
- `front:tab:emptied` L1065-1080: tab:emptied (FL 16S, tau 0.987): n 7,215/8,815/8,762; share of queries by hits deleted (0-19: 84.4/44.0/37.3%; 20-49; 50-99; emptied 4.7/22.9/33.4%) and uncover
- `front:tab:threepanel-genecap` L1101-1137: tab:threepanel 16S gene-capture block: FL 0.926/0.767/0.682, V4 0.855/0.665/0.576, V7-V9 0.870/0.681/0.594; Delta coverage -0.243/-0.279/-0.275; Delta cond F1 -
- `front:tab:threepanel-panelref` L1101-1137: tab:threepanel panel-as-reference block: cond gene-capture F1 0.809/0.788/0.783 (100%), 0.806/0.786/0.789 (30%); Delta -0.026/-0.017.
- `val16s:emptied_vs_uncovered_panels` L742: Full-length 16S emptied-record rate 4.7% (original) vs 22.9% (k-means) (8.3% vs 29.4% pooled over 8 regions); uncovered survivors 3.4% vs 2.8% (3.1% vs 2.4% poo
- `mag:abstract_mag_margin` L63: Abstract: the ANI front end captures more of the true gene complement than every one of the eight amplicon regions wherever it answers (+0.10 to +0.21 F1 at ful
- `mag:magval_completeness_null` L823; tab:magval 896-919 resolve-rate and Cond. F1 columns: Under self-exclusion the resolve rate is 0.306-0.308 and conditional gene-capture F1 0.805-0.807 from 100% to 30%; tab:magval resolve 0.308/0.307/0.309/0.295/0.
- `mag:magval_paired_main` L824; tab:magval 896-919 (MAG/16S, delta, n columns); tab:magval caption 897-900 (V4 pairing 1,893-2,311): Paired within organism (depth-100 union amplicon cells), the MAG arm leads full-length 16S by 0.103 and V4 by 0.212 F1, with the same margin at 30% as at 100%; 
- `mag:magval_no_crossover` L826-830: At 1% completeness (37 kb, fifty-fold below MIMAG) the MAG arm still leads FL by 0.077 and V4 by 0.190 F1, winning 1,129 organisms to 473 (p~1e-79); degrading t
- `mag:magval_vs_every_region` L832-835: Among amplicons both arms resolve, MAG captures more than every region: +0.103 to +0.116 F1 vs FL, V1-V2, V1-V3; +0.168 to +0.212 vs the five V3-onward regions 
- `mag:tab_magvs16s` L842-868: tab:magvs16s: per region n, conditional MAG/16S F1, delta (median delta), coverage MAG/16S, population MAG/16S (e.g. FL 1,916 0.813 0.710 +0.103 (+0.056) 0.304 
- `mag:magval_population_view` L837-840: Population view reverses the ranking: MAG answers ~30% of organisms vs 86-93% for an amplicon; 0.605 vs 0.248 (FL), 0.483 vs 0.253 (V4).
- `mag:magval_af_gate_counterfactual` L878-880: Relaxing the ANI floor 0.83 -> 0 moves the share of scoreable organisms with a qualifying surviving hit 0.311 -> 0.316; relaxing the AF gate 15% -> 2% moves it 
- `mag:fig_magcomp_magcross` L921-938: F14 (A: AF_ref vs AF_query medians over the 2,000-organism GTDB subsample; B: resolve rate flat in completeness; C: MAG leads both amplicons at every completene
- `mag:twopanel_mag_panelarm` L1018-1022; tab:twopanel MAG rows 1045-1048; tab:threepanel panel-arm rows 1132-1134: [sec:twopanel, cross-reference] MAG conditional gene-capture F1 within 0.016 (random) / 0.010 (k-means) of the original at every completeness, Welch p<1e-3; MAG
- `mag:twopanel_mag_abstention_strata` L1085-1087: [sec:twopanel, cross-reference] Among organisms that are the only member of their species, MAG abstention 0.6895 (k-means, n=3,942) vs 0.6883 (original, n=9,580
- `panels:panels.condF1_crude_bounds` L1018-1019: Crude conditional gene-capture F1: random below original in every amplicon region by at most 0.029 (V7-V9) and within 0.016 at every MAG completeness; k-means w
- `panels:panels.coverage_disagrees` L1021-1023: Alternative panels resolve more MAGs (0.342 and 0.428 vs 0.309 crude at full completeness) and fewer amplicons (0.750 and 0.650 vs 0.921, full-length 16S).
- `panels:panels.tab_twopanel_16S` L1025-1043: tab:twopanel 16S block (crude coverage | conditional gene-capture F1, original/k-means/random): FL 0.921/0.750/0.650 | 0.657/0.645/0.641; V1-V3 0.927/0.757/0.66
- `panels:panels.tab_twopanel_MAG` L1044-1048: tab:twopanel MAG block (resolve | cond F1, original/k-means/random): 100% 0.309/0.342/0.428 | 0.807/0.797/0.791; 50% 0.309/0.343/0.428 | 0.806/0.798/0.792; 30% 
- `panels:panels.fl_gap_standardized` L1055-1056: Once a PGFam cache gap in the k-means run was closed both panels are scored against the same cache state; phylum standardization leaves most of the FL coverage 
- `panels:panels.surviving_uncovered` L1057-1060: Surviving full-length records uncovered: 3.4% / 2.8% / 2.3% (original/k-means/random); 1.5-3.7% in every FL stratum; 3.1% / 2.4% / 1.8% pooled over eight region
- `panels:panels.emptied_arithmetic` L1061-1063: Emptied by construction: uncovered 1.000; cross-panel coverage difference = difference in how often the record empties: 33.4% random, 22.9% k-means, 4.7% origin
- `panels:panels.surviving_pred_empty` L1064: Uncovered surviving full-length records include cells whose selected references carry no gene families: 107 / 76 / 62 (original / k-means / random).
- `panels:panels.tab_emptied_shares` L1066-1081: tab:emptied population and share columns: n = 7,215 / 8,815 / 8,762; hits deleted 0-19: 84.4% / 44.0% / 37.3%; 20-49: 7.2% / 17.9% / 14.8%; 50-99: 3.8% / 15.2% 
- `panels:panels.tab_emptied_uncovered` L1066-1081: tab:emptied uncovered columns: 0-19: 0.037 / 0.028 / 0.023; 20-49: 0.017 / 0.022 / 0.021; 50-99: 0.015 / 0.032 / 0.028; emptied 1.000 / 1.000 / 1.000.
- `panels:panels.std_gap_share` L1095-1096: In every region, and for gene-capture coverage, the standardized gap keeps 76-92% of the crude one and its interval excludes zero.
- `panels:panels.cond_f1_random_gap` L1097-1098: Conditional gene-capture F1 lower on the random panel by 0.013 to 0.028 (phylum-standardized), an order of magnitude less than the coverage gap.
- `panels:panels.tab_threepanel_genecap_block` L1120-1124: tab:threepanel gene-capture block (phylum-std coverage original/k-means/random; Delta coverage [CI]; Delta conditional F1 [CI]): FL 0.926/0.767/0.682, -0.243 [-
- `panels:panels.tab_threepanel_magpanel_block` L1130-1133: tab:threepanel MAG panel-as-reference block (phylum-std conditional gene-capture F1; Delta random-original [CI]): 100% 0.809/0.788/0.783, -0.026 [-0.032,-0.020]
- `panels:panels.taxonomy_bootstrap_CIs` L1094, 1095, 1097, 1115-1119: Bootstrap intervals and 'interval excludes zero' verdicts on taxonomy endpoints: V4 genus-coverage Delta (-0.276 to -0.238); tab:threepanel genus-block coverage
- `panels:panels.selfincl_genecap_penalty` L1149-1150: Standardized self-inclusion gene-capture F1 falls from 0.837 to 0.730 on V4 (-0.107, -0.114 to -0.100) and from 0.906 to 0.805 on V7-V9; on full-length 16S gene
- `panels:panels.what_this_licenses` L1160-1163: Standardized self-exclusion genus and gene-capture coverage falls by 0.20 to 0.29 against the stratified panel in every region; dense genome-scale reference mov

### R2.4 (22 items)

- `front:abs-picrust2-asshipped` L74: As shipped, prFBA self-excluded conditional EC F1 0.802 vs PICRUSt2 0.727 on amplicons both answer; no EC prediction for 11.1% of EC-scoreable amplicons (9.6% a
- `front:fig:picrust2` L570-579: F13_prfba_vs_picrust2.pdf: (A) PICRUSt2 precision 0.735 -> 0.235 across NSTI bins, prFBA to 0.606; (B) below the 0.865 family floor prFBA abstains, PICRUSt2 pre
- `val16s:ec_coverage_bins_and_noanswer` L74 (11.1%; 9.6%), 715 (EC coverage 0.982 and 0.829), 1179 (answers 89%), 1181 (11%; 9.6%; 1.5%): EC coverage of the self-excluded arm 0.982 (0.945-0.987 bin) and 0.829 (>=0.995); prFBA answers 89% of self-excluded amplicons; 11.1% of EC-scoreable amplicons 
- `picrust:AS-04` L516; 544 (n=67,038); 1179 (67,038): Unscoreable amplicons: 953 amplicons (from 161 source genomes) whose source genome has no gene-family annotation in the BV-BRC cache are excluded; scoreable n =
- `picrust:AS-05` L517: Cache-repair disclosure: empty family lists behind 675 of the 953 unscoreable amplicons and behind the references of 818 of the 1,060 self-excluded amplicons sc
- `picrust:AS-08` L521; 1179: PICRUSt2 coverage: only 26 of 67,038 scoreable amplicons (27 of all 67,991) fall below its alignment threshold.
- `picrust:AS-09` L521; 525; 544-547; 555-557 (tab:picrust2 'All amplicons' rows); 1179 ('89%'); 1187 ('-0.016 ... 0.713 against 0.729'); 1216 ('from -0.016'); 733 ('0.713, still below PICRUSt2's 0.729'): Headline as shipped: PICRUSt2 coverage 1.000, recall 0.849, precision 0.669, F1 0.729, population F1 0.729 (0.729 over all 67,012 it answers); prFBA exclude-sel
- `picrust:AS-11` L74; 1181: Non-answers: prFBA gives no EC prediction for 11.1% of EC-scoreable amplicons (11% in Discussion): 9.6% abstain outright, the other 1.5% select references that 
- `picrust:AS-12` L521; 576 (fig C): Per-region paired Wilcoxon on per-amplicon F1: advantage in all eight regions, median DeltaF1 +0.072 to +0.095, Holm p ~ 0; prFBA wins between three and seven a
- `picrust:AS-13` L81 ('+0.046 to -0.202'); 523; 1217: Annotation dialect: on the 19,686 amplicons whose source genome is in both databases, scoring against PICRUSt2's own annotation instead of BV-BRC's moves condit
- `picrust:AS-14` L82 ('0.768 for full-length 16S against 0.668 for V4'); 526: Population EC F1 by region: mean reverses only in the five regions beginning at V3 or later; prFBA higher for full-length 16S, V1-V2, V1-V3 (0.768, 0.760, 0.766
- `picrust:AS-15` L527: ModelSEED reaction space tracks EC space within 0.01 overall (conditional F1 prFBA 0.805, PICRUSt2 0.721), in every region and both source-in-reference strata; 
- `picrust:AS-16` L530; 538 ('0.235 at NSTI >= 0.5'); 574 (fig A); 1179: NSTI bins: PICRUSt2 precision falls 0.735 -> 0.235 while recall moves 0.870 -> 0.769; prFBA precision 0.822 -> 0.606; prFBA answers 86.0% at NSTI >= 0.5 vs 92.1
- `picrust:AS-17` L531: For 248 of the 603 high-NSTI answered amplicons prFBA's database holds a different organism at >= 0.995 identity; all 603 answers rest on references below 0.987
- `picrust:AS-18` L532; 575 (fig B: 0.32; 7% = 547 of 7,422; 145 further; 692 = 10.8% of 6,407); 1179 ('none at all below its family floor'): Below the 0.865 family floor prFBA abstains on all 547 amplicons whose nearest different organism sits below it, where PICRUSt2 predicts at precision 0.316 (0.3
- `picrust:AS-20` L537 (0.868 recall, 0.735 precision, median NSTI 0.010); 538 (0.642 source absent); 559-563 (tab:picrust2 strata rows): Strata: source in PICRUSt2 ref (n=19,686): PICRUSt2 1.000/0.868/0.735/0.792/0.792, prFBA exclude 0.932/0.871/0.828/0.838/0.781; source absent (n=47,352): PICRUS
- `picrust:AS-21` L539: Include-self prFBA reproduces the truth genome (EC Jaccard 0.94 against truth); PICRUSt2's predictions overlap include-self prFBA's on only Jaccard 0.56-0.61 of
- `picrust:AS-23` L540 (second clause): prFBA's self-excluded selections have a median size of one reference.
- `picrust:AS-24` L572-577 (fig:picrust2, images/prfba/F13_prfba_vs_picrust2.pdf): Figure F13 panels A (NSTI precision decay), B (family-floor bin), C (per-region conditional F1 advantage); caption numbers are inventoried in AS-12, AS-16, AS-1
- `picrust:XR-01` L715 ('the EC coverage of fig:picrust2B ... is 0.982 and 0.829'): sec:depth cites the as-shipped EC coverage by nearest-non-self identity bin: 0.982 (0.945-0.987) and 0.829 (>= 0.995).
- `panels:disc.picrust2_failure_modes` L1179: PICRUSt2 returns a gene complement for every amplicon apart from 26 of the 67,038 scoreable amplicons removed by its input filters; precision decays 0.735 -> 0.
- `panels:disc.prfba_answer_rate` L1179-1181: prFBA answers 89% of self-excluded amplicons and none below its family floor; 11% non-answer (9.6% of EC-scoreable amplicons abstain outright; 1.5% select refer

### R2.5 (11 items)

- `front:abs-depth-ec` L67: Conditional EC F1 0.802 (depth 100) vs 0.798 (depth 20); depth-500 recovered answers EC F1 0.783 vs 0.798, precision 0.747 vs 0.790 for amplicons that already a
- `front:tab:picrust2` L542-568: tab:picrust2 (n=67,038 scoreable): All: PICRUSt2 1.000/0.849/0.669/0.729/0.729; prFBA excl 0.889/0.847/0.787/0.802/0.713; prFBA incl 1.000/0.999/0.936/0.959/0.9
- `val16s:depth_probe_quality` L67 (EC F1 0.783 vs 0.798, precision 0.747 vs 0.790, genus 0.60), 728 (EC recall 0.854/0.747/0.783 vs 0.835/0.790/0.798; 80.2% family, 60.1% genus), 731 (0.601 vs 0.651 genus; 0.802 vs 0.800 family), 1204: Probe-recovered amplicons: EC recall 0.854, precision 0.747, F1 0.783 (vs 0.835/0.790/0.798 for the depth-20 answered set); anchors 80.2% family-correct, 60.1% 
- `val16s:conditional_ec_d20_vs_d100` L67 (0.802 vs 0.798), 732 (F1 0.802 vs 0.798; precision 0.787 vs 0.790; recall 0.847 vs 0.835), 1204: On the as-shipped EC axis the union's conditional F1 is 0.802 at depth 100 vs 0.798 at depth 20 (precision 0.787 vs 0.790, recall 0.847 vs 0.835)
- `val16s:population_ec_and_projections` L733 (0.600 -> 0.713 vs PICRUSt2 0.729), 734-735 (projections: coverage ~0.96, population genus ~0.49, population EC F1 ~0.75): As-shipped population EC F1 of prFBA rises 0.600 -> 0.713 at depth 100, still below PICRUSt2 0.729; projecting the probe's per-region depth-500 rates onto the 1
- `picrust:AS-03` L515: Shared vocabulary: metrics on the 2,132 EC numbers both annotation systems can express and the 7,042 ModelSEED reactions those reach.
- `picrust:AS-07` L518 (2,689 scoreable; 19,686 amplicons); 524 (29% of amplicons); 523; 536; 548; 560 (n=19,686); 563 (n=47,352); 81; 1217: In-reference stratum sizes: 2,689 scoreable in-reference genomes, 19,686 scoreable in-reference amplicons (29% of amplicons), 47,352 source-absent amplicons.
- `picrust:AS-19` L536; 548; 565 (tab:picrust2 ceiling row); 1189 ('90.4%'): Annotation ceiling: for the 2,689 scoreable in-reference genomes PICRUSt2's own annotation recovers 90.4% of the BV-BRC truth ECs at 74.3% precision (Jaccard 0.
- `picrust:XR-02` L733 ('population EC F1 rises from 0.600 to 0.713, still below PICRUSt2's 0.729'); 1187 ('-0.129 at the earlier record depth of 20'): Depth-20 vs depth-100 as-shipped comparison: prFBA population EC F1 0.600 (depth 20) -> 0.713 (depth 100) vs PICRUSt2 0.729; population gap -0.129 at depth 20.
- `panels:lim.depth_probe_quality` L1204: Depth probe recovered calls: 60.1% genus-correct, 80.2% family-correct, EC F1 0.783, at or below the answered set.
- `panels:lim.depth20_vs_100` L1204: Raising record depth 20 -> 100 left conditional figures roughly unchanged: as-shipped conditional EC F1 0.802 vs 0.798; under the legacy reducer conditional gen

### R3 (23 items)

- `front:abs-controlled-genus` L76: Controlled: PICRUSt2 population EC F1 0.779 vs prFBA 0.687 (Delta -0.092 [-0.107,-0.081]); exact-tie 0.685 (Delta -0.094 [-0.109,-0.083]); forced 0.747 (Delta -
- `front:abs-dialect-density` L81: On 19,686 amplicons whose source is in PICRUSt2's reference, scoring against PICRUSt2's own annotation swings the conditional difference from +0.046 to -0.202; 
- `front:abs-region-penalty` L82: As shipped, population EC F1 0.768 full-length vs 0.668 V4; controlled 0.692 vs 0.710 (exact-tie 0.690 vs 0.708), coverage 0.890 vs 0.933.
- `front:tab:sharedref` L681-707: tab:sharedref (shared 7,175-genome reference, n=63,871 cells, matched 57,021; random abstain 50,946, family donor 41,257): coverage/pop F1/cond F1/skill for PIC
- `methods:sharedref_selection_sizes` L219; 303 ("86.8% of answered amplicons on the controlled benchmark"): Controlled self-excluded benchmark at full density, shipped exact-tie selector: fold 0 one genome for 9,846 of 11,751 answered amplicons (83.8%), two for 1,354 
- `methods:sharedref_ngenomes_vs_identity` L228: Across novelty bins on the controlled benchmark, P(n_genomes >= 2) rises monotonically with rank-1 identity.
- `val16s:controlled_forced_recovery` L736-737: Forcing prFBA to answer 99.95% of cells recovers +0.061 of its 0.092 population-F1 deficit (+0.061 of 0.094 exact-tie); residual -0.031 [-0.037, -0.028] (exact-
- `val16s:controlled_selection_sizes` L219, 228, 232 (6.4%; 35.8%), 352 (86.8% of 50,842; 61.7% legacy), 1170 (44,132 of 50,842; 1.199; 61.7%), 1171 (6.4%), 1229, 1236 (4.8% of multi-genome selections; 6.4% of unique tie sets): Controlled benchmark, exact-tie: fold 0 single 9,846 of 11,751 (83.8%), two 1,354 (11.5%), >=3 551 (4.7%); pooled 86.8% (44,132 of 50,842), mean 1.199; 61.7% un
- `picrust:SR-02` L77 ('picked legacy in every fold'); 340 (10/10, 0.001-0.003); 588-590 (legacy in all five folds at both ranks by 0.001-0.003 F1; exact-tie 0.002 lower at identical coverage): Leave-one-fold-out policy selection chose legacy in every fold at both ranks (10/10), by margins of 0.001-0.003 F1; shipped exact-tie scores 0.002 lower in popu
- `picrust:SR-03` L76; 596-597; 1184; 1187 ('-0.092 (-0.094 for the exact-tie rule)'); 1216 ('-0.092, or to -0.094'): Headline, leave-genus-out: PICRUSt2 population EC F1 0.779 vs prFBA 0.687 at coverage 0.893, Delta -0.092 [-0.107, -0.081], beyond the pre-registered 0.02 SESOI
- `picrust:SR-04` L76 ('the deficit holds in all eight regions'); 598: Regions and domains: abstaining -0.059 to -0.182, forced -0.029 to -0.034, every BH q <= 0.01; exact-tie -0.060 to -0.184 and -0.031 to -0.036; domains abstaini
- `picrust:SR-05` L598 (AURC 0.812 vs 0.788, exact-tie 0.787); 599 (cross-knob 0.813 and 0.788, exact-tie 0.786; family swap ~0.007): Risk-coverage: PICRUSt2 leads at every operating point, AURC 0.812 vs 0.788 (exact-tie 0.787; prFBA swept on rank-1 identity, PICRUSt2 on NSTI); exchanging sign
- `picrust:SR-06` L82 ('0.692 vs. 0.710; exact-tie 0.690 vs. 0.708 ... coverage 0.890 against 0.933'); 600: Within prFBA, full-length 16S does not outscore V4: population EC F1 0.692 vs 0.710 (exact-tie 0.690 vs 0.708); coverage 0.890 vs 0.933 under either rule.
- `picrust:SR-07` L601-602: Reaction space: PICRUSt2 population F1 0.780 vs prFBA 0.691, Delta -0.089 [-0.105, -0.078] (exact-tie -0.092 [-0.107, -0.081]); forced -0.028 (exact-tie -0.031)
- `picrust:SR-08` L77 ('copy-the-nearest-reference baseline ... -0.058 [-0.069, -0.051]; exact-tie -0.060'; '0.747 vs. 0.745; exact-tie 0.746'); 605-608; 1184: Reference arms: skill scale from constant (population F1 0.646) to clade-max oracle (0.790); PICRUSt2 covers 0.92, prFBA forced 0.70 (exact-tie forced 0.69), nn
- `picrust:SR-09` L78; 611-616; 1190 ('+0.021 ... +0.019 for the exact-tie rule'); 1192 ('two thirds of the deficit'): Abstention value: declining the 10.7% of cells prFBA declines is worth +0.021 [0.014, 0.030] (q=0.002) over a random mask of the same size; legacy-masked forced
- `picrust:SR-10` L620-621: Novelty bins (nearest surviving reference identity): PICRUSt2 precision declines gently, 0.862 in the 0.945-0.987 band to 0.680 below 0.865; prFBA precision low
- `picrust:SR-11` L683-706 (tab:sharedref, all rows and caption: matched 57,021; random abstain 50,946; family donor 41,257); 613 (0.769 vs 0.767); 1187 (matched -0.027 = 0.769 vs 0.796; exact-tie -0.029, 0.767); 1216 (0.769, or 0.767, against 0.796): tab:sharedref: PICRUSt2 1.000/0.779/0.796/0.92; NSTI-routed 1.000/0.775/0.793/0.90; prFBA 0.893/0.687/0.769/0.28; forced 0.9995/0.747/0.769/0.70; exact-tie 0.89
- `picrust:SR-12` L81 ('+0.107, or +0.106 ... against +0.030 ... fourfold sweep'); 523 ('25, 50 and 100% subsamples ... 3.5 times as steeply (0.107, or 0.106 ... versus 0.030)'); 1189; 1217 ('quartering the shared reference costs prFBA 0.064 ... against PICRUSt2's 0.018'): Reference-density sweep: per tenfold more reference tips prFBA gains +0.107 population EC F1 (exact-tie +0.106) vs PICRUSt2 +0.030 (~3.5x as steep); quartering 
- `picrust:XR-03` L219; 232 (6.4% of multi-genome tie sets span phyla; 35.8% share a genus); 303; 352 (86.8% of 50,842; 61.7% legacy); 1170 (44,132 of 50,842; mean 1.199; 61.7%); 1171 (6.4% of tie sets on the controlled benchmark); 1229: ADJACENT (Methods/Validation/Discussion/Limitations, derived from the controlled benchmark): exact-tie selects one genome for 86.8% of 50,842 answered amplicons
- `panels:disc.selection_size_controlled` L1170, 1229: Shipped exact-tie selector returns a single genome for 86.8% of answered amplicons on the controlled benchmark (44,132 of 50,842; mean selection size 1.199; 61.
- `panels:disc.margin_reversal` L1184, 1187-1188, 1215-1216: As-shipped conditional margin +0.075 on matched amplicons (0.802 vs 0.727) reverses to -0.027 on matched cells (0.769 vs 0.796; -0.029, 0.767 exact-tie); popula
- `panels:disc.vocab_breadth` L1189, 1217: PICRUSt2 reproduces only 90.4% of the BV-BRC truth ECs for genomes it already contains; prFBA's density slope 0.107 (0.106 exact-tie) vs 0.030 population F1 per

### R3F (4 items)

- `front:abs-controlled-family-nn1` L77: Leave-family-out widens the gap to -0.123 [-0.150,-0.099] (exact-tie -0.126 [-0.152,-0.102]); nested CV picked legacy in every fold; abstaining arm vs nn1 -0.05
- `front:abs-controlled-random-mask` L78: Refusals better targeted than a random mask: +0.021 population F1 (+0.018 [-0.004,0.040] at family rank; exact-tie +0.019 and +0.015); answering everywhere reco
- `picrust:FR-02` L77 ('-0.123 ([-0.150, -0.099]) ... -0.126 ([-0.152, -0.102])'); 624: Leave-family-out headline: PICRUSt2 0.749 vs prFBA 0.626 at coverage 0.856, Delta -0.123 [-0.150, -0.099] (exact-tie 0.623, -0.126 [-0.152, -0.102]); forced prF
- `picrust:FR-03` L78 ('+0.018 ([-0.004, 0.040]) at family rank ... +0.015 for the exact-tie rule'); 599 (swap ~0.007); 625-628; 1192: Leave-family-out secondary: every region further against prFBA (-0.076 to -0.245, all q=0.002; exact-tie -0.077 to -0.247); forced-arm equivalence with nn1 pers

### R4P (10 items)

- `front:abs-phenotype-readout` L79: Phenotype readout: no mean difference (PICRUSt2 ahead on 2 of 22 media); n=755 ATP-informative genomes (747 V4), EC F1 differ by -0.028; prFBA drafts mismatch t
- `picrust:PH-01` L79 ('n=755 ... (747 of them V4) ... differ by only -0.028 ... 0.13 more media ... p=0.49'); 632 (1,000-genome subsample, V4 where available, exact-tie; -0.094 vs -0.092); 634; 1185; 1193 (755; -0.028); 1222 (755): Primary ATP readout: over 755 paired ATP-informative genomes (747 V4) of a 1,000-genome subsample, prFBA's exact-tie arm mismatches the truth ATP profile on 0.1
- `picrust:PH-02` L79 ('two of 22 media'); 631 (Ac.NO and Etho.O2, q=0.002); 639; 1185: Per-medium McNemar favours PICRUSt2 on two of the 22 testable media (Ac.NO and Etho.O2, q=0.002).
- `picrust:PH-03` L634 (correlation -0.14); 1190 (Spearman rho = -0.14, p = 2x10^-4): Per-genome correlation between the methods' EC-F1 difference and phenotype difference is only -0.14 (Spearman, p = 2e-4).
- `picrust:PH-04` L635-637; 1186: Against nn1 the phenotype favours prFBA: 0.146 fewer mismatched media (n=755, p=0.017, q=0.028; p=0.28 over all paired genomes); on those genomes the two tie on
- `picrust:PH-05` L79 ('only 22 of 54 media carry ATP signal'); 638 (53 media + empty control; 22 carry signal); 1222: The ATP vector spans 53 media plus an empty-medium control; only 22 ever carry signal for the truth or compared arms under this template route.
- `panels:disc.atp_vs_picrust2` L1185: ATP-production profiles not detectably further from truth: prFBA mismatches truth on 0.13 more media than PICRUSt2 (Wilcoxon p=0.49, q=0.49, n=755); per-medium 
- `panels:disc.atp_vs_nn1` L1186: Against nn1, prFBA's profiles sit closer to truth over 755 ATP-informative genomes (-0.146 media, p=0.017, q=0.028; p=0.28 over all paired); on those genomes EC
- `panels:disc.refusal_credit` L1190-1193: Choosing which amplicons to decline is worth +0.021 population F1 over declining as many at random (+0.019 exact-tie); EC-F1 vs ATP-distance difference Spearman
- `panels:lim.phenotype_scope` L1222: Phenotype: null at 755 paired genomes; ATP vector carries signal in only 22 of 54 media; gapfilled growth/essentiality on a 200-genome gapfilled subsample of th

### R4S2 (4 items)

- `front:abs-phenotype-controls` L80: Pre-registered C12 failed (5 of 200 gapfilled truth models gated); after A-C12 and A-C13 the 1,000-genome ATP pass reproduced every ATP vector (43,140 of 43,140
- `picrust:PH-06` L641-646; 1222 (n=200): Structure/gapfilled readouts (A-C13; unpaired per-arm medians): draft reaction Jaccard PICRUSt2 0.737 vs prFBA 0.707 (nn1 0.705) over 1,000 genomes; essentialit
- `picrust:AM-02` L80 ('43,140 of 43,140 rows ... 0 of 8,614 rows, against 49 ... valid under the amendments'); 667-673; 676-679 (seven-reaction draft correction: 865 -> 95 of 1,181 drafts; medians move <= 0.004; burden <= 1 reaction, 1.5 as shipped; 138 of 200 truth drafts make ATP from nothing; 43,140 ATP-pass rows with no empty-medium ATP); 1193 (ATP pass reproduces every vector): Properties of the amended runs: the 1,000-genome ATP pass reproduces every ATP vector (43,140 of 43,140 rows) under A-C12 and A-C13; no gapfilled model of any a
- `picrust:AM-03` L648-652: Controls: nine fold-scoring controls (C3, C4, C6, C7, C8, C10, C13, C14, C15) evaluated with no hard failure; C2 self-inclusion ceiling run separately passes: p

### R4H (11 items)

- `front:abs-objective-posthoc` L61: prFBA objective (post-hoc in-silico only): as released moves little flux and breaks no more ties than shuffled weights; with probabilities from each query's 100
- `methods:wobj_results` L307-311: p_func lowered truth-absent flux at the face midpoint by 0.0061 [0.0025, 0.0097] (4.7% of two-tier 0.130; below SESOI 0.01); face minimum -0.0025 [-0.0053, +0.0
- `methods:msprob_released_equivalence` L285, 313-314: Released single-model path fails; corrected copy and released source with two lines fixed gave identical objective values on 360 of 360; fix commits b6857a0, 09
- `methods:stageO_design_and_split` L317-318: Four factors: probability source, cost shape, gapfilled probability (0, 0.084 = measured support rate, 0.2), exchange weight (100, 10, 1); 90 objectives per amp
- `methods:stageO_results` L319-324: Selected objective (estimate, linear, G0, X1) lowered truth-absent flux vs released by 0.0172 [0.0118, 0.0221]; face minimum 0.0076 [0.0047, 0.0133]; beat shuff
- `methods:stageO_wrapper_equivalence` L325: The shipped recommended call (exchange weight 1, gapfilled p 0, biomass floored) reproduces the tuned objective bitwise on eight amplicons spanning every select
- `val16s:pfunc_calibration` L225 (modal fractional weight 0.5), 353, 1172, 1225 (slope 0.78; informative only for >=3 refs): p_func ranks functions correctly but is overconfident: probability-scale calibration slope 0.78 [0.76, 0.79]; skill over a constant only for selections of >=3 r
- `val16s:stage_e_estimator` L354-360, 1172 (calibration), 1225 ('tuned on a benchmark reference of about 7,000 genomes'): Neighbour-and-clade estimate (E3: top-100 hits, tau 0.02, family prior, pseudo-count 2; no prior for 13% of cells; same settings in all 5 folds) over 57,021 cel
- `panels:disc.objective_test1` L1172, 1225: First post-hoc objective test (300 multi-genome + 60 single-genome amplicons): weighting cuts gene-backed reactions free to vary on the optimal face by 21%, shu
- `panels:disc.objective_test2` L1172, 1225: Second test: configuration chosen on 121 of the 201 holdout groups; on the other 80, truth-absent flux fell by 0.0172 [0.0118, 0.0221] vs the released objective
- `panels:lim.reconstruction_tests` L1225: Eq. objective evaluated only on the same 300 multi-genome and 60 single-genome amplicons; second test's configuration chosen on 121 of 201 holdout groups; Eq. f

### R5A (9 items)

- `methods:sharedref_tieset_coherence` L232 (first clause): Among multi-genome tie sets of the controlled benchmark, 6.4% unite genomes from different phyla and only 35.8% share a genus.
- `methods:recon_modal_fractional_weight` L225: Where the weighting is active the modal fractional weight is 0.5.
- `picrust:AS-10` L74; 521; 1187 ('+0.075 on matched amplicons (0.802 against 0.727)'); 1216 ('matched 0.802 against 0.727 as shipped'): Matched both-answer comparison: on the 59,607 amplicons both methods answer, conditional EC F1 0.802 (prFBA) vs 0.727 (PICRUSt2), +0.075.
- `mag:magval_identity_bins` L874: Binned by the amplicon arm's best non-self identity, at 16S identity 0.945-0.987 the amplicon answers 98% of the time while whole-genome ANI answers 8-20%.
- `mag:magval_scope_claims` L882-885: The MAG advantage is 0.08-0.21 F1 at any completeness >= 1%; the conditional comparison covers 31% of organisms; on the other 69% the amplicon answers at F1 0.6
- `mag:magval_tau99_paired` L892-895: At tau=99 the paired advantage over FL widens +0.103 -> +0.127 and population F1 rises 0.248 -> 0.301 vs the amplicon's 0.605; at matched stringency FL reaches 
- `panels:panels.welch_tests` L1020: Original vs k-means conditional F1 difference significant in all eight regions (Welch p<0.01; p<1e-8 for V4, V4-V5, V7-V9), full-length +0.011 (p=2.3e-4), and a
- `panels:panels.mag_singleton_abstention` L1083: MAG arm, organisms that are the only member of their species in the panel: abstention 0.6895 on k-means (n=3,942) vs 0.6883 on the original (n=9,580), differenc
- `panels:panels.mag_multiplicity` L1084-1085: Higher multiplicities only on the k-means panel: abstention 0.559 for species with >=20 panel members; that stratum accounts for k-means' lower overall MAG abst

### D1 (6 items)

- `front:abs-sharedref-size` L75: Controlled design: species-dereplicated shared reference of 6,788-7,175 genomes per fold (median 7,175).
- `picrust:SR-01` L75 (6,788-7,175; median 7,175); 584; 587 (8,872 eligible; 63,871 cells); 588 (2,376 groups, n_eff 6,716, ICC 0.33); 683 (7,175; n=63,871); 1215; 1219 (7,175, 'forty-fourth'); 1220 (2,376; 6,716): Controlled design sizes: per-fold shared reference 6,788-7,175 genomes (median 7,175); 8,872 eligible genomes; 63,871 scored cells; cluster bootstrap over 2,376
- `picrust:FR-01` L624 (2,059 held-out genomes; 14,692 cells; median 8,468-genome reference); 627 (171 groups vs 2,376; n_eff 598 vs 6,716); 1220: Leave-family-out design sizes: 2,059 held-out genomes, 14,692 cells, median 8,468-genome reference, 171 holdout groups, effective n 598.
- `picrust:FR-04` L627 (18.6%, 1,653 of 8,872); 1221: Family blocking imperfect for 18.6% (1,653 of 8,872) of eligible panel genomes whose family label is an unclassified placeholder.
- `panels:lim.database_sizes` L1214-1215, 1219: Databases differ roughly tenfold (~315,000 BV-BRC genomes behind 459,301 unique >=200 bp 16S vs 27,870 PICRUSt2 reference genomes); 7,175-genome shared referenc
- `panels:lim.holdout_power` L1220-1221: Leave-family-out resolves into 171 holdout groups vs leave-genus-out's 2,376 (effective n 598 vs 6,716); family blocking imperfect for 18.6% of the 8,872 eligib

### D7 (1 items)

- `val16s:backend_concordance` L384-386: Identical representative genome set for 99.94-99.97% of ASVs both runs built (3,724 EmilyKin, 3,219 bioreactor); divergences are equal-score ties; the two missi

### NO (124 items)

- `front:abs-16S-database` L61: BV-BRC 16S database: 1,428,909 raw records; 460,117 MD5-deduplicated unique sequences; 459,301 retained at >=200 bp.
- `front:abs-embedding-abandoned` L62: Embedding retrieval abandoned: on a 97,624-insert V4-V5 store cosine #1 was the highest-identity hit among cosine's top 20 for only 54.7% of ASVs; cosine top-10
- `front:abs-gpu-aligner` L63: 57.3 GCUPS exhaustive scoring; 7.8x over a 60-process Biopython CPU baseline on a 32-core/64-thread host; bit-exact on 400/400 control pairs and on all 6.0 mill
- `front:abs-panel-size` L63: Deterministic 10,000-organism, eight-region, 67,991-amplicon in-silico validation.
- `front:abs-selfincl-genus` L64: Source in database: genus accuracy 0.89-1.00.
- `front:abs-selfexcl-genus` L65: Self-exclusion (0.987): 9.5% of amplicons abstain; conditional genus accuracy 0.58-0.73; whole-population genus accuracy 0.37-0.54.
- `front:abs-tau0995-taxonomy` L66: Excluding only hits >=0.995: abstention 6.0%; conditional genus accuracy 0.63-0.75; whole-population genus accuracy 0.44-0.59.
- `front:abs-depth-selection` L67: 86.7% of abstentions arise because exclusion empties a record truncated to 100 hits; raising depth 20->100 cut abstention 23.4%->9.5%; re-aligning at depth 500 
- `front:abs-threepanel-taxonomy-quality` L69: Conditional family accuracy differs by at most 0.052; conditional genus accuracy higher on the k-means panel in every region (up to 0.07) and on the random pane
- `front:abs-threepanel-genus-coverage` L70: Self-exclusion genus coverage on the random panel falls by 0.20 to 0.26 in every region (part of 'genus and gene-capture coverage falls by 0.20 to 0.29').
- `front:abs-gtdb-reversal` L72: Searched as synthetic MAGs against GTDB (species exclusion empties a 100-hit record in <1% of cases), random-panel organisms resolve more often than the stratif
- `front:abs-family-vs-genus-accuracy` L73: Family-level assignment 0.76-0.85 vs genus 0.58-0.73 conditional under self-exclusion; short amplicons lower family-rank accuracy, by less than genus.
- `front:intro-no-numbers` L90-123: No quantitative claims; the only numeric tokens are the ModelSEEDpy commit c2e31a0 (released prFBA function fails on both paths) and branch fix-prfba-single-mod
- `front:tab:main-taxonomy` L420-442: tab:main genus accuracy incl. (0.891-0.996), excl. cond/pop (0.576/0.371 .. 0.726/0.538) and family excl. cond (0.759-0.848) per region; caption: anchor predict
- `front:tab:tau-taxonomy` L444-467: tab:tau genus cond, genus coverage, genus pop, family cond and abstained at 0.987 / 0.995 per region (e.g. V4 0.576/0.628, abstained 11.8%/8.0%). [Hit-removal c
- `front:fig:accuracy` L469-475: F1_accuracy_heatmap.pdf: per-region, per-rank concordance (anchor/consensus) and anchor coverage under self-inclusion; caption: anchor genus 0.89-1.00, species 
- `front:fig:leakage-genus` L485-491: F11_self_inclusion_exclusion_delta.pdf: genus leakage gap 0.25-0.27 (FL, V1-V3, V1-V2), 0.31-0.32 other five regions.
- `front:fig:perrank` L493-499: F2_perrank_concordance.pdf (self-exclusion): family-level assignment 0.76-0.85.
- `front:fig:embedding` L501-506: fig_pairwise_hist.png: random 16S pairs average cosine ~0.90.
- `front:tab:mag` L771-806: tab:mag: 291 real EBPR MAGs vs GTDB r226 (143,614 species reps) and a 2,160-genome BV-BRC subset: novelty classes (105/79/14/62/31; 3/7/5/98/178), AF-qualified 
- `front:tab:panelcomp` L964-988: tab:panelcomp: Archaea 1,289 (12.9%)/494 (4.9%)/392 (3.9%); phyla 206/142/130; family bins 1,412/861/753 (named 883/564/478); genus bins 3,364/1,996/1,699 (name
- `front:fig:panelspace` L996-1007: F16_panel_space.pdf: 210,779 clusterable 16S sequences (140,197 organisms); density = cosine distance to 21st neighbour in 1024-d (210,779/10,000 ~ 21); Spearma
- `front:tab:threepanel-genus` L1101-1137: tab:threepanel 16S genus block: phylum-standardized genus coverage (FL 0.700/0.544/0.479 ... V7-V9 0.640/0.457/0.382), Delta coverage (-0.220 .. -0.257 with CIs
- `front:tab:threepanel-gtdb` L1101-1137: tab:threepanel GTDB block: resolve rate 0.577/0.675/0.715 (100%), 0.579/0.677/0.717 (50%), 0.568/0.671/0.713 (30%); Delta +0.139/+0.137/+0.145 with CIs.
- `front:tab:threepanel-caption` L1101-1110: tab:threepanel caption: 113 phyla common to all three panels; GTDB block over 68 phyla present in all three 20% subsamples, which hold 93.9%, 98.5%, 99.7% of th
- `methods:refdb_counts` L134-136: Raw BV-BRC 16S pull 1,428,909 header->sequence records; MD5 dedup to 460,117 unique; >=200 bp filter leaves 459,301 unique references.
- `methods:align_prefilter_fraction` L141; 146 ("<0.11%"): edlib prefilter eliminates ~99.9% of candidates; the 500-candidate shortlist scores <0.11% of the database.
- `methods:align_tie_cluster_width` L146: On the 10,000-organism panel 208 of 60,461 distinct amplicon sequences (0.34%) have >=200 shortlisted references tied at the top score and 95 (0.16%) have all 5
- `methods:align_scheme_and_kernel_constants` L143-144, 148, 152, 155-158: Top-500 edlib k-band heap; SW scheme +2/-3/-5/-2; Stage 3 aligns top 100 (--prefilter-k 500 --k2 100), fallback re-score over --gpu-cand 200; sm_120 RTX 5070 Ti
- `methods:align_gpu_exactness` L152 (no ASV needed the fallback); 159; 333 (60,461 records byte for byte, no fallback): GPU score bit-exact to Biopython on all 6,046,100 hits of the depth-100 records; ranking by GPU score reproduces all 60,461 records byte for byte with 0 fallbac
- `methods:hardware_stage3_timings` L333: Stage 3 shortens from 7,128 s to 3,425 s on 40 workers; re-score ~4.3 ms per candidate, traceback ~16 ms; edlib prefilter took 70% of depth-100 wall time, ~four
- `methods:mag_gtdb_database` L168: GTDB r226 skani database of 143,614 species-representative genomes; lineages joined from 732,475 bac120+ar53 accessions.
- `methods:mag_af_gate_and_ani_knobs` L171-173, 176: A spurious reference displaced the best hit for roughly two dozen MAGs in development; AF gate 15%; novelty bins >=95 / 85-95 / <85, no_hit below ~80 ANI; ANI k
- `methods:mag_union_size_one` L178-179: At genome scale the union was of size one for every MAG resolved in sec:magreal and for 99.50% of the 49,985 synthetic MAGs of sec:magval, against the amplicon 
- `methods:select_thresholds` L183-184: Yarza/Stackebrandt bands t_species 0.987, t_genus 0.945, t_family 0.865; reliability gate coverage <0.90 or |SW - edlib identity| > 0.02.
- `methods:select_union_mean_size` L187 ("mean 2.88 genomes"): With the source present the 100-hit exact-tie union has mean 2.88 genomes on the 10,000-organism panel.
- `methods:recon_bvbrc_row_limit` L196: Every BV-BRC query pages past the 25,000-row response limit.
- `methods:recon_literature_calibration` L212: Species-boundary genomes share 65-90% of genes; core ~80% of a genome; single conspecific reference supports p ~0.7-0.9; up to 10-35% of true content missing.
- `methods:deployment_selection_sizes` L220; 253 ("1,410 multi-genome synthetic genomes of the 3,950-ASV deployment"); 303 ("62.5% of answered ASVs"): Real-amplicon deployment (EmilyKin, depth 100): 2,345 of 3,755 answered ASVs at n=1 (62.5%), 597 at two, 813 at three+, mean 2.289, 195 abstaining; 3,950 ASVs; 
- `methods:refdb_genome_count` L220, 264 ("~315,000 BV-BRC genomes"); 265 ("roughly ten times"): The 459,301 references come from ~315,000 BV-BRC genomes; roughly ten times PICRUSt2's 27,870.
- `methods:e3_release_design_constants` L251, 253: E3: identity-weighted share over the first 100 hits passing gate and family floor, one per genome (production: one per species), shrunk to anchor-family prevale
- `methods:prior_art_reference_sizes` L265: PICRUSt2 2.6.3 default references hold 27,870 genomes; PanFP 2,400.
- `methods:prfba_objective_constants` L269, 277, 285, 286, 294, 299, 313: ele_limit 100; kinCoef 1000; biomass floor >=95% of max; exchange weight 100; min_prob 0.01; released function at ModelSEEDpy c2e31a0.
- `methods:wobj_design` L305-306, 312: Declared post-hoc test sampled 360 amplicons of the leave-genus-out benchmark (100 each with two, three and four+ genomes, 60 with one), five weight sets on the
- `val16s:panel_design` L63, 388-394, 1197: 10,000 deterministic taxonomy-stratified BV-BRC genomes; 8 regions (V1-V2, V1-V3, V3-V4, V4, V4-V5, V4-V5_944R, V7-V9, FL); 67,991 amplicons; Caporaso/Parada-Ap
- `val16s:region_amplification_archaea` L414, 1198: 292-bp V4 amplicon (0.987 allows at most 3 differences); FL/V1-V2/V1-V3/V3-V4 amplify 35-76 archaeal amplicons each (<=6% extraction); V4, V4-V5, V4-V5_944R 67-
- `val16s:unnamed_genus_denominator` L399: 22-25% of source organisms have no named NCBI genus and are excluded from both denominators
- `val16s:tabmain_taxonomy` L64-65, 73, 401-403, 408, 432-439 (genus incl., genus excl. cond./pop., family excl. cond. columns), 473, 497, 1175 (family 0.76-0.80 cond., 54-59% pop. family, 37-42% / 0.51-0.54 pop. genus, 'roughly a third' no genus call): Self-inclusion anchor genus 0.89-1.00; self-exclusion conditional genus 0.58-0.73, population genus 0.37-0.54; conditional family 0.76-0.85; per-region tab:main
- `val16s:taxonomy_derived_gaps` L404 (genus gap), 408 (11-18 points), 410 (0.16-0.18), 473 (species 0.65-0.88, FL 0.99), 489 (fig:leakage caption genus gap): Genus leakage gap 0.25-0.27 (FL, V1-V3, V1-V2) vs 0.31-0.32 (five V3-onward regions); family 11-18 points above genus; family->genus gap 0.16-0.18 for V3-onward
- `val16s:taxonomy_denominators` L65, 397 (taxonomic part), 1175 ('no genus call for roughly a third'): 9.5% of all amplicons receive no taxonomic call under self-exclusion (10.3% of genus-scoreable); a further 22% of genus-scoreable amplicons resolve to a referen
- `val16s:tau0995_taxonomy` L66, 394, 415-416, 456-463 (tab:tau genus cond., coverage, pop., family cond., abstained columns): At 0.995 the exclusion removes 881,501 hits instead of 1,267,002 and leaves 4,101 amplicons without a call instead of 6,482; abstention 9.5% -> 6.0%; conditiona
- `val16s:mcnemar` L396, 405, 1197: Paired McNemar, Holm across 28 pairs: 26 of 28 differ at genus under self-exclusion; n=5,533-7,238 per pair, 5,481 common to all eight; n.s. V1-V2/V1-V3 (Holm p
- `val16s:tau_rationale_deployment_bands` L418: In the 3,950-ASV deployment 47% of ASVs have no hit >=0.987, 19% best hit in 0.987-0.995, 34% >=0.995 (52/14/34 for the 3,276-ASV bioreactor deployment)
- `val16s:union_selection_sizes_10k` L178 (34.0%), 187 (mean 2.88), 352 (68.1% of 61,509; 61,534 / 68.1%): Amplicon exact-tie union holds >1 genome for 34.0% of self-inclusion selections, mean 2.88 genomes; under self-exclusion one genome for 68.1% of the 61,509 answ
- `val16s:abstention_decomposition_d100` L67 (86.7%; 23.4% -> 9.5%), 398 (13.3%, 6,457), 415 (6,457 vs 6,482, 25 unexplained), 712-713, 738, 1203, 1205: Union self-exclusion abstains on 6,457 (9.5%) at depth 100 vs 15,888 (23.4%) at depth 20; 5,598 (86.7%, vs 94.6%) are emptied records; of the other 859 (13.3%),
- `val16s:coverage_by_identity_bin` L715 (0.998 vs 0.843; 30.8 of 100 removed; 0.604, 12.1 of 20 at depth 20), 716 (include-self 1.000 every bin): Exclude-self selection coverage 0.998 for amplicons whose nearest other organism is at 0.945-0.987 vs 0.843 at >=0.995 (mean 30.8 of 100 hits removed); 0.604 wi
- `val16s:depth_probe_counts` L67 (88.4%, 283/320), 718-727, 728 (281 of 283 genus-tier; median 0.9856): 320 emptied-at-depth-20 amplicons (40/region, 316 unique) re-aligned at depth 500: 283 (88.4%) resolve; 37 have all 500 hits >=tau (<=6 source copies; 170-495 o
- `val16s:d100_realized_recovery` L727 (9,424 of 15,026 = 62.7%; 5,598 still emptied), 730 (3,601 of 60,461 = 6.0%), 733 (coverage 0.766 -> 0.905; probe-implied 0.912): Depth 100 answers 9,424 of the 15,026 depth-20-emptied amplicons (62.7%), 5,598 still emptied; for 3,601 of 60,461 distinct sequences (6.0%) the 100th hit is st
- `val16s:legacy_tabmain_genus_d20_vs_d100` L732 (conditional genus <0.01), 733 (population genus 0.375 -> 0.472), 734 (base 0.375), 1204 (genus within 0.01): Under the legacy reducer conditional genus accuracy moves by <0.01 in every region from depth 20 to 100 (-0.0068..+0.0091); pooled population genus accuracy ris
- `val16s:alignment_cost_d100` L739-740: edlib prefilter takes 17,470 of the depth-100 run's 24,949 s (70%); Biopython stage 1,147 s on 60 workers (depth 20) -> 7,128 s on 40 (4.1x worker-seconds; ~3/5
- `val16s:union_self_recovery` L1201-1202: Residual union self-recovery gap of the seven sub-full-length regions: 1,145 of 60,693 amplicons; 525 (46%) source not in the top-100 (byte-identical to >=100 o
- `val16s:all_tied_records` L146 (208 of 60,461 = 0.34% with >=200 shortlisted refs tied at top score; 95 = 0.16% all 500 tied), 1206 (453 of 60,461 = 0.75%; 2.2% of 67,991; median 19 species, 22 per record): Shortlist and record tie statistics on the 10k panel
- `val16s:gpu_bitexact_and_benchmarks` L63 (6.0 million hits; 57.3 GCUPS; 7.8x; 400/400), 159, 365-367: GPU SW bit-exact to Biopython on 400/400 pairs and on all 6,046,100 depth-100 hits; exhaustive 20 x 459,301 sweep 65.9 s at 57.3 GCUPS, 7.8x over the 511 s 60-p
- `val16s:deployment_selection_sizes` L220, 1170 (2,345 of 3,755, mean 2.289), 1208 (62.5%; 2.034 -> 2.289; largest 19 -> 77), 1229 (62.5%): EmilyKin depth-100 deployment (3,950 ASVs vs 459,301 references, ~315,000 genomes): 2,345 of 3,755 answered at n=1 (62.5%), 597 at n=2, 813 at >=3, mean 2.289, 
- `val16s:deployment_truncation` L1207: In 25 of the 3,950 deployment ASVs the 100th retained hit ties the best on identity and score; those selections (7 to 77 genomes) are lower bounds
- `val16s:edlib_prefilter_3950` L336 (99.47%), 370-373: edlib top-500 prefilter recovers the exhaustive top-1 for 99.47% of the 3,950 ASVs, mean top-20 overlap 97.6%; 219 ASVs (5.54%) lost >=1 top-20 hit scoring abov
- `picrust:AS-01` L510; 511; 265; 582; 1214; 1219: Benchmark size and PICRUSt2 reference: 60,461 MD5-unique amplicons underlying 67,991 amplicon-organism cells; PICRUSt2 2.6.3 default bacterial+archaeal referenc
- `picrust:AS-02` L511; 582; 1214; 1219: prFBA database size: ~315,000 BV-BRC genomes behind 459,301 unique >=200 bp 16S sequences, roughly tenfold PICRUSt2's 27,870; the 7,175-genome shared reference 
- `picrust:AS-06` L518 (2,720 of 9,746); 1213: Source assembly is in PICRUSt2's (GTDB-based, species-dereplicated) reference for 2,720 of 9,746 benchmark source genomes.
- `picrust:AS-22` L540 (first clause): PICRUSt2 routes 10.2% of archaeal amplicons (467 of 4,587) to its bacterial tree.
- `picrust:AM-01` L80 ('five of 200 gapfilled truth models gated out where zero are allowed'); 653-664; 666 (166 of 166 pinned rows); 668 (180 genomes outside the validation); 669 ('136, in 50 genomes, did not before'); 674-675 (four models; 72 of 72; commit e1254ba); 1222 (five truth models, 10 gated rows): Record of the pre-registered failure and the fix's validation: 5 of 200 truth models gated (10 rows); fraction leg 49 gated rows of 8,614 vs a 10% ceiling; exit
- `picrust:XR-04` L1187 (last clause: re-aligns each fold, retains 200 hits per query, --k2 200); 589 (5-fold grouping); 588 (SESOI 0.02 pre-registered): Design constants: 5-fold grouping, pre-registered SESOI 0.02 raw DeltaF1, prFBA fold runs retain 200 hits (--k2 200), tau 0.987 leakage edges.
- `mag:methods_ani_db_counts` L168: GTDB r226 skani sketch of 143,614 species-representative genomes; lineages joined from 732,475 accessions of the bac120+ar53 taxonomies.
- `mag:methods_af_guard_anecdote` L171-173: A single spurious high-ANI/near-zero-AF reference displaced the correct best hit for roughly two dozen MAGs in development; AF_query gate 15%; novelty classes >
- `mag:methods_ani_thresholds` L175-176: Recalibrated selector thresholds: t_species 0.95, t_genus 0.85, abstention floor 0.83, coverage gate 0.15 (replacing 0.987 / 0.865 / 0.90).
- `mag:methods_union_size_mag` L178-179; 1173-1174 (Discussion, qualitative restatement): At genome scale the union is of size one for every resolved real MAG and for 99.50% of the 49,985 synthetic MAGs of sec:magval; exact ANI ties occur nearly two 
- `mag:magreal_assessment_set` L751-753: 291 dereplicated EBPR MAGs (5 hybrid metagenomes, MetaBAT2, dRep 95% ANI); CheckM2 mean completeness 87.9%, mean contamination 1.63%, all MIMAG medium or better
- `mag:magreal_gtdb_headline` L754-755; tab:mag 771-806 (GTDB r226 column: novelty classes 105/79/14/62/31, 198 (68%), 95.66/72.5%, 119 s): Against GTDB r226, 198 of 291 MAGs (68%) obtain an AF-qualified best reference at median ANI 95.66 over median 72.5% of the query; search 119 s on 60 threads; 1
- `mag:magreal_concordance` L756-757; tab:mag concordance block (GTDB column): Concordance with GTDB-Tk: domain 198/198, phylum 178/198, class 186/198, order 196/198, family 197/198, genus 176/197, species 87/97 (90%); identical accession 
- `mag:magreal_rename_attribution` L759-763: All 20 phylum disagreements are renames (15 Bacteroidota->Bacteroidota_A, 5 Patescibacteria->Patescibacteriota); all 12 class disagreements renames (9 Verrucomi
- `mag:magreal_bvbrc_arm` L765-768; tab:mag BV-BRC column (3/7/5/98/178, 15 (5%), 87.43/34.0%, 6.3 s, 15/15 15/15 14/15 14/15 10/15 15/15 3/8): BV-BRC 2,160-genome subset: 15 MAGs reach an AF-qualified best hit, 178 no hit, median qualified AF 34% (ANI 87.43); 14 of the 15 clear the abstention floor and
- `mag:magreal_selection` L769-770: Selector on both records: a reference genome for the 198 GTDB-qualified MAGs and abstention for 93; abstention causes low_AF_unresolved 62, nothing found 31.
- `mag:magval_construction` L812-820: Contigs retrieved for 9,997 of 10,000 panel genomes (3 hold no sequence); ~50 kb mean synthetic contigs; 19.8% of panel genomes closed, 40.4% <10 contigs, 1 to 
- `mag:magval_afref_completeness` L829; fig:magcomp caption 921-929 panel A; 1210 (Limitations): AF_ref ~ 0.99 x completeness (an essentially unbiased completeness estimator the pipeline discards); Limitations: about 0.99x when the source is in the referenc
- `mag:magval_marker_reach` L871-873: Among the 6,011 panel organisms with a same-genus congener in the reference, the best surviving hit aligns over a median 20.3% of the query and 45.1% fall below
- `mag:magval_gtdb_density` L875-876: A ~14-fold denser reference (143,614 GTDB r226 species reps) nearly doubles the resolve rate to 0.56 with completeness-invariance intact; 0.556 at full complete
- `mag:magval_universe_asymmetry` L885-887: The arms search reference universes differing ~thirtyfold (~315,000 BV-BRC genomes behind 459,301 unique 16S sequences vs 9,997 panel genomes); the 95-ANI rule 
- `mag:magval_tau99_mag_only` L888-890: Re-scoring at tau=99 retains 1,285 references in the 95-99 band; genus accuracy slips 0.910 -> 0.895.
- `mag:twopanel_mag_depth_probe` L1014-1016: [sec:twopanel, cross-reference] MAG depth probe (447 capped + 200 control queries at -n 2000: 11 selections changed, all abstention -> resolution; depth 200 con
- `mag:twopanel_mag_gtdb` L72 (Abstract); 1127-1130 (tab:threepanel GTDB rows); 1151-1158: [sec:twopanel, cross-reference] GTDB arm on ~20% phylum-proportional subsamples (2,000/1,995/1,989): species exclusion empties a depth-capped record in under 1%
- `panels:panels.design_scope` L941-945: Original panel = 10,000 organisms, taxonomy-stratified, one per species; replication on two further 10,000-organism panels: 16S all eight regions, self-inclusio
- `panels:panels.kmeans_construction` L947-951: k-means++ k=10,000, seed 1729, 15 iterations on L2-normalized 1024-d NT embeddings of 210,779 full-length (>=1400 bp, non-atypical) 16S of 140,197 organisms; co
- `panels:panels.random_construction` L953-957: Uniform-random panel: 10,000 organisms without replacement, seed 1729, from the same 140,197 admissible organisms; representative 16S = longest full-length non-
- `panels:panels.tab_panelcomp` L959-960, 966-987: tab:panelcomp and prose: Archaea 1,289 / 494 / 392; phyla 206 / 142 / 130; family bins 1,412 / 861 / 753 (named 883 / 564 / 478); genus bins 3,364 / 1,996 / 1,6
- `panels:panels.archaea_shares` L961: Archaea 12.9% (stratified), 4.9% (k-means), 3.9% (random), against 3.5% (4,945 of 140,197) of the admissible pool.
- `panels:panels.overlaps` L962: Overlap 804 original-k-means, 524 original-random, 701 k-means-random, 49 in all three; 7,292 of the original's organisms belong to the admissible pool.
- `panels:panels.space_deciles` L989-993: Two sparsest deciles hold 19.2% of sequences and supply 49.1% of the original panel's on-frame picks vs 36.6% of k-means (excess ~1.7x); densest decile holds 12
- `panels:panels.offframe_counts` L994, 1004: 903 of the original panel's organisms are off frame: 743 (603 shorter than 1,400 bp, 140 atypical) and 160 lineage-gate failures (140 index rows without taxon i
- `panels:panels.fig16_caption` L996-1007: F16: 210,779 clusterable sequences; density = cosine distance to 21st NN (k=210,779/10,000~21); ranking stable in k (Spearman 0.964, k=21 vs 50); deciles hold 7
- `panels:panels.mag_random_sketch` L1009-1010: Random panel's panel-as-sketch MAG arm mirrors k-means: same synthesizer and seed, completeness {100,90,70,50,30}%, 500-genome batches; sketch reproduces k-mean
- `panels:panels.mag_depth_probe` L1014: MAG depth probe: depth 100 vs 2,000 on 447 capped + 200 control queries; 11 selections changed, all abstention -> resolution; depth 200 already converged.
- `panels:panels.mag_capped_counts` L1015-1016: Only 0.16% of MAG queries (156 of 99,970 records incl. contaminated grid; 93 of the 49,985 uncontaminated 30-100% records, 0.19%) were capped and emptied at dep
- `panels:panels.ladder_meta` L1089-1092: Common-support rung leaves 9,265 / 9,963 / 9,987 organisms (full-length best 16S); pooled mix of the 113 phyla present in all three panels; 1,000 phylum-stratif
- `panels:panels.v4_genus_coverage` L1094: Standardized V4 genus coverage 0.615 / 0.431 / 0.358 (original / k-means / random); random-minus-original -0.257 vs -0.330 crude.
- `panels:panels.cond_accuracy_taxonomy` L1097: Conditional genus accuracy under self-exclusion differs random vs original by -0.016 to +0.062 (higher on random in six of eight regions); k-means higher than o
- `panels:panels.tab_threepanel_genus_block` L1101-1119: tab:threepanel genus block, phylum-std genus coverage original/k-means/random and Delta random-original coverage | conditional genus accuracy: FL 0.700/0.544/0.
- `panels:panels.tab_threepanel_gtdb_block` L1103-1104, 1125-1129: tab:threepanel GTDB block (phylum-std resolve original/k-means/random; Delta random-original [CI]): 100% 0.577/0.675/0.715, +0.139 [+0.107,+0.170]; 50% 0.579/0.
- `panels:panels.selfincl_recovery` L1141-1142: Self-inclusion: on the original panel the source genome is among its amplicon's 100 retained hits for 99.1% of amplicons (99.2%, 67,466 of 67,991, if a same-tax
- `panels:panels.selfincl_miss_causes` L1143-1147: Misses: in 8.2% (k-means) and 15.1% (random) of amplicons the 100 retained hits are all exact ties at identity 1.0 (tie eviction); in 18.0% and 16.9% the source
- `panels:panels.selfincl_genus_penalty` L1149-1150: Standardized self-inclusion genus accuracy falls by 0.07 to 0.08 (V4 -0.079, V7-V9 -0.065); on full-length 16S genus accuracy barely at all (-0.006).
- `panels:panels.selfincl_fulllength_credit` L1150-1151: On full-length 16S the tie eviction vanishes (no amplicon affected) but for 23% the source's full-length sequence is credited to another organism; on the random
- `panels:panels.gtdb_design` L1153-1156: GTDB arm: ~20% phylum-proportional subsample of each panel (2,000, 1,995 and 1,989 scored organisms; round(0.2n) per phylum) against 143,614 species-representat
- `panels:panels.gtdb_resolve` L1157: Against GTDB the alternative panels resolve more novel organisms: 0.715 and 0.675 vs 0.577 at full completeness after standardization; random-minus-original +0.
- `panels:panels.gtdb_survivors` L1158-1159: Median best surviving relative 92.1 and 91.4 ANI (random, k-means) vs 89.5 (original); original more often has no qualifying hit (16.1% vs 7.6%) or only a low-A
- `panels:disc.cosine_vs_identity` L1168: 16S cosine ordering within a near-identical cloud correlates only ~0.55 with true identity.
- `panels:disc.region_tradeoff_taxonomy` L1175: Novel organisms, V3-onward amplicons: right at family level 0.76-0.80 conditional; roughly a third receive no genus call; 54-59% correct family and 37-42% corre
- `panels:lim.run_scale` L1197: Conclusions rest on the full 10,000-genome, eight-region, 67,991-amplicon run with Wilson 95% intervals and paired McNemar comparisons.
- `panels:lim.archaeal_amplification` L1198: Four Bacteria-scoped regions amplify almost no archaea (35-76 archaeal amplicons each, <=6% extraction); V4, V4-V5, V4-V5_944R amplify 67-95% of archaea; V7-V9 
- `panels:lim.selfrecovery_gap` L1201-1202: Residual self-recovery gap of the seven sub-full-length regions under the union: 1,145 of 60,693 amplicons; 525 (46%) byte-identical to >=100 other genomes' cop
- `panels:lim.emptied_share` L1203, 1205: Depth ceiling accounts for 86.7% of abstentions under the union at depth 100 (5,598 of 6,457); 5,598 amplicons still emptied.
- `panels:lim.all_tied_records` L1206: For 0.75% of distinct amplicon sequences (453 of 60,461; 2.2% of the 67,991 records) all 100 retained hits tie at the best score, spanning a median of 19 specie
- `panels:lim.af_completeness` L1210: Reference aligned fraction scales ~0.99 x completeness when the source is in the reference, roughly half under self-exclusion; the 15% AF gate bounds MAG covera
- `panels:lim.modelseedpy_bug` L1226-1227: The prFBA function of ModelSEEDpy c2e31a0 raises an error on every call; fix on branch fix-prfba-single-model.
- `panels:lim.external_citations` L1230-1232: Panaroo: accessory genome nearly tenfold larger without error correction (Tonkin-Hill 2020); 95-99% soft-core convention (Roary, Page 2015).
