# Calibration of prFBA's p = k/n weighting

`scripts/calibrate_pfunc.py` (sha256 dccb91749094), 2026-10-20 (commit b9fa5c7a7c); read-only on existing benchmark artefacts. Question: prFBA gives each function x of a selected reference set of n genomes the weight p(x) = k_x/n. Is p a calibrated probability that x is present in the organism the amplicon came from?

## Key results

Population: controlled genus-blocked benchmark, exact_tie d100, EC, multi-genome cells (n >= 2; 7,973 cells, 3,851 source genomes, 1168 holdout groups), union denominator, unless stated. **95% CIs are a cluster bootstrap over holdout groups** (score_folds' pre-specified unit, 2000 draws, seed 1729); genus and source-genome intervals are sensitivities. The largest group holds 32% of the multi-cell points, so the group interval mostly reflects how often that one cluster is drawn; the genus interval is the middle ground.

1. **p is overconfident.** The Murphy reliability term (the part of p's Brier score that perfect recalibration would remove) is 0.0210 [0.0194, 0.0242] (genus [0.0191, 0.0230]; source [0.0199, 0.0224]) on a Brier of 0.1846 [0.1801, 0.1931]; out-of-fold isotonic recalibration recovers 0.0209 [0.0193, 0.0239] of it. The probability-scale calibration slope (weighted least squares of observed presence on p over the distinct-p bins, weights = points; 1 = calibrated, < 1 = too extreme) is 0.78 [0.76, 0.79] (genus [0.76, 0.79]; source [0.77, 0.79]), intercept 0.06 [0.04, 0.08]. The overconfidence is at p >= 1/2 (n = 2, p = 1/2 is present 0.449 [0.437, 0.461]; p = 1 from n >= 2 references 0.837 [0.826, 0.843]); low p is slightly under-confident (n = 4, p = 1/4: 0.286 [0.241, 0.312]). The logistic slope on logit p is set by where p = 1 is clipped and is secondary: 0.41 [0.40, 0.42] at eps 0.01, 0.28 [0.26, 0.28] at eps 0.001.
2. **p ranks functions correctly.** Within n = 2, 3 and 4 observed presence rises strictly with p, and adjacent bins' group CIs do not overlap; n = 3: p=1/3: 0.317 [0.301, 0.332]; p=2/3: 0.560 [0.523, 0.580]; p=3/3: 0.854 [0.834, 0.866].
3. **Skill over a constant confidence comes from n >= 3 only.** Pooled multi Brier skill vs the constant c = pooled precision is 0.061 [0.034, 0.074] (genus [0.046, 0.075]; source [0.051, 0.071]) (> 0). By bucket, **for two-way ties (n = 2) p does not beat a constant confidence**: skill 0.010 [-0.009, 0.024] (genus [-0.006, 0.025]; source [-0.001, 0.020]), CI covers 0. n = 2 is 71% of multi cells and 70% of their points. n = 3: 0.132 [0.079, 0.165]; n >= 4: 0.201 [0.171, 0.228].
4. **p = 1 is more often right when several references agree.** p = 1 from n >= 2 references is present 0.837 [0.826, 0.843] of the time, p = 1 from a single (cross-genus) reference 0.787 [0.772, 0.796]; paired difference 0.050 [0.043, 0.061] (genus [0.043, 0.057]; source [0.045, 0.055]) (> 0). Functions outside the union (p = 0) are present 0.039 [0.034, 0.043].
5. **Thresholding p barely moves the benchmark.** Against the pre-registered SESOI of 0.02 F1 (fold_plan meta.frozen.sesoi_delta_f1; score_folds' rule: exceeded iff |dF1| >= SESOI and the CI excludes 0), a majority vote (p >= 0.5) vs the union changes macro F1 by 0.012 [0.008, 0.015] at n = 3 (2.3% of the 57,021 scored cells; micro 0.016 [0.012, 0.018]): not exceeded (non-zero, below the SESOI); and by 0.025 [0.018, 0.030] (genus [0.018, 0.031]; source [0.022, 0.028]) at n >= 4 (1.7% of cells; micro 0.028 [0.022, 0.034]): exceeded (but the CI bound nearer 0, +0.018, is inside the SESOI). At n = 2 it is the union. Over all scored cells the rule changes the benchmark's macro F1 by 0.0007 [0.0005, 0.0009]: not exceeded (non-zero, below the SESOI). The intersection (p = 1; at n = 2 also p > 0.5) loses to the union in every n bucket (macro dF1 at n = 2 -0.022 [-0.027, -0.018], n >= 4 -0.024 [-0.038, -0.016]).
6. **EC is marginally worse calibrated than either reaction layer (by 0.0010 and 0.0008 on a reliability term of 0.021); the two reaction layers are indistinguishable** (reliability term ec 0.021 [0.019, 0.024]; rxn_ec 0.020 [0.018, 0.023]; rxn_role 0.020 [0.019, 0.023]; pgfam 0.098 [0.093, 0.103]; product 0.044 [0.041, 0.051]). Paired differences from the shared draws: ec - rxn_ec 0.0010 [0.0004, 0.0016]; ec - rxn_role 0.0008 [0.0004, 0.0015]; rxn_ec - rxn_role -0.0002 [-0.0005, 0.0004]. PGFam and the product-string proxy are worse calibrated (ec - pgfam -0.077 [-0.080, -0.073], ec - product -0.023 [-0.027, -0.021]). Skill vs c: ec 0.061 [0.034, 0.074]; rxn_ec 0.055 [0.029, 0.069]; rxn_role 0.036 [0.014, 0.048]; pgfam -0.119 [-0.166, -0.095]; product 0.017 [-0.026, 0.040] (covers 0).
7. **Skill exists only when the nearest reference is close**: rank-1 identity 0.865-0.945 (3,532 cells) -0.007 [-0.053, 0.018] (genus [-0.033, 0.014]; source [-0.023, 0.008]), CI covers 0: not distinguishable from 0; rank-1 identity 0.945-0.987 (4,441 cells) 0.111 [0.096, 0.125] (genus [0.096, 0.126]; source [0.100, 0.121]), > 0. Forced answers below the family floor (exact_tie_forced, rank-1 identity < 0.865, 537 multi cells): -0.438 [-0.540, -0.353], worse than a constant.
8. **Other runs** (EC multi skill vs c): legacy|d100 0.092 [0.060, 0.111]; graduated|d100 0.214 [0.193, 0.226]; exact_tie_guards|d100 0.021 [-0.007, 0.034] **covers 0** (excludes 0 under genus, source); exact_tie|d050 0.056 [0.043, 0.074]; exact_tie|d025 0.051 [0.013, 0.076]; exact_tie_forced|d100 0.045 [0.016, 0.061].
9. **As shipped (secondary; CIs over source genus)**: include-self multi skill 0.624 [0.600, 0.646] is a ceiling (the source genome votes); exclude_0.987 0.264 [0.248, 0.278], WLS slope 0.88 [0.86, 0.90], a post-hoc exclusion on truncated top-20 records.
10. **Exploratory, not pre-registered**: top-20 identity-weighted prevalence (tau 0.01) thresholded at 0.5 on all 57,021 exact_tie cells changes macro F1 vs the shipped union by 0.0252 [0.0229, 0.0275] (genus [0.0234, 0.0274]; source [0.0245, 0.0260]) (SESOI rule, were it pre-registered: exceeded); universe Brier 0.0627 [0.0615, 0.0647] vs 0.0889 [0.0870, 0.0911] for k/n. On non-saturated records (40,686 cells: the 200th coverage-qualified hit is already below the family floor; every record holds 200 hits, so this too assumes identity tracks align_score beyond rank 200) the macro F1 change is 0.0212 [0.0194, 0.0226]; that is the only truncation evidence (the 'exposed' flag in scores.json is a heuristic: the cutoff is by align_score, not identity). A hypothesis for a pre-registered test, not a result.

## Tables

exact_tie d100, multi-genome cells, union denominator. CIs over holdout groups. REL = Murphy reliability term; WLS = probability-scale slope; logit slope at eps 0.01 / 0.001.

| unit | REL | WLS slope | skill vs c | OOF isotonic gap | logit slope .01 / .001 |
|---|---|---|---|---|---|
| ec | 0.0210 [0.0194, 0.0242] | 0.78 [0.76, 0.79] | 0.061 [0.034, 0.074] | 0.0209 [0.0193, 0.0239] | 0.41 / 0.28 |
| rxn_ec | 0.0200 [0.0184, 0.0232] | 0.76 [0.74, 0.77] | 0.055 [0.029, 0.069] | 0.0199 [0.0182, 0.0229] | 0.40 / 0.27 |
| rxn_role | 0.0202 [0.0186, 0.0232] | 0.71 [0.69, 0.73] | 0.036 [0.014, 0.048] | 0.0200 [0.0184, 0.0228] | 0.38 / 0.25 |
| pgfam | 0.0980 [0.0932, 0.1029] | 0.89 [0.86, 0.91] | -0.119 [-0.166, -0.095] | 0.0977 [0.0928, 0.1025] | 0.52 / 0.35 |
| product | 0.0440 [0.0406, 0.0509] | 0.82 [0.80, 0.84] | 0.017 [-0.026, 0.040] | 0.0437 [0.0402, 0.0504] | 0.42 / 0.28 |

By n bucket (EC):

| bucket | cells (% of scored) | % of multi points | skill vs c | REL | WLS slope |
|---|---|---|---|---|---|
| n2 | 5,693 (10.0%) | 70% | 0.010 [-0.009, 0.024] | 0.0233 [0.0214, 0.0264] | 0.76 [0.73, 0.78] |
| n3 | 1,313 (2.3%) | 17% | 0.132 [0.079, 0.165] | 0.0162 [0.0133, 0.0218] | 0.82 [0.79, 0.84] |
| n4+ | 967 (1.7%) | 13% | 0.201 [0.171, 0.228] | 0.0163 [0.0132, 0.0211] | 0.81 [0.78, 0.85] |
| multi | 7,973 (14.0%) | 100% | 0.061 [0.034, 0.074] | 0.0210 [0.0194, 0.0242] | 0.78 [0.76, 0.79] |

Binary rules (EC), dF1 vs the union; SESOI 0.02:

| bucket | macro dF1 p>=0.5 | micro dF1 p>=0.5 | macro dF1 p=1 | verdict (p>=0.5, macro) |
|---|---|---|---|---|
| n2 | 0.0000 [0.0000, 0.0000] | 0.0000 [0.0000, 0.0000] | -0.0215 [-0.0271, -0.0181] | identical to the union |
| n3 | 0.0122 [0.0084, 0.0147] | 0.0159 [0.0120, 0.0184] | -0.0298 [-0.0416, -0.0167] | not exceeded (non-zero, below the SESOI) |
| n4+ | 0.0248 [0.0184, 0.0300] | 0.0278 [0.0225, 0.0342] | -0.0239 [-0.0383, -0.0162] | exceeded (but the CI bound nearer 0, +0.018, is inside the SESOI) |
| multi | 0.0050 [0.0037, 0.0058] | 0.0059 [0.0045, 0.0066] | -0.0232 [-0.0282, -0.0198] | not exceeded (non-zero, below the SESOI) |
| all | 0.0007 [0.0005, 0.0009] | 0.0008 [0.0006, 0.0010] | -0.0032 [-0.0037, -0.0026] | not exceeded (non-zero, below the SESOI) |

Observed presence at each exact p (EC; n >= 5 and the decile rows are in reliability.csv):

| pop | bin | points | cells | observed [group CI] | genus CI |
|---|---|---|---|---|---|
| single | n=1,p=1/1 | 22,910,365 | 49,048 | 0.787 [0.772, 0.796] | [0.781, 0.792] |
| multi | p=0 (outside union) | 13,365,951 | 7,973 | 0.039 [0.034, 0.043] | [0.037, 0.041] |
| n2 | n=2,p=1/2 | 704,323 | 5,620 | 0.449 [0.437, 0.461] | [0.437, 0.462] |
| n2 | n=2,p=2/2 | 2,398,260 | 5,693 | 0.828 [0.817, 0.835] | [0.821, 0.836] |
| n3 | n=3,p=1/3 | 127,679 | 1,311 | 0.317 [0.301, 0.332] | [0.300, 0.332] |
| n3 | n=3,p=2/3 | 120,989 | 1,312 | 0.560 [0.523, 0.580] | [0.535, 0.584] |
| n3 | n=3,p=3/3 | 529,365 | 1,313 | 0.854 [0.834, 0.866] | [0.841, 0.866] |
| n4+ | n=4,p=1/4 | 39,202 | 520 | 0.286 [0.241, 0.312] | [0.261, 0.314] |
| n4+ | n=4,p=2/4 | 28,192 | 519 | 0.420 [0.389, 0.440] | [0.389, 0.448] |
| n4+ | n=4,p=3/4 | 33,181 | 520 | 0.599 [0.566, 0.641] | [0.557, 0.639] |
| n4+ | n=4,p=4/4 | 193,961 | 520 | 0.869 [0.853, 0.893] | [0.846, 0.893] |

Every primary run (EC, multi, union):

| run | multi cells | skill vs c [group] | genus | source | REL |
|---|---|---|---|---|---|
| exact_tie|d100 | 7,973 | 0.061 [0.034, 0.074] | [0.046, 0.075] | [0.051, 0.071] | 0.0210 [0.0194, 0.0242] |
| legacy|d100 | 22,388 | 0.092 [0.060, 0.111] | [0.079, 0.105] | [0.085, 0.099] | 0.0189 [0.0167, 0.0231] |
| graduated|d100 | 30,747 | 0.214 [0.193, 0.226] | [0.204, 0.223] | [0.208, 0.219] | 0.0110 [0.0094, 0.0140] |
| exact_tie_guards|d100 | 6,867 | 0.021 [-0.007, 0.034] | [0.004, 0.037] | [0.010, 0.031] | 0.0245 [0.0226, 0.0281] |
| exact_tie|d050 | 5,850 | 0.056 [0.043, 0.074] | [0.038, 0.074] | [0.045, 0.066] | 0.0219 [0.0197, 0.0243] |
| exact_tie|d025 | 3,783 | 0.051 [0.013, 0.076] | [0.033, 0.069] | [0.038, 0.063] | 0.0252 [0.0216, 0.0319] |
| exact_tie_forced|d100 | 8,513 | 0.045 [0.016, 0.061] | [0.029, 0.060] | [0.035, 0.056] | 0.0236 [0.0217, 0.0276] |

Exploratory top-k (EC, all exact_tie cells; dF1 = macro F1 of q >= 0.5 minus the shipped union):

| k | universe Brier | d vs k/n | dF1, all records | dF1, non-saturated records |
|---|---|---|---|---|
| 5 | 0.0674 [0.0662, 0.0691] | -0.0215 [-0.0223, -0.0205] | 0.0177 [0.0162, 0.0199] | 0.0153 [0.0142, 0.0166] |
| 10 | 0.0643 [0.0631, 0.0661] | -0.0246 [-0.0255, -0.0233] | 0.0225 [0.0206, 0.0250] | 0.0194 [0.0179, 0.0207] |
| 20 | 0.0627 [0.0615, 0.0647] | -0.0262 [-0.0273, -0.0247] | 0.0252 [0.0229, 0.0275] | 0.0212 [0.0194, 0.0226] |
| 50 | 0.0618 [0.0605, 0.0639] | -0.0271 [-0.0283, -0.0254] | 0.0268 [0.0242, 0.0291] | 0.0222 [0.0202, 0.0237] |
| 100 | 0.0616 [0.0603, 0.0637] | -0.0273 [-0.0285, -0.0256] | 0.0271 [0.0245, 0.0294] | 0.0224 [0.0204, 0.0238] |

## Method

- **Population**: score_folds' scored cells (`cells.csv.gz`, space ec, outcome scored) of exact_tie, legacy, graduated, exact_tie_guards (d100), exact_tie d050 / d025 and exact_tie forced. Reps = `predictions.json pred[mode][md5].reps`; truth = the cell's `src_genome_id` (`_cache_panel_ec` through norm_ec); p is recomputed from the reps' function sets. Abstentions have no p. Strata: all, single, multi, n2, n3, n4+, crossed with novelty_bin, rank1_rank, region and a twin-free one-cell-per-md5 population.
- **Units**: `ec` (2,235 scored ECs); `rxn_ec` (score_folds' rxn space, w(R) = max over member ECs); `rxn_role` (PGFam role route, w(R) = max over member PGFams); `pgfam`; `product` (PGFam family_product strings minus empty and generic 'hypothetical protein' variants, 2919 distinct strings dropped; a proxy for the RAST strings the merge counts). Union denominator = the union's own functions; universe = the fixed vocabulary with p = 0 as a bin (ec, rxn_ec, rxn_role). The max-over-genes reaction weight differs from genome-level prevalence on 0.13% (rxn_ec) and 1.73% (rxn_role) of multi-cell union reactions.
- **Scores**: Brier of p, of the union (= 1 - precision) and of the constant c = pooled precision (in-sample; an out-of-fold c is also reported); skill = 1 - Brier(p)/Brier(c); Murphy decomposition over distinct p (exact; residual reported); WLS calibration line; logistic slope on logit(clip(p, eps)); ECE; out-of-fold isotonic / Platt recalibration (4 folds fit, 5th scored; refitting noise is not in the CI); binary rules union, p >= 0.5, p > 0.5, p = 1, micro and macro (= the benchmark's per-cell mean) F1.
- **Uncertainty**: score_folds.cluster_bootstrap's draw verbatim (`default_rng(seed).integers(0, K, size=(B, K))` over the stratum's sorted cluster keys), percentile 95%. Primary unit = holdout group (`group_id`); sensitivities = genus (`ci_genus`) and source genome (`ci_src`) in scores.json and reliability.csv. Sources nest in genera and genera in groups (asserted). One set of draws per stratum serves every unit, so unit contrasts (`_unit_contrasts`) are paired. The as-shipped secondary has no holdout groups; there the primary unit is the source genus.
- **Leakage** (asserted, every cell of every run): the source genome is never a rep or a top-k neighbour, and no rep or neighbour shares its genus, species, taxon_id or holdout group.

## Caveats

- Genus-blocked folds: every primary number is cross-genus transfer and pessimistic for deployment; single-reference p = 1 is not the manuscript's conspecific estimate. The panel is species-dereplicated, so multi sets are multi-species.
- One holdout group (grp00000, all of fold 0) carries 32% of the multi points, so group CIs mostly reflect how often that group is drawn; across the headline metrics of the five units they are 1.3-3.2x as wide as source-genome CIs.
- Absence = absence of BV-BRC PGFam annotation on both sides; partial ECs dropped.
- The as-shipped secondary uses pseudo-folds for OOF recalibration and predates the 20 -> 100 retained-hit rerun; include-self is a ceiling.
- Top-k is exploratory; its truncation evidence is the non-saturated-records population only, and that check is itself a heuristic, not a bound: it assumes identity follows align_score order beyond rank 200, which fails near the floor in a large share of those cells.

## Outputs and assertions

`scores.json` (per run x unit x stratum; `lo`/`hi` = primary unit, `ci_genus`, `ci_src`; `clusters` gives the cluster counts and the largest cluster's share), `reliability.csv` (long format, one stratification at a time), `assertions.json` (every join, invariant and leakage check; all required checks zero: **True**).
