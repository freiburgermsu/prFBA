# Objective tuning (Stage O): results

> Filled by `scripts/objective_tuning_results.py` from `/home/freiburger/Documents/prFBA/region_validation/shared_reference/data/objective_tuning` on 2026-11-16T14:51:15-0600. Check every number by an independent re-derivation before quoting it.



Design: `PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md`, section "Stage O" (declared 2026-11-02, post hoc, binding). Harness: `scripts/objective_tuning.py`. Selection and confirmation: `scripts/objective_tuning_results.py`. Outputs: `data/objective_tuning/{per_amplicon.jsonl, meta.json, summary.json}`. Summary generated 2026-11-16T14:51:15-0600.

## 1. What was run

- **Sample.** The weighted-objective test's amplicons (`data/weighted_objective/population.json`): 360 rows, status {"ok": 360}. The gapfilled model of each selection's union EC set was rebuilt as before: ModelSEEDpy c2e31a0, templates 21a8084, reversals applied, complete medium, biomass >= 0.95 x max.
- **Split.** A holdout group goes to the tune half if sha256("1729|" + group_id) is even, else to the confirm half. Rows by half: {"tune": 167, "confirm": 193}.
- **Probability sources.** `kn` is the shipped k/n. `est` is Stage E's estimator at role level, from `/home/freiburger/Documents/prFBA/region_validation/shared_reference/data/pfunc_estimators/stage_o_weights_E3.json` (sha256 ead0d11eeb94): E3. `est_shuffled` is est permuted among the gene-backed reactions, using the prior shuffled arm's permutation, one per amplicon. `two_tier` sets gene-backed p = 1. `oracle` sets gene-backed p = 1 if truth-supported, else 0.
- **Factors.**
  - C, the cost on rxn reactions: `lin` = max(0.01, 1 − p) (shipped) or `log` = max(0.01, ln max(p, 0.01) / ln 0.01).
  - G, the p of gapfilled reactions, in every arm: 0 (shipped), g_meas = 0.083538 (1021 / 12222 over 167 tune-half rows of the prior run), or 0.2.
  - X, the exchange weight: 100 (shipped), 10 or 1.
  - In total: 5 P × 2 C × 3 G × 3 X = 90 arms per amplicon, each on a fresh copy of the model.
- **Measure.** The total truth-absent gene-backed flux on the optimal face (objective <= opt + 1e-6), as its minimum, midpoint and maximum, each divided by D. D is the gene-backed flux of the two-tier, linear, ex-100 arm, one constant per amplicon. The guard metrics are total exchange flux, gene-backed flux, flux precision and gapfilled flux.
- **Pooling.** A (primary): the n >= 2 strata weighted n2 0.6975, n3 0.1869, n4+ 0.1156. B (reported, not used to select): all answered amplicons, n1 0.8680, n2 0.0921, n3 0.0247, n4+ 0.0153.

## 2. Checks

| check | result |
|---|---|
| rows ok / total | 360 / 360 |
| rebuilt model = prior row (n_new_rxn, growth 1e-6, n_gene_backed, n_gapfilled) | 358 / 360 |
| D (two_tier/lin/G0/X100 gene-backed flux) bitwise = prior | 357 / 360 |
| kn/lin/G0/X100 objective bitwise = prior `real` | 357 / 360 |
| kn/lin/G0/X100 returned vertex (gene-backed, truth-absent flux) bitwise = prior | 357 / 360, 357 / 360 |
| kn/lin/G0/X100 face min/mid/max within 1e-9 of prior `real` | 356 / 360 (max abs diff mid 0.0007726599323445715) |
| two_tier/lin/G0/X100 objective bitwise = prior `two_tier` | 357 / 360 |
| two_tier/lin/G0/X100 face within 1e-9 of prior | 357 / 360 (max abs diff mid 0.0011824526945120872) |
| coefficient map of the new function = prfba_released_copy (shipped arm) | 360 / 360 |
| released_check pass; = prior objective_released bitwise | 360 / 360; 357 / 360 |
| determinism pairs identical (same LP) | 3240 / 3240 |
| prior arms re-run through the prior code, bitwise (repro) | not run |
| arm solves optimal | 32400 / 32400 |
| biomass at 0.95 x max (min, max ratio; off floor > 1e-6) | 0.9499999999999998, 0.9500000000000001; 0 |
| face LP failures (after one cold retry); retried | 0; 0 |
| gapfilled reactions that are gene-backed | 0 |
| rows with alpha > 0 and a null clade (literal formula) | 0 |
| wall mean / max per amplicon; largest worker RSS | 72.2 s / 109.5 s; 2.211 GiB |

## 3. Selection (tune half)

Tune rows by stratum: {"n2": 46, "n3": 49, "n4+": 42, "n1": 30}.

**Eligibility.** A configuration is eligible if its median total exchange flux and its median gene-backed flux are each within 10% of the shipped arm's. 36 of 36 configurations are eligible. By exchange weight: {"X100": 12, "X10": 12, "X1": 12}.

| configuration | exchange flux ratio | gene-backed flux ratio | precision ratio | eligible | pooled A | pooled B |
|---|---|---|---|---|---|---|
| kn/lin/G0/X100 | 1.000 | 1.000 | 1.000 | yes | 0.1362 | 0.1883 |
| kn/lin/G0/X10 | 1.000 | 1.000 | 0.997 | yes | 0.1351 | 0.1881 |
| kn/lin/G0/X1 | 1.005 | 0.992 | 1.002 | yes | 0.1339 | 0.1883 |
| kn/lin/Gmeas/X100 | 1.000 | 1.000 | 1.000 | yes | 0.1362 | 0.1883 |
| kn/lin/Gmeas/X10 | 1.000 | 1.000 | 0.997 | yes | 0.1351 | 0.1881 |
| kn/lin/Gmeas/X1 | 1.005 | 0.992 | 1.002 | yes | 0.1340 | 0.1883 |
| kn/lin/G0.2/X100 | 1.000 | 1.000 | 0.997 | yes | 0.1362 | 0.1883 |
| kn/lin/G0.2/X10 | 1.000 | 1.000 | 0.997 | yes | 0.1351 | 0.1881 |
| kn/lin/G0.2/X1 | 1.005 | 0.994 | 1.004 | yes | 0.1339 | 0.1883 |
| kn/log/G0/X100 | 1.000 | 1.000 | 1.002 | yes | 0.1364 | 0.1883 |
| kn/log/G0/X10 | 1.004 | 0.997 | 0.997 | yes | 0.1370 | 0.1884 |
| kn/log/G0/X1 | 1.004 | 0.994 | 0.997 | yes | 0.1354 | 0.1885 |
| kn/log/Gmeas/X100 | 1.000 | 1.000 | 1.002 | yes | 0.1363 | 0.1883 |
| kn/log/Gmeas/X10 | 1.004 | 0.997 | 0.997 | yes | 0.1369 | 0.1884 |
| kn/log/Gmeas/X1 | 1.004 | 0.997 | 0.997 | yes | 0.1355 | 0.1885 |
| kn/log/G0.2/X100 | 1.000 | 1.000 | 0.997 | yes | 0.1363 | 0.1883 |
| kn/log/G0.2/X10 | 1.004 | 0.997 | 0.997 | yes | 0.1369 | 0.1884 |
| kn/log/G0.2/X1 | 1.004 | 0.997 | 1.003 | yes | 0.1356 | 0.1886 |
| est/lin/G0/X100 | 1.000 | 1.018 | 1.017 | yes | 0.1234 | 0.1777 |
| est/lin/G0/X10 | 1.000 | 1.018 | 1.023 | yes | 0.1222 | 0.1775 |
| est/lin/G0/X1 | 1.004 | 1.011 | 1.029 | yes | 0.1196 | 0.1753 |
| est/lin/Gmeas/X100 | 1.000 | 1.018 | 1.017 | yes | 0.1234 | 0.1777 |
| est/lin/Gmeas/X10 | 1.000 | 1.018 | 1.023 | yes | 0.1223 | 0.1776 |
| est/lin/Gmeas/X1 | 1.004 | 1.011 | 1.029 | yes | 0.1196 | 0.1753 |
| est/lin/G0.2/X100 | 1.000 | 1.018 | 1.023 | yes | 0.1234 | 0.1777 |
| est/lin/G0.2/X10 | 1.000 | 1.018 | 1.023 | yes | 0.1223 | 0.1776 |
| est/lin/G0.2/X1 | 1.004 | 1.011 | 1.029 | yes | 0.1196 | 0.1753 |
| est/log/G0/X100 | 1.000 | 1.030 | 1.004 | yes | 0.1257 | 0.1789 |
| est/log/G0/X10 | 1.004 | 1.030 | 1.017 | yes | 0.1254 | 0.1789 |
| est/log/G0/X1 | 1.004 | 1.016 | 1.015 | yes | 0.1249 | 0.1791 |
| est/log/Gmeas/X100 | 1.000 | 1.030 | 1.004 | yes | 0.1257 | 0.1789 |
| est/log/Gmeas/X10 | 1.004 | 1.030 | 1.004 | yes | 0.1255 | 0.1789 |
| est/log/Gmeas/X1 | 1.004 | 1.016 | 1.016 | yes | 0.1250 | 0.1791 |
| est/log/G0.2/X100 | 1.000 | 1.030 | 1.015 | yes | 0.1257 | 0.1789 |
| est/log/G0.2/X10 | 1.004 | 1.030 | 1.015 | yes | 0.1255 | 0.1789 |
| est/log/G0.2/X1 | 1.004 | 1.016 | 1.016 | yes | 0.1251 | 0.1791 |

**The G factor.** G is read only on rows where the gapfilled flux varies across G's levels:

| family | eligible G | rows where gapfilled flux varies (n2/n3/n4+) | choice | why |
|---|---|---|---|---|
| kn/lin/*/X100 | G0, Gmeas, G0.2 | 3/4/2 | kn/lin/G0/X100 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| kn/lin/*/X10 | G0, Gmeas, G0.2 | 3/4/2 | kn/lin/G0/X10 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| kn/lin/*/X1 | G0, Gmeas, G0.2 | 6/6/5 | kn/lin/G0/X1 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| kn/log/*/X100 | G0, Gmeas, G0.2 | 4/3/4 | kn/log/G0/X100 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| kn/log/*/X10 | G0, Gmeas, G0.2 | 6/3/5 | kn/log/G0/X10 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| kn/log/*/X1 | G0, Gmeas, G0.2 | 11/3/6 | kn/log/G0/X1 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| est/lin/*/X100 | G0, Gmeas, G0.2 | 1/1/0 | est/lin/G0/X100 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| est/lin/*/X10 | G0, Gmeas, G0.2 | 2/1/0 | est/lin/G0/X10 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| est/lin/*/X1 | G0, Gmeas, G0.2 | 3/1/0 | est/lin/G0/X1 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| est/log/*/X100 | G0, Gmeas, G0.2 | 2/4/1 | est/log/G0/X100 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| est/log/*/X10 | G0, Gmeas, G0.2 | 4/4/2 | est/log/G0/X10 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| est/log/*/X1 | G0, Gmeas, G0.2 | 7/5/3 | est/log/G0/X1 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |

**Selection.** Selected: **est/lin/G0/X1** (pooled A 0.1196 against the shipped 0.1362; lowest pooled A 0.1196 at est/lin/G0/X1; tie set within 0.002: est/lin/G0/X1; 2 departure(s) from the shipped configuration). The joint selection over all eligible configurations on all rows (sensitivity) picks est/lin/G0/X1 (agrees).

- Ranking under pooling A (eligible configurations): est/lin/G0/X1 (0.1196), est/lin/Gmeas/X1 (0.1196), est/lin/G0.2/X1 (0.1196), est/lin/G0/X10 (0.1222), est/lin/G0.2/X10 (0.1223), est/lin/Gmeas/X10 (0.1223), est/lin/G0.2/X100 (0.1234), est/lin/Gmeas/X100 (0.1234), est/lin/G0/X100 (0.1234), est/log/G0/X1 (0.1249)
- Ranking under pooling B (reported only): est/lin/Gmeas/X1 (0.1753), est/lin/G0/X1 (0.1753), est/lin/G0.2/X1 (0.1753), est/lin/G0/X10 (0.1775), est/lin/G0.2/X10 (0.1776), est/lin/Gmeas/X10 (0.1776), est/lin/G0.2/X100 (0.1777), est/lin/Gmeas/X100 (0.1777), est/lin/G0/X100 (0.1777), est/log/G0/X10 (0.1789)
- Oracle headroom captured on the tune half: 0.679 of the oracle's headroom at its own cell; 0.910 in the shipped cell's frame (the shipped configuration: 0.259)

## 4. Confirmation (confirm half)

Confirm rows by stratum: {"n2": 54, "n3": 51, "n4+": 58, "n1": 30}. Intervals are 95% percentile cluster bootstraps over holdout groups (2,000 draws, seed 1729 + 1000 × stratum index). The SESOI is 0.01 of D.

| contrast | pooling | face mid | face min | face max | as returned (taf_gb) | label |
|---|---|---|---|---|---|---|
| selected-shipped (est/lin/G0/X1 - kn/lin/G0/X100) | A | -0.0172 [-0.0221, -0.0118] | -0.0076 [-0.0133, -0.0047] | -0.0268 [-0.0323, -0.0178] | -0.0149 [-0.0203, -0.0100] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| selected-shipped (est/lin/G0/X1 - kn/lin/G0/X100) | B | -0.0065 [-0.0119, +0.0032] | -0.0014 [-0.0075, +0.0085] | -0.0116 [-0.0170, -0.0014] | -0.0106 [-0.0159, -0.0003] | inconclusive: the midpoint CI crosses 0 and a SESOI bound |
| selected-two_tier (est/lin/G0/X1 - two_tier/lin/G0/X1) | A | -0.0222 [-0.0268, -0.0159] | -0.0100 [-0.0149, -0.0063] | -0.0345 [-0.0404, -0.0244] | -0.0221 [-0.0276, -0.0170] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| selected-two_tier (est/lin/G0/X1 - two_tier/lin/G0/X1) | B | -0.0070 [-0.0124, +0.0025] | -0.0017 [-0.0077, +0.0082] | -0.0124 [-0.0179, -0.0024] | -0.0109 [-0.0156, -0.0016] | inconclusive: the midpoint CI crosses 0 and a SESOI bound |
| selected-shuffled (est/lin/G0/X1 - est_shuffled/lin/G0/X1) | A | -0.0283 [-0.0382, -0.0195] | -0.0279 [-0.0379, -0.0192] | -0.0286 [-0.0386, -0.0198] | -0.0235 [-0.0322, -0.0083] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| selected-shuffled (est/lin/G0/X1 - est_shuffled/lin/G0/X1) | B | -0.0200 [-0.0282, -0.0069] | -0.0210 [-0.0292, -0.0068] | -0.0189 [-0.0272, -0.0070] | -0.0182 [-0.0253, -0.0056] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| oracle-two_tier@selected_cell (oracle/lin/G0/X1 - two_tier/lin/G0/X1) | A | -0.0324 [-0.0372, -0.0242] | -0.0201 [-0.0253, -0.0143] | -0.0446 [-0.0507, -0.0327] | -0.0302 [-0.0356, -0.0224] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| oracle-two_tier@selected_cell (oracle/lin/G0/X1 - two_tier/lin/G0/X1) | B | -0.0200 [-0.0260, -0.0098] | -0.0136 [-0.0193, -0.0050] | -0.0265 [-0.0334, -0.0139] | -0.0207 [-0.0267, -0.0098] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| shipped-two_tier@shipped_cell (kn/lin/G0/X100 - two_tier/lin/G0/X100) | A | -0.0057 [-0.0074, -0.0036] | -0.0031 [-0.0042, -0.0012] | -0.0083 [-0.0116, -0.0054] | -0.0074 [-0.0110, -0.0048] | detectable but below SESOI: the midpoint CI excludes 0, |mean| < 0.01, and the CIs do not lie inside +-0.01 |
| shipped-two_tier@shipped_cell (kn/lin/G0/X100 - two_tier/lin/G0/X100) | B | -0.0008 [-0.0010, -0.0005] | -0.0004 [-0.0005, -0.0002] | -0.0011 [-0.0015, -0.0007] | -0.0010 [-0.0014, -0.0006] | equivalent within SESOI: the midpoint and both endpoint CIs lie inside +-0.01, nonzero (midpoint CI excludes 0) |
| oracle-two_tier@shipped_cell (oracle/lin/G0/X100 - two_tier/lin/G0/X100) | A | -0.0245 [-0.0288, -0.0183] | -0.0123 [-0.0170, -0.0085] | -0.0367 [-0.0426, -0.0266] | -0.0236 [-0.0290, -0.0176] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| oracle-two_tier@shipped_cell (oracle/lin/G0/X100 - two_tier/lin/G0/X100) | B | -0.0161 [-0.0207, -0.0078] | -0.0096 [-0.0131, -0.0033] | -0.0227 [-0.0285, -0.0115] | -0.0182 [-0.0243, -0.0075] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |

Shuffled comparator: P = est: est_shuffled is its own shuffled arm.

Per stratum:

| contrast | stratum | n (groups) | face mid | face min | face max |
|---|---|---|---|---|---|
| selected-shipped | n2 | 54 (26) | -0.0166 [-0.0228, -0.0098] | -0.0086 [-0.0160, -0.0046] | -0.0247 [-0.0318, -0.0137] |
| selected-shipped | n3 | 51 (26) | -0.0183 [-0.0324, -0.0099] | -0.0052 [-0.0191, +0.0012] | -0.0313 [-0.0476, -0.0195] |
| selected-shipped | n4+ | 58 (21) | -0.0188 [-0.0232, -0.0047] | -0.0056 [-0.0147, -0.0019] | -0.0320 [-0.0401, -0.0071] |
| selected-shipped | n1 | 30 (20) | -0.0048 [-0.0113, +0.0063] | -0.0004 [-0.0075, +0.0108] | -0.0093 [-0.0157, +0.0023] |
| selected-two_tier | n2 | 54 (26) | -0.0204 [-0.0268, -0.0125] | -0.0105 [-0.0173, -0.0059] | -0.0303 [-0.0385, -0.0175] |
| selected-two_tier | n3 | 51 (26) | -0.0252 [-0.0342, -0.0159] | -0.0092 [-0.0200, +0.0000] | -0.0411 [-0.0528, -0.0279] |
| selected-two_tier | n4+ | 58 (21) | -0.0288 [-0.0346, -0.0089] | -0.0086 [-0.0137, +0.0020] | -0.0490 [-0.0589, -0.0180] |
| selected-two_tier | n1 | 30 (20) | -0.0047 [-0.0111, +0.0063] | -0.0004 [-0.0075, +0.0108] | -0.0090 [-0.0155, +0.0022] |
| selected-shuffled | n2 | 54 (26) | -0.0257 [-0.0400, -0.0142] | -0.0256 [-0.0399, -0.0142] | -0.0258 [-0.0401, -0.0143] |
| selected-shuffled | n3 | 51 (26) | -0.0378 [-0.0504, -0.0245] | -0.0363 [-0.0495, -0.0236] | -0.0393 [-0.0513, -0.0254] |
| selected-shuffled | n4+ | 58 (21) | -0.0284 [-0.0351, -0.0103] | -0.0284 [-0.0351, -0.0103] | -0.0284 [-0.0351, -0.0103] |
| selected-shuffled | n1 | 30 (20) | -0.0187 [-0.0281, -0.0035] | -0.0199 [-0.0293, -0.0034] | -0.0175 [-0.0271, -0.0035] |

Guard metrics, selected − shipped:

| contrast | pooling | exchange flux | gene-backed flux | flux precision | gapfilled flux |
|---|---|---|---|---|---|
| selected-shipped | A | +2.5045 [+1.2036, +3.3873] | +1.1801 [-12.9610, +7.9732] | +0.0149 [+0.0100, +0.0203] | -0.2501 [-0.7399, -0.0900] |
| selected-shipped | B | +0.9219 [+0.4818, +1.5527] | +8.8665 [+1.6885, +13.3366] | +0.0106 [+0.0003, +0.0159] | -0.2961 [-1.0766, -0.0431] |

Oracle headroom captured on the confirm half (pooling A): selected (own cell) 0.687; shipped 0.233.

**Degeneracy** (FVA on the optimal face, finalists only, second pass):

not run (second pass: objective_tuning.py --half confirm --arms <finalists> --fva-arms <finalists> --outdir <new dir>; then --fva-dir here)

## 5. Interpretation choices (fixed before the full run)

- **Eligibility.** The guard compares medians over all tune-half rows, n1 included: the ratio of the configuration's median to the shipped arm's, for total exchange flux and for gene-backed flux. An eligible configuration must also be optimal on every row, with bio1 = 0.95 × max on every row.
- **The G factor.** G is read per (P, C, X) family, on the n >= 2 rows where the face midpoint of the gapfilled flux differs across the three G levels by more than 1e-6 × max(1, D). On those rows the stratum weights are renormalised over the strata present. Ties within 0.002 go to G0. A family with no such row keeps G0. As a sensitivity, the joint selection over all 36 configurations on all rows is reported beside it.
- **Selection across families.** The lowest pooled A over all tune rows of the n >= 2 strata wins. Ties within 0.002 go to the configuration with the fewest departures from kn/lin/G0/X100, then to the lowest pooled A.
- **Headroom.** Captured headroom is (selected − two_tier) / (oracle − two_tier), using the controls at the selected cell. The shipped-cell version is given beside it.
- **Shuffled comparator.** The design has no kn-shuffled arm. If kn is selected, the shuffled comparator is est_shuffled at the same cell, and this is flagged.
- **Null clade.** If Stage E gives an amplicon a null clade and α > 0, the literal formula applies: α stays in the denominator. The count of such rows is in the checks.

## 6. Limitations

- The design's caveats apply: transfer from the ~7,000-genome benchmark reference to the ~315,000-genome production reference, and the out-of-scope expansion of support.
- The weighted-objective test's limitations also apply (`WEIGHTED_OBJECTIVE_RESULTS.md`, section 10). Truth is an annotation-level proxy, the EC route stands in for function strings, there is one medium, and so on.
- Pooling B rests on the n1 stratum: about 30 rows per half.
