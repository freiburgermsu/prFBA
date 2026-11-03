# Objective tuning (Stage O): results

Run 2026-11-03 (unit `otune-run`, 08:50 to 10:27; FVA pass `otune-fva`). Design:
`PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md`, section "Stage O" (declared 2026-11-02, post hoc).
**Independently re-derived on 2026-11-03** (separate code from the raw rows,
`data/objective_tuning/verification/VERIFICATION.md`): all 972 confirm-half means and interval
bounds, the guard table, the rankings, the headroom, the checks and the FVA degeneracy match
exactly. Its nine corrections are applied below. Sections 1 to 6 were filled by
`scripts/objective_tuning_results.py`; the answer that follows is written from them.

## The answer

**Selected on the tune half (167 amplicons): the Stage E estimator (E3), the shipped linear
cost, gapfilled reactions left at p = 0, exchange weight 1.** On the held-out confirm half (193
amplicons, n >= 2 pooling):

- **It moves flux off reactions the true organism lacks, past the declared SESOI.** Against
  the shipped objective the truth-absent flux falls by −0.0172 [−0.0221, −0.0118] of D at the
  face midpoint, with both face endpoints below 0 (minimum −0.0076, maximum −0.0268). The tune
  half gave −0.0165, so there is no winner's curse.
- **The gain is the probability information.** It beats the same costs on shuffled
  probabilities by −0.0282 [−0.0382, −0.0195], and its own two-tier arm by −0.0222. It captures
  0.69 of the oracle's headroom; the shipped objective captures 0.23. The mechanism is a
  precision gain: truth-supported gene-backed flux rises, gene-less and gapfilled flux fall.
- **Where it comes from (exploratory).** The estimator alone, at the shipped exchange weight,
  gives −0.0128 [−0.0186, −0.0088]; lowering the exchange weight to 1 adds −0.0043. The
  normalised log cost is slightly worse than linear (+0.0039). Pricing gapfilled reactions at
  the measured 0.084 or at 0.2 instead of 0 changes nothing (+0.0000): they carry flux because
  the model needs them to grow, so their price cannot steer.
- **Degeneracy.** The selected objective leaves 13.6 gene-backed reactions free on the optimal
  face against 35.4 shipped, but shuffled costs leave 2.1: any gradation breaks ties, so this
  is not evidence of the probabilities' information.
- **Not established for single-genome selections**, where prFBA usually runs (86.8% of
  answered amplicons): −0.0048 [−0.0113, +0.0063] on the confirm half (30 rows), after −0.0124
  on the tune half. Pooled at the population's weights the effect is inconclusive.
- **Caveats** (section 6): the confirm half was held out from selection, not from the design
  (its amplicons were in the weighted-objective test, whose exchange-weight-1 arm had been
  seen); one holdout group holds 96 of the 193 confirm rows (without it −0.0164
  [−0.0236, −0.0100]); uptake composition under exchange weight 1 was not recorded; truth is an
  annotation-level proxy; the benchmark reference is ~7,000 genomes against ~315,000 in
  production.

**Recommendation.** If prFBA adopts a change, adopt the estimator first: it carries most of
the gain at the shipped exchange weight, and it gives single-genome selections graded
probabilities. Report it next to the shipped objective with these confirm-half effect sizes.
Lowering the exchange weight is a smaller, separate change that alters uptake, and it should be
tested on measured exchange data before it ships. A log cost and a gapfill price are not
supported by this test.

Design: `PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md`, section "Stage O" (declared 2026-11-02, post hoc, binding). Harness: `scripts/objective_tuning.py`. Selection and confirmation: `scripts/objective_tuning_results.py`. Outputs: `data/objective_tuning/{per_amplicon.jsonl, meta.json, summary.json}`. Summary generated 2026-11-03T11:27:48-0600.

## 1. What was run

- **Sample.** The weighted-objective test's amplicons (`data/weighted_objective/population.json`): 360 rows, status {"ok": 360}. The gapfilled model of each selection's union EC set was rebuilt as before (the prior LP exactly on 357 / 360 rows; see the checks): ModelSEEDpy c2e31a0, templates 21a8084, reversals applied, complete medium, biomass >= 0.95 x max.
- **Split.** A holdout group goes to the tune half if sha256("1729|" + group_id) is even, else to the confirm half. Rows by half: {"tune": 167, "confirm": 193}.
- **Probability sources.** `kn` is the shipped k/n. `est` is Stage E's estimator at role level, from `/home/freiburger/Documents/prFBA/region_validation/shared_reference/data/pfunc_estimators/stage_o_weights_E3.json` (sha256 330750c4c28b): E3. `est_shuffled` is est permuted among the gene-backed reactions, using the prior shuffled arm's permutation, one per amplicon. `two_tier` sets gene-backed p = 1. `oracle` sets gene-backed p = 1 if truth-supported, else 0.
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
| rebuilt model = prior row (n_new_rxn, growth 1e-6, n_gene_backed, n_gapfilled) | 358 / 360 (fail: 1674805.4__V3-V4 tune n4+; 1798454.3__V4 confirm n4+). A third row, 1017270.21__V4-V5 (confirm n4+), passes these counts but rebuilds a different LP (D 667.78 vs 656.27; shipped and released objectives differ), so 357 / 360 rows rebuild the prior LP. Excluding all three changes neither the selection nor any conclusion (selected − shipped, A: −0.0171 [−0.0220, −0.0119]). |
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
| wall mean / max per amplicon; largest worker RSS | 131.5 s / 291.1 s; 2.141 GiB |

## 3. Selection (tune half)

Tune rows by stratum: {"n2": 46, "n3": 49, "n4+": 42, "n1": 30}.

**Eligibility.** A configuration is eligible if its median total exchange flux and its median gene-backed flux are each within 10% of the shipped arm's. 36 of 36 configurations are eligible. By exchange weight: {"X100": 12, "X10": 12, "X1": 12}. The guard did not bind: no candidate's median ratio departs from 1 by more than 0.45% (exchange flux) or 3.0% (gene-backed flux). Per row it is wider. For the selected arm, the gene-backed flux ratio to the shipped arm spans 0.65–1.45 on the tune half (29 / 167 rows outside ±10%) and 0.80–1.27 on the confirm half (12 / 193). The exchange flux ratio spans 1.000–1.12 (1 / 167) and 1.000–1.18 (4 / 193), and it never falls below 1.

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
| est/lin/G0/X1 | 1.004 | 1.016 | 1.029 | yes | 0.1196 | 0.1753 |
| est/lin/Gmeas/X100 | 1.000 | 1.018 | 1.017 | yes | 0.1234 | 0.1777 |
| est/lin/Gmeas/X10 | 1.000 | 1.018 | 1.023 | yes | 0.1223 | 0.1776 |
| est/lin/Gmeas/X1 | 1.004 | 1.016 | 1.029 | yes | 0.1196 | 0.1753 |
| est/lin/G0.2/X100 | 1.000 | 1.018 | 1.023 | yes | 0.1234 | 0.1777 |
| est/lin/G0.2/X10 | 1.000 | 1.018 | 1.023 | yes | 0.1223 | 0.1776 |
| est/lin/G0.2/X1 | 1.004 | 1.016 | 1.029 | yes | 0.1197 | 0.1753 |
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
| est/lin/*/X1 | G0, Gmeas, G0.2 | 2/1/0 | est/lin/G0/X1 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| est/log/*/X100 | G0, Gmeas, G0.2 | 2/4/1 | est/log/G0/X100 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| est/log/*/X10 | G0, Gmeas, G0.2 | 4/4/2 | est/log/G0/X10 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |
| est/log/*/X1 | G0, Gmeas, G0.2 | 7/5/3 | est/log/G0/X1 | lowest pooled A on the G-varying rows; ties within 0.002 -> fewer departures |

**Selection.** Selected: **est/lin/G0/X1** (pooled A 0.1196 against the shipped 0.1362; lowest pooled A 0.1196 at est/lin/G0/X1; tie set within 0.002: est/lin/G0/X1; 2 departure(s) from the shipped configuration). The joint selection over all eligible configurations on all rows (sensitivity) picks est/lin/G0/X1 (agrees).

- Ranking under pooling A (eligible configurations): est/lin/G0/X1 (0.1196), est/lin/Gmeas/X1 (0.1196), est/lin/G0.2/X1 (0.1197), est/lin/G0/X10 (0.1222), est/lin/G0.2/X10 (0.1223), est/lin/Gmeas/X10 (0.1223), est/lin/G0.2/X100 (0.1234), est/lin/Gmeas/X100 (0.1234), est/lin/G0/X100 (0.1234), est/log/G0/X1 (0.1249)
- Ranking under pooling B (reported only): est/lin/Gmeas/X1 (0.1753), est/lin/G0/X1 (0.1753), est/lin/G0.2/X1 (0.1753), est/lin/G0/X10 (0.1775), est/lin/G0.2/X10 (0.1776), est/lin/Gmeas/X10 (0.1776), est/lin/G0.2/X100 (0.1777), est/lin/Gmeas/X100 (0.1777), est/lin/G0/X100 (0.1777), est/log/G0/X10 (0.1789)
- Oracle headroom captured on the tune half: 0.677 of the oracle's headroom at its own cell; 0.908 in the shipped cell's frame (the shipped configuration: 0.259)

## 4. Confirmation (confirm half)

Confirm rows by stratum: {"n2": 54, "n3": 51, "n4+": 58, "n1": 30}. Intervals are 95% percentile cluster bootstraps over holdout groups (2,000 draws, seed 1729 + 1000 × stratum index). The SESOI is 0.01 of D.

| contrast | pooling | face mid | face min | face max | as returned (taf_gb = truth-absent / the arm's own gene-backed flux = 1 − flux precision; not / D) | label |
|---|---|---|---|---|---|---|
| selected-shipped (est/lin/G0/X1 - kn/lin/G0/X100) | A | -0.0172 [-0.0221, -0.0118] | -0.0076 [-0.0133, -0.0047] | -0.0268 [-0.0323, -0.0178] | -0.0149 [-0.0203, -0.0100] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| selected-shipped (est/lin/G0/X1 - kn/lin/G0/X100) | B | -0.0065 [-0.0119, +0.0032] | -0.0014 [-0.0075, +0.0085] | -0.0116 [-0.0170, -0.0014] | -0.0106 [-0.0159, -0.0003] | inconclusive: the midpoint CI crosses 0 and a SESOI bound |
| selected-two_tier (est/lin/G0/X1 - two_tier/lin/G0/X1) | A | -0.0222 [-0.0268, -0.0159] | -0.0100 [-0.0149, -0.0062] | -0.0344 [-0.0403, -0.0244] | -0.0221 [-0.0275, -0.0169] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| selected-two_tier (est/lin/G0/X1 - two_tier/lin/G0/X1) | B | -0.0070 [-0.0124, +0.0025] | -0.0017 [-0.0077, +0.0082] | -0.0123 [-0.0179, -0.0024] | -0.0109 [-0.0156, -0.0016] | inconclusive: the midpoint CI crosses 0 and a SESOI bound |
| selected-shuffled (est/lin/G0/X1 - est_shuffled/lin/G0/X1) | A | -0.0282 [-0.0382, -0.0195] | -0.0279 [-0.0379, -0.0192] | -0.0286 [-0.0386, -0.0198] | -0.0235 [-0.0322, -0.0082] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| selected-shuffled (est/lin/G0/X1 - est_shuffled/lin/G0/X1) | B | -0.0199 [-0.0282, -0.0069] | -0.0210 [-0.0292, -0.0068] | -0.0189 [-0.0272, -0.0070] | -0.0182 [-0.0253, -0.0056] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| oracle-two_tier@selected_cell (oracle/lin/G0/X1 - two_tier/lin/G0/X1) | A | -0.0324 [-0.0372, -0.0242] | -0.0201 [-0.0253, -0.0143] | -0.0446 [-0.0507, -0.0327] | -0.0302 [-0.0356, -0.0224] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| oracle-two_tier@selected_cell (oracle/lin/G0/X1 - two_tier/lin/G0/X1) | B | -0.0200 [-0.0260, -0.0098] | -0.0136 [-0.0193, -0.0050] | -0.0265 [-0.0334, -0.0139] | -0.0207 [-0.0267, -0.0098] | lowers truth-absent flux by >= SESOI: midpoint CI below 0 with mean <= -0.01, and both face endpoints' CIs below 0 |
| shipped-two_tier@shipped_cell (kn/lin/G0/X100 - two_tier/lin/G0/X100) | A | -0.0056 [-0.0074, -0.0036] | -0.0030 [-0.0041, -0.0012] | -0.0083 [-0.0115, -0.0054] | -0.0073 [-0.0110, -0.0048] | detectable but below SESOI: the midpoint CI excludes 0, |mean| < 0.01, and the CIs do not lie inside +-0.01 |
| shipped-two_tier@shipped_cell (kn/lin/G0/X100 - two_tier/lin/G0/X100) | B | -0.0007 [-0.0010, -0.0005] | -0.0004 [-0.0005, -0.0002] | -0.0011 [-0.0015, -0.0007] | -0.0010 [-0.0014, -0.0006] | equivalent within SESOI: the midpoint and both endpoint CIs lie inside +-0.01, nonzero (midpoint CI excludes 0) |
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
| selected-two_tier | n3 | 51 (26) | -0.0249 [-0.0340, -0.0159] | -0.0089 [-0.0200, +0.0000] | -0.0408 [-0.0527, -0.0279] |
| selected-two_tier | n4+ | 58 (21) | -0.0288 [-0.0346, -0.0089] | -0.0086 [-0.0137, +0.0020] | -0.0490 [-0.0589, -0.0180] |
| selected-two_tier | n1 | 30 (20) | -0.0047 [-0.0111, +0.0063] | -0.0004 [-0.0075, +0.0108] | -0.0090 [-0.0155, +0.0022] |
| selected-shuffled | n2 | 54 (26) | -0.0257 [-0.0400, -0.0142] | -0.0256 [-0.0399, -0.0142] | -0.0258 [-0.0401, -0.0143] |
| selected-shuffled | n3 | 51 (26) | -0.0375 [-0.0504, -0.0244] | -0.0360 [-0.0495, -0.0236] | -0.0390 [-0.0511, -0.0254] |
| selected-shuffled | n4+ | 58 (21) | -0.0284 [-0.0351, -0.0103] | -0.0284 [-0.0351, -0.0103] | -0.0284 [-0.0351, -0.0103] |
| selected-shuffled | n1 | 30 (20) | -0.0187 [-0.0281, -0.0035] | -0.0199 [-0.0293, -0.0034] | -0.0175 [-0.0271, -0.0035] |

Guard metrics, selected − shipped:

| contrast | pooling | exchange flux | gene-backed flux | flux precision | gapfilled flux |
|---|---|---|---|---|---|
| selected-shipped | A | +2.5045 [+1.2036, +3.3873] | +1.1801 [-12.9610, +7.9732] | +0.0149 [+0.0100, +0.0203] | -0.2501 [-0.7399, -0.0900] |
| selected-shipped | B | +0.9219 [+0.4818, +1.5527] | +8.8665 [+1.6885, +13.3366] | +0.0106 [+0.0003, +0.0159] | -0.2961 [-1.0766, -0.0431] |

The flux-precision column is the taf_gb column with its sign flipped (taf_gb = 1 − precision), not separate evidence. At n = 1, where prFBA usually runs (86.8% of answered amplicons), selected − shipped is −0.0048 [−0.0113, +0.0063] (30 rows, 20 groups): not detected. On the tune half it was −0.0124 [−0.0273, −0.0030], so it did not replicate, and pooling B is inconclusive (−0.0065 [−0.0119, +0.0032]). The result that meets the SESOI holds under pooling A (n ≥ 2) only. At the best vertex of each face (face min) the gain is −0.0076, below the SESOI. The solver's returned vertex gains −0.0146 [−0.0201, −0.0097] (taf_point_D).

**Where the gain comes from (exploratory; contrasts not declared).** At X100, est − kn is −0.0128 [−0.0186, −0.0088] at the face midpoint. It re-routes flux among gene-backed reactions: truth-supported flux +0.034 D, with exchange, gene-less and gapfilled flux unchanged. It also narrows the TAF face by 0.019 D. With est, X1 − X100 adds −0.0043 [−0.0061, −0.0017] as a uniform shift of the face, by a different route: exchange flux rises (+0.0034 D) and every internal class falls (truth-supported −0.013 D, gene-less −0.002 D, gapfilled −0.0004 D). About a third of that is less gene-backed flux against the fixed D; the rest is a precision gain of +0.0029. The two parts add up (−0.0169 against −0.0172). Over the selected configuration as a whole, the drop at the returned vertex is all selectivity (volume term +0.0002 [−0.0016, +0.0017]), and gene-less and gapfilled flux fall. Neither run records per-exchange fluxes, so which nutrients X1 imports more of is not known.

Oracle headroom captured on the confirm half (pooling A): selected (own cell) 0.685; shipped 0.231; selected in the shipped cell's frame 0.933 (tune: 0.908).

**Degeneracy** (FVA on the optimal face, finalists only, second pass; the change in the number of gene-backed reactions with FVA range > 1e-3, pooling A; arm levels: selected 13.6, shipped 35.4, two_tier at X1 41.5, est_shuffled at X1 2.1. The shuffled arm is the least degenerate and the worst on TAF, so breaking ties does not by itself lower TAF):

- selected-shipped: opt face -21.8900 [-24.6501, -14.8117]; 1e-7 face -22.4329 [-25.0767, -15.9331]
- selected-two_tier: opt face -27.9726 [-30.7292, -20.2375]; 1e-7 face -28.4936 [-31.1012, -21.5362]
- selected-shuffled: opt face +11.4578 [+9.9520, +13.9285]; 1e-7 face +11.4961 [+9.6409, +13.4811]

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
- The confirm half is held out from selection, not from the design. Its 193 amplicons were in the weighted-objective test, whose exploratory exchange-weight-1 arm (identical to kn/lin/G0/X1 on 356 / 360 rows; −0.0090 [−0.0111, −0.0050] against two_tier_ex1 on this half) had been seen when the X factor was declared.
- One holdout group (grp00000, fold 0, 96 source genomes) holds 96 of the 193 confirm rows (43–60% of each stratum). Without it, selected − shipped is −0.0164 [−0.0236, −0.0100] (A); within it, −0.0175. No single-group deletion moves the pooled-A mean outside [−0.0174, −0.0162].
- The guard is on total exchange flux only. Uptake composition under X1 was not recorded, and the selected arm imports more on 57% of rows (up to +18%).
- scripts/objective_tuning.py, scripts/objective_tuning_results.py and OBJECTIVE_TUNING_RESULTS.md are untracked in git. "Fixed before the full run" rests on file times (2026-11-02 23:05 / 23:15; the run started 2026-11-03 08:50) and on the hashes recorded in meta.json and summary.json.
