# Objective tuning (Stage O): results

> Filled by `scripts/objective_tuning_results.py` from `/home/freiburger/Documents/prFBA/region_validation/shared_reference/data/objective_tuning/smoke` on 2026-11-02T23:06:18-0600. Check every number by an independent re-derivation before quoting it.

**PLUMBING / SMOKE OUTPUT: NOT RESULTS.** Weights file is a STUB (meta.stub), not Stage E's estimator. 

Design: `PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md`, section "Stage O" (declared 2026-11-02, post hoc, binding). Harness: `scripts/objective_tuning.py`. Selection and confirmation: `scripts/objective_tuning_results.py`. Outputs: `data/objective_tuning/{per_amplicon.jsonl, meta.json, summary.json}`. Summary generated 2026-11-02T23:06:18-0600.

## 1. What was run

- **Sample.** The weighted-objective test's amplicons (`data/weighted_objective/population.json`): 6 rows, status {"ok": 6}. The gapfilled model of each selection's union EC set was rebuilt as before: ModelSEEDpy c2e31a0, templates 21a8084, reversals applied, complete medium, biomass >= 0.95 x max.
- **Split.** A holdout group goes to the tune half if sha256("1729|" + group_id) is even, else to the confirm half. Rows by half: {"tune": 4, "confirm": 2}.
- **Probability sources.** `kn` is the shipped k/n. `est` is Stage E's estimator at role level, from `/home/freiburger/Documents/prFBA/region_validation/shared_reference/data/objective_tuning/smoke/stage_o_weights.STUB.json` (sha256 a92922436494): E1 hits (stub): K=100, tau=0.01, dedup=genome, first 100 hits, stage 0 = calibrate_pfunc.stage0_ranked; clade = anchor (reps[0]) genus in the fold reference, family if < 3; alpha=0.0. `est_shuffled` is est permuted among the gene-backed reactions, using the prior shuffled arm's permutation, one per amplicon. `two_tier` sets gene-backed p = 1. `oracle` sets gene-backed p = 1 if truth-supported, else 0.
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
| rows ok / total | 6 / 6 |
| rebuilt model = prior row (n_new_rxn, growth 1e-6, n_gene_backed, n_gapfilled) | 6 / 6 |
| D (two_tier/lin/G0/X100 gene-backed flux) bitwise = prior | 6 / 6 |
| kn/lin/G0/X100 objective bitwise = prior `real` | 6 / 6 |
| kn/lin/G0/X100 returned vertex (gene-backed, truth-absent flux) bitwise = prior | 6 / 6, 6 / 6 |
| kn/lin/G0/X100 face min/mid/max within 1e-9 of prior `real` | 5 / 6 (max abs diff mid 3.5806965741769746e-05) |
| two_tier/lin/G0/X100 objective bitwise = prior `two_tier` | 6 / 6 |
| two_tier/lin/G0/X100 face within 1e-9 of prior | 6 / 6 (max abs diff mid 2.1669521732547992e-10) |
| coefficient map of the new function = prfba_released_copy (shipped arm) | 6 / 6 |
| released_check pass; = prior objective_released bitwise | 6 / 6; 6 / 6 |
| determinism pairs identical (same LP) | 72 / 72 |
| prior arms re-run through the prior code, bitwise (repro) | 6 / 6 |
| arm solves optimal | 540 / 540 |
| biomass at 0.95 x max (min, max ratio; off floor > 1e-6) | 0.9499999999999998, 0.95; 0 |
| face LP failures (after one cold retry); retried | 0; 0 |
| gapfilled reactions that are gene-backed | 0 |
| rows with alpha > 0 and a null clade (literal formula) | 0 |
| wall mean / max per amplicon; largest worker RSS | 118.5 s / 189.2 s; 2.786 GiB |

## 3. Selection (tune half)

Tune rows by stratum: {"n2": 1, "n3": 1, "n4+": 1, "n1": 1}.

**Eligibility.** A configuration is eligible if its median total exchange flux and its median gene-backed flux are each within 10% of the shipped arm's. 36 of 36 configurations are eligible. By exchange weight: {"X100": 12, "X10": 12, "X1": 12}.

| configuration | exchange flux ratio | gene-backed flux ratio | precision ratio | eligible | pooled A | pooled B |
|---|---|---|---|---|---|---|
| kn/lin/G0/X100 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/lin/G0/X10 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/lin/G0/X1 | 1.000 | 0.998 | 0.998 | yes | 0.2523 | 0.5468 |
| kn/lin/Gmeas/X100 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/lin/Gmeas/X10 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/lin/Gmeas/X1 | 1.000 | 0.998 | 0.998 | yes | 0.2523 | 0.5468 |
| kn/lin/G0.2/X100 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/lin/G0.2/X10 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/lin/G0.2/X1 | 1.000 | 0.998 | 0.998 | yes | 0.2523 | 0.5468 |
| kn/log/G0/X100 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/log/G0/X10 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/log/G0/X1 | 1.000 | 0.998 | 1.000 | yes | 0.2523 | 0.5468 |
| kn/log/Gmeas/X100 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/log/Gmeas/X10 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/log/Gmeas/X1 | 1.000 | 0.998 | 1.000 | yes | 0.2523 | 0.5468 |
| kn/log/G0.2/X100 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/log/G0.2/X10 | 1.000 | 1.000 | 1.000 | yes | 0.2523 | 0.5465 |
| kn/log/G0.2/X1 | 1.000 | 0.998 | 1.000 | yes | 0.2523 | 0.5468 |
| est/lin/G0/X100 | 1.000 | 1.002 | 1.049 | yes | 0.2263 | 0.5430 |
| est/lin/G0/X10 | 1.000 | 1.002 | 1.049 | yes | 0.2263 | 0.5430 |
| est/lin/G0/X1 | 1.000 | 1.000 | 1.049 | yes | 0.2263 | 0.5433 |
| est/lin/Gmeas/X100 | 1.000 | 1.002 | 1.049 | yes | 0.2263 | 0.5430 |
| est/lin/Gmeas/X10 | 1.000 | 1.002 | 1.049 | yes | 0.2263 | 0.5430 |
| est/lin/Gmeas/X1 | 1.000 | 1.000 | 1.049 | yes | 0.2263 | 0.5433 |
| est/lin/G0.2/X100 | 1.000 | 1.002 | 1.049 | yes | 0.2263 | 0.5430 |
| est/lin/G0.2/X10 | 1.000 | 1.002 | 1.049 | yes | 0.2263 | 0.5430 |
| est/lin/G0.2/X1 | 1.000 | 1.000 | 1.049 | yes | 0.2263 | 0.5433 |
| est/log/G0/X100 | 1.000 | 1.000 | 1.033 | yes | 0.2355 | 0.5443 |
| est/log/G0/X10 | 1.000 | 1.000 | 1.033 | yes | 0.2355 | 0.5443 |
| est/log/G0/X1 | 1.000 | 0.998 | 1.033 | yes | 0.2355 | 0.5446 |
| est/log/Gmeas/X100 | 1.000 | 1.000 | 1.033 | yes | 0.2355 | 0.5443 |
| est/log/Gmeas/X10 | 1.000 | 1.000 | 1.033 | yes | 0.2355 | 0.5443 |
| est/log/Gmeas/X1 | 1.000 | 0.998 | 1.033 | yes | 0.2355 | 0.5446 |
| est/log/G0.2/X100 | 1.000 | 1.000 | 1.033 | yes | 0.2355 | 0.5443 |
| est/log/G0.2/X10 | 1.000 | 1.000 | 1.033 | yes | 0.2355 | 0.5443 |
| est/log/G0.2/X1 | 1.000 | 0.998 | 1.033 | yes | 0.2355 | 0.5446 |

**The G factor.** G is read only on rows where the gapfilled flux varies across G's levels:

| family | eligible G | rows where gapfilled flux varies (n2/n3/n4+) | choice | why |
|---|---|---|---|---|
| kn/lin/*/X100 | G0, Gmeas, G0.2 | 0/0/0 | kn/lin/G0/X100 | gapfilled flux varies on no row: G not read |
| kn/lin/*/X10 | G0, Gmeas, G0.2 | 0/0/0 | kn/lin/G0/X10 | gapfilled flux varies on no row: G not read |
| kn/lin/*/X1 | G0, Gmeas, G0.2 | 0/0/0 | kn/lin/G0/X1 | gapfilled flux varies on no row: G not read |
| kn/log/*/X100 | G0, Gmeas, G0.2 | 0/0/0 | kn/log/G0/X100 | gapfilled flux varies on no row: G not read |
| kn/log/*/X10 | G0, Gmeas, G0.2 | 0/0/0 | kn/log/G0/X10 | gapfilled flux varies on no row: G not read |
| kn/log/*/X1 | G0, Gmeas, G0.2 | 0/0/0 | kn/log/G0/X1 | gapfilled flux varies on no row: G not read |
| est/lin/*/X100 | G0, Gmeas, G0.2 | 0/0/0 | est/lin/G0/X100 | gapfilled flux varies on no row: G not read |
| est/lin/*/X10 | G0, Gmeas, G0.2 | 0/0/0 | est/lin/G0/X10 | gapfilled flux varies on no row: G not read |
| est/lin/*/X1 | G0, Gmeas, G0.2 | 0/0/0 | est/lin/G0/X1 | gapfilled flux varies on no row: G not read |
| est/log/*/X100 | G0, Gmeas, G0.2 | 0/0/0 | est/log/G0/X100 | gapfilled flux varies on no row: G not read |
| est/log/*/X10 | G0, Gmeas, G0.2 | 0/0/0 | est/log/G0/X10 | gapfilled flux varies on no row: G not read |
| est/log/*/X1 | G0, Gmeas, G0.2 | 0/0/0 | est/log/G0/X1 | gapfilled flux varies on no row: G not read |

**Selection.** Selected: **est/lin/G0/X100** (pooled A 0.2263 against the shipped 0.2523; lowest pooled A 0.2263 at est/lin/G0/X1; tie set within 0.002: est/lin/G0/X100, est/lin/G0/X10, est/lin/G0/X1; 1 departure(s) from the shipped configuration). The joint selection over all eligible configurations on all rows (sensitivity) picks est/lin/G0/X100 (agrees).

- Ranking under pooling A (eligible configurations): est/lin/G0.2/X100 (0.2263), est/lin/Gmeas/X1 (0.2263), est/lin/G0.2/X1 (0.2263), est/lin/G0.2/X10 (0.2263), est/lin/G0/X1 (0.2263), est/lin/Gmeas/X10 (0.2263), est/lin/G0/X10 (0.2263), est/lin/Gmeas/X100 (0.2263), est/lin/G0/X100 (0.2263), est/log/G0/X1 (0.2355)
- Ranking under pooling B (reported only): est/lin/G0/X10 (0.5430), est/lin/Gmeas/X10 (0.5430), est/lin/G0.2/X10 (0.5430), est/lin/G0.2/X100 (0.5430), est/lin/Gmeas/X100 (0.5430), est/lin/G0/X100 (0.5430), est/lin/G0.2/X1 (0.5433), est/lin/Gmeas/X1 (0.5433), est/lin/G0/X1 (0.5433), est/log/G0/X10 (0.5443)
- Oracle headroom captured on the tune half: 0.966 of the oracle's headroom at its own cell; 0.966 in the shipped cell's frame (the shipped configuration: 0.110)

## 4. Confirmation (confirm half)

Confirm rows by stratum: {"n2": 1, "n3": 0, "n4+": 0, "n1": 1}. Intervals are 95% percentile cluster bootstraps over holdout groups (2,000 draws, seed 1729 + 1000 × stratum index). The SESOI is 0.01 of D.

| contrast | pooling | face mid | face min | face max | as returned (taf_gb) | label |
|---|---|---|---|---|---|---|
| selected-shipped (est/lin/G0/X100 - kn/lin/G0/X100) | A | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| selected-shipped (est/lin/G0/X100 - kn/lin/G0/X100) | B | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| selected-two_tier (est/lin/G0/X100 - two_tier/lin/G0/X100) | A | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| selected-two_tier (est/lin/G0/X100 - two_tier/lin/G0/X100) | B | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| selected-shuffled (est/lin/G0/X100 - est_shuffled/lin/G0/X100) | A | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| selected-shuffled (est/lin/G0/X100 - est_shuffled/lin/G0/X100) | B | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| oracle-two_tier@selected_cell (oracle/lin/G0/X100 - two_tier/lin/G0/X100) | A | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| oracle-two_tier@selected_cell (oracle/lin/G0/X100 - two_tier/lin/G0/X100) | B | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| shipped-two_tier@shipped_cell (kn/lin/G0/X100 - two_tier/lin/G0/X100) | A | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| shipped-two_tier@shipped_cell (kn/lin/G0/X100 - two_tier/lin/G0/X100) | B | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| oracle-two_tier@shipped_cell (oracle/lin/G0/X100 - two_tier/lin/G0/X100) | A | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |
| oracle-two_tier@shipped_cell (oracle/lin/G0/X100 - two_tier/lin/G0/X100) | B | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | not evaluable |

Shuffled comparator: P = est: est_shuffled is its own shuffled arm.

Per stratum:

| contrast | stratum | n (groups) | face mid | face min | face max |
|---|---|---|---|---|---|
| selected-shipped | n2 | 1 (1) | -0.0156 [-0.0156, -0.0156] | -0.0103 [-0.0103, -0.0103] | -0.0210 [-0.0210, -0.0210] |
| selected-shipped | n1 | 1 (1) | +0.0001 [+0.0001, +0.0001] | +0.0002 [+0.0002, +0.0002] | +0.0000 [+0.0000, +0.0000] |
| selected-two_tier | n2 | 1 (1) | -0.0302 [-0.0302, -0.0302] | -0.0242 [-0.0242, -0.0242] | -0.0362 [-0.0362, -0.0362] |
| selected-two_tier | n1 | 1 (1) | +0.0001 [+0.0001, +0.0001] | +0.0002 [+0.0002, +0.0002] | +0.0000 [+0.0000, +0.0000] |
| selected-shuffled | n2 | 1 (1) | -0.0242 [-0.0242, -0.0242] | -0.0242 [-0.0242, -0.0242] | -0.0242 [-0.0242, -0.0242] |
| selected-shuffled | n1 | 1 (1) | +0.0002 [+0.0002, +0.0002] | +0.0002 [+0.0002, +0.0002] | +0.0002 [+0.0002, +0.0002] |

Guard metrics, selected − shipped:

| contrast | pooling | exchange flux | gene-backed flux | flux precision | gapfilled flux |
|---|---|---|---|---|---|
| selected-shipped | A | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ |
| selected-shipped | B | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ | no rows in stratum n3,n4+ |

Oracle headroom captured on the confirm half (pooling A): selected (own cell) n/a; shipped n/a.

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
