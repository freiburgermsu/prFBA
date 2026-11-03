# Objective tuning (Stage O): results

<!-- TEMPLATE-BANNER-START -->
> **TEMPLATE — no results yet.** This file is the draft that the full Stage O run fills. `scripts/objective_tuning_results.py` reads it and writes the filled copy to `data/objective_tuning/OBJECTIVE_TUNING_RESULTS.filled.md`. The `{{...}}` fields below come from `data/objective_tuning/summary.json`. This file is never overwritten. The prose around the fields is a draft. Before any number is quoted, check it by an independent re-derivation (design, "Order of work", step 4).
<!-- TEMPLATE-BANNER-END -->

{{status_banner}}

Design: `PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md`, section "Stage O" (declared 2026-11-02, post hoc, binding). Harness: `scripts/objective_tuning.py`. Selection and confirmation: `scripts/objective_tuning_results.py`. Outputs: `data/objective_tuning/{per_amplicon.jsonl, meta.json, summary.json}`. Summary generated {{generated}}.

## 1. What was run

- **Sample.** The weighted-objective test's amplicons (`data/weighted_objective/population.json`): {{n_rows}} rows, status {{status_counts}}. The gapfilled model of each selection's union EC set was rebuilt as before: ModelSEEDpy c2e31a0, templates 21a8084, reversals applied, complete medium, biomass >= 0.95 x max.
- **Split.** A holdout group goes to the tune half if sha256("1729|" + group_id) is even, else to the confirm half. Rows by half: {{halves}}.
- **Probability sources.** `kn` is the shipped k/n. `est` is Stage E's estimator at role level, from `{{weights_path}}` (sha256 {{weights_sha256}}): {{estimator}}. `est_shuffled` is est permuted among the gene-backed reactions, using the prior shuffled arm's permutation, one per amplicon. `two_tier` sets gene-backed p = 1. `oracle` sets gene-backed p = 1 if truth-supported, else 0.
- **Factors.**
  - C, the cost on rxn reactions: `lin` = max(0.01, 1 − p) (shipped) or `log` = max(0.01, ln max(p, 0.01) / ln 0.01).
  - G, the p of gapfilled reactions, in every arm: 0 (shipped), g_meas = {{g_meas}} ({{g_meas_info}}), or 0.2.
  - X, the exchange weight: 100 (shipped), 10 or 1.
  - In total: 5 P × 2 C × 3 G × 3 X = 90 arms per amplicon, each on a fresh copy of the model.
- **Measure.** The total truth-absent gene-backed flux on the optimal face (objective <= opt + 1e-6), as its minimum, midpoint and maximum, each divided by D. D is the gene-backed flux of the two-tier, linear, ex-100 arm, one constant per amplicon. The guard metrics are total exchange flux, gene-backed flux, flux precision and gapfilled flux.
- **Pooling.** A (primary): the n >= 2 strata weighted {{weights_A}}. B (reported, not used to select): all answered amplicons, {{weights_B}}.

## 2. Checks

{{checks_table}}

## 3. Selection (tune half)

Tune rows by stratum: {{tune_n}}.

**Eligibility.** A configuration is eligible if its median total exchange flux and its median gene-backed flux are each within 10% of the shipped arm's. {{n_eligible}} of 36 configurations are eligible. By exchange weight: {{eligible_by_X}}.

{{guard_table}}

**The G factor.** G is read only on rows where the gapfilled flux varies across G's levels:

{{g_reading_table}}

**Selection.** {{selection_paragraph}}

- Ranking under pooling A (eligible configurations): {{ranking_A}}
- Ranking under pooling B (reported only): {{ranking_B}}
- Oracle headroom captured on the tune half: {{headroom_tune}}

## 4. Confirmation (confirm half)

Confirm rows by stratum: {{confirm_n}}. Intervals are 95% percentile cluster bootstraps over holdout groups (2,000 draws, seed 1729 + 1000 × stratum index). The SESOI is 0.01 of D.

{{confirm_table}}

Shuffled comparator: {{shuffled_note}}

Per stratum:

{{confirm_strata_table}}

Guard metrics, selected − shipped:

{{confirm_guard_table}}

Oracle headroom captured on the confirm half (pooling A): {{headroom_confirm}}.

**Degeneracy** (FVA on the optimal face, finalists only, second pass):

{{degeneracy}}

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
