# Fix the directions of six reactions that let gapfilled models make ATP from nothing

Branch `fix-irreversible-reaction-directions` (commit `7d104ac`), against `main` (`2294299`).

## Summary

Each of these six template reactions closes an energy-generating cycle in gapfilled
models: the model makes ATP from ADP + Pi with no net exchange, so it produces ATP on a
medium of water and protons alone. Five contradict the ModelSEED database's direction.
The sixth, rxn10966, has no database direction; it is fixed through the reaction it
couples to.

This PR sets each reaction's direction in the v7.0 templates (GramNeg, GramPos, Archaea,
and Core for catalase). For each reaction, `direction` and `GapfillDirection` are set,
and `lower_bound`, `upper_bound`, `maxforflux` and `maxrevflux` follow `direction` where
they differ. Nothing else in the files changes; this was checked by parsing every file
before and after.

| Reaction | Name | Before (direction / GapfillDirection) | After | Database | dG'0, forward (kcal/mol) |
|---|---|---|---|---|---|
| rxn00006 | catalase (2 H2O2 => 2 H2O + O2) | `<` / `<` (all four templates, Core included) | `>` | `>` | -45.2 (GOLD) |
| rxn05119 | L-aspartate:NAD+ oxidoreductase (deaminating) | `>` / `=` | `<` | `<` | +13.3 |
| rxn05117 | L-aspartate:NADP+ oxidoreductase (deaminating) | `>` / `=` | `<` | `<` | NADP analogue of rxn05119 |
| rxn10966 | nitrate reductase (thioredoxin), NTRARz | `=` / `=`, reverse_penalty -2.94 | `>` | `?` (coupled rxn00572 is `<`) | reversed with thioredoxin reductase: +34.8 |
| rxn00381 | adenylylsulfate sulfohydrolase | `=` / `=`, reverse_penalty -1.08 | `>` | `>` | -19.6 (GOLD) |
| rxn00084 | NADPH:hydrogen-peroxide oxidoreductase | `=` / `=`, reverse_penalty -6.71 | `>` | `>` | -76.5 (SILVER) |

dG'0 values are from the ModelSEEDDatabase thermodynamics grades (eQuilibrator where
available). The Gram templates are single-line JSON, so the diff is easiest to read with
`git diff --word-diff main...fix-irreversible-reaction-directions`.

## How the errors show up

We found these while debugging a validity gate in a benchmark that builds ModelSEED models
from predicted EC sets: MSATPCorrection, then MSGapfill on complete medium. The gate
requires no ATP or biomass on water and protons.

We rebuilt 218 gapfilled models from 20 genomes, reproducing the original runs' gapfill
solutions and gate verdicts exactly. 44 of them make ATP from nothing (ATPM flux
150-1000, no biomass). In every one the leak is a closed internal loop, and in each model
a single reaction added by MSGapfill closes it:

- **rxn10966 added in reverse (14 models).** With NADPH:thioredoxin reductase run backwards,
  it is NADP-dependent nitrite oxidation. The respiratory chain turns the resulting NADPH
  into a proton gradient and ATP.
- **rxn09180 (anaerobic protoporphyrinogen oxidase), 15 models.** The loop runs through
  catalase, which every template allows only as 2 H2O + O2 -> 2 H2O2, together with the
  reversible O2-dependent protoporphyrinogen oxidase rxn20627.
- **rxn08167 added in reverse (8 models).** The uphill step is L-aspartate dehydrogenase
  (rxn05119, or rxn05117 when rxn05119 is blocked) running forward, against the database.
- **rxn00381 added in reverse (6 models)** and **rxn00084 added in reverse (1 model)**: each
  reaction is itself the infeasible step.

In the 29 models where one of the four was added in reverse, the template's negative
`reverse_penalty` had made that the cheapest direction. The gapfilling cost of a direction
is (2 - gene score) x |base_cost + penalty|, so the reverse costs 3.4-6.1 against 8-12
forward. When we blocked their reverse directions, the gapfiller added none of the
four forward. Negating the negative reverse penalties alone, on the original
templates, closes 27 of the 29 leaks. In the other 2 (Archaea, rxn00381) the reverse is
still chosen at the higher cost.

## Validation

We rebuilt all 218 models from the patched template files, with the gapfilling procedure
otherwise unchanged and the gapfiller's reversals integrated, as
`MSGapfill.integrate_gapfill_solution` does:

- All 44 leaking models now pass. No model that passed before starts leaking.
- Every model grows on complete medium, including 10 that did not grow before.
- The passing models' growth is unchanged at the median (ratio 1.000), but 13 of 164
  moved by more than 10% (0.65x to 4.5x).
- Gapfill solutions differ from the original in 190 of the 218 models (median Jaccard
  0.94), mostly because catalase's draft direction changes in 172 of the 218 drafts.

Our pipeline's ModelSEEDpy fork integrates only the gapfiller's new reactions, not its
reversals. With that fork, 43 of the 44 grow, and the remaining one and 7 models that
already failed to grow lack only a dropped reversal. `integrate_gapfill_solution`
applies reversals, so this does not affect this PR.

Removing the gapfilled culprit reactions instead of fixing directions also closes the
leaks, but leaves the direction errors in place for every other use of the templates.
Setting every database-irreversible reaction to its database direction (over 1,000
reactions) also closes them, at a large cost in growth; see the companion direction
audit. So this PR is limited to the six reactions with direct evidence.

## Not in this PR (for the maintainers to decide)

- **Older template versions.** The v5.0 and v6.0 templates carry the same six entries. They
  are left unchanged here so that published models stay reproducible.
- **Other direction mismatches.** Five more reactions that the templates allow only
  backwards are in a companion PR (`fix-opposite-template-directions`). An audit of all
  the template directions against the database is in a separate issue. It finds that,
  for most mismatches, the database direction is the wrong one.
- **Negative reverse penalties.** A companion PR (`fix-reverse-penalty-sign`) flips the
  sign of the negative reverse penalties on reactions whose forward direction is
  favoured.

## Reproduce

Analysis code and per-model results are in freiburgermsu/prFBA (branch
`hardware-aware-alignment-dispatch`), under `region_validation/shared_reference/`:
- `scripts/egc_rebuild.py`: rebuild with `--template-dir <checkout>/templates/v7.0
  --variant reversed_only` (or `none` to drop reversals);
- `scripts/template_validation_summary.py`;
- `data/template_direction_audit/`: `validation_summary.json` and
  `validation_per_model.csv`, rows `baseline`, `PR1`, `PR1+rev`, `mainPEN`;
- `scripts/egc_analyze.py`;
- `EGC_C12_REPORT.md` and `EGC_C12_gated_models.csv`.

The full model files are not in the repository.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
