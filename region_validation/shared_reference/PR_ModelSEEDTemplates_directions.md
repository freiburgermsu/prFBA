# Fix the directions of six reactions that let gapfilled models make ATP from nothing

Branch `fix-irreversible-reaction-directions` (commit `fcfe00d`), against `main` (`2294299`).

## Summary

Six template reactions have directions that contradict the ModelSEED database and the
thermodynamics. In gapfilled models, each closes an energy-generating cycle: the model
makes ATP from ADP + Pi with no net exchange, so it produces ATP on a medium of water and
protons alone. This PR sets each reaction to its database direction in the v7.0 templates
(GramNeg, GramPos, Archaea, and Core for catalase). For each reaction, `direction`,
`GapfillDirection`, `lower_bound`, `upper_bound`, `maxforflux` and `maxrevflux` change
together. Nothing else in the files changes; this was checked by parsing every file
before and after.

| Reaction | Name | Before (direction / GapfillDirection) | After | Database | dG'0, forward (kcal/mol) |
|---|---|---|---|---|---|
| rxn00006 | catalase (2 H2O2 => 2 H2O + O2) | `<` / `<` (all four templates, Core included) | `>` | `>` | -45.2 (GOLD) |
| rxn05119 | L-aspartate:NAD+ oxidoreductase (deaminating) | `>` / `=` | `<` | `<` | +13.3 |
| rxn05117 | L-aspartate:NADP+ oxidoreductase (deaminating) | `>` / `=` | `<` | `<` | NADP analogue of rxn05119 |
| rxn10966 | nitrate reductase (thioredoxin), NTRARz | `=` / `=`, reverse_penalty -2.94 | `>` | `?` (coupled rxn00572 is irreversible, `<`) | reversed with thioredoxin reductase: +34.8 |
| rxn00381 | adenylylsulfate sulfohydrolase | `=` / `=`, reverse_penalty -1.08 | `>` | `>` | -19.6 (GOLD) |
| rxn00084 | NADPH:hydrogen-peroxide oxidoreductase | `=` / `=`, reverse_penalty -6.71 | `>` | `>` | -76.5 (SILVER) |

dG'0 values are from the ModelSEEDDatabase thermodynamics grades (eQuilibrator where
available).

## How the errors show up

We found these while debugging a validity gate in a benchmark that builds ModelSEED models
from predicted EC sets (MSATPCorrection + MSGapfill on complete medium, as in
`MSBuilder.build_and_gapfill`). The gate requires no ATP or biomass on water and protons.
We rebuilt 218 gapfilled models from 20 genomes, reproducing the original runs' gapfill
solutions and gate verdicts exactly. 44 of them make ATP from nothing (ATPM flux 150-1000, no biomass), and in every
one the leak is a closed internal loop. In each model a single reaction added by
MSGapfill closes the loop:

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

The template's negative `reverse_penalty` decided the direction MSGapfill chose in 29 of
the 44 models. The gapfilling objective costs a direction at 2 x |base_cost + penalty|, so
a negative reverse penalty makes the reverse direction the cheaper one.

## Validation

Rebuilt all 218 models from the patched template files, with the gapfilling procedure
otherwise unchanged:

- All 44 leaking models now pass. No model that passed before starts leaking.
- The other models' gapfill solutions and growth on complete medium are unchanged in
  practice: the median growth ratio is 1.000. The gapfilling MILP is degenerate, so a few
  solutions move to equally cheap alternatives.
- Under stock ModelSEEDpy, 43 of the 44 grow. The remaining one, and 7 models that already
  failed to grow, need a separate ModelSEEDpy fix (below). With that fix applied, all 218
  models grow and none leaks.

We also tested two alternatives:
- Removing the reactions outright cleared the leaks, but left 13 models unable to grow.
- Setting every database-irreversible reaction to its database direction (more than 1,000
  reactions) cleared the leaks, but stopped 11 models from growing and lowered median
  growth by 25%.

So this PR is limited to the six reactions with direct evidence.

## Not in this PR (for the maintainers to decide)

- **Older template versions.** The v5.0 and v6.0 templates carry the same six entries. They
  are left unchanged here so that published models stay reproducible.
- **Broader direction mismatches.** On `GapfillDirection`, 833 template reactions are
  reversible where the database says `>`, and 336 where it says `<` (Archaea: 835 and 338).
- **Negative penalties.** 1,594 reactions carry a negative `reverse_penalty`. Because of the
  `abs()` in `GapfillingPkg.build_gapfilling_objective_function`, this makes the reverse
  direction cheaper than the forward one for 1,190 of them.
- **ModelSEEDpy drops reversals.** `MSBuilder.build_and_gapfill` integrates only the
  gapfiller's `new` reactions and drops its `reversed` set. In some models this removes
  the growth the gapfill solution was built for. This is a separate ModelSEEDpy issue.

## Reproduce

Analysis code and per-model results are in freiburgermsu/prFBA (branch
`hardware-aware-alignment-dispatch`):
- `region_validation/shared_reference/scripts/egc_rebuild.py`, including `--template-dir`;
- `region_validation/shared_reference/scripts/egc_analyze.py`;
- `region_validation/shared_reference/EGC_C12_REPORT.md`;
- `region_validation/shared_reference/EGC_C12_gated_models.csv`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
