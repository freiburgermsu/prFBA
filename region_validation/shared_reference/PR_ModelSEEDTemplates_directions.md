# Fix the directions of five reactions that let gapfilled models make ATP from nothing

Branch `fix-irreversible-reaction-directions` (commit `21a8084`), against `main` (`2294299`).

## Summary

Each of these five template reactions is a step in an energy-generating cycle in
gapfilled models. The model makes ATP from ADP + Pi with no net exchange, so it produces
ATP on a medium of water and protons alone.

This PR restricts each reaction to its physiological direction in the v7.0 templates
(GramNeg, GramPos, Archaea, and Core for catalase). For each one, `direction` and
`GapfillDirection` are set, and `lower_bound`, `upper_bound`, `maxforflux` and
`maxrevflux` follow `direction` where they differ. Nothing else in the files changes;
this was checked by parsing every file before and after.

| Reaction | Name | Before (direction / GapfillDirection) | After | Database | dG'0, forward (kcal/mol) |
|---|---|---|---|---|---|
| rxn00006 | catalase (2 H2O2 => 2 H2O + O2) | `<` / `<` (all four templates, Core included) | `>` | `>` | -45.2 (GOLD) |
| rxn10966 | nitrate reductase (thioredoxin), NTRARz | `=` / `=`, reverse_penalty -2.94 | `>` | `?` | reversed with thioredoxin reductase: +34.8 |
| rxn00381 | adenylylsulfate sulfohydrolase | `=` / `=`, reverse_penalty -1.08 | `>` | `>` | -19.6 (GOLD) |
| rxn00084 | NADPH:hydrogen-peroxide oxidoreductase | `=` / `=`, reverse_penalty -6.71 | `>` | `>` | -76.5 (SILVER) |
| rxn08167 | L-aspartate oxidase (menaquinone-8) | `=` / `=`, reverse_penalty -3.32 | `>` | `<` | see below |

**Three follow the database directly** (rxn00006, rxn00381, rxn00084). rxn10966 has no
database direction. Its direction follows from the database's `<` for rxn00572, the net
reaction of rxn10966 and thioredoxin reductase.

**rxn08167 is set against the database.** Aspartate oxidation is L-aspartate oxidase's
(NadB's) physiological direction. The database's `<` rests on a BRONZE eQuilibrator
value (+12.6 kcal/mol). That value is inconsistent with the database's own value for the
NAD-linked L-aspartate dehydrogenase (NadX, rxn05119, +13.3): the two reactions differ
only in the electron acceptor, and menaquinone is about 11 kcal/mol the better one.
- **Thermodynamics alone does not settle it.** A consistent estimate puts rxn08167 near
  equilibrium (about +2 kcal/mol). Group contribution and dGPredictor make no call
  (-7.1 ± 13.1 and -17.1 ± 51.5, both `?` in the database).
- **The direction rests on physiology.** NadB oxidises aspartate, pulled by quinolinate
  synthase (NadA).

dG'0 values are from the ModelSEEDDatabase thermodynamics grades (eQuilibrator where
available). The Gram templates are single-line JSON, so the diff is easiest to read with
`git diff --word-diff main...fix-irreversible-reaction-directions`.

## How the errors show up

We found these while debugging a validity gate in a benchmark that builds ModelSEED models
from predicted EC sets: MSATPCorrection, then MSGapfill on complete medium. The gate
requires no ATP or biomass on water and protons.

We rebuilt 218 gapfilled models from 20 genomes, reproducing the original runs' gapfill
solutions and gate verdicts exactly. 44 of them make ATP from nothing (ATPM flux
150-1000, no biomass). In every one the leak is a closed internal loop, and a single
reaction added by MSGapfill closes it:

- **rxn10966 added in reverse (14 models).** With NADPH:thioredoxin reductase run backwards,
  it is NADP-dependent nitrite oxidation. The respiratory chain turns the resulting NADPH
  into a proton gradient and ATP.
- **rxn09180 (anaerobic protoporphyrinogen oxidase), 15 models.** The loop runs through
  catalase, which every template allows only as 2 H2O + O2 -> 2 H2O2, together with the
  reversible O2-dependent protoporphyrinogen oxidase rxn20627.
- **rxn08167 added in reverse (8 models).** Iminoaspartate + menaquinol -> L-aspartate +
  menaquinone. Together with NadX (rxn05119 or rxn05117) running forward, this makes
  NAD(P)H from menaquinol, about 11 kcal/mol uphill.
- **rxn00381 added in reverse (6 models)** and **rxn00084 added in reverse (1 model)**: each
  reaction is itself the infeasible step.

In 29 of the 44 leaking models, one of rxn10966, rxn08167, rxn00381 and rxn00084 was
added in reverse. In each, the template's negative `reverse_penalty` had made that the
cheaper direction. The steering is wider than the leaks: rxn08167 alone is added in
reverse in 72 of the 218 models. The gapfilling cost of a direction is (2 - gene score) x |base_cost + penalty|,
so the reverse costs 3.4-6.1 against 8-12 forward. Two experiments separate the penalty
from the direction:
- **Blocking the reverse directions.** In those 29 models the gapfiller then added none of
  the four, in either direction, and all 29 pass.
- **Negating the negative reverse penalties alone, on the original templates.** This
  prevents 27 of the 29 leaks. In the other 2 (Archaea, rxn00381) the reverse is still
  chosen at the higher cost. Of the 142 models we rebuilt that way, it also opens one new
  leak in a model that passed before.

## Validation

We rebuilt all 218 models from the patched template files, with the gapfilling procedure
otherwise unchanged. Every row applies the gapfiller's reversals, as
`MSGapfill.integrate_gapfill_solution` does (opened to ±1000; the rebuild does not
prune):

| Templates | Leaking models | New leaks | Models that do not grow | Passing models' growth, median ratio | Passing models moved > 10% (of 174) |
|---|---|---|---|---|---|
| original v7.0 (reference) | 44 | — | 0 | — | — |
| **this PR** | **0** | **0** | **0** | **1.000** | **18** (0.56x to 4.5x) |

**The reference row is a fresh rebuild.** It uses the original templates, with reversals
applied. It matches the original runs' gate verdicts in 218 of 218 models, and their
gapfill solutions in 211. In the other 7, the ATP-correction expansion differs. Our
ModelSEEDpy fork's biochemistry and thermodynamics modules changed between the original
runs and the rebuild, which probably explains it.

**The gapfill solutions shift.** 6 of 218 are identical to the reference (8 against the
original runs; median Jaccard 0.92). The changed directions alter the gapfilling problem
in every model, and the gapfilling MILP has many near-equal optima.

**Our pipeline's ModelSEEDpy fork drops the gapfiller's reversals.** It integrates only the
gapfiller's new reactions, through the fork's `MSBuilder.build_and_gapfill`. That
explains the 10 models that did not grow in our original runs. With that fork, this PR
still closes all 44 leaks, but 3 models lose growth to dropped reversals.
`MSGapfill.integrate_gapfill_solution` applies reversals, so this does not affect this
PR.

**The alternatives fall short.**
- **NadX instead of rxn08167.** This is the same change, but with rxn05119 and rxn05117 set
  to the database direction `<` and rxn08167 left `=` (branch
  `six-reaction-fix-nadx-variant`, `7d104ac`). It also closes every leak, with nearly
  the same growth: 4 models differ by more than 10%. But it blocks NadX's physiological
  forward reaction, the first step of NAD synthesis in organisms without NadB.
- **Removing the gapfilled culprits.** This also closes the leaks, but leaves the direction
  errors in place for every other use of the templates.
- **Every database direction.** Setting every database-irreversible reaction to its
  database direction (over 1,000 reactions, with rxn10966 forward as here) also closes
  them. But 2 count as closed only because gapfilling becomes infeasible, 11 models stop
  growing, and median growth falls by a quarter. See the companion direction audit.

So this PR is limited to the five reactions with direct evidence.

## Not in this PR (for the maintainers to decide)

- **Older template versions.** The v5.0 and v6.0 templates carry the same five entries.
  They are left unchanged here so that published models stay reproducible.
- **Other direction mismatches.**
  - **The companion PR** (`fix-opposite-template-directions`) fixes five more reactions
    whose template direction fields oppose the database. Four are allowed only in the
    forbidden direction. For rxn11937, `GapfillDirection` contradicts its draft
    direction.
  - **The separate audit issue** compares all the template directions with the database.
    For most of the 41 reactions the templates allow only backwards (in at least one
    field), it finds that the database direction is the wrong one. The other ~1,150
    reactions the templates leave reversible were not reviewed one by one.
- **Negative reverse penalties.** A companion PR (`fix-reverse-penalty-sign`) flips the
  sign of the negative reverse penalties on reactions whose forward direction is
  favoured.

## Reproduce

Analysis code and per-model results are in freiburgermsu/prFBA (branch
`hardware-aware-alignment-dispatch`), under `region_validation/shared_reference/`:
- `scripts/egc_rebuild.py`: rebuild with `--template-dir <checkout>/templates/v7.0
  --variant reversed_only` (or `none` to drop reversals);
- `scripts/template_validation_summary.py` writes
  `data/template_direction_audit/validation_summary.json` and `validation_per_model.csv`.
  The relevant entries are:
  - `base+rev`, `PR1v2+rev` and `PR1v2`;
  - `PR1+rev`, the NadX alternative; its `vs_five_reaction_fix` field lists the 4 models;
  - `mainPEN`;
- `scripts/egc_analyze.py`;
- `EGC_C12_REPORT.md` and `EGC_C12_gated_models.csv`.

The rebuild inputs (the benchmark's EC sets and the model files) are not in the
repository.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
