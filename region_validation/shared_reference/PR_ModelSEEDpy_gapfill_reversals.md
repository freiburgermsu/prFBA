## Summary

`MSGapfill.run_gapfilling` returns two kinds of change:
- `"new"`: template reactions to add;
- `"reversed"`: reactions the model already has, which the solution needs in a direction
  their bounds forbid.

`MSBuilder.build_and_gapfill` built its solution from `"new"` alone and integrated it with
the static `MSBuilder.integrate_gapfill_solution`, which only adds reactions. Every
reversal was dropped. The returned model lacks part of its own gapfill solution, and can
fail to grow on the very medium it was gapfilled on.

The package's other integration paths already apply reversals: `MSGapfill.integrate_gapfill_solution`
and `MSBuilder.gapfill` both do, and `GapfillingPkg.run_test_conditions` includes them.
`build_and_gapfill` was the one path that did not.

## Change

- **New `MSBuilder.apply_gapfill_reversals(model, reversed_reactions)`.**
  - It opens each reversed reaction in the direction the gapfiller chose.
  - The bound is the one `get_reaction_constraints_from_direction` gives that direction:
    ±1000, the bound `build_and_gapfill` already gives a new reaction.
  - The other bound is left as it was, and a bound that is already wider is kept.
  - A reaction absent from the model is skipped with a warning.
  - It returns the ids it opened.
- **`build_and_gapfill` calls it** after `integrate_gapfill_solution`. Its return value is
  unchanged.

## Tests

`tests/core/test_msbuilder_gapfill_reversals.py` has 5 tests:
- the helper on small cobra models: reverse and forward openings, only the reversed side
  changing, wider bounds kept, absent reactions and empty input;
- a wiring test that runs `build_and_gapfill` with ATP correction and gapfilling stubbed,
  and checks that a reversed reaction comes back open while the base model is untouched.

On this branch the new tests and the existing `tests/core/test_msbuilder.py` pass (13
tests). On c2e31a0 the wiring test fails exactly where it should: the reversed reaction
keeps (0, 1000) instead of opening to (−1000, 1000).

## Validation against a pipeline that applies reversals itself

Four gapfilled models from prFBA's benchmark were built two ways in one process, from the
same EC set, template (21a8084) and complete medium: once through the package's
`build_and_gapfill`, and once through a harness that inlines it and applies the reversals
itself. The four are two Gram-positive and two Gram-negative, each needing one reversal.

| Model | Reversal | Growth, unfixed c2e31a0 | Growth, this branch | Growth, harness |
|---|---|---|---|---|
| 2582646.3 truth | rxn05289_c0 `<` | **0.0** | 14.587 | 14.587 |
| 1160095.3 family donor | rxn05289_c0 `<` | **0.0** | 16.723 | 16.723 |
| 1797343.3 truth | rxn00260_c0 `<` | **0.0** | 4.605 | 4.605 |
| 157688.7 PICRUSt2 | rxn00260_c0 `<` | **0.0** | 11.274 | 11.274 |

- On c2e31a0 the gapfill solution's new reactions are the same, but the reversed reaction
  keeps its original bound, so **none of the four models grows**.
- On this branch all four match the harness exactly: the same new reactions, the reversed
  reaction opened to (−1000, 1000), the same growth, and the same ATP sanity gate (pass).
- The branch run used the working tree that was then committed unchanged as e1254ba.

Script: prFBA `region_validation/shared_reference/scripts/validate_builder_reversals.py`.
Results: `data/builder_reversal_validation_{fixed,unfixed}.json`.

## Effect in a 200-genome benchmark

prFBA's phenotype benchmark gapfills 2,362 models from 200 genomes on complete medium
(templates 21a8084, ModelSEEDpy d884a50). Its harness inlines `build_and_gapfill` and
applies the reversals itself.
- 107 of the 2,362 models need a reversal; with the reversals applied, every model grows.
- An earlier run on the original templates dropped the reversals and left 136 models unable
  to grow. 70 of those are among the 107.

Details: freiburgermsu/prFBA, `region_validation/shared_reference/C12_SCOPE.md`.
