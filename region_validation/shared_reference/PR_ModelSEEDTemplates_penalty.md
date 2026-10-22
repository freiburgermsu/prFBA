# Flip negative reverse penalties on reactions whose forward direction the database favours

Branch `fix-reverse-penalty-sign` (commit `09f6bab`), stacked on
`fix-irreversible-reaction-directions` (`21a8084`, the five-reaction fix). Until that PR
merges, this one shows both commits; the change here is the second. It is an optional
clean-up: the five-reaction fix already closes every leak we traced.

## Summary

1,594 reactions in each genome-scale v7.0 template, and 26 in Core, carry a negative
`reverse_penalty`. The v5.0 and v6.0 templates carry the same ones. This PR changes
only v7.0; the ModelSEEDpy alternative below would cover all versions.

**A negative reverse penalty makes the reverse the cheaper direction.** `GapfillingPkg`
costs a direction at (2 - gene score) x |base_cost + penalty|. A negative reverse penalty
lowers base_cost + penalty, so the reverse becomes cheaper than the forward direction for
1,190 of the 1,594 (all 26 in Core). The abs() only keeps the cost non-negative.

**For most of them the database says the forward direction is favoured.** Across all
4,031 such entries with a graded ModelSEED dG estimate. The grades are the `best_source`
values in ModelSEEDDatabase dev's
`Biochemistry/Thermodynamics/SourceGrading/results/thermo_grades/reaction_grades.tsv`,
which master does not have:
- 92% are on reactions whose forward direction is favoured;
- the penalty tracks the forward dG (Spearman 0.84), as if its sign were flipped.

**Only a minority can affect gapfilling.** The gapfiller builds a reverse candidate, and
so uses the reverse penalty, only where the template's lower bound is negative and
`GapfillDirection` is not `>`. After the five-reaction fix (`21a8084`), that leaves 407
per genome-scale template and 4 in Core (412 and 5 on `main`). Among those with a dG
estimate (313 per template), the relationship is much weaker:
- 74% (233) are forward-favoured;
- the Spearman correlation is 0.18;
- 80 are on reactions whose reverse direction is favoured.

**This PR flips the sign only where the database's graded dG favours the forward
direction.** That is 1,234 reactions per genome-scale template and 22 in Core. For 859 of
them the reverse is the cheaper direction. The gapfiller can add 233 of them in reverse
(4 in Core), and for 227 of those it is steered towards the disfavoured direction. The
PR leaves alone the 101 per template whose reverse direction is favoured, and the 259 per
template without a dG estimate. Nothing else in the files changes. No template has a
negative `forward_penalty`, and every `base_cost` is at least 1.

## How the bias showed up

In 44 of 218 gapfilled models we traced, the model made ATP from nothing. In 29 of them,
one of four reactions (rxn10966, rxn08167, rxn00381, rxn00084) was added in reverse. The
negative penalty had made that the cheaper direction, at 3.4-6.1 against 8-12 forward.
Removing the bias alone, on the original templates and in its blanket form, prevents 27
of those 29 leaks. In the other 2 (Archaea, rxn00381) the reverse is chosen even at the
higher cost.

This PR is not what closes those leaks:
- **The five-reaction fix closes all 44.** After it, none of the four has a reverse
  candidate left.
- **This PR's own rule would have changed the culprit's penalty in only 7 of the 29.**
  rxn10966 (14 models) has no dG estimate, so it is left alone. rxn08167 (8 models) has a
  graded dG favouring the reverse: an eQuilibrator value of +12.6 kcal/mol, graded
  BRONZE. Group contribution and dGPredictor favour the forward direction for it.

The case for the PR is the bias itself, on the reactions that have not been seen to leak
yet.

## Validation

We rebuilt 218 gapfilled models from 20 genomes. Every row applies the gapfiller's
reversals:

| Templates | Leaking models | New leaks | Models that do not grow | Passing models' growth, median ratio to original | Passing models moved > 10% vs original (of 174) |
|---|---|---|---|---|---|
| original v7.0 | 44 | — | 0 | — | — |
| five-reaction fix | 0 | 0 | 0 | 1.000 | 18 |
| **five-reaction fix + this PR** | **0** | **0** | **0** | **1.000** | **18** |
| five-reaction fix + every negative penalty negated (alternative below) | 0 | 0 | 0 | 1.000 | 21 |

- **Gapfill solutions change in 203 of 218 models** (median Jaccard 0.93 against the
  five-reaction fix), and 19 of the 174 passing models move by more than 10%.
  - Negation only ever raises a candidate's cost.
  - In 189 of the 203, the earlier solution used none of the reactions whose penalty
    changed, and every stage before the final gapfill is identical. Those are alternative
    optima: the gapfilling MILP has many equal-cost solutions, and any change to
    candidate costs moves it.
- **The blanket form changes 199 solutions** (median Jaccard 0.94) and moves 9 models.
- **Neither form shifts growth more than the other in any consistent way.**
- **This validation cannot show that the PR closes leaks,** since the five-reaction fix
  already closes all 44.

## Alternatives

- **Negate every negative reverse penalty.** Branch `negate-all-negative-reverse-penalties`
  (`6430128`) does this. It also flips the 80 reachable reverse-favoured penalties per
  template, which then favour the thermodynamically disfavoured direction. We also
  applied it to the original templates, without the five-reaction fix. We rebuilt 142
  of the 174 models of the 16 genomes behind the 29 penalty-driven leaks (all arms but
  base_constant), with reversals not integrated. It opened one new leak among them
  (2065.12, base_nn1, controlled).
- **Fix it in ModelSEEDpy instead.** Costing each direction as base_cost + |penalty| has the
  same effect as the blanket negation, given no negative forward penalties and base costs
  of at least 1. It would go where the sum is formed: in
  `GapfillingPkg.extend_model_with_template_for_gapfilling`, and in the copies in
  `core/fbabuilder.py` and `core/gapfillinghelper.py`. It would then apply to every
  template. Which layer to fix is the maintainers' call.

## Reproduce

In freiburgermsu/prFBA (branch `hardware-aware-alignment-dispatch`), under
`region_validation/shared_reference/`:
- `scripts/template_direction_audit.py` (`--db <dev Biochemistry dir>` selects the
  database):
  - `apply --negate-reverse-penalty-where forward_favoured`;
  - `--negate-negative-reverse-penalty` for the blanket form;
  - `audit --templates <21a8084 checkout>/templates/v7.0` for the statistics, in the
    `negative_reverse_penalty` block of `direction_audit_summary.json`;
- `scripts/egc_rebuild.py --template-dir <checkout>/templates/v7.0 --variant reversed_only`;
- the counterfactual on `main` (`mainPEN`):
  - `apply --negate-negative-reverse-penalty` on main's templates;
  - then `egc_rebuild.py --template-dir <dir> --variant none --genomes <the 16 genomes>
    --arms base_family_donor,base_nn1,picrust2,prfba_exact_tie,truth` (all arms but
    base_constant);
- `scripts/template_validation_summary.py`: entries `PENcur2+rev`, `PEN2+rev`,
  `PR1v2+rev`, `base+rev`, `mainPEN`;
- `data/template_direction_audit/`: `validation_summary.json`, `validation_per_model.csv`.

The rebuild inputs (the benchmark's EC sets and the model files) are not in the
repository.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
