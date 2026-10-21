# Flip negative reverse penalties on reactions whose forward direction the database favours

Branch `fix-reverse-penalty-sign` (commit `4a39553`), stacked on
`fix-irreversible-reaction-directions` (`7d104ac`, the six-reaction fix). Until that PR
merges, this one shows both commits; the change here is the second.

## Summary

1,594 reactions in each genome-scale v7.0 template, and 26 in Core, carry a negative
`reverse_penalty`.

**A negative reverse penalty makes the reverse the cheaper direction.** `GapfillingPkg`
costs a direction at (2 - gene score) x |base_cost + penalty|. A negative reverse penalty
lowers base_cost + penalty, so the reverse becomes cheaper than the forward direction
for 1,190 of the 1,594 (all 26 in Core). The abs() only keeps the cost non-negative.

**For most of them the database says the forward direction is favoured.** Across all
4,031 such entries with a ModelSEED dG estimate:
- 92% are on reactions whose forward direction is favoured;
- the penalty tracks the forward dG (Spearman 0.84), as if its sign were flipped.

**Only a minority can affect gapfilling.** The gapfiller builds a reverse candidate, and
so uses the reverse penalty, only where the template's lower bound is negative and
`GapfillDirection` is not `>`. That leaves 410 per genome-scale template and 4 in Core.
Among those with a dG estimate (316 per template), the relationship is much weaker:
- 74% (233) are forward-favoured;
- the Spearman correlation is 0.16;
- 83 are on reactions whose reverse direction is favoured. For these, the negative
  penalty points the thermodynamically right way.

**This PR flips the sign only where the database's dG favours the forward direction:**
1,234 reactions per genome-scale template (233 of them reachable by the gapfiller) and
22 in Core. It leaves alone the 101 per template whose reverse direction is favoured,
and the 259 per template without a dG estimate. Nothing else in the files changes. No
template has a negative `forward_penalty`, and every `base_cost` is at least 1.

## Why it matters

In 44 of 218 gapfilled models we traced, the model made ATP from nothing:
- **The bias chose the leaking direction in 29 of them.** One of four reactions (rxn10966,
  rxn08167, rxn00381, rxn00084) was added in reverse, at 3.4-6.1 against 8-12 forward.
- **The bias alone is enough to cause most of those leaks.** On the original templates,
  negating every negative reverse penalty closes 27 of the 29.

The six-reaction fix closes all 44 by itself, so this PR is not needed for those leaks.
It removes the bias for the other reactions. Of the four culprits, only rxn08167 still
has a reverse candidate after the six-reaction fix. Its database dG favours the reverse
direction, so this PR leaves its penalty alone.

## Validation

We rebuilt 218 gapfilled models from 20 genomes, with the gapfiller's reversals
integrated:

| Templates | Leaks | New leaks | Models that do not grow | Passing models' growth (median ratio to original) | Passing models moved > 10% (of 164) |
|---|---|---|---|---|---|
| six-reaction fix | 0 | 0 | 0 | 1.000 | 13 |
| **six-reaction fix + this PR** | **0** | **0** | **0** | **1.000** | **18** |
| six-reaction fix + every negative penalty negated (alternative below) | 0 | 0 | 0 | 1.000 | 22 |

- Gapfill solutions shift, because the cheapest direction changes for some candidates:
  - 21 of 218 are identical to the six-reaction fix's (median Jaccard 0.94);
  - 12 of the 164 passing models moved by more than 10% against the six-reaction fix.
- The blanket alternative shifts more: 7 identical (median Jaccard 0.92), and 16 moved by
  more than 10%.
- This validation cannot show that the PR closes leaks, since the six-reaction fix
  already closes all 44.

## Alternatives

- **Negate every negative reverse penalty.** Branch `negate-all-negative-reverse-penalties`
  (`a82ced5`) does this. It also flips the 83 reachable reverse-favoured penalties per
  template, which then favour the thermodynamically disfavoured direction. We also
  applied it to the original templates, without the six-reaction fix, rebuilding the 142
  models of the 16 genomes behind the leaks. It opened one new leak there (2065.12,
  base_nn1, controlled).
- **Fix it in ModelSEEDpy instead.** Costing each direction as base_cost + |penalty| has the
  same effect as the blanket negation, given no negative forward penalties and base
  costs of at least 1. It would go where the sum is formed: in
  `GapfillingPkg.extend_model_with_template_for_gapfilling`, and in the copies in
  `core/fbabuilder.py` and `core/gapfillinghelper.py`. It would then apply to every
  template. Which layer to fix is the maintainers' call.

## Reproduce

In freiburgermsu/prFBA (branch `hardware-aware-alignment-dispatch`), under
`region_validation/shared_reference/`:
- `scripts/template_direction_audit.py`:
  - `apply --negate-reverse-penalty-where forward_favoured`;
  - `--negate-negative-reverse-penalty` for the blanket form;
  - `audit` for the statistics;
- `scripts/template_validation_summary.py`;
- `data/template_direction_audit/`: `direction_audit_summary.json` (the
  `negative_reverse_penalty` block), `validation_summary.json`, `validation_per_model.csv`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
