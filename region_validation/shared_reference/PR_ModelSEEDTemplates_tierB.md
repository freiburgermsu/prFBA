# Fix five reactions the templates allow only in the direction the database forbids

Branch `fix-opposite-template-directions` (commit `1f45108`), stacked on
`fix-irreversible-reaction-directions` (`7d104ac`, the six-reaction fix). Until that PR
merges, this one shows both commits; the change here is the second.

## Summary

The templates allow these five reactions only in the direction that the ModelSEED
database marks forbidden. In each case the database direction is also the physiological
one:

| Reaction | Name | Templates | Field: before → after | Why the database direction is right |
|---|---|---|---|---|
| rxn00206 | superoxide dismutase | GramNeg, GramPos, Archaea, Core | `direction` and `GapfillDirection`: `<` → `>` | the templates allow only O2 + H2O2 → 2 superoxide, the reverse of dismutation |
| rxn00380 | ADP-sulfurylase (sulfate adenylyltransferase, ADP) | GramNeg, GramPos, Archaea | `direction` and `GapfillDirection`: `>` → `<` | the enzyme runs APS + Pi → ADP + sulfate, in sulfite oxidation; the database call is GOLD, and all three estimators agree |
| rxn09566 | pyridoxamine oxidase | GramNeg, GramPos, Archaea | `direction` and `GapfillDirection`: `>` → `<` | the template direction makes pyridoxamine from pyridoxal, NH3 and O2 |
| rxn11937 | N2O reductase (NAD) | GramNeg, GramPos, Archaea, Core | `GapfillDirection`: `>` → `<` (`direction` already matches) | gapfilling could add N2 + NAD+ → N2O + NADH, the reverse of N2O reduction |
| rxn16828 | aminoacrylate peracid reductase (RutC) | GramNeg, GramPos, Archaea | `direction`: `<` → `>` (`GapfillDirection` stays `=`) | the template allows only the reverse of the peracid reduction |

The draft bounds and max fluxes follow `direction`. Each reaction is written the same way
round as in the database (shared-compound signs, protons ignored), and none has
non-default bounds. A full before/after parse confirms nothing else in the files
changes. `git diff --word-diff` shows the changes in the single-line Gram templates.

Two of the database entries also need attention, which this PR does not attempt:
- **rxn09566**'s stoichiometry is wrong: it consumes 2 H2O2 and releases 0.5 O2. It
  should be pyridoxamine + O2 + H2O → pyridoxal + NH3 + H2O2.
- **rxn16828** is mass-imbalanced in the database (H, O).

## How these five were chosen

An audit of every v7.0 template direction against the database found 38 reactions per
genome-scale template that the templates allow only in the database-forbidden direction.
We reviewed each against known physiology. Only these five turned out to be template
errors:
- for 22, the database direction is the wrong one, e.g. guanine deaminase, GMP reductase,
  heme d synthase, PPi-releasing adenylyltransferases, phospho-glycoside hydrolases, and
  two coenzyme M biosynthesis steps;
- 7 run in both directions in cells, so the database should say `=`: the shared
  succinate dehydrogenase / fumarate reductase half-reaction, ATP sulfurylase and two
  hydrogenases among them;
- 2 are mis-specified reactions that should be removed or blocked rather than flipped
  (rxn11944, rxn13150);
- 2 are unresolved.

The per-reaction readings, and a branch that applies all 38 for reference, are described
in the companion direction-audit issue.

## Validation

We rebuilt 218 gapfilled models (20 genomes; MSATPCorrection, then MSGapfill on complete
medium, with the gapfiller's reversals integrated) from the patched templates:

| Templates | Energy-generating leaks | New leaks | Models that do not grow | Passing models' growth (median ratio to original) | Passing models moved > 10% (of 164) |
|---|---|---|---|---|---|
| original v7.0 | 44 | — | 10 | — | — |
| six-reaction fix | 0 | 0 | 0 | 1.000 | 13 |
| **six-reaction fix + this PR** | **0** | **0** | **0** | **1.000** | **11** |

- Of the five reactions, superoxide dismutase is in 182 of the 218 drafts and rxn11937 in
  14. The other three occur in none of them, so this validation does not exercise them.
- Changing superoxide dismutase's draft direction shifts most gapfill solutions: 25 of 218
  are identical to the six-reaction fix's, median Jaccard 0.94. Growth barely moves:
  13 of the 164 passing models moved by more than 10% against the six-reaction fix, and
  the median ratio is 1.000.
- The validation uses complete medium, where the gapfiller can route around a blocked
  pathway.

## Reproduce

In freiburgermsu/prFBA (branch `hardware-aware-alignment-dispatch`), under
`region_validation/shared_reference/`:
- `scripts/template_direction_audit.py`:
  - `audit`: run against ModelSEEDDatabase `dev` (`078a395f`), it classifies every
    reaction the same way;
  - `apply --tiers B --reactions rxn00206,rxn00380,rxn09566,rxn11937,rxn16828`;
- `scripts/egc_rebuild.py --template-dir ... --variant reversed_only`;
- `scripts/template_validation_summary.py`;
- `data/template_direction_audit/`: `direction_audit.tsv`, `tierB_reactions.tsv`,
  `tierB_readings.tsv`, `validation_summary.json`, `validation_per_model.csv`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
