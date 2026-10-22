# Fix five reactions whose template direction opposes the database

Branch `fix-opposite-template-directions` (commit `0805094`), stacked on
`fix-irreversible-reaction-directions` (`cf63a8f`, the five-reaction fix). Until that PR
merges, this one shows both commits; the change here is the second.

## Summary

In each of these five reactions, a template direction field opposes the ModelSEED
database, and the database direction is the physiological one. For four of them the
templates allow only the direction the database forbids. For rxn11937 the
`GapfillDirection` contradicts its own draft direction, so gapfilling cannot add the
reaction at all.

| Reaction | Name | Templates | Field: before → after | Why the database direction is right |
|---|---|---|---|---|
| rxn00206 | superoxide dismutase | GramNeg, GramPos, Archaea, Core | `direction` and `GapfillDirection`: `<` → `>` | the templates allow only O2 + H2O2 → 2 superoxide, the reverse of dismutation |
| rxn00380 | ADP-sulfurylase (sulfate adenylyltransferase, ADP) | GramNeg, GramPos, Archaea | `direction` and `GapfillDirection`: `>` → `<` | the enzyme runs APS + Pi → ADP + sulfate in sulfite oxidation; the database call is GOLD, and all three estimators agree |
| rxn09566 | pyridoxamine oxidase | GramNeg, GramPos, Archaea | `direction` and `GapfillDirection`: `>` → `<` | the template direction makes pyridoxamine from pyridoxal, NH3 and O2 |
| rxn11937 | N2O reductase (NAD) | GramNeg, GramPos, Archaea, Core | `GapfillDirection`: `>` → `<` | `direction` already allows only N2O reduction (`<`). `GapfillDirection` `>` contradicts it, which leaves the gapfill candidate bounded 0/0 in ModelSEEDpy. The fix makes N2O reduction available to gapfilling |
| rxn16828 | aminoacrylate peracid reduction (annotated to RutC) | GramNeg, GramPos, Archaea | `direction`: `<` → `>` (`GapfillDirection` stays `=`) | the draft allowed only the reverse of the peracid reduction |

The draft bounds and max fluxes follow `direction`. Each reaction is written the same way
round as in the database (shared-compound signs, protons ignored), and none has
non-default bounds. A full before/after parse confirms nothing else in the files
changes. `git diff --word-diff` shows the changes in the single-line Gram templates.

**Database entries that need attention.** This PR does not attempt these:
- **rxn09566:** the stoichiometry is wrong. It consumes 2 H2O2 and releases 0.5 O2, and
  should be pyridoxamine + O2 + H2O → pyridoxal + NH3 + H2O2.
- **rxn16828:** the reaction is mass-imbalanced (H, O).

**Three more reactions have the rxn11937 conflict.** Pyruvate kinase (rxn00148),
pyruvate formate-lyase (rxn00157) and ADP-forming acetate:CoA ligase (rxn00172) carry
`direction` `<` with `GapfillDirection` `>` in all four templates, so gapfilling cannot
add them either. The database gives them `?` or `=`, so this PR leaves them out of
scope.

## How these five were chosen

An audit of every v7.0 template direction against the database (on top of the
five-reaction fix) found 41 reactions per genome-scale template that the templates allow
only in the database-forbidden direction, in at least one field. We reviewed each against
known physiology. Only these five turned out to be template errors. Of the rest:
- **24:** the database direction is the wrong one. Examples include:
  - guanine deaminase, GMP reductase and heme d synthase;
  - NadX, and the L-aspartate oxidase that the five-reaction fix sets forward;
  - PPi-releasing adenylyltransferases and phospho-glycoside hydrolases;
  - two coenzyme M biosynthesis steps.
- **8:** reversible or near equilibrium, so the database should say `=`. Among them are the
  succinate dehydrogenase / fumarate reductase half-reaction, ATP sulfurylase and two
  hydrogenases. The two HQT entries are reversible only in vitro and are read with low
  confidence.
- **2:** mis-specified reactions that should be removed or blocked rather than flipped
  (rxn11944, rxn13150).
- **2:** unresolved.

The per-reaction readings, and a branch that applies all 41 for reference, are described
in the companion direction-audit issue.

## Validation

We rebuilt 218 gapfilled models (20 genomes; MSATPCorrection, then MSGapfill on complete
medium) from the patched templates. Every row applies the gapfiller's reversals:

| Templates | Leaking models | New leaks | Models that do not grow | Passing models' growth, median ratio to original | Passing models moved > 10% vs original (of 174) |
|---|---|---|---|---|---|
| original v7.0 | 44 | — | 0 | — | — |
| five-reaction fix | 0 | 0 | 0 | 1.000 | 18 |
| **five-reaction fix + this PR** | **0** | **0** | **0** | **1.000** | **19** |

- **Against the five-reaction fix alone,** 14 of the 174 passing models moved by more
  than 10%. 18 of the 218 gapfill solutions are identical to the five-reaction fix's
  (median Jaccard 0.94).
- **The solutions shift even in drafts that carry none of the five reactions.** All 36
  such drafts are identical between the two runs, yet every one gets a different
  solution (median Jaccard 0.92). Any change to the gapfill candidate pool moves the MILP
  to another near-equal optimum.
- **Only superoxide dismutase is exercised.** It is in 182 of the 218 drafts.
  - rxn00380, rxn09566 and rxn16828 are in no draft.
  - rxn11937 is in 14 drafts, but its change is to `GapfillDirection` only.
  - The gapfiller adds none of the four.
- **The validation uses complete medium,** where the gapfiller can route around a blocked
  pathway.

## Reproduce

In freiburgermsu/prFBA (branch `hardware-aware-alignment-dispatch`), under
`region_validation/shared_reference/`:
- `scripts/template_direction_audit.py`:
  - `audit`: the committed audit used the ModelSEEDDatabase fork at `e8446e3e`, which is
    `dev` `078a395f` plus one unmerged stoichiometry-correction commit. Against `dev`
    itself, every reaction is classified the same way;
  - `apply --tiers B --reactions rxn00206,rxn00380,rxn09566,rxn11937,rxn16828`;
- `scripts/egc_rebuild.py --template-dir <checkout>/templates/v7.0 --variant reversed_only`;
- `scripts/template_validation_summary.py`: entries `Bcur2+rev`, `PR1v2+rev`, `base+rev`;
- `data/template_direction_audit/`: `direction_audit.tsv`, `tierB_reactions.tsv`,
  `tierB_readings.tsv`, `validation_summary.json`, `validation_per_model.csv`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
