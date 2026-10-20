# Why the C12 validity gate fails: the gapfilled reactions that make ATP from nothing

2026-10-20. Scripts: `scripts/egc_rebuild.py` (rebuild with stage provenance and
counterfactual gapfill variants), `scripts/egc_analyze.py` (leak, stage replay, minimum
supports, minimum cuts). Per-model table: `EGC_C12_gated_models.csv`. Run artefacts
(models, per-model JSON, verifier scratch) are under `data/egc_analysis/` (git-ignored).
Every mechanism below was re-derived by an independent verifier with its own code.

## Scope and reproduction

C12 fails because 5 of the 200 gapfilled truth models (10 rows, one build shared by both
pipelines) pass ATP or biomass on water and protons alone (`phenotype_arm.sanity_gate`).
The same gate excludes 39 method-arm rows. All 49 gated rows sit in 20 genomes. Every
gapfilled model of those 20 genomes was rebuilt (218 models: truth once, five arms x two
pipelines), truth first exactly as `run_genome()` orders them, through a statement-for-
statement copy of `gapfill_block()` that snapshots the model after each stage (draft;
media evaluation; ATP-correction gapfill; genome-scale expansion; final MSGapfill).

The rebuild reproduces the stored rows exactly: 238 of 238 (`gapfilled.new_rxn` and
`sanity_gate_ok`), including all 49 gated rows.

## What leaks

- All 44 gated models make ATP (ATPM 150-1000); none makes biomass. None of the 174
  passing models leaks.
- Every leak is a **closed internal cycle**: zero net exchange, and it persists with every
  exchange closed, water and protons included. Internally it is ADP + Pi -> ATP + H2O from
  nothing; "ATP from water" is the gate's medium, not the mechanism.
- Stage replay: 190 of the 218 raw drafts leak (the known template artefact loops); the
  ATP-correction expansion step removes them by narrowing a median of 3 reactions; the
  final MSGapfill step re-opens a leak in exactly the 44 gated models (42 re-open a
  draft-level loop, 2 create a new one).
- In every gated model exactly one reaction added by MSGapfill is a single-reaction cut
  (removing it abolishes the leak without moving it elsewhere). Draft reactions in the
  same loop are also single cuts (F1-ATPase in 38/44, catalase 15, ...), so the
  gapfilled reaction is "the added reaction that closed the loop", not always "the
  reaction that is wrong".

## The five loop families

| Gapfilled reaction (as added) | Models | Templates | Truth models | Thermodynamically wrong step |
|---|---|---|---|---|
| rxn10966 nitrate reductase (thioredoxin), **reversed** | 14 | 11 GP, 3 GN | 586733.26, 3166642.3 | rxn10966 rev + rxn05289 rev = NADP-dependent nitrite oxidation, +34.8 kcal/mol (database rxn00572 '<') |
| rxn09180 anaerobic protoporphyrinogen oxidase, forward | 13 (+2) | mostly GN | — | draft **catalase rxn00006 run as 2 H2O + O2 -> 2 H2O2**, +45.2 kcal/mol (GOLD); every template allows only this direction ('<'), the database says '>' |
| rxn08167 L-aspartate oxidase (MK8), reversed | 8 | GP | 2862822.4, 3132284.3 | draft aspartate dehydrogenase rxn05119 (NAD) / rxn05117 (NADP) forward, +13.3 kcal/mol; template '>' vs database '<' |
| rxn00381 APS sulfohydrolase, **reversed** | 6 | 4 GP, 2 archaea | — | rxn00381 reversed, +19.6 kcal/mol (GOLD) |
| rxn00084 NADPH peroxidase, **reversed** | 1 | GN | 1797949.3 | rxn00084 reversed (2 H2O + NADP -> NADPH + H2O2), +76.5 kcal/mol (SILVER) |

In two models of genome 2045206.3 (base_family_donor), reactions added by the ATP-
correction gapfill (rxn60859, rxn14426) are alternative cuts of the rxn09180 loop.

## Root causes, in the ModelSEED stack

1. **Template direction errors.** The gapfiller builds candidates from `GapfillDirection`,
   on which 833 template reactions are reversible where the database says '>' and 336
   where it says '<' (GramNeg/GramPos v7; 835/338 Archaea v6). Catalase rxn00006 is '<'
   in all four templates, Core included. These supply the uphill step in 23 of the 44
   models (catalase, aspartate dehydrogenase).
2. **Negative reverse penalties.** 1,594 template reactions carry a negative
   `reverse_penalty`; the gapfilling objective costs a direction as
   2 x |base_cost + penalty|, which makes the reverse direction cheaper than the forward
   one for 1,190 of them. It decided the direction MSGapfill chose in 29 of 44 models
   (rxn10966, rxn08167, rxn00381, rxn00084). In 21 of those the chosen direction is
   itself thermodynamically infeasible.
3. **Dropped reversals.** `MSBuilder.build_and_gapfill` (and `phenotype_fold.gapfill_block`,
   which inlines it) integrates only the gapfiller's `new` reactions and drops its
   `reversed` set. 23 of the 218 rebuilt models have a dropped reversal. In the rxn10966
   family the shipped models grow only through the leak; once the leak is blocked, the
   gapfill solution needs rxn05289 reopened forward, which is dropped, so growth is zero
   until the reversal is applied (then 1.1-16.3, gate passing).

## Counterfactual gapfills (all 218 models of the 20 genomes)

| Variant | Gated models now passing | New leaks | Effect on growth |
|---|---|---|---|
| blacklist5: remove the five gapfilled reactions | 44/44 | 0 | 13 gated models stop growing (dropped reversals) |
| dirfix: block only the infeasible direction of four, blacklist rxn09180 | 44/44 | 0 | as blacklist5 |
| principled: database directions for rxn00006, rxn05119, rxn10966, rxn00381, rxn00084; integrate reversals | 36/44 | 0 | every model grows, including 10 passing models that had zero growth |
| **principled2: principled + rxn05117 (NADP aspartate dehydrogenase)** | **44/44** | **0** | **every model grows; passing models' growth unchanged (median ratio 1.000)** |
| dbdir: database direction for every database-irreversible template reaction | 44/44 | 0 | disruptive: 11 models stop growing, passing models' growth median ratio 0.75 |

Culprit-free gapfill solutions change slightly under any variant (the gapfilling MILP is
degenerate: 63/99 identical size, median Jaccard 0.93, growth ratio 1.00), so the
variants are comparable alternative optima, not different reconstructions.

## What this means for C12

The failure is a reconstruction-pipeline artefact that touches every arm (10 truth rows,
39 method-arm rows: prFBA 14, nn1 11, PICRUSt2 7, family donor 7), not a property of
any method. The minimal fix that clears it in all 20 affected genomes without changing
any passing model's growth is **principled2**: six template direction corrections
(rxn00006 '>', rxn05119 '<', rxn05117 '<', rxn10966 '>', rxn00381 '>', rxn00084 '>') plus
integrating the gapfiller's reversals. It is untested on the other 180 genomes of the
gapfilled pass; adopting it means a SPEC amendment and a re-run of that pass (the ATP
pass does not gapfill and is unaffected), after which C12 would be re-evaluated. The
integration fix also changes the growth readout of models with dropped reversals, so the
gapfilled phenotype numbers would move. That decision is the authors'.
