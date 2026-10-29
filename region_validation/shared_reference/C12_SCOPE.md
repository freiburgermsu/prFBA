# How far the ATP-from-nothing problem reaches, and what amendment A-C12 changes

2026-10-27/28. Companion to `EGC_C12_REPORT.md` (why C12 fails) and SPEC section 14
(amendment A-C12).

**Provenance.** Every number here comes from runs on the following code and data:
- ModelSEEDpy `d884a50`, the commit the pre-registered run imported, pinned through
  `PYTHONPATH` (SPEC section 14);
- templates main `2294299` (original) and `21a8084` (the five-reaction fix);
- ModelSEEDDatabase `e8446e3e`.

**Scripts.**
- Draft-level scan: `scripts/egc_scope_scan.py` and `scripts/egc_draft_loops.py`, run by
  `scripts/egc_scope_chain_pinned.sh`;
- probe templates: `scripts/egc_template_variants.py`;
- before/after comparison: `scripts/c12_compare.py`;
- the comparisons and the genus rescoring: `scripts/c12fix_finish.sh`.

**Data.** The following are committed, with `.meta.json` provenance sidecars where they
apply. Row-level phenotype files and models stay git-ignored.
- `data/egc_analysis/scope_scan_{200,fixonly,variants}.jsonl`;
- `data/egc_analysis/draft_loops_{200,fixonly,variants}.jsonl`;
- `data/egc_analysis/templates_scope/manifest.json`;
- `data/c12fix_compare_{gf,atp}.json`;
- `data/phenotype_{gf,atp}_c12fix/phenotype_summary.json`;
- the rescored `data/summary_genus.json` and `data/controls.json`. The pre-amendment pair
  is kept beside them as `*.pre-c12fix` and in git history.

## Summary

- **C12 passes under amendment A-C12.**
  - The gapfilled pass gates 0 of 8,614 rows, against 49 before (10 of them truth rows).
  - There is no new leak in any of the 200 genomes, including the 180 outside the
    amendment's validation.
  - Rescored, the leave-genus-out artefact is valid (exit 0, all 15 controls passing).
  - Every EC number, the answer-key fingerprint and the scorer fingerprint are unchanged.
- **The amendment moves no ATP-pass number.** The ATP pass reproduces the pre-registered
  pass row for row.
- **The gapfilled readouts change slightly, with the same pattern.**
  - Every gapfilled model now grows; 136 did not.
  - All 200 genomes are scored, where the pre-registered run dropped 5.
- **Most raw drafts make ATP from nothing** (72%), in all 200 genomes. ATP correction
  removes these loops before gapfilling, so none reaches a readout.
  - Three template-direction families account for them: GAPN, maltodextrin glucosidase and
    cytochrome c peroxidase.
  - Seven direction changes remove 89% of the draft leaks. This was tested at the draft
    level only.
- **The five-reaction fix exposes one latent draft loop** (rxn00014 reversed), in 14
  drafts of 7 genomes. It survives into no final model.
- **The first launch imported a moved ModelSEEDpy.** It was caught by a 1–4-reaction
  difference in the drafts. Both passes were rerun from zero, pinned to the pre-registered
  commit, and the pin reproduces every stored row rebuilt.
- **Amendment A-C13 reran the arm on the updated ModelSEEDpy (`c2e31a0`, section 7).**
  - That commit fixes demand reactions built as free sources and an ATP-correction step
    that undid its own loop filter.
  - C12 passes again; every ATP-pass number is unchanged, and the gapfilled medians move by less
    than 0.005.
  - The fork's dropped reversals are fixed in a branch whose builder reproduces this
    harness's models.

## 1. Most raw drafts make ATP from nothing

The scan covers every (genome, arm) draft of the gapfilled pass: the same 200 genomes,
arms, EC sets and seed. It builds each draft exactly as `gapfill_block` does, adds ATPM,
and applies the C12 gate's medium: every uptake closed, water and protons free. Nothing
is gapfilled.

| | Original templates | Five-reaction fix |
|---|---|---|
| leaking drafts | 851 of 1,181 (72.1%) | 865 of 1,181 (73.2%) |
| leaking truth drafts | 136 of 200 | 138 of 200 |
| genomes with a leaking draft | 200 of 200 | 200 of 200 |
| leaks at the 1,000 ATPM bound | 833 | 833 |
| drafts making biomass from nothing | 0 | 0 |

**By template:**
- Archaea 98/142;
- GramNeg 373/515, or 387 with the fix;
- GramPos 380/524.

**By arm, on the original templates:**
- constant baseline 200/200 (one large EC set for every genome);
- truth 136, prFBA 136, PICRUSt2 136, nn1 133;
- family donor 110 of 181.

**This is the drafts' normal state, not a property of the 20 genomes C12 flagged.** It
does not reach the readouts. The ATP readout is the core-model ATP vector, and its
empty-medium gate passed in all 43,140 pre-registered rows (one per genome, arm and
pipeline). In the gapfilled pass,
MSATPCorrection's expansion step closes the draft loops before MSGapfill runs
(`EGC_C12_REPORT.md`, stage replay). C12 is about what MSGapfill re-opens.

**The leak verdicts do not depend on the ModelSEEDpy commit.** The same scan on the moved
fork (c2e31a0) gave identical leak flags and empty-medium ATP for all 2,362 drafts; only
the demand-reaction counts differed.

## 2. Three loop families in the raw drafts

The scan enumerated up to five minimum-cardinality supports of ATPM for each of the 136
leaking truth drafts on the original templates. A draft can hold several families.

| Family | Reactions | The thermodynamically wrong step | Template `direction` / database | Truth drafts (of 136) |
|---|---|---|---|---|
| **GAPN** | rxn23850 (non-phosphorylating G3P dehydrogenase) reversed, GAPDH rxn00781, PGK rxn01100 | GAPN run as 3PG + NADH → G3P + NAD + H2O, uphill by the database's 9.75 kcal/mol (about 41 kJ/mol) | `=` in all three genome-scale templates / `>` (ΔG'° −9.75 kcal/mol, GOLD) | **130**, always the first (smallest) support |
| **Maltodextrin** | maltodextrin glucosidase rxn08934/rxn08935 reversed, maltodextrin phosphorylase rxn08941/rxn08942, phosphoglucomutase rxn00704, glucokinase rxn00216 reversed; or the alpha-anomer route rxn15249 + rxn20583 bridged by lactose hydrolysis (rxn00816, rxn15230) | glucosidase and glucokinase both run backwards | glucosidase `=` / `>`; glucokinase `direction` `=` while its `GapfillDirection` is `>` / database `?` | 24 |
| **Cytochrome c peroxidase** | rxn00014 reversed, cbb3 oxidase (rxn14426, rxn60859), F1-ATPase rxn08173; on the original templates protoporphyrinogen oxidase (rxn20627, rxn02304) stands in for catalase | rxn00014 run as 2 ferricytochrome c + 2 H2O → 2 ferrocytochrome c + H2O2, about +214 kJ/mol (E'° of H2O2/H2O +1.35 V against cytochrome c +0.25 V) | `=` / `?` (no database direction) | 7 |

GAPN alone opens a loop in 65% of all truth drafts (130 of 200).
- Neither GAPN nor rxn00014 is among the five reactions the fix corrects. Neither loop
  appears in any of the 44 gated final models: ATP correction removes both.
- The direction audit (issue ModelSEED/ModelSEEDTemplates#10) lists rxn23850 only in bulk,
  among the 1,166 reversible-in-template reactions (tier C).
- It cannot see rxn00014 at all, since the database gives that reaction no direction.

## 2b. Seven direction changes remove most draft-level leaks

Probe templates (`scripts/egc_template_variants.py`) restrict more reactions to their forward
direction on top of 21a8084, edited the way 21a8084 edits a reaction. Scanned over the same
1,181 drafts (`data/egc_analysis/scope_scan_variants.jsonl`):

| Templates | Leaking drafts | Leaking truth drafts | At the 1,000 bound |
|---|---|---|---|
| original (2294299) | 851 (72.1%) | 136 | 833 |
| five-reaction fix (21a8084) | 865 (73.2%) | 138 | 833 |
| + rxn00014 `>` | 849 (71.9%) | 136 | 833 |
| + rxn00014, GAPN rxn23850 and the maltodextrin glucosidases (rxn05746, rxn08933–rxn08936), all `>` | **95 (8.0%)** | **21** | 53 |

- **rxn00014 forward-only closes the loop the fix opens:** all 14 fix-only drafts, plus 2
  drafts whose only loop was the original templates' rxn00014 route.
- **With GAPN and the glucosidases also forward-only, 754 of the 849 remaining leaking drafts
  stop leaking (89%).** Every constant-baseline draft is among them, although all 200
  leaked before.
- **The 95 drafts that still leak hold a long tail of a different kind.** The same
  conversion is catalysed with two electron carriers, both reversible: run forward with one
  carrier and backward with the other, it interconverts them, and a proton-pumping step then
  makes ATP (F1-ATPase in 19 of the 21 truth drafts). Most of them leak below the bound
  (11 of 21). Examples:
  - the Mg-protoporphyrin cyclase steps with a generic acceptor (rxn11665–rxn11667) against
    NAD (rxn11940–rxn11942);
  - D-lactate dehydrogenase with cytochrome c (rxn00146) against ubiquinone (rxn08783),
    through cytochrome bc1 (rxn13689);
  - glutamate synthase with NADP (rxn00085) against ferredoxin (rxn12822), with the
    ATP-driven protochlorophyllide reductase (rxn45223) run backwards;
  - hydrogen cycling (rxn05887, rxn08736, rxn10124).
- **These probes are draft-level only.** None was run through ATP correction or gapfilling,
  so their effect on growth is untested. They are not part of amendment A-C12, which
  pre-commits to no further template changes.

## 3. The fix opens one loop at the draft level, in 14 drafts

The five-reaction fix closes no draft-level loop: 0 drafts go from leaking to not leaking.
It opens one, in 14 drafts that did not leak before.
- **Where:** 7 GramNeg genomes (157466.6, 1903292.3, 1912957.4, 2580566.3, 2715253.317,
  2818388.4, 2949308.3).
- **Which arms:** nn1 6, prFBA 5, truth 2 (1903292.3 and 2715253.317), family donor 1.

All 14 drafts have the same single support, at ATPM 500:
- catalase rxn00006, now forward;
- rxn00014 reversed;
- cbb3 oxidase rxn14426 and rxn60859;
- F1-ATPase rxn08173.

Catalase forward is the physiological direction, so the fix itself is right. On the
original templates catalase ran only as 2 H2O + O2 → 2 H2O2, and that blocked this route.
The fix exposes the latent error, rxn00014's reversibility. Restricting rxn00014 to `>`
stops all 14 fix-only drafts leaking, and 2 more drafts whose only loop was the original
templates' rxn00014 route (section 2b).

**The loop does not survive into a final model.** In the pinned gapfilled pass, all 14
fix-only drafts, in both pipelines, gapfill to models that pass the gate and grow (growth 2.0
to 12.3). Every other model of the 7 genomes passes too. MSATPCorrection closes the loop,
and MSGapfill does not reopen it.

The exposure is therefore latent in this pipeline. It still matters to anyone who uses the
drafts directly, or gapfills without ATP correction.

## 4. The final models under the amendment

Pinned gapfilled pass, 200 genomes, 8,614 rows matched one to one with the pre-registered
pass (`data/c12fix_compare_gf.json`).

| Gated rows | truth | prFBA | nn1 | PICRUSt2 | family donor | all |
|---|---|---|---|---|---|---|
| pre-registered | 10 | 14 | 11 | 7 | 7 | 49 |
| amended | 0 | 0 | 0 | 0 | 0 | **0** |

**Both legs of C12 pass.** Leg (b) has 0 gated truth models. Leg (c) gates 0 of 8,614 rows,
against a 10% ceiling. All 49 leaks are fixed, and no model leaks that did not leak before.
That holds in the 20 genomes of the validation and in the other 180.

**C11 still passes.** Two builds from one EC set give identical reaction sets, ATP vectors
and gapfill solutions.

**Every gapfilled model grows.** Of the 2,362 main-arm gapfilled models, 136 in 50 genomes
did not grow before; all grow now.
- 107 models apply a reversal the fork would have dropped; the harness applies every
  reversal the gapfiller proposes.
- 70 of the former non-growers are among them.
- The other 66 grow on a gapfill solution the corrected templates changed.

**The models that grew before keep their growth:** median after/before ratio 1.000, with
274 of 2,177 moving by more than 10%. The gapfill burden changes in 1,304 models: by a mean
of +0.21 reactions and a median of 0.

**The 7 fix-only genomes pass.** Every model passes the gate and grows (section 3).

## 5. What the amendment does to the readouts

**ATP pass (1,000 genomes): unchanged.**
- The rerun on the fixed templates reproduces the pre-registered ATP pass exactly: all 43,140
  rows match in every field but wall time, including the core ATP vectors, the draft
  reaction counts and every vs-truth readout.
- The summary's arms, contrasts, ladders and pipeline-effect blocks are identical, so every
  ATP-pass number reported stands.
- **Why the fix leaves it unchanged:** the five corrections change reaction directions, not
  which reactions a draft holds, so the structural readouts cannot move. None of them changes
  which media the core model makes ATP on.
- **The pre-registered file held 1,208 duplicate rows:** 31 genomes, identical but for wall
  time. The summarizer already keeps one row per (genome, arm, pipeline), so only its
  `meta.n_models` was inflated (44,348 against 43,140); no statistic is affected.

**Gapfilled pass: the same pattern on all 200 genomes.**
- Values are medians from the controlled pipeline; the pre-registered value is in
  parentheses.
- n is now 200 for every arm but the family donor (181), against 170–195 before.

| Arm | Essentiality disagreement | Growth log-ratio to truth | Gapfill burden |
|---|---|---|---|
| truth | 0 | 0 | 81.5 (83) |
| prFBA | 0.069 (0.068) | 0.004 (0.002) | 79 (80.5) |
| PICRUSt2 | 0.066 (0.066) | −0.007 (−0.009) | 83 (85) |
| nn1 | 0.070 (0.070) | 0.006 (0.004) | 83 (84) |
| family donor | 0.068 (0.064) | 0.004 (0.000) | 83 (83) |
| constant | 0.091 (0.092) | −0.149 (−0.067) | 72 (71) |

- **The ordering survives.** PICRUSt2 is marginally ahead on essentiality, and prFBA is
  indistinguishable from nn1.
- **The means move more than the medians.** prFBA's mean growth log-ratio falls from 0.255
  to 0.044, because the leaking and non-growing models are gone.

**Genus rescoring.**
- `score_folds.py` ran in fold mode against the amended gapfilled summary, with the pinned
  main 2294299 templates and the pinned ModelSEEDpy. The scorer imports it.
- **It reused the cells table:** the fingerprints matched.
- **C11, C12 and every other control pass** (C1 and C2 come from their own modes), so the
  artefact is valid.
- **Every EC number is unchanged.** The arms, contrasts and nested-CV blocks are identical;
  only validity flags and provenance differ.
  - Under the amendment, prFBA reaches population F1 0.687 against PICRUSt2's 0.779
    (−0.092, 95% CI −0.107 to −0.081).
- **The family summary is not rescored.** It was scored without the phenotype arm, so it
  carries no C11 or C12 and was valid already. SPEC section 14's "the genus and family
  summaries are rescored the same way" was wrong for it.

## 6. Reach beyond prFBA

1. **ModelSEED templates.** The draft-level loops are template direction errors that every
   ModelSEED reconstruction inherits, not just prFBA's.
   - GAPN rxn23850 is the single highest-yield fix: one GOLD-graded direction (database
     `>`), and it opens the loop in 96% of the leaking truth drafts. It stays invisible in
     a pipeline that runs ATP correction.
   - rxn00014 needs a thermodynamic argument, because the database has no direction for
     it.
   - Both belong in a follow-up to the five-reaction PR and to issue #10; see the note
     below.
2. **The fork's dropped reversals.** `MSBuilder.build_and_gapfill` in the fork drops the
   gapfiller's `reversed` set (`EGC_C12_REPORT.md`, root cause 3). That affects every user
   of the fork's builder, not only this harness. It is fixed on the branch
   `fix-gapfill-reversals` (section 7), and the PR is ready.
3. **The moving ModelSEEDpy.** f395ae6 (ATP-expansion filter) and 8ac600b (demand
   reactions) change gapfilled results: all 80 rows of the two smoke genomes differ.
   Amendment A-C13 adopted them, pinned at c2e31a0; see section 7.
4. **D11.** SPEC D11 dropped the minimal-media growth panel partly because
   minimal-medium gapfilling "produces models that fail a thermodynamic gate". Those are
   the same leaks: gapfilling closes them against a medium. Whether the amended templates
   rescue a minimal-media panel is untested here.

## 7. Amendment A-C13: the same arm on the updated ModelSEEDpy (2026-10-29)

**What changed.** Both passes were rerun with ModelSEEDpy `c2e31a0` instead of `d884a50`
(SPEC section 15); everything else is as in A-C12. Two build-path fixes separate the
commits:
- `8ac600b`: demand reactions are no longer built as free sources of six drain compounds,
  which the C12 gate could not see because it closes only exchanges;
- `f395ae6`: `build_tests()` no longer undoes the ATP-correction expansion filter when no
  test medium is selected.

**ATP pass: unchanged.** Against A-C12, all 43,140 rows match except their draft reaction
counts (the demand reactions). Every ATP vector and every statistic in the summary is
identical, so every ATP-pass number stands.

**Gapfilled pass.**
- **Gate:** 0 of 8,614 rows gated (A-C12: 0; pre-registered: 49), and no new leak
  in any genome.
- **Growth:** every model grows, against 136 non-growing models pre-registered.
- **Reversals:** 114 models apply one, against 107 under A-C12.
- **Against A-C12:** the median growth ratio is 1.000, with 87 of 2,362 models
  moving by more than 10%. The gapfill burden changes in 1,068 models, by a mean of
  +0.016 reactions.
- **C11:** it passes, with identical reaction sets, ATP vectors and gapfill solutions.

Controlled medians, A-C13 / A-C12 / pre-registered (the pre-registered run drops the five
gated genomes):

| Arm | Essentiality disagreement | Growth log-ratio to truth | Gapfill burden |
|---|---|---|---|
| truth | 0 / 0 / 0 | 0 / 0 / 0 | 82.5 / 81.5 / 83 |
| prFBA | 0.068 / 0.069 / 0.068 | 0.008 / 0.004 / 0.002 | 80 / 79 / 80.5 |
| PICRUSt2 | 0.065 / 0.066 / 0.066 | −0.006 / −0.007 / −0.009 | 84 / 83 / 85 |
| nn1 | 0.069 / 0.070 / 0.070 | 0.008 / 0.006 / 0.004 | 82.5 / 83 / 84 |
| family donor | 0.066 / 0.068 / 0.064 | 0.006 / 0.004 / 0 | 83 / 83 / 83 |
| constant | 0.090 / 0.091 / 0.092 | −0.134 / −0.149 / −0.067 | 72 / 72 / 71 |

- The method arms' essentiality and growth medians move by less than 0.005, and their burden
  by at most one reaction.
- None of the comparisons the paper draws changes: PICRUSt2 is marginally ahead on
  essentiality, prFBA is indistinguishable from nn1, and the burden is about the truth
  model's.

**Genus rescoring.**
- Against `data/phenotype_gf_c13/phenotype_summary.json`, C11 is PASS and C12 is PASS.
  The leave-genus-out artefact is valid (exit 0).
- The answer-key fingerprint is unchanged (a9dc024d); the EC numbers are those of A-C12 and the
  pre-registered run.
- The A-C12-scored summary and controls are kept as `*.a-c12`.

**The reversal fix in ModelSEEDpy itself.** The fork's `MSBuilder.build_and_gapfill` still
drops the gapfiller's reversals in `c2e31a0`.
- **The fix:** branch `fix-gapfill-reversals` (e1254ba), adding
  `MSBuilder.apply_gapfill_reversals`. PR text: `PR_ModelSEEDpy_gapfill_reversals.md`.
- **Validation:** four models needing a reversal were built through the package and
  through this harness (`scripts/validate_builder_reversals.py`). On the branch, all four
  are identical to the harness. On `c2e31a0`, the package returns all four unable to grow.
- **This harness does not call `build_and_gapfill`.** On the branch it reproduces A-C13's
  rows exactly (72 of 72), so A-C13 stands for the fixed package as well.

## Note for the five-reaction PR (draft, for the authors to post)

> **A follow-up from a 200-genome scan: one latent loop this PR exposes, and the commonest
> draft-level loop**
>
> We scanned the raw (un-gapfilled) drafts of 200 genomes from our benchmark: six EC sets
> each, 1,181 drafts, with every uptake closed and water and protons free. We used the
> original templates and this PR.
>
> **The scan:**
> - **Without gapfilling, most drafts make ATP from nothing:** 72% on the original templates,
>   73% with this PR.
> - In our pipeline MSATPCorrection removes these loops before gapfilling. So they reach no
>   final model, and with this PR none of our 8,614 gapfilled models leaks.
> - They do matter to anyone who uses drafts directly, or gapfills without ATP correction.
>
> **1. This PR exposes one latent loop (14 drafts in 7 Gram-negative genomes).**
> - Catalase now runs forward, which is correct.
> - Cytochrome c peroxidase (rxn00014) is reversible in all three genome-scale templates; the
>   database gives it no direction. Reversed, it makes H2O2 from water and ferricytochrome c.
> - Catalase turns that H2O2 into O2, and cbb3 oxidase (rxn14426, rxn60859) uses the O2 to
>   pump protons for ATP synthase.
> - The reversal is uphill by about 214 kJ/mol (E'° of H2O2/H2O is +1.35 V; of cytochrome c,
>   +0.25 V).
> - On the original templates, protoporphyrinogen oxidase (rxn20627 with rxn02304) stands in
>   for catalase in 7 more drafts.
> - **Restricting rxn00014 to `>` on top of this PR stops all 14 drafts leaking,** and 2 more
>   drafts whose only loop was the original-template route. In our runs the loop never
>   survived ATP correction; we suggest `>` as hygiene.
>
> **2. The commonest draft loop is GAPN (rxn23850).**
> - It is reversible in all three genome-scale templates, although the database says `>`
>   (ΔG'° −9.75 kcal/mol, GOLD).
> - Reversed, with GAPDH (rxn00781) and PGK (rxn01100), it makes ATP from ADP and phosphate.
>   It is the smallest loop in 130 of the 136 leaking truth drafts.
> - **Five maltodextrin glucosidases** (rxn05746, rxn08933–rxn08936) are reversible too.
>   Hydrolysis is irreversible, and the database says `>` for three of them. Run backwards
>   with maltodextrin phosphorylase and glucokinase, they close a second loop.
>   - Glucokinase rxn00216 has `direction` `=` but `GapfillDirection` `>`.
>
> **3. Together these remove most draft-level leaks.** Restricting rxn00014, rxn23850 and the
> five glucosidases to `>`, on top of this PR, cuts the leaking drafts from 849 to 95 (89%).
> The remainder is a long tail of carrier-interconversion loops: one conversion with two
> cofactors, both reversible.
>
> **Caveat:** these seven changes were tested at the draft level only, not through ATP
> correction and gapfilling, so their effect on growth is untested. We are happy to validate
> them on the 218 models we used for this PR and open a follow-up PR. GAPN and the
> glucosidases are among the 1,166 reversible-in-template reactions of issue #10. rxn00014
> is not, because the database has no direction for it.
>
> Data and scripts: freiburgermsu/prFBA, region_validation/shared_reference/C12_SCOPE.md.
