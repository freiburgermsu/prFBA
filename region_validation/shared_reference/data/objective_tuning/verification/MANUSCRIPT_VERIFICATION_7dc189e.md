# Verification of the uncommitted Stage E / Stage O manuscript edits

2026-11-03. Target: `/home/freiburger/Documents/prFBA_paper`, branch `depth100-standard`, working tree on `62fbaf1` (`prFBA.tex`, `prFBA_andrew_voice.tex`). Nothing in any repository was edited or committed. Every script, build and output is in this folder.

## How it was checked

- **The edits themselves.** HEAD plus the edit script's 11 replacements reproduces each working file byte for byte, in both variants (`replay.py`).
- **Numbers.** Every number in the new text was traced to a results file. The Stage O contrasts, the G effect, the half split and the gapfilled-active share were re-computed from the raw rows (`check2.py`, `check3.py`). Stage E figures were read from `scores.json`, `selected.json` and `assertions.json`.
- **Code.** For the code claims I read c2e31a0's `prFBA`, MSCommunity's constructor and `MSModelUtil.__init__`, and the 09856f5 diff.
- **Template report.** The figures in `REPORT.md` were checked against `scan/scan_summary.json` and `gf_pass/compare_vs_c13.json`.
- **Validity of the fixes.** All 23 proposed corrections were checked with `corrections.py`: each old string occurs exactly once in each variant, including when they are applied in order. Both variants were compiled before the corrections (`build/`) and after them (`build2/`).

## Verdict per edit

| # | Edit | Verdict |
|---|---|---|
| 1 | Abstract | **Not confirmed as worded.** It drops the exchange-weight change from the result and overstates the truth proxy (C5). |
| 2 | Community kinetic caveat | **Not confirmed as worded.** On c2e31a0 the community path returns no solution at all (C4). |
| 3 | ModelSEEDpy fix commits (09856f5) | Confirmed. |
| 4 | Stage O paragraph | Every number confirmed. One misstatement (C6). Recommended: R1, R2, R3, R7a, R8, R13. |
| 5 | Scope: "the exceptions, two post-hoc in-silico tests" | Confirmed. |
| 6 | "Apart from those tests" | Confirmed. |
| 7 | "evaluated only by the post-hoc in-silico tests ..., against an alternative estimate" | Confirmed. |
| 8 | Neighbour-and-clade passage after the calibration sentence | **Two errors** (C1, C2); the rest is confirmed. Recommended: R4a, R4b, R5, R11. |
| 9 | Discussion sentence | Numbers confirmed, but it **contradicts the unchanged opening of the same sentence** (C3). Recommended: R7b, R9. |
| 10 | Limitations | Confirmed. Recommended: R7c. |
| 11 | Seven-reaction template sentence | Confirmed. Recommended: R10. |

## 1. Abstract (edit 1): not confirmed as worded

- **Attribution.** The text attributes the SESOI-passing result to "probabilities estimated from each query's 100-hit record and its nearest available family" alone.
  - The declared contrast that passed is the *selected configuration* (estimate **and** exchange weight 1) against the released objective: −0.0172 [−0.0221, −0.0118].
  - The estimate alone at the released exchange weight is an undeclared, exploratory contrast: −0.0128 [−0.0186, −0.0088]. Its interval reaches inside the SESOI.
  - The Discussion sentence (edit 9) does name "a lower exchange weight"; the abstract should match it.
- **Truth proxy.** "Reactions the source organism lacks" overstates the truth proxy. The measure is flux through gene-backed reactions *the source genome's draft does not support*, which the body text says.
- **Wording.** "Declared smallest effect" should read "declared smallest effect of interest".
- **Confirmed.** "as released it moves little flux and breaks no more ties than shuffled weights do" and "in multi-genome selections only": at n = 1 the gain is −0.0048 [−0.0113, +0.0063], and the population-weighted pooling B is inconclusive, −0.0065 [−0.0119, +0.0032].

## 2. Community kinetic caveat (edit 2): not confirmed as worded

- **The structural claim holds.** Sources: c2e31a0 `modelseedpy/core/msprobability.py` l.165–174 (md5 `8d40e57a…`), MSCommunity 03bd80c, INVESTIGATION.md N6.
  - When `not commkinetics` and there is more than one model, l.172 constructs `MSCommunity(mdlUtil.model, …)`.
  - That constructor runs `self.util = MSModelUtil(model, True, …)` (mscommsim.py:230). The second argument is `copy` (`MSModelUtil.__init__(self, model, copy=False, …)`, msmodelutl.py:372), so `model.copy()` is constrained.
  - `CommKineticPkg` is then built on `self.util.model` (mscommsim.py:291–293), i.e. on the copy.
  - `commkinetics=True` skips MSCommunity, and `build_from_species_models` takes `commkinetics` only in its signature (commhelper.py:45).
  - N6 measured the objective at 763.8 for both kinetics_coef = 1,000 and kinetics_coef = 1.
- **But the released function never reaches a solution on this path.** l.173–174 build `ElementUptakePkg` a second time, after the constructor, and raise `PackageError`.
  - `repro_c2e31a0.log` FP3b: "RAISED PackageError: Package with name ElementUptakePkg already in model!" on a two-member community.
  - WEIGHTED_OBJECTIVE_RESULTS.md l.198: "On c2e31a0 a two-member community also raised the `PackageError` … The branch runs it, but MSCommunity installs its kinetic rows on its own copy … That remains unfixed."
  - The no-effect measurement was therefore taken on b6857a0, or on c2e31a0 with the duplicate build removed.
  - Saying "neither … affects the solution" of c2e31a0 implies that the released community path runs. It does not.
- **Typography.** `kinCoef` is bare text; the equation two lines up uses `$\text{kinCoef}$`.
- **Confirmed.** "no result here uses the community path": WEIGHTED_OBJECTIVE_RESULTS.md l.172, "the community branch of the released function was not run".

## 3. Fix commits (edit 3): confirmed

- **The commits.** `git -C ModelSEEDpy-prfba log --oneline -3` lists 09856f5, b6857a0 and c2e31a0, on branch `fix-prfba-single-model`.
- **What 09856f5 does.** It adds `min_growth=1e-6`. If the floored objective's maximum is NaN (infeasible) or below 1e-6, it raises `ValueError`, because the 95% floor would be about 0 and the minimisation would "return the all-zero flux vector". The text matches this.
- **Notes, no change needed.**
  - It tests the *floored objective*, which is growth only when biomass is that objective. The preceding sentence explains that the floor goes on whatever objective the model carries.
  - It also refuses an infeasible model.
  - `git branch -vv` shows no upstream for `fix-prfba-single-model`, so both commits are local. "Committed" is literally true, but readers cannot fetch them. This was already the case for b6857a0.

## 4. The second post-hoc test (edit 4): numbers confirmed, one misstatement

Contrast values are from `summary.json` and `VERIFICATION.md` (Stage O, confirm half, pooling A, face midpoint unless stated). Where my own recompute applies, it agrees.

| Manuscript | Source | OK |
|---|---|---|
| "declared before it ran" | Design committed in b46598b, 2026-11-02 22:06:05 −0600. Its sha256 `27dc89a9…` equals `meta.json` `design_sha256`, and the file is unchanged since. The run went from 2026-11-03 08:50:28 to 10:27:57. | yes (see note a) |
| same 360 amplicons; 90 objectives | 360 rows ok; 5 P × 2 C × 3 G × 3 X = 90 | yes |
| 0.084 = rate at which gapfilled reactions are draft-supported | g_meas = 1021 / 12222 = 0.08354, from the tune half of the first test (0.0844 over all rows) | yes |
| log cost "normalised to the same values at p = 0 and p = 1" | max(0.01, ln max(p, 0.01) / ln 0.01): 1 at p = 0 and 0.01 at p = 1, as for the linear cost | yes |
| 167 / 193 | n_by_half tune 167, confirm 193 | yes (wording: R7a) |
| lowered by 0.0172 [0.0118, 0.0221] | −0.0172 [−0.0221, −0.0118]; recompute −0.01719 | yes |
| "past the smallest effect of interest, with both face endpoints lowered" | SESOI 0.01. Rule: "midpoint CI below 0 with mean ≤ −0.01, and both face endpoints' CIs below 0". Label: "lowers … by ≥ SESOI". The midpoint interval also lies wholly beyond 0.01. | yes (R1) |
| minimum 0.0076 [0.0047, 0.0133] | face min −0.0076 [−0.0133, −0.0047]; recompute −0.00762. The source adds: "below the SESOI". | yes, but R1 |
| beat shuffled by 0.0282 [0.0195, 0.0382] | −0.0282 [−0.0382, −0.0195] | yes (R2) |
| 69% against 23% of the oracle's headroom | headroom_confirm.A: 0.6854 (own cell) and 0.2305. Each is measured against its own cell's two-tier and oracle arms. In the released cell's frame the selected arm captures 0.933, so the stated comparison is the conservative one. | yes |
| exploratory: estimate alone 0.0128 [0.0088, 0.0186] | est/lin/G0/X100 − kn/lin/G0/X100 = −0.0128 [−0.0186, −0.0088]; recompute −0.01285 | yes |
| lower exchange weight a further 0.0043 [0.0017, 0.0061] | est/lin/G0/X1 − est/lin/G0/X100 = −0.0043 [−0.0061, −0.0017]. The two parts telescope to the headline. | yes |
| log worse by 0.0039 [0.0011, 0.0053] | est/log/G0/X1 − est/lin/G0/X1 = +0.0039 [+0.0011, +0.0053] | yes |
| "pricing gapfilled reactions at 0.084 or 0.2 … changed nothing" | G sets the gapfilled *probability*, so the price is 1 − p (0.916, 0.8). Gmeas − G0 = +0.000038 [−0.000003, +0.000105]. Over all 12 families and both halves the largest effect is 0.000202. | **no** (C6) |
| "a median 100% of them carry flux at the optimum" | First test, `arms.real.n_active_gapfilled / model.n_gapfilled`: median 1.000 over 360 rows (193 rows all active; mean 0.988; min 0.904; pooled 26,414 / 26,789 = 0.986) | yes |
| "the model needs them to grow" | Interpretation. It is reasonable: they are added for growth and carry flux even at full cost. | acceptable |
| 13.6 against 35.4; shuffled 2.1 | FVA levels: selected 13.6, shipped 35.4, est_shuffled at X1 2.1; contrast −21.89 [−24.65, −14.81] | yes |
| single-genome 0.0048 [−0.0063, 0.0113], 30 amplicons, after 0.0124 | n1 confirm −0.0048 [−0.0113, +0.0063] (30 rows, 20 groups); tune −0.0124 | yes |
| 96 of 193; without it 0.0164 [0.0100, 0.0236] | grp00000 holds 96 of 193; without it −0.0164 [−0.0236, −0.0100] | yes |
| held out from selection, not from design | Results §6 and VERIFICATION, "Design-level exposure" | yes (R13) |
| "The released objective is unchanged." | Nothing in ModelSEEDpy or the pipeline changed the objective | yes |

**Notes on edit 4.**

- (a) **What "declared before it ran" rests on.** It holds for the design: factors, split, selection rule, SESOI and controls. Two finer points rest on files, not commits; both are disclosed in the sources:
  - The interpretation rules (results §5) and the scripts are untracked, so they are fixed only by mtime (2026-11-02 23:05) and recorded hashes.
  - The estimator carried was E3, not the Stage E script's winner E4(E3). That choice was made at 08:50:15, after Stage E's results and 13 s before Stage O started.
- (b) **"Weight" for gapfilled reactions.** Elsewhere the manuscript calls p a "weight" ("five sets of weights", "an oracle weight of 1"), so "the weight of gapfilled reactions (0 …)" is internally consistent. But the same sentence uses "exchange weight" for a cost coefficient. "Gapfilled weight 0 and exchange weight 1" can therefore be read as free gapfilled flux, when p = 0 is in fact the full cost (R8). "Pricing … at 0.084" is wrong under either reading (C6).
- (c) **Omitted context (R2, R3).**
  - The shuffled comparator is itself worse than the released objective: est_shuffled − shipped = +0.0110 [+0.0031, +0.0222]. The 0.0282 margin is therefore partly a harmful comparator. The first test's paragraph disclosed the analogous fact.
  - Exchange weight 1 raises total exchange flux: +2.50 [+1.20, +3.39] units against a median of about 165 (it rises on 110 of 193 rows, by up to +18%). Which exchanges carry the increase was not recorded, and the source recommends testing on measured exchange data before shipping.
- (d) **The split.** It is the sha256 parity of the group id: 121 groups (167 amplicons) tune, 80 groups (193 amplicons) confirm, 201 groups in all (recomputed; 0 hash mismatches). 167 is not "one half of the amplicons", and 121 of 201 is not "half of the holdout groups" (R7a–c).

## 5–7. Scope statement (edits 5–7): confirmed

- The two tests, the "Apart from those tests" sentence, and "evaluated only for calibration, against an alternative estimate" all match the sources.
- One optional point, with no correction proposed: the second test also solved *variants* of eq:objective (log cost, exchange weights 10 and 1, gapfilled p). "That solve \Cref{eq:objective}" is loose but not misleading.

## 8. Neighbour-and-clade passage (edit 8): two errors

| Manuscript | Source | OK |
|---|---|---|
| hits among the first 100, weight exp(−(id₁ − id_j)/0.02) | `hits_used` 100; K = 100 means every stage-0 survivor of the 100 hits; τ = 0.02 (`weights.meta.tau`, `selected.json`) | yes |
| "coverage-qualified hits" | Stage 0 has three filters: coverage ≥ 0.90, \|SW − edlib identity\| ≤ 0.02, and identity ≥ 0.865. Hits are then deduplicated by genome. The manuscript calls the first two the "reliability gate". | imprecise (R4a) |
| shrunk toward the anchor's family prevalence among the fold's reference genomes, pseudo-count 2, order when the family is unnamed or has < 3 genomes | α = 2. Rule: "family → order … back off one rank while the anchor's clade has < 3 genomes in the fold reference … or its panel label is a placeholder". | yes, but omits the third case: no prior when the order also fails, 7,591 of 57,021 cells (13%) and 21 of the 360 Stage O amplicons (R4b) |
| nested CV chose the same settings in all five folds | `selected.json` E3 per_fold, 0–4: `E3\|K=100\|tau=0.02\|genome\|family\|alpha=2` | yes |
| no neighbour or clade member shares genus, species, taxon or holdout group with the source | `assertions.json` `required.f_*` and `i_*` are all 0 at d100, d050 and d025 | yes |
| universe Brier 0.0889 → 0.0578, −0.0311 [−0.0322, −0.0297] | E0 0.08891; E3 0.05779; contrast −0.03113 [−0.03215, −0.02972] | yes (sign style: R11) |
| **union Brier 0.209 → 0.130** | E0 0.2087; **E3 0.1267**. The 0.130 is **E1's** union Brier (0.1301), quoted by PFUNC_ESTIMATORS_RESULTS.md §1 for E1. | **no** (C1) |
| union AUROC 0.543 → 0.816 | 0.5426 → 0.8160 | yes |
| n = 1 AUROC 0.805 | 0.8054 | yes |
| recalibration lowers it by 0.0411; further gain 0.0408 [0.0384, 0.0424]; "about half" | E4(E0) − E0 union −0.0411 [−0.0466, −0.0373]; E3 − E4(E0) −0.0408 [−0.0424, −0.0384]; 0.0411 / 0.0820 = 50% | yes |
| prevalence ranking 0.739; estimate adds 0.067 [0.062, 0.072] | Stage E VERIFICATION: 0.739 [0.725, 0.751] at n = 1; E3 − prevalence +0.067 [0.062, 0.072] | yes |
| macro-F1 gain +0.032, mostly outside the union | E3 d_macro_F1_vs_union 0.0317 [0.0289, 0.0339]; "about three quarters … is support expansion" | yes |
| within the union +0.0085 [+0.0072, +0.0096], below the SESOI | +0.0085 [0.0072, 0.0096]; SESOI 0.02; "not exceeded" | yes |
| **"(+0.025 at selections of two or more)"** | n ≥ 2 within-union gain +0.0249 [0.0228, 0.0265], **"exceeded" (lower bound above 0.02)**. Placed after "below the smallest effect of interest", the parenthetical reads as the SESOI's value, and it hides that the n ≥ 2 gain passes it. | **no** (C2) |
| "not part of the released pipeline" | No tracked pipeline file outside `region_validation/` references the estimator; Stage O ran it from a study weights file | yes |
| (population) | Stage E scores all 57,021 cells (49,048 at n = 1), while the preceding sentence is about multi-reference cells | unstated (R5) |

## 9. Discussion (edit 9): contradiction with the unchanged opening

- **The contradiction.** The sentence still opens "…and \Cref{eq:objective} only in one post-hoc in-silico test on selections of the controlled benchmark … :" and then, in the same paragraph, the new text says "A second post-hoc test shows …" (C3).
- **Confirmed.** 0.0172 [0.0118, 0.0221] is the confirm half against the released objective. "Specific to the estimates" holds (selected − shuffled −0.0282 [−0.0382, −0.0195]; see R2), as does "only in multi-genome selections".
- **Wording.** "Chosen on half of the holdout groups … on the other half" (R7b). "More of the reference set's signal" should be the alignment record's: the estimate uses all 100 hits and a clade prior, not the selected reference set (R9).

## 10. Limitations (edit 10): confirmed

- **Confirmed.** "Two post-hoc in-silico tests on the same 300 multi-genome and 60 single-genome amplicons" holds. So does "tuned on a benchmark reference of about 7,000 genomes": 6,788 (fold 0) to 7,175 (folds 1–4), against about 315,000 in production. "Post hoc" and "not released" hold.
- **Wording.** Only "chosen on half of the holdout groups" needs changing (R7c).

## 11. Seven-reaction sentence (edit 11): confirmed

- **The reactions.** There are seven: rxn23850 (GAPN, the non-phosphorylating G3P dehydrogenase), rxn05746 and rxn08933–rxn08936 (five maltodextrin glucosidases), and rxn00014 (cytochrome c peroxidase). Each goes from `=` to `>`, and each is one step of an energy-generating loop when reversed (REPORT §1).
- **Draft leaks.** Leaking drafts fall from 865 to 95 of 1,181, with 770 stopping and 0 starting (`scan_summary.json`; REPORT §2 corrects the commit message's "849").
- **Through gapfilling, on all 200 genomes.**
  - Gating: `gated_after` 0 of 8,614 rows.
  - Growth: `growing_after` 2,362 of 2,362.
- **Medians (REPORT §4).**
  - Controlled pipeline: essentiality moves by at most 0.0005, growth by at most 0.0035, burden by at most 1 reaction.
  - As shipped: essentiality by at most 0.0024 and growth by at most 0.0016, but burden by up to **1.5** reactions (PICRUSt2, 84 → 82.5).
  - So "at most one reaction" is the controlled pipeline, which is the one the paper reports (R10).
- **Amendment status.** "Not part of any amendment; every figure stays on A-C13" matches REPORT §5: "The paper stays on A-C13".

## Check 3: consistency with the unchanged text (both files)

- **L1102, Discussion.** "only in one post-hoc in-silico test" contradicts the appended second-test sentence (C3).
- **L339.** "used as the weight of \Cref{eq:objective} in one post-hoc in-silico test": p_func was also the released arm of the second test, and that arm agrees with the first (shipped − two_tier on the confirm half −0.0056 [−0.0074, −0.0036], below the SESOI). "The first" is safer (R6).
- **L338.** "no number below should be read as evidence for or against the weighting" now sits directly above a passage that compares p_func with an alternative estimate. The tension is pre-existing, because the calibration sentence was already below, but the new passage sharpens it (R12).
- **L310.** "whose exploratory exchange-weight-1 arm": the description of the first test (L292) lists five weight sets and never introduces this arm (R13).
- **Abstract against Discussion.** The Discussion includes the exchange weight in the result and the abstract does not (C5).
- **No other stale "one test" statements remain.** The only other `one post-hoc` hits are L339 and L1102. L332, L336 and L1153 were updated by the edits.

## Check 4: the two variants

- **Before the edits.** `diff prFBA.tex prFBA_andrew_voice.tex` shows only lines 121–124 and 239, identically at HEAD and in the working tree. HEAD plus the 11 edits reproduces both working files exactly.
- **After the corrections.** All 23 proposed old strings occur exactly once in each file, and applying them leaves the variant diff unchanged.

## Check 5: compile

- **Working files** (scratch copy, `build/`). `latexmk -pdf` exits 0 for both variants: 90 pages each, and `grep -c undefined` on each log gives 0.
- **Corrected files** (all 23 corrections, `build2/`). Exit 0, 90 pages each, 0 undefined.

## Corrections (exact strings; each old string occurs exactly once in each .tex file, also when applied in this order)

The same pairs, in machine-readable form, are in `corrections.py`. Run `python corrections.py <dir>` to check them, and add `--apply` to write them; use a scratch copy only. Old strings are the current working-tree text. C6 uses "probability", which is correct whether or not R8 is adopted.

### Required (errors, contradictions, misattributions)

**C1**: Stage E union Brier is E1's value, not the estimate's (E3 0.1267).

old:
```latex
and over the union's functions from 0.209 to 0.130, and it raises
```
new:
```latex
and over the union's functions from 0.209 to 0.127, and it raises
```

**C2**: the n>=2 within-union F1 gain exceeds the 0.02 SESOI; the parenthetical reads as the SESOI.

old:
```latex
below the smallest effect of interest ($+0.025$ at selections of two or more).
```
new:
```latex
below the smallest effect of interest of 0.02, although at selections of two or more it is $+0.0249$ [$+0.0228$, $+0.0265$], above it.
```

**C3**: Discussion still says 'only in one post-hoc in-silico test' before 'A second post-hoc test'.

old:
```latex
and \Cref{eq:objective} only in one post-hoc in-silico test on selections of the controlled benchmark (the \nameref{sec:prfba} section): the weighting cuts
```
new:
```latex
and \Cref{eq:objective} only in two post-hoc in-silico tests on selections of the controlled benchmark (the \nameref{sec:prfba} section). In the first, the weighting cuts
```

**C4**: community kinetics: on c2e31a0 the community path raises before any solve (FP3b); the private-copy result is on the fix branch, which leaves it unfixed.

old:
```latex
In the released \texttt{prFBA} function (ModelSEEDpy c2e31a0), however, the community package installs this row on a private copy of the community model rather than on the model that is solved, so neither kinCoef nor the \texttt{commkinetics} option affects the solution; no result here uses the community path.
```
new:
```latex
In the released \texttt{prFBA} function (ModelSEEDpy c2e31a0), however, the community path, too, fails at the carbon-uptake constraint built twice (below), and where that is corrected the community package installs this row on a private copy of the community model rather than on the model that is solved, so neither $\text{kinCoef}$ nor the \texttt{commkinetics} option affects the solution; the fix described below does not change this, and no result here uses the community path.
```

**C5**: abstract: the SESOI-passing configuration includes exchange weight 1; truth is the source genome's draft, not the organism.

old:
```latex
whereas with probabilities estimated from each query's 100-hit record and its nearest available family it moves flux off reactions the source organism lacks by more than the declared smallest effect, in multi-genome selections only).
```
new:
```latex
whereas with probabilities estimated from each query's 100-hit record and its nearest available family, and a lower exchange weight, it moves flux off reactions the source genome's draft does not support by more than the declared smallest effect of interest, in multi-genome selections only).
```

**C6**: G is a probability (price 1 - p), and the measured effect is up to 0.0002, not nothing.

old:
```latex
and pricing gapfilled reactions at 0.084 or 0.2 instead of 0 changed nothing, because
```
new:
```latex
and giving gapfilled reactions a probability of 0.084 or 0.2 instead of 0 changed truth-absent flux by at most 0.0002, because
```

### Recommended (overstatement by omission, imprecision, stale wording)

**R1**: the face-minimum gain is below the SESOI (source: 'below the SESOI').

old:
```latex
with both face endpoints lowered (minimum 0.0076 [0.0047, 0.0133]);
```
new:
```latex
with both face endpoints lowered (the minimum by 0.0076 [0.0047, 0.0133], below the smallest effect of interest);
```

**R2**: the shuffled comparator is itself worse than the released objective (+0.0110 [+0.0031, +0.0222]).

old:
```latex
it beat the same costs on shuffled estimates by 0.0282 [0.0195, 0.0382],
```
new:
```latex
it beat the same costs on shuffled estimates by 0.0282 [0.0195, 0.0382], partly because shuffled estimates raise truth-absent flux (by 0.0110 [0.0031, 0.0222] against the released objective),
```

**R3**: omitted caution: exchange weight 1 raises total exchange flux, composition not recorded.

old:
```latex
so they were held out from the selection but not from the design.
```
new:
```latex
so they were held out from the selection but not from the design. The lower exchange weight also raises total exchange flux (by 2.5 [1.2, 3.4] units, against a median of about 165), and which exchanges carry the increase was not recorded.
```

**R4a**: the filter is the reliability gate plus the family floor, one hit per genome, not coverage alone.

old:
```latex
the identity-weighted share of the coverage-qualified hits among a query's first 100 that carry it
```
new:
```latex
the identity-weighted share of the hits among a query's first 100 that pass the reliability gate and the family floor (one per genome) and carry it
```

**R4b**: the clade rule has a third case: no prior when the order also fails (7,591 of 57,021 cells).

old:
```latex
(pseudo-count 2; the order when the family is unnamed or has fewer than three reference genomes);
```
new:
```latex
(pseudo-count 2; the order when the family is unnamed or has fewer than three reference genomes, and no prior, for 13\% of cells, when the order fails too);
```

**R5**: Stage E numbers are over all 57,021 scored cells, not the multi-reference cells of the preceding sentence.

old:
```latex
Against $p_{\text{func}}$ it lowers the Brier score over the whole EC vocabulary
```
new:
```latex
Against $p_{\text{func}}$, over all 57{,}021 scored cells (49{,}048 of them single-genome selections), it lowers the Brier score over the whole EC vocabulary
```

**R6**: p_func was also the released arm of the second test.

old:
```latex
used as the weight of \Cref{eq:objective} in one post-hoc in-silico test, it lowers
```
new:
```latex
used as the weight of \Cref{eq:objective} in the first post-hoc in-silico test, it lowers
```

**R7a**: 167 of 360 is not half the amplicons; groups split 121/80 by hash.

old:
```latex
A configuration was selected on the holdout groups of one half of the amplicons (167) and tested on the other (193).
```
new:
```latex
The holdout groups were split in two by a hash of their identifiers, and a configuration was selected on one half (121 groups, 167 amplicons) and tested on the other (80 groups, 193 amplicons).
```

**R7b**: Discussion: give the realised split.

old:
```latex
chosen on half of the holdout groups, truth-absent flux fell by 0.0172 [0.0118, 0.0221] against the released objective on the other half,
```
new:
```latex
chosen on 121 of the 201 holdout groups, truth-absent flux fell by 0.0172 [0.0118, 0.0221] against the released objective on the other 80,
```

**R7c**: Limitations: give the realised split.

old:
```latex
and the second's configuration chosen on half of the holdout groups,
```
new:
```latex
and the second's configuration chosen on 121 of the 201 holdout groups,
```

**R8a**: 'weight 0' for gapfilled reactions reads as free flux; it is p = 0, cost 1.

old:
```latex
the weight of gapfilled reactions (0 as released, 0.084,
```
new:
```latex
the probability given to gapfilled reactions (0 as released, so full cost; 0.084,
```

**R8b**: same.

old:
```latex
gapfilled weight 0 and exchange weight 1
```
new:
```latex
gapfilled probability 0 and exchange weight 1
```

**R9**: the estimate uses the whole 100-hit record and a clade prior, not only the reference set.

old:
```latex
the objective can carry more of the reference set's signal:
```
new:
```latex
the objective can carry more of the alignment record's signal:
```

**R10**: burden 'at most one reaction' is the controlled pipeline; as shipped it is 1.5.

old:
```latex
and moved the method arms' medians by at most 0.004 and their gapfill burden by at most one reaction.
```
new:
```latex
and moved the method arms' medians by at most 0.004 and, in the controlled pipeline, their gapfill burden by at most one reaction (1.5 as shipped).
```

**R11**: sign convention: every other reduction in the edits is a positive magnitude.

old:
```latex
($-0.0311$ [$-0.0322$, $-0.0297$])
```
new:
```latex
(by 0.0311 [0.0297, 0.0322])
```

**R13**: the first test's description lists five weight sets; the exchange-weight-1 arm is never introduced.

old:
```latex
whose exploratory exchange-weight-1 arm had been seen when the factors were chosen,
```
new:
```latex
whose exploratory exchange-weight-1 arm, not reported above, had been seen when the factors were chosen,
```

**R12**: pre-existing, now sharper: the new passage below compares the weighting with an alternative.

old:
```latex
and no number below should be read as evidence for or against the weighting.
```
new:
```latex
and no selector-arm number below should be read as evidence for or against the weighting.
```
