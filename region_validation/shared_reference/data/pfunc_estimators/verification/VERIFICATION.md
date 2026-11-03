# Stage E (function-probability estimators): independent verification

Verified 2026-11-03 against `scripts/pfunc_estimators.py` (sha256 `47d5b6c13217…`, matching the run), its outputs in `data/pfunc_estimators/`, and the draft `PFUNC_ESTIMATORS_RESULTS.md`. The design (`PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md`, "Stage E", sha256 `27dc89a9b340…`) is unchanged since b46598b (committed 22:06). Every Stage E output and dev directory in the author's scratchpad is timestamped 22:34 or later.

**Bottom line.** Every number I checked reproduces. No leakage was found. The nested selection is genuinely nested. There are no numerical errors in the draft. The problems are in the framing:

- Three statements about Stage O are stale.
- One phrase ("declared rule") is wrong.
- One "within 1e-6" claim has a 1.3e-6 exception.
- One reference size is misquoted.
- Three headline comparisons need a baseline next to them, or they overstate what the hit list adds:
  - recalibration vs new information;
  - the n = 1 AUROC against a no-information ranker;
  - how much of the decision-metric F1 gain comes from outside the union.

## Method

I wrote my own code from the design and the raw inputs. It does not import `pfunc_estimators`, `calibrate_pfunc` or `select_references`. Stage 0 was re-implemented from `select_references.select_representatives`:

- coverage and |SW − edlib| are rounded to 4 decimals before comparison;
- identity is compared raw;
- survivors are sorted stably by (−identity, −align_score), then deduplicated by genome.

The code, all in this directory:

- `vctx.py`: panel, EC universe through norm_ec, truth matrix, cells filter.
- `p0_sanity.py`: population and fold-reference checks.
- `p1_stream.py`: ijson stream of the 5 d100 `predictions.json` and hits files. It builds the 100-hit and 200-hit lists in 44 s.
- `p2_pass.py`: 72 variants per cell. These are the declared grids (species variants are omitted after checking they are a no-op) plus out-of-grid probes and baselines. It also runs the nested selection, the composites, my own out-of-fold isotonic (10,000 bins, against the author's 4,000), and leakage on every cell.
- `p3_stage_o.py`: all 360 Stage O entries of all three weight files.
- `p4_analyze.py`: metrics and the cluster bootstrap. It uses the full **2,000 draws**, seed 1729, and score_folds' draw over the stratum's sorted holdout-group keys, so the CIs are directly comparable.
- `p5_compare.py`: side-by-side with `scores.json`.
- `p6_extra.py`: per-cell AUROC and truncation.
- `p7_k10.py`: the K = 10 F1 discrepancy.
- `p8_leaky.py`: counterfactual leaky prior.

Results are in `p0_sanity.json`, `p2_summary.json`, `p3_stage_o.json`, `p4_results.json`, `p5_compare.json`, `p6_extra.json` and `p8_leaky.json`.

The population checks reproduce: 57,021 cells; 8,459 sources; 2,292 groups; n = 1/2/3/≥4 = 49,048 / 5,693 / 1,313 / 967.

On every cell:

- n_truth equals `cells.csv.gz`;
- union F1 equals `cells.f1` within 1e-3;
- the union of the reps' ECs equals `predictions.json` `ec`;
- the 100-hit and 200-hit tie sets equal the reps.

Each fold's `_cache_ec_subset.json` EC rows equal the panel truth. Its key set equals the fold plan's reference list (6,788 / 7,175 / 7,175 / 7,175 / 7,175) and equals the panel minus that fold's test genomes. That set shares no genome, genus or holdout group with the fold's test set.

## Verdicts

### Claim 1: CONFIRMED

| quantity | mine | study / calibration | abs diff |
|---|---|---|---|
| E0 universe Brier | 0.0889144246 | 0.08891442 | 4.6e-9 |
| E0 union macro F1 | 0.7672147746 | 0.7672148 | 2.5e-8 |
| E0 log-loss / union micro F1 / union-support Brier / p ≥ 0.5 macro F1 | 0.41445486 / 0.78426500 / 0.20868144 / 0.76791519 | 0.4144549 / 0.784265 / 0.2086814 / 0.7679152 | ≤ 4e-8 |
| E1 K20 τ0.01 genome, **200 hits**, universe Brier | 0.0627383198 | 0.06273832 | 2.3e-10 |
| E1 K20 τ0.01 genome, **100 hits**, universe Brier | 0.0627383362 | 0.06273834 | ~0 |
| E1 200-hit K = 5/10/50/100 universe Brier | 0.06744084 / 0.06434223 / 0.06184415 / 0.06161656 | same | ≤ 5e-9 |
| E1 200-hit K = 10, q ≥ 0.5 macro F1 | 0.7897647 (= the study) | calibration 0.7897634 | **1.3e-6** |

The K = 10 F1 discrepancy belongs to calibrate_pfunc, not to this study. At K = 10, 18,388 EC points sit at q = 0.5 exactly. A dot product (mine and the study's) gives 0.7897647; calibrate's cumulative sums give 0.7897634 (`p7_k10.py`). So the q ≥ 0.5 rule at exact ties depends on summation order. The effect is immaterial, but it falsifies the draft's "all within 1e-6" (see correction C13).

The design's literal E1 target "0.06270" is a rounding typo for 0.06274. The script asserted the correct calibration value. No change is needed in the results file.

### Claim 2: CONFIRMED

Nested CV pools the outer folds over all 57,021 cells, contrasted against E0 with holdout-group CIs from 2,000 draws. Mine and the study's agree to every printed digit:

| | mine | study |
|---|---|---|
| E3 − E0 universe Brier | −0.03112881 [−0.03215481, −0.02971987] | identical |
| E3 − E0 universe log-loss | −0.2126739 [−0.2184941, −0.2052910] | −0.2126739 [−0.2184941, −0.2052910] |
| E3 − E0 union-support Brier | −0.08196353 [−0.08700278, −0.07890226] | identical |
| E3 union WLS slope | 0.903493 [0.891836, 0.922438] | 0.9034928 [0.8918363, 0.9224377] |
| E3 / E0 union AUROC | 0.816007 / 0.542571 | 0.816007 / 0.5425707 |
| E3 macro dF1 vs union | +0.0316921 [0.0288705, 0.0339469] | identical |
| n = 1 (49,048 cells): E3 − E0 universe Brier | −0.0327677 [−0.0338055, −0.0312031] | identical |
| n = 1: E3 union AUROC | 0.805398 [0.791890, 0.817423] | 0.8053983 [0.7918904, 0.8174226] |

Over E0, E1, E2 and E3, every metric (universe Brier, log-loss, union Brier, WLS slope, AUROC, macro and micro dF1) in all six strata agrees with `scores.json` to ≤ 5e-8. This holds for estimates and for contrasts against E0. CI bounds agree to ≤ 9e-6; the largest differences are in the binned AUROC replicates.

E4(E0) and E4(E3) were refitted with my own out-of-fold isotonic, using a different binning. They agree to ≤ 6e-5 (estimates) and ≤ 1.6e-4 (CI bounds). Examples: E4(E3) − E3 universe Brier −0.0003 [−0.0007, 0.0001]; union Brier −0.0005 [−0.0024, 0.0010]; macro dF1 −0.0042 [−0.0066, −0.0017], against the study's [−0.0067, −0.0017].

Not re-derived:

- E4(E1) and E4(E2);
- ECE;
- the universe OLS slope;
- the genus and source-genome sensitivity CIs;
- all d050 and d025 rows.

### Claim 3: CONFIRMED

My nested selection picks the study's variant in every outer fold:

| family | fold 0 | fold 1 | fold 2 | fold 3 | fold 4 |
|---|---|---|---|---|---|
| E1 | K50 τ.05 | K100 τ.02 | K50 τ.05 | K100 τ.02 | K50 τ.05 |
| E2 | family α4 | family α4 | family α8 | family α8 | family α4 |
| E3 | K100 τ.02 genome family α2 | same | same | same | same |

- **Log-loss values.** My training-fold universe log-loss for all 49 distinct grid variants × 5 outer folds matches `selected.json` `train_universe_logloss_by_variant_and_outer_fold` to ≤ 5e-8. The outer-fold log-loss matches too.
- **E3 is a clear optimum in all five folds.** It beats the runner-up (α = 1) by 3.0–3.7e-4.
- **E2's α = 4 vs α = 8 is essentially a tie** in folds 3 and 4 (margins 2.7e-6 and 2.1e-6).
- **The family start beats the genus start in every fold.** The margin is 0.014–0.017 for E2 and 0.0011–0.0014 for E3.
- **Species dedup changes 0 of 50,842 records.** Raw panel species labels are unique (8,872). One pair collides after the Candidatus normalization: "Candidatus Sumerlaeia bacterium" and "Sumerlaeia bacterium", both in grp00960. That pair never alters a list.

### Claim 4: CONFIRMED (no leakage)

I checked every cell, both clade starts, and the 100-hit, 200-hit and rep lists. All counts are zero:

- the source is never a neighbour, rep or clade member;
- no neighbour, rep or clade member shares the source's genus, panel species (raw or normalised), taxon_id or holdout group;
- every one is in the cell's own fold reference;
- every clade is drawn from the cell's fold;
- the anchor's genus never equals the source's.

All 200 stored hits of every record are fold-reference genomes. Clade rank counts reproduce the README: family start 40,280 family / 9,150 order / 7,591 none; genus start 19,626 / 20,911 / 9,057 / 7,427.

### Claim 5: CONFIRMED

All 360 entries were recomputed, in all three files: `stage_o_weights.json`, `stage_o_weights_E1.json` (K50 τ.05) and `stage_o_weights_E3.json`. The 40 random keys (seed 20261103) are reported separately in `p3_stage_o.json`.

- **What matches exactly:** neighbour ids and order (first 100 hits → stage 0 → rank → genome dedup → top K); clade rank and member set; α (0 exactly when the clade is null); reps; md5; fold.
- **Weights match to ≤ 5e-8.** The files store weights to 7 significant digits, so a consumer's q differs from the estimator's by ≤ 4.5e-8.
- **The draft's counts reproduce:** family 305 / order 34 / none 21; median 100 neighbours; median 18 members among the 339 with a clade.
- **The E3 file is the main file without the map.** Its `amplicons` are byte-for-byte equal to `stage_o_weights.json`'s, with `isotonic: null`.
- **The population join is clean.** md5, fold, source, group and reps all match the scored cells.

## The requested checks (a)–(f)

**(a) Clade prior or hits drawing on the held-out genus: none.**

- **Members.** Clade members come only from the cell's fold reference, which contains no test-fold genome, genus or group.
- **Counterfactual.** I recomputed E3 with the prior taken over the *full panel*, which brings the held-out genomes back in. Universe Brier would be 0.05565 and union Brier 0.1224. The reported values are 0.05778561 and 0.1267, which match my fold-reference prior to 1e-8. That leak is therefore ruled out; in 27,315 cells the source itself would have been a member.
- **Hits.** All 200 hits per record are fold-reference genomes.
- **Family prior (a non-finding).** It does include the source's own *family*: the anchor's family equals the source's in 44.7% of cells, and its order in 66.4%. That is the design's "nearest available clade". The genus is blocked, not the family, so this is not leakage.

**(b) Nested selection excluded the outer fold: yes.**

- **Selection.** My selection uses only the four training folds' log-loss, and it reproduces the choices and the training-loss values exactly. The study's composite rows by fold equal those choices.
- **E4 maps.** My out-of-fold maps for E4(E0) and E4(E3), each fitted on the other four folds, reproduce E4 within 2e-4.
- **Two caveats that are not leakage:**
  - The training folds' records were searched against references that include the outer fold's genomes, as references only. That is inherent to the benchmark.
  - The CIs hold the selected variants and the fitted maps fixed. They therefore exclude selection and refitting noise, which is small for E3 because all five folds chose the same variant.

**(c) Unsupported ECs.**

- **Scoring is consistent.** Every estimator is scored on all 2,235 ECs against the same truth, and log-loss clips every point to [0.01, 0.99].
- **Exact decomposition.** E0's universe squared error splits into union-support squared error plus the count of truth ECs outside the union, with a residual of −9e-10.
- **Who puts mass outside the union.** E2 and E3 put α·π/(den + α) there, and E4 puts r(0) ≈ 0.008 on every q = 0.
- **The asymmetry this creates.** Universe metrics reward support expansion that Stage O cannot use; the draft says so. The same expansion also drives most of the decision-metric gain, which the draft does not say (see the findings below).

**(d) Fairness of the n = 1 comparison.**

- **AUROC.** E0's 0.5 at n = 1 is by construction (q ≡ 1; CI [0.5, 0.5]). E3's 0.805 is real, but most of the gain over 0.5 needs no hit-list information:
  - ranking the anchor's ECs by their prevalence in the fold reference, the same score for every amplicon, gives a pooled AUROC of 0.739 [0.725, 0.751] at n = 1;
  - the family prior alone gives 0.722.
- **What the hit list adds.** E3 − prevalence AUROC is +0.067 [0.062, 0.072] at n = 1 (+0.073 [0.068, 0.078] over all cells). As per-cell mean AUROC: E3 0.845, prevalence 0.772, E0 0.5.
- **Pooled AUROC also mixes in between-cell differences.** E4(E0), whose q at n = 1 is a per-fold constant, has a pooled AUROC of 0.479 at n = 1.
- **Brier.** At n = 1, E0's union-support Brier (0.2134) is 1 − precision. The recalibrated constant, E4(E0) with q = r(1) ≈ 0.79, gives 0.1680. E3's information gain at n = 1 is therefore E3 − E4(E0) = −0.0415 [−0.0431, −0.0391] on union-support Brier, out of the −0.0869 reported against E0. On universe Brier it is −0.0206 [−0.0219, −0.0187], against −0.0328.

**(e) Grid edges.** I ran an out-of-grid probe. It is exploratory and was not used for any choice.

- **E1.** Unweighted or τ = 0.1, 0.2 (K = 50, 100) is worse on every outer fold's training folds, by +0.0018 to +0.0022 universe log-loss over the grid optimum.
- **E3.** Unweighted or τ = 0.05, 0.1 (family, α = 1, 2, 4) is worse in every fold, by +0.0014 to +0.0019. Pooled over the outer folds, E3 at τ = 0.05 vs 0.02 gives universe Brier +0.0011 [0.0008, 0.0013] and union Brier +0.0016 [0.0004, 0.0025]. Unweighted E3 gives +0.0030 [0.0026, 0.0033] and +0.0052 [0.0032, 0.0065].
- **Implication.** The edge hides no better configuration at these values. The E3 optimum is bracketed by 0.01 and 0.05, and E1's by 0.01 and 0.1. Identity weighting matters: flat weights are clearly worse.
- **The clade prior lowers the optimal τ.** E3's optimum (0.02) sits below E1's modal choice (0.05).
- **What this does to "E3 vs E1".** Part of the "E3 vs E1" difference is that τ difference. At matched hits (E1 at K = 100, τ = 0.02):

  | metric | E3 − E1 at matched τ |
  |---|---|
  | universe Brier | −0.0012 [−0.0014, −0.0011] |
  | union Brier | −0.0030 [−0.0039, −0.0025] |
  | union slope | 0.88 → 0.90 |
  | macro dF1 | +0.0007 [0.0000, 0.0012] |

  The draft's 0.84 → 0.90 compares against E1's composite, which uses τ = 0.05 in three folds; that variant's slope is 0.826.

**(f) The isotonic map for Stage O is in-sample, as disclosed.**

- **My refit reproduces the file's map.** Fitting on all five folds with the carried-forward E3 gives r(0) = 0.008323, r(0.5) = 0.4831 and r(1) = 0.943375. The file has 0.008323, 0.4831 and 0.943375. The maximum gap on a 2,001-point grid is 0.0057, which comes from the different binning. The out-of-fold maps give r(1) = 0.936–0.947.
- **The in-sample share is negligible.** The 360 amplicons are 0.63% of cells.
- **The material points:**
  - the map caps p at 0.943, which puts the cost floor at 0.057, not 0.01;
  - it is fitted per EC but would be applied at role level.
- **Stage O does not use it.** Stage O was launched with `stage_o_weights_E3.json`, which has no map: `data/objective_tuning/meta.json` gives `args.weights` = that file, `weights.meta.estimator` = "E3" and `isotonic` = None. The draft's §1, §4 and §6 still describe the E4 file as the Stage O file (corrections C4, C10–C12, C18).

## Further findings (overstatement checks)

1. **About half of the improvement over k/n is recalibration, which k/n itself can have.**
   - E4(E0) − E0 = −0.0113 universe / −0.0411 union.
   - E3 − E4(E0), the information gain, is −0.0198 [−0.0211, −0.0180] universe Brier, −0.0408 [−0.0424, −0.0384] union-support Brier, +0.291 [0.265, 0.301] union AUROC and +0.0332 [0.0303, 0.0356] macro dF1.
   - The design's "beats k/n" rule is met by E4(E0) alone, whose AUROC is 0.525 and macro dF1 −0.0015. So "all seven alternatives beat k/n" shows better calibration; it does not by itself show more information.
2. **Most of the decision-metric F1 gain comes from outside the union.**
   - E3's macro dF1 vs the union is +0.0317 [0.0289, 0.0339].
   - Restricted to the union support (thresholding only the ECs the draft contains), it is **+0.0085 [0.0072, 0.0096]**. That is below the 0.02 SESOI, so "not exceeded".
   - The rest comes from 48.7 ECs per cell outside the union with q ≥ 0.5, of which 29.0 are true.
   - By stratum: at n = 1, 0.0317 against +0.0058 [0.0046, 0.0072], not exceeded; at n ≥ 2, 0.0317 against +0.0249 [0.0228, 0.0265], exceeded (lower bound above 0.02).
   - E1 behaves the same way: 0.0295 against +0.0066.
   - The design defines the decision metric over the universe, so the reported number is per design. But the draft's support-expansion caveat claims "Universe metrics include this" and leaves the F1 rule out, although that rule is computed over the universe too.
3. **The draft's statements check out on everything else:**
   - section 1's E1, E2, E3 and E4(E3) numbers;
   - all of section 3's E0–E3 table cells and the n = 1 and n ≥ 2 tables;
   - section 4's selection;
   - the 7,591 (13.3%) no-clade cells;
   - the 43.5% (22,107 of 50,842) of d100 records whose last coverage-qualified hit among the first 100 is still above the family floor;
   - the 88 required assertions, all of which pass.

## Proposed corrections to `PFUNC_ESTIMATORS_RESULTS.md` (old → new)

**C1** (line 3, banner):
- OLD: `**DRAFT, not yet independently re-derived.** The design (`PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md`, "Order of work" item 4) requires an independent re-derivation before any manuscript or deck text changes. That has not been done. Everything below comes from one script and the checks it runs on itself, plus one spot check.`
- NEW: `**Independently re-derived on 2026-11-03** (separate code from the raw inputs; `pfunc_verify/VERIFICATION.md`). Every d100 Brier, log-loss, slope, AUROC and F1 number for E0–E3 below reproduced to within 5e-8, with the same 2,000 bootstrap draws. E4(E0) and E4(E3) reproduced to within 2e-4 with a separately fitted isotonic map. Not re-derived: the density rows, ECE, the universe OLS slope, E4(E1), E4(E2), and the genus and source-genome intervals.`

**C2** (line 5):
- OLD: `` `stage_o_weights.json`, `stage_o_weights_E1.json` and `run.log`. ``
- NEW: `` `stage_o_weights.json`, `stage_o_weights_E1.json` and `run.log`. `stage_o_weights_E3.json`, written 2026-11-03 08:50, holds the same 360 entries without the isotonic map; it is the file Stage O was launched with. ``

**C3** (§1):
- OLD: `the nested-CV contrasts against E0 have CIs below 0 for both universe Brier and union-support Brier.`
- NEW: `the nested-CV contrasts against E0 have CIs below 0 for both universe Brier and union-support Brier. Recalibrating k/n alone (E4(E0): union AUROC 0.525, macro dF1 −0.0015) also meets the rule. So the rule shows better calibration; it does not on its own show more information. The contrast with E4(E0), below, isolates the information.`

**C4** (§1):
- OLD: `It cuts universe Brier from 0.0889 to 0.0593 and union-support Brier from 0.209 to 0.130.`
- NEW: `It cuts universe Brier from 0.0889 to 0.0593 and union-support Brier from 0.209 to 0.130. About half of the union-support gain is calibration that k/n can have without new information: the out-of-fold isotonic map of k/n alone, E4(E0), reaches 0.1675 (−0.0411). Against E4(E0), E3 changes universe Brier by −0.0198 [−0.0211, −0.0180], union-support Brier by −0.0408 [−0.0424, −0.0384] and union AUROC by +0.291 [0.265, 0.301].`

**C5** (§1, last bullet of the E3-vs-E1 list):
- OLD: `- the union-support calibration slope from 0.84 to 0.90.`
- NEW: `- the union-support calibration slope from 0.84 to 0.90. Most of that slope change is E1's flatter τ = 0.05 in three folds. Against E1 at the same hits (K = 100, τ = 0.02), the clade prior changes universe Brier by −0.0012 [−0.0014, −0.0011], union-support Brier by −0.0030 [−0.0039, −0.0025], the slope from 0.88 to 0.90, and macro dF1 by +0.0007 [0.0000, 0.0012].`

**C6** (§1):
- OLD: `**Winner by the declared rule: E4(E3), the out-of-fold isotonic recalibration of E3.**`
- NEW: `**Winner by the script's winner rule: E4(E3), the out-of-fold isotonic recalibration of E3.** The design names no winner rule. The script takes the lowest nested-CV universe log-loss among the estimators that beat k/n.`

**C7** (§1):
- OLD: `` In practice E3 and E4(E3) are the same estimator on the scores that decide. `stage_o_weights.json` carries E3 in full plus the isotonic map. Drop the map and Stage O runs E3. ``
- NEW: `` In practice E3 and E4(E3) are the same estimator on the scores that decide. Stage O therefore runs E3. It was launched with `stage_o_weights_E3.json`, which holds the same 360 entries as `stage_o_weights.json` without the isotonic map (`data/objective_tuning/meta.json`, `args.weights`). ``

**C8** (§1, n = 1 paragraph):
- OLD: `There every union EC gets p = 1, so k/n cannot rank anything (AUROC 0.5). E3 ranks them, with a union-support AUROC of 0.805, and lowers union-support Brier from 0.213 to 0.127.`
- NEW: `There every union EC gets p = 1, so k/n cannot rank anything (AUROC 0.5 by construction). E3 ranks them, with a union-support AUROC of 0.805, and lowers union-support Brier from 0.213 to 0.127. Two baselines put this in proportion. Ranking the anchor's ECs by their prevalence in the fold reference, which carries no information about the amplicon, already gives an AUROC of 0.739 [0.725, 0.751]; E3 adds +0.067 [0.062, 0.072] to that. The recalibrated constant (E4(E0), q ≈ 0.79) already gives a union-support Brier of 0.168. E3's information gain at n = 1 is therefore −0.0415 [−0.0431, −0.0391] of the −0.087.`

**C9** (§1, decision-metric bullets):
- OLD: `- E3: +0.032 [0.029, 0.034], exceeded, with both bounds above the SESOI;`
- NEW: `- E3: +0.032 [0.029, 0.034], exceeded, with both bounds above the SESOI. About three quarters of this is support expansion: E3 puts q ≥ 0.5 on 48.7 ECs per cell outside the union, 29.0 of them true. Thresholding only the union's ECs gives +0.0085 [0.0072, 0.0096], which is not exceeded (+0.0058 at n = 1; +0.0249 [0.0228, 0.0265] at n ≥ 2, which is exceeded; E1: +0.0066);`

**C10** (§4, grid edges):
- OLD: `- **Grid edges.** E1's τ = 0.05 and E3's τ = 0.02 are the largest values in their declared grids, so the optimum may lie at larger τ. This study does not look further.`
- NEW: `- **Grid edges.** E1's τ = 0.05 and E3's τ = 0.02 are the largest values in their declared grids. The independent re-derivation probed beyond them; the probe is exploratory and informed no choice. On every outer fold's training folds, E1 at τ = 0.1, 0.2 or unweighted (K = 50, 100) has universe log-loss 0.0018–0.0022 above E1's grid optimum, and E3 at τ = 0.05, 0.1 or unweighted (family, α = 1, 2, 4) has 0.0014–0.0019 above E3's. Pooled, E3 at τ = 0.05 is worse than at 0.02: universe Brier +0.0011 [0.0008, 0.0013], union-support Brier +0.0016 [0.0004, 0.0025]. The E3 optimum is therefore bracketed by 0.01 and 0.05, and E1's by 0.01 and 0.1.`

**C11** (§4, Stage O file):
- OLD: `- **Stage O file** (`stage_o_weights.json`, 360 amplicons): estimator E4, K = 100, τ = 0.02, dedup genome, clade start family (back-off to order), α = 2. It includes a 747-knot isotonic map, with r(0) = 0.008, r(0.5) = 0.48 and r(1) = 0.94; the out-of-fold maps give r(1) between 0.936 and 0.947.`
- NEW: `` - **Stage O files** (360 amplicons). `stage_o_weights.json` holds estimator E4: K = 100, τ = 0.02, dedup genome, clade start family (back-off to order), α = 2. It includes a 747-knot isotonic map with r(0) = 0.008, r(0.5) = 0.48 and r(1) = 0.943; the out-of-fold maps give r(1) between 0.936 and 0.947. The map caps p at 0.943, so no reaction reaches the cost floor, and it is in-sample. Stage O was therefore launched with `stage_o_weights_E3.json`: the same entries without the map (estimator E3). ``

**C12** (§4 sub-bullets):
- OLD: `  - Each amplicon has a median 100 neighbours and a median 18 clade members.`
- NEW: `  - Each amplicon has a median 100 neighbours; the 339 with a clade have a median 18 members.`
- OLD: `  - The map is fitted on all five folds, so it is in-sample for these 360 benchmark amplicons (0.6% of cells).`
- NEW: `  - The map is fitted on all five folds, so it is in-sample for these 360 benchmark amplicons (0.6% of cells). It does not reach Stage O, which runs the E3 file.`

**C13** (§5):
- OLD: `` - **Reproductions.** All agree with `data/calibration/scores.json` to within 1e-6 or better: ``
- NEW: `` - **Reproductions.** All agree with `data/calibration/scores.json` to within 1e-6 or better, with one exception. The 200-hit K = 10 q ≥ 0.5 macro F1 is 0.7897647 here and 0.7897634 there (1.3e-6): 18,388 points sit at q = 0.5 exactly, and calibrate_pfunc's cumulative sums put some of them just below 0.5. ``

**C14** (§5, append to the spot-check bullet):
- OLD: `All 30 matched.`
- NEW: `All 30 matched. The independent re-derivation then recomputed all 360 entries of all three files. Neighbour ids and order, clade rank and members, α and reps match exactly. Weights match to 5e-8, because the files store 7 significant digits.`

**C15** (§6, transfer):
- OLD: `The benchmark reference is about 6,800 species-dereplicated genomes per fold.`
- NEW: `The benchmark reference is 6,788 (fold 0) to 7,175 (folds 1–4) species-dereplicated genomes per fold.`

**C16** (§6, clade coverage):
- OLD: `Their anchor's order is a placeholder or has fewer than 3 fold-reference genomes. For those cells E3 equals E1.`
- NEW: `Their anchor's order is a placeholder (4,551 cells) or has fewer than 3 fold-reference genomes (3,040). For those cells E3 equals E1 at K = 100, τ = 0.02.`

**C17** (§6, truncation):
- OLD: `the 100th coverage-qualified hit is still above the family floor`
- NEW: `the last coverage-qualified hit among the first 100 is still above the family floor`

**C18** (§6, E4 map):
- OLD: `Stage O would apply it to role-level q on the union support, which this study does not validate.`
- NEW: `Stage O would apply it to role-level q on the union support, which this study does not validate. Stage O was therefore launched without it (`stage_o_weights_E3.json`).`

**C19** (§6, bootstrap details):
- OLD: `the binned and exact point estimates agree to within about 1e-5.`
- NEW: `the binned and exact point estimates agree to within about 1e-5. The intervals hold the selected variants and the fitted isotonic maps fixed. They therefore exclude selection and refitting noise, which is small for E3 because all five folds chose the same variant.`

**Not for the draft but worth fixing (script-generated README).** The README's "Cells without an admissible clade up to order … by reason" dict counts the *genus*-start "none" set:

- what it says: {order_placeholder 4,387, order_lt3 3,040} = 7,427;
- for the carried-forward family start (7,591 cells) the split is 4,551 / 3,040.

The README's "Stage O weights" section also predates `stage_o_weights_E3.json`.

## Files written

All in `/tmp/claude-1000/-home-freiburger-Documents-prFBA/0c6f5bc6-8224-5a34-9b20-7d5643487901/scratchpad/pfunc_verify/`:

- this `VERIFICATION.md`;
- the code: `vctx.py` and `p0_sanity.py`–`p8_leaky.py`;
- the results: `p0_sanity.json`, `p2_summary.json`, `p3_stage_o.json`, `p4_results.json`, `p5_compare.json`, `p6_extra.json`, `p8_leaky.json`;
- the caches: `cache/stream_f{0..4}.pkl`, `cache/passA.npz`, `cache/passB.pkl`, `cache/union_quantities.npz`.

No repository file was modified and nothing was committed.
