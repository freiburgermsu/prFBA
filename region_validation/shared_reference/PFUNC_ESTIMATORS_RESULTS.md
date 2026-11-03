# Function-probability estimators (Stage E): results

**Independently re-derived on 2026-11-03** (separate code from the raw inputs; `data/pfunc_estimators/verification/VERIFICATION.md`). Every d100 Brier, log-loss, slope, AUROC and F1 number for E0–E3 below reproduced to within 5e-8, with the same 2,000 bootstrap draws. E4(E0) and E4(E3) reproduced to within 2e-4 with a separately fitted isotonic map. Not re-derived: the density rows, ECE, the universe OLS slope, E4(E1), E4(E2), and the genus and source-genome intervals.

Run 2026-11-02, 22:59 to 23:06 (systemd unit, 12 workers), `scripts/pfunc_estimators.py` (sha256 `47d5b6c13217…`). Design: `PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md`, section "Stage E", committed in b46598b. Outputs are in `data/pfunc_estimators/`: `scores.json`, `assertions.json`, `selected.json`, `README.md` (all tables), `stage_o_weights.json`, `stage_o_weights_E1.json` and `run.log`. `stage_o_weights_E3.json`, written 2026-11-03 08:50, holds the same 360 entries without the isotonic map; it is the file Stage O was launched with. The study is **post hoc**. It was declared after the calibration study's exploratory top-k analysis had been seen, and the estimator grids were informed by that analysis. Nested selection keeps the tuning out of the reported numbers. It does not undo the choice of what to tune.

Intervals are 95% percentile cluster bootstraps over holdout groups (2,000 draws, seed 1729, score_folds' draw). All estimators share the same draws, so the contrasts are paired. Genus and source-genome intervals are in `scores.json` and the README. They are narrower than the group intervals and change no conclusion below.

## 1. The answer

**Can p be better than k/n?** Yes, by a wide margin, and every family does it. All seven alternatives to k/n beat it by the design's rule (E1–E3, and E4 of each family including k/n itself): the nested-CV contrasts against E0 have CIs below 0 for both universe Brier and union-support Brier. Recalibrating k/n alone (E4(E0): union AUROC 0.525, macro dF1 −0.0015) also meets the rule. So the rule shows better calibration; it does not on its own show more information. The contrast with E4(E0), below, isolates the information. Most of the gain comes from using more of the hit list than the exact-tie set. E1 is identity-weighted prevalence over the stage-0 survivors. Its nested composite used the top 50 at τ = 0.05 in three outer folds (the carried-forward configuration) and the top 100 at τ = 0.02 in the other two. It cuts universe Brier from 0.0889 to 0.0593 and union-support Brier from 0.209 to 0.130. About half of the union-support gain is calibration that k/n can have without new information: the out-of-fold isotonic map of k/n alone, E4(E0), reaches 0.1675 (−0.0411). Against E4(E0), E3 changes universe Brier by −0.0198 [−0.0211, −0.0180], union-support Brier by −0.0408 [−0.0424, −0.0384] and union AUROC by +0.291 [0.265, 0.301].

**Does the clade prior add anything?** A little, and the gain is distinguishable from 0. E3 is the hit list shrunk toward the anchor's family prior (K = 100, τ = 0.02, α = 2). Compared with E1 it changes:

- universe Brier by −0.0016 [−0.0018, −0.0013];
- union-support Brier by −0.0034 [−0.0040, −0.0028];
- the union-support calibration slope from 0.84 to 0.90. Most of that slope change is E1's flatter τ = 0.05 in three folds. Against E1 at the same hits (K = 100, τ = 0.02), the clade prior changes universe Brier by −0.0012 [−0.0014, −0.0011], union-support Brier by −0.0030 [−0.0039, −0.0025], the slope from 0.88 to 0.90, and macro dF1 by +0.0007 [0.0000, 0.0012].

The clade prior on its own (E2, k/n shrunk to the family prior) is the weakest family. Its union-support Brier is 0.149, against 0.130 for E1.

**Winner by the script's winner rule: E4(E3), the out-of-fold isotonic recalibration of E3.** The design names no winner rule. The script takes the lowest nested-CV universe log-loss among the estimators that beat k/n. It ranks first by universe log-loss, universe Brier and union Brier alike. However, it is **not distinguishable from E3** on either Brier score:

- universe Brier: E4(E3) − E3 = −0.0003 [−0.0007, +0.0001];
- union Brier: −0.0005 [−0.0024, +0.0010].

It is better on log-loss, −0.0041 [−0.0068, −0.0020], and on calibration slope: 0.94 against 0.90 on the union support, and 1.00 against 0.94 on the universe. It loses macro F1 at q ≥ 0.5, by −0.0042 [−0.0067, −0.0017]. In practice E3 and E4(E3) are the same estimator on the scores that decide. Stage O therefore runs E3. It was launched with `stage_o_weights_E3.json`, which holds the same 360 entries as `stage_o_weights.json` without the isotonic map (`data/objective_tuning/meta.json`, `args.weights`).

**Where k/n fails most** is where prFBA usually runs: one selected genome, 86% of scored cells. There every union EC gets p = 1, so k/n cannot rank anything (AUROC 0.5 by construction). E3 ranks them, with a union-support AUROC of 0.805, and lowers union-support Brier from 0.213 to 0.127. Two baselines put this in proportion. Ranking the anchor's ECs by their prevalence in the fold reference, which carries no information about the amplicon, already gives an AUROC of 0.739 [0.725, 0.751]; E3 adds +0.067 [0.062, 0.072] to that. The recalibrated constant (E4(E0), q ≈ 0.79) already gives a union-support Brier of 0.168. E3's information gain at n = 1 is therefore −0.0415 [−0.0431, −0.0391] of the −0.087. The gain at n ≥ 2 is smaller but clear: union-support Brier −0.057 [−0.059, −0.053]. There the calibration slope goes from 0.78 to 0.97.

**Decision metric (q ≥ 0.5 against the union, SESOI 0.02, score_folds' rule).** Thresholding gives these macro F1 changes:

- E3: +0.032 [0.029, 0.034], exceeded, with both bounds above the SESOI. About three quarters of this is support expansion: E3 puts q ≥ 0.5 on 48.7 ECs per cell outside the union, 29.0 of them true. Thresholding only the union's ECs gives +0.0085 [0.0072, 0.0096], which is not exceeded (+0.0058 at n = 1; +0.0249 [0.0228, 0.0265] at n ≥ 2, which is exceeded; E1: +0.0066);
- E4(E3): +0.027 [0.025, 0.031], exceeded;
- E1: +0.030 [0.027, 0.032], exceeded;
- E2: +0.014 [0.011, 0.017], not exceeded (non-zero, below the SESOI).

For comparison, k/n's own p ≥ 0.5 rule changes it by only +0.0007.

**At lower reference density** the carried-forward E3 still beats k/n by the same margin, without re-tuning. Universe Brier changes by −0.033 at d050 and −0.033 at d025; union Brier by −0.087 and −0.089. The isotonic step gains slightly more at lower density: E4(E3*) − E3* union Brier is −0.0028 [−0.0052, −0.0009] at d025.

## 2. What was run

- **Population.** As calibrate_pfunc defined it: the controlled genus-blocked benchmark, `cells.csv.gz` arm `prfba_exact_tie`, space ec, outcome scored, density 1.0. That is 57,021 cells, 8,459 source genomes and 2,292 holdout groups. By selection size: n = 1 49,048; n = 2 5,693; n = 3 1,313; n ≥ 4 967. Truth is the source genome's EC set on the 2,235-EC universe. Reps come from the fold run's `predictions.json`.
- **Records.** Each fold's own hits, searched against that fold's reference only. Only the first 100 of the 200 stored hits were used, in align-score rank order (asserted). Then come select_references stage 0, ranking by identity and then align score, and genome dedup. "Top K" is the first K survivors; K = 100 means every survivor of the 100 hits. The 200-hit lists were kept only for the reproduction check.
- **Estimators.** The design's E0–E4, 88 tuned variants in all:
  - E1: 30 variants (K, τ, dedup);
  - E2: 10 (clade start genus or family, α);
  - E3: 48;
  - E4: isotonic of each family's selected variant.

  The anchor is `reps[0]` (prFBA's own anchor). Clades use panel_eligible lineage over the fold reference (`_cache_ec_subset.json` keys). A clade backs off genus → family → order when it has fewer than 3 reference genomes or its label is a panel placeholder (`unclassified_*`).
- **Nested selection.** For each outer fold, each family's variant minimises the pooled universe log-loss (ε 0.01) of the other four folds. E4's map for a fold is fitted on the other four folds, using the variant selected for that fold.
- **Metrics.** As designed. Universe Brier, log-loss and decile ECE. Union-support Brier, WLS slope and AUROC; the union support is E0's union for every estimator. Macro and micro F1 of q ≥ 0.5 against the union. Each is given overall and for n = 1, 2, 3, ≥ 4 and ≥ 2.

## 3. Results (nested CV, d100)

All 57,021 cells. Each cell shows the estimate and, below it, the contrast against E0.

| | E0 (k/n) | E1 | E2 | E3 | E4(E3) |
|---|---|---|---|---|---|
| universe Brier | 0.0889 [0.0870, 0.0911] | 0.0593 [0.0583, 0.0610] | 0.0672 [0.0652, 0.0704] | 0.0578 [0.0567, 0.0596] | 0.0575 [0.0564, 0.0590] |
| Δ vs E0 | | −0.0296 [−0.0306, −0.0282] | −0.0217 [−0.0233, −0.0197] | −0.0311 [−0.0322, −0.0297] | −0.0314 [−0.0324, −0.0303] |
| universe log-loss | 0.4145 [0.4057, 0.4247] | 0.2113 [0.2072, 0.2178] | 0.2571 [0.2459, 0.2740] | 0.2018 [0.1975, 0.2087] | 0.1977 [0.1943, 0.2026] |
| Δ vs E0 | | −0.2031 [−0.2092, −0.1958] | −0.1573 [−0.1672, −0.1449] | −0.2127 [−0.2185, −0.2053] | −0.2168 [−0.2231, −0.2107] |
| universe ECE (deciles) | 0.0872 | 0.0193 | 0.0380 | 0.0155 | 0.0005 |
| union-support Brier | 0.2087 [0.2001, 0.2214] | 0.1301 [0.1243, 0.1384] | 0.1494 [0.1408, 0.1614] | 0.1267 [0.1209, 0.1351] | 0.1262 [0.1218, 0.1328] |
| Δ vs E0 | | −0.0785 [−0.0836, −0.0754] | −0.0593 [−0.0615, −0.0567] | −0.0820 [−0.0870, −0.0789] | −0.0825 [−0.0893, −0.0781] |
| union WLS slope | 0.699 [0.669, 0.716] | 0.841 [0.828, 0.862] | 0.801 [0.767, 0.827] | 0.903 [0.892, 0.922] | 0.940 [0.928, 0.970] |
| union AUROC | 0.543 [0.535, 0.549] | 0.799 [0.783, 0.812] | 0.734 [0.715, 0.751] | 0.816 [0.802, 0.828] | 0.815 [0.805, 0.825] |
| macro dF1 vs union (q ≥ 0.5) | 0.0007 [0.0005, 0.0009] | 0.0295 [0.0270, 0.0320] | 0.0142 [0.0108, 0.0165] | 0.0317 [0.0289, 0.0339] | 0.0274 [0.0252, 0.0310] |
| micro dF1 vs union | 0.0008 | 0.0275 [0.0252, 0.0294] | 0.0139 [0.0110, 0.0159] | 0.0300 [0.0275, 0.0318] | 0.0259 [0.0238, 0.0291] |
| ECs with q > 0, per cell | 480 | 972 | 770 | 1,065 | 2,235 |

ECE is a sum of absolute bin errors and is biased upward under resampling. For E4(E3) its percentile interval, [0.0006, 0.0073], lies above the estimate, so it is given without an interval. E4(E0), E4(E1) and E4(E2) are in the README. E4 of k/n alone also beats k/n: universe Brier −0.0113, union Brier −0.0411. The isotonic step improves calibration. It does not improve ranking.

**n = 1** (49,048 cells; contrast against E0):

| | E1 | E3 | E4(E3) |
|---|---|---|---|
| universe Brier | −0.0311 [−0.0321, −0.0295] | −0.0328 [−0.0338, −0.0312] | −0.0331 [−0.0341, −0.0318] |
| universe log-loss | −0.2138 | −0.2244 [−0.2304, −0.2156] | −0.2290 [−0.2356, −0.2220] |
| union Brier (E0 0.2134) | −0.0830 [−0.0889, −0.0795] | −0.0869 [−0.0925, −0.0835] | −0.0875 [−0.0949, −0.0828] |
| union AUROC (E0 0.500) | 0.785 | 0.805 [0.792, 0.817] | 0.804 [0.795, 0.814] |
| union WLS slope (E0 undefined) | 0.818 | 0.887 [0.873, 0.910] | 0.922 [0.906, 0.957] |
| macro dF1 vs union | +0.0295 [0.0271, 0.0325] | +0.0317 [0.0290, 0.0344] | +0.0274 [0.0246, 0.0317] |

**n ≥ 2** (7,973 cells; contrast against E0):

| | E1 | E3 | E4(E3) |
|---|---|---|---|
| universe Brier (E0 0.0755) | −0.0202 [−0.0219, −0.0175] | −0.0210 [−0.0227, −0.0185] | −0.0211 [−0.0224, −0.0189] |
| universe log-loss | −0.1372 | −0.1405 [−0.1495, −0.1262] | −0.1418 [−0.1491, −0.1296] |
| union Brier (E0 0.1846) | −0.0556 [−0.0580, −0.0526] | −0.0567 [−0.0589, −0.0533] | −0.0567 [−0.0602, −0.0542] |
| union WLS slope (E0 0.778) | 0.929 | 0.965 [0.955, 0.978] | 1.004 [0.996, 1.026] |
| union AUROC (E0 0.705) | 0.854 | 0.857 [0.845, 0.865] | 0.855 [0.846, 0.863] |
| macro dF1 vs union (E0 0.0050) | 0.0293 [0.0234, 0.0329] | 0.0317 [0.0264, 0.0351] | 0.0277 [0.0245, 0.0296] |

In the separate buckets, E3 − E0 universe Brier is −0.0227 at n = 2, −0.0173 at n = 3 and −0.0162 at n ≥ 4. Union Brier is −0.0616, −0.0467 and −0.0439. Every one has a CI below 0. The clade prior's extra gain over E1 sits mostly at n = 1, where E3 − E1 union Brier is −0.0039. At n ≥ 2 it is −0.0011 [−0.0015, −0.0002].

**Density sensitivity** (carried-forward configurations, no re-tuning; E4 refitted out of fold within each density):

| | d050 (54,808 cells) | d025 (52,260 cells) |
|---|---|---|
| E0 universe / union Brier | 0.0934 / 0.2188 | 0.0963 / 0.2260 |
| E3* − E0, universe Brier | −0.0328 [−0.0340, −0.0310] | −0.0330 [−0.0344, −0.0309] |
| E3* − E0, union Brier | −0.0872 [−0.0904, −0.0847] | −0.0891 [−0.0918, −0.0851] |
| E3* − E1*, union Brier | −0.0040 | −0.0042 |
| E4(E3*) − E3*, union Brier | −0.0016 [−0.0036, 0.0000] | −0.0028 [−0.0052, −0.0009] |

## 4. Selected hyperparameters (`selected.json`)

| family | outer folds 0 / 1 / 2 / 3 / 4 | carried forward |
|---|---|---|
| E1 | K50 τ0.05 / K100 τ0.02 / K50 τ0.05 / K100 τ0.02 / K50 τ0.05 (genome) | K = 50, τ = 0.05, genome |
| E2 | family α4 / α4 / α8 / α8 / α4 | family, α = 4 |
| E3 | K100 τ0.02 genome family α2, all five folds | K = 100, τ = 0.02, genome, family, α = 2 |

- The family clade start beat the genus start in every fold, for both E2 and E3. Genome dedup was chosen only by the tie rule: the species variants are identical to the genome variants, because species dedup changed 0 of 50,842 records (the panel holds one genome per species).
- **Grid edges.** E1's τ = 0.05 and E3's τ = 0.02 are the largest values in their declared grids. The independent re-derivation probed beyond them; the probe is exploratory and informed no choice. On every outer fold's training folds, E1 at τ = 0.1, 0.2 or unweighted (K = 50, 100) has universe log-loss 0.0018–0.0022 above E1's grid optimum, and E3 at τ = 0.05, 0.1 or unweighted (family, α = 1, 2, 4) has 0.0014–0.0019 above E3's. Pooled, E3 at τ = 0.05 is worse than at 0.02: universe Brier +0.0011 [0.0008, 0.0013], union-support Brier +0.0016 [0.0004, 0.0025]. The E3 optimum is therefore bracketed by 0.01 and 0.05, and E1's by 0.01 and 0.1.
- **Stage O files** (360 amplicons). `stage_o_weights.json` holds estimator E4: K = 100, τ = 0.02, dedup genome, clade start family (back-off to order), α = 2. It includes a 747-knot isotonic map with r(0) = 0.008, r(0.5) = 0.48 and r(1) = 0.943; the out-of-fold maps give r(1) between 0.936 and 0.947. The map caps p at 0.943, so no reaction reaches the cost floor, and it is in-sample. Stage O was therefore launched with `stage_o_weights_E3.json`: the same entries without the map (estimator E3).
  - Clades: 305 amplicons have a family-level clade, 34 order-level, 21 none.
  - Each amplicon has a median 100 neighbours; the 339 with a clade have a median 18 members.
  - `stage_o_weights_E1.json` holds the best E1 (K = 50, τ = 0.05, genome).
  - The map is fitted on all five folds, so it is in-sample for these 360 benchmark amplicons (0.6% of cells). It does not reach Stage O, which runs the E3 file.

## 5. Checks

- **Reproductions.** All agree with `data/calibration/scores.json` to within 1e-6 or better, with one exception. The 200-hit K = 10 q ≥ 0.5 macro F1 is 0.7897647 here and 0.7897634 there (1.3e-6): 18,388 points sit at q = 0.5 exactly, and calibrate_pfunc's cumulative sums put some of them just below 0.5.
  - E0 universe Brier 0.08891442 and union macro F1 0.7672148;
  - E0 log-loss, union micro F1, union-support Brier and p ≥ 0.5 macro F1;
  - E1 at 200 hits, τ = 0.01, K = 5/10/20/50/100; the required K = 20 value is 0.06273832, exact;
  - E0 universe Brier at d050 and d025.

  At 100 hits, E1 K = 20 gives 0.06273834.
- **88 required assertions, all zero or true.**
  - Leakage, on every cell, for every neighbour list and both clade starts: the source is never a neighbour or clade member. No neighbour or member shares its genus, species, taxon_id or holdout group. All of them are in the fold reference, and the anchor's genus is never the source's.
  - The joins match calibrate_pfunc's: union equals the predictions' EC set, F1 and n_truth equal `cells.csv.gz`, the EC cache equals the truth matrix, and the 200-hit tie set equals the reps.
  - The 360 population amplicons match the scored cells.
  - Pass B equals pass A.
  - A consumer recomputation from `stage_o_weights.json` equals the estimator.
- **Spot check.** Separate scratch code recomputed 30 Stage O entries directly from the raw hits, caches and panel: neighbour ids and weights, tie set equal to the reps, clade rank and members, and leakage. All 30 matched. The independent re-derivation then recomputed all 360 entries of all three files. Neighbour ids and order, clade rank and members, α and reps match exactly. Weights match to 5e-8, because the files store 7 significant digits.
- **Determinism.** Two full runs gave identical estimates.

## 6. Caveats

- **Post hoc** (see the top). The decision rule held for every family, so the headline does not hang on one tuned choice. Even so, this is not a pre-registered result.
- **Transfer is the main open question.** The benchmark reference is 6,788 (fold 0) to 7,175 (folds 1–4) species-dereplicated genomes per fold. Production searches about 315,000 genomes, where the first 100 hits of a well-sequenced species can all be that one species. E1 and E3 would then average one species' strains rather than a neighbourhood. Nothing here measures that:
  - The species-dedup variant was meant to address it, but it is a no-op on this panel.
  - Density sensitivity thins the reference but adds no conspecific redundancy, so it speaks to the question only partly.
- **Cross-genus transfer.** The source's genus is never in the reference, so the clade prior is the nearest available clade, never the source's own. Absolute levels, and probably the gains, differ when a congeneric or conspecific reference exists.
- **Support expansion.** E1 and E3 give q > 0 to about twice as many ECs as the union, and E4 to all 2,235. Universe metrics include this. Stage O uses only the union support, so for Stage O the union-support rows are the relevant ones. On those rows E3 and E4(E3) are tied.
- **Clade coverage.** 7,591 cells (13%) have no admissible clade up to order with the family start. Their anchor's order is a placeholder (4,551 cells) or has fewer than 3 fold-reference genomes (3,040). For those cells E3 equals E1 at K = 100, τ = 0.02. The schema's ranks stop at order, so the back-off was not extended.
- **Truncation at 100 hits** leaves the tie set and the K = 20 numbers unchanged. But in 43% of d100 records the last coverage-qualified hit among the first 100 is still above the family floor (28% at d050, 11% at d025). For those records, K = 100 is cut short by the record depth.
- **E4 map in Stage O.** The map was fitted on EC-level q over the universe. Stage O would apply it to role-level q on the union support, which this study does not validate. Stage O was therefore launched without it (`stage_o_weights_E3.json`). The in-tree consumer (`scripts/objective_tuning.py`, read only, not modified) matches the file on three points. It reads the per-amplicon `alpha`, which is 0 for the 21 amplicons with a null clade, so they get no prior term. It applies `meta.isotonic` by clipped linear interpolation, as this script does. It uses `meta.estimator` ("E4") only as a label.
- **Bootstrap details.** AUROC replicates use 1,001-bin histograms; the binned and exact point estimates agree to within about 1e-5. The intervals hold the selected variants and the fitted isotonic maps fixed. They therefore exclude selection and refitting noise, which is small for E3 because all five folds chose the same variant. The holdout-group intervals are dominated by the largest group, as in the calibration study.

## 7. Runtime and files

- **Runtime.** A full run takes 411 s. Streaming the 15 hits and predictions files takes 45 s, pass A over 94 variants at three densities 155 s, and pass B plus bootstrap analysis about 170 s.
- **Memory.** Peak RSS was 7.9 GB in the parent and 9.9 GB in the largest child; the child figure counts pages shared with the parent through fork. The unit's cgroup peak was 14 to 21 GB, including page cache.
- **Files.** The script is `scripts/pfunc_estimators.py`. The outputs are in `data/pfunc_estimators/`; that directory is gitignored like the other `data/` outputs, so committing it needs `git add -f`. The stream cache (scratch only) is `/tmp/claude-1000/-home-freiburger-Documents-prFBA/0c6f5bc6-8224-5a34-9b20-7d5643487901/scratchpad/pfunc_est/cache_full/`.
