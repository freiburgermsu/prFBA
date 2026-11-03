# Stage O (objective tuning): independent verification

Re-derived 2026-11-03 with my own code from the raw rows. I imported none of the author's scripts. The bootstrap follows the stated convention: per stratum, sorted group keys, `default_rng(seed).integers(0, K, (2000, K))`, ratio of sums, seed 1729 + 1000 × index (n2, n3, n4+, n1). Pooled intervals are percentiles of the weighted sum of the strata's draws.

- **Inputs:** `data/objective_tuning/{per_amplicon.jsonl, meta.json, summary.json, OBJECTIVE_TUNING_RESULTS.filled.md}` (the 11:27 version, with the FVA section); `data/objective_tuning_fva/per_amplicon.jsonl`, read only after unit otune-fva had exited with EXIT=0 (193 rows); `data/weighted_objective/{per_amplicon.jsonl, population.json, recomputed_contrasts.json}`.
- **Weights:** A is 2474/663/410 over 3547 (n2 0.69749, n3 0.18692, n4+ 0.11559). B is n1 0.868, with the remaining 0.132 split in A's proportions.
- No repository file was modified, nothing was committed, and no unit was touched.

## Bottom line

- **Numbers.** Every number in the filled document that I checked reproduces: 972 confirm-half means and CI bounds with a maximum absolute difference of 0.0, and all 36 guard-table rows at the printed precision. The selection, rankings, headroom, checks table and degeneracy section also reproduce.
- **Leakage.** No confirm-half row reaches the selection, and g_meas comes from tune rows only.
- **Robustness.** The headline (selected − shipped, pooling A, −0.0172 [−0.0221, −0.0118]) survives:
  - leaving out any one group: the mean stays in [−0.0174, −0.0162];
  - removing the dominant confirm group: −0.0164 [−0.0236, −0.0100];
  - other seeds, a joint bootstrap and 20,000 draws;
  - excluding the rows whose rebuilt model differs from the prior run.
- **What I would change** is wording and missing context, not numbers:
  1. **n = 1.** At n = 1, where prFBA usually runs, nothing is detected on the confirm half, and the tune-half n1 effect did not replicate.
  2. **X1 mechanism.** The X1 component works by uptake substitution (exchange flux up, every internal flux class down), not by re-routing among reactions. About a third of its TAF drop is lower total gene-backed flux against the fixed D; two thirds is a precision gain.
  3. **taf_gb.** The "as returned (taf_gb)" column is not divided by D. It is 1 − flux precision, so it duplicates the precision guard.
  4. **Rebuilt LPs.** Three rows, not two, rebuild a different LP from the prior run.
  5. **Design exposure.** The confirm half is held out from selection, not from the design. The X1 level followed the prior study's exploratory exchange-weight-1 arm, which was run on all 360 amplicons.
  6. **Dominant group.** One holdout group holds half of the confirm rows.

## Claim by claim

### 1. Split: CONFIRMED
- I recomputed the hash from `group_id` for all 360 rows. It matches the recorded `half` field on every row: 167 tune and 193 confirm.
- There are 201 groups: 121 tune and 80 confirm. No group spans the halves, and every `group_id` equals the one in population.json.
- No genus, source genome or amplicon md5 occurs in both halves.
- Rows (groups) by stratum:
  - tune: n2 46 (40), n3 49 (40), n4+ 42 (33), n1 30 (28);
  - confirm: n2 54 (26), n3 51 (26), n4+ 58 (21), n1 30 (20).
- **Structural asymmetry, not in the document.** Group `grp00000` (fold 0; 96 distinct source genomes from 57 genera, led by Streptococcus with 18) holds **96 of the 193 confirm rows**. That is 28/54 n2, 22/51 n3, 35/58 n4+ and 11/30 n1, or 43–60% of each stratum. The largest tune group has 8 rows. All 96 fold-0 rows are confirm rows.

### 2. Eligibility: CONFIRMED (36 of 36), but the guard never bound
- My ratio-of-medians table reproduces the document's guard table exactly. Every candidate is optimal on all 167 tune rows, with biomass at 0.95 × max.
- The largest deviation of any candidate's median ratio from 1 is 0.45% for exchange flux and 3.0% for gene-backed flux, against the 10% guard. Medians of per-row ratios deviate by at most 0.02% and 0.85%.
- Per row, for the selected arm against the shipped arm (min / 5% / median / 95% / max):

| half | ratio | min / 5% / median / 95% / max | rows outside ±10% |
|---|---|---|---|
| tune | gene-backed flux | 0.645 / 0.838 / 1.000 / 1.116 / 1.453 | 29 / 167 (17%); 23% outside ±5% |
| tune | exchange flux | 1.000 / 1.000 / 1.000 / 1.065 / 1.119 | 1 / 167 |
| confirm | gene-backed flux | 0.800 / 0.930 / 1.001 / 1.073 / 1.272 | 12 / 193 (6%) |
| confirm | exchange flux | 1.000 / 1.000 / 1.000 / 1.060 / 1.184 | 4 / 193 |

- Exchange flux never falls below the shipped arm's (minimum ratio 1.00000000). It rises on 95/167 tune rows and 110/193 confirm rows. The median rise is 0.03 and 0.014 flux units, the 95th percentile 13.7 and 9.4, and the maximum 19.8 and 30.0. The shipped arm's median exchange flux is about 165 units.
- So the median guard says nothing about large per-row changes in the internal flux distribution.
- **Uptake composition is not recoverable.** Neither run records per-exchange fluxes. The arm fields hold only the total Σ|v_ex|, which counts secretion as well as uptake, and the FVA pass adds only gene-backed FVA counts. Whether X1 changes *which* nutrients are taken up cannot be answered from the rows. See the mechanism section for what the totals do show.

### 3. Selection: CONFIRMED, with one imprecision in the claim's wording
- **Selected arm.** est/lin/G0/X1 has pooled A 0.119648, against the shipped arm's 0.136160. Its pooled B is 0.175319, against 0.188290. The tie set within 0.002 is {est/lin/G0/X1}, with 2 departures from the shipped configuration.
- **Next family.** est/lin/G0/X10 is 0.122160, a margin of 0.0025. est/lin/G0/X100 is 0.123364, 0.0037 away and with 1 departure. A tie rule wider than 0.0037 would therefore have picked X100.
- **Sensitivities.** The joint (naive) selection over all 36 configurations agrees. So would selection by pooling B, where G0, Gmeas and G0.2 at X1 tie.
- **Log never wins.** At every (P, G, X), log ≥ lin on the tune half. On the confirm half, est/log/G0/X1 − est/lin/G0/X1 = +0.0039 [+0.0011, +0.0053].
- **The G factor never decides.** G0 is chosen in all 12 families, and my counts of G-varying rows equal the document's table. The brief's parenthetical, "gapfilled flux does not vary across G levels on the rows that would decide it", is **imprecise**:
  - Gapfilled flux *does* vary on 2 to 20 n ≥ 2 tune rows per family, including 3 for the selected family (2 n2, 1 n3).
  - On those rows G0 is lowest, or within 0.002 of the lowest. In est/lin/*/X10, G0 is lowest by 0.0045. In est/lin/*/X1, Gmeas is lower than G0 by 7e-9, a tie broken toward G0.
  - Over all tune rows, G moves pooled A by at most 0.0002 in any family.
- **Stability (new).** I resampled the tune groups (2,000 draws, per stratum and also jointly) and re-selected each time. est wins 100% of draws. est/lin/G0/X1 wins 93.3% (per stratum) and 91.8% (joint); est/lin/G0/X100 wins the rest.
- **Tune-half headroom.** 0.677 at the selected cell, 0.908 in the shipped cell's frame, and 0.259 for the shipped arm. All confirmed.
- **No winner's curse under pooling A.** Tune selected − shipped is −0.0165 [−0.0242, −0.0100] (A), against −0.0172 on the confirm half.

### 4. Confirmation: CONFIRMED (exact)
All means and CIs below equal summary.json to the last bit.

| contrast (confirm half) | pooling A, face mid | other values |
|---|---|---|
| selected − shipped | −0.0172 [−0.0221, −0.0118] | face min −0.0076 [−0.0133, −0.0047]; face max −0.0268 [−0.0323, −0.0178]; B −0.0065 [−0.0119, +0.0032]; n1 stratum −0.0048 [−0.0113, +0.0063] |
| selected − est_shuffled/lin/G0/X1 | −0.0282 [−0.0382, −0.0195] | |
| selected − two_tier/lin/G0/X1 | −0.0222 [−0.0268, −0.0159] | |
| oracle − two_tier at the selected cell | −0.0324 [−0.0372, −0.0242] | |
| shipped − two_tier at the shipped cell | −0.0056 [−0.0074, −0.0036] | |

- **Headroom captured** (ratio of pooled-A means): selected 0.6854 at its own cell; shipped 0.2305. In the shipped cell's frame the selected arm captures 0.933; the document does not print this for the confirm half (tune: 0.908).
- **Stronger agreement with the prior study than "in kind".**
  - Over all 360 Stage O rows, shipped − two_tier is **−0.0061 [−0.0097, −0.0025]**, identical to the prior study's.
  - The prior study's own rows restricted to the confirm half give −0.0056 [−0.0074, −0.0036], identical to Stage O's.
  - The prior tune half gives −0.0066 [−0.0128, +0.0005].
  - Stage O's kn/lin/G0/X1 − two_tier/lin/G0/X1 over all rows (−0.0086 [−0.0123, −0.0045]) equals the prior study's exploratory real_ex1 − two_tier_ex1. The two arm pairs agree to within 1e-9 on 356/360 rows.
- **Guards** (selected − shipped, A):
  - exchange flux +2.5045 [+1.2036, +3.3873];
  - gene-backed flux +1.1801 [−12.9610, +7.9732];
  - flux precision +0.0149 [+0.0100, +0.0203];
  - gapfilled flux −0.2501 [−0.7399, −0.0900].
- **Signs per row** (selected − shipped, face mid, negative / all): n2 50/54, n3 38/51, n4+ 49/58, n1 22/30.
- **Labels.** Re-applying the SESOI rules gives the same label for every contrast.

### 5. Decomposition (confirm half, pooling A): CONFIRMED. These are undeclared, exploratory contrasts
They are not in summary.json or in the document. If they are added, they must be labelled exploratory.

| contrast | face mid | face min | face max |
|---|---|---|---|
| est alone at X100: est/lin/G0/X100 − kn/lin/G0/X100 | −0.0128 [−0.0186, −0.0088] | −0.0033 [−0.0100, −0.0010] | −0.0224 [−0.0285, −0.0148] |
| X1 alone with kn: kn/lin/G0/X1 − kn/lin/G0/X100 | −0.0040 [−0.0064, −0.0012] | −0.0040 | −0.0040 (a uniform shift) |
| X1 given est: est/lin/G0/X1 − est/lin/G0/X100 | −0.0043 [−0.0061, −0.0017] | −0.0043 | −0.0043 |
| est given X1: est/lin/G0/X1 − kn/lin/G0/X1 | −0.0132 [−0.0175, −0.0096] | | |
| log vs lin, est at X1 | +0.0039 [+0.0011, +0.0053] | | |
| Gmeas − G0, est/lin/X1 | +0.000038 [−0.000003, +0.000105] | | |
| G0.2 − G0, est/lin/X1 | +0.000038 [−0.000003, +0.000105] | | |

- **Additivity.** est at X100 plus X1 with kn gives −0.0169, against the joint −0.0172; the interaction is −0.0003.
- **X1 under the other probability sources.** two_tier X1 − X100 = −0.0006 [−0.0019, −0.0001]; oracle X1 − X100 = −0.0085 [−0.0124, −0.0044].

### 6. Health: CONFIRMED, plus one row the count check misses
- **Run health.**
  - All 360 rows are ok, each with 90 arms, and all 32,400 solves are optimal.
  - Biomass ratio lies in [0.9499999999999998, 0.9500000000000001]; no arm is off the floor by more than 1e-6.
  - There are no face-LP failures or retries.
  - All 3,240 determinism pairs are identical.
  - No gapfilled reaction is gene-backed, and no row has α > 0 with a null clade.
  - Internal consistency: taf_face_mid equals (min + max) / 2 / D exactly; exchange_term / X equals flux_exchange to 2e-15; gene-less flux equals the internal-term split to 3e-13; D equals the two_tier/lin/G0/X100 gene-backed flux exactly.
- **`prior_match.ok` false (2 rows):**
  - `1674805.4__V3-V4`: tune, n4+, grp00037. n_new_rxn 50 vs 56, growth 13.947 vs 13.541, n_gene_backed 1158 vs 1100, n_gapfilled 50 vs 56; D 965.08 vs 956.24.
  - `1798454.3__V4`: confirm, n4+, grp00115. n_new_rxn 58 vs 59, growth 4.667 vs 4.796, n_gapfilled 58 vs 59; D 836.29 vs 827.59.
- **A third row the count check misses.** `1017270.21__V4-V5` (confirm, n4+, grp00089) passes every count, but rebuilds a different LP. Its D is 667.78 against 656.27. Its shipped objective is 16758.36 against 16775.29, and its released objective also differs from the prior one. These three rows are exactly the 3/360 failures of the "D bitwise", "objective bitwise" and "released = prior bitwise" checks. A fourth row, `2053306.99__V4` (confirm), has bitwise-equal D and shipped objective, but its shipped face midpoint differs by 3.6e-5.
- **Excluding them changes nothing.** Excluding 2, 3 or 4 of these rows leaves the selection at est/lin/G0/X1, with tune pooled A 0.11967. Selected − shipped (A) becomes −0.0171 [−0.0221, −0.0118] (2 excluded) or −0.0171 [−0.0220, −0.0119] (3 or 4). Selected − two_tier is −0.0222 to −0.0223, selected − shuffled −0.0282 to −0.0283, and shipped − two_tier −0.0057.

## Leakage and method

- **Execution: no leakage found.**
  - `objective_tuning_results.py` selects from `r["half"] == "tune"` only, and that field equals the recomputed hash on every row.
  - g_meas is recomputed from the prior run's tune-half rows as 1021 / 12222 = 0.0835379, equal to meta.json. The confirm half would give 0.0851, and G does not matter anyway.
  - The smoke run used stub weights (`stage_o_weights.STUB.json`, 6 rows, 2 of them confirm).
  - The E3 weights file was written at 08:50:15. The run started at 08:50:28.
  - The SHA-256 hashes of the analysis script, harness, design and weights equal those recorded in summary.json and meta.json. The analysis script and template were last modified on 2026-11-02 at 23:05, before the run. **They are untracked in git**, so "fixed before the full run" rests on file mtimes and the recorded hashes, not on a commit.
- **Design-level exposure of the confirm half.** These are caveats, not defects in Stage O:
  - All 193 confirm amplicons were in the weighted-objective test, which reported results on all 360 before this design was declared. Among them was the exploratory exchange-weight-1 arm (real_ex1, identical to kn/lin/G0/X1 on 356/360 rows). On the confirm half that arm gave −0.0090 [−0.0111, −0.0050] against two_tier_ex1.
  - The X1 level, which supplies about 25% of the headline, was therefore put on the menu with the confirm half's data already seen. The confirm half is held out from *selection*, not from the *design*.
  - Stage E chose E3's hyperparameters on all five benchmark folds, which contain the confirm amplicons' cells (0.6% of cells). All five outer folds picked the same variant, so the influence is negligible.
- **Bootstrap independence.** The per-stratum bootstrap resamples a group independently in each stratum. Ten confirm groups have rows in several strata, including grp00000 in all four. A joint resample (each group drawn once across strata) gives nearly identical intervals:
  - A: [−0.0221, −0.0116], against [−0.0221, −0.0118];
  - B: [−0.0123, +0.0042], against [−0.0119, +0.0032].
- **Single permutation.** est_shuffled uses one permutation per amplicon, so the intervals carry no permutation variance. I did not verify the permutation itself, because per-reaction values are not stored.

## The TAF mechanism: is it steering?

I split each confirm row's returned vertex into flux classes, in units of D, and pooled with A:
- gene-backed truth-absent (`flux_absent`);
- gene-backed truth-supported (`flux_gene_backed − flux_absent`);
- gapfilled;
- gene-less draft (`flux_rxn − gene_backed − gapfilled`);
- exchange.

| contrast | absent | supported (gene-backed) | gapfilled | gene-less draft | exchange | TAF face width |
|---|---|---|---|---|---|---|
| selected − shipped | −0.0146 [−0.0201, −0.0097] | +0.0207 [+0.0038, +0.0294] | −0.0004 [−0.0011, −0.0001] | −0.0019 [−0.0034, −0.0004] | +0.0034 [+0.0016, +0.0044] | −0.0191 [−0.0234, −0.0100] |
| est at X100 − shipped | −0.0103 [−0.0165, −0.0067] | +0.0341 [+0.0246, +0.0443] | 0.0000 | −0.0000 | −0.0000 | −0.0191 |
| X1 given est (est X1 − est X100) | −0.0043 [−0.0061, −0.0017] | −0.0134 [−0.0307, −0.0078] | −0.0004 | −0.0019 | +0.0034 | +0.0000 |
| X1 with kn (kn X1 − kn X100) | −0.0041 | −0.0041 [−0.0156, +0.0002] | −0.0003 | −0.0019 | +0.0022 | −0.0001 |

I also split the absent-flux change exactly into a volume term and a selectivity term. With s = absent / gene-backed for an arm, Δabsent = s_b·Δgene-backed + gene-backed_a·Δs.

| contrast | volume | selectivity |
|---|---|---|
| selected − shipped | +0.0002 [−0.0016, +0.0017] | −0.0148 [−0.0197, −0.0100] |
| est at X100 − shipped | +0.0023 | −0.0125 |
| X1 given est | −0.0015 [−0.0025, −0.0008] | −0.0028 [−0.0039, −0.0008] |
| X1 with kn | −0.0013 | −0.0028 |

**Verdict:**
- **For the headline contrast, the TAF drop is steering, not displacement.** Truth-supported gene-backed flux rises. Gene-less and gapfilled flux fall rather than rise. The vertex drop is entirely selectivity: the volume term is +0.0002. Flux precision rises by +0.0149.
- **The two components work differently.**
  - *The estimator* is pure re-routing among gene-backed reactions. At X100, exchange, gene-less and gapfilled flux are unchanged, truth-absent flux falls and truth-supported flux rises. It also narrows the TAF face by 0.019 D: most of its midpoint gain sits at the face maximum (−0.0224 against −0.0033 at the face minimum).
  - *X1* shifts the whole face by a uniform amount, with no change in width. It raises exchange flux and lowers **every** internal class: truth-supported flux falls three times as much as truth-absent flux. This is uptake substitution, importing medium components rather than synthesising them. About a third of X1's TAF drop (−0.0015 of −0.0043) is volume, less gene-backed flux against the fixed D. The other two thirds is a precision gain (+0.0029 [+0.0008, +0.0041]).
- **Row level** (163 confirm rows with n ≥ 2, selected − shipped):
  - absent flux falls on 112 rows;
  - on 41 of those, truth-supported flux falls as well;
  - on 36, gene-less flux rises;
  - the largest gene-backed drop is −0.20 D (551275.6__V1-V2), almost all of it attributable to X1.
- **Face minimum.** At the best vertex of each face the gain is only −0.0076, below the SESOI. The declared midpoint rule and the as-returned vertex (taf_point_D −0.0146 [−0.0201, −0.0097]) both pass.

## est scale (E3 against k/n, all 360 rows)
- **Few reactions reach the top.** Under E3, 36,460 of 312,243 gene-backed reactions (11.7%; per-row median 7.8%) reach exactly p = 1. 48,052 (15.4%) have p ≥ 0.995 and are certainly at the linear cost floor. About 75,556 (24%) have p ≥ 0.985. 11 rows have no reaction at p = 1. Under k/n, 241,237 (77.3%) are at p = 1, and 100% at n1.
- **Mean linear cost of gene-backed reactions.**

| stratum | E3 | k/n |
|---|---|---|
| n1 | 0.181 | 0.010 (18× lower) |
| n2 | 0.205 | 0.111 |
| n3 | 0.204 | 0.135 |
| n4+ | 0.234 | 0.161 |
| all rows | 0.209 | 0.115 |

  So E3 roughly doubles gene-backed costs overall and narrows their gap to gene-less cost (1).
- **Log cost.** The log cost reaches its floor at p ≥ 0.955. That holds about 39% of gene-backed reactions under E3, against 78% under k/n. The log shape therefore flattens the top, where most of the flux is; this is consistent with log never winning.
- **Within-row discrimination.** Mean AUROC of truth-supported against truth-absent reactions is 0.872 for E3 and 0.674 for k/n.
- **Scale or ranking?** The est_shuffled control keeps E3's scale and destroys its ranking, and it is *worse* than the shipped arm:
  - est_shuffled − shipped: +0.0110 [+0.0031, +0.0222] (face minimum +0.0203);
  - est_shuffled − two_tier at X1: +0.0060 [−0.0011, +0.0174].
  So the ranking, not the scale, carries the gain. It also means that selected − shuffled (−0.0282) partly reflects a harmful comparator. Measured from a neutral base, the gain is the −0.0172 against shipped.

## Robustness
- **Leave one group out** (confirm half, pooling A):
  - selected − shipped: [−0.0174, −0.0162]. The largest influences are grp00640 (+0.0010) and grp00000 (+0.0008). No single deletion moves the mean above −0.01.
  - selected − two_tier: [−0.0225, −0.0199]; selected − shuffled: [−0.0288, −0.0259]; est at X100 − shipped: [−0.0132, −0.0118].
  - Pooling B for selected − shipped: [−0.0092, −0.0043].
- **grp00000** (selected − shipped, pooling A):
  - Without it: face mid −0.0164 [−0.0236, −0.0100], face min −0.0094 [−0.0157, −0.0039], face max −0.0233 [−0.0331, −0.0147].
  - Pooling B without it: −0.0043 [−0.0127, +0.0067].
  - Within it alone: −0.0175.
  - The effect is homogeneous in and out of the dominant group. Bootstrap draws hold 0 copies of it 34–37% of the time per stratum.
- **Seeds** (pooling A CI for selected − shipped, face mid):

| seed | CI |
|---|---|
| 1729 | [−0.0221, −0.0118] |
| 424242 (wobj_results' ALT_SEED) | [−0.0217, −0.0120] |
| 1 | [−0.0221, −0.0121] |
| 2026 | [−0.0224, −0.0120] |
| 99991 | [−0.0220, −0.0119] |
| 1729 with B = 20,000 | [−0.0221, −0.0119] |

  The pooling B CI stays across 0 under every seed, its upper end between +0.0027 and +0.0040. The n1 CI stays across 0, with its upper end between +0.0057 and +0.0071.
- **n = 1 did not replicate.** The tune n1 stratum gave −0.0124 [−0.0273, −0.0030]; the confirm n1 stratum gives −0.0048 [−0.0113, +0.0063]. Pooling B went from −0.0130 [−0.0258, −0.0047] on the tune half to −0.0065 [−0.0119, +0.0032] on the confirm half.

## Degeneracy (FVA second pass): reproduced
- The four FVA-pass arms on all 193 confirm rows equal the main run bitwise on objective, face min/mid/max, absent, gene-backed and exchange flux (0 mismatches).
- Pooled-A contrasts of the count of gene-backed reactions with FVA range > 1e-3 reproduce the document:

| contrast | optimal face | 1e-7 face |
|---|---|---|
| selected − shipped | −21.89 [−24.65, −14.81] | −22.43 [−25.08, −15.93] |
| selected − two_tier | −27.97 [−30.73, −20.24] | −28.49 [−31.10, −21.54] |
| selected − shuffled | +11.46 [+9.95, +13.93] | +11.50 [+9.64, +13.48] |

- Pooled-A levels on the optimal face: selected 13.6, shipped 35.4, two_tier at X1 41.5, est_shuffled at X1 2.1. The shuffled arm is the *least* degenerate and the worst on TAF, so tie-breaking alone does not explain the gain.

## Proposed corrections to `data/objective_tuning/OBJECTIVE_TUNING_RESULTS.filled.md`
No number in the document is wrong. Each item below gives the exact old text and the new text.

1. **Checks table, line 27** (the count check misses the third LP difference).
   - Old: `| rebuilt model = prior row (n_new_rxn, growth 1e-6, n_gene_backed, n_gapfilled) | 358 / 360 |`
   - New: `| rebuilt model = prior row (n_new_rxn, growth 1e-6, n_gene_backed, n_gapfilled) | 358 / 360 (fail: 1674805.4__V3-V4 tune n4+; 1798454.3__V4 confirm n4+). A third row, 1017270.21__V4-V5 (confirm n4+), passes these counts but rebuilds a different LP (D 667.78 vs 656.27; shipped and released objectives differ), so 357 / 360 rows rebuild the prior LP. Excluding all three changes neither the selection nor any conclusion (selected − shipped, A: −0.0171 [−0.0220, −0.0119]). |`

2. **Section 1, Sample** (same point).
   - Old: `The gapfilled model of each selection's union EC set was rebuilt as before:`
   - New: `The gapfilled model of each selection's union EC set was rebuilt as before (the prior LP exactly on 357 / 360 rows; see the checks):`

3. **Confirm table header** (taf_gb is not on the D scale).
   - Old: `| contrast | pooling | face mid | face min | face max | as returned (taf_gb) | label |`
   - New: `| contrast | pooling | face mid | face min | face max | as returned (taf_gb = truth-absent / the arm's own gene-backed flux = 1 − flux precision; not / D) | label |`

4. **After the guard table** (insert after `| selected-shipped | B | +0.9219 [+0.4818, +1.5527] | +8.8665 [+1.6885, +13.3366] | +0.0106 [+0.0003, +0.0159] | -0.2961 [-1.0766, -0.0431] |`):
   - New: `The flux-precision column is the taf_gb column with its sign flipped (taf_gb = 1 − precision), not separate evidence. At n = 1, where prFBA usually runs (86.8% of answered amplicons), selected − shipped is −0.0048 [−0.0113, +0.0063] (30 rows, 20 groups): not detected. On the tune half it was −0.0124 [−0.0273, −0.0030], so it did not replicate, and pooling B is inconclusive (−0.0065 [−0.0119, +0.0032]). The result that meets the SESOI holds under pooling A (n ≥ 2) only. At the best vertex of each face (face min) the gain is −0.0076, below the SESOI. The solver's returned vertex gains −0.0146 [−0.0201, −0.0097] (taf_point_D).`

5. **New paragraph after item 4** (mechanism; exploratory, undeclared contrasts).
   - New: `**Where the gain comes from (exploratory; contrasts not declared).** At X100, est − kn is −0.0128 [−0.0186, −0.0088] at the face midpoint. It re-routes flux among gene-backed reactions: truth-supported flux +0.034 D, with exchange, gene-less and gapfilled flux unchanged. It also narrows the TAF face by 0.019 D. With est, X1 − X100 adds −0.0043 [−0.0061, −0.0017] as a uniform shift of the face, by a different route: exchange flux rises (+0.0034 D) and every internal class falls (truth-supported −0.013 D, gene-less −0.002 D, gapfilled −0.0004 D). About a third of that is less gene-backed flux against the fixed D; the rest is a precision gain of +0.0029. The two parts add up (−0.0169 against −0.0172). Over the selected configuration as a whole, the drop at the returned vertex is all selectivity (volume term +0.0002 [−0.0016, +0.0017]), and gene-less and gapfilled flux fall. Neither run records per-exchange fluxes, so which nutrients X1 imports more of is not known.`

6. **Section 3, eligibility** (the guard never bound).
   - Old: `36 of 36 configurations are eligible. By exchange weight: {"X100": 12, "X10": 12, "X1": 12}.`
   - New: `36 of 36 configurations are eligible. By exchange weight: {"X100": 12, "X10": 12, "X1": 12}. The guard did not bind: no candidate's median ratio departs from 1 by more than 0.45% (exchange flux) or 3.0% (gene-backed flux). Per row it is wider. For the selected arm, the gene-backed flux ratio to the shipped arm spans 0.65–1.45 on the tune half (29 / 167 rows outside ±10%) and 0.80–1.27 on the confirm half (12 / 193). The exchange flux ratio spans 1.000–1.12 (1 / 167) and 1.000–1.18 (4 / 193), and it never falls below 1.`

7. **Section 4, headroom line** (comparability).
   - Old: `Oracle headroom captured on the confirm half (pooling A): selected (own cell) 0.685; shipped 0.231.`
   - New: `Oracle headroom captured on the confirm half (pooling A): selected (own cell) 0.685; shipped 0.231; selected in the shipped cell's frame 0.933 (tune: 0.908).`

8. **Degeneracy bullets** (units and pooling).
   - Old: `**Degeneracy** (FVA on the optimal face, finalists only, second pass):`
   - New: `**Degeneracy** (FVA on the optimal face, finalists only, second pass; the change in the number of gene-backed reactions with FVA range > 1e-3, pooling A; arm levels: selected 13.6, shipped 35.4, two_tier at X1 41.5, est_shuffled at X1 2.1. The shuffled arm is the least degenerate and the worst on TAF, so breaking ties does not by itself lower TAF):`

9. **Section 6, Limitations** (append after `- Pooling B rests on the n1 stratum: about 30 rows per half.`).
   - New: `- The confirm half is held out from selection, not from the design. Its 193 amplicons were in the weighted-objective test, whose exploratory exchange-weight-1 arm (identical to kn/lin/G0/X1 on 356 / 360 rows; −0.0090 [−0.0111, −0.0050] against two_tier_ex1 on this half) had been seen when the X factor was declared.`
   - New: `- One holdout group (grp00000, fold 0, 96 source genomes) holds 96 of the 193 confirm rows (43–60% of each stratum). Without it, selected − shipped is −0.0164 [−0.0236, −0.0100] (A); within it, −0.0175. No single-group deletion moves the pooled-A mean outside [−0.0174, −0.0162].`
   - New: `- The guard is on total exchange flux only. Uptake composition under X1 was not recorded, and the selected arm imports more on 57% of rows (up to +18%).`
   - New: `- scripts/objective_tuning.py, scripts/objective_tuning_results.py and OBJECTIVE_TUNING_RESULTS.md are untracked in git. "Fixed before the full run" rests on file times (2026-11-02 23:05 / 23:15; the run started 2026-11-03 08:50) and on the hashes recorded in meta.json and summary.json.`

## Files (all in this directory)
- `extract.py`: streams the raw rows into `rows.pkl` and `prior.pkl`.
- `vlib.py`: my own helpers (half hash, weights, cluster bootstrap, joint bootstrap).
- `c1_split_health.py` → `out_c1_c6.json`: claims 1 and 6, and row consistency.
- `c2_c3_selection.py` → `out_c2_c3.json`: claims 2 and 3, per-row guard ratios, G reading.
- `c4_c5_confirm.py` → `out_c4_c5.json`: claims 4 and 5, alternative seed and joint bootstrap, a comparison of 972 numbers with summary.json.
- `c7_mechanism.py` → `out_c7.json`: flux-class decomposition.
- `c8_robust.py` → `out_c8.json`: precision decomposition, leave-one-group-out, grp00000, seeds, exclusions.
- `c9_est_scale.py` → `out_c9.json`: E3 and k/n placement on the cost scale.
- `c10_prior_fva.py` → `out_c10.json`: prior-study cross-checks, g_meas, FVA pass.
- `c11_tables.py`: guard table and checks table against the document.
- `c12_selection_stability.py` → `out_c12.json`: tune-half selection bootstrap, tune against confirm by stratum.
- `c13_volume_split.py` → `out_c13.json`: volume and selectivity split.
