# Weighted-objective test: results

Run 2026-10-29, 20:41 to 22:25 (population stage, then 360 amplicons on 8 workers). Design: `WEIGHTED_OBJECTIVE_TEST.md` (declared post hoc; decision rules amended before launch, section 12). Outputs: `data/weighted_objective/{population.json, per_amplicon.jsonl, summary.json}`; log `logs/wobj_run.log`. Every number below was recomputed from `per_amplicon.jsonl` by `scripts/wobj_results.py`. That script writes `data/weighted_objective/recomputed_contrasts.json` and checks it against `summary.json`. It compares 645 values (the mean and both CI bounds of 215 contrast-by-measure-by-stratum cells, n1 to pooled, for 7 contrasts and 7 measures); all agree, the largest absolute difference being 5.2e-18. It also re-applies the section-10 decision rules to its own numbers and reaches the same C3, H1, H1s and H2 labels as `summary.json` `verdict`. The script also draws the two figures.

Independent re-derivation (2026-10-30, separate code; scripts `rederive.py`, `rederive2.py` and `altseed.py` in the session scratch directory `wobj_verify/`): it recomputed 705 contrast values from `per_amplicon.jsonl` (largest difference from `summary.json` 1.8e-15), the arm levels, H2 counts, concentration, composition, checks, Wilcoxon/BH and the alternative-seed shifts. It confirmed the numbers here and corrected these statements: the union-only role count (section 3; `wobj_results.py` had reported 0), the all-arm vertex re-solve maximum (section 3), the SESOI sentence and the oracle-remainder explanation (section 1), and the fix branch and community-path status (sections 1 and 11).

Intervals are 95% percentile cluster bootstraps over holdout groups (2,000 draws), as declared. "Pooled" means the n >= 2 strata weighted by eligible source genomes: n2 0.6975, n3 0.1869, n4+ 0.1156 (2,474 / 663 / 410). Flux measures are fractions of `D`, the two-tier arm's gene-backed flux, so 0.01 is one percentage point of it. That is the declared SESOI.

## 1. The answer

**Does the weighting steer flux toward the true organism's functions?** A little, and not robustly. The question is asked only where the weighting is active, on amplicons with two or more selected genomes. There, the real-p objective's optimal face carries 0.6 percentage points less truth-absent flux than the two-tier objective at its midpoint: −0.0061 [−0.0097, −0.0025]. That is a 4.7% relative reduction of two-tier's 0.130. The reduction is −0.0096 [−0.0146, −0.0054] at the face maximum and −0.0025 [−0.0053, +0.0011] at the face minimum. So the sign is not established at the low end of the face. The midpoint mean and its whole CI lie short of the 0.01 SESOI; among the face measures only the face-maximum CI reaches past it (to −0.0146), as does the secondary as-returned TAF_gb CI (to −0.0117). The mechanical verdict is **"detectable but below SESOI"**.

The effect is concentrated. The median amplicon moves essentially nothing in n2 and n3. On 81 of 300 multi-genome amplicons real p lowers the midpoint by at least one point; on 18 it raises it by at least one point. The 30 largest reductions carry 73% of the net.

**Is the flux effect specific to p?** Yes, by the declared rule. Real beats the shuffled control, which places the same coefficients at random, by −0.0120 [−0.0179, −0.0069]. About half of that gap is shuffling doing harm, not p doing good. Shuffled carries more truth-absent flux than two-tier: +0.0060 [+0.0028, +0.0106].

**Does the weighting reduce degeneracy?** Yes, but not specifically (next paragraph). On the optimal face, real leaves 21% fewer gene-backed reactions free to vary than two-tier: 32.7 against 41.5, a difference of −8.9 [−10.3, −7.5]. A 1e-7 face confirms it.

**Is the degeneracy effect specific to p?** No. Shuffled coefficients shrink the face more: 27.9 variable reactions, a 33% reduction, and real − shuffled is +4.8 [+3.4, +6.0]. Plain parsimony (uniform) leaves 41.8, the same as two-tier. Any gradation of the gene-backed coefficients breaks ties. The gene-support tier alone does not. The information in p breaks fewer ties than a random gradation does.

**How much could any weighting achieve?** An oracle that prices only the truth-absent gene-backed reactions lowers the midpoint by −0.0249 [−0.0291, −0.0210]. So the objective can express a reweighting (control C3 passes). Real p captures 24% of that headroom. Even the oracle removes only 19% of two-tier's truth-absent flux. The oracle already prices every truth-absent gene-backed reaction at coefficient 1, as it prices a gene-less reaction, so the flux it leaves there is flux the constraints (95% of maximum biomass at the minimum exchange flux) gave it no cheaper way to reroute. That remainder was not decomposed further.

**Why so small?** The exchange term is 98.9% of the objective (median). It is identical across the six exchange-weight-100 arms on 357 of 360 amplicons, so the objective acts lexicographically. It first minimizes exchange flux; the internal weights then only choose among those solutions. The weighted gene-backed p < 1 term is a median 0.11% of the objective value.

**Where the paper says prFBA usually runs** (one selected genome, 86.8% of answered amplicons on this benchmark), real and two-tier are the same objective. Their returned solutions were bitwise identical on 60 of 60 single-genome amplicons.

**The released single-model path of `MSProbability.prFBA` does not run** (section 11). The test ran a documented faithful copy and re-solved the real arm through the released source on every amplicon; the two agreed exactly on 360 of 360. Fix: ModelSEEDpy branch `fix-prfba-single-model`, commit b6857a0 (parent c2e31a0), in the local worktree `/home/freiburger/Documents/ModelSEEDpy-prfba`. It is not pushed.

## 2. What was run

- **Population** (`population.json` `meta`): the controlled shared-reference benchmark (leave-genus-out folds, shipped exact-tie rule, full reference density, EC space). It has 54,243 candidate cells with truth >= 200 ECs. One amplicon per source genome and per md5 was drawn, excluding md5s scored for more than one source genome. The sample is n2 100 (66 holdout groups), n3 100 (66), n4+ 100 (54; 4 to 13 references) and n1 60 (48). Templates: grampos 193, gramneg 126, archaea 41.
- **Model per amplicon**: the gapfilled ModelSEED model of the union EC set (ModelSEEDpy c2e31a0, templates 21a8084, gapfiller reversals applied). It runs complete medium at uptake 100, a carbon cap of 100 and biomass >= 0.95 × max. The constraint set is the same in every arm; only the objective coefficients differ.
- **p** follows the production rule on the EC route: p(role) is the share of selected genomes whose EC set covers the role, and p(reaction) is the maximum over its genes. Gene-less and gapfilled reactions get 0, so their coefficient is 1.
- **Arms**:
  - `two_tier` (gene-backed p = 1): what production gives at n = 1.
  - `real`.
  - `uniform`: coefficient 0.01 everywhere.
  - `shuffled`: real p permuted among gene-backed reactions.
  - `oracle`: p = 1 if the reaction is truth-supported, else 0.
  - `real_exp2`: prob_exp 2.
  - `real_ex1` and `two_tier_ex1`: exchange weight 1.
- **Truth support** comes from the truth genome's own draft, built by the same route. A gene-backed reaction it does not support is "truth-absent".
- **Measures**:
  - Optimal face: objective <= opt + 1e-6, one absolute slack for every arm.
  - `taf_face_min` / `max` / `mid`: total truth-absent gene-backed flux on that face, divided by `D`.
  - Degeneracy: the number of gene-backed reactions whose range on that face exceeds 1e-3.
- **Provenance**:
  - Harness sha256 `d44e7fa3…`, identical to the committed file and to the second smoke's.
  - prFBA repository at f8205be.
  - Libraries: cobra 0.29.1, gurobipy 13.0.1, optlang 1.9.1, cobrakbase 5e4eedb, mscommunity 176b6e1 (`summary.json` `meta.libraries`).

## 3. Checks (the gate and the run's health)

| check | result |
|---|---|
| amplicons ok / gated / unscoreable | 360 / 0 / 0 |
| arm solves optimal | 2,880 / 2,880 |
| C2 `released_check` (copy vs released source, real arm) | 360 / 360 pass, max relative difference 0.0, one released-source sha256 (`6eeecfb7…`) |
| C1 n1 returned solution, real vs two_tier | bitwise identical 60 / 60 |
| n1 face values within 1e-9 | 59 / 60 (largest difference 7.2e-5, 2053306.99__V4: LP warm-start tolerance, three orders below any decision threshold) |
| face-LP failures | 0, after one cold retry of 157 LPs (9 on 1e-6 faces, 148 on 1e-7 faces) |
| feasible-scan LP failures | 0 |
| growth | every arm solve at 0.95 × max biomass (ratio 0.95 ± 2e-16); max biomass identical across arms |
| gapfill sanity gate | 360 / 360 pass; median 72.5 gapfilled reactions; reversals applied in 57 models (64 reactions) |
| `D` face-robust | largest `D_face_rel_width` 2.0e-5 |
| opt7 (1e-7) count differs from opt by > 5% | two_tier 34, real 38, uniform 41, shuffled 49 of 360 rows; the H2 label is the same at 1e-7 |
| vertex re-solves (Method 0/1/2) | all reach the same objective (max relative difference 1.3e-13 for real, 2.0e-12 over all eight arms, the largest on the exchange-weight-1 pair); TAF as returned depends on the vertex on 203 / 360 rows for real and 269 / 360 for two_tier (largest spread 0.196 and 0.197); the oracle's on 1 / 360 |
| genes outside template roles; union-only roles | 0; 15 roles on 6 of 360 rows (all multi: 2 n2, 2 n3, 2 n4+; 1 to 5 per row), each given p = 0, but no gene-backed reaction ends at p = 0 (section 10) |
| time, memory | wall mean 124.6 s per amplicon (max 267.3 s); largest worker 1.804 GiB (the plan estimated 1.75), parent 0.743 GiB |

C1 and C2 hold, so `verdict()` reads the contrasts.

## 4. H1: truth-absent flux on the optimal face

Mean paired difference, fraction of `D`, 95% group-cluster CI.

| contrast | measure | n2 | n3 | n4+ | pooled n >= 2 |
|---|---|---|---|---|---|
| real − two_tier | face mid | −0.0054 [−0.0102, −0.0004] | −0.0061 [−0.0092, −0.0022] | −0.0100 [−0.0137, −0.0052] | **−0.0061 [−0.0097, −0.0025]** |
| real − two_tier | face min | −0.0023 [−0.0062, +0.0029] | −0.0024 [−0.0050, +0.0012] | −0.0039 [−0.0072, +0.0005] | −0.0025 [−0.0053, +0.0011] |
| real − two_tier | face max | −0.0085 [−0.0152, −0.0030] | −0.0099 [−0.0147, −0.0043] | −0.0160 [−0.0217, −0.0091] | −0.0096 [−0.0146, −0.0054] |
| real − two_tier | as returned (TAF_gb) | −0.0074 [−0.0125, −0.0039] | −0.0085 [−0.0123, −0.0038] | −0.0118 [−0.0183, −0.0072] | −0.0081 [−0.0117, −0.0054] |
| real − shuffled | face mid | −0.0101 [−0.0177, −0.0038] | −0.0091 [−0.0143, −0.0027] | −0.0282 [−0.0525, −0.0166] | **−0.0120 [−0.0179, −0.0069]** |
| shuffled − two_tier | face mid | +0.0048 [+0.0011, +0.0102] | +0.0030 [−0.0013, +0.0063] | +0.0183 [+0.0067, +0.0435] | +0.0060 [+0.0028, +0.0106] |
| real − uniform | face mid | −0.0039 [−0.0095, +0.0023] | −0.0049 [−0.0084, −0.0000] | −0.0108 [−0.0143, −0.0061] | −0.0049 [−0.0088, −0.0002] |
| oracle − two_tier | face mid | −0.0228 [−0.0285, −0.0182] | −0.0297 [−0.0362, −0.0213] | −0.0297 [−0.0363, −0.0189] | **−0.0249 [−0.0291, −0.0210]** |

- **Arm levels** (pooled face mid): two_tier 0.1302, uniform 0.1290, real 0.1241, shuffled 0.1361, oracle 0.1053. Real's optimal face is also narrower than two-tier's: mean width 0.0133 against 0.0195 in n2, 0.0198 against 0.0272 in n3, and 0.0207 against 0.0328 in n4+.
- **Per-stratum tests** (Wilcoxon signed-rank, Pratt; BH over the declared family of 6):
  - real − two_tier face mid: q 2.3e-7 (n2), 7.8e-5 (n3) and 3.1e-7 (n4+); rank-biserial −0.61, −0.46 and −0.60.
  - real − shuffled face mid: q at most 2e-6 in every stratum.
- **Concentration.** Of the 300 multi-genome amplicons (unweighted), 81 have real − two_tier <= −0.01 at the face midpoint and 18 have >= +0.01. The median difference is −6e-7 in n2, −8e-8 in n3 and −0.0035 in n4+. The quantiles are 5%: −0.052, 25%: −0.011, 50%: −0.0003, 75%: 0.0000 and 95%: +0.011. The 30 largest reductions sum to 73% of the net reduction.
- **Why.** The weighting can move only truth-absent flux that sits on p < 1 reactions. In the two-tier arm that addressable share averages 0.385 (n2), 0.474 (n3) and 0.616 (n4+). The median model has 913.5 gene-backed reactions, of which 191.5 have p < 1 (n2 119.5, n3 202, n4+ 236). Among gene-backed reactions the p < 1 ones are truth-supported 0.467 of the time and the p = 1 ones 0.867 (static precision, multi mean). This matches the calibration finding that p ranks well.
- **Verdict** (`summary.json` `verdict.H1`): "detectable but below SESOI: the midpoint CI excludes 0, |mean| < SESOI, and the CIs do not lie inside +-SESOI". The midpoint mean misses the SESOI, the face-minimum CI includes 0 and the face-maximum CI extends past −0.01. `H1s_specific_to_p` is true.

## 5. The oracle ceiling

The positive control passes (`positive_control_oracle_expressible` true). The oracle's face-midpoint CI lies below −0.01, and its face-minimum and face-maximum CIs lie below 0: min −0.0138 [−0.0179, −0.0109], max −0.0360 [−0.0416, −0.0301]. A reweighting of gene-backed reactions can therefore move truth-absent flux under this objective. The H1 shortfall is a fact about p, not a structural null.

Real p captures 0.243 of the oracle's pooled face-midpoint effect: 0.236 in n2, 0.207 in n3 and 0.336 in n4+ (0.327 on the as-returned TAF). The oracle itself removes 0.0249 of two-tier's 0.1302, i.e. 19%. Even at n = 1, where p carries no information, the oracle lowers the midpoint by −0.0170 [−0.0250, −0.0091]. A perfect annotation signal would help there too; p = 1 cannot supply it.

## 6. H2: degeneracy

The measure is the number of gene-backed reactions with range > 1e-3 on the optimal face (slack 1e-6). The table gives each arm's mean and the paired difference.

| | n2 | n3 | n4+ | pooled n >= 2 |
|---|---|---|---|---|
| two_tier mean | 39.9 | 43.3 | 48.7 | 41.5 |
| real mean | 32.5 | 32.9 | 33.6 | 32.7 |
| shuffled mean | 29.0 | 25.4 | 25.1 | 27.9 |
| uniform mean | 40.3 | 42.8 | 48.8 | 41.8 |
| real − two_tier | −7.4 [−9.3, −5.6] (19%) | −10.4 [−12.5, −8.0] (24%) | −15.2 [−18.0, −12.7] (31%) | **−8.9 [−10.3, −7.5] (21% [18%, 25%])** |
| real − two_tier, 1e-7 face | −7.6 [−9.7, −5.8] | −10.5 [−12.8, −8.1] | −15.2 [−18.1, −12.8] | −9.0 [−10.5, −7.7] (22%) |
| shuffled − two_tier | −10.9 [−13.4, −8.4] | −17.8 [−21.3, −13.3] | −23.6 [−25.9, −20.0] | −13.7 [−15.5, −11.7] (33% [28%, 37%]) |
| real − shuffled | +3.5 [+1.9, +5.1] | +7.4 [+4.3, +9.9] | +8.5 [+5.6, +10.0] | **+4.8 [+3.4, +6.0]** |
| real − uniform | −7.8 [−10.0, −6.1] | −10.0 [−12.2, −7.7] | −15.2 [−18.2, −12.8] | −9.1 [−10.6, −7.8] |

- Real has fewer variable reactions than two-tier on 247 of 300 multi-genome amplicons. The opt-face sum of ranges falls too: −20.9 [−32.3, −2.4] against a two-tier 161.8.
- **Verdict** (`verdict.H2`, the same at the 1e-7 face): "reduces degeneracy: >= 10% fewer variable gene-backed reactions on the optimal face; NOT specific to p: shuffled costs do as well or better".
- **Noise scale.** On n1, real and two-tier solve the same LP, yet their opt-face counts agree on only 38 of 60 amplicons. The largest difference is 4 reactions and the mean +0.20 [0.00, +0.42]. The difference is exactly 0 at the 1e-7 face. The pooled −8.9 is far outside this.
- **Tolerance faces carry no arm signal.** On the face at 0.1% of the two-tier internal term, the count is 238.4 for real, two_tier and shuffled alike (multi mean). On the face at 0.1% of the whole objective, the real arm's sum of ranges is a mean 52,832 against 147 on the optimal face.

## 7. Objective composition (H3)

The real arm's objective on multi-genome amplicons:

- **Overall split.** The internal weighted term is a median 1.05% of the objective (range 0.53% to 3.44%). The exchange term is the rest (at least 96.6% on every row). The two-tier arm's internal share is a median 0.95%, and the uniform arm's 0.045%. `internal + exchange = objective` holds to a relative 2.0e-12 on every arm and row.
- **Inside the internal term** (median of per-row shares): gene-less draft reactions at coefficient 1 take 80.4%. Gene-backed p < 1 reactions take 11.1%, gene-backed p = 1 reactions 3.3% and gapfilled reactions 2.6%. The means are 77.4%, 14.7%, 3.2% and 4.7%. The gene-backed p < 1 term, the only part the weighting changes, is a median 0.11% of the objective value (max 2.5%).
- **Lexicographic order.** The exchange term is the same in the six exchange-weight-100 arms on 357 of 360 amplicons: relative spread <= 1e-6. The other three are 1777136.3__V7-V9 (1.24e-4, the maximum), 2697503.3__V4-V5 (5.3e-5) and 984206.3__V1-V3 (1.7e-5). The arms first agree on the exchange flux; the internal coefficients then choose among those solutions. The weighting therefore acts at the optimum itself. But any tolerance on the objective value is swamped by the exchange term.
- **Exchange weight 1 (exploratory).** The internal share rises to a median 50.9% (34.4% to 73.3%). real_ex1 − two_tier_ex1 at the face midpoint is −0.0086 [−0.0123, −0.0045], with the face minimum −0.0050 [−0.0080, −0.0006]. That is somewhat larger than the released objective's −0.0061, and still below the SESOI.
- **prob_exp 2.** real_exp2 − two_tier is −0.0061 [−0.0098, −0.0026], indistinguishable from prob_exp 1.

## 8. The degenerate n1 control

On 60 single-genome amplicons, every gene-backed reaction has p = 1 (p < 1 reactions: 0 on every row). Real and two-tier are therefore the same objective:

- Returned solutions: bitwise identical, 60 of 60.
- Face values: within 1e-9 on 59 of 60 (largest difference 7.2e-5).
- Opt-face counts: differ by FVA noise only (section 6).
- Flux and degeneracy contrasts: real − two_tier face mid is −6e-7 [−2e-6, +4e-12].

Against plain parsimony, the prFBA objective at n = 1 (two-tier plus the exchange term) carries slightly more truth-absent flux: real − uniform face mid +0.0008 [+0.0001, +0.0018]. That is far below the SESOI. Its degeneracy is no different: +0.37 [−0.18, +0.95].

This is the regime of 86.8% of answered amplicons on the controlled benchmark (manuscript sec:recon). On those, the weighting does nothing by construction.

## 9. How post hoc this is

- **Declared after seeing data.** The test was declared after a probe of the released code, one orientation model and a 10-amplicon smoke (design, header). It is not a pre-registration.
- **Decision rules amended after the first smoke** (design, sections 11 and 12, R1 to R3). This was done after that smoke's point contrasts had been seen, before the full run was launched:
  - **R1.** H1, H1s and C3 moved from the solver's returned vertex to the optimal face.
  - **R2.** The optimal face got one absolute slack of 1e-6, with a 1e-7 confirmation face.
  - **R3.** A slack-dominated interval was dropped.
  - **Withdrawn.** The first smoke's printed verdict ("works"; "reduces degeneracy, NOT specific to p") was withdrawn.
- **The amendment did not decide the H1 label.**
  - **Old rule.** The pre-review rule decided H1 on the point TAF_gb with the same 0.01 SESOI. That version is in the scratch copy `wobj_v2/WEIGHTED_OBJECTIVE_TEST.md`, section 10.
  - **Full-run result under it.** The point contrast is −0.0081 [−0.0117, −0.0054]: its CI excludes 0, its mean misses the SESOI, and its CI does not lie inside ±0.01. The old rule also reads "detectable but below SESOI".
  - **H2.** The H2 label does not depend on the slack between 1e-7 and 1e-6. The full run did not compute the pre-review arm-specific slack, so H2 under the old rule is unknown; on one review amplicon that slack flipped the sign of real − shuffled (design, section 12, R2).
- **No change after launch.** The full run used the declared harness unchanged (sha256 `d44e7fa3…`). Nothing in it was changed after the run started. `verdict()` applied the rules mechanically.
- **Stale design text.** The design document's header and section 14 still say the full run "has not been launched". This file reports that run.

## 10. Limitations

- **Truth is an annotation-level proxy.** A truth-absent reaction is a gene-backed reaction that the source genome's own BV-BRC-based draft does not support. It is not demonstrated absent biochemistry.
- **The EC route stands in for the function-string merge.** p is computed by the production rule on template roles granted by EC sets, not on the RAST function strings of eq:frequency. The drafts use `allow_all_non_grp_reactions=False` where production uses True (design, section 15).
- **Sample and scope.**
  - 360 amplicons (300 with n >= 2) from 54,243 candidates of one benchmark: the controlled leave-genus-out benchmark at full reference density. The deployment database's multi-genome merges are larger and more often cross-phylum.
  - One medium (complete), one carbon cap, one biomass floor, one shuffled permutation per amplicon.
  - md5s scored for more than one source genome are excluded.
- **Selections of one genome dominate.** n >= 2 covers 13.2% of answered amplicons on this benchmark. The pooled effect is an effect on that minority. On the majority the objective is two-tier exactly.
- **Union-only roles.** On 6 of 300 multi rows, 15 roles are granted by the union's EC set but by no single member's. The EC-route rule gives such a role p = 0. A reaction backed only by such roles would then cost 1 in real against 0.01 in two-tier: an artefact of the EC route, not prevalence information. On these 6 rows no gene-backed reaction carries p = 0 (`model.p_values_gene_backed` has no 0 bin), so every reaction with a union-only gene also has a member-supported gene, and the artefact changes no coefficient here.
- **The face interval bounds the total truth-absent flux, not each reaction.** What a user receives is one vertex. On 203 of 360 rows the real arm's returned TAF depends on the vertex, by up to 0.196.
- **Not tested.** Accuracy of growth, exchange or flux predictions against measurements; community models (the community branch of the released function was not run); expression data.

## 11. The released single-model path of `MSProbability.prFBA` is broken

ModelSEEDpy c2e31a0, `modelseedpy/core/msprobability.py` l.155-201. The source is byte-identical at e1254ba (design, section 2). On a gapfilled single-genome model:

1. **Bare model.** `prFBA(model, ...)` raises `TypeError: object of type 'Model' has no len()` (l.165).
2. **One-element list.** `prFBA([model], ...)` passes the list itself to `MSModelUtil`: `AttributeError: 'list' object has no attribute 'id'`.
3. **No environment.** With `environment=None` (the default), `add_medium(None)` raises `AttributeError`.
4. **Duplicate package.** With the model unwrapped, a second `ElementUptakePkg` on the same model raises `PackageError: Package with name ElementUptakePkg already in model!` (l.174), on every call.
5. **Constraint print.** The function prints the full constraint-name list to stdout on every call (l.194).
6. **Gapfilled reactions have no note.** They carry no `notes['probability']`, so the coefficient loop raises `KeyError: 'probability'` (l.187). Production writes the note on the draft, before gapfilling.
7. **Wrong objective floored.** A gapfilled model as built carries `ATPM_c0` as its objective, and the function floors whatever objective it finds. The 95% floor then lands on ATP maintenance (bio1 = 0, ATPM 314.69 on the probe model).

Patches written at step 1 (applied to neither pinned checkout; the harness ran its own copy, and the fix branch below supersedes them):

- `msprobability_single_model.patch` makes four changes. It unwraps a bare model or a one-element list, skips `add_medium` when `environment is None`, deletes the duplicate `ElementUptakePkg`, and moves the print under `printLP`.
- `msprobability_gapfilled_default.patch` (optional) gives an un-noted reaction weight 0.

Both are copied, with the test script and its log, to `data/weighted_objective/patches/` (originals in the session scratch directory `wobj/patch/`). With both, a scratch copy gives objective 18,799.709290 and bio1 1.404069 on the probe model in all three call forms (`data/weighted_objective/patches/test_patched.log`). Hazard 7 is the caller's: set `bio1` as the objective before the call.

**Fix branch:** `fix-prfba-single-model`, commit b6857a0 on c2e31a0, in the worktree `/home/freiburger/Documents/ModelSEEDpy-prfba`. It is committed but not pushed; the PR text is `PR_ModelSEEDpy_prfba_single_model.md`. According to its author's report and test logs (scratch `prfba_fix/`):

- It fixes items 1 to 5.
- Items 6 and 7 become explicit keywords whose defaults keep the released behavior: `missing_probability=None` still raises KeyError, but before the model is changed; `biomass_id=None` still floors the objective found on the model.
- On three rebuilt amplicons, called with `biomass_id='bio1'`, it reproduces the harness copy's objective exactly (relative difference 0.0).
- On c2e31a0 a two-member community also raised the `PackageError`, so the community path was broken too. The branch runs it, but MSCommunity installs its kinetic rows on its own copy of the model, not on the model prFBA solves. That remains unfixed.

**Harness compatibility.** `released_prfba()` in `scripts/weighted_objective_test.py` asserts that the c2e31a0 source text occurs verbatim. On a ModelSEEDpy carrying this fix, the harness therefore stops at startup, even with `--no-check-released`. A rerun on a fixed library needs the harness changed first.

## 12. Figures and reproduction

- `data/weighted_objective/fig_wobj_flux.{png,pdf}`: real − two_tier, real − shuffled and oracle − two_tier at the face midpoint, per stratum and pooled, with CIs and the SESOI band.
- `data/weighted_objective/fig_wobj_degeneracy.{png,pdf}`: variable gene-backed reactions on the optimal face, two_tier / uniform / real / shuffled, per stratum and pooled.
- `scripts/wobj_results.py` recomputes every contrast from `per_amplicon.jsonl` (the declared bootstrap draw, plus an independent seed as a sensitivity check: on the face measures its CI bounds move by at most 0.0006 pooled and 0.0011 within a stratum, n4+ real − shuffled face max; on the opt-face count by at most 0.12 pooled), reports its agreement with `summary.json` under `agreement_with_summary`, writes `recomputed_contrasts.json` and draws both figures. Run it with `~/Documents/py_venv/bin/python scripts/wobj_results.py`. It runs in seconds and needs no solver.
