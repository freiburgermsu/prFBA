# In-silico test of prFBA's weighted objective

Declared 2026-10-29. This is a **declared, post-hoc experiment, not a pre-registration**. Before the hypotheses and decision rules below were written down, the author had already seen three things: (i) the probe of the released code (section 2); (ii) one orientation model (genome 100053.8, truth EC set: the objective's composition and the FVA cost); (iii) a 10-amplicon plumbing smoke (section 13.1). That smoke led to the measurement changes listed in section 11. An independent review of the design and the harness (section 12) then showed that the primary flux measure and the "optimal face", as first declared, depended on the solver's choice of vertex and on an arm-specific slack. The amendments it required were made **before launch**, after the first smoke's point contrasts had been seen; so they too are post hoc. They change what is measured and how the decision rules read it; they reverse no hypothesis. A second 10-amplicon smoke on the amended harness is in section 13.2. The full run (section 14) has not been launched.

Harness: `scripts/weighted_objective_test.py`. Launch script (not executed): `scripts/launch_wobj.sh` (systemd unit `wobj-run`; it supersedes `scripts/launch_weighted_objective_test.sh`). Outputs of the full run: `data/weighted_objective/{population.json, per_amplicon.jsonl, summary.json}`; log `logs/wobj_run.log`.

## 1. Question

prFBA.tex `eq:objective` (l.278-283) minimizes

    sum_r  max(min_prob, 1 - p_r^prob_exp) |v_r|  +  sum_ex 100 |v_ex|,    biomass >= 0.95 max,

with `p` the per-function prevalence `k/n` over the `n` selected references (`eq:frequency`, l.192-197). The paper makes three claims. The selector returns one genome for 86.8% of answered amplicons on the controlled benchmark (l.210). At `n = 1` every gene-backed reaction gets the same floored coefficient, so prFBA reduces to two-level parsimony (l.212-213). And (l.291): *"Whether the weighting reduces solution degeneracy where it is active has not been measured."* This test measures it where it is active (`n >= 2`). It asks two questions:

1. Does the fractional `p` move flux off reactions the source organism lacks?
2. Does it shrink the set of alternative optima?

## 2. Does the released single-model path run? (step 1)

The released function is `MSProbability.prFBA`, ModelSEEDpy c2e31a0, `modelseedpy/core/msprobability.py` l.155-201. Its source is byte-identical in the `fix-gapfill-reversals` worktree (e1254ba; `diff` against `git show HEAD:` is empty). It was called on one gapfilled ModelSEED model: the truth EC set of 100053.8 (515 ECs, gramneg), built by the phenotype harness's `gapfill_block` with the A-C13 settings (61 gapfilled reactions, growth 4.909 on complete medium, sanity gate passed). Probe scripts and logs are in the scratch directory `wobj/step1/`.

| call | what happens |
|---|---|
| `prFBA(model, environment=complete)` | `TypeError: object of type 'Model' has no len()` at l.165 (`len(model_s_)`) |
| `prFBA([model], environment=complete)` | `AttributeError: 'list' object has no attribute 'id'`: `MSModelUtil` receives the list itself (l.165 passes `model_s_`, not `model_s_[0]`), so `msmodelutl.py` l.375 `self.id = model.id` fails |
| `prFBA([model])` (default `environment=None`) | the same AttributeError. With the list unwrapped, `add_medium(None)` raises `AttributeError: 'NoneType' object has no attribute 'get_media_constraints'` (`fbahelper.py` l.292) |
| list unwrapped (edit 1 only), `environment=complete` | `PackageError: Package with name ElementUptakePkg already in model!`: the second `ElementUptakePkg(...)` (l.174) re-registers a package already registered at l.168 |
| edits 1 + 2 (second ElementUptakePkg dropped), bio1 objective set | runs: optimal, objective 18,799.709, bio1 1.4041 (95% of the 1.4780 maximum under the carbon cap). Prints the full constraint-name list (19,348 characters) to stdout (l.194) |
| edits 1 + 2, gapfilled reactions without `notes['probability']` (production order, section 3) | `KeyError: 'probability'` at l.187 |
| edits 1 + 2, the model's objective as built | runs, but the gapfilled model carries `ATPM_c0` as its objective (`1.0*ATPM_c0 - 1.0*ATPM_c0_reverse`), so the 95% floor is put on ATP maintenance: bio1 = 0.0000, ATPM 314.69 |

So the released single-model path does not run. In both call forms it fails at l.165, and with that fixed it fails at l.174 on every call. Two further hazards make a fixed call meaningless unless the caller guards against them:

- **`notes['probability']` is missing on gapfilled reactions.** Production attaches the note at draft build, before gapfilling.
- **The function floors whatever objective the model carries.** A gapfilled model as built carries ATPM_c0, not bio1.

A further medium detail: an `MSMedia` passed as `environment` goes through `FBAHelper.convert_kbase_media`, whose default `uniform_uptake=1000` sets every medium compound's uptake to 1000. A `{cpd: uptake}` dict keeps its own bounds; the harness passes the gapfilling medium's compounds at 100.

**Patches (not applied, not committed).** Both are in the scratch directory `wobj/patch/`. Each applies with `patch -p1` to `ModelSEEDpy-reversals` (branch fix-gapfill-reversals). Both were checked with `patch --dry-run` and then applied to a copy; the result is byte-identical to the edited file.

- `msprobability_single_model.patch` is the minimal fix. It (a) accepts a bare Model or a one-element list and unwraps it, (b) skips `add_medium` when `environment is None`, (c) deletes the duplicate ElementUptakePkg build, (d) moves the constraint-name print under `printLP`. The community branch is not touched (`build_from_species_models`, `MSCommunity`, `commkinetics`). That branch was not tested. The second ElementUptakePkg build sits outside the community condition, so the same PackageError is expected there too.
- `msprobability_gapfilled_default.patch` is optional and applies on top of the first. It reads `rxn.notes.get("probability", 0)`, giving a reaction without a note weight 0: the value the builders give every gene-less reaction, and the paper's "for example a weight of 0" (l.215). It is a policy choice, so it is kept separate.

With both patches, a scratch copy of the module runs on the probe model three ways: bare Model, `[Model]`, and un-noted gapfilled reactions. All three give objective 18,799.709290 and bio1 1.404069 (`wobj/patch/test_patched.log`). With `environment=None` it returns 0 / 0, because a gapfilled model's own exchange bounds admit no uptake.

**What the harness runs instead.** `prfba_released_copy()` in the harness is a literal copy of l.155-201. It keeps every bound, constraint and coefficient line verbatim, with only edit (1), edit (2) and the dropped print. The caller sets `bio1` as the objective before the call; this deviation from the model as built is stated in the code. On every amplicon the harness also re-solves the real arm through the **released source itself**. That source is read from the imported module and exec'd with edits (1) and (2), each asserted to occur exactly once; its sha256 is recorded. The two objective values must agree to a relative 1e-6. In the first smoke they agreed exactly on 10 of 10 amplicons (`rel_diff` 0.0); the second smoke is in section 13.2.

## 3. How `notes['probability']` is populated in production (step 2)

**Per function** (`codiffusion_bioreactor/bvbrc_to_kbase_genome.py`, `create_synthetic_genome`, l.1108-1197): a function seen for the first time creates a feature with `'probability': 1` (l.1137). A later genome carrying the same function increments it once (`elif function not in genome_functions: functions[function]['probability'] += 1`, l.1177-1180). `genome_functions = {}` is reset for each genome (l.1108). Finally `functions[function]['probability'] /= num_genomes` (l.1188), and the result is written onto the feature and its CDS. This is `p = k/n`.

**Per reaction**, `codiffusion_bioreactor/build_prfba_community_models.py` `_build_one` (l.356-362):

```python
probs = {f.get("id"): f.get("probability", 0) for f in g.get("features", [])}
m = MSBuilder(msg, _W["tpl"][name], name=rid).build_base_model(rid, "0", True, False)
...
for r in m.reactions:
    p = max((probs.get(gg.id, 0) for gg in r.genes), default=0.0) if r.genes else 0.0
    r.notes["probability"] = str(p)
```

`load_genomes_to_modelseedpy.py` applies the same rule (l.280 `feature.get('probability', 0)`, l.343 `probability = max(gene_probs)`, l.348 `probability = 0.0` for reactions without genes, stored as `str`). The deck states the same rule: `presentations/make_deck.py` l.1174, "Reaction weight = max over its genes". The note is written on the **draft** (the `build` stage, before `mediagf`), so gapfilled reactions carry none. No production code in these repositories calls `MSProbability.prFBA` (grep for `prFBA(` / `MSProbability` in codiffusion_bioreactor finds only an import in `util.py`).

**The rule in this test.** The harness's features are template roles: `phenotype_arm.genome_from_ec` makes one feature per granted role, with the role id as feature id. The smoke confirmed that every model gene id is a template role id (0 exceptions over 10 models). The EC-route analogue of "genome `g` carries function `x`" is "every EC in role `x`'s name is in `EC(g)`". So:

- `p(role) = #{members g : ECs(role) ⊆ EC(g)} / n`
- `p(R) = max over R's genes of p(gene)`, and `0.0` for a reaction without genes. This is stored as `str` in `rxn.notes['probability']`.

The rule is applied after gapfilling, so gapfilled reactions get `0.0`: the value the production rule gives every gene-less reaction, and the paper's "for example a weight of 0". Their coefficient is therefore `max(0.01, 1 - 0) = 1`. A role granted by the union but by no single member would get `p = 0`; this happened 0 times in the smoke (`n_union_only_roles`).

## 4. Hypotheses

The flux hypotheses are decided on each arm's **optimal face**, not at the vertex the solver happens to return (section 12, R1). For one arm on one amplicon:

- `F_lo`, `F_hi`: the minimum and maximum of the total truth-absent gene-backed flux over the arm's optimal face (objective <= opt + 1e-6, the same absolute slack for every arm).
- `D`: the two-tier arm's gene-backed flux at its returned optimum. It is one constant per amplicon, the same for every arm, and it is itself face-robust (section 8).
- `TAF_face_min = F_lo / D`, `TAF_face_max = F_hi / D`, `TAF_face_mid = (F_lo + F_hi) / 2D`.

Hypotheses:

- **H1 (primary).** Where `n >= 2`, the real-p objective's optimal face carries less truth-absent flux than the two-tier objective's on the same model: real − two_tier < 0 on `TAF_face_mid`, with the sign holding at both face endpoints (`TAF_face_min` and `TAF_face_max`).
- **H1s (specificity).** Real also beats shuffled, which has the same coefficient multiset placed at random, on `TAF_face_mid`. So the effect comes from the information in `p`, not from the gradation alone.
- **H2 (degeneracy).** On the optimal face, real leaves fewer gene-backed reactions free to vary (range > 1e-3) than two-tier. A 1e-7 face confirms the count. Specificity is checked the same way, real vs shuffled.
- **H3 (descriptive, the diagnostic the paper needs).** The objective's composition at the optimum: the internal weighted term against the exchange term, and whether the exchange term is the same in every arm (then the objective acts lexicographically: exchange first, internal weights second). Also the addressable share of truth-absent flux (section 8).
- **Controls.**
  - C1, determinism: at `n = 1`, real and two-tier have identical coefficients, so their returned objectives and fluxes must be bitwise equal (their face values are checked too).
  - C2, faithfulness: `released_check` (section 2).
  - C3, positive control: the oracle arm, on the same face measures. It shows whether *any* reweighting of gene-backed reactions can move truth-absent flux under this objective.

## 5. Population

**Source.** `data/cells.csv.gz`, arm `prfba_exact_tie`, space `ec`, density 1.0, outcome `scored`: the controlled shared-reference benchmark. It uses leave-genus-out folds (`fold_plan.json`), the shipped exact-tie rule, and full density. The fold runs are `data/prfba/f0{0..4}_d100/predictions.json`, mode `fold`, not forced. They hold the per-cell selections at `pred.exact_tie[md5].reps`, with the union EC set at `.ec`. The member and truth EC sets come from `data/_cache_panel_ec.json` through `norm_ec`, intersected with `data/panel_ec_universe.json` (2,235 ECs). These are the same sources `calibrate_pfunc.py` uses (`layerA_job`: `rec["reps"]`, `E[ri]`). Templates come from the truth genome's `panel_eligible.json:template`.

Counts below were computed by the harness's population stage, written to `population.json` `meta`.

| | n1 | n2 | n3 | n4+ | total |
|---|---|---|---|---|---|
| scored cells | 49,048 | 5,693 | 1,313 | 967 | 57,021 |
| candidates (truth >= 200 ECs) | 46,548 | 5,494 | 1,272 | 929 | 54,243 |
| eligible source genomes (twin md5 excluded) | 7,521 | 2,474 | 663 | 410 | |
| **sampled** | **60** | **100** | **100** | **100** | **360** |

**Exclusions.**

- 2,778 scored cells whose truth genome has fewer than 200 ECs (`min_truth_ec`, the phenotype harness's threshold).
- Cells of the 3,621 md5s scored for more than one source genome ("twins"). A twin's model would be identical while its truth differs. This restricts the **target population**; it is not needed as a correction, since the sampler already takes at most one amplicon per md5. The review (section 12, O1) measured what it removes, from `cells.csv.gz` with truth >= 200 ECs: 632 of 3,106 n2, 218 of 881 n3 and 188 of 598 n4+ source genomes (20%, 25%, 31%). Those totals agree with the table above (3,106 − 632 = 2,474; 881 − 218 = 663; 598 − 188 = 410). Twin and non-twin mean F1 are close (n2 0.796 vs 0.791, n3 0.798 vs 0.793, n4+ 0.780 vs 0.776), so the restriction probably biases little. Including twins would move the pooled weights from 0.6975 / 0.1869 / 0.1156 to 0.677 / 0.192 / 0.130.

**Integrity**, over all 54,243 candidates: the union of the members' EC sets equals the prediction's EC set in every case (0 mismatches). The recomputed union-vs-truth F1 matches `cells.csv` f1 within 1e-3 in every case (0 disagreements). The truth EC count equals `cells.csv` `n_truth` in every case (0). Leakage (source genome among the reps) occurred 0 times.

**Sampling.**

- At most one amplicon per source genome and per md5 across the whole sample.
- Strata are drawn rarest first (n4+, n3, n2, n1), under one global one-per-source-genome rule. A source genome drawn in n4+ or n3 is then not available to n2, which slightly under-represents, in n2, source genomes that also have n3 or n4+ cells.
- Within a stratum, the sorted source genomes are shuffled with `random.Random("1729:wobj:<stratum>")`.
- Each source genome contributes its V4 cell if it has one, else the lexicographically first region.
- Jobs are interleaved n2, n3, n4+, n1, so `--limit k` takes a stratified prefix.

**Size.** The population is sized to the budget, about 2 h on 8 workers as first planned (section 14 has the amended estimate). n1 is the degenerate control: it tests C1, and it measures real vs uniform where the paper says prFBA runs almost always.

## 6. Model construction (identical for every arm of an amplicon)

1. Build the gapfilled model of the **union** EC set with `phenotype_fold.gapfill_block(m, union, template, id, complete_media, cpds, "as_shipped", None, need_ess=False)`. This is the phenotype harness's gapfilled pass verbatim: draft → ATP correction on its own selected media → MSGapfill on complete medium → gapfiller reversals applied (`--integrate-reversals`, A-C12/A-C13). Templates: `ModelSEEDTemplates-c12run/templates/v7.0` (commit 21a8084). ModelSEEDpy: the pinned clone c2e31a0 first on PYTHONPATH, enforced by `--expect-modelseedpy c2e31a0`.
2. Gates. A model that fails `phenotype_arm.sanity_gate` (biomass or ATP from water and protons) is `gated`. A model with no gapfill solution, or no growth on complete medium, is `unscoreable`. Neither enters any statistic; both are counted.
3. Truth support: the truth genome's base model, with the same template and the same `build_base_model(..., allow_all_non_grp_reactions=False, annotate_with_rast=False)` call, no gapfilling. A prediction reaction is classed by that model:
   - **truth-supported**: gene-associated there.
   - **truth-absent**: gene-backed in the prediction but not supported.
   - **unannotated**: gene-less in the prediction (gapfilled or spontaneous) and not supported. These are never counted as truth-absent.
4. The prFBA constraint set, through `prfba_released_copy`, on a fresh `model.copy()` per arm:
   - the complete-medium compounds at uptake 100 (dict environment);
   - one `ElementUptakePkg` carbon cap of 100;
   - `bio1` set as objective, then `min_value` >= 0.95 × max.

   The constraint set is the same for every arm; only the coefficients differ.

## 7. Arms

| arm | p per rxn reaction | prob_exp | ex_weight | role |
|---|---|---|---|---|
| `two_tier` | gene-backed 1, gene-less 0 | 1 | 100 | what production gives at n = 1; the primary comparator |
| `real` | production rule (section 3) | 1 | 100 | the method |
| `uniform` | all 1 (coefficient 0.01 everywhere: plain parsimony plus the exchange term) | 1 | 100 | the task's "uniform"; conflates the gene-support tier with the gradation |
| `shuffled` | real p permuted among the model's gene-backed reactions (`rng_for(1729, key, "shuffle")`) | 1 | 100 | specificity |
| `oracle` | gene-backed: 1 if truth-supported else 0; gene-less 0 | 1 | 100 | positive control (headroom) |
| `real_exp2` | real p | 2 | 100 | sensitivity |
| `real_ex1`, `two_tier_ex1` | real / two-tier | 1 | 1 | exploratory: the weighting with the exchange term no longer dominant |

## 8. Metrics (per arm)

`|v|` is the forward plus the reverse primal: the quantity the objective prices.

**Primary, vertex-robust (every arm; 2 LPs each).**

- The **optimal face** is objective <= opt + 1e-6, with one absolute slack (`OPT_SLACK`) for every arm and every amplicon. On it, `F_lo` and `F_hi` are the minimum and maximum total truth-absent gene-backed flux.
- `D` is the two-tier arm's gene-backed flux at its returned optimum: one per-amplicon constant, the same for every arm. Its own min and max over the two-tier optimal face are recorded (`D_face_rel_width`), so its face-robustness is checked on every row.
- `taf_face_min = F_lo / D`, `taf_face_max = F_hi / D`, `taf_face_mid = (F_lo + F_hi) / 2D` (the decision measure), `taf_face_width = (F_hi − F_lo) / D`.
- At slack `s`, futile forward-plus-reverse flux on the cheapest truth-absent reaction (coefficient `c_min`) raises `F_hi` by `s / c_min` at no other cost. That increase is recorded per face as `futile_offset_bound` (it is a lower bound on what the slack adds to `F_hi`). On the optimal face it is at most 1e-6 / 0.01 = 1e-4 flux units.

**Secondary, "as returned" (the solver's default vertex; what a user of prFBA receives).**

- **TAF_gb**: flux through truth-absent gene-backed reactions ÷ flux through all gene-backed reactions, both of this arm. `taf_point_D`: the same numerator ÷ `D`. TAF_all: the numerator ÷ flux through all rxn-prefixed reactions. `flux_absent` itself is reported too.
- **Vertex spread.** The arm's LP is re-solved from scratch (`gurobipy.Model.reset`) with Method 0 (primal simplex), 1 (dual simplex) and 2 (barrier with crossover). Per re-solve: objective, gene-backed and truth-absent flux, TAF_gb, `taf_point_D`. `taf_gb_vertex_spread` is max − min over the default vertex and the three re-solves; `objective_vertex_rel_diff_max` checks that every re-solve reached the same optimum. Method is reset to automatic before any face LP.
- Flux precision, gene-backed on both sides: (gene-backed − truth-absent flux) ÷ gene-backed flux, i.e. 1 − TAF_gb. Count precision: truth-supported ÷ all active gene-backed reactions (|v| > 1e-6). Gene-less truth-supported reactions (gapfilled, spontaneous) are on neither side.
- `n_active`, split into all, gene-backed, truth-absent and gapfilled.
- Growth: bio1 at the optimum against the maximum. It must equal 0.95 × max; this is a sanity check.

**Objective composition (H3).** The internal term and the exchange term, and their shares. The internal term is further split into gene-backed p = 1, gene-backed p < 1, gapfilled and gene-less draft reactions; each part uses that arm's own p. Per row, `exchange_term_rel_spread_ex100` is (max − min) ÷ max of the exchange term over the six ex_weight-100 arms.

**Addressable share** (the interpretation key): truth-absent flux on p < 1 reactions ÷ all truth-absent flux, in the two-tier arm at its returned vertex. The weighting can only move flux off p < 1 reactions. Where this share is small, the weighting is inert whatever p's calibration. Also reported: the static precision of p < 1 against p = 1 reactions.

**Degeneracy** (FVA arms two_tier, real, uniform, shuffled):

- **Feasible scan** (once per amplicon, on the shared constraint set): min/max of every gene-backed reaction. A reaction fixed there is fixed on every arm's face, so it is skipped.
- **opt face** (1e-6, as above): the number of gene-backed reactions with range > 1e-3 (**H2 primary**), the sum of ranges, and the counts among p < 1 and truth-absent reactions. The count at range > 1e-6 is stored but never read: it does not plateau in the slack (section 12, R2).
- **opt7 face** (1e-7): the same FVA, as confirmation. A row is flagged (`opt7_count_flag`) when |count(opt) − count(opt7)| ÷ max(1, count(opt)) > 0.05; `summary.checks.opt7_count_flags` counts them per arm.
- **int1e-3 face**: one absolute slack per amplicon, 0.1% of the two-tier arm's internal term, for every FVA arm. Descriptive: sum of ranges. Its truth-absent interval is stored but not read, because it is mostly the futile offset (section 12, R3).
- **obj1e-3 face**: 0.1% of the whole objective, real arm only. It is the diagnostic that a tolerance on the objective is swamped by the exchange term.
- The other arms (oracle, real_exp2, real_ex1, two_tier_ex1) get the optimal-face interval only, no FVA.
- Descriptors: n_reps, the p-value histogram of gene-backed reactions, model sizes, the gapfill block, and the L1 flux distance between arms at their returned vertices.

## 9. Statistics

- Paired by amplicon: real − two_tier (primary), real − shuffled, real − uniform, oracle − two_tier, real_exp2 − two_tier, real_ex1 − two_tier_ex1, shuffled − two_tier.
- Per stratum (n2, n3, n4+; n1 as the control): the mean and median paired difference; the fractions < 0, = 0 and > 0; Wilcoxon signed-rank with `zero_method="pratt"` (exact zeros, i.e. identical optima, are expected); and the matched-pairs rank-biserial correlation over the non-zero differences.
- **CIs** are a 95% percentile cluster bootstrap, 2,000 draws, with score_folds' draw (`default_rng(seed).integers(0, K, (B, K))` over the sorted cluster keys). The seed is 1729 + 1000 × stratum index (n2 0, n3 1, n4+ 2, n1 3), so the strata draw independently. Within a stratum, every contrast and metric shares the same draws. The primary unit is the holdout group (`group_id`), the benchmark's pre-specified convention. The source genome is the sensitivity unit; with one amplicon per source genome it is an ordinary bootstrap.
- **Pooled multi (n >= 2)** = Σ_s w_s × mean_s, with w_s each stratum's share of eligible source genomes: n2 0.6975, n3 0.1869, n4+ 0.1156 (2,474 / 663 / 410). Its CI combines independent within-stratum cluster draws with the same weights. The comparator mean that H2's relative reduction divides by is pooled with the same weights (`mean_comparator_weighted`). The Wilcoxon test on the pooled sample is unweighted, and is labelled so.
- **Multiplicity.** Each verdict is one pooled test per hypothesis. H1 "works" requires the midpoint condition, both endpoint conditions and the specificity condition at once: an intersection-union test, each component at the 95% level, which needs no correction. BH serves the per-stratum reporting only: `taf_face_mid` × {real − two_tier, real − shuffled} × {n2, n3, n4+}, 6 Wilcoxon p-values (`bh_q_per_stratum_family`).

## 10. Decision rules (applied mechanically by `verdict()` in the harness)

**SESOI.** 0.01 absolute in TAF_face, i.e. one percentage point of `D`, the two-tier gene-backed flux. For degeneracy, a 10% relative reduction in the opt-face count.

**Gate.** C1 (every n1 row's real and two-tier returned solution identical) and C2 (every `released_check` passes) must hold. If either fails, or `released_check` was not run, `verdict()` returns `harness_defect` with the reason and reads no other number.

All rules below use the pooled multi contrast and its group-cluster 95% CI.

- **C3, positive control first.** Oracle − two_tier: is the `taf_face_mid` CI upper bound < 0 with mean <= −0.01, and are the `taf_face_min` and `taf_face_max` CI upper bounds < 0? If not, the objective cannot express a reweighting of gene-backed reactions. A null on H1 is then a *structural* null about the objective, not evidence about `p`.
- **H1 "the weighting works"**: real − two_tier `taf_face_mid` CI upper bound < 0 and mean <= −SESOI; **and** the `taf_face_min` and `taf_face_max` CI upper bounds < 0 (the sign is not a vertex choice); **and** real − shuffled `taf_face_mid` CI upper bound < 0 (H1s).
- **H1 "harmful"**: the mirror image. Midpoint CI lower bound > 0 and mean >= +SESOI, and both endpoint CI lower bounds > 0.
- **H1 "does nothing" (equivalence)**: the midpoint and both endpoint CIs lie entirely inside (−0.01, +0.01), and the midpoint CI contains 0. If they lie inside but the midpoint CI excludes 0, the label is "equivalent within SESOI, nonzero". A non-significant p-value is not evidence of no effect.
- Otherwise, in this order:
  - "lowers truth-absent flux by >= SESOI at the face midpoint, but" not at both face endpoints (vertex-dependent), and/or not specific to p;
  - "raises ... at the face midpoint, but not at both face endpoints";
  - "detectable but below SESOI": the midpoint CI excludes 0, |mean| < SESOI, and the CIs do not lie inside ±SESOI;
  - "inconclusive": the midpoint CI crosses 0 and a SESOI bound.
- **H2**: real − two_tier on the opt-face count (`faces.opt.n_var_1e-3`). "Reduces degeneracy" if the CI upper bound < 0 and the relative reduction (−mean ÷ the stratum-weighted two-tier mean) >= 10%; "below the 10% SESOI" if only the CI condition holds; "increases" if the CI lower bound > 0; otherwise "no detectable change". In every case the label says whether real is also below shuffled ("specific to p") or not. The same rule is applied to the opt7 count (`H2_at_1e-7_face`); `H2_confirmed_at_1e-7_face` is false when the two labels differ, and the H2 result is then reported as slack-sensitive.
- Per-stratum results are reported beside the pooled verdict. The stratum results and the oracle headroom captured (Δ(real − two_tier) ÷ Δ(oracle − two_tier) on `taf_face_mid`; the point-TAF version beside it) say *where* and *how much*. The verdict says *whether*.
- Reported but outside the verdict: the point ("as returned") contrasts and their vertex spread, real − uniform (whole objective vs plain parsimony), real_exp2, the ex_weight-1 pair, and objective composition.

## 11. Changes made after the first smoke (measurement validity, not direction)

1. **Near-optimal slack.** It was 0.1% of each arm's *own* internal term. That handed the real arm, whose internal term is larger, a looser face. It is now one absolute slack per amplicon, 0.1% of the two-tier arm's internal term (the int1e-3 face).
2. **Degeneracy primary.** It moved to the count on the numerically optimal face. On any face with a positive slack, every reaction that can move near the optimum moves, so the count of variable reactions was **identical for every arm on 10 of 10 smoke amplicons** (e.g. 365/365/365/365 for two_tier/real/uniform/shuffled). On those faces only the sum of ranges carries information. *(Amended by section 12, R2: the "numerically optimal face" was then defined with an arm-specific slack, max(1e-6, 1e-6 × the arm's internal term). It is now one absolute slack of 1e-6 for every arm.)*
3. **Count threshold.** The count threshold is range > 1e-3, not > 1e-6. One count at 1e-6 changed by 1 between two otherwise identical runs whose only difference was the order of the FVA LPs.
4. **Cost.** A feasible-region scan now skips fixed reactions. Face measurements are limited to the arms the decision rules read: two_tier, real, uniform and shuffled; obj1e-3 on real only. Point metrics were unchanged by this (verified on 5 amplicons).

## 12. Review (pre-launch amendments, 2026-10-29)

An independent reviewer re-built two smoke amplicons with the harness's own steps (2021369.7__V7-V9, n3; 1871331.13__V4, n4+) and re-solved each arm from scratch. Probe script and outputs: scratch `wobj_review/slack_probe.py`, `probe_2021369.json`, `probe_1871331.json`. The verdict was "needs changes", with four required changes. Each is listed below with its evidence and with what was done. The amended harness reproduces the probe on 2021369.7 (section 13.2, test row `wobj_v2/t1`).

### R1. The flux hypotheses were decided at one arbitrary vertex

**Finding.** H1, H1s and C3 read TAF_gb at the vertex the solver returned, and that vertex is arbitrary among alternative optima.

- 2021369.7: two-tier TAF_gb is 0.0166 (truth-absent flux 12.92) at the harness's default solve. Under primal, dual and barrier it is 0.00521, 0.00444 and 0.00495, with the objective identical (15732.056167). That amplicon's real − two_tier contrast goes from −0.0139 to about −0.002.
- 1871331.13: real TAF_gb is 0.365 by default and 0.243 under primal simplex (same objective, 20859.214531). real − two_tier goes from −0.0012 to −0.123.
- In the first smoke's rows, the two-tier vertex sits at 0.89–1.00 of its own optimal-face truth-absent interval on 5 of 8 multi rows (0.905, 0.893, 0.924, 0.991, 1.000; recomputed from `wobj/smoke3/per_amplicon.jsonl`). So the point contrast is biased by the vertex choice.
- These per-amplicon shifts are the size of the SESOI (0.01) and of the first smoke's pooled effect (−0.0139): the verdict could be set by the solver path.

**Done.**

- H1, H1s and C3 now read the truth-absent flux interval on the optimal face, `F_lo` and `F_hi`, divided by the arm-independent constant `D` (section 4, 8). `verdict()` decides on `taf_face_mid` and requires the sign at both endpoints (section 10).
- The oracle and the three other non-FVA arms now get the face interval too (2 LPs each, no FVA).
- The point TAF is kept as the secondary "as returned" measure, with its spread over Method 0/1/2 re-solves (`vertex_resolves`, 3 LPs per arm; about 0.2 s per arm in the test row).
- Reproduction on 2021369.7 (`wobj_v2/t1/per_amplicon.jsonl`): the two-tier re-solves give TAF_gb 0.00521 / 0.00444 / 0.00495 against 0.0166 at the default vertex (spread 0.0122), as in the probe. On the face, two-tier is [0.00419, 0.01786] (mid 0.01103) and real [0.00271, 0.00347] (mid 0.00309), so real − two_tier is −0.0079 at the midpoint, −0.0015 at the minimum and −0.0144 at the maximum. `D_face_rel_width` is 1.3e-7: `D` does not depend on the vertex.

### R2. The "opt face" was neither optimal nor common to the arms

**Finding.** Its slack was max(1e-6, 1e-6 × the arm's own internal term): 1.5e-4 for two_tier, 1.85e-4 for real, 2.85e-4 for shuffled and 9.2e-6 for uniform on 2021369.7. With a 1e-3 range threshold, a slack `s` counts any reroute whose marginal cost is below s / 1e-3, about 0.15 per unit of flux: that is up to 15 extra p = 1 steps, not alternative optima. The count of gene-backed reactions with range > 1e-3, at slack (arm-specific | absolute 1e-8 | 1e-7 | 1e-6 | 1e-5 | 1e-4):

| amplicon | arm | arm-specific | 1e-8 | 1e-7 | 1e-6 | 1e-5 | 1e-4 |
|---|---|---|---|---|---|---|---|
| 2021369.7 | two_tier | 224 | 76 | 74 | 75 | 114 | 223 |
| | real | 127 | 29 | 28 | 29 | 48 | 119 |
| | shuffled | 114 | 29 | 29 | 29 | 59 | 102 |
| 1871331.13 | two_tier | 211 | 52 | 51 | 51 | 93 | 198 |
| | real | 72 | 25 | 24 | 25 | 28 | 48 |
| | shuffled | 127 | 17 | 17 | 17 | 20 | 54 |

- The counts plateau between 1e-8 and 1e-6. The harness slack sat above the plateau and inflated the counts 3–7 times.
- On 1871331.13, real − shuffled is −55 at the harness slack but +8 on the plateau: the H2 specificity call flips sign.
- The count at range > 1e-6 never plateaus (two-tier on 2021369.7: 111 / 223 / 242 at 1e-8 / 1e-7 / 1e-6).
- 1e-8 gave 2 LP failures on 1871331.13.
- The opt-face sum of ranges and truth-absent interval are insensitive to the slack (two-tier on 2021369.7: sum 119.5–122.1; interval [3.25–3.265, 13.91–13.94]).

**Done.**

- The optimal face is now objective <= opt + 1e-6 (`OPT_SLACK`), one absolute slack for every arm.
- A 1e-7 confirmation face (`opt7`) is added for the four FVA arms. Rows whose counts differ by more than 5% are flagged, and the H2 rule is also applied to the opt7 count (section 10).
- H2 stays on the range > 1e-3 count at the 1e-6 face. The range > 1e-6 count is stored but is in no decision and in no summary contrast.
- Test row on 2021369.7: counts at 1e-6 | 1e-7 are two_tier 75 | 76, real 28 | 29, uniform 74 | 75, shuffled 29 | 29 (the probe at 1e-6: 75, 29, 75, 29). This is within the ±1 count noise of section 11 item 3. 0 face-LP failures.

### R3. `absent_face_width.int1e-3` was mostly the slack

**Finding.** At slack `s` the LP can add futile forward-plus-reverse flux on the cheapest truth-absent reaction (coefficient 0.01), which adds s / 0.01 of "absent flux".

- Real on 2021369.7: opt-face max 2.7184 vs int1e-3 max 17.7456, a difference of 15.03, against 0.15046 / 0.01 = 15.05.
- Real on 1871331.13: 358.958 vs 378.296, a difference of 19.34, against 0.19363 / 0.01 = 19.36.
- two_tier and shuffled show the same offset. The level is therefore mostly a function of the slack, not the near-optimal freedom of the absent flux.

**Done.**

- `absent_face_width.int1e-3` is removed from the metrics and the summary. The int1e-3 face remains for its sum of ranges only.
- Every face now records `futile_offset_bound = s / c_min`, with `c_min` the arm's cheapest truth-absent coefficient. In the test row, real's int1e-3 maximum is 17.7456 with an offset of 15.0457. On the 1e-6 face the offset is 1e-4 flux units, against `D` of 778.8.
- The doc no longer describes any positive-slack interval as "the width of the truth-absent flux range".

### R4. Declare, implement, re-smoke, withdraw

**Done.**

- R1–R3 are declared here and implemented in `face_block`, `vertex_resolves` (new), `arm_metrics`, `FACES`, `METRICS` and `verdict()`.
- Re-smoke: section 13.2 (`--limit 10`, new scratch outdir `wobj_v2/smoke_final`, same 10 amplicons as the first smoke).
- Withdrawn: the first smoke's readout "on the smoke, then, the shrinkage of the optimal face comes from the gradation, not from `p`'s information", and its printed verdict "works" / "reduces degeneracy, NOT specific to p". Both rested on the measures R1 and R2 invalidate.

### Optional improvements

| # | suggestion | decision |
|---|---|---|
| O1 | twin-md5 exclusion not needed (the sampler already takes one amplicon per md5); keep twins with a seeded pick, or state the restriction | **stated, not changed.** The exclusion is now declared part of the target population (section 5), with the review's measure of what it removes. Changing it would have changed the sampled population after the first smoke. |
| O2 | precision_flux / precision_count asymmetric | **adopted**: both are gene-backed on both sides (section 8) |
| O3 | TAF_gb's denominator depends on the arm; report flux_absent | **adopted** as part of R1 (`D`; `taf_point_D`; `flux_absent` in the metrics) |
| O4 | `verdict()` does not gate on C1/C2 | **adopted**: `harness_defect` (section 10) |
| O5 | equivalence branch precedes "detectable but below SESOI" | **adopted**: "equivalent within SESOI, nonzero" |
| O6 | H2 relative reduction used the unweighted base | **adopted**: `mean_comparator_weighted` |
| O7 | BH vs bootstrap verdict | **adopted**: section 9, multiplicity; BH now covers the 6 per-stratum tests only |
| O8 | 3–5 shuffled permutations | **declined**: one permutation per amplicon; its noise averages over the about 300 multi amplicons and enters the bootstrap CI. Stated in section 15. |
| O9 | H3 framing: exchange term identical across the ex_weight-100 arms | **adopted**: `exchange_term_rel_spread_ex100` per row and in `checks` (section 8, 13.2) |
| O10 | production drafts use `allow_all_non_grp_reactions=True` | **documented** (section 15) |
| O11 | rarest-first drawing under the global one-per-source rule | **documented** (section 5) |
| O12 | cost of the fixes | measured on the second smoke (section 13.2, 14) |

## 13. Smoke (plumbing, 10 amplicons; these are not results)

Neither smoke is a result. Each has 10 amplicons (n1 2, n2 3, n3 3, n4+ 2), so 8 multi rows and at most 3 holdout-group clusters per stratum. They show that the harness runs and that the readouts have dynamic range.

### 13.1 First smoke (harness before the review; superseded)

Run: `--limit 10 --workers 5`, two batches of 5 under PYTHONPATH=c2e31a0, in the scratch directory `wobj/smoke3/`. These plumbing facts do not depend on the measures the review changed, and the second smoke reproduces them:

- All 10 `ok`. None gated (sanity gate passed on all). Every `released_check` passed (`max_rel_diff` 0.0). C1 holds: 2 of 2 n1 rows identical. 0 genes outside the template roles, 0 union-only roles.
- **Objective composition.** The internal weighted term is **0.89-1.34%** of the real arm's objective; the exchange term is the rest. Two-tier: 0.87-1.21%. Uniform: 0.02-0.06%. With ex_weight 1, the internal share rises to 0.47-0.57%. Within the real arm's internal term, gene-less draft reactions at coefficient 1 contribute 51-92%. Gene-backed p = 1 reactions contribute 1-5%, and gene-backed p < 1 reactions 0-35%.
- **The 0.1%-of-objective tolerance is vacuous.** On the real arm, the obj1e-3 face has a sum of ranges of 19,399-110,660 against 0.9-635 on the opt face as then defined.
- **Addressable share** (two-tier, multi rows): 0.42-1.00 on 7 of 8 multi amplicons. It is 0.0 on one, 84109.4, which has only 20 p < 1 reactions. Mean 0.66. The share is 0 by construction on n1.
- **Static precision of gene-backed reactions**: p < 1 mean 0.48 against p = 1 mean 0.90 (multi rows). This is consistent with the EC calibration (`data/calibration/README.md`: p = 1/2 present 0.449, p = 2/2 present 0.828).
- **Point contrasts (now secondary).** real − two_tier TAF_gb at the returned vertex: weighted pooled −0.0139, CI [−0.0212, −0.0024]. The second smoke returns the same vertices and the same value.
- **Withdrawn** (section 12, R4): the readout "on the smoke, then, the shrinkage of the optimal face comes from the gradation, not from `p`'s information", and the printed verdict "works" / "reduces degeneracy, NOT specific to p". Both rested on the vertex-dependent point TAF (R1) and on the arm-specific opt-face slack (R2).
- Time and memory of that harness: wall mean 103.8 s per amplicon (57.5-156.8), of which gapfill 76.7 s and arms 24.1 s; 1.595 GiB per worker; process tree peak 8.42 GiB with 5 workers.

### 13.2 Second smoke (amended harness)

**Run.** The same 10 amplicons: the first smoke's `population.json` was copied into the new outdir, and the harness's spec check accepted it. PYTHONPATH=`/home/freiburger/Documents/ModelSEEDpy-c2e31a0`, `--expect-modelseedpy c2e31a0` (clean checkout at c2e31a0c). Two batches, `--limit 5 --workers 5` and then `--limit 10 --workers 5`, in the scratch directory `wobj_v2/smoke_final/`. Harness sha256 `d44e7fa3...` (`harness.sha256`); the second batch's `summary.json` covers all 10 rows.

**Run before it.** The amended harness without the LP retry (sha256 `39aa834d...`) ran the same 10 amplicons first (`wobj_v2/t1`, `wobj_v2/g`, `wobj_v2/smoke`). 4 single FVA LPs failed there, one on a 1e-6 face and three on 1e-7 faces. That led to one cold retry per failed LP (`_solve_retry`). Otherwise the two runs agree: the largest absolute difference in objective, TAF_gb and face values over all arms and rows is 1.0e-10, and no opt or opt7 count differs by more than 1.

**Checks** (`summary.json` `checks`):

- 10/10 `ok`, 0 gated, 0 unscoreable.
- C2: `released_check` 10/10, `max_rel_diff` 0.0, one released-source sha256 (`6eeecfb7...`).
- C1: 2/2 n1 rows bitwise identical at the returned vertex. Their face values agree within 1e-9 (max abs difference 1.1e-11).
- Face LPs: 0 failures after the retry; 4 LPs needed it.
- Vertex re-solves: all 8 arms on all 10 rows (0 skipped); every re-solve reached the same objective (max relative difference 3.6e-14). The returned TAF_gb depends on the vertex on 8 of 10 rows for both real and two_tier (spread > 1e-6); the largest spread is 0.122 for real and 0.045 for two_tier. The oracle's TAF_gb never moves (0 of 10).
- `D` is face-robust: the largest `D_face_rel_width` is 9.0e-7.
- opt7 confirmation: 2 of 40 FVA arm-rows flagged, uniform on 76935.14 (90 vs 84) and shuffled on 3125783.3 (22 vs 24). Real and two_tier: 0 flagged. The H2 label is the same at 1e-7.
- **Exchange term (H3).** It is the same in all six ex_weight-100 arms on 10 of 10 rows (largest relative spread 3.5e-14). So the objective acts lexicographically: the arms first agree on the exchange term, and the internal weights only choose among the solutions that minimise it. The internal term is 0.89-1.34% of the real arm's objective. That share matters for any tolerance-based face, not for the optimum itself.

**Per amplicon** (multi rows; face interval [`taf_face_min`, `taf_face_max`]; counts = range > 1e-3 on the opt face):

| amplicon | stratum | two_tier face | real face | shuffled face | oracle mid | TAF_gb as returned, two_tier → real | count two_tier / real / shuffled / uniform |
|---|---|---|---|---|---|---|---|
| 2021369.7__V7-V9 | n3 | [0.0042, 0.0179] | [0.0027, 0.0035] | [0.0040, 0.0233] | 0.0027 | 0.0166 → 0.0026 | 75 / 28 / 29 / 74 |
| 1235557.6__V3-V4 | n2 | [0.0572, 0.0692] | [0.0432, 0.0540] | [0.0585, 0.0696] | 0.0318 | 0.0680 → 0.0423 | 40 / 33 / 27 / 41 |
| 76935.14__V7-V9 | n4+ | [0.0121, 0.0273] | [0.0283, 0.0344] | [0.0132, 0.0245] | 0.0059 | 0.0262 → 0.0304 | 86 / 36 / 38 / 90 |
| 2703885.11__V7-V9 | n2 | [0.2332, 0.2502] | [0.2138, 0.2220] | [0.2385, 0.2482] | 0.2132 | 0.2341 → 0.2086 | 79 / 62 / 42 / 78 |
| 3125783.3__V4-V5_944R | n3 | [0.0217, 0.0300] | [0.0003, 0.0004] | [0.0217, 0.0222] | 0.0002 | 0.0218 → 0.0003 | 27 / 20 / 22 / 26 |
| 2715752.3__V4 | n3 | [0.1737, 0.1880] | [0.1737, 0.1737] | [0.1744, 0.1744] | 0.1737 | 0.1737 → 0.1737 | 20 / 18 / 12 / 20 |
| 84109.4__V3-V4 | n2 | [0.3338, 0.3790] | [0.3340, 0.3790] | [0.3338, 0.3790] | 0.3219 | 0.3787 → 0.3790 | 41 / 36 / 38 / 42 |
| 1871331.13__V4 | n4+ | [0.2553, 0.3664] | [0.2455, 0.3690] | [0.4266, 0.4266] | 0.2453 | 0.3664 → 0.3651 | 51 / 25 / 17 / 52 |

**Pooled multi contrasts** (stratum-weighted mean, group-cluster 95% CI; plumbing, not estimates):

| contrast | measure | mean | 95% CI |
|---|---|---|---|
| real − two_tier | `taf_face_mid` | −0.0110 | [−0.0186, −0.0029] |
| real − two_tier | `taf_face_min` | −0.0088 | [−0.0149, −0.0014] |
| real − two_tier | `taf_face_max` | −0.0132 | [−0.0225, −0.0033] |
| real − two_tier | TAF_gb as returned | −0.0139 | [−0.0212, −0.0024] |
| real − shuffled | `taf_face_mid` | −0.0177 | [−0.0306, −0.0041] |
| oracle − two_tier | `taf_face_mid` | −0.0291 | [−0.0339, −0.0243] |
| real − two_tier | opt-face count | −14.6 | [−21.4, −8.9] (weighted two-tier base 52.7) |
| real − two_tier | opt7-face count | −15.3 | [−21.7, −9.9] |
| real − shuffled | opt-face count | +6.1 | [−0.8, 14.4] |
| shuffled − two_tier | opt-face count | −20.7 | [−34.2, −10.1] |

- Per stratum, real − two_tier `taf_face_mid` is −0.0128 (n2), −0.0136 (n3) and +0.0040 (n4+, 2 rows); BH q 0.75 / 0.75 / 1.0.
- Oracle headroom captured: 0.377 on `taf_face_mid` (0.354 on the point TAF).
- On the corrected face, both real and shuffled leave fewer reactions free than two-tier, and shuffled at least as few as real.
- The mechanical `verdict()` on these 10 rows: C3 expressible; H1 "works"; H2 "reduces degeneracy ..., NOT specific to p", the same at 1e-7. With 8 multi amplicons this shows only that the amended rules execute.

**Time and memory.**

- Wall per amplicon, mean 131.6 s (56.3-256.2), of which gapfill 66.0 s and arms 63.4 s. The first smoke spent 24.1 s on the arms of the same amplicons: the four opt7 FVA passes are new, and FVA on the tight faces is slower per LP. The vertex re-solves cost about 0.2 s per arm.
- Peak resident memory (`ru_maxrss`; the harness labels it `gb`, the unit is GiB):
  - pool workers, one task each: 1.462-1.592 GiB;
  - parent (pool mode, including the summary's bootstrap): 0.731-0.735 GiB;
  - process tree (5 workers + parent, polled every 2 s from `/proc/<pid>/status`): 8,819,508 kB = 8.41 GiB and 8,687,064 kB = 8.28 GiB for the two batches.
- Growth over tasks in one process: one process (`--workers 1`, the run in `wobj_v2/growth4/`) ran 4 tasks in a row, 2212470.277 (gramneg), 3125783.3 (archaea), 84109.4 (grampos), 2715752.3 (archaea); its running peak after each was 1.542, 1.688, 1.906 and 2.028 GiB. Part of the rise is the template cache (a new template on each of the first three tasks), but the fourth task, on an already loaded template, still added 0.122 GiB. An earlier pair in one process (`wobj_v2/g`: 3125783.3 then 2715752.3) gave 1.512 then 1.651 GiB, against 1.472 GiB for 2715752.3 in a fresh pool worker. The full run therefore restarts workers every 2 tasks (section 14).

## 14. Full run plan (not launched)

- **Population**: 360 amplicons (n2 = n3 = n4+ = 100, n1 = 60), the harness default. Its first 10 jobs are the smokes' 10.
- **Launch**: `bash scripts/launch_wobj.sh`. It starts the transient user unit `wobj-run` with PYTHONPATH=`/home/freiburger/Documents/ModelSEEDpy-c2e31a0`, `--expect-modelseedpy c2e31a0`, 8 workers, `--maxtasksperchild 2` and MemoryMax 20G. Outputs go to `data/weighted_objective/` (a new directory), the log to `logs/wobj_run.log`. The script refuses to start if MemAvailable < 35 GiB, if `wobj-run` is active, or if the output directory or the log exists.
- **Stages**: (1) `--population-only` builds `population.json` alone (4.17 GiB peak, about 2 min, first smoke); (2) the pool runs every amplicon and writes `summary.json`.
- **Memory**: estimated peak 8 × 1.75 GiB + 0.74 GiB = **14.7 GiB**.
  - 1.75 GiB per worker is the largest single-task worker peak of the second smoke (1.604 GiB, 2703885.11 in the run before the retry; 1.592 GiB in the final run) plus the measured growth from a second task in the same process (+0.146 GiB).
  - 0.74 GiB is the parent's peak (0.731-0.735 GiB).
  - `--maxtasksperchild 2` caps the growth: over 4 tasks one process grew 1.542 → 1.688 → 1.906 → 2.028 GiB (section 13.2), so the harness default of 10 tasks could reach about 3 GiB per worker.
  - MemoryMax = ceil(1.3 × 14.74 GiB) = **20G**. The smoke's tree peak (8.41 GiB with 5 workers) scales to 13.0 GiB at 8, within the estimate.
  - A worse case: a worker whose second task is a large model on a template it has not loaded yet. In the 4-task run, a task on a new template added 0.218 GiB (1.688 → 1.906), so such a worker could reach 1.604 + 0.218 = 1.82 GiB, and 8 of them 15.3 GiB. That still leaves MemoryMax 20G a margin of 1.31×.
- **Failure mode to watch.** If the cgroup's OOM killer kills a pool worker, `multiprocessing.Pool` replaces it, but the task it held is lost and `imap_unordered` waits for it forever. The unit then stays active while `logs/wobj_run.log` gets no new `[wobj] <key>` rows. The longest smoke amplicon took 256 s, so no new row for about 15 min means a stall, not slowness. A stopped run resumes safely (finished keys are skipped): run stage 2 by hand with a new log name.
- **Duration**: 360 × 131.6 s ÷ 8 = 5,922 s = 1.65 h at the second smoke's load. Worker restarts add 180 × at most 20 s ÷ 8 = 0.13 h, and the population stage about 0.04 h. Under the heavier load the first smoke met (147.3 s against 103.8 s, a ratio of 1.42), the amplicon time would be 2.34 h. So about **1.8-2.5 h**. The second smoke ran with 5 concurrent workers; 8 workers and other sessions' jobs on the machine may slow each worker.

## 15. Limitations and open points

- **Truth support is an annotation-level proxy.** It is the gene-associated reactions of the truth genome's own draft, built by the same EC→role route. Flux through a truth-absent reaction means flux through a reaction the source's BV-BRC annotation does not support, not demonstrated biochemistry.
- **The EC route stands in for the production function-string route.** Production counts RAST function strings, while this test counts template roles granted by EC sets. `calibrate_pfunc.py` found EC-level and reaction-level calibration close (reliability 0.021 against 0.020). The `product` proxy for function strings is worse calibrated (0.044).
- **Draft-building fidelity gap.** Production drafts are built with `allow_all_non_grp_reactions=True` (`build_prfba_community_models.py` l.357, positional `build_base_model(rid, "0", True, False)`; `load_genomes_to_modelseedpy.py` l.329). `gapfill_block` uses False. On the two review probe unions this adds 12 and 14 universal or spontaneous reactions to the draft (1,177 → 1,189 and 1,140 → 1,154), of which the gapfilled model already had 1 and 0 (`probe_*.json`). Small, but these models are not byte-for-byte production drafts.
- **One medium, one carbon cap, one floor.** All arms use complete medium, `C <= 100` and 95% biomass. Both the composition and the addressable share depend on the medium.
- **The face interval is of the total truth-absent flux.** `F_lo` and `F_hi` bound the sum over truth-absent reactions, not each reaction. The point metrics are what a user of prFBA receives, but they are secondary, because the vertex is arbitrary among alternative optima (section 12, R1). Their vertex spread is reported beside them.
- **One shuffled permutation per amplicon** (review O8, declined). Its per-amplicon noise enters the bootstrap CI of real − shuffled; more permutations would narrow that CI, not move its centre.
- **Target population.** Amplicons whose md5 is scored for one source genome only (section 5). Strata are drawn rarest first under one global one-per-source-genome rule.
- **The ex_weight-1 arms are exploratory**: no claim about the released objective rests on them.
- **The mscommunity checkout moves.** `/home/freiburger/Documents/MSCommunity` is another session's checkout; it moved from 375572d to 176b6e1 during the first session. The released prFBA imports it, but never calls it on the single-model path. Its commit is recorded in `meta.libraries.mscommunity`.
