# PR: MSProbability.prFBA, make the single-model path run and make its two policies explicit

- Repository: freiburgermsu/ModelSEEDpy
- Base: `atp-safe-and-drain-fixes` (c2e31a0, the same base as the reversal fix, PR #19)
- Head: `fix-prfba-single-model`, commits b6857a0 and 09856f5. Two commits on c2e31a0; it does not depend on `fix-gapfill-reversals`.
- Worktree: `/home/freiburger/Documents/ModelSEEDpy-prfba`
- Files: `modelseedpy/core/msprobability.py` (+95 / -10) and the new `tests/core/test_msprobability_prfba.py` (277 lines, 12 tests)
- Investigation behind it (reproductions on both trees, root-cause commits, callers): `prFBA/region_validation/shared_reference/data/weighted_objective/msprob_investigation/INVESTIGATION.md`

## Submit (not yet pushed)

```bash
git -C /home/freiburger/Documents/ModelSEEDpy-prfba push origin fix-prfba-single-model
cd /home/freiburger/Documents/ModelSEEDpy-prfba
gh pr create --repo freiburgermsu/ModelSEEDpy --base atp-safe-and-drain-fixes \
  --head fix-prfba-single-model \
  --title "MSProbability.prFBA: make the single-model path run, and make its two policies explicit" \
  --body-file /home/freiburger/Documents/prFBA/region_validation/shared_reference/PR_ModelSEEDpy_prfba_single_model.md
```

(`--body-file` posts this whole file. To post only the description, cut it at the "PR body" heading below.)

---

## PR body

### Summary

On c2e31a0, `MSProbability.prFBA` fails on every call with one model, and on every call with several. The single-model path has not run in any version since the first draft (f351cb4, 2024-02-13), and no test exercised `prFBA` before this PR. This PR makes the smallest changes that let it run, refuses a model that cannot grow instead of returning an all-zero "optimal" answer, and turns two silent assumptions into documented keywords that default to the current behavior.

### What was broken, and the fix

| # | Released behavior (c2e31a0) | Fix |
|---|---|---|
| 1 | `MSModelUtil(model_s_ if len(model_s_) == 1 else ...)`: a bare `cobra.Model` raises `TypeError: object of type 'Model' has no len()`. A one-item list reaches `MSModelUtil([model])`, which raises `AttributeError` (a list has no `.id`). | A bare model is wrapped, and a one-item list is unwrapped, so `MSModelUtil` gets the model itself. An empty list raises `ValueError`. |
| 2 | `mdlUtil.add_medium(environment)` is called even when `environment=None`, the default, and raises (`AttributeError` in `FBAHelper.convert_kbase_media`). Root cause: prFBA's own guard was dropped in the merge 937829d (2024-08-26), harmless until 5149cc3 (2024-09-13) removed `add_medium`'s `if media is None` guard. | `add_medium` is skipped when `environment is None`, so the model keeps its own exchange bounds. `MSModelUtil.add_medium(None)` itself still raises; every other caller guards it, so this PR fixes only prFBA's call site. |
| 3 | `ElementUptakePkg(...).build_package({"C": 100})` runs twice. The second call raises `PackageError: Package with name ElementUptakePkg already in model!`. The duplicate is outside the `if not commkinetics and len(model_s_) > 1` block, so the community path raises it too (shown on a two-member fixture below). | The package is built once, at its first position. |
| 4 | `print([cons.name for cons in model.constraints])` runs on every call. | The print is moved under `printLP`, next to the LP file it goes with. |
| 5 | When the medium admits no growth, the 95% floor is 0 and the minimization returns the all-zero flux vector with status `optimal`, silently. With fix 2 this is what `environment=None` returns on a gapfilled model whose own bounds admit no growth. An infeasible model reaches the solver as a NaN bound (an opaque `GurobiError`). | New keyword `min_growth` (default 1e-6): below it, or on NaN, raise `ValueError` naming the objective. `min_growth=None` restores the unchecked behavior. (09856f5) |

No bound, constraint or objective coefficient changes. The coefficient is still `max(min_prob, 1 - p**prob_exp)` on both directions of every `rxn` reaction and `ex_weight` on both directions of every `EX_` exchange. The floor is still 95% of the entry objective's maximum under `C <= 100`. The comment above the coefficient loop used to say `(1-p^prob_exp)*flux + min_prob`, which is not what the code computes; it now matches the code.

### Two policy decisions, kept separate

**(a) Reactions without `notes["probability"]`. Decision: keep raising by default, and add an opt-in keyword.**
`missing_probability=None` (the default) raises `KeyError`, as before. Two things change: the check runs before the model is modified, and the message names the reactions (count plus the first five). `missing_probability=0.0` treats them as having no genetic evidence, giving coefficient `max(min_prob, 1 - 0) = 1`.
Why it matters: a caller that writes the probabilities on the draft and then gapfills ends up with gapfilled reactions that carry no note. Treating them as 0 is the rule the prFBA pipeline applies (p(reaction) = max over its genes, and 0.0 when it has no genes). That rule belongs to the caller, so the library does not apply it silently.

**(b) The objective on entry. Decision: keep the current behavior by default, document the hazard, and add an opt-in keyword.**
`prFBA` floors whatever objective the model carries. `MSATPCorrection.build_tests` sets the objective to ATP hydrolysis (`msatpcorrection.py` l.606), so a model built through ATP correction and gapfilling can arrive with `ATPM_c0` as its objective. In that case the 95% floor holds ATP maintenance, not growth. All three gapfilled models below arrived with objective `1.0*ATPM_c0 - 1.0*ATPM_c0_reverse_...`. On one earlier probe model, the released call (once it was made to run) returned bio1 = 0.0000 and ATPM 314.69 (`WEIGHTED_OBJECTIVE_TEST.md` section 2). `biomass_id=None` (the default) keeps the current behavior. `biomass_id="bio1"` sets the objective before the floor is computed.

The docstring is new. It states the objective, that the model is modified in place (pass `model.copy()`; a second call on the same object raises `PackageError`), how each `environment` form is handled (I checked this in `MSModelUtil.add_medium` and `FBAHelper.convert_kbase_media`: an `MSMedia` goes through `convert_kbase_media`, whose default `uniform_uptake=1000` sets every listed compound's uptake to 1000), and the two keywords.

### Tests

`tests/core/test_msprobability_prfba.py` uses `e_coli_core.json`, with compartment suffixes added as in `test_msgapfill.py`. Every internal reaction is renamed `rxnNNNNN_c0` and given a probability from {0, 0.25, 0.5, 0.75, 1}. Biomass and ATPM keep their ids, so their coefficient is 0. The tests check that:

1. `prFBA(model)` and `prFBA([model])` both run, and their objectives agree to a relative 1e-9.
2. `environment=None` matches a dict of the model's own medium; a lower glucose uptake lowers the floor.
3. Biomass is held at 95% of its maximum under `C <= 100`: the floor's `lb`, the floor's variables (the biomass reaction's), and biomass flux at the optimum.
4. Every objective coefficient is `max(min_prob, 1 - p**prob_exp)` for `rxn`, `ex_weight` for `EX_` and 0 otherwise, with `min_prob=0.05, prob_exp=2, ex_weight=7`. The solver's objective equals the priced flux.
5. The call prints nothing and writes no `prFBA.lp` unless `printLP=True`.
6. A missing note raises `KeyError` naming the reaction, and the model is left unchanged (no `min_value`, no `C_elements`, same objective).
7. `missing_probability=0.0` matches annotating the reaction with `"0"`.
8. With `ATPM_c0` as the entry objective, the floor is on ATPM. `biomass_id` moves it back to biomass and reproduces the biomass run.
9. The two-member community fixture (from MSCommunity's `tests/test_f_biomass_drains.py`) runs to `optimal`. This checks only that the `PackageError` is gone.
10. A model that cannot grow (`environment={}` with ATP maintenance released) raises `ValueError` instead of returning all-zero fluxes.
11. An infeasible model (`environment={}` with ATP maintenance at its lower bound) raises `ValueError` naming infeasibility, not a solver error.
12. `min_growth=None` keeps the unchecked behavior (status optimal, every flux 0).

Results (`~/Documents/py_venv/bin/python -m pytest -p no:cacheprovider`, `PYTHONPATH` set to each tree). The new file ran from a scratch directory for both trees. For `tests/core`, c2e31a0 ran from a scratch copy of its `tests/`, with a `../ModelSEEDDatabase` symlink; the branch ran from its own worktree:

| Tree | New test file | `tests/core` |
|---|---|---|
| c2e31a0 | **9 failed** (`TypeError: object of type 'Model' has no len()`, `PackageError ... already in model!`, `unexpected keyword argument 'missing_probability'`); tests 10-12 postdate that run | 1 failed, 69 passed |
| this branch (09856f5) | **12 passed** (gurobi; also 12 passed with `cobra.Configuration().solver = "glpk"`) | 1 failed, 81 passed |

The one `tests/core` failure on both trees is the known, unrelated `tests/core/test_msgapfill.py::test_ms_gap_fill2`. The 12 extra passes on the branch are the new file. To run `test_msatpcorreption.py` from outside a checkout, `../ModelSEEDDatabase` must resolve: it is a cwd-relative path in `modelseedpy/biochem/modelseed_biochem.py`.

### Equivalence with the faithful copy used in the prFBA weighted-objective test

The prFBA weighted-objective test (`prFBA/region_validation/shared_reference/WEIGHTED_OBJECTIVE_TEST.md`, harness `scripts/weighted_objective_test.py`) could not call the released function. It used `prfba_released_copy()` instead, a literal copy with the harness's edits (1)-(3): the unwrap, dropping the duplicate package, and dropping the print, `printLP` and the community branch.

For the check, three amplicons were rebuilt in one process with the harness's own functions, one per multi-genome stratum (the fastest to gapfill in each). The functions were `gapfill_block` (as_shipped, complete medium, reversals applied), `role_probabilities`, the max-over-genes rule and `set_probs`, on templates `ModelSEEDTemplates-c12run` v7.0, with `PYTHONPATH` set to this branch. On each, the real arm was solved four ways:

- (a) the harness copy, with `bio1` set as the objective by the caller, as the harness does;
- (b) the fixed function on the bare model as built (objective ATPM_c0), with `biomass_id="bio1"`;
- (c) the fixed function on `[model]`, with the notes stripped from every gapfilled reaction and `missing_probability=0.0`;
- (d) the fixed function on (c)'s model with the default `missing_probability`.

| Amplicon | Stratum | Reactions | Gapfilled | Objective (a) copy | (b) fixed | rel. diff (b) vs (a) | (c) list + missing 0 | rel. diff (c) vs (a) | bio1 (a) = (b) | (d) default | Stored real-arm objective in `per_amplicon.jsonl` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 298656.3__V7-V9 | n2 | 1011 | 97 | 17132.478474090174 | 17132.478474090174 | 0.0 | 17132.478474090174 | 0.0 | 1.3648274532 | KeyError | 17132.478474090174 |
| 47303.8__V4-V5 | n3 | 974 | 84 | 20223.754412977967 | 20223.754412977967 | 0.0 | 20223.754412977967 | 0.0 | 1.1284120847 | KeyError | 20223.754412977967 |
| 3101495.5__V4 | n4+ | 856 | 105 | 19818.088318215458 | 19818.088318215458 | 0.0 | 19818.088318215458 | 0.0 | 1.2565526157 | KeyError | 19818.088318215458 |

All three agree exactly (relative difference 0.0, well under the 1e-9 target). (b) printed nothing to stdout (0 characters). The rebuilt models also matched the stored run: the same number of gapfilled reactions (97, 84, 105) and the same real-arm objective. That three-way match is informational, since gapfilling is rebuilt rather than loaded. The deliverable is the in-process agreement between the copy and the fixed function. Script and output: `/tmp/claude-1000/-home-freiburger-Documents-prFBA/0c6f5bc6-8224-5a34-9b20-7d5643487901/scratchpad/prfba_fix/equiv/` (`equiv.py`, `equiv.json`, `equiv.log`; scratch, not durable).

### Not changed, and worth knowing

- **Community path, kinetics: a prFBA call-site regression, not fixed here.** With two models, the fixed function runs (status optimal, objective 763.8 on the two-member fixture). But the community kinetic rows never reach the solved model:
  - prFBA constructs `MSCommunity(mdlUtil.model, ...)`, discards it, and solves `mdlUtil.model`. MSCommunity has always worked on `MSModelUtil(model, copy=True)`, so its `*_commKin` rows land on its copy: the solved model has 19 constraints and no kinetic rows, the copy 20, including `memA_commKin` and `memB_commKin` (MSCommunity 176b6e1; reproduced at 03bd80c).
  - Until fdaf0a2 (2024-03-13) prFBA built `CommKineticPkg` through the solved model's package manager; that commit commented the two lines out. The comment added in 0050e06, "the MSCommunity constructor installs the community kinetic rows", is true only of the copy.
  - Consequences: `kinetics_coef` has no effect (objective 763.8 at 1,000 and at 1); `commkinetics=True` builds no kinetics anywhere, because `build_from_species_models` never reads the flag; `abundances` are not forwarded to MSCommunity, whose copy is renormalized to 1/n.
  - A follow-up should solve the model MSCommunity constrains, or rebuild the rows on the solved model; MSCommunity's copy also adds `rxn00062_c1/_c2` without probability notes, so that change needs `missing_probability`. The new test asserts nothing about kinetics.
- **Downstream effect on the prFBA harness.** `weighted_objective_test.py`'s `released_prfba()` exec-edits the c2e31a0 source. It asserts that `mdlUtil = MSModelUtil(model_s_ if len(model_s_) == 1 else` occurs exactly once, and `worker_init` calls it unconditionally. On any ModelSEEDpy with this fix, the harness therefore stops at startup with `RuntimeError: released prFBA source changed`, even with `--no-check-released`. The equivalence check above stubbed it out for that reason. The remedy belongs in prFBA, not here: call `MSProbability.prFBA(..., biomass_id="bio1")` directly as the released check, or skip the exec-edit when the library is not c2e31a0.
- The new tests `importorskip` `mscommunity` and `cobrakbase`, because `msprobability.py` imports both at module level. On a CI without them, all 9 tests skip.
- `abundances`, `commkinetics`, `kinetics_coef` and `expression` are unchanged. `expression` is unused, and `commkinetics` and `kinetics_coef` are functionally unused too (previous bullet).
- Only reactions whose id starts with `rxn` are priced. On a model with non-ModelSEED ids (BiGG, e.g. e_coli_core) no internal reaction is priced and the probabilities have no effect; the tests rename e_coli_core's reactions for this reason.
- Each call leaves the model in two class-level caches, and `msprobability.py` imports `mscommunity` at module level, so the in-function import is redundant.
- For reviewers: the open cshenry PR #48 would carry failures 1-4 into cshenry/main, which today lacks 2 and 3. freiburgermsu `main` fails earlier still, at the deleted `modelseedpy.community.commhelper` import.
- The function still modifies its argument in place and registers `ElementUptakePkg` on it, so it cannot be called twice on the same object. The docstring says so.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
