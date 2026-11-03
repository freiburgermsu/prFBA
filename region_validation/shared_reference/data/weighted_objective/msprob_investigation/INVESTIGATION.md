# `MSProbability.prFBA`: failure points, root causes, history and impact

2026-11-02. Read-only investigation. No ModelSEEDpy tree was checked out, reset or modified, and nothing was committed.

**Code examined**

- Released source: ModelSEEDpy c2e31a0, `modelseedpy/core/msprobability.py`, `prFBA` at l.154-201 (l.154 is `@staticmethod`). The file is byte-identical (md5 `8d40e57a…`) at e1254ba (`ModelSEEDpy-reversals`), in the shared clone (`ModelSEEDpy`, also at c2e31a0) and in EmilyKin's two pinned copies (`community_modeling/_pinned/ModelSEEDpy-{c2e31a0,d884a50}`).
- Fix: `fix-prfba-single-model`, b6857a0 (worktree `ModelSEEDpy-prfba`, md5 `1014bf2e…`). On b6857a0, `prFBA` is at l.154-268.
- Environment: MSCommunity 03bd80c (editable install, `~/Documents/MSCommunity`), cobra 0.29.1, optlang 1.9.1, Gurobi.

Line numbers below are c2e31a0's unless marked "b68".

## Files in this folder

| File | What |
|---|---|
| `repro.py` | Self-contained reproduction. Run it once per tree (commands in its docstring). It takes about 20 s and peaks at 0.9 GB RSS. |
| `repro_c2e31a0.log` | Run with `PYTHONPATH=…/ModelSEEDpy-c2e31a0`. The header proves the tree: `modelseedpy.__file__`, the source file of `prFBA`, its md5, and HEAD c2e31a0. |
| `repro_b6857a0.log` | Run with `PYTHONPATH=…/ModelSEEDpy-prfba` (HEAD b6857a0). |
| `fixture_g100053.8_gf.json.gz` (+ `.meta.json`) | The gapfilled probe model of WEIGHTED_OBJECTIVE_TEST.md section 2, copied from scratch `wobj/step1/`: genome 100053.8, gramneg template, 1,261 reactions, 61 of them gapfilled. Its objective as built is `ATPM_c0`. |
| `fixture_kbase_complete_media.txt` | The 347 compound ids of `ModelSEEDDatabase/Media/KBaseMedia.cpd`. `ENV = {cpd: 100.0}`. |
| `aux_branch_probe.py` (+ `aux_branch_probe_{cshenry_main,origin_main}.log`) | Runs the same calls on two other branches, extracted read-only with `git archive` into scratch. |

`repro.py` has four parts:

- **Part A** runs the released call forms on both trees.
- **Part B** runs on c2e31a0 only. It reaches defects that an earlier crash hides, using *line-count-preserving* edits of the released module source, compiled under the pseudo file name `msprobability.py[c2e31a0+F1(+F2)]`. F1 unwraps one model on l.165; F2 removes the duplicate build on l.173-174. Every traceback therefore cites the released file's own line numbers. The edits exist only to get past earlier defects.
- **Part C** exercises b6857a0's keywords.
- **Part D** covers the further problems on both trees. On c2e31a0 it runs through the F1+F2 variant; every number agrees with b6857a0.

Case ids (`FP…`, `C…`, `N…`) refer to the logs.

## Summary

| # | Failure point | c2e31a0 | b6857a0 | Introduced by |
|---|---|---|---|---|
| FP1 | One model, bare or as `[model]` | **Confirmed**: `TypeError` at l.165; `AttributeError` at msmodelutl.py:375 | **Fixed** (FP1a/FP1b: 18,799.709290, bio1 1.404069, nothing printed) | f351cb4 (2024-02-13), the first draft; rewritten into the released ternary at d91d096 (2024-02-21). The single-model path has not run in any version examined. |
| FP2 | `environment=None` (default) | **Confirmed**: `AttributeError` in `convert_kbase_media` (fbahelper.py:292) via l.166 | **No longer raises, but returns a silent all-zero "optimal" solution** on the fixture (FP1c, C6). `MSModelUtil.add_medium(None)` itself still raises (FP2b). | Two steps. 937829d (2024-08-26) dropped prFBA's own guard, which was harmless then. 5149cc3 (2024-09-13) removed `add_medium`'s `if media is None` guard, and from that commit the call crashes. |
| FP3 | `ElementUptakePkg` built twice | **Confirmed** on the single path (FP3a) and the community path (FP3b): `PackageError` at l.174 | **Fixed** (FP1a, FP3b) | 937829d (2024-08-26), a merge-conflict splice of three lines from `community_simulation` |
| FP4 | Constraint names printed on every call | **Confirmed**: 19,348 characters per call on the fixture (FP4) | **Fixed**: moved under `printLP` (C5) | 04e940d (2024-03-10) |
| FP5 | Gapfilled reactions carry no `notes['probability']` (hazard) | **Confirmed**: `KeyError` at l.187, after the caller's model has already been changed (FP5) | Default still raises, but before any change and naming the reactions (C1). `missing_probability=0.0` runs (C2). | Since the first draft f351cb4; the current lines date from d91d096 |
| FP6 | The 95% floor goes on whatever objective the model carries (hazard) | **Confirmed**: floor on ATPM_c0, bio1 = 0, ATPM 314.6875 (FP6) | Default unchanged (C4). `biomass_id='bio1'` moves the floor to biomass (C3). | `slim_optimize` floor since f351cb4. The ATPM objective comes from `MSATPCorrection.build_tests` (56021362, Christopher Henry, 2022-07-12). |

New problems, detailed in section 7:

- **N1**: a silent zero solution whenever the medium admits no growth. This is the fix's new default behavior on a gapfilled model as built.
- **N2**: an infeasible model crashes with an opaque Gurobi error.
- **N3**: the comment's objective, read literally, makes fully supported reactions free (757,864.6 units of loop flux on the optimal face, against 588.9).
- **N4**: unpriced reaction classes, and non-ModelSEED ids leave every internal reaction unpriced.
- **N5**: the exchange term is 99.06% of the objective, and half of it prices secretion.
- **N6**: on b6857a0 the community kinetic rows land on MSCommunity's private copy, so `kinetics_coef`, `commkinetics` and `abundances` do not reach the solved model.
- **N7**: in-place mutation, process-lifetime caches, and a `PackageError` on a second call.
- **N8**: probability values are not validated.
- **N9**: hard module-level imports.
- **N10**: from 2024-09-05 to 2026-09-01 every call raised at the function's first line. It still does on freiburgermsu `main`.
- **N11**: two different definitions of "exchange".
- **N12**: no test exercised `prFBA` before b6857a0.

**On the refactor hypothesis.** The refactor that moved community modeling into the standalone `mscommunity` package (0050e06, 2026-09-01) introduced none of FP1-FP6. Since fc33219 (2024-09-05), every call had failed at the function's first statement, an import of the deleted `modelseedpy.community.commhelper` (N10). 0050e06 repointed that import and so *exposed* FP1-FP4, which date from f351cb4 (FP1, 2024-02-13), 04e940d (FP4, 2024-03-10), the merge 937829d (FP2's call-site guard and FP3, 2024-08-26) and 5149cc3 (FP2's `add_medium` guard, 2024-09-13). 0050e06 is responsible for one false comment, "the MSCommunity constructor installs the community kinetic rows" (l.171; see N6), and for the "lazy import" claim in its commit message (N9).

No code in prFBA, codiffusion_bioreactor or EmilyKin calls the released function (section 9). The weighted-objective results used a private copy. That copy was cross-checked on 360 of 360 rows against the released source with the two failing lines corrected.

---

## 1. FP1: one model, bare or as a one-element list

**Symptom.** `prFBA(model, …)` raises `TypeError: object of type 'Model' has no len()` at l.165. `prFBA([model], …)` raises `AttributeError: 'list' object has no attribute 'id'` at msmodelutl.py:375 (`self.id = model.id`). `environment=None` makes no difference (FP1c).

**Repro.** FP1a, FP1b and FP1c in `repro_c2e31a0.log`:

```
<tree>/modelseedpy/core/msprobability.py:165 in prFBA
    mdlUtil = MSModelUtil(model_s_ if len(model_s_) == 1 else build_from_species_models(...))
<tree>/modelseedpy/core/msmodelutl.py:375 in __init__        (list form only)
    self.id = model.id
```

**Root cause.** l.165 tests `len(model_s_) == 1` and then passes `model_s_` itself, not `model_s_[0]`. Neither `cobra.Model` nor cobrakbase's `FBAModel` defines `__len__` (checked), so a bare model fails at `len()`. A one-element list passes the length test, and then `MSModelUtil` receives the list.

**Introducing commit.** The first draft already had the defect: f351cb4, 2024-02-13, "Andrew's WorkMacbook <freiburgermsu@aol.com>", "prFBA code draft". It read `if len(model_s_) > 1: model_s_ = build_from_species_models(model_s_)` followed by `mdlUtil = MSModelUtil(model_s_)`. `MSModelUtil.__init__` already did `self.id = model.id` there (msmodelutl.py:109). d91d096 (2024-02-21, same author, "Commhelper and MSProbability updates") rewrote it into the released ternary, and `git blame` attributes l.165 to 04e940d (2024-03-10), which only appended `commkinetics=`. **The single-model path has not run in any version examined.** That covers the 18 historical versions of the file extracted with `git show <commit>:…` (f351cb4 through c2e31a0, with `self.id = model.id` checked at f351cb4, d91d096, 937829d, 0050e06 and c2e31a0) and every branch in section 8. cshenry/main fails differently; see `aux_branch_probe_cshenry_main.log`.

**b6857a0.** Fixed at b68 l.215-221: a bare model is wrapped, a one-element list or tuple is unwrapped, and an empty list raises `ValueError` (C8). FP1a and FP1b give objective 18,799.709290 and bio1 1.404069. These are the values the patched probe (`patches/test_patched.log`) and the harness copy gave.

**Impact.** Nothing could call the single-model path, so no result came from it. The weighted-objective harness wrote its own copy for this reason (section 9).

## 2. FP2: `environment=None`, the default

**Symptom.** With FP1 bypassed, `add_medium(None)` raises `AttributeError: 'NoneType' object has no attribute 'get_media_constraints'` (FP2):

```
msprobability.py[c2e31a0+F1]:166 in prFBA            mdlUtil.add_medium(environment)
<tree>/modelseedpy/core/msmodelutl.py:2284 in add_medium
    if not hasattr(media, "items"):  media = FBAHelper.convert_kbase_media(media)
<tree>/modelseedpy/core/fbahelper.py:292 in convert_kbase_media
    return {"EX_"+exID: uniform_uptake for exID in kbase_media.get_media_constraints().keys()}
```

FP2b shows the same traceback from `MSModelUtil(model).add_medium(None)` alone, on **both** trees.

**Root cause.** Two guards were removed, eight months apart:

1. **prFBA's guard.** Every version before the merge had `if environment is not None:  mdlUtil.add_medium(environment)` (f351cb4, d91d096, 04e940d…35d3a40). It was removed in the merge-conflict resolution 937829d (2024-08-26 10:57:45 -0500, "Freiburgermsu <freiburgermsu@aol.com>", "merge conflict resolution", parents f7ca721 and 35d3a40). In `git show 937829d` the line `mdlUtil.add_medium(environment)` appears as `++` (in neither parent) in place of the guarded line. It is copied from f7ca721's `community_simulation`, which called `add_medium(environment)` unguarded.
2. **add_medium's own guard.** At 937829d this was still harmless, because `MSModelUtil.add_medium` began with `if media is None:  return self.model.medium`. That guard was added at 85156f0 (2023-12-11, "Andrew Freiburger", "error catching of MSUtils"). It was **removed at 5149cc3** (2024-09-13 18:55:02 +0000, "Freiburgermsu_Sequoia <freiburgermsu@aol.com>", "CommScores associated edits"), the same diff that added the cpd-to-`EX_…_e0` mapping (`git log -L '/def add_medium/,+10:modelseedpy/core/msmodelutl.py'`).

**The crash has been live since 5149cc3.** It was masked by N10's import failure until 0050e06.

**b6857a0.** prFBA no longer calls `add_medium` when `environment is None` (b68 l.232-233). The guard matches every other caller in the package (msminimalmedia.py:108, 211, 328, 344 all test `if environment:`) and `MSModelUtil.__init__` (l.374). `add_medium(None)` itself is not fixed (FP2b), so its root cause remains in `MSModelUtil`.

The fix exposes a worse outcome. "Keep the model's own bounds" on the gapfilled fixture means bounds under which bio1 cannot grow. They are the bounds the last ATP-correction or gapfilling medium left: 332 uptakes are open, but not O2, NH3, SO4, K+ or Ca2+, and bio1 max is 0.0. The floor becomes 0.95 × 0 = 0, and the call returns `status=optimal objective=0 bio1=0` with every flux zero and no warning (FP1c and C6 in `repro_b6857a0.log`). This is problem N1.

**Impact.** No caller used the default. The harness always passes `environment` (weighted_objective_test.py l.388-392).

## 3. FP3: `ElementUptakePkg` built twice

**Symptom.** With FP1 bypassed, every single-model call raises `PackageError: Package with name ElementUptakePkg already in model!` (FP3a). The community path raises it unmodified (FP3b), after `build_from_species_models` and the `MSCommunity` constructor (l.172) have already run:

```
<tree>/modelseedpy/core/msprobability.py:174 in prFBA
    elepkg = ElementUptakePkg(mdlUtil.model)  ;  elepkg.build_package({"C": 100})
<tree>/modelseedpy/fbapkg/basefbapkg.py:46 in __init__         self.pkgmgr.addpkgobj(self)
<tree>/modelseedpy/fbapkg/mspackagemanager.py:57 in addpkgobj
    raise PackageError("Package with name " + classname + " already in model!")
```

**Root cause.** l.168 registers an `ElementUptakePkg` on `mdlUtil.model`. l.173-174 construct a second one on the same model, outside the `if not commkinetics and len(model_s_) > 1` block. `MSPackageManager.get_pkg_mgr` is keyed by the model object, and `addpkgobj` refuses a second package of the same class. That rule has existed since e261105 (2021-07-12, "Christopher Henry", "Adding package manager to avoid overlapping packages"), so the package manager did not change.

**Introducing commit.** 937829d, the same merge as FP2. In `git show 937829d -- modelseedpy/core/msprobability.py`, these three lines are `++` lines, present in neither parent:

```
++        # constrain carbon consumption
++        elepkg = ElementUptakePkg(mdlUtil.model)  ;  elepkg.build_package({"C": 100})
++        # evaluate the metabolomics data over time
```

They are f7ca721's `community_simulation` body (`# constrain carbon consumption`, the same `ElementUptakePkg` line, `#  the metabolomics data over time`), spliced into prFBA by the conflict resolution. Every earlier version built the package once.

**b6857a0.** Fixed: built once, at b68 l.235 (FP1a, FP3b).

**Impact.** From 937829d onward the community path could not run either. This is moot in practice, because of N10.

## 4. FP4: constraint names printed on every call

**Symptom.** l.194 `print([cons.name for cons in mdlUtil.model.constraints])` writes 19,348 characters per call on the fixture (FP4, FP6, N1, N3…N8 in `repro_c2e31a0.log`; 77,392 characters for N7's four calls).

**Introducing commit.** 04e940d (2024-03-10 00:46:43 -0600, "Andrew's WorkMacbook <freiburgermsu@aol.com>", "MSProbability and associated edits"). The first draft f351cb4 printed the whole coefficient dict (`print(coef)`), and d91d096 removed that print.

**b6857a0.** Fixed: printed only under `printLP=True` (b68 l.262-263). C5 shows 0 characters and no `prFBA.lp` by default; with `printLP=True` it prints 19,348 characters and writes `prFBA.lp`. `printLP` still writes `./prFBA.lp` into the current directory, overwriting any existing file.

**Impact.** Cosmetic, but the harness had to drop this line (edit 3 in `prfba_released_copy`) and redirect stdout around `released_prfba()` (l.846).

## 5. FP5 (hazard): reactions without `notes['probability']`

**Symptom.** If the gapfilled reactions carry no note (the production order: notes are written on the draft, then the model is gapfilled), l.187 raises `KeyError: 'probability'` (FP5). By that point the caller's model has already been changed: `C_elements` and `min_value` have been added and the medium replaced (332 → 347 entries). The objective is not yet replaced.

**Root cause.** `float(rxn.notes["probability"])` with no default. Reactions added after annotation have no note: gapfilled reactions, and also the `rxn00062_c1/_c2` that MSCommunity adds to its copy (N6).

**Introducing commit.** This is a design gap present since the first draft f351cb4 (`rxn.notes["probability"]` in a comprehension over `internal_list()`). The current two lines, both directions over `model.reactions`, are d91d096 (2024-02-21).

**b6857a0.** Explicit. `missing_probability=None` raises `KeyError` naming the count and the first five ids, before the model is touched (b68 l.222-227; C1 shows the model unchanged). `missing_probability=0.0` gives coefficient 1, identical to annotating `"0.0"` (C2: relative difference 0).

**Impact.** The harness annotates every `rxn` reaction after gapfilling (l.784-786 production rule, `set_probs` l.457-460), gapfilled ones with 0.0. The manuscript states this (prFBA.tex l.215).

## 6. FP6 (hazard): the floor goes on whatever objective the model carries

**Symptom.** On the fixture as built (objective `1.0*ATPM_c0 - 1.0*ATPM_c0_reverse_76e6c`), the 95% floor goes on ATP maintenance (`min_value` holds ATPM_c0, lb 314.6875). The result is `objective=31762.734375 bio1=0.000000 ATPM_c0=314.687500` (FP6 on c2e31a0, C4 on b6857a0).

**Root cause.** l.178 `maxBioSol = mdlUtil.model.slim_optimize()` maximizes the entry objective, and l.179 `add_minimal_objective_cons(maxBioSol*.95)` floors `self.model.objective.expression` (msmodelutl.py:1039-1042). `MSATPCorrection.build_tests` sets `self.model.objective = self.atp_hydrolysis.id` (msatpcorrection.py:606; blame 56021362, Christopher Henry, 2022-07-12), so a model built through ATP correction and gapfilling arrives with ATPM as its objective.

The floor is computed after the carbon cap (l.168), so it is 95% of the *carbon-capped* maximum. Nothing checks the direction: an entry objective in `min` direction would be "floored" at 95% of its minimum.

**Introducing commit.** The `slim_optimize` floor dates from f351cb4. The `add_minimal_objective_cons(maxBioSol*.95)` idiom predates prFBA (77f1ce4, 2024-01-10, in `community_simulation`).

**b6857a0.** Default unchanged (C4). `biomass_id='bio1'` sets the objective before the floor (b68 l.229-230): floor on bio1 at lb 1.404069, objective 18,799.709290 (C3). That is the value the harness gets by setting bio1 itself.

**Impact.** The harness sets `M.objective = "bio1"` before every call (l.829).

---

## 7. Further problems in the same function

### N1. A silent zero solution when the medium admits no growth (new default behavior on b6857a0)

When the medium admits no biomass, `maxBioSol` is 0, the floor is 0, and the minimization returns the all-zero vector with `status=optimal`:

- N1 (both trees): `environment={}` closes every uptake and gives objective 0 and every flux zero.
- FP1c and C6 (b6857a0): the default `environment=None` on the fixture as built gives the same.

Nothing warns. b6857a0 makes `environment=None` "work", so a user who omits the medium on a gapfilled model now gets a plausible-looking empty answer instead of a crash.

A related inconsistency: msminimalmedia skips `add_medium` for any falsy environment (`if environment:`), but b6857a0 tests `is not None`, so an empty dict closes all uptake.

Suggested guard: raise, or at least warn, when `maxBioSol` is NaN or ≤ a tolerance.

### N2. An infeasible model fails with an opaque solver error

On e_coli_core with `environment={}` (ATPM lb 8.39 with no uptake), `slim_optimize()` returns NaN. `add_minimal_objective_cons(nan)` then makes optlang hand Gurobi a NaN bound: `GurobiError: Invalid argument to Model.addLConstr` at msmodelutl.py:1042, reached via l.179 (b68 l.246). This happens identically on both trees. Other solvers may behave differently. Either way the user gets no "infeasible medium" message.

### N3. The comment's objective against the code's: which is intended, and what difference it makes

- **Which is intended.** The comment at l.182 reads `sum_r ((1-p^prob_exp)*flux_r + min_prob) + …`, but l.187-188 compute `max(min_prob, 1-p^prob_exp)`. The two have disagreed since the first draft f351cb4. The manuscript's eq:objective (prFBA.tex l.276-283) uses the `max` form and says the floor "keeps a minimum parsimony penalty on fully supported reactions so that none becomes cost-free" (l.286). The code is the intended form; the comment is wrong. b6857a0 corrects the comment (b68 l.249).
- **Literal reading.** `+ min_prob` outside the product is a constant (|R|·min_prob) with no effect on the argmin, so the objective reduces to `(1-p^e)·|v|` and every reaction with p = 1 is free. On the fixture (N3, same constraints), the solver returned the same vertex. But the optimal face then admits **757,864.6** units of flux through the 810 p = 1 reactions (zero-cost loops), against **588.9** under the code's coefficients. With graded p the figures are 211,891.4 against 199.6. The floor exists to remove exactly this degeneracy.
- **Additive reading.** The coefficient `(1-p^e) + min_prob` differs from `max` only by +min_prob on reactions with p ≤ 1 - min_prob. On the fixture it returns the same flux vector (L1 0.0000) in both the two-tier and graded settings. Only the objective value shifts, by 0.01 × Σ|v| over p < 1 reactions: 18,801.414475 against 18,799.709290. So the additive form would change little.

The `max` form has one side effect: all p ≥ 1 - min_prob (p ≥ 0.99 at the default) tie at min_prob. N8 shows p = 0.999 receiving the same coefficient as p = 1.

### N4. Reaction classes left at zero cost

Only ids starting with `"rxn"` (l.186) and `"EX_"` (l.189) are priced.

**On the gapfilled fixture,** 5 of 1,261 reactions are unpriced (N4): `bio1`, `ATPM_c0`, `SK_cpd02701_c0`, `SK_cpd03091_c0` and `SK_cpd11416_c0`. There are no `DM_` reactions. All five are one-way, with bounds (0, 1000), so they cannot form a zero-cost cycle.

- **ATPM_c0** is free ATP dissipation. Its lb is 0, so no maintenance demand is imposed.
- **The two non-biomass sinks** are free outlets for `cpd02701` and `cpd03091`. The carbon cap does not count them, because ElementUptakePkg counts only `'EX_' in id` (msmodelutl.py:680-681; N4: `C_elements` counts no `SK_` variable).
- **Measured: no free route is exploited on this model.** Both sinks carry 0 flux. Pricing them like an exchange leaves the objective unchanged (18,799.709290 → 18,799.709290), and the same holds under uptake-only pricing (N5).
- **The community merge adds more unpriced reactions.** N6 lists `bio1`, `bio2`, `bio3`, `SK_cpd11416_c0/_c1/_c2`.

**On non-ModelSEED ids, nothing internal is priced.** On e_coli_core (BiGG ids, every reaction annotated), 0 of 75 internal reactions are priced. The objective is identical with every p = 1 and every p = 0 (10,427.422669; identical flux vectors), so **the probabilities have no effect at all**. Internal Σ|v| is 390.34 at the returned vertex but up to 79,883.15 on the optimal face. The function silently degenerates to "minimize exchange flux". b6857a0's own test renames e_coli_core's reactions to `rxnNNNNN_c0` for this reason, but the docstring does not warn about it.

### N5. Exchanges priced `ex_weight` in both directions

l.189-191 price both the forward variable (secretion, for `cpd_e0 <=>`) and the reverse variable (uptake) at 100. On the fixture with bio1 and ENV (N5):

- **The exchange term dominates.** It is 18,623.301851 of the 18,799.709290 objective (**99.06%**), consistent with the manuscript's 98.9% median (l.297). The internal, probability-weighted term is 176.41.
- **Half of it prices secretion.** Uptake is 94.09 and secretion 92.14, so secretion is 49.48% of the exchange term. The largest exchanges are CO2 out (+49.39), H2O out (+42.73), pyruvate in (-31.73), O2 in (-29.15), H+ in (-22.31) and NH3 in (-9.61). **H+ and H2O alone are 34.92% of the exchange term.** The term mostly prices obligatory byproducts and proton/water bookkeeping, which depends on protonation conventions.
- **Consequence 1: the optimal flux distribution changes.** These are measured deltas only. With uptake-only pricing on the same constraints, bio1 stays at its floor (1.404069) but the flux vector moves by L1 14.37:

  | | code (both directions priced) | uptake-only pricing |
  |---|---|---|
  | total uptake | 94.09 | 93.60 |
  | total secretion | 92.14 | 92.71 |
  | internal Σ\|v\| | 759.41 | 765.16 |
  | active `rxn` reactions | 362 | 363 |
  | CO2 secretion | 49.39 | 50.10 |
  | H2O secretion | 42.73 | 42.58 |

  The pricing of secretion is therefore not neutral, even on a single genome.
- **Consequence 2: the internal weights act only as a tie-break.** At ex_weight = 100, internal weights choose only among distributions of equal exchange flux (the manuscript, l.297).
- **Consequence 3: free outlets win over secretion.** Any unpriced outlet (N4) beats secretion at 100 per unit. That does not happen on this fixture, but it is structurally possible.
- **Consequence 4: the community rationale depends on the kinetic rows.** The manuscript's rationale for the term ("maximizing intracellular transport and thereby interactions among community members", l.287) applies to communities. There, secretion to the environment competes with transfer through the shared pool, and that transfer is priced by probability. But per N6 the community kinetic rows never reach the solved model.

### N6. Community path: the kinetic rows, `kinetics_coef`, `commkinetics` and `abundances` never reach the solved model

This runs on b6857a0, or on c2e31a0 with F2 (FP3c). It uses the two-member fixture of b6857a0's test, with a spy subclass of `mscommunity.mscommsim.MSCommunity` (N6).

- **The rows go to the copy.** `MSCommunity` is handed the solved model (True) but constrains a different object (False). It does `self.util = MSModelUtil(model, True, …)`, i.e. `copy=True` (mscommsim.py:230), and installs `CommKineticPkg` on `self.util.model` (mscommsim.py:291-293). The solved model has 19 constraints and **no** `*_commKin` rows; the copy has 20, including `memA_commKin` and `memB_commKin`. `MSCommObj` is never used after l.172.
- **`kinetics_coef` has no effect.** The objective is 763.800000 at both 1,000 and 1.
- **`commkinetics=True` builds no kinetics anywhere.** A truthy value skips MSCommunity (`if not commkinetics`, l.170) and delegates kinetics to `build_from_species_models`, which accepts `commkinetics` but never reads it: the name occurs once in its source, in the signature (commhelper.py:45). The commented l.164 shows the intended default. Neither setting therefore puts kinetics on the solved model.
- **`abundances` diverge between the two models.** They reach the solved model only through `build_from_species_models`: community bio1 on the solved model has c1 -0.9 and c2 -0.1. They are not forwarded to `MSCommunity(…)`, which assumes 1/n (mscommsim.py:263) and rewrites *its copy's* bio1 to -0.5/-0.5 (`set_abundance`, mscommsim.py:300-311). The two models hold different community biomass definitions. `close_member_drains` would also act on the copy only.
- **MSCommunity adds `rxn00062_c1/_c2` to its copy only, without notes.** If the call site were changed to solve MSCommunity's model, the default `missing_probability` would raise on them.

**History.** At 04e940d (2024-03-10), prFBA explicitly built the rows through the *solved* model's package manager:

```
pkgmgr = MSPackageManager.get_pkg_mgr(mdlUtil.model)
pkgmgr.getpkg("CommKineticPkg").build_package(kinetics_coef, MSCommObj)
```

The in-package `MSCommunity` already copied then (`self.util = MSModelUtil(model, True)`), so its members indexed the copy. Whether those 2024 rows were correct was not tested here. **fdaf0a2** (2024-03-13 21:06:22 -0500, "Andrew's WorkMacbook", "MSprobability and CommKinetics updates") commented those two lines out and left only the constructor call. From then on, no kinetic row has been built on the solved model in this lineage.

0050e06 (2026-09-01, "Andrew Freiburger", "Remove the community subpackage…") switched to the standalone package and added `kinetic_coeff=kinetics_coef`, with the comment "the MSCommunity constructor installs the community kinetic rows" (l.171). That is true only of MSCommunity's copy. The standalone MSCommunity has copied since its initial commit (6246555, 2024-08-30; the current line is from a9eef4ca, 2025-02-27).

cshenry/main still has the explicit on-model build, but cannot run (section 8). The defect therefore sits at the ModelSEEDpy call site; it is not internal to MSCommunity. One remedy is to solve the model MSCommunity constrains.

**Manuscript impact.** The manuscript (prFBA.tex l.264-272) describes a `CommKinetics` row "in a multi-member community model" at kinCoef = 1000. The released function never applies one to the model it solves. No reported result uses the community path (WEIGHTED_OBJECTIVE_RESULTS.md section 10: "the community branch of the released function was not run"), but the methods text overstates the code.

### N7. In-place mutation and process-lifetime caches

- **The caller's model is changed in place:** medium replaced, `C_elements` and `min_value` added, objective replaced by the minimization (N7: direction `min`, medium 347 entries).
- **A second call on the same object raises** `PackageError` on both trees (C7 at b68 l.235; N7).
- **Every call leaks.** Each call leaves one entry in each of `MSPackageManager.pkgmgrs` and `MSModelUtil.mdlutls`, which are class-level dicts keyed by model and hold strong references (N7: (0, 0) → (3, 3) after three calls). The harness needed `purge_caches()` (l.198-204) to avoid holding every copy. b6857a0's docstring documents the in-place change and the `PackageError`, but not the leak.

### N8. Probability values are not validated

The coefficient is `max(min_prob, 1 - float(p)**prob_exp)` (N8):

- `"nan"` → 0.01: a NaN becomes "certain", because `max(0.01, nan)` returns 0.01.
- `"1.5"` → 0.01.
- `"-0.5"` → 1.5.
- `"0.999"` → 0.01, tied with p = 1.
- `prob_exp=0` would make every coefficient min_prob.

### N9. Hard import-time dependencies

The module imports `from mscommunity import MSCommunity` at l.3 (added at 954b5f6, 2025-09-30, "Freiburgermsu", "updating the probability code") and cobrakbase at l.1 and l.4. The lazy imports inside `prFBA` (l.158-159) are therefore redundant, and 0050e06's commit message ("lazy, so no hard dependency at import time") is wrong.

`import modelseedpy.core.msprobability` fails without mscommunity. That is why b6857a0's tests `importorskip` both packages, and why codiffusion's `util.py:19` import drags them in.

### N10. Historical: from 2024-09-05 to 2026-09-01 every call raised at its first statement

The function's first statement imports `build_from_species_models` from `modelseedpy.community.commhelper`. fc33219 (2024-09-05, "Freiburgermsu_Sequoia", "CommScores updates") deleted `modelseedpy/community/commhelper.py` along with most of the subpackage. Later merges brought back a `community/__init__.py` that imports the missing `commkineticpkg`.

Until 0050e06 repointed the import to `mscommunity.commhelper`, every call in the c2e31a0 lineage raised `ModuleNotFoundError` before reaching FP1-FP4. **This is still the state of freiburgermsu `main`** (971df5d, its live tip). All three call forms raise `ModuleNotFoundError: No module named 'modelseedpy.community.commkineticpkg'` at msprobability.py:156 → community/__init__.py:9 (`aux_branch_probe_origin_main.log`).

### N11. Two definitions of "exchange"

The objective prices ids with the *prefix* `EX_` (l.189). `ElementUptakePkg` and `add_medium` use `MSModelUtil.exchange_list()`, which matches the *substring* `'EX_' in rxn.id` (msmodelutl.py:680-681). They agree on ModelSEED models, but any id containing `EX_` after its first character would be capped and given a medium bound while being priced as an internal reaction, or left unpriced.

### N12. No test ever exercised `prFBA` before b6857a0

Pickaxe searches over every ref in the shared clone, including the origin, cshenry and upstream remotes as fetched, find no test, notebook or example that calls it:

- `git log --all -S'prFBA' -- tests/`
- `-S'MSProbability'`
- `-S'msprobability'`
- `-S'prFBA' -- . ':!modelseedpy/core/msprobability.py'`

The only hits are b6857a0 itself and an HTML table header in `msfbareport`'s template (92a59ff).

### Minor

- `expression` is unused; the manuscript says so (l.288).
- With one model, `abundances` is silently ignored.
- `min_prob` is a floor on the *cost*, not on the probability.
- `solution.objective_value` is the weighted cost, not growth.
- ElementUptakePkg's `build_package` has a mutable default (`exception_reactions=[]`, appended to at elementuptakepkg.py:29). It is harmless here because prFBA passes no exception compounds.

---

## 8. Which branches carry the code

`git branch -r` was used, together with `git ls-remote` (live tips) and `gh api` (read-only), and nothing was fetched. Upstream ModelSEED/ModelSEEDpy has **no** `msprobability.py` on `main` (a5ea486, 2022-12-15, the live tip and the local ref) or on `dev` (062714b, 2026-02-11, `gh api` 404). `MSProbability` exists only in the freiburgermsu and cshenry forks.

| Branch (tip) | prFBA text (function md5) | FP1 | FP2 | FP3 | FP4 | Runs? |
|---|---|---|---|---|---|---|
| atp-safe-and-drain-fixes = origin (c2e31a0), element-closure-pkg (164377c), fix-gapfill-reversals (e1254ba), remove-community-module (d884a50), thermo-strict-margin-and-units (fb30087), and their origin/ counterparts | c2e31a0 (`2cdef284…`) | yes | yes | yes | yes | no (Part A) |
| fix-prfba-single-model (b6857a0, local only, not pushed per `git ls-remote`) | fixed (`547e2a39…`) | fixed | silent zero (N1) | fixed | fixed | yes, with N1-N8 remaining |
| origin/main (971df5d, the live freiburgermsu default branch), main (c347511), fix/stale-msexpression-import = origin (8a0b617) | pre-0050e06 (`4372aea4…`) | yes | yes (5149cc3 is an ancestor) | yes | yes | no: ModuleNotFoundError first (N10) |
| origin/cshenry-main (cca64e5, 2024-08-29) | same text | yes | no (5149cc3 is not an ancestor; add_medium keeps its None guard) | yes | yes | no |
| **upstream-pr/sync-core-msprobability** (eee9b2d; **cshenry/ModelSEEDpy PR #48, OPEN**, "Sync core/msprobability.py from freiburgermsu fork") | same text, onto cshenry/main | yes | yes | yes | yes | no. Merging it would bring FP2 and FP3 into cshenry/main, which lacks them today. |
| cshenry/main (local b758b44; live tip ac02db7, whose `msprobability.py` is identical per `gh api`), and the ~50 `upstream-pr/sync-*` branches based on it | 13274cc (2024-03-17, Christopher Henry) (`2368ddd9…`) | yes | no (guarded) | no (one build) | yes | **No.** cshenry's `MSModelUtil` has `__init__(self, model)` only and no `add_medium`, `add_minimal_objective_cons` or `add_objective`. A bare model gives `TypeError` (no len) at l.186. `[model]` gives `TypeError: unhashable type: 'list'` (mspackagemanager.py:19). Two models give `TypeError: MSModelUtil.__init__() got an unexpected keyword argument 'copy'` (community/commhelper.py:83). See `aux_branch_probe_cshenry_main.log`. It keeps the explicit on-model `CommKineticPkg` build (N6). |
| upstream-pr/remove-community (a0f397b; cshenry PR #31, OPEN) | cshenry text with the MSCommunity import moved into the branch (`7ff2def4…`) | yes | no | no | yes | no (same API mismatch) |
| origin/to_merge (994a851, 2024-03-11) | early text (`f7c38956…`) | yes | no | no | yes | no |

Ancestry was checked with `git merge-base --is-ancestor <commit> <branch>` for f351cb4, 04e940d, fdaf0a2, 937829d, fc33219, 5149cc3, 954b5f6 and 0050e06.

---

## 9. Callers of `MSProbability.prFBA`, or of a copy of its logic

The search covered tracked files (`git grep`) and untracked or ignored code files (`git ls-files -o [-i]`, `*.py`, `*.ipynb`, `*.sh`, `*.R`) in all three repos. The patterns were `prFBA(`, `MSProbability`, `msprobability`, `max(min_prob`, `notes['probability']`, `prob_exp`, `ex_weight`, `ElementUptakePkg` and `add_minimal_objective_cons`. Nothing was modified.

| Repo | File:line | What it does | Goes through the released function? | Does a reported result depend on it? |
|---|---|---|---|---|
| prFBA | `region_validation/shared_reference/scripts/weighted_objective_test.py:125-163` `prfba_released_copy()`, called at l.830 for every arm of every amplicon | **Private copy.** It unwraps one model (edit 1), drops the duplicate build (edit 2) and drops the print and the community branch (edit 3); every bound and coefficient line is verbatim. It avoids FP1-FP4 by construction. The caller sets bio1 first (l.829, avoiding FP6), annotates every `rxn` reaction with gapfilled ones at 0.0 (l.784-786 and `set_probs` l.457-460, avoiding FP5), and always passes `environment` (avoiding FP2). | No | **Yes**: WEIGHTED_OBJECTIVE_RESULTS.md and the manuscript's in-silico test (prFBA.tex l.291-301) |
| prFBA | same file, l.168-195 `released_prfba()`, called at l.393 (`worker_init`), l.1425 (`main`) and l.847 (per-row `released_check`) | Exec's the imported module's **released source** with two asserted edits (F1, F2) as a cross-check | Released source, edited | It validates the copy. In `per_amplicon.jsonl`, 360 of 360 rows have `released_check.pass` with rel_diff 0.0. Source sha256 `6eeecfb7…` equals c2e31a0's `prFBA` (recomputed). On b6857a0 the edit target is absent, so `worker_init` and `main` raise `RuntimeError` even with `--no-check-released`, which gates only the per-row check (l.842, 1395). |
| prFBA | `region_validation/shared_reference/data/weighted_objective/patches/test_patched.py` | Loads a *patched scratch copy* of the module (probe) | No | No (probe log only) |
| prFBA | `presentations/make_deck.py:3016`; prFBA.tex l.264-273, 288, 300-301 | Text describing the function | — | Text only. l.264-272 describes community kinetics that the code never applies (N6). |
| codiffusion_bioreactor | `util.py:19` | `from modelseedpy.core.msprobability import MSProbability`, never used. The import pulls in mscommunity and cobrakbase (N9). | Import only | No |
| codiffusion_bioreactor | `build_prfba_community_models.py:362`, `load_genomes_to_modelseedpy.py:355` (+ `.ipynb_checkpoints/` copies) | **Write** `notes['probability']` (producers) | No | No |
| codiffusion_bioreactor | `compartment_aware_gapfill.py:76` | Reads `notes.get("probability", 1.0)` for its own gapfill scoring. Its default of 1.0 for a missing note is the opposite of the prFBA pipeline's 0.0. | No | No |
| EmilyKin | `community_modeling/ek_pipeline.py:885`, `load_genomes_to_modelseedpy.py:355`, `codiffusion_copies/{build_prfba_community_models.py:355, load_genomes_to_modelseedpy.py:355}` | Producers of the note | No | No |
| EmilyKin | `community_modeling/codiffusion_copies/util.py:19` | Import only | Import only | No |
| EmilyKin | `community_modeling/_pinned/ModelSEEDpy-{c2e31a0,d884a50}/…/msprobability.py` (untracked) | Library copies, byte-identical to c2e31a0 | — | No. `plan_evidence/round11/codif_port/modelseedpy.md:209` records that EmilyKin never calls `MSProbability.prFBA`. |

**Conclusion.** No code anywhere calls the released `prFBA`, so none could have hit FP1-FP6 in production. The only solved prFBA objective behind a reported number is the harness's private copy. That copy guards against FP5 and FP6 explicitly, and it was shown equal to the released source with the two failing lines corrected (360 of 360). The cross-check edits FP2's line only implicitly, because it always passes an environment.

The manuscript's remaining exposure is textual: the community `CommKinetics` description (N6), and the "min_biomass" wording for a floor that the code puts on the entry objective (FP6).

---

## 10. Corrections and omissions in `PR_ModelSEEDpy_prfba_single_model.md`

The PR's factual claims about the four fixed defects, the keyword behavior, the file statistics (+77/-10, 248 lines, 9 tests), the base (PR #19 is `fix-gapfill-reversals → atp-safe-and-drain-fixes`), msatpcorrection.py l.606, and the 19-versus-20 constraint observation all check out, the last also at MSCommunity 03bd80c. These items need changing:

1. **Community kinetics, "lives in another package, so it is not fixed here"** ("Not changed", first bullet). The attribution is inaccurate. The rows go to a copy because the ModelSEEDpy call site (l.170-172; b68 l.237-239) constructs `MSCommunity`, discards it and solves `mdlUtil.model`, while MSCommunity has always worked on `MSModelUtil(model, copy=True)`. Until fdaf0a2 (2024-03-13), prFBA explicitly built `CommKineticPkg` through the solved model's package manager, and cshenry/main still has that line. 0050e06 added the comment "the MSCommunity constructor installs the community kinetic rows", which is false for the solved model. This is a prFBA-side regression and could be fixed in this PR or a follow-up, for example by solving the model MSCommunity constrains. The bullet also omits three things (N6): `kinetics_coef` has no effect (763.8 at 1,000 and at 1); `commkinetics=True` builds no kinetics anywhere, because `build_from_species_models` never reads the flag; and `abundances` are not forwarded to MSCommunity, whose copy is renormalized to 1/n.
2. **Omission: `environment=None` now returns a silent zero** (N1). On the gapfilled models the PR targets, "keeps the model's own exchange bounds" means bounds that admit no growth (fixture: bio1 max 0). The fixed call then returns `status=optimal`, objective 0, all fluxes zero (C6, FP1c). The docstring should warn, or the function should raise when the floor is 0 or NaN. An infeasible medium gives a raw `GurobiError` (N2). The same imprecision is in WEIGHTED_OBJECTIVE_TEST.md section 2, whose "own exchange bounds admit no uptake" should read "admit no growth": 332 uptakes are open on the fixture.
3. **Omission: the root cause of item 2.** The table describes the symptom only. The crash comes from `MSModelUtil.add_medium` losing its `if media is None` guard at 5149cc3 (2024-09-13); prFBA's own guard was dropped at 937829d (2024-08-26). The PR fixes the call site only, and `add_medium(None)` still raises on b6857a0 (FP2b). That is consistent with the other callers, which all guard, but worth one sentence.
4. **"`abundances`, `commkinetics`, `kinetics_coef` and `expression` are unchanged. `expression` is still unused."** This is true but incomplete. `commkinetics` and `kinetics_coef` are functionally unused as well: `commkinetics=True` delegates to a path that does nothing, and `kinetics_coef` reaches only MSCommunity's copy (N6).
5. **Equivalence section, "a literal copy with edits 1 and 3 and the duplicate package dropped".** This garbles the harness's own numbering, where (1) is the unwrap, (2) is dropping the duplicate and (3) is dropping the print, `printLP` and the community branch. It should read "edits (1)-(3)".
6. **"(MSCommunity 176b6e1)".** MSCommunity HEAD is now 03bd80c (2026-10-31). The observation reproduces there, and the copy line dates from a9eef4ca (2025-02-27). Fine as a pin, but say so.
7. **Worth adding for reviewers.** The single-model path has not run in any version examined since the first draft f351cb4 (FP1), and no test ever exercised `prFBA` before this PR (N12). freiburgermsu `main` fails even earlier, at the deleted `modelseedpy.community.commhelper` import (N10). The open cshenry PR #48 would propagate FP1-FP4 into cshenry/main, which today lacks FP2 and FP3.
8. **Minor docstring omissions.** Non-ModelSEED reaction ids leave every internal reaction unpriced, so the probabilities have no effect (N4, e_coli_core). Each call leaves the model in two class-level caches (N7). The module-level `from mscommunity import MSCommunity` (l.3) makes the in-function import redundant and contradicts 0050e06's "lazy" note (N9).

## 11. Commands used (all read-only)

```
git -C ~/Documents/ModelSEEDpy log --full-history --format='%h %ad %an | %s' --date=iso c2e31a0 -- modelseedpy/core/msprobability.py
git -C ~/Documents/ModelSEEDpy blame -L 154,202 --date=short c2e31a0 -- modelseedpy/core/msprobability.py
git -C ~/Documents/ModelSEEDpy show 937829d -- modelseedpy/core/msprobability.py      # the ++ lines (FP2, FP3)
git -C ~/Documents/ModelSEEDpy log -L '/def add_medium/,+10:modelseedpy/core/msmodelutl.py' c2e31a0   # 85156f0 / 5149cc3
git -C ~/Documents/ModelSEEDpy show <commit>:modelseedpy/core/msprobability.py       # each version of prFBA
git -C ~/Documents/ModelSEEDpy log --all -S'prFBA' -- tests/                           # N12
git -C ~/Documents/ModelSEEDpy merge-base --is-ancestor <commit> <branch>              # section 8
git -C ~/Documents/ModelSEEDpy branch -r ; git ls-remote <repo> refs/heads/main ...    # section 8
gh api repos/ModelSEED/ModelSEEDpy/contents/modelseedpy/core?ref=dev ; gh pr list --repo cshenry/ModelSEEDpy
git -C ~/Documents/MSCommunity log -S'MSModelUtil(model, True' -- mscommunity/mscommsim.py
```
