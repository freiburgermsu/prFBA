# Run 4, TF lane: R4.3 post-hoc studies, 2026-11-16

Scope: calibrate_pfunc, Stage E (`pfunc_estimators.py`) and its verification, the producer of
`stage_o_weights_E3.json`, Stage O (`objective_tuning.py` + `objective_tuning_results.py`), the FVA variant, the
objective wrapper equivalence, and Stage O's verification. The weighted-objective test (wobj) was the other Run 4
agent's: it ran on FM.

Executor: Run 4 TF agent. Progress log with timestamps: `logs/RUN4_TF_PROGRESS.md`. TF unit logs and memory
samples: `logs/tf/r43/tf_logs/`; TF unit records: `logs/tf/r43/units_tf.jsonl`. Commit messages:
`logs/commit_msgs/r4tf_*.txt`.

**Outcome.**
- Every stage ran on TF through `tf_unit`/`prfba-ns`, and every control passed.
- **No selected configuration and no verdict changes.**
  - Stage E: same winner (E4(E3)), same per-fold variants and rankings, all SESOI verdicts unchanged.
  - Stage O: same selection (est/lin/G0/X1), the −0.0172 headline unchanged at print, every label unchanged.
- Three commits of mine, path-limited, not pushed: `95675ac`, `bb5dc5b`, `06251c8`. `2520f4f`, in the same range,
  is the other agent's wobj commit.
- The cache was never written. `stamp_state --compare R3.final`: CACHE UNCHANGED at the start
  (`stamps/R4TF.0.json`) and at the end (`stamps/R4TF.final.json`).
- TF's clone is at FM HEAD `06251c8` with a clean tracked tree.

---------------------------------------------------------------------------------------------------------------

## Preflight and inputs

- **Stamp.** `stamp_state.py --compare R3.final` → `stamps/R4TF.0.json`: CACHE UNCHANGED (`05f03bcd…`), catalogue
  current. d20 baseline 15/32: the failing 17 are exactly the plan's R5.5 list.
- **Code.** TF was already at FM HEAD `1a64f47`. It was then synced to `95675ac` (the producer commit) and, at the
  end, to `06251c8`. Before the last sync, TF's 35 locally modified tracked outputs were reverted after checking
  them against FM HEAD's blobs: 34 were equal; the 35th is the Stage E README, which on FM has the corrections
  section re-appended.
- **Pushes, all verified (sha256 + size + mtime):**
  - `push_R4.3_posthoc.txt`: 68 files VERIFY OK. TF read `cells.csv.gz` d57709e3, `.fp` 224f90fb and
    `_cache_panel_ec.json` 8eb296b1.
  - The 161 static SR inputs (list C minus wobj/otune) are sha256-identical. Only the mtimes of the 40 files that
    R3.4 rewrote byte-identically differ.
  - `population.json` 5d13cf1e, after the marker `R4_WOBJ_POPULATION_DONE`.
  - wobj `per_amplicon.jsonl` d44a27e7, after the FM wobj main run ended at 13:54. Stage O needs it as the
    "prior" for `prior_match`, `released_check`, the est_shuffled permutation and g_meas.
- **TF set-aside.** TF's stale copies were moved to TF `$W/pre_repair/` before any stage, keeping the `*.py`, the
  VERIFICATION docs, the template and `smoke/`. This covers the pre-re-score `weighted_objective/{population,
  per_amplicon}`, `objective_tuning/per_amplicon.jsonl` (Stage O resumes finished keys) and
  `calibration/scores.json`.
- **D11 pins.**
  - Stage O, FVA and the wrapper read ModelSEEDDatabase `5389554a`, passed as `--msd …-5389554a`, so
    `meta.json args.msd` records the commit. The default-path symlink also points there.
  - TF's tree hashes of Biochemistry `b0717c45`, Media `93b22cac`, Objects `56de5184`, Libs `ba60ff07`, README and
    LICENSE were recomputed today: they equal FM's git objects at 5389554a.
  - ModelSEEDpy-c2e31a0, ModelSEEDpy-prfba `09856f5` and Templates-c12run `21a8084` are clean at their commits.
  - wobj (FM, other agent) recorded `--msd …-e8446e3e`.
  - calibrate and Stage E use no ModelSEED tree; they ran with their committed `meta.args`, and fresh caches
    `$W/calib_cache_r4` and `$W/pfe_cache_r4`.
- **Code-path equivalence.** No R4.3 script, and nothing they import, changed since `cc8014f`, the code of the
  runbook §7 TF equivalence runs: `git diff --stat cc8014f HEAD` on the SR scripts shows only the new
  `refresh_panel_ec_cache.py`. Runbook §7 (a)–(e) therefore remains the "today's code on pre-repair inputs
  reproduces the committed numbers" control for every stage here.

## The missing producer: `scripts/stage_o_weights_e3.py` (commit `95675ac`)

**What it does.**
- Copies the `amplicons` unchanged.
- Sets `meta.estimator` to the winner's base family (E4 → E3) and `meta.isotonic` to null.
- Adds `meta.derived_from {file, sha256 of the source bytes, change}` and `meta.generated_e3`, keeping every
  other meta key in order.
- Serialises with `json.dumps` (default separators, no indent, no trailing newline), the form of the 11-03 file.

**The change sentence** reads its numbers from `scores.json`, rather than carrying them as text: the E4(E3) − E3
pairwise contrasts (d100 "all": universe Brier, union Brier, macro F1 at q ≥ 0.5), the amplicon count and the
map's maximum.

**It refuses** unless the source and `scores.json` come from one Stage E run with an E3-based winner. It checks:
- the base family;
- the estimator label;
- `script_sha256`;
- K, τ, dedup and α against the carried-forward variant;
- the isotonic map against the winner's `isotonic_all_folds`.

**Control.** Run on the committed 11-02 Stage E files with `--generated-e3 2026-11-03T08:50:15-06:00`, the output is
**byte-identical** to the committed `stage_o_weights_E3.json` (sha256 330750c4…, 931,791 B). No serialisation
difference.

**Tests.** `tests/test_stage_o_weights_e3.py`, 8 pass: derivation and serialisation, text that follows the numbers,
four refusals, the isotonic mismatch, and the committed-files consistency check. They were re-run after `06251c8`
and still pass, so the committed re-scored trio is self-consistent: `stage_o_weights.json` 54b781fc, `scores.json`
and `stage_o_weights_E3.json` ead0d11e.

## Per stage

| stage | where | unit wall | cgroup peak / max RSS | MemoryMax | result |
|---|---|---|---|---|---|
| calibrate_pfunc (4 workers) | TF | 391 s | 15.3 / 4.9 GB | 20G | rc 0; `all_required_zero` true |
| Stage E (12 workers) + E3 producer | TF | 350 s | 14.4 / 9.2 GB | 20G | rc 0; 88/88 required assertions |
| Stage E verification p0–p8 | TF | 289 s | 9.4 / 7.0 GB | 40G | all rc 0 |
| Stage O (10 workers, maxtasksperchild 1) | TF | 3,354 s | 15.9 / 2.2 GB | 34G | rc 0; validate-only ok |
| FVA variant (10 workers) | TF | 1,442 s | 13.2 / 1.7 GB | 34G | rc 0 |
| wrapper equivalence | TF | 292 s | 2.0 / 2.2 GB | 8G | rc 0, all_pass |
| results with `--fva-dir` + Stage O verification | TF | 5 s | 0.2 GB | 16G | all rc 0 |

- FVA and the wrapper ran concurrently (42G of MemoryMax). Every other unit ran alone, so at most 42G of MemoryMax
  was in use at once.
- TF wall time: 11:54–15:16. Idle from 12:12 to 13:55, waiting for the FM wobj main run (7,094 s on FM).

### Controls (with evidence)

- **calibrate, Stage E, Stage O, FVA, wrapper.** Runbook §7 (a)–(e): on the pre-repair inputs, TF with today's
  code reproduced FM's committed outputs bit for bit (Stage O 360/360 arm objectives on the subset; FVA 15/15;
  wrapper all fields; calibrate and Stage E every number). The code has not changed since.
- **Stage O results script, on FM today.** `objective_tuning_results.py` on the 11-03 rows reproduces the
  committed `summary.json` exactly. Only `generated`, `script.sha256` and `outdir` differ: the script was edited
  after 11-03 without any numeric effect.
- **Doc-number extractors.** Both are controlled on the 11-02/11-03 outputs:
  - `$W/tools/pfe_doc_numbers.py`: every one of the 422 mapped tokens of `PFUNC_ESTIMATORS_RESULTS.md` reproduces
    from the old outputs;
  - `$W/tools/otune_doc_numbers.py`: 0 failures, and it reproduces `OBJECTIVE_TUNING_RESULTS.md` byte-identically
    from the old outputs.
- **Exchange-rise count.** `$W/tools/otune_exchange_rise.py` reproduces 95/167 and 110/193 (the "57%"), which the
  11-03 verifier had counted ad hoc.
- **Built-in checks.**
  - Stage E: 88/88 assertions. One caveat: the assertion named `b_…_reproduces_0.06273832` passes against the
    **re-scored** calibration target 0.06275733 (|d| 4.2e-9); only its name carries the 11-02 value.
  - Stage O: validate-only ok; `released_check` 360/360; prior_match failures are the same 2 rows as on 11-03;
    `coef map = prfba_released_copy` 360/360; determinism pairs pass.
- **Independent verifications** (re-run, exit 0):
  - Stage E: E0–E3 reproduce within 5.0e-8, E4 within 1.6e-4, leakage 0, selection agrees in every fold.
  - Stage O: 972 confirm-half numbers equal `summary.json` (max diff 0.0); 36/36 guard rows.

### Headline numbers, old → new

**Calibration** (`data/calibration/README.md`, regenerated).
- Unchanged at print: EC reliability 0.0210 [0.0194, 0.0242]; WLS slope 0.78 [0.76, 0.79]; multi skill 0.061
  [0.034, 0.074]; n = 2 skill 0.010 (covers 0); thresholding +0.0007.
- The as-shipped secondary moves: exclude_0.987 skill 0.264 → 0.277 and slope 0.88 → 0.91; include-self
  0.624 → 0.638.

**Stage E** (`PFUNC_ESTIMATORS_RESULTS.md` §3; nested CV, d100, 57,021 cells). The 360 Stage O weight entries are
identical; the isotonic map was refitted (747 → 725 knots; r(0) 0.008, r(0.5) 0.48, r(1) 0.943).

| | E0 | E1 | E2 | E3 | E4(E3) |
|---|---|---|---|---|---|
| universe Brier | 0.0889 | 0.0593 → **0.0594** | 0.0672 | 0.0578 | 0.0575 |
| universe log-loss | 0.4145 → **0.4146** | 0.2113 → **0.2114** | 0.2571 → **0.2572** | 0.2018 | 0.1977 |
| union-support Brier | 0.2087 | 0.1301 | 0.1494 | 0.1267 | 0.1262 |
| union WLS slope | 0.699 | 0.841 | 0.801 → **0.802** | 0.903 | 0.940 |
| union AUROC | 0.543 | 0.799 | 0.734 | 0.816 | 0.815 → **0.814** |
| macro dF1 vs union | 0.0007 | 0.0295 | 0.0142 | 0.0317 | 0.0274 |

- Contrasts against E0: E1 union −0.0785 → **−0.0786**. For E3, universe −0.0311 [−0.0322, −0.0297] and union
  −0.0820 [−0.0871, −0.0789]; only the lower bound of the union CI moved, from −0.0870.
- Several CI bounds in the table move by one unit as well (`PFUNC_ESTIMATORS_RESULTS.md` §3).
- E4(E3) − E3: universe −0.0003 [−0.0007, +0.0001]; union −0.0005 → **−0.0006** [−0.0024, +0.0010]; macro F1
  −0.0042 → **−0.0043**.
- n = 1: E3 union AUROC 0.805; information gain over E4(E0) −0.0415.
- Within-union macro dF1: +0.0085 [0.0072, 0.0096 → **0.0097**].
- In all, 63 printed values moved. Apart from the 7–8-digit reproduction values and the knot count, each moved by
  at most two units in the last printed digit.

**Stage O** (`OBJECTIVE_TUNING_RESULTS.md`).
- Selected: **est/lin/G0/X1, unchanged**. Pooled A 0.1196 against the shipped 0.1362; tie set {est/lin/G0/X1};
  naive selection agrees.
- **Confirm-half selected − shipped, pooling A: −0.0172 [−0.0221, −0.0118], unchanged** (−0.0171875015 against
  −0.0171875014). Face min −0.0076 [−0.0133, −0.0047]; face max −0.0268.
- **Decomposition, unchanged:** est alone at X100 −0.0128 [−0.0186, −0.0088]; X1 adds −0.0043 [−0.0061, −0.0017];
  additivity −0.0169 vs −0.0172; log +0.0039 [0.0011, 0.0053]; G +0.0000; volume term +0.0002 [−0.0016, +0.0017].
- **n = 1, unchanged:** −0.0048 [−0.0113, +0.0063] on the confirm half; −0.0124 on the tune half. Pooling B −0.0065.
- Moved:

  | quantity | old → new |
  |---|---|
  | tune half | −0.0165 → −0.0166 |
  | shuffled comparator | −0.0282 → −0.0283 |
  | headroom | 0.685 / 0.231 / 0.933 → 0.687 / 0.233 / 0.935 (tune 0.908 → 0.910) |
  | selection-bootstrap share | 0.932 → 0.942 |
  | selected arm's median gene-backed flux ratio to shipped | 1.016 → 1.011 |
  | G-varying n2 rows (est/lin/*/X1) | 2 → 3 |

  47 printed values moved, each effect by at most 0.0003.
- Robustness, unchanged: leave-one-group-out [−0.0174, −0.0162]; without grp00000 −0.0164 [−0.0236, −0.0100];
  excluding the 3 rebuilt-LP rows −0.0171 [−0.0220, −0.0119].
- Row-level change: 53/360 rows and 1,827/32,400 arm solves changed (est's EC sets via the 17 repaired genomes);
  the 2 changed population jobs changed their truth counts.
- FVA: 20/193 rows changed. Degeneracy arm levels 13.6 / 35.4 / 41.5 / 2.1, unchanged.

**Wrapper equivalence.** all_pass; tuned and released bitwise on 8/8. Every number equals the 11-11 file. Only
provenance differs: TF library paths, HEAD, `prfba_objective.py` sha (6e2c2eb, as the runbook notes), the weights
sha and the commit-named MSD path. All 8 amplicons are unchanged rows.

**`collect.py`** over the 11 R4H items (`control/collect_run4_tf.{json,md}`):
- 49 values; every old value equals the print; 0 errors; 3 changed at print precision.
- Mine: Stage O "beat shuffled" 0.0282 → 0.0283 (L319-321).
- The other agent's wobj values: "beat shuffled" 0.0120 → 0.0121; "shuffled free reactions" 27.9 → 27.8.

## Docs updated

All have a dated re-score note; declared-criteria text is unchanged.

| doc | how it was updated |
|---|---|
| `PFUNC_ESTIMATORS_RESULTS.md` | explicit token→path map |
| `data/pfunc_estimators/README.md` | regenerated, and the 11-03 corrections section re-appended (its counts are unchanged) |
| `data/pfunc_estimators/verification/VERIFICATION.md` | re-score section with the verification's key results, old → new; the body is kept as the 11-03 record |
| `OBJECTIVE_TUNING_RESULTS.md` | generated lines realigned from the new filled document; hand lines mapped token by token |
| `data/objective_tuning/verification/VERIFICATION.md` | re-score section, as for Stage E |
| `data/calibration/README.md` | regenerated by the script |

`data/objective_tuning/OBJECTIVE_TUNING_RESULTS.filled.md` is now the script's own output. The 11-03 version had the
hand corrections; those live in `OBJECTIVE_TUNING_RESULTS.md`.

**Judgment call.** Only the two results docs had their bodies re-numbered. Both `VERIFICATION.md` bodies are kept
as the 11-03 record. Each gets a dated section with the verification's own key results, old → new, and points to
the results docs for the current study numbers.

Other provenance notes:
- The Stage O §2 wall-time row now shows TF's timing.
- `wrapper_equivalence.json` records TF's library paths.

## Commits (prFBA, branch `hardware-aware-alignment-dispatch`, not pushed)

| commit | content |
|---|---|
| `95675ac` | `scripts/stage_o_weights_e3.py`, `tests/test_stage_o_weights_e3.py` |
| `bb5dc5b` | calibration (README, assertions, reliability.csv); `pfunc_estimators/*` and its verification outputs and VERIFICATION.md; `PFUNC_ESTIMATORS_RESULTS.md`. Git-ignored `calibration/scores.json` sha256 395c42b7. `pfunc_estimators/run.log` is reconstructed: the `[pfunc]` lines of the TF unit log plus `EXIT=0`, 116 lines like the 11-02 file, which was the run's stdout |
| `06251c8` | `objective_tuning/{meta,summary,wrapper_equivalence}.json`, `per_amplicon.jsonl.gz` (gzip -9 -n, the 11-03 form), filled.md, `verification/out_*.json` + VERIFICATION.md; `objective_tuning_fva/{meta.json,per_amplicon.jsonl.gz}`; `OBJECTIVE_TUNING_RESULTS.md`. Git-ignored rows: d12f7a36 (Stage O), 7ef13c3d (FVA) |

The 11-02/11-03 outputs (and copies of the edited docs) are in
`$W/pre_repair/region_validation/shared_reference/{PFUNC_ESTIMATORS_RESULTS.md,OBJECTIVE_TUNING_RESULTS.md,data/{calibration,pfunc_estimators,objective_tuning,objective_tuning_fva}}`.
The move logs are `logs/r4tf/moveaside_*.log`.

## What Run 5 needs

1. **Manuscript (R5.7), R4.3 values that change at print:**

   | line | old → new |
   |---|---|
   | L319-321 | "beat … shuffled estimates by 0.0282 [0.0195, 0.0382]" → **0.0283** [0.0195, 0.0382] |
   | L319-321 | "by 0.0110 [0.0031, 0.0222] against the released objective" → **0.0111 [0.0031, 0.0223]** |
   | L358 | "an out-of-fold isotonic recalibration … lowers it by 0.0411" → **0.0412** |
   | L360 | "within the union the gain is +0.0085 [+0.0072, +0.0096]" → [+0.0072, **+0.0097**] |

   - Everything else at L317-325 and L353-360 is unchanged at print: calibration slope 0.78 [0.76, 0.79]; +0.0007;
     0.0889 → 0.0578 (0.0311 [0.0297, 0.0322]); 0.209 → 0.127; 0.543 → 0.816; 0.805; 0.0408 [0.0384, 0.0424];
     0.739; 0.067 [0.062, 0.072]; +0.032; +0.0249 [+0.0228, +0.0265]; 13%; −0.0172 [0.0118, 0.0221]; 0.0076
     [0.0047, 0.0133]; 69%/23%; 0.0128 / 0.0043 / 0.0039; 13.6 / 35.4 / 2.1; n = 1 0.0048 / 0.0124; 96 of 193,
     0.0164 [0.0100, 0.0236]; 2.5 [1.2, 3.4].
   - The wobj values (0.0120 → 0.0121; 27.9 → 27.8) are the other agent's.
   - `val16s:pfunc_calibration` and `val16s:stage_e_estimator` remain "manual" in `claims_spec.json`. The check
     above was done by hand against the new outputs; a selector could now be added.
2. **Decks (R5.3).** `pfunc_estimators/{scores,selected,assertions}.json`, `calibration/{scores.json,reliability.csv}`
   and `objective_tuning/{summary,wrapper_equivalence,verification/out_c8}.json` are regenerated.
   - `make_deck.py` L424 `PFE_WITHIN_UNION_F1["all"]` = (0.0085, 0.0072, 0.0096): the last value becomes **0.0097**.
     `["n2+"]` = (0.0249, 0.0228, 0.0265) is unchanged.
   - `OTUNE_LOG_CI` = (0.0011, 0.0053), L420, is unchanged.
3. **WEIGHTED_OBJECTIVE_RESULTS.md** and `wobj_results.py` belong to the other agent (Stage O read only its
   `population.json` and `per_amplicon.jsonl`).
4. **Disk and TF state.**
   - `$W/from_tf/pull_R4.3_tf_*` (≈ 0.2 GB) and TF's `$W/{calib_cache_r4,pfe_cache_r4,r43_prelim}` can be deleted
     after review. TF's `$W/pre_repair` holds TF's old copies.
   - TF's canonical git-ignored outputs are now the new ones: `objective_tuning{,_fva}/per_amplicon.jsonl`,
     `calibration/scores.json`, `weighted_objective/{population.json,per_amplicon.jsonl}`. A future TF re-run of
     these stages must move them aside first.
5. **Runbook corrections.**
   - Stage O needs the wobj **main-run rows**, not only `population.json`, so Stage O waits for the full wobj run.
   - The Stage E README carries a hand-appended section that every re-run drops; it was re-appended here.
   - The c8/vlib constants (`SELECTED`, `EXC2`) are hard-coded, so check them whenever the selection or the
     prior-match failures change. Both still hold.
