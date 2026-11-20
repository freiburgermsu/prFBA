# Run 2, TF lane: R2.1 (E3 gate + review) and R2.2 (MAG panel arms), 2026-11-15

Executor: the Run 2 TF agent. Scope: plan §0–4 (Run 2) and §8, and TF_RUNBOOK.md. Every prFBA command on TF
ran through `tf_unit`/`prfba-ns`. Progress log with timestamps: `logs/RUN2_TF_PROGRESS.md`. TF unit records:
`logs/tf/` (pulled from TF's `$W/logs`), TF `$W/stamps/units_tf.jsonl`.

**Outcome.**
- Every control passed. R2.2 and R2.1 were both produced, verified and committed (not pushed):
  - `abdfca6` (R2.2);
  - `e722aaa` (R2.1).
- `logs/R2.2_DONE` was written at 22:37, right after the R2.2 commit, so R2.3 could start.
- **The E3 verdict is unchanged: PASS.**
- **Nothing is blocking.**

---------------------------------------------------------------------------------------------------------------

## Preflight

**FM stamp.** `stamp_state.py --compare R1.6` → `stamps/R2.0_tf.json`:
- CACHE UNCHANGED (`05f03bcd…`, 159,753 entries); catalogue current.
- d20 baseline 22/32. The 10 failures are the expected R1.4 files.

**`tf_code_sync.sh`, first use: works.**
- First sync: TF cc8014f → 3ffc97d (37 s), tracked tree clean.
- Final sync: TF → `e722aaa` (FM HEAD). Before it, TF's 120 locally modified tracked outputs were reverted with
  `git checkout --`. Their blobs had first been checked equal to FM HEAD's, so the fast-forward could apply.
- TF's clone is now at FM HEAD with a clean tracked tree, ready for R4.3.

**Pushes.** `push_R2.1_e3gate.txt` + `push_R2.2_mag_indices_if_built_on_FM.txt`:
- 0.79 GB on the wire, 70 s;
- `tf_verify`: 8 files identical on FM and TF (sha256, size, mtime). This includes the cache, so TF read the R1.6 freeze.

**TF fixes and set-asides.**
- `mag_validation{,_km,_rnd}/logs` did not exist on TF. Every driver redirects into them with no `mkdir`, so every
  batch would have FAILED. They were created.
- Old outputs on TF moved to TF `$W/pre_repair`: the 3 old indices, and the 4 gate outputs. Tracked old batch
  dirs were deleted (their old versions are in git).

---------------------------------------------------------------------------------------------------------------

## R2.2 MAG panel arms (commit `abdfca6`)

### Indices (W1)

The FM `.new` indices were built 11-15 14:18–14:23 on cache `05f03bcd` (= R1.6). They still match their stamps.
They are installed at the canonical paths on FM and TF; the old ones are in `$W/pre_repair/` on both machines.

| index | new sha | old (= W1 `--restore-old-values` control) sha | repaired genomes |
|---|---|---|---|
| `mag_validation/pgfam_panel.json` | `cd5c3f49` | `dcde696b` | 140 |
| `kmeans_panel/pgfam_panel_kmeans.json` | `22bb5115` | `b597a00f` | 20 |
| `random_panel/pgfam_panel_random.json` | `7986e108` | `64c2ea64` | 38 |

- The `*.stamp.json` sidecars are committed. The indices themselves are git-ignored.
- **Verify on TF (`r22-w1-verify`): PASS.**
  - `--check-only` passes for all 3 against the pushed cache.
  - Independent TF rebuilds are byte-identical to the installed indices.
  - The `--restore-old-values` rebuilds are byte-identical to the old ones.

### Control (`r22-mag-control`): PASS

The W1 control indices were run with HEAD `3ffc97d` code on **all six arms + smoke** (121 batches, `-P 6`).
- **121/121 batch tables** (`mag_cells.csv` and `mag_summary.json`) are byte-identical to the old tables:
  - tracked: tau99, n500, k-means, smoke;
  - git-ignored `pre_repair` copies: main, lowc, random.
- Concatenated: **`mag_cells_all.csv` md5 a57df30b8346fcb72f14f7f33708669a, `mag_cells_tau99.csv` md5
  3cb0e487aa31c2f30f1699fcbdf88c92: MATCH** (199,940 rows each).

### Re-score (`r22-mag-rescore`)

- The six drivers ran as committed, plus the smoke set: **120 batches, 0 FAILED**.
- Concatenation, then `aggregate.py` ×4.
- The unit's rc = 141 is SIGPIPE from the final `git status | head` listing only. Every step reported rc 0.
- New md5: `mag_cells_all.csv` `74fa46de…`, `mag_cells_tau99.csv` `05ac68a5…`.

### Pull and install on FM (`tools/r22_install.sh`)

- 368 files pulled. Then:
  - FM's git-ignored `mag_validation/{out_batches,lowc_out_batches}` and `mag_validation_rnd/out_batches` moved to
    `$W/pre_repair/`;
  - the 120 driver logs that the pull overwrites were backed up to `$W/pre_repair/*/logs`;
  - the tracked batch dirs were removed, then the pulled tree was installed.
- **`tf_verify`: 248 output files + 120 logs identical on FM and TF.**

### Checks

**Rows.** 199,940 (main), 199,940 (τ99), 99,970 (n500), 20 lowc batches (179,946), k-means 99,960, random 99,980.

**W8 (`tools/r22_w8.sh` → `control/w8/mag/attrib_*.json`): 7/7 arms PASS.**
- 0 one-sided rows. **0 identity-column changes:** `n_selected`, anchor, tier and the rank cells are unchanged,
  so no selection changed.
- **0 unattributed.**

| arm | changed rows | source repaired | anchor repaired | multi-pick (n_sel ≥ 2) |
|---|---|---|---|---|
| main | 2,935 | 2,518 | 396 | 21 |
| τ99 | 3,082 | 2,530 | 531 | 21 |
| n500 | 1,456 | 1,260 | 192 | 4 |
| lowc | 2,639 | 2,265 | 357 | 17 |
| k-means | 212 | 123 | 89 | 0 |
| random | 645 | 328 | 314 | 3 |
| smoke | 27 | 20 | 5 | 2 |

- **All 68 multi-pick rows were checked as well.** Each has a repaired genome among its MAG's skani hit
  references (`control/w8/mag/multi_pick_check.json`).

**Resolve-rate prediction (`expected/resolve_impact.json`): all 38 points match exactly**, old and new: τ95 at 14
completeness levels × 2 regimes, and τ99 at 5 × 2. The outcome transitions match too (`control/r22_checks.json`).
- "Resolved" means ok with a non-empty prediction, as in `aggregate.py`.
- A first pass that counted `pred_empty` as resolved reported false DIFFs. That was the definition, not the data.

### Old → new (self-exclusion, contamination 0; full tables in `control/r22_headlines.txt`)

| arm | resolve 100% | resolve 1% | conditional F1 100% |
|---|---|---|---|
| main τ95 | 0.3075 → **0.3134** | 0.2533 → **0.2581** | 0.8068 → 0.8066 |
| τ99 | 0.3656 → **0.3722** | – | 0.8230 → 0.8220 |
| depth 500 | 0.3086 → 0.3145 | – | 0.8065 → 0.8063 |
| k-means | 0.3422 → 0.3424 | – | 0.7971 → 0.7969 |
| random (pooled from cells) | 0.4276 → 0.4277 | – | 0.7906 → 0.7896 |

- Population F1 (main, 100%): 0.2481 → 0.2528.
- Genus accuracy is unchanged everywhere. Self-inclusion stays at 1.000.

### Manuscript items (`collect.py --only` the 4 R2.2 items → `control/collect_run2_tf_r22.{json,md}`)

- 16 values: 10 changed, **7 at print precision**. All old values match the printed numbers.
- L889-891 / tab:magval:
  - τ95 resolve 0.308 → **0.313**;
  - τ99 0.366 → **0.372**;
  - +19% unchanged;
  - conditional F1 τ99 0.823 → **0.822**.
- L822: unchanged.
- **Depth 500 (L943-945, L1011-1013).**
  - The depth-100 side is 0.3071 → **0.3127** (conditional F1 0.8066 → **0.8063**).
  - The depth-500 side reads `two_panel/mag_panel_comparison.json`. That file is **R2.3's output** and is not
    regenerated yet.
  - The "ad hoc" part was re-derived with `tools/mag_depth500.py`, which reproduces the printed 49,185 / 88 / 87 /
    1 / 0 on the old tables. New: **49,795** scoreable cells per regime; still 88 outcome changes (87 resolve, 1
    gene-less reference) and 0 selection switches (`control/mag_depth500_{old,new}.json`).

---------------------------------------------------------------------------------------------------------------

## R2.1 E3 gate and review (commit `e722aaa`)

### Controls (all PASS)

**C-E3 (`r21-e3-control`).**
- Inputs:
  - the pre-repair cache (`f0a84bf1`) + the pre-repair catalogue (`46488fa0`), bound over `bvbrc_cache/`;
  - the cc8014f reference cells (`db5f6e1c` / `5a84dd24`), bound over the canonical paths.
- Outputs were bound to `$W/control_e3/out`, so the clone's tracked files were never written.
- Result: `results`, `per_region_n2+`, verdict, `verdict_detail` and diagnostics (timing masked) are
  **IDENTICAL** to the committed `pfunc_e3_gate.json`.
  - `_meta` store, reproduction (135,982 match), neighbours, keys and provider misses are identical too.
  - `per_amplicon.csv.gz`: 135,982 rows × 55 non-timing columns, **0 differences**. `run_meta` is identical
    apart from times and args.
  - Report: `control/e3/c_e3_compare.json`.

**Review control (`r21-e3-review-control`).**
- Inputs: the same pre-repair cache + catalogue, and the committed `per_amplicon.csv.gz`.
- Against the committed `review_recompute.json`: all counts exact. 41 float leaves differ at ≤ 1.2e-15 relative
  (summation order under `imap_unordered`); 0 differ at 1e-12.

**Determinism.**
- Two smokes (`--limit 500 --workers 8 --verify-q`) on the repaired inputs: identical in every block and in all
  1,000 rows × 57 columns (`control/e3/smoke_determinism.json`).
- verify-q max |q − independent q| = 3e-15.

### Run (`r21-e3-run`)

- Smokes: 92 / 91 s.
- **Gate: 944 s at 10 workers, cgroup peak 14.7 GB.** Review: 205 s.
- **Checks:**
  - `reproduction_vs_gene_capture_cells` match **135,982/135,982** against the R1.3/R1.5 cells;
  - `catalogue.is_current` true;
  - 0 provider misses;
  - run valid; exit 0.
- Pulled and installed on FM; `tf_verify` 4/4 identical. The review dir is at `$W/e3_review/`, and
  `review_recompute.json` equals `e3_review/indep_recompute.json`.
- Old outputs: FM `$W/pre_repair/region_validation_d100/data/pfunc_e3_gate*` (and in git).

### Old → new

**Verdict: PASS → PASS.**

exclude_0.987:
- Scored amplicons 59,571 → 61,047; `pred_empty` 1,060 → 243; `truth_empty` 953 → 278.
  - The plan expected about 675 more per mode. The actual +1,476 equals the R1.3 gene-capture outcome shift exactly.
- **n ≥ 2: 19,158 → 19,542 amplicons; 6,776 → 6,891 sources; 68.7M → 70.5M cells.**

| n ≥ 2 | k/n | E3 (species) | difference [95% CI] |
|---|---|---|---|
| Brier | 0.1547 → 0.1536 | 0.1017 → 0.1016 | **−0.0530 [−0.0547, −0.0513] → −0.0520 [−0.0536, −0.0505]** |
| AUROC | 0.8745 → 0.8781 | 0.9283 → 0.9282 | **+0.0538 [+0.0517, +0.0560] → +0.0501 [+0.0484, +0.0519]** |
| slope | 0.44 → 0.45 | 1.10 → 1.09 [1.07, 1.12] | |

- exclude_0.995, n ≥ 2: Brier difference −0.0428 → −0.0419; AUROC difference +0.0445 → +0.0419.
- n = 1, E3: Brier 0.139 (unchanged); AUROC 0.867 → 0.868.

**Review** (new `review_recompute.json`):

| quantity | old → new |
|---|---|
| recalibration share of the Brier gain | 44.5% → **47.9%** |
| saturated k/n Brier | 0.1311 → 0.1287 |
| saturated k/n AUROC | 0.8806 → 0.8840 |
| E3 vs saturated k/n | −0.0294 / +0.0477 → **−0.0271 / +0.0442** |
| Brier difference by n (n = 2 / 3 / ≥ 4) | −0.1007 / −0.0623 / −0.0240 → −0.0998 / −0.0617 / −0.0229 |
| re-entry | 13.5% → 13.6% of amplicons; 1.19% → 1.17% of weight |
| clean subset | 16,445 → 16,753 amplicons; −0.0597 → −0.0587 |

- Hard-leakage counts: still 0.
- Agreement with the gate CSV: 0 mismatches.

### PFUNC_E3_GATE.md

- All numbers were re-scored.
- A re-score paragraph was added: date, machine, inputs, controls, what moved.
- **The criterion paragraph is byte-identical.**

Caveat numbers that neither the gate nor the review stores now have a producer:
- **Which numbers.** The by-n CIs, n ≥ 10, the family/order clustering widening, and the 11,557 → 13,623
  no-lost-weight subset.
- **New script.** `region_validation_d100/scripts/pfunc_e3_gate_caveats.py` writes
  `data/pfunc_e3_gate/caveats_recompute.json`, using the gate's own bootstrap. Both files are in `e722aaa`; the
  JSON was added with `-f` because that directory is git-ignored except for named files.
- **Control on the 2026-11-12 inputs.**
  - Exact: the gate's pooled CI, the 11,557 subset, the n counts, every point estimate, and 1.9× for family.
  - Within 0.0003: the other CI bounds, compared with the reviewer's ad hoc ones.
  - The order widening came out 2.6× vs the printed 2.7×. That is recorded in the docstring.
- **New values:**
  - family widening 1.9× → 2.0×;
  - order 2.7× → 2.8×;
  - n ≥ 10: −0.0083 → −0.0076;
  - missing neighbours: 32,021 → 32,157 distinct genomes.

### Manuscript items (`collect.py --only` the 3 R2.1 items → `control/collect_run2_tf_r21.{json,md}`)

- 30 values: all changed, **26 at print precision**.
- Old values match the printed ones, except the known rounding of "about 45%" (unrounded 44.5%).

| item (L254-255; also L361, 1225) | old → new |
|---|---|
| Brier | 0.155 → **0.154**, and 0.102 |
| Brier difference | −0.053 [−0.055, −0.051] → **−0.052 [−0.054, −0.050]** |
| AUROC | 0.874 → **0.878**, and 0.928 |
| AUROC difference | +0.054 [+0.052, +0.056] → **+0.050 [+0.048, +0.052]** |
| amplicons | 19,158 → **19,542** |
| sources | 6,776 → **6,891** |
| recalibration share | about 45% → **about 48%** |
| E3 vs saturated k/n | 0.029 / 0.048 → **0.027 / 0.044** |

---------------------------------------------------------------------------------------------------------------

## TF timings and memory (cgroup `memory.peak` read inside the unit; GNU time max RSS)

| unit | wall | cgroup peak | max RSS | MemoryMax |
|---|---|---|---|---|
| r22-w1-verify | 212 s | 9.3 GiB | 3.3 GiB | 16G |
| r22-mag-control (121 batches, P6) | 569 s | 28.8 GiB | 4.9 GiB | 36G |
| r22-mag-rescore (120 batches + smoke + agg) | 652 s | 25.0 GiB | 4.9 GiB | 40G |
| r21-e3-control (gate, 10 workers) | 936 s | 20.9 GiB | 10.6 GiB | 32G |
| r21-e3-review-control (8 workers) | 203 s | 14.4 GiB | 10.7 GiB | 32G |
| r21-e3-run (2 smokes + gate + review) | 1,333 s | 14.6 GiB | 10.9 GiB | 32G |

- Never more than 52G of MemoryMax at once: w1-verify ‖ mag-control. All other units ran alone.
- The TF lane's wall time was 22:14–23:17, about 63 min.

## Files committed (prFBA, branch `hardware-aware-alignment-dispatch`, not pushed; path-limited)

- **`abdfca6`** (119 files):
  - `mag_validation/out_batches_{tau99,n500}/`, `mag_validation/smoke/out/`, `mag_validation_km/out_batches/`;
  - `mag_validation/out/{mag_cells_all,mag_cells_tau99}.csv`, `{mag_summary_full,mag_summary_tau99,summary_n500}.json`;
  - `mag_validation_km/out/summary_km.json`;
  - 3 new `*.stamp.json` index sidecars.
- **`e722aaa`** (7 files):
  - `region_validation_d100/data/pfunc_e3_gate.json`;
  - `data/pfunc_e3_gate/{per_amplicon.csv.gz,run_meta.json,review_recompute.json}`;
  - new `data/pfunc_e3_gate/caveats_recompute.json`;
  - new `scripts/pfunc_e3_gate_caveats.py`;
  - `PFUNC_E3_GATE.md`.
- Not touched: the other agent's `picrust2_comparison` changes, `presentations/`, KBase, and the ModelSEED,
  codiffusion and EmilyKin checkouts.

## Notes for later stages

1. **R2.3 (other agent).** Everything it reads from R2.2 is installed on FM.
   - Not regenerated here: `mag_validation/out/{paired_*,crossover*,af_gate_counterfactual}.json`, `F_mag_*`,
     `three_panel/*`, `two_panel/*`.
   - The depth-500 claim's d500 side comes from `two_panel/mag_panel_comparison.json`.
2. **R2.5 step 3 (P2C MAG arm).** `mag_prfba_ec.py` reads the new `pgfam_panel.json` and checks against the new
   `out_batches`/`lowc_out_batches`. Both are installed, so it can run.
3. **`$W/pre_repair` now also holds** the old MAG indices, the git-ignored MAG batch tables and logs, and the old
   gate outputs. `$W/from_tf/pull_R2.{1_e3gate,2_mag}` hold the pulled copies and can be deleted after review.
4. **New tools in `$W/tools`:** `tf/{unit_run_tf,mag_control,w1_verify,mag_rescore,e3_control,e3_review_control,e3_run}.sh`
   (copies on TF), `r22_install.sh`, `r22_w8.sh`, `r22_checks.py`, `mag_depth500.py`, `cmp_e3_full.py`,
   `cmp_json.py`, `e3_caveats.py`.
   - `tf_unit` + `unit_run_tf.sh` give the TF cgroup memory peak, which the journal does not.
5. **Runbook corrections.**
   - Create `mag_validation{,_km,_rnd}/logs` on TF before any driver.
   - Before `tf_code_sync.sh`, revert TF-local tracked outputs with `git checkout --`, after checking they equal the
     FM commit.
