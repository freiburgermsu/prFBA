# Run 1 report: cache writers and controls (R1.0–R1.6), 2026-11-15

Executor: Run 1 agent. Plan: `plan/PLAN.md` §0–4 (Run 1) and §8. Decisions D1–D11 as fixed by the orchestrator;
R1.K was not run here (D5: the KBase upload was a separate process, which committed its uploader as `7adf85c` at
17:09 and finished during this run).

**Outcome.** Every Run 1 stage ran and every control passed or was diagnosed:
- The six depth-100 trees are re-scored on the repaired cache: gene capture in all six, and selection, taxonomy
  and McNemar in the three legacy trees.
- The d20 original tree is re-selected and re-scored with today's code (D2), with the C-D20 control. The d20
  k-means and random gene capture is re-scored.
- E1 (0.995) is re-scored.
- No BV-BRC fetch happened anywhere, so the shared cache is byte-identical to its R1.0 state. The E3 catalogue
  stays current (no rebuild), and the cache is frozen at R1.6 (20:47:47).
- W8: every changed gene-capture cell in the 6 d100 trees and d20 k-means has a repaired source or rep.
  15 d20 cells are not the repair. All are explained by cache history from before 11-11 (none by code):
  - 14 d20-original cells on the four SNAP20J genomes;
  - 1 d20-random cell whose rep entered the cache on 10-21.
- Results are committed in three path-limited commits (not pushed).

Progress log with timestamps: `logs/RUN1_PROGRESS.md`. Every unit's memory and wall time: `stamps/units.jsonl`.

---------------------------------------------------------------------------------------------------------------

## Deviations from the plan text, and why

1. **Memory peaks.** The plan's `U()` relies on the journal's end-of-unit "Memory peak". On this machine
   (systemd 255) that value is wrong: it reported 256K for a 1.5 GB allocation, because it reads an emptied cgroup.
   - Every unit here runs through `tools/U.sh` + `tools/unit_run.sh`. The wrapper samples the unit's own cgroup
     every 30 s (`logs/<unit>.mem.log`) and records `memory.peak` from inside the cgroup, plus GNU time's max
     RSS, in `stamps/units.jsonl`. Self-test: 1.48 GiB recorded.
   - The cgroup peak includes page cache. The max RSS is the largest single process.
   - `U.sh` also adds `LimitNOFILE=1048576` (the driver warns below 65,536; the random panel once died on
     Errno 24) and `PYTHONUNBUFFERED=1`.
2. **R1.0 step 5 (uploader commit) and R1.K were skipped (D5).** The tag `rescore-2026-11-start` already existed (cc8014f).
3. **C-L4 used SNAP20J.** On SNAP20, all 17 substantive differences landed on the 4 SNAP20J genomes, which is
   the condition of handoff item 3. The SNAP20 result is kept as `control/d20_code_delta*_snap20.json`.
4. **C-D20 was run although the selection code delta turned out to be schema-only.** The plan allows dropping it in that case.
   The gap-close postdates the July gene capture, and the July cache state is not SNAP20J, so the control was
   needed for an exact three-point attribution. It cost 12 min.
5. **R1.4 additions.**
   - `prefetch_gene_families.py --select-mode legacy` on the d20 hits runs before `select_references`, as the d100
     driver does, because each lazy miss is a whole-cache rewrite. It found 28,639 cached and 0 to fetch.
   - The flock is taken with `flock -n`, so the run fails instead of waiting.
   - The plan's single-quoted `bash -c` block was written as a script, `tools/r14_d20.sh`.
6. **W8 needs reps dumps.** `tools/dump_genecap_reps.py` gained `--legacy`; without it, legacy trees would be
   re-selected under the union rule.
   - Every dump runs through `ctl_run` on a read-only cache copy, with `--keep-tables`, so each dump also
     re-proves its tree's gene-capture tables. Byte-for-byte reproduction:
     - old side: the 5 d100 trees run here, the original union (helper run), and d20 k-means;
     - new side: all 6 d100 trees;
     - `reps_d20.sh`: all 4 dumps (the control tree, and the new d20 original, k-means and random).
   - Two old-side dumps do NOT reproduce their committed tables:
     - d20 original on SNAP20 with the June selection: 35 cells, because the July cache was not SNAP20;
     - d20 random on SNAP: 1 cell, the committed `pgfam_missing` whose rep was cached on 10-21.
     Their reps are used only as attribution inputs.
7. **Tree figures.** The git-ignored `figures/` of the six d100 trees were moved to `pre_repair/`, so that the
   renderers, which skip existing files, re-rendered them. All 44 figures per tree were re-rendered.
8. **Incident (no effect).** `kmeans_panel/close_exclude_rep_gap.py --help` is not argparse: it started a
   gap-close on `region_validation_kmeans`.
   - It was killed while loading the selection and hits, before any cache access. The cache stamp was unchanged,
     there was no `.tmp`, and no tree file was written.
   - From then on, sources were read rather than probed.

New tools in `$W/tools`:
- `U.sh`, `unit_run.sh`
- `stamp_state.py`
- `cmp_selection.py` (record-level selection diff and the touched set)
- `r13_moveaside.sh`, `r13_checks.sh`, `r14_d20.sh`, `r15_e1.sh`, `r16_final.sh`
- `reps_new_d100.sh`, `reps_d20.sh`
- `genecap_headline.py`, `taxacc_headline.py`, `post_stage_mem.py`, `d20_delta_detail.py`

Other tool changes:
- `dump_genecap_reps.py` gained `--legacy`.
- `tools/test_tools.py` passes 9/9.

---------------------------------------------------------------------------------------------------------------

## R1.0 Preflight

- **Writer check:** `pgrep` showed only the KBase upload (allowed), and there was no `$CACHE.tmp`.
- **`stamps/R1.0.json`:**
  - cache: 5,118,271,515 B, 2026-11-13T10:38:50, sha256 `05f03bcd…`, key-set sha256 `4ef0950c…`, 159,753 entries;
  - catalogue `catalogue_is_current` true;
  - CORRECTIONS sha256 `429c7b7a…`;
  - prFBA `23773d3`, paper `0214a14`;
  - d20 baseline 32/32 OK.
  - The full key list and the 1,006 repaired values are saved alongside, for R1.6.
- **`stamps/ignored_before.sha256`:** 442 git-ignored files (selections, taxacc JSONL, genecap long, MAG batches,
  SR predictions).
- **`expected/`:** the 11-13 measurements were already present.

## R1.1 Snapshots

These were built in the helper run and re-verified here:
- SNAP `f0a84bf1…` (= the logged pre-repair file);
- SNAP20 `86cfdbe1…`;
- SNAP20J `e57a129c…`.

C-L ran on copies (`cl_snap_copy.json`, `d20ctl/snap20{,j}_copy.json`), never on the shared snapshots.

## R1.2 Controls

| control | result | evidence |
|---|---|---|
| **C-U** (D8) | **PASS** | `score_taxacc --force` + `mcnemar --force` in `region_validation_d100`: git diff clean on all 9 tracked outputs; `taxacc_per_record.jsonl` sha256 `091cc592…` == R1.0. 263 s, 22.3 GiB |
| **C-L1** (gate) | **PASS** | `select_references --select-mode legacy --gene-cache <SNAP copy>` on the d100 legacy hits: 67,991/67,991 records identical, `_meta` identical with `generated_from` masked, 0 lazy misses (copy stat unchanged). Touched records 20,644 (= plan). `control/cl1_*.json` |
| **C-L2** | **PASS** | k-means legacy 74,704/74,704 identical; touched 13,265 (= plan) |
| **C-L3** | **PASS** | random legacy 75,328/75,328 identical; touched 12,994 (= plan) |
| **C-L4** (D2) | code delta = **output schema only** | On SNAP20: 17 records differ in substance, all `no_gene_set → selected` on the 4 SNAP20J genomes. On SNAP20J: 0 substantive differences, `n_selected` identical in 67,991/67,991, `_meta.summary` identical. 67,427 records differ only by fields added since June (`selected[].equivalent_genomes`, `selected[].gene_set_empty`, `provenance.n_equivalent_genomes`), and knobs gain `select_all_*` and `select_mode`. `control/d20_code_delta.json`, `_substantive.json` |
| **C-G** | PASS (helper run) | k-means/random d100 union `score_genecap` on SNAP: byte-identical. Re-confirmed here: the old-side reps dumps of both trees reproduce cells, summary and long table byte for byte |
| **C-D20** (D2) | done | Today's code on SNAP20J in a control tree (`control/region_validation_ctl`), every step via `ctl_run`; gap-close 0 missing; snapshot copy unchanged. **Committed → control:** taxacc tables and `taxacc_per_record.jsonl` byte-identical (summary `_meta` paths differ); McNemar pairs, self-recovery diagnosis and failures identical (`mcnemar_regions.json`: paths and `correct_definition` text); gene capture 17/135,982 cells differ (exclude_0.987 only: reps absent from the July cache, filled by the October gap-close). 699 s, 42.6 GiB |

Old-side reps dumps (W8 inputs), checked against the committed tables:
- **Byte-identical reproduction:** 5 d100 trees, plus d20 k-means.
- **d20 random:** differs in 1 cell. It is committed `pgfam_missing`; the rep was added to the cache after 10-18.
- **d20 original (SNAP20 + June selection):** differs in 35 cells, so the July cache was not SNAP20.

## R1.3 Six-tree d100 POST

Pre-launch:
- MemAvailable was 209 GiB (≥ 170 required).
- No other job above 10 GB was running beside the POST: only the user's codif units, and the KBase unit had
  already ended. The `cmp_selection` units run beside it peaked below 1 GiB.
- Move-aside: 66 items to `pre_repair/<repo path>` (`logs/r13_moveaside.log`).

The run:
- Unit `rv-rescore-post`, MemoryMax=150G: rc 0, **8,844 s (2 h 27 min)**.
- Cgroup peak **112.0 GiB**, including page cache, which builds up tree by tree. Largest process **62.5 GiB RSS**.
- Per tree and per stage: `control/r13/post_stage_mem.json`.

| tree | wall | cgroup max | heaviest stage |
|---|---|---|---|
| original union | 16.7 min | 61.2 GiB | P7 |
| original legacy | 31.8 | 72.6 | P7 72.6, P5 70.3, P2 67.1 |
| k-means union | 17.1 | 85.1 | P7 |
| k-means legacy | 32.2 | 95.1 | P7 95.1, P5 90.8, P2 87.2 |
| random union | 16.8 | 99.5 | P7 |
| random legacy | 32.6 | 111.0 | P7 111.0, P5 107.9, P2 104.3 |

Checks (`logs/r13_checks.log`, all PASS):
- **Check 5:** passes in all six trees (`select_mode`, `pgfam_missing == 0`, gap close).
- **Fetches:** 0 everywhere (every P1/P3/P4 `to_fetch=0`; gap-close `0 missing`). The post-R1.3 cache sha256 equals R1.0's.
- **Predicted output:** the original union `gene_capture_{cells,by_region,long}` are **byte-identical** to `expected/repro_0987`.
- **Union trees:** `selection.json` and `taxacc_per_record.jsonl` sha256 equal R1.0; git shows only the 2 gene-capture files changed.
- **Legacy selections:** every changed record is a touched record.
  - Original: 1,506 of 20,644.
  - k-means: 1,103 of 13,265.
  - random: 1,822 of 12,994.
  - The main stage moves are `no_gene_set → marginal_gain / included / tier_cap`.
- **W8 (6/6 PASS, 0 unattributed, 0 one-sided, 0 identity changes):**

| tree | cells | changed | rep_repaired | source_repaired |
|---|---|---|---|---|
| original union | 135,982 | 3,732 | 2,122 | 1,610 (= 11-13 measurement) |
| original legacy | 135,982 | 3,866 | 2,256 | 1,610 |
| k-means union | 149,408 | 1,996 | 1,744 | 252 |
| k-means legacy | 149,408 | 1,730 | 1,478 | 252 |
| random union | 150,656 | 2,502 | 1,998 | 504 |
| random legacy | 150,656 | 2,024 | 1,520 | 504 |

## R1.4 d20 trees (D2)

- **Unit and lock:** `rescore-r14-d20` under `flock -n` on the cache-writer lock; 1,721 s; cgroup 48.3 GiB, RSS 43.3 GiB.
- **Fetches:** candidate prefetch 0, lazy misses 0, gap-close 0 missing. The cache stat was unchanged throughout.
- **d20 union trees:** `selection.json` and `taxacc_per_record.jsonl` sha256 equal R1.0; only the 2 gene-capture files changed.
- **Selection, control → new (data step, `--d20`):** 1,417 records changed.
  - 1,411 are touched by repaired or poisoned genomes.
  - The other 6 all fall on the SNAP20J reset genomes (33953.3 ×5, 176233.20 ×1, `no_gene_set → selected`). These 4
    are neither repaired nor poisoned: they hold PGFams in the pre-repair cache. So the 6 are June→11-11 cache
    history, carried by the control's June emulation.
  - Mean `n_selected` 1.298275 → 1.309217; abstained 564 → 167; `no_gene_set` exclusions 4,278 → 828.
- **Self-recovery diagnosis (d20):** failures 11,300 → 10,821; terminal stage `no_gene_set` 974 → 262.
- **Three-point attribution (committed → control → new):**
  - taxonomy: identical → changes (include mode only);
  - gene capture: 17 cells (code and cache history) → see W8 below.
- **W8:** see the "d20 W8" section below.

## d20 W8 (gene capture; `tools/reps_d20.sh`; 1,516 s, RSS 43.4 GiB)

Reps dumps reproduce their tables byte for byte in 4 of 4 cases: the control tree on SNAP20J, and the new d20
original, k-means and random trees on a copy of the frozen cache.

| tree / step | cells | changed | rep_repaired | source_repaired | not repair | verdict |
|---|---|---|---|---|---|---|
| d20 original, committed → control (code and cache history) | 135,982 | 17 | – | – | 17: exclude_0.987 reps absent from the July cache, filled by the October gap-close | code step (as designed) |
| d20 original, control → new (data) | 135,982 | 3,549 | 1,925 | 1,610 | **14, all on the SNAP20J genomes**: source 33953.3 ×10 (truth `[]` in the June emulation), rep 294671.3 ×3, rep 176233.20 ×1 | PASS after explanation |
| d20 k-means, committed → new | 149,408 | 1,494 | 1,242 | 252 | 0 | **PASS** |
| d20 random, committed → new | 150,656 | 1,895 | 1,390 | 504 | **1**: 2026724.2804 V1-V3 exclude, committed `pgfam_missing`; rep 2026734.1643 was first cached by the d100 random gap-close on 10-21 (cache key #157,329) | PASS after explanation |

- The four SNAP20J genomes are not among the 1,006 repaired genomes or the 20 poisoned keys. They hold PGFams
  in the 11-11 cache, which is how today's d20 selection picks them.
- Their cells are therefore June→11-11 cache history, introduced by the control's June emulation. They are not
  the 11-13 repair and not code.
- Per-cell evidence: `control/w8/d20_unattributed_explained.json` (15 explained, 0 unexplained).
- The informative committed → new run (`attrib_region_validation_committed_to_new.json`) has 18 non-repair
  cells: the same 14, plus 4 that lie inside the committed → control step.


## R1.5 E1 (0.995)

- Unit `e1-rescore`: smoke 600/600 rows matched; full run 615 s (unit 703 s); RSS 9.3 GiB.
- 0 provider misses and an empty supplement, so there was no `--merge-supplement` and no cache write.
- **Built-in 0.987 reproduction:** byte-identical against the new P7 files (3/3).
- **`gene_capture_tau_compare.json`:** equals `genecap_repaired_impact/gene_capture_tau_compare.json` in content.
  The only difference is `_meta.repro_0987_byte_identical`, which was false at 11-13 (before P7) and is true now.
- The tau0995 cells, summary and long table are byte-identical to the 11-13 prediction.

## R1.6 Cache-final checkpoint

- `stamps/R1.6.json`: cache UNCHANGED versus R1.0 (size, mtime, sha256, key set); all 1,006 repaired values
  unchanged; catalogue current.
- No catalogue rebuild was needed. The deployments' catalogue stamp stays current.
- **FREEZE at 2026-11-15T20:47:47.**

---------------------------------------------------------------------------------------------------------------

## Old → new headline numbers

### d100 gene capture: 8-region mean macro recall %, resolved (OK cells) / population (abstain → 0)

| tree | mode | resolved | population | outcome changes |
|---|---|---|---|---|
| original union | include | 99.460 → 99.466 | 99.420 → 99.448 | ok 67,010→67,700; pred_empty 28→13; truth_empty 953→278 |
| original union | excl 0.987 | 67.793 → 67.783 | **60.433 → 61.287** | ok 59,571→61,047; abstained 6,407→6,423; pred_empty 1,060→243 |
| k-means union | include | 96.731 → 96.743 | 96.560 → 96.639 | ok 74,255→74,347; pred_empty 131→80; truth_empty 318→277 |
| k-means union | excl 0.987 | 67.099 → 67.111 | 46.427 → 46.735 | ok 51,262→51,625; abstained 22,527→22,528; pred_empty 597→274 |
| random union | include | 96.564 → 96.580 | 96.427 → 96.495 | ok 74,710→74,932; pred_empty 104→65; truth_empty 514→331 |
| random union | excl 0.987 | 67.792 → 67.824 | 39.730 → 40.006 | ok 43,594→43,986; abstained 30,804→30,829; pred_empty 416→182 |
| original legacy | include | 96.939 → 96.936 | 96.874 → 96.926 | ok 66,991→67,706; abstained 47→7 |
| original legacy | excl 0.987 | 70.777 → 70.798 | 63.396 → 64.106 | ok 59,863→61,137; abstained 7,175→6,576 |
| k-means legacy | include | 94.442 → 94.465 | 94.298 → 94.391 | ok 74,272→74,369; abstained 114→58 |
| k-means legacy | excl 0.987 | 70.006 → 70.023 | 48.622 → 48.859 | ok 51,459→51,729; abstained 22,927→22,698 |
| random legacy | include | 93.021 → 93.034 | 92.915 → 92.989 | ok 74,730→74,962; abstained 84→35 |
| random legacy | excl 0.987 | 70.847 → 70.894 | 41.690 → 41.889 | ok 43,773→44,063; abstained 31,041→30,934 |

### d20 gene capture

| tree | mode | resolved | population |
|---|---|---|---|
| original (legacy) | include | 96.918 → 96.913 | 96.842 → 96.903 |
| original (legacy) | excl 0.987 | 69.785 → 69.783 | 52.968 → 53.707 |
| k-means | include | 96.293 → 96.310 | 96.123 → 96.206 |
| k-means | excl 0.987 | 65.664 → 65.683 | 33.235 → 33.516 |
| random | include | 95.758 → 95.780 | 95.622 → 95.695 |
| random | excl 0.987 | 66.014 → 66.039 | 27.581 → 27.808 |

- Original: truth_empty 958 → 278 and abstained (excl) 16,538 → 15,980.
- Random: `pgfam_missing` 1 → 0.
- Control → new (original, excl 0.987): population 52.974 → 53.707.

### Legacy taxonomy

Anchor predictor, 8-region mean % `correct_population`:
- Include-mode Genus:
  - d100 original 91.914 → 92.488
  - k-means 92.307 → 92.647
  - random 91.151 → 91.830
  - d20 original 91.915 → 92.512
- Include Species moves by +0.2 to +0.6 points.
- **Exclude-mode taxonomy is unchanged in every legacy tree.** Its re-selection reads no gene provider.
  `control/taxacc_headline_run1.json`.

### Legacy selections

| tree | mean `n_selected` | abstained | `no_gene_set` exclusions |
|---|---|---|---|
| d100 original | 1.307129 → 1.317630 | 549 → 167 | 5,098 → 1,053 |
| k-means | 1.426122 → 1.428424 | 211 → 129 | 4,718 → 1,728 |
| random | 1.469878 → 1.473038 | 284 → 129 | 6,693 → 1,922 |

**Caution on the definition of `n_selected`.** These means are over all selection records. The plan's "legacy
include-mode `n_selected` mean was 1.3176 before" is a different quantity: it is close to the mean over
non-abstained records (1.317769 before → 1.320875 after for the d100 original). The manuscript's "1.32 for the
ordered reducer" (collect.py `select_legacy_mean_size`: mean over include-mode `ok` gene-capture cells) is
unchanged at print precision. Do not read the new all-records mean, 1.317630, as the plan's "before" value.

### McNemar (the 4 regenerated tables: 3 d100 legacy trees and the d20 original)

- In every table, 112 of the 224 pair-tests changed, **all of them include-mode**. Exclude-mode McNemar is
  unchanged, consistent with the unchanged exclude-mode taxonomy.
- Holm-significant region pairs, anchor|Genus include, out of 28:

  | tree | old → new |
  |---|---|
  | d100 original legacy | 25 → 28 |
  | k-means legacy | 24 → 25 |
  | random legacy | 22 → 24 |
  | d20 original | 25 → 26 |

- Anchor|Genus exclude_0.987 is unchanged: 24, 27, 26 and 28 of 28.

### E1 (0.995)

Population 0.987 → 0.995, by region:
- FL 64.0 → 69.2
- V4 59.3 → 64.0

Across regions:
- 0.995 population range 62.5–69.0 → **63.4–69.7**
- 0.987 range 57.6–64.4 → **58.6–65.2**

### Manuscript numbers (`manuscript_numbers/collect.py --only` the 21 Run 1 items)

Outputs: `control/collect_run1.{json,md}`; the repo's `number_diff.*` is untouched.
- 101 values in total; all old values match the printed numbers; 0 errors.
- 92 changed, **38 of them at print precision**. Among those:
  - **R1.3U:**
    - include capture min 99.1 → 99.2
    - exclusion population range 57.6–64.4 → 58.6–65.2
    - V3-onward 57.6–59.2 → 58.6–60.3
    - FL/V1 63.3–64.4 → 64.0–65.2
    - V4 58.1 → 59.3, FL 63.3 → 64.0
    - leakage gap V4-V5 32.9 → 33.0
    - share with `n_selected`>1: 34.0 → 34.1
    - abstract 58–64 → 59–65
  - **R1.5:**
    - tab:tau population 0.995: 62.5–69.0 → 63.4–69.7
    - score-zero at 0.987: 11.1 → 9.9
    - `pred_empty`: 1.6 → 0.4 (0.987) and 1.4 → 0.4 (0.995)
    - V4 0.995: 63.1 → 64.0; FL 0.995: 68.7 → 69.2
    - population gain max 5.5 → 5.4
  - **R1.3L:** union V4 include precision 0.795 → 0.794 (the legacy 0.856 is unchanged).
  - **R1.4:** d20-vs-d100 conditional F1 maximum drop 1.2 → 1.1.
- The legacy self-recovery items (12,602 / 67,466, about 19%) read the d100-legacy `self_recovery_diagnosis.json`.
  That file is regenerated by Run 2 D5, so they are not final yet.

---------------------------------------------------------------------------------------------------------------

## Files committed (prFBA, branch `hardware-aware-alignment-dispatch`; not pushed)

| commit | stage | content |
|---|---|---|
| `81e9b34` | R1.3 | 33 files: the 6 d100 trees' `gene_capture_{by_region.json,cells.csv}`, and the 3 legacy trees' `taxacc_region_summary.json`, `{asv_region_concordance,concordance_long,confusion_genus,self_recovery,mcnemar_region_pairs}.csv`, `mcnemar_regions.json` |
| `000fcaa` | R1.5 | `region_validation_d100/data/{gene_capture_by_region_tau0995.json,gene_capture_cells_tau0995.csv,gene_capture_tau_compare.json,genecap_tau0995/run_meta.json}` |
| `3ffc97d` | R1.4 | 19 files: `region_validation/data/` taxacc, gene capture (incl. `gene_capture_long.csv`), McNemar, self-recovery (6 + 6 files); `region_validation_{kmeans,random}/data/gene_capture_{by_region.json,cells.csv}`; new `region_validation/data/control_currentcode_precache/{selection_meta,taxacc_region_summary,gene_capture_by_region}.json` (C-D20, D2) |

- All three commits sit on top of `7adf85c`, the KBase lane's commit.
- No code fixes were needed in prFBA.
- The work-dir tools are outside the repo.
- Left untouched: `presentations/prFBA_updates_2026-10.user-edits-1029.pptx`, `upload_kbase_genomes.py`
  and the KBase files, and the ModelSEED, codiffusion and EmilyKin checkouts.

## What Run 2 needs to know (blockers and notes)

1. **Cache frozen at R1.6, byte-identical to R1.0** (`05f03bcd…`). The catalogue is current and unchanged.
   - R2.1 (E3 gate) can use the existing catalogue.
   - Its `reproduction_vs_gene_capture_cells` must now match the new cells (hazard H7).
   - Start each Run 2 stage with `stamp_state.py --compare $W/stamps/R1.6.json`.
2. **All six `post.done` are back**, so the gates of `run_picrust2_d100.sh` and `run_downstream_d100.sh` are open.
3. **The d100 legacy `self_recovery_diagnosis.json` and `self_recovery_failures.csv` are stale** against the new
   legacy selection until R2.3 D5. `diagnose_self_recovery` is in `run_downstream_d100.sh`.
4. **Handoff item 6 is resolved.** The W1 `.new` panel indices record cache sha256 `05f03bcd…`, which equals the
   R1.6 stamp, so no rebuild is needed for R2.2.
5. **SNAP20J, not SNAP20, is the d20 pre-repair state that reproduces the June selection.** Use it for the R2.5
   `prfba_selection` control (handoff item 3).
6. **The d20 baseline now fails check 6 by design.** Ten of the 32 files changed (`region_validation/data/{taxacc_region_summary,gene_capture_by_region}.json`,
   `{asv_region_concordance,gene_capture_cells,mcnemar_region_pairs}.csv`, `selection.json`;
   `region_validation_{kmeans,random}/data/gene_capture_{by_region.json,cells.csv}`).
   - Re-stamp it in R5.5.
   - Until then use `check_d100.py post`, not `all`.
7. **Not regenerated, and they need a decision or a later stage:**
   - **The 44 tracked d20 tree figures** (`region_validation/figures/`). They were last rendered in June, which already
     predates their July/September data. Re-render them in R5.2 with `--force` and `SOURCE_DATE_EPOCH`.
   - **The June d20 variant-selection files** (`*_et/_etg/_gd/_gr/_ib/_ir`). They are a historical variant study.
8. **D10 / R2.5:** the d20 as-shipped chain reads the d20 `selection.json`, which is now today's selection on the
   repaired cache. The old one is in `pre_repair/region_validation/data/selection.json`.
9. **Disk:** `$W/control` holds about 35 GB of cache and snapshot copies (`cl_snap_copy`, `cache_post_r13`,
   `cache_post_r14`, `d20ctl/*`). They can be deleted after Run 2.
10. **KBase (D5):** the other process uploaded both builds and committed `7adf85c`. The R1.K post-upload items
    (README §11 is in that commit; deck constants; `bvbrc_rebuild_summary.json`) belong to that lane or to R5.3.
