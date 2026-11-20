# Run 2 report, local lane (R2.4, R2.5, R2.3), 2026-11-15/16

Executor: Run 2 local agent. Plan: `plan/PLAN.md` §0–4 (Run 2) and §8. R2.1 (E3 gate) and R2.2 (MAG re-score)
ran on the TF machine (other agent; `logs/R2.2_DONE` 22:37, commits `abdfca6`, `e722aaa`).
Progress log with timestamps: `logs/RUN2_LOCAL_PROGRESS.md`. Units: `stamps/units.jsonl`.

**Outcome.** Every local stage ran, every control passed or is diagnosed to a documented cause, and W8 attributes
every changed row of every regenerated cell table. The cache was never written (stamps `R2L.0`, `R2L.1`, `R2L.2`,
`R2L.final`: unchanged vs `R1.6`; every selection `_meta` has 0 reps/sources missing from the cache).
Five prFBA commits, not pushed.

---------------------------------------------------------------------------------------------------------------

## Deviations from the plan text, and why

1. **Code: the fetcher (commit `772832a`, before any fetch).** `fetch_bvbrc_annotations.py` (both byte-identical
   copies) now sends `prefetch_gene_families`' User-Agent, spaces requests by `--min-interval` (default 1 s), and
   takes `--selection` as a comma-separated list that defaults to **both** selection files, as
   `genome_function_sets.py` projects both. So D10 cannot recur when the d20 README pipeline is run as written.
   The runs still fetched once per file, as the plan and the d100 driver do. Tests: `tests/test_fetch_bvbrc_annotations.py`.
2. **Code: the AF-gate levels (commit `7b14dc3`).** `af_gate_counterfactual.py` hard-coded the three FL 16S levels
   (0.6912 / 0.71 / 0.7092). They moved (0.692 / 0.7116 / 0.7097), so the script now reads them from the paired
   comparisons. **Control:** run on the committed pairings with the pre-repair panel index, it reproduces the committed JSON
   exactly (`mags` and `summary` equal, plus a new provenance block). Test: `tests/test_af_gate_levels.py`.
3. **Selection controls needed a cache-history emulation.** On the full snapshots, the union selections differ from
   the committed ones by P3 tie-break swaps to genomes cached *after* the committed run. Cache key order is
   insertion order, so `tools/truncate_snapshot.py` cuts a snapshot to the cache size the run's log reports
   (148,300 at 10-20; 111,789 at 09-17). On the truncated snapshots both modes are **record-identical** (below).
4. **The d20 control state is SNAP20, not SNAP20J.** SNAP20J resets `33953.3` to `[]` (its June state), but the
   09-17 legacy selection selected it. SNAP20J[:111,789] therefore fails the legacy arm (5 records, all `33953.3`),
   while SNAP20[:111,789] passes both arms. SNAP20J remains the right state for the June `select_references` run (Run 1).
5. **Worker counts and memory caps.**
   - **d100 `prfba_selection`:** 12 workers hit a cgroup OOM at 90G (`rescore-r24-p2d`, 22:29; anon 95.7 GB, from
     fork copy-on-write; nothing was written). The re-run used `SEL_WORKERS=8` under 150G and peaked at 117 GiB.
   - **d20:** 12 workers instead of the plan's 24.
   - **`mag_prfba_ec.py`:** `Pool(40)` is hard-coded and OOMed at 60G after `mag_direct_feed.py` had finished
     (`rescore-r25-magarm`, 23:24). The re-run under 180G peaked at 107.5 GiB.
   - Results do not depend on any worker count; the selection controls compare content, not key order.
6. **Sequencing.** The P2C MAG arm (R2.5 step 3) ran after `R2.2_DONE`: `mag_prfba_ec.py` reads the R2.2
   `pgfam_panel.json` and `out_batches`. Its installed sha256 `cd5c3f49…` equals its stamp, which records cache `05f03bcd…`.
7. **W8 extensions (outside the repo).**
   - `tools/attrib_diff.py` gained `--extra-genomes` / `--extra-keys`: documented non-repair reasons. The Run 1
     version is kept as `attrib_diff.py.run1`; `test_tools.py` 9/9.
   - New `tools/attrib_rows.py` covers the row tables `attrib_diff` has no profile for.
   - New `tools/cmp_prfba_selection.py`, `gfs_diff.py`, `as_headline.py`.
   - New scripts `r2_*.sh`, `r24_p2d.sh`, `r25_*.sh`, `r23a.sh`, `r23b.sh`; commit helper `pcommit.sh`.
8. **The function-map pin control was not re-run.** It is W5 (`logs/w5_control.log`, EQUAL CONTENT for both pins).
   Both products files were unchanged since: sha256 recorded before the fetch in `stamps/R2L.products_before.sha256`.

---------------------------------------------------------------------------------------------------------------

## R2.4: d100 as-shipped PICRUSt2 (`rescore-r24-p2d-b`: 1,155 s, cgroup 117.1 GiB, RSS 48.7 GiB)

**Controls** (`control/r2/`):

| control | result |
|---|---|
| fmaps pin e8446e3e | W5: EQUAL CONTENT |
| `prfba_selection --legacy` on SNAP | **records identical** (67,991) |
| `prfba_selection` (default) on SNAP | 305 records differ, reps only. Every one swaps an early-cached genome out and a genome cached after the 10-20 run in (`cmp_d100.json`) |
| default on SNAP[:148,300] (`snap_d100_1020.json`) | **records identical**, `_meta` equal |

**Checks** (`control/r2/d100/`):
- **Inputs unchanged:**
  - `fetch genomes` to_fetch 0, and `bvbrc_source_genomes.json` is byte-equal;
  - `domain_routing.csv` is unchanged.
- **Fetch:** 39,229 + 2,183 new products.
- **Maps:** +39 `pgfam_ec` and +40 `pgfam_rxn` entries, 0 changed; `ec_rxn` identical.
- **EC vocabulary:** unchanged (2,277 BV-BRC; 2,132 shared).
- **GFS:** 729 changed, all repaired; 0 unexplained; 267 new and 160 dropped genomes (reps entering or leaving).
- **Selections old → new:**
  - legacy: 2,503 changed, all touched;
  - default: 755 changed = 305 cache history + 473 changed SNAP→new, 0 untouched.
- **W8 `comparison_cells.csv.gz`:** **PASS**. 24,827 changed = 11,270 source_repaired + 13,176 rep_repaired + 381
  cache_history_1020.
  - 0 identity-column changes; PICRUSt2 rows change only through their (repaired) source's truth.
- **F13:** `manuscript_figures/p2/render_comparison_d100.py` → `manuscript_figures/out/F13_prfba_vs_picrust2.{pdf,png}`.
  - Not copied into the paper (R5).
  - Suptitle n 67,038 → 67,713. Median ΔF1 0.0716–0.0951 → 0.0727–0.0952; the label "+0.07 to +0.10" is unchanged.
  - The stock figure was rendered by the driver.

**Old → new (EC, all regions):**
- scoreable n 67,038 → **67,713**
- prfba_exclude:
  - coverage 0.889 → **0.902**
  - conditional F1 0.8017 → 0.8016 (0.802)
  - population F1 0.713 → **0.723**
- PICRUSt2 conditional F1 0.7290 → 0.7293
- no EC prediction 11.1% → **9.8%** (abstain 9.6 → 9.5; gene-less reps 1.5 → 0.4)
- dialect amplicons 19,686 → 19,843; conditional ΔF1 vs BV-BRC truth +0.046 → +0.047

## R2.5: d20 as-shipped chain, ceiling, dialect, depth probe, MAG arm

Units:
- `rescore-r25-p2c`: 1,071 s, cgroup 66.5 GiB.
- MAG arm: `rescore-r25-magarm` (direct feed), then `-b` (239 s, 107.5 GiB).
- `rescore-r25-checks`: 363 s.

**Controls** (`control/r2/`):

| control | result |
|---|---|
| fmaps pin ac619c6e | W5: EQUAL CONTENT |
| default / legacy on SNAP20J | 3 / 20 records differ (late-cached genomes and `33953.3`) |
| default / legacy on **SNAP20[:111,789]** | **records identical** (both). `_meta.n_reps_missing` is 0 → 23 only because the 20 poisoned keys, and 3 others re-inserted later, sit past the cut; an absent key and `[]` are equivalent for the cache-only provider |
| default / legacy on SNAP20J[:111,789] | identical / 5 records (`33953.3`) |
| depth-probe sample keys | **identical** to committed, so the k500 hits are reused (no GPU) |

**Three-point split** (`d20/cmp_*`):
- committed → SNAP20 (cache history): default 11 / legacy 35 records;
- SNAP20 → new (data): default 367 / legacy 2,301 records, **0 untouched** in either.

**Checks:**
- **Fetch:** 40,895 (default) + **350,476 (legacy, D10)** new products.
- **Maps:** +305 `pgfam_ec` and +342 `pgfam_rxn` entries, 0 changed.
- **EC vocabulary:** BV-BRC 2,274 → 2,275 (+5.1.1.11); shared unchanged at 2,132.
- **GFS:** 801 changed = 709 repaired + 92 whose newly fetched products add ECs (D10); 0 unexplained.
  - `n_with_empty_ec` 876 → 242.
- **W8 `comparison_cells.csv.gz`:** **PASS**. 22,555 changed = 11,270 source_repaired + 11,048 rep_repaired + 202
  rep_newproducts + 35 cache_history_0917.
  - The D10 cells are all `prfba_*_legacy`: 151 exclude, 51 include.
  - PICRUSt2 changes only through a repaired source (1,610 cells).
- **W8 rows** (`w8_rows.json`, `w8_rows_magarm.json`, `d20/depth_probe_attrib.json`): all PASS.
  - `annotation_ceiling` per-genome rows: 31 changed, all repaired.
  - Depth probe: 63 changed rows = 17 repaired source or rep + 4 newly mapped products + 42 others. The 42 rows are
    identical between the SNAP20 re-evaluation (`rescore-r25-dpctl`, current code and maps on the pre-repair cache) and
    the new run, so their change lies entirely in committed (09-22) → SNAP20. That step is cache values differing from
    the 09-22 state; each of the 42 has a rep or source past the 09-17 key cut (insertion index >= 111,789; the 09-22
    size is unknown). Every SNAP20 → new change touches a repaired genome.
  - MAG arm selections: 279,916 rows, 0 changed (all match the R2.2 scorer's recorded anchor and size).
  - MAG arm cells: 11,676 source_repaired + 557 pick_repaired + 75,349 `base_constant` rows. The
    panel-prevalence floor is recomputed from all panel truths, which include the 140 repaired genomes.
- `mag_direct_feed_summary.json` regenerated byte-identical.

**Old → new:**
- prfba_exclude:
  - conditional F1 0.7985 → 0.7980 (0.798)
  - population F1 **0.600 → 0.609**
  - coverage 0.751 → 0.763
- annotation ceiling:
  - n **2,689 → 2,719**
  - shared recall/precision 0.904/0.743 unchanged; Jaccard 0.688 → 0.687
- depth probe:
  - recovery 283/320 unchanged
  - recovered EC F1 **0.783 → 0.829**, precision 0.747 → 0.791, recall 0.854 → 0.904
- d20 population gap vs PICRUSt2 −0.129 → −0.120

## R2.3: downstream (`rv-rescore-downstream`: 1,241 s, 36.5 GiB; `rescore-r23b`: 312 s)

**Controls:**
- **Taxonomy `est` unchanged:** in `three_panel_comparison_d100.json` every taxonomy crude, phylum-std and delta `est`
  is unchanged (96/96 per regime). The CIs move.
  - 5 `excludes_zero` verdicts flip in `_d100` and 16 in `_d100_legacy` (14 of them include mode, where legacy
    taxonomy moved in R1.3). All are borderline; listed in `control/r2/r23/three_panel_ci_flips.json` for R5.
- **Hash-seed proof:** `compare_three_panels` is byte-identical under PYTHONHASHSEED 0 and 777, for `_d100` and `_d100_legacy`.
- **Union D5:** `region_validation_d100/data/self_recovery_*` show no diff against HEAD.
- **check_d100 all:** only check 6 fails. A fresh `sha256sum -c` after R2.3 part B (`stamps/R2L.final.json`) gives
  17 of 32 failing, every one a re-score output:
  - R1.4's 10 d20 tree files;
  - R2.5's 5 P2C files (comparison cells and summary, both selections, GFS);
  - R2.3's `three_panel/three_panel_comparison.json` and `three_panel_mag_panelarm.json`.

  `space_deciles_three.json`, `three_panel_mag_gtdb.json`, the hits and the k-means/random selections and taxacc still
  pass. Re-stamp in R5.5.

**Checks:**
- **Emptied strata:**
  - original N 7,215 → 7,268;
  - emptied 336 and 2,019 unchanged;
  - random 2,926 → 2,927, from `129817.44`, a repaired genome that is newly scoreable. Emptiness reads only hits; the
    population moved.
- **Coverage:** original FL gene-capture coverage 0.921 → **0.931**, V4 0.863 → **0.881** (as predicted).
- **W8 rows:**
  - legacy `self_recovery_failures.csv`: 836 changed rows, all on records whose legacy selection changed in R1.3;
  - AF-gate per-MAG outcomes: 207 changed, all with a repaired source or skani reference.

**Old → new:**
- legacy self-recovery failures 12,602 → **12,127** of 67,466 (19% → 18%)
- MAG vs FL 16S at d100, 100% complete:
  - n 1,916 → 1,972
  - MAG 0.813 → 0.814, 16S 0.710 → 0.712, Δ +0.103 → +0.102
  - V4 Δ +0.212 → +0.209
- AF gate:
  - newly admitted 1,367 → 1,391; mean F1 0.599 → 0.600; pooled 0.742 → 0.743
  - share above the union FL level 26.7% → 26.5%; coverage gain +45.2% → +44.6%
- new `paired_tau99_d100.json`: FL conditional Δ **+0.127**; MAG population F1 0.304 (printed 0.301) against the
  re-scored 16S population 0.612 (printed "unchanged 0.605")

## Manuscript numbers (`collect.py --only`, outputs in `control/`, repo `number_diff.*` untouched)

- **`collect_run2_local.{json,md}`** (75 items: R2.3 42, R2.4 22, R2.5 11):
  - 761 values; every old value matches the printed one; 0 errors;
  - 439 changed at print precision;
  - 56 pending, all W6 derived numbers (Run 5).
- **`collect_run2_selfrec`** (the two R1.3L legacy self-recovery items D5 feeds): 12,602 → 12,127; 19% → 18%.
- **`collect_r24`** (R2.4 only, run before R2.5) is superseded by the above.
- **`collect_run2_r22.{json,md}`** covers the four R2.2 items, whose depth-500 side reads the R2.3
  `two_panel/mag_panel_comparison.json`. 16 values, 0 pending, 11 changed at print precision:
  - resolve τ95 0.307 → 0.313, τ99 0.366 → 0.372; conditional F1 τ99 0.823 → 0.822;
  - pooled resolve, depth 100 / 500: 0.3071 / 0.3089 → 0.3127 / 0.3144;
  - conditional F1, depth 100 / 500: 0.8066 / 0.8061 → 0.8063 / 0.8059.
- **`mag_panel_comparison.json` at 4 dp** (`control/r2/r23/mag_panel_comparison_4dp.json`): 91 of 166 leaves move.
  Scoreable cells rise (original 98,370 → 99,590; k-means 99,510 → 99,600), as the R2.2 data predicts. Its byte
  reproducibility on unchanged inputs was not re-tested, because the inputs changed.
- **Abstract:** the abstract-level R2.3 numbers do not move at print precision. L74 moves: 11.1% → 9.8% no EC
  prediction, 9.6% → 9.5% abstain. L67 moves: probe EC F1 0.783 → 0.829.

## Commits (prFBA, branch `hardware-aware-alignment-dispatch`; not pushed; path-limited)

| commit | content |
|---|---|
| `772832a` | fetcher: both selections by default, prefetch User-Agent, rate limit; README D10 note; tests |
| `c21b029` | R2.4: P2D `comparison_summary.json`, `annotation_dialect.json`, stock figure |
| `d354551` | R2.5: P2C ceiling, dialect, comparison summary, depth-probe results, `mag_score_summary.json`, figure; README D10 note with the measured effect |
| `7b14dc3` | `af_gate_counterfactual.py` reads its 16S levels; test |
| `ec1d0c3` | R2.3: three_panel, two_panel, MAG pairings, crossover, AF gate, d100-legacy self-recovery, MAG figures, new `paired_tau99_d100.json` |

Untouched: `presentations/`, ModelSEED checkouts, KBase, codiffusion/EmilyKin, paper repo.

## For Run 3 and later (blockers and notes)

1. **The shared-reference answer key has changed.** These inputs to `score_folds.py` and the rest of R3 are regenerated:
   - `P2C/data/genome_function_sets.json` (sha256 changed from `6f4f817c…`)
   - `function_maps.json` (ac619c6e pin; +305 `pgfam_ec`)
   - `annotation_ceiling.json` (n 2,719)

   `answer_key_fp` `a9dc024df800` will move. W3 (`refresh_panel_ec_cache.py --functions`) must read the new GFS.
   Keep the hazard order: regenerate `_cache_panel_ec.json` before moving `kfold_audit/panel_ec_sets.json` aside.
2. **The pre-repair copies W6/W7 read are in `pre_repair/`:**
   - git-ignored: P2C and P2D `comparison_cells.csv.gz`, both selections, maps and GFS (moved aside);
   - copies: both `pgfam_products.json` (`pre_repair/`, sha256 in `stamps/R2L.products_before.sha256`), the MAG-arm
     `panel_ec.json`, `cells.csv.gz` and `prfba_mag_selections.csv.gz`.
3. **R5.5:** re-stamp `d20_baseline.sha256`. It fails 17 of 32 files: R1.4's 10, R2.5's 5 P2C files, and R2.3's 2
   three_panel files. All are re-score outputs.
4. **R5.4 docs:**
   - P2C/P2D READMEs still quote n = 67,038 and the old results tables.
   - `mag_validation/out/RESULTS_*.md` were not regenerated.
5. **Memory for later runs.** `prfba_selection` at d100 needs about 120 GiB with 8 workers; 12 workers OOMed at 90G.
   `mag_prfba_ec.py` peaks near 108 GiB (40 processes). Plan Run 3 caps accordingly.
6. **Not run, by design:**
   - `three_panel/compare_three_panels_mag.py` (`three_panel_mag_gtdb.json`), whose GTDB arm was not re-scored (map deferred);
   - the `mag_validation/out/RESULTS_*.md` text.
7. **Disk:** about 15 GB of truncated-snapshot controls (`control/r2/snap*.json`) on top of Run 1's 35 GB. They can
   go after Run 2 is accepted.
