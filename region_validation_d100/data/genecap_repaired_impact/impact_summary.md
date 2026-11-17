# Depth-100 10k-panel gene capture: re-scored on the repaired PGFam cache

**2026-11-16 (re-score):** the committed files now equal this measurement's *new* pass: the R1.3 POST re-score reproduced `repro_0987` byte-identically, and R1.5 reproduced `gene_capture_tau_compare.json` (in content; only `_meta.repro_0987_byte_identical`, now true, differs). Everything below is the 2026-11-13 record.

2026-11-13. Measurement only: no committed file, cache, repo or other project was written. The analysis covers `region_validation_d100` (67,991 amplicons in 8 regions). The committed outputs are `gene_capture_cells.csv` and `gene_capture_by_region.json` (2026-10-20), plus `gene_capture_*_tau0995` and `gene_capture_tau_compare.json` (2026-11-11).

## How it was run

- `region_validation_d100/scripts/genecap_tau0995.py` was imported unchanged. Only its module-level output paths were redirected to this directory.
- It was scored twice against one loaded provider:
  - **Old pass:** the 1,006 repaired shared-cache entries were set back to their `backup_old_values` from `CORRECTIONS_2026-11-13.json`. It reproduces **all 7 committed files byte for byte**: the cells, long-form and by-region outputs at 0.987, the same three at 0.995, and `tau_compare`. So every difference below comes from the repair and nothing else.
  - **New pass:** the repaired cache on disk (size 5,118,271,515; sha256 `05f03bcd…`, unchanged during the run).
- The script's own 0.987 byte-identity check now fails, by design. The tables keep the same keys, but rows differ: 1,389 include rows, 2,343 at exclude_0.987 and 2,171 at exclude_0.995.
- 781 repaired genomes are in the scoring universe: 674 were empty lists (`[]`) when they should not have been ("stale"), and 107 were truncated.
- No BV-BRC request was made: 0 genomes were missing, and fetching was disabled.
- **The selected representatives (reps) did not change in any cell**, in any mode. The repair changed only PGFam sets, never selection.
- Each cell is classed by the repaired genomes among its source and its reps (old and new).
  - Cells with no repaired genome: 0 changed.
  - The 240 metric values recomputed from the cells match the summary JSONs exactly.
- The `foreign_patric_id` amendment (kbase files 294.202 and 44139.8) does not affect this result: gene capture reads only the PGFam cache.

## Headline (pooled over the 8 regions)

| mode | population recall % | population F1 % | resolved recall % | truth_empty | pred_empty | OK |
|---|---|---|---|---|---|---|
| include | 99.40 → 99.43 (+0.03) | 91.24 → 91.22 (−0.01) | 99.44 → 99.44 (+0.01) | 953 → 278 | 28 → 13 | 67,010 → 67,700 |
| exclude_0.987 | 60.17 → 61.04 (**+0.86**) | 53.73 → 54.51 (+0.77) | 67.72 → 67.70 (−0.02) | 953 → 278 | 1,060 → 243 | 59,571 → 61,047 |
| exclude_0.995 | 65.21 → 65.96 (**+0.75**) | 58.10 → 58.79 (+0.68) | 70.44 → 70.42 (−0.02) | 953 → 278 | 953 → 248 | 62,054 → 63,419 |

**Population denominators.** They grow from 67,038 to 67,713 cells in every mode, because 675 truth_empty cells per mode become scorable.

**Abstentions.** Abstentions rise by 16 (0.987) and 15 (0.995). These are newly scorable stale sources whose self-excluded reselection abstains.

**Per region, exclude_0.987 population recall** gains between +0.67 pp (V3-V4) and +1.14 pp (V4). At 0.995 the range is +0.56 pp (FullLength16S) to +0.94 pp (V4). Resolved recall moves by at most 0.05 pp in any region or mode (largest: V7-V9 exclude_0.995, −0.046 pp).

## Attribution: stale `[]` entries, not truncation

- **Stale entries drive the population metrics.** Using the exact additive decomposition, 99–101% of the population recall and population F1 change in both exclude modes comes from cells that involve a stale entry. Truncated-only cells contribute ≤ 0.01 pp (about 1%), and cells that involve both kinds about 0%.
- **Include mode is too small to split.** Every include change is ≤ 0.03 pp, so its shares are not meaningful; population F1 (−0.01 pp), for example, splits 35% truncated and 65% stale.
- **Counterfactual check.** Repairing only the truncated entries leaves exclude_0.987 population recall at 60.18%. Repairing only the stale entries gives 61.03%. The full repair gives 61.04%.
- **This re-score depends on keeping the 839 stale-entry repairs.** If that decision were reversed and only the proven truncations kept, the result would fall back to within 0.01 pp of the committed tables, pooled, in every mode (60.17 → 60.18 at 0.987; 65.21 → 65.21 at 0.995). The largest move in any single region would be 0.011 pp.
- **What the stale repair does:**
  - 675 truth_empty cells per mode, all with a stale source, become scored. At exclude_0.987: 658 become OK (the 653 stale-only ones have a mean recall of 0.688), 16 abstain and 1 is pred_empty.
  - Another 4 cells stay truth_empty but count as changed. Their source, 2916410.3, is legitimately `[]`; a stale rep changed only pred_size.
  - 818 exclude_0.987 pred_empty cells become OK: 815 stale-only (mean recall 0.626) plus 3 in the both class.
  - 482 OK cells with a stale rep gain recall: +4.5 pp on average, with precision −6.2 pp.
- **Why truncation barely registers:** it touches only about 0.5% of cells (351–362 per mode), and its effects cancel.
  - At exclude_0.987, the 102 cells with a truncated source lose 1.2 pp of recall on average (worst −5.3 pp), because truth grows by about 670 PGFams.
  - The 237 cells with a truncated rep gain 1.7 pp on average (best +12.2 pp).
- **Resolved recall barely moves (−0.02 pp).** Its decomposition into shares (stale 121–136%, truncated −20 to −35%) is not meaningful for a change this small.
- **The 278 cells that remain truth_empty** come from 38 sources whose lists are still `[]`. None is a repaired genome; they are among the 963 entries the audit judged legitimately empty.

## Files

This directory keeps `impact.json` (every number below, with checks and attribution per region and
mode), this summary, `impact_tables.md`, and the re-score's `gene_capture_tau_compare.json` and
`run_meta.json`. The per-cell outputs (the re-scored cell tables, the byte-identical
reproduction of the committed files on the pre-repair values, and the pickled cells) were made in
the session scratchpad and are not kept; re-create them with `scripts/genecap_tau0995.py`, its
output paths redirected, once on the repaired cache and once with the corrections log's
`backup_old_values` restored in memory. The committed d100 gene-capture files are unchanged.

## Per region and mode (old → new; Δ in percentage points; trunc / stale / both = shares of the population-recall change)

### include

| region | OK old→new | truth_empty | pred_empty | abstained | resolved recall % | population recall % | population F1 % | Δpop recall: trunc / stale / both |
|---|---|---|---|---|---|---|---|---|
| pooled (all regions) | 67010→67700 | 953→278 | 28→13 | 0→0 | 99.44→99.44 (+0.01) | 99.40→99.43 (+0.03) | 91.24→91.22 (-0.01) | -0% / 100% / 0% |
| FullLength16S | 7215→7268 | 83→30 | 0→0 | 0→0 | 99.98→99.98 (+0.00) | 99.98→99.98 (+0.00) | 99.78→99.78 (-0.01) | 0% / 100% / 0% |
| V1-V2 | 7446→7499 | 83→32 | 3→1 | 0→0 | 99.51→99.51 (+0.01) | 99.47→99.50 (+0.03) | 94.13→94.13 (+0.00) | -0% / 100% / 0% |
| V1-V3 | 7505→7559 | 85→33 | 2→0 | 0→0 | 99.65→99.65 (+0.00) | 99.62→99.65 (+0.03) | 96.16→96.17 (+0.02) | -0% / 100% / 0% |
| V3-V4 | 8249→8312 | 94→34 | 4→1 | 0→0 | 99.36→99.36 (+0.00) | 99.31→99.35 (+0.04) | 90.22→90.22 (+0.00) | -0% / 100% / 0% |
| V4 | 9443→9566 | 158→38 | 7→4 | 0→0 | 99.15→99.16 (+0.01) | 99.07→99.12 (+0.04) | 84.70→84.70 (-0.00) | -0% / 100% / 0% |
| V4-V5 | 9406→9528 | 158→38 | 5→3 | 0→0 | 99.31→99.32 (+0.01) | 99.26→99.29 (+0.03) | 88.36→88.35 (-0.01) | -0% / 100% / 0% |
| V4-V5_944R | 9001→9106 | 142→38 | 4→3 | 0→0 | 99.33→99.34 (+0.01) | 99.29→99.30 (+0.02) | 88.59→88.58 (-0.01) | -0% / 100% / 0% |
| V7-V9 | 8745→8862 | 150→35 | 3→1 | 0→0 | 99.40→99.41 (+0.01) | 99.36→99.39 (+0.03) | 91.34→91.32 (-0.02) | -2% / 102% / 0% |

### exclude_0.987

| region | OK old→new | truth_empty | pred_empty | abstained | resolved recall % | population recall % | population F1 % | Δpop recall: trunc / stale / both |
|---|---|---|---|---|---|---|---|---|
| pooled (all regions) | 59571→61047 | 953→278 | 1060→243 | 6407→6423 | 67.72→67.70 (-0.02) | 60.17→61.04 (+0.86) | 53.73→54.51 (+0.77) | 1% / 99% / 0% |
| FullLength16S | 6643→6766 | 83→30 | 107→36 | 465→466 | 68.72→68.73 (+0.01) | 63.27→63.99 (+0.71) | 60.46→61.16 (+0.70) | 1% / 99% / 0% |
| V1-V2 | 6853→6978 | 83→32 | 125→51 | 471→471 | 69.90→69.88 (-0.02) | 64.31→65.02 (+0.71) | 59.27→59.92 (+0.65) | 1% / 99% / 0% |
| V1-V3 | 6959→7088 | 85→33 | 119→42 | 429→429 | 69.52→69.52 (+0.00) | 64.44→65.19 (+0.75) | 60.18→60.87 (+0.69) | 1% / 99% / 0% |
| V3-V4 | 7224→7359 | 94→34 | 112→35 | 917→919 | 66.68→66.69 (+0.01) | 58.37→59.04 (+0.67) | 51.62→52.25 (+0.63) | 1% / 99% / 0% |
| V4 | 8160→8430 | 158→38 | 172→19 | 1118→1121 | 67.30→67.27 (-0.03) | 58.12→59.26 (+1.14) | 48.30→49.28 (+0.98) | 0% / 99% / 1% |
| V4-V5 | 8172→8415 | 158→38 | 151→24 | 1088→1092 | 66.38→66.36 (-0.02) | 57.64→58.59 (+0.95) | 50.35→51.18 (+0.83) | 0% / 99% / 1% |
| V4-V5_944R | 7863→8067 | 142→38 | 125→21 | 1017→1021 | 66.54→66.52 (-0.02) | 58.10→58.91 (+0.81) | 50.98→51.68 (+0.70) | 0% / 99% / 1% |
| V7-V9 | 7697→7944 | 150→35 | 149→15 | 902→904 | 67.30→67.29 (-0.01) | 59.21→60.31 (+1.09) | 52.29→53.27 (+0.98) | 0% / 99% / 1% |

### exclude_0.995

| region | OK old→new | truth_empty | pred_empty | abstained | resolved recall % | population recall % | population F1 % | Δpop recall: trunc / stale / both |
|---|---|---|---|---|---|---|---|---|
| pooled (all regions) | 62054→63419 | 953→278 | 953→248 | 4031→4046 | 70.44→70.42 (-0.02) | 65.21→65.96 (+0.75) | 58.10→58.79 (+0.68) | 0% / 99% / 0% |
| FullLength16S | 6925→7033 | 83→30 | 91→35 | 199→200 | 71.53→71.52 (-0.01) | 68.66→69.21 (+0.56) | 65.51→66.06 (+0.55) | 1% / 99% / 0% |
| V1-V2 | 7053→7169 | 83→32 | 108→43 | 288→288 | 72.91→72.88 (-0.03) | 69.04→69.67 (+0.63) | 62.76→63.34 (+0.58) | 1% / 99% / 0% |
| V1-V3 | 7132→7258 | 85→33 | 114→40 | 261→261 | 72.18→72.16 (-0.02) | 68.58→69.29 (+0.71) | 63.35→64.05 (+0.69) | 0% / 100% / 0% |
| V3-V4 | 7603→7728 | 94→34 | 103→36 | 547→549 | 69.37→69.36 (-0.01) | 63.91→64.48 (+0.58) | 56.49→57.04 (+0.56) | 0% / 100% / 0% |
| V4 | 8554→8791 | 158→38 | 141→22 | 755→757 | 69.70→69.71 (+0.00) | 63.10→64.03 (+0.94) | 52.63→53.44 (+0.81) | 0% / 99% / 0% |
| V4-V5 | 8519→8753 | 158→38 | 148→30 | 744→748 | 69.05→69.01 (-0.04) | 62.51→63.38 (+0.87) | 54.75→55.54 (+0.79) | 0% / 99% / 0% |
| V4-V5_944R | 8215→8405 | 142→38 | 118→28 | 672→676 | 69.33→69.35 (+0.02) | 63.25→63.99 (+0.74) | 55.61→56.28 (+0.67) | 0% / 99% / 0% |
| V7-V9 | 8053→8282 | 150→35 | 130→14 | 565→567 | 70.20→70.16 (-0.05) | 64.62→65.56 (+0.93) | 57.13→57.96 (+0.83) | 1% / 99% / 0% |

### Share of each pooled change by cell class (exact additive; percentage points)

| mode | metric | Δ total | truncated-only | stale-only | both | none |
|---|---|---|---|---|---|---|
| include | recall_resolved | +0.01 | -0.00 (-2%) | +0.01 (102%) | +0.00 (0%) | +0.00 |
| include | recall_population | +0.03 | -0.00 (-0%) | +0.03 (100%) | +0.00 (0%) | +0.00 |
| include | f1_population | -0.01 | -0.00 (35%) | -0.01 (65%) | -0.00 (0%) | +0.00 |
| exclude_0.987 | recall_resolved | -0.02 | +0.01 (-35%) | -0.02 (136%) | +0.00 (-1%) | +0.00 |
| exclude_0.987 | recall_population | +0.86 | +0.00 (1%) | +0.85 (99%) | +0.00 (0%) | +0.00 |
| exclude_0.987 | f1_population | +0.77 | -0.01 (-1%) | +0.78 (101%) | +0.00 (0%) | +0.00 |
| exclude_0.995 | recall_resolved | -0.02 | +0.00 (-20%) | -0.02 (121%) | +0.00 (-2%) | +0.00 |
| exclude_0.995 | recall_population | +0.75 | +0.00 (0%) | +0.75 (99%) | +0.00 (0%) | +0.00 |
| exclude_0.995 | f1_population | +0.68 | -0.00 (-1%) | +0.69 (101%) | +0.00 (0%) | +0.00 |

### Cells by class (pooled): cells / changed / rep set changed

| mode | truncated-only | stale-only | both | none |
|---|---|---|---|---|
| include | 351 / 351 / 0 | 1032 / 1032 / 0 | 6 / 6 / 0 | 66602 / 0 / 0 |
| exclude_0.987 | 362 / 362 / 0 | 1971 / 1971 / 0 | 10 / 10 / 0 | 65648 / 0 / 0 |
| exclude_0.995 | 344 / 344 / 0 | 1821 / 1821 / 0 | 6 / 6 / 0 | 65820 / 0 / 0 |

### Outcome transitions of changed cells (pooled)

* **include**: truncated_only: ok->ok 351; stale_only: ok->ok 342, pred_empty->ok 15, truth_empty->ok 675; both: ok->ok 6
* **exclude_0.987**: truncated_only: abstained->abstained 8, ok->ok 354; stale_only: ok->ok 482, pred_empty->ok 815, truth_empty->abstained 16, truth_empty->ok 653, truth_empty->pred_empty 1, truth_empty->truth_empty 4; both: ok->ok 2, pred_empty->ok 3, truth_empty->ok 5
* **exclude_0.995**: truncated_only: abstained->abstained 5, ok->ok 339; stale_only: ok->ok 442, pred_empty->ok 705, truth_empty->abstained 15, truth_empty->ok 656, truth_empty->pred_empty 1, truth_empty->truth_empty 2; both: ok->ok 2, pred_empty->ok 1, truth_empty->ok 3

### Counterfactual (pooled): population recall % with only one class repaired

| mode | old | truncated entries only | stale entries only | full repair |
|---|---|---|---|---|
| include | 99.40 | 99.40 | 99.43 | 99.43 |
| exclude_0.987 | 60.17 | 60.18 | 61.03 | 61.04 |
| exclude_0.995 | 65.21 | 65.21 | 65.95 | 65.96 |

### Other committed results that read the pre-repair PGFam data

* **k-means d100 gene capture** (region_validation_kmeans_d100/data/gene_capture_{cells.csv,by_region.json,long.csv}). Reads: shared PGFam cache via score_genecap.py (pre-repair). Expected: Same mechanism as the d100 re-score measured here (pooled exclude_0.987: population recall +0.0086, resolved recall -0.0002; include: resolved +0.0001). Expect population recall to rise in both modes as stale-rep pred_empty cells become OK. Stale-source truth_empty cells joining the denominators move resolved recall in either direction. Here the stale exposure is smaller for sources (41 of 318 truth_empty) and 323 of 597 exclude pred_empty cells have a stale rep.
* **random d100 gene capture** (region_validation_random_d100/data/gene_capture_{cells.csv,by_region.json,long.csv}). Reads: shared PGFam cache via score_genecap.py (pre-repair). Expected: Same mechanism as the d100 re-score measured here (pooled exclude_0.987: population recall +0.0086, resolved recall -0.0002; include: resolved +0.0001). Expect population recall to rise in both modes as stale-rep pred_empty cells become OK. Stale-source truth_empty cells joining the denominators move resolved recall in either direction. 183 of 514 truth_empty cells have a stale source and 234 of 416 exclude pred_empty cells a stale rep.
* **E3 validation gate (d100, exclude_0.987 / 0.995)** (region_validation_d100/data/pfunc_e3_gate.json, pfunc_e3_gate/). Reads: shared PGFam cache (truth, reps, neighbours) + e3_clade_catalogue.json built from it (the catalogue excluded all 1,802 [] lists; the rebuilt one has 839 more genomes). Expected: small; the k/n baseline gains most (stale reps were counted as gene-less in n, which pulled k/n down), so the E3-minus-k/n margin should shrink slightly; a 0.053 Brier margin with a CI half-width of about 0.002 makes a verdict flip very unlikely. About 675 more amplicons per mode become scorable. Needs a re-run against the rebuilt catalogue to state it as a number.
* **shared_reference k-fold benchmark (PICRUSt2 0.779 vs prFBA 0.687 macro EC F1)** (region_validation/picrust2_comparison/data/genome_function_sets.json, shared_reference per-fold _cache_pgfam_subset.json, comparison cells). Reads: frozen copies of the pre-repair lists (a cache repair does not change them); 17 truncated panel genomes, 0 stale. Expected: negligible (about 0.2% of cells per arm). Truth EC sets grow for 17 genomes, so recall falls a little on those cells in every arm. prFBA also gains predicted ECs on about 0.3% of cells with a truncated rep. The PICRUSt2-over-prFBA ranking cannot plausibly change.
* **MAG validation (original, k-means and random panels)** (mag_validation/out/{mag_cells_all,mag_cells_tau99}.csv, mag_summary_*.json; mag_validation_km/out, mag_validation_rnd/out). Reads: NOT the shared cache: frozen int-encoded panel indices extracted from it on 2026-10-08 (mag_validation/pgfam_panel.json, kmeans_panel/pgfam_panel_kmeans.json, random_panel/pgfam_panel_random.json). Not in the audit census; checked here: every repaired panel genome holds its OLD length. A cache repair does not reach these indices; they need re-extraction. Expected: 1,220 of the 1,600 truth_empty MAG cells per regime have a stale source and would become scorable. Exclude-regime OK cells with a stale anchor (287) gain predicted families, so recall rises. Truncated sources (170 cells per regime) lower recall slightly. Net: more scorable cells and a small recall rise, but only after the panel indices are rebuilt; as committed they cannot change.
