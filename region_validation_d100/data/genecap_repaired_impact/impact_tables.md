> **Re-score 2026-11-16: the numbers in this document are historical.** These tables are the 2026-11-13 measurement of what the cache repair does (old pass vs new pass); the committed gene-capture files now equal the 'new' pass (re-score R1.3/R1.5). Current values: `manuscript_numbers/number_diff.md` (every manuscript claim, old -> new, with its source artefact and tex lines) and the re-scored artefacts themselves (re-score runs R1-R4 on the 2026-11-13 repaired BV-BRC PGFam cache; prFBA commits 81e9b34 through bb41947).


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
* **MAG validation (original, k-means and random panels)** (mag_validation/out/{mag_cells_all,mag_cells_tau99}.csv, mag_summary_*.json; mag_validation_km/out, mag_validation_rnd/out). Reads: NOT the shared cache: frozen int-encoded panel indices extracted from it on 2026-10-08 (mag_validation/pgfam_panel.json, kmeans_panel/pgfam_panel_kmeans.json, random_panel/pgfam_panel_random.json). Not in the audit census; checked here: every repaired panel genome holds its OLD length. A cache repair does not reach these indices; they need re-extraction.. Expected: 1,220 of the 1,600 truth_empty MAG cells per regime have a stale source and would become scorable. Exclude-regime OK cells with a stale anchor (287) gain predicted families, so recall rises. Truncated sources (170 cells per regime) lower recall slightly. Net: more scorable cells and a small recall rise, but only after the panel indices are rebuilt; as committed they cannot change.
