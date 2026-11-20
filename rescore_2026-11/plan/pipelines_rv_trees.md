# Rescore map: region-validation trees, panel comparisons, depth probe, manuscript figures

Read-only survey, 2026-11-15. prFBA HEAD `cc8014f` (branch `hardware-aware-alignment-dispatch`), paper HEAD `0214a14`.
Nothing was written outside this scratch directory. Helper scripts and their outputs are in
`rescore/rv_trees_work/` (`stage_times.py`, `legacy_exposure.py`, `legacy_exposure_*.json`).

`PY=/home/freiburger/Documents/py_venv/bin/python`, `R=/home/freiburger/Documents/prFBA`,
`CACHE=$R/bvbrc_cache/genome_gene_families.json` (repaired 2026-11-13; 5,118,271,515 bytes;
sha256 `05f03bcd…`; 159,753 keys, key-set sha256 unchanged by the repair).

---------------------------------------------------------------------------------------------

## 0. Findings that decide the plan

1. **The scoring code has not changed since the committed runs; only the data and the fetchers have.**
   Last commits: `select_references.py` 7153311 (2026-10-20 09:43, before the d100 POST at 22:01);
   every d100 `scripts/score_*.py`, `prefetch_universe.py`, `mcnemar_regions.py`,
   `diagnose_self_recovery.py`, `kmeans_panel/close_exclude_rep_gap.py`: b73d768 (2026-10-20 setup);
   `three_panel/emptied_strata.py` b73d768; `three_panel/compare_three_panels.py` 5d2a204 (10-29, the
   outputs were regenerated in that commit); `two_panel/compare_panels.py` b4f9db1 (10-18);
   `genecap_tau0995.py` f5abb6d; `pfunc_e3_panel_gate.py` 4f7871d. Changed since: `prefetch_gene_families.py`
   and the nine `region_validation*/scripts/prefetch_pgfams.py` copies (paging, `eq(annotation,PATRIC)`,
   raise instead of partial), which only matter for genomes **not** already cached; `pfunc_estimator.py`
   (817f20b, 552e23f, 8cfc297) — 817f20b's "keep k/n when E3 is 0" sits in `reweight_asv` (build path),
   not in `e3_probabilities`, which the gate calls, so the gate's formula is unchanged;
   `align_hits.py` (d0a926f, dbbc8f7, 15f3caa: GPU picks the top 100, Biopython only traces back), which
   re-ran stage 3 of the original d100 panel from its saved shortlist **byte-identically**
   (GPU_CAND_200.md). So "latest code" adds nothing to re-run beyond what the cache repair forces.
2. **Union trees vs legacy trees behave differently.**
   - Under the exact-tie union (`select_all=True`, every `*_d100`, d20 kmeans, d20 random),
     `select_references.select_representatives` never reads the gene provider when choosing members
     (the stage-7 block that reads it runs only for `select_all=False`; `select_references.main` passes
     no `gene_count_fn`, so the P3 tie-break is off and `est_genes` is the 3000 default). Measured:
     the 2026-11-13 census (scratch `bvbrc_trunc/census_*.json`) re-selected every amplicon that
     touches a repaired genome on the repaired cache and found `include_selection_reproduced`
     20,644/20,644 (original), 13,268/13,268 (k-means), 13,005/13,005 (random); the repaired-cache
     re-score (`region_validation_d100/data/genecap_repaired_impact/`) found "rep set changed" 0 in every
     mode. Hence `selection.json`, all `score_taxacc` outputs, McNemar, region meta and the
     self-recovery diagnosis of the union trees are **reusable**; only gene capture and its consumers move.
   - Under the legacy reducer (`*_d100_legacy`, d20 original), stage 7 asks the provider for each
     candidate's PGFam set and **drops a genome whose set is empty as `no_gene_set`**, and uses set
     sizes for the marginal-gain stop. The 839 stale `[]` entries were therefore excluded as gene-less.
     Measured here by streaming each legacy `selection.json` (`rv_trees_work/legacy_exposure_*.json`):

     | legacy selection | records | records touching a repaired genome | records with a repaired genome excluded at `no_gene_set` | records with a repaired (truncated) genome selected |
     |---|---|---|---|---|
     | region_validation_d100_legacy | 67,991 | 20,644 | **1,317** | 187 |
     | region_validation_kmeans_d100_legacy | 74,704 | 13,265 | **933** | 162 |
     | region_validation_random_d100_legacy | 75,328 | 12,994 | **1,674** | 134 |
     | region_validation (d20 original, legacy) | 67,991 | 7,606 | **1,234** | 157 |

     So every legacy tree must be re-run **from `select_references` onward**, taxonomy included
     (include-mode anchors come from `selection.json`), plus `diagnose_self_recovery.py` (its cited
     `failure_terminal_stage_overall.no_gene_set` = 961 for d100 legacy is exactly this effect).
3. **No BV-BRC fetch is expected.** Every prefetch enumerates its genome set from the hits and
   knobs only (`prefetch_gene_families.candidate_genome_ids` runs `build_audit` with no provider), the
   repair kept the key set (`keys_sha256_unchanged`), and every one of those sets was fetched on
   2026-10-20/21. The only possible misses are exclude-mode reps the legacy re-selection newly names
   in `close_exclude_rep_gap.py` (1–2 per tree at the time). Network is a fallback, not a step.
4. **POST peak memory has never been recorded**, and one POST was OOM-killed (below). Plan for it.
5. **The d20 trees cannot be re-run through their own drivers** without silently changing the cited
   rule or the alignment (section 6, H5). Re-scoring them at all is a decision, not a default
   (section 7, Q1).

---------------------------------------------------------------------------------------------

## 1. Trees and what cites them

| tree | panel | record depth (k2) | rule (selection / re-selection) | aligns | cited by | scripts vs `region_validation_d100/scripts` |
|---|---|---|---|---|---|---|
| region_validation_d100 | original (taxonomy-stratified) | 100 | exact-tie union | yes | manuscript throughout (taxonomy, gene capture, 0.995, McNemar, E3 gate, F1/F2/F4b/F11), three-panel, two-panel, MAG-vs-16S, deck | reference; also holds `genecap_tau0995.py`, `pfunc_e3_panel_gate.py`, `pfunc_e3_gate_review_recompute.py` |
| region_validation_d100_legacy | original | 100 | legacy | no (hits copied) | precision 0.856, legacy self-recovery, d20-vs-d100 legacy comparisons | `_common.BASE` only |
| region_validation_kmeans_d100 | k-means | 100 | union | yes | three-panel, two-panel tables, emptied strata | `_common.BASE`, `build_truth.py` (manifest taxid patch) |
| region_validation_kmeans_d100_legacy | k-means | 100 | legacy | no | `_d100_legacy` comparisons | same |
| region_validation_random_d100 | random | 100 | union | yes | three-panel, emptied strata | same |
| region_validation_random_d100_legacy | random | 100 | legacy | no | `_d100_legacy` comparisons | same |
| region_validation | original | 20 | legacy (selection June 10; genecap `--legacy` 07-22; taxacc `--legacy` 09-17) | done | d20-vs-d100 comparisons (inventory items `legacy_tabmain_genus_d20_vs_d100`, `legacy_genecap_d20_vs_d100`, `lim.depth20_vs_100`), depth probe (truth, fan-out), PICRUSt2 d20 comparison, d20 baseline | older `prefetch_universe.py` (no `--select-mode`), `prefetch_pgfams.py` (repaged 8cfc297), `render_accuracy.py`, `diagnose_self_recovery.py` |
| region_validation_kmeans | k-means | 20 | union | done | d20 baseline; D100_VS_D20.md | as above + `build_truth.py` |
| region_validation_random | random | 20 | union | done | d20 baseline; D100_VS_D20.md | as above + `build_truth.py` |

Committed (git-tracked) outputs per d100 tree: `data/{ambiguity,amplification_by_phylum,asv_region_concordance,concordance_long,confusion_genus,mcnemar_region_pairs,self_recovery}.csv`,
`data/{gene_capture_by_region,mcnemar_regions,taxacc_region_summary,benchmark_truth,selection_10k}.json`,
`data/gene_capture_cells.csv`, `data/hits/raw/{asv_alignment_summary.csv,method_selection.json}`;
plus `self_recovery_diagnosis.json`/`self_recovery_failures.csv` in the two original trees; plus in
`region_validation_d100`: `gene_capture_{by_region,cells}_tau0995.*`, `gene_capture_tau_compare.json`,
`genecap_tau0995/run_meta.json`, `pfunc_e3_gate.json`, `pfunc_e3_gate/{per_amplicon.csv.gz,run_meta.json,review_recompute.json,production_ek_retro.json}`,
`genecap_repaired_impact/*`. Big intermediates (hits 3.3–4.4 GB, `selection.json` 3.2–4.0 GB,
`taxacc_per_record.jsonl`, `gene_capture_long.csv`, shortlist/GPU scores) are on disk, git-ignored.
Tree `figures/` are not committed.

Paper figures (`prFBA_paper/images/prfba/`, md5-matched): F1, F11, F2 (= `F2_perrank_concordance_exclude_self_ALT`), F4b, F14 (= `F_mag_completeness`), F15 (= `F_mag_crossover`) are byte copies of `manuscript_figures/out/`; F16_panel_space = `two_panel/F_panel_space.pdf`. F13 is PICRUSt2 (other report).

---------------------------------------------------------------------------------------------

## 2. Command sequences (as the committed files were produced)

Dependency class per step: **N** = reads no BV-BRC gene data (16S, taxonomy, hits only);
**S** = reads the PGFam cache only through the legacy selection (cache-dependent in legacy trees,
cache-independent in union trees); **C** = reads PGFam sets directly.

### 2.0 Upstream of everything (N, reuse)

Front-load (June 10 original; Oct 9 k-means; Oct 16 random): `select_10k.py` / `kmeans_panel/select_kmeans.py` /
`random_panel/select_random.py` → `insilico_pcr.py --emax 3` → `dedup_amplicons.py` → `build_truth.py`; the
seven files are byte copies in the d100 trees (`run_validation.sh` checks they exist and are not
symlinks). Alignment (unit `rv-d100-align`, `run_align_all_d100.sh`, 2026-10-20/21):

```
$PY $R/align_hits.py --fasta <T>/data/combined_amplicons.unique.fasta --outdir <T>/data/hits/raw \
  --prefilter-k 500 --k2 100 --gpu-cand 500 --workers 40 --min-ref-len 200 --force-backend cuda --save-shortlist
```
Legacy trees receive copies via `seed_legacy_hits_d100.sh <union tree>` (cp + cmp + link-count check,
run by the POST driver each time). None of this reads PGFam data: **reuse**.

### 2.1 POST, one d100 tree `T` (unit `rv-d100-post` → `run_post_all_d100.sh` → `STAGE=post NO_TMUX=1 POST_LOCK_HELD=1 bash $R/$T/scripts/run_validation.sh full`)

Working directory `$R/$T/scripts`. `MODE=exact-tie`, `LEG=""` for union; `MODE=legacy`, `LEG=--legacy` for legacy.
Checkpointed stages (`run_stage`) skip when their output exists; the prefetch stages always run.

| # | command | reads | writes | class |
|---|---|---|---|---|
| P0 | inline python: raw md5-keyed hits × `amplicon_md5_fanout.json` (checkpoint `data/hits/asv_top20_alignment_hits.json`) | hits/raw, fan-out | fanned hits (4.0–4.4 GB) | N |
| P1 | `$PY $R/prefetch_gene_families.py --hits $T/data/hits/asv_top20_alignment_hits.json --cache $CACHE --workers 16 --select-mode $MODE` (always runs) | hits, cache (json.load) | cache (only if something was fetched; save every 2000) | C (fetch only) |
| P2 | `$PY $R/select_references.py --hits … --out data/selection.json.tmp --gene-provider auto --gene-cache $CACHE --select-mode $MODE && mv` (checkpoint `selection.json`) | hits, cache (all lists → sets) | `selection.json` (3.2–4.0 GB) | union: N in content; legacy: **C** |
| P3 | `$PY $T/scripts/prefetch_pgfams.py` (truth = `selection_10k.json`, reps = `selection.json`; always runs) | cache, selection | cache (only if fetched) | C (fetch only) |
| P4 | `$PY $T/scripts/prefetch_universe.py --select-mode $MODE` (always runs) | hits, cache | cache (only if fetched) | C (fetch only) |
| P5 | `$PY $R/kmeans_panel/close_exclude_rep_gap.py --tree $T $LEG` (checkpoint `data/close_exclude_rep_gap.done`, log `.log`) | selection, hits, cache (loaded per round; reloaded whole to merge) | cache (only if fetched) | union: N; legacy: C |
| P6 | `$PY $T/scripts/score_taxacc.py $LEG` (checkpoint `taxacc_region_summary.json`) | hits, selection, truth; exclude re-selection has **no** provider | `taxacc_per_record.jsonl`, `concordance_long.csv`, `asv_region_concordance.csv`, `self_recovery.csv`, `confusion_genus.csv`, `ambiguity.csv`, `taxacc_region_summary.json` | union N; legacy S |
| P7 | `$PY $T/scripts/score_genecap.py $LEG` (checkpoint `gene_capture_by_region.json`) | selection, hits, cache via `bvbrc_gene_provider` (whole cache as sets); exclude re-selection with provider | `gene_capture_cells.csv`, `gene_capture_long.csv`, `gene_capture_by_region.json` | **C** |
| P8 | `$PY $T/scripts/mcnemar_regions.py` (checkpoint `mcnemar_region_pairs.csv`) | `taxacc_per_record.jsonl` | `mcnemar_region_pairs.csv`, `mcnemar_regions.json` | union N; legacy S |
| P9 | `$PY $T/scripts/build_region_meta.py` (checkpoint `region_meta.csv`) | stats, manifest, truth, fan-out | `region_meta.csv`, `amplification_by_phylum.csv` | N |
| P10 | `render_accuracy.py`, `render_genecap.py` (no `--force` → **skip existing figures**) | the CSVs above | `$T/figures/*` (not committed) | follows inputs |
| P11 | driver: `$PY $R/region_validation_d100/check_d100.py post $T` → `data/post.done` | | | check |

### 2.2 Region_validation_d100-only steps

| # | command (as run) | inputs | outputs | class |
|---|---|---|---|---|
| E1 | `cd $R/region_validation_d100 && systemd-run --user … $PY -u scripts/genecap_tau0995.py > logs/genecap_tau0995.log` (smoke first: `--limit 300`); then `scripts/genecap_tau0995.py --merge-supplement` if it fetched | selection, hits, cache streamed (ijson) | `data/gene_capture_{cells,long}_tau0995.csv`, `data/gene_capture_by_region_tau0995.json`, `data/gene_capture_tau_compare.json`, `data/genecap_tau0995/{run_meta.json,pgfam_supplement.json,repro_0987/}`; merge rewrites `$CACHE` | C |
| E2 | `$PY $R/pfunc_estimator.py build-catalogue --cache $CACHE --out $R/bvbrc_cache/e3_clade_catalogue.json` (only when the cache changed; done 2026-11-13 10:43) | cache streamed | catalogue | C |
| E3 | `$PY region_validation_d100/scripts/pfunc_e3_panel_gate.py --workers 12` (cwd `$R`; smoke `--limit 500 --workers 8 --verify-q`) | selection, hits (ijson), cache (int store), catalogue; refuses a stale catalogue | `data/pfunc_e3_gate.json`, `data/pfunc_e3_gate/{per_amplicon.csv.gz,run_meta.json}` | C |
| E4 | `PFUNC_GATE_REVIEW_OUT=<scratch> WORKERS=8 $PY region_validation_d100/scripts/pfunc_e3_gate_review_recompute.py`, then copy `indep_recompute.json` → `data/pfunc_e3_gate/review_recompute.json` | gate CSV, hits, cache, catalogue | review JSON (cited caveats in PFUNC_E3_GATE.md) | C |
| E5 | measurement only: `genecap_repaired_impact/` (scratch driver `bvbrc_foreign/genecap_scripts/rescore_repaired.py` importing `genecap_tau0995.py` with redirected paths) | | `impact*.json/md`, `run_meta.json` | superseded by the rescore |
| E6 | `production_ek_retro.json`: retro E3 over the old mixed EmilyKin build (scratch tool) | | | deployment scope; superseded by the BV-BRC-only rebuild |

### 2.3 Cross-panel downstream (unit `rv-d100-downstream` → `run_downstream_d100.sh`; needs `post.done` in all six trees; for `S` in `_d100`, `_d100_legacy`; `O=$R/region_validation$S`, `K=$R/region_validation_kmeans$S`)

| # | command | reads | output | class |
|---|---|---|---|---|
| D1 | `$PY three_panel/compare_three_panels.py --suffix $S` | 3 trees' `selection_10k.json`, `benchmark_truth.json`, `asv_region_concordance.csv`, `gene_capture_cells.csv` | `three_panel/three_panel_comparison$S.json` (seeded bootstrap B=1000, seed 1729) | genecap endpoints C; taxonomy endpoints N (union) / S (legacy) |
| D2 | `$PY three_panel/emptied_strata.py --suffix $S` | 3 trees' `gene_capture_cells.csv` (FL, exclude_0.987, outcome != truth_empty) + `method_selection.json` | `three_panel/emptied_strata$S.json` (deck EMPTIED, tab:emptied) | C (population = non-truth_empty) |
| D3 | `$PY two_panel/compare_panels.py --orig-truth $O/data/benchmark_truth.json --km-truth $K/data/benchmark_truth.json --orig-genecap $O/data/gene_capture_cells.csv --km-genecap $K/data/gene_capture_cells.csv --out two_panel/panel_comparison$S.json` | | `two_panel/panel_comparison$S.json` | C |
| D4 | `$PY two_panel/mechanism_depth.py --depth 100 --original-root $O --kmeans-root $K --out two_panel/mechanism_depth$S.json` | truths, FL genecap cells, MAG cells `mag_validation/out_batches_n500/b*/mag_cells.csv`, `mag_validation_km/out_batches/b*/mag_cells.csv` | `two_panel/mechanism_depth$S.json` | C (both arms) |
| D5 | `$PY $O/scripts/diagnose_self_recovery.py` | `$O` selection.json, truth | `$O/data/self_recovery_diagnosis.json`, `self_recovery_failures.csv` | union N; legacy S |
| D6 | `$PY mag_validation/compare_arms.py --amp $O/data/gene_capture_cells.csv --mag-cells mag_validation/out/mag_cells_all.csv --out mag_validation/out/paired_comparison$S.json` | | → F14 | C (both arms) |
| D7 | `$PY mag_validation/crossover.py --amp $O/data/gene_capture_cells.csv --out mag_validation/out/crossover$S.json mag_validation/out/mag_cells_all.csv mag_validation/lowc_out_batches/b*/mag_cells.csv` | | → F15 | C (both arms) |
| D8 | `$PY region_validation_d100/check_d100.py all` | | checks 1–6 | check |

RUNBOOK §5 warns that the `--mag-cells` / crossover inputs of D6–D7 were inferred from the cited
outputs' tier sets, not from a recorded command.

### 2.4 Other cross-panel outputs in scope

| output | producer | inputs | class |
|---|---|---|---|
| `three_panel/three_panel_mag_panelarm.json` | `$PY three_panel/compare_three_panels_magpanel.py` | MAG panel-arm cells `mag_validation/out_batches_n500/b*`, `mag_validation_km/out_batches/b*`, `mag_validation_rnd/out_batches/b*` (`mag_cells.csv`, scored by `mag_validation/score_mag.py --pgfams <panel index>`), d20 truths via `PANELS` | C (via panel indices) |
| `two_panel/mag_panel_comparison.json` | `$PY two_panel/compare_mag_panels.py` | same MAG cells (original, k-means) | C (via panel indices) |
| `three_panel/three_panel_mag_gtdb.json` | `$PY three_panel/compare_three_panels_mag.py` | `mag_validation{,_km,_rnd}/out_gtdb/gtdb_cells*.csv` (ANI/taxonomy only), d20 truths | N: reuse |
| `three_panel/space_deciles_three{,_plot}.json`, `F_panel_space_three*` | `three_panel/plot_space_three.py [--slide]` | `two_panel/figdata_space.npz`, `knn_density.npz`, selections | N: reuse |
| `two_panel/{space_xy.npy,space_meta.parquet}`, `knn_density.*`, `figdata_space.*`, `F_panel_space.*` (= paper F16), `space_deciles.json` | `project_space.py` → `knn_density.py` → `figdata_space.py` → `plot_space.py` | embeddings, index.parquet, selections | N: reuse (and do not re-run `project_space.py`: UMAP is unseeded by design) |
| d20 `three_panel/three_panel_comparison.json`, `two_panel/panel_comparison.json`, `two_panel/mechanism_depth.json` | same scripts, default (d20) roots | d20 trees | C; historical (d20 baseline / RESULTS.md) |

### 2.5 Panel PGFam indices (kmeans_panel/, random_panel/, mag_validation/)

`mag_validation/pgfam_panel.json` (2026-10-08, 106 MB), `kmeans_panel/pgfam_panel_kmeans.json` (10-09, 141 MB),
`random_panel/pgfam_panel_random.json` (10-18, 153 MB): format `{"n_vocab": N, "genomes": {gid: [int ids]}}`,
each file its own first-appearance vocabulary. **No generator is in the repo** (only `kmeans_panel/pgfam_index.log`:
"index: 10000 genomes, 1879463 distinct families"). Consumers: `mag_validation/score_mag.py --pgfams`
(`score_km.sh`, `score_rnd.sh`, the original n500 scoring) → MAG cells → section 2.4 rows and D4/D6/D7.
Repaired genomes inside them (counted here; every one still holds its **old** length):

| index | genomes | repaired genomes in index | stale `[]` | truncated | empty lists in index |
|---|---|---|---|---|---|
| original | 10,000 | 140 | 123 | 17 | 161 |
| k-means | 10,000 | 20 | 9 | 11 | 45 |
| random | 10,000 | 38 | 29 | 9 | 72 |

### 2.6 Depth probe / retention analyses

- `region_validation/picrust2_comparison/scripts/depth_probe_sample.py` (seed 1729; N) →
  `align_hits.py --k2 500 --gpu-cand 500` into `data/depth_probe/hits_k500/` (N; 248 s wall, 5.3 GB RSS) →
  `$PY region_validation/picrust2_comparison/scripts/depth_probe_eval.py` (defaults) →
  `data/depth_probe/depth_probe_results.json` (2026-09-22, b72b71c). **C**: it loads the whole PGFam cache
  (orjson) for the P3 gene-count tie-break **and** the EC projection (`function_maps.json` `pgfam_ec`), and uses
  the d20 original tree's truth/fan-out through `region_validation/scripts/_common`. Only the eval re-runs.
- Retention mechanism: `two_panel/mechanism_depth.py` (D4), `three_panel/emptied_strata.py` (D2),
  `check_d100.py` check 4 (truncation at depth 100; N), `diagnose_self_recovery.py` (D5).
- D100_VS_D20.md / UNION_VERIFICATION.md / UNION_THROUGHOUT_CHANGES.md quote ad hoc computations
  (selfrec.py, Welch tests, emptied counts); their gene-capture-derived numbers need recomputation by hand.

### 2.7 manuscript_figures/ (run from each script's directory; output `manuscript_figures/out/`, git-ignored; then copy)

| script | output → paper file | data | rerun? |
|---|---|---|---|
| `cd rv && $PY render_genecap.py --force --only F4b` | `F4b_percent_genecapture_by_region.pdf` → same name | `region_validation_d100/data/gene_capture_cells.csv` | **yes** |
| `cd rv && $PY render_accuracy.py --force --only F1,F11` | F1, F11 | union taxacc CSVs | data unchanged (optional) |
| `cd rv && $PY run_exclude_alt.py` | `F2_perrank_concordance_exclude_self_ALT.pdf` → `F2_perrank_concordance.pdf` | union taxacc | data unchanged (optional; it passes force=True) |
| `cd mag && $PY plot_results_d100.py` | `F_mag_completeness.pdf` → `F14_mag_completeness.pdf` | `mag_validation/out/mag_summary_full.json`, `out_gtdb/gtdb_summary_2000.json`, `out/paired_comparison_d100.json` | **yes**, after D6 and the MAG re-score |
| `cd mag && $PY plot_crossover_d100.py` | `F_mag_crossover.pdf` → `F15_mag_crossover.pdf` | `mag_validation/out/crossover_d100.json` | **yes**, after D7 and the MAG re-score |

The renderers skip existing outputs unless `--force` (except `run_exclude_alt.py`).

### 2.8 How the d20 trees' committed files were actually made (not by their drivers)

- `region_validation` (legacy): `selection.json` 2026-06-10 (June `select_references`, no `select_mode` key);
  `gene_capture_*` 2026-07-22 (54d1150, `score_genecap.py --legacy` by hand); `taxacc_*`, `mcnemar_*`
  2026-09-17 (ed2a2b7, `score_taxacc.py --legacy`); `self_recovery_diagnosis.json` 06-11.
- `region_validation_kmeans` (union): selection 10-09 23:16; taxacc 10-10 19:31; genecap re-scored 10-11 10:38
  after `close_exclude_rep_gap.py --tree region_validation_kmeans` (kmeans_panel/close_gap.log, `rescore.done`).
- `region_validation_random` (union): units `prfba-rv-random{,2,3}` 10-16..10-18.
- The d20 `scripts/run_validation.sh` now aligns with `--k2 100 --gpu-cand 200` (dbbc8f7/d0a926f edited it),
  passes no `--legacy` and no `--select-mode`, and accepts `--force`.

---------------------------------------------------------------------------------------------

## 3. Re-run vs reuse

### 3.1 Matrix

| step | union d100 ×3 | legacy d100 ×3 | d20 original (legacy) | d20 k-means/random (union) |
|---|---|---|---|---|
| alignment, fan-out (P0), front-load | reuse | reuse (copies) | reuse | reuse |
| P1/P3/P4 prefetch | run (driver always runs them; expected 0 fetches) | run (expected 0) | n/a by hand | n/a |
| P2 select_references | **reuse** (cache-independent; census-verified) | **re-run** | re-run if d20 is rescored (see H5) | reuse |
| P5 gap close | reuse (`.done`; union rep set unchanged) | **re-run** | re-run (d20 copy of prefetch_universe has no `--select-mode`) | reuse |
| P6 score_taxacc (+P8 McNemar) | reuse (optionally re-run once into scratch as a byte-reproduction check) | **re-run** | re-run if rescored | reuse |
| P7 score_genecap | **re-run** | **re-run** | re-run if rescored | re-run if rescored |
| P9 region meta | reuse | reuse | reuse | reuse |
| P10 tree figures | re-render with `--force` (not committed) | same | same | same |
| D5 diagnose_self_recovery | reuse (original union) | **re-run** (original legacy; d20 if rescored) | | |
| E1 tau0995 | **re-run** (original union only) | | | |
| E2 catalogue | rebuild only if the cache changed after 2026-11-13 10:38 | | | |
| E3/E4 E3 gate + review | **re-run** (inputs changed: 839 genomes enter the store/catalogue) | | | |
| D1–D4, D6–D7, D8 | **re-run** (both suffixes) | | | d20 versions only if d20 is rescored |
| panel indices → MAG cells → magpanel/two-panel MAG | **regenerate** (other workstream owns MAG scoring) | | | |
| depth_probe_eval | **re-run** (hits_k500 reused) | | | |
| GTDB arm, panel-space figures | reuse | | | |
| manuscript F4b, F14, F15 | **re-render** | | | |

### 3.2 Exact delete lists (d100 trees refuse `--force`; `post.done` makes the driver skip a tree)

Union tree (`region_validation_d100`, `_kmeans_d100`, `_random_d100`):
```
rm data/post.done data/gene_capture_by_region.json      # P7 re-runs; P1/P3/P4 run (0 fetch); P2/P5/P6/P8/P9 skip
# optional, to refresh the tree figures: rm figures/F4*.p* figures/F5b* figures/F8* figures/F9* figures/F12*
```
Legacy tree (`*_d100_legacy`):
```
rm data/post.done data/selection.json data/close_exclude_rep_gap.done data/taxacc_region_summary.json \
   data/gene_capture_by_region.json data/mcnemar_region_pairs.csv
# keep data/hits/ (copies + fan-out), the front-load, region_meta.csv (N)
```
Then `systemd-run --user --unit=rv-rescore-post --collect -p LimitNOFILE=1048576 --working-directory=$R /bin/bash -c 'bash $R/region_validation_d100/run_post_all_d100.sh >> $R/region_validation_d100/logs/post_all_rescore.log 2>&1'`
(add `-p MemoryMax=…`, see H3). For the original legacy tree also re-run D5 (`$PY $R/region_validation_d100_legacy/scripts/diagnose_self_recovery.py`, it is in `run_downstream_d100.sh`).

Deleting `post.done` also closes two other gates: `run_picrust2_d100.sh` (PICRUSt2 workstream) refuses
until both original trees have it again, and `run_downstream_d100.sh` until all six do. The PICRUSt2
d100 step must not start against an original tree mid-rescore (RUNBOOK §4: `prfba_selection.py` is
cache-only and relies on the legacy POST having cached the legacy reducer's candidates).

### 3.3 Order

1. Pre-flight: no cache writer running (`pgrep -af 'prefetch_|select_references|close_exclude_rep_gap|score_genecap|prfba_selection|genecap_tau0995|pfunc_e3'`; no `$CACHE.tmp`); record cache size/mtime/sha256 and catalogue meta; confirm `sha256sum -c region_validation_d100/d20_baseline.sha256` (all 32 OK today).
2. POST of the six d100 trees (driver order; flock). Gate on `check_d100.py post` per tree.
3. `genecap_tau0995.py` (E1) with the old `pgfam_supplement.json` moved aside → its built-in 0.987 check must now be byte-identical to the **new** `gene_capture_*` from step 2 (a free cross-check). If it fetched anything, run `--merge-supplement` (a cache writer).
4. After the **last** cache writer (steps 2–3): if the cache size/mtime differs from the 2026-11-13 10:38 stamp, rebuild the catalogue (E2). Only then:
5. E3 gate (E3) → its `reproduction_vs_gene_capture_cells` must match the new cells → review (E4).
6. Panel indices + MAG re-score (MAG workstream), then `run_downstream_d100.sh` (D1–D8) — or run D1–D3, D5 now and D4/D6/D7 after the MAG cells are new.
7. `depth_probe_eval.py`; `compare_three_panels_magpanel.py`, `compare_mag_panels.py`.
8. manuscript_figures F4b/F14/F15 with `SOURCE_DATE_EPOCH` set; copy to the paper; update hand-entered numbers (H11).

---------------------------------------------------------------------------------------------

## 4. Runtime, memory, network (measured where a log exists)

### 4.1 POST stage walls, 2026-10-20/21 (from `region_validation_d100/logs/post_all.log`, `rv_trees_work/stage_times.py`)

| stage | orig union | orig legacy | k-means union | k-means legacy | random union | random legacy |
|---|---|---|---|---|---|---|
| P0 fan-out | 4.0 | 4.2 | 4.4 | 4.4 | 3.0 | 4.1 |
| P1 candidate prefetch | 6.0 (1,529 fetched, 181 s) | 5.8 (490) | 7.6 (2,218) | 7.7 (2,596) | 4.1 (242) | 7.5 (1,266); 2.9 at resume (0) |
| P2 select_references | 6.2 | 7.5 | 7.8 | 9.1 | 5.4 | 11.6 |
| P3 prefetch_pgfams | 1.5 (0) | 1.5 (0) | 1.7 (0) | 1.6 (0) | 1.2 (0) | 1.5; 1.2 at resume |
| P4 prefetch_universe | 9.5 (1,029, 297 s) | 7.7 (896) | 8.9 (1,553) | 9.2 (1,862) | 6.0 (555) | 7.3 (977); 3.7 at resume (0) |
| P5 gap close | 11.4 (1 fetched) | 6.4 | 5.9 | 18.0 (2) | 8.0 (1) | OOM-killed at 1.8 min; 11.7 at resume (1) |
| P6 score_taxacc | 5.2 | 5.5 | 5.3 | 4.5 | 3.5 | 3.7 |
| P7 score_genecap | 7.3 | 7.2 | 8.4 | 7.1 | 6.8 | 6.9 |
| P8 McNemar | 0.5 | 0.5 | 0.5 | 0.5 | 0.3 | 0.3 |
| P9 region meta | <0.1 | <0.1 | <0.1 | <0.1 | <0.1 | <0.1 |
| P10 renders | ~1.5 | ~1–1.5 | ~1–1.5 | ~1–1.5 | ~1–1.5 | ~1–1.5 |
| **tree total** | **53** | **48** | **51** | **64** | **39** | **32 (resumed) + 34 before the OOM** |

Minutes. Unit CPU: `rv-d100-post` 4 h 46 min until the OOM kill (2026-10-21 12:30:39), then 32 min.
Fetch rate then: 3–10 genomes/s at 16 workers, 0 errors.
**Rescore estimate**: union tree ≈ P1+P3+P4 (≈ 12–17 min of cache/hits loading, 0 fetch) + P7 (7–8) + renders
≈ 20–27 min; legacy tree ≈ 45–60 min; all six ≈ 3.5–4.5 h serial (the flock forbids parallel POST).

### 4.2 Other steps

| step | wall | memory | source |
|---|---|---|---|
| alignment per panel (not re-run) | 24,949 / 26,271 / 21,579 s (edlib 17,470 / 18,800 / 14,479 s; GPU 321 / 345 / 320 s; rescore 7,128 / 7,098 / 6,756 s), 40 workers | smoke 4.9 GB RSS | `method_selection.json`; unit 1 month 1 d CPU |
| optional stage-3 re-run from saved shortlist (latest align code, verification only) | 3,425 s (original panel) | not recorded | GPU_CAND_200.md |
| E1 genecap_tau0995 | 662 s (2026-11-11; incl. 180 fetches in 15 s); merge 75 s | maxrss 9.27 GB | `genecap_tau0995/run_meta.json`, journal |
| repaired-cache re-score (E5 driver, two passes) | 1,355 s | maxrss 9.74 GB; cgroup peak 12.3 GB | `genecap_repaired_impact/run_meta.json`, journal |
| E2 catalogue rebuild | 64.5 s | — | CORRECTIONS log `derived` |
| E3 E3 gate | 1,595 s, 12 workers | main 10.6 GB; cgroup peak 21.2 GB | `pfunc_e3_gate/run_meta.json` |
| E4 review recompute | 550 s, 8 workers | not recorded | `review_recompute.json` `elapsed_s` |
| D1 compare_three_panels | 3.6 min per suffix | not recorded | downstream.log |
| D2–D7 | ≈1 min per suffix together | not recorded | downstream.log |
| D8 / check_d100 all | 4.3 min | **22.2 GB peak** | journal `rv-d100-checkall` |
| whole downstream unit | 13.5 min (14 min 36 s CPU) | not recorded | journal |
| depth_probe_eval | not recorded (results written ~5 min after alignment, 2026-09-22) | not recorded; loads the whole cache with orjson | file mtimes |
| manuscript renders | ~1–2 min | small | out/ mtimes |

**POST memory was never recorded.** P1–P5 and P7 each `json.load` the whole 5.1 GB cache (P2/P5/P7
also hold it as Python sets via `bvbrc_gene_provider`; P5 and P7 additionally re-read it for the key set
via `load_cache_keys`, P5 once more per round to merge; `score_mag.py`'s docstring put the
111,789-genome cache at ~25 GB of string objects, it now has 159,753 genomes). P1/P2/P4/P5/P6/P7 also
load the 4.0–4.4 GB fanned hits whole and P2/P3/P5/P6/P7 the 3.2–4.0 GB selection whole. One
OOM kill is on record (P5, random legacy). The streaming re-scorers (`genecap_tau0995.py`, the gate)
stay at 9–21 GB. Today: 251 GB RAM, ~204 GB available, swap 59/59 GB used, and other user services
running (`codif-v7fpost-*`).

### 4.3 Network

BV-BRC is needed only for genomes absent from the cache. Expected: 0 (section 0, item 3). All fetch
paths now page (`fetch_patric_cds_rows`, `eq(annotation,PATRIC)`, raise `IncompleteFetch` rather than
cache a partial or empty set): `prefetch_gene_families.fetch_pgfams` (P1, P4, P5, E1, the provider's
lazy miss) and `prefetch_pgfams.fetch_pgfams_capaware` (P3). No region-validation, panel or
manuscript-figure script reads `bvbrc_cache/kbase_genome_cache` (grep), so the 15 replaced genome dicts
do not reach this report's scope.

---------------------------------------------------------------------------------------------

## 5. Stale copies of BV-BRC gene data in scope, and how to regenerate them

| copy | state | regenerate |
|---|---|---|
| `bvbrc_cache/genome_gene_families.json` | repaired 2026-11-13 (1,006 entries; log with old values `CORRECTIONS_2026-11-13.json`) | — (source of truth) |
| `bvbrc_cache/e3_clade_catalogue.json` | rebuilt 10:43 on the repaired cache (`catalogue_is_current` true) | E2 again after any later cache write |
| `region_validation_d100/data/genecap_tau0995/pgfam_supplement.json` | 180 entries fetched 2026-11-11 by the old unpaged query; 2 truncated (1609099.3, 1733.9553); merged into the cache, which is now repaired. `genecap_tau0995.py` uses a supplement entry only when the cache lacks the genome, so it is inert today | move aside before E1; E1 rebuilds it from the paged fetcher (expected empty) |
| panel indices (section 2.5) | frozen 10-08/10-09/10-18; 140 / 20 / 38 repaired genomes hold old lengths | write an extractor (none in repo): stream `$CACHE` with ijson, keep the panel's genome ids (`region_validation/data/selection_10k.json`, `kmeans_panel/selection_kmeans.json`, `random_panel/selection_random.json`), int-encode, write `{"n_vocab","genomes"}`; check every list length equals the cache's set size; then re-score MAG batches (`score_km.sh`, `score_rnd.sh`, the original n500 scoring) |
| legacy `selection.json` ×3 (and d20 original) | 933–1,674 records per tree excluded a repaired genome as `no_gene_set` | P2 |
| every `gene_capture_*` (6 d100 trees, tau0995) | scored on the pre-repair cache | P7, E1 |
| `pfunc_e3_gate*` | pre-repair store; catalogue excluded 1,802 `[]` genomes (now 963) | E3, E4 |
| `depth_probe_results.json` | pre-repair cache (P3 tie-break + EC projection) | `depth_probe_eval.py` |
| d20 trees' `gene_capture_*` (and legacy taxacc/selection for the original) | pre-repair | only if Q1 says so |
| Out of scope here (PICRUSt2/shared_reference report): `region_validation/picrust2_comparison/data/genome_function_sets.json` and its d100 copy (frozen PGFam lists; 17 truncated panel genomes), `prfba_selection*.json` (cache-only provider; the legacy variant is cache-dependent), shared_reference per-fold `_cache_pgfam_subset.json` / `_cache_ec_subset.json`. `fetch_bvbrc_annotations.py` queries the genome and protein-family tables by id batches with `limit(len(chunk)+10)` — not CDS rows, so not affected by the 25,000-row window. | | |

---------------------------------------------------------------------------------------------

## 6. Hazards

- **H1 Shared cache writers.** The cache is rewritten whole (load → mutate → `os.replace`); two writers lose
  each other's fetches. The flock (`region_validation_d100/.post_cache_writer.lock`) covers only
  `run_post_all_d100.sh` / the d100 `run_validation.sh`. Outside it: `genecap_tau0995.py --merge-supplement`,
  the provider's lazy miss inside `select_references` / `score_genecap` (a whole-file rewrite **per miss**),
  the d20 trees' scripts, `bvbrc_cache/complete_kmeans_pgfams.py`, any deployment build that fetches.
  Readers are safe (atomic replace). Check for `$CACHE.tmp` and running writers first.
- **H2 Cache change ⇒ catalogue stale.** Any write (even one gap-close fetch) changes size/mtime, so
  `pfunc_e3_panel_gate.py` refuses (`catalogue is stale … rebuild first`) and `reweight_build_e3.py`'s
  preflight flags the deployments' catalogue as stale. The BV-BRC-only builds the user wants uploaded were
  built against the 10:43 catalogue; a later cache write does not change them, but their recorded catalogue
  stamp stops being current. If any fetch happens (P1–P5, the provider's lazy miss, or E1's
  `--merge-supplement`), rebuild E2 after the last write and before E3. A 0-fetch prefetch never
  touches the file (`populate_cache` / `run_pass` / `prefetch_universe` return before saving).
- **H3 Memory.** One OOM kill (P5, 2026-10-21 12:30, while other jobs ran); swap is full today. Run one
  tree at a time (the flock does), under `systemd-run --user` with memory accounting
  (`-p MemoryAccounting=yes`, a `MemoryMax` with headroom) so the journal records `memory peak`
  — that would also give the first POST memory numbers. The streaming path (`genecap_tau0995.py`) is
  9–10 GB but exists only in `region_validation_d100/scripts`, is pinned to that tree by its own `_common`
  import, and has no `--legacy` switch (a wrapper that pre-imports another tree's `_common`/`score_genecap`
  and sets `S.LEGACY_SELECT` would work, as the E5 driver did; that is new code).
- **H4 Resumability.** POST is resumable stage by stage, but checkpoints are "output exists": stale outputs
  are silently kept unless deleted (section 3.2). `--force` is refused in d100 trees. P1/P3/P4 always run
  and reload the cache (~3 loads per tree even with 0 fetches). A failed `check_d100.py post` deletes
  `close_exclude_rep_gap.done` and `gene_capture_by_region.json` and writes no `post.done`. Tree renderers
  skip existing figures (no `--force` in the driver). Prefetch saves every 2,000 fetches; a crash loses
  at most that many (retried next run). `run_downstream_d100.sh` reruns everything (no checkpoints).
- **H5 d20 trees.** Their `run_validation.sh` would align at `--k2 100 --gpu-cand 200` if `align.done` were
  missing, scores without `--legacy`, prefetches without `--select-mode`, and still accepts `--force`
  (which re-runs `select_10k.py`). Running it would switch the d20 original from the cited legacy rule
  to the union. Re-score d20 trees, if at all, by hand with explicit flags
  (`select_references.py … --select-mode legacy`, `score_taxacc.py --legacy --force`,
  `score_genecap.py --legacy --force`, `mcnemar_regions.py --force`). The June selection was made by older
  `select_references` code (54d1150, f6b3be7, 7153311 came later), so a reproduction on unchanged inputs
  is unverified. Rewriting any of the 32 files in `d20_baseline.sha256` fails check 6 by design; the
  baseline must be re-stamped deliberately, as on 2026-10-29.
- **H6 Determinism / byte-for-byte.**
  - `score_genecap`: reproduced byte-for-byte twice (2026-11-11 tau0995 run vs the 10-20 files; 2026-11-13
    control pass on the logged old values reproduced all 7 committed files).
  - Union `select_references`: record-for-record on 20,644 / 13,268 / 13,005 touched records (census);
    whole-file byte identity not checked (JSON `indent=1`, no timestamps, dict order deterministic).
  - `score_taxacc`, `mcnemar_regions`, `diagnose_self_recovery`, `compare_panels`, `mechanism_depth`:
    no recorded reproduction; outputs carry no timestamps. Re-running union taxacc once into a scratch
    copy and `cmp`-ing is a cheap check before trusting "reuse".
  - `compare_three_panels*`: byte-identical under PYTHONHASHSEED 0/777/12345 since 5d2a204; `emptied_strata`
    has no RNG.
  - `pfunc_e3_panel_gate`: fork pool + seeded bootstrap (seed 1729, B 2000); run-to-run identity unverified.
  - Alignment: stage-3 byte-identical for the original panel only; k-means/random unchecked.
  - Cache file bytes change with insertion order on every rewrite; compare entries (or the key-set
    sha256), not the file sha256.
  - PDFs embed CreationDate: re-renders are not byte-identical (the paper copies match `out/` only because
    they were copied). Set `SOURCE_DATE_EPOCH` or compare PNG pixels.
  - `two_panel/project_space.py` UMAP is deliberately unseeded: never re-run it; reuse `space_xy.npy`.
- **H7 Built-in reproduction checks invert meaning after the repair.** `genecap_tau0995.py` compares its
  0.987 cells with `data/gene_capture_*`; the gate compares with `gene_capture_cells.csv`. Both fail on
  the repaired cache until P7 has been re-run, and must pass afterwards (use them as cross-checks; order
  P7 → E1 → E3).
- **H8** `pfunc_e3_gate_review_recompute.py` writes into the repo (`data/pfunc_e3_gate/review`) unless
  `PFUNC_GATE_REVIEW_OUT` is set.
- **H9 Alignment guard.** `run_validation.sh` (d100) requires `method_selection.json` params
  k2=100/gpu_cand=500/prefilter_k=500; re-aligning into a tree with today's defaults (gpu_cand 200) would
  fail that check. Any stage-3 verification belongs in a scratch outdir with `--shortlist-in`.
- **H10 Cross-workstream order.** D4, D6, D7, F14, F15, `three_panel_mag_panelarm.json` and
  `mag_panel_comparison.json` need the MAG cells re-scored on rebuilt panel indices; running
  `run_downstream_d100.sh` before that mixes new 16S cells with stale MAG cells. D6/D7's MAG inputs were
  inferred (RUNBOOK §5).
- **H11 Hand-carried numbers.** `presentations/deck_constants_d100.json` EMPTIED (4.7/22.9/33.4) comes from
  `emptied_strata_d100.json`, whose N (outcome != truth_empty) grows when stale truths become scorable
  (original FL truth_empty 83 → 30); PFUNC_E3_GATE.md, D100_VS_D20.md, UNION_VERIFICATION.md,
  UNION_THROUGHOUT_CHANGES.md, two_panel/RESULTS.md, `genecap_repaired_impact/impact_summary.md` and the
  manuscript quote values that move; the Welch tests (UNION_THROUGHOUT_CHANGES 1.72) have no script.
- **H12 Disk.** 2.7 TB free; each cache rewrite needs a 5.1 GB tmp file.

### Validation checks that exist

- `region_validation_d100/check_d100.py`: `align` (checks 1–4: 100 hits/record and params; top-20 prefix vs
  depth 20 with GPU-score tie analysis; 316 probe md5s = `hits_k500[:100]`; truncation report),
  `copy` (legacy hit copies), `post` (check 5: `select_mode` in selection/taxacc/genecap matches the tree,
  `n_empty_after_exclude` reported, `pgfam_missing == 0`, gap-close sentinel), `all` (+ check 6: the 32-file
  d20 baseline). Last full run 2026-10-21: all PASS.
- `run_validation.sh` guards: front-load present and not symlinked; tree/rule match (exit 2); alignment
  params; `_common.BASE` isolation; GPU-share guard; flock.
- `close_exclude_rep_gap.py` fixed point (no missing exclude-mode reps); never caches `[]`.
- `genecap_tau0995.py` 0.987 byte-identity vs `data/gene_capture_*`; `pfunc_e3_panel_gate.py`
  `reproduction_vs_gene_capture_cells` (135,982 cells matched) and stale-catalogue refusal.
- The 2026-11-13 controls (scratch `bvbrc_trunc/census_*.json`, `genecap_repaired_impact/run_meta.json`):
  old-value reproduction and per-cell attribution — reusable as the before/after harness.
- `tests/` (pytest: `test_bvbrc_paging.py`, `test_select_method.py`, `test_pfunc_estimator.py`, …).
- Expected after the rescore (from the impact measurement, original union): exclude_0.987 population
  recall 60.17 → 61.04 %, resolved 67.72 → 67.70 %, truth_empty 953 → 278 (include and both exclude
  modes), pred_empty 1,060 → 243, abstained 6,407 → 6,423, include rows changed 1,389.

---------------------------------------------------------------------------------------------

## 7. Open questions for the lead

- **Q1** Re-score the d20 trees? Their gene-capture and legacy-taxonomy numbers are cited only in
  d20-vs-d100 comparisons and in the 32-file baseline. Cost: original legacy ≈ 45–60 min by hand plus
  re-stamping `d20_baseline.sha256`; k-means/random ≈ 10 min each (genecap only). Risk: the June
  `select_references` is not today's, so a re-selection may differ for reasons other than the cache.
- **Q2** Re-run union `score_taxacc` anyway (≈ 5 min/tree) to prove byte reproduction, or rely on the
  census argument?
- **Q3** Who owns the panel-index extractor and the MAG re-score (D4/D6/D7/F14/F15 wait on it)?
- **Q4** If a legacy gap-close fetch adds a genome to the cache, accept the catalogue rebuild (and the
  deployments' catalogue stamp going stale), or route such fetches to a side file?
