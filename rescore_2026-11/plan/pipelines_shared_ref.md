# Rescore plan: `region_validation/shared_reference` and `region_validation/picrust2_comparison`

Read-only investigation, 2026-11-15, prFBA HEAD `cc8014f` (branch
`hardware-aware-alignment-dispatch`). Nothing in the repo was edited, run or written. Every
number marked "measured" comes from a read-only script I ran; scripts and outputs are in
`/tmp/claude-1000/-home-freiburger-Documents-prFBA/400ec395-f172-59da-a98d-3dbe6b3a4256/scratchpad/rescore/scripts/`
and `/tmp/claude-1000/-home-freiburger-Documents-prFBA/400ec395-f172-59da-a98d-3dbe6b3a4256/scratchpad/rescore/work/`
(my files there: `repaired.json`, `gfs_staleness*.json`, `changed_sets.json`,
`new_candidates_free_filters.json`, `fold_pgfam_staleness.json`, `universe_check.json`,
`unknown_pgfam_classes.json`; `panel_len_check*` and `tieset_check*` belong to another worker).
Paths below are relative to `~/Documents/prFBA` unless absolute. `SR=region_validation/shared_reference`,
`P2C=region_validation/picrust2_comparison`.

---

## 0. Summary

1. **Both studies read BV-BRC gene data through one file:**
   `P2C/data/genome_function_sets.json`, built 2026-09-17 15:22 from the PGFam cache as it was
   then. The shared-reference campaign never reads the PGFam cache for truth or traits. It
   reads it only for gene counts in prFBA selection, through per-fold `_cache_pgfam_subset.json`.
   The as-shipped benchmark reads the cache directly in `prfba_selection.py`,
   `fetch_bvbrc_annotations.py` and `genome_function_sets.py`.
   **The kbase genome cache and the E3 clade catalogue are read by nothing in either study**
   (grep of both `scripts/` trees).
2. **The d20 `genome_function_sets.json` is stale for 713 of its 49,109 genomes** (measured:
   the current cache re-projected through the stored `function_maps.json`).
   - 693 are 2026-11-13 repairs: 614 empty-list poison, 40 PATRIC-cap truncations, 39 RefSeq
     truncations.
   - **20 are a second, older stale class: the 2026-10-09 poisoned-key repair.** All 20 are in
     `bvbrc_cache/poisoned_keys.json`, and all had `n_pgfam = 0` in the 09-17 build. The task
     text does not name this class; none of the 20 is a 10k-panel genome.
   - The d100 copy (`region_validation_d100/picrust2_comparison/data/genome_function_sets.json`,
     built 10-20) is stale for 733 of 59,697 genomes, all from the 11-13 class (622 / 66 / 45).
3. **The repair reaches the shared-reference design, not only its scores.**
   - 140 of the 10,000 panel source genomes changed.
   - 17 of them are eligible (16 Bacteria, 1 Archaea). Each gains 9–114 ECs and loses none,
     counting only PGFams that can be projected now (measured).
   - **123 of the 161 genomes dropped as `empty_ec` now have ECs.** 122 of those also pass the
     length and charset filters: 65 Archaea (48 Methanobacteriota) and 57 Bacteria, and 104 of
     them share a genus with an eligible genome. HMM occupancy, the last filter, was not
     computed. At the panel's measured 5.7% occupancy drop rate, expect about 115 new eligible
     genomes.
   - Re-running `plan_folds.py --phase panel` on regenerated data therefore changes
     `n_eligible`, from 8,872 to about 8,990. Through the genus edges it changes the component
     structure and the k=5 bin packing, and so every fold. Section 4 costs both paths: Path A
     keeps the panel and fold plan; Path B re-plans.
4. **Path A reaches the following** (measured):
   - 23 of the 30 genus fold references carry a changed genome. All 45 per-fold
     `_cache_pgfam_subset.json` files (data, data_family, data_c2) differ from the current cache
     in exactly these 17 genomes.
   - 130 test amplicon cells (fold 0: 24, 1: 47, 2: 24, 3: 31, 4: 4) have changed truth.
   - 3 of the 1,000 ATP-pass phenotype genomes and 0 of the 200 gapfilled-pass genomes are
     among the 17.
   - No EC that can be projected now falls outside the 2,235-EC `panel_ec_universe.json`, for
     either the 17 or the 122. PGFams missing from `pgfam_products.json` cannot be projected
     until they are fetched (next point).
   - **These counts hold in EC space only.** If `function_maps.json` is regenerated against
     the moved ModelSEED trees (point 7), `ec_rxn` and `pgfam_rxn` change for every genome, so
     reaction-space truth moves for every cell, not only these 130.
     `function_maps.json:_meta` records paths, not commits, so the 09-17 ModelSEED state cannot
     be reconstructed. The closest state that can be pinned is the d100 build (10-20), which
     differs from the d20 maps in 2 `ec_rxn` values.
5. **A network step is unavoidable.** The repaired genomes carry PGFam ids that are not in
   `P2C/data/pgfam_products.json` (measured below), so `fetch_bvbrc_annotations.py pgfams`
   (BV-BRC `protein_family_ref`, resumable, incremental) must run before `function_maps.py` and
   `genome_function_sets.py`. Without it, `genome_function_sets.py` silently drops those
   families' ECs (`pgfam_ec.get(p, ())`).

   | Genome class | Genomes carrying unprojectable PGFams | Distinct PGFams |
   |---|---|---|
   | sources | 120 | 15,975 |
   | default-selection reps | 484 | 25,896 |
   | legacy-only reps | 5,442 of 6,508 | 352,922 |

   **The legacy-only row is a defect present since the 09-17 build.** The fetcher reads only
   `--selection prfba_selection.json`, but `genome_function_sets.py` also projects the reps in
   `prfba_selection_legacy.json`. The d20 `prfba_*_legacy` as-shipped arms were therefore
   under-projected from the start. The d100 copy does not have this gap: its 588 genomes with
   unknown PGFams are all repaired genomes.
6. **Everything downstream reuses its own cached copy whenever the file exists, without
   checking where it came from.** These copies are `_cache_panel_ec.json`, the per-fold
   `_cache_ec_subset.json` and `_cache_pgfam_subset.json`, `kfold_audit/panel_ec_sets.json` and
   `mag_arm/panel_ec.json`. Regenerating `genome_function_sets.json` alone changes no
   shared-reference number. The copies must be deleted in the order given in section 5. One
   ordering trap: if `data/_cache_panel_ec.json` is deleted before it is regenerated,
   `build_reference.py` silently falls back to the stale
   `P2C/data/kfold_audit/panel_ec_sets.json`.
7. **The "latest code" has a tension.** The aligner changed after the campaign's commit
   `b72b71cf`, in `dbbc8f7`, `d0a926f` and `15f3caa`.
   - Those commits report byte-identical records at `--gpu-cand >= --k2`. But `d0a926f` also
     replaced an unstable argsort pick that "would have changed 2,180 records at 200".
   - Reusing `hits/` keeps 10-14 alignments, which the align stamp accepts because it hashes
     inputs, not code.
   - Re-aligning is the campaign's longest step, at about 4,000 s per prFBA fold run, two at a
     time.
   - The ModelSEED inputs of `function_maps.py` also moved: ModelSEEDTemplates HEAD `e7fc410`
     (10-21, rxn39860 removed) and ModelSEEDDatabase reaction records (09-11 to 09-29). The d20
     and d100 maps already differ in 2 `ec_rxn` values. Regenerating `function_maps.json`
     therefore moves reaction-space numbers for reasons unrelated to the cache.

---

## 1. Where BV-BRC gene data enters

| BV-BRC artefact | First reader in these studies | Derived copies (all must be regenerated or deleted) |
|---|---|---|
| `bvbrc_cache/genome_gene_families.json` (shared PGFam cache; 5.12 GB; repaired 11-13 10:38) | `P2C/scripts/fetch_bvbrc_annotations.py pgfams` (PGFam id set), `P2C/scripts/genome_function_sets.py` (projection), `P2C/scripts/prfba_selection.py` (P3 gene-count tie-break), `P2C/scripts/depth_probe_eval.py` (gene counts), `SR/scripts/run_prfba_fold.py` (`--pgfam-cache`, gene counts) | `SR/data*/prfba/*/_cache_pgfam_subset.json` (45 files). Indirect: everything that reads `genome_function_sets.json` below |
| BV-BRC `protein_family_ref` (live API) | `fetch_bvbrc_annotations.py pgfams` → `P2C/data/pgfam_products.json` (09-17 13:44, 4,174,289 families) | `P2C/data/function_maps.json` (09-17 15:12) |
| BV-BRC `genome` (live API) | `fetch_bvbrc_annotations.py genomes` → `P2C/data/bvbrc_source_genomes.json` (assembly accessions; unaffected by the repair) | — |
| `P2C/data/genome_function_sets.json` (1.71 GB, 09-17 15:22, sha256 `6f4f817c…`) | `SR/scripts/plan_folds.py`, `build_reference.py`, `run_prfba_fold.py`, `score_folds.py`, `phenotype_fold.py`, `calibrate_pfunc.py`, `design_audit_baselines.py`; `P2C/scripts/compare.py`, `annotation_ceiling.py`, `annotation_dialect_test.py`, `mag_common.py` | `SR/data/_cache_panel_ec.json` and `SR/data_family/_cache_panel_ec.json` (10-13/10-14); `SR/data/panel_ec_universe.json` (hard link of `SR/data_family/panel_ec_universe.json`); `SR/data/panel_eligible.json` (`n_ec`, `dropped`; hard link of data_family's); 30+30 `refs/*/traits_EC.txt.gz` + `manifest.json` (+2 `data_c2/refs/incself_*`); 45 `prfba/*/_cache_ec_subset.json`; `P2C/data/kfold_audit/panel_ec_sets.json` (10-13, fallback source for `build_reference.py`); `P2C/data/mag_arm/panel_ec.json` (10-16) |
| `mag_validation/pgfam_panel.json` (10-08; predates both repairs; out of this scope) | `P2C/scripts/mag_prfba_ec.py` via `mag_validation/score_mag.load_pgfams` (gene provider for the MAG arm's selection) | `P2C/data/mag_arm/prfba_mag_selections.csv.gz` |
| `bvbrc_cache/kbase_genome_cache/`, `bvbrc_cache/e3_clade_catalogue.json` | **not read by either study** | — |

`score_folds.py` treats `--functions`, `--maps`, `--ec-universe`, `--fanout`,
`--annotation-ceiling` and `--templates` as the ANSWER KEY. Their bytes are fingerprinted
(`answer_key_fp`, currently `a9dc024df800`). Any regeneration therefore changes `answer_key_fp`
and `inventory_fp`, and voids the reuse of `cells.csv.gz`.

---

## 2. Measured staleness (read-only checks)

| Check | Result |
|---|---|
| d20 `genome_function_sets.json` vs current cache re-projected through stored `function_maps.json` (`gfs_staleness.json`, 108 s, about 1.5 GB RSS) | 48,396 identical, 713 differ. 614 `empty_list_poison`, 40 `truncated_by_patric_cap`, 39 `truncated_by_refseq` (all 11-13), and 20 from the 10-09 poisoned-key repair. In every one of the 713, EC, `rxn_ec` and `rxn_role` all change. No genome is missing from the cache |
| d100 copy, same test (`gfs_staleness_d100.json`) | 58,964 identical, 733 differ (622 / 66 / 45, all 11-13); 679 genomes overlap with the d20 diff set |
| Panel genomes among the changed (`changed_sets.json`) | 140 of 10,000 panel source genomes: 17 eligible, 123 `empty_ec` drops. All 140 are truth sources of the as-shipped benchmark |
| The 123 `empty_ec` drops (`new_candidates_free_filters.json`) | 122 pass length [1200, 1800] and the ACGTN filter, 1 fails length. 65 Archaea, 57 Bacteria; 104 share a genus with an eligible genome. Occupancy not computed. 38 `empty_ec` drops stay empty |
| The 17 eligible genomes (`universe_check.json`) | Each gains between 9 and 114 ECs from projectable PGFams and loses 0. 813.301 has 8,640 → 14,057 PGFams, 2,256 of them not yet in `pgfam_products.json`. 0 projectable ECs fall outside the 2,235-EC universe, for the 17 and for the 122 |
| Per-fold PGFam subsets vs current cache (`fold_pgfam_staleness.json`, 76 s) | All 45 `_cache_pgfam_subset.json` (data 20, data_family 20, data_c2 5) differ in exactly the same 17 genomes (from 1–17 per subset) and nowhere else. The 10-14 subsets already carried the 10-09 repair |
| Where the 17 sit (genus campaign) | 23 of 30 refs (unaffected: f02/f03/f04 arc d025/d050 and f04_arc_d100). Test genomes per fold 3/6/3/4/1 → 24/47/24/31/4 amplicon cells. Phenotype: 3 of 1,000 ATP-pass genomes, 0 of 200 gapfilled-pass genomes |
| `annotation_ceiling.json` | 1 of its 2,689 per-genome rows is a changed source. 31 genomes were skipped for empty truth, and some may now be scoreable |
| Unprojectable PGFams by class (`unknown_pgfam_classes.json`, 40 s) | Sources: 120 genomes, 15,975 PGFams. Default reps: 484 genomes, 25,896. Legacy-only reps: 5,442 of 6,508 genomes, 352,922. Distinct total: 389,145 |
| `function_maps.py` inputs | ModelSEEDTemplates shared tree at `e7fc410` (10-21, rxn39860 removed). ModelSEEDDatabase commits on 09-11, 09-15, 09-24, 09-25, 09-28 and 09-29 touch reactions. d20 vs d100 maps: `pgfam_ec` and `pgfam_rxn` identical on common keys, `ec_rxn` differs in 2 values |

---

## 3. Command sequences that produced each committed result file

Conventions: `PY=~/Documents/py_venv/bin/python`; `SR=$HOME/Documents/prFBA/region_validation/shared_reference`.
PICRUSt2 is called by the scripts through `bash -lc 'source ~/opt/picrust2/activate.sh && …'`.
Runtimes come from `SR/data/ledger.jsonl` (`wall_s`), the logs named, or `JOBS_AUDIT_2026-11-02.md`.
"Unmeasured" means no log records it.

### 3A. As-shipped PICRUSt2 benchmark, d20 (`P2C`, README "Pipeline", commit `8dcc58c`, 09-17)

Several committed files feed the manuscript and the shared-reference study:
`comparison_summary.json` (manuscript line 733, "0.600"), `annotation_ceiling.json` (0.904/0.743,
n = 2,689; also an answer-key input of `score_folds.py`) and `depth_probe/depth_probe_results.json`
(manuscript lines 722–728). The manuscript's as-shipped section itself
(0.802 vs 0.727, coverage 88.9%) comes from the **d100 copy**
(`region_validation_d100/picrust2_comparison`), which is outside this scope but has the same
script chain.

| # | Command (cwd `P2C`) | Output | BV-BRC data? | Runtime / memory | Network | Path A or B |
|---|---|---|---|---|---|---|
| a1 | `python scripts/prepare_inputs.py --outdir data/full_inputs --batches 8` | `data/full_inputs/` | no | unmeasured (seconds) | no | **reuse** |
| a2 | `bash scripts/run_picrust2.sh data/full_inputs data/full_picrust2 2 30` | `data/full_picrust2/` (stock PICRUSt2 2.6.3) | no | 13:27→15:09 = 1 h 42 m (`full_picrust2_driver.log`) | no | **reuse** |
| a3 (run before a5: the fetcher takes its genome set from the selection on disk; a3 reads only hits, truth and cache, none of which the rescore regenerates) | `python scripts/prfba_selection.py --out data/prfba_selection.json --workers 24` and `… --legacy --out data/prfba_selection_legacy.json` | selections, include + exclude_0.987, on `region_validation/data/hits/asv_top20_alignment_hits.json` (depth 20, 06-10) | **yes**: whole PGFam cache via orjson, for P3 gene counts | 74 s wall, **25.8 GB max RSS** when the cache was 111,789 genomes (`prfba_selection.log`). Now 159,753 genomes / 5.12 GB: expect about 40 GB (unmeasured) | no | **re-run** (gene counts changed; a former `[]` genome can now win a tie) |
| a4 | `python scripts/fetch_bvbrc_annotations.py genomes` | `data/bvbrc_source_genomes.json` | no gene data | unmeasured | **yes** | **reuse** |
| a5 | `python scripts/fetch_bvbrc_annotations.py pgfams` | `data/pgfam_products.json` (resumable: fetches only missing ids) | **yes**: loads the whole cache (orjson) to collect ids | 958 s for 4.17M families at 5,000 per batch (`fetch_pgfams.log`). Incremental: about 42k (sources + default reps) or 395k (with legacy) ids, roughly 10–90 s of API time plus the cache load (about 40 GB, unmeasured) | **yes** | **re-run after a3**, once per selection file: `fetch_bvbrc_annotations.py pgfams` (default) and `fetch_bvbrc_annotations.py pgfams --selection data/prfba_selection_legacy.json` (hazard H8) |
| a6 | `python scripts/function_maps.py` | `data/function_maps.json` | products only | unmeasured (about 1 min) | no | **re-run** (new PGFams). Reads the moving ModelSEEDTemplates `v7.0` and ModelSEEDDatabase: hazard H9 |
| a7 | `python scripts/genome_function_sets.py` | `data/genome_function_sets.json` | **yes**: whole cache (orjson) + maps + both selections | unmeasured; about 40 GB expected | no | **re-run** |
| a8 | `python scripts/annotation_ceiling.py` | `data/annotation_ceiling.json` (committed) | via a7 | unmeasured | no | **re-run** |
| a9 | `python scripts/compare.py` | `data/comparison_cells.csv.gz`, `data/comparison_summary.json` (committed), `data/domain_routing.csv` (committed; PICRUSt2-only, unchanged) | via a3, a6, a7 | unmeasured | no | **re-run** |
| a10 | `python scripts/render_comparison.py` | `figures/F1_prfba_vs_picrust2.{png,pdf}` (committed) | via a9 | seconds | no | **re-run** |
| a11 | `python scripts/annotation_dialect_test.py` | `data/annotation_dialect.json` (committed; the manuscript uses the d100 copy's) | via a3, a7 | unmeasured | no | **re-run** |
| a12 | `python scripts/depth_probe_sample.py` then the k500 re-alignment (README: `--k2 500 --gpu-cand 500`, about 4 min) then `python scripts/depth_probe_eval.py` | `data/depth_probe/{wiped_sample.keys.json, depth_probe_results.json}` (committed) | `depth_probe_eval.py` loads the whole PGFam cache (orjson) + `function_maps.json` | unmeasured | no | **k500 hits: reuse.** **Sample: re-run** `depth_probe_sample.py` (seconds; it reads `--selection prfba_selection.json`, which a3 regenerates) and diff the result against the committed `wiped_sample.keys.json`. If the sample is identical, the k500 hits stand; if not, re-align the new keys (about 4 min). **eval: re-run** |
| a13 | MAG/16S study (commit `1c03234`): `mag_16s_loci.py`, `mag_extract_16s.py`, `mag_16s_retention.py`, `mag_gtdb_ssu.py`, `src16s_hmm_check.py`, `mag_direct_feed.py`, `mag_prfba_ec.py`, `mag_score.py` | `data/mag_*_summary.json`, `data/src16s_hmm_check_summary.json` (committed; read by the deck, not the manuscript) | Truth: `mag_common.load_panel_ec()` reuses `data/mag_arm/panel_ec.json` if it exists, else derives it from a7. `mag_prfba_ec.py` reads `mag_validation/pgfam_panel.json` | PICRUSt2 panel run 12:03→12:49 (`mag_arm/logs/panel_picrust2.log`), sample run 11:48→12:09; the rest unmeasured | `mag_gtdb_ssu.py` reads local GTDB r232 metadata | **reuse**: loci, retention, extraction, PICRUSt2 runs on extracted 16S, GTDB SSU, hmm check (no gene data). **re-run**: `mag_direct_feed.py` (its F1 uses the truth), `mag_prfba_ec.py` (only if `pgfam_panel.json` is regenerated upstream), `mag_score.py`. **Delete `data/mag_arm/panel_ec.json` first** |

### 3B. Shared-reference leave-genus-out campaign (`SR/data`, commit `b72b71cf`, 10-13 to 10-14)

The driver run is recorded in `SR/data/ledger.jsonl` (`run_campaign.py --root data --k 5
--folds 5 --build-jobs 12`; the first launch with k = 10 aborted at plan, exit 2, because the
largest component, 2,084 genomes, is more than 1.5·n/k). Committed results:
`data/panel_eligible.json`, `data/fold_plan.json`, `data/summary_genus.json`,
`data/controls.json`, `data/summary_recipe.json`, the phenotype summaries and the post-hoc
studies in 3D–3F.

| Stage | Command (as recorded) | Outputs | Reads BV-BRC data? | Measured runtime | Network | Path A | Path B |
|---|---|---|---|---|---|---|---|
| 1 panel | `$PY $SR/scripts/plan_folds.py --phase panel --outdir $SR/data --threads 32` | `panel_eligible.json`, `panel_{bac,arc}.{fna,afa}`, `panel_md5_seq.json`, `panel_md5_id.json`, `panel_ec_universe.json`, `panel_16s_queries.fna`, `midas_empty.csv`, caches `_cache_panel_ec.json` and `_cache_panel_md5_seq.json`, `_hmmalign/` (content-addressed) | **yes**: streams `genome_function_sets.json` into `_cache_panel_ec.json`; eligibility filters 1–2 | 190 s | no | **There is no CLI for this.** `--phase panel` re-decides eligibility. Path A needs a one-off that regenerates `_cache_panel_ec.json` only, e.g. `plan_folds._cache_panel_ec(functions, panel_ids, cache, force=True)` over the 10,000 panel ids. **`panel_ec_universe.json` and `panel_eligible.json` can stay byte-identical:** the EC list is unchanged as measured, `score_folds.load_ec_universe`'s prevalence return value is discarded (`_prev`, line 4259), nothing else reads `prevalence`, and nothing reads `panel_eligible.json:n_ec` (grep of `SR/scripts`). Keeping them identical keeps their sha256 in the phenotype binding and in `answer_key_fp`. The 123 `dropped: empty_ec` reasons and the 17 `n_ec` values stay stale but are cosmetic | **re-run** |
| 2 prepass | `$PY $SR/scripts/plan_folds.py --phase prepass --outdir $SR/data --k2 100 --gpu-cand 100 --force-backend cuda --threads 32` | `prepass_full16s.json` (521 MB), `prepass_amplicon.json` (6.7 GB), `neighbors.json`, `hits_full16s/`, `hits_amplicon/` | no (16S only) | **19,231 s (5.3 h)** | no | **reuse** | **re-run** (fingerprint-gated on the panel md5 set, so it rebuilds itself) |
| 3 plan | `$PY $SR/scripts/plan_folds.py --phase plan --outdir $SR/data --k 5 --rank genus` | `fold_plan.json` | no (genus, md5, identity edges) | 12 s | no | **reuse** | **re-run** (new components and new packing) |
| 4 refs | 30× `$PY $SR/scripts/build_reference.py --mode fold --fold F --domain {bac,arc} --density {0.25,0.5,1.0} --fold-plan … --panel … --outdir $SR/data/refs --threads 3` (12 concurrent, one serial warm-up) | `refs/fNN_dom_dDDD/{.fna,.hmm,.tre,.model,.probe.tre,traits_EC.txt.gz,traits_16S.txt.gz,manifest.json}` | **yes**: EC sets from `_cache_panel_ec.json`, falling back to `panel_ec_sets.json` and then the stale `kfold_audit/panel_ec_sets.json` | 30 ok jobs, 8,634 s summed; median 148 s; bac_d100 762–940 s (FastTree is the long pole). Wall 04:09→04:40 including the polytomy retry | no | **Partial.** Delete `manifest.json` (and `traits_EC.txt.gz`) in the **23 affected dirs only**. Keep `.fna`, `.hmm`, `.tre`, `.fasttree.log`, `.fasttree.stderr` and `.probe.tre`, and re-run **without `--force`**: FastTree, hmmbuild and the probe placement are reused, and only the crash screen (`hsp.py` on the probe tree) and the manifest are redone (cost unmeasured, about a minute each). The new manifests record `fasttree_seconds: 0.0` on this path (cosmetic; it shows in the C13 descriptives). Do not refresh the 7 unaffected arc refs "for uniformity": `write_traits` gzips with the current mtime, so their `traits_EC.txt.gz` bytes change even with identical content, their `stage_place.json` stamps stop matching, and their 7 PICRUSt2 runs would have to re-run. `--force` would redo hmmalign, FastTree and the probe placement | **re-run all** (`--force` or a fresh root) |
| 5 PICRUSt2 | 30× `run_picrust2_fold.py --mode fold --routing forced --fold F --domain dom --density d --ref-dir … --fold-plan … --panel … --outdir $SR/data/picrust2/fNN_dom_dDDD --processes 14` and 5× `--routing nsti --density 1.0 --ref-dir …bac_d100 --ref-dir …arc_d100 --outdir …/fNN_nsti` (4 concurrent) | `predictions.json`, `combined_EC_predicted.tsv.gz`, `placed.tre`, `stage_place.json`, `stage_hsp.json` | via `traits_EC.txt.gz` | 35 runs, 4,408 s summed, median 159 s, max 344 s; wall 04:40→05:02 | no | **Re-run every run whose reference's traits changed** (23 forced + 5 nsti; the 7 unaffected arc runs can stay, since their stamps still match). `stage_place.json` hashes `traits_sha256`, so placement re-runs too although epa-ng never reads traits (a false dependency; cheap here) | **re-run all** |
| 6 prFBA | 15× `run_prfba_fold.py --mode fold --fold F --density d --fold-plan … --panel … --ref-manifest …bac… --ref-manifest …arc… --outdir $SR/data/prfba/fNN_dDDD --select-modes legacy,exact_tie,graduated,exact_tie_guards --k2 200 --gpu-cand 400 --prefilter-k 1500 --force-backend cuda --workers 30` (2 concurrent), then 5× the same at d100 with `--forced-answer --hits-from …/fNN_d100 --outdir …/fNN_d100_forced` | `ref_db.json`, `ref_id_db.json`, `queries.fna`, `hits/asv_top20_alignment_hits.json`, `align_stamp.json`, `_cache_ec_subset.json`, `_cache_pgfam_subset.json`, `predictions.json` (549 MB at f00_d100) | **yes**: `_cache_ec_subset.json` from `genome_function_sets.json`; `_cache_pgfam_subset.json` from the PGFam cache (gene counts, about 90k `gene_count_fn` calls per run, so it is live, not inert) | Base runs 3,102–4,396 s each (alignment 4,073 s of 4,198 at f00_d100; the two streaming loads 12 s + 4 s; selection about 100 s). Forced runs 156–186 s with `--hits-from`. 23 ok runs, 53,186 s summed; wall 06:15→12:50 | no (cache only) | **Re-select on reused alignments**: per run dir delete `predictions.json`, `_cache_ec_subset.json` and `_cache_pgfam_subset.json`, keep `hits/` and `align_stamp.json`, and run the same command **without `--force`**, which takes the `[skip] reusing hits` branch. About 3 min per run, 20 runs | **re-align all** (3,100–4,400 s per base run, 2 GPU-concurrent: about 7–8 h) |
| 7a C1 recipe control | `build_reference.py --mode recipe_control --n-stock 960 --n-holdout 500 …`; 3× `run_picrust2_fold.py --mode recipe_control --ref-dir data/refs/recipe_{A,B,C} --queries data/refs/recipe_control_queries.fna --truth data/refs/recipe_control_truth.json --traits ~/opt/picrust2/…/bacteria/ec.txt.gz --outdir data/picrust2/recipe_X --processes 12`; `score_folds.py --mode recipe_control … --cells data/cells_recipe.csv.gz --summary data/summary_recipe.json --controls data/controls.json` (`logs/c1_*.log`) | `summary_recipe.json`, the C1 entry in `controls.json` | **no** (stock PICRUSt2 genomes and traits) | builds about 1 min; place 8 s per arm | no | **reuse** (optionally re-score so C1's provenance carries the new `answer_key_fp`) | reuse (independent of the panel) |
| 7b score, pass 1 | `score_folds.py --mode fold --fold-plan … --panel … --picrust2 data/picrust2 --prfba data/prfba --refs data/refs --neighbors data/neighbors.json --cells data/cells.csv.gz --summary data/summary_genus.json --controls data/controls.json --rank-blocked genus --bootstrap 2000 --expect-densities 0.25,0.5,1.0 --expect-domains bac,arc --expect-nsti primary --expect-forced-answer primary --force` (driver-printed; `score_genus.log`) | `cells.csv.gz` (+ `.fp`; 4,216,740 rows), `summary_genus.json`, `controls.json` | **yes**: full `jload` of `--functions` and `--maps` | **5,596 s** (`score_genus.log`), 5,658 s (`score_genus_pheno.log`). RSS not measured in any log (a full parse of the 1.7 GB file) | no | **re-run** | **re-run** |
| 7c phenotype | see 3D | | | | | | |
| 7d score, pass 2 (bound) | `scripts/c13_finish.sh` (A-C13) | `summary_genus.json`, `controls.json`; `.a-c12` backups; `c13_compare_*.json`, `c13_vs_prereg_compare_*.json` | as 7b | 4,975 s in the reuse path (`logs/score_genus_c13.log`); A-C12 4,572 s | no | **re-run** after 3D | **re-run** |

### 3C. C2 include-self (`SR/data_c2`, `scripts/run_c2_include_self.sh`, commit `6b62212`)

Steps:
1. `build_reference.py --mode include_self --domain {arc,bac} --density 1.0 … --outdir data_c2/refs --threads 16`
2. 10× `run_picrust2_fold.py --mode include_self --routing forced --fold F --domain dom --ref-dir data_c2/refs/incself_dom_d100 … --outdir data_c2/picrust2/incself_fNN_dom_d100 --processes 14`
3. 5× `run_prfba_fold.py --mode include_self … --ref-manifest data_c2/refs/incself_{bac,arc}_d100/manifest.json --outdir data_c2/prfba/incself_fNN_d100 --select-modes … --k2 200 --gpu-cand 400 …`
4. `score_folds.py --mode include_self … --cells data_c2/cells_incself.csv.gz --summary data_c2/summary_incself.json --controls data/controls.json --bootstrap 2000 --force`

Outputs are the committed `data_c2/summary_incself.json` and the C2 entry merged into
`data/controls.json` (manuscript line 649: 52,979 cells, 0.899).

Runtime: prFBA about 3,349 s alignment and 3,472 s total per fold (5 folds); scoring 3,059 s.

**Path A:**
- Refresh the two incself refs' traits (delete `manifest.json`; keep the tree).
- Re-run the 10 PICRUSt2 runs.
- Re-select the 5 prFBA runs on their kept `hits/`. Delete `predictions.json` and both
  `_cache_*_subset.json` first: the script skips any run whose `predictions.json` exists.
- Re-score.

**Path B:** everything.

### 3D. Phenotype arm (`SR/data/phenotype_*`; SPEC §9, §14, §15)

Committed: `phenotype_{atp,gf}/` (pre-registered, 10-18), `phenotype_{atp,gf}_c12fix/` (A-C12),
`phenotype_{atp,gf}_c13/` (A-C13, current), `data/phenotype_summary.json`,
`c12fix_compare_*.json` and `c13_*compare_*.json`.

Current (A-C13) commands, reconstructed from `meta`: `script`, `arms_requested`, `readouts`,
`templates`, `libraries` and `integrate_reversals`. The literal argv is not logged.

```
PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-c2e31a0 $PY $SR/scripts/phenotype_fold.py \
  --cells $SR/data/cells.csv.gz --fold-plan $SR/data/fold_plan.json --panel $SR/data/panel_eligible.json \
  --prfba $SR/data/prfba --picrust2 $SR/data/picrust2 --select-mode exact_tie \
  --templates /home/freiburger/Documents/ModelSEEDTemplates-c12run/templates/v7.0/ \
  --integrate-reversals --expect-modelseedpy c2e31a0 --maxtasksperchild N \
  --n-genomes 1000 --readouts atp,structure --workers 10 --outdir $SR/data/phenotype_atp_c13
# gapfilled: --n-genomes 200 --readouts atp,growth,essentiality,structure --workers 16 --outdir $SR/data/phenotype_gf_c13
```

The templates checkout `ModelSEEDTemplates-c12run` is at `21a8084`.

Runtime and memory:
- ATP pass: 10-29 10:15→12:27, **2 h 12 m**, 10 workers.
- Gapfilled pass: 10:15→13:09, **2 h 55 m**, 16 workers. Per genome 78–90 s (ATP) and about
  750 s (gapfilled).
- An earlier, unpinned A-C12 gapfilled pass was **oom-killed at 125.0 GB peak** (JOBS_AUDIT).
  The parent `jload`s all of `genome_function_sets.json`.

BV-BRC data: through `--functions` (truth EC sets), `cells.csv.gz` and the predictions.

**Re-run in both paths:**
- `meta.binding` stores the sha256 of `cells.csv.gz`, `panel_eligible.json`,
  `panel_ec_universe.json`, `fold_plan.json` and `--functions`. Any change voids C11/C12 under
  `--require-phenotype-binding`.
- The genome sample is drawn from the cells: per-genome ΔF1 quintiles, prFBA answered,
  `n_truth_ec >= 200`. The 1,000 and 200 genomes can therefore change even under Path A.

**The pre-registered and A-C12 passes are historical records.** SPEC §14 and §15 keep them as
the record. Re-running them is a choice, not a requirement; if they are re-run, use their own
pins (`ModelSEEDpy-d884a50`; templates `2294299` for pre-registered, `21a8084` for A-C12).

### 3E. Leave-family-out (`SR/data_family`, `chain_family.sh` + `score_family.sh`)

Stages 1–2 were not run here. `data_family/` holds **hard links** to `data/`'s
`panel_eligible.json`, `panel_ec_universe.json`, `neighbors.json`, `prepass_*.json` and
`hits_amplicon/asv_top20_alignment_hits.json` (link count 4). `_cache_panel_ec.json` is a
separate copy dated 10-14 15:08.

Steps:
1. `run_campaign.py --root data_family --rank family --folds 5 --k 5 --build-jobs 12 --stages 3`
   writes `fold_plan_family.json` (11 s).
2. Stage 4: 30 builds, 13,201 s summed, max 1,216 s.
3. `chain_family.sh`, which runs `run_campaign.py … --stages 5-6 --p2-jobs 3 --p2-processes 12 --prfba-jobs 1`:
   - PICRUSt2: 35 runs, 1,525 s summed.
   - prFBA: 20 runs, 10,478 s summed. A base run takes about 600 s (564 s alignment); a forced
     run about 45 s.
4. `score_family.sh` (452 s).

Committed: `data_family/{summary_family.json, controls.json, fold_plan_family.json, panel_eligible.json}`.

**Path A:** the same partial refresh as genus (refs traits, PICRUSt2, prFBA re-select, score).
First delete and regenerate `data_family/_cache_panel_ec.json`, which build_reference reads
beside `--panel`.

**Path B:** everything, after re-linking (or copying) the new stage 1–2 artefacts from `data/`.
`plan_folds.write_json` writes through tmp + `os.replace`, which breaks the hard link (hazard H6).

### 3F. Post-hoc studies on the genus campaign (all read `SR/data` outputs)

| Study | Command | Reads | Measured runtime / memory | Re-run? |
|---|---|---|---|---|
| p(function) calibration (10-20, `b9fa5c7`) | `$PY scripts/calibrate_pfunc.py` (→ `data/calibration/{README.md,assertions.json,reliability.csv,scores.json}`) | `cells.csv.gz`, `fold_plan.json`, `panel_*`, `_cache_panel_ec.json`, per-fold `_cache_pgfam_subset.json` (pgfam unit), `P2C/data/{function_maps,pgfam_products,prfba_selection,genome_function_sets}.json`, `comparison_cells.csv.gz` (as-shipped secondary) | unmeasured (units `prfba-calib-pfunc{,2}`, `calib-fix-full`) | **re-run** after 3A and 3B |
| Weighted-objective test (`c91cf19`) | `bash scripts/launch_wobj.sh` (systemd-run `wobj-run`, MemoryMax 20G, PYTHONPATH ModelSEEDpy-c2e31a0: `weighted_objective_test.py --population-only`, then `--workers 8 --maxtasksperchild 2`), then `scripts/wobj_results.py` | `cells.csv.gz`, `fold_plan`, `panel`, `_cache_panel_ec.json`, `panel_ec_universe.json`, `prfba/*/predictions.json` | 20:40→22:25 (**1.75 h**); parent 0.74 GB, worker 1.80 GB | **re-run** (the launcher refuses to run if `data/weighted_objective` exists: move it aside) |
| Stage E estimators (`f2ec19a`) | `$PY scripts/pfunc_estimators.py --workers 12 --bootstrap 2000 --seed 1729 --densities 1.0,0.5,0.25 --outdir data/pfunc_estimators --cache-dir <scratch>` | `cells.csv.gz`, `prfba/f*_dDDD/{hits,predictions.json,_cache_ec_subset.json}`, `_cache_panel_ec.json` | **411 s**; parent RSS 7.94 GB, largest child 9.87 GB (`run.log`) | **re-run** |
| Stage O objective tuning (`b52984f`, shipped `6e2c2eb`) | `OTUNE_WEIGHTS=…/stage_o_weights_E3.json bash scripts/launch_objective_tuning.sh` (`otune-run`, MemoryMax 34G, 10 workers, `--maxtasksperchild 1`), then `$PY scripts/objective_tuning_results.py`; FVA variant to `data/objective_tuning_fva` (193 rows) | `weighted_objective/population.json`, `stage_o_weights*.json`, per-fold `_cache_ec_subset.json` | started 11-03 08:50; 131.5 s mean and 291 s max per amplicon; worker 2.14 GiB; 360 rows (`OBJECTIVE_TUNING_RESULTS.md`). End time unmeasured; the estimate is 1.3–1.9 h | **re-run** (the launcher refuses if `per_amplicon.jsonl` exists) |
| EGC / template-direction audits (`egc_*.py`, `template_direction_audit.py`, `template_validation_summary.py`, `alignment_growth_probe.py`, `validate_builder_reversals.py`) | see `logs/egc_*`, `logs/tmplval_*` | the **pre-registered** phenotype rows and cells | — | **no**: diagnostic record of the pre-registered C12 failure, fixed by amendment |
| `design_audit_*.py` | standalone CLIs; `score_folds.py` imports their baseline and ICC code | `genome_function_sets.json`, `comparison_cells.csv.gz` | — | no committed output; nothing to re-run |

---

## 4. The two paths

### Path A: frozen design (keep `panel_eligible.json` eligibility, fold plans, refs' trees, placements' 16S, alignments)

As-shipped chain, in dependency order:
1. `prfba_selection.py` (default and `--legacy`). Local; reads hits, truth and the cache.
2. **Network:** `P2C/scripts/fetch_bvbrc_annotations.py pgfams`, once with the default
   selection and once with `--selection data/prfba_selection_legacy.json` (H8). It takes its
   genome set from the selections step 1 wrote.
3. `function_maps.py` (pin the templates and database: H9)
4. `genome_function_sets.py`
5. `annotation_ceiling.py`, `compare.py`, `render_comparison.py`, `annotation_dialect_test.py`,
   `depth_probe_sample.py` (diff against the committed keys), `depth_probe_eval.py`
6. MAG: delete `mag_arm/panel_ec.json`; run `mag_direct_feed.py` and `mag_score.py`

Shared reference:

7. **Genus campaign** (`SR/data`):
   - Regenerate `data/_cache_panel_ec.json`, and only that file: `panel_ec_universe.json` and
     `panel_eligible.json` can stay byte-identical (3B stage 1). No CLI exists for this
     one-off.
   - Refresh the 23 affected refs (traits only).
   - Re-run the 28 PICRUSt2 runs.
   - Re-select the 20 prFBA runs.
   - Score pass 1 (`--force`).
8. **C2** (`SR/data_c2`): refs traits, PICRUSt2, prFBA re-select, `score_folds --mode include_self`.
9. **Phenotype:** ATP and gapfilled passes (A-C13 configuration).
10. **Score pass 2:** `--require-phenotype-binding --phenotype-summary data/phenotype_gf_c13/phenotype_summary.json`,
    templates pinned at `2294299`, `PYTHONPATH` ModelSEEDpy-d884a50 (as in `c13_finish.sh`).
11. **Family** (`SR/data_family`): regenerate `_cache_panel_ec.json`; then refs traits,
    PICRUSt2, prFBA re-select, `score_family.sh`.
12. **Post-hoc studies:** `calibrate_pfunc`, `pfunc_estimators`, `launch_wobj.sh`, then
    `launch_objective_tuning.sh` (+ results, FVA).

Wall-clock estimate, from measured parts:

| Block | Estimate |
|---|---|
| As-shipped (steps 1–6) | about 0.5–1 h, mostly cache loads (about 40 GB each) |
| Refs | unmeasured; about 0.5 h at 12 concurrent |
| PICRUSt2 | about 20 min |
| prFBA re-select | about 30 min |
| Score pass 1 | 1.6 h |
| C2 | about 1.5 h |
| Phenotype (both passes concurrently) | about 3 h |
| Score pass 2 | 1.4 h |
| Family | about 1 h |
| Post-hoc studies | about 4–5 h (wobj, then otune) |
| **Total** | **about 14–16 h, mostly sequential** |

What Path A leaves untrue:
- `panel_eligible.json` keeps 123 `empty_ec` drop reasons that are no longer true.
- The panel under-represents those 65 Archaea relative to a fresh run of the same code.
- The prFBA alignments come from the 10-14 aligner code.

### Path B: re-plan (`plan_folds.py --phase panel` from regenerated data)

Steps 1–6 as in Path A, then the whole of 3B–3F, from stage 1 and with `--force` or a fresh
root:
- prepass about 5.3 h (old aligner code; the current one may differ);
- 30 builds about 31 min wall;
- 35 PICRUSt2 runs about 22 min;
- 15 base prFBA alignments, 3,100–4,400 s each at 2 concurrent (about 7–8 h), plus 5 forced
  runs at about 3 min each;
- score 1.6 h;
- C2 (5 × 3,400 s alignments plus 51 min scoring);
- phenotype about 3 h;
- score 1.4 h;
- family (plan, 30 builds, 13,201 s summed; 20 prFBA runs, 10,478 s summed; score);
- post-hoc 4–5 h.

**About 1.5–2 days.**

It is an amendment, not a rescore. SPEC §1 freezes only the `frozen` block, but `n_eligible`,
the fold composition and the median reference size of 7,175 are stamped into every artefact
and into the manuscript sentence. Expect the median reference to grow by about 90 genomes and
archaea by about 5%.

---

## 5. Every stale copy of BV-BRC-derived data, and how to regenerate it (in this order)

| # | Stale copy | Built from | How to regenerate | Reused silently if left? |
|---|---|---|---|---|
| 1 | `P2C/data/prfba_selection{,_legacy}.json` | PGFam cache (gene counts) | `prfba_selection.py [--legacy]` (first: #2 takes its genome set from these) | — |
| 2 | `P2C/data/pgfam_products.json` | live API, 09-17 | `fetch_bvbrc_annotations.py pgfams`, once per selection file (network; incremental) | yes: missing ids are projected as no EC |
| 3 | `P2C/data/function_maps.json` | #2 + ModelSEED templates and database | `function_maps.py` | — |
| 4 | `P2C/data/genome_function_sets.json` | PGFam cache + #1 + #3 | `genome_function_sets.py` | — |
| 5 | `region_validation_d100/picrust2_comparison/data/{pgfam_products,function_maps,genome_function_sets,prfba_selection*}.json` | same chain, 10-20 | same scripts in that tree (out of scope) | — |
| 6 | `SR/data/_cache_panel_ec.json`, `SR/data_family/_cache_panel_ec.json` | #4 | Path A: one-off call of `plan_folds._cache_panel_ec(..., force=True)`. Path B: `plan_folds.py --phase panel --force` | **yes**: `plan_folds.fresh()` and `build_reference.load_ec_sets()` reuse it whenever it exists |
| 7 | `P2C/data/kfold_audit/panel_ec_sets.json` (10-13) | #4 | regenerate, or delete it | **yes**: `build_reference.load_ec_sets()` falls back to it when #6 is missing (covers all panel genomes, so it passes the coverage check) |
| 8 | `SR/data/panel_ec_universe.json` (= `data_family/…`, hard link) | #6 | Path B: `--phase panel`. Path A: leave it. The EC list is unchanged as measured, and the `prevalence` block is read by nothing | — |
| 9 | `SR/data/panel_eligible.json` (= `data_family/…`) | #6 | Path B: `--phase panel`. Path A: leave it. `n_ec` is read by nothing; the 123 `empty_ec` reasons stay stale | — |
| 10 | `SR/data{,_family}/refs/f*/traits_EC.txt.gz` + `manifest.json`; `SR/data_c2/refs/incself_*/…` | #6 | `build_reference.py` after deleting `manifest.json` | **yes**: a dir with `manifest.json` is skipped |
| 11 | `SR/data*/prfba/*/_cache_ec_subset.json` (45) | #4 | delete; `run_prfba_fold.py` re-streams it | **yes**: `restricted_map()` reuses it if it covers the keys (unless `--force`) |
| 12 | `SR/data*/prfba/*/_cache_pgfam_subset.json` (45) | PGFam cache | delete, as #11 | **yes**, as #11 |
| 13 | `SR/data*/picrust2/*/predictions.json` | #10 | re-run (stamps detect the traits change) | the driver's `p2_done` checks meta only: **yes** via `run_campaign.py` without `--force`; the runner itself re-runs on a stamp mismatch |
| 14 | `SR/data*/prfba/*/predictions.json` | #11, #12 | delete, then re-run | **yes**: the runner returns `[skip]` when the stamp and meta match |
| 15 | `SR/data/cells.csv.gz` (+ `.fp`), `data_c2/cells_incself.csv.gz`, `data_family/cells.csv.gz` | everything above | `score_folds.py --force` | no: the `.fp` sidecar checks the answer-key bytes |
| 16 | `P2C/data/mag_arm/panel_ec.json` | #4 | delete; `mag_common.load_panel_ec()` rebuilds it (orjson load of #4) | **yes** |
| 17 | `mag_validation/pgfam_panel.json` (10-08) | PGFam cache | MAG-arm extraction step (out of scope) | yes, for `mag_prfba_ec.py` |
| 18 | `P2C/data/comparison_cells.csv.gz` | #3, #4 | `compare.py` | read by `calibrate_pfunc.py` and `design_audit_*` |

Hard links to remember: `data_family/{panel_eligible,panel_ec_universe,neighbors,prepass_*}.json`
point at `data/`'s inodes, and `prepass_amplicon.json` and `hits_amplicon/asv_top20_alignment_hits.json`
share an inode four ways.

---

## 6. Hazards

- **H1. Caches are reused when they merely exist.** Stale copies #6, #7, #10, #11, #12, #14 and
  #16 are reused without any provenance check. Delete them in the order of section 5. Missing
  #6 while #7 exists gives stale traits silently.
- **H2. `--force` means re-align.**
  - `run_prfba_fold.py --force` re-runs the alignment, about 4,000 s per base run.
  - `build_reference.py --force` re-runs hmmalign, FastTree and the probe placement.
  - `run_campaign.py --force` passes `--force` to every job.
  - For Path A, delete the outputs named in 3B and run without `--force`.
  - The driver's completion predicates (`p2_done`, `pf_done`, `manifest.json` exists) read only
    meta. `stale_dependents()` warns about what a `--force` leaves stale.
- **H3. Shared cache writer.** `bvbrc_cache/genome_gene_families.json` is shared with production
  code, which can write it (the 10-09 postcheck tested "provider wrote no poison"), and with
  other sessions. The pending KBase upload of the new builds is one such session.
  - Writes are atomic (tmp + `os.replace`), so a reader sees one consistent version. Two steps
    run hours apart can still see different versions.
  - Record the cache sha256 at the start (11-13 value: `05f03bcd…`) and check it at the end.
  - Only `run_prfba_fold.py` takes `--pgfam-cache`. `prfba_selection.py`, `genome_function_sets.py`,
    `fetch_bvbrc_annotations.py` and `depth_probe_eval.py` use the `_common.PATHS.pgfam_cache`
    constant, so they cannot be pointed at a snapshot without a code edit.
- **H4. Memory.**
  - `prfba_selection.py` ran at 25.8 GB RSS on the smaller 09-17 cache. Expect about 40 GB now
    for each orjson load of the 5.12 GB cache (`prfba_selection`, `fetch … pgfams`,
    `genome_function_sets`, `depth_probe_eval`). This is not measured in any log.
  - `score_folds.load_functions` and `phenotype_fold` `jload` all of the 1.7 GB file; their RSS
    is not measured in any log. The comment in `run_prfba_fold.restricted_map` puts a
    whole-file parse of both large files at 30–40 GB.
  - The phenotype gapfilled pass was oom-killed at 125 GB (unpinned, A-C12). The pinned runs
    used 10 and 16 workers with `--maxtasksperchild`.
  - The prepass reads a 6.7 GB `prepass_amplicon.json` by streaming.
  - Launch long steps as `systemd-run --user` units (memory note); wobj and otune already do.
- **H5. Resumable or not.**
  - `plan_folds` phases: fingerprint- or existence-gated.
  - `build_reference`: per-dir manifest; FastTree is skipped when the tree, log and stderr exist.
  - `run_picrust2_fold`: stamp-gated per stage.
  - `run_prfba_fold`: align stamp plus meta.
  - `run_campaign`: ledger plus predicates. Stages 7+ are not in the driver.
  - `phenotype_fold`: skips finished genomes in `per_genome.jsonl` unless `--force`. A changed
    cells table should get a **new outdir**, because otherwise rows built for the old genome set
    are kept.
  - `launch_wobj.sh` and `launch_objective_tuning.sh` refuse to overwrite existing outputs.
  - `score_folds`: the `.fp` sidecar.
- **H6. Hard links between `data/` and `data_family/`.** Re-running `plan_folds` in `data/`
  breaks the links (tmp + replace) and leaves `data_family/` on the old inodes, with nothing to
  say so. Re-link explicitly, or rebuild data_family from `--stages 1`.
- **H7. Determinism and byte identity.**
  - Every JSON meta embeds `generated` and `git_commit`, and `cells.csv.gz.fp` embeds
    timestamps, so **no committed file reproduces byte for byte**.
  - What is reproduced: content-level numbers. The A-C12 and A-C13 rescorings reused
    `cells.csv.gz` by fingerprint and reproduced arms, contrasts and nested CV exactly (seeded
    bootstrap). The A-C12 ATP pass reproduced all 43,140 pre-registered rows. Phenotype rows
    are independent of worker history (84/84, 42/42). The aligner commits report byte-identical
    `asv_top20_alignment_hits.json` from saved shortlists.
  - `prfba_selection` and `select_references` are deterministic.
  - FastTree and epa-ng determinism are **not verified in any log**. To check, compare
    `manifest.tip_labels_sha256`, `model_string` and the tree bytes.
  - `resolve_nonroot_polytomies` is idempotent.
  - `phenotype_fold` and `pfunc_estimators` are seeded (1729).
- **H8. The legacy-reps fetch gap** (section 0, point 5). Fix it by fetching PGFams for both
  selection files.
- **H9. Moving ModelSEED inputs.**
  - `function_maps.py` hard-codes `MST=/home/freiburger/Documents/ModelSEEDTemplates/templates/v7.0`
    (moved, at `e7fc410`) and `MSD=/home/freiburger/Documents/ModelSEEDDatabase` (moved). It has
    no CLI flag for either.
  - `score_folds.py --templates` defaults to the same moving tree. The c12 and c13 rescorings
    pinned `ModelSEEDTemplates-main2294299` and ModelSEEDpy `d884a50`. Keep both pins, or
    declare the change: it moves `answer_key_fp`.
  - The phenotype passes need `ModelSEEDTemplates-c12run` (`21a8084`) and ModelSEEDpy `c2e31a0`
    pinned through PYTHONPATH and `--expect-modelseedpy`.
- **H10. Aligner code changed after `b72b71cf`** (section 0, point 7). Reusing hits and
  re-aligning are different claims; write down which one the rescore uses.
- **H11. `controls.json` is a merge across modes** (recipe_control C1, include_self C2, fold).
  - A stale C1 shows up as `from_other_modes` with `answer_key_matches=false`. That is a printed
    warning, not a failure.
  - Re-run in this order: C2, then fold pass 1, phenotype, fold pass 2. Optionally re-score C1
    first.
- **H12. Path A has no CLI** for refreshing the answer-key caches while keeping eligibility
  (3B stage 1). The one-off should regenerate only `_cache_panel_ec.json` (in `data/` and
  `data_family/`). Do not run `plan_folds.py --phase panel`: it would re-decide eligibility,
  and through the new md5 set change `panel_fingerprint`, after which prepass and plan rebuild,
  which is Path B. The phenotype binding moves anyway, through `cells.csv.gz` and `--functions`.

---

## 7. Validation checks that exist

| Check | Where |
|---|---|
| `plan_folds`: leakage preflight, all counters 0 (exit 2); k2-truncation guard (re-run with doubled k2, ×3, then abort); `neighbors.json` panel-fingerprint check; md5 collision; header round trip | `plan_folds.py` |
| `build_reference`: §5.1 geometry abort (gap > 0.10, width > 2,000); C4 set equality (`manifest.set_equality_ok`); crash screen (one retry, then abort) | `build_reference.py` |
| PICRUSt2 runner: C3 (preflight and post-hoc; exit 3), C4, C7 accounting | `run_picrust2_fold.py` |
| prFBA runner: C6a–d, C8 canaries (exit 3); forced-answer invariant (0 mismatches on answered cells); `reps_outside_reference == 0` | `run_prfba_fold.py` |
| Driver: `campaign_inventory.json`, i.e. expected vs discovered runs plus a ledger-vs-disk re-fingerprint (exit 1) | `run_campaign.py` |
| Scorer: C1–C15; exit 0/1/2/3; `.fp` reuse contract; `answer_key_fp`, `scorer_code_fp`; `inputs_changed_during_scoring`; D22 denominator tripwire | `score_folds.py` |
| Phenotype: `--expect-modelseedpy`; `meta.binding` (contract 1); C11 identity control; C12 gate (0 gated truth, ≤ 10% gated) | `phenotype_fold.py`, `score_folds.py --require-phenotype-binding` |
| Self-tests: `_selftest_run_picrust2_fold.py`, `_selftest_score_folds_inventory.py`, `_selftest_score_folds_provenance.py`, `_selftest_fixture_build_reference.py`, `test_ref_db_override.py` | `SR/scripts/` |
| Post-hoc: `pfunc_estimators/assertions.json` ("all required zero/true"); `objective_tuning.py --validate-only`; `prior_match` / `released_check` counters; `c12_compare.py` before/after diffs | per study |
| As-shipped: `prfba_selection` meta `n_reps_missing_from_pgfam_cache` (0 at 09-17); `genome_function_sets` `_meta.n_with_empty_ec` (876 at 09-17) | `P2C` |
