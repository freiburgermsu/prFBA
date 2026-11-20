# Rescore plan: MAG validation, real MAGs, depth-100 deployments, E3 gate, downstream_e3

Written 2026-11-15 by a read-only investigation (prFBA HEAD `cc8014f`, branch
`hardware-aware-alignment-dispatch`; working tree: `deployment_d100/README.md` and
`upload_kbase_genomes.py` modified, `deployment_d100/kbase_upload_bvbrc_d100.sh` untracked).
Nothing in the repo, the caches or the deployments was written. The two measurement scripts
and their outputs used below are in
`/tmp/claude-1000/-home-freiburger-Documents-prFBA/400ec395-f172-59da-a98d-3dbe6b3a4256/scratchpad/rescore/work/`
(`panel_len_check.{py,json}`, `tieset_check.{py,json}`). Paths are absolute unless they start
with a repo-relative directory name, which means `/home/freiburger/Documents/prFBA/<path>`.

User requests that frame this: the full rescore ("all of the values in the manuscript and in
the repo [are] the result of the most recent code", 2026-11-15) and "yes also upload all of the
new builds to KBase". The upload is covered under E.6 below as a step with ordering
constraints; this report does not run it.

---------------------------------------------------------------------------------------------

## 0. Conclusions

1. **The MAG panel arms are stale, but only in gene capture, and only through the three
   int-encoded panel indices.** `score_mag.py` never reads the shared cache. It reads
   `mag_validation/pgfam_panel.json`, `kmeans_panel/pgfam_panel_kmeans.json` and
   `random_panel/pgfam_panel_random.json` (git-ignored; built 2026-10-08 13:18, 10-09 14:00,
   10-18 17:26). Streaming the current shared cache (measured, 23 s, 1.0 GB) shows every panel
   genome is present (30,000/30,000), and the per-genome list lengths differ in exactly
   **140 / 20 / 38** genomes. Every one of them is a 2026-11-13 repaired entry
   (original: 123 empty_list_poison + 16 truncated_by_refseq + 1 truncated_by_patric_cap;
   k-means 9 + 8 + 3; random 29 + 5 + 4). No other drift. Empty lists fall 161 -> 38,
   45 -> 36 and 72 -> 43.
2. **Per-cell selection is invariant to the repair in the MAG arm.** Anchor, picks, tier,
   n_selected, ani_class and the six per-rank taxonomy cells cannot change, and neither can
   the per-rank taxonomy summaries (their scoreable set comes from the truth lineage).
   `score_mag.py` passes `gene_provider` but not `gene_count_fn`; `ANI_KNOBS` does not
   override `select_all`, so the union block of `select_representatives`
   (`select_references.py` lines 327-483) runs, and it reads only `gene_count_fn`. The
   `gene_provider` reads are at lines 560-610, in the legacy reducer only.
   **The aggregate figures do change.** A cell is labelled `truth_empty` whenever its
   source's PGFam set is empty, whatever was picked. `aggregate.py` computes
   `coverage_answered = n_ok / (n_ok + n_abstained)`; it drops truth_empty cells and counts
   pred_empty as unanswered. After the repair:
   * the 1,220 stale-source cells per regime enter that denominator as ok or abstained;
   * stale-pick pred_empty cells become ok.
   So these move: the "resolve rate 0.306-0.308", the coverage column of `tab:magvs16s`,
   `n_ok` / `n_abstained` / `n_truth_empty`, conditional and population gene capture, and
   the paired n (1,916 ...). The paired n moves again when the amplicon side's 675
   stale-source cells per mode become scorable. Only the GTDB arm (no truth lineage
   lookup, no PGFams) is invariant in every number.
3. **The index builder is not in the repo**, and the only project transcript that mentions
   `n_vocab` (99b1c5b4, 2026-10-16) shows reads of it only. It has to be
   rewritten (spec in A.4). The int encoding is arbitrary for `score_mag.py` (it only
   intersects and unions sets within one index), so the index file cannot and need not be
   byte-identical. The cell tables can be, and that is the control (A.5).
4. **Every score driver skips a batch whose `mag_cells.csv` exists**, so re-running
   `score_all.sh` and the others in place is a silent no-op. Move the `out_batches*`
   directories aside first, or write to new directories.
5. **Hard ordering edge: the amplicon gene-capture tables of `region_validation_d100`
   (other report's scope) must be re-scored before (a) the E3 gate and (b) every paired
   MAG-vs-16S output.** The gate's built-in reproduction check compares each amplicon's
   (outcome, n_reps, truth_size, pred_size, intersection) with the committed
   `gene_capture_cells.csv` and `gene_capture_cells_tau0995.csv`. On the repaired cache
   1,389 / 2,343 / 2,171 rows differ (include / 0.987 / 0.995, impact study), so a gate run
   before that re-score is marked "INVALID RUN" and exits 3.
6. **Reusable unchanged**: every contig file, both donor sets, the three panel sketches,
   every hits TSV and synth manifest (main grid, tau99 (same hits), n500, lowc, GTDB x3,
   k-means, random), the whole GTDB arm's scored outputs, the pilot, the 291 real-MAG
   outputs, the deployments' depth-100 hits, and the two BV-BRC-only deployment builds
   (already built from the repaired caches).
7. **The deployments are already on repaired data.** The open questions are code-stamp
   purity for codiffusion (its build predates `8cfc297`; equivalent by scan, not by
   rebuild), whether to refresh the tie-set caches, and the KBase upload, which must come
   after any rebuild decision.
8. **New finding, not in any audit: the contig fetches are capped at 25,000 rows too.**
   `fetch_contigs*.py` sends `limit(50000)`, but two random-panel genomes (1733.9197,
   1733.9535, batch b004) and two `~/Documents/sludge_genomes` genomes (1733.9522,
   1733.9680) hold exactly 25,000 contigs. All four genome ids share the `1733.` prefix. The two
   random-panel ones are *M. tuberculosis* (taxon 1773 in `random_panel/selection_random.json`)
   and hold 15.7 and 14.1 Mb of contig sequence. The sludge files were counted for headers
   only. Probably truncated; the impact is tiny
   (A.6/B, D).

---------------------------------------------------------------------------------------------

## A. MAG panel arm, original panel (`mag_validation/`)

### A.1 What produces each committed file

Earliest BV-BRC step: the contig fetch (sequence, not gene data). The first gene-data step
is the panel-index extraction.

| # | step (command as run) | output | reads BV-BRC data? |
|---|---|---|---|
| 0a | `python mag_validation/fetch_contigs.py` (NWORK=16; reads `region_validation/data/selection_10k.json`) | `mag_validation/contigs/<gid>.fna` (9,997; 3 empty) | genome_sequence API, network |
| 0b | `bin/skani sketch -l mag_validation/ref_files.txt -o mag_validation/skani_db_panel -t 61` (form as in `run_mag_rnd.sh`) | `skani_db_panel/` (3.7 GB) | contigs |
| 0c | donor selection (contaminated tiers only; withdrawn) | `donor_candidates.json`, `donors.json`, `donor_assignment.json`, `donor_contigs/` | contigs |
| **1** | **panel-index extraction (script lost; see A.4)** | **`mag_validation/pgfam_panel.json`** (`{"n_vocab":1862396,"genomes":{gid:[int]}}`) | **shared PGFam cache** |
| 2 | `bash mag_validation/run_full.sh` (or `run_batch.sh` x 8 concurrent, `panel_par.log`): `synth_mags.py` -> `skani search -n 100 --min-af 0` -> delete MAGs | `hits/bNNN.tsv`, `manifests/bNNN.json` (20 batches, 10 tiers, 99,970 MAGs) | contigs, sketch |
| 3 | `bash mag_validation/score_all.sh` (`score_mag.py --hits --manifest --outdir out_batches/bNNN`; defaults `--pgfams mag_validation/pgfam_panel.json`, `--truth region_validation/data/benchmark_truth.json`; xargs -P 6) | `out_batches/bNNN/{mag_cells.csv,mag_summary.json}` (git-ignored) | **panel index** |
| 4 | concatenation, no script: header of b000 + `tail -n+2` of each `out_batches/b*/mag_cells.csv` in name order (verified: md5 `a57df30b8346fcb72f14f7f33708669a` both ways) | `out/mag_cells_all.csv` | - |
| 5 | `python mag_validation/aggregate.py` (defaults) | `out/mag_summary_full.json` | - |
| 6 | `python mag_validation/compare_arms.py --mag-cells out/mag_cells_all.csv --out out/paired_comparison.json` (amp default = depth-20 `region_validation/data/gene_capture_cells.csv`) | `out/paired_comparison.json` | amplicon cells |
| 7 | `region_validation_d100/run_downstream_d100.sh` MAG block, for S in `_d100`, `_d100_legacy`: `compare_arms.py --amp region_validation$S/data/gene_capture_cells.csv --mag-cells out/mag_cells_all.csv --out out/paired_comparison$S.json`; `crossover.py --amp region_validation$S/data/gene_capture_cells.csv --out out/crossover$S.json out/mag_cells_all.csv lowc_out_batches/b*/mag_cells.csv` | `out/paired_comparison_d100{,_legacy}.json`, `out/crossover_d100{,_legacy}.json` | amplicon cells |
| 8 | `python mag_validation/crossover.py out/mag_cells_all.csv lowc_out_batches/b*/mag_cells.csv` (depth-20 amp default) | `out/crossover.json` | amplicon cells |
| 9 | tau99: `bash mag_validation/score_tau.sh 99` -> `out_batches_tau99/` (tracked); concat (md5 `3cb0e487aa31c2f30f1699fcbdf88c92` verified) -> `out/mag_cells_tau99.csv`; `aggregate.py --batches out_batches_tau99 --out out/mag_summary_tau99.json`; `compare_arms.py --mag-cells out/mag_cells_tau99.csv --out out/paired_tau99.json` | as named | panel index |
| 10 | n500: `bash mag_validation/run_mag_n500.sh` (5 tiers, -n 500) -> `hits_n500/`, `manifests_n500/`; `bash mag_validation/score_n500.sh` (-P 5) -> `out_batches_n500/` (tracked); `aggregate.py --batches out_batches_n500 --out out/summary_n500.json` | `out/summary_n500.json`, `out_batches_n500/b*/` | panel index |
| 11 | lowc sweep (9 tiers 25%..1%): synth+search driver not committed (`run_batch.sh` with a GRID env and `lowc_*` dirs; manifests carry the post-fix `completeness_vs_pool` field) -> `lowc_hits/`, `lowc_manifests/`; `bash mag_validation/score_lowc.sh` -> `lowc_out_batches/` (git-ignored, but read by step 7/8) | `lowc_out_batches/b*/mag_cells.csv` | panel index |
| 12 | `python mag_validation/af_gate_counterfactual.py` (no CLI; `SM.load_pgfams()` default path; reads `manifests/b000..b019`, `hits/b000..b019`; mutates `ANI_KNOBS['cov_min']` in-process; hard-coded 16S levels 0.6912 / 0.71 / 0.7092) | `out/af_gate_counterfactual.json` | panel index |
| 13 | smoke: `score_mag.py --hits smoke/hits.tsv --manifest smoke/manifest.json --outdir smoke/out` | `smoke/out/{mag_cells.csv,mag_summary.json}` | panel index |
| 14 | figures: `manuscript_figures/mag/plot_results_d100.py` (reads `out/mag_summary_full.json`, `out_gtdb/gtdb_summary_2000.json`, `out/paired_comparison_d100.json`) -> F14; `manuscript_figures/mag/plot_crossover_d100.py` (reads `out/crossover_d100.json`) -> F15; copy to `~/Documents/prFBA_paper/images/prfba/F14_mag_completeness.pdf`, `F15_mag_crossover.pdf`. The old `mag_validation/plot_results.py` / `plot_crossover.py` -> `out/F_mag_*.{pdf,png}` read the depth-20 files | figures | - |
| 15 | cross-panel: `two_panel/compare_mag_panels.py` (reads `out_batches_n500` + `mag_validation_km/out_batches`) -> `two_panel/mag_panel_comparison.json`; `three_panel/compare_three_panels_magpanel.py` (reads `out_batches_n500`, km and rnd `out_batches`) -> `three_panel/three_panel_mag_panelarm.json` | as named | - |

Text artefacts (`out/RESULTS_v2.md`, `RESULTS_tau99.md`, `CORRECTIONS.md`, `DEPTH_PROBE.md`)
quote numbers from these files by hand.

### A.2 Re-run or reuse

* **Reuse** steps 0a-0c and 2 (hits and manifests do not depend on gene data), the
  lowc and n500 hits, and the tau99 batches' inputs (tau99 re-scores the same `hits/`).
  Do **not** re-run synthesis or search: `synth_mags.py` was fixed after the main run
  (contamination dose and `completeness_actual`, `out/CORRECTIONS.md`). A re-synthesis would
  not reproduce the committed contaminated tiers or the manifests' completeness fields, and
  it buys nothing for the rescore.
* **Re-run**: step 1 (new index), then 3-13 and 15, then 14 and the hand-written text.
  Steps 6-8 and 12 also need the re-scored amplicon cells (section H).
* Expected change (impact study's exposure counts for `out/mag_cells_all.csv`, per regime):
  1,220 of the 1,600 truth_empty cells have a stale source and become scorable; 287
  exclude-regime ok cells have a stale anchor (377 at tau99); 170 cells have a truncated
  source; 107 exclude ok cells a truncated anchor. Non-anchor picks are not recorded in the
  cell table, so a few more cells can change through a repaired second pick. In the
  summaries, `coverage_answered` (the resolve rate), `n_ok`, `n_abstained`,
  `n_truth_empty`, `n_pred_empty` and every gene-capture block move. The per-rank
  `n_scoreable` / `n_called` / `correct` / `correct_population` blocks do not move.

### A.3 Runtime, memory, network

| step | measured | memory | network |
|---|---|---|---|
| 0a contig fetch | 488 s, 16 workers, 20.5 genomes/s (`logs/fetch.log`) | small | yes (not needed) |
| 2 synth+search, -n 100 | 507 s for 20 batches, 8 concurrent x 6 threads; per batch synth 23-29 s, search 92-178 s (`logs/panel_par.log`) | not logged | no |
| n500 synth+search | 176 s (`logs/run_n500.log`) | not logged | no |
| lowc synth+search | 65 s; 12 batches lost 1-3 queries each (`!!MISSING` in `logs/lowc.log`) | not logged | no |
| **1 index extraction (new)** | one ijson pass over the 5.12 GB cache: 23-26 s (measured twice here, yajl2_c) | ~1 GB for all three panels | no: all 30,000 genomes are cached |
| **3 score, main grid** | 153 s for 20 batches at -P 6 (span of `logs/score_b0*.log` mtimes) | not logged; measured here: index load 1.54 GB (2.7 s), taxopy `TaxDb` 2.7 GB (13 s), so about 4.5-5.5 GB per process with truth and hits, about 30 GB at -P 6 | no |
| tau99 score | 549 s | same | no |
| n500 score | 104 s (-P 5) | same | no |
| lowc score | not logged (same order) | same | no |
| 5 aggregate, 6-8 compare/crossover, 12 af-gate | seconds to a minute each (not logged) | < 2 GB (af-gate loads the index: ~2 GB) | no |

### A.4 Stale copies and how to regenerate them

| file | stale how | regenerate |
|---|---|---|
| `mag_validation/pgfam_panel.json` (106 MB, git-ignored, built 2026-10-08) | 140 repaired genomes hold their old lists (`len_eq_old` 140, `len_eq_new` 0) | rebuild from the shared cache (spec below) |
| `mag_validation/out_batches/`, `out_batches_tau99/`, `out_batches_n500/`, `lowc_out_batches/`, `smoke/out/` | scored on the old index | re-score |
| `out/mag_cells_all.csv`, `mag_cells_tau99.csv`, `mag_summary_full.json`, `mag_summary_tau99.json`, `summary_n500.json`, `paired_*.json`, `crossover*.json`, `af_gate_counterfactual.json` | derived | re-derive |

There is no committed or surviving builder of the index. Regeneration spec (proposed script,
not in the repo; the same for all three panels):

```python
# build_pgfam_panel_index.py --old-index OLD --out NEW [--restore-old-values]
# Read-only on the shared cache. Genome set and key order come from OLD (the original panel's
# order matches neither selection_10k.json nor benchmark_truth.json, so it cannot be re-derived;
# k-means and random follow their selection files). Ints are assigned by first appearance over
# genomes in that order, PGFam strings sorted within a genome; lists are sorted unique ints.
import json, os, hashlib, ijson
CACHE = "/home/freiburger/Documents/prFBA/bvbrc_cache/genome_gene_families.json"
CORR  = "/home/freiburger/Documents/prFBA/bvbrc_cache/CORRECTIONS_2026-11-13.json"
keys = list(json.load(open(OLD))["genomes"]); want = set(keys)
st0 = os.stat(CACHE)
vals = {}
with open(CACHE, "rb") as fh:
    for g, v in ijson.kvitems(fh, ""):
        if g in want:
            vals[g] = sorted(set(v))
if RESTORE_OLD:                                   # the control index: pre-repair values
    bk = json.load(open(CORR))["backup_old_values"][CACHE]
    for g in want & set(bk):
        vals[g] = sorted(set(bk[g]))
assert want <= set(vals), "panel genome missing from the cache"   # 0 missing today
vocab, genomes = {}, {}
for g in keys:
    genomes[g] = sorted(vocab.setdefault(f, len(vocab)) for f in vals[g])
tmp = NEW + ".tmp"
json.dump({"n_vocab": len(vocab), "genomes": genomes}, open(tmp, "w")); os.replace(tmp, NEW)
# + a sidecar NEW + ".stamp.json": cache size/mtime/sha256 before and after (must be equal),
#   corrections-log sha256, n_genomes, n_empty, list of genomes whose length differs from OLD.
```

Checks on the new index: same key set and order as the old; `len(new[g]) == len(set(cache[g]))`
for every g; lengths equal to the old index for every genome outside the corrections log;
old length == log `old_size` and new length == log `new_size` for the 140 / 20 / 38;
`|A & B|` equal between old and new indices for a random sample of genome pairs that avoid
repaired genomes (proves the encoding is consistent).

Where to put it: keep `score_mag.py`'s default path working by renaming the old file to
`mag_validation/pgfam_panel.pre_repair_2026-11-13.json` and writing the new one at
`mag_validation/pgfam_panel.json` (the same for `kmeans_panel/` and `random_panel/`).
`af_gate_counterfactual.py` has no `--pgfams` flag and reads only the default path.

### A.5 Hazards and validation

* **Silent no-op**: every driver (`score_all.sh`, `score_tau.sh`, `score_n500.sh`,
  `score_lowc.sh`, `score_km.sh`, `score_rnd.sh`) prints "cached" and skips any batch whose
  `mag_cells.csv` exists. Move each `out_batches*` aside (or call `score_mag.py` directly
  into new directories and swap). `mag_validation/out_batches/` and `lowc_out_batches/` are
  git-ignored, so their old contents exist only on disk: keep a copy until the new run is
  verified.
* **Byte-identity control (do first)**: build a control index with `--restore-old-values`
  and score the main grid with it into scratch. The concatenated cells must reproduce
  `out/mag_cells_all.csv` byte for byte (md5 `a57df30b...`), and the tau99 cells
  `out/mag_cells_tau99.csv` (`3cb0e487...`). This is expected (`score_mag.py` iterates
  `sorted(meta.items())`, hit sorting is stable, set operations are invariant to the int
  relabelling), but not yet demonstrated. It also proves that the current
  `select_references.py` (changed since the scoring only in `bvbrc_gene_provider`'s fetch
  path, `7153311`) still scores identically.
* **Expected-diff check on the real run**: the columns `mag, genome_id, completeness,
  contamination, regime, ani_class, tier, best_ani, n_hits, n_selected, anchor` and the six
  rank cells must be identical to the committed tables in every row. Rows whose
  `outcome/recall/precision/f1/jaccard` change must have the source genome or the anchor
  among the 140 repaired genomes, or `n_selected >= 2` (a repaired non-anchor pick).
  Row counts: 199,940 (main, tau99), 99,970 (n500), 20 lowc batches.
* **Memory**: about 5 GB per scorer process, about 30 GB at the drivers' -P 6, all
  unmeasured as a whole. Run under `systemd-run --user ... -p MemoryMax=48G` (memory note:
  nohup jobs die on container replacement) or lower -P.
* **Shared-cache writers**: the extraction only reads the cache. Stamp size, mtime and
  sha256 before and after, and do not run anything that writes the cache meanwhile
  (`select_references.py` CLI with its default `--gene-provider auto`,
  `prefetch_gene_families.py` with its default `--cache`, region-validation prefetchers).
  Such a write would also stale the E3 catalogue (section F).
* **Determinism of the consumers**: `aggregate.py`, `compare_arms.py` (scipy Wilcoxon) and
  `crossover.py` have no RNG. The three-panel bootstraps are seeded and hash-seed-proof since
  `5d2a204`. `two_panel/compare_mag_panels.py` builds its phylum/domain weights from sets
  (no RNG, but float summation order follows the hash seed): run it with
  `PYTHONHASHSEED=0` and compare at 4 dp.
* **Dependency on the other report's amplicon re-score**: steps 6-8 and 12 pair against
  `region_validation{,_d100,_d100_legacy}/data/gene_capture_cells.csv`. The
  `af_gate_counterfactual.json` thresholds (0.71 union full-length 16S level at depth 100,
  0.6912 legacy, 0.7092 depth 20) are hard-coded and must be re-read from the re-scored
  amplicon summaries. The impact study puts resolved recall within 0.05 pts, so 0.71 likely
  stands, but check it.
* **No automatic checker** covers the MAG outputs (`region_validation_d100/check_d100.py`
  has no MAG check). The checks above are the validation.

### A.6 Contig-level caveat (original panel): none

The largest original-panel contig file has 15,664 contigs (813.301), under the 25,000 cap.
Note: the manuscript (line ~815) says the panel spans "1 to 5,265 contigs", while 813.301's
file holds 15,664 (5,265 is the second largest, 114683.4). Worth checking against whatever
the text counted.

---------------------------------------------------------------------------------------------

## B. MAG k-means and random panel arms (`mag_validation_km/`, `mag_validation_rnd/`)

### B.1 Commands

| arm | step | command | output |
|---|---|---|---|
| k-means | contigs | `python mag_validation_km/fetch_contigs_km.py` (reads `mag_validation_km/panel_manifest.json`) | `contigs/` (9,996; 1,032 s, 12 workers) |
| k-means | **index** | lost builder; 2026-10-09 14:00, after the poisoning repair | **`kmeans_panel/pgfam_panel_kmeans.json`** (141 MB) |
| k-means | synth+search | `bash mag_validation_km/run_mag_km.sh` (5 tiers, -n 500, 8 x `run_batch.sh` with `SCRIPTS=mag_validation PANEL=mag_validation_km`) | `hits/`, `manifests/` (1,107 s; b017 search 955 s) |
| k-means | score | `bash mag_validation_km/score_km.sh` (`--truth region_validation_kmeans/data/benchmark_truth.json --pgfams kmeans_panel/pgfam_panel_kmeans.json`, -P 5) | `out_batches/` (tracked, 170 s) |
| k-means | summary | `aggregate.py --batches mag_validation_km/out_batches --truth region_validation_kmeans/data/benchmark_truth.json --out mag_validation_km/out/summary_km.json` (n_cells 99,960, 20 batches) | `out/summary_km.json` |
| random | contigs | `python mag_validation_rnd/fetch_contigs_rnd.py` then `fetch_contigs_rnd_full.py` (520 s for 8,011) | `contigs/` (9,998) |
| random | **index** | lost builder; 2026-10-18 17:26 | **`random_panel/pgfam_panel_random.json`** (153 MB) |
| random | sketch+synth+search | `bash mag_validation_rnd/run_mag_rnd.sh` (sketch 12 s if absent; 5 tiers, -n 500) | `skani_db_panel/`, `hits/`, `manifests/` (1,834 s) |
| random | score | `bash mag_validation_rnd/score_rnd.sh` (`--truth region_validation_random/data/benchmark_truth.json --pgfams random_panel/pgfam_panel_random.json`) | `out_batches/` (118 s; git-ignored: `mag_validation_rnd/*` is ignored except scripts and `out_gtdb`) |
| random | summary | none (`mag_validation_rnd/out/` is empty) | - |
| both | cross-panel | `three_panel/compare_three_panels_magpanel.py` -> `three_panel/three_panel_mag_panelarm.json`; `two_panel/compare_mag_panels.py` (k-means only) -> `two_panel/mag_panel_comparison.json` | tracked |

### B.2 Re-run or reuse

Reuse contigs, sketches, hits, manifests. Rebuild both indices (A.4 spec, `--old-index` =
the current file; genome order = the selection file's), move `out_batches/` aside, re-score
with `score_km.sh` / `score_rnd.sh`, re-run `aggregate.py` for k-means, then the two
cross-panel scripts. Exposure: 20 k-means and 38 random repaired genomes. The random arm's
committed record is only `three_panel_mag_panelarm.json` (its batches are git-ignored), so
keep the old `out_batches/` until the new JSON is checked.

### B.3 Stale copies

* the two indices (as A.4);
* **random-panel contigs of 1733.9197 and 1733.9535**: exactly 25,000 contigs each
  (15.7 Mb and 14.1 Mb of sequence), the server's row cap, although the request says
  `limit(50000)`. Both genomes' PGFam lists were also repaired as truncated_by_patric_cap
  (42,063 and 45,949 PATRIC CDS rows), and 1733.9535's genome-cache file too. They are queries
  of batch b004 (10 synthetic MAGs) and references in the random sketch. Neither is in the
  random GTDB subsample. The k-means maximum is 24,673 (1733.9588, complete).
  Regenerating them means a paged contig fetch (`limit(25000,start)` with a Content-Range
  check, as `prefetch_gene_families.fetch_patric_cds_rows` does), a re-sketch of the
  random panel and a re-search of all 20 random batches (they all search against the
  sketch). Measured cost: about 31 min of search plus the sketch. Recommendation: confirm
  the truncation with one count query per genome, and record it as a caveat rather than
  re-run (2 of 10,000 genomes, 10 of 49,990 queries).

### B.4 Hazards

As A.5. In addition: `score_mag.py` defaults to the ORIGINAL panel's truth and index, and
silently scores a k-means/random batch against the wrong organisms if `--truth/--pgfams` are
left off (the docstring warns). Use the `score_km.sh` / `score_rnd.sh` variables.

---------------------------------------------------------------------------------------------

## C. GTDB dense-reference arm (all three panels)

Commands: `bash mag_validation/run_gtdb.sh` (batches b000-b003, contamination grid) and the
contamination-free rest (b004-b007, `logs/run_gtdb_rest.log`); `score_gtdb.py` ->
`out_gtdb/gtdb_{cells.csv,summary.json}` (b000-b003) and `out_gtdb_batches/b004..b007/`;
`python mag_validation/combine_gtdb.py` -> `out_gtdb/gtdb_{cells,summary}_2000.*`; the
k-means and random arms via their `gtdb_arm` runs -> `mag_validation_{km,rnd}/out_gtdb/`;
`three_panel/compare_three_panels_mag.py` -> `three_panel/three_panel_mag_gtdb.json`.

`score_gtdb.py` reads only ANI, AF_query and AF_ref (no taxonomy mapping, no PGFams).
**Nothing to re-run.** Search cost if ever needed: about 150-180 s per 1,250-MAG batch
(GTDB index load included), 20-40 min per panel (logs).

---------------------------------------------------------------------------------------------

## D. The 291 real MAGs (EmilyKin; `mag_skani_hits.py`)

Commands (`run_emilykin_mag_skani.sh`): `mag_skani_hits.py sketch` (the sludge DB, 2 s) and
`mag_skani_hits.py search ... --db-kind bvbrc|gtdb --topk 20 --resume` (GTDB: 119 s search
plus about 104 s marker load and indexing, 60 threads; `run_stats.gtdb_r226.json`) ->
`~/Documents/EmilyKin/mag_skani_hits/` (hits JSON, summary, concordance, raw TSV). The
selections `mag_reference_selection.{gtdb_r226,bvbrc_sludge}.json` (2026-10-16) carry
`gene_data: "estimated"`, `select_mode: exact-tie`. In union mode the selection reads no gene
data, so they would be identical with or without a provider.

**Nothing reads BV-BRC gene data. Nothing to re-run.** If a determinism re-run is wanted:
* write to a scratch `--outdir`, not `~/Documents/EmilyKin/mag_skani_hits/` (outside prFBA;
  `--resume` may skip work);
* re-select with `select_references.py --hits ... --out <scratch> --knobs-json ani_knobs.json
  --gene-provider estimate`. The CLI default `--gene-provider auto` would, run from the prFBA
  directory, load the 5 GB shared cache whole (~28-40 GB) and rewrite it atomically on any
  fetch.

Stale copy (minor): `~/Documents/sludge_genomes/1733.9522.fna` and `1733.9680.fna` hold
exactly 25,000 contigs (same cap). The real MAGs have 0 hits to them
(`mag_vs_bvbrc_sludge_skani_raw.tsv`), so no reported number is affected.

Related, unconsumed: `reference_maps/gtdb_r226_acc_to_bvbrc.*` (2026-10-20) used "non-empty
PGFam entry" as tie-break step 5 and was built before the repair. Nothing reported reads it
(the manuscript says the map is not yet integrated). It is stale but costly to rebuild
(BV-BRC queries). Leave it and note it.

---------------------------------------------------------------------------------------------

## E. Depth-100 deployments (`deployment_d100/`)

### E.1 Builds in play

| build | dir | status |
|---|---|---|
| shipped (mixed PATRIC+RefSeq, k/n; 2026-10-29/30) | `~/Documents/{codiffusion_bioreactor,EmilyKin}/synthetic_genomes_d100/` | superseded; in KBase as version 1 (GenomeSets 261948/6467/1, 275568/3757/1); reproducible only at or before `3ce12d0`; `upload_kbase_genomes.py` refuses it |
| BV-BRC-only E3 (2026-11-13) | `.../synthetic_genomes_bvbrc_d100/` (16.7 / 19.4 GiB) | current: from the repaired caches, re-selected from the shipped d100 hits; not uploaded |
| BV-BRC-only k/n arms | `deployment_d100/downstream_e3/kn_builds/{codiffusion,emilykin}` (git-ignored) | comparison arms for downstream_e3 |
| mixed E3 re-weightings | `synthetic_genomes_d100_e3` | deleted 2026-11-13 |

### E.2 Command sequence (per deployment; `d100_pipeline.sh` via `run_{codif,emilykin}_d100.sh`)

| # | stage | command (as run) | output | BV-BRC data |
|---|---|---|---|---|
| 1 | align | `align_hits.py --force-backend cuda --k2 100 --prefilter-k 500 --gpu-cand 200 --workers 12 --save-shortlist --fasta F --taxonomy-csv T --outdir OUT` | `OUT/asv_top20_alignment_hits.json` (+ shortlist, `method_selection.json`) | 16S reference set only |
| 2 | tieset | `deployment_d100/tieset_gene_cache.py --hits H --out OUT/gene_families_tieset.json --force` (skipped when the tie-set file is newer than the hits) | small PGFam cache + `.tieset.json` sidecar | **shared PGFam cache** (streamed) |
| 3 | prefetch | `prefetch_gene_families.py --hits H --cache OUT/gene_families_tieset.json --workers 8` | adds missing selected genomes to the small cache | network for missing |
| 4 | build | `deployment_d100/build_d100.py --tieset GC.tieset.json --hits H --gene-cache GC --genome-cache-dir prFBA/bvbrc_cache/kbase_genome_cache --pfunc e3 --out-dir SYN --selection-out SYN/asv_reference_selection_d100.json --procs 8 --io-workers 8` (PYTHONHASHSEED=0, systemd-run MemoryMax=40G MemorySwapMax=0) | SYN (genomes, provenance, manifest, `pfunc/pgfam_store.json`, `pfunc/pgfam_fetched.json`, `pfunc/source_digests.json`, audit) | tie-set cache (selection), **shared PGFam cache + `e3_clade_catalogue.json`** (E3 store; streamed, never written), **kbase_genome_cache** (sources; missing ones fetched and added), network for PGFam sets not in the shared cache |
| 5 | stats | `deployment_d100/deployment_stats.py --selection SEL --hits H --provenance-dir SYN/provenance --out SYN/deployment_stats_d100.json --force` | stats | - |
| 6 | provenance | `write_mapping_provenance.py --hits H --selection SEL --genome-dir SYN --asv-fasta ... --taxonomy-csv ... --alignment-stats OUT/method_selection.json --gene-cache GC` (dataset names in `bvbrc_rebuild_summary.json` `deployments.*.mapping_provenance`) | `SYN/mapping_provenance.json`, `asv_mapping_summary.csv`, `abstained_asvs.{csv,md}` | - |
| 7 | k/n arm | `build_d100.py --tieset ... --selection SYN/asv_reference_selection_d100.json --gene-cache GC --genome-cache-dir ... --pfunc kn --out-dir downstream_e3/kn_builds/<dep> --procs 8 --io-workers 8` | k/n arm | genome cache |
| 8 | upload | `deployment_d100/kbase_upload_bvbrc_d100.sh` (DRY_RUN=1 first) | new object versions + `*_d100_bvbrc` GenomeSets in ws 261948 / 275568; `SYN/kbase_upload_log.json`; `deployment_d100/logs/kbase_upload_bvbrc_d100.log` | KBase network |

Committed records: `deployment_d100/bvbrc_rebuild_summary.json` (976b2c4; every number and
command of the 2026-11-13 builds), `deployment_d100/stats_d100_{codif,emilykin}.json`
(c91cf19; **the shipped builds**), `deployment_d100/compare/*.json` (`compare_d20_d100.py`,
shipped builds), `deployment_d100/e3_reweight_summary.json` (d9b05a1; mixed re-weighting,
superseded), `stats_d20/` (20-hit builds, historical).

### E.3 Re-run or reuse

* **Reuse** stage 1 (hits; align cost 2,379 s codif / 3,024 s EmilyKin on the GPU).
* **The BV-BRC-only builds already reflect the repaired data**: tie-set caches repaired
  in place (22 / 35 entries), shared PGFam cache, catalogue and genome cache repaired before
  they were built; `pfunc/pgfam_fetched.json` was fetched with the paged fetcher; the
  uploader's feature gate passes both.
* **Code**: the EmilyKin build's recorded hashes of `pfunc_estimator.py`,
  `build_synthetic_genomes.py`, `build_d100.py`, `build_synthetic_genomes_parallel.py`
  all equal HEAD (checked here). The fifth module in the chain is `prefetch_gene_families.py`;
  the E3 store's fetches use it. Its recorded hash (`7a70dfba...`, `how.prefetch_gene_families_py_sha256`)
  is not HEAD's (`1911983c...`). The recorded version survives as
  `/tmp/claude-1000/-home-freiburger-Documents-prFBA/400ec395-f172-59da-a98d-3dbe6b3a4256/scratchpad/paging_mut/backup/prefetch_gene_families.py`,
  and the diff against HEAD is docstring-only (29 diff lines, two docstrings), so the code
  is the same. Codiffusion's first three differ from HEAD: it predates
  the own-genome `fig|<gid>.` rule and the `load_source_genome` refactor (all in `8cfc297`).
  A scan found no foreign id in it, so a HEAD rebuild should give identical genomes and
  differ only in `_meta` stamps. That is equivalence by scan, not by rebuild.
  Decision for the executor:
  - (a) keep it (recommended); or
  - (b) rebuild into a NEW SYN, because the existing SYN is refused, exit 2. This needs
    either repointing `run_codif_d100.sh`, `kbase_upload_bvbrc_d100.sh` and
    `downstream_e3/common.py` at the new SYN, or moving the old SYN aside first, as was done
    for EmilyKin. Then diff every per-ASV file. Cost: about 9 min for E3 and 5 min for k/n,
    17 GB disk, and the E3 store's 266 PGFam fetches from BV-BRC (hold them to 5/s or
    less). If any genome changes, downstream_e3 parts 2-4 must be re-run.
* **Committed stats are the shipped builds'**: `stats_d100_{codif,emilykin}.json` differ
  from `SYN/deployment_stats_d100.json` only in `selection_rule.p3_rep_is_lexicographic_smallest`
  (545 -> 547, 590 -> 591), the selection and provenance paths, and wall time (checked here).
  Every manuscript deployment number is unchanged. Read from the committed
  `stats_d100_emilykin.json`:
  * answered 3,755, abstained 195;
  * n_genomes 1 / 2 / 3+ = 2,345 / 597 / 813, mean 2.2895;
  * cross-phylum 145/454 (31.9%) and 488/784 (62.2%);
  * all_multi resolved 1,238 of 1,410, bacteria_archaea 20. To make the repo
  current, re-run step 5 on SYN (about 30 s each) into the two committed paths, or add
  `stats_d100_*_bvbrc.json` beside them.
* **Tie-set caches vs the current shared cache** (measured here): every entry present in
  both is equal (EmilyKin 6,709, codiffusion 7,418). But 7 / 5 tie-set genomes absent from
  the tie-set files are now in the shared cache (150054.82, 3039803.773, 316.1475,
  39486.180, 39486.313, 85698.504, 91750.3255; codiffusion: 3039803.773, 316.1475,
  85698.504, 91750.3255, 91750.3960). And 611 / 702 tie-set genomes are absent from both;
  P3 counts them as 0 genes. EmilyKin's tie-set file also holds 19 entries that the shared
  cache lacks: prefetch wrote them into the small cache on 2026-10-29, and the repair audit
  covered them. `open_issues` in `bvbrc_rebuild_summary.json`: adding the 7 / 5 changes no
  selected genome, only `est_genes` / P2 order in 8 / 6 records. A fresh extraction is
  optional. If done:
  - move `OUT/gene_families_tieset.json` aside (the pipeline reuses it when newer than the
    hits);
  - expect the 19 tie-set-only EmilyKin entries to be dropped and partly re-fetched by
    prefetch;
  - build into a new SYN.
  Fetching the 611 / 702 never-cached genomes would be new data acquisition (P3 inputs),
  not a repair.
* `deployment_d100/compare/` (`compare_d20_d100.py --out-dir deployment_d100/compare
  --scratch <dir>`; offline, network disabled): its control selections read the tie-set
  cache, now repaired, so their P3-dependent fields would shift by the repaired entries. The
  manuscript uses only its hit-concordance numbers (219 of 3,950 ASVs, 1,696 hits), which
  do not depend on gene data. The script is pinned to the shipped build dirs for
  validation. Re-run is optional.

### E.4 Runtime, memory, network (measured)

| stage | codiffusion | EmilyKin | memory | network |
|---|---|---|---|---|
| align (2026-10-29) | 2,379 s | 3,024 s | cgroup peak 5.6 GB / 26.6 GB (chain) | no |
| tieset | 20 s | 21 s | ~0.9 GB (smoke) | no |
| prefetch | 7 s (shipped) / 10 s (2026-11-13, 0 to fetch) | 16 s / 10 s | ~2 GB cgroup | only for missing |
| build E3 | 483 s stage; 526.5 s unit | 480 s stage; 522 s unit | anon 10.8 / 10.0 GB (2 s samples); cgroup incl. page cache 20.9 / 31.8 GB; MemoryMax=40G | 3 genomes fetched to the genome cache (codif), 0 (EmilyKin); 266 / 306 PGFam sets for the E3 store |
| build k/n | 285 s | 290 s | anon 2.5 / 2.2 GB; cgroup 12.4 / 21.5 GB | 0 |
| stats | 33 s | 32 s | small | no |
| catalogue rebuild (if the cache changes) | 45-65 s, 2.7 GB | | | no |
| KBase upload (2026-11-04, mixed builds) | 205 min, 15.9 GB, 172 calls | 247 min, 19.5 GB, 205 calls | `MemoryMax=16G` planned | KBase |

The BV-BRC-only per-ASV JSONs hold about half the features, so the upload should take less
than the mixed one (an estimate; the bvbrc builds have not been uploaded).

### E.5 Stale copies of BV-BRC data (deployments)

| copy | state | action |
|---|---|---|
| `~/Documents/*/bvbrc_alignment_hits_d100/gene_families_tieset.json` | repaired in place 2026-11-13 (old values in the corrections log); 7 / 5 newly cached genomes not added | optional fresh extraction (E.3) |
| `.../gene_families_tieset.json.tieset.json` | still stamps the pre-repair source cache (size 5,085,756,050, 159,573 entries) plus a `corrections_2026-11-13` block | rewritten by a fresh extraction |
| `~/Documents/*/bvbrc_alignment_hits_d100/asv_reference_selection_d100.json` + `deployment_stats_d100.json` | shipped audit/stats from pre-repair caches; protected (the shipped manifest points at them) | leave |
| `synthetic_genomes_d100/` (both) + KBase v1 objects | mixed annotation | superseded; leave |
| `deployment_d100/stats_d100_*.json` | shipped-build stats | refresh from SYN (E.3) |
| `deployment_d100/compare/*.json` | shipped builds, pre-repair tie-set | optional |
| `deployment_d100/e3_reweight_summary.json`, `region_validation_d100/data/pfunc_e3_gate/production_ek_retro.json` | mixed builds (55.2% features kept k/n) | historical records; superseded by `bvbrc_rebuild_summary.json` (manuscript's "11%" = 415,207 / 3,841,468) |
| `SYN/pfunc/pgfam_store.json`, `pgfam_fetched.json` | made from repaired cache / paged fetcher | current |

### E.6 Hazards

* **Refusals by design**: the build stage exits 2 if SYN exists and is not empty, so
  `run_*_d100.sh` (SYN = `synthetic_genomes_bvbrc_d100`) are now refused. A rebuild needs a
  new SYN ending in `_d100`, or the old one moved aside first.
* **The catalogue gate**: an E3 build stops before selection if
  `bvbrc_cache/e3_clade_catalogue.json` does not match the shared cache's size and mtime.
  Any write to the shared cache (other sessions' prefetchers, `select_references.py
  --gene-provider auto`) forces a rebuild:
  `python pfunc_estimator.py rebuild-command` (45-65 s), and that changes E3 everywhere.
* **Shared genome cache**: the build adds new files to `bvbrc_cache/kbase_genome_cache/`
  (atomic `os.replace` via `build_d100.py`). That is additive, but concurrent builds share it.
* **BV-BRC rate**: the E3 store's fetch pool (8 threads, `--pfunc-fetch-workers`, not exposed
  by `d100_pipeline.sh`) has no rate limit of its own. The EmilyKin rebuild was held to 5/s
  or less only by a scratch `sitecustomize.py` shim.
* **Determinism**: rebuilding EmilyKin with an unchanged selection reproduced every other
  ASV's four files byte for byte, as well as the audit, the PGFam store and the fetched side
  file (`foreign_patric_id_rebuild`). What differs between runs: the manifest's run stamps
  (code hashes, `seconds_prepare`, the audit mtime in `selection_stamp`), and the stats'
  `wall_seconds`. Run with `PYTHONHASHSEED=0`, as the builds were.
* **Uncommitted uploader**: `upload_kbase_genomes.py` carries 324 insertions and 35
  deletions of uncommitted changes (the `--new-version`, `--prior-upload-log` and
  `--dataset-label` support the new script needs). `kbase_upload_bvbrc_d100.sh` is
  untracked, and it writes "Uploaded with prFBA $HEAD" (`git rev-parse --short HEAD`) into
  both GenomeSet descriptions. Commit both before the real upload, or the recorded commit
  misdescribes the code that uploaded.
* **KBase upload ordering**: `--new-version` stamps each object version with
  `genome_json_sha256`, and the GenomeSet references exact versions. Make every rebuild
  decision (codiffusion code stamp, tie-set refresh) before uploading. Otherwise each
  changed genome needs a third version and a new set. Run `DRY_RUN=1` first (it reads the
  token and both workspaces, saves nothing). The uploader's feature gate (8 processes,
  about 28 s per deployment) refuses any foreign or missing `PATRIC_id`. Its notes say
  version 1 of the EmilyKin objects carries a "codiffusion bioreactor, V4-V5" notes text;
  the new versions are fixed by `--dataset-label`. The decks
  (`presentations/make_deck.py` line ~719: "They are not yet in KBase") and
  `deployment_d100/README.md` section 11 need updating after the upload.
* **Repo hygiene**: the two SYN dirs are untracked and not git-ignored in the codiffusion and
  EmilyKin repos (`open_issues`). A `git add -A` there would stage 36 GB.

---------------------------------------------------------------------------------------------

## F. E3 panel gate (`region_validation_d100/scripts/pfunc_e3_panel_gate.py`)

### F.1 Command and outputs

```
systemd-run --user --unit=e3-gate-$(date +%s) --collect -p MemoryMax=40G -p MemorySwapMax=0 \
  -p WorkingDirectory=/home/freiburger/Documents/prFBA --wait --pipe \
  /home/freiburger/Documents/py_venv/bin/python \
  /home/freiburger/Documents/prFBA/region_validation_d100/scripts/pfunc_e3_panel_gate.py --workers 12
```
(defaults: `--bootstrap 2000 --seed 1729`; `--limit N` writes a smoke run to
`data/pfunc_e3_gate/smoke/`). Outputs: `region_validation_d100/data/pfunc_e3_gate.json`,
`data/pfunc_e3_gate/per_amplicon.csv.gz`, `data/pfunc_e3_gate/run_meta.json`. Then
`PFUNC_GATE_REVIEW_OUT=<dir> python region_validation_d100/scripts/pfunc_e3_gate_review_recompute.py`
(WORKERS=8; 550 s; reads `per_amplicon.csv.gz`, the hits, the shared cache, the catalogue) ->
copy its `indep_recompute.json` to `data/pfunc_e3_gate/review_recompute.json`. The text in
`region_validation_d100/PFUNC_E3_GATE.md`, the manuscript and the decks quotes these numbers.

Inputs: `region_validation_d100/data/selection.json` (key set and include n_selected),
`data/hits/asv_top20_alignment_hits.json`, the **shared PGFam cache** (streamed into an int32
store), **`bvbrc_cache/e3_clade_catalogue.json`**, and as reference cells
`data/gene_capture_cells.csv` (0.987) and `data/gene_capture_cells_tau0995.csv` (0.995).

### F.2 Re-run: yes, after the amplicon gene-capture re-score

The committed run (2026-11-12) used the pre-repair cache (size 5,092,217,865) and catalogue
(1,802 empty lists excluded; the rebuilt one excludes 963). The catalogue on disk is already
current against the repaired cache (rebuilt 2026-11-13 10:43), so no rebuild is needed now.
Expected effect (impact study): k/n gains most, because stale reps counted as gene-less in
n; the margin should shrink slightly from -0.0530 [-0.0547, -0.0513] Brier; a verdict flip is
very unlikely. About 675 more amplicons per mode become scorable.

### F.3 Runtime, memory, network

1,595 s on 12 workers; main maxrss 10.6 GB; cgroup peak 21.2 GB under MemoryMax=40G;
store 318 M ids, 6.85 M distinct PGFams. No network: neighbours missing from the cache are
only counted (32,021 genomes, `neighbours_absent_from_cache`). The review recompute takes 550 s.

### F.4 Hazards and validation

* **Reproduction check (hard dependency)**: `load_ref_cells` plus the `repro` counter. Any
  mismatch with the committed gene-capture cells makes the verdict "... (INVALID RUN ...)"
  and exits 3. On the repaired cache the old tables differ in 2,343 (0.987) and 2,171
  (0.995) rows, so re-score `region_validation_d100` gene capture first and commit or write
  those tables. A valid run then doubles as a cross-check of that re-score
  (`reproduction_vs_gene_capture_cells` must read `match` only, 135,982 rows).
* **Catalogue gate**: exits if `catalogue_is_current` is false. Keep the shared cache frozen
  from catalogue build to gate end.
* **Determinism**: the rows of `per_amplicon.csv.gz` are sorted and the bootstrap clusters
  re-indexed by sorted name, with a fixed seed. So the JSON `results` blocks should
  reproduce on unchanged inputs. `_meta` (timestamps, rss) differs, and the gzip header
  holds an mtime, so compare decompressed content, not file bytes.
* `production_ek_retro.json` is a record of the mixed build; do not regenerate it.
  `reweight_build_e3.py` now refuses unstamped builds.

---------------------------------------------------------------------------------------------

## G. `deployment_d100/downstream_e3`

### G.1 Commands (`downstream_e3/README.md`, "Reproduce")

```
systemd-run --user --unit=ds-e3-<name>-$(date +%s) --collect -p MemoryMax=40G -p MemorySwapMax=0 \
  -p WorkingDirectory=/home/freiburger/Documents/prFBA -E PYTHONDONTWRITEBYTECODE=1 -E PYTHONHASHSEED=0 \
  /home/freiburger/Documents/py_venv/bin/python -W ignore \
  /home/freiburger/Documents/prFBA/deployment_d100/downstream_e3/<script> --builds bvbrc <args>
```

| script | output (`results_bvbrc/`, `models_bvbrc/`) | time (BV-BRC-only run) |
|---|---|---|
| `census.py --builds bvbrc --scan` | `census.json`, `census_sweep.json` | 32 s |
| `penalties.py --builds bvbrc --workers 8` | `penalties.json`, `part2_work/` | 1,190 s |
| `flux.py --builds bvbrc --workers 8` | `flux.json`, `part3/{deployment}/{asv}.json` | 2,157 s; peak worker RSS 1.91 GB (mixed run 3.99 GB) |
| `vocab.py --builds bvbrc --against mixed --workers 8` | `vocab.json`, `part4/` | 252 s |

Inputs: the E3 builds (`synthetic_genomes_bvbrc_d100`), the k/n arms
(`downstream_e3/kn_builds/*`), and pinned `git archive` exports of ModelSEEDpy 09856f5,
MSCommunity 03bd80c, ModelSEEDTemplates b6fe86b / 21a8084 and ModelSEEDDatabase 5389554a.
The exports live in `PINS` (default: this session's scratchpad `downstream/pins`), and
`common.ensure_pins()` recreates them from those commits.

### G.2 Re-run or reuse

Reuse. `results_bvbrc/` was re-run (2026-11-13 11:48-12:30) on the builds from the repaired
caches. The later EmilyKin foreign-id rebuild changed only ASV 7ac9047b, which no part
models, and census `--scan` was re-run at 15:17. Re-run only if a deployment is rebuilt
(E.3 b, or a tie-set refresh), and then all four parts.

### G.3 Hazards

* Not byte-reproducible: GPR term order varies run to run (order-only GPR strings
  3,679 -> 3,846 in part 2, 7,737 -> 7,956 in part 3, for byte-identical genomes). Compare
  the summaries, not the files.
* The census table has 21 anchor failures and 3 unclassified sweep files, all caused by
  another session's EmilyKin edits. It needs refreshing; that is a documentation task, not
  data.
* `--builds mixed` cannot be re-run: its E3 dirs were deleted.
* The pins depend on the session scratchpad, which is not durable. Pass `DS_E3_PINS=<durable
  dir>`, or let `ensure_pins()` rebuild them; that needs the library repos to hold those
  commits.

---------------------------------------------------------------------------------------------

## H. Run order and dependency edges (this report's scope)

```
[other report] region_validation_d100 (+ _legacy, + depth-20 if in scope) gene capture
               re-scored on the repaired cache -> gene_capture_cells*.csv
        |                                   |
        v                                   v
F  E3 gate (needs the new cells for its      A7/A8/A12  paired_comparison_d100*, crossover_d100*,
   reproduction check) -> review recompute              af_gate_counterfactual (+ depth-20 twins)
                                                        ^
A1/B  rebuild 3 panel indices (+ control index) --------+
      -> control re-score (byte-identity) -> real re-score of out_batches, tau99, n500, lowc,
         km, rnd, smoke -> concat -> aggregate -> paired/crossover/af-gate -> two_panel /
         three_panel MAG JSONs -> F14/F15 -> RESULTS_v2.md, manuscript sec:magval, decks
E  deployments: decide (a)/(b) and the tie-set refresh -> stats refresh -> [rebuild =>
   downstream_e3 parts 1-4] -> KBase upload (kbase_upload_bvbrc_d100.sh, DRY_RUN=1 first)
   -> README section 11, decks
C, D: nothing to run.
```

Steps A1/B, E and F use different resources and can overlap. The shared cache must stay
frozen while A1, F and any E3 build run.

Wall-clock estimate, excluding the amplicon re-score:
* panel indices: about 2 min;
* control plus real MAG re-scores, all sets: about 30-40 min at -P 5-6;
* comparisons and figures: minutes;
* E3 gate plus review: about 36 min;
* deployment stats: 1 min;
* optional codiffusion rebuild: about 15 min, plus downstream_e3 at about 60 min;
* KBase upload: probably 3-7 h for both deployments, by the 2026-11-04 rates.

---------------------------------------------------------------------------------------------

## I. Every stale copy of BV-BRC data in this scope

| location | content | stale because | fix | needed for the manuscript? |
|---|---|---|---|---|
| `mag_validation/pgfam_panel.json` | int PGFam sets, 10,000 genomes | 140 repaired genomes old | re-extract (A.4) | yes (MAG gene capture) |
| `kmeans_panel/pgfam_panel_kmeans.json` | same | 20 | re-extract | yes (three/two-panel MAG) |
| `random_panel/pgfam_panel_random.json` | same | 38 | re-extract | yes (three-panel MAG) |
| MAG batch cell tables and every derived MAG JSON/CSV | derived from the above | - | re-score | yes |
| `mag_validation_rnd/contigs/1733.9197.fna`, `1733.9535.fna` | contigs | 25,000-row cap (probable) | paged re-fetch, re-sketch, re-search (about 30 min) | no; caveat |
| `~/Documents/sludge_genomes/1733.9522.fna`, `1733.9680.fna` | contigs | same | same, for the sludge DB | no (0 hits) |
| deployment tie-set caches | PGFam lists | repaired; 7 / 5 newly cached genomes missing | optional fresh extraction | no (no selected genome changes) |
| `deployment_d100/stats_d100_*.json` | stats of the shipped builds | pre-repair selection | copy/re-run from SYN | numbers identical |
| `deployment_d100/compare/*.json` | control selections on the old tie-set | pre-repair tie-set | optional re-run | no |
| `region_validation_d100/data/pfunc_e3_gate*` | gate on pre-repair cache and catalogue | pre-repair | re-run after the amplicon re-score | yes |
| `reference_maps/gtdb_r226_acc_to_bvbrc.*` | accession map, PGFam tie-break | pre-repair | leave (unconsumed) | no |
| `e3_reweight_summary.json`, `production_ek_retro.json`, `synthetic_genomes_d100/`, KBase v1 objects | mixed builds | superseded | leave | no |

Not stale (checked): `bvbrc_cache/e3_clade_catalogue.json` (current against the repaired
cache), both `synthetic_genomes_bvbrc_d100/` (built after the repair; EmilyKin after the
foreign-id fix), their `pfunc/` stores, `downstream_e3/results_bvbrc/`, the GTDB arm, and
the real-MAG outputs.

---------------------------------------------------------------------------------------------

## J. Open questions for the orchestrator

1. Re-score the depth-20 pairings (`paired_comparison.json`, `paired_tau99.json`,
   `crossover.json`, `out/F_mag_*`)? They pair against
   `region_validation/data/gene_capture_cells.csv` (depth 20), which the manuscript no longer
   cites. Same question for the `_legacy` pairings.
2. Rebuild codiffusion so its manifest's code stamp is HEAD's (E.3 b), or accept
   equivalence by scan?
3. Refresh the deployments' tie-set caches from the current shared cache (7 / 5 genomes;
   only `est_genes` / P2 order change)?
4. Keep the 839 empty-list (poison) repairs? They drive most of the MAG and amplicon
   gene-capture change, and 2 of the 3 deployment selection changes
   (`bvbrc_rebuild_summary.json` `open_issues`).
5. Re-fetch the two random-panel genomes' contigs (25,000-row cap), or record a caveat?
6. Upload timing: after decisions 2-3, `kbase_upload_bvbrc_d100.sh`. Upload the k/n
   arms too, or only the two E3 builds? The script uploads only the two
   `synthetic_genomes_bvbrc_d100` builds.
