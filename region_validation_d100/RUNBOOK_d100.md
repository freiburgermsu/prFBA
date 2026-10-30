# Depth-100 rerun of the amplicon arm: runbook

This runbook re-runs the 16S amplicon benchmark with 100 retained hits per record (k2 = 100) instead of 20, on all three panels and under both selection rules. It does not change any depth-20 artefact the manuscript cites. The setup is done: trees built, code changes applied, smoke test passed. Nothing below has been launched.

Throughout, `PY=/home/freiburger/Documents/py_venv/bin/python` and `R=/home/freiburger/Documents/prFBA`.

## Decisions this follows

| | Decision |
|---|---|
| D1 | Records keep 100 hits, using `--k2 100 --gpu-cand 500 --prefilter-k 500`. The whole 500-candidate shortlist is Biopython-rescored, so a record depends only on its sequence. |
| D2 | All three panels are scored under both selection rules: the exact-tie union (`*_d100`) and the legacy reducer (`*_d100_legacy`). |
| D3 | Each panel is aligned separately: three sequential alignments. |
| D4 | The edlib shortlist and the GPU scores are saved during alignment (`--save-shortlist`). Any later depth then needs only the rescore stage. |

## Trees

Every tree is a copy of its depth-20 source:

- `scripts/` is copied without `__pycache__`, `*.bak*`, `*.orig` or `run_full_at_515.sh`. That last script would `mv` the *original* tree's `data/`.
- `region_panel{,_final}.json` is copied.
- `data/` holds byte copies (not links) of the seven front-load files: `selection_10k.json`, `src16s_10k.fasta`, `combined_amplicons.fasta`, `combined_amplicons.stats.json`, `combined_amplicons.unique.fasta`, `amplicon_md5_fanout.json` and `benchmark_truth.json`.
- `_common.BASE` points at the tree itself.

| Tree | Panel | Rule | Aligns? |
|---|---|---|---|
| `region_validation_d100` | original (taxonomy-stratified) | exact-tie union | yes |
| `region_validation_d100_legacy` | original | legacy | no: copies hits from `region_validation_d100` |
| `region_validation_kmeans_d100` | k-means | exact-tie union | yes |
| `region_validation_kmeans_d100_legacy` | k-means | legacy | no: copies hits |
| `region_validation_random_d100` | uniform random | exact-tie union | yes |
| `region_validation_random_d100_legacy` | uniform random | legacy | no: copies hits |

**What `run_validation.sh` does in every tree**

Environment knobs (see the script header):

- `STAGE=all|align|post`
- `NO_TMUX=1`: run in the foreground under systemd, with no tmux.
- `ALIGN_WORKERS` (default 40).
- `SELECT_MODE`: fixed per tree. A mismatch exits with code 2.
- `ALLOW_GPU_SHARE=1`

Built-in guards:

- **Refuses `--force`.**
- **Front-load is check-only.** A missing front-load file or a symlinked one aborts the run; it is never rebuilt.
- **Legacy trees never align.**
- **Aborts (exit 3) when `align_hits.py` or `run_prfba_fold.py` is already running.** The C2 campaign uses the GPU in its prFBA step.
- **Raises the fd limit.**
- **Checks `method_selection.json` before writing `align.done`.** The check is k2=100, gpu_cand=500, prefilter_k=500.
- **Resumes from a saved shortlist.** If a run dies after the edlib stage, the rerun uses `--shortlist-in`.

The selection rule is threaded through POST as follows:

| Stage | exact-tie trees | legacy trees |
|---|---|---|
| `prefetch_gene_families.py` | `--select-mode exact-tie` | `--select-mode legacy` |
| `select_references.py` | `--select-mode exact-tie` | `--select-mode legacy` |
| `prefetch_universe.py` | `--select-mode exact-tie` | `--select-mode legacy` |
| `close_exclude_rep_gap.py --tree <tree>` | (no flag) | `--legacy` |
| `score_taxacc.py` / `score_genecap.py` | (no flag) | `--legacy` |

`close_exclude_rep_gap.py` now runs after the universe prefetch. It writes a sentinel, `data/close_exclude_rep_gap.done`.

## Smoke test result (setup, 2026-10-20)

The setup ran the 316-md5 depth-probe FASTA with the new flags and `--save-shortlist`, using 32 workers at load average about 43:

```
$PY $R/align_hits.py --fasta $R/region_validation/picrust2_comparison/data/depth_probe/wiped_sample.fna \
  --outdir <scratch>/k100probe --prefilter-k 500 --k2 100 --gpu-cand 500 --workers 32 \
  --min-ref-len 200 --force-backend cuda --save-shortlist
```

- **Fresh run:** 207 s wall. Of that, edlib took 125.2 s, GPU 1.9 s and rescore 75.4 s. Peak RSS was 4.9 GB.
- **Comparison with the depth-500 probe:** all 316 of 316 records are identical to `hits_k500[:100]` in every field (score, identity, md5, genome, edlib distance, coordinates, lineage). There were no swaps inside ties: with gpu_cand 500 = k1, the final `(-score, dist, idx)` sort is deterministic.
- **Reload run:** `--shortlist-in <scratch>/k100probe` skipped edlib and GPU, took 76 s, and produced a byte-identical hits JSON (same sha256).
- **Resume test:** a 24-sequence scratch tree exercised `STAGE=align NO_TMUX=1` through fresh, skip (`align.done`) and resume (`--shortlist-in` picked up automatically). All three gave byte-identical output.

**Cost model from the probe:** edlib costs 12.7 worker-seconds per ASV under the current load. Rescore at (k2 100, gc 500) costs 7.6 worker-seconds per ASV. That rescore figure is an upper bound: the probe sequences sit in dense identity clusters, and the plan's model gives about 5.8.

At 40 workers, each panel takes about 5–5.5 h of edlib plus 1.5–3.2 h of rescore plus about 6 min of GPU, so **about 7–9 h per panel and about 21–27 h for all three**. Each panel's shortlist plus GPU scores take about 0.4 GB.

## 0. Before launch (lead)

```bash
PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA
CACHE=$R/bvbrc_cache/genome_gene_families.json
cd $R && sha256sum -c --quiet region_validation_d100/d20_baseline.sha256 && echo BASELINE_OK
# (1) commit the setup's code (see "Files changed" at the end) for provenance
# (2) snapshot the shared PGFam cache while NO writer runs (4.4 GB)
pgrep -af 'prefetch_|select_references|close_exclude_rep_gap|score_genecap|prfba_selection' ; ls $CACHE.tmp 2>/dev/null
cp -a --reflink=never $CACHE $R/bvbrc_cache/genome_gene_families.SNAPSHOT_before_d100.json
# (3) optional: write-protect exactly the 32 cited depth-20 files (not whole dirs:
#     picrust2_comparison/data/mag_arm is still written by other work)
awk '{print $2}' region_validation_d100/d20_baseline.sha256 | xargs chmod a-w
# (4) GPU must be free of the C2 campaign's prFBA folds (they use --force-backend cuda)
pgrep -af 'align_hits\.py|run_prfba_fold\.py' || echo GPU_FREE
```

Check two more things before launch:

- **Driver mismatch.** `nvidia-smi` reports a driver/library mismatch: module 595.84, library 595.99.02. CuPy works, and was verified inside a `systemd-run --user` unit. A pending reboot would kill the campaign, so either reboot first or not until it finishes.
- **Snapshot cleanup.** Old cache snapshots in `bvbrc_cache/` total about 22 GB. Clean them only if disk gets tight; 855 GB is free.

## 1. Alignments: three panels, sequential, GPU/CPU exclusive

There is one unit for all three alignments. It skips any panel whose `data/hits/raw/align.done` exists, and a panel that died after edlib resumes from its shortlist. Re-launching the same command is always safe.

```bash
systemd-run --user --unit=rv-d100-align --collect -p LimitNOFILE=1048576 --working-directory=$R \
  /bin/bash -c 'ALIGN_WORKERS=40 bash /home/freiburger/Documents/prFBA/region_validation_d100/run_align_all_d100.sh \
    >> /home/freiburger/Documents/prFBA/region_validation_d100/logs/align_all.log 2>&1'
```

Monitor it with any of:

- `systemctl --user status rv-d100-align`
- `tail -f $R/region_validation{,_kmeans,_random}_d100/data/hits/raw/align.log`
- the log modification times.

Each panel goes through these stages:

1. `prefilter n/N` (edlib, 5+ h)
2. `saved edlib shortlist` (from this point a crash costs only the rescore)
3. GPU stage, then `saved GPU scores`
4. `rescore n/N`
5. `[ok] alignment params k2=100 ...`, then `align.done`

Check each panel as soon as its `align.done` appears. These checks are read-only and use about 20 GB of RAM:

```bash
$PY $R/region_validation_d100/check_d100.py align region_validation_d100          # checks 1-4 (+ probe, check 3)
$PY $R/region_validation_d100/check_d100.py align region_validation_kmeans_d100
$PY $R/region_validation_d100/check_d100.py align region_validation_random_d100
```

## 2. Copy hits into the legacy trees

This happens automatically: `run_post_all_d100.sh` (step 3) runs it for each legacy tree. To run it by hand, it is idempotent, copies only (never links), writes the sentinel last, and verifies `cmp` and the link count:

```bash
bash $R/region_validation_d100/seed_legacy_hits_d100.sh            # every aligned union tree
$PY $R/region_validation_d100/check_d100.py copy region_validation_d100_legacy
```

## 3. POST stages: six trees, strictly one cache writer at a time

There is one unit for this step. It takes an `flock` (`region_validation_d100/.post_cache_writer.lock`), so a second POST driver refuses to start. The trees run in this order: original union, original legacy, k-means union, k-means legacy, random union, random legacy.

Each tree waits for its panel's `align.done`, polling every 5 min. The unit can therefore start as soon as the original panel has aligned, and its POST overlaps the next panel's alignment. Alignment never writes the cache.

After a tree's POST, the unit runs `check_d100.py post <tree>`:

- **The checks pass:** it writes `data/post.done`, and the tree is skipped on re-launch.
- **The checks fail:** it writes no `post.done`, carries on, and exits 1 at the end.

```bash
systemd-run --user --unit=rv-d100-post --collect -p LimitNOFILE=1048576 --working-directory=$R \
  /bin/bash -c 'bash /home/freiburger/Documents/prFBA/region_validation_d100/run_post_all_d100.sh \
    >> /home/freiburger/Documents/prFBA/region_validation_d100/logs/post_all.log 2>&1'
```

Per tree, POST runs these stages in order:

1. fan-out (about 4 GB written)
2. candidate prefetch (rule-matched)
3. `select_references`
4. `prefetch_pgfams` (truth + reps)
5. `prefetch_universe` (rule-matched)
6. gap close (rule-matched; fixed point)
7. `score_taxacc`
8. `score_genecap`
9. `mcnemar`
10. `build_region_meta`
11. renders

Expect the CPU stages to take about 1.5–2 h per tree at depth 100, plus BV-BRC fetching. The plan estimates 10k–30k new genomes in total at 2.5–3.9 genomes/s; the save cadence is now every 2000 fetches.

Do not run any other job that writes `bvbrc_cache/genome_gene_families.json` while this unit is active. That means another prefetch, a `--gene-provider bvbrc` selection, or another tree's `score_genecap`.

To restrict the list: `POST_TREES="region_validation_d100 region_validation_d100_legacy" bash .../run_post_all_d100.sh`.

## 4. PICRUSt2 as-shipped comparison, original panel, both rules

This step needs `post.done` in `region_validation_d100` and `region_validation_d100_legacy`. The legacy POST is what caches the legacy reducer's candidates, and `prfba_selection.py` is cache-only.

It writes only to `region_validation_d100/picrust2_comparison/{data,figures}`, and reads the depth-independent PICRUSt2 run read-only from `region_validation/picrust2_comparison/data/full_picrust2`.

`bvbrc_source_genomes.json` and `pgfam_products.json` were seeded by copy, so the BV-BRC fetches only fill the gap. This step reads the PGFam cache and never writes it, so it may overlap the remaining trees' POST.

```bash
systemd-run --user --unit=rv-d100-picrust2 --collect -p LimitNOFILE=1048576 --working-directory=$R \
  /bin/bash -c 'SEL_WORKERS=12 bash /home/freiburger/Documents/prFBA/region_validation_d100/run_picrust2_d100.sh \
    >> /home/freiburger/Documents/prFBA/region_validation_d100/logs/picrust2.log 2>&1'
```

The script runs these steps in order:

1. `prfba_selection.py --out prfba_selection.json`
2. `prfba_selection.py --legacy --out prfba_selection_legacy.json` (fork pool of 12 workers, about 15 GB of hits per worker copy-on-write)
3. `fetch_bvbrc_annotations.py genomes`
4. `fetch_bvbrc_annotations.py pgfams` for each selection
5. `function_maps.py`
6. `genome_function_sets.py --selection both`
7. `compare.py --picrust2 <original full_picrust2> --selection ... --legacy-selection ...`
8. `render_comparison.py`

## 5. Three-panel and the other cross-panel comparisons, both rules

This step needs `post.done` in all six trees. Every output is a new `_d100` / `_d100_legacy` file, and the cited depth-20 JSONs are never written.

```bash
systemd-run --user --unit=rv-d100-downstream --collect -p LimitNOFILE=1048576 --working-directory=$R \
  /bin/bash -c 'bash /home/freiburger/Documents/prFBA/region_validation_d100/run_downstream_d100.sh \
    >> /home/freiburger/Documents/prFBA/region_validation_d100/logs/downstream.log 2>&1'
```

It runs these, for S in `_d100` and `_d100_legacy`:

- `three_panel/compare_three_panels.py --suffix S` → `three_panel/three_panel_comparison{S}.json`
- `three_panel/emptied_strata.py --suffix S` → `three_panel/emptied_strata{S}.json`. Its `emptied_pct` is the deck's EMPTIED. Strata are 0–19 / 20–49 / 50–99 / 100.
- `two_panel/compare_panels.py` (path flags) → `two_panel/panel_comparison{S}.json`
- `two_panel/mechanism_depth.py --depth 100` → `two_panel/mechanism_depth{S}.json`. Its bins extend to 20–49 / 50–99 / 100+.
- `diagnose_self_recovery.py` in the original tree (`min(tie, 100)`) → that tree's `data/self_recovery_*`
- `mag_validation/compare_arms.py --amp <tree>/data/gene_capture_cells.csv` → `mag_validation/out/paired_comparison{S}.json`
- `mag_validation/crossover.py --amp ... --out mag_validation/out/crossover{S}.json <mag cells>`

**Verify before trusting the MAG-vs-16S rows.** The `--mag-cells` input (`mag_validation/out/mag_cells_all.csv`) and the crossover inputs (`mag_cells_all.csv` plus `lowc_out_batches/b*/mag_cells.csv`) were inferred from the tier sets in the cited `paired_comparison.json` and `crossover.json`, not from a recorded command. The cited depth-20 MAG comparisons and two-panel comparisons paired the original panel's legacy cells with k-means union cells. The `_d100_legacy` and `_d100` outputs are each single-rule.

## 6. Checks (section 8 of the plan)

```bash
$PY $R/region_validation_d100/check_d100.py all        # every check whose inputs exist; exit 1 on any FAIL
cd $R && sha256sum -c region_validation_d100/d20_baseline.sha256     # check 6, also run at the very end
```

1. Every record has 100 hits, and `params.k2=100`, `gpu_cand=500`, `prefilter_k=500`.
2. **Top-20 prefix versus the depth-20 record.**
   - The score vectors must be equal for 100% of md5s.
   - A hit may differ only where its score's tie group crosses GPU rank 100, i.e. more than 100 shortlist candidates score at or above it. The depth-20 run rescored only the GPU top 100, so such a tie was sampled arbitrarily. This is verified against the saved `gpu_scores.npy`.
3. The 316 probe md5s inside the original panel's full-scale record equal `hits_k500[:100]`.
4. Records whose 100th hit is still at identity ≥ 0.987. This uses the same definition as `plan_folds._truncation_violations` and is reported, not gated.
5. **POST checks.**
   - `_meta.n_empty_after_exclude` is reported.
   - `pgfam_missing == 0` in `gene_capture_cells.csv`.
   - `select_mode` in `selection.json`, the taxacc summary and the genecap summary must match the tree's rule (`exact-tie` / `default_union` versus `legacy`).
   - The gap-close sentinel must be present.
6. The depth-20 baseline sha256 values are unchanged. The baseline is a copy of the lead's file, at `region_validation_d100/d20_baseline.sha256`.
   - Amended 2026-10-29: three baseline entries were replaced on purpose:
     `three_panel/three_panel_comparison.json`, `three_panel_mag_gtdb.json` and
     `three_panel_mag_panelarm.json`.
   - The three-panel bootstraps had drawn phylum strata in hash-seed order, so their
     intervals moved between runs (UNION_VERIFICATION.md, note N2).
   - The scripts now sort the strata, and the files were regenerated from them. Rates and
     point estimates are byte-for-byte unchanged; only interval bounds moved.
   - The old files remain in git history.

## 7. Deck and manuscript (last)

- `presentations/make_deck.py --depth 100 [--rule legacy|union]` reads three inputs from the depth-100 trees:
  - taxacc/gene-capture from `region_validation_d100[_legacy]`;
  - `comparison_summary.json` from `region_validation_d100/picrust2_comparison/data`;
  - `three_panel_comparison_d100[_legacy].json`.

  The s04 alignment-hours stat comes from `region_validation_d100/data/hits/raw/method_selection.json`. The output is `prFBA_updates_2026-10_d100.pptx`.
- **Hand-entered deck constants.** The deck refuses to build until `presentations/deck_constants_d100.json` has every entry of `SELF_TOP20`, `SELF_SAT`, `SELF_CRED` and `EMPTIED` filled. The file says how to compute each: EMPTIED comes from `emptied_strata.py`; SELF_TOP20 is the `self_in_top20` share in `asv_region_concordance.csv` (include mode, one row per amplicon). SELF_SAT and SELF_CRED have no generator in the repo, and the three must sum to 100.
- **Slide prose still describes depth 20.** `make_deck` prints the list of slides to rewrite.
- `--depth 20`, the default, still produces the cited deck. The setup verified this: all 155 members of the pptx zip are identical except the embedded workbooks' `docProps/core.xml` timestamps.
- **Column names keep their depth-20 spelling.** `self_in_top20`, `asv_top20_alignment_hits.json` and the `top20` key keep their names by design (more than 50 consumers). At depth 100, `self_in_top20` means "in the top 100". The depth is recorded in `hits/raw/method_selection.json` `params.k2`, and `render_accuracy.py` reads it for its F6 labels.

## Files changed by the setup (uncommitted)

**Repo root and shared scripts (backward compatible: default behaviour unchanged):**

- `align_hits.py` (C1): adds `--save-shortlist` and `--shortlist-in DIR`.
  - The saved files are `edlib_top{k1}_shortlist.npz` (idx/dist int32, padded with -1, plus `n`), `gpu_scores.npy` and `shortlist_meta.json`. The meta file holds sha256 values of the ASV list and of the reference md5 order, plus the edlib timing.
  - A content mismatch on reload exits with code 2, as does using either flag on the CPU path.
  - The exit-2 guard for `--gpu-cand < --k2` is kept.
  - `method_selection.json` params gain `save_shortlist` / `shortlist_in`.
- `prefetch_gene_families.py` (C4/C6): adds `knobs_for_select_mode()`; `populate_cache(select_mode=...)`; `--select-mode exact-tie|legacy` (default exact-tie); `--save-every`. The default save cadence goes from 200 to 2000.
- `kmeans_panel/close_exclude_rep_gap.py` (C5): adds `--legacy`, which sets `score_genecap.LEGACY_SELECT = True`. It now validates `--tree` and logs the mode.
- `three_panel/compare_three_panels.py` (C7): adds `--suffix` / `--roots`. `--out` defaults to `three_panel_comparison{suffix}.json`, and `--roots` without `--suffix` requires `--out`. Roots are recorded in meta only on non-default runs.
- `mag_validation/compare_arms.py` and `mag_validation/crossover.py` (C7): add `--amp`. `crossover.py` also gains `--out` and now uses argparse; the MAG-cells files remain positional, so existing invocations are unchanged.
- `two_panel/mechanism_depth.py` (C7): adds `--original-root`, `--kmeans-root`, `--depth` (which sets the bins via `bins_for`; depth 20 gives the original bins) and `--out`. The default output is `mechanism_depth.json` at depth 20 and `mechanism_depth_d<N>.json` otherwise.
- `presentations/make_deck.py` (C7): adds `--depth`, `--rule` and `--constants`. Note that a concurrent session's commit `6fc6bd2` swept these edits in.
- **New files:**
  - `three_panel/emptied_strata.py` (C8: em2.py promoted; it reproduces em2.py exactly at depth 20, with EMPTIED 14.3 / 40.4 / 49.9);
  - `presentations/deck_constants_d100.json`;
  - this directory's `run_align_all_d100.sh`, `seed_legacy_hits_d100.sh`, `run_post_all_d100.sh`, `run_picrust2_d100.sh`, `run_downstream_d100.sh`, `check_d100.py` and `d20_baseline.sha256`.

**Per tree:**

- `scripts/_common.py`: `BASE` (C3).
- `scripts/run_validation.sh`: C2, generated from one template.
- `scripts/prefetch_universe.py`: `--select-mode` and `SAVE_EVERY = 2000` (C4/C6).
- `scripts/prefetch_pgfams.py`: `SAVE_EVERY = 2000` (C6).
- `scripts/diagnose_self_recovery.py`: `min(tie, 100)` (C8).
- `scripts/render_accuracy.py`: F6 labels follow `params.k2` (C8).

The depth-20 trees' scripts are untouched.

**Not changed:** `depth_probe_eval.py` already has `--hits/--keys/--out` (C7 needs nothing more), and it is probe-specific: its bins run to 500. `compare_variants.py`'s "top-20" labels are left as they are; that script belongs to the optional June variant studies.
