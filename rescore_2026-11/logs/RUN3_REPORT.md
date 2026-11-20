# Run 3 report: shared reference, Path A part 1 (R3.1–R3.5), 2026-11-16

Executor: Run 3 agent. Plan: `plan/PLAN.md` §0–4 (Run 3), §6, §8 and `plan/pipelines_shared_ref.md`.
Decisions: D1 Path A, D3/D11 pins per stage, D9 (historical phenotype passes not re-run).
Progress log with timestamps: `logs/RUN3_PROGRESS.md`. Every unit's wall time and memory: `stamps/units.jsonl`.
Controls, checks and W8 reports: `control/r3/`. Per-run logs: `logs/r3/`. Moved-aside and copied pre-repair
outputs: `pre_repair/region_validation/shared_reference/` (and `pre_repair/.../kfold_audit/panel_ec_sets.json`).

**Outcome.**
- Every R3.1–R3.5 stage ran on FM, wall time 07:47–11:30.
- Every control the plan lists passed except one pair. **Genus pass 1 exits 3 because C11/C12 fail**: they are bound
  to the pre-repair phenotype summary (`data/phenotype_summary.json`, 10-16). This was predicted before the run,
  is by construction, and is diagnosed below.
  - That branch stops at pass 1, as the plan intends; R4.1/R4.2 clear it.
  - Every genus pass-1 number in this report is **provisional, pending pass 2**. Pass 2 must reproduce them, and
    it reuses the cells table by fingerprint.
- One W8 failure on the first pass was a gap in the attribution tool (the accessory vocabulary), not in the data.
  It was diagnosed and fixed (§ "Deviations", item 5).
- The cache was never written. `stamp_state.py --compare R1.6`: CACHE UNCHANGED at the start (`stamps/R3.0.json`)
  and at the end (`stamps/R3.final.json`).
- TF ran no stage. The measured local runtime of R3.4 is under the 45-minute threshold (arithmetic in
  "Machine split"). TF's code was fast-forwarded to HEAD at the end, ready for R4.3.
- Commits are path-limited and not pushed (see "Commits").

**Headlines** (genus provisional). Leave-genus-out prfba_nestedcv 0.6867 vs PICRUSt2 0.7787, Δ −0.092
[−0.1073, −0.0812]: unchanged at 4 dp. Leave-family-out Δ −0.1235 → −0.1234. C2 Δ 0.0928 → 0.0927.
- `collect.py` over the 27 R3/R3F items: 360 values, every old value matches the print, 0 errors, 34 pending
  (W6, Run 5), 17 changed at print precision.
- Of those 17, only 7 are shared-reference values: SR-04 −0.019 → −0.020, SR-10 0.782 → 0.783 and 0.803 → 0.804,
  family FR-02/FR-03 ×4. The other 10 are R2.4/R2.5 as-shipped/dialect values these items also cite.
- **Optional aligner proof.** Today's aligner re-aligned fold f00_d100 to a hits file byte-identical to the
  10-14 one.

---------------------------------------------------------------------------------------------------------------

## Machine split (task rule: R3.4 on TF only if its local runtime exceeds ~45 min and the transfer is well under that)

- **Measured.** The f00_d100 control fold (the largest: 6,788 refs, 12,773 queries) took 134 s wall at 10.0 GB RSS
  on the full streaming path. Selection is single-threaded; `--workers 30` only feeds `align_hits`, which is skipped.
- **Extrapolated.** About 65 min serial for the 45 runs, scaling by hits size:
  - genus: 20 runs over 26 GB of hits;
  - family: 20 runs over 6 GB;
  - C2: 5 runs over 6.5 GB.
  That is about 32 min at the plan's 2-way and 15–20 min at 5-way. FM had ~200 GB free, so 5-way fits in memory.
- **TF route cost.** TF lacked:
  - the family and C2 hits (12.5 GB, ~0.9 GB zstd);
  - the new GFS;
  - all ref manifests and `panel_md5_{seq,id}.json`.
  The predictions and subsets coming back (~15 GB, ~5 GB zstd at the measured 3.2× for predictions) would need
  a checksummed pull, plus a one-fold equivalence proof.
- **Decision.** The transfer alone (~10 min) is not "well under" a 15–30 min local runtime, so R3.4 ran on FM
  beside R3.3: 40 genus+family runs at 5-way (1,110 s) and 5 C2 runs at 2-way (575 s), concurrent with PICRUSt2.
- **TF code.** `tf_code_sync.sh` was run at the end, so TF's clone is at FM HEAD for R4.3 (see "Run 4").

---------------------------------------------------------------------------------------------------------------

## Stages, controls and checks

### R3.1 W3 answer-key cache (`rescore-r31-w3`, 65 s, 0.9 GiB)

| check | result |
|---|---|
| **W3 control** on the OLD GFS (`pre_repair/…/genome_function_sets.json`, 6f4f817c) | `data/` and `data_family/` caches **equal in content, key order and bytes** to the on-disk `_cache_panel_ec.json` (185ca60f) |
| real refresh from the new GFS (9463ed94) | 9,762 keys (set unchanged). **140 changed = `changed_sets.json` changed_panel**; 17 eligible = changed_eligible; 0 lost any EC; 0 new ECs outside the 2,235-EC universe (eligible or not) |
| guarded files | `panel_eligible.json` / `panel_ec_universe.json` sha256 **and inode** unchanged (data_family hard links intact) |
| gained ECs | 9–120 (plan: 9–114). The difference is one genome: 1503470.4 gains 112 → 120. Its 183 PGFams that could not be projected on 11-15 were fetched in R2.5. Every other genome equals `expected/rescore_work/universe_check.json` exactly |
| order hazard (#14/#15) | `kfold_audit/panel_ec_sets.json` moved to `pre_repair` **after** both caches were regenerated |
| data vs data_family | the new caches are identical (8eb296b1) |

### R3.2 references (`rescore-r32-refs`, 1,815 s wall, 12 concurrent, cgroup 17.8 GiB)

- **Affected references** (`control/r3/r32_affected_refs.json`): genus **23 of 30**, exactly the measured set. The 7
  untouched ones are f02/f03/f04 arc d025/d050 and f04_arc_d100. Family: 25 of 30 (untouched: f00 arc ×3, f03 arc
  d025/d050). C2: 2 of 2.
- **Procedure.** `manifest.json` and `traits_EC.txt.gz` were moved to `pre_repair`; the other rewritten files (marker
  table, `.model`, probe FASTA, `hsp_screen/`) were copied there. Then `build_reference.py` was re-run with the
  ledger argv, without `--force` (C2: the `run_c2_include_self.sh` argv).
- **Build times.** Median 45 s, max 1,196 s; the bac_d100 descriptives are the long pole.
- **Checks: `r32_checks.py`, 50/50 OK** (`control/r3/r32_checks.json`):
  - changed files only `traits_EC.txt.gz`, `manifest.json`, `traits_16S.txt.gz` (gzip mtime only; decompressed
    equal) and `hsp_screen/*`;
  - `.tre`, `.fna`, `.hmm`, `.probe.tre`, `.model` and `probe_queries.fna` **byte-identical**; FastTree reused in
    every build;
  - `tip_labels_sha256`, `model_string`, `reference_genomes`, `n_tips`, `msa_width` and `set_equality_ok` (C4)
    equal; tree descriptives (MPD, monophyly) equal;
  - decompressed trait rows differ **only for the 17 genomes**: 284 rows (genus 123, family 144, C2 17);
  - 5 refs gained one trait column through the prevalence screen (2.5.1.114 ×4, 3.2.1.39 ×1); no crash-screen
    retry anywhere;
  - every build logged "EC sets from <panel dir>/_cache_panel_ec.json". None derived EC sets or fell through, and no
    `panel_ec_sets.json` was written;
  - the unaffected refs are untouched (0 files changed).

### R3.3 PICRUSt2 (68 runs: genus 23 forced + 5 nsti, family 25 + 5, C2 10; 4×14 and 2×14)

Units: `rescore-r33-p2-gf` (1,926 s, 38.5 GiB) and `rescore-r33-p2-c2` (932 s, 19.1 GiB). Argv: the ledger stage-5
rows and `run_c2_include_self.sh`, without `--force`.

**`r33_checks.py`: 68/68 OK** (`control/r3/r33_checks.json`):
- C3, C3_preflight and C4 pass in every run; C7 reported.
- **`placed.tre` is byte-identical to the pre-repair copy in every domain workdir.** Placement re-ran in all 68,
  because the stage stamp hashes `traits_sha256`, so epa-ng determinism is now shown on 68 re-placements.
- `stage_place.json` differs only in `traits_sha256`.
- Unaffected runs: 0 files changed (stamp taken before any run vs after).
- Query EC predictions changed: genus 5,003 queries (21 of 28 runs), family 2,051 (18 of 30), C2 628 (6 of 10).

### R3.4 prFBA re-select (45 runs; ledger stage-6 / `run_c2` argv, no `--force`)

**Control first, f00_d100 on the PRE-REPAIR inputs** (`rescore-r34-ctl`, 155 s; `control/r3/prfba_ctl/`). It ran in
scratch outdirs with copies of the hits and the pre-repair manifests, so the canonical directory was untouched.

| variant | inputs | result |
|---|---|---|
| (a) plan's literal control | old `_cache_ec_subset` / `_cache_pgfam_subset` kept (reuse path) | predictions **identical in content** to the committed run (all 4 arms × 12,773 cells; meta equal minus generated/git_commit/wall/align seconds) |
| (b) streaming on the pre-repair snapshot | 09-17 GFS (`pre_repair`) + SNAP f0a84bf1, no subsets | predictions **identical**; the streamed subsets **equal the 10-14 subsets** |

**Real runs.**
- `r34_one.sh`, per run:
  - pre-checks that `align_stamp.json` matches `ref_db.json`/`queries.fna`;
  - moves `predictions.json` + both subsets to `pre_repair`;
  - makes `hits/` read-only, so an align-stamp mismatch fails instead of re-aligning;
  - requires `[skip] reusing` and no `[align]`.
- Units: `rescore-r34-prfba-gf` (40 runs, 5-way, 1,110 s, 54.6 GiB) and `rescore-r34-prfba-c2` (5 runs, 2-way, 575 s).

**`r34_checks.py`: 45/45 OK** (`control/r3/r34_checks.json`, per-run reports in `control/r3/r34_cmp/`):
- EXIT 0; hits reused in every run (no re-alignment).
- The rewritten `ref_db`, `ref_id_db`, `queries`, `midas` and `align_stamp` files are byte-identical.
- C6/C8 canaries clean in every fold run; `hits_outside_reference` 0; `reps_outside_reference` 0.
- Forced-answer contract held in all 10 forced runs.
- C2 hard gate: exact_tie source recovery 1.0 in all 5 folds. Legacy and graduated recovery are identical to the
  10-20 run.
- Meta equal minus run fields; assertions (`gene_count_fn_calls`, k2-saturation, C6 counters) unchanged.
- **Subsets differ only in the changed genomes, and the per-run count equals `fold_pgfam_staleness.json` in all 45**
  (1–17 per run). The EC subsets differ in the same genomes.
- **Predictions W8: 0 unattributed.** Every changed cell has a changed genome among its old or new reps.
  - Changed prediction cells: genus 3,172, family 1,979, C2 538.
  - ECs: 191,336 gained, 851 lost; the losses come from rep swaps, all attributed.

### R3.5 scoring (pins of `c13_finish.sh` on every scorer: `--templates ModelSEEDTemplates-main2294299`, `PYTHONPATH=ModelSEEDpy-d884a50`)

Before scoring, `r35_backup.sh` copied every scorer output to `pre_repair`: the three cells tables and their
`.fp`, the summaries, both `controls.json`, `summary_recipe` and `cells_recipe`. These are the W8 input, and
they include the old SR `cells.csv.gz` (Run 1 handoff).

| step | unit | wall / peak | result |
|---|---|---|---|
| C1 argv reconstruction (scratch) | `rescore-r35-c1ctl` | 62 s | The argv was never logged. The reconstruction `--mode recipe_control --picrust2 data/picrust2 --refs data/refs --truth data/refs/recipe_control_truth.json` reproduces `cells_recipe` (decompressed) and every arm and contrast of `summary_recipe.json`. Its accounting shows the pinned templates are **byte-identical** to what the 10-14 scorer read (only the path differs), and the answer key moves only through functions/maps/ceiling |
| C1 (real) | `rescore-r35-c1` | 56 s / 5.1 GiB | C1 pass; numbers identical; provenance now carries answer key `c2b3c809cf33` |
| C2 | `rescore-r35-c2` | 3,437 s / 26.9 GiB | rc 0; C2 pass; C3–C8, C10, C13–C15 pass; C1 seen as fresh; valid |
| family (`score_family.sh` argv + pins) | `rescore-r35-family` | 470 s / 23.5 GiB | rc 0; C3, C4, C6, C7, C8, C10, C13, C14, C15 pass; valid |
| genus pass 1 (7b argv + `--force` + pins) | `rescore-r35-genus` | 6,456 s / 41.5 GiB | **exit 3: C11, C12 fail** (defects `phenotype_bound_to_a_different_cells`, `…_functions` against `data/phenotype_summary.json`). C1–C10 and C13–C15 pass on the new answer key; C14 55/55 runs; `inputs_changed_during_scoring` empty; D22 tripwire gap 0.0 (tolerance 0.0005); inner nested-CV winner **legacy in 5 of 5** folds, as before; summary stamped INVALID by the scorer (provisional sentence kept) |

**Answer key.** `a9dc024df800` → **`c2b3c809cf33`** (functions 6f4f817c → 9463ed94, maps e9a8854f → 122f5669, ceiling
938c7d58 → 460d866e; templates, EC universe and fanout unchanged). Scorer code fp `89a47bdfb384` is unchanged.

**The `.fp` answer-key check** (`tools/fp_check.py` calls `score_folds.cells_reuse_verdict` with the inventory the
new summary records; `control/r3/fp_check_*.json`):
- the pre-repair cells table + sidecar is **refused** ("inventory fingerprint … ANSWER KEY a9dc… vs c2b3… changed");
- the new table is **reusable**;
- result for all three tables (C2, family, genus): old refused, new reusable.

**W8 on the cells tables** (`tools/attrib_sr_cells.py`; `control/r3/w8_cells_<root>.json`):

| table | rows | changed | unattributed | reasons |
|---|---|---|---|---|
| C2 `cells_incself.csv.gz` | 1,150,020 | 4,596 | **0** | source 2,340; picrust2 ref 1,194; prFBA prediction 640; constant/donor/clade ref set 320; nn1 102 |
| family `cells.csv.gz` | 969,936 | 32,020 | **0** | accessory vocab 12,459; constant/donor/clade ref set 7,759; picrust2 ref 4,402; prFBA prediction 4,154; source 2,376; nn1 550; random-abstain 320 |
| genus `cells.csv.gz` (pass 1) | 4,216,740 | 111,014 | **0** | accessory vocab 82,904 (moved in 10 of 30 fold × density × space); picrust2 ref 10,858; source 8,580; prFBA prediction 6,724; nn1 1,132; constant/donor/clade ref set 460; random-abstain 356 |

Strict invariants:
- PICRUSt2 rows on untouched refs must be identical unless the source changed: 0 violations.
- prFBA rows with unchanged predictions must be identical unless the source changed, apart from the
  accessory-vocabulary columns: 0 violations.

---------------------------------------------------------------------------------------------------------------

## Deviations from the plan text, and why

1. **R3.4 ran on FM, not TF.** See "Machine split" for the arithmetic.
2. **R3.3/R3.4 were split by root.** Genus+family started at 08:23, as soon as their 48 refs passed `r32_checks`.
   C2 started at 08:24, after its 2 refs. This kept C2 off the critical path (C2 must be scored before genus pass 1).
3. **Scoring order: C2 → C1 (during C2) → family (in parallel) → genus pass 1.**
   - The plan orders "optionally C1, then C2, then pass 1". C1 fingerprints every `data/picrust2` prediction, so it
     waited for genus R3.3. It ran (56 s) while the C2 scorer was still computing.
   - `score_folds` reads `controls.json` only at its final read-modify-write merge, so C2's merge saw the new C1,
     and C2's summary lists C1 as fresh.
   - C1 and C2 write disjoint control ids. Genus pass 1 started only after both.
4. **C1 was re-scored**, although the plan marks it optional. The reconstructed argv was first proven in scratch
   (above), so the C1 entry now carries the new answer key instead of `answer_key_matches=false`.
5. **W8 tool fix: the accessory vocabulary.**
   - The first family W8 pass left 8,013 rows unattributed and flagged 6,443 "prFBA rows changed with unchanged
     predictions". In every one of them, only `n_truth_accessory` and/or `f1_accessory` changed.
   - Cause: `score_folds.FoldContext` defines the accessory vocabulary as EC prevalence in [0.05, 0.95] over each
     fold × density reference set, which contains the 17 changed genomes.
   - Fix: `attrib_sr_cells.py` gained the reason `accessory_vocab_moved`. It is allowed only when those two columns
     are the only change **and** the recomputed old/new vocabulary for that fold × density × space actually differs.
     The recomputation uses the old/new `_cache_panel_ec.json`, the manifests and `ec_rxn`, with the scorer's own
     semantics.
   - On re-run: 0 unattributed, 0 strict violations. The vocabulary moved in 8 of 30 family cells.
6. **Pass 1 sees a phenotype summary that the 10-14 pass 1 did not.** `score_folds` defaults `--phenotype-summary`
   to the file beside `--cells`. `data/phenotype_summary.json` (10-16, bound to the pre-repair cells
   `ffc9d3a7…`) did not exist at the 10-14 pass 1. With the new cells, C11/C12 therefore fail on binding defects, and
   the summary is stamped INVALID with the provisional sentence kept.
   - The plan's literal command was run **without hiding the file**. Hiding it would have left the A-C13 C11/C12
     "pass" entries in `controls.json`, carrying the old answer key and the old cells. That is the misleading state.
   - Treated under the task rule "if a control fails, stop that branch, diagnose, report": nothing in Run 3 proceeds
     from the genus summary's validity.
   - The pass-1 summary was committed alone (`1a64f47`), so it can be reverted by itself if HEAD should keep the
     pre-repair summary until pass 2.
7. **Optional aligner proof (plan R3.5, optional).** It ran on the idle GPU beside pass 1 (`rescore-r3-alignreplay`,
   2,449 s):
   - command: `run_prfba_fold.py … --force` for f00_d100 into `$W/align_replay/sr_f00_d100`, with today's
     `align_hits.py` (last changed `15f3caa`, HEAD `4a2b378`);
   - `asv_top20_alignment_hits.json` is **byte-identical** to the 10-14 alignment, and `align_stamp.json` is
     identical;
   - the predictions equal the new canonical f00_d100 predictions;
   - alignment took 2,317 s vs 4,073 s on 10-14.
   This turns H10 ("the SR hits come from the 10-14 aligner") from an assumption into a demonstrated identity for
   that fold.

---------------------------------------------------------------------------------------------------------------

## Old → new headline numbers (`collect.py --only` the 27 R3/R3F items; `control/collect_run3.{json,md}`; the repo's `number_diff.*` is untouched)

- **Totals.** 27 items, 360 values: every old value equals the printed number (8 at a rounding boundary), 0 errors,
  34 pending (all W6 `derived:sr_*`/`as10_*`, Run 5), 99 changed, **17 changed at print precision**.
- **The headline sentence.** It now begins with `INVALID`. No extractor reads it, so nothing is affected; that
  changes when pass 2 is written.

Changed at print precision:

| item (tex lines) | value | printed | old → new | source stage |
|---|---|---|---|---|
| SR-04 (L76) | prfba_exact_tie_forced Archaea Δ | −0.020 | −0.019 → −0.020 (old at rounding boundary) | R3 genus (provisional) |
| SR-10 (L620-621) | prFBA precision, identity 0.945–0.987 | 0.783 | 0.782 → 0.783 (boundary) | R3 genus (provisional) |
| SR-10 (L620-621) | exact-tie precision, 0.945–0.987 | 0.803 | 0.803 → 0.804 | R3 genus (provisional) |
| FR-02 (L77, 624) | family PICRUSt2 pop F1 | 0.749 | 0.749 → 0.750 | R3F |
| FR-03 (L78, 599, 625-628, 1192) | family forced pop F1 | 0.7095 | 0.7095 → 0.7096 | R3F |
| FR-03 | family nn1 pop F1 | 0.7070 | 0.7070 → 0.7072 | R3F |
| FR-03 | family exact-tie forced pop F1 | 0.7072 | 0.7072 → 0.7074 | R3F |
| abs-dialect-density (L81), disc.vocab_breadth (L1189, 1217) | dialect ΔF1 BV-BRC / PICRUSt2 truth; amplicons | +0.046 / −0.202; 19,686 | +0.047 / −0.203; 19,843 | R2.4 (already in Run 2) |
| abs-region-penalty (L82) | as-shipped FL / V4 pop F1 | 0.768 / 0.668 | 0.776 / 0.681 | R2.4 |
| disc.margin_reversal (L1184-1216) | as-shipped deficit; prFBA pop; d20 deficit | −0.016; 0.713; −0.129 | −0.006; 0.723; −0.120 | R2.4/R2.5 |

**Unchanged at print precision.** The genus headline (L76-78: 0.687 vs 0.779, −0.092 [−0.107, −0.081]), the tab:sharedref
cells and the family headline (−0.123 [−0.150, −0.099]). The family CI hi moves −0.0991 → −0.099 inside the
same print; the C2 numbers (0.9914 / 0.8986) are unchanged.

---------------------------------------------------------------------------------------------------------------

## Timings and memory (cgroup `memory.peak` read inside the unit / GNU time max RSS)

| unit | stage | wall | cgroup peak | max RSS | MemoryMax |
|---|---|---|---|---|---|
| rescore-r31-w3 | R3.1 | 65 s | 0.9 GiB | 0.9 GiB | 16G |
| rescore-r34-ctl | R3.4 | 155 s | 23.6 GiB | 9.5 GiB | 60G |
| rescore-r35-c1ctl | R3.5 | 62 s | 6.6 GiB | 5.1 GiB | 20G |
| rescore-r32-refs | R3.2 | 1,815 s | 17.8 GiB | 0.5 GiB | 80G |
| rescore-r34-prfba-c2 | R3.4 | 575 s | 23.3 GiB | 9.5 GiB | 40G |
| rescore-r33-p2-c2 | R3.3 | 932 s | 19.1 GiB | 13.9 GiB | 60G |
| rescore-r34-prfba-gf | R3.4 | 1,110 s | 54.6 GiB | 9.5 GiB | 80G |
| rescore-r34-checks | R3.4 | 596 s | 29.7 GiB | 6.3 GiB | 40G |
| rescore-r33-p2-gf | R3.3 | 1,926 s | 38.5 GiB | 11.5 GiB | 100G |
| rescore-r35-c1 | R3.5 | 56 s | 5.1 GiB | 5.1 GiB | 20G |
| rescore-r35-family | R3.5 | 470 s | 23.5 GiB | 23.5 GiB | 60G |
| rescore-r35-w8-family | R3.5 | 33 s | 1.6 GiB | 1.6 GiB | 40G |
| rescore-r35-w8-family2 | R3.5 | 92 s | 2.3 GiB | 2.2 GiB | 40G |
| rescore-r35-c2 | R3.5 | 3,437 s | 26.9 GiB | 24.4 GiB | 60G |
| rescore-r35-w8-c2 | R3.5 | 75 s | 5.3 GiB | 5.2 GiB | 40G |
| rescore-r3-alignreplay | opt. | 2,449 s | 11.0 GiB | 9.5 GiB | 60G |
| rescore-r35-w8-genus | R3.5 | 213 s | 7.6 GiB | 7.6 GiB | 60G |
| rescore-r35-genus | R3.5 | 6,456 s | 41.5 GiB | 41.4 GiB | 100G (rc 3) |

- Concurrency: R3.3 genus+family (4×14) ran beside R3.4 genus+family (5-way) and the C2 lane (08:23–08:42). The
  four units' cgroup peaks sum to 135.5 GiB, an upper bound because the peaks were not simultaneous. The lowest
  MemAvailable sampled by any Run 3 unit was 134 GB.
- The user's `codif-v7fpost-*` units kept running throughout.
- The critical path was the include-self bac reference (1,193 s), then C2 scoring (3,437 s), then genus pass 1
  (6,456 s).

---------------------------------------------------------------------------------------------------------------

## Commits (prFBA, branch `hardware-aware-alignment-dispatch`, path-limited, not pushed)

| commit | stage | content |
|---|---|---|
| `4a2b378` | R3.1–R3.5 (C1, C2, family) | `data/summary_recipe.json`, `data_c2/summary_incself.json`, `data_family/summary_family.json`, `data_family/controls.json` |
| `1a64f47` | R3.5 genus pass 1 (provisional) | `data/summary_genus.json`, `data/controls.json`; stamps of the git-ignored `cells.csv.gz` (d57709e3), `.fp` (224f90fb), inventory fp (a32df47c) and `_cache_panel_ec.json` (8eb296b1) in the message |

Commit messages: `logs/commit_msgs/r3_c1_c2_family.txt`, `r3_genus_pass1.txt`. TF's clone was then fast-forwarded
to `1a64f47` (`tf_code_sync.sh`; its tracked tree was clean).

No code change in prFBA was needed. Work-dir tools are new in `$W/tools`:
- `r31_w3.sh`, `r32_refs.sh`, `r32_one.sh`, `r32_checks.py`
- `r33_p2.sh`, `r33_one.sh`, `r33_checks.py`
- `r34_ctl.sh`, `r34_one.sh`, `r34_all.sh`, `r34_checks.py`, `cmp_prfba_pred.py`
- `r35_backup.sh`, `r35_w8.sh`, `attrib_sr_cells.py`, `fp_check.py`

Untouched: `presentations/` (including the user-edits deck), ModelSEED checkouts, KBase, codiffusion/EmilyKin, the
paper repo.

---------------------------------------------------------------------------------------------------------------

## What Run 4 needs

1. **The genus cells table is final.**
   - Files: `data/cells.csv.gz` d57709e3, `.fp` 224f90fb, answer key `c2b3c809cf33`, inventory a32df47c.
   - `c13_finish.sh`'s pass-2 scorer has no `--force`, so pass 2 reuses the table by fingerprint and its arms must
     equal pass 1's.
   - R4.3 can therefore read it now, so the TF push can go before R4.2.
2. **R4.1 must bind to the new inputs.** The phenotype passes need `--cells data/cells.csv.gz` (d57709e3) and the
   new `--functions` (9463ed94), in fresh canonical directories (move `phenotype_{atp,gf}_c13` aside).
   **R4.2** must pass `--phenotype-summary data/phenotype_gf_c13/phenotype_summary.json` explicitly, as
   `c13_finish.sh` does; that is what clears C11/C12. Run only its `score_folds` line.
3. **`data/phenotype_summary.json` must be regenerated or replaced.** The file is from 10-16 (outdir
   `data/phenotype`, sha256 957f284b); plan Stage-0 #20 covers it. Until then, every scorer invocation that does
   not name a phenotype summary auto-reads it and fails C11/C12, and the deck (R5.3) reads it too.
4. **TF for R4.3.** TF's clone is at FM HEAD `1a64f47`. Push `tf/lists/push_R4.3_posthoc.txt` (`tf_push.sh`, then
   `tf_verify.sh`); it overwrites the pre-re-score copies of list D:
   - `SR/data/cells.csv.gz` + `.fp` (final, above);
   - `SR/data/_cache_panel_ec.json` (8eb296b1);
   - all 20 genus `data/prfba/*/predictions.json`, `_cache_ec_subset.json` and `_cache_pgfam_subset.json`
     (re-selected in R3.4; hits and stamps already on TF and unchanged);
   - R2.5's `P2C/data/{function_maps,pgfam_products,prfba_selection,genome_function_sets}.json` and
     `comparison_cells.csv.gz`.
   Wire estimate: about 3–4 GB compressed. Predictions compress only about 3.2× with zstd-3, not the 14× in the
   runbook manifest.
5. **Blockers for R4.3, still open:**
   - **`stage_o_weights_E3.json` has no producer** in the repo (only consumers: `launch_objective_tuning.sh`,
     `objective_wrapper_equivalence.py`). It needs a small scripted W-item after Stage E.
   - **D11 pins per stage:** wobj on MSD `e8446e3e`; Stage O, FVA and the wrapper on `5389554a` (pass `--msd` or
     switch TF's symlink). All three trees are on TF.
   - **Order correction** (runbook §4): calibrate ‖ wobj `--population-only`, then Stage E ‖ wobj main, then the E3
     weights, then Stage O, FVA, results, wrapper and verification.
   - wobj's barrier re-solve is CPU-dependent at 1e-14. Run the wobj main run on FM if bit-identity is required.
6. **Moved aside in `pre_repair/region_validation/shared_reference/`** (needed by R4/R5 diffs and W6 `--state pre`):
   - old manifests and traits;
   - PICRUSt2 run dirs, minus `intermediate/`;
   - prFBA predictions and subsets;
   - cells tables, summaries and controls;
   - `_cache_panel_ec.json` ×2;
   - `kfold_audit/panel_ec_sets.json`.
7. **Disk.** `$W/align_replay` (2.1 GB) and `$W/control/r3/prfba_ctl` (4.2 GB) can be deleted after review. `pre_repair/region_validation/shared_reference` is 24 GB.
