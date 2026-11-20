# TF runbook: running parts of the prFBA re-score on the second machine

Written 2026-11-15 by the TF preparation step. Scope assigned to TF by the re-score plan
(`plan/PLAN.md`): **R2.1** (E3 gate and review), the **MAG re-scoring of R2.2-R2.3**, and **R4.3**
(post-hoc studies: calibrate_pfunc, Stage E, weighted-objective test, Stage O with FVA, wrapper
equivalence, verification). Everything else stays on FM. Nothing of the re-score itself was run here.

FM = this machine (`frankensteinsmonster`, Threadripper 3970X, 32C/64T, Zen 2, AVX2, 251 GB).
TF = `ssh TF` (`andrew-MS-7E70`, Ryzen 9 9900X, 12C/24T, Zen 5, AVX-512, 91 GB, RTX 4060 Ti 8 GB,
user `andrew`, home `/home/andrew`, 678 GB free).

Helper scripts named below live in two places:
- FM side: `~/Documents/prFBA_rescore_2026-11/tf/bin/` (`tf_push.sh`, `tf_verify.sh`, `tf_pull.sh`,
  `tf_code_sync.sh`) and the file lists in `tf/lists/`.
- TF side: `~/tf_runs/prfba_rescore/bin/` (`prfba-ns`, `tf_env.sh`, `tf_unit`, the `eq_*.sh` equivalence
  scripts). Copies are in `tf/tf_side_bin/` here.

---------------------------------------------------------------------------------------------------

## 0. Read this first: TF already mounts FM's live `~/Documents`

TF has **FM's `/home/freiburger/Documents` mounted read-write over sshfs at
`/home/freiburger/Documents`** (`FM:/home/freiburger/Documents fuse.sshfs rw`), and read-only at
`/home/andrew/FM`. Every hard-coded `/home/freiburger/Documents/...` path in prFBA therefore resolves
on TF to **FM's live repo over the network**. A prFBA script run plainly on TF would read FM's files at
about 12 MB/s and **write its outputs into FM's working tree**.

So every prFBA command on TF runs through **`~/tf_runs/prfba_rescore/bin/prfba-ns`**. It is a private
user+mount namespace (`unshare --user --map-current-user --mount`, no root needed) in which
`/home/freiburger` is bind-mounted from the TF-local shadow tree
`~/tf_runs/prfba_rescore/home_freiburger`. Inside it:
- every hard-coded path resolves exactly as on FM, and so does `os.getcwd()`. A symlink tree would not
  give that: `getcwd()` returns the physical path, and `abspath()` of a relative path would then differ
  from FM's, which breaks, for example, `objective_wrapper_equivalence.py`'s weights-path check;
- FM's sshfs mount is hidden, so nothing can read or write FM by accident;
- the uid stays `andrew` (TF's Gurobi licence is a named-user licence), and capabilities are dropped
  before the command runs;
- nothing changes outside the namespace. No prFBA code was edited and nothing was committed.

**Two rules:**
1. Never run prFBA code on TF outside `prfba-ns` (or `tf_unit`, which wraps it).
2. **Never run a launcher that itself calls `systemd-run`** (`scripts/launch_wobj.sh`,
   `scripts/launch_objective_tuning.sh`), not even inside `prfba-ns`. The unit it creates starts in
   TF's normal namespace, where `/home/freiburger/Documents` is FM's live tree. Section 6 gives the
   equivalent `tf_unit` commands. The same applies to any FM driver that calls `systemd-run`.

The namespace refuses to start if the shadow tree is not in place. It checks that the marker
`Documents/.TF_SHADOW_TREE` exists and that `Documents` is on the shadow tree's device.

---------------------------------------------------------------------------------------------------

## 1. What is on TF, and where

Physical paths are under `/home/andrew/tf_runs/prfba_rescore/` (TFR). Inside `prfba-ns`, `TFR/home_freiburger`
*is* `/home/freiburger`.

| inside prfba-ns | physical (TFR/...) | content | provenance / check |
|---|---|---|---|
| `/home/freiburger/Documents/prFBA` | `home_freiburger/Documents/prFBA` | git clone of FM's prFBA, branch `hardware-aware-alignment-dispatch` at **cc8014fcee40** (= FM HEAD = tag `rescore-2026-11-start`), clean | from `git bundle` (220 MB); `git rev-parse HEAD` equal on both machines |
| `.../prFBA/<untracked inputs>` | same | re-score-invariant inputs, sections 1a-1c | `rsync --checksum`; sha256+size+mtime verified for all 441 files (section 7) |
| `/home/freiburger/Documents/ModelSEEDpy-c2e31a0` | `home_freiburger/Documents/ModelSEEDpy-c2e31a0` | clone at **c2e31a0c775e**, tree `9bb05aeb1f14`, clean | same commit and tree hash as FM's `~/Documents/ModelSEEDpy-c2e31a0` (clean, no ignored files other than `__pycache__`). TF's own `~/Documents/ModelSEEDpy-c2e31a0` is also at c2e31a0 and clean; the run clone is a separate copy so that nothing else on TF can move it |
| `/home/freiburger/Documents/ModelSEEDpy-prfba` | `home_freiburger/Documents/ModelSEEDpy-prfba` | clone at **09856f5f018a** (branch fix-prfba-single-model), tree `5e7358d4d0df`, clean, history contains c2e31a0 | same commit and tree as FM's worktree `~/Documents/ModelSEEDpy-prfba` |
| `/home/freiburger/Documents/ModelSEEDTemplates-c12run` | same under `home_freiburger/Documents` | git checkout at **21a808428560** | bundle of FM's ModelSEEDTemplates; clean |
| `/home/freiburger/Documents/ModelSEEDTemplates-main2294299` | " | git checkout at **22942996dacb** | not read by any TF stage (FM's R2.4/R2.5/R3 pin); here for completeness |
| `/home/freiburger/Documents/ModelSEEDDatabase-5389554a` | " | `git archive` of Biochemistry, Media, Objects, Libs, README.md, LICENSE at **5389554abfa8** | git tree hashes recomputed on TF equal FM's: Biochemistry `b0717c45a4f5`, Media `93b22cac2512`, Objects, Libs, README.md, LICENSE |
| `/home/freiburger/Documents/ModelSEEDDatabase-e8446e3e` | " | the same at **e8446e3ec205** | Biochemistry `318eda230997`, the rest identical |
| `/home/freiburger/Documents/ModelSEEDDatabase-41345b04` | " | the same at **41345b041f8c** (FM's current HEAD, branch tautomer-adenine-n1_20261113) | Biochemistry `c828fbb7dc58`, the rest identical |
| `/home/freiburger/Documents/ModelSEEDDatabase` | symlink, now **-> ModelSEEDDatabase-5389554a** | the path every TF stage reads by default | **decision needed, section 4** |
| `/home/freiburger/Documents/py_venv` | symlink -> `/home/andrew/Documents/py_venv` | TF's venv, same versions as FM's (section 3) | `sys.executable` reads `/home/freiburger/Documents/py_venv/bin/python`, as on FM |
| `/home/freiburger/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files/{bacteria,archaea}/ec.txt.gz` | `home_freiburger/opt/...` | the only PICRUSt2 files any TF stage reads (`calibrate_pfunc.P2_DEFAULTS`) | sha256 `cbd50a85...`, `09fd99e5...`, equal to FM's |
| `/home/freiburger/gurobi.lic` | symlink -> `/home/andrew/gurobi.lic` | TF's academic licence (named user andrew, host andrew-MS-7E70, expires 2027-11-04) | `GRB_LICENSE_FILE` points at it too |
| `/home/freiburger/Documents/prFBA_rescore_2026-11` | `home_freiburger/Documents/prFBA_rescore_2026-11` | TF's `$W`: `equivalence/`, `pre_repair/`, `logs/` | |
| (outside the namespace) | `control/bvbrc_cache_pre_repair/` | pre-repair shared PGFam cache (sha256 `f0a84bf198c3...`, 5,092,217,865 B, mtime 2026-11-11 11:30:22.7086 -0600) and pre-repair E3 catalogue (sha256 `46488fa0df4a...`) | rebuilt from `CORRECTIONS_2026-11-13.json` backup_old_values. Same bytes as FM's `$W/control/genome_gene_families.pre_repair.json`. Used for the R2.1 C-E3 control |
| (outside) | `bin/`, `xfer/`, `logs/` | wrappers, bundles, file lists | |

Two empty directories were created in TF's clone:
- `bvbrc_cache/`: the bind target for the control cache; the real cache goes here at R1.6;
- `mag_validation_rnd/out_batches/`: the bind target used by the MAG downstream equivalence.

The equivalence runs left the clone's tracked files unchanged: `git status` shows only these
directories, the pushed untracked inputs, and the smoke/`$W` outputs under `prFBA_rescore_2026-11`.

Symlink tree outside `TFR`: none. The only thing outside `~/tf_runs/prfba_rescore` that this setup
uses is TF's own venv `~/Documents/py_venv` (read only). `cobrakbase` and `MSCommunity` resolve through
that venv's `fm_repos.pth` to `/home/andrew/Documents/{cobrakbase,MSCommunity}`.

### 1a. E3 gate inputs (pushed)
`region_validation_d100/data/selection.json` (3.31 GB) and
`region_validation_d100/data/hits/asv_top20_alignment_hits.json` (4.00 GB). Both are union-tree
outputs, unchanged by the re-score. Also `nodes.dmp`, `names.dmp` and `merged.dmp`.

### 1b. MAG inputs (pushed)
- `mag_validation/{hits,manifests,hits_n500,manifests_n500,lowc_hits,lowc_manifests}/`,
  `mag_validation_km/{hits,manifests}/` and `mag_validation_rnd/{hits,manifests}/`: 3.0 GB, unchanged
  by the re-score.
- The three **current (pre-repair) panel indices** at their canonical paths:
  `mag_validation/pgfam_panel.json`, `kmeans_panel/pgfam_panel_kmeans.json` and
  `random_panel/pgfam_panel_random.json`. They are the W1 `--old-index` inputs, and they are stale
  (section 6, R2.2).
- FM's git-ignored old batch tables, **not** at their canonical paths. They were put straight into
  TF's `$W/pre_repair/mag_validation/{out_batches,lowc_out_batches}` and
  `$W/pre_repair/mag_validation_rnd/out_batches`, so the drivers' "skip existing batch" cannot reuse
  them. They are the W8 attribution input.
- Tracked old outputs **are** at their canonical paths, because they came with the clone:
  `mag_validation/out_batches_tau99`, `out_batches_n500`, `smoke/out` and `mag_validation_km/out_batches`.
  They must be deleted before re-scoring, as on FM (plan R2.2 step 3).

### 1c. Post-hoc (R4.3) inputs (pushed)
- Invariant: `region_validation/shared_reference/data/panel_ec_universe.json`, and, for all 20
  `data/prfba/fNN_dXXX/`, `hits/` (26.7 GB), `align_stamp.json`, `queries.fna`, `ref_db.json`,
  `ref_id_db.json` and `midas_empty.csv`.
- **Pre-re-score versions, pushed only for the equivalence runs. They must be replaced before R4.3**
  (list `tf/lists/pushed_2026-11-15_D_sr_PRE_RESCORE_replace_before_R4.3.txt`):
  - `cells.csv.gz(.fp)`, `_cache_panel_ec.json` and `calibration/scores.json`;
  - every `prfba/*/{predictions.json,_cache_ec_subset.json,_cache_pgfam_subset.json}`;
  - P2C `function_maps`, `pgfam_products`, `prfba_selection`, `genome_function_sets` and
    `comparison_cells`;
  - `weighted_objective/{population.json,per_amplicon.jsonl}` and `objective_tuning/per_amplicon.jsonl`.

  Pushing the regenerated versions over them with `tf_push.sh` replaces them. rsync only sends the
  changed blocks.

---------------------------------------------------------------------------------------------------

## 2. Reproducibility environment (`bin/tf_env.sh`, sourced by `prfba-ns`)

| variable | value | why |
|---|---|---|
| `HOME` | `/home/freiburger` | `score_mag.py`, `aggregate.py`, `compare_*`, `crossover.py` and the MAG drivers use `~/Documents/prFBA` and `~/Documents/py_venv` |
| `GRB_LICENSE_FILE` | `/home/andrew/gurobi.lic` | the licence is found independently of `HOME` |
| `OPENBLAS_CORETYPE` | `Haswell` | FM's OpenBLAS (numpy 0.3.31, scipy 0.3.30) picks Haswell kernels; TF's would pick SkylakeX (AVX-512) |
| `NPY_DISABLE_CPU_FEATURES` | `X86_V4 AVX512_ICL` | FM's numpy dispatches X86_V3 only; TF's adds AVX-512 loops. With this, TF's `numpy.show_config()` SIMD block equals FM's exactly |
| `TZ` | `America/Chicago` | both machines already use it; the E3 catalogue currency check formats the cache mtime in local time |
| `PATH` | the venv first | |
| `PYTHONPATH` | unset | each pinned run sets it on its own command line |

`tf_unit` adds the plan's `U()` environment: `PYTHONHASHSEED=0` and `SOURCE_DATE_EPOCH=1794700800`.

The CPU pins are kept as belt and braces. The E3 smoke was **also identical without them**
(`PRFBA_NS_NO_CPU_PINS=1`, a diagnostic-only switch). Gurobi has no such switch; the Stage O, wobj and
wrapper runs below show its objectives are bitwise equal on both CPUs anyway.

---------------------------------------------------------------------------------------------------

## 3. Versions on both machines (`tf/versions_FM.json`, `tf/versions_TF.json`)

Identical:
- Python 3.12.3, glibc 2.39, kernel 7.0.0-42.
- numpy 2.4.6 (scipy-openblas), scipy 1.17.1, pandas 3.0.3, cobra 0.29.1, optlang 1.9.1.
- gurobipy 13.0.1 (Gurobi 13.0.1), highspy, swiglpk, ijson 3.5.0 (yajl2_c), taxopy 0.14.0, orjson 3.11.7.
- matplotlib 3.10.3, statsmodels 0.14.6, sklearn 1.7.0, pyarrow 23.0.1, numba 0.65.1, Bio 1.85,
  torch 2.11.0+cu130.
- cobrakbase 5e4eedb (clean), MSCommunity 5616699 (clean on TF).
- The default solver is optlang.gurobi_interface on both.

Differences:
- **The CPU**, and with it the SIMD dispatch, which section 2 pins.
- The core count (64 vs 24 threads).
- The Gurobi licence file (each machine's own academic licence; FM's expires 2027-04-14).
- The paths of the venv's editable repos. On FM, `~/Documents/ModelSEEDpy` etc. are found through
  `__editable__` finders; on TF, through `fm_repos.pth` -> `/home/andrew/Documents/...`. Every TF stage
  that imports modelseedpy pins it with `PYTHONPATH`, and `--expect-modelseedpy` / `check_library`
  enforce the pin. The wrapper and Stage O runs below show the pinned paths resolve exactly to
  `/home/freiburger/Documents/ModelSEEDpy-*`.
- MSCommunity is imported, but it is not on the single-model path. Stage O ran with 03bd80c and the
  wrapper with 0532529-dirty; TF's 5616699 reproduces both bitwise.

Two departures from the preparation brief:
- **MSCommunity was not synced at the stages' commits.** TF's venv copy, 5616699, is used. It is the
  same commit as FM's current checkout, and it is not on the computational path, as the bitwise results
  show. TF's own repo contains 03bd80c if a run must record it.
- **`ModelSEEDpy-d884a50` was not synced.** It is the R3/R4.2 `score_folds` pin, and those stages stay on
  FM. TF's `~/Documents/ModelSEEDpy-d884a50` exists but was not checked.

---------------------------------------------------------------------------------------------------

## 4. DECISION NEEDED: which ModelSEEDDatabase the R4.3 model-building stages read

The wobj, Stage O, FVA and wrapper stages call
`phenotype_fold.worker_init(args.msd, ...)` -> `phenotype_arm.preload_biochemistry(msd)`. The `--msd`
default is `/home/freiburger/Documents/ModelSEEDDatabase`. The media
(`phenotype_arm.KB_MEDIA`, `KB_COMPLETE_CPD`) are always read from that path, but `Media/` is identical
at all three commits. None of the launchers passes `--msd`, and none of the outputs records the
database **commit**; only the path is recorded. FM's tree has moved under these studies:

| when | FM's `~/Documents/ModelSEEDDatabase` HEAD (reflog) | used by |
|---|---|---|
| 2026-10-14 to 10-30 16:42 | `e8446e3e` | **wobj** (run 10-29 20:40-22:25), the A-C12 scope (`C12_SCOPE.md`) |
| 2026-10-30 16:45 to 11-13 00:25 | `5389554a` (cpd11620 ferredoxin charge 4 -> 5; Fd doubled in 43 reactions; 95 files) | **Stage O** (11-03), **FVA** (11-03), **wrapper equivalence** (11-11) |
| since 2026-11-13 12:19 | `41345b04` (+ untracked `._.DS_Store`): Marvin/RDKit protonation and tautomer changes, 56 Biochemistry files, ~220k lines | nothing committed yet |

A re-run with the default path on FM today reads `41345b04`. That would move model-building numbers
for reasons unrelated to the cache, as D3 found for the function maps. It also explains the committed
Stage O `prior_match_fail` (2 keys) and the wrapper's `BAD_KEYS`: those are the rows whose LP changed
between `e8446e3e` (wobj) and `5389554a` (Stage O).

**Recommendation (D11): historical pins, so the re-score is a data-only delta.** Run wobj on `e8446e3e`,
and Stage O, FVA and the wrapper on `5389554a`. Either:
- switch TF's symlink between the stages:
  `ln -sfn ModelSEEDDatabase-e8446e3e ~/tf_runs/prfba_rescore/home_freiburger/Documents/ModelSEEDDatabase`;
  or
- pass `--msd /home/freiburger/Documents/ModelSEEDDatabase-<commit>` explicitly. `args.msd` then
  records the commit-named path, and the wrapper re-reads it from Stage O's `meta.json`.

FM's R4.1 phenotype passes face the same question (SPEC §9 reads the default path; `C12_SCOPE.md` names
`e8446e3e`). Whatever is decided, FM and TF must use the same commit for the same stage; all three
trees are on TF.

Another gap for R4.3: **`stage_o_weights_E3.json` has no producer in the repo.** It was derived by hand
on 2026-11-03 08:50 from `stage_o_weights.json`. The differences are `meta.estimator` "E4" -> "E3",
`meta.isotonic` -> null, and added `meta.derived_from {file, sha256, change}` and
`meta.generated_e3`; the 360 amplicon entries are identical. Stage E writes only `stage_o_weights.json`
and `_E1.json`. After Stage E is re-run, `launch_objective_tuning.sh` would stop with "missing" unless
this derivation is scripted (a small W-item) and re-run.

**Dependency correction to plan R4.3** ("calibrate ‖ Stage E ‖ wobj"): `pfunc_estimators.py` reads
`weighted_objective/population.json` (hard, L1452) and `calibration/scores.json` (optional, L1513).
The order is therefore:
1. `calibrate` ‖ wobj `--population-only`;
2. then Stage E ‖ the wobj main run;
3. then the E3 weights derivation;
4. then Stage O, FVA, results, wrapper and verification.

---------------------------------------------------------------------------------------------------

## 5. Moving data

### 5a. Code: FM commit -> TF
Commits only travel as git objects. **Run this after every FM commit that a TF stage needs**: the
R1.0 commits (W1 `build_pgfam_panel_index.py`, the uploader), any edit of `af_gate_counterfactual.py`'s
16S levels (R2.3 step 4), and the E3-weights derivation.

```bash
~/Documents/prFBA_rescore_2026-11/tf/bin/tf_code_sync.sh
#  bundles FM's commits that TF lacks (TF_HEAD..branch), fetches and fast-forwards on TF, and checks that
#  both HEADs are equal. It never writes FM's tree.
```

`tf_code_sync.sh` has not been exercised yet: at the time of writing, FM had no commit after cc8014f.
Its first use is the real test, and it stops with a message if the HEADs do not end up equal.
`tf_push.sh` (a no-op re-push, then VERIFY OK), `tf_pull.sh` (files and whole directories, checked
byte-identical) and `tf_verify.sh` were tested on 2026-11-15.

If TF's clone has local edits to tracked files (pushed regenerated outputs), the fast-forward can
refuse. In that case run `git -C <TF clone> stash` or `checkout -- <paths>` first: the commit carries
the same content. Uncommitted FM edits do not travel this way; push them as files and note them.

### 5b. Inputs: FM -> TF (checksummed)
```bash
T=~/Documents/prFBA_rescore_2026-11/tf
$T/bin/tf_push.sh $T/lists/<list>.txt   # rsync -a -R --checksum -z(zstd-3) --partial, then tf_verify.sh
$T/bin/tf_verify.sh $T/lists/<list>.txt # sha256 + size + mtime on both sides; exit 1 on any difference
```

`mtime` matters: the E3 gate refuses to run unless the cache's size and mtime (to the second) match
the catalogue's `_meta`. rsync `-a` preserves it, and `tf_verify.sh` checks it.

**Push-later manifest.** These inputs are regenerated by FM stages, so they are not on TF yet, except as
the pre-re-score versions of section 1c. Sizes are FM's current files. Wire estimates use measured
zstd-3 ratios and 11 MB/s; the measured link is 10.9-12.0 MB/s.

| list file (`tf/lists/`) | paths | size now | on the wire | min | regenerated on FM by | needed on TF by |
|---|---|---|---|---|---|---|
| `push_R2.1_e3gate.txt` | `bvbrc_cache/genome_gene_families.json` | 5.12 GB | 0.70 GB (x7.3) | 1.1 | R1.6 freeze | R2.1, W1 if built on TF |
| " | `bvbrc_cache/e3_clade_catalogue.json` | 10.7 MB | 1.5 MB | 0 | R1.6 (only if the cache moved) | R2.1 |
| " | `bvbrc_cache/CORRECTIONS_2026-11-13.json` | 16.8 MB | 2.4 MB | 0 | static | W1 `--restore-old-values` if on TF |
| " | `region_validation_d100/data/gene_capture_cells{,_tau0995}.csv` (tracked) | 2 x 20.6 MB | 8 MB | 0 | R1.3 P7, R1.5 E1 | R2.1 reference cells |
| `push_R2.2_mag_indices_if_built_on_FM.txt` | the 3 W1 panel indices | 0.40 GB | 0.10 GB (x4.0) | 0.2 | R2.2 W1 | R2.2 |
| `push_R2.3_mag_downstream.txt` | `region_validation{_d100,_d100_legacy,}/data/gene_capture_cells.csv` + `gene_capture_by_region.json` (tracked) | 62 MB | 13 MB | 0 | R1.3, R1.4 | R2.3 MAG pairing, af_gate levels |
| `push_R4.3_posthoc.txt` | `SR/data/cells.csv.gz(.fp)` | 56 MB | 56 MB | 0.1 | R3.5 pass 1; final after R4.2 (below) | all of R4.3 |
| " | `SR/data/_cache_panel_ec.json` | 49 MB | 6.6 MB | 0 | R3.1 W3 | all of R4.3 |
| " | `SR/data/prfba/*/predictions.json` (20) | 9.23 GB | 0.65 GB (x14.2) | 1.0 | R3.4 | calibrate, Stage E, wobj population |
| " | `SR/data/prfba/*/_cache_ec_subset.json` (20) | 0.46 GB | 0.06 GB | 0.1 | R3.4 | Stage E, Stage O, wrapper, calibrate |
| " | `SR/data/prfba/*/_cache_pgfam_subset.json` (20) | 2.43 GB | 0.47 GB (x5.2) | 0.7 | R3.4 | calibrate |
| " | `P2C/data/{function_maps,pgfam_products,prfba_selection,genome_function_sets}.json`, `comparison_cells.csv.gz` | 1.98 GB | 0.20 GB | 0.3 | R2.5 | calibrate |
| **total** | | **19.8 GB** | **≈ 2.3 GB** | **≈ 3-4** | | |

Generated by `tf/make_manifest.py` (`tf/push_later_manifest.json`). If the W8/attribution step wants the
pre-repair MAG index on TF too, it is already there (section 1b).

`cells.csv.gz` is rewritten by R4.2 (pass 2), and the plan starts R4.3 after R3.5, in parallel with R4.1.
If R4.3 is to read the pass-1 file, push it once after R3.5. Pass 2's arms must equal pass 1's, so the
cells used are the same. Record the sha256 that TF read: `tf_verify.sh` prints it.

### 5c. Outputs: TF -> FM
```bash
$T/bin/tf_pull.sh $T/lists/pull_<stage>.txt      # -> ~/Documents/prFBA_rescore_2026-11/from_tf/pull_<stage>/
# inspect / diff there, then install:
rsync -a --checksum ~/Documents/prFBA_rescore_2026-11/from_tf/pull_<stage>/ ~/Documents/prFBA/
```

The pull lists: `pull_R2.1_e3gate.txt`, `pull_R2.2_mag.txt`,
`pull_R2.2_mag_indices_if_built_on_TF.txt`, `pull_R2.3_mag_downstream.txt` and `pull_R4.3_posthoc.txt`.

TF's `$W` outputs, such as the E3 review dir, are pulled directly:
```bash
rsync -a --checksum TF:/home/andrew/tf_runs/prfba_rescore/home_freiburger/Documents/prFBA_rescore_2026-11/e3_review/ \
      ~/Documents/prFBA_rescore_2026-11/e3_review/
```

Pull R2.2 **before** FM runs `run_downstream_d100.sh` (R2.3 step 1): that script's MAG block reads
`mag_validation/out/mag_cells_all.csv` and `lowc_out_batches`.

---------------------------------------------------------------------------------------------------

## 6. Launching a stage on TF

```bash
ssh TF
TFR=~/tf_runs/prfba_rescore; U=$TFR/bin/tf_unit
R=/home/freiburger/Documents/prFBA; SR=$R/region_validation/shared_reference
W=/home/freiburger/Documents/prFBA_rescore_2026-11; L=$W/logs       # paths as seen INSIDE prfba-ns
# tf_unit UNIT MEMORYMAX LOG 'COMMAND': transient user unit, MemoryMax, MemorySwapMax=0, PYTHONHASHSEED=0,
# SOURCE_DATE_EPOCH=1794700800, COMMAND run in prfba-ns with cwd $R; "EXIT=<rc>" appended to LOG.
systemctl --user status <unit>; journalctl --user -u <unit> -o cat | tail   # the memory peak is printed on stop
```

**Interactive or short commands:** `$TFR/bin/prfba-ns <cmd>` (cwd `$R`; set `PRFBA_NS_CWD` to change it).

**Linger.** `andrew` has `Linger=no`. User units survive while any andrew session exists; today that is
the seat0 desktop session. If the desktop logs out with no ssh session open, systemd kills running
units. Either keep the desktop logged in, or ask the user to run `loginctl enable-linger andrew` on TF
once. That changes TF's configuration, so it is their call.

**Memory budget.** TF shows `MemAvailable` ≈ 69 GB with the desktop and its apps running (about
20 GB used). Keep ≈ 8 GB for the desktop and run at most **about 60 GB of `MemoryMax` at once**.

**Pre-flight before every stage:**
- `tf_code_sync.sh` reports TF HEAD == FM HEAD;
- `tf_verify.sh` passes on that stage's push list(s);
- the stage's old outputs on TF have been moved to `$W/pre_repair/` (TF's `$W`), as on FM;
- `readlink /home/freiburger/Documents/ModelSEEDDatabase` inside `prfba-ns` is the agreed pin (R4.3).

### R2.1 E3 gate + review (after R1.6 and `push_R2.1_e3gate.txt`)
```bash
# (optional) C-E3 control: the 2026-11-12 inputs = pre-repair cache + old catalogue + the COMMITTED (cc8014f)
# reference cells. The regenerated cells arrive with push_R2.1 or tf_code_sync, so bind the committed copies
# (TFR/control/ref_cells_cc8014f/, sha256 db5f6e1c... and 5a84dd24..., = cc8014f) over the canonical paths.
# The control may then run before or after those pushes. It writes the canonical outputs in TF's clone, so run
# it BEFORE moving those outputs aside, and compare / restore after:
CC=$TFR/control/ref_cells_cc8014f; D=$R/region_validation_d100/data
PRFBA_NS_BINDS="$TFR/control/bvbrc_cache_pre_repair=$R/bvbrc_cache;$CC/gene_capture_cells.csv=$D/gene_capture_cells.csv;$CC/gene_capture_cells_tau0995.csv=$D/gene_capture_cells_tau0995.csv" \
  $U e3-control 32G $L/e3_control.log "python -u region_validation_d100/scripts/pfunc_e3_panel_gate.py --workers 10"
#   (file binds need the target files to exist: they do once the clone is checked out.) Then compare the 'results'
#   blocks with git show cc8014f:region_validation_d100/data/pfunc_e3_gate.json, copy the control outputs to $W/control_e3/,
#   and git checkout -- region_validation_d100/data/pfunc_e3_gate.json region_validation_d100/data/pfunc_e3_gate/
# move the old outputs aside (TF clone), then the smoke twice (determinism), then the gate:
mkdir -p $TFR/home_freiburger/Documents/prFBA_rescore_2026-11/pre_repair/region_validation_d100/data
#   mv .../prFBA/region_validation_d100/data/pfunc_e3_gate{.json,/per_amplicon.csv.gz,/run_meta.json,/review_recompute.json} -> pre_repair
$U e3-smoke1 20G $L/e3_smoke1.log "python -u region_validation_d100/scripts/pfunc_e3_panel_gate.py --limit 500 --workers 8 --verify-q && mv region_validation_d100/data/pfunc_e3_gate/smoke $W/e3_smoke1"
$U e3-smoke2 20G $L/e3_smoke2.log "python -u region_validation_d100/scripts/pfunc_e3_panel_gate.py --limit 500 --workers 8 --verify-q && mv region_validation_d100/data/pfunc_e3_gate/smoke $W/e3_smoke2"
$U e3-gate-rescore 32G $L/e3_gate.log "python -u region_validation_d100/scripts/pfunc_e3_panel_gate.py --workers 10"
$U e3-review 32G $L/e3_review.log "PFUNC_GATE_REVIEW_OUT=$W/e3_review WORKERS=8 python -u region_validation_d100/scripts/pfunc_e3_gate_review_recompute.py && cp $W/e3_review/indep_recompute.json region_validation_d100/data/pfunc_e3_gate/review_recompute.json"
```
The checks are the plan's: `reproduction_vs_gene_capture_cells` reads `match` for all 135,982 rows,
and `catalogue_is_current` is true. Then pull `pull_R2.1_e3gate.txt` and `e3_review/`.

### R2.2 MAG re-score (after `tf_code_sync.sh` brings W1)
```bash
# 1. W1 indices, here (needs push_R2.1's cache + CORRECTIONS json) or on FM (then push_R2.2_...):
#    the plan's R2.2 step 1 commands, run via prfba-ns; keep the old indices in $W/pre_repair.
# 2. Control (md5 a57df30b... / 3cb0e487...): the plan's step 2 with --pgfams $W/control/pgfam_panel.pre_repair.json.
# 3. Clear the tracked old batches in TF's clone (the git-ignored ones are already only in $W/pre_repair):
#    rm -rf $R/mag_validation/{out_batches_tau99,out_batches_n500,smoke/out} $R/mag_validation_km/out_batches
#    The rnd dir mag_validation_rnd/out_batches exists, empty (an equivalence bind target); leave it.
# 4. Re-score: the drivers run -P 6 (main, tau, lowc) and -P 5 (n500, km, rnd); ~4.3-4.8 GB per scorer.
$U mag-rescore 40G $L/mag_rescore.log "bash mag_validation/score_all.sh; bash mag_validation/score_tau.sh 99; bash mag_validation/score_n500.sh; bash mag_validation/score_lowc.sh; bash mag_validation_km/score_km.sh; bash mag_validation_rnd/score_rnd.sh; python mag_validation/score_mag.py --hits mag_validation/smoke/hits.tsv --manifest mag_validation/smoke/manifest.json --outdir mag_validation/smoke/out"
# 5. Concatenate and aggregate: the plan's R2.2 step 5, via prfba-ns. The concatenation recipe is in
#    bin/eq_mag_agg.sh (cat_cells).
```
Don't run the MAG re-score at the same time as the E3 gate: about 29 + 21 GB plus page cache.

### R2.3 MAG parts (after `push_R2.3_mag_downstream.txt`; R2.2 done)
Run the plan's R2.3 steps 2-5 via `prfba-ns`, with `PYTHONHASHSEED=0` (tf_unit sets it):
- `compare_three_panels_magpanel.py`;
- `compare_mag_panels.py`;
- the new `paired_tau99_d100.json`;
- `af_gate_counterfactual.py`, after `tf_code_sync.sh` if its 16S levels were edited on FM;
- the d20 pairings, crossover and plots.

Then pull `pull_R2.3_mag_downstream.txt`.

### R4.3 post-hoc (after the push lists, the D11 decision and the E3-weights producer)
Each command is equivalent to its FM launcher with the same python command lines, pins and memory caps.
The launchers' `MemAvailable` checks become: check `free -g` first.

Before starting, move to TF's `$W/pre_repair/`:
- `SR/data/{calibration,pfunc_estimators,weighted_objective,objective_tuning,objective_tuning_fva}`
  (keep the `verification/*.py` scripts);
- the pre-re-score `cells.csv.gz` etc., which the pushes overwrite anyway.

```bash
export PRFBA_NS_CWD=$SR
PIN="PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-c2e31a0"
# (1) calibrate || wobj population
$U calib-pfunc 20G $SR/logs/tf_calibrate.log "python $SR/scripts/calibrate_pfunc.py --forced --workers 4 --cache-dir $W/calib_cache"
#     (= the committed run's args: scores.json meta.args; --modes/--densities/--spaces/--bootstrap/--seed/--eps are the defaults)
$U wobj-pop 12G $SR/logs/wobj_run.log "$PIN python -W ignore $SR/scripts/weighted_objective_test.py --outdir $SR/data/weighted_objective --expect-modelseedpy c2e31a0 [--msd <D11 pin>] --population-only"
# (2) Stage E || wobj main run
$U pfe-run 20G $SR/logs/tf_pfunc_estimators.log "python $SR/scripts/pfunc_estimators.py --workers 12 --bootstrap 2000 --seed 1729 --densities 1.0,0.5,0.25 --outdir data/pfunc_estimators --cache-dir $W/pfe_cache"
$U wobj-run 20G $SR/logs/wobj_run.main.log "$PIN python -W ignore $SR/scripts/weighted_objective_test.py --outdir $SR/data/weighted_objective --expect-modelseedpy c2e31a0 [--msd <D11 pin>] --workers 8 --maxtasksperchild 2"
#     then: python scripts/wobj_results.py      (prfba-ns, cwd $SR)
# (3) E3 weights from the new stage_o_weights.json (the W-item in section 4)
# (4) Stage O: validate, run, results, FVA, wrapper, verification
$TFR/bin/prfba-ns bash -c "cd $SR && $PIN python -W ignore scripts/objective_tuning.py --outdir $SR/data/objective_tuning --weights $SR/data/pfunc_estimators/stage_o_weights_E3.json --expect-modelseedpy c2e31a0 [--msd <D11 pin>] --validate-only"
$U otune-run 34G $SR/logs/otune_run.log "$PIN python -W ignore $SR/scripts/objective_tuning.py --outdir $SR/data/objective_tuning --weights $SR/data/pfunc_estimators/stage_o_weights_E3.json --expect-modelseedpy c2e31a0 [--msd <D11 pin>] --workers 10 --maxtasksperchild 1"
$TFR/bin/prfba-ns bash -c "cd $SR && python scripts/objective_tuning_results.py"
$U otune-fva 34G $SR/logs/otune_fva_run.log "$PIN python -W ignore $SR/scripts/objective_tuning.py --outdir data/objective_tuning_fva --weights data/pfunc_estimators/stage_o_weights_E3.json --half confirm --arms two_tier/lin/G0/X100,kn/lin/G0/X100,est/lin/G0/X1,est_shuffled/lin/G0/X1,two_tier/lin/G0/X1 --fva-arms est/lin/G0/X1,kn/lin/G0/X100,two_tier/lin/G0/X1,est_shuffled/lin/G0/X1 --expect-modelseedpy c2e31a0 [--msd <D11 pin>] --workers 10 --maxtasksperchild 1"
#     (the FVA args are those in data/objective_tuning_fva/meta.json, b52984f; relative paths as recorded)
$TFR/bin/prfba-ns bash -c "cd $SR && PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-prfba python -W ignore scripts/objective_wrapper_equivalence.py"
#     verification: python data/pfunc_estimators/verification/p*.py and data/objective_tuning/verification/{extract,c1..c13}.py
#     in prfba-ns, cwd $SR, in the order of their VERIFICATION.md files (extract.py first).
```
Then pull `pull_R4.3_posthoc.txt`.

---------------------------------------------------------------------------------------------------

## 7. Equivalence evidence (2026-11-15, all on TF through prfba-ns, CURRENT committed inputs)

Every run below read the inputs that produced FM's committed results, pre-repair where the committed
run was pre-repair. TF reproduced FM in every computed value. The only differences are timings,
timestamps, paths deliberately redirected to `$W/equivalence`, the git commit, and fields that FM's own
re-run also changes. Each of those is explained in the table. The outputs are in TF's
`$W/equivalence/` (physical `TFR/home_freiburger/Documents/prFBA_rescore_2026-11/equivalence/`); copies
of the comparison files are in `tf/equivalence/`.

| # | what (pins) | compared with | result |
|---|---|---|---|
| (a) | `objective_wrapper_equivalence.py`: 8 amplicons rebuilt (gapfill + 2 harness + 3 wrapper solves each); ModelSEEDpy-prfba 09856f5, templates 21a8084, MSD 5389554a | committed `wrapper_equivalence.json` (FM, 2026-11-11) | `all_pass`; `rebuild_matches_stage_o_all`; `harness_vs_stage_o_bitwise_all` true (tuned and released); every objective, flux and primal bitwise equal to FM's file. 1 field differs: `wrapper.sha256` of `prfba_objective.py`, which changed in 6e2c2eb (11-12) after FM's run. FM HEAD has the same new sha. So the shipped wrapper also reproduces bitwise |
| (a') | Stage O subset: `objective_tuning.py --only` 4 amplicons (n2, n3, n4+, and the slowest n3), 90 arms each; ModelSEEDpy c2e31a0, templates 21a8084, MSD 5389554a, E3 weights | committed `objective_tuning/per_amplicon.jsonl` rows | **360/360 arm objectives bitwise equal; 0 differing fields** in all 4 rows (every arm metric, `model`, `gapfill`, `prior_match` = True, `released_check` = True, `determinism`); timing and memory excluded |
| (a'') | Stage O **FVA variant** subset: `objective_tuning.py` with the FVA run's args (`--half confirm`, 5 arms, 4 FVA arms), 3 confirm-half amplicons (n1, n2, n4+); MSD 5389554a | committed `objective_tuning_fva/per_amplicon.jsonl.gz` rows | **all identical**: 15/15 arm objectives bitwise equal, 0 differing fields, FVA faces included |
| (a''') | wobj subset: `weighted_objective_test.py` rebuilds the population from FM's current inputs, then 3 amplicons (n2, n3, n4+) x 8 arms; MSD e8446e3e | committed `weighted_objective/{population.json,per_amplicon.jsonl}`, and an FM re-run of the same 3 amplicons today | population: `jobs`, `spec`, `assertions`, candidate counts **identical** (only `generated`, `git_commit` differ). Rows: every arm's primary solve, metrics and optimal-face (FVA) values **bitwise identical**, and the same holds for the primal and dual simplex re-solves (`vertex_resolves` "0", "1"). **The only difference is the Gurobi barrier re-solve** (`vertex_resolves` "2", a degeneracy diagnostic): 69 fields, max relative difference **1.4e-14**, plus the spreads derived from it (`*_vertex_spread`, `objective_vertex_rel_diff_max`). The FM re-run reproduces the committed rows bitwise, so the barrier is deterministic on each machine but CPU-dependent (Zen 5 vs Zen 2). See the wobj recommendation in section 8 |
| (b) | E3 gate smoke `--limit 500 --workers 3 --verify-q` on the pre-repair cache (sha `f0a84bf1`) + old catalogue (`46488fa0`) + committed reference cells | FM's smoke of 2026-11-12 (same args) and the committed full `per_amplicon.csv.gz` | 1,000 (key, mode) rows: **0 differing rows** vs FM's smoke (all 57 columns except `t_select_ms`, `t_e3_ms`). **0 differing rows** vs the committed full run (55 shared columns). `results`, `verdict`, `verdict_detail` and `per_region_n2+` (bootstrap CIs included) are **identical**. `_meta.store` (159,753 genomes, 318,108,350 ids), `reproduction_vs_gene_capture_cells` {match: 1000} and `catalogue.is_current` true. Only `diagnostics.*median_t_e3_ms*` differ (timing). Also identical **without the CPU pins and with 4 workers** |
| (c) | MAG `score_mag.py` with the CURRENT (pre-repair) panel indices, batch b000 of 4 arms: main, τ99, k-means, random | FM's on-disk `out_batches/b000`, tracked `out_batches_tau99/b000`, tracked `mag_validation_km/out_batches/b000`, on-disk `mag_validation_rnd/out_batches/b000` | **`mag_cells.csv` and `mag_summary.json` byte-identical in all 4** (cells md5 c93e2432 / 1abb9600 / c3e2c85e / 896d0e11) |
| (c') | MAG downstream from FM's batch tables: concatenation, `aggregate.py` x4, `crossover.py` x3, `compare_three_panels_magpanel.py` | committed `out/mag_cells_all.csv` (md5 **a57df30b**), `mag_cells_tau99.csv` (**3cb0e487**), `mag_summary_full/tau99`, `summary_n500`, `summary_km`, `crossover{,_d100,_d100_legacy}`, `three_panel_mag_panelarm.json` | **all byte-identical** |
| (c'') | the hash-seed-dependent ones, `PYTHONHASHSEED=0` on both machines: `compare_arms.py` x4, `af_gate_counterfactual.py`, `compare_mag_panels.py` | FM re-run of the same commands into scratch | **6/6 byte-identical FM vs TF**. Against the committed files (random hash seed at the time), `paired_*` differ only in the order of keys of equal-completeness tiers. The values are identical (dict equality). `compare_mag_panels` agrees at 4 dp (summation order; the plan expects this). `af_gate_counterfactual.json` has an identical summary and all 9,997 per-MAG records; only its `_about` text was edited after generation |
| (d) | Stage E `pfunc_estimators.py --workers 12 --bootstrap 2000 --seed 1729` on FM's current inputs | committed `data/pfunc_estimators/*` | `selected.json`, `assertions.json` **byte-identical**. `stage_o_weights.json` and `_E1.json` differ only in `meta.generated` (all 360 entries and the 747-knot isotonic map identical). `scores.json`: every number identical. The diffs are timings, `args.outdir/cache_dir` (redirected), `git_commit`, and the completion order of `meta.stream_jobs` (imap_unordered). README: the same apart from the header provenance, the runtime line and FM's hand-appended corrections section |
| (e) | `calibrate_pfunc.py` (committed args: `--forced --workers 4`, defaults otherwise) on FM's current inputs | committed `data/calibration/*` (`scores.json` untracked on FM) | `reliability.csv` (72,188 rows) **byte-identical**. `scores.json`: **every number identical**; the only diffs are in `meta` (generated, git_commit, args.outdir/cache_dir, layer-A job completion order and per-job seconds/RSS). `assertions.json` differs only in `secondary.seconds` / `peak_rss_gb`. README identical apart from its date/commit header |

**Pushed files.** `tf_verify.sh` over the four lists pushed on 2026-11-15 reports **441 files (45 GB)
identical: sha256, size and mtime**. The 120 MAG `pre_repair` copies are sha256-identical too. The MSD
trees are verified by git tree hashes, the pinned git checkouts by commit and tree hash, and the
control cache by sha256 (all in section 1).

**What this proves.** These code paths reproduce FM bit for bit on TF:
- E3 gate (incl. bootstrap CIs);
- MAG scoring and every MAG downstream script;
- calibrate and Stage E (numpy- and bootstrap-heavy);
- model building, gapfilling and every LP solve of Stage O, FVA, the wrapper and wobj.

There is one exception: wobj's Gurobi-barrier re-solve diagnostic (≤ 1.4e-14 relative). The
namespace's paths, the venv, the pins and the CPU environment are therefore sufficient. The only
residual risk is inputs, which `tf_verify.sh` checks before each stage.

Not exercised on TF:
- `pfunc_e3_gate_review_recompute.py`. It uses the same IntStore/E3 code as the gate (550 s on FM);
  run it once on the control cache if wanted.
- The W1 index builder. It is not written yet; its own control (md5 `a57df30b` / `3cb0e487`) covers it.
- The verification scripts and `objective_tuning_results.py` / `wobj_results.py`. These are
  single-process, solver-free post-processing of the rows proven identical above.

---------------------------------------------------------------------------------------------------

## 8. Throughput: TF vs FM on the same unit, and recommended settings

| unit | FM | TF | TF/FM speed |
|---|---|---|---|
| MAG `score_mag.py`, 1 batch, 1 process (4.3-4.8 GB) | 33.2 s (b000, today, load ≈ 9) | 21.0 s (b000); 20.9 τ99, 23.4 km, 24.3 rnd | 1.6x |
| E3 smoke, 500 amplicons, 3 workers (incl. streaming the 5.1 GB cache into the int32 store, 78 s on TF) | 229 s, cgroup peak 5.15 GB | 112 s (pinned), 104 s (unpinned, 4 workers), peak 5.12 GB | 2.0x; per-amplicon E3 median 35 vs 75 ms (2.1x) |
| wrapper equivalence, 8 amplicons, 1 process | 612 s (gapfill 42-60 s each) | 284 s (gapfill 17-27 s each) | 2.2x |
| Stage O amplicon (gapfill + 90 arms), 4 concurrent on TF | 125-291 s (10 concurrent) | 69-89 s (worst 89 vs 291) | 1.8x median |
| Stage E, full, 12 workers | 411 s; RSS parent 7.94, child 9.87 GB | 365 s; **cgroup peak 8.8 GB**; pass A slower (212 vs 155 s) | 1.1x; 12 busy workers saturate TF's 12 cores |
| calibrate, full, 4 workers | ≈ 735 s, 26.8 min CPU, 11.4 GB peak (10-20 journal) | 406 s, 15.3 min CPU | 1.8x |
| wobj amplicon (gapfill + 8 arms + FVA faces + 3 re-solves), 3 concurrent | 125 / 135 / 160 s (FM re-run today) | 73 / 85 / 99 s | 1.6-1.7x |
| Stage O FVA amplicon (5 arms, 4 with FVA), 3 concurrent | 117-138 s (10 concurrent) | 45-55 s | 2.5x |
| `aggregate.py` x4 | - | 12.6 s | - |

**The scheduling fact: with 12 busy workers TF is only about FM's speed.** Stage E: 365 vs 411 s, with
pass A slower on TF (212 vs 155 s). TF has 12 cores and dual-channel memory. Single-process work is
1.6-2.5x faster on TF. Do not plan for more than about 10-12 concurrent CPU-bound workers on TF, and
expect parallel stages to finish at FM's speed or slower.

The FM timings are biased in TF's favour. They come from committed runs under unknown load, or from
re-runs today at load average 7-9 (the KBase upload and two codif units were running), so the
single-process ratios are probably somewhat optimistic.

| stage on TF | workers | MemoryMax | expected wall on TF | basis |
|---|---|---|---|---|
| R2.1 gate | `--workers 10` | 32G | ≈ 15-20 min | FM 1,595 s at 12 workers, peak 21.2 GB; TF 2.1x per worker; cache load 80 s |
| R2.1 review recompute | `WORKERS=8` | 32G | ≈ 5-8 min | FM 550 s at 8 workers |
| R2.2 MAG re-score | drivers as committed (-P 6 / -P 5) | 40G | ≈ 10 min for all six drivers + smoke | 21-24 s per batch at P1; FM 153 s (main), 549 s (τ99), 104 s (n500), 170 s (km), 118 s (rnd) |
| R2.2/R2.3 aggregate, compare, crossover, af_gate | 1 | 8G | ≈ 1-2 min | measured |
| R4.3 calibrate | `--workers 4` | 20G | ≈ 8-12 min | FM ≈ 12 min, 11.4 GB peak (journal of `calib-fix-full`) |
| R4.3 Stage E | `--workers 12` | 20G | ≈ 6 min | measured 365 s, 8.8 GB |
| R4.3 wobj | `--workers 8 --maxtasksperchild 2` (launcher) | 20G | ≈ 60-70 min + 2 min population | FM 1.75 h; TF ≈ 1.65x per amplicon; **see the barrier note below** |
| R4.3 Stage O | `--workers 10 --maxtasksperchild 1` | 34G | ≈ 45-55 min | 360 x ≈ 75 s / 10 |
| R4.3 FVA (193 rows) | 10 | 34G | ≈ 30-35 min | FM 3,375 s |
| R4.3 wrapper | 1 | 8G | ≈ 5 min | measured 284 s |

Concurrency on TF, staying within about 60 GB:
- the E3 gate (32G) alone, or with the R4.3 single-process steps;
- the MAG re-score (40G) alone;
- R4.3: calibrate (20G) ‖ wobj population (12G); then Stage E (20G) ‖ wobj (20G); then Stage O (34G),
  then FVA (34G).

The R2 and R4 lanes are hours apart in the plan, so they do not compete.

**wobj and the barrier note.** On TF, wobj's per-row Gurobi barrier re-solve (`vertex_resolves["2"]`)
differs from FM's by up to 1.4e-14 relative. So do the three spreads derived from it, which
`summary.json` `checks.vertex_resolves` and `WEIGHTED_OBJECTIVE_RESULTS.md` §3 report (max relative
objective difference across methods 1.3e-13 / 2.0e-12; rows whose TAF depends on the vertex, thresholded
at 1e-6). Gurobi exposes no CPU-dispatch switch. Everything else in wobj is bitwise identical, and Stage O
reads nothing from these fields.
- **If the repository must be bit-identical to what FM would produce, run the wobj main run on FM**
  (≈ 1.75 h at 8 workers). The population step and Stage E/O/FVA/wrapper can stay on TF.
- Otherwise run it on TF and record that the barrier diagnostic is machine-dependent at the 1e-14 level.

Everything else assigned to TF is bitwise reproducible there.
