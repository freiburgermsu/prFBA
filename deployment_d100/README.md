# Depth-100 rerun of the two real-amplicon deployments

**Status (2026-11-12): run and uploaded.** Both depth-100 builds were made with k/n feature
probabilities (2026-10-29/30; RESULTS.md) and uploaded to KBase on 2026-11-04 by
`kbase_upload_d100.sh` (user-approved), 0 failures (`logs/kbase_upload_d100.log`):
codiffusion 3,225 genomes -> GenomeSet 261948/6467/1 (beside the 20-hit build in that workspace),
EmilyKin 3,755 genomes -> GenomeSet 275568/3757/1 in a new narrative. Since 2026-11-12 the
builders default to E3 for multi-genome selections (section 9). Since 2026-11-13 every build
uses BV-BRC (PATRIC) annotation only (section 10).
- The shipped builds and their uploads are k/n builds that also merged BV-BRC's imported RefSeq
  CDSs. They are reproducible only from git at or before `3ce12d0`.
- Re-running `d100_pipeline.sh` (default `PFUNC=e3`) makes a NEW BV-BRC-only build into a NEW
  `SYN`, and refuses to write into an existing one.

The sections below record the setup as it was planned on 2026-10-29.

The benchmark moved to 100-hit records (`align_hits.py` defaults `--k2 100 --prefilter-k 500
--gpu-cand 200`; CUDA path: edlib top-500 prefilter, then GPU Smith-Waterman ranks all 500,
then Biopython traceback of the top 100). The two real-amplicon deployments were made on 20-hit
records:

| | ASVs | 20-hit hits | 20-hit build selection | 20-hit synthetic genomes |
|---|---|---|---|---|
| A codiffusion bioreactor (V4-V5) | 3,276 | `codiffusion_bioreactor/bvbrc_alignment_hits/asv_top20_alignment_hits_gpu.json` | `asv_reference_selection_gpu.54d1150.json` (KBase build) | `synthetic_genomes_54d1150/` (uploaded to ws 261948) |
| B EmilyKin | 3,950 | `EmilyKin/bvbrc_alignment_hits/asv_top20_alignment_hits_gpu.json` | `asv_reference_selection_gpu.json` | `EmilyKin/synthetic_genomes/` |

None of these files are touched. The depth-100 outputs go to NEW directories only:

| | depth-100 hits, selection, stats | depth-100 synthetic genomes |
|---|---|---|
| A | `/home/freiburger/Documents/codiffusion_bioreactor/bvbrc_alignment_hits_d100/` | `/home/freiburger/Documents/codiffusion_bioreactor/synthetic_genomes_d100/` |
| B | `/home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits_d100/` | `/home/freiburger/Documents/EmilyKin/synthetic_genomes_d100/` |

Inside each `bvbrc_alignment_hits_d100/`: `asv_top20_alignment_hits.json` (the pipeline's fixed
file name and `top20` key; each record holds up to 100 hits), `method_selection.json`,
`edlib_top500_shortlist.npz` + `gpu_scores.npy` + `shortlist_meta.json` (resume cache),
`gene_families_tieset.json` (+ `.tieset.json`), `asv_reference_selection_d100.json`,
`deployment_stats_d100.json`.

## 1. Which 20-hit file and counting method give the published numbers

`deployment_stats.py` (this directory) is the counting method. Run on
`EmilyKin/bvbrc_alignment_hits/asv_reference_selection_gpu.json` with the default
`--taxon-source index`, it reproduces every published 20-hit deployment number exactly
(`stats_d20/emilykin_gpu.json`):

| quantity | published (prFBA.tex l.211, 223, 290, 1072, 1134; deck) | reproduced |
|---|---|---|
| ASVs / answered / abstaining | 3,950 / 3,753 / 197 | 3,950 / 3,753 / 197 (139 below family, 58 reliability) |
| n_genomes = 1 | 2,344 (62.5%) | 2,344 (62.5% of answered) |
| n_genomes = 2 / >= 3 | 599 / 810 | 599 / 810 |
| mean n_genomes | 2.034 | 2.0341, **over answered ASVs** (the selection `_meta.summary.mean_handful` 1.933 is over all 3,950) |
| cross-phylum, two-genome | 31.7% (143/451) | 143/451 = 31.7% |
| cross-phylum, larger | 59.0% (462/783) | 462/783 = 59.0% |
| resolved multi-genome merges | 1,234 of 1,409 | 1,234 of 1,409 |
| Bacteria + Archaea merges | 18 | 18 |
| (None-counted variant, superseded) | 40.6% / 66.2% (old text) | 183/451 = 40.6% / 518/783 = 66.2% |

How the number is counted (mirrors the audit's re-run script that traced the published
figures, `scratchpad/v3/rerun.py`, which ran over the build's `provenance/*.json`):
* `n_genomes` = `len(selected)`; equals the build provenance's `genomes` list for all 3,753
  answered ASVs (checked: 0 mismatches).
* cross-phylum: for each merge with >= 2 genomes, genome_id -> NCBI taxon id from
  `prFBA/index.parquet` (rows with NaN taxon_id dropped; this index holds one representative
  genome per unique 16S sequence, so a selected genome that is no sequence's representative
  stays unresolved), lineage by taxopy on `prFBA/nodes.dmp` + `names.dmp`; a merge needs >= 2
  resolved members (else it leaves the denominator); it is bucketed by its full member count
  (2 vs >= 3); it spans phyla when its resolved members carry > 1 distinct non-None phylum.
  `--taxon-source selection` (the taxon_id recorded per selected genome) resolves 1,399 of
  1,409 and gives 171/589 = 29.0% and 480/810 = 59.3%; it is not the published method.

Every variant tried (`stats_d20/*.json`):

| selection file | answered | n1 | n2 | n>=3 | mean | x-phylum 2 | x-phylum >=3 | rule |
|---|---|---|---|---|---|---|---|---|
| **EmilyKin `_gpu.json`** (= `_gpu.preP3.json`, identical md5) | 3,753 | 2,344 | 599 | 810 | 2.0341 | **143/451** | **462/783** | union, P3 off |
| EmilyKin `_gpu.p3.json` | 3,753 | 2,344 | 599 | 810 | 2.0341 | 145/452 | 462/777 | union, P3 on |
| EmilyKin `_gpu.legacy.json` | 3,724 | 2,132 | 1,442 | 150 | 1.4901 | 359/943 | 28/139 | legacy |
| EmilyKin `asv_reference_selection.json` (edlib hits) | 3,727 | 2,135 | 1,442 | 150 | 1.4897 | 359/943 | 28/139 | legacy |
| EmilyKin `.estimated.json` (edlib hits) | 3,753 | 2,036 | 1,543 | 174 | 1.5337 | 364/984 | 38/162 | legacy, estimated genes |
| **codif `_gpu.54d1150.json` (KBase build)** | 3,224 | 2,188 | 525 | 511 | 1.8955 | 116/397 | 191/484 | union, P3 on |
| codif `_gpu.json` (= `.preP3`) | 3,224 | 2,188 | 525 | 511 | 1.8955 | 115/393 | 200/494 | union, P3 off |
| codif `_gpu.p3.json` | 3,224 | 2,188 | 525 | 511 | 1.8955 | 115/395 | 191/484 | union, P3 on |
| codif `_gpu.legacy.json` | 3,219 | 2,007 | 1,070 | 142 | 1.4467 | 152/677 | 30/135 | legacy |
| codif `_edlib.json` | 3,219 | 2,008 | 1,069 | 142 | 1.4464 | 152/676 | 30/135 | legacy |

Only EmilyKin `_gpu.json` / `.preP3.json` (the selection the current `EmilyKin/synthetic_genomes`
manifest points to) matches every published number; `.p3.json` matches the n_genomes figures
but not the cross-phylum ones.

**Deployment A, 20-hit KBase build** (`asv_reference_selection_gpu.54d1150.json`,
`stats_d20/codif_gpu.54d1150.json`): 3,276 ASVs, 3,224 answered, 52 abstaining (38 below family,
14 reliability); n_genomes = 1: 2,188 (67.9% of answered), 2: 525, >= 3: 511; mean 1.8955 over
answered (1.8654 over all); cross-phylum 116/397 = 29.2% (two-genome) and 191/484 = 39.5%
(larger), 881 of 1,036 multi-genome merges resolved, 33 Bacteria + Archaea; provenance check
against `synthetic_genomes_54d1150/provenance`: 0 mismatches.

## 2. Selection rule

* The depth-100 build applies the **shipped exact-tie union** (select_references
  `DEFAULT_KNOBS`: `select_all=True`, `select_all_window=0.0`, species dedup on) with the **P3
  gene-richest-conspecific tie-break on** (build_synthetic_genomes_parallel's default; the
  paper's l.180 "keeping the highest-identity assembly, gene-richest on a further tie"). The
  smoke selections confirm it: every answered ASV carries `select_all:window=0.0`, the selected
  identities never spread, and `deployment_stats.py` reports `rule=exact-tie union, P3=on`.
* Old builds: the published EmilyKin 20-hit build (`_gpu.json`) was the exact-tie union with
  **P3 off** -- in all 1,046 co-optimal conspecific cases the representative is the
  lexicographically smallest genome_id (in `.p3.json`, 569 of 1,042 are not). The codif KBase
  build (`.54d1150`) was the exact-tie union with **P3 on** (its log: `P3 gene-tiebreak=on`; 539
  of 983 representatives not the smallest id). The `.legacy` / edlib files are the legacy ordered
  reducer. `prFBA/emilykin_gpu_rebuild.log` records a legacy build (3,724 built / 226 abstain)
  that the later union build replaced in `EmilyKin/synthetic_genomes/` (left untouched).
* Consequence for the paper: its 20-hit cross-phylum figures come from a P3-off selection while
  its method text describes P3 on; under P3 on (`.p3.json`) they would be 145/452 = 32.1% and
  462/777 = 59.5% (1,229 of 1,409 resolved). The n_genomes figures do not depend on P3.

## 3. The pipeline (`d100_pipeline.sh`, per deployment; resumable)

1. **align**: `align_hits.py --force-backend cuda --k2 100 --prefilter-k 500 --gpu-cand 200
   --workers 12 --save-shortlist --fasta F --taxonomy-csv T --outdir <hits_d100>`. There is no
   `--backend` flag; `--force-backend cuda` is the equivalent and aborts if no GPU. `--outdir`
   is always explicit (the default is the 20-hit `EmilyKin/bvbrc_alignment_hits`). A crash after
   the edlib stage resumes with `--shortlist-in` (tested: byte-identical hits). The stage ends
   with a check that every record holds up to 100 hits, ranks 1..n.
2. **tieset**: `tieset_gene_cache.py` writes a NEW small PGFam cache with exactly the genomes an
   exact-tie union can consult, streamed read-only from the shared 5.09 GB
   `bvbrc_cache/genome_gene_families.json` (159,573 genomes). Loading the shared cache whole
   would take ~28 GB (lists, prefetch) to ~40 GB (sets, the provider) -- measured on a 500 MB
   prefix (2.72 GB / 3.88 GB) and scaled by bytes -- and the prefetch would rewrite the shared
   file another session may be using. Exactness: in union mode the gene cache reaches selection
   only through `gene_count_fn` (P3 tie-break, est_genes, gene_set_empty), which is queried only
   for the exact-tie pool (selected + species-deduplicated genomes), a set that does not depend
   on gene counts; copying those entries verbatim gives identical answers (and None where the
   shared cache has none). `build_d100.py --tieset` re-checks it at run time and aborts on any
   query outside the tie set (smoke: 481 and 495 queries, 0 outside).
3. **prefetch**: `prefetch_gene_families.py --hits H --cache <small cache> --workers 8` fetches
   the selected genomes the shared cache lacks (BV-BRC read-only GETs) into the small cache.
4. **build**: `build_d100.py --tieset ... --hits H --gene-cache <small cache> --out-dir
   <synthetic_genomes_d100> --selection-out <hits_d100>/asv_reference_selection_d100.json
   --procs 8 --io-workers 8`. It runs `build_synthetic_genomes_parallel.main()` unchanged except
   for four run-time patches:
   * **BV-BRC User-Agent.** `build_synthetic_genomes._api_json` (the `fast` genome builder used
     by the warm step) sends urllib's bare User-Agent, which BV-BRC answers with HTTP 403
     (re-tested 2026-10-29: bare 403, with prefetch_gene_families' header 200). Without the
     patch every selected genome missing from `bvbrc_cache/kbase_genome_cache` fails to fetch
     and is silently dropped from its merge (or the ASV fails). In the codif smoke 61 of 130
     selected genomes had to be fetched (all succeeded with the patch). The permanent fix is to
     send the same header in `build_synthetic_genomes._api_json`, as
     `prefetch_gene_families.fetch_pgfams` does (not edited here).
   * the tie-set guard (above), and `gc.freeze()` before the merge pool forks (keeps the
     workers' GC from copy-on-write duplicating the parent's hits/audit; output unaffected).
   * **atomic genome-cache writes**: `build_synthetic_genomes.build_genome` writes
     `kbase_genome_cache/{gid}.json` in place, so a kill mid-write (e.g. by MemoryMax) would
     leave a truncated file in that shared cache that every later build fails to parse and
     silently drops from its merge; the wrapper writes a temp file and `os.replace`s it
     (tested from an empty scratch cache: 29 genomes fetched, 0 temp files left, selection
     identical).
   Since 2026-11-12 the builders' feature probability defaults to E3 (`--pfunc e3`; section 9).
   `d100_pipeline.sh` passes `PFUNC` as `--pfunc "$PFUNC"`, default `e3` since 2026-11-13. It used
   to pin `kn` to reproduce the shipped builds; no setting does that now (section 10).
   The warm step adds NEW per-genome files to the shared `bvbrc_cache/kbase_genome_cache/`
   (additive, as the 20-hit builds did). No ModelSEEDpy import anywhere in the chain
   (`modelseedpy imported: False` in every build log), so no PYTHONPATH pin is needed. The
   merge module is `codiffusion_bioreactor/bvbrc_to_kbase_genome.py` (d822b7f, clean); the
   pipeline logs its commit at start.
5. **stats**: `deployment_stats.py` -> `<hits_d100>/deployment_stats_d100.json`.

Caveat on P3 (pre-existing pipeline behaviour, unchanged): the tie-break counts genes only for
genomes present in the gene cache (an uncached conspecific counts 0), and the prefetch fetches
only the P3-off representatives. Which conspecific represents a species therefore depends on
the cache contents at run time; n_genomes and abstentions do not. The `.tieset.json` sidecar
records the shared cache's size, mtime and entry count at extraction.

Old entry points now forward to the new scripts (their 20-hit commands are in git history at
a90dd17): `run_newdata_gpu_align.sh` -> `run_codif_d100.sh` STAGE=align;
`run_newdata_build.sh` -> `run_codif_d100.sh` STAGE=build; `run_emilykin_gpu_rebuild.sh` ->
`run_emilykin_d100.sh` (all stages).

## 4. Smoke tests (scratch dirs; all passed)

Memory: every process tree sampled every 0.5 s (`smoke/mem_*.json`; PSS sum = physical
memory, RSS sum double-counts the fork-shared reference blob).

| test | n ASVs | wall | peak PSS sum | peak RSS sum | max 1-proc RSS | GPU |
|---|---|---|---|---|---|---|
| align codif (12 workers) | 24 | 74.6 s | 5.03 GB | 23.33 GB | 5.08 GB | 869 MiB |
| align EmilyKin (12 workers) | 48 | 111.2 s | 5.02 GB | 23.38 GB | 5.08 GB | 869 MiB |
| tieset codif / EmilyKin | 24 / 24 | 14.1 / 12.4 s | 0.86 / 0.80 GB | | | |
| build codif (8 procs; 61 genomes fetched) | 24 | 78.9 s | 1.34 GB | 3.15 GB | 0.59 GB | |
| build EmilyKin | 24 | 4.1 s | 0.25 GB | 0.53 GB | 0.09 GB | |
| whole `d100_pipeline.sh` codif | 20 | 181.9 s | 4.83 GB | 23.35 GB | 4.89 GB | 869 MiB |
| whole `d100_pipeline.sh` EmilyKin | 20 | 132.5 s | 5.02 GB | 23.35 GB | 5.10 GB | 869 MiB |

Checks: every record holds exactly 100 hits with ranks 1..100 (24 + 48 + 20 + 20 records);
0 GPU re-score fallbacks; the 20-ASV and 24-ASV codif runs give identical records and
selections for the shared ASVs; resume from the shortlist is byte-identical; selection = union
P3 on. Smoke stats (`smoke/stats_*20_d100.json`, the FASTA's first 20 ASVs, which are the most
abundant): codif mean n_genomes 6.25 (the same ASVs averaged far less at 20 hits; e.g. 12 -> 29
and 6 -> 17 genomes where the exact-tie cloud exceeded 20); EmilyKin mean 1.41, 3 abstaining.
At 20 hits the tie cloud was cut at the 20th hit in 148 codif and 89 EmilyKin records.

## 5. Resources (estimates for the orchestrator)

Alignment (per-ASV costs measured on the smokes, under the other session's ~32-core load):
edlib 14.8 worker-s/ASV (codif) and 16.9 (EmilyKin); traceback 2-4 worker-s/ASV; ~1 min fixed.
At 12 workers: **A 1.3-1.45 h, B 1.75-1.9 h.** The depth-100 region-validation run
(`region_validation_d100/data/hits/raw/method_selection.json`: 17,470 s edlib for 60,461
amplicons of mean 542 bp at 40 workers) scales to ~8.0 worker-s/ASV at the ASVs' 374-377 bp
mean, about half (47-54%) of the smoke rate, so on a less loaded machine the alignment may take
only A ~0.8 h, B ~0.95 h.
Peak memory: 5.0 GB PSS measured; full run est. **~10 GB** (+ ~2-3 GB of per-ASV
shortlists/records at 3,950 ASVs and copy-on-write in the 12 traceback workers).

Build: merge time is superlinear in the union size (profiled: 1 genome 0.2-0.8 s; 29 genomes
45.7 s); est. merge 0.2-0.5 h per deployment at 8 procs, plus fetching the selected genomes
missing from `kbase_genome_cache` at the measured 2.6-6.6 genomes/s with 8 threads (est. 1-3k
genomes, 3-20 min),
plus prefetch and tie-set extraction (minutes): **A ~0.5-1 h, B ~0.3-0.8 h.** Peak memory: 1.34
GB PSS measured on 24 heavy ASVs; full run est. **~8 GB typical, 20 GB worst case** (parent:
100-hit hits + audit ~1.1 GB, scaled from 0.21 GB measured for the 79,000-hit 20-hit EmilyKin
file; tie-set gene sets ~2-4 GB at 7.76 GB per GB of cache JSON; 8 merge workers at 0.6 GB
measured for a 29-genome union, up to ~2 GB for a 100-genome union).

**Chain: est. 2.5-5.5 h; peak est. 20 GB -> MemoryMax=26G** (x1.3), CPUQuota=1200%, GPU.
Disk: est. 30-60 GB of synthetic genomes per deployment (20-hit: 24 GB codif, 31 GB EmilyKin;
the smoke's heavy ASVs averaged 15 MB vs 7.4 MB per genome at 20 hits).

## 6. Launch (orchestrator)

```
bash /home/freiburger/Documents/prFBA/deployment_d100/launch_dep100.sh
```
runs (after preflight: MemAvailable >= 40 GB, GPU free >= 4 GB, disk >= 150 GB, unit absent):
```
systemd-run --user --unit=dep100-chain --collect -p MemoryMax=26G -p CPUQuota=1200% \
  -p WorkingDirectory=/home/freiburger/Documents/prFBA \
  -E OMP_NUM_THREADS=1 -E OPENBLAS_NUM_THREADS=1 -E MKL_NUM_THREADS=1 -E PYTHONUNBUFFERED=1 \
  /bin/bash -c "/bin/bash /home/freiburger/Documents/prFBA/deployment_d100/dep100_chain.sh > /home/freiburger/Documents/prFBA/deployment_d100/logs/dep100-chain.log 2>&1; echo EXIT=\$? >> /home/freiburger/Documents/prFBA/deployment_d100/logs/dep100-chain.log"
```
Logs: `logs/dep100-chain.log`, `logs/codif_d100.log`, `logs/emilykin_d100.log` (each stage logs
its wall time and the unit's cgroup memory.peak). That counter includes page cache charged to
the unit (the 5 GB gene-cache stream, the reference JSONs, written genomes), so it will read
higher than the PSS-based estimates in section 5 and is not an RSS overrun; MemoryMax reclaims
page cache before it kills anything. The chain also waits before each heavy stage
until MemAvailable >= 32 GB. Relaunching after a failure resumes (finished alignments are
skipped; a partial one restarts from its saved shortlist; the build reruns).

## 7. Statistics after the run (analyst)

The chain already writes `<hits_d100>/deployment_stats_d100.json`; to recompute:
```
cd /home/freiburger/Documents/prFBA
PY=~/Documents/py_venv/bin/python
E=/home/freiburger/Documents/EmilyKin; C=/home/freiburger/Documents/codiffusion_bioreactor
$PY deployment_d100/deployment_stats.py --selection $E/bvbrc_alignment_hits_d100/asv_reference_selection_d100.json \
    --hits $E/bvbrc_alignment_hits_d100/asv_top20_alignment_hits.json \
    --provenance-dir $E/synthetic_genomes_d100/provenance --out deployment_d100/stats_d100_emilykin.json
$PY deployment_d100/deployment_stats.py --selection $C/bvbrc_alignment_hits_d100/asv_reference_selection_d100.json \
    --hits $C/bvbrc_alignment_hits_d100/asv_top20_alignment_hits.json \
    --provenance-dir $C/synthetic_genomes_d100/provenance --out deployment_d100/stats_d100_codif.json
```
Fields: `answered`, `abstained` (+ `abstained_by_flag`), `n_genomes_1/2/3plus` (+ `_pct_of_answered`),
`mean_n_genomes_answered` (the paper's mean), `cross_phylum.{n=2,n>=3}.{cross_phylum,resolved}`
(the paper's numerators/denominators), `cross_phylum.all_multi.{resolved,multi_genome_merges,
bacteria_archaea}` (the "1,234 of 1,409" and "18"), `hits_audit` (hits per record), `selection_rule`,
`provenance_check` (must be 0 mismatched). Compare EmilyKin d100 against `stats_d20/emilykin_gpu.json`
(published) and codif d100 against `stats_d20/codif_gpu.54d1150.json` (KBase build). Places that
quote the 20-hit numbers: prFBA.tex and prFBA_andrew_voice.tex lines 211, 223, 290, 1072, 1127,
1134; `presentations/make_deck.py` DEPLOY_N1, XPHYLUM_DEPLOY, KBASE_ASVS/BUILT/ABSTAIN/WS,
DEPLOY_ASVS and the '20-hit records' footers/notes (lines 249-280, 815-816, 904, 1171-1190,
2882-2905).

## 8. KBase upload (done 2026-11-04 with `kbase_upload_d100.sh`; the plan below predates it)

Deployment A only (the uploader's GenomeSet description is hard-coded for the codiffusion V4-V5
ASVs). Run after the chain, by hand:
```
cd /home/freiburger/Documents/prFBA
~/Documents/py_venv/bin/python upload_kbase_genomes.py \
    --genome-dir /home/freiburger/Documents/codiffusion_bioreactor/synthetic_genomes_d100 \
    --ws-id 261948 --genomeset-name codiffusion_synthetic_genomes_prFBA_d100 --dry-run
# then the same command without --dry-run
```
**This writes into workspace 261948** (freiburgermsu:narrative_1780783998282, which holds the
20-hit build: 3,224 Genomes + GenomeSet 261948/3226/1 per `newdata_upload_54d1150.log`).
Before running it, note (from the local logs and code; the workspace itself was not queried):
* The uploader treats a Genome as "already present" by object NAME and names new objects by ASV
  md5. Per `synthetic_genomes_54d1150/kbase_upload_log.json` / `kbase_rename_log.json` the 20-hit
  objects were renamed to iterativeIDs (e.g. `Methanobacterium.1`), so the d100 genomes would
  not collide: all ~3.2k built d100 genomes would be saved as NEW objects (md5 names) next to the
  20-hit set (the 20-hit upload was 14.83 GB / 184 min). Any 20-hit object still carrying its
  md5 name would instead be "adopted", i.e. the d100 GenomeSet would point at the 20-hit genome.
  A fresh narrative (new `--ws-id`) avoids both.
* `--selection` defaults to the manifest's `_meta.selection` (the d100 selection). The optional
  `mapping_provenance.json` (commit + inputs for the object metadata) is written by
  `write_mapping_provenance.py`, which still hard-codes the 20-hit alignment description
  ("exhaustive GPU Smith-Waterman (gpu_align.py --full) ... top-20", l.334) and the codiffusion
  dataset name (l.325); fix those before generating it for a d100 upload:
  `write_mapping_provenance.py --hits <hits_d100>/asv_top20_alignment_hits.json --selection
  <hits_d100>/asv_reference_selection_d100.json --genome-dir <synthetic_genomes_d100> --asv-fasta
  .../dna-sequences_codif_all.fasta --taxonomy-csv .../taxonomy_codif_all.csv --alignment-stats
  <hits_d100>/method_selection.json --gene-cache <hits_d100>/gene_families_tieset.json`
  (tested on the smoke build: runs, 0.6 s).

## 9. Feature probability: E3 is the build default (2026-11-12)

The top-level `README.md` ("Feature probability of the synthetic genomes: E3") has the full
description. What it means for these deployments:

* **What changed.** `build_d100.py` passes everything through to the parallel builder, whose
  `--pfunc` now defaults to **e3**. `d100_pipeline.sh` passes the env variable `PFUNC` to the
  build stage as `--pfunc "$PFUNC"`. Until 2026-11-13 it pinned `PFUNC=kn` to reproduce the shipped
  k/n builds. Since the BV-BRC-only rule (section 10) no setting reproduces them, so the default
  is `e3`, the builders' default. The bullets below were written for the E3 build of 2026-11-12;
  section 10 lists what changed since.
  * **`PFUNC=kn`**: k/n everywhere. The manifest `_meta` gains
    `pfunc: {estimator: kn, estimator_id: k/n, gating_rule, min_n}` and the two input stamps
    (`selection_stamp`, `source_digests`, with `SYN/pfunc/source_digests.json`; top-level README).
  * **`PFUNC=e3` (the default) with a NEW `SYN`** gives an E3 build:
    * For ASVs with >= 2 built source genomes, the feature (and CDS) probabilities become the E3
      estimate, a 16S-neighbour vote in PGFam space plus an NCBI family/order prior.
    * Features without a PGFam keep k/n. Before the BV-BRC-only rule that was about 55% of
      features, because of the RefSeq/PATRIC vocabulary mix in `fast_build_genome`'s CDS query.
      Now it is about 17% (20-ASV smoke, section 10).
    * n = 1 genomes are unchanged.
    * The build stage **refuses to start** (exit 2, nothing written) when `SYN` exists and is not
      empty, for any `PFUNC` since 2026-11-13. The shipped `synthetic_genomes_d100/` is such a
      directory, and so is a half-written build after a crash (delete it first). `SYN` must still
      end in `_d100` (the pipeline's guard), e.g. `synthetic_genomes_bvbrc_d100`.
    * The selection audit and the stats go to `$SYN/asv_reference_selection_d100.json` and
      `$SYN/deployment_stats_d100.json`, not `$OUT` (every `PFUNC` since 2026-11-13). So the
      shipped audit, which the shipped manifest points at, is not rewritten.
    * The build stage re-runs selection from `$OUT`'s hits and tie-set cache. Its prefetch adds any
      genome still missing from that cache, as every re-run does. So if a genome that failed to
      fetch before now succeeds, the build can select differently from the shipped one.
    * `run_codif_d100.sh` / `run_emilykin_d100.sh` export the shipped `SYN`, so their build and
      all stages are always refused. Call the pipeline directly with a NEW `SYN` instead, e.g. for
      EmilyKin:
      ```
      systemd-run --user --unit=e3-build-ek-$(date +%s) --collect -p MemoryMax=40G -p MemorySwapMax=0 \
        -p WorkingDirectory=/home/freiburger/Documents/prFBA --wait --pipe \
        -E NAME=emilykin_bvbrc_d100 -E STAGE=build \
        -E FASTA=/home/freiburger/Documents/EmilyKin/asvs.fasta \
        -E TAXCSV=/home/freiburger/Documents/EmilyKin/bvbrc_embedding_hits/inputs/taxonomy.csv \
        -E OUT=/home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits_d100 \
        -E SYN=/home/freiburger/Documents/EmilyKin/synthetic_genomes_bvbrc_d100 \
        /bin/bash /home/freiburger/Documents/prFBA/deployment_d100/d100_pipeline.sh
      ```
      This command form (then with `-E PFUNC=e3` and an `_e3_d100` SYN) was tested on 2026-11-12
      with codiffusion inputs, a scratch `OUT` holding a 37-ASV hits subset, a scratch `SYN` and a
      scratch `GENOME_CACHE_DIR` (plus `-E PYTHONHASHSEED=0 -E PROCS=4 -E IOW=4`), not on the
      deployment. Its output was byte-identical (33/33 ASVs, all four per-ASV files) to
      `reweight_build_e3.py` applied to the same subset's k/n build. That was before the
      BV-BRC-only rule; the pipeline itself was not re-run on a deployment since.
* **Why n >= 2 only.** Stage O established the flux gain only for multi-genome selections:
  -0.0172 [-0.0221, -0.0118] of D, against -0.0048 [-0.0113, +0.0063] at n = 1. n = 1 is 62.5%
  (EmilyKin) and 67.9% (codiffusion) of the built genomes.
* **Gate.** The production estimator on the 10k panel (exclude_0.987, n >= 2), E3(species) minus
  k/n: Brier -0.0530 [-0.0547, -0.0513] and AUROC +0.0538 [+0.0517, +0.0560]
  (`region_validation_d100/data/pfunc_e3_gate.json`).
* **Inputs an e3 build needs.**
  * The shared PGFam cache `bvbrc_cache/genome_gene_families.json` (`--pfunc-pgfam-cache`). It is
    streamed and never written. It is NOT the tie-set `--gene-cache`.
  * The clade catalogue `bvbrc_cache/e3_clade_catalogue.json`. It must match the cache. Otherwise the
    run stops before selection and prints the rebuild command:
    `/home/freiburger/Documents/py_venv/bin/python /home/freiburger/Documents/prFBA/pfunc_estimator.py rebuild-command`;
    the rebuild takes ~45 s and 2.7 GB.
  * The pipeline's `--gene-cache` is the tie-set file, so its prefetch never stales the catalogue.
  * PGFam sets missing from the cache are fetched into `OUT/pfunc/pgfam_fetched.json`
    (`--pfunc-no-fetch` skips this). With every needed genome cached, no network is used.
* **Outputs added:**
  * `OUT/pfunc/pgfam_store.json` and `pgfam_fetched.json`;
  * for E3 genomes, `provenance/{asv}.json` keys `probability_kn` (the replaced k/n) and `pfunc`
    (diagnostics);
  * a `pfunc` summary per manifest entry;
  * `_meta.pfunc`: estimator, parameters, counts, and `code_provenance`, which holds the sha256 of
    `pfunc_estimator.py`, the builders and `build_d100.py`, plus git HEAD and a dirty flag;
  * in every mode (k/n too): `_meta.selection_stamp` (sha256, size and mtime_ns of the audit) and
    `SYN/pfunc/source_digests.json` (`{genome_id: source_digest}` of every built source genome,
    recorded in `_meta.source_digests`). The retro tool checks both.
* **Tie-set guard with `--selection`.** No selection runs, so there are no gene_count_fn queries.
  `build_d100.py` therefore does not install the guard and logs "tie-set guard not applicable with
  --selection". Before this change it logged "tie-set guard on" even though the guard did nothing.
* **The shipped builds, retroactively: done 2026-11-12, now refused.** `reweight_build_e3.py`
  wrote `synthetic_genomes_d100_e3` from each shipped build. Those re-weightings, like the
  shipped builds, merge PATRIC and RefSeq features. Since 2026-11-13 the tool refuses any source
  build without the BV-BRC-only annotation stamp (section 10), so the command below now exits
  with that refusal and writes nothing. It is kept as the record of how the re-weightings were
  made (copy-on-write; it wrote a NEW directory and never touched the shipped one):
  ```
  systemd-run --user --unit=e3-retro-ek-$(date +%s) --collect -p MemoryMax=40G -p MemorySwapMax=0 \
    -p WorkingDirectory=/home/freiburger/Documents/prFBA --wait --pipe \
    /home/freiburger/Documents/py_venv/bin/python /home/freiburger/Documents/prFBA/reweight_build_e3.py \
      --build-dir /home/freiburger/Documents/EmilyKin/synthetic_genomes_d100 \
      --out-dir   /home/freiburger/Documents/EmilyKin/synthetic_genomes_d100_e3 --procs 8
  ```
  * Measured on all 1,410 EmilyKin n >= 2 ASVs (outputs since deleted): 354 s. Parent RSS 3.1 GB;
    each worker's PSS 0.6-0.8 GB. The store was 0.69 GB on disk with 20,697 of 21,003 genomes.
  * About 21 GB of genomes were written. The cgroup's memory.peak reached 40G through page cache,
    not process memory.
  * What the tool checks, in order (details in the top-level README):
    1. Run level, the selection audit. A stamped build (`_meta.selection_stamp`): the audit's sha256
       must equal the stamp, else every built ASV in scope (n = 1 too) is copied unchanged (k/n)
       with `inputs drifted since build: selection audit`. A legacy build: the run stops before
       writing anything if the audit's mtime is newer than the source manifest's, unless
       `--allow-audit-newer`.
    2. Per ASV with n >= 2, each failure copying the ASV unchanged (k/n):
       * `selected genomes`: the record's selected ids equal the built plus failed genomes;
       * every built source genome is in `--genome-cache-dir` (else `source genome(s) missing from
         the genome cache`);
       * `function_prevalence` of the cached source genomes equals the provenance's;
       * `source digest`: the cached source genomes' digests equal the build's
         `pfunc/source_digests.json`. For a legacy build without that file, the check is instead
         `source genome cache newer than build`: no cached source-genome file may be newer than
         the source manifest, unless `--allow-cache-newer`.
    The audit's `decisions[]` (the neighbour list) are covered only by check 1.
  * The shipped builds predate the input stamps, so the two legacy rules apply to them. Neither
    fires on them today (checked 2026-11-12, and every check above passes on all 1,410 EmilyKin
    and 1,034 codiffusion n >= 2 ASVs). Both audits are older than their manifests: EmilyKin
    22:24 vs 22:32, codiffusion 21:21 vs 21:32, both on 2026-10-29. No n >= 2 ASV uses a cached
    source genome newer than its manifest. A copy of the audit or of the genome cache that does
    not preserve mtimes will trip them.
  * Re-run `write_mapping_provenance.py` on the new directory.
* **Build memory.** An e3 build adds the PGFam store to the parent before the merge pool forks,
  about 3-3.5 GB for a whole deployment. `gc.freeze()` covers it.
  * The store holds the clade members of every selected rep, because the prior is anchored on the
    first rep that builds. Over the n >= 2 ASVs that is EmilyKin 23,709 genomes and codiffusion
    23,262 (~13% more than the retro tool's exact scope), all in the shared cache.
  * The full `--hits` e3 path was not measured at deployment scale.
* **KBase.** For E3 builds, the uploader records the estimator in each object's provenance
  (`method_params[0].feature_probability`) and ws metadata (`probability_estimator`); k/n uploads
  are unchanged. Objects are named by ASV md5 and an existing name is skipped or adopted. So
  uploading an E3 build into a workspace that already holds the same ASVs' k/n genomes would keep
  the k/n objects. Use a **new workspace** (`--ws-id`) or new object names (section 8).

## 10. Source annotation: BV-BRC (PATRIC) only (2026-11-13)

The top-level `README.md` ("Source annotation of the synthetic genomes") has the full
description. BV-BRC serves two annotations per genome: its own PATRIC (RASTtk) one, and the
RefSeq (NCBI PGAP) one it imports. Up to `3ce12d0` the fast builder took both. Every build now
takes the PATRIC CDSs only, with no option to include RefSeq:
- `fast_build_genome` queries `eq(annotation,PATRIC)`.
- Every cached or fetched source genome passes through `build_synthetic_genomes.bvbrc_only`.
  That includes `build_d100.py`'s atomic cache patch, which now calls the same
  `load_source_genome` / `fetch_source_genome` helpers.
- The shared `bvbrc_cache/kbase_genome_cache/` is never rewritten. Its older mixed files are
  filtered on read.
- Every manifest records `_meta.annotation`.

What this means here:
- **The shipped builds and their uploads mix both vocabularies.** These are
  `synthetic_genomes_d100/` (both deployments), their KBase GenomeSets 261948/6467/1 and
  275568/3757/1, and the `synthetic_genomes_d100_e3/` re-weightings. They are reproducible only
  from git at or before `3ce12d0`.
- **Re-running `d100_pipeline.sh` makes a NEW BV-BRC-only build into a NEW `SYN`.** This holds for
  any `PFUNC`; the default is now `e3`. The audit and the stats go into `$SYN`. An existing
  non-empty `SYN` is refused (exit 2, nothing written), and so are the build and all stages of
  `run_codif_d100.sh` / `run_emilykin_d100.sh`. Selection reads PGFam sets only, so it is
  unaffected by the rule.
- **`reweight_build_e3.py` refuses the shipped builds.** Rebuild them instead.
- **Smoke test.** 20 EmilyKin d100 audit records (`build_d100.py --selection`, the same selected
  genomes as the shipped build; 2026-11-13; scratch output):
  - every feature and CDS carries a non-empty `PATRIC_id` alias;
  - features fell from 63,051 to 34,864 (median 0.556 of the shipped count per genome, range
    0.35-1.0);
  - the share of features resolving to a PGFam rose from 45.5% to 82.2%, and is 82.8% over the
    12 E3-reweighted n >= 2 genomes;
  - all 8 n = 1 genomes' function sets equal their source genome's PATRIC products.

## Files

| file | role |
|---|---|
| `deployment_stats.py` | selection (+hits, +provenance) -> stats JSON; the published counting method |
| `tieset_gene_cache.py` | small exact PGFam cache for an exact-tie union build |
| `build_d100.py` | build_synthetic_genomes_parallel + User-Agent fix + tie-set guard + gc.freeze + atomic cache writes; `--pfunc` default e3 passes through (section 9); BV-BRC PATRIC annotation only (section 10) |
| `d100_pipeline.sh` | the 5-stage pipeline (env-configured) |
| `run_codif_d100.sh`, `run_emilykin_d100.sh` | deployment A / B entry points (STAGE=align/build/all) |
| `dep100_chain.sh` | A then B (the unit's body) |
| `launch_dep100.sh` | preflight + the systemd-run command (unit `dep100-chain`) |
| `stats_d20/` | the 20-hit stats for every selection variant |
| `smoke/` | smoke-test memory samples, logs and stats |
