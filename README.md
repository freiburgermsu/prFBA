# prFBA — 16S rRNA embedding space

A dense vector representation of every unique 16S rRNA gene in BV-BRC, built with
Nucleotide-Transformer-v2-500m on an RTX 5070 Ti. **Read `METHODS.md` for the full write-up**
(model, parameters, density/sparsity geometry, and how to reference new sequences).

## Contents

| | |
|---|---|
| `embeddings.f16.npy` | (460,116 × 1024) L2-normalized float16 — the embedding space (cosine = dot) |
| `index.parquet` | per-unique-sequence metadata (row-aligned to the matrix) |
| `members.parquet` | every input sequence → its row (reverse lookup) |
| `manifest.json` | provenance: model, params, counts, dedup, timings, environment |
| `density_stats.json` | geometry metrics (kNN density, PCA, intrinsic dim, strata, …) |
| `fig_*.png` | density / PCA / length-strata / UMAP figures |
| `nt_embed.py` | shared encoder (used by both scripts below) |
| `embed_16s.py` | build the store (resumable, GPU) |
| `analyze_space.py` | density/sparsity analysis (+ `--min-len/--max-len/--exclude-atypical`, `--tag`) |
| `query.py` | reference new sequences by cosine similarity (+ same length filters) |
| `watch_and_run.sh` | waits for the dataset, then runs embed → analyze |
| `tests/` | correctness validators (analysis + query), all passing |
| **`REPORT_amplicon_hits.md`** | **findings: experimental V4–V5 ASVs vs the embedding space** |
| `insilico_pcr.py` | excise a 16S sub-region (V4–V5, 515F/926R) from references for region-matching |
| `hit_amplicons.py` | embed ASVs → top-20 cosine hits + BV-BRC metadata (+ naive comparison) |
| `concordance.py` | multi-rank taxonomic concordance vs MiDAS (uses NCBI taxdump via taxopy) |
| `align_hits.py` | **hardware-aware alignment dispatch** (GPU `gpusw` / edlib→Biopython / full Biopython) — see below |
| `edlib_biopython_hits.py` | two-stage CPU pipeline (edlib prefilter → Biopython local SW); reusable stages + enrichment |
| `gpu_align.py` | original inline-CUDA-kernel GPU Smith-Waterman (superseded as a tier by the `gpusw` package) |
| `asv_top20_hits.json` | **per-ASV top-20 BV-BRC hits**: cosine + organism/genome_name/taxon_id/genome_id/feature_id/n_genomes |
| `asv_summary.csv` / `asv_concordance.csv` / `concordance_by_rank.json` / `findings_stats.json` | per-ASV + aggregate results |

> **Data files are git-ignored.** `embeddings.f16.npy` (≈0.9 GB), `*.parquet`, `v4v5_refs.json`, the
> NCBI taxdump (`*.dmp`), and `umap_*.npy` exceed GitHub's file limits, so they are not committed.
> They are fully regenerable — `python embed_16s.py` rebuilds the store, `python analyze_space.py` the
> analysis. The committed `manifest.json`, `density_stats.json`, `fig_*.png`, the `asv_*` result files,
> and `REPORT_amplicon_hits.md` capture the results.

## Experimental ASV analysis (V4–V5)

Matching the study's 3,276 experimental 16S ASVs against this space — **region-aware**, since the
amplicons are ~373 bp V4–V5 (515F/926R), not full genes. See **`REPORT_amplicon_hits.md`** for the
full findings (median best cosine 0.9985; the dominant *Methanobacterium* methanogenic core matched at
cosine 1.0; abundance-weighted genus concordance with MiDAS ≈ 0.89).

```bash
# 1. excise the V4–V5 region from every reference (in-silico PCR)
python insilico_pcr.py --fasta .../model_inputs/BV_BRC_16S.json --out v4v5_refs.json
# 2. embed the region-matched references
python embed_16s.py --fasta v4v5_refs.json --metadata .../model_inputs/16S_metadata.json --outdir v4v5_store
# 3. match the ASVs (region-matched + naive comparison) -> asv_top20_hits.json, asv_summary.csv
python hit_amplicons.py --amplicons .../new_data/dna-sequences_codif_all.fasta \
    --store v4v5_store --naive-store . --asv-ids .../asv_IDs.csv --taxonomy .../taxonomy.csv --topk 20 --outdir .
# 4. taxonomic concordance vs MiDAS
python concordance.py
```

## Alignment-based hits (edlib top-500 prefilter → hardware-accelerated scoring)

For exact Smith-Waterman hits of the ASVs against the BV-BRC references (rather than
embedding cosine), `align_hits.py` uses a fixed, hardware-**independent** prefilter and a
hardware-**dependent** scoring backend — same local-SW scheme (match +2 / mismatch −3 /
gap_open −5 / gap_extend −2) everywhere, so the outputs are directly comparable:

1. **edlib top-500 prefilter (always, every machine).** An edlib HW edit-distance scan over
   every unique reference keeps the 500 lowest-distance candidates per ASV — discarding the
   ~99.9% of obviously-irrelevant references before any Smith-Waterman work.
2. **Smith-Waterman scoring of just those 500 candidates**, on:

| probe | scoring backend |
|---|---|
| CUDA present | GPU SW kernel (`gpu_align.py`, CuPy/NVRTC) scoring only each ASV's 500 candidates |
| Apple Metal/MPS present | detected; no Metal SW kernel yet → falls back to CPU (framework in place) |
| otherwise | Biopython CPU local SW on the 500 |

3. **Biopython alignment of the top 100** for exact % identity / aligned length / reference
   coordinates (`--k2 100`, the standard record depth since the depth-100 rerun,
   `region_validation_d100/RUNBOOK_d100.md`). On the GPU path the GPU score, bit-exact to
   Biopython's, ranks the shortlist under the rescore's own tie-break and only its top 100
   are aligned; each alignment's score is checked against the GPU's, and an ASV whose
   scores disagree falls back to a Biopython re-score of the GPU's best 200 (`--gpu-cand
   200`). The record is the one a re-score of all 500 gives
   (`region_validation_d100/GPU_CAND_200.md`). On the CPU path Biopython scores all 500.
   The output file keeps its historical name `asv_top20_alignment_hits.json`.

Why top-500 (not exhaustive): a deep prefilter captures the best representative(s) per ASV —
including the large equal-SW-score tie clusters (empirically up to ~200 genomes tie at the top
score) that decide which near-identical genome is named the representative — **without** paying
to Smith-Waterman the whole 459k-reference DB. Validated: edlib-500 → GPU/CPU reproduces the
GPU-exhaustive best score exactly, and the exhaustive top-20 up to equal-score tie arbitration.

```bash
python align_hits.py --explain               # print the hardware probe + chosen backend, run nothing
python align_hits.py                          # auto-detect backend, run on all ASVs
python align_hits.py --fasta X.fasta --taxonomy-csv X.csv --outdir DIR
python align_hits.py --force-backend cpu      # override the scoring backend (cuda|metal|cpu)
python align_hits.py --prefilter-k 500        # edlib shortlist depth (default 500)
python align_hits.py --k2 20 --gpu-cand 100   # the earlier depth-20 record
```

Every path emits `asv_top20_alignment_hits.json` + `asv_alignment_summary.csv` (identical
schema) plus `method_selection.json` (the probe, decision, and run stats). `gpu_align.py` holds
the CUDA SW kernel (also runnable standalone for an exhaustive all-reference mapping). The
backend decision is unit-tested in `tests/test_select_method.py`.

## MAG whole-genome hits — skani ANI vs reference databases (GTDB)

`mag_skani_hits.py` extends prFBA from 16S amplicons to **whole MAGs**: instead of referencing a
373 bp amplicon into BV-BRC by 16S identity, it references a genome into a sketched database by
**skani average nucleotide identity** (skani v0.3.2 static binary in `bin/skani`, git-ignored).
Assessed on the 291 EmilyKin dereplicated MAGs against two databases:

| DB | contents | taxonomy join |
|---|---|---|
| `skani_db/skani_gtdb_r226-v0.3` | GTDB **R226** pre-sketched skani DB (skani authors, 38.7 GB tarball) | `skani_db/gtdb_r226/bac120+ar53_taxonomy_r226.tsv.gz` (accession → lineage) |
| `skani_db/bvbrc_sludge` | 2,160 BV-BRC `sludge_genomes` sketched locally | `~/Documents/sludge/gtdb_taxonomy.json` (genome_id → ranks) |

```bash
./run_emilykin_mag_skani.sh          # sketch (if needed) + both searches + concordance
# or the stages individually:
python mag_skani_hits.py sketch --fasta-dir DIR --out skani_db/NAME
python mag_skani_hits.py search --mags DIR --db skani_db/NAME --db-name NAME --db-kind gtdb|bvbrc \
    --taxonomy TSV_OR_JSON --gtdbtk gtdbtk.bac120.summary.tsv --outdir OUT --topk 20
python mag_skani_hits.py concordance --hits OUT/mag_top20_skani_hits.NAME.json --gtdbtk TSV --outdir OUT
```

Outputs per DB (in `~/Documents/EmilyKin/mag_skani_hits/`): `mag_top20_skani_hits.<db>.json`
(**schema-compatible with `asv_top20_alignment_hits.json`** — `identity` = ANI/100, `align_score`
= ANI, `aligned_len` = AF_query·genome_len — so `select_references.py` consumes it unchanged with
the ANI-recalibrated `ANI_KNOBS`), `mag_skani_summary.<db>.csv` (AF-aware best hit + ANI novelty
class per MAG), `mag_skani_concordance.<db>.csv` + `concordance_by_rank.<db>.json` (per-rank
best-hit lineage vs the GTDB-Tk call), `run_stats.<db>.json`.

The **best hit is AF-aware** (highest ANI with `AF_query ≥ 15%`, else `low_AF_unresolved`) — the
guard from the EmilyKin server handoff, where one spurious high-ANI/near-zero-AF reference
outranked the true match for ~24 MAGs. ANI classes follow Table S6: `>=95` same species,
`85-95` putatively species-distinct, `<85` divergent, `low_AF_unresolved`, `no_hit`.
BV-BRC hits carry real BV-BRC genome_ids, so the synthetic-genome merge
(`select_references.py` → `build_synthetic_genomes.py`, PGFam gene provider) works end-to-end.

## Feature probability of the synthetic genomes: E3 (default since 2026-11-12)

The synthetic-genome builders (`build_synthetic_genomes.py`, `build_synthetic_genomes_parallel.py`,
`deployment_d100/build_d100.py`) give every merged feature a probability `p`, which prFBA's objective
turns into a reaction cost. The merge's own value is **k/n**, the share of the ASV's n built source
genomes carrying the feature's BV-BRC product string. **E3** (`pfunc_estimator.py`, selected in Stage E
and tuned in Stage O, `region_validation/shared_reference/`) replaces it with an identity-weighted vote
of the ASV's 16S neighbours. The neighbours are the selected references, then the stage-0 survivors of
the first 100 audit hits, one per BV-BRC species, K = 100, w = exp(-(id1 - id)/0.02). The vote asks
whether each neighbour's BV-BRC PGFam set contains the feature's PGFam. It is shrunk toward the PGFam
prevalence of the anchor's NCBI family, else order (alpha = 2; at least 3 member species, at most 100,
one genome per species). The anchor is the first selected reference that built. A reference that failed
to build stays a neighbour.

**Scope: n >= 2 only.** E3 replaces k/n only for ASVs with at least 2 built source genomes. At n = 1
the probability stays k/n (= 1.0), because Stage O established the flux gain only there:
- n >= 2: truth-absent flux falls by 0.0172 [0.0118, 0.0221] of D (held-out confirm half).
- n = 1: confirm half -0.0048 [-0.0113, +0.0063]; the tune half's gain did not replicate.

n = 1 is 62.5% (EmilyKin) / 67.9% (codiffusion) of the deployed genomes, so E3 reaches about a third of
them.

**Validation gate (passed; `region_validation_d100/data/pfunc_e3_gate.json`).** This is the production
implementation on the 10k panel, in its own PGFam namespace. Setting: self-exclusion at 0.987,
n >= 2; 19,158 amplicons, 6,776 source genomes. E3 (species dedup) minus k/n:
- Brier -0.0530 [-0.0547, -0.0513]
- AUROC +0.0538 [+0.0517, +0.0560]

CIs are 95%, by bootstrap over source genomes. Two caveats from the gate review:
- The gain shrinks with n: dBrier -0.101 at n = 2, -0.062 at n = 3, -0.024 at n >= 4.
- About 45% of the Brier gain is calibration that a recalibrated k/n would also get. The AUROC gain
  survives recalibration.

**Namespace and reach.** E3 works on PGFams. A merged feature's family comes from its own `PGFAM`
alias, else from the PGFams of source-genome CDSs with the same product string. Features with neither
keep k/n, and that is **about 55%** of features (55.2% over all 1,410 EmilyKin d100 n >= 2 ASVs). The
cause is the RefSeq/PATRIC vocabulary mix in `fast_build_genome`'s CDS query:
- About half of the cached genomes carry both PATRIC CDSs (with `pgfam_id`) and RefSeq CDSs (with no
  PATRIC id and no `pgfam_id`).
- A merged feature that a RefSeq CDS introduced has no PGFAM alias.
- The product-string match recovers only about 2.4% of features.

The query was deliberately left unchanged, because changing it would change the shipped gene union. A
reaction's p is the maximum over its genes, so E3 and k/n values mix there. With `--builder full`, no
feature should resolve to a PGFam, so E3 would keep k/n everywhere. This is inferred from the converter's
code and was not run; the builder warns.

**Flags** (all three builders; `build_d100.py` passes them through):

| flag | meaning |
|---|---|
| `--pfunc {e3,kn}` | default **e3**; `kn` reproduces the earlier builds (genome, provenance and taxonomy files and the manifest entries byte-identical; the manifest `_meta` gains `pfunc: {estimator: kn, estimator_id: k/n, gating_rule, min_n}` plus the two input stamps below) |
| `--pfunc-min-n 2` | E3 only for ASVs with at least this many built source genomes |
| `--pfunc-dedup {species,genome}` | neighbour dedup (default species) |
| `--pfunc-pgfam-cache` | shared PGFam cache `bvbrc_cache/genome_gene_families.json`, streamed with ijson, never loaded whole or written (not `--gene-cache`, which a d100 build points at its small tie-set file) |
| `--pfunc-catalogue` | E3 clade catalogue `bvbrc_cache/e3_clade_catalogue.json` |
| `--pfunc-allow-stale-catalogue` | continue with a catalogue that does not match the PGFam cache |
| `--pfunc-no-fetch` | do not fetch PGFam sets missing from the cache from BV-BRC |
| `--pfunc-fetch-workers 8` | fetch threads |

The catalogue check runs **before selection**. An e3 run stops, printing the rebuild command, if the
catalogue is missing or no longer matches the PGFam cache (size and mtime). The rebuild takes about
45 s and 2.7 GB. It writes `bvbrc_cache/e3_clade_catalogue.json` (the `--out` default):
```bash
/home/freiburger/Documents/py_venv/bin/python /home/freiburger/Documents/prFBA/pfunc_estimator.py rebuild-command   # prints the line below, with explicit --cache/--out
systemd-run --user --unit=e3-catalogue-$(date +%s) --collect -p MemoryMax=8G -p MemorySwapMax=0 \
  -p WorkingDirectory=/home/freiburger/Documents/prFBA --wait --pipe \
  /home/freiburger/Documents/py_venv/bin/python /home/freiburger/Documents/prFBA/pfunc_estimator.py build-catalogue
```
Sometimes the run's own selection or prefetch writes new genomes into the PGFam cache. That happens
with a `--gene-cache` equal to the shared cache, the serial builder's default. The run then warns and
records `staled_by_this_run`. Rebuild the catalogue before the next run.

With every needed genome cached, an e3 build makes no network call. A genome that cannot be fetched is
left out of the neighbour list and the prior, and is counted in `_meta.pfunc.pgfam_store`.

**Outputs.**
- `OUT/pfunc/pgfam_store.json` holds the slim PGFam sets the run used.
- `OUT/pfunc/pgfam_fetched.json` holds the genomes fetched from BV-BRC. It is never written to the
  shared cache, and a failed fetch is never stored.
- For E3 genomes, `provenance/{asv}.json` gains:
  - `probability_kn`: the replaced k/n values;
  - `pfunc`: the full diagnostics, including the clade rank, neighbour and member coverage, the anchor,
    `reps_not_built`, `anchor_built` and per-source `source_digests`.
- Each manifest entry gets a `pfunc` summary, with the reason whenever E3 was not applied.
- `_meta.pfunc` holds the estimator, its parameters, the gating rule, `code_provenance` (sha256 of
  `pfunc_estimator.py` and the builder files, git HEAD, dirty flag), the catalogue and store
  provenance, and counts.
- If E3 raises for an ASV, that genome keeps k/n with the reason `error: ...`. It is never lost.

**Input stamps (every build, `--pfunc kn` too).** They let the retro tool below tell whether a
build's inputs changed since. No per-ASV file changes.
- `_meta.selection_stamp`: `{sha256, size, mtime_ns}` of the selection audit, taken when the build
  reads it (`--selection`) or writes it (`--hits`). The parallel builder's `--hits --limit` run writes
  no audit and stamps `null`.
- `OUT/pfunc/source_digests.json`: `{genome_id: source_digest}` for every built source genome
  (`pfunc_estimator.source_digest`, a sha256 over the genome's sorted (function, PGFam) pairs),
  computed on the genome dicts as loaded, before the merge. `_meta.source_digests` records the file,
  the count and any genome whose digest differed between two ASVs of one run (stored as `null`).
- Cost: the digests add about 10% builder CPU. Measured on a serial `--pfunc kn` build of the first
  300 EmilyKin d100 ASVs: 210.8 s user time before, 232.5 s after. The per-ASV files stayed
  byte-identical 286/286.

**Applying E3 to an existing k/n build** (copy-on-write: only a NEW directory is written; the same
code as at build time):
```bash
systemd-run --user --unit=e3-retro-$(date +%s) --collect -p MemoryMax=40G -p MemorySwapMax=0 \
  -p WorkingDirectory=/home/freiburger/Documents/prFBA --wait --pipe \
  /home/freiburger/Documents/py_venv/bin/python /home/freiburger/Documents/prFBA/reweight_build_e3.py \
    --build-dir /home/freiburger/Documents/EmilyKin/synthetic_genomes_d100 \
    --out-dir   /home/freiburger/Documents/EmilyKin/synthetic_genomes_d100_e3 --procs 8   # [--limit N] [--asvs LIST]
```
The retro tool applies E3 only where the build's inputs are unchanged. Exactly these are checked,
in this order:
1. **The selection audit (whole run).** For a build with `_meta.selection_stamp`: the audit's
   sha256 now must equal the stamp. Otherwise every built ASV in scope, n = 1 included, is copied
   unchanged (k/n) with `inputs drifted since build: selection audit`, and nothing is streamed or
   fetched. For a legacy build without the stamp (every build made before the stamps were added on
   2026-11-12, including the shipped d100 builds): the run **stops before writing anything** if the
   audit's mtime is newer than the source manifest's, unless `--allow-audit-newer` (a warning,
   recorded in `_meta.pfunc.input_checks.selection_audit`). This is the only check on the audit's content
   beyond the selected genomes. A legacy build whose audit was rewritten but kept an older mtime is
   not caught.
2. **Per ASV with n >= 2**, each failure copying that ASV unchanged (k/n):
   - `selected genomes`: the audit record's selected genome ids equal the built plus failed ones;
   - every built source genome is in `--genome-cache-dir` (else `source genome(s) missing from the
     genome cache`);
   - `function_prevalence` of the cached source genomes equals the provenance's;
   - `source digest`: with the build's `pfunc/source_digests.json`, each cached source genome's digest
     equals the recorded one (this also catches a changed PGFam alias alone). A legacy build without
     that file instead fails the ASV with `source genome cache newer than build` when any of its
     cached source-genome files has an mtime newer than the source manifest's, unless
     `--allow-cache-newer`. A scratch copy of the cache made without preserving mtimes trips this,
     so check such a copy and then pass the flag.
   The reason is `inputs drifted since build: <which>`, and the ids involved go to the entry's
   `pfunc.drift_genomes`.

Other cases:
- A record missing from the audit gets `record not in selection audit`.
- An ASV that raises is copied unchanged with `error: ...`.
- After every ASV, its `{asv}.json`, `genome_objects/` and `provenance/` outputs are parsed. If one is
  missing or unparseable (a source build that was already broken), the entry becomes `failed`, with
  the reason in `error`. The exit code is then 1, as it is when an ASV raises.

The tool does not write `mapping_provenance.json` or the abstention reports. Re-run
`write_mapping_provenance.py` on the new directory.

**Memory.** Measured on the full EmilyKin d100 retro: 1,410 ASVs, 20,697 genomes, store file 0.69 GB,
354 s at `--procs 8`. That was before the per-ASV output parse was added. The parse costs about 94 ms
per built ASV (measured serially on 300 EmilyKin ASVs), about 6 CPU-minutes, or roughly 45 s at
`--procs 8`, for all 3,755 built EmilyKin genomes. The parent's RSS was 3.1 GB and each worker's PSS 0.6-0.8 GB; the workers share
the store copy-on-write through fork. The cgroup's `memory.peak` reached a 40G MemoryMax through page
cache (the 5 GB cache stream and the written genomes), not process memory. Run full builds and retros
under `systemd-run --user ... -p MemoryMax=40G -p MemorySwapMax=0`.

A build cannot know before the merge which references will build, so its store also holds the clade
members of every selected reference. That is about 13% more genomes than the retro tool's exact scope
(EmilyKin 23,709 vs 21,003; codiffusion 23,262 vs 20,774), all already in the cache. An e3 `--hits`
build therefore adds about 3-3.5 GB for the store to the parent; that full path was not measured at
deployment scale. The NCBI taxdump (several GB) is loaded only when an anchor is outside the catalogue.

**KBase.** `upload_kbase_genomes.py` records the estimator for E3 builds only, as provenance
`method_params[0].feature_probability` (the `_meta.pfunc` stamp) and ws metadata
`probability_estimator`. k/n uploads are unchanged. Objects are named by ASV and an existing name is
skipped, so an E3 rebuild of an already uploaded build needs a **new workspace or new object names**.

**The tuned objective.** `prfba_objective.run_prfba(model, environment)` is prFBA's recommended prFBA
call, with the Stage O configuration `est/lin/G0/X1`:
- `TUNED_PRFBA_KWARGS`: ex_weight=1, min_prob=0.01, prob_exp=1, missing_probability=0.0, biomass bio1.
- `set_reaction_probabilities(model, feature_probabilities)` puts the genome's p onto the reactions.

It needs the ModelSEEDpy branch `fix-prfba-single-model` (`PYTHONPATH=~/Documents/ModelSEEDpy-prfba`).
The released library raises a clear ImportError.

Caveat on ex_weight=1: an import then costs the same as a gene-less reaction, so total exchange flux
rises, by +2.50 [+1.20, +3.39] units against a median of about 165. Which nutrients are taken up more
was not recorded, so check it against measured exchange data before relying on it. `tuned=False` or
`ex_weight=100` keeps the shipped exchange term. ModelSEEDpy's library default (ex_weight=100) is
unchanged, because other callers (community studies) use it.

Downstream users of p were not re-run with E3 (codiffusion's community gapfill penalty 2 - p,
ek_pipeline evidence filters). `METHODS.md` covers the embedding space and does not describe p(function).

## Quick start

```bash
# reference a new sequence (or a FASTA of many) against the space
python query.py --fasta new.fasta --topk 10
python query.py --fasta new.fasta --topk 10 --min-len 1400 --max-len 1600   # full-length refs only

# re-run / extend the analysis (full space, or a length-restricted view)
PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True python analyze_space.py --store .
PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True python analyze_space.py --store . \
    --min-len 1400 --max-len 1600 --tag fulllen      # writes *_fulllen.* alongside the full run
```

Environment: uv venv at `~/Documents/py_venv`, torch 2.11+cu130, **transformers 4.44.2** (pinned —
NT-v2's gated FFN is incompatible with transformers ≥ 5).

## Large files (> 100 MB, git-ignored)

These exceed GitHub's 100 MB limit and are **git-ignored** (regenerable). Expected to be absent from a clone:

| file | size | regenerate |
|---|---|---|
| `embeddings.f16.npy` | 900 MB | `python embed_16s.py` |
| `v4v5_refs.json` | 534 MB | `python insilico_pcr.py` |
| `names.dmp` / `nodes.dmp` | 277 MB / 206 MB | NCBI taxdump auto-downloaded by `taxopy` (concordance.py) |
| `v4v5_store/embeddings.f16.npy` | 192 MB | `python embed_16s.py --fasta v4v5_refs.json --outdir v4v5_store` |
| `old_asv_comparison/asv_top100_hits.json` | ~100 MB | `python hit_amplicons.py --topk 100` |
| `bvbrc_cache/genome_gene_families.json` | ~1.5 GB | per-genome PGFam cache, rebuilt by `prefetch_gene_families.py` (whole `bvbrc_cache/` is git-ignored) |
| `region_validation/data/hits/asv_top20_alignment_hits.json` | 774 MB | `region_validation/scripts/run_validation.sh` (align stage; fanned md5→genome top-20 hits) |
| `region_validation/data/hits/raw/asv_top20_alignment_hits.json` | 690 MB | `region_validation/scripts/run_validation.sh` (align stage; raw md5-keyed top-20 hits) |
| `region_validation/data/selection.json` | 666 MB | `region_validation/scripts/run_validation.sh` (select_references over 67,991 amplicons) |
| `region_validation/data/taxacc_per_record.jsonl` | 270 MB | `region_validation/scripts/run_validation.sh` (score_taxacc per-record output) |

(Threshold = GitHub's 100 MB hard limit; no file is near 100 GB. The result-summary JSONs,
figure-input CSVs, figures, and the manuscript under `region_validation/` are all < 100 MB and **are** tracked.)
