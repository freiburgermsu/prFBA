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

3. **Biopython re-alignment** of the surviving top candidates for exact % identity / aligned
   length / reference coordinates; emit the top-20.

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
```

Every path emits `asv_top20_alignment_hits.json` + `asv_alignment_summary.csv` (identical
schema) plus `method_selection.json` (the probe, decision, and run stats). `gpu_align.py` holds
the CUDA SW kernel (also runnable standalone for an exhaustive all-reference mapping). The
backend decision is unit-tested in `tests/test_select_method.py`.

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
