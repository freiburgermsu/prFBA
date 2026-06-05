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
