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

> **Data files are git-ignored.** `embeddings.f16.npy` (≈0.9 GB), `*.parquet`, and the `umap_*.npy`
> exceed GitHub's file limits, so they are not committed. They are fully regenerable —
> `python embed_16s.py` rebuilds the store, `python analyze_space.py` the analysis. The committed
> `manifest.json`, `density_stats.json`, and `fig_*.png` capture the results.

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
