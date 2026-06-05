# A 16S rRNA embedding space from BV-BRC — methods, parameters, and usage

This directory (`../prFBA`) contains a dense vector representation ("embedding space") of every
unique 16S ribosomal-RNA gene sequence collected from BV-BRC, together with the code that
produced it, an analysis of the space's geometry, and a tool to reference new sequences against
it by cosine similarity.

> Status note: the embedding store and the density/sparsity numbers in §6 are written by the
> automated pipeline once the upstream dataset (`model_inputs/BV_BRC_16S.json`) finishes
> regenerating. §1–§5, §7 and §8 are final; §6 is populated from `density_stats.json` at the end
> of the run.

---

## 1. Inputs and provenance

| File (in `codiffusion_bioreactor/model_inputs/`) | Content |
|---|---|
| `BV_BRC_16S.json` | `{ fasta_header : nucleotide_sequence }` — 16S rRNA gene sequences exported from the BV-BRC `genome_feature` API (`annotation=PATRIC, feature_type=rRNA, product=16S`), DNA-FASTA format. |
| `16S_metadata.json` | `{ feature_id : record }` — the matching feature records (`genome_id`, `taxon_id`, `genome_name`, `na_length`, `na_sequence_md5`, …) from the same query. |

Both files are produced by `BV_BRC-Copy1.ipynb`: the API is paged in 25 000-record batches into
shard files, which are then merged (cell 12) into the two JSON files above. The header of each
sequence encodes its join key and organism, e.g.

```
fig|1505.7.rna.1| 16s_rRNA [[Clostridium] sordellii strain JGS6382 | 1505.7]
     └─ feature_id ┘                └────── organism ──────┘ └ genome_id ┘
```

The **join key** between the two files is `header.split("|")[1]` (the feature id, e.g.
`1505.7.rna.1`), which is exactly the metadata key. `genome_id` and `organism` are recovered from
the header even when a metadata record is absent, so every embedded sequence retains a taxonomic
label.

---

## 2. Embedding model and representation

**Model — `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species`** (Dalla-Torre et al.,
*Nature Methods* 2024). A 498 M-parameter, ESM-style transformer pretrained by masked-language
modelling on genomes from 850+ species (bacteria, archaea, eukaryotes), making it well matched to
the multi-domain, mostly-bacterial 16S content of BV-BRC. It is, to date, one of the strongest
general-purpose nucleotide foundation models for producing taxonomically structured sequence
embeddings.

**Tokenization.** Non-overlapping **6-mers** (≈ 6 bp/token), 4 107-token vocabulary, 2 048-token
context (≈ 12 kb). A full-length 16S gene (~1 500 bp) maps to ~260 tokens and is therefore embedded
**without truncation**; only atypical multi-kb mis-exports are clipped (see `max_tokens`).

**Sequence embedding.** We take the **final hidden layer** and **mean-pool over the real
(non-padding) tokens** using the attention mask, then **L2-normalize** the 1 024-dim vector.
Mean-pooling (rather than the CLS token) is the representation InstaDeep recommends for NT and is
the standard choice for sequence-level similarity. Because every stored vector is unit-norm,

> **cosine similarity = dot product**, and cosine distance = 1 − cosine.

**Why this exact recipe on this hardware.** The target GPU is an **RTX 5070 Ti** (Blackwell,
compute capability `sm_120`, 16 GB) with PyTorch 2.11 / CUDA 13. Two compatibility facts shaped the
recipe:

1. NT-v2 uses a **gated FFN** and ships custom modelling code. transformers ≥ 5 forces the
   in-library `Esm*` classes (plain FFN), which shape-mismatch the gated checkpoint, so the model
   is loaded under **transformers 4.44.2** with `trust_remote_code=True` (its own architecture).
2. NT-v2's attention upcasts the softmax to fp32 while values remain in storage dtype, so loading
   weights directly in bf16 raises a matmul dtype clash on this stack. We therefore keep **fp32
   master weights** and run inference inside **`torch.autocast(bf16)`**, which reconciles the
   operand dtypes and still executes the matmuls on bf16 tensor cores. Measured throughput
   ≈ 155 sequences/s at 5–11 GB GPU memory — compute-bound, i.e. the GPU is well utilized.

---

## 3. Pipeline (`embed_16s.py`)

1. **Load** `BV_BRC_16S.json` + `16S_metadata.json` (188 GB RAM host; loaded fully).
2. **Exact dedup.** 16S genes are massively redundant across genomes, so each distinct nucleotide
   string (md5 key) is embedded **once**, carrying a `multiplicity` (how many input records share
   it) and `n_genomes`. This cuts GPU work by the dedup factor and — importantly — makes the
   geometry analysis meaningful (otherwise "density" would just count identical gene copies).
   Sequences shorter than `min_len` (20 bp) are recorded but not embedded.
3. **Length-sorted batching.** Unique sequences are processed shortest→longest so each padded batch
   wastes minimal compute; results are scattered back to a stable `row_id`.
4. **Embed** with the §2 recipe; vectors stream into a float16 memmap.
5. **Resumable.** `progress.json` records how many length-sorted items are done; re-running
   continues from there (inputs + sort order are deterministic).
6. **Write** the store (§4) and a full `manifest.json` provenance record.

### Parameters

| Parameter | Value | Meaning |
|---|---|---|
| `model` | `nucleotide-transformer-v2-500m-multi-species` | encoder |
| `embed_dim` | 1024 | vector dimensionality |
| pooling | masked mean of last hidden layer | sequence → vector |
| normalization | L2 (unit norm) | cosine = dot product |
| `batch_size` | 128 | sequences per forward pass |
| `max_tokens` | 512 (~3 kb) | token cap; clips only atypical mis-exports, never real 16S |
| `min_len` | 20 bp | below this: recorded, not embedded |
| `atypical_len` | > 2000 or < 200 bp | flag for non-16S-like length |
| precision | fp32 weights + bf16 autocast | Blackwell-safe mixed precision |
| dtype stored | float16 | compact; cosine unaffected |
| transformers / torch | 4.44.2 / 2.11+cu130 | pinned for NT-v2 remote code + sm_120 |

---

## 4. The embedding store (files written here)

| File | Content |
|---|---|
| `embeddings.f16.npy` | `(N_unique, 1024)` L2-normalized **float16** — the embedding space. |
| `index.parquet` | one row per unique sequence: `row_id, md5, seq_len, multiplicity, n_genomes, feature_id, organism, genome_id, taxon_id, genome_name, atypical_len`. Row order matches the embedding matrix. |
| `members.parquet` | every input sequence → its `row_id` (`feature_id, row_id, seq_len`); reverse lookup from any BV-BRC feature to its vector. |
| `sort_order.npy` | length-sorted processing order (for resume). |
| `manifest.json` | model, parameters, counts, dedup factor, length/multiplicity stats, environment, timings. |
| `density_stats.json` + `fig_*.png` + `umap_*.npy` | geometry analysis (§6). |

Load the space in three lines:

```python
import numpy as np, pandas as pd
E   = np.load("embeddings.f16.npy", mmap_mode="r")   # (N, 1024) unit vectors
idx = pd.read_parquet("index.parquet")               # row_id-aligned metadata
# cosine(i, j) = float(E[i].astype('f4') @ E[j].astype('f4'))
```

---

## 5. Density / sparsity analysis methodology (`analyze_space.py`)

All metrics are computed by **exact** (not approximate) cosine kNN on the GPU, chunked so the score
block stays within memory at any N. We characterize the cloud on the unit sphere with:

- **Packing tightness** — cosine to the 1st nearest neighbour (NN1) and mean over the k = 10
  nearest. High → the point sits among near-identical/closely-related sequences (dense); low → it
  is isolated (sparse).
- **Local density counts** — number of neighbours within cosine radii {0.99, 0.97, 0.95, 0.90,
  0.80}. The distribution separates crowded cores from sparse frontiers.
- **Global background** — cosine of millions of random pairs: where the bulk of mass sits and how
  concentrated the space is overall.
- **Effective dimensionality** — PCA of the full covariance: variance in PC1/top-10, number of
  components for 50/90/95/99 % variance, and the **participation ratio** (Σλ)²/Σλ².
- **Intrinsic dimensionality** — the **TwoNN** estimator (Facco et al. 2017), maximum-likelihood
  form `d = n / Σ ln(r₂/r₁)` on chord distances, discarding the top 5 % of ratios (near-duplicate /
  outlier tail): the dimensionality of the manifold the data actually occupy, typically far below
  1024.
- **Taxonomic coherence** — fraction of each point's k neighbours sharing its genus: does proximity
  in the embedding correspond to taxonomy (a validity check on the space).
- **Isolation / novelty frontier** — fraction of points whose NN1 < 0.80 (no close reference).
- **UMAP** (cosine metric, sub-sampled) colored by genus and by local density, for visualization.

**Validation** (`tests/validate_analysis.py`, `tests/validate_query.py`; all pass). The analysis
internals are checked against independent references: GPU exact-kNN reproduces a CPU brute-force
cosine computation (NN1 / mean-kNN max Δ < 1e-3, radius counts within ±1); chunked scoring is
stable to fp16 matmul-tiling noise; the PCA participation ratio and 90 %-variance dimension match a
numpy float64 eigendecomposition; and TwoNN recovers the known intrinsic dimension of synthetic
S^d manifolds (d = 3, 5, 10, 20) to ~10 % and monotonically. The query tool's reported cosines
match a brute-force `qᵀE` to < 3e-3 with correct top-k ordering, exact sequences are detected by md5
(their cosine is ≈ 0.999–1.0; it is not forced to exactly 1.0 because the bulk and query passes pad
to different batch shapes), novel sequences are flagged sparse (cos < 0.8), and batch == single,
`--min-sim` filtering, and FASTA/JSON parsing all behave correctly.

---

## 6. Results — geometry of the space

Built from **1,428,909** BV-BRC 16S sequences → **460,116 unique** (3.11× exact redundancy;
724 below 20 bp skipped; 1,465 = 0.32 % flagged atypical-length). Embedding: **50.1 min** on the
RTX 5070 Ti (~153 seq/s averaged over the length-sorted run; peak 277 seq/s; GPU 100 %).
All figures are in this directory; raw numbers in `density_stats.json` / `manifest.json`.

### 6.1 The space is globally concentrated and extremely dense

16S rRNA is *the* conserved phylogenetic marker, and the embeddings reflect that: even **random
pairs** of sequences average **cosine 0.90** (median 0.937, p99 0.987). The entire 460k-point cloud
occupies a small cap of the 1023-sphere rather than spreading over it.

Locally it is denser still. The mean **nearest-neighbour cosine is 0.995 (median 0.998)**, and an
average sequence has **~192,600 neighbours within cosine 0.95** — i.e. ~42 % of the whole database
sits within 0.95 of any given point:

| neighbours within… | mean count (of 460,116) |
|---|---|
| cos ≥ 0.99 | 1,390 |
| cos ≥ 0.97 | 111,396 |
| cos ≥ 0.95 | 192,565 |
| cos ≥ 0.90 | 313,329 |
| cos ≥ 0.80 | 411,240 |

**85.1 %** of unique sequences have a near-duplicate (NN1 > 0.99) and **37.8 %** have NN1 > 0.999 —
deep redundancy (the same species sequenced across many genomes) survives even exact dedup.

### 6.2 There is essentially no sparse / novel frontier *within* the database

The flip side of that density: almost nothing is isolated. The **isolated fraction (NN1 < 0.80) is
0.0000**, and only **0.09 %** of sequences have *no* neighbour above cosine 0.95. Every 16S gene in
BV-BRC has close relatives already in BV-BRC. (Sparsity therefore matters for *new external*
queries — see §7 — not for the reference set itself.)

### 6.3 Dimensionality — a low-dim manifold, densely filled

- **PCA (effective dim):** PC1 alone explains 38 % of variance; 2 PCs reach 50 %, **44 PCs reach
  90 %**, 111 reach 95 %; **participation ratio ≈ 5.0**. So the dominant taxonomic variation lives
  in a handful of directions, with a ~tens-of-dimensions tail — far below the ambient 1024.
- **Intrinsic dim (TwoNN):** on the raw set ≈ **1.1** — *not* a manifold dimension but a direct
  signature of the near-duplicate filaments (85 % of points have a near-twin, so the finest-scale
  structure is 1-D chains). After thinning near-duplicates (keep NN1 < 0.99 → 68,482 distinct
  representatives), the **taxonomic manifold dimension rises to ≈ 24** (≈ 26 at NN1 < 0.97),
  consistent with the PCA tail. The two estimators together say: a **~24-D manifold of taxonomic
  variation, blanketed in dense near-duplicate clusters**.

### 6.4 Proximity tracks taxonomy (validity check)

The 10 nearest neighbours of a sequence share its genus **50 %** of the time overall, rising to
**67 % for full-length (1400–1600 bp)** sequences — neighbourhoods are taxonomically coherent, most
cleanly so where the full gene (all 9 variable regions) is present. Genus is a crude first-token
label, so this is a conservative floor.

### 6.5 Length stratification (answering the partial-vs-full question with data)

| length band | n | NN1 cos | isolated | genus purity |
|---|---|---|---|---|
| <800 bp | 164,486 | 0.993 | 0.000 | 0.355 |
| 800–1200 | 56,538 | 0.991 | 0.000 | 0.350 |
| 1200–1400 | 19,871 | 0.992 | 0.000 | 0.417 |
| **1400–1600 (full)** | **215,941** | **0.997** | **0.000** | **0.666** |
| ≥1600 | 3,280 | 0.989 | 0.000 | 0.315 |

Full-length sequences form the **tightest, most taxonomically-coherent core** (highest NN1, 2× the
genus purity of partials). Crucially, the partials are **not** stranded: of the 240,895 partial
(<1400 bp) sequences, **88.9 % have a full-length 16S within cosine 0.95** (median best 0.979) and
only **0.0008 %** are "orphans" (best < 0.80). So length heterogeneity *softens* taxonomic
resolution (lower purity) but does **not** fragment the space — partials embed onto the same
manifold as their full-length relatives (`fig_partial_vs_fulllength.png`). This empirically
supports keeping all lengths with flags, and optionally restricting to ≥1400 bp (`--min-len 1400`)
when maximum taxonomic resolution is wanted.

### 6.6 Figures
`fig_nn1_hist.png` (packing tightness) · `fig_local_density.png` (neighbours ≥0.95, log y) ·
`fig_pairwise_hist.png` (global background) · `fig_pca_spectrum.png` (scree + cumulative variance) ·
`fig_length_strata.png` (NN1 + isolation by length band) · `fig_partial_vs_fulllength.png` ·
`fig_umap_genus.png` / `fig_umap_density.png` (UMAP, 60k sample).

---

## 7. Referencing new sequences by cosine similarity (`query.py`)

A new 16S sequence is placed in the space by embedding it with the **identical encoder** (so it is
directly comparable), L2-normalizing, and taking the dot product against every stored vector — one
GPU matmul. The nearest references (highest cosine) are returned with their taxonomy, and an exact
md5 match is flagged.

```bash
# one sequence
python query.py --seq ACGT...        --topk 10

# many at once (FASTA or {header: seq} JSON), save a table, report only confident hits
python query.py --fasta new_seqs.fasta --topk 5 --min-sim 0.90 --out hits.csv
```

The tool returns, per query, the top-k `cosine`, whether it is an exact md5 match, and the
reference `organism / genome_name / taxon_id / seq_len / multiplicity`.

**Interpreting the score** — *calibrated to this space* (§6). Because 16S is highly conserved,
absolute cosines run high: **random sequence pairs already average 0.90**, so the informative,
taxonomy-resolving band is compressed near the top. Read scores relative to that 0.90 background,
not from zero:

| cosine to nearest reference | reading |
|---|---|
| ≥ ~0.999 | near-identical 16S — usually an exact/near-exact DB match (check the md5 flag) |
| ~0.99–0.999 | same species / very close relative |
| ~0.97–0.99 | same genus / family neighbourhood |
| ~0.90–0.97 | conserved-core background — weak/ambiguous signal, *not* a confident relative |
| < ~0.85 | genuinely divergent or non-16S — below the 5th percentile of random pairs; a true sparse/novel region |

A practical consequence: a hit at 0.92 is *not* a good match here — it is essentially background.
Treat **≥ 0.97** as the threshold for a meaningful taxonomic neighbour, and **< 0.85** as a real
novelty flag.

Programmatic use (any framework, no GPU required for the search itself):

```python
import numpy as np, pandas as pd, torch, nt_embed
E   = torch.from_numpy(np.load("embeddings.f16.npy")).cuda()      # (N,1024)
idx = pd.read_parquet("index.parquet")
tok, model = nt_embed.load_model()
q   = nt_embed.embed_batch(["ACGT..."], tok, model)               # (1,1024), unit norm
sims = (q.half() @ E.T).squeeze(0)                                # cosine to every reference
top  = torch.topk(sims.float(), 10)
print(idx.iloc[top.indices.cpu().numpy()].assign(cosine=top.values.cpu().numpy()))
```

For very large or repeated query workloads, the same unit vectors can be dropped into an ANN index
(FAISS `IndexFlatIP`, HNSW, ScaNN); inner product on these vectors *is* cosine, so no reconfiguration
is needed.

---

## 8. Reproducibility

```
GPU            NVIDIA RTX 5070 Ti (Blackwell, sm_120, 16 GB)
python         3.12  (uv venv at ~/Documents/py_venv)
torch          2.11.0+cu130
transformers   4.44.2   (pinned: NT-v2 remote code; >=5 mismatches the gated FFN)
tokenizers     0.19.1
model          InstaDeepAI/nucleotide-transformer-v2-500m-multi-species
```

End-to-end:

```bash
python embed_16s.py        # build embeddings.f16.npy + index/members/manifest
python analyze_space.py    # density_stats.json + figures
python query.py --fasta new.fasta --topk 10   # reference new sequences
```

`watch_and_run.sh` wires the first two steps to fire automatically once
`model_inputs/BV_BRC_16S.json` finishes regenerating.
