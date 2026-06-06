# Comprehensive comparison — alignment vs 16S embedding vs k-mers (all iterations)

Task: for each experimental V4-V5 ASV, find its **closest reference** among 97,624 region-matched BV-BRC
V4-V5 inserts. Ground truth = ranking by **Biopython global %identity over the whole DB**. All accuracy
figures below are measured against that ground truth.

## 1. Master table — closest-hit accuracy & cost

*(closest-hit metrics: n = 200 ASVs vs the full-DB true ranking, unless noted)*

| Method / iteration | exact top-1¹ | recall@100 #1² | recall@500 #1 | all top-5 @100³ | all top-5 @500 | speed / query⁴ | 1-time build | ML | GPU | role |
|---|---|---|---|---|---|---|---|---|---|---|
| **Brute-force alignment** (Biopython, full DB) — *ground truth* | 100%* | 100% | 100% | 100% | 100% | **357 s** | — | no | no | exact |
| edlib edit-distance (full DB)⁵ | ~exact | ~100% | ~100% | ~100% | ~100% | **1.4 s** | — | no | no | exact-ish |
| **NT-v2 embedding cosine** — region-matched | **50.5%** | **84.5%** | **88.5%** | **37%** | **48%** | **0.06 ms** | ~50 min embed | yes | yes | retriever (coarse ranker) |
| NT-v2 cosine — naive (full-length refs)⁶ | — | — | — | — | — | 0.06 ms | ~50 min | yes | yes | retriever (worse) |
| **k-mer 8** (TF-IDF cosine) | **75.5%** | **99.0%** | **99.5%** | **96%** | **98%** | **16 ms** | 12 s | **no** | **no** | retriever **+ ranker** |
| **k-mer 12** (TF-IDF cosine) | 72.5% | 98.0% | 98.5% | 94.5% | 98% | 12 ms | 18 s | no | no | retriever + ranker |
| RRF( cosine ⊕ k-mer 8 )⁷ | 57.0% | 99.0% | 99.5% | 95% | 97.5% | 16 ms + 0.06 ms | both | yes | yes | fused retriever |
| RRF( cosine ⊕ k-mer 12 ) | 56.5% | 99.0% | 99.0% | 92.5% | 97% | 12 ms + 0.06 ms | both | yes | yes | fused retriever |

¹ **exact top-1** = the method's rank-1 reference *is* the single highest-%identity reference in the whole DB.
  100% for alignment by definition. (Cosine's *within-top-20-candidate* top-1 was 54.7%, Spearman ρ=0.55; the
  50.5% here is the stricter full-DB figure.)
² **recall@k #1** = the true #1 is somewhere in the method's top-k = the ceiling an exact alignment-rerank of
  that top-k could reach.
³ **all top-5 @k** = all five of the true top-5 (by %identity) are inside the method's top-k.
⁴ speed: cosine = GPU brute-force over 97,624 refs; k-mer = CPU sparse dot incl. full scan; alignment = full
  pairwise DP vs every ref. (Reranking only a top-k *shortlist* costs k × per-pair, e.g. Biopython ≈ 3.7 ms ×
  k, or sub-second with VSEARCH.)
⁵ edlib is a fast alignment-free *edit distance*; it was the prefilter used to compute the true ranking, so its
  standalone recall is ~100% by construction (not independently benchmarked for exact top-1).
⁶ naive full-length refs give strictly lower cosines than region-matched (median 0.992 vs 0.998; region-matched
  beats naive for 98.3% of ASVs, +0.006 mean) and reach max 0.998 (never 1.0); not separately benchmarked on the
  recall metrics.
⁷ RRF (k=60) fuses the cosine and k-mer rank lists; needs the embedding (ML + GPU). It matches k-mer's recall but
  its **top-1 is worse** — folding the weaker cosine order back in dilutes k-mer's better ordering.

**Headline:** the classical **k-mer-8** retriever (no ML, no GPU) beats the NT-v2 embedding on every closest-hit
metric — 99% vs 84.5% recall@100 of the true best hit, 96% vs 37% all-top-5, 75.5% vs 50.5% exact top-1 — because
one substitution flips up to *k* k-mers, so k-mer overlap is a near-surrogate for nucleotide identity that
mean-pooling averages away. RRF fusion adds nothing (k-mer recall is already ~99%) and costs ordering.

## 2. Hierarchical embedding tree (greedy coarse-to-fine, K=10)

*(different metric — does greedy "cluster → align reps → descend" reach the true #1's cluster; n = 150)*

| variant | reaches true-#1's cluster | alignments / query |
|---|---|---|
| greedy 1-level | 19.3% | 10 |
| beam-2 (keep top-2 reps) | 45.3% | 10 |
| beam-3 | 62.0% | 10 |
| **greedy 2-level (recurse once)** | **8.7%** | 20 |

vs flat cosine recall@500 of the true #1 = **88.5%**. Greedy descent on the embedding routes wrong ~80% of the
time (a single representative can't summarize a ~9,762-member cluster) and recursion compounds the loss —
catastrophically worse than flat search, which is already 0.06 ms/query (so the tree optimizes a non-problem).

## 3. Supporting metrics measured along the way

**3a. Embedding vs the prior ASV→genome mapping** (agreement at cosine top-k; n = 1,736 representative ASVs):

| cosine top-k | genome | taxon | genus |
|---|---|---|---|
| 1 | 71.1% | 76.8% | 64.6% |
| 5 | 84.6% | 89.4% | 83.2% |
| 20 | 89.5% | 93.7% | 90.0% |
| 50 | 91.7% | 95.2% | 93.5% |
| 100 | 93.3% | 96.1% | 95.6% |

**3b. Top-5 hit sequence %identity — old all-vs-all method vs embedding** (n = 1,697):

| method (top-5) | best (mean) | best (median) | mean-of-5 |
|---|---|---|---|
| old all-vs-all (prior study) | 97.48% | 98.93% | 96.65% |
| NT-v2 embedding cosine | 97.51% | 99.19% | 95.18% |

embedding best-hit ≥ old best-hit for **80%** of ASVs (mean +0.03 pp).

**3c. Embedding match quality & cosine-vs-identity** (n = 3,215 / 200 / 149):
median best cosine **0.998**; 99.9% of ASVs match ≥0.97; cosine-vs-%identity Spearman within candidates **0.55**;
brute cosine **6×10⁶×** faster than full-DB Biopython alignment; region-match gain +0.006 over naive.

## 4. Verdict

For **closest-hit by identity** in this dense near-identical V4-V5 regime, ranked best→worst:
1. **Exact alignment** — perfect but 357 s/query un-indexed (seconds with VSEARCH/edlib indexing).
2. **k-mer-8 TF-IDF** — ~99% recall, 75.5% exact top-1, **CPU-only, no ML, 16 ms/query** → the right *retriever*
   (and a decent ranker); pair with an alignment rerank of its top-k for ~99% exact top-1.
3. **NT-v2 embedding cosine** — fastest (0.06 ms) and a fine *recaller* for coarse/taxon work, but a coarse
   identity ranker (84.5% recall@100, 50.5% top-1) — the wrong tool for fine identity ranking here.
4. **RRF fusion** — no benefit over k-mer alone (and hurts ordering).
5. **Hierarchical embedding tree** — counter-productive (lowers recall, solves a non-problem at this scale).

**Recommended production pipeline (fully classical):** `k-mer-8 TF-IDF retrieve top-k → exact alignment rerank`.
