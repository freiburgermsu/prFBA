# k-mer-augmented search — prototype results (no ML beyond the existing embedding)

We prototyped the classical pieces of "Option 2" — a **TF-IDF k-mer retriever** (sklearn char k-grams;
deterministic, training-free) and **Reciprocal Rank Fusion** (RRF, k=60) with the NT-v2 cosine ranking —
and measured them against the **true** ranking by Biopython %identity over the whole 97,624-insert DB
(edlib prefilter → Biopython rescore), n = 200 ASVs. No learned reranker (omitted per request).

## Results

| metric | cosine (NT-v2) | **k-mer 8** | k-mer 12 | RRF(cos⊕k8) | RRF(cos⊕k12) |
|---|---|---|---|---|---|
| top-1 correct (method #1 == true #1) | 50.5% | **75.5%** | 72.5% | 57.0% | 56.5% |
| recall@100 of true #1 | 84.5% | **99.0%** | 98.0% | 99.0% | 99.0% |
| recall@200 of true #1 | 87.0% | 99.5% | 98.0% | 99.0% | 99.0% |
| recall@500 of true #1 | 88.5% | 99.5% | 98.5% | 99.5% | 99.0% |
| all true top-5 present @100 | 37.0% | **96.0%** | 94.5% | 95.0% | 92.5% |
| all true top-5 present @200 | 40.0% | 96.0% | 95.5% | 96.5% | 96.0% |
| all true top-5 present @500 | 48.0% | 98.0% | 98.0% | 97.5% | 97.0% |

## What it shows

1. **The classical k-mer retriever decisively beats the NT-v2 embedding at this task — with no ML, no GPU.**
   k-mer-8 TF-IDF cosine recalls the true best hit in its top-100 for **99%** of ASVs (vs the embedding's
   84.5%), contains *all* of the true top-5 for **96%** (vs 37%), and even ranks the exact #1 first **75.5%**
   of the time (vs 50.5%). This is the predicted mechanism realized: a single nucleotide substitution flips
   up to *k* k-mers, so k-mer overlap is a near-surrogate for nucleotide identity, whereas mean-pooling
   averages that signal away.

2. **k = 8 ≥ k = 12 here.** Shorter k-mers share more between near-identical V4-V5 sequences (more robust to
   the few mismatches), giving marginally higher recall and top-1.

3. **RRF fusion helps *recall* but hurts *ordering*.** Fused recall@100 (99%) matches k-mer alone — but
   top-1 correctness *drops* to 57% (from k-mer's 75.5%) because folding the weaker cosine ranking back in
   dilutes the better k-mer order. Fusion only pays off when two retrievers are *complementary with gaps*;
   here k-mer already has ~99% recall, so there is nothing for cosine to add and its noise only costs ordering.

## Recommendation (revised by the data)

For closest-hit 16S mapping in this near-identical regime, the simplest, most accurate, **fully classical**
pipeline is:

> **k-mer-8 TF-IDF cosine retrieve (top-k) → exact alignment rerank (VSEARCH `usearch_global`).**

- The embedding is **not needed** for this task (it underperforms k-mer on every metric), and the RRF fusion
  is unnecessary (k-mer recall is already ~99%). Drop both for closest-hit.
- The k-mer retriever's recall ceiling is **~99% @100** (vs the embedding's 84.5%), so an exact rerank of its
  top-k makes the top-1 essentially exact (~99%, bounded by that recall) — without the embedding's 15-pp recall hole.
- k-mer-8 retrieval alone already gets the exact closest hit 75.5% of the time with **no alignment at all**.
- The embedding remains useful for *other* goals (remote homology, length/indel robustness, semantic
  clustering, cross-region search) — just not for fine identity-ranking against a dense near-identical DB.

Reproduce: `python kmer_augmented_search.py` → `kmer_augmented_search.json`. Classical stack only
(scikit-learn TF-IDF + numpy RRF + edlib/Biopython for the ground truth); the embedding is loaded only for
the head-to-head baseline.
