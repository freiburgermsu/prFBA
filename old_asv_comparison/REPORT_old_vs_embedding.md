# Embedding vs prior ASV→genome mappings, and cosine vs alternatives at scale

## Executive summary

We re-mapped the same 16S V4-V5 ASVs by NT-v2 embedding cosine against a region-matched BV-BRC space (97,624 unique inserts) and asked two questions: how well does the embedding reproduce the study's prior curated `ASV_genomeIDs.json` mapping, and is cosine the right similarity method? On the first, the embedding recovers the prior mapping well at the neighborhood/taxon level (taxon top-20 = 0.937, genome top-20 = 0.895, 93.7% of ASVs concordant at top-20) but poorly at rank 1 (genome top-1 = 0.714), and genome Jaccard is tiny (mean 0.049, median 0.030) because V4-V5 cannot resolve to a single genome — the embedding returns the whole degenerate con-generic neighborhood while the old curation is narrow (median 2, mean 5.2 genomes/ASV). On the second, cosine is an excellent, massively scalable retriever — the globally highest-identity reference falls inside the cosine top-20 for 98.3% of ASVs (1/60 missed) at 0.06 ms/query, a ~6×10⁶× speedup over exhaustive pairwise alignment (357 s/query) — but only a coarse re-ranker, equalling the true identity-#1 hit just 54.7% of the time (Spearman ρ = 0.547 within the top-20; cosine-#1 averages 95.2% identity vs 96.1% for the best candidate, a 0.9-pp gap). Nothing "beats" cosine as the embedding metric (on L2-normalized vectors it is monotone with dot-product and Euclidean distance), and nothing beats alignment for true identity, so the answer is to use both in series. The recommendation is a hybrid retrieve-then-rerank workflow: cosine ANN to generate a high-recall top-k, then alignment (MMseqs2/VSEARCH) to rerank those few candidates by true % identity.

## Methods recap

The codiffusion anaerobic-digester study had assigned each V4-V5 ASV to a set of BV-BRC genomes (`modeling_files/ASV_genomeIDs.json`; ASVs with redundant genome sets were combined, collapsing 3,215 ASVs to 1,736 representatives). We re-mapped the **same** ASVs by NT-v2 embedding cosine against a region-matched V4-V5 BV-BRC embedding space of 97,624 unique inserts, taking each ASV's top-20 nearest references (and the "expanded" set of all genomes sharing a matched insert). We then compared the embedding hits against the old curated sets at three levels — genome ID, taxon (organism name), and genus — at top-1 and top-20, plus set-overlap Jaccard. Separately (Part 2), we tested whether cosine actually retrieves the highest-sequence-identity references by re-scoring candidates with true % identity via `Bio.Align` (within-candidate ranking on n=149 ASVs over their cosine top-20; a global miss-check on n=60 ASVs against a 4,000-reference random pool), and benchmarked brute-force cosine vs pairwise alignment timing.

## 1. Embedding vs the prior ASV→genome mapping

We asked how often the embedding's hit(s) recover the study's prior curated assignment (n = 1,736 representative ASVs; "expanded" = every genome sharing the matched V4-V5 insert).

**Agreement (embedding hit ∈ old set).**

| Level | top-1 | top-20 |
|---|---|---|
| Genome | 0.714 | 0.895 |
| Taxon (organism name) | 0.771 | 0.937 |
| Genus (n = 1,146 evaluable) | 0.637 | 0.900 |

At the neighborhood level the methods are nearly interchangeable: **93.7%** of ASVs have `genome_top20 OR taxon_top20 = True` (1,626/1,736). Genome-top20 (89.5%) is a strict subset of taxon-top20 here (every genome hit implies a taxon hit), so the OR equals the taxon rate exactly. Only **6.3%** (110 ASVs) are fully discordant (neither genome nor taxon recovered in the top-20).

**Why top-1 ≪ top-20.** The gap is rank degeneracy, not retrieval failure. V4-V5 is too short to separate con-generic genomes, so many references tie at cosine ≈ 1.0 and the rank-1 genome is essentially a random draw from a degenerate cluster. Concretely, **18.1%** of ASVs (314/1,736) miss at genome-top-1 yet hit at genome-top-20 — the right genome is present but not at slot 1. The drop is smallest at the taxon level (0.771→0.937, a coarser label is harder to mis-rank) and largest at genome (0.714→0.895), exactly the resolution where the ties bite.

**Why genome Jaccard is tiny (≈0.05) despite 0.90 containment.** Containment is high because the old set is usually *inside* the embedding's returned set; the sets are not the same size. The old curation is narrow (old genomes/ASV median = 2, mean = 5.25), while one matched V4-V5 insert expands to a far broader genome set, so `|old ∩ emb| / |old ∪ emb|` stays small even when old ⊂ emb. Mean Jaccard(top20) = **0.049**, median = **0.030**. Broad-insert cases make this explicit: ASVs whose insert maps to 47–53 genomes have Jaccard 0.018–0.036, e.g. the 53-genome `uncultured Gammaproteobacteria` insert (Jaccard 0.036). V4-V5 simply cannot resolve to a single genome — the embedding returns the whole degenerate neighborhood, inflating the union and deflating Jaccard.

**Concrete discordant examples** (`genome_top20 = taxon_top20 = False`, despite strong matches):
- `79ba36cf…`: old set = 10 genomes/10 taxa; embedding #1 = *Sphingobium yanoikuyae A-TP*.
- `a3d67e84…`: old set = 29 genomes/19 taxa; embedding #1 = *Chlamydiales bacterium SCN18_25_1_16_R1_B_38_15*.
- `2780b1a1…`: old set = 25 genomes/13 taxa; embedding #1 = *Pseudomonas putida PP2323*.
- `da380d24…`: old set = 21 genomes/17 taxa; embedding #1 = *Parachlamydiaceae bacterium SSA5*.

These are not low-quality matches (best cosines for fully-discordant ASVs run ≈ 0.997–0.998); they are cases where the embedding and the prior curation landed on different members/clades of an indistinguishable V4-V5 set, plus references absent from the older BV-BRC snapshot.

**Implication.** The embedding reproduces the prior mapping at the neighborhood/taxon level (93.7% top-20, 0.937 taxon, 0.895 genome) and does so as a fast first pass, but the low top-1 and tiny Jaccard show it cannot replace genome-level resolution. V4-V5 cosine yields the correct *neighborhood*; collapsing that to a single curated genome still requires identity-based scoring (Part 2).

> **Caveat — not a symmetric ground-truth benchmark.** The embedding side returns the *broad* "expanded" set of all genomes sharing a matched insert, whereas the old side is a *small, curated* set (median 2, mean 5.2). With Jaccard(top20) mean ≈0.05, the high top-1/top-20 "agreement" largely reflects that the wide embedding net catches the few curated genomes — not that two independent methods converge on the same answer. Genus agreement also rests on only 1,146 of 1,736 evaluable rows, and the old mapping is itself an imperfect 16S→genome assignment, not a gold standard.

## 2. Does cosine find the closest hits?

We tested whether NT-v2 embedding cosine actually retrieves the highest-sequence-identity references, scoring true % identity with `Bio.Align`.

**Within-candidate ranking (n = 149 ASVs, cosine top-20 re-scored by identity; "best candidate" = best within the top-20, not the global best).** Cosine and % identity are positively but only moderately correlated: Spearman ρ = **0.547**. The cosine-#1 hit equals the identity-#1 hit only **54.7%** of the time. The cost of this imperfect ordering is small: the cosine-#1 candidate averages **95.2%** identity versus **96.1%** for the best-identity candidate in the same top-20 — a gap of just **0.9 percentage points**.

**Miss check — does the best hit fall inside top-20? (n = 60 ASVs vs a 4,000-reference random pool).** The globally highest-identity reference is contained in the cosine top-20 for **98.3%** of ASVs (only **1 / 60** missed; recall@20 = 1 − 1/60). Recall@20 is therefore high even when the #1 slot is mis-ordered.

**Speed.** Brute-force cosine costs **0.06 ms/query** (0.20 s for the full 3,215 × 97,624 comparison), versus **3.7 ms/pair** for pairwise alignment = **357 s/query** against the full DB — a **~5.9 × 10⁶×** speedup (the stored value 5,864,546× is computed from un-rounded per-query times; the rounded inputs give 5.95 × 10⁶×).

| Metric | Value |
|---|---|
| Spearman(cosine, %id), within top-20 | 0.547 |
| cosine-#1 == identity-#1 | 54.7% |
| cosine-#1 mean %id vs best-candidate mean %id | 95.2% vs 96.1% (−0.9 pp) |
| top-20 recall of global highest-identity ref | 98.3% (1/60 missed) |
| cosine vs alignment per query | 0.06 ms vs 357 s (~5.9 × 10⁶×) |

**Interpretation.** Cosine is an excellent, massively scalable **retriever** — its top-20 recovers the globally best-identity reference 98.3% of the time at ~6 million-fold lower cost than exhaustive alignment — but a **coarse re-ranker**: with ρ = 0.547 it does not perfectly order candidates by identity, picking the exact best hit barely above half the time. The practical workflow is thus to use cosine for cheap top-k retrieval, then re-rank the short list by true alignment identity to recover the final ~0.9 pp. (Identity figures rest on small samples — within-candidate n=149, miss-check n=60 vs a 4,000-ref pool rather than the full 97,624-ref DB.)

## 3. Is cosine the best method, or is there something better at scale?

Short answer: cosine is the correct *retrieval* metric and is essentially free (0.06 ms/query, a 5.9×10⁶× speedup over pairwise alignment), but it is *not* a correct *ranking* metric for sequence identity — our own numbers place the true highest-identity reference at rank 1 only 54.7% of the time. The right design is not "cosine vs. alignment" but a **hybrid: ANN retrieval, then alignment re-rank**.

### (a) Cosine on L2-normalized embeddings is already the optimal embedding metric — the lever is the index, not the metric

For unit-normalized NT-v2 embeddings, cosine similarity equals the dot product and is monotone with squared Euclidean distance (‖a−b‖² = 2 − 2·a·b), so cosine, inner-product, and L2 retrieval return the identical ranking. There is no "better similarity function" to swap in; on the embedding side the only remaining question is how to avoid the brute-force `3,215 × 97,624` scan (here 0.20 s — trivial, but linear in DB size). Approximate-NN indexing makes this sublinear and billion-scalable while preserving cosine:
- **FAISS** — `IndexFlatIP`/normalized-L2 for exact; `IVF` (inverted-file coarse quantization) + **PQ** (product quantization) for memory-compressed billion-vector search; or **HNSW** (graph) for highest recall at fixed latency; GPU-accelerated. IVF-PQ is the standard memory/speed/recall balance; HNSW trades RAM for recall.
- **ScaNN** (Google) — anisotropic quantization tuned specifically for maximum-inner-product search (MIPS), typically the recall/latency leader on dense retrieval.
- **DiskANN** — graph index that spills to SSD for single-node billion-scale at low RAM.

For our 97,624-insert region-matched space these are unnecessary (brute force is already 0.2 s); they matter only if the reference space grows to the full BV-BRC (or cross-region) scale. So: keep cosine, add an index *only when the DB outgrows brute force*.

### (b) Sequence-native alternatives that return *actual* identity at scale

Embedding cosine is a learned proxy for identity (Spearman = 0.547 within our top-20). If the deliverable is a real % identity or a species/strain call, use a sequence-native method:
- **MMseqs2** — k-mer-prefiltered, vectorized Smith–Waterman; ~10,000× faster than BLAST and near-BLAST sensitivity at ~100×, scales across cores/nodes, and ships a 16S/18S/SSU-vs-SILVA taxonomy workflow. Best when you want alignment-grade identity with throughput.
- **VSEARCH / USEARCH** — k-mer prefilter + global alignment; the *de facto* 16S standard (OTU clustering, chimera detection, `usearch_global` identity assignment). Directly gives the identity threshold-based calls (e.g., 97% / 98.7% species cutoffs) that amplicon pipelines expect.
- **MinHash / FracMinHash sketching — Mash, sourmash** — estimate ANI from k-mer sketches in near-constant time per comparison; **sourmash Branchwater** searches petabyte-scale collections. Ideal for massive prefiltering / containment / ANI estimation, though for short V4–V5 amplicons the small k-mer set makes alignment (VSEARCH/MMseqs2) the more reliable identity caller.

For 16S taxonomy specifically, alignment/k-mer identity (VSEARCH-class, SILVA/MiDAS references) remains the established gold standard; embeddings are not yet a validated replacement for threshold-based species calls.

### (c) Recommended pattern for this data: hybrid retrieve-then-rerank

Our two results define the design exactly:
1. **Recall is excellent**: the globally highest-identity reference (vs. a 4,000-ref random pool) lands inside the cosine top-20 for **98.3%** of ASVs (1 miss in 60). Cosine is a near-perfect *candidate generator*.
2. **Top-1 ranking is weak**: cosine's #1 equals the identity #1 only **54.7%** of the time; cosine-#1 averages **95.2%** identity vs. **96.1%** for the best top-20 candidate (a ~0.9-point identity gap, recoverable).

So:

> **Stage 1 — retrieve (cosine ANN):** take the cosine top-k (here k=20) per ASV. Cost ≈ 0.06 ms/query; recall ≈ 98.3%.
> **Stage 2 — rerank (alignment):** align only those ~20 candidates with MMseqs2 or VSEARCH (`Bio.Align` here costs 3.7 ms/pair, so 20 pairs ≈ 0.07 s/query) and rank by true % identity.

This captures the ~6×10⁶× speed advantage of cosine (you never align against the full DB — only 20 candidates instead of 97,624) *and* alignment's precision, converting the 54.7% top-1 hit rate into ~98.3% (bounded by stage-1 recall, tunable upward with larger k). This is the canonical retrieve-then-rerank pattern and is the recommendation for production ASV mapping.

### When to prefer each (decision summary)

| Use case | Method |
|---|---|
| Fast high-recall candidate generation; billion-scale DBs | Embedding cosine + FAISS/ScaNN/DiskANN |
| Remote homology / divergent refs; length & indel robustness; downstream ML features | Embedding (cosine) |
| Exact identity value; 97% / 98.7% species or strain calls; interpretable alignments | Alignment (MMseqs2 / VSEARCH) |
| Massive prefilter / ANI / containment across petabyte collections | FracMinHash (sourmash, Mash) |
| **Our ASV→genome mapping (precision *and* scale)** | **Hybrid: cosine ANN retrieve → MMseqs2/VSEARCH rerank** |

Bottom line: nothing "beats" cosine as the embedding metric, and nothing beats alignment for true identity — so use both in series. Pure embedding is preferable only when you want speed / remote-homology / ML features and can tolerate the 54.7% top-1 imprecision; pure alignment is preferable when you need exact thresholds but can be made affordable at scale by letting cosine prune the search space first.

## Recommendation

Adopt a **two-stage retrieve-then-rerank pipeline** for ASV→genome mapping:

1. **Retrieve (cosine ANN, high recall):** L2-normalize NT-v2 embeddings and retrieve the cosine top-k. Knobs:
   - **k:** start at k=20 (98.3% recall against a 4,000-ref pool); raise k to lift recall on larger/full-DB reference spaces.
   - **Index:** brute-force `IndexFlatIP` while the reference space is ≤~10⁵ inserts (our 97,624-insert space scans in 0.20 s). Switch to an ANN index **only when the DB outgrows brute force** — FAISS **IVF-PQ** for the memory/speed/recall balance, **HNSW** for max recall at fixed latency (RAM-heavy), **ScaNN** for MIPS-optimal recall/latency, **DiskANN** for billion-scale at low RAM. All preserve the cosine ranking.
2. **Rerank (alignment, high precision):** align only the ~k retrieved candidates per ASV with **MMseqs2** or **VSEARCH** and rank by true % identity; apply the amplicon thresholds (97% / 98.7%) for species/strain calls. Cost is ~k alignments/query (≈0.07 s for k=20), never the full DB.

This delivers cosine's ~6×10⁶× speed advantage *and* alignment-grade identity, turning the 54.7% raw top-1 hit rate into ~98.3% precision-at-rank-1 bounded only by stage-1 recall.

## Outputs

| File | Contents |
|---|---|
| `comparison_stats.json` | Part-1 aggregate agreement stats (genome/taxon/genus top1/top20, Jaccard) over n=1,736 |
| `old_vs_embedding.csv` | Per-ASV comparison rows (genome/taxon/genus top1/top20, jaccard, best_hit, n_old_genomes) |
| `findings_stats.json` | Old-ASV embedding match-quality stats (median best cosine 0.998; 99.9% ≥0.97) |
| `asv_top20_hits.json` | Per-ASV cosine top-20 reference hits (inserts/genomes/taxa/cosines) |
| `part2_results.json` | Part-2 cosine-vs-identity results (Spearman, top-1 agreement, recall@20, timings/speedup) |
| `compare_old_mappings.py` | Script: Part-1 embedding vs prior-mapping comparison |
| `part2_cosine_vs_identity.py` | Script: Part-2 cosine vs true % identity (`Bio.Align`) and speed benchmark |

All files are in `/home/freiburger/Documents/prFBA/old_asv_comparison/`.

---

# Appendix: adversarial verification

All numbers verified. The stored speedup (5,864,546) differs slightly from the naive 357/0.06e-3 = 5,950,000 because the JSON value is computed from the precise (un-rounded) per-query times, but both are ~6e6x as claimed. All other figures match exactly.

## Verification

**Claim 1 — per-level agreement means (n=1,736). VERDICT: CONFIRMED.**
Recomputed directly from `old_vs_embedding.csv`, exactly matching `comparison_stats.json`:
- genome top1 = 0.7143, top20 = 0.8952
- taxon top1 = 0.7713, top20 = 0.9366
- genus top1 = 0.6370, top20 = 0.8997, **n_evaluable = 1,146** (the 590 non-evaluable rows have blank genus columns and are correctly excluded)
- (bonus) genome Jaccard(top20): mean = 0.0492, median = 0.0304 — also confirmed.

**Claim 2 — coverage union and old-genome count. VERDICT: CONFIRMED (with one note).**
- Fraction with (genome_top20 OR taxon_top20) True = **0.9366**. Note: this exactly equals taxon_top20 alone — i.e. genome_top20 is a strict subset of taxon_top20 (every genome hit implies a taxon hit), so the union adds nothing.
- median n_old_genomes = **2.0**; mean = **5.25** (≈5.2). Confirmed.

**Claim 3 — Part 2 consistency, recall@20, and speedup. VERDICT: CONFIRMED.**
- Miss-check internally consistent: 1/60 missed = 1.67% (stored missed_pct 1.7), recall@20 = 1 − 1/60 = **98.33%** (stored 98.3%). Confirmed.
- Speedup: 357 s/query ÷ 0.06 ms/query = **5,950,000×** from the rounded inputs; the JSON stores **5,864,546×**, computed from the precise un-rounded per-query times. Both are ~**6×10⁶**, so the order-of-magnitude claim is CONFIRMED. The component checks also reconcile: 0.20 s ÷ 3,215 queries = 0.062 ms/query (≈0.06); 3.7 ms/pair × 97,624 refs = 361 s (≈357, the stored value uses the precise per-pair time).
- within_candidate: Spearman = 0.547, cosine-#1 == identity-#1 = 54.7%, cosine-#1 mean %id = 95.2 vs best-candidate 96.1 (0.9-pt gap). Confirmed. **Caveat on n:** this within-candidate analysis was run on **n=149** ASVs (not 1,736), and "best_candidate" = best within the top-20 cosine set, not the global best.

**Claim 4 — single biggest caveat. VERDICT: stated.**
The comparison is **not symmetric / not a clean ground-truth benchmark**. The embedding side returns the *broad* "expanded" set of all BV-BRC genomes sharing a matched V4-V5 insert, whereas the old side is a *small, curated* set (median 2, mean 5.2 genomes/ASV). With Jaccard(top20) mean ≈0.05 (median 0.03), the high top-1/top-20 "agreement" largely reflects that the wide embedding net catches the few curated genomes, not that the two methods independently converge on the same answer. Secondary caveats: genus agreement rests on only **n_eval=1,146** of 1,736 (590 unevaluable); Part 2's identity analysis uses tiny samples (within-candidate n=149, miss-check n=60 vs a 4,000-ref random pool, not the full 97,624-ref DB); and "old↔embedding agreement" treats the old mapping as a reference even though it too was an imperfect 16S→genome assignment.
