# Option 2 components + feasibility of a hierarchical cluster-descend search

**Setting.** Closest-hit 16S V4-V5 mapping: ASVs against 97,624 region-matched BV-BRC inserts, using NT-v2 mean-pooled embedding cosine. Measured behavior of that embedding: it is a strong **recaller** but a coarse **re-ranker** — cosine-#1 equals alignment-identity-#1 only **54.7%** of the time; the true single best-%identity reference is in cosine top-100 / 200 / 500 only **84.5% / 87% / 88.5%** of the time; brute cosine is **0.06 ms/query** over 97k refs, so flat search is already instant at this scale.

---

## Q1 — What Option 2 is made of: ML vs classical

Option 2 was: **k-mer-sparse retrieval → Reciprocal Rank Fusion with the embedding → optional learned %identity reranker → exact alignment.**

**Net answer up front.** Only **two** components are machine learning, and only **one is strictly required**. The NT-v2 embedding (the recaller you already have) is the only *required* ML piece. The learned-%identity reranker (Identity/FASTCAR) is the *only other* ML piece, and it is optional. Everything else — k-mer sparse retrieval, RRF, exact alignment — is classical, deterministic, training-free information-retrieval / alignment-free sequence math. And in this near-identical V4-V5 regime the k-mer retriever is not merely a "complement" to the embedding: it can be the **primary** retriever, because k-mer overlap tracks single-nucleotide %identity *better* than a mean-pooled embedding does (mechanism below).

### (a) k-mer sparse retrieval — NOT ML; classical IR / alignment-free (a technique, several tools)

No training, no learned parameters, no objective minimized — deterministic feature counting plus a fixed similarity measure. Two classical flavors:

- **k-mer count / TF-IDF sparse vectors.** Decompose each sequence into its multiset of length-k substrings, build a sparse vector (counts, presence/absence, or TF-IDF-weighted), score candidates by cosine / dot-product / Jaccard. Textbook bag-of-words IR applied to DNA — a *technique*; tooling is generic (any sparse-vector / inverted-index library, or a few lines of numpy/scipy), and it is exactly the k-mer seeding/prefilter stage inside BLAST-family and USEARCH/VSEARCH/CD-HIT.
- **MinHash / FracMinHash sketches.** Store a small bottom-/scaled hash sketch and estimate Jaccard (**MinHash**, via **Mash**) or Jaccard-containment (**FracMinHash**, via **sourmash**). Still classical — hashing + set-overlap estimation, no learning. (Note: Mash uses LSH/MinHash, *not* TF-IDF; the TF-IDF sparse-vector flavor and the MinHash-sketch flavor are two distinct classical approaches, not the same thing.)

**Why k-mer overlap tracks single-nucleotide %identity better than a mean-pooled embedding here — the load-bearing point.** There is a direct, closed-form statistical link between k-mer overlap and nucleotide identity. Under a Poisson model of point substitutions, **Mash** converts MinHash Jaccard *j* into a **Mash distance** D = −(1/k)·ln(2j/(1+j)) that is a *monotone* function of the per-base mutation rate and correlates tightly with ANI/%identity — i.e. **D ≈ 1 − ANI** (a bounded-error approximation, reported RMSE 0.00274 at k=21, s=1000; it is a tight monotone correlation, *not* a literal identity). A single nucleotide substitution destroys up to *k* overlapping k-mers, so k-mer overlap is sharply and monotonically sensitive to exactly the quantity you re-rank on. A mean-pooled embedding, by contrast, averages per-position representations into one fixed vector, which (i) discards positional information and (ii) was trained for a semantic/functional objective, not for resolving 1-2 nt differences — hence a great recaller but a coarse re-ranker (the measured 54.7%). k-mer overlap has no such averaging loss; it is a near-surrogate for the alignment metric itself. That is why the k-mer stage can be **primary**, not just complementary, in this regime.

### (b) Reciprocal Rank Fusion (RRF) — NOT ML; pure arithmetic (a technique, ~one line)

No learned parameters (one fixed constant), operates only on *ranks* not scores — which is precisely why it is robust to the embedding's cosine and the k-mer score living on incommensurable scales. For a reference *d* appearing across a set of rankings *R*, with rank_r(d) its rank in ranking *r*:

> **RRFscore(d) = Σ_{r ∈ R} 1 / (k + rank_r(d))**, standard constant **k = 60**.

Origin: **Cormack, Clarke & Büttcher, "Reciprocal Rank Fusion Outperforms Condorcet and Individual Rank Learning Methods," SIGIR 2009.** (The title itself frames RRF as the *non-learned* alternative that beats learning-to-rank — direct support for "not ML.") The k=60 default comes from that paper and is now the default hybrid-search fuser in Elasticsearch/OpenSearch, Azure AI Search, Weaviate, etc. It is the deterministic glue combining the embedding ranking and the k-mer ranking.

### (c) the learned %identity reranker (Identity / FASTCAR) — YES, ML; the only other ML piece, and optional

Confirmed: **FASTCAR** (James, Luczak & Girgis, 2018) predicts the *global alignment identity score* of two DNA sequences using **self-supervised general linear models trained on a small number of alignment-free k-mer statistics** (training labels generated by mutating sequences to known identities), in linear time/space (reported ~40% faster than BLAST, 6-10× faster than USEARCH). That learned regression is what makes it ML — as opposed to the deterministic Mash formula in (a). Its role is to re-rank the fused shortlist by *predicted* %identity before paying for exact alignment: a learned approximation of the final objective.

### (d) exact alignment (final stage) — NOT ML; classical dynamic programming

Smith-Waterman / Needleman-Wunsch (or the alignment inside BLAST/USEARCH/VSEARCH) is deterministic DP. It is the ground-truth %identity you ultimately rank by. No learning.

### Tally: ML vs classical in Option 2

| Component | ML? | Required? | Type | Tool / origin |
|---|---|---|---|---|
| NT-v2 embedding (recaller) | **Yes** | **Required** (only required ML) | learned representation | your existing embedding |
| k-mer sparse retrieval | No | can be **primary** retriever | classical IR / alignment-free | TF-IDF sparse vectors; Mash (MinHash); sourmash (FracMinHash) |
| Reciprocal Rank Fusion | No | optional glue | rank arithmetic, k=60 | Cormack, Clarke & Büttcher, SIGIR 2009 |
| learned %identity reranker | **Yes** | **Optional** (only other ML) | self-supervised GLM on k-mer stats | FASTCAR (Girgis et al.) |
| exact alignment | No | required for ground truth | dynamic programming | Smith-Waterman / BLAST / USEARCH / VSEARCH |

**Bottom line.** Option 2 is mostly **classical IR and alignment-free sequence math**, not a stack of ML models. The embedding is the only *required* learned component; the FASTCAR-style reranker is the only *other* (and optional) one. The k-mer retriever and RRF are deterministic and training-free, and — because k-mer overlap is a bounded-error surrogate for single-nucleotide %identity (Mash: D ≈ 1−ANI) — the k-mer stage can serve as a **primary** retriever in this near-identical regime, not a mere complement.

---

## Q2 — The hierarchical cluster-descend scheme: feasibility & accuracy

### What the scheme actually is (it is not novel — it is a known family)

The proposal — "k-means the embeddings into ~10 clusters, pick a representative per cluster, align the query to the 10 reps, descend greedily into the best rep's cluster, sub-cluster into ~10, recurse to the leaves" — is a textbook **hierarchical k-means search tree**, the structure made famous as the **vocabulary tree** of Nistér & Stewénius (CVPR 2006): cluster all vectors into K centers, recurse k-means with K inside each partition, giving Kⁿ leaves at depth n, then route a query top-down through the centroids. With K=10 and ~97.6k refs you reach the leaves in ~5 levels (10⁵ = 100k). Mapped to the methods you asked me to cite:

- **Vocabulary tree / hierarchical k-means** (Nistér–Stewénius) — the exact data structure proposed.
- **FAISS IVF coarse quantizer + `nprobe`** — a *one-level* version: k-means centroids (the coarse quantizer) define inverted lists; at query time you probe the `nprobe` nearest centroids. "Descend into the single best rep" is **IVF with `nprobe` = 1** — the lowest-recall, fastest setting. FAISS's own guidance is that `nprobe` is the recall↔speed knob, raised to 16-64 to recover 90-95% recall; the proposal hard-codes the worst end of that knob.
- **ScaNN** — Google's k-means *partitioning tree* + anisotropic quantization + late re-ranking; same partition-then-rerank skeleton, and notably it keeps an explicit re-ranking stage rather than trusting the routing.
- **Metric trees — VP-tree / ball-tree / cover-tree / M-tree** — the *principled* "descend a tree of representatives" with branch-and-bound **backtracking** and a triangle-inequality pruning bound that lets them revisit branches. The proposal as written has **no backtracking** — it is the greedy, bound-free degenerate case, known in the ANN literature as **"defeatist search"** (descending via per-level decision boundaries without backtracking).
- **CD-HIT / UCLUST / mmseqs2 linclust** — greedy incremental **centroid clustering** for sequences; "pick a representative per cluster" is exactly their centroid abstraction, and they are routinely used to *build* such reference sets — but on alignment/k-mer identity, not embeddings.
- **Phylogenetic placement — EPA-ng / SEPP / DEPP** — the genuinely principled "descend a reference tree with backtracking" method for this domain. EPA-ng does a cheap **preplacement** to nominate a *set* of promising candidate branches, then thoroughly optimizes likelihood over that set — it never commits to a single greedy branch.

So nothing here is new; it is IVF-`nprobe`-1 / a vocabulary tree, and the right baselines are metric trees with backtracking and phylogenetic placement.

### Why it hurts in *our* dense, near-identical regime

The decisive fact: routing is done **on the embedding metric** — the same metric measured to be a strong recaller but a coarse re-ranker (cosine-#1 = identity-#1 only **54.7%**; true best in cosine top-100/200/500 only **84.5% / 87% / 88.5%**). Every failure mode below is what happens when you move that unreliable metric from a *wide flat top-k* into a *narrow sequence of hard commitments*.

1. **Greedy ("defeatist") descent is unrecoverable.** Pick the wrong rep at level 1 and the true best hit's entire subtree is never visited — there is no triangle-inequality bound to force a backtrack (unlike a real VP/ball/cover tree, whose branch-and-bound *can* revisit). A single bad routing decision out of ~5 caps recall permanently, and with only 54.7% top-1 agreement on the embedding metric, errors compound multiplicatively across depth.

2. **A representative is a lossy summary; the true best hit can live in a cluster whose rep does NOT win.** "Best mean-pooled centroid" ≠ "contains the best single-nucleotide match." Mean-pooling already discards the per-position information %identity is computed from; averaging a cluster into one rep discards it *again* (hard quantization "inevitably incurs severe quantization loss"). The query can be closest to centroid A while its true best reference sits in cluster B whose centroid ranks 2nd or 3rd — the routing analogue of the 54.7% miss, now *worse* because the rep is an average of many sequences.

3. **Near-identical sequences are split across cluster boundaries.** In a dense regime references differ by a handful of nucleotides; k-means cuts straight through these tight neighborhoods, so a query's top few true hits routinely land in *different* leaves (points near a partition boundary get separated from their neighbors). Greedy descent reaches at most one. This boundary-split error is most severe exactly when points are densely packed — our case.

4. **A tree prunes harder than flat top-k, so recall can only go down.** Brute cosine examines all ~97.6k refs and you keep a *wide* net (top-100/200/500) to reach 84.5-88.5% recall. Greedy descent examines only ~K·depth ≈ 10·5 = **~50 candidates total** and commits to a single leaf — a far narrower funnel than even top-100 flat. Since the routing metric is the *same* unreliable cosine, narrowing the candidate set cannot raise recall above flat cosine; it can only lower it. To claw recall back you must widen the beam / add backtracking (≈ raising `nprobe`, multi-branch descent), and as you widen toward "examine as many candidates as flat top-k," you erase the speedup and converge back to flat search.

5. **It solves a non-problem.** Trees exist to avoid scanning everything. But at 97.6k refs brute cosine is **0.06 ms/query** — flat search is already instant. There is no latency budget to recover, so you trade a measurable recall loss for an unmeasurable speed gain. The hierarchical scheme is the wrong tool for a problem that is not speed-bound.

### Verdict

For this regime the hierarchical cluster-descend scheme is a **known method** (vocabulary tree / hierarchical k-means = IVF with `nprobe`=1, run as greedy "defeatist search") and an **accuracy regression**: because it routes on the same coarse embedding metric but through a far narrower, unrecoverable funnel, it can only **lower recall versus flat cosine**, while delivering a speedup you don't need (flat search is already 0.06 ms/query). It solves a non-problem at 97k-ref scale. The principled "descend a reference tree with backtracking" version of this idea is **phylogenetic placement** (EPA-ng/SEPP/DEPP) or a tree built on **alignment/k-mer distance with a wide beam / backtracking** — not a greedy descent over an embedding k-means tree.

`**Measured on our data** (greedy descent on the 97,624 embeddings, K=10, representative = member nearest each centroid; true #1 by Biopython %identity over the whole DB; n=150 ASVs):

| metric | result |
|---|---|
| greedy 1-level reaches the true-#1's cluster | **19.3%** |
| beam-2 (keep top-2 representatives' clusters) | 45.3% |
| beam-3 | 62.0% |
| greedy **2-level** (one recursion) | **8.7%** |
| *flat cosine recall@500 of the true #1 (reference)* | *88.5%* |

Greedy descent reaches the right cluster only **19.3%** of the time (barely above the 10% chance baseline) and, after one recursion, retains the true #1 only **8.7%** of the time — versus **88.5%** for flat cosine top-500. With ~9,762 members per cluster, a single representative cannot indicate whether its cluster holds the query's near-identical best hit; greedy descent then discards 9 of 10 clusters irrecoverably and recursion compounds the loss. This confirms the analysis: the tree prunes *harder* than flat top-k on a metric that is unreliable for fine ordering, so it lowers the recall you cannot afford — while solving a non-problem (flat search is 0.06 ms/query at this scale).`

(The empirical greedy-routing test running separately is expected to show **greedy-descent recall < flat-cosine recall**, consistent with the ANN literature on defeatist search and with failure mode 4 above.)

### When to use it / when not

**Use a cluster-descend tree when:**
- You actually have a **latency/scale bottleneck** (millions+ of refs where flat search is too slow) — not the case at 97k.
- You route on a **reliable metric**: build the tree from **alignment or k-mer/MinHash distance**, not mean-pooled cosine. k-mer overlap tracks single-nucleotide %identity by construction (Mash converts MinHash Jaccard to a bounded-error Mash distance ≈ 1−ANI; every substitution knocks out up to *k* k-mers), whereas a mean-pooled vector has averaged that signal away.
- The clusters correspond to **real, well-separated taxa** (coarse taxonomic routing), where a single mis-route doesn't cost the best hit — not for resolving among near-identical references at the leaves.
- You **never descend greedily**: use backtracking / a wide beam (a real VP/ball/cover metric tree with branch-and-bound, IVF with `nprobe` ≫ 1, or beam search over the top-b reps per level) — accepting that this restores recall at the cost of the speedup.

**Do NOT use it when (our case):**
- Flat search is already instant (0.06 ms/query over 97k) — there is no speed problem to solve.
- Routing would be done on the **embedding metric**, whose 54.7% top-1 / 84.5-88.5% top-k unreliability is amplified by hard greedy commitments.
- References are **dense and near-identical**, so cluster boundaries split true neighbors and a single lossy representative cannot summarize them.
- The correct architecture is the **flat-recall-then-rerank** pipeline of Option 2: keep the wide flat-cosine net (optionally RRF-fused with k-mer/MinHash to lift recall), then spend compute on an **exact-alignment re-rank** of the top-k — never on a greedy hierarchical router that can only throw recall away.

---

### Sources

- Cormack, Clarke & Büttcher, *Reciprocal Rank Fusion Outperforms Condorcet and Individual Rank Learning Methods*, SIGIR 2009 — https://cormack.uwaterloo.ca/cormacksigir09-rrf.pdf
- RRF formula / k=60 (Elasticsearch docs) — https://www.elastic.co/docs/reference/elasticsearch/rest-apis/reciprocal-rank-fusion
- Ondov et al., *Mash: fast genome and metagenome distance estimation using MinHash*, Genome Biology 2016 (Mash distance ≈ 1−ANI) — https://pmc.ncbi.nlm.nih.gov/articles/PMC4915045/
- sourmash / FracMinHash ANI — https://sourmash.readthedocs.io/en/latest/classifying-signatures.html
- James, Luczak & Girgis, *FASTCAR: alignment-free prediction of sequence alignment identity scores*, bioRxiv 2018 — https://www.biorxiv.org/content/10.1101/380824v2.full
- Nistér & Stewénius, *Scalable Recognition with a Vocabulary Tree*, CVPR 2006 — https://web.ece.ucsb.edu/~manj/ece181b-Winter2010/lectures/12%20bagofwords.pdf
- FAISS indexes wiki (IVF coarse quantizer, `nprobe`) — https://github.com/facebookresearch/faiss/wiki/Faiss-indexes ; Douze et al., *The Faiss Library*, 2024 — https://arxiv.org/abs/2401.08281
- ScaNN algorithms (Google Research) — https://github.com/google-research/google-research/blob/master/scann/docs/algorithms.md
- Bruch, *Foundations of Vector Retrieval*, arXiv 2401.09350 (defeatist search; boundary-split & quantization loss) — https://arxiv.org/pdf/2401.09350
- Metric trees / branch-and-bound / cover trees — https://en.wikipedia.org/wiki/Nearest_neighbor_search ; Beygelzimer et al., *Cover Trees for Nearest Neighbor* — https://www.cs.princeton.edu/courses/archive/spr05/cos598E/bib/covertree.pdf
- Steinegger & Söding, *Linclust: clustering huge protein sequence sets in linear time*, Nat Commun 2018 (CD-HIT/UCLUST greedy centroid clustering) — https://www.nature.com/articles/s41467-018-04964-5
- Barbera et al., *EPA-ng*, Systematic Biology 2019 (preplacement → candidate-set, not greedy) — https://academic.oup.com/sysbio/article/68/2/365/5079844 ; DEPP scaling — https://pmc.ncbi.nlm.nih.gov/articles/PMC11193062/

---

# Appendix: verification

Strong confirmation across the board. The "Foundations of Vector Retrieval" (arxiv 2401.09350) and spill-tree results directly confirm: (1) defeatist search = greedy routing without backtracking, (2) points near partition boundaries get separated from their neighbors (the boundary-split failure mode), and (3) hard quantization "inevitably incurs severe quantization loss" (the lossy-representative failure mode). I have everything needed.

## Verification

I adversarially checked all four claim clusters via WebSearch/WebFetch. Every load-bearing claim is corroborated by primary or authoritative sources. Details below, with the few overstatements/caveats flagged.

**(1) k-mer/MinHash and RRF are NOT machine learning — CONFIRMED.**
- **RRF is deterministic rank arithmetic, not ML.** Confirmed it is Cormack, Clarke & Büttcher, SIGIR 2009 ([IR Anthology](https://ir.webis.de/anthology/2009.sigirconf_conference-2009.146/), [dblp](https://dblp.org/rec/conf/sigir/CormackCB09.html), [author PDF](https://cormack.uwaterloo.ca/cormacksigir09-rrf.pdf)). Formula confirmed: `RRFscore(d) = Σ_r 1/(k + rank_r(d))`, k=60 default ([Elasticsearch docs](https://www.elastic.co/docs/reference/elasticsearch/rest-apis/reciprocal-rank-fusion), [Serghei blog](https://blog.serghei.pl/posts/reciprocal-rank-fusion-explained/)). It operates on ranks with a single fixed constant — no trained parameters. Note: the paper's title literally contrasts RRF with "individual rank *learning* methods," i.e., RRF is presented as the *non-learned* alternative that beats learning-to-rank — direct support for "not ML."
- **Mash/sourmash are deterministic sketching, not ML.** Confirmed Mash = MinHash locality-sensitive hashing + a closed-form Poisson mutation model; "Mash distance strongly correlates with ANI" ([Genome Biology 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4915045/), [PubMed](https://pubmed.ncbi.nlm.nih.gov/27323842/)). Hashing + set-overlap estimation, no training. One precision caveat: the section states "**D ≈ 1−ANI** (RMSE 0.00274 at k=21, s=1000)." The RMSE figure is reported by Mash, but the relationship is **monotone**, not a literal identity `D = 1−ANI`; the exact Mash formula is D = −(1/k)·ln(2j/(1+j)) from Jaccard j. "D ≈ 1−ANI" is a fair shorthand and the source supports the tight correlation, but it is an approximation, not an equation — keep the "≈". Also note Mash uses LSH/MinHash (it does not literally use TF-IDF); the section's framing of TF-IDF sparse vectors vs. MinHash sketches as two *distinct* classical flavors is correct and does not conflate them.

**(2) Hierarchical scheme = IVF coarse quantizer / hierarchical k-means / vocabulary tree / CD-HIT centroids; EPA-ng/DEPP = principled tree descent — CONFIRMED.**
- FAISS IVF: coarse quantizer is k-means centroids; `nprobe` probes the most-promising clusters; "descend into single best rep" = nprobe=1 (lowest-recall, fastest) ([FAISS indexes wiki](https://github.com/facebookresearch/faiss/wiki/Faiss-indexes), [The Faiss Library](https://arxiv.org/html/2401.08281v4)). Confirmed.
- Vocabulary tree = recursive k-means routing (Nistér–Stewénius); the proposal is exactly this structure. Confirmed via multiple retrieval surveys.
- EPA-ng: confirmed two-phase — cheap **preplacement** nominates a candidate *set* of branches (default reduces thousands to often <10), then thorough likelihood optimization; "never commits to a single greedy branch" is accurate ([Systematic Biology 2019](https://academic.oup.com/sysbio/article/68/2/365/5079844)). DEPP confirmed as deep-learning distance-into-tree-space placement ([Syst Biol 2023](https://academic.oup.com/sysbio/article/72/1/17/6575921), [scaling paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11193062/)).

**(3) Greedy tree descent without backtracking is a known recall-loser — STRONGLY CONFIRMED.** The canonical term is **"defeatist search"**: descending the tree using decision boundaries at each level *without backtracking*. Sources confirm it trades recall for speed, and that spill-trees/backtracking exist specifically to recover the recall it loses ([Foundations of Vector Retrieval, arXiv 2401.09350](https://arxiv.org/pdf/2401.09350); [spatial-trees NN paper](https://arxiv.org/pdf/1507.03338)). This is a stronger, named confirmation than the section claimed — you may add the term "defeatist search."

**(4) "Representative is a lossy summary" and "prune-harder-than-flat-top-k → lower recall" — SOUND.**
- Lossy-summary / boundary-split: directly confirmed — "points near the decision boundary of splits can be separated from their neighbors," and hard quantization "inevitably incurs severe quantization loss" (Foundations of Vector Retrieval; spill-tree literature). Failure modes (ii) and (iii) are textbook, not speculative.
- Prune-harder → lower recall: logically sound and consistent with the nprobe=1 vs. nprobe≫1 recall knob (FAISS guidance that you raise nprobe to recover 90–95% recall). The argument that a *same-metric* narrower funnel cannot exceed flat top-k recall holds.

**(5) FASTCAR is the only-other ML piece — CONFIRMED (with a wording nuance).** Confirmed FASTCAR (James, Luczak & Girgis, 2018) uses **self-supervised general linear models** trained on **a small number of alignment-free k-mer statistics** to predict global alignment identity in linear time, with training data generated by **mutating sequences** to known identities ([bioRxiv v2](https://www.biorxiv.org/content/10.1101/380824v2.full)). All accurate.

**Overstatements to flag (minor):**
1. `D ≈ 1−ANI` — keep the "≈"; it is a tight *monotone correlation*, not the exact Mash equation. Do not present it as `D = 1−ANI`.
2. "**provably** tracks single-nucleotide %identity" (Q2) — Mash gives a closed-form *statistical estimator* with bounded RMSE, which is strong, but "provably" overstates it slightly; "by construction / with bounded error" is more defensible than "provably."
3. "FASTCAR ... ~40% faster than BLAST, 6–10x faster than USEARCH" — confirmed verbatim from the abstract, fine to keep.
4. The 54.7% / 84.5% / 87% / 88.5% and 0.06 ms/query figures are **your own internal measurements** — not web-verifiable and I did not attempt to; the analysis's *use* of them (greedy descent recall < flat cosine recall) is logically consistent with the verified ANN literature.

Net: no claim was contradicted. Strengthen by (a) adding the term "defeatist search" for greedy-no-backtrack descent, and (b) softening `D ≈ 1−ANI` and "provably" to reflect that Mash is a bounded-error monotone estimator, not an identity.
