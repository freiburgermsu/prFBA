# Old ASV→genome mapping vs the embedding mapping (+ cosine vs alternatives)

Two questions, answered in **`REPORT_old_vs_embedding.md`**:

1. **How well does the NT-v2 embedding reproduce the codiffusion study's prior `ASV_genomeIDs.json`
   mapping?** The same V4–V5 ASVs (`codiffusion_bioreactor/model_inputs/dna-sequences.fasta`, 3,215;
   1,736 redundancy-collapsed representatives) were re-mapped by cosine against `prFBA/v4v5_store`.
   → top-20 agreement **0.937 taxon / 0.895 genome / 0.90 genus**; top-1 much lower (0.71 genome)
   because V4–V5 ties con-generic genomes; genome Jaccard ~0.05 (embedding returns a far broader set
   than the old curated median of 2). Reproduces the *neighborhood*, not single-genome resolution.

2. **Is cosine the best similarity method, or is there something better at scale?** Empirically,
   cosine is a near-perfect **retriever** (the true highest-%identity reference is in the cosine
   top-20 for **98.3%** of ASVs, at 0.06 ms/query — ~6×10⁶× faster than full-DB alignment) but a
   **coarse re-ranker** (cosine-#1 = identity-#1 only 54.7% of the time). Recommendation: **hybrid
   retrieve-then-rerank** — cosine ANN (FAISS/ScaNN) for high-recall candidates, then alignment
   (MMseqs2/VSEARCH) to rank the top-k by true identity.

## Files
| | |
|---|---|
| `REPORT_old_vs_embedding.md` | full report + adversarial-verification appendix |
| `compare_old_mappings.py` | Part 1: embedding hits vs `ASV_genomeIDs.json` (genome/taxon/genus, top-1/top-20) |
| `part2_cosine_vs_identity.py` | Part 2: cosine vs true %identity (`Bio.Align`) + speed benchmark |
| `expand_comparison.py` | top-k agreement curve (k=1→100) + top-5 identity comparison (old vs embedding) |
| `comparison_stats.json` / `old_vs_embedding.csv` | Part-1 aggregate + per-ASV |
| `comparison_stats_topk.json` / `old_vs_embedding_topk.csv` | top-k capture at k=1,5,20,50,100 |
| `identity_top5_comparison.json` / `identity_top5_per_asv.csv` | top-5 %identity: old all-vs-all vs cosine |
| `asv_top100_genome_mapping.json` | acquired top-100 mapping: per ASV → 100 × [genome_id, cosine] |
| `part2_results.json` | Part-2 metrics |
| `asv_top20_hits.json` / `asv_summary.csv` / `findings_stats.json` | old-ASV embedding hits + match stats |

**Top-k capture:** at the top-100 cosine hits, **93.3% genome / 96.1% taxon / 95.6% genus** of ASVs
recover their prior mapping (71/77/65% at top-1). **Top-5 identity:** the embedding's top-5 are as
sequence-identical as the old all-vs-all method's (best-hit 97.5% for both; embedding ≥ old for 80%).

**Beyond NT-v2 cosine** — `FRAMEWORKS_beyond_cosine.md` (web-researched + fact-checked) surveys
alternative frameworks for the fine identity-ordering problem: the recommended fix is **retrieve
(cosine) → rerank with exhaustive alignment** (`vsearch usearch_global --maxaccepts 0 --maxrejects 0`);
alternatives include learned %identity rerankers (Identity/FASTCAR), contrastive metric-learning
(Scorpio/DNABERT-S recipe), and — for taxonomy/novelty — phylogenetic placement (DEPP/C-DEPP, EPA-ng)
and calibrated classifiers (IDTAXA). No published DNA embedding is shown to beat alignment at fine 16S
%identity ranking. Includes a decision table and a ranked "what to prototype next".

## Reproduce
```bash
# from prFBA/ : embed the old ASVs, then compare + validate
python hit_amplicons.py --amplicons ../codiffusion_bioreactor/model_inputs/dna-sequences.fasta \
    --store v4v5_store --naive-store . --topk 20 --outdir old_asv_comparison
cd old_asv_comparison
python compare_old_mappings.py
python part2_cosine_vs_identity.py
```
