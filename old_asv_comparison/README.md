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
| `comparison_stats.json` / `old_vs_embedding.csv` | Part-1 aggregate + per-ASV |
| `part2_results.json` | Part-2 metrics |
| `asv_top20_hits.json` / `asv_summary.csv` / `findings_stats.json` | old-ASV embedding hits + match stats |

## Reproduce
```bash
# from prFBA/ : embed the old ASVs, then compare + validate
python hit_amplicons.py --amplicons ../codiffusion_bioreactor/model_inputs/dna-sequences.fasta \
    --store v4v5_store --naive-store . --topk 20 --outdir old_asv_comparison
cd old_asv_comparison
python compare_old_mappings.py
python part2_cosine_vs_identity.py
```
