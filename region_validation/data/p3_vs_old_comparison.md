# Old (P2) vs P3 (gene-richest-conspecific) synthetic-genome comparison

P3 keeps the gene-richest assembly as each species' representative; the anchor (taxonomic
call) is held stable by a deterministic genome_id tie-break. Layout: `synthetic_genomes/`
= old, `synthetic_genomes_p3/` = P3.

## EmilyKin (3950 ASVs)
- genome SET changed (gene-richer within-species rep): **494**
- anchor genome_id changed (same-species, gene-richer assembly): 306
- anchor taxonomy changed -> species: 6, genus: 2, family: 2  (stable by design)
- mean synthetic-genome features: old 3196 -> P3 3202  (+6); larger in P3: 316 ASVs, smaller: 172

## codiffusion (3276 ASVs)
- genome SET changed (gene-richer within-species rep): **451**
- anchor genome_id changed (same-species, gene-richer assembly): 237
- anchor taxonomy changed -> species: 12, genus: 2, family: 0  (stable by design)
- mean synthetic-genome features: old 2928 -> P3 2928  (-0); larger in P3: 267 ASVs, smaller: 175

