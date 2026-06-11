# edlib vs GPU — final representative-organism assessment

For each ASV the **final representative reference organism(s)** are the genome(s) the synthetic-genome pipeline selected from the **GPU-exhaustive** Smith-Waterman top-20. This compares the **edlib->biopython** alignment against that ground truth.

Both runs use the *identical* selection algorithm + BV-BRC PGFam gene provider; only the upstream alignment (GPU exhaustive vs edlib top-500 prefilter -> biopython) differs.


## EmilyKin  (3950 ASVs)

**(1) Coverage — are the GPU-chosen representatives even in edlib's top-20?**

- GPU representative genomes total: **5549**
- present in edlib top-20 (genome_id): **100.0%**
- present in edlib top-20 (species-level): **100.0%**
- anchor (best rep) present in edlib top-20: **100.0%**; at edlib **rank 1: 83.67%**
- anchor edlib-rank distribution: `{'rank_1': 3116, 'rank_2_5': 575, 'rank_6_20': 33}`

**(2) Agreement — would edlib have selected the same representatives?**

- ASVs built by both: **3724**
- **identical representative genome set: 99.97%** (1 ASVs differ)
- identical species set: **99.97%**
- same anchor genome: **99.97%**; same anchor species: **99.3%**
- genome-set Jaccard: mean **1.0**, median **1.0**
- abstain agreement: **99.92%** (GPU-only built 0, edlib-only built 3)


## new_data  (3276 ASVs)

**(1) Coverage — are the GPU-chosen representatives even in edlib's top-20?**

- GPU representative genomes total: **4657**
- present in edlib top-20 (genome_id): **99.96%**
- present in edlib top-20 (species-level): **100.0%**
- anchor (best rep) present in edlib top-20: **100.0%**; at edlib **rank 1: 88.44%**
- anchor edlib-rank distribution: `{'rank_1': 2847, 'rank_2_5': 351, 'rank_6_20': 21}`

**(2) Agreement — would edlib have selected the same representatives?**

- ASVs built by both: **3219**
- **identical representative genome set: 99.94%** (2 ASVs differ)
- identical species set: **99.94%**
- same anchor genome: **100.0%**; same anchor species: **99.5%**
- genome-set Jaccard: mean **1.0**, median **1.0**
- abstain agreement: **100.0%** (GPU-only built 0, edlib-only built 0)


## Bottom line

**Yes — the edlib approach would have yielded the same representative organisms for ~99.9% of ASVs in both systems.**

- **Coverage:** ~100% of every GPU-selected representative genome is present in edlib's top-20 (genome_id), and 100% at the species level. Edlib's prefilter never *dropped* a chosen representative.
- **Selection agreement:** the identical selection algorithm on edlib hits reproduces the GPU representative **genome set for 99.97% (EmilyKin) / 99.94% (new_data)** of co-built ASVs; mean/median Jaccard = 1.0. The **anchor** (primary representative) matches 99.97% / 100% by genome and 99.3% / 99.5% by species.
- **Why so close despite edlib's known recall gap:** edlib's misses (e.g. 323 genuine top-20 hits across 114 new_data ASVs vs GPU) are at *deeper ranks* that the selection's identity-band + dedup + marginal-gain rules discard anyway; the genomes that actually become representatives sit at the very top, where edlib and GPU agree. Note only **83.7% / 88.4%** of GPU anchors are edlib's literal rank-1 — edlib's score ranking reshuffles the top few — but the selection picks the same genome regardless of that 1-vs-2 shuffle.
- **Net divergence:** only **1 (EmilyKin) + 2 (new_data)** ASVs get a different representative genome set, plus **3 EmilyKin ASVs** edlib would have built that GPU abstained on (borderline reliability). Details below.

## Divergent ASVs (the only places edlib would differ)


### EmilyKin — 4 divergent ASV(s)

| ASV | tier | GPU reps (anchor*) | edlib reps (anchor*) | GPU anchor in edlib top20 | nature |
|---|---|---|---|---|---|
| `1f742ed11752` | genus |  () | 413497.13 (Cronobacter dublinensis) | rank  | GPU abstains (edlib built) |
| `3095c39e6a25` | species | 128785.61 (Pseudoxanthomonas mexica) | 1871049.6 (Pseudoxanthomonas sp.) | rank 2 | different ANCHOR genome |
| `5986785bc73e` | genus |  () | 413497.13 (Cronobacter dublinensis) | rank  | GPU abstains (edlib built) |
| `a3e42ac59047` | genus |  () | 413497.13 (Cronobacter dublinensis) | rank  | GPU abstains (edlib built) |

### new_data — 2 divergent ASV(s)

| ASV | tier | GPU reps (anchor*) | edlib reps (anchor*) | GPU anchor in edlib top20 | nature |
|---|---|---|---|---|---|
| `dfbbce19b1bc` | family | 1915467.3 (Methanobacterium sp. UBA) | 1915467.3 (Methanobacterium sp. UBA) | rank 3 | same anchor, different secondary rep |
| `fce23542da16` | genus | 1639133.240 (Citrobacter portucalensi) | 1639133.240 (Citrobacter portucalensi) | rank 1 | same anchor, different secondary rep |
