# Two-panel comparison: k-means vs taxonomy-stratified 10,000-organism panels

All numbers from this session. Where a number supersedes one reported earlier in the
session, the superseded value is given too, with the reason.

## 1. What was compared

| | original panel | k-means panel |
|---|---|---|
| selection | taxonomy-stratified, sqrt-flattened availability, rare phyla and archaea up-weighted, species-dereplicated | k-means++ (k = 10,000, seed 1729, 15 iters) over L2-normalised 1024-d nucleotide-transformer embeddings; the organism nearest each centroid |
| organisms | 10,000 | 10,000 |
| distinct species | 10,000 (one representative per species) | 4,935 (NOT dereplicated) |
| phyla | 206 | 142 |
| 16S arm | in-silico PCR, 8 regions, top-20 alignment hits, self-exclusion at 16S identity tau = 0.987 | identical pipeline, isolated tree (`region_validation_kmeans`) |
| MAG arm | synthetic MAGs, skani vs the 10,000-genome panel, self-exclusion at ANI tau = 95 | identical |

Both MAG arms are compared at skani retention depth `-n 500` and contamination 0.

## 2. ANI depth: the original panel at -n 100 vs -n 500

The original panel's main MAG run used `-n 100`; the k-means run used `-n 500`. The
original was therefore re-run at `-n 500` so the cross-panel comparison is depth-matched.
Depth barely matters (contamination 0, 49,985 cells per regime):

| regime | resolve n100 | resolve n500 | cond. F1 n100 | cond. F1 n500 |
|---|---|---|---|---|
| include | 1.0000 | 1.0000 | 0.9998 | 0.9998 |
| exclude | 0.3071 | 0.3089 | 0.8066 | 0.8061 |

A five-fold deeper retention moves the self-exclusion resolve rate by +0.0018 and
conditional F1 by -0.0005. This confirms the depth-probe finding that the truncation
residual has converged well before depth 200, and it means the manuscript's existing
`-n 100` tables need no revision.

## 3. MAG arm, cross-panel (depth 500, contamination 0, tau_ANI = 95)

Self-exclusion regime, coverage on the standardisation ladder:

| completeness | coverage crude (km / orig) | domain-std | phylum-std | cond. F1 (km / orig) |
|---|---|---|---|---|
| 1.00 | 0.3422 / 0.3086 | 0.3291 / 0.3064 | 0.3021 / 0.2958 | 0.7971 / 0.8065 |
| 0.90 | 0.3413 / 0.3089 | 0.3282 / 0.3066 | 0.3014 / 0.2960 | 0.7976 / 0.8069 |
| 0.70 | 0.3425 / 0.3079 | 0.3293 / 0.3060 | 0.3026 / 0.2957 | 0.7981 / 0.8066 |
| 0.50 | 0.3426 / 0.3089 | 0.3293 / 0.3070 | 0.3025 / 0.2963 | 0.7980 / 0.8065 |
| 0.30 | 0.3454 / 0.3101 | 0.3318 / 0.3078 | 0.3059 / 0.2972 | 0.7962 / 0.8042 |

Two results:

1. **Conditional gene-capture F1 is panel-invariant** (0.796-0.798 vs 0.804-0.807).
   Whatever the selection scheme, once whole-genome ANI finds a qualifying reference the
   metabolism it recovers is the same.
2. **The k-means panel has HIGHER MAG coverage** (+3.4 points crude). Most of that is
   taxon mix: phylum-standardized the gap falls to +0.6 points. The residual is in the
   expected direction -- the k-means panel is not species-dereplicated, so after removing
   same-species hits at ANI >= 95 it still retains sub-genus relatives that pass the
   aligned-fraction gate, whereas the dereplicated original panel is left with
   genus-level relatives near skani's detection floor.

Completeness-invariance holds on both panels from 100% down to 30%.

## 4. GTDB dense-reference arm, now over the full 2,000

The dense-reference arm had completed only 4 of 8 batches (1,000 of the intended 2,000
genomes). `gtdb_sub_ids.txt` is not randomly ordered in phylum -- the first 1,000 ids span
85 phyla, the unrun 1,000 span 123 -- so the completed half over-weighted common phyla.
The remaining four batches were run (contamination-free) and the arm recombined over all
2,000 at contamination 0:

| completeness | resolve (2,000) | resolve (biased 1,000) |
|---|---|---|
| 1.00 | 0.5560 | 0.5670 |
| 0.90 | 0.5545 | 0.5630 |
| 0.70 | 0.5540 | 0.5630 |
| 0.50 | 0.5570 | 0.5690 |
| 0.30 | 0.5500 | 0.5640 |

The manuscript's claim -- a 73-fold denser reference nearly doubles the resolve rate to
~0.56 with completeness-invariance intact -- survives at full size; the biased half had
overstated it by about 1.4 points.

## 5. Where each panel samples the embedding space (figure `F_panel_space`)

Local density is measured in the full 1024-d embedding (cosine distance to the 21st
nearest neighbour among the 210,779 clusterable sequences), NOT in the UMAP projection,
which preserves neighbourhoods but not density. The ranking is stable in k
(Spearman 0.964 between k = 21 and k = 50).

Share of each panel's on-frame picks per density decile (1 = sparsest):

| decile | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| k-means | 18.9 | 17.7 | 11.9 | 14.1 | 9.7 | 8.2 | 7.7 | 4.9 | 2.7 | 4.1 |
| original | 24.6 | 24.5 | 15.7 | 14.3 | 7.1 | 5.1 | 4.1 | 2.7 | 1.1 | 0.8 |

Both schemes over-weight the sparse part of the space relative to a population-uniform
10% per decile, which is what they are designed to do. The taxonomy-stratified panel does
it about twice as hard: it draws 49.1% of its picks from the two sparsest deciles against
36.6% for k-means, and reaches the densest decile five times less often (0.8% vs 4.1%).
The dense deciles are where a few species carry thousands of near-identical operons;
species-dereplication collapses each of them to one pick.

**An asymmetry worth stating:** 903 of the original panel's 10,000 organisms have no
position in this space at all -- 743 because their best 16S is shorter than 1,400 bp or
flagged atypical, and 160 because their lineage fails the domain+phylum gate. The k-means
panel is on-frame by construction (10,000/10,000), since the gate was applied before
clustering. Panel C is therefore plotted as shares of each panel's on-frame picks.

## 6. Corrections made to this session's own earlier numbers

- **Contamination**: the k-means MAG arm was synthesized contamination-free, so all
  cross-panel MAG comparisons are restricted to the original's contamination-0 cells.
- **GTDB arm**: see section 4 -- previously quoted from a phylum-biased half.
- **16S cross-panel coverage**: the first pass counted `pgfam_missing` cells as
  not-covered. See section 7.

## 7. 16S arm, cross-panel -- after closing a PGFam cache gap

**The gap.** `run_validation.sh` has a stage, `prefetch_universe.py`, that warms the
PGFam cache for *every* top-20 hit genome. The k-means run never executed it (there is
no `prefetch_universe.log` in its data directory); it ran only the narrower candidate
prefetch, which covers the genomes named in the as-is, include-mode selection -- 39,060
of them, all cached. But exclude mode does not reuse that selection: `score_genecap`
rebuilds each amplicon's hit record with self-hits removed and re-runs
`select_representatives` over it, which reaches deeper into the top-20 and names
genomes the include-mode selection never did. **4,954 such reps were absent from the
cache**, and every cell that selected one scored `pgfam_missing` -- counted as
not-covered, against 0 such cells in the original panel's run.

The 4,954 were fetched (0 errors, 47 genuinely gene-less) and the enumeration iterated
to a fixed point. Two checks:

- The exclude-mode rep set is **invariant to cache state** -- 27,489 reps before the
  fetch and 27,489 after. Selection does not depend on what is cached (`select_all` is
  on with coherence and contamination scoring off), so the original panel's July
  scoring is *not* stale and did not need redoing.
- Running the same enumeration on the original tree finds 27,683 reps and **0 missing**,
  confirming both panels are now scored against one cache state.

After re-scoring, `pgfam_missing` is 0 in both modes across all 8 regions.

**Corrected coverage** (self-exclusion; the pre-fix k-means figures reported earlier in
this session are given for reference, and were too low by 9-12 points). Denominators exclude `truth_empty` / `pgfam_missing` / `no_hits`, matching `score_genecap`; `pred_empty` is kept and counted as not-covered:

| region | coverage crude (km / orig) | phylum-std (km / orig) | cond. F1 (km / orig) | km crude, pre-fix |
|---|---|---|---|---|
| FullLength16S | 0.5780 / 0.8288 | 0.6243 / 0.8363 | 0.6428 / 0.6450 | 0.484 |
| V1-V2 | 0.5596 / 0.8070 | 0.5961 / 0.8143 | 0.6401 / 0.6498 | -- |
| V1-V3 | 0.5804 / 0.8242 | 0.6113 / 0.8351 | 0.6429 / 0.6493 | 0.452 |
| V3-V4 | 0.4833 / 0.7379 | 0.5438 / 0.7495 | 0.5916 / 0.5977 | -- |
| V4 | 0.4545 / 0.7062 | 0.4985 / 0.6897 | 0.5629 / 0.5850 | 0.353 |
| V4-V5 | 0.4572 / 0.7111 | 0.5007 / 0.6941 | 0.5739 / 0.5903 | -- |
| V4-V5_944R | 0.4504 / 0.7155 | 0.5049 / 0.7073 | 0.5813 / 0.5936 | -- |
| V7-V9 | 0.4747 / 0.7342 | 0.5251 / 0.7176 | 0.5788 / 0.5986 | 0.379 |

**Conditional gene-capture F1 is panel-invariant** (within 0.002-0.022 everywhere, with
the original marginally ahead). Coverage is not -- and section 8 shows why.

## 8. The 16S coverage gap is a retention-depth artifact, not a panel property

The 16S arm keeps the top 20 alignment hits per amplicon. Self-exclusion then removes
every hit with identity >= 0.987 or a matching genome/taxon id. When it removes all 20,
the query abstains **by construction**, with no evidence either way about whether a
usable reference sits at rank 21. `n_removed_self` records exactly this, per cell
(FullLength16S, self-exclusion):

| hits removed from the top 20 | k-means: share / abstain | original: share / abstain |
|---|---|---|
| 0-4 | 20.8% / 0.024 | 61.3% / 0.037 |
| 5-9 | 17.4% / 0.014 | 15.3% / 0.028 |
| 10-14 | 11.6% / 0.011 | 6.0% / 0.012 |
| 15-19 | 9.9% / 0.028 | 3.1% / 0.013 |
| **20 (record emptied)** | **40.4% / 1.000** | **14.3% / 1.000** |

Where the record survives, abstention is 1-4% on both, the original marginally more often (3.30% against 1.90% pooled; CMH p = 0.002) and both far below the emptied-record rate.
Where it is emptied, abstention is 1.0 on both, by construction. The entire cross-panel
coverage difference is therefore the difference in how often the record is emptied:
40.4% vs 14.3%. The arithmetic closes -- 1 - 0.404 = 0.596 against an observed 0.578,
and 1 - 0.143 = 0.857 against an observed 0.829, the small residual being the few
percent of surviving records that still abstain plus `pred_empty` cells.

The k-means panel is hit 2.8x harder because it is not species-dereplicated and samples
the dense core of the embedding space (section 5): a query there has far more than 20
near-identical sequences in BV-BRC, so depth 20 cannot hold a single non-self relative.

The `n_removed_self` table above is the whole of the depth evidence; the MAG arm is a
separate line of support and not a second measurement of the same mechanism (the two
arms also differ seventy-fold in reference universe, and the MAG arm's coverage edge was
attributed to surviving sub-genus relatives, not to depth). What the MAG arm does show
is that the deficit is not a property of the k-means panel: conditional on the same
same-species-mate count, MAG abstention is identical across panels -- 0.6885 (k-means)
vs 0.6883 (original) for organisms with exactly one panel-mate. Once multiplicity is
controlled the MAG arm is panel-invariant, and at depth 500 the k-means panel is in fact
slightly ahead on coverage (0.342 vs 0.309).

**Consequence for the manuscript:** the self-exclusion coverage figures are lower bounds
set by retention depth, and the bound binds unequally across panels. A cross-panel
coverage comparison at fixed depth measures the panels' redundancy structure, not their
reference quality. The defensible cross-panel statement is the conditional one:
*gene-capture F1, given that a reference is found, does not depend on the selection
scheme* -- in either arm.

## 9. Superseded numbers from earlier in this session

| quantity | reported earlier | corrected | why |
|---|---|---|---|
| k-means 16S coverage, FullLength16S crude | 0.484 | 0.5780 | 4,954 exclude-mode reps missing from the PGFam cache were counted as not-covered |
| k-means 16S coverage, V1-V3 / V4 / V7-V9 crude | 0.452 / 0.353 / 0.379 | 0.580 / 0.455 / 0.475 | same |
| GTDB dense-reference resolve rate | 0.564 (c0.3) | 0.5500 (c0.3) | computed over a phylum-biased half of the subsample; now all 2,000 |
| cause of the k-means coverage deficit | "not species-dereplicated, so self-exclusion empties top-20 records" | correct in direction but the covariate was wrong: panel-mate count does not govern the 16S arm, whose reference universe is all of BV-BRC. The governing quantity is `n_removed_self` (section 8). | |

## 10. Definitions used for the MAG arm (second correction)

The cross-panel MAG numbers above use `aggregate.py`'s exact definitions, after a first
pass got both wrong:

- A cell with `outcome == "ok"` and a **blank precision** is a `pred_empty`: a reference
  was selected but contributed no gene families. It is excluded from the answered set.
  Counting it as resolved reproduces the bug the `score_mag.py` audit fixed, and
  deflated conditional F1 to 0.797 on the original panel against the 0.807 that
  `aggregate.py` -- and the manuscript's existing table -- report.
- The coverage denominator is `outcome in {ok, abstained}`, not every cell. The original
  panel has 161 gene-less source genomes against the k-means panel's 45, so keeping them
  penalises the original about three times as hard.

The corrected script reproduces `aggregate.py` on the original panel exactly (coverage
0.3086 / 0.3089 / 0.3101 and conditional F1 0.8065 / 0.8065 / 0.8042 at completeness
1.00 / 0.50 / 0.30), which is the check that it is now measuring the same thing the
manuscript does. The k-means coverage edge survives the correction (phylum-standardised
0.302 against 0.296) and conditional F1 now differs by up to 0.010 rather than 0.003,
with the original panel marginally ahead.
