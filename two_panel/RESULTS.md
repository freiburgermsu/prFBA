# Two-panel comparison: k-means vs taxonomy-stratified 10,000-organism panels

> **Re-score 2026-11-16.** The numbers below were re-derived from the artefacts re-scored on the repaired
> BV-BRC PGFam cache (2026-11-13 repair; re-score runs R2.2/R2.3; the depth-20 16S cells of sections 7-8 were
> re-scored in R1.4). Headline changes: original-panel MAG self-exclusion coverage at depth 500, 100%
> completeness, 0.3086 -> 0.3145 (k-means 0.3422 -> 0.3424), so the k-means crude edge narrows from +3.4 to
> +2.8 points and the phylum-standardised one from +0.6 to +0.4; original-panel full-length 16S coverage
> 0.8288 -> 0.8359; gene-less source genomes 161 -> 38 (original) and 45 -> 36 (k-means); in section 8 the
> surviving-record abstention difference between the panels is no longer significant once stratified (CMH
> p 0.002 -> 0.71). Unchanged: the GTDB arm (section 4), the embedding-space analysis (section 5) and the
> `n_removed_self` strata shares. Left historical: section 9 (this session's superseded numbers, all pre-repair),
> the pre-fix k-means figures and the "47 genuinely gene-less" count in section 7, and the first-pass bug
> value in section 10. The MAG singleton-abstention pair in section 8 now uses the manuscript's recipe
> (`derived_numbers.py` `panels_mag_single_*`, depth 500): the old k-means 0.6885 matched no recipe (that
> recipe gave 0.6895).
> Added: a random-panel contig-cap caveat (section 11, decision D6).

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
| exclude | 0.3127 | 0.3144 | 0.8063 | 0.8059 |

A five-fold deeper retention moves the self-exclusion resolve rate by +0.0017 and
conditional F1 by -0.0004. This confirms the depth-probe finding that the truncation
residual has converged well before depth 200, and it means the manuscript's existing
`-n 100` tables need no revision.

## 3. MAG arm, cross-panel (depth 500, contamination 0, tau_ANI = 95)

Self-exclusion regime, coverage on the standardisation ladder:

| completeness | coverage crude (km / orig) | domain-std | phylum-std | cond. F1 (km / orig) |
|---|---|---|---|---|
| 1.00 | 0.3424 / 0.3145 | 0.3297 / 0.3107 | 0.3030 / 0.2990 | 0.7969 / 0.8063 |
| 0.90 | 0.3415 / 0.3146 | 0.3289 / 0.3108 | 0.3023 / 0.2992 | 0.7975 / 0.8067 |
| 0.70 | 0.3426 / 0.3134 | 0.3297 / 0.3101 | 0.3032 / 0.2987 | 0.7980 / 0.8065 |
| 0.50 | 0.3428 / 0.3139 | 0.3299 / 0.3105 | 0.3033 / 0.2989 | 0.7979 / 0.8061 |
| 0.30 | 0.3455 / 0.3159 | 0.3322 / 0.3121 | 0.3065 / 0.3005 | 0.7961 / 0.8040 |

Two results:

1. **Conditional gene-capture F1 is panel-invariant** (0.796-0.798 vs 0.804-0.807).
   Whatever the selection scheme, once whole-genome ANI finds a qualifying reference the
   metabolism it recovers is the same.
2. **The k-means panel has HIGHER MAG coverage** (+2.8 points crude). Most of that is
   taxon mix: phylum-standardized the gap falls to +0.4 points. The residual is in the
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

The 4,954 were fetched (0 errors, 47 genuinely gene-less *(pre-repair, historical: counted on
the cache before the 2026-11-13 repair, which found stored failed fetches among such empties)*) and the enumeration iterated
to a fixed point. Two checks:

- The exclude-mode rep set is **invariant to cache state** -- 27,489 reps before the
  fetch and 27,489 after. Selection does not depend on what is cached (`select_all` is
  on with coherence and contamination scoring off), so the original panel's July
  scoring is *not* stale and did not need redoing.
- Running the same enumeration on the original tree finds 27,683 reps and **0 missing**,
  confirming both panels are now scored against one cache state.

After re-scoring, `pgfam_missing` is 0 in both modes across all 8 regions.

**Corrected coverage** (self-exclusion; the pre-fix k-means figures reported earlier in
this session are given for reference, and were too low by 9-12 points *(pre-repair, historical:
measured against the pre-repair corrected values)*). Denominators exclude `truth_empty` / `pgfam_missing` / `no_hits`, matching `score_genecap`; `pred_empty` is kept and counted as not-covered:

| region | coverage crude (km / orig) | phylum-std (km / orig) | cond. F1 (km / orig) | km crude, pre-fix |
|---|---|---|---|---|
| FullLength16S | 0.5814 / 0.8359 | 0.6285 / 0.8444 | 0.6428 / 0.6446 | 0.484 |
| V1-V2 | 0.5627 / 0.8175 | 0.5997 / 0.8235 | 0.6399 / 0.6492 | -- |
| V1-V3 | 0.5840 / 0.8348 | 0.6154 / 0.8439 | 0.6429 / 0.6491 | 0.452 |
| V3-V4 | 0.4875 / 0.7456 | 0.5489 / 0.7571 | 0.5917 / 0.5980 | -- |
| V4 | 0.4590 / 0.7209 | 0.5042 / 0.7025 | 0.5624 / 0.5838 | 0.353 |
| V4-V5 | 0.4625 / 0.7227 | 0.5082 / 0.7042 | 0.5736 / 0.5893 | -- |
| V4-V5_944R | 0.4544 / 0.7243 | 0.5109 / 0.7139 | 0.5813 / 0.5927 | -- |
| V7-V9 | 0.4800 / 0.7485 | 0.5340 / 0.7324 | 0.5789 / 0.5977 | 0.379 |

**Conditional gene-capture F1 is panel-invariant** (within 0.002-0.021 everywhere, with
the original marginally ahead). Coverage is not -- and section 8 shows why.

## 8. The 16S coverage gap is a retention-depth artifact, not a panel property

The 16S arm keeps the top 20 alignment hits per amplicon. Self-exclusion then removes
every hit with identity >= 0.987 or a matching genome/taxon id. When it removes all 20,
the query abstains **by construction**, with no evidence either way about whether a
usable reference sits at rank 21. `n_removed_self` records exactly this, per cell
(FullLength16S, self-exclusion):

| hits removed from the top 20 | k-means: share / abstain | original: share / abstain |
|---|---|---|
| 0-4 | 20.8% / 0.024 | 61.2% / 0.029 |
| 5-9 | 17.4% / 0.014 | 15.3% / 0.015 |
| 10-14 | 11.6% / 0.011 | 6.0% / 0.007 |
| 15-19 | 9.9% / 0.028 | 3.1% / 0.009 |
| **20 (record emptied)** | **40.4% / 1.000** | **14.4% / 1.000** |

Where the record survives, abstention is 1-3% on both, the original more often in the pooled rate but not once stratified (2.41% against 1.90% pooled; CMH p = 0.71), and both far below the emptied-record rate.
Where it is emptied, abstention is 1.0 on both, by construction. The entire cross-panel
coverage difference is therefore the difference in how often the record is emptied:
40.4% vs 14.4%. The arithmetic closes -- 1 - 0.404 = 0.596 against an observed 0.581,
and 1 - 0.144 = 0.856 against an observed 0.836, the small residual being the few
percent of surviving records that still abstain plus `pred_empty` cells.

The k-means panel is hit 2.8x harder because it is not species-dereplicated and samples
the dense core of the embedding space (section 5): a query there has far more than 20
near-identical sequences in BV-BRC, so depth 20 cannot hold a single non-self relative.

The `n_removed_self` table above is the whole of the depth evidence; the MAG arm is a
separate line of support and not a second measurement of the same mechanism (the two
arms also differ seventy-fold in reference universe, and the MAG arm's coverage edge was
attributed to surviving sub-genus relatives, not to depth). What the MAG arm does show
is that the deficit is not a property of the k-means panel: conditional on the same
same-species-mate count, MAG abstention is the same within sampling error across panels --
0.6889 (k-means) vs 0.6823 (original), difference +0.007 (95% CI -0.011 to +0.024), for
organisms with exactly one panel-mate. Once multiplicity is
controlled the MAG arm is panel-invariant, and at depth 500 the k-means panel is in fact
slightly ahead on coverage (0.342 vs 0.314).

**Consequence for the manuscript:** the self-exclusion coverage figures are lower bounds
set by retention depth, and the bound binds unequally across panels. A cross-panel
coverage comparison at fixed depth measures the panels' redundancy structure, not their
reference quality. The defensible cross-panel statement is the conditional one:
*gene-capture F1, given that a reference is found, does not depend on the selection
scheme* -- in either arm.

## 9. Superseded numbers from earlier in this session

*(Session record, pre-repair, historical: the "corrected" column holds the values as corrected at the
time, on the cache before the 2026-11-13 repair. The current values are in sections 4 and 7.)*

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
  `aggregate.py` -- and the manuscript's existing table -- report *(pre-repair, historical:
  the repair cut the original panel's `pred_empty` cells from 30-39 to 8-9 per completeness
  tier, so the first-pass bug would now deflate F1 far less)*.
- The coverage denominator is `outcome in {ok, abstained}`, not every cell. The original
  panel has 38 gene-less source genomes against the k-means panel's 36, so keeping them
  would now penalise the two about equally (before the 2026-11-13 cache repair, which
  showed most of the original's empties to be stored failed fetches, it was 161 against
  45, about three times as hard on the original).

The corrected script reproduces `aggregate.py` on the original panel exactly (coverage
0.3145 / 0.3139 / 0.3159 and conditional F1 0.8063 / 0.8061 / 0.8040 at completeness
1.00 / 0.50 / 0.30), which is the check that it is now measuring the same thing the
manuscript does. The k-means coverage edge survives the correction (phylum-standardised
0.303 against 0.299) and conditional F1 now differs by up to 0.010 rather than 0.003,
with the original panel marginally ahead.

## 11. Random panel: contig-fetch cap (caveat added 2026-11-16, decision D6)

No random-panel note exists elsewhere (`mag_validation_rnd/` holds none), so the caveat is recorded here.
The BV-BRC contig fetch behind the MAG arms sent one request per genome, and the server caps a response at
25,000 rows. Two uniform-random-panel genomes, 1733.9197 and 1733.9535, have exactly 25,000 contigs in
`mag_validation_rnd/contigs/` (as do two sludge genomes), so their contig sets may be truncated
(`fetch_contigs_rnd_full.py`, paged since `d2a37a4`, records the two as truncated by that cap) and their
synthetic MAGs drawn from an incomplete genome. This touches 10 of the random panel's 49,990 MAG queries
(2 genomes x 5 completeness tiers). Decision D6 of the 2026-11 re-score: a caveat, not a re-fetch (a paged
re-fetch, re-sketch and re-search of the 20 random-panel batches would take about 35 min); no random-panel
number was regenerated for it.
