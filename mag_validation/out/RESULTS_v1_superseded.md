# MAG-matching validation on the 10,000-genome panel — results (panel arm)

Run 2026-10-08. Reference = skani sketch of the 9,997 panel genomes themselves.
99,970 synthetic MAGs (9,997 genomes x 5 completeness x 2 contamination), 2 regimes,
199,940 scored cells, 1,428,212 skani hit rows.

## Headline

Completeness is a NULL axis from 100% down to 30%, in both regimes, and 10%
contamination costs ~1 point. Under self-inclusion the MAG arm is at ceiling
(coverage 1.000, genus 1.000, gene capture 1.000) at every tier. Under
self-exclusion coverage is ~0.31 and FLAT across completeness (0.303-0.318).

Paired against the amplicon arm on the same organisms, conditional gene-capture F1:

| regime | MAG tier | vs region | n | MAG | 16S | median delta | MAG wins | 16S wins | p |
|---|---|---|---|---|---|---|---|---|---|
| exclude | 100% | FullLength16S | 1519 | 0.809 | 0.710 | +0.063 | 1212 | 303 | 9e-147 |
| exclude | 30%  | FullLength16S | 1520 | 0.805 | 0.710 | +0.062 | 1206 | 311 | 1e-140 |
| exclude | 100% | V4            | 1591 | 0.802 | 0.641 | +0.128 | 1409 | 181 | 1e-208 |
| exclude | 30%  | V4            | 1575 | 0.801 | 0.639 | +0.134 | 1396 | 178 | 9e-206 |
| include | 100% | FullLength16S | 7212 | 1.000 | 0.978 | - | 1112 | 11 | 5e-112 |
| include | 30%  | FullLength16S | 7212 | 1.000 | 0.978 | - | 1111 | 19 | 3e-110 |

**The question the proposal posed -- "the completeness at which a MAG stops beating a
full-length 16S amplicon" -- has the answer "it does not, down to 30% completeness."**
A 30%-complete MAG carrying 10% contamination recovers the same metabolism as a
complete genome, and beats a full-length amplicon by the same margin.

## Per-tier detail (panel reference)

self-INCLUDE: coverage 1.000, genus cond/pop 1.000/1.000, family 1.000/1.000,
gene capture cond/pop 1.000/1.000 at 100/90/70/50/30% completeness, 0% contamination.
At 10% contamination the 30% tier dips to genus 0.995, gene capture 0.997.

self-EXCLUDE (tau = 95 ANI; source genome, source taxon and every >=95-ANI hit removed):

| compl | cont | n | answered | genus cond/pop | family cond/pop | genecap cond/pop |
|---|---|---|---|---|---|---|
| 100% | 0%  | 9997 | 0.311 | 0.910 / 0.333 | 0.960 / 0.349 | 0.806 / 0.251 |
| 90%  | 0%  | 9997 | 0.312 | 0.910 / 0.333 | 0.957 / 0.349 | 0.807 / 0.252 |
| 70%  | 0%  | 9997 | 0.309 | 0.910 / 0.330 | 0.959 / 0.347 | 0.807 / 0.250 |
| 50%  | 0%  | 9997 | 0.310 | 0.911 / 0.331 | 0.959 / 0.348 | 0.809 / 0.251 |
| 30%  | 0%  | 9997 | 0.311 | 0.912 / 0.333 | 0.958 / 0.348 | 0.805 / 0.250 |
| 30%  | 10% | 9997 | 0.318 | 0.889 / 0.322 | 0.935 / 0.339 | 0.788 / 0.250 |

Conditional accuracy far exceeds the amplicon arm (genus 0.91 vs 0.60-0.73;
family 0.96 vs 0.79-0.86; gene capture 0.81 vs 0.68-0.72). Coverage is far lower
(0.31 vs 0.77), which is the reference-universe artifact below.

## Mechanism

* 100% of synthetic MAGs recover their own source genome among the hits (self-hit
  ANI min 98.07, median 100.00). The amplicon arm's retention-depth ceiling does
  not exist here: a whole genome has no equal-score tie clouds.
* Every self-hit is >=95 ANI, so exclusion removes it cleanly (0 of 5,000 survive).
* After exclusion the median MAG has 2 surviving hits and 28.5% have none at all.
* ANI is invariant to incompleteness because a clean contig subset still aligns
  almost entirely: AF_query stays ~95% at every tier while AF_ref ~= 0.95 x
  completeness. The selector gates on AF_query, so completeness is invisible to it.
  AF_ref is therefore a free completeness estimator the pipeline currently discards.

## What this cannot support

* **Reference-universe asymmetry.** The MAG arm searched 9,997 genomes; the amplicon
  arm searched 459,301 16S references drawn from ~112,000 genomes. The panel is
  species-dereplicated (0 of 9,997 genomes share a species; 61.7% share a genus), so
  after exclusion the nearest surviving relative is a genus-level one near skani's
  ~83 ANI detection floor. The MAG arm's 0.31 coverage and its population numbers
  (genus 0.33, gene capture 0.25) are therefore LOWER BOUNDS set by reference
  sparsity, not measurements of ANI's reach. Population-level claims in the other
  direction (16S 0.535 vs MAG 0.248 population F1) are not attributable to the
  marker. The conditional comparison is the one the design supports.
* **In-silico degradation is optimistic.** Dropping whole fragments cannot produce
  chimeric contigs, strain-mixture inflation, or coverage-biased mis-binning, which
  is how real MAGs actually fail -- the same way in-silico PCR makes the amplicon
  arm optimistic. Both arms are upper bounds, which is what keeps them comparable.
* **Contamination is spliced from out-of-panel donors** (108 representative genomes
  of a different phylum, held OUT of the reference sketch) precisely so splicing
  could not manufacture the exact ANI-100 ties a panel-internal donor produces.
  Real contamination is not a clean graft of a single foreign genome.
* **Fragmentation is synthetic.** Deposited contigs were first cut to ~50 kb mean
  because 19.8% of panel genomes are closed and 40.4% have <10 contigs, so dropping
  deposited contigs could not reach a 10%-granular grid and would have made
  degradation severity depend on each reference's assembly quality (1 to 5,265
  contigs across the panel). Contiguity therefore varies with the source assembly
  and is not itself a controlled axis.

## Comparability

Truth, metric definitions and name resolution are the amplicon arm's, not
reimplementations: truth_lineage from benchmark_truth.json (NCBI via taxopy, null at
unnamed ranks -> unscoreable, never wrong); predicted lineage re-resolved through the
same _common.lineage_for() from taxon_id = int(genome_id.split('.')[0]); names
compared after _common._norm(); three-state cells; conditional vs
population(anti-inflation) denominators; Wilson 95% CI for proportions and
normal-approx mean CI for means of proportions. PGFam truth sets are byte-identical
to the amplicon arm's (spot-checked: 3760/3288/1591 for 100053.8/1001822.3/100225.5).
25.2% of panel organisms have no named NCBI genus and are excluded from genus
denominators, matching the amplicon arm's 22-25%.

## Selection effect on the conditional comparison (checked directly)

The conditional comparison restricts to organisms BOTH arms answer. That subset is
not representative: it is the 31% of organisms for which whole-genome ANI finds a
qualifying reference, and those organisms are easier for the amplicon arm too.
Measured on the amplicon arm's own cells (exclude regime, conditional F1):

| region | 16S over all organisms it answers | over the subset the MAG arm also answers | over organisms the MAG arm abstained on |
|---|---|---|---|
| FullLength16S | 0.6450 (n=5978) | **0.7096** (n=1519) | 0.6230 (n=4459) |
| V4            | 0.5851 (n=6671) | **0.6414** (n=1591) | 0.5674 (n=5080) |

So the subset is +0.065 / +0.056 easier for 16S than its own overall conditional
performance. Because the comparison is PAIRED within that subset -- the same
organisms scored both ways -- the subset's easiness does not bias the DIFFERENCE,
and the MAG advantage of +0.063 (FullLength) and +0.128 (V4) stands as a paired
estimate. What it does bound is the population the claim covers.

**Correct statement of the headline:** on the organisms where whole-genome ANI finds
a qualifying reference, a MAG beats a 16S amplicon on gene capture by +0.06
(full-length) to +0.13 (V4) F1, and does so just as well at 30% completeness with
10% contamination as at 100%. This is NOT a claim that a MAG beats an amplicon in
general: on the other 69% of organisms the MAG arm abstains entirely while the
amplicon arm still answers at F1 0.62 (full-length) / 0.57 (V4). With this
reference universe, the marker that answers more often is the amplicon.

## Self-inclusion is near-tautological

Under self-inclusion the source genome is in the reference set, so ANI identifies it
exactly and prediction equals truth by construction -- hence coverage 1.000 and gene
capture 1.000 at every tier. The informative content of that regime is narrow: it
shows ANI recovers the right genome even from a 30%-complete, 10%-contaminated draft
(which is the completeness-invariance result), not that the method is accurate. The
include-regime comparison against the amplicon arm (1.000 vs 0.978) should be read as
"16S cannot always resolve its own source genome, ANI can", not as a performance gap.
All substantive comparison belongs to the self-exclusion regime.

## Dense-reference arm (GTDB r226) — the sparsity caveat, measured

The panel arm's 0.31 self-exclusion coverage could be reference sparsity or a real
limit of ANI. To separate them, a phylum-stratified 1,000-genome subsample (10,000
synthetic MAGs, 9,887 returning hits) was searched against GTDB r226 — 732,475
accessions, a 73x larger reference universe. Exclusion here is ANI-only (>= 95),
since GTDB accessions carry no BV-BRC taxon id; 95 ANI is the species boundary, so
the criterion stands on its own.

| regime | compl | resolve rate | 95% CI | med best ANI | med AF_query | med AF_ref |
|---|---|---|---|---|---|---|
| include | 100% | 0.941 | 0.925-0.954 | 99.99 | 99.89 | 99.27 |
| include | 70%  | 0.936 | 0.919-0.950 | 99.99 | 99.88 | 69.30 |
| include | 30%  | 0.937 | 0.920-0.950 | 99.99 | 99.82 | 29.69 |
| exclude | 100% | 0.567 | 0.536-0.597 | 90.87 | 53.16 | 55.02 |
| exclude | 70%  | 0.563 | 0.532-0.593 | 90.93 | 54.12 | 39.38 |
| exclude | 50%  | 0.569 | 0.538-0.599 | 90.89 | 53.06 | 27.56 |
| exclude | 30%  | 0.564 | 0.533-0.594 | 90.64 | 52.41 | 16.49 |

Three results:

1. **A dense reference nearly doubles self-exclusion coverage, 0.31 -> 0.56.** The
   panel arm's low coverage was indeed reference sparsity, so its population-level
   numbers are lower bounds, as stated. Coverage is a property of the reference
   universe, not of the marker.
2. **Completeness-invariance replicates on a 73x larger database.** Resolve rate is
   0.563-0.569 across 100/90/70/50/30% completeness, every CI overlapping. The null
   axis is not an artifact of the sparse panel.
3. **The mechanism is now explicit in the numbers.** AF_ref ~= 0.99 x completeness
   (99.27, 88.87, 69.30, 49.54, 29.69) while AF_query is flat (99.89 -> 99.82) and
   median best ANI is flat (99.99 include, ~90.9 exclude). The selector gates on
   AF_query, so completeness is invisible to it by construction. 10% contamination
   moves AF_query proportionally (99.9 -> ~90.8) and nothing else.
   **AF_ref is a free, nearly unbiased completeness estimator that the pipeline
   currently computes and discards.**

Self-include resolve rate is 0.94 rather than 1.00 here because ~6% of panel
organisms have no >= 83-ANI match anywhere in GTDB, unlike the panel reference which
contains every source genome by construction.

Caveat on this arm: the 4 batches were finalized from a run whose driver was stopped
for memory safety, so 113 of 10,000 queries (1.13%) are absent. They are not a
truncation artifact — every absent query is a 0%-contamination one, spread evenly
across completeness tiers (21/24/25/22/21), i.e. ~2.2% of clean MAGs have no GTDB
hit at all while their contaminated siblings always hit the spliced donor. A
truncated file would have lost a contiguous tail of the sort order instead.

## Resource notes

A single skani search against the 57 GB GTDB sketch peaks at 33-38 GB RSS; eight
concurrent searches OOM-killed the box (confirmed in dmesg), which is why this arm
ran at 2-4 way concurrency and was scoped to 1,000 genomes. The panel arm (3.7 GB
sketch) ran 10 batches concurrently at 6 threads each, ~10x the throughput of one
61-thread search, because a single skani search saturates only ~5 cores regardless
of -t. This machine also carries four unrelated long-running jobs of the user's
(~60 GB, 13 days), which bounded the safe memory envelope throughout.

## Per-rank detail (panel arm, 0% contamination)

self-exclude, conditional / population:

| compl | answered | Phylum | Class | Order | Family | Genus | Species |
|---|---|---|---|---|---|---|---|
| 100% | 0.311 | 0.983/0.312 | 0.982/0.342 | 0.971/0.343 | 0.960/0.349 | 0.910/0.333 | 0.000/0.000 |
| 90%  | 0.312 | 0.982/0.312 | 0.980/0.342 | 0.970/0.343 | 0.957/0.349 | 0.910/0.333 | 0.000/0.000 |
| 70%  | 0.309 | 0.981/0.309 | 0.979/0.339 | 0.969/0.340 | 0.959/0.347 | 0.910/0.330 | 0.000/0.000 |
| 50%  | 0.310 | 0.980/0.310 | 0.979/0.339 | 0.968/0.341 | 0.959/0.348 | 0.911/0.331 | 0.000/0.000 |
| 30%  | 0.311 | 0.982/0.311 | 0.980/0.340 | 0.968/0.341 | 0.958/0.348 | 0.912/0.333 | 0.000/0.000 |

self-include is 1.000 at every rank and tier except Species (0.999).

**Species accuracy of 0.000 under exclusion is structural, not a failure.** Exclusion
removes every hit at >= 95 ANI, and 95 ANI *is* the species boundary, so no surviving
reference can be conspecific with the source. The amplicon arm has the same property
at tau = 0.987. Species should not be reported as a scored rank for either arm's
exclusion regime.

The conditional per-rank gradient is much shallower than the amplicon arm's: phylum
0.983 -> genus 0.910 is a 7-point decline over five ranks, where the amplicon arm
falls from family 0.79-0.86 to genus 0.60-0.73 (a 12-19 point family->genus gap
alone). When ANI finds a qualifying reference it is usually the right genus.

## Abstention composition (panel arm, self-exclude, 100% completeness)

n_ok 3,064 | n_abstained 6,773 | n_truth_empty 160 (excluded as undefined)

| cause | n | share of the 9,837 scoreable |
|---|---|---|
| `no_hit` — nothing survives exclusion at all | 3,436 | 35% |
| `low_AF_unresolved` — hits exist but none clears AF_query >= 15% | 3,379 | 34% |
| `<85` ANI — AF-qualified but below the 83 selector floor | 419 | 4% |
| resolved (85-95 ANI band, genus or family tier) | 3,064 | 31% |

This differs in kind from the amplicon arm, where 94.6% of self-exclusion
abstentions were records the filter emptied after truncation to the top 20 hits.
Nothing here is truncated: persistence depth was 100 from the outset and the median
MAG retains only 2 post-exclusion hits, so these records are genuinely empty or
genuinely low-AF. The MAG arm's self-exclusion numbers are therefore NOT bounded by
a retention-depth ceiling the way the amplicon arm's are.

## Leakage check (direct)

How often is the SELECTED anchor the source genome itself, by regime and tier
(0% contamination, n=9,997 each)?

| regime | 100% | 90% | 70% | 50% | 30% |
|---|---|---|---|---|---|
| self-include | 0.9991 | 0.9990 | 0.9987 | 0.9986 | 0.9988 |
| self-exclude | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

Two things follow. The exclusion filter is **airtight**: not one of 9,997 x 5
exclude-regime selections recovered the source genome, so no self-leakage inflates
the exclusion numbers. And self-inclusion is confirmed **near-tautological**: in
99.9% of cases the method retrieves the exact source genome, so gene capture of
1.000 is prediction-equals-truth by construction, not an accuracy measurement. Its
real content is that the retrieval stays at 99.9% even from a 30%-complete,
10%-contaminated draft.
