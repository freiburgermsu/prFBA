# Corrections after adversarial audit

Six auditors attacked RESULTS_v1.md and each objection was independently
adjudicated. Everything below was confirmed against the data before being acted on.
RESULTS_v1.md is superseded by RESULTS_v2.md.

## WITHDRAWN — every contamination claim (invalidating)

`take_to_fraction()` cannot skip the first fragment in shuffled order
(`if acc + len(f) > target and kept: continue`), so the minimum donor graft is one
whole fragment however small the intended dose. Verified over all 49,985
contaminated MAGs: contamination_actual min 0.047, median 0.0904 (on target),
**max 0.9676**; contam_bp max 1,592,099; 8.7% exceed 0.115, 5.3% exceed 0.15, and
**116 MAGs hold more donor than source sequence**.

Overshoot frequency rises monotonically with incompleteness — 1.7 / 2.3 / 3.4 / 5.7
/ **13.5%** at 100/90/70/50/30% — which is precisely the shape v1 read as a
"contamination x completeness interaction".

Decomposed, the effect vanishes on-dose: restricted to organisms at
contamination_actual <= 0.115 the deltas are **+0.0001 genus and +0.0029 recall**,
while the 2,070 over-dosed organisms give -0.139 genus and -0.099 recall. Dose-response
is flat to 0.15 and entirely above it.

The mechanism is also informative: for over-dosed MAGs the answered rate RISES
(0.242 -> 0.300) while conditional phylum accuracy COLLAPSES (0.968 -> 0.711) — the
donor graft aligns to the donor's own relatives, converting abstentions into
confidently wrong calls. That is a real finding about heavy contamination, but it is
not a measurement of a 10% dose.

**So: at a controlled 10% dose contamination has no measurable effect; the
synthesizer's dose control fails at low completeness and the v1 contamination
numbers measure that failure.** This does not touch the completeness-null headline
or the paired MAG-vs-amplicon table, whose rows are 0%-contamination cells.

## CORRECTED — conditional gene capture was understated

`score_genecap.py` classifies a selection whose references carry no PGFams as
`pred_empty` and keeps it out of the conditional pool; `score_mag.py` had no such
branch and logged 377 of them as `ok` with recall 0. All 43 anchor genomes involved
are present in the PGFam index with empty family sets (161 of 10,000 panel genomes
have none), so they are genuinely pred_empty, not missing. Applying the amplicon
arm's own rule to both sides:

| quantity | v1 | corrected |
|---|---|---|
| conditional gene-capture F1, 100% compl, exclude | 0.7966 | **0.8068** |
| paired MAG vs FullLength16S (100%) | 0.8088 / 0.7096 | **0.8136 / 0.7092** |
| paired MAG vs V4 (100%) | 0.8018 / 0.6414 | **0.8109 / 0.6413** |
| median paired delta vs FullLength16S | +0.063 | **+0.064** |
| median paired delta vs V4 | +0.128 | **+0.133** |

The correction favours the MAG arm; v1 was conservative.

## CORRECTED — numbers that did not reproduce

| claim in v1 | measured |
|---|---|
| "28.5% have no surviving hits" (was one batch) | **34.4%** at 100% completeness, rising to **36.1%** at 30% |
| "self-hit ANI min 98.07" (was one batch) | min **90.85** over all 99,970 queries; median 100.00 |
| "every self-hit is >=95 ANI, 0 of 5,000 survive" | **2 of 99,970** self-hits fall below 95 and would survive exclusion |
| "AF_query stays ~95%, AF_ref ~= 0.95 x completeness" (pilot carry-overs) | AF_query ~99.9, AF_ref ~= **0.99** x completeness |
| "gene capture cond/pop" column label | those columns are **recall**, not F1 |
| "PGFam truth byte-identical to the amplicon arm" | true for 9,745 of 9,746; genome 33953.3 has truth_size 0 |

The no_hit fraction is NOT flat in completeness (34.4 -> 36.1%). The resolve rate is
flat not because the failure modes compensate but because the organisms moving from
low_AF to no_hit were unresolved under either label.

## NOTED — scope limits v1 did not state

* **254 panel organisms are absent from every taxonomic denominator.** The MAG arm
  reuses `benchmark_truth.json`, which exists only for the 9,746 genomes that
  yielded an amplicon. Taxonomic denominators are therefore 9,746-based, not 9,997.
* **The paired comparison covers 7,295 of 9,997 organisms** (FullLength16S has no
  cell for 2,702 panel genomes), and the conditional comparison 1,510.
* **`fragment()` discards every deposited contig shorter than 2,500 bp** and
  computes `total` after the loss, so `completeness_actual` reads 1.0000 for a MAG
  that is already missing its short contigs.
* **The retention rule is a size-biased greedy fill, not a random drop** (the same
  `and kept: continue` branch), so contiguity is mildly confounded with completeness.
* **238 of 640 stored p-values underflowed to literal 0.0.** They are now reported
  as "< 1e-300"; `median_delta` now averages the two middle values for even n.
* **95 ANI and 0.987 16S identity remove very different amounts.** The MAG rule
  removes a median of 1 hit (67.4% remove exactly the self genome); the amplicon
  rule removes a median of 4. Matching on stringency instead (exclude_0.995) shrinks
  the genus margin by ~37%, though the MAG arm still wins at both. This is the
  reference-density asymmetry already conceded, surfacing in the exclusion rule.

## Both synthesizer bugs are now FIXED in synth_mags.py

`take_to_fraction()` gained an `allow_empty` mode used for the contamination graft:
the first shuffled fragment is now skippable and an overshooting fragment is
truncated rather than taken whole. Verified: realized contamination is now exactly
0.0909 at 100 / 70 / 30 / 10 / 3% completeness, where it previously reached 0.97.
`synth_one()` now records `completeness_actual` against the DEPOSITED genome length
rather than the post-`MIN_FRAG` fragment pool (and keeps the old ratio as
`completeness_vs_pool`). Measured exhaustively over all 9,997 genomes, the retained
pool is the full source for the median genome but falls below 95% of it for **517
genomes (5.2%)**, below 75% for **114 (1.1%)**, below 50% for **37 (0.4%)** and
below 25% for 9; the worst retains **6.2%**. So 517 of the MAGs labelled "100%
complete" hold under 95% of their source, and the loss scales with the source
assembly's fragmentation, which is the confound the fragmentation step was meant to
remove. Two things limit the damage: the loss is a constant per-genome offset
applied identically to every tier, so it cannot create or destroy a completeness
trend; and its direction is conservative -- the 517 lossy genomes answer at 0.193 vs
0.313 and reach genus 0.812 vs 0.913, so the defect drags the MAG arm down rather
than flattering it.

`run_full.sh` / `run_gtdb.sh` therefore reproduce a correctly-dosed grid. The
results in RESULTS_v2.md come from the pre-fix run and are 0%-contamination only,
which the contamination bug does not touch.
