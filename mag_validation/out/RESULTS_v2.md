# MAG-matching validation on the 10,000-genome panel — results

Run 2026-10-08 on 9,997 BV-BRC genomes of the 10,000-genome panel used for the 16S
amplicon validation. 99,970 synthetic MAGs, 199,940 scored cells, 1,428,212 skani
hit rows, plus a 10,000-MAG dense-reference arm against GTDB r226.

Supersedes RESULTS_v1.md. Every number here survived an adversarial audit; what the
audit changed or withdrew is in CORRECTIONS.md and summarised at the end.

**All contamination tiers are excluded from every result below.** The synthesizer's
dose control fails at low completeness (see CORRECTIONS.md), so the contaminated
cells measure that failure rather than a 10% dose. Everything reported here is the
0%-contamination grid.

## Headline

**Completeness is a null axis from 100% down to 30%.** Nothing the pipeline decides
moves with it — not which reference is selected, not whether one is selected at all,
not how much metabolism is recovered.

| self-exclude, 0% contamination | 100% | 90% | 70% | 50% | 30% |
|---|---|---|---|---|---|
| resolve rate | 0.307 | 0.308 | 0.306 | 0.307 | 0.307 |
| conditional gene-capture F1 | 0.8068 | 0.8070 | 0.8071 | 0.8067 | 0.8052 |
| conditional genus accuracy | 0.910 | 0.910 | 0.910 | 0.911 | 0.912 |
| conditional family accuracy | 0.960 | 0.957 | 0.959 | 0.959 | 0.958 |

Paired within organism against the amplicon arm's own cells on the same genomes:

| regime | MAG tier | vs | n | MAG | 16S | median delta | MAG wins | 16S wins |
|---|---|---|---|---|---|---|---|---|
| exclude | 100% | FullLength16S | 1510 | 0.8136 | 0.7092 | **+0.064** | 1212 | 294 |
| exclude | 30%  | FullLength16S | 1510 | 0.8107 | 0.7092 | **+0.063** | 1206 | 301 |
| exclude | 100% | V4            | 1573 | 0.8109 | 0.6413 | **+0.133** | 1409 | 163 |
| exclude | 30%  | V4            | 1558 | 0.8094 | 0.6388 | **+0.137** | 1396 | 161 |

Wilcoxon signed-rank p < 1e-140 throughout (many underflow double precision).

**The question the proposal posed — "the completeness at which a MAG stops beating a
full-length 16S amplicon" — has the answer "it does not, down to 30% completeness."**
A 30%-complete MAG wins by the same margin as a complete genome.

## Mechanism

* **100% of 99,970 synthetic MAGs recover their own source genome** among the hits
  (median self-hit ANI 100.00, min 90.85; only 2 fall below the 95 exclusion cut).
  The amplicon arm's retention-depth ceiling does not exist here — a whole genome
  has no equal-score tie clouds, and persistence depth was 100 from the outset.
* **AF_query is flat (99.89 -> 99.82) while AF_ref ~= 0.99 x completeness**
  (99.27 / 88.87 / 69.30 / 49.54 / 29.69). **AF_ref is a free, nearly unbiased
  completeness estimator the pipeline already computes and discards.**
* The gate is blind but the hit screen is not. Completeness is invisible to the AF
  gate, which reads AF_query; it is NOT invisible to skani's screening step, where
  the post-exclusion candidate pool shrinks from a mean of 11.60 to 10.80 hits per
  query between 100% and 30% completeness (median 2 at every tier, so the loss is in
  the tail). The net effect on the resolve rate is nil, but "invisible by
  construction" would overstate it: incompleteness costs marginal candidates, it just
  does not cost the ones that were going to be selected.
* Exclusion is airtight: **0 of 9,997 x 5 exclude-regime selections** recovered the
  source genome. Self-inclusion is correspondingly near-tautological — 99.9% of
  include-regime anchors ARE the source genome, so its 1.000 gene capture is
  prediction-equals-truth by construction. Its real content is that retrieval stays
  at 99.9% from a 30%-complete draft.

## Abstention composition (self-exclude)

| cause | 100% compl | 30% compl |
|---|---|---|
| `no_hit` — nothing survives exclusion | 34.4% | 36.1% |
| `low_AF_unresolved` — hits exist, none clears AF_query >= 15% | 33.8% | 32.0% |
| `<85` ANI, below the 83 selector floor | ~4% | ~4% |
| resolved | 30.7% | 30.7% |

no_hit rises 34.4 -> 36.1% while low_AF falls 33.8 -> 32.0%; the resolved fraction
is unchanged because the organisms that lose their marginal-AF relatives as
completeness drops were not being resolved anyway. Nothing here is
a truncation artifact, unlike the amplicon arm where 94.6% of self-exclusion
abstentions were records emptied after truncation to the top 20 hits.

## Per-rank (self-exclude, conditional / population)

| compl | Phylum | Class | Order | Family | Genus | Species |
|---|---|---|---|---|---|---|
| 100% | 0.983/0.312 | 0.982/0.342 | 0.971/0.343 | 0.960/0.349 | 0.910/0.333 | 0.000/0.000 |
| 30%  | 0.982/0.311 | 0.980/0.340 | 0.968/0.341 | 0.958/0.348 | 0.912/0.333 | 0.000/0.000 |

**Species 0.000 is structural, not a failure**: exclusion removes every hit at >= 95
ANI and 95 ANI *is* the species boundary, so no surviving reference can be
conspecific. Species should not be reported as a scored rank for either arm's
exclusion regime.

## Dense-reference arm (GTDB r226, 1,000-genome subsample)

| regime | compl | resolve rate | 95% CI | med ANI | med AF_query | med AF_ref |
|---|---|---|---|---|---|---|
| include | 100% | 0.941 | 0.925-0.954 | 99.99 | 99.89 | 99.27 |
| include | 30%  | 0.937 | 0.920-0.950 | 99.99 | 99.82 | 29.69 |
| exclude | 100% | 0.567 | 0.536-0.597 | 90.87 | 53.16 | 55.02 |
| exclude | 30%  | 0.564 | 0.533-0.594 | 90.64 | 52.41 | 16.49 |

**A dense reference nearly doubles self-exclusion coverage, 0.31 -> 0.56**, and
completeness-invariance replicates on a 73x larger database (every CI overlapping).
Coverage is a property of the reference universe, not of the marker.

## What this does and does not support

**Supported.** On organisms where whole-genome ANI finds a qualifying reference, a
MAG beats a 16S amplicon on gene capture by +0.064 (full-length) to +0.133 (V4) F1,
and does so just as well at 30% completeness as at 100%. Reference selection by ANI
is insensitive to MAG completeness over the whole MIMAG medium-quality-and-above
range.

**Not supported.**
* *Not* a claim that a MAG beats an amplicon in general. The conditional comparison
  covers the 31% of organisms the MAG arm answers; on the other 69% it abstains
  entirely while the amplicon arm still answers at F1 0.62 (full-length) / 0.57 (V4).
  Scored over the whole population the amplicon arm wins (0.535 vs 0.248 F1).
* The MAG-answered subset is **+0.065 easier for the amplicon arm too** (16S
  conditional F1 0.7096 there vs 0.6450 over everything it answers). The pairing
  means this does not bias the difference, but it bounds the population covered.
* **Reference-universe asymmetry.** The MAG arm searched 9,997 genomes; the amplicon
  arm 459,301 16S references from ~112,000 genomes. The panel is species-dereplicated
  (0 of 9,997 share a species; 61.7% share a genus). Population-level MAG numbers are
  lower bounds.
* **Exclusion stringency is not matched.** The 95-ANI rule removes a median of 1 hit;
  the amplicon arm's 0.987 removes a median of 4. Matching stringency instead
  (exclude_0.995) shrinks the genus margin ~37%; the MAG arm still wins at both.
* **In-silico degradation is optimistic**, exactly as in-silico PCR is for the
  amplicon arm: dropping fragments cannot produce chimeric contigs, strain-mixture
  inflation or coverage-biased mis-binning. Both arms are upper bounds, which is what
  keeps them comparable.
* **Contamination is unmeasured.** The grid exists but its dose control is broken.

## Comparability

Truth, metric definitions and name resolution are the amplicon arm's, not
reimplementations: `truth_lineage` from `benchmark_truth.json` (NCBI via taxopy, null
at unnamed ranks -> unscoreable, never wrong); predicted lineage re-resolved through
the same `_common.lineage_for()`; names compared after `_common._norm()`; three-state
cells; conditional vs population (anti-inflation) denominators; `pred_empty` kept out
of the conditional pool and folded into the population pool per `score_genecap.py`;
Wilson 95% CI for proportions, normal-approx mean CI for means of proportions. PGFam
truth sets match the amplicon arm's for 9,745 of 9,746 genomes. Taxonomic
denominators are 9,746-based: 254 panel organisms have no entry in
`benchmark_truth.json`.

## Audit outcome

One invalidating finding (all contamination claims withdrawn — synthesizer dose
overshoot), one material correction that **raised** the MAG numbers (`pred_empty`
handling, conditional F1 0.7966 -> 0.8068), four material corrections to
non-reproducing secondary numbers, and several cosmetic fixes (median of even-length
samples, p-value underflow, a recall/F1 column label). The completeness-null headline
and the paired MAG-vs-amplicon result survived unchanged. Details in CORRECTIONS.md.

## Low-completeness sweep: where does the MAG arm fall to 16S?

Nine further tiers (25/20/15/10/7/5/3/2/1% completeness, 0% contamination,
9,997 genomes each = 89,973 additional synthetic MAGs) were run to find the
crossover. Realized completeness tracks target to within 0.1% at the median of every
tier; at 1% the median MAG holds 37 kb of source sequence. Only 26 of 89,973 queries
found no reference at all, and they concentrate exactly where expected (0 at
25-10%, 1 at 7%, rising to 17 at 1%).

| completeness | resolve rate | conditional F1 | MAG / 16S full-length | delta | MAG / V4 | delta |
|---|---|---|---|---|---|---|
| 100% | 0.3075 | 0.8068 | 0.8136 / 0.7092 | **+0.104** | 0.8109 / 0.6413 | **+0.170** |
| 50%  | 0.3074 | 0.8067 | 0.8121 / 0.7087 | +0.103 | 0.8100 / 0.6405 | +0.170 |
| 25%  | 0.3092 | 0.8062 | 0.8106 / 0.7069 | +0.104 | 0.8099 / 0.6391 | +0.171 |
| 15%  | 0.3009 | 0.8049 | 0.8097 / 0.7082 | +0.102 | 0.8086 / 0.6405 | +0.168 |
| 10%  | 0.2953 | 0.8027 | 0.8056 / 0.7042 | +0.101 | 0.8045 / 0.6359 | +0.169 |
| 5%   | 0.2917 | 0.7965 | 0.7962 / 0.7061 | +0.090 | 0.7986 / 0.6383 | +0.160 |
| 3%   | 0.2879 | 0.7926 | 0.7932 / 0.7032 | +0.090 | 0.7957 / 0.6398 | +0.156 |
| 2%   | 0.2826 | 0.7913 | 0.7922 / 0.7060 | +0.086 | 0.7903 / 0.6364 | +0.154 |
| 1%   | 0.2533 | 0.7830 | 0.7806 / 0.7009 | **+0.080** | 0.7813 / 0.6399 | **+0.141** |

n = 1,260-1,580 paired organisms per row. These are differences of paired means; the
median per-organism paired difference is smaller (+0.064 at 100%, +0.050 at 1%)
because the difference distribution is right-skewed -- many near-ties plus a tail of
large MAG wins. At 1% the MAG arm still wins 928 organisms to 331, p = 9e-81. Both
statistics are in `out/crossover.json`.

**There is no crossover.** At 1% completeness -- a 37 kb fragment, 50x below MIMAG's
medium-quality floor and not a MAG any study would retain -- whole-genome ANI still
beats a full-length 16S amplicon by +0.080 F1 and V4 by +0.141.

The reason is structural: the MAG arm's prediction is the gene content of the
*selected reference genome*, not of the MAG itself. MAG completeness therefore
matters only insofar as it changes which reference gets selected, and ANI is robust
enough that it rarely does. Degrading the query by two orders of magnitude costs
0.024 of conditional F1 (0.8068 -> 0.7830).

**What does degrade is coverage**: the resolve rate falls 0.308 -> 0.253 (-18%
relative), flat to ~15% completeness and then declining, with the sharpest drop in
the last step (0.283 at 2% -> 0.253 at 1%) as the remaining sequence stops clearing
skani's screen. The practical limit on a low-completeness MAG is whether a reference
is found at all, not how good the answer is when one is.

Figure: `out/F_mag_crossover.pdf` / `.png`.

## Why the MAG arm abstains on 69% of organisms

Four stacked causes, measured at 100% completeness, 0% contamination, self-exclusion.

**1. Marker reach is the dominant cause.** ANI is a whole-genome measure, and at genus
distance two genomes retain high local identity over only a small fraction of their
length. Among the 6,011 organisms that HAVE a same-genus congener in the panel, the
best surviving hit has a median AF_query of just **20.3%** (p25 = 5.1%), so **45.1%
fall below the 15% AF gate** even though a genuine genus-mate is present. A 16S
amplicon has no such problem: the whole amplicon aligns by construction, so its
aligned fraction is ~100% regardless of distance, and its identity stays 0.90-0.99
across genus and family. Binned by the amplicon arm's own `best_nonself_identity`:

| 16S best non-self identity | MAG resolves | 16S answers | median post-exclusion ANI |
|---|---|---|---|
| >= 0.995 | 0.484 | 0.642 | 89.8 |
| 0.987-0.995 | 0.430 | 0.879 | 87.4 |
| 0.97-0.987 | 0.202 | **0.983** | 85.1 |
| 0.945-0.97 | 0.082 | **0.988** | 85.7 |
| 0.90-0.945 | 0.049 | **0.986** | 89.7 |
| < 0.90 | 0.045 | 0.653 | 86.7 |

At 16S identity 0.945-0.987 -- the genus-to-family range -- the amplicon answers
98-99% of the time and whole-genome ANI answers 8-20%. (The amplicon's own dip to
0.642 in the closest bin is the retention-depth artifact that bounds its
self-exclusion figures, not a reach failure.)

**2. The detection floor.** skani stops reporting below ~80 ANI. 15.9% of organisms
with a same-genus congener present have no surviving hit at all; that rises to 51%
when only a family relative exists and 72% when nothing at order level or below does.

**3. Reference sparsity amplifies both.** The panel is species-dereplicated, so
exclusion leaves genus-level relatives at best. A 73x denser reference (GTDB r226)
nearly doubles the resolve rate, 0.31 -> 0.56.

**4. The ANI floor is almost irrelevant; the AF gate is almost everything.**
Counterfactually, holding one threshold and sweeping the other:

| change | resolve rate |
|---|---|
| as shipped (ANI >= 83, AF >= 15%) | 0.319 |
| ANI floor 83 -> 0, AF gate unchanged | 0.323 (**+0.004**) |
| AF gate 15% -> 5%, floor unchanged | 0.388 |
| AF gate 15% -> 2%, floor unchanged | 0.462 (**+0.143**) |
| both relaxed (ANI >= 80, AF >= 2%) | 0.579 |
| any surviving hit at all (ceiling) | 0.667 |

**The gate is not a defect -- it is the method working.** The references a looser gate
would admit are genuinely worse: the 1,367 organisms that newly resolve at AF >= 2%
have mean gene-capture F1 **0.599**, against 0.807 for those already resolving, and
only 26.9% of them beat the full-length 16S constant of 0.709. Relaxing the gate
would raise coverage by 45% relative while pooling mean F1 down from 0.807 to 0.743,
and would do so by adding references that are, on average, worse than simply using
the cheaper marker. prFBA abstains precisely where answering would be worse than a
16S amplicon -- which is the behaviour the method claims and, here, the behaviour it
delivers.
