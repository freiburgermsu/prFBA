# Run 5, second half (R5.3 decks, R5.7 manuscript), 2026-11-16

Executor: Run 5B agent (the decks were done by a subagent it briefed; its report is folded in below). Inputs:
`RUN5A_REPORT.md`, `prFBA/manuscript_numbers/number_diff.{json,md}` (178 claims), `plan/PLAN.md` §1 and §5, and the
RUN3 / RUN4_LOCAL / RUN4_TF reports.

**Outcome.**
- **Paper:** commits **`75d6f03`** and **`2387ac2`** on `main` of `~/Documents/prFBA_paper` (not pushed; origin/main was
  re-fetched just before the first commit and had no new commits). `2387ac2` narrows the provenance sentence to
  "every result that reads the gene-family cache", adds the legacy as-shipped figure to the §select list of
  non-default-rule results (and a `\label{sec:discussion}` for `\nameref`), and changes "below" to "in this section". One file edited, `prFBA.tex` (1,241 → 1,251 lines), plus the four
  re-rendered figures. Messages: `logs/commit_msgs/r5b_paper.txt`, `r5b_paper_2.txt`.
- **Decks:** commit **`a6812be`** on `hardware-aware-alignment-dispatch` of `~/Documents/prFBA` (path-limited, not
  pushed). The four `*_2026-11.pptx` decks were rebuilt; no `*_2026-10*` deck and not the untracked user-edits file
  were touched.
- **Build:** `pdflatex prFBA && bibtex prFBA && pdflatex prFBA && pdflatex prFBA` → 98 pages, **0 undefined references
  or citations** (bibtex: 0 warnings), and **only the two pre-existing overfull boxes** (measured on a baseline build
  of the unmodified tex, `scratchpad/r5b/build0`): `\hbox 23.69pt` at the reference-database paragraph (L133–138) and
  `\hbox 16.96pt` in `tab:picrust2` (now L554–571). The edits briefly added a 0.34 pt overfull `\vbox` on the float
  page holding `tab:picrust2` + F13. Including F13 at `0.98\linewidth` removes it; no text was changed for layout.
- **Numbers:** about **480 printed numbers change in place** (sequence alignment of all numeric tokens, old vs new tex:
  477 replaced, plus 85 removed and 156 added inside rewritten sentences). Of these, 234 are table cells in 58 rows,
  regenerated from the source JSON or cells; the rest are prose values.
- **Meaning:** every one of the 24 meaning-change items was rewritten. Old and new text is given below. Three further
  claims changed meaning and are rewritten too (the "three fifths" bound, the 0.995 denominators, the legacy
  as-shipped lead).

---------------------------------------------------------------------------------------------------------------

## Where each step ran

Everything ran on this machine (FM). LaTeX, git and text editing stay here by the user's rule. The small checks were:
- the unrounded MAG paired means (`compare_arms.py`/`crossover.py` logic, imported read-only, over a 22 MB and a 20 MB
  CSV plus the low-completeness batches);
- three cells-table reads (as-shipped 10 MB gz, shared-reference family 13 MB gz);
- JSON reads.

Each ran in under a minute and under 1 GB, so none was worth a TF push. The deck subagent built the decks here in
`systemd-run --user` units with a memory cap; none ran longer than about a minute. TF does not hold the git-ignored
depth-100 inputs the generators read. Nothing ran on TF and nothing wrote the cache.

---------------------------------------------------------------------------------------------------------------

## Manuscript: how the numbers were applied

Line numbers in this section are those of the old file (`0214a14`) unless marked "new".

1. **Tables (`scratchpad/r5b/tables.py`).** Every numeric cell of the nine re-scored tables was regenerated:
   - `tab:main` and `tab:tau` gene-capture columns come from `gene_capture_by_region.json` and
     `gene_capture_tau_compare.json`; their taxonomy columns are unchanged by construction.
   - `tab:picrust2`, `tab:sharedref`, `tab:twopanel` and `tab:emptied` come from number_diff labels.
   - `tab:magvs16s` and `tab:magval` come from an unrounded recomputation (see "Rounding").
   - `tab:threepanel` comes from `three_panel_comparison_d100.json` and `three_panel_mag_panelarm.json`; the GTDB
     block's source is byte-unchanged.

   For every regenerated cell the old value was checked against the print. The only mismatches are the known
   old-side rounding boundaries (`69.05→69.1`, AS-20 0.8405 / 0.6415 / 0.7855, magvs16s 0.8155 / 0.5035 / 0.1575),
   which were not "repaired".
2. **Prose (`locate.py`, `apply.py`).** Each changed value was located at its claim's tex lines, with one replacement
   per position and right-to-left within a line. Multi-hit positions were checked by hand. Excluded positions were
   either historical or handled by hand:
   - L646, the pre-registered essentiality, stays 0.068;
   - the McNemar counts and q values;
   - the A-C13 row counts.

   Two locator conflicts were fixed by hand: L823 had to read the 100–30% range high values 0.314 / 0.807, and L1095
   the V4 CI 0.239. One locator bug turned "PICRUSt2" into "PICRUSt4" at L79; it was found by a letter-adjacency check
   and fixed.
3. **Sweeps.** I then grepped the whole file for every changed value's old string and classified each survivor. This
   found repeats that no claim's line range covered; all are now fixed:
   - L402 gene capture "58--64%" → 59–65%;
   - L1162 emptied rates 4.7 / 33.4% → 4.6 / 33.3%;
   - L521 "67,038" → 67,713;
   - L891 τ = 99 population 0.248 / 0.301 / 0.605 → 0.250 / 0.304 / 0.612;
   - L899 V4 n 2,299 → 2,441;
   - L1062 arithmetic inputs 0.334 / 0.047 → 0.333 / 0.046.

   A second check covered values whose `changed_printed` flag is false although the print differs from the new value,
   because the old value sat on a rounding boundary:
   - family exact-tie 0.623 → 0.624 and forced 0.709 → 0.710 (L624);
   - PICRUSt2 rxn conditional F1 0.721 → 0.722 (L527);
   - PICRUSt2 recall by NSTI bin 0.870 / 0.769 → 0.871 / 0.768 (L530, L1179).

   Reverse check: every changed value's new print is present near its claim (the 5 exceptions are values written as
   words, e.g. "four media", or the deliberately historical 43,140).
4. **Hand-carried values (number_diff `manual_values`).**
   - Stage O shuffled vs released: +0.0110 [0.0031, 0.0222] → 0.0111 [0.0031, 0.0223] (L319).
   - Stage E isotonic recalibration: 0.0411 → 0.0412 (L358).
   - Stage E within-union ΔF1 CI: +0.0096 → +0.0097 (L360).
   - p-values come from the unrounded `new` field: Spearman 1.47e-7 → 1.5×10⁻⁷ (L634, L1190); Welch FL 1.47e-4 →
     1.5×10⁻⁴ (L1020); MAG 1% Wilcoxon 3.4e-83 → 3×10⁻⁸³ (L826). The "p < 10⁻¹⁶¹" and "p < 10⁻²⁹" bounds still hold
     (9.4e-167, 1.1e-30). All bound-type claims ("> 99", "≥ 0", "< 0.01", "within 0.010", "original ahead at every
     completeness") still hold except "within 0.016", now 0.017 (rewritten).

### Rounding (policy, and every boundary case)

collect.py's `new_fmt` formats the stored value with Python's rounding. For sources stored at 4 dp this is unreliable
when the 4th decimal is 5. The policy:
- **MAG paired / crossover cells** (`tab:magvs16s`, `tab:magval`): recomputed unrounded from the cells with the
  producers' own logic (`scratchpad/r5b/mag_unrounded.json`). This reproduces every stored 4-dp value. Δ is the
  difference of the unrounded means, as the `tab:magvs16s` caption says.
  - Three cells differ from `new_fmt`: V1–V3 16S 0.704464 → **0.704** (`new_fmt` 0.705), V4 MAG 0.808503 → **0.809**
    (0.808), V7–V9 median Δ 0.125462 → **+0.125** (+0.126).
  - In `tab:magval` the same applies at 1%: FL MAG 0.784522 → **0.785** (`new_fmt` 0.784).
- **Three-panel bootstrap CIs** (4 dp stored, no finer source): rounded half up, the paper's convention for this
  table. This gives genus conditional-accuracy V3–V4 CI low −0.0385 → **−0.039** (`new_fmt` −0.038).
- **Point estimates recomputed from the unrounded phystd fields:** panel-as-reference Δ at 100% −0.026520 → **−0.027**
  (`new_fmt` −0.026).
- **`tab:emptied` shares come from the counts:** random 50–99 share 1,269/8,780 = 0.14453 → **14.5%** (`new_fmt` 14.4).
- **As-shipped, from `comparison_cells.csv.gz`:**
  - prFBA FL population EC F1 0.776505 → **0.777** (`new_fmt` 0.776; L82, L526);
  - the source-absent PICRUSt2 coverage stays **0.999** (0.999478; `new_fmt` 1.000).

  The family-rank PICRUSt2 population F1 (0.749517 → 0.750) and the exact-tie value (0.623665 → 0.624) were confirmed
  from the family cells.
- **Raw value unchanged:** where `new_fmt` differs from the print only because the 4-dp value is unchanged
  (SR-05 0.788, SR-07 −0.092, SR-08 0.030, FR-03 CI −0.053, the depth probe's 80.2%, the deployment's 2.289), the print
  was kept.
- **Float noise:** FR-03 "exact-tie vs nn1" (−0.0835 as a difference of two 4-dp values) keeps −0.083.
- **Everything else:** `new_fmt`.

### Figures (step 4)

Copied from `prFBA/manuscript_figures/out/`; the md5s match the Run 5A report:

| tex include | source | md5 |
|---|---|---|
| `F4b_percent_genecapture_by_region.pdf` | same name | eee40cef1c3ea6bb4fd6afee16e2302c |
| `F13_prfba_vs_picrust2.pdf` | same name | f7dfa4a58721f09612fb099c481dc2c6 |
| `F14_mag_completeness.pdf` | `F_mag_completeness.pdf` | 2f9b18cccfc2596ac1b69a7931f75ef3 |
| `F15_mag_crossover.pdf` | `F_mag_crossover.pdf` | 03ca35710ed5e1ea807385b42d9bf50b |

- F1, F2, F11 and F16 are pixel-identical (Run 5A) and were not copied.
- Captions agree with the figures: the F14 note reads "+0.102 / +0.209", the F15 note "0.313 → 0.258", and
  F13 n = 67,713. Each matches the text.
- The fig:genecap and fig:picrust2 captions were re-numbered (59–65%; 8%, 548 of 6,660, 147, 695, 6,423, 0.234).

---------------------------------------------------------------------------------------------------------------

## Meaning changes: old and new text

The line numbers below are from the old file (`0214a14`) and the new file (`75d6f03`; `2387ac2` adds no lines). The provenance
sentence of item 9 was then narrowed in `2387ac2` to: "Every result here that reads the gene-family cache was recomputed on the
repaired cache, … the exceptions, … are the backend concordance check, the shared-reference panel's eligibility, and the
phenotype passes and side studies run on the earlier phenotype sample."

#### 1. Phenotype vs nn1 (sec:sharedref)

Old:

> L635: Against the nn1 baseline the phenotype endpoint favors prFBA: its models mismatch the truth profile on 0.146 fewer media than nn1's ($n = 755$, $p = 0.017$, $q = 0.028$).
>
> L636: That does not reverse the EC result against nn1.
>

New:

> L642: Against the nn1 baseline the phenotype endpoint detects no difference either: prFBA's models mismatch the truth profile on 0.030 fewer media than nn1's ($n = 759$, Wilcoxon $p = 0.41$, $q = 0.51$), and over all 1{,}000 paired genomes the point estimate favors nn1 instead (by 0.06 media, $p = 0.46$).
>
> L643: That matches the EC result against nn1.
>

#### 1b. Phenotype vs nn1 (discussion)

Old:

> L1186: Against the nn1 baseline, prFBA's profiles sit slightly closer to truth over the 755 ATP-informative genomes it answers ($-0.146$ media, $p = 0.017$, $q = 0.028$; $p = 0.28$ over all paired genomes), but this is not a reversal of the EC result: on those same genomes the two tie on ECs (EC F1 0.782 against 0.781, $p = 0.64$), and prFBA's EC shortfall against nn1 comes entirely from abstentions scored as zero.
>

New:

> L1196: Against the nn1 baseline the profiles show no difference either ($-0.030$ media for prFBA over the 759 ATP-informative genomes it answers, $p = 0.41$, $q = 0.51$; $+0.06$ over all 1{,}000 paired genomes, $p = 0.46$), matching the EC result: on those same genomes the two tie on ECs (EC F1 0.782 against 0.780, Wilcoxon signed-rank $p = 0.84$), and prFBA's EC shortfall against nn1 comes entirely from abstentions scored as zero.
>

#### 2. McNemar media (abstract)

Old:

> L79: The phenotype readout detects no mean difference (per medium, PICRUSt2 is ahead on two of 22 media), but it cannot test most of the gap: on $n=755$ ATP-informative genomes (747 of them V4) where prFBA's exact-tie arm answered and the two methods' EC F1 differ by only $-0.028$, ATP-production profiles of draft models built from each method's predicted ECs mismatch truth on 0.13 more media for prFBA (Wilcoxon $p=0.49$; only 22 of 54 media carry ATP signal for these arms) --- a null, not a demonstrated equivalence.
>

New:

> L79: The phenotype readout detects no mean difference (per medium, PICRUSt2 is ahead on four of 22 media), but it cannot test most of the gap: on $n=759$ ATP-informative genomes (750 of them V4) where prFBA's exact-tie arm answered and the two methods' EC F1 differ by only $-0.031$, ATP-production profiles of draft models built from each method's predicted ECs mismatch truth on 0.07 more media for prFBA (Wilcoxon $p=0.79$; only 22 of 54 media carry ATP signal for these arms) --- a null, not a demonstrated equivalence.
>

#### 2b. McNemar media (sec:sharedref)

Old:

> L631: On this endpoint's mean the EC-space deficit is not detected, although the assay's per-medium McNemar tests favor PICRUSt2 on two of the 22 testable media (Ac.NO and Etho.O2, $q = 0.002$).
>
> L639: Taken with the EC result, the reading is that prFBA's predictions are measurably worse as annotation sets while showing no detectable mean difference as inputs to a metabolic model at this assay's resolution (per medium, PICRUSt2 is ahead on two media).
>

New:

> L636: On this endpoint's mean the EC-space deficit is not detected, although the assay's per-medium McNemar tests favor PICRUSt2 on four of the 22 testable media (Ac.NO, Ac.O2, Etho.O2 and For.O2, $q = 0.029$--$0.046$).
>
> L646: Taken with the EC result, the reading is that prFBA's predictions are measurably worse as annotation sets while showing no detectable mean difference as inputs to a metabolic model at this assay's resolution (per medium, PICRUSt2 is ahead on four media).
>

#### 2c. McNemar media (discussion)

Old:

> L1185: That qualifier earns its place, because the one endpoint built on top of those predictions shows no significant overall gap: ATP-production profiles from draft models are not detectably further from truth (prFBA mismatches truth on 0.13 more media than PICRUSt2, Wilcoxon $p = 0.49$, $q = 0.49$, $n = 755$), although the assay's per-medium McNemar tests detect PICRUSt2 advantages on two media (Ac.NO and Etho.O2, $q = 0.002$).
>

New:

> L1195: That qualifier earns its place, because the one endpoint built on top of those predictions shows no significant overall gap: ATP-production profiles from draft models are not detectably further from truth (prFBA mismatches truth on 0.07 more media than PICRUSt2, Wilcoxon $p = 0.79$, $q = 0.79$, $n = 759$), although the assay's per-medium McNemar tests detect PICRUSt2 advantages on four media (Ac.NO, Ac.O2, Etho.O2 and For.O2, $q = 0.029$--$0.046$).
>

#### 3. Phenotype sample, readout and Spearman

Old:

> L634: Over 755 paired ATP-informative genomes this arm mismatches the truth profile on 0.13 more media than PICRUSt2 does, which is not distinguishable from zero (Wilcoxon $p = 0.49$, $q = 0.49$), and the per-genome correlation between the two methods' EC-F1 difference and their phenotype difference is only $-0.14$.
>
> L1193: The phenotype result is a null at 755 genomes, not an equivalence; on those genomes prFBA's mean EC-F1 gap to PICRUSt2 is only $-0.028$, so the assay tests the propagation of a small deficit, and the ATP pass, re-run on drafts built with the template-direction correction (amendment A-C12) and again on the later ModelSEEDpy (A-C13), reproduces every ATP vector and summary statistic of the pre-registered pass (the \nameref{sec:sharedref} section).
>
> L1190: What survives the controls is narrower but still real: the choice of which amplicons to decline carries information, being worth $+0.021$ population F1 over declining as many at random ($+0.019$ for the exact-tie rule against its own masked run); a prFBA prediction names the genomes it came from, which an ancestral-state reconstruction cannot; and the EC deficit leaves no detectable mean difference in the ATP phenotype, although per genome it carries through weakly (Spearman $\rho = -0.14$ between the EC-F1 and ATP-distance differences, $p = 2\times10^{-4}$).
>

New:

> L641: Over 759 paired ATP-informative genomes this arm mismatches the truth profile on 0.07 more media than PICRUSt2 does, which is not distinguishable from zero (Wilcoxon $p = 0.79$, $q = 0.79$), and the per-genome correlation between the two methods' EC-F1 difference and their phenotype difference is weak, though distinguishable from zero (Spearman $\rho = -0.19$, $p = 1.5\times10^{-7}$).
>
> L1203: The phenotype result is a null at 759 genomes, not an equivalence; on those genomes prFBA's mean EC-F1 gap to PICRUSt2 is only $-0.031$, so the assay tests the propagation of a small deficit, and on the sample of the time the ATP pass, re-run on drafts built with the template-direction correction (amendment A-C12) and again on the later ModelSEEDpy (A-C13), reproduced every ATP vector and summary statistic of the pre-registered pass, while the figures reported here come from a re-run of the A-C13 configuration on the repaired data and a re-drawn sample (the \nameref{sec:sharedref} section).
>
> L1200: What survives the controls is narrower but still real: the choice of which amplicons to decline carries information, being worth $+0.021$ population F1 over declining as many at random ($+0.019$ for the exact-tie rule against its own masked run); a prFBA prediction names the genomes it came from, which an ancestral-state reconstruction cannot; and the EC deficit leaves no detectable mean difference in the ATP phenotype, although per genome it carries through weakly (Spearman $\rho = -0.19$ between the EC-F1 and ATP-distance differences, $p = 1.5\times10^{-7}$).
>

#### 4. A-C13 counts are a re-run on a new sample (abstract)

Old:

> L80: One pre-registered phenotype control failed (five of 200 gapfilled truth models gated out where zero are allowed), which under the design's own rule marked the leave-genus-out run invalid; after two disclosed post-hoc amendments made with knowledge of that failure (a template-direction correction with the gapfiller's reversals applied, A-C12, and a later ModelSEEDpy commit that fixes two further build defects, A-C13), both phenotype passes were re-run, the 1{,}000-genome ATP pass reproduced every ATP vector (43{,}140 of 43{,}140 rows), no gapfilled model of any arm is gated (0 of 8{,}614 rows, against 49), and the run is valid under the amendments, with no EC number changed.
>

New:

> L80: One pre-registered phenotype control failed (five of 200 gapfilled truth models gated out where zero are allowed), which under the design's own rule marked the leave-genus-out run invalid; after two disclosed post-hoc amendments made with knowledge of that failure (a template-direction correction with the gapfiller's reversals applied, A-C12, and a later ModelSEEDpy commit that fixes two further build defects, A-C13), both phenotype passes were re-run and the 1{,}000-genome ATP pass reproduced every ATP vector of the pre-registered pass (43{,}140 of 43{,}140 rows); in the A-C13 configuration, re-run on the repaired gene-family data with a re-drawn sample, no gapfilled model of any arm is gated (0 of 8{,}660 rows, against 49 in the pre-registered run), and the run is valid under the amendments, which change no EC number.
>

#### 4b. A-C13 counts (sec:sharedref)

Old:

> L667: The 1{,}000-genome ATP pass reproduces the pre-registered pass exactly (43{,}140 of 43{,}140 rows), so every ATP figure above stands: the corrections change reaction directions, not which reactions a draft holds or on which media its core model makes ATP.
>
> L668: In the gapfilled pass no model of any arm is gated (0 of 8{,}614 rows, against 49), including the 180 genomes outside the validation.
>
> L672: The ATP pass again reproduces every ATP vector and summary statistic; no gapfilled model is gated, every model grows, and 114 apply a reversal (107 under A-C12); the method arms' essentiality and growth medians move by less than 0.005 and their gapfill burden by at most one reaction, and none of the comparisons drawn above changes; and, rescored against this summary, the artifact remains valid.
>
> L673: The gapfilled figures above are A-C13's.
>
> L679: They reach no reported readout: the ATP vector is that of each draft's core model, whose empty-medium control makes no ATP in any of the 43{,}140 rows of the ATP pass, and in the gapfilled pass ATP correction closes the draft loops before gapfilling, so that under both amendments no gapfilled model makes ATP or biomass from water and protons alone (0 of 8{,}614 rows gated).
>

New:

> L674: On the sample of the time, the 1{,}000-genome ATP pass reproduced the pre-registered pass exactly (43{,}140 of 43{,}140 rows): the corrections change reaction directions, not which reactions a draft holds or on which media its core model makes ATP.
>
> L675: In the gapfilled pass no model of any arm was gated (0 of 8{,}614 rows, against 49), including the 180 genomes outside the validation.
>
> L679: On that sample the ATP pass again reproduced every ATP vector and summary statistic, no gapfilled model was gated, every model grew and 114 applied a reversal (107 under A-C12), the method arms' essentiality and growth medians moved by less than 0.005 and their gapfill burden by at most one reaction, and, rescored against that summary, the artifact remained valid.
>
> L680: After the gene-family cache was repaired, the A-C13 configuration alone was re-run, on the re-drawn sample: its ATP pass has 43{,}174 rows; in its gapfilled pass no model of any arm is gated (0 of 8{,}660 rows), all 2{,}372 gapfilled models grow and 107 apply a reversal; and, rescored against this summary, the leave-genus-out artifact is valid.
>
> L681: The phenotype figures above are this re-run's.
>
> L687: They reach no reported readout: the ATP vector is that of each draft's core model, whose empty-medium control makes no ATP in any of the 43{,}174 rows of the re-run ATP pass, and in the gapfilled pass ATP correction closes the draft loops before gapfilling, so that under both amendments, and in the re-run, no gapfilled model makes ATP or biomass from water and protons alone (0 of 8{,}660 rows gated in the re-run).
>

#### 4c. Gapfilled readouts (PH-06)

Old:

> L644: On every one of these readouts prFBA is indistinguishable from copying the nearest surviving reference: reaction Jaccard 0.707 against 0.705, essentiality disagreement 0.068 against 0.069, and median growth log-ratio against the truth model 0.008 against 0.008 (PICRUSt2 $-0.006$).
>
> L646: These gapfilled figures come from the re-run under amendments A-C12 and A-C13 (below), in which no model fails its sanity control, so all 200 genomes are scored; the pre-registered run, which excluded the five genomes whose truth model failed it, gave the same pattern (essentiality disagreement 0.066 against 0.068, growth log-ratio 0.002 against 0.004).
>

New:

> L651: On these readouts prFBA stays close to copying the nearest surviving reference (reaction Jaccard 0.704 against 0.700, essentiality disagreement 0.063 against 0.063), while its median growth log-ratio against the truth model is 0.020 against nn1's 0.011 (PICRUSt2 $-0.003$).
>
> L653: These gapfilled figures come from the re-run of amendment A-C13 on the re-drawn sample (below), in which no model fails its sanity control, so all 200 genomes are scored; the pre-registered run, which excluded the five genomes whose truth model failed it, gave the same essentiality ordering (0.066 against 0.068) and growth log-ratios near zero (0.002 against 0.004).
>

#### 5. Depth probe (abstract)

Old:

> L67: The coverage and whole-population self-exclusion figures remain lower bounds: under the shipped exact-tie union, 86.7\% of the abstentions arise because the exclusion filter empties a hit record already truncated to its top 100 hits; raising the record depth from 20 to 100 cut the abstention rate from 23.4\% to 9.5\% while leaving the conditional figures roughly unchanged (conditional EC F1 0.802 against 0.798); and re-aligning at persistence depth 500 a sample of records emptied at depth 20 recovers an answer for 88.4\% (283/320) of them at comparable quality (EC F1 0.783 vs.\ 0.798 for amplicons that already answered at depth 20, precision 0.747 vs.\ 0.790, genus accuracy 0.60).
>

New:

> L67: The coverage and whole-population self-exclusion figures remain lower bounds: under the shipped exact-tie union, 86.7\% of the abstentions arise because the exclusion filter empties a hit record already truncated to its top 100 hits; raising the record depth from 20 to 100 cut the abstention rate from 23.4\% to 9.5\% while leaving the conditional figures roughly unchanged (conditional EC F1 0.802 against 0.798); and re-aligning at persistence depth 500 a sample of records emptied at depth 20 recovers an answer for 88.4\% (283/320) of them, and in EC space those answers are no worse than the amplicons that already answered at depth 20 (EC F1 0.829 vs.\ 0.798, precision 0.791 vs.\ 0.790; genus accuracy 0.60).
>

#### 5b. Depth probe (sec:depth)

Old:

> L728: The recovered references are genus-tier (281 of 283) at a median surviving identity of 0.9856, just below the cutoff, and --- the point that decides whether the recovery is worth having --- they are of ordinary quality rather than dregs: EC recall 0.854, precision 0.747 and F1 0.783, against 0.835, 0.790 and 0.798 for the amplicons that already answered at depth 20, with 80.2\% family-correct and 60.1\% genus-correct anchors.
>
> L731: The conditional figures are not: under the same shipped selector, the probe's recovered amplicons score below the depth-20 answered set on EC precision (0.747 versus 0.790), F1 (0.783 versus 0.798) and genus anchor accuracy (0.601 versus 0.651 pooled), match it on family anchor accuracy (0.802 versus 0.800) and exceed it on EC recall (0.854 versus 0.835), which predicted conditional figures roughly unchanged or slightly lower.
>

New:

> L736: The recovered references are genus-tier (281 of 283) at a median surviving identity of 0.9856, just below the cutoff, and --- the point that decides whether the recovery is worth having --- they are not dregs: in EC space they score at or above the amplicons that already answered at depth 20 (recall 0.904, precision 0.791 and F1 0.829, against 0.834, 0.790 and 0.798), and their anchors are 80.2\% family-correct and 60.1\% genus-correct.
>
> L739: The conditional figures are not: under the same shipped selector, the probe's recovered amplicons score at or above the depth-20 answered set in EC space (precision 0.791 versus 0.790, F1 0.829 versus 0.798, recall 0.904 versus 0.834), match it on family anchor accuracy (0.802 versus 0.800) and fall below it on genus anchor accuracy (0.601 versus 0.651 pooled), which predicted conditional figures roughly unchanged: no lower in EC space, slightly lower at genus rank.
>

#### 5c. Depth probe (limitations)

Old:

> L1204: Conditional accuracies are not bounded this way: in the depth probe the recovered calls are 60.1\% genus-correct and 80.2\% family-correct with EC F1 0.783, at or below the answered set, and raising the record depth from 20 to 100 left the conditional figures roughly unchanged (as-shipped conditional EC F1 0.802 against 0.798; under the legacy reducer, the only rule with a depth-20 \Cref{tab:main} scoring, conditional genus accuracy within 0.01 in every region, and conditional gene capture, which is macro recall, up 0.7 to 1.3 points while conditional gene-capture F1 fell by 0.1 to 1.2 points).
>

New:

> L1214: Conditional accuracies are not bounded this way: in the depth probe the recovered calls are 60.1\% genus-correct and 80.2\% family-correct, below and level with the answered set, but score above it in EC space (EC F1 0.829 against 0.798); and raising the record depth from 20 to 100 left the conditional figures roughly unchanged (as-shipped conditional EC F1 0.802 against 0.798; under the legacy reducer, the only rule with a depth-20 \Cref{tab:main} scoring, conditional genus accuracy within 0.01 in every region, and conditional gene capture, which is macro recall, up 0.7 to 1.3 points while conditional gene-capture F1 fell by 0.1 to 1.1 points).
>

#### 6. Surviving-uncovered and CMH

Old:

> L742: A second, independently selected panel converts this argument into a direct measurement: on a panel built to be denser in near-identical relatives, with both panels scored under the exact-tie union at depth 100, the full-length 16S emptied-record rate rises from 4.7\% to 22.9\% (8.3\% to 29.4\% pooled over all eight regions) while the share of full-length records that survive exclusion yet go uncovered stays near 3\% on both (3.4\% and 2.8\%; 3.1\% and 2.4\% pooled over all eight regions), which accounts for the whole coverage difference between them (the \nameref{sec:twopanel} section).
>
> L1059: Where the record survives, only a few percent of each panel's queries go uncovered (3.4\%, 2.8\% and 2.3\% of surviving full-length records for the original, $k$-means and random panels, and 1.5--3.7\% in every full-length stratum; 3.1\%, 2.4\% and 1.8\% pooled over all eight regions), some 29- to 43-fold below the emptied-record rate on full-length 16S.
>
> L1060: These full-length surviving-record rates differ by about a point (1.1 points, and 1.3 pooled over all eight regions) and in the direction opposite to the coverage gap: the original panel, the best covered, leaves the most surviving records uncovered; stratified by the number of deleted hits, only the original-versus-random difference is significant (Cochran--Mantel--Haenszel $p=0.004$; original versus $k$-means $p=0.15$).
>
> L1062: The cross-panel coverage difference is, to within two points, the difference in how often that happens: 33.4\% of random-panel queries, 22.9\% of $k$-means and 4.7\% of the original's. The arithmetic closes --- $1-0.334=0.666$ against an observed 0.650, $1-0.229=0.771$ against 0.750, and $1-0.047=0.953$ against 0.921, the residual being the few percent of surviving records that still go uncovered.
>
> L1063: (The ``uncovered'' columns report $1-{}$coverage: gene-capture abstentions plus cells whose selected references carry no gene families, of which the original, $k$-means and random panels have 107, 76 and 62 among their surviving full-length records.)
>

New:

> L750: A second, independently selected panel converts this argument into a direct measurement: on a panel built to be denser in near-identical relatives, with both panels scored under the exact-tie union at depth 100, the full-length 16S emptied-record rate rises from 4.6\% to 22.9\% (8.2\% to 29.4\% pooled over all eight regions) while the share of full-length records that survive exclusion yet go uncovered stays near 2\% on both (2.4\% and 2.3\%; 1.8\% and 1.8\% pooled over all eight regions), which accounts for the whole coverage difference between them (the \nameref{sec:twopanel} section).
>
> L1068: Where the record survives, only about 2\% of each panel's queries go uncovered (2.4\%, 2.3\% and 1.8\% of surviving full-length records for the original, $k$-means and random panels, and at most 3.0\% in any full-length stratum, down to none of the original panel's 274 records in the 50--99 stratum; 1.8\%, 1.8\% and 1.3\% pooled over all eight regions), some 42- to 55-fold below the emptied-record rate on full-length 16S.
>
> L1069: These full-length surviving-record rates differ by about half a point (0.6 points, and 0.5 pooled over all eight regions) and in the direction opposite to the coverage gap: the original panel, the best covered, leaves the largest share of surviving full-length records uncovered, by a fraction of a point (pooled over regions it ties the $k$-means panel); stratified by the number of deleted hits, neither difference is significant (Cochran--Mantel--Haenszel $p=0.117$ for original versus random, $p=0.99$ for original versus $k$-means).
>
> L1071: The cross-panel coverage difference is, to within 2.3 points, the difference in how often that happens: 33.3\% of random-panel queries, 22.9\% of $k$-means and 4.6\% of the original's.
>
> L1072: The arithmetic closes --- $1-0.333=0.667$ against an observed 0.655, $1-0.229=0.771$ against 0.754, and $1-0.046=0.954$ against 0.931, the residual being the roughly 2\% of surviving records that still go uncovered.
>
> L1073: (The ``uncovered'' columns report $1-{}$coverage: gene-capture abstentions plus cells whose selected references carry no gene families, of which the original, $k$-means and random panels have 36, 42 and 31 among their surviving full-length records.)
>

#### 7. Random within 0.016 -> 0.017

Old:

> L1019: In crude terms the random panel sits below the original in every amplicon region, by at most 0.029 (V7--V9), and within 0.016 at every MAG completeness; the $k$-means panel sits within 0.023 of the original in the amplicon arm and 0.010 in the MAG arm, the original ahead in all eight regions and at all five completeness levels.
>

New:

> L1028: In crude terms the random panel sits below the original in every amplicon region, by at most 0.030 (V7--V9), and within 0.017 at every MAG completeness; the $k$-means panel sits within 0.023 of the original in the amplicon arm and 0.010 in the MAG arm, the original ahead in all eight regions and at all five completeness levels.
>

#### 8. MAG multiplicity strata

Old:

> L1087: Abstention there falls to 0.559 for organisms whose species has 20 or more members in the panel, and that stratum accounts for the $k$-means panel's lower overall MAG abstention (0.658 against 0.688 over organisms with a species-level truth lineage); without it the $k$-means panel abstains at 0.687.
>

New:

> L1097: Abstention there falls to 0.559 for organisms whose species has 20 or more members in the panel, and that stratum accounts for the $k$-means panel's lower overall MAG abstention (0.657 against 0.682 over organisms with a species-level truth lineage); without it the $k$-means panel abstains at 0.686, slightly above the original.
>

#### 9. Cache-repair disclosure (L516-517)

Old:

> L516: The 953 amplicons (from 161 source genomes) whose source genome has no gene-family annotation in the BV-BRC cache, and therefore no truth EC set, are unscoreable and excluded.
>
> L517: A later audit of that cache found that the empty family lists of most of these source genomes (behind 675 of the 953 amplicons), and of the references behind 818 of the 1{,}060 self-excluded amplicons scored as empty predictions, were stored failed fetches rather than BV-BRC's annotation, and that 167 genomes' lists had been cut short by the 25{,}000-row response limit; re-scored on the repaired cache, the 10{,}000-organism panel's whole-population gene capture rises by 0.75--0.86 points under self-exclusion, and capture over resolved amplicons moves by at most 0.05 points, so the gene-capture figures reported here are slightly conservative.
>

New:

> L516: The 278 amplicons (from 38 source genomes) whose source genome has no gene-family annotation in the BV-BRC cache, and therefore no truth EC set, are unscoreable and excluded.
>
> L517: Before the cache was repaired on 2026-11-13, this exclusion held 953 amplicons from 161 source genomes: an audit found that the empty family lists behind 675 of those 953 amplicons, and behind 818 of the 1{,}060 self-excluded amplicons then scored as empty predictions, were stored failed fetches rather than BV-BRC's annotation, and that 167 genomes' lists had been cut short by the 25{,}000-row response limit.
>
> L518: Every result reported here was recomputed on the repaired cache, which holds only BV-BRC's own annotation fetched past that limit (the \nameref{sec:recon} section); the exceptions, which predate the repair and are labeled where they appear, are the backend concordance check (the \nameref{sec:concordance} section), the shared-reference panel's eligibility (the \nameref{sec:sharedref} section) and the pre-registered and first amended phenotype passes.
>
> L519: In this comparison the repair made 675 more amplicons scoreable (67{,}713) and turned most of those empty predictions into answers (243 remain), raising prFBA's coverage from 0.889 to 0.902 and its population EC F1 from 0.713 to 0.723 while leaving its conditional EC F1 at 0.802; on the 10{,}000-organism panel it raised whole-population gene capture under self-exclusion by 0.7--1.1 points by region and moved capture over resolved amplicons by at most 0.03 points.
>

#### 10. Gene-capture denominators (L397, L419)

Old:

> L397: We report two denominators for the self-exclusion figures: a \emph{conditional} accuracy over amplicons the pipeline still resolves to a reference named at that rank (for gene capture, to any reference that carries gene families), and a \emph{whole-population} accuracy in which every other amplicon counts as unresolved: the amplicons that receive no taxonomic call once their species-level relatives are removed (9.5\% of all amplicons, 10.3\% of genus-scoreable ones) and, at genus, a further 22\% of genus-scoreable amplicons whose selected reference has no named genus (genus coverage 0.64--0.74 by region); in gene capture, the 11.1\% of scoreable amplicons for which the selection keeps no genome (9.6\%) or keeps only references that carry no gene families (1.6\%) score zero.
>
> L419: Gene capture moves the same way: the share of scoreable amplicons for which the selection keeps no genome falls from 9.6\% to 6.0\% (and those keeping only references without gene families from 1.6\% to 1.4\%), so whole-population capture rises by 4.1--5.5 points in every region (to 62.5--69.0\%; V4 from 58.1\% to 63.1\%) and capture over resolved amplicons by 2.4--3.0 points (to 69.1--72.9\%).
>

New:

> L397: We report two denominators for the self-exclusion figures: a \emph{conditional} accuracy over amplicons the pipeline still resolves to a reference named at that rank (for gene capture, to any reference that carries gene families), and a \emph{whole-population} accuracy in which every other amplicon counts as unresolved: the amplicons that receive no taxonomic call once their species-level relatives are removed (9.5\% of all amplicons, 10.3\% of genus-scoreable ones) and, at genus, a further 22\% of genus-scoreable amplicons whose selected reference has no named genus (genus coverage 0.64--0.74 by region); in gene capture, the 9.8\% of scoreable amplicons for which the selection keeps no genome (9.5\%) or keeps only references that carry no gene families (0.4\%) score zero.
>
> L419: Gene capture moves the same way: the share of scoreable amplicons for which the selection keeps no genome falls from 9.5\% to 6.0\% (the share keeping only references without gene families stays at 0.4\%), so whole-population capture rises by 4.1--5.4 points in every region (to 63.4--69.7\%; V4 from 59.3\% to 64.0\%) and capture over resolved amplicons by 2.4--3.0 points (to 69.0--72.9\%).
>

#### 11. As-shipped deficit magnitude (+ legacy)

Old:

> L525: Scored over the whole population, where prFBA's abstentions count as zero, the ranking reverses on the mean, though only narrowly (F1 0.713 versus 0.729), while the median amplicon still favors prFBA --- the honest statement of the tradeoff, not a contradiction.
>
> L1187: Sharing the reference and the vocabulary and blocking holdouts at genus rank reverses the as-shipped conditional margin, from $+0.075$ on matched amplicons (0.802 against 0.727) to $-0.027$ on matched cells (0.769 against 0.796; $-0.029$, 0.767, for the exact-tie rule), and widens the population deficit, from $-0.016$ under the shipped exact-tie union (0.713 against 0.729; $-0.129$ at the earlier record depth of 20) to $-0.092$ ($-0.094$ for the exact-tie rule, the same rule at both ends); the controlled design also re-aligns each fold against its own restricted reference and retains 200 hits per query (\texttt{-{}-k2 200}).
>
> L1216: Under it prFBA's conditional advantage reverses (matched conditional EC F1 0.769, or 0.767 for the exact-tie rule, against 0.796, from a matched 0.802 against 0.727 as shipped), while PICRUSt2's population lead persists and widens ($\Delta$ from $-0.016$ to $-0.092$, or to $-0.094$ for the exact-tie rule).
>

New:

> L527: Scored over the whole population, where prFBA's abstentions count as zero, the ranking reverses on the mean under the shipped exact-tie union, though only narrowly (F1 0.723 versus 0.729), while the median amplicon still favors prFBA --- the honest statement of the tradeoff, not a contradiction.
>
> L528: Under the legacy reducer, the policy the controlled comparison below selects, the mean does not reverse (population F1 0.734).
>
> L1197: Sharing the reference and the vocabulary and blocking holdouts at genus rank reverses the as-shipped conditional margin, from $+0.075$ on matched amplicons (0.802 against 0.727) to $-0.027$ on matched cells (0.769 against 0.796; $-0.029$, 0.767, for the exact-tie rule), and widens the population deficit, from $-0.006$ under the shipped exact-tie union (0.723 against 0.729; $-0.120$ at the earlier record depth of 20) to $-0.092$ ($-0.094$ for the exact-tie rule, the same rule at both ends; for the legacy policy at both ends, from a $+0.005$ lead as shipped to $-0.092$); the controlled design also re-aligns each fold against its own restricted reference and retains 200 hits per query (\texttt{-{}-k2 200}).
>
> L1226: Under it prFBA's conditional advantage reverses (matched conditional EC F1 0.769, or 0.767 for the exact-tie rule, against 0.796, from a matched 0.802 against 0.727 as shipped), while PICRUSt2's population lead under the shipped exact-tie rule persists and widens ($\Delta$ from $-0.006$ to $-0.094$), and the legacy policy goes from a $+0.005$ lead as shipped to $-0.092$.
>

#### 12. Three fifths (discussion)

Old:

> L1175: Third, the 10{,}000-organism validation exposes a region-versus-generalization tradeoff --- one that follows the variable region an amplicon spans rather than its length, since V1--V2, the second-shortest amplicon, generalizes about as well as V1--V3 --- and resolves it at two taxonomic ranks: when the source organism itself is in the database, even a short marker recovers nearly all of its metabolism (over 99\% gene capture, at macro precision 0.79--0.88 for the regions from V3 onward and 0.92--0.94 for V1--V2 and V1--V3), whereas for genuinely novel organisms an amplicon from V3 onward (V4, V4--V5, V3--V4, V7--V9) is right at the family level more than three times in four when it makes a call (0.76--0.80 conditional) but, over the self-excluded population, receives no genus call for roughly a third (abstentions plus answered amplicons whose reference has no named genus), captures less than three fifths of the gene complement (57.6--59.2\%) and places 54--59\% of organisms in the correct family and 37--42\% in the correct genus --- a gap that full-length 16S and the two V1-anchored regions (V1--V2, V1--V3) measurably narrow, to 63.3--64.4\% gene capture and 0.51--0.54 population genus accuracy.
>

New:

> L1185: Third, the 10{,}000-organism validation exposes a region-versus-generalization tradeoff --- one that follows the variable region an amplicon spans rather than its length, since V1--V2, the second-shortest amplicon, generalizes about as well as V1--V3 --- and resolves it at two taxonomic ranks: when the source organism itself is in the database, even a short marker recovers nearly all of its metabolism (over 99\% gene capture, at macro precision 0.79--0.88 for the regions from V3 onward and 0.92--0.94 for V1--V2 and V1--V3), whereas for genuinely novel organisms an amplicon from V3 onward (V4, V4--V5, V3--V4, V7--V9) is right at the family level more than three times in four when it makes a call (0.76--0.80 conditional) but, over the self-excluded population, receives no genus call for roughly a third (abstentions plus answered amplicons whose reference has no named genus), captures about three fifths of the gene complement (58.6--60.3\%) and places 54--59\% of organisms in the correct family and 37--42\% in the correct genus --- a gap that full-length 16S and the two V1-anchored regions (V1--V2, V1--V3) measurably narrow, to 64.0--65.2\% gene capture and 0.51--0.54 population genus accuracy.
>

#### 13. Limitations: amendments

Old:

> L1222: On the phenotype side, the absence of a detectable mean ATP difference is a null at 755 paired genomes rather than a demonstration of equivalence, the ATP vector carries signal in only 22 of its 54 media under this template route (in the truth models, in both methods' arms and in the nearest-neighbor and family-donor baselines; no calibration rung switches on any other medium), the assay is restricted to the subset prFBA answers, and the gapfilled growth and essentiality readouts are reported only as unpaired per-arm summaries with no paired test: they are scored on a 200-genome gapfilled subsample of the 1{,}000 (controlled pipeline: growth and essentiality $n=200$ for both methods); and the validity of the leave-genus-out artifact rests on post-hoc amendments: the pre-registered gapfill sanity control failed its zero-tolerance leg, with five truth models gated out (one gapfill per genome, shared by both pipelines, so 10 gated rows), and it passes only after the template-direction correction diagnosed in the \nameref{sec:sharedref} section was adopted and both phenotype passes were re-run, and again after a later ModelSEEDpy fixed two build defects (A-C13), changes made with knowledge of the failure, although they move no EC number and no ATP figure.
>

New:

> L1232: On the phenotype side, the absence of a detectable mean ATP difference is a null at 759 paired genomes rather than a demonstration of equivalence, the ATP vector carries signal in only 22 of its 54 media under this template route (in the truth models, in both methods' arms and in the nearest-neighbor and family-donor baselines; no calibration rung switches on any other medium), the assay is restricted to the subset prFBA answers, and the gapfilled growth and essentiality readouts are reported only as unpaired per-arm summaries with no paired test: they are scored on a 200-genome gapfilled subsample of the 1{,}000 (controlled pipeline: growth and essentiality $n=200$ for both methods); and the validity of the leave-genus-out artifact rests on post-hoc amendments: the pre-registered gapfill sanity control failed its zero-tolerance leg, with five truth models gated out (one gapfill per genome, shared by both pipelines, so 10 gated rows), and it passes only after the template-direction correction diagnosed in the \nameref{sec:sharedref} section was adopted and both phenotype passes were re-run, and again after a later ModelSEEDpy fixed two build defects (A-C13), changes made with knowledge of the failure, although they move no EC number and no ATP figure; the A-C13 configuration was then re-run on the repaired data with a re-drawn sample, and it still passes.
>

---------------------------------------------------------------------------------------------------------------

## Disclosures added or rewritten (step 3)

- **Path A (D1)** after the "primary analysis" sentence (new L591–592): panel eligibility was decided 2026-10-13 from
  the cache of the time.
  - 123 of the 161 empty-EC exclusions were stored fetch failures (122 pass the filters, 65 of those archaea).
  - Panel, folds, trees, placements and alignments are kept; re-planning would have admitted about 115 more genomes
    and changed every fold.
  - Truth and traits are recomputed. This changes the EC sets of 17 eligible genomes, which reach 23 of 30 genus
    references and 130 of 63,871 test cells. The 130 were re-counted from `fold_plan.json` × `changed_sets.json`:
    24 / 47 / 24 / 31 / 4 by fold. The 23 of 30 is from RUN3.
  - n_eff 6,716 → 6,715 at all three places.
- **Cache-repair disclosure** (old L516–517 → new L516–519). The old "slightly conservative / +0.75–0.86 / ≤ 0.05 points"
  passage is replaced by:
  - a past-tense correction note (953 amplicons from 161 genomes before the repair; 675 and 818 stored fetch failures;
    167 truncated lists);
  - a statement that every result that reads the gene-family cache was recomputed on the repaired, BV-BRC-only
    cache, naming the exceptions that predate it: the D7 concordance check, the panel's eligibility (D1), and the
    phenotype passes and side studies run on the earlier phenotype sample (wording as of `2387ac2`);
  - what the repair changed: in the as-shipped comparison +675 scoreable amplicons, empty predictions 1,060 → 243,
    coverage 0.889 → 0.902, population EC F1 0.713 → 0.723, conditional F1 still 0.802. On the panel, population gene
    capture rises 0.7 to 1.1 points by region and resolved capture moves at most 0.03 points (computed from the old
    and new `gene_capture_by_region.json`).
- **Phenotype re-draw** (new L638–639, after the subsample sentence):
  - the subsample is stratified by domain × template × quintile of the EC-F1 difference;
  - 200 of 7,855 candidates change, moving 13 of 15 strata;
  - 442 of 1,000 ATP genomes and 62 of 200 gapfilled genomes are kept;
  - a kept genome's rows are identical unless its own inputs changed (106 of 19,072 ATP rows; 0 of 2,684 gapfilled
    rows).

  Source: RUN4_LOCAL §3 and C12_SCOPE §7b.
- **D6** (new L1019, after the random panel's panel-as-sketch sentence):
  - 1733.9197 and 1733.9535 returned exactly 25,000 contigs: 10 of the arm's 49,990 MAG queries;
  - the arm was not re-fetched;
  - two genomes of the real-MAG arm's BV-BRC subset (1733.9522, 1733.9680) are capped too, but no real MAG hits them
    (`pipelines_mag_deploy.md`). That is the plan's "4 genomes".
- **D7** (L385): "This check predates the shipped rule, the record depth and the repair of the gene-family cache: …
  on the gene-family cache of the time".

---------------------------------------------------------------------------------------------------------------

## Not recomputed / unresolved (kept as printed, or dropped)

| claim | where | decision |
|---|---|---|
| AF-gate "1,386 newly admitted incl. 19 gene-less" | L879 | **dropped** (not recomputed; the store has only the non-empty count 1,391). The sentence now starts "Of the organisms a 2% gate would newly admit, the 1,391 that receive a non-empty prediction …"; 26.7% → 26.5% and 0.710 → 0.712 (Run 2 values in number_diff's manual note) |
| depth probe's answered-set anchors 0.651 (genus) / 0.800 (family) | L731 | kept as printed; number_diff marks them "not recomputed, expected unchanged (d20 exclude-mode taxonomy unchanged)" |
| depth-500 projections ~0.96 / ~0.49 / ~0.75 | L734 | ~0.96 and ~0.49 kept (taxonomy-based, unchanged inputs). The EC projection was updated **~0.75 → ~0.77**, the value the re-numbered P2C README gives; the paper's own recipe (base + 0.221 × 0.884 × recovered F1) reproduces both the old 0.753 and the new 0.771 |
| AS-15 rxn-vs-EC divergences 0.026 / 0.014 / 0.05 / 0.03 | L527 | kept (not recomputed) |
| SR-02 inner margins 0.001–0.003 | L340, L590 | kept (not recomputed) |
| wobj "weighted term median 0.11%" | L310 | kept (not recomputed; the wobj rows moved in 2 of 360 jobs) |
| historical figures: the pre-registered pair 0.066 / 0.068 and 0.002 / 0.004; 49 gated rows of 8,614; 136 models in 50 genomes; 107 reversals under A-C12; 166 of 166; 865 → 95 of 1,181; 138 of 200; 218 models of 20 genomes; 72 of 72 | L646–683 | kept, now labelled as the pre-registered run or the earlier sample |

---------------------------------------------------------------------------------------------------------------

## Decks (R5.3, subagent; commit `a6812be`)

**Constants (old → new, source):**
- `AF_SWEEP`: 0.311/0.316/0.379/0.452/0.567/0.655 → 0.314/0.319/0.381/0.454/0.569/0.656 (`deck_AF_SWEEP_*`).
- `BAND_MAG`: 0.082/0.202 → 0.081/0.200. The depth-100 bins are used and a comment says so; both print as 8–20%.
- `MAG_REACH` and `XPHYLUM_BENCH`: unchanged.
- `EMPTIED` (`deck_constants_d100.json`): 4.7/22.9/33.4 → 4.6/22.9/33.3.
- `PHENO_EC_GAP`: −0.0284 → −0.03133.
- `AF2_NEW` / `AF2_BEAT_FL`: 1,367 / 0.599 / 26.7% → 1,391 / 0.600 / 26.5%, now read from
  `af_gate_counterfactual.json`. The old table had no entry for the new 0.7116 level, so the build failed until this
  was changed.
- `PFE_WITHIN_UNION_F1`: upper bound 0.0096 → 0.0097, now read from `p4_results.json`.
- `OTUNE_LOG_CI` and `INCL_N_SELECTED`: unchanged; both are now read from files.
- `GAP_CI_D100`: computed at build time with a seeded bootstrap. Union [−0.0110, −0.0015]; legacy [+0.0006, +0.0101].
- E3 "about 45%" → 48%.
- `KBASE_D100`: → 261948/6468/1 and 275568/3758/1, dated 2026-11-15.
- `INTEGRATION_STATUS`: updated.

**Sentences rewritten:**
- Slide 25 (condensed 13, brief 8), phenotype:
  - nn1 "−0.15, q = 0.028" → "−0.03, q = 0.51, n.s." (no longer in the win colour);
  - McNemar "2 of 22 (Ac.NO, Etho.O2; q = 0.002)" → "4 of 22 (Ac.NO, Ac.O2, Etho.O2, For.O2; q 0.029–0.046)";
  - re-drawn sample and re-run counts;
  - "49 in the pre-registered run".
- Slide 27: the ATP pass now reads "re-run", and the 10-29 comparison is labelled.
- Slide 29: the seven-reaction study is labelled "pre-repair sample".
- Slides 19 and 30: the legacy gap is no longer "within resampling noise".
- Slide 31: the legacy cells now read "(significant)", and a union-costs-population-F1 sentence was added.
- Slide 42: CMH not significant (p 0.117 / 0.99).
- KBase sentences on slides 1, 2, 6 and 45 and in the condensed and brief decks.
- E3 gate: "passed again 2026-11-15".

**Rendering and files:**
- All four decks were rendered (soffice → PDF → PNG) and every changed slide was viewed: full deck 5–45, condensed
  8 and 10–20, brief 7–10, brief-panels 6–10. No overflow or overlap.
- The commit also contains `presentations/condensed/F_heatmap.png` and `figure_stats.json`. The condensed decks embed
  that heatmap, so `make_condensed_figures.py` was re-run; the other condensed figures came out byte-identical.

**Deck items not re-scored:**
- The slide-6 notes' deployment set changes (149/130/129), a 10-30 comparison not recomputed after the 11-13
  rebuild;
- hit-level MAG/GTDB/source-recovery shares (not gene-family dependent);
- slides 9–10 notes still say "no pull request opened yet" while the face says PR #20 is open. This predates the
  re-score.

**Judgment calls the user may overrule:**
- the closing slide's Done/Open items ("Re-scored on repaired BV-BRC data"; "Manuscript: re-scored revision" now
  open);
- one re-score sentence in each title-slide's notes.

I spot-checked the rebuilt decks' text and notes for the stale strings (0.146, two media, q = 0.002, near 3%, 0.783,
"slightly conservative", 43,140, 0.713, 67,038, 0.308, "within 0.016"). The only survivors are labelled historical
statements: the 10-29 C12 comparison's 8,614, the pre-repair seven-reaction study, and deployment counts that merely
look similar (3,755).

---------------------------------------------------------------------------------------------------------------

## Flags for the user

1. **Under the legacy reducer prFBA now leads PICRUSt2 as shipped.** The deck agent found this; number_diff does not
   track it.
   - Population EC F1 0.734 vs 0.729 (+0.005). The deck build's bootstrap CI is +0.0006 to +0.0101; the paper cites
     only the point values. It was −0.002 [−0.007, +0.002] before the repair.
   - The manuscript's as-shipped arm is the exact-tie union (−0.006), so its claims hold under that rule.
   - Because the controlled comparison's primary arm is the legacy policy, I added a short clause in three places
     (new L527–528, L1197, L1226): under legacy the as-shipped mean does not reverse, and it goes from a +0.005 lead to
     −0.092. Only the point values 0.734 and +0.005 are cited, read from `comparison_summary.json`
     (`prfba_exclude_legacy`). The §select enumeration of non-default-rule results now lists it (`2387ac2`).
2. **Repeats outside the claim inventory.** The inventory's line ranges missed L402, L1162, L521, L891, L899, L624
   and L527 (listed above). All are fixed. Any other prose repeat that no inventory item lists, and that my old-value
   sweeps did not catch, would still be stale. The sweeps covered every changed value's old string over the whole file.
3. **Plan errata still open** (plan §5, not re-score matters, not edited):
   - "460,117" (L61, L135; the planner counts 460,116);
   - "1 to 5,265 contigs" (L814; 813.301 holds 15,664);
   - marker reach 6,011 (W6 still counts 6,011);
   - F16 vs three-panel deciles.

   Each needs its own check before editing.
4. **Pre-existing doc errors from Run 5A** (not in the paper): the method paper's Table 6 V7–V9 Δ (+0.193 vs
   +0.167), "~112k" genomes with cached PGFams (the cache holds 159,753), and RESULTS_v2 "73x" vs the deck's 14x.
5. **Not committed in the paper repo:** the untracked `k_md5.txt` and `*.bak-*` files (pre-existing). The build
   products (`prFBA.pdf`, `.aux`, `.log`, `.bbl`) are git-ignored and were left in place.

Work files: `scratchpad/r5b/` in the session scratchpad. It holds `tables.py`, `tables_mag.py`, `mag_unrounded.json`,
`locate.py`, `apply.py`, `applied.json`, `rewrite.py`, `rewrites.json`, `sweep*.txt`, `meaning_pairs.md` and the
baseline build in `build0/`.
