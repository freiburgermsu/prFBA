# Exact-tie union throughout: change list for prFBA.tex

Written 2026-10-27. Read-only on both repositories; this file is the only one written into a repo.
Manuscript: `/home/freiburger/Documents/prFBA_paper/prFBA.tex` (branch `depth100-standard`, HEAD bcc65d6, 1,135 lines).

Commit bcc65d6 (17:04) landed during this work. It rewrote alignment and hardware sentences (L63, L146, L153, L297–298, L321, L639), none of which is rule-dependent. All substring checks below were run after it, against the current files.

## How to use this list

- Every item gives the current line, the exact current substring, the replacement, and the source of each new number with its unrounded value.
- Apply by substring, not by line number. Two replacements add a line (1.9 and 3.26), which shifts later line numbers by one.
- Apply section 1 in order, then section 3. Section 1b is optional; apply it last.
- Repo paths are relative to `/home/freiburger/Documents/prFBA` unless absolute.
- **U100**: the depth-100 exact-tie union trees (`region_validation_d100`, `region_validation_kmeans_d100`, `region_validation_random_d100`).
- **L100**: the same trees with `_legacy`.

## How it was checked

- **Values recomputed.** Every value was recomputed with `~/Documents/py_venv/bin/python` from the U100 files.
  - Before reading a U100 value, I reproduced the L100 number the manuscript now prints from the same key path in the `_legacy` file.
  - All of them reproduced: tab:main, tab:magvs16s, tab:magval, tab:twopanel, tab:threepanel, tab:emptied, McNemar, Welch, the tier capture, the single-genome shares, τ=99, the AF-gate 31.8%, and the controlled arms.
  - This catches wrong key paths before they propagate.
- **Rounding from the source.** Where a 4-dp summary value sits on a half (for example 0.7455, 0.8035, −0.2475, −0.0195, −0.1235), the value was recomputed from the cells.
- **Substrings tested by script** (`verify.py`, reproduction section). The tested substrings are the ones printed here, because this file was generated from the same data. Each "Current" substring:
  - occurs exactly once in `prFBA.tex`;
  - starts on the stated line;
  - contains no newline;
  - occurs exactly once in `prFBA_andrew_voice.tex`. The two files differ only at lines 121–124 and 238, and no item touches those lines, so every edit applies unchanged to the voiced variant.
- **Applied in sequence.** All section-1 and section-3 replacements apply in order without overlapping, and the optional items still match afterwards.
- **Compiled.** The fully edited copy (sections 1, 3 and 1b applied) compiles with pdflatex + bibtex:
  - 0 errors and no undefined references;
  - 84 pages against 82 for the current file;
  - one new 13 pt overfull vbox, from float placement.

## Counts

| Part | Count |
|---|---|
| Section 1: text, table and caption changes | 108 items, 22 of them qualitative |
| Section 1b: optional | 4 items |
| Section 2: figures | 6 must be re-rendered from union inputs; 3 need no change |
| Section 3: add the exact-tie arm alongside the controlled benchmark | 33 items |
| Section 4: unresolved | 9 points |

## Qualitative changes in brief

- **Self-inclusion.**
  - Genus accuracy falls to 0.89–1.00 and gene-capture recall rises to 99–100%.
  - Gene-capture precision falls on short amplicons (V4 0.795 against 0.856 under legacy).
  - "Near-ceiling, falling only at species rank" no longer holds for short regions: V4 anchor phylum is 0.955, family 0.922 and genus 0.891.
- **Region ordering in gene capture changes.**
  - V1–V2 and V1–V3 are now best, not full-length 16S; V4–V5 is now worst, not V4.
  - The gene-capture leakage gap is no longer near-uniform: 29.6–32.9 points, against 25.9–26.5.
- **McNemar.**
  - 26 of 28 pairs now differ, not 24.
  - The non-differing pairs are V1–V2 vs V1–V3 and V3–V4 vs V7–V9.
  - Full-length vs V1–V3 now differs, and so does the V4–V5 primer pair, which is no longer the closest pair by accuracy.
- **Single-genome share under self-exclusion** rises to 68.1%, from 46.8% (taxonomic) and 51.0% (gene content).
- **References that carry no gene families.**
  - The union creates 1,060 such exclude cells (`pred_empty`); legacy created none.
  - The conditional and population capture definitions must say how they are counted.
  - The tab:emptied note "none" becomes 107, 76 and 62.
- **Three panels.**
  - Conditional genus accuracy is no longer higher on the random panel in every region: the difference is −0.016 to +0.062, and the CI excludes zero only for full-length 16S and V1–V3.
  - The family bound rises from 0.04 to 0.05, and the gene-capture F1 bound from 0.027 to 0.028.
  - No panel has a uniform family advantage: the family-rank conditional differences are mixed in sign. The CI excludes zero only for V1–V3, where the alternative panels are higher (+0.042 random, +0.052 k-means), and V3–V4, where the original is higher (−0.025, −0.032); the other six regions cross zero.
  - The self-inclusion panel penalty more than doubles (V4 −0.107; genus 0.07–0.08).
  - "Below 3%" and "under a point" for the surviving records fail (3.4%; 1.1 points).
- **MAG vs 16S.**
  - The advantage is +0.10 to +0.21. The smallest region is now full-length 16S, not V1–V2.
  - At 1% completeness the full-length lead is +0.077, with 1,129 wins to 473.
  - The amplicon arm leads the population view by more (0.605 against 0.248).
- **Retention depth.** Three legacy depth-20 → depth-100 clauses have no union counterpart and are deleted (section 4).
- **PICRUSt2 as shipped.**
  - This section was already union; only the legacy parenthetical is deleted.
  - The regions where prFBA leads PICRUSt2 in population F1 are unchanged (full-length 16S, V1–V2, V1–V3).
- **Controlled benchmark.**
  - The exact-tie arm sits 0.002 below the nested-CV arm at the same coverage, and adding it changes no conclusion.
  - The nested-CV family headline −0.124 is a double rounding of −0.123459. It is flagged in section 4, not changed.

## 1. Text changes: legacy-sourced claims and their union replacements

108 items, in manuscript order. Apply in this order. Items marked QUALITATIVE change what a sentence claims, not just a number.

### 1.1. L63 (Abstract)

Current:

~~~latex
($+0.11$ to $+0.19$ F1 at full completeness, narrowing to $+0.08$ to $+0.17$ at 1\% completeness
~~~

Replacement:

~~~latex
($+0.10$ to $+0.21$ F1 at full completeness, narrowing to $+0.08$ to $+0.19$ at 1\% completeness
~~~

Source and values:

- Paired MAG-vs-16S conditional gene-capture F1, union amplicon cells (`region_validation_d100/data/gene_capture_cells.csv`, self_mode exclude_0.987) x MAG cells (`mag_validation/out/mag_cells_all.csv` + `lowc_out_batches/b*/mag_cells.csv`, regime exclude, contamination 0), compare_arms.py/crossover.py pairing logic recomputed at full precision.
- 100%: min FullLength16S +0.102997, max V4 +0.212059 (8 regions). `paired_comparison_d100.json` `exclude|c1_x0|FullLength16S|f1`.conditional mag 0.813 - amp 0.71 = 0.103; `exclude|c1_x0|V4|f1` 0.8101 - 0.598 = 0.2121.
- 1%: min FullLength16S +0.076708, max V4 +0.190273; these are the two regions `crossover_d100.json` carries (row completeness 0.01: `FullLength16S_delta` 0.0767, `V4_delta` 0.1903). Positive in every region at every one of the 14 levels (min over all levels +0.0767, max +0.2145).

Note: L100 check: 100% +0.1053 (V1-V2) to +0.1917 (V4); 1% +0.0815 (V1-V2) to +0.1692 (V4) reproduce the current text. Under the union the smallest gap moves from V1-V2 to full-length 16S.

### 1.2. L63 (Abstract)

Current:

~~~latex
scored with the earlier ordered (legacy) reducer except for the as-shipped figures, which use the shipped exact-tie union.
~~~

Replacement:

~~~latex
scored with the shipped exact-tie union selection rule.
~~~

Source and values:

- Rule statement; every 10,000-organism figure now comes from the union trees (`region_validation_d100`, `_kmeans_d100`, `_random_d100`; `_meta.select_mode` = `default_union`).

### 1.3. L64 (Abstract) — QUALITATIVE

Current:

~~~latex
When the source organism is in the database, genus accuracy is 0.94--0.99 and gene capture 95--99\%.
~~~

Replacement:

~~~latex
When the source organism is in the database, genus accuracy is 0.89--1.00 and gene capture 99--100\%.
~~~

Source and values:

- `region_validation_d100/data/taxacc_region_summary.json` regions[r].anchor.include.all.all.per_rank.Genus.correct: min V4 0.890933, max FullLength16S 0.996484.
- `region_validation_d100/data/gene_capture_by_region.json` _meta.headline_macro_recall_pct.include: min V4 99.148, max FullLength16S 99.979.

Note: QUALITATIVE: recall rises but the larger include-mode union costs precision: include macro precision is 0.7945 on V4 (0.8565 under legacy), 0.9967 on full-length (regions[r].include.macro.precision.value). See optional item O1.

### 1.4. L65 (Abstract) — QUALITATIVE

Current:

~~~latex
among amplicons resolved to a reference named at genus, genus accuracy is 0.60--0.74, and among amplicons resolved to any reference, gene capture is 69--73\%, while over the whole novel-organism population (for genus accuracy, abstentions and references with no named genus scored as unresolved; for gene capture, abstentions scored as zero) genus accuracy is 0.40--0.56 and gene capture 60--68\%.
~~~

Replacement:

~~~latex
among amplicons resolved to a reference named at genus, genus accuracy is 0.58--0.73, and among amplicons resolved to a reference that carries gene families, gene capture is 66--70\%, while over the whole novel-organism population (for genus accuracy, abstentions and references with no named genus scored as unresolved; for gene capture, abstentions and references without gene families scored as zero) genus accuracy is 0.37--0.54 and gene capture 58--64\%.
~~~

Source and values:

- taxacc_region_summary.json (U100) anchor.exclude_0.987.all.all.per_rank.Genus.correct: 0.576329 (V4) to 0.726323 (FullLength16S); .correct_population: 0.370554 (V4) to 0.537691 (FullLength16S).
- gene_capture_by_region.json (U100) _meta.headline_macro_recall_pct.exclude_0.987: 66.383 (V4-V5) to 69.903 (V1-V2); _meta.headline_macro_recall_pct_population.exclude_0.987: 57.643 (V4-V5) to 64.443 (V1-V3).

Note: QUALITATIVE (definition): under the union 1,060 exclude cells (1.6% of 67,038 scoreable) select references that carry no gene families (outcome pred_empty; 0 under legacy). The conditional capture excludes them (macro recall over outcome==ok only, per _meta.headline_note) and the population figure scores them as zero, so the two parenthetical definitions must say so.

### 1.5. L68 (Abstract) — QUALITATIVE

Current:

~~~latex
phylum-standardized conditional gene-capture F1 differs from the stratified panel by at most 0.027 and conditional family accuracy by at most 0.04, while conditional genus accuracy is higher on both alternative panels in every region, by up to 0.07.
~~~

Replacement:

~~~latex
phylum-standardized conditional gene-capture F1 differs from the stratified panel by at most 0.028 and conditional family accuracy by at most 0.05, while conditional genus accuracy is higher on the $k$-means panel in every region, by up to 0.07, and on the random panel in six of eight ($-0.016$ to $+0.062$).
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987]:
- endpoint genecap_cond_f1: max |delta_random_minus_original.est|, |delta_kmeans_minus_original.est| = 0.0278 (random-original, V7-V9); MAG arm (`three_panel/three_panel_mag_panelarm.json`, rule-independent) max 0.0262.
- endpoint family_cond_acc: max |Δ vs original| = 0.0523 (kmeans-original, V1-V3).
- endpoint genus_cond_acc: delta_kmeans_minus_original.est +0.0027 (V3-V4) to +0.0700 (V1-V3), all 8 positive; delta_random_minus_original.est -0.0161 (V7-V9) to +0.0622 (V1-V3), negative in V3-V4 (-0.0100) and V7-V9 (-0.0161).

Note: QUALITATIVE: under legacy the alternative panels were higher in all 8 regions on both panels (random-original +0.0016..+0.0737); under the union the random panel is lower in V3-V4 and V7-V9, so 'higher on both alternative panels in every region' fails. Family rank (no uniform 'family advantage'): random-original / k-means-original family_cond_acc are mixed in sign; excludes_zero only for V1-V3 (+0.0419 / +0.0523, alternative panels higher) and V3-V4 (-0.0251 / -0.0324, original higher); the other six regions cross zero.

### 1.6. L69 (Abstract)

Current:

~~~latex
on the random panel its genus and gene-capture coverage falls by 0.22 to 0.29 in every region
~~~

Replacement:

~~~latex
on the random panel its genus and gene-capture coverage falls by 0.20 to 0.29 in every region
~~~

Source and values:

- three_panel_comparison_d100.json exclude_0.987 genus_coverage delta_random_minus_original.est -0.2037 (V1-V2) to -0.2576 (V4-V5_944R); genecap_coverage -0.2418 (V1-V3) to -0.2899 (V4-V5_944R). Every CI excludes zero (16/16).

### 1.7. L70 (Abstract) — QUALITATIVE

Current:

~~~latex
falling on the random panel by up to 0.05 in short regions and 0.02 on full-length 16S
~~~

Replacement:

~~~latex
falling on the random panel by up to 0.11 in short regions and 0.02 on full-length 16S
~~~

Source and values:

- three_panel_comparison_d100.json cells[regime=include, endpoint=genecap_cond_f1].delta_random_minus_original.est: V4 -0.1069 (largest), V4-V5 -0.1059, V4-V5_944R -0.1048, V7-V9 -0.1014, V3-V4 -0.0910, V1-V2 -0.0617, V1-V3 -0.0504, FullLength16S -0.0159.

Note: QUALITATIVE (magnitude): the self-inclusion panel penalty roughly doubles under the union (legacy max -0.0497), because the union's include-mode selections grow with tie-cluster size (mean n_selected 2.88 vs 1.32).

### 1.8. L72 (Abstract)

Current:

~~~latex
(0.80--0.86 against 0.60--0.74 conditional under self-exclusion), but that robustness does not carry over to metabolism: answers resolved only at family tier recover 61\% of the true gene complement (58--63\% by region) against 72\% for genus-tier answers
~~~

Replacement:

~~~latex
(0.76--0.85 against 0.58--0.73 conditional under self-exclusion), but that robustness does not carry over to metabolism: answers resolved only at family tier recover 58\% of the true gene complement (56--61\% by region) against 69\% for genus-tier answers
~~~

Source and values:

- taxacc_region_summary.json (U100) anchor.exclude_0.987.all.all.per_rank.Family.correct 0.759021 (V4) to 0.847912 (FullLength16S); Genus.correct 0.576329 to 0.726323.
- `region_validation_d100/data/gene_capture_cells.csv`, self_mode exclude_0.987, outcome ok: mean recall tier==family 0.582178 (n 8,451), by region 0.55663 (V4-V5_944R) to 0.60668 (V1-V2); tier==genus 0.692858 (n 51,120). L100 check: 0.610469 (0.58324-0.63039) vs 0.722909 reproduces 61% (58-63%) vs 72%.

Note: The clause 'short amplicons lower family-rank accuracy too, though by less than at genus rank' still holds (family FL-V4 0.089 < genus FL-V4 0.150).

### 1.9. L183 (Methods)

Current:

~~~latex
(12{,}602 of the 67{,}466 self-inclusion amplicons whose source is in the 100-hit record).
~~~

Replacement:

~~~latex
(12{,}602 of the 67{,}466 self-inclusion amplicons whose source is in the 100-hit record).
The ordered reducer remains available as an option (\texttt{select\_references.py -{}-select-mode legacy}), but every result in this paper uses the exact-tie union except the pre-registered primary prFBA arm of the controlled benchmark, whose policy was chosen by leave-one-fold-out selection (the \nameref{sec:sharedref} section).
~~~

Source and values:

- New sentence (the single Methods mention of the alternative). Option name from `/home/freiburger/Documents/prFBA/select_references.py` line 797: `--select-mode {exact-tie,graduated,legacy}`, default exact-tie.

Note: Insertion; the replacement keeps the existing sentence verbatim and appends one line. The numbers in L181/L183 that compare the union with the ordered reducer are retained (see the retained-mentions list).

### 1.10. L316 (Validation preamble) — QUALITATIVE

Current:

~~~latex
46.8\% of the 61{,}509 answered amplicons under self-exclusion on the 10{,}000-organism panel (the legacy reducer's taxonomic scoring in \Cref{tab:main}; its gene-content selection, which feeds the gene-capture columns, answers 60{,}647 and selects one genome for 51.0\%)
~~~

Replacement:

~~~latex
68.1\% of the 61{,}509 answered amplicons under self-exclusion on the 10{,}000-organism panel (the taxonomic scoring of \Cref{tab:main}; its gene-content selection, which feeds the gene-capture columns, answers 61{,}534 and selects one genome for 68.1\% of them)
~~~

Source and values:

- `region_validation_d100/data/taxacc_per_record.jsonl`: self_modes.exclude_0.987, scores.anchor.abstained false for 61,509 of 67,991; n_selected==1 for 41,894 (0.681104).
- `region_validation_d100/data/gene_capture_cells.csv`, self_mode exclude_0.987: n_selected>0 for 61,534 (= 67,991 - 6,457, the prfba_selection abstention count); n_selected==1 for 41,916 (0.681184).
- L100 check: 28,775/61,509 = 0.467818 and 30,941/60,647 = 0.510182 reproduce 46.8% / 51.0%.

Note: QUALITATIVE: under self-exclusion the union is the MORE often single-genome rule (68.1% vs 46.8%/51.0% under legacy); the taxonomic and gene-content selections now agree to 0.01 points.

### 1.11. L342 (10,000-organism validation)

Current:

~~~latex
and score its selection under the legacy reducer (\Cref{tab:main}) against
~~~

Replacement:

~~~latex
and score its selection under the shipped exact-tie union (\Cref{tab:main}) against
~~~

Source and values:

- Rule statement.

### 1.12. L348 (10,000-organism validation) — QUALITATIVE

Current:

~~~latex
still resolves to a reference named at that rank (for gene capture, to any reference), and
~~~

Replacement:

~~~latex
still resolves to a reference named at that rank (for gene capture, to any reference that carries gene families), and
~~~

Source and values:

- Definition: the conditional capture is macro recall over outcome==ok cells (`gene_capture_by_region.json` _meta.headline_note); under the union 1,060 answered cells are pred_empty (_meta.outcome_counts_total.exclude_0.987.pred_empty).

Note: QUALITATIVE (definition), see item for L65.

### 1.13. L348 (10,000-organism validation)

Current:

~~~latex
and, at genus, a further 19\% of genus-scoreable amplicons whose selected reference has no named genus (genus coverage 0.67--0.76 by region); in gene capture, the 10.7\% of scoreable amplicons for which the legacy selection keeps no genome score zero.
~~~

Replacement:

~~~latex
and, at genus, a further 22\% of genus-scoreable amplicons whose selected reference has no named genus (genus coverage 0.64--0.74 by region); in gene capture, the 11.1\% of scoreable amplicons for which the selection keeps no genome (9.6\%) or keeps only references that carry no gene families (1.6\%) score zero.
~~~

Source and values:

- `region_validation_d100/data/taxacc_per_record.jsonl`, exclude_0.987, records with truth_lineage.Genus named (52,026): abstained 5,365 (0.103122, unchanged from legacy; the abstention set is identical under both rules); answered with no named anchor genus 11,230 (0.215854; legacy 9,713 = 0.186695).
- taxacc_region_summary.json (U100) anchor.exclude_0.987.all.all.per_rank.Genus.coverage 0.642956 (V4) to 0.740292 (FullLength16S).
- gene_capture_by_region.json (U100) _meta.outcome_counts_total.exclude_0.987: abstained 6,407, pred_empty 1,060, truth_empty 953; scoreable 67,991-953 = 67,038; 6,407/67,038 = 0.095573, 1,060/67,038 = 0.015812, sum 0.111385. L100: 7,175/67,038 = 0.107029 (reproduces 10.7%).

Note: The '9.5% of all amplicons, 10.3% of genus-scoreable ones' earlier in the sentence is rule-invariant (taxacc 6,482/67,991 = 0.095336) and stays.

### 1.14. L352 (10,000-organism validation) — QUALITATIVE

Current:

~~~latex
Under realistic self-inclusion, anchor genus accuracy is 0.94--0.99 and gene capture (macro recall) is 95--99\% across all eight regions
~~~

Replacement:

~~~latex
Under realistic self-inclusion, anchor genus accuracy is 0.89--1.00 and gene capture (macro recall) is 99--100\% across all eight regions
~~~

Source and values:

- As item L64 (Genus.correct include 0.890933-0.996484; headline_macro_recall_pct.include 99.148-99.979).

Note: QUALITATIVE: include precision drops on short regions (V4 0.7945); see optional item O2.

### 1.15. L353 (10,000-organism validation) — QUALITATIVE

Current:

~~~latex
Under honest self-exclusion the numbers fall, least for full-length 16S and the two V1-anchored regions and most for V4 (\Cref{tab:main}): among resolved amplicons, genus accuracy is 0.60--0.74 and gene capture 69--73\%, while over the whole novel-organism population genus accuracy is 0.40--0.56 and gene capture 60--68\%.
~~~

Replacement:

~~~latex
Under honest self-exclusion the numbers fall, least for full-length 16S and the two V1-anchored regions and most for the five regions that start at V3 or later (\Cref{tab:main}): among resolved amplicons, genus accuracy is 0.58--0.73 and gene capture 66--70\%, while over the whole novel-organism population genus accuracy is 0.37--0.54 and gene capture 58--64\%.
~~~

Source and values:

- Ranges as item L65. Region ordering per column (U100): genus cond lowest V4 0.5763; genus pop lowest V4 0.3706; gene capture cond lowest V4-V5 66.383 (V4 67.303 is 4th lowest); gene capture pop lowest V4-V5 57.643 (V4 58.115). FL and both V1-anchored regions are the top three in every column.

Note: QUALITATIVE: 'most for V4' holds for genus accuracy only; for gene capture the lowest region is V4-V5.

### 1.16. L354 (10,000-organism validation) — QUALITATIVE

Current:

~~~latex
Full-length 16S and the long V1--V3 generalize best; short V4 generalizes worst, and V4--V5 and V7--V9 sit between.
~~~

Replacement:

~~~latex
In genus accuracy full-length 16S and the long V1--V3 generalize best, short V4 worst, and V4--V5 and V7--V9 sit between; gene capture orders the regions differently, highest for V1--V2 and V1--V3 (69.9\% and 69.5\% conditional, against 68.7\% for full-length 16S) and lowest for V4--V5 (66.4\%).
~~~

Source and values:

- gene_capture_by_region.json (U100) _meta.headline_macro_recall_pct.exclude_0.987: V1-V2 69.903, V1-V3 69.518, FullLength16S 68.723, V4-V5 66.383 (lowest), V4 67.303. Genus (taxacc U100 Genus.correct / correct_population): FL 0.7263/0.5377, V1-V3 0.7063/0.5128, V1-V2 0.7026/0.5074, V4 0.5763/0.3706 (lowest), V4-V5 0.6015/0.3909, V7-V9 0.6260/0.4234.

Note: QUALITATIVE: under legacy full-length led gene capture (72.9%) and V4 was lowest (69.0%); under the union V1-V2 leads and V4-V5 is lowest.

### 1.17. L355 (10,000-organism validation) — QUALITATIVE

Current:

~~~latex
The gene-capture leakage gap is near-uniform at $\sim$26 points (conditional) across every region, whereas the genus-accuracy gap widens from 0.24--0.25 for full-length 16S, V1--V3 and V1--V2 to 0.30--0.33 for the five regions that start at V3 or later
~~~

Replacement:

~~~latex
The gene-capture leakage gap varies less across regions, 30--33 points (conditional; 29.6 for V1--V2 to 32.9 for V4--V5), whereas the genus-accuracy gap widens from 0.25--0.27 for full-length 16S, V1--V3 and V1--V2 to 0.31--0.32 for the five regions that start at V3 or later
~~~

Source and values:

- Gene-capture gap = headline_macro_recall_pct.include - .exclude_0.987 (= 100 x regions[r].leakage_gap_macro_recall): FL 31.256, V1-V2 29.607, V1-V3 30.127, V3-V4 32.683, V4 31.845, V4-V5 32.928, V4-V5_944R 32.794, V7-V9 32.097.
- Genus gap = taxacc anchor Genus.correct include - exclude_0.987 (same quantity F11 plots from concordance_long.csv): FL 0.270161, V1-V2 0.247016, V1-V3 0.262514; V3-V4 0.307632, V4 0.314604, V4-V5 0.316421, V4-V5_944R 0.308151, V7-V9 0.312709.
- L100 check: gene-capture 25.883-26.526; genus 0.244702-0.249772 and 0.304635-0.333260.

Note: QUALITATIVE: the gene-capture gap is no longer near-uniform (spread 3.3 points vs 0.6 under legacy) and the genus contrast between the two region groups narrows; the conclusion 'region costs taxonomic resolution more than raw functional recovery' still holds, more weakly.

### 1.18. L356 (10,000-organism validation) — QUALITATIVE

Current:

~~~latex
Under the same legacy reducer, the paired McNemar test finds 24 of the 28 region pairs significantly different at the genus rank under self-exclusion (Holm $p<0.05$; each pair tested on its common-amplified organisms with a named genus, $n=5{,}533$--$7{,}238$, of which $5{,}481$ are common to all eight regions); the four that do not differ are full-length 16S against V1--V3, V1--V2 against V1--V3, V3--V4 against V7--V9 (Holm $p=0.06$--$0.13$ across the four), and the two near-identical V4--V5 primer variants, still the closest pair.
~~~

Replacement:

~~~latex
Under the same exact-tie union, the paired McNemar test finds 26 of the 28 region pairs significantly different at the genus rank under self-exclusion (Holm $p<0.05$; each pair tested on its common-amplified organisms with a named genus, $n=5{,}533$--$7{,}238$, of which $5{,}481$ are common to all eight regions); the two that do not differ are V1--V2 against V1--V3 and V3--V4 against V7--V9 (Holm $p=0.14$ and $0.09$), while the two near-identical V4--V5 primer variants, which disagree on the fewest organisms of any pair (176), differ significantly (Holm $p=6\times10^{-5}$).
~~~

Source and values:

- `region_validation_d100/data/mcnemar_regions.json` by["anchor|Genus|exclude_0.987"]: 26 of 28 pairs significant; n_common 5,533-7,238; n_common_all8 5,481 (both unchanged).
- Non-significant: V1-V2/V1-V3 p_holm 0.141238 (b/c 386/429), V3-V4/V7-V9 p_holm 0.088866 (853/771). V4-V5/V4-V5_944R p_holm 6.3e-05, b/c 58/118 (176 discordant, the fewest; accuracy gap 0.0085, second-closest after V1-V2/V1-V3 0.0075). FullLength16S/V1-V3 p_holm 0.006004 (now significant). Same values in `mcnemar_region_pairs.csv`.

Note: QUALITATIVE: the set of non-differing pairs changes (FL vs V1-V3 and the V4-V5 primer pair become significant) and the primer pair is no longer the closest by accuracy.

### 1.19. L359 (10,000-organism validation)

Current:

~~~latex
conditional family accuracy holds at 0.80--0.86 --- 12--20 points above the genus figure for the same amplicon.
~~~

Replacement:

~~~latex
conditional family accuracy holds at 0.76--0.85 --- 11--18 points above the genus figure for the same amplicon.
~~~

Source and values:

- taxacc_region_summary.json (U100) anchor.exclude_0.987 Family.correct minus Genus.correct: FL 0.12159, V1-V2 0.11444, V1-V3 0.11397, V3-V4 0.16233, V4 0.18269, V4-V5 0.18032, V4-V5_944R 0.18218, V7-V9 0.17071. Family range 0.759021-0.847912.

### 1.20. L360 (10,000-organism validation)

Current:

~~~latex
Answers resolved only at family tier still recover 61\% of the true gene complement under self-exclusion (58--63\% by region), against 72\% for genus-tier answers.
~~~

Replacement:

~~~latex
Answers resolved only at family tier still recover 58\% of the true gene complement under self-exclusion (56--61\% by region), against 69\% for genus-tier answers.
~~~

Source and values:

- As item L72 (gene_capture_cells.csv U100, exclude ok cells: family tier 0.582178, by region 0.55663-0.60668; genus tier 0.692858).

### 1.21. L361 (10,000-organism validation)

Current:

~~~latex
the family$\rightarrow$genus gap is 0.17--0.20, i.e.\ roughly one in five novel organisms
~~~

Replacement:

~~~latex
the family$\rightarrow$genus gap is 0.16--0.18, i.e.\ roughly one in six novel organisms
~~~

Source and values:

- Family-minus-genus conditional accuracy, V3-V4 0.16233, V4 0.18269, V4-V5 0.18032, V4-V5_944R 0.18218, V7-V9 0.17071 (taxacc U100).

Note: Wording: 0.16-0.18 is about one in six rather than one in five.

### 1.22. L366 (10,000-organism validation)

Current:

~~~latex
\caption{Full 10{,}000-organism validation (anchor predictor under the legacy reducer, 100-hit records;
~~~

Replacement:

~~~latex
\caption{Full 10{,}000-organism validation (anchor predictor under the shipped exact-tie union, 100-hit records;
~~~

Source and values:

- Rule statement for tab:main.

### 1.23. L367 (10,000-organism validation) — QUALITATIVE

Current:

~~~latex
capture over amplicons resolved to any reference) and over the whole novel-organism population (for accuracy, abstentions and references unnamed at that rank scored as unresolved; for capture, abstentions scored as zero).}
~~~

Replacement:

~~~latex
capture over amplicons resolved to a reference that carries gene families) and over the whole novel-organism population (for accuracy, abstentions and references unnamed at that rank scored as unresolved; for capture, abstentions and references without gene families scored as zero).}
~~~

Source and values:

- Definition change, as item L65 (1,060 pred_empty exclude cells under the union, 0 under legacy).

Note: QUALITATIVE (definition).

### 1.24. L376 (10,000-organism validation)

Current:

~~~latex
FullLength16S & 0.986 & 0.739 / 0.563 & 99.4\% & 72.9\% / 67.5\% & 0.865 \\
~~~

Replacement:

~~~latex
FullLength16S & 0.996 & 0.726 / 0.538 & 100.0\% & 68.7\% / 63.3\% & 0.848 \\
~~~

Source and values:

- Columns: genus incl. = taxacc anchor.include Genus.correct; excl. cond./pop. = anchor.exclude_0.987 Genus.correct / correct_population; gene capture incl. = gene_capture_by_region.json _meta.headline_macro_recall_pct.include; excl. cond./pop. = headline_macro_recall_pct.exclude_0.987 / headline_macro_recall_pct_population.exclude_0.987; fam. excl. (cond.) = anchor.exclude_0.987 Family.correct. All from `region_validation_d100/data/` (U100).
- Unrounded (incl.; cond / pop; capture incl.; cond / pop; family): 0.996484; 0.726323 / 0.537691; 99.979; 68.723 / 63.275; 0.847912

Note: 99.979 rounds to 100.0\% at one decimal; if a ceiling-looking value is unwanted, write 99.98\% in this one cell.

### 1.25. L377 (10,000-organism validation)

Current:

~~~latex
V1--V2        & 0.967 & 0.722 / 0.537 & 97.5\% & 71.6\% / 66.1\% & 0.838 \\
~~~

Replacement:

~~~latex
V1--V2        & 0.950 & 0.703 / 0.507 & 99.5\% & 69.9\% / 64.3\% & 0.817 \\
~~~

Source and values:

- Columns: genus incl. = taxacc anchor.include Genus.correct; excl. cond./pop. = anchor.exclude_0.987 Genus.correct / correct_population; gene capture incl. = gene_capture_by_region.json _meta.headline_macro_recall_pct.include; excl. cond./pop. = headline_macro_recall_pct.exclude_0.987 / headline_macro_recall_pct_population.exclude_0.987; fam. excl. (cond.) = anchor.exclude_0.987 Family.correct. All from `region_validation_d100/data/` (U100).
- Unrounded (incl.; cond / pop; capture incl.; cond / pop; family): 0.949595; 0.702579 / 0.507398; 99.510; 69.903 / 64.310; 0.817022

### 1.26. L378 (10,000-organism validation)

Current:

~~~latex
V1--V3        & 0.972 & 0.723 / 0.545 & 98.1\% & 72.0\% / 66.9\% & 0.846 \\
~~~

Replacement:

~~~latex
V1--V3        & 0.969 & 0.706 / 0.513 & 99.6\% & 69.5\% / 64.4\% & 0.820 \\
~~~

Source and values:

- Columns: genus incl. = taxacc anchor.include Genus.correct; excl. cond./pop. = anchor.exclude_0.987 Genus.correct / correct_population; gene capture incl. = gene_capture_by_region.json _meta.headline_macro_recall_pct.include; excl. cond./pop. = headline_macro_recall_pct.exclude_0.987 / headline_macro_recall_pct_population.exclude_0.987; fam. excl. (cond.) = anchor.exclude_0.987 Family.correct. All from `region_validation_d100/data/` (U100).
- Unrounded (incl.; cond / pop; capture incl.; cond / pop; family): 0.968789; 0.706275 / 0.512821; 99.645; 69.518 / 64.443; 0.820248

### 1.27. L379 (10,000-organism validation)

Current:

~~~latex
V3--V4        & 0.953 & 0.648 / 0.447 & 96.5\% & 70.0\% / 61.6\% & 0.818 \\
~~~

Replacement:

~~~latex
V3--V4        & 0.922 & 0.614 / 0.406 & 99.4\% & 66.7\% / 58.4\% & 0.777 \\
~~~

Source and values:

- Columns: genus incl. = taxacc anchor.include Genus.correct; excl. cond./pop. = anchor.exclude_0.987 Genus.correct / correct_population; gene capture incl. = gene_capture_by_region.json _meta.headline_macro_recall_pct.include; excl. cond./pop. = headline_macro_recall_pct.exclude_0.987 / headline_macro_recall_pct_population.exclude_0.987; fam. excl. (cond.) = anchor.exclude_0.987 Family.correct. All from `region_validation_d100/data/` (U100).
- Unrounded (incl.; cond / pop; capture incl.; cond / pop; family): 0.921914; 0.614282 / 0.405550 (2,572/6,342 = 0.4055503); 99.362; 66.679 / 58.365; 0.776615

### 1.28. L380 (10,000-organism validation)

Current:

~~~latex
V4            & 0.935 & 0.602 / 0.401 & 95.1\% & 69.0\% / 59.8\% & 0.800 \\
~~~

Replacement:

~~~latex
V4            & 0.891 & 0.576 / 0.371 & 99.1\% & 67.3\% / 58.1\% & 0.759 \\
~~~

Source and values:

- Columns: genus incl. = taxacc anchor.include Genus.correct; excl. cond./pop. = anchor.exclude_0.987 Genus.correct / correct_population; gene capture incl. = gene_capture_by_region.json _meta.headline_macro_recall_pct.include; excl. cond./pop. = headline_macro_recall_pct.exclude_0.987 / headline_macro_recall_pct_population.exclude_0.987; fam. excl. (cond.) = anchor.exclude_0.987 Family.correct. All from `region_validation_d100/data/` (U100).
- Unrounded (incl.; cond / pop; capture incl.; cond / pop; family): 0.890933; 0.576329 / 0.370554; 99.148; 67.303 / 58.115; 0.759021

### 1.29. L381 (10,000-organism validation)

Current:

~~~latex
V4--V5        & 0.950 & 0.628 / 0.431 & 96.0\% & 69.7\% / 60.9\% & 0.825 \\
~~~

Replacement:

~~~latex
V4--V5        & 0.918 & 0.601 / 0.391 & 99.3\% & 66.4\% / 57.6\% & 0.782 \\
~~~

Source and values:

- Columns: genus incl. = taxacc anchor.include Genus.correct; excl. cond./pop. = anchor.exclude_0.987 Genus.correct / correct_population; gene capture incl. = gene_capture_by_region.json _meta.headline_macro_recall_pct.include; excl. cond./pop. = headline_macro_recall_pct.exclude_0.987 / headline_macro_recall_pct_population.exclude_0.987; fam. excl. (cond.) = anchor.exclude_0.987 Family.correct. All from `region_validation_d100/data/` (U100).
- Unrounded (incl.; cond / pop; capture incl.; cond / pop; family): 0.917909; 0.601488 / 0.390884; 99.311; 66.383 / 57.643; 0.781812

### 1.30. L382 (10,000-organism validation)

Current:

~~~latex
V4--V5\_944R  & 0.950 & 0.630 / 0.435 & 96.1\% & 70.1\% / 61.6\% & 0.829 \\
~~~

Replacement:

~~~latex
V4--V5\_944R  & 0.917 & 0.609 / 0.400 & 99.3\% & 66.5\% / 58.1\% & 0.791 \\
~~~

Source and values:

- Columns: genus incl. = taxacc anchor.include Genus.correct; excl. cond./pop. = anchor.exclude_0.987 Genus.correct / correct_population; gene capture incl. = gene_capture_by_region.json _meta.headline_macro_recall_pct.include; excl. cond./pop. = headline_macro_recall_pct.exclude_0.987 / headline_macro_recall_pct_population.exclude_0.987; fam. excl. (cond.) = anchor.exclude_0.987 Family.correct. All from `region_validation_d100/data/` (U100).
- Unrounded (incl.; cond / pop; capture incl.; cond / pop; family): 0.917336; 0.609185 / 0.400342; 99.331; 66.537 / 58.099; 0.791363

### 1.31. L383 (10,000-organism validation)

Current:

~~~latex
V7--V9        & 0.958 & 0.646 / 0.461 & 96.8\% & 70.9\% / 62.7\% & 0.831 \\
~~~

Replacement:

~~~latex
V7--V9        & 0.939 & 0.626 / 0.423 & 99.4\% & 67.3\% / 59.2\% & 0.797 \\
~~~

Source and values:

- Columns: genus incl. = taxacc anchor.include Genus.correct; excl. cond./pop. = anchor.exclude_0.987 Genus.correct / correct_population; gene capture incl. = gene_capture_by_region.json _meta.headline_macro_recall_pct.include; excl. cond./pop. = headline_macro_recall_pct.exclude_0.987 / headline_macro_recall_pct_population.exclude_0.987; fam. excl. (cond.) = anchor.exclude_0.987 Family.correct. All from `region_validation_d100/data/` (U100).
- Unrounded (incl.; cond / pop; capture incl.; cond / pop; family): 0.938718; 0.626009 / 0.423398; 99.396; 67.299 / 59.214; 0.796716

### 1.32. L392 (10,000-organism validation) — QUALITATIVE

Current:

~~~latex
Accuracy is near-ceiling with the source present, falling only at species rank (anchor 0.66--0.87 for every region except full-length 16S, 0.97)
~~~

Replacement:

~~~latex
Accuracy is high with the source present (anchor genus 0.89--1.00), falling most at species rank (anchor 0.65--0.88 for every region except full-length 16S, 0.99)
~~~

Source and values:

- taxacc_region_summary.json (U100) anchor.include.all.all.per_rank.Species.correct: V4 0.6543 (min) ... V1-V3 0.8773 (max of the seven sub-full-length regions); FullLength16S 0.9916. Genus.correct 0.8909-0.9965.
- L100 check: species 0.6574-0.8721, FL 0.9703 (reproduces 0.66-0.87, 0.97).

Note: QUALITATIVE: under the union the short-region anchor loses accuracy at every rank with the source present (V4 anchor phylum 0.9552, family 0.9218, genus 0.8909, against 0.9753/0.9625/0.9352 under legacy), so 'near-ceiling, falling only at species rank' no longer holds for short regions.

### 1.33. L400 (10,000-organism validation) — QUALITATIVE

Current:

~~~latex
Include-self capture is 95--99\%; under honest self-exclusion the columns are population means with abstentions scored as zero capture (60--68\%), and among resolved amplicons capture is 69--73\%
~~~

Replacement:

~~~latex
Include-self capture is 99--100\%; under honest self-exclusion the columns are population means with abstentions and references without gene families scored as zero capture (58--64\%), and among resolved amplicons capture is 66--70\%
~~~

Source and values:

- headline_macro_recall_pct.include 99.148-99.979; headline_macro_recall_pct_population.exclude_0.987 57.643-64.443; headline_macro_recall_pct.exclude_0.987 66.383-69.903 (U100).
- F4b's bars (render_genecap.py fig_F4b) keep abstained and pred_empty cells (both recall 0 in gene_capture_cells.csv; only truth_empty/pgfam_missing/no_amplicon are dropped), so the union bars are exactly these population means (recomputed: FL 63.27 ... V4-V5 57.64).

Note: QUALITATIVE (definition): pred_empty cells now enter the bars as zero.

### 1.34. L408 (10,000-organism validation) — QUALITATIVE

Current:

~~~latex
The genus gap is 0.24--0.25 for full-length 16S, V1--V3 and V1--V2 and 0.30--0.33 for the other five regions (V3--V4, V4, both V4--V5 variants and V7--V9); the gene-capture gap, not shown, is near-uniform at $\sim$26 points
~~~

Replacement:

~~~latex
The genus gap is 0.25--0.27 for full-length 16S, V1--V3 and V1--V2 and 0.31--0.32 for the other five regions (V3--V4, V4, both V4--V5 variants and V7--V9); the gene-capture gap, not shown, is 30--33 points
~~~

Source and values:

- As item L355 (genus gap 0.247016-0.270161 and 0.307632-0.316421; gene-capture gap 29.607-32.928).

Note: QUALITATIVE: no longer near-uniform.

### 1.35. L416 (10,000-organism validation)

Current:

~~~latex
Family-level assignment remains robust (0.80--0.86) where
~~~

Replacement:

~~~latex
Family-level assignment remains robust (0.76--0.85) where
~~~

Source and values:

- taxacc U100 anchor.exclude_0.987 Family.correct 0.759021-0.847912 (the quantity F2's exclude-self panel plots).

### 1.36. L443 (PICRUSt2 as shipped)

Current:

~~~latex
(F1 0.713 versus 0.729; under the legacy reducer 0.726, a gap within resampling noise)
~~~

Replacement:

~~~latex
(F1 0.713 versus 0.729)
~~~

Source and values:

- Deletes the legacy-key comparison (`region_validation_d100/picrust2_comparison/data/comparison_summary.json` by_region.ec.all.all.prfba_exclude_legacy.population.f1 = 0.7265; cells 0.72648). The union values stay: prfba_exclude 0.713, picrust2 0.7287.

Note: Optional replacement for the deleted significance statement: the genome-clustered bootstrap of D100_VS_D20.md (scratchpad p2cells.json, not a repo output) gives the union gap -0.016 [-0.021, -0.011].

### 1.37. L631 (Retention depth)

Current:

~~~latex
recall 0.847 against 0.835), and under the legacy reducer conditional genus accuracy moves by less than 0.01 in every region, while conditional gene capture rises slightly rather than falling, by 0.7 to 1.3 points.
~~~

Replacement:

~~~latex
recall 0.847 against 0.835).
~~~

Source and values:

- Deletion. No like-for-like depth-20 union scoring of tab:main exists (region_validation/data/taxacc_region_summary.json is select_mode legacy; the June `*_et.json` exact-tie variants predate the corrected scorer: no select_mode, no correct_population, no unscoreable-rank exclusion).

Note: UNRESOLVED for a union replacement; see section 4.

### 1.38. L632 (Retention depth)

Current:

~~~latex
; population genus accuracy under the legacy reducer, the rule of \Cref{tab:main}, rises from 0.375 to 0.472; and on the as-shipped PICRUSt2 axis prFBA's population EC F1 rises from 0.600 to 0.713, still below PICRUSt2's 0.729 (prFBA 0.726 under the legacy reducer).
~~~

Replacement:

~~~latex
; and on the as-shipped PICRUSt2 axis prFBA's population EC F1 rises from 0.600 to 0.713, still below PICRUSt2's 0.729.
~~~

Source and values:

- Deletion of two legacy clauses. The union depth-100 pooled population genus accuracy is 0.437916 (taxacc U100: sum n_correct / sum n over Genus.exclude_0.987, 22,783/52,026 = 0.4379), but there is no depth-20 union base to pair it with (section 4).

Note: UNRESOLVED for a d20->d100 union genus statement.

### 1.39. L633 (Retention depth)

Current:

~~~latex
gives coverage $\sim$0.96, population genus accuracy $\sim$0.49 (the legacy-reducer depth-20 base of 0.375, with the recovered amplicons credited at the shipped selector's probe accuracy) and as-shipped population EC F1 $\sim$0.75.
~~~

Replacement:

~~~latex
gives coverage $\sim$0.96 and as-shipped population EC F1 $\sim$0.75.
~~~

Source and values:

- Deletion of the legacy-based genus projection (its base 0.375 is the legacy depth-20 tab:main figure). The coverage (0.766 -> ~0.96) and EC F1 (0.600 -> ~0.75) projections are union-based and stay.

Note: UNRESOLVED for a union genus projection.

### 1.40. L641 (Retention depth) — QUALITATIVE

Current:

~~~latex
with both panels scored under the legacy reducer at depth 100, the full-length 16S emptied-record rate rises from 4.7\% to 22.9\% (8.3\% to 29.4\% pooled over all eight regions) while the share of full-length records that survive exclusion yet go uncovered stays below 3\% on both (2.8\% and 2.3\%; 2.6\% and 2.1\% pooled over all eight regions)
~~~

Replacement:

~~~latex
with both panels scored under the exact-tie union at depth 100, the full-length 16S emptied-record rate rises from 4.7\% to 22.9\% (8.3\% to 29.4\% pooled over all eight regions) while the share of full-length records that survive exclusion yet go uncovered stays near 3\% on both (3.4\% and 2.8\%; 3.1\% and 2.4\% pooled over all eight regions)
~~~

Source and values:

- `three_panel/emptied_strata_d100.json` panels.original.survivors non_ok/n = 236/6,879 = 0.034307; panels.kmeans.survivors 188/6,796 = 0.027663; emptied_pct 4.7 / 22.9 (rule-invariant: identical in _d100_legacy).
- Pooled over 8 regions from each tree's gene_capture_cells.csv (exclude_0.987, outcome != truth_empty, emptied = n_removed_self >= 100, surviving non-ok = outcome != ok): original 0.031033, k-means 0.024318; emptied 0.08292 / 0.29368.

Note: QUALITATIVE (wording): 'stays below 3%' fails for the original panel (3.4%); 'near 3%' holds.

### 1.41. L723 (MAG versus 16S)

Current:

~~~latex
(record depth 100, scored under the legacy reducer as in \Cref{tab:main}), the MAG arm leads full-length 16S by 0.122 and V4 by 0.192 in gene-capture F1
~~~

Replacement:

~~~latex
(record depth 100, scored under the exact-tie union as in \Cref{tab:main}), the MAG arm leads full-length 16S by 0.103 and V4 by 0.212 in gene-capture F1
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`.
- 100%: FL Δ +0.102997 (n 1,916), V4 Δ +0.212059 (n 2,299). 30%: FL +0.100978, V4 +0.212686, so 'the same margin at 30% as at 100%' still holds.

### 1.42. L725 (MAG versus 16S)

Current:

~~~latex
the MAG arm still leads full-length 16S by 0.096 F1 and V4 by 0.169, winning 1{,}240 organisms to 373 ($p\approx 3\times10^{-121}$; \Cref{fig:magcross}).
~~~

Replacement:

~~~latex
the MAG arm still leads full-length 16S by 0.077 F1 and V4 by 0.190, winning 1{,}129 organisms to 473 ($p\approx 1\times10^{-79}$; \Cref{fig:magcross}).
~~~

Source and values:

- `mag_validation/out/crossover_d100.json` row completeness 0.01: FullLength16S_delta 0.0767 (full precision 0.076708), V4_delta 0.1903 (0.190273), FullLength16S_mag_wins 1129, FullLength16S_amp_wins 473, FullLength16S_p 1.3798e-79 (Wilcoxon, zsplit).

Note: The FL win ratio weakens (1,129:473 vs 1,240:373) and p rises by 42 orders of magnitude, but the lead is positive everywhere.

### 1.43. L732 (MAG versus 16S)

Current:

~~~latex
it is $+0.105$ to $+0.122$ F1 against full-length 16S and the two V1-anchored regions (V1--V2, V1--V3) and $+0.167$ to $+0.192$ against the five regions that start at V3 or later, with V4 the widest gap.
~~~

Replacement:

~~~latex
it is $+0.103$ to $+0.116$ F1 against full-length 16S and the two V1-anchored regions (V1--V2, V1--V3) and $+0.168$ to $+0.212$ against the five regions that start at V3 or later, with V4 the widest gap.
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`.
- 100% Δ: FL +0.102997, V1-V2 +0.116243, V1-V3 +0.109723; V3-V4 +0.186035, V4 +0.212059, V4-V5 +0.192489, V4-V5_944R +0.188331, V7-V9 +0.167893.

Note: The smallest gap is now full-length 16S (under legacy V1-V2); V4 is still the widest.

### 1.44. L733 (MAG versus 16S)

Current:

~~~latex
(MAG 0.816--0.823 against 0.728--0.784)
~~~

Replacement:

~~~latex
(MAG 0.816--0.824 against 0.703--0.768)
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Metric recall.
- Conditional recall at 100%: MAG 0.816420 (V7-V9) to 0.823717 (V3-V4); 16S 0.703338 (V4-V5) to 0.767589 (V1-V2).

### 1.45. L734 (MAG versus 16S)

Current:

~~~latex
significant at $p < 10^{-179}$.
~~~

Replacement:

~~~latex
significant at $p < 10^{-161}$.
~~~

Source and values:

- `paired_comparison_d100.json` exclude|c1_x0|<region>|f1 conditional.p: largest V1-V3 6.2835e-162 (FL 4.0017e-165, V1-V2 1.5371e-179, V3-V4 2.8424e-260, V7-V9 4.1902e-303; V4, V4-V5, V4-V5_944R stored null = underflow below ~1e-300).

### 1.46. L737 (MAG versus 16S)

Current:

~~~latex
where an amplicon answers for 87--93\%, scoring an abstention as zero recovered metabolism puts the amplicon ahead everywhere --- 0.594 against 0.248 for full-length 16S, and 0.497 against 0.253 for V4.
~~~

Replacement:

~~~latex
where an amplicon answers for 86--93\%, scoring an abstention as zero recovered metabolism puts the amplicon ahead everywhere --- 0.605 against 0.248 for full-length 16S, and 0.483 against 0.253 for V4.
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`.
- Population 16S coverage (share of paired organisms with amplicon ok): 0.863463 (V4) to 0.926982 (V1-V3). Population F1: FL MAG 0.247726 / 16S 0.604584; V4 0.252973 / 0.482968.

Note: 'The gap narrows for the five regions from V3 on' (next sentence) still holds: population gap FL 0.357, V1-V2 0.351, V1-V3 0.359 vs V3+ 0.230-0.274.

### 1.47. L747 (MAG versus 16S)

Current:

~~~latex
Amplicon cells come from the depth-100 records, scored under the legacy reducer as in \Cref{tab:main}.}
~~~

Replacement:

~~~latex
Amplicon cells come from the depth-100 records, scored under the exact-tie union as in \Cref{tab:main}.}
~~~

Source and values:

- Rule statement (tab:magvs16s caption).

Note: The caption's 'ordering is unchanged at every completeness down to 1%' still holds: under the union the region order of Δ (V4 > V4-V5 > V4-V5_944R > V3-V4 > V7-V9 > V1-V2 > V1-V3 > FL) is identical at all 14 levels.

### 1.48. L757 (MAG versus 16S)

Current:

~~~latex
Full-length & 1{,}929 & 0.813 & 0.691 & $+0.122$ ($+0.087$) & 0.304 & 0.926 & 0.248 & \textbf{0.594} \\
~~~

Replacement:

~~~latex
Full-length & 1{,}916 & 0.813 & 0.710 & $+0.103$ ($+0.056$) & 0.304 & 0.921 & 0.248 & \textbf{0.605} \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Completeness 100%.
- Unrounded: n 1916; MAG 0.813011; 16S 0.710014; Δ 0.102997; median 0.055620; cov MAG 0.304173 / 16S 0.920699; pop MAG 0.247726 / 16S 0.604584

### 1.49. L758 (MAG versus 16S)

Current:

~~~latex
V1--V2        & 1{,}999 & 0.814 & 0.709 & $+0.105$ ($+0.062$) & 0.298 & 0.922 & 0.242 & \textbf{0.598} \\
~~~

Replacement:

~~~latex
V1--V2        & 1{,}992 & 0.814 & 0.698 & $+0.116$ ($+0.059$) & 0.298 & 0.920 & 0.242 & \textbf{0.593} \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Completeness 100%.
- Unrounded: n 1992; 0.814043; 0.697800; 0.116243; 0.059069; 0.297704 / 0.919968; 0.242474 / 0.592708

### 1.50. L759 (MAG versus 16S)

Current:

~~~latex
V1--V3        & 2{,}026 & 0.813 & 0.703 & $+0.110$ ($+0.063$) & 0.298 & 0.930 & 0.243 & \textbf{0.601} \\
~~~

Replacement:

~~~latex
V1--V3        & 2{,}020 & 0.813 & 0.704 & $+0.110$ ($+0.052$) & 0.298 & 0.927 & 0.243 & \textbf{0.602} \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Completeness 100%.
- Unrounded: n 2020; 0.813247; 0.703525; 0.109723; 0.052100; 0.297935 / 0.926982; 0.242559 / 0.601768

### 1.51. L760 (MAG versus 16S)

Current:

~~~latex
V3--V4        & 1{,}902 & 0.815 & 0.638 & $+0.178$ ($+0.145$) & 0.298 & 0.880 & 0.242 & \textbf{0.519} \\
~~~

Replacement:

~~~latex
V3--V4        & 1{,}885 & 0.815 & 0.629 & $+0.186$ ($+0.148$) & 0.298 & 0.875 & 0.242 & \textbf{0.516} \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Completeness 100%.
- Unrounded: n 1885; 0.815472; 0.629436; 0.186035; 0.147995; 0.297661 / 0.875288; 0.242070 / 0.516159

### 1.52. L761 (MAG versus 16S)

Current:

~~~latex
V4            & 2{,}309 & 0.809 & 0.617 & $+0.192$ ($+0.163$) & 0.313 & 0.867 & 0.253 & \textbf{0.497} \\
~~~

Replacement:

~~~latex
V4            & 2{,}299 & 0.810 & 0.598 & $+0.212$ ($+0.171$) & 0.313 & 0.863 & 0.253 & \textbf{0.483} \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Completeness 100%.
- Unrounded: n 2299; 0.810101; 0.598042; 0.212059; 0.170766; 0.313294 / 0.863463; 0.252973 / 0.482968

### 1.53. L762 (MAG versus 16S)

Current:

~~~latex
V4--V5        & 2{,}342 & 0.809 & 0.619 & $+0.189$ ($+0.164$) & 0.314 & 0.874 & 0.254 & \textbf{0.507} \\
~~~

Replacement:

~~~latex
V4--V5        & 2{,}334 & 0.809 & 0.616 & $+0.192$ ($+0.157$) & 0.314 & 0.868 & 0.254 & \textbf{0.504} \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Completeness 100%.
- Unrounded: n 2334; 0.808570; 0.616081; 0.192489; 0.157451; 0.314274 / 0.868318; 0.253839 / 0.503519

### 1.54. L763 (MAG versus 16S)

Current:

~~~latex
V4--V5 944R & 2{,}269 & 0.811 & 0.622 & $+0.189$ ($+0.164$) & 0.318 & 0.879 & 0.257 & \textbf{0.512} \\
~~~

Replacement:

~~~latex
V4--V5 944R & 2{,}255 & 0.811 & 0.623 & $+0.188$ ($+0.157$) & 0.318 & 0.873 & 0.257 & \textbf{0.510} \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Completeness 100%.
- Unrounded: n 2255; 0.811428; 0.623097; 0.188331; 0.156866; 0.317894 / 0.873153; 0.257432 / 0.509743

### 1.55. L764 (MAG versus 16S)

Current:

~~~latex
V7--V9        & 2{,}321 & 0.808 & 0.641 & $+0.167$ ($+0.141$) & 0.326 & 0.885 & 0.264 & \textbf{0.522} \\
~~~

Replacement:

~~~latex
V7--V9        & 2{,}301 & 0.808 & 0.641 & $+0.168$ ($+0.126$) & 0.326 & 0.880 & 0.264 & \textbf{0.523} \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Completeness 100%.
- Unrounded: n 2301; 0.808428; 0.640535; 0.167893; 0.126156; 0.326092 / 0.879831; 0.263601 / 0.522857

### 1.56. L772 (MAG versus 16S)

Current:

~~~latex
the amplicon answers 98--99\% of the time
~~~

Replacement:

~~~latex
the amplicon answers 98\% of the time
~~~

Source and values:

- U100 gene_capture_cells.csv exclude_0.987, outcomes ok/abstained/pred_empty (pred_empty counted as no gene-content answer), binned by best_nonself_identity, restricted to organisms in the MAG population: FL 0.97-0.987 0.9774 (n 1,504), 0.945-0.97 0.9819 (n 993); all regions 0.9825 (n 10,127) and 0.9810 (n 6,537). Same populations under L100 give 0.9845-0.9889 (reproduces 98-99%). MAG side (8-20%) is rule-independent.

Note: If pred_empty cells are instead dropped from the denominator the rate is 0.995-0.999; the replacement keeps the L100 population definition.

### 1.57. L778 (MAG versus 16S)

Current:

~~~latex
and only 31.8\% of them exceed the full-length 16S level of 0.691
~~~

Replacement:

~~~latex
and only 26.7\% of them exceed the full-length 16S level of 0.710
~~~

Source and values:

- Comparator = union FL 16S conditional F1 at 100% in tab:magvs16s (0.710014). Of the 1,367 newly admitted 2%-AF-gate MAGs with a non-empty prediction (mean F1 0.59946), 365 have F1 > 0.71001: 365/1,367 = 0.267008.
- MAG-side cells: scratchpad pickle `verify/agent_scratch2/sk3/afgate2.pkl` (15%- and 2%-gate MAG cells; not a repo file). L100 check: 435/1,367 = 0.318215 at 0.69117 reproduces 31.8%.

Note: Provenance caveat: the 2%-gate MAG cells exist only in the session scratchpad.

### 1.58. L780 (MAG versus 16S)

Current:

~~~latex
by 0.08--0.19 F1, at any completeness at or above 1\%.
~~~

Replacement:

~~~latex
by 0.08--0.21 F1, at any completeness at or above 1\%.
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`.
- Over all 14 levels x 8 regions: min +0.076708 (FL, 1%), max +0.214472 (V4, 90%). L100: +0.0815 to +0.1938.

### 1.59. L782 (MAG versus 16S)

Current:

~~~latex
the amplicon arm still answers at F1 0.62 (full-length) and 0.56 (V4), so scored over the whole population the amplicon arm wins, 0.594 against 0.248.
~~~

Replacement:

~~~latex
the amplicon arm still answers at F1 0.64 (full-length) and 0.54 (V4), so scored over the whole population the amplicon arm wins, 0.605 against 0.248.
~~~

Source and values:

- U100 amplicon ok-cell mean F1 on organisms the MAG arm (100%, exclude, x0) does not answer: FL 0.635016 (n 4,725), V4 0.544154 (n 5,859). Population FL 0.604584 vs 0.247726. L100 check 0.62126 / 0.55548.

### 1.60. L783 (MAG versus 16S)

Current:

~~~latex
are also 0.050 easier for the amplicon arm
~~~

Replacement:

~~~latex
are also 0.053 easier for the amplicon arm
~~~

Source and values:

- U100 FL amplicon ok-cell mean F1 on MAG-answered organisms 0.710014 (n 1,916) minus all FL ok cells 0.656627 (n 6,643) = +0.053387. L100 check +0.04975.

### 1.61. L790 (MAG versus 16S)

Current:

~~~latex
widens from $+0.122$ to $+0.145$ and population F1 rises from 0.248 to 0.301 against the amplicon's unchanged 0.594
~~~

Replacement:

~~~latex
widens from $+0.103$ to $+0.127$ and population F1 rises from 0.248 to 0.301 against the amplicon's unchanged 0.605
~~~

Source and values:

- `mag_validation/out/mag_cells_tau99.csv` (exclude, completeness 1.0, contamination 0) paired with U100 FL cells: conditional n 2,296, MAG 0.829406, 16S 0.701969, Δ +0.127437; population MAG 0.300818 vs 16S 0.604584. L100 check: Δ +0.144954, 16S pop 0.594221 (reproduces +0.145 / 0.594).

### 1.62. L797 (MAG versus 16S)

Current:

~~~latex
(depth-100 records, legacy reducer); $\Delta$ is the difference
~~~

Replacement:

~~~latex
(depth-100 records, exact-tie union); $\Delta$ is the difference
~~~

Source and values:

- Rule statement (tab:magval caption).

### 1.63. L808 (MAG versus 16S)

Current:

~~~latex
100\% & 0.308 & 0.807 & 0.813 / 0.691 ($+0.122$) & 0.809 / 0.617 ($+0.192$) & 1{,}929 \\
~~~

Replacement:

~~~latex
100\% & 0.308 & 0.807 & 0.813 / 0.710 ($+0.103$) & 0.810 / 0.598 ($+0.212$) & 1{,}916 \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Resolve rate and MAG cond. F1 columns are MAG-only and unchanged.
- Unrounded: FL 0.813011/0.710014 (+0.102997), V4 0.810101/0.598042 (+0.212059), n(FL) 1916 (crossover_d100.json at 4 dp agrees).

### 1.64. L809 (MAG versus 16S)

Current:

~~~latex
50\%  & 0.307 & 0.807 & 0.812 / 0.692 ($+0.121$) & 0.809 / 0.616 ($+0.193$) & 1{,}940 \\
~~~

Replacement:

~~~latex
50\%  & 0.307 & 0.807 & 0.812 / 0.710 ($+0.102$) & 0.809 / 0.596 ($+0.213$) & 1{,}927 \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Resolve rate and MAG cond. F1 columns are MAG-only and unchanged.
- Unrounded: FL 0.81210/0.71001 (+0.10209), V4 0.80946/0.59623 (+0.21323), n 1927 (crossover_d100.json at 4 dp agrees).

### 1.65. L810 (MAG versus 16S)

Current:

~~~latex
25\%  & 0.309 & 0.806 & 0.811 / 0.690 ($+0.122$) & 0.808 / 0.615 ($+0.192$) & 1{,}946 \\
~~~

Replacement:

~~~latex
25\%  & 0.309 & 0.806 & 0.811 / 0.708 ($+0.103$) & 0.808 / 0.595 ($+0.213$) & 1{,}932 \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Resolve rate and MAG cond. F1 columns are MAG-only and unchanged.
- Unrounded: FL 0.81110/0.70779 (+0.10331), V4 0.80847/0.59516 (+0.21330), n 1932 (crossover_d100.json at 4 dp agrees).

### 1.66. L811 (MAG versus 16S)

Current:

~~~latex
10\%  & 0.295 & 0.803 & 0.807 / 0.689 ($+0.117$) & 0.804 / 0.614 ($+0.190$) & 1{,}857 \\
~~~

Replacement:

~~~latex
10\%  & 0.295 & 0.803 & 0.807 / 0.708 ($+0.099$) & 0.805 / 0.594 ($+0.211$) & 1{,}845 \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Resolve rate and MAG cond. F1 columns are MAG-only and unchanged.
- Unrounded: FL 0.80664/0.70759 (+0.09905), V4 0.80512/0.59415 (+0.21096), n 1845 (crossover_d100.json at 4 dp agrees).

### 1.67. L812 (MAG versus 16S)

Current:

~~~latex
5\%   & 0.292 & 0.797 & 0.798 / 0.689 ($+0.109$) & 0.798 / 0.616 ($+0.181$) & 1{,}834 \\
~~~

Replacement:

~~~latex
5\%   & 0.292 & 0.797 & 0.798 / 0.707 ($+0.091$) & 0.798 / 0.596 ($+0.202$) & 1{,}821 \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Resolve rate and MAG cond. F1 columns are MAG-only and unchanged.
- Unrounded: FL 0.79787/0.70662 (+0.09125), V4 0.79799/0.59614 (+0.20185), n 1821 (crossover_d100.json at 4 dp agrees).

### 1.68. L813 (MAG versus 16S)

Current:

~~~latex
2\%   & 0.283 & 0.791 & 0.794 / 0.689 ($+0.105$) & 0.792 / 0.615 ($+0.177$) & 1{,}804 \\
~~~

Replacement:

~~~latex
2\%   & 0.283 & 0.791 & 0.794 / 0.707 ($+0.087$) & 0.792 / 0.593 ($+0.199$) & 1{,}791 \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Resolve rate and MAG cond. F1 columns are MAG-only and unchanged.
- Unrounded: FL 0.79352/0.70694 (+0.08658), V4 0.79234/0.59308 (+0.19926), n 1791 (crossover_d100.json at 4 dp agrees).

### 1.69. L814 (MAG versus 16S)

Current:

~~~latex
1\%   & 0.253 & 0.783 & 0.783 / 0.687 ($+0.096$) & 0.783 / 0.614 ($+0.169$) & 1{,}614 \\
~~~

Replacement:

~~~latex
1\%   & 0.253 & 0.783 & 0.783 / 0.706 ($+0.077$) & 0.783 / 0.593 ($+0.190$) & 1{,}603 \\
~~~

Source and values:

- Paired within organism with compare_arms.py logic: MAG cells `mag_validation/out/mag_cells_all.csv` (+ `lowc_out_batches/b*/mag_cells.csv` for <=25%), regime exclude, contamination 0, MAG ok = outcome ok with non-empty precision; amplicon cells `region_validation_d100/data/gene_capture_cells.csv` exclude_0.987, ok = outcome ok; conditional = both ok; population = both scoreable, non-ok scored 0. Cross-checked against `mag_validation/out/paired_comparison_d100.json` (keys `exclude|c1_x0|<region>|f1`, 4-dp) and `crossover_d100.json`. Resolve rate and MAG cond. F1 columns are MAG-only and unchanged.
- Unrounded: FL 0.78294/0.70624 (+0.07671), V4 0.78312/0.59285 (+0.19027), n 1603 (crossover_d100.json at 4 dp agrees).

### 1.70. L841 (Panel composition)

Current:

~~~latex
under one selection rule, the legacy reducer)
~~~

Replacement:

~~~latex
under one selection rule, the shipped exact-tie union)
~~~

Source and values:

- Rule statement; union trees region_validation_d100, region_validation_kmeans_d100, region_validation_random_d100 (`three_panel_comparison_d100.json` meta.roots).

### 1.71. L917 (Panel composition)

Current:

~~~latex
by at most 0.031 (V7--V9), and within 0.016 at every MAG completeness
~~~

Replacement:

~~~latex
by at most 0.029 (V7--V9), and within 0.016 at every MAG completeness
~~~

Source and values:

- `three_panel_comparison_d100.json` exclude_0.987 genecap_cond_f1 crude: random-original V7-V9 -0.029470 (largest), k-means-original V7-V9 -0.022918 (largest; 'within 0.023' unchanged); original ahead of both panels in all 8 regions.

### 1.72. L918 (Panel composition)

Current:

~~~latex
(Welch $p\leq 0.005$, and $p<10^{-7}$ for V4, V4--V5, V4--V5\_944R and V7--V9), full-length 16S included ($+0.014$, $p=1.4\times10^{-6}$)
~~~

Replacement:

~~~latex
(Welch $p<0.01$, and $p<10^{-8}$ for V4, V4--V5 and V7--V9), full-length 16S included ($+0.011$, $p=2.3\times10^{-4}$)
~~~

Source and values:

- Welch t-test (scipy ttest_ind, equal_var=False) of exclude_0.987 ok-cell F1, `region_validation_d100` vs `region_validation_kmeans_d100` gene_capture_cells.csv: FL +0.011234 p 2.3469e-4; V1-V2 p 5.08e-4; V1-V3 p 8.27e-3 (largest); V3-V4 2.76e-4; V4 2.91e-9; V4-V5 1.91e-9; V4-V5_944R 5.87e-6; V7-V9 5.91e-13. L100 check reproduces p<=0.005, the four p<1e-7 and FL +0.014, p 1.41e-6.

Note: V4-V5_944R drops out of the p<10^-8 group; all eight remain significant (uncorrected; 7/8 under Bonferroni).

### 1.73. L919 (Panel composition)

Current:

~~~latex
(0.753 and 0.654 against 0.926 for full-length 16S)
~~~

Replacement:

~~~latex
(0.750 and 0.650 against 0.921 for full-length 16S)
~~~

Source and values:

- `three_panel_comparison_d100.json` FullLength16S exclude_0.987 genecap_coverage crude: original 0.920721, kmeans 0.749631, random 0.650422.

### 1.74. L926 (Panel composition)

Current:

~~~latex
Crude values; self-exclusion throughout; 16S arm under the legacy reducer;
~~~

Replacement:

~~~latex
Crude values; self-exclusion throughout; 16S arm under the exact-tie union;
~~~

Source and values:

- Rule statement (tab:twopanel caption).

### 1.75. L938 (Panel composition)

Current:

~~~latex
full-length & 0.926 & 0.753 & 0.654 & 0.641 & 0.627 & 0.620\\
~~~

Replacement:

~~~latex
full-length & 0.921 & 0.750 & 0.650 & 0.657 & 0.645 & 0.641\\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987, endpoint=genecap_coverage / genecap_cond_f1].{original,kmeans,random}.crude (the paper's crude column reproduces each panel's scorer output).
- Unrounded original/k-means/random: cov 0.920721/0.749631/0.650422; cond F1 0.656630/0.645396/0.641129

### 1.76. L939 (Panel composition)

Current:

~~~latex
V1--V3      & 0.930 & 0.759 & 0.667 & 0.646 & 0.638 & 0.631\\
~~~

Replacement:

~~~latex
V1--V3      & 0.927 & 0.757 & 0.664 & 0.649 & 0.641 & 0.636\\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987, endpoint=genecap_coverage / genecap_cond_f1].{original,kmeans,random}.crude (the paper's crude column reproduces each panel's scorer output).
- Unrounded original/k-means/random: cov 0.927001/0.756902/0.664318; cond F1 0.649184/0.640975/0.636156

### 1.77. L940 (Panel composition)

Current:

~~~latex
V4          & 0.867 & 0.648 & 0.537 & 0.573 & 0.554 & 0.549\\
~~~

Replacement:

~~~latex
V4          & 0.863 & 0.647 & 0.536 & 0.559 & 0.541 & 0.536\\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987, endpoint=genecap_coverage / genecap_cond_f1].{original,kmeans,random}.crude (the paper's crude column reproduces each panel's scorer output).
- Unrounded original/k-means/random: cov 0.863492/0.647484/0.535591; cond F1 0.559324/0.540624/0.535636

### 1.78. L941 (Panel composition)

Current:

~~~latex
V7--V9      & 0.885 & 0.661 & 0.550 & 0.590 & 0.567 & 0.559\\
~~~

Replacement:

~~~latex
V7--V9      & 0.880 & 0.659 & 0.547 & 0.594 & 0.571 & 0.565\\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987, endpoint=genecap_coverage / genecap_cond_f1].{original,kmeans,random}.crude (the paper's crude column reproduces each panel's scorer output).
- Unrounded original/k-means/random: cov 0.879858/0.658700/0.547493; cond F1 0.594277/0.571359/0.564807

### 1.79. L954 (Panel composition)

Current:

~~~latex
gap of 0.161 against 0.173 crude
~~~

Replacement:

~~~latex
gap of 0.159 against 0.171 crude
~~~

Source and values:

- `three_panel_comparison_d100.json` FullLength16S exclude_0.987 genecap_coverage: delta_kmeans_minus_original.est -0.1588; crude 0.749631 - 0.920721 = -0.171090.

### 1.80. L957 (Panel composition)

Current:

~~~latex
(2.8\%, 2.3\% and 1.9\% of surviving full-length records for the original, $k$-means and random panels, and 0.4--3.1\% in every full-length stratum; 2.6\%, 2.1\% and 1.4\% pooled over all eight regions), some 35- to 53-fold below
~~~

Replacement:

~~~latex
(3.4\%, 2.8\% and 2.3\% of surviving full-length records for the original, $k$-means and random panels, and 1.5--3.7\% in every full-length stratum; 3.1\%, 2.4\% and 1.8\% pooled over all eight regions), some 29- to 43-fold below
~~~

Source and values:

- `three_panel/emptied_strata_d100.json` panels.<p>.survivors: original 236/6,879 = 0.034307, kmeans 188/6,796 = 0.027663, random 137/5,836 = 0.023475.
- Per-stratum (0-19 / 20-49 / 50-99) uncovered from the cells: original 0.03664 / 0.01731 / 0.01471; k-means 0.02835 / 0.02219 / 0.03211; random 0.02264 / 0.02080 / 0.02837 -> range 1.47-3.66%.
- Pooled 8 regions (gene_capture_cells.csv, as item L641): 0.031033 / 0.024318 / 0.018350.
- Fold below the emptied-record uncovered rate of 1.000: 1/0.034307 = 29.1, 1/0.023475 = 42.6.

### 1.81. L958 (Panel composition) — QUALITATIVE

Current:

~~~latex
differ by under a point (1.2 points pooled over all eight regions)
~~~

Replacement:

~~~latex
differ by about a point (1.1 points, and 1.3 pooled over all eight regions)
~~~

Source and values:

- FL survivors 3.4307% - 2.3475% = 1.083 points; pooled 3.1033% - 1.8350% = 1.268 points.

Note: QUALITATIVE (wording): 'under a point' fails (1.08). The direction claim (original panel leaves the most surviving records uncovered) still holds.

### 1.82. L958 (Panel composition)

Current:

~~~latex
(Cochran--Mantel--Haenszel $p=0.018$; original versus $k$-means $p=0.28$)
~~~

Replacement:

~~~latex
(Cochran--Mantel--Haenszel $p=0.004$; original versus $k$-means $p=0.15$)
~~~

Source and values:

- `emptied_strata_d100.json` tests.cmh_survivor_strata.original_vs_random.p 0.004066; .original_vs_kmeans.p 0.148370 (L100 file: 0.01772 / 0.28443).

Note: 'Only the original-versus-random difference is significant' still holds (stratified test).

### 1.83. L960 (Panel composition)

Current:

~~~latex
to within a point and a half, the difference
~~~

Replacement:

~~~latex
to within two points, the difference
~~~

Source and values:

- FL coverage differences vs emptied-rate differences (union): original-random 0.92072-0.65042 = 0.27030 vs 0.33394-0.04657 = 0.28737 (1.7 points); original-k-means 0.17109 vs 0.18247 (1.1); k-means-random 0.09921 vs 0.10490 (0.6). L100 max was 1.4 points.

### 1.84. L960 (Panel composition)

Current:

~~~latex
against an observed 0.654, $1-0.229=0.771$ against 0.753, and $1-0.047=0.953$ against 0.926, the residual being the few percent of surviving records that still abstain.
~~~

Replacement:

~~~latex
against an observed 0.650, $1-0.229=0.771$ against 0.750, and $1-0.047=0.953$ against 0.921, the residual being the few percent of surviving records that still go uncovered.
~~~

Source and values:

- Observed FL gene-capture coverage (crude): random 0.650422, k-means 0.749631, original 0.920721 (three_panel_comparison_d100.json); emptied shares rule-invariant.

Note: Wording: under the union part of the residual is surviving records that select references without gene families (107/76/62), not abstentions.

### 1.85. L961 (Panel composition) — QUALITATIVE

Current:

~~~latex
plus any cells whose selected references carry no gene families, of which there are none on any of the three panels under the legacy reducer.)
~~~

Replacement:

~~~latex
plus cells whose selected references carry no gene families, of which the original, $k$-means and random panels have 107, 76 and 62 among their surviving full-length records.)
~~~

Source and values:

- `emptied_strata_d100.json` panels.<p>.survivors.non_ok - .strict: 236-129 = 107, 188-112 = 76, 137-75 = 62 (= pred_empty survivor cells in each tree's gene_capture_cells.csv, FullLength16S exclude_0.987). Over all 8 regions: 1,060 / 597 / 416.

Note: QUALITATIVE: 'none ... under the legacy reducer' becomes a non-zero count under the union.

### 1.86. L966 (Panel composition)

Current:

~~~latex
uncovered rates under the legacy reducer.}
~~~

Replacement:

~~~latex
uncovered rates under the exact-tie union.}
~~~

Source and values:

- Rule statement (tab:emptied caption).

### 1.87. L975 (Panel composition)

Current:

~~~latex
0--19                      & 84.4\% & 0.031 & 44.0\% & 0.023 & 37.3\% & 0.019\\
~~~

Replacement:

~~~latex
0--19                      & 84.4\% & 0.037 & 44.0\% & 0.028 & 37.3\% & 0.023\\
~~~

Source and values:

- `three_panel/emptied_strata_d100.json` panels.<p>.strata[stratum].non_ok_rate (4 dp), recomputed exactly from each union tree's gene_capture_cells.csv (FullLength16S, exclude_0.987, outcome != truth_empty, binned on n_removed_self). Shares are rule-invariant and unchanged; the emptied row (1.000) is unchanged.
- Unrounded: uncovered original 223/6,087 = 0.03664; k-means 110/3,880 = 0.02835; random 74/3,269 = 0.02264

### 1.88. L976 (Panel composition)

Current:

~~~latex
20--49                     & 7.2\%  & 0.012 & 17.9\% & 0.018 & 14.8\% & 0.017\\
~~~

Replacement:

~~~latex
20--49                     & 7.2\%  & 0.017 & 17.9\% & 0.022 & 14.8\% & 0.021\\
~~~

Source and values:

- `three_panel/emptied_strata_d100.json` panels.<p>.strata[stratum].non_ok_rate (4 dp), recomputed exactly from each union tree's gene_capture_cells.csv (FullLength16S, exclude_0.987, outcome != truth_empty, binned on n_removed_self). Shares are rule-invariant and unchanged; the emptied row (1.000) is unchanged.
- Unrounded: 9/520 = 0.01731; 35/1,577 = 0.02219; 27/1,298 = 0.02080

### 1.89. L977 (Panel composition)

Current:

~~~latex
50--99                     & 3.8\%  & 0.004 & 15.2\% & 0.028 & 14.5\% & 0.020\\
~~~

Replacement:

~~~latex
50--99                     & 3.8\%  & 0.015 & 15.2\% & 0.032 & 14.5\% & 0.028\\
~~~

Source and values:

- `three_panel/emptied_strata_d100.json` panels.<p>.strata[stratum].non_ok_rate (4 dp), recomputed exactly from each union tree's gene_capture_cells.csv (FullLength16S, exclude_0.987, outcome != truth_empty, binned on n_removed_self). Shares are rule-invariant and unchanged; the emptied row (1.000) is unchanged.
- Unrounded: 4/272 = 0.01471; 43/1,339 = 0.03211; 36/1,269 = 0.02837

### 1.90. L993 (Panel composition)

Current:

~~~latex
genus coverage on V4 is 0.641, 0.441 and 0.368 on the original, $k$-means and random panels, a random-minus-original difference of $-0.273$ ($-0.289$ to $-0.254$) against $-0.342$ crude.
~~~

Replacement:

~~~latex
genus coverage on V4 is 0.615, 0.431 and 0.358 on the original, $k$-means and random panels, a random-minus-original difference of $-0.257$ ($-0.276$ to $-0.239$) against $-0.330$ crude.
~~~

Source and values:

- `three_panel_comparison_d100.json` V4 exclude_0.987 genus_coverage phystd 0.615292 / 0.430820 / 0.358027; delta_random_minus_original est -0.2573, ci [-0.2758, -0.2393]; crude random - original -0.329636.

### 1.91. L994 (Panel composition)

Current:

~~~latex
keeps 76--94\% of the crude one
~~~

Replacement:

~~~latex
keeps 76--92\% of the crude one
~~~

Source and values:

- Ratio delta_random_minus_original.est / (random.crude - original.crude), exclude_0.987: genus_coverage 0.764 (V4-V5) to 0.906 (V1-V3); genecap_coverage 0.828 (V7-V9) to 0.921 (V1-V3). All 16 CIs exclude zero.

### 1.92. L996 (Panel composition) — QUALITATIVE

Current:

~~~latex
conditional genus accuracy under self-exclusion is higher on the random panel than on the original in every region, by 0.002 to 0.074, with intervals that exclude zero in five of eight (V1--V2 and V4--V5 only marginally), and conditional gene-capture F1 is lower on the random panel by 0.014 to 0.027
~~~

Replacement:

~~~latex
conditional genus accuracy under self-exclusion differs between the random panel and the original by $-0.016$ to $+0.062$ (higher on the random panel in six of eight regions), with intervals that exclude zero only for full-length 16S and V1--V3, and conditional gene-capture F1 is lower on the random panel by 0.013 to 0.028
~~~

Source and values:

- `three_panel_comparison_d100.json` exclude_0.987 genus_cond_acc delta_random_minus_original.est: FL +0.0377 [0.0125, 0.0609], V1-V2 +0.0156, V1-V3 +0.0622 [0.0387, 0.0856], V3-V4 -0.0100, V4 +0.0161, V4-V5 +0.0149, V4-V5_944R +0.0245, V7-V9 -0.0161; excludes_zero true only for FL and V1-V3.
- genecap_cond_f1 delta_random_minus_original.est -0.0128 (V1-V2) to -0.0278 (V7-V9).

Note: QUALITATIVE: the uniform-sign statement ('higher ... in every region') fails under the union; the 'accuracy moves far less than coverage' conclusion holds.

### 1.93. L1004 (Panel composition)

Current:

~~~latex
The 16S blocks score all three panels at record depth 100 under the legacy reducer.}
~~~

Replacement:

~~~latex
The 16S blocks score all three panels at record depth 100 under the exact-tie union.}
~~~

Source and values:

- Rule statement (tab:threepanel caption).

### 1.94. L1014 (Panel composition)

Current:

~~~latex
full-length  & 0.728 & 0.558 & 0.493 & $-0.235$ \scriptsize[$-0.265$, $-0.203$] & $+0.045$ \scriptsize[$+0.022$, $+0.066$] \\
~~~

Replacement:

~~~latex
full-length  & 0.700 & 0.544 & 0.479 & $-0.220$ \scriptsize[$-0.251$, $-0.191$] & $+0.038$ \scriptsize[$+0.013$, $+0.061$] \\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987]: columns 2-4 = {original,kmeans,random}.phystd of genus_coverage (genus block) or genecap_coverage (gene-capture block); column 5 = delta_random_minus_original.est [ci] of the same endpoint; column 6 = delta_random_minus_original.est [ci] of genus_cond_acc / genecap_cond_f1. CI bounds are stored at 4 dp; half-way bounds (marked *) are rounded away from zero, the convention the current table uses; their third decimal is Monte Carlo noise (PYTHONHASHSEED-dependent phylum order).
- Unrounded: 0.699818 / 0.543770 / 0.479476; -0.2203 [-0.2507, -0.1905*]; +0.0377 [0.0125*, 0.0609]

### 1.95. L1015 (Panel composition)

Current:

~~~latex
V1--V3       & 0.726 & 0.584 & 0.488 & $-0.239$ \scriptsize[$-0.261$, $-0.216$] & $+0.074$ \scriptsize[$+0.050$, $+0.095$] \\
~~~

Replacement:

~~~latex
V1--V3       & 0.688 & 0.546 & 0.473 & $-0.215$ \scriptsize[$-0.237$, $-0.194$] & $+0.062$ \scriptsize[$+0.039$, $+0.086$] \\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987]: columns 2-4 = {original,kmeans,random}.phystd of genus_coverage (genus block) or genecap_coverage (gene-capture block); column 5 = delta_random_minus_original.est [ci] of the same endpoint; column 6 = delta_random_minus_original.est [ci] of genus_cond_acc / genecap_cond_f1. CI bounds are stored at 4 dp; half-way bounds (marked *) are rounded away from zero, the convention the current table uses; their third decimal is Monte Carlo noise (PYTHONHASHSEED-dependent phylum order).
- Unrounded: 0.688417 / 0.545752 / 0.473250; -0.2152 [-0.2368, -0.1936]; +0.0622 [0.0387, 0.0856]

### 1.96. L1016 (Panel composition)

Current:

~~~latex
V3--V4       & 0.675 & 0.486 & 0.417 & $-0.258$ \scriptsize[$-0.285$, $-0.226$] & $+0.002$ \scriptsize[$-0.029$, $+0.029$] \\
~~~

Replacement:

~~~latex
V3--V4       & 0.643 & 0.465 & 0.397 & $-0.246$ \scriptsize[$-0.274$, $-0.217$] & $-0.010$ \scriptsize[$-0.039$, $+0.019$] \\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987]: columns 2-4 = {original,kmeans,random}.phystd of genus_coverage (genus block) or genecap_coverage (gene-capture block); column 5 = delta_random_minus_original.est [ci] of the same endpoint; column 6 = delta_random_minus_original.est [ci] of genus_cond_acc / genecap_cond_f1. CI bounds are stored at 4 dp; half-way bounds (marked *) are rounded away from zero, the convention the current table uses; their third decimal is Monte Carlo noise (PYTHONHASHSEED-dependent phylum order).
- Unrounded: 0.643326 / 0.465323 / 0.397470; -0.2459 [-0.274, -0.2168]; -0.0100 [-0.0385*, 0.0188]

### 1.97. L1017 (Panel composition)

Current:

~~~latex
V4           & 0.641 & 0.441 & 0.368 & $-0.273$ \scriptsize[$-0.289$, $-0.254$] & $+0.029$ \scriptsize[$-0.003$, $+0.058$] \\
~~~

Replacement:

~~~latex
V4           & 0.615 & 0.431 & 0.358 & $-0.257$ \scriptsize[$-0.276$, $-0.239$] & $+0.016$ \scriptsize[$-0.015$, $+0.045$] \\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987]: columns 2-4 = {original,kmeans,random}.phystd of genus_coverage (genus block) or genecap_coverage (gene-capture block); column 5 = delta_random_minus_original.est [ci] of the same endpoint; column 6 = delta_random_minus_original.est [ci] of genus_cond_acc / genecap_cond_f1. CI bounds are stored at 4 dp; half-way bounds (marked *) are rounded away from zero, the convention the current table uses; their third decimal is Monte Carlo noise (PYTHONHASHSEED-dependent phylum order).
- Unrounded: 0.615292 / 0.430820 / 0.358027; -0.2573 [-0.2758, -0.2393]; +0.0161 [-0.0146, 0.0452]

### 1.98. L1018 (Panel composition)

Current:

~~~latex
V7--V9       & 0.678 & 0.486 & 0.405 & $-0.272$ \scriptsize[$-0.291$, $-0.252$] & $+0.014$ \scriptsize[$-0.015$, $+0.043$] \\
~~~

Replacement:

~~~latex
V7--V9       & 0.640 & 0.457 & 0.382 & $-0.257$ \scriptsize[$-0.276$, $-0.237$] & $-0.016$ \scriptsize[$-0.043$, $+0.013$] \\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987]: columns 2-4 = {original,kmeans,random}.phystd of genus_coverage (genus block) or genecap_coverage (gene-capture block); column 5 = delta_random_minus_original.est [ci] of the same endpoint; column 6 = delta_random_minus_original.est [ci] of genus_cond_acc / genecap_cond_f1. CI bounds are stored at 4 dp; half-way bounds (marked *) are rounded away from zero, the convention the current table uses; their third decimal is Monte Carlo noise (PYTHONHASHSEED-dependent phylum order).
- Unrounded: 0.639869 / 0.456636 / 0.382442; -0.2574 [-0.2759, -0.237]; -0.0161 [-0.0427, 0.0125*]

### 1.99. L1021 (Panel composition)

Current:

~~~latex
full-length  & 0.931 & 0.770 & 0.686 & $-0.245$ \scriptsize[$-0.256$, $-0.235$] & $-0.021$ \scriptsize[$-0.033$, $-0.007$] \\
~~~

Replacement:

~~~latex
full-length  & 0.926 & 0.767 & 0.682 & $-0.243$ \scriptsize[$-0.254$, $-0.233$] & $-0.019$ \scriptsize[$-0.031$, $-0.005$] \\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987]: columns 2-4 = {original,kmeans,random}.phystd of genus_coverage (genus block) or genecap_coverage (gene-capture block); column 5 = delta_random_minus_original.est [ci] of the same endpoint; column 6 = delta_random_minus_original.est [ci] of genus_cond_acc / genecap_cond_f1. CI bounds are stored at 4 dp; half-way bounds (marked *) are rounded away from zero, the convention the current table uses; their third decimal is Monte Carlo noise (PYTHONHASHSEED-dependent phylum order).
- Unrounded: 0.925712 / 0.766921 / 0.682440; -0.2433 [-0.2541, -0.2329]; -0.0189 [-0.0309, -0.005]

### 1.100. L1022 (Panel composition)

Current:

~~~latex
V4           & 0.857 & 0.666 & 0.579 & $-0.278$ \scriptsize[$-0.289$, $-0.267$] & $-0.024$ \scriptsize[$-0.031$, $-0.017$] \\
~~~

Replacement:

~~~latex
V4           & 0.855 & 0.665 & 0.576 & $-0.279$ \scriptsize[$-0.291$, $-0.268$] & $-0.024$ \scriptsize[$-0.031$, $-0.017$] \\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987]: columns 2-4 = {original,kmeans,random}.phystd of genus_coverage (genus block) or genecap_coverage (gene-capture block); column 5 = delta_random_minus_original.est [ci] of the same endpoint; column 6 = delta_random_minus_original.est [ci] of genus_cond_acc / genecap_cond_f1. CI bounds are stored at 4 dp; half-way bounds (marked *) are rounded away from zero, the convention the current table uses; their third decimal is Monte Carlo noise (PYTHONHASHSEED-dependent phylum order).
- Unrounded: 0.855448 / 0.665352 / 0.576364; -0.2791 [-0.2905*, -0.268]; -0.0243 [-0.0313, -0.0168]

### 1.101. L1023 (Panel composition)

Current:

~~~latex
V7--V9       & 0.873 & 0.685 & 0.598 & $-0.275$ \scriptsize[$-0.286$, $-0.263$] & $-0.027$ \scriptsize[$-0.033$, $-0.020$] \\
~~~

Replacement:

~~~latex
V7--V9       & 0.870 & 0.681 & 0.594 & $-0.275$ \scriptsize[$-0.287$, $-0.263$] & $-0.028$ \scriptsize[$-0.035$, $-0.021$] \\
~~~

Source and values:

- `three_panel/three_panel_comparison_d100.json` cells[regime=exclude_0.987]: columns 2-4 = {original,kmeans,random}.phystd of genus_coverage (genus block) or genecap_coverage (gene-capture block); column 5 = delta_random_minus_original.est [ci] of the same endpoint; column 6 = delta_random_minus_original.est [ci] of genus_cond_acc / genecap_cond_f1. CI bounds are stored at 4 dp; half-way bounds (marked *) are rounded away from zero, the convention the current table uses; their third decimal is Monte Carlo noise (PYTHONHASHSEED-dependent phylum order).
- Unrounded: 0.869580 / 0.680991 / 0.594261; -0.2753 [-0.2869, -0.2631]; -0.0278 [-0.0351, -0.0209]

### 1.102. L1045 (Panel composition) — QUALITATIVE

Current:

~~~latex
Standardized self-inclusion gene-capture F1 falls from 0.887 to 0.839 on V4 ($-0.048$, $-0.053$ to $-0.043$) and from 0.915 to 0.866 on V7--V9, and genus accuracy by 0.01.
~~~

Replacement:

~~~latex
Standardized self-inclusion gene-capture F1 falls from 0.837 to 0.730 on V4 ($-0.107$, $-0.115$ to $-0.100$) and from 0.906 to 0.805 on V7--V9, and genus accuracy by 0.07 to 0.08.
~~~

Source and values:

- `three_panel_comparison_d100.json` cells[regime=include]: genecap_cond_f1 V4 original.phystd 0.836842, random.phystd 0.729982, delta_random_minus_original -0.1069 [-0.1145*, -0.0998]; V7-V9 0.906473 -> 0.805086 (-0.1014); genus_cond_acc delta_random_minus_original V4 -0.0792, V7-V9 -0.0653.
- L100 check: V4 0.887142 -> 0.839496 (-0.0476 [-0.0528, -0.043]); V7-V9 0.915285 -> 0.865565; genus -0.0121 / -0.0106.

Note: QUALITATIVE (magnitude): the self-inclusion panel penalty more than doubles; the genus penalty rises from 0.01 to 0.07-0.08.

### 1.103. L1046 (Panel composition) — QUALITATIVE

Current:

~~~latex
On full-length 16S gene-capture F1 falls by less than half as much ($-0.021$) and genus accuracy not at all ($+0.004$)
~~~

Replacement:

~~~latex
On full-length 16S gene-capture F1 falls by less than a sixth as much ($-0.016$) and genus accuracy barely at all ($-0.006$)
~~~

Source and values:

- include FullLength16S genecap_cond_f1 delta_random_minus_original -0.0159 [-0.0176, -0.014] (0.0159/0.1069 = 0.149); genus_cond_acc -0.0064 [-0.0089, -0.0037] (CI excludes zero).

Note: QUALITATIVE: FL genus accuracy now falls slightly (CI excludes zero) rather than 'not at all'.

### 1.104. L1059 (Panel composition)

Current:

~~~latex
genus and gene-capture coverage falls by 0.22 to 0.29 against the stratified panel
~~~

Replacement:

~~~latex
genus and gene-capture coverage falls by 0.20 to 0.29 against the stratified panel
~~~

Source and values:

- As item L69.

### 1.105. L1061 (Panel composition)

Current:

~~~latex
differs by at most 0.027 between any two panels in either arm
~~~

Replacement:

~~~latex
differs by at most 0.028 between any two panels in either arm
~~~

Source and values:

- Max |pairwise Δ| of exclude_0.987 genecap_cond_f1 phystd over RO/KO/RK and 8 regions: 0.0278 (RO V7-V9); MAG panel arm (three_panel_mag_panelarm.json, exclude) max 0.0262. Still an order of magnitude below the 0.20-0.29 coverage gap.

### 1.106. L1073 (Discussion) — QUALITATIVE

Current:

~~~latex
even a short marker recovers nearly all of its metabolism (95--97\% gene capture), whereas for genuinely novel organisms an amplicon from V3 onward (V4, V4--V5, V3--V4, V7--V9) is right at the family level about four times in five when it makes a call (0.80--0.83 conditional) but, over the self-excluded population, receives no genus call for roughly 30\% (abstentions plus answered amplicons whose reference has no named genus), captures only about three fifths of the gene complement (59.8--62.7\%) and places 59--65\% of organisms in the correct family and 40--46\% in the correct genus --- a gap that full-length 16S and the two V1-anchored regions (V1--V2, V1--V3) measurably narrow, to 66.1--67.5\% gene capture and 0.54--0.56 population genus accuracy.
~~~

Replacement:

~~~latex
even a short marker recovers nearly all of its metabolism (over 99\% gene capture, at macro precision 0.79--0.88), whereas for genuinely novel organisms an amplicon from V3 onward (V4, V4--V5, V3--V4, V7--V9) is right at the family level more than three times in four when it makes a call (0.76--0.80 conditional) but, over the self-excluded population, receives no genus call for roughly a third (abstentions plus answered amplicons whose reference has no named genus), captures less than three fifths of the gene complement (57.6--59.2\%) and places 54--59\% of organisms in the correct family and 37--42\% in the correct genus --- a gap that full-length 16S and the two V1-anchored regions (V1--V2, V1--V3) measurably narrow, to 63.3--64.4\% gene capture and 0.51--0.54 population genus accuracy.
~~~

Source and values:

- V3-onward regions = V3-V4, V4, V4-V5, V4-V5_944R, V7-V9 (U100). headline_macro_recall_pct.include 99.148-99.396; regions[r].include.macro.precision 0.7945 (V4)-0.8799 (V7-V9).
- taxacc anchor.exclude_0.987: Family.correct 0.759021-0.796716; Genus.coverage 0.642956-0.676344 (uncovered 32.4-35.7%); Family.correct_population 0.53587-0.58690; Genus.correct_population 0.370554-0.423398.
- headline_macro_recall_pct_population.exclude_0.987 57.643-59.214 (V3 onward); FL/V1-V2/V1-V3 63.275/64.310/64.443; Genus.correct_population 0.537691/0.507398/0.512821.

Note: QUALITATIVE: the 'nearly all of its metabolism' clause now needs the precision qualifier; family 'four times in five' becomes 'more than three in four'; genus non-call 'roughly 30%' becomes 'roughly a third'.

### 1.107. L1101 (Limitations)

Current:

~~~latex
which makes every self-exclusion coverage and population-level figure here, as-shipped or legacy-scored, a lower bound.
~~~

Replacement:

~~~latex
which makes every self-exclusion coverage and population-level figure here a lower bound.
~~~

Source and values:

- Removes the 'legacy-scored' category (no longer used).

### 1.108. L1102 (Limitations)

Current:

~~~latex
(as-shipped conditional EC F1 0.802 against 0.798; conditional genus accuracy under the legacy reducer within 0.01 in every region, and conditional gene capture, which is macro recall, up 0.7 to 1.3 points while conditional gene-capture F1 fell by 0.1 to 1.2 points).
~~~

Replacement:

~~~latex
(conditional EC F1 0.802 against 0.798).
~~~

Source and values:

- Keeps only the union-at-both-depths comparison (`region_validation{,_d100}/picrust2_comparison/data/comparison_summary.json` by_region.ec.all.all.prfba_exclude.conditional.f1: 0.7980 -> 0.8017).

Note: UNRESOLVED for a union d20->d100 tab:main statement (section 4).

## 1b. Optional items

Apply only if wanted, and only after sections 1 and 3 (item O1 edits the text that item 1 for L64 produces). Each old substring was checked to occur exactly once in the edited text.

### O1. L64 (Abstract)

Current:

~~~latex
When the source organism is in the database, genus accuracy is 0.89--1.00 and gene capture 99--100\%.
~~~

Replacement:

~~~latex
When the source organism is in the database, genus accuracy is 0.89--1.00 and gene capture 99--100\%, although the larger tie sets lower gene-capture precision on short amplicons (0.79 on V4).
~~~

Source and values:

- Apply AFTER section-1 item for L64 (old = that item's new text). gene_capture_by_region.json (U100) regions.V4.include.macro.precision.value 0.794532.

Note: OPTIONAL qualitative disclosure.

### O2. L352 (10,000-organism validation)

Current:

~~~latex
even a short amplicon recovers nearly all of its metabolism (\Cref{fig:accuracy}, \Cref{fig:genecap}).
~~~

Replacement:

~~~latex
even a short amplicon recovers nearly all of its metabolism (\Cref{fig:accuracy}, \Cref{fig:genecap}), although the larger exact-tie selections lower its gene-capture precision (macro precision 0.79 on V4, against 1.00 on full-length 16S).
~~~

Source and values:

- regions[r].include.macro.precision.value: V4 0.794532, FullLength16S 0.996716 (U100).

Note: OPTIONAL qualitative disclosure.

### O3. L349 (10,000-organism validation)

Current:

~~~latex
Under the exact-tie union only 13.3\% of the abstentions
~~~

Replacement:

~~~latex
Only 13.3\% of the abstentions
~~~

Source and values:

- Redundant once the whole paper is union-scored.

Note: OPTIONAL cleanup.

### O4. L611 (Retention depth)

Current:

~~~latex
Under the shipped exact-tie union, which every as-shipped statistic in this section uses unless a rule is named, 6{,}457 amplicons
~~~

Replacement:

~~~latex
Under the shipped exact-tie union, 6{,}457 amplicons
~~~

Source and values:

- After the section-1 deletions no other rule is named in sec:depth.

Note: OPTIONAL cleanup.

## 1c. Legacy mentions that stay unchanged

**Controlled-benchmark mentions** (the authors' exception). Section 3 adds the exact-tie value next to each number. They are at:

- L76 ("which picked legacy in every fold");
- L305;
- L316 ("61.7% of the same amplicons under the legacy reducer that leave-one-fold-out selection selects there");
- L506–508;
- L527–528;
- L541 ("+0.016 like-for-like against the legacy policy's own masked forced run");
- L547, L1068 and L1123.

**Methods contrasts with the alternative** (L181 and L183). These are not headline or table numbers. Keep them as the "available alternative" context, or cut them if the authors want only the new one-sentence mention (item 1.9).

- L181, "keeps gene-capture precision under self-exclusion at or above the ordered reducer's". Verified: U100 ≥ L100 conditional precision in all 8 regions (V4 0.5360 against 0.5343, the narrowest).
- L181, "(mean 2.88 genomes against 1.32 for the ordered reducer ...)". Verified: `gene_capture_cells.csv`, include, outcome ok, mean n_selected 2.8836 (U100) against 1.3176 (L100).
- L181, "(V4 0.795 against 0.856)". Verified: `regions.V4.include.macro.precision.value` 0.7945 against 0.8565.
- L183, "∼19% of self-recovery round-trips (12,602 of the 67,466 ...)". This is the legacy self-recovery loss that motivated the union.

**Already union, verified, and unchanged:**

- L1082, "a baseline its shipped, abstaining arm does not reach". This now literally names the exact-tie arm, and the claim holds for it (0.685 against nn1 0.745), so no edit is needed.
- L172: the 34.0% share. U100 include share with n_selected > 1 among answered cells is 0.3404.
- Abstract L66, L73–74 and L78.
- L180, L209 and L289.
- L349: the 13.3%, from `prfba_selection.json`.
- sec:picrust2 L428–458, tab:picrust2 and the fig:picrust2 caption. All come from the `prfba_exclude` key.
- sec:depth L609–630 and L634–640.
- L1077–1079 and L1099–1101.

**Rule-independent and unchanged:**

- the self-recovery shares at L1038–1047 and L1062 (computed from the hits);
- the emptied shares and the emptied row of tab:emptied;
- every MAG-only and GTDB number;
- tab:mag and tab:panelcomp;
- fig:panelspace and fig:embedding.

## 2. Figures

There are nine `\includegraphics` in prFBA.tex. Each paper file was matched to its render by md5, not by size. The current renders were made from patched copies of the repo scripts in this session's scratchpad. Scratchpad root:

`/tmp/claude-1000/-home-freiburger-Documents-prFBA/0c6f5bc6-8224-5a34-9b20-7d5643487901/scratchpad/`

Paths below that start with `ms_figs_` are under it. That directory is not durable, so the union re-renders should be committed as real scripts.

### 2.1 Six figures to re-render from union inputs

| Figure (label, line) | Paper file | Rendered by | Reads now (legacy) | Switch to (union) |
|---|---|---|---|---|
| fig:accuracy, L390 | `images/prfba/F1_accuracy_heatmap.pdf` | `ms_figs_scripts/rv_legacy/render_accuracy.py` (render_F1, include-self), via patched `ms_figs_scripts/rv_legacy/_common.py`; md5 53bb5e0d = `ms_figs_scripts/raw_out/rv_legacy/F1_accuracy_heatmap.pdf` | `BASE = .../region_validation_d100_legacy` → `data/concordance_long.csv` (fallback `asv_region_concordance.csv`) | `BASE = /home/freiburger/Documents/prFBA/region_validation_d100` |
| fig:genecap, L398 | `images/prfba/F4b_percent_genecapture_by_region.pdf` | `ms_figs_scripts/rv_legacy/render_genecap.py` (fig_F4b) with two scratch patches; md5 87fe2367 = `ms_figs_scripts/raw_out/rv_legacy/F4b_percent_genecapture_by_region.pdf` | `BASE` legacy → `data/gene_capture_cells.csv` | `BASE` = `region_validation_d100`. Keep both patches |
| fig:leakage, L406 | `images/prfba/F11_self_inclusion_exclusion_delta.pdf` | `ms_figs_scripts/rv_legacy/render_accuracy.py` (F11: include minus exclude of `correct`); md5 49b14fa0 = `ms_figs_scripts/raw_out/rv_legacy/F11_self_inclusion_exclusion_delta.pdf` | `BASE` legacy → `data/concordance_long.csv` | `BASE` = `region_validation_d100` |
| fig:perrank, L414 | `images/prfba/F2_perrank_concordance.pdf` | `ms_figs_scripts/rv_legacy/render_accuracy_exclude_alt.py` (render_F2 with `_filter_selfmode(raw, "exclude")`), run by `ms_figs_scripts/rv_legacy/run_exclude_alt.py`; md5 7a247418 = `ms_figs_alternates/F2_perrank_concordance_exclude_self_ALT.pdf` | `BASE` legacy → `data/concordance_long.csv`, exclude rows | `BASE` = `region_validation_d100` |
| fig:magcomp, L821 | `images/prfba/F14_mag_completeness.pdf` | `ms_figs_scripts/mag/plot_results_d100_legacy.py` (patched copy of `mag_validation/plot_results.py`); md5 19d9caf4 = `ms_figs_scripts/raw_out/mag/F_mag_completeness.pdf` | panel C: `mag_validation/out/paired_comparison_d100_legacy.json` (line 24 of the script) | `mag_validation/out/paired_comparison_d100.json` |
| fig:magcross, L831 | `images/prfba/F15_mag_crossover.pdf` | `ms_figs_scripts/mag/plot_crossover_d100_legacy.py` (patched copy of `mag_validation/plot_crossover.py`); md5 f4862da3 = `ms_figs_scripts/raw_out/mag/F_mag_crossover.pdf` | `mag_validation/out/crossover_d100_legacy.json` (line 21) | `mag_validation/out/crossover_d100.json` |

What each figure needs beyond the input switch:

- **fig:accuracy.**
  - A current union render already exists: `region_validation_d100/figures/F1_accuracy_heatmap.pdf`.
  - It was rendered at 2026-10-20 22:53 by `region_validation_d100/scripts/render_accuracy.py`, which is identical code (only `FIGS` differs). That is after its input (`concordance_long.csv`, 22:44) and after the gap close (22:40).
  - The caption change is item 1.32.
- **fig:genecap.**
  - Keep both scratch patches: `ax.set_ylim(0, min(122, top + 24))` and `legend(loc="upper center", ncol=2)`. The include-self ticks now sit at 99–100%, so the headroom matters more.
  - `region_validation_d100/figures/F4b_percent_genecapture_by_region.pdf` is unpatched; don't copy it.
  - The bars count abstained and `pred_empty` cells as 0, so they show the 58–64% population means.
  - The caption change is item 1.33.
- **fig:leakage.**
  - A current union render exists at `region_validation_d100/figures/F11_self_inclusion_exclusion_delta.pdf` (same code).
  - The caption change is item 1.34.
- **fig:perrank.**
  - No union render of this exclude-self variant exists.
  - `region_validation_d100/figures/F2_perrank_concordance.pdf` is the include-self version and must not be used.
  - The caption change is item 1.35.
- **fig:magcomp.** Two hard-coded strings must change:
  - footer: `(depth-100 records, legacy reducer)` → `(depth-100 records, exact-tie union)`;
  - panel C annotation: `$p<10^{-179}$` → `$p<10^{-164}$`. The union full-length p is 4.0e-165; V4 underflows.
  - "+0.05 easier" stays (the union value is +0.053).
  - The computed annotation becomes "+0.103 vs full-length 16S, +0.212 vs V4", and the computed n range becomes 1,912–2,299.
  - The panel C y-limits (0.55, 0.93) still contain the union V4 series (0.593–0.598).
  - Panels A and B are MAG-only and unchanged. The caption needs no change.
- **fig:magcross.**
  - No rule text is hard-coded.
  - The computed annotation becomes "still +0.077 at 1% complete (+0.103 at 100%)", and the footer n range becomes 1,603–2,311.
  - There is still no crossover: every Δ is positive.
  - The caption needs no change.

### 2.2 Three figures that need no change

| Figure (label, line) | Paper file | Why it is unchanged |
|---|---|---|
| fig:picrust2, L490 | `images/prfba/F13_prfba_vs_picrust2.pdf` | md5 c94ec692 = `ms_figs_scripts/raw_out/p2/F1_prfba_vs_picrust2.pdf`, from `ms_figs_scripts/p2/render_comparison_d100.py`, which reads `comparison_summary.json` key `prfba_exclude` (the union). The median-ΔF1 title is computed. |
| fig:panelspace, L896 | `images/prfba/F16_panel_space.pdf` | md5 e3a53986 = `two_panel/F_panel_space.pdf` (`two_panel/plot_space.py`). Embedding density only; rule-independent. |
| fig:embedding, L422 | `images/prfba/fig_pairwise_hist.png` | md5 f30d2ad8 = `/home/freiburger/Documents/prFBA/fig_pairwise_hist.png`. Embedding only; no selection. |

### 2.3 Re-render recipe

1. Copy `ms_figs_scripts/rv_legacy/` and set two lines in `_common.py`:
   - `BASE = "/home/freiburger/Documents/prFBA/region_validation_d100"`;
   - `FIGS` = a scratch output directory.
2. Render F1 and F11 with `render_accuracy.py` and F4b with `render_genecap.py`, the way the current renders were made.
3. Render F2 with `run_exclude_alt.py`, after pointing its `C.FIGS` at the same scratch directory.
4. Copy `mag/plot_results_d100_legacy.py` and `mag/plot_crossover_d100_legacy.py`. In each, drop `_legacy` from the input JSON name and make the two string edits above.
5. Copy the PDFs to `images/prfba/` under the manuscript's file names.

## 3. Controlled benchmark: add the exact-tie arm alongside

### What stays

The pre-registered headline stays. Its primary prFBA arm is chosen by leave-one-fold-out ("nested CV") selection (`arms.prfba_nestedcv`). That selection picked legacy in all 5 folds at both ranks, so the arm equals `arms.prfba_legacy` to the last digit. Every item below keeps that number and adds the shipped exact-tie value (`arms.prfba_exact_tie`, `arms.prfba_exact_tie_forced`) next to it.

### Sources and method

- Genus rank: `region_validation/shared_reference/data/summary_genus.json`.
- Family rank: `region_validation/shared_reference/data_family/summary_family.json`.
- Both summaries store values at 4 dp. Every 3-dp value below was recomputed from the cells to avoid double rounding: `data/cells.csv.gz` and `data_family/cells.csv.gz`, space `ec`, density 1.0, `unscoreable_min_align` dropped, population F1 = f1 if outcome `scored` else 0. For example, exact-tie forced 0.7455 is 0.745538 → 0.746, not 0.745.
- The summaries have no exact-tie contrast against nn1, random abstain, constant or family donor. Those differences are point estimates from the cells, labelled "no interval computed".

### Side-by-side values

| Quantity | nested-CV (kept) | exact-tie (added) |
|---|---|---|
| **Genus rank** | | |
| population EC F1 | 0.687 (0.686676) | 0.685 (0.684933) |
| coverage | 0.893 | 0.893 (identical) |
| Δ vs PICRUSt2 [95% CI] | −0.092 [−0.107, −0.081] | −0.094 [−0.109, −0.083] |
| forced: population F1 | 0.747 (0.747258) | 0.746 (0.745538) |
| forced: Δ vs PICRUSt2 [95% CI] | −0.031 [−0.037, −0.028] | −0.033 [−0.039, −0.029] |
| conditional F1, matched 57,021 cells | 0.769 | 0.767 (0.767215) |
| skill, abstaining / forced | 0.28 / 0.70 | 0.27 / 0.69 |
| AURC / cross-knob | 0.788 / 0.788 | 0.787 / 0.786 |
| vs nn1 | −0.058 [−0.069, −0.051] | −0.060 (no interval) |
| vs random abstain | +0.021 [0.014, 0.030] | +0.019 (no interval) |
| regions, abstaining / forced | −0.059…−0.182 / −0.029…−0.034 | −0.060…−0.184 / −0.031…−0.036 |
| domains, abstaining (Bacteria, Archaea) | −0.095, −0.052 | −0.097, −0.053 |
| domains, forced (Bacteria, Archaea) | −0.032, −0.019 | −0.034, −0.020 |
| density slope per log10 tips | 0.107 | 0.106 |
| quartering cost | 0.064 | 0.064 |
| full-length / V4 population F1 | 0.692 / 0.710 | 0.690 / 0.708 |
| band below 0.865, forced population F1 | 0.567 | 0.567 |
| precision in the 0.945–0.987 / 0.865–0.945 bands | 0.783 / 0.740 | 0.803 / 0.751 |
| **Family rank** | | |
| population EC F1 | 0.626 | 0.623 |
| coverage | 0.856 | 0.856 |
| Δ vs PICRUSt2 [95% CI] | −0.124 [−0.150, −0.099] | −0.126 [−0.152, −0.102] |
| forced: population F1 | 0.709 | 0.707 |
| forced: Δ vs PICRUSt2 [95% CI] | −0.040 [−0.048, −0.031] | −0.042 [−0.050, −0.033] |
| regions | −0.076…−0.245 | −0.077…−0.247 |
| forced, vs nn1 (0.7070) | 0.7095 | 0.7072 |
| vs nn1 | −0.081 [−0.107, −0.059] | −0.083 (no interval) |
| vs constant (prevalence floor) | −0.017 [−0.053, 0.018] | −0.020 (no interval) |
| skill | −0.14 | −0.16 |
| vs random abstain | +0.018 [−0.004, 0.040] | +0.015 (no interval) |
| vs family donor | +0.107 [0.044, 0.179] | +0.104 (no interval) |

### Places that already report the exact-tie value

These need nothing:

- L78, L547–560, L1083–1084 and L1091: the phenotype arm is the exact-tie arm.
- L209, L289, L316, L1068, L1123 and L1130: the selection-size and tie-set shares are the exact-tie values.
- L508: "the shipped exact-tie default scores 0.002 lower in population F1 at identical coverage". Verified: 0.686676 against 0.684933, coverage 0.892753 for both.
- L509: the random-abstain arm is built from exact-tie forced predictions.
- L527–528: already give +0.019 for the exact-tie policy.
- L564, the C2 control: a selection-level pass with no prFBA-arm number.

### Items


### 3.1. L75 (Abstract)

Current:

~~~latex
against prFBA's 0.687 ($\Delta -0.092$, 95\% CI $[-0.107, -0.081]$), and against 0.747 for prFBA forced to answer everywhere ($\Delta -0.031$); both intervals lie wholly beyond
~~~

Replacement:

~~~latex
against prFBA's 0.687 ($\Delta -0.092$, 95\% CI $[-0.107, -0.081]$; shipped exact-tie rule 0.685, $\Delta -0.094$, $[-0.109, -0.083]$), and against 0.747 for prFBA forced to answer everywhere ($\Delta -0.031$; exact-tie 0.746, $\Delta -0.033$); each interval lies wholly beyond
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie.ec: coverage 0.8928, population.f1.value 0.6849 (cells 0.684933); contrasts.prfba_exact_tie__vs__picrust2.ec: delta_f1_population -0.0937 (cells -0.093718), ci95 [-0.109, -0.0829], q_bh 0.002022
- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie_forced.ec: coverage 0.9995, population.f1.value 0.7455 (cells 0.745538, so 0.746 not 0.745); contrasts.prfba_exact_tie_forced__vs__picrust2.ec: delta -0.0331 (cells -0.033112), ci95 [-0.0392, -0.0291]
- Exact-tie forced CI [-0.039, -0.029] also lies beyond the 0.02 SESOI (sesoi_verdict exceeded), hence 'each interval'.

### 3.2. L76 (Abstract)

Current:

~~~latex
it widens to $-0.124$ $([-0.150, -0.099])$.
~~~

Replacement:

~~~latex
it widens to $-0.124$ $([-0.150, -0.099])$, and to $-0.126$ $([-0.152, -0.102])$ for the exact-tie rule.
~~~

Source and values:

- `region_validation/shared_reference/data_family/summary_family.json` arms.prfba_exact_tie.ec: coverage 0.8563, population.f1.value 0.6235 (cells 0.623497 -> 0.623); contrasts.prfba_exact_tie__vs__picrust2.ec delta -0.1259 (cells -0.125930), ci95 [-0.1522, -0.1017]

Note: Pre-registered nested-CV value -0.124 kept as is (but see section 4: cells give -0.123459).

### 3.3. L76 (Abstract)

Current:

~~~latex
registers as its bar ($-0.058$, $[-0.069, -0.051]$);
~~~

Replacement:

~~~latex
registers as its bar ($-0.058$, $[-0.069, -0.051]$; exact-tie $-0.060$, no interval computed);
~~~

Source and values:

- No contrast prfba_exact_tie__vs__base_nn1 exists in `region_validation/shared_reference/data/summary_genus.json`; point difference from `region_validation/shared_reference/data/cells.csv.gz` (space ec, density 1.0, unscoreable_min_align dropped; population F1 = f1 if outcome scored else 0): exact_tie - base_nn1 = -0.060135 (0.684933 - 0.745068).

### 3.4. L76 (Abstract)

Current:

~~~latex
(0.747 vs.\ 0.745) and no further.
~~~

Replacement:

~~~latex
(0.747 vs.\ 0.745; exact-tie 0.746) and no further.
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie_forced.ec: coverage 0.9995, population.f1.value 0.7455 (cells 0.745538, so 0.746 not 0.745); contrasts.prfba_exact_tie_forced__vs__picrust2.ec: delta -0.0331 (cells -0.033112), ci95 [-0.0392, -0.0291]
- base_nn1 0.745068; exact-tie forced - nn1 = +0.000470.

### 3.5. L77 (Abstract)

Current:

~~~latex
($+0.021$ population F1, though only $+0.018$ $([-0.004, 0.040])$ at family rank)
~~~

Replacement:

~~~latex
($+0.021$ population F1, though only $+0.018$ $([-0.004, 0.040])$ at family rank; $+0.019$ and $+0.015$ for the exact-tie rule against its own masked run, no intervals computed)
~~~

Source and values:

- base_random_abstain is the exact-tie forced run masked to prFBA's coverage. Genus: exact_tie - base_random_abstain = +0.018881 (0.684933 - 0.666052); family: +0.015277 (0.623497 - 0.608220). Computed from `region_validation/shared_reference/data/cells.csv.gz` (space ec, density 1.0, unscoreable_min_align dropped; population F1 = f1 if outcome scored else 0) and `region_validation/shared_reference/data_family/cells.csv.gz` (same filters); no contrast exists in the summaries.
- 'answering everywhere recovers two thirds of the deficit' holds for exact-tie too: +0.060605 of 0.093718 = 0.647.

### 3.6. L79 (Abstract)

Current:

~~~latex
($+0.107$ against $+0.030$ population EC F1 per tenfold more reference tips over a fourfold sweep)
~~~

Replacement:

~~~latex
($+0.107$, or $+0.106$ for the exact-tie rule, against $+0.030$ population EC F1 per tenfold more reference tips over a fourfold sweep)
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_nestedcv.ec.density_slope_f1_per_log10_tips 0.1068; arms.prfba_exact_tie.ec.density_slope_f1_per_log10_tips 0.106; arms.picrust2 0.0301. by_density f1 (0.25/0.5/1.0): exact-tie 0.6211/0.6527/0.6849.

### 3.7. L80 (Abstract)

Current:

~~~latex
does not appear in the controlled comparison (0.692 vs.\ 0.710), where full-length queries abstain more often (coverage 0.890 against 0.933).
~~~

Replacement:

~~~latex
does not appear in the controlled comparison (0.692 vs.\ 0.710; exact-tie 0.690 vs.\ 0.708), where full-length queries abstain more often (coverage 0.890 against 0.933 under either rule).
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie.ec.by_region: FullLength16S population.f1 0.6901, coverage 0.8898; V4 0.7076, 0.9329 (nested-CV: 0.6924/0.8898, 0.7097/0.9329).

### 3.8. L441 (PICRUSt2 as shipped)

Current:

~~~latex
(0.107 versus 0.030 per tenfold)
~~~

Replacement:

~~~latex
(0.107, or 0.106 for the exact-tie rule, versus 0.030 per tenfold)
~~~

Source and values:

- As item L79 (density_slope_f1_per_log10_tips). Ratio 0.106/0.0301 = 3.5, so 'about 3.5 times' holds.

### 3.9. L514 (Shared reference (controlled benchmark))

Current:

~~~latex
PICRUSt2 reaches population EC F1 0.779, prFBA 0.687 at coverage 0.893 --- a difference of $-0.092$ whose 95\% CI, $[-0.107, -0.081]$, lies entirely beyond the pre-registered 0.02 threshold.
~~~

Replacement:

~~~latex
PICRUSt2 reaches population EC F1 0.779, prFBA 0.687 at coverage 0.893 --- a difference of $-0.092$ whose 95\% CI, $[-0.107, -0.081]$, lies entirely beyond the pre-registered 0.02 threshold; the shipped exact-tie rule scores 0.685 at the same coverage, a difference of $-0.094$ $([-0.109, -0.083])$.
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie.ec: coverage 0.8928, population.f1.value 0.6849 (cells 0.684933); contrasts.prfba_exact_tie__vs__picrust2.ec: delta_f1_population -0.0937 (cells -0.093718), ci95 [-0.109, -0.0829], q_bh 0.002022

### 3.10. L515 (Shared reference (controlled benchmark))

Current:

~~~latex
narrows the deficit to $-0.031$ ($[-0.037, -0.028]$), still wholly beyond threshold.
~~~

Replacement:

~~~latex
narrows the deficit to $-0.031$ ($[-0.037, -0.028]$), still wholly beyond threshold (exact-tie: 0.746, $-0.033$, $[-0.039, -0.029]$).
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie_forced.ec: coverage 0.9995, population.f1.value 0.7455 (cells 0.745538, so 0.746 not 0.745); contrasts.prfba_exact_tie_forced__vs__picrust2.ec: delta -0.0331 (cells -0.033112), ci95 [-0.0392, -0.0291]

### 3.11. L516 (Shared reference (controlled benchmark))

Current:

~~~latex
(abstaining $-0.059$ to $-0.182$; forced $-0.029$ to $-0.034$; every BH-adjusted $q \leq 0.01$) and in both domains (abstaining, Bacteria $-0.095$ and Archaea $-0.052$; forced, $-0.032$ and $-0.019$)
~~~

Replacement:

~~~latex
(abstaining $-0.059$ to $-0.182$; forced $-0.029$ to $-0.034$; every BH-adjusted $q \leq 0.01$; exact-tie $-0.060$ to $-0.184$ and $-0.031$ to $-0.036$) and in both domains (abstaining, Bacteria $-0.095$ and Archaea $-0.052$; forced, $-0.032$ and $-0.019$; exact-tie $-0.097$ and $-0.053$, forced $-0.034$ and $-0.020$)
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` contrasts.prfba_exact_tie__vs__picrust2.ec.by_region delta_f1_population: -0.0598 (V4-V5_944R; cells -0.05985) to -0.1838 (V1-V2; cells -0.18377), all q_bh 0.002022; contrasts.prfba_exact_tie_forced__vs__picrust2.ec.by_region: -0.0307 (V7-V9) to -0.0364 (V3-V4).
- by_domain: exact-tie Bacteria -0.0966 (cells -0.096562), Archaea -0.0528 (-0.052791); forced Bacteria -0.0341 (-0.034055), Archaea -0.0195 (cells -0.019542, so -0.020).

### 3.12. L516 (Shared reference (controlled benchmark))

Current:

~~~latex
(AURC 0.812 versus 0.788;
~~~

Replacement:

~~~latex
(AURC 0.812 versus 0.788, exact-tie 0.787;
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie.ec.aurc 0.7868 (nested-CV 0.7881, picrust2 0.8125). PICRUSt2 leads exact-tie at every stored operating point too (f1_at_coverage 0.5/0.75/0.9/1.0: 0.8149/0.8065/0.7960/0.7787 vs 0.7911/0.7793/0.7660/0.7455).

### 3.13. L517 (Shared reference (controlled benchmark))

Current:

~~~latex
(0.813 and 0.788 cross-knob)
~~~

Replacement:

~~~latex
(0.813 and 0.788 cross-knob; exact-tie 0.786)
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie.ec.aurc_cross_knob 0.7862 (nested-CV 0.7875, picrust2 0.8133). Family rank exact-tie 0.7493 -> 0.7422 cross-knob, the same ~0.007 shift as the nested-CV arm.

### 3.14. L520 (Shared reference (controlled benchmark))

Current:

~~~latex
PICRUSt2 covers 0.92 of the interval, prFBA forced 0.70, and nn1 0.69.
~~~

Replacement:

~~~latex
PICRUSt2 covers 0.92 of the interval, prFBA forced 0.70 (exact-tie forced 0.69), and nn1 0.69.
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie_forced.ec.skill_score 0.6922 (nested-CV forced 0.7047); abstaining exact-tie 0.2702.

### 3.15. L521 (Shared reference (controlled benchmark))

Current:

~~~latex
(0.747 versus 0.745)
~~~

Replacement:

~~~latex
(0.747, or 0.746 for the exact-tie rule, versus 0.745)
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie_forced.ec: coverage 0.9995, population.f1.value 0.7455 (cells 0.745538, so 0.746 not 0.745); contrasts.prfba_exact_tie_forced__vs__picrust2.ec: delta -0.0331 (cells -0.033112), ci95 [-0.0392, -0.0291]

### 3.16. L523 (Shared reference (controlled benchmark))

Current:

~~~latex
($[-0.069, -0.051]$, $q = 0.002$), while PICRUSt2 clears it
~~~

Replacement:

~~~latex
($[-0.069, -0.051]$, $q = 0.002$; exact-tie $-0.060$, no interval computed), while PICRUSt2 clears it
~~~

Source and values:

- exact_tie - base_nn1 = -0.060135 from `region_validation/shared_reference/data/cells.csv.gz` (space ec, density 1.0, unscoreable_min_align dropped; population F1 = f1 if outcome scored else 0); no such contrast in the summary.

### 3.17. L530 (Shared reference (controlled benchmark))

Current:

~~~latex
answering everywhere gains $+0.061$ population F1 (0.687 to 0.747), two thirds of the 0.092 deficit
~~~

Replacement:

~~~latex
answering everywhere gains $+0.061$ population F1 (0.687 to 0.747; exact-tie 0.685 to 0.746, also $+0.061$), two thirds of the 0.092 deficit
~~~

Source and values:

- From `region_validation/shared_reference/data/cells.csv.gz` (space ec, density 1.0, unscoreable_min_align dropped; population F1 = f1 if outcome scored else 0): legacy(=nested-CV) forced - abstaining +0.060582; exact-tie +0.060605 (0.745538 - 0.684933), 65% of its 0.0937 deficit.

### 3.18. L531 (Shared reference (controlled benchmark))

Current:

~~~latex
its forced variant still reaches population F1 0.567 against the zero
~~~

Replacement:

~~~latex
its forced variant still reaches population F1 0.567 (0.567 for the exact-tie rule as well) against the zero
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie_forced.ec.by_novelty_bin['<0.865'].population.f1.value 0.5669 (nested-CV forced 0.5666; PICRUSt2 0.6338).

### 3.19. L535 (Shared reference (controlled benchmark))

Current:

~~~latex
(0.783 and 0.740)
~~~

Replacement:

~~~latex
(0.783 and 0.740; exact-tie 0.803 and 0.751)
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie.ec.by_novelty_bin['0.945-0.987'].population.precision 0.8035 and ['0.865-0.945'] 0.7511; from `region_validation/shared_reference/data/cells.csv.gz` (space ec, density 1.0, unscoreable_min_align dropped; population F1 = f1 if outcome scored else 0) (precision with abstentions 0): 0.803484 and 0.751083, so 0.803 (the 4-dp 0.8035 would double-round to 0.804). Still below PICRUSt2 (0.862143, 0.816910).

### 3.20. L539 (Shared reference (controlled benchmark))

Current:

~~~latex
a difference of $-0.124$ $([-0.150, -0.099])$; forced to answer, prFBA reaches 0.709 and still trails by $-0.040$ $([-0.048, -0.031])$.
~~~

Replacement:

~~~latex
a difference of $-0.124$ $([-0.150, -0.099])$ (exact-tie 0.623, $-0.126$, $[-0.152, -0.102]$); forced to answer, prFBA reaches 0.709 and still trails by $-0.040$ $([-0.048, -0.031])$ (exact-tie 0.707, $-0.042$, $[-0.050, -0.033]$).
~~~

Source and values:

- `region_validation/shared_reference/data_family/summary_family.json` arms.prfba_exact_tie.ec: coverage 0.8563, population.f1.value 0.6235 (cells 0.623497 -> 0.623); contrasts.prfba_exact_tie__vs__picrust2.ec delta -0.1259 (cells -0.125930), ci95 [-0.1522, -0.1017]
- `region_validation/shared_reference/data_family/summary_family.json` arms.prfba_exact_tie_forced.ec: population.f1.value 0.7072 (cells 0.707243); contrasts.prfba_exact_tie_forced__vs__picrust2.ec delta -0.0422 (cells -0.042184), ci95 [-0.0503, -0.0334]

### 3.21. L540 (Shared reference (controlled benchmark))

Current:

~~~latex
($-0.076$ to $-0.245$, all $q = 0.002$), the forced-arm equivalence with nn1 persists (0.7095 versus 0.7070), and the abstaining arm's shortfall against nn1 grows to $-0.081$ $([-0.107, -0.059])$.
~~~

Replacement:

~~~latex
($-0.076$ to $-0.245$, all $q = 0.002$; exact-tie $-0.077$ to $-0.247$), the forced-arm equivalence with nn1 persists (0.7095 versus 0.7070; exact-tie 0.7072), and the abstaining arm's shortfall against nn1 grows to $-0.081$ $([-0.107, -0.059])$ (exact-tie $-0.083$, no interval computed).
~~~

Source and values:

- `region_validation/shared_reference/data_family/summary_family.json` contrasts.prfba_exact_tie__vs__picrust2.ec.by_region: -0.0773 (V7-V9; cells -0.07734) to -0.2475 (V1-V2; cells -0.24747, so -0.247), all q_bh 0.00226.
- From `region_validation/shared_reference/data_family/cells.csv.gz` (same filters): exact_tie_forced 0.707243 vs base_nn1 0.706985; exact_tie - base_nn1 = -0.083488.

### 3.22. L541 (Shared reference (controlled benchmark))

Current:

~~~latex
becomes a point estimate of $-0.017$ with an interval spanning zero $([-0.053, 0.018])$, which drives its skill score negative ($-0.14$),
~~~

Replacement:

~~~latex
becomes a point estimate of $-0.017$ with an interval spanning zero $([-0.053, 0.018])$ (exact-tie $-0.020$), which drives its skill score negative ($-0.14$; exact-tie $-0.16$),
~~~

Source and values:

- From `region_validation/shared_reference/data_family/cells.csv.gz` (same filters): exact_tie - base_constant = -0.019799 (0.623497 - 0.643296); `region_validation/shared_reference/data_family/summary_family.json` arms.prfba_exact_tie.ec.skill_score -0.1606.

### 3.23. L541 (Shared reference (controlled benchmark))

Current:

~~~latex
($+0.016$ like-for-like against the legacy policy's own masked forced run).
~~~

Replacement:

~~~latex
($+0.016$ like-for-like against the legacy policy's own masked forced run; $+0.015$ for the exact-tie policy against its own).
~~~

Source and values:

- From `region_validation/shared_reference/data_family/cells.csv.gz` (same filters): exact_tie - base_random_abstain = +0.015277.

### 3.24. L543 (Shared reference (controlled benchmark))

Current:

~~~latex
by $+0.107$ $([0.044, 0.179])$, $q = 0.002$.
~~~

Replacement:

~~~latex
by $+0.107$ $([0.044, 0.179])$, $q = 0.002$ (exact-tie $+0.104$, no interval computed).
~~~

Source and values:

- From `region_validation/shared_reference/data_family/cells.csv.gz` (same filters): exact_tie - base_family_donor = +0.104175.

### 3.25. L587 (Shared reference (controlled benchmark))

Current:

~~~latex
Clade-max is an oracle bound, not an achievable method.}
~~~

Replacement:

~~~latex
Clade-max is an oracle bound, not an achievable method. The \emph{prFBA} rows use the policy leave-one-fold-out selection chose (legacy in every fold); the shipped exact-tie rows are shown alongside.}
~~~

Source and values:

- Caption sentence for the added rows.

### 3.26. L597 (Shared reference (controlled benchmark))

Current:

~~~latex
prFBA (forced to answer)     & 0.9995 & 0.747 & 0.769 & 0.70 \\
~~~

Replacement:

~~~latex
prFBA (forced to answer)     & 0.9995 & 0.747 & 0.769 & 0.70 \\
prFBA, shipped exact-tie     & 0.893 & 0.685 & 0.767 & 0.27 \\
exact-tie (forced to answer) & 0.9995 & 0.746 & 0.767 & 0.69 \\
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie.ec: coverage 0.8928, population.f1 0.6849, conditional_matched.f1 0.7672 (n 57,021; cells 0.767215), skill_score 0.2702; arms.prfba_exact_tie_forced.ec: 0.9995, 0.7455 (cells 0.745538), 0.7672, 0.6922.

Note: Two rows inserted after the existing forced row; the existing rows are untouched.

### 3.27. L635 (Retention depth)

Current:

~~~latex
recovers $+0.061$ population F1 of its $0.092$ deficit, so two thirds of that deficit does sit in the refusals, but the residual $-0.031$ ($[-0.037, -0.028]$) persists at full coverage.
~~~

Replacement:

~~~latex
recovers $+0.061$ population F1 of its $0.092$ deficit ($+0.061$ of $0.094$ for the exact-tie rule), so two thirds of that deficit does sit in the refusals, but the residual $-0.031$ ($[-0.037, -0.028]$; exact-tie $-0.033$, $[-0.039, -0.029]$) persists at full coverage.
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie.ec: coverage 0.8928, population.f1.value 0.6849 (cells 0.684933); contrasts.prfba_exact_tie__vs__picrust2.ec: delta_f1_population -0.0937 (cells -0.093718), ci95 [-0.109, -0.0829], q_bh 0.002022
- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie_forced.ec: coverage 0.9995, population.f1.value 0.7455 (cells 0.745538, so 0.746 not 0.745); contrasts.prfba_exact_tie_forced__vs__picrust2.ec: delta -0.0331 (cells -0.033112), ci95 [-0.0392, -0.0291]

### 3.28. L1085 (Discussion)

Current:

~~~latex
to $-0.027$ on matched cells (0.769 against 0.796)
~~~

Replacement:

~~~latex
to $-0.027$ on matched cells (0.769 against 0.796; $-0.029$, 0.767, for the exact-tie rule)
~~~

Source and values:

- From `region_validation/shared_reference/data/cells.csv.gz` (space ec, density 1.0, unscoreable_min_align dropped; population F1 = f1 if outcome scored else 0), cells all three of exact_tie, legacy and picrust2 answer (57,021): exact_tie 0.767215, legacy 0.769167, picrust2 0.795978; exact_tie - picrust2 = -0.028763.

### 3.29. L1085 (Discussion)

Current:

~~~latex
(0.713 against 0.729; $-0.129$ at the earlier record depth of 20) to $-0.092$;
~~~

Replacement:

~~~latex
(0.713 against 0.729; $-0.129$ at the earlier record depth of 20) to $-0.092$ ($-0.094$ for the exact-tie rule, the same rule at both ends);
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` arms.prfba_exact_tie.ec: coverage 0.8928, population.f1.value 0.6849 (cells 0.684933); contrasts.prfba_exact_tie__vs__picrust2.ec: delta_f1_population -0.0937 (cells -0.093718), ci95 [-0.109, -0.0829], q_bh 0.002022

Note: With the exact-tie arm both ends of this comparison use the union.

### 3.30. L1087 (Discussion)

Current:

~~~latex
(0.107 against 0.030 population F1 per $\log_{10}$ reference genomes)
~~~

Replacement:

~~~latex
(0.107, or 0.106 for the exact-tie rule, against 0.030 population F1 per $\log_{10}$ reference genomes)
~~~

Source and values:

- As item L79.

### 3.31. L1088 (Discussion)

Current:

~~~latex
being worth $+0.021$ population F1 over declining as many at random;
~~~

Replacement:

~~~latex
being worth $+0.021$ population F1 over declining as many at random ($+0.019$ for the exact-tie rule against its own masked run);
~~~

Source and values:

- As item L77 (+0.018881).

### 3.32. L1112 (Limitations)

Current:

~~~latex
(matched conditional EC F1 0.769 against 0.796, from a matched 0.802 against 0.727 as shipped), while PICRUSt2's population lead persists and widens ($\Delta$ from $-0.016$ to $-0.092$).
~~~

Replacement:

~~~latex
(matched conditional EC F1 0.769, or 0.767 for the exact-tie rule, against 0.796, from a matched 0.802 against 0.727 as shipped), while PICRUSt2's population lead persists and widens ($\Delta$ from $-0.016$ to $-0.092$, or to $-0.094$ for the exact-tie rule).
~~~

Source and values:

- As items L1085 (0.767215) and L75 (-0.093718).

### 3.33. L1113 (Limitations)

Current:

~~~latex
and quartering the shared reference costs prFBA 0.064 population F1 against PICRUSt2's 0.018
~~~

Replacement:

~~~latex
and quartering the shared reference costs prFBA 0.064 population F1 under either rule against PICRUSt2's 0.018
~~~

Source and values:

- `region_validation/shared_reference/data/summary_genus.json` by_density f1 at 1.0 minus 0.25: nested-CV 0.6867-0.6224 = 0.0643; exact-tie 0.6849-0.6211 = 0.0638; picrust2 0.7787-0.7606 = 0.0181.

## 4. Unresolved points and caveats

1. **No like-for-like depth-20 union baseline for tab:main's scorer.**
   - Affected items: 1.37, 1.38, 1.39 and 1.108 (L631, L632, L633, L1102).
   - These sentences compare tab:main's figures at depth 20 and depth 100 under the legacy reducer (conditional genus within 0.01; gene capture +0.7 to +1.3 points; population genus 0.375 → 0.472; the ~0.49 projection).
   - The original panel's only depth-20 tab:main scoring is legacy (`region_validation/data/taxacc_region_summary.json`, `select_mode` legacy). The June exact-tie variants (`taxacc_region_summary_et.json`, `gene_capture_by_region_et.json`) predate the corrected scorer: no `select_mode`, no `correct_population`, no exclusion of unscoreable ranks.
   - The edits therefore delete these clauses rather than replace them. What remains is union at both depths (coverage 0.766 → 0.905; conditional EC F1 0.798 → 0.802; population EC F1 0.600 → 0.713).
   - The union depth-100 pooled population genus accuracy is 0.4379 (22,783 of 52,026), if a level without a comparison is wanted.
   - Producing the baseline: re-score the depth-20 union selection with the current `score_taxacc.py` / `score_genecap.py` into a scratch tree, with a patched `_common.py` so nothing is written into `region_validation/`.
     - Option (a): re-score `region_validation/data/selection_exacttie.json`, after confirming that its reps match `region_validation/picrust2_comparison/data/prfba_selection.json`.
     - Option (b): first re-select the depth-20 records with `--select-mode exact-tie`, then re-score.
     - Not done here, because the scorer writes into its tree.
2. **The nested-CV family headline is double-rounded (pre-registered, so not changed).**
   - `summary_family.json` stores `contrasts.prfba_nestedcv__vs__picrust2.ec.delta_f1_population` = −0.1235.
   - The cells give −0.123459 (0.625968 − 0.749427), which rounds to −0.123.
   - The manuscript prints −0.124 at L76 and L539. The authors should decide whether to correct it.
3. **The AF-gate share (item 1.57, 26.7%) rests on MAG cells that are only in the scratchpad.**
   - The 2%-AF-gate cells are in `verify/agent_scratch2/sk3/afgate2.pkl`, not in the repo.
   - The same pickle reproduces the current 31.8% at the legacy comparator, so the logic is the cited one, but the provenance is weak.
4. **Half-way confidence bounds.**
   - In tab:threepanel and item 1.102, some stored 4-dp bounds sit exactly half-way: −0.1905, +0.0125 (twice), −0.0385, −0.2905, −0.1145.
   - They were rounded away from zero, the convention the current table uses.
   - Their third decimal cannot be recovered from the 4-dp JSON, and it is Monte Carlo noise in any case: the phylum iteration order depends on PYTHONHASHSEED.
5. **tab:main full-length include-self gene capture is 99.979%.** At the table's one decimal it prints as 100.0\%. Use 99.98\% in that one cell if a ceiling-looking value is unwanted.
6. **The PICRUSt2 population-gap interval is not a repo output.**
   - Item 1.36 deletes the legacy comparison "(under the legacy reducer 0.726, a gap within resampling noise)".
   - The only interval for the union gap, −0.016 [−0.021, −0.011], is the genome-clustered bootstrap in D100_VS_D20.md (scratchpad `p2cells.json`). It is not a repo output, so it is not proposed as text.
7. **Not re-derived.**
   - The deployment-run shares at L210 and L289 (62.5%, 2,344 of 3,753). They come from a separate depth-20 run whose selection rule I could not find recorded in a repo file, and they are not a legacy-tree result.
   - The backend-concordance result at L337–338 ("identical representative genome set for 99.94--99.97% of ASVs"). It is a development-era selection result whose rule is not recorded in the manuscript or in a file I found. It is not a tab:main or cross-panel result, so it is not a change item.
   - The real-MAG numbers (tab:mag), which are rule-independent.
8. **Outside the manuscript.**
   - The voiced variant takes every edit unchanged (verified).
   - `presentations/make_deck.py` already defaults to `--rule union` (line 3659), so a deck rebuilt now matches the union text.
   - `region_validation_d100/D100_VS_D20.md` keeps its L100 framing.
9. **Figures.** The six figure re-renders in section 2 have not been run, because only this report may be written. Until they are, captions edited by items 1.32–1.35 will describe union numbers that the legacy images do not show. fig:magcomp and fig:magcross have no number-bearing caption text, but their images carry legacy values.

## Reproduction

All scripts are under:

`/tmp/claude-1000/-home-freiburger-Documents-prFBA/0c6f5bc6-8224-5a34-9b20-7d5643487901/scratchpad/union_chg/`

They are read-only on both repositories.

| Script | What it computes |
|---|---|
| `edits_s1.py`, `edits_s3.py` | the edit data this file was generated from (old, new, sources) |
| `verify.py` | uniqueness, line and sequential-application checks; writes the edited copies `prFBA_union_applied*.tex` |
| `render.py` | this report |
| `build/` | baseline and edited pdflatex builds |
| `tabmain.py`, `abst.py`, `abst2.py`, `tier.py` | tab:main; abstention, coverage and single-genome shares; tier capture |
| `mag.py`, `mag2.py` (output `mag_out.json`) | tab:magvs16s, tab:magval, ranges, τ=99, easier-by, identity bins |
| `tp3.py`, `tp3b.py`, `welch.py`, `empt8.py` | three-panel and two-panel values, Welch tests, pooled emptied and survivor rates |
| `sr*.py`, `srcells*.py`, `srcond.py` | controlled-benchmark arms, contrasts and the cell-level recomputations |
