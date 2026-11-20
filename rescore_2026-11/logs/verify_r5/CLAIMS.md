# Independent verification of the R5 manuscript rewrite: claims, disclosures, build, style

Verifier: independent agent, 2026-11-16. Scope: `~/Documents/prFBA_paper/prFBA.tex` at **2387ac2**
(diff `0214a14..2387ac2`, commits 75d6f03 and 2387ac2). Report only: no tracked file in any repo was edited.
Line numbers (L) are those of prFBA.tex at 2387ac2 (identical to the working tree).
Scratch scripts and renders: `/tmp/claude-1000/-home-freiburger-Documents-prFBA/400ec395-f172-59da-a98d-3dbe6b3a4256/scratchpad/verify_r5/claims/`
(`gap_ci.py`, `as15.py`, `nn1_wilcoxon.py`, `cnt130.py`, `style.py`, `build/`, `f13page-44.png`, `f13fig-44.png`).
RUN5B_REPORT.md was used only to locate items.

## Verdict

**No blocking issues.** Every one of the 24 meaning changes is stated correctly against the source artefacts,
and none is overstated. No sentence still asserts an old conclusion. The paper builds with 0 undefined
references or citations and only the two pre-existing overfull boxes. The disclosures match PLAN.md and the
run reports, and the provenance exception list is complete as scoped.

There are about a dozen minor or cosmetic items, listed below. The ones most worth fixing before submission:
1. L741: "still below PICRUSt2's 0.729" is not tied to a selection rule, in a sentence that also cites the
   legacy reducer. Under the legacy reducer the depth-100 figure is 0.734, above PICRUSt2.
2. L643: "That matches the EC result against nn1." The EC result against nn1 in the controlled comparison is a
   significant population deficit (−0.058). The phenotype null matches only the conditional tie.
3. L659: "rescored the cells table and reproduced it byte for byte". In the re-run that this text now
   describes, score_folds **reused** the cells table by fingerprint and did not rescore it
   (`SR/logs/score_genus_c13.log` line 15: "reusing existing … cells.csv.gz: fingerprints match").
4. Balance. The paper reports the nn1 all-genome sensitivity (+0.06, p 0.46) but not the PICRUSt2 one, which is
   now +0.284 media against prFBA, p = 0.057 (was 0.221, p 0.11).
5. F13 panel B: the orange annotation now overprints the 0.98/0.99/0.97 bar labels.

**Bootstrap CI (legacy reducer, as-shipped).** prFBA legacy minus PICRUSt2 population EC F1 is **+0.0054**, with
a 95% genome-clustered CI of **[+0.0006, +0.0101]**. That is 10,000 replicates, seed 20261116, resampling the
9,708 source genomes over 67,713 amplicons, using exactly make_deck.py's `gap_ci_d100()`. The interval excludes
zero, but narrowly: 1.26% of replicates are ≤ 0. For the exact-tie union the gap is −0.0062, CI [−0.0110, −0.0015].
Both match the deck's GAP_CI_D100. The paper states the legacy lead as a point estimate only, so nothing is
overstated.

---------------------------------------------------------------------------------------------------------------

## 1. Meaning changes against source

| claim (tex) | source and value | verdict |
|---|---|---|
| ATP vs nn1 no longer significant: −0.030 media, Wilcoxon p 0.41, q 0.51; all 1,000 genomes +0.06, p 0.46 (L642, L1196) | `SR/data/phenotype_atp_c13/phenotype_summary.json` contrasts.prfba_exact_tie__vs__base_nn1.controlled: Δ −0.030303, p 0.4075, q 0.5094, n 759; all-genome Δ +0.06, p 0.4561, n 1000 | **correct**; the sign reading ("favors nn1") is right |
| McNemar favours PICRUSt2 on four of 22 testable media: Ac.NO, Ac.O2, Etho.O2, For.O2, q 0.029–0.046 (L79, L636, L646, L1195) | same summary, mcnemar_per_medium: 22 media with p (BH over 22); q 0.0289 / 0.0396 / 0.0461 / 0.0461, c > b in all four, none with b > c | **correct** |
| ATP vs PICRUSt2: 0.07 more media, p 0.79, n 759 (750 V4), ΔEC-F1 −0.031, ρ −0.19, p 1.5e-7 | summary Δ 0.068511, p 0.7909; derived_numbers.post.json 750, −0.03133, −0.18938, 1.47e-7 | **correct**. See 2.d for the unreported all-genome sensitivity |
| Depth probe recovered answers at or above the answered set in EC space: F1 0.829 vs 0.798, precision 0.791 vs 0.790, recall 0.904 vs 0.834 (L67, L736, L739, L1214) | `region_validation/picrust2_comparison/data/depth_probe/depth_probe_results.json` 0.8288 / 0.7914 / 0.904; d20 `comparison_summary.json` prfba_exclude conditional 0.798 / 0.7898 / 0.8344 | **correct**. Minor wording at L67 and L1214, see 2 |
| Surviving-uncovered rates 2.4/2.3/1.8% FL; 1.8/1.8/1.3% pooled; strata 0–3.0%, original 50–99 = 0 of 274; 42–55-fold; spreads 0.6 (0.5); CMH p 0.117 and 0.99 (L1068–1069) | `three_panel/emptied_strata_d100.json`: survivors non_ok_rate 0.0239/0.0227/0.0181; 50–99 n = 274, rate 0.0; CMH p 0.11667 (vs random), 0.99475 (vs k-means); folds 41.8/55.2 | **correct**. Direction sentence, see 2.e |
| Coverage arithmetic "to within 2.3 points" (L1070–1071) | 0.667−0.655, 0.771−0.754, 0.954−0.931 = 1.2, 1.7, 2.3 points | **correct**. The old "within two points" was in fact 3.2 points, so this also fixes a pre-existing error |
| As-shipped population deficit −0.006 (0.723 vs 0.729); depth 20 −0.120 (L527, L1197, L1226) | d100 `comparison_summary.json` by_region.ec.all.all population f1: prfba_exclude 0.7228, picrust2 0.7290; d20 0.6088 | **correct** |
| EC F1 vs nn1 0.782 vs 0.780, Wilcoxon signed-rank p 0.84 (L1196) | reproduced independently (`nn1_wilcoxon.py`): 759 pairs, 0.7820 vs 0.7802, scipy default (zero_method "wilcox") p = 0.836 | **correct**. Caveat in section 5 |
| A-C13 counts as a re-run on a re-drawn sample: 43,174 ATP rows, 0 of 8,660 gated, 2,372 models all grow, 107 reversals (L80, L680, L687) | RUN4_LOCAL §2; derived_numbers.post.json pheno_atp_rows 43174, pheno_gf_rows 8660, gated 0, models 2372, growing 2372, reversal 107 | **correct**. The 43,140 / 8,614 / 114 figures are now labelled "on the sample of the time" (L674–679) |
| Cache-repair statement (L516–519): 278 amplicons from 38 genomes; 953 / 161 before; 675 of 953; 818 of 1,060; 167 truncated; 675 more scoreable (67,713); 243 pred_empty remain; coverage 0.889 → 0.902; population 0.713 → 0.723; conditional 0.802; 10k gene capture +0.7–1.1 points population, ≤ 0.03 points conditional | comparison_cells outcome counts (truth_empty 278); gene_capture_by_region.json (pred_empty 243); pre-repair vs new headline: population +0.672 to +1.145 points, conditional |Δ| ≤ 0.029 points | **correct** |
| Denominators L397 / L419: 9.8% = 9.5% + 0.4%; 9.5% → 6.0%; 0.4% stays 0.4% | number_diff: 0.09844, 0.09486, 0.00359; 0.995: 0.05975, 0.00366 | **correct**. Cosmetic: the printed parts add up to 9.9% (rounding) |
| "about three fifths (58.6–60.3%)" (L1185) | tab:main population capture for V3–V4, V4, V4–V5, V4–V5_944R, V7–V9: 59.0, 59.3, 58.6, 58.9, 60.3 | **correct**. "less than three fifths" correctly became "about" |
| "within 0.017 at every MAG completeness"; "at most 0.030 (V7–V9)" (L1028) | number_diff: max gap 0.016718; V7–V9 0.029659 | **correct** |
| Legacy reducer as shipped 0.734 vs 0.729 (+0.005) (L528, L1197, L1226) | comparison_summary prfba_exclude_legacy 0.7343 vs picrust2 0.7290 | **correct**. See section 3 |

---------------------------------------------------------------------------------------------------------------

## 2. Search by phrasing for surviving old conclusions

I searched the whole tex for: favor(s) prFBA, phenotype advantage, outperform, two media, two of (the) 22,
below the answered set, at or below, resampling noise, (slightly) conservative, reproduce\*, significant\*,
closer to truth, advantage, 43,140. I also searched for the old values: 755, 747, 0.13 more, p = 0.49, 0.146,
67,038, 19,686, 2,689, 0.713, 0.600, 11.1%, 9.6%, 1.6%, 3.4%, 29- to, p = 0.004, 107/76/62, "two points",
"few percent". The paper has no conclusions section and no supplement; the abstract, introduction, results,
discussion, limitations, captions and table notes were all covered.

**No survivor asserts an old conclusion.** Every remaining hit is either a valid current statement or an
explicitly historical one: 43,140 at L80, L674 and L1203 is labelled "pre-registered pass" or "sample of the
time"; 8,614 and 49 at L661 and L675 are pre-registered or A-C12; 0.713 at L519 is the "from" value.
Residual wording issues:

| L | text | problem | suggested fix |
|---|---|---|---|
| 741 | "on the as-shipped PICRUSt2 axis prFBA's population EC F1 rises from 0.609 to 0.723, still below PICRUSt2's 0.729" | The as-shipped mean is stated without a selection rule, in a sentence whose previous clause is "compared under the legacy reducer". Under legacy the depth-100 value is 0.734 (d20 0.616), above PICRUSt2. | "… under the shipped exact-tie union prFBA's population EC F1 rises from 0.609 to 0.723, still below PICRUSt2's 0.729 (the legacy reducer's 0.734 is above it)", or at least insert "under the shipped exact-tie union" |
| 643 | "That matches the EC result against nn1." | The EC result against nn1 in the controlled comparison is a significant population deficit (−0.058, L613). Only the conditional EC tie matches. L644 explains it, but the sentence on its own over-reads. | "That matches the EC result against nn1 on the cells prFBA answers." |
| 659 | "The scorer found the summary bound to this run, rescored the cells table and reproduced it byte for byte, …" | In the re-score pass 2, which produced the "final gapfilled summary" this sentence now refers to, score_folds reused cells.csv.gz by fingerprint and did not rescore it (score_genus_c13.log L15; RUN4_LOCAL §3). | "… found the summary bound to this run, reused the cells table by its fingerprint (byte-unchanged), left every arm and contrast unchanged, …" |
| 642 / 1196 vs 641 | nn1 all-genome sensitivity reported (+0.06, p 0.46); PICRUSt2 all-genome sensitivity not | Asymmetric reporting. PICRUSt2's all-genome sensitivity is now prFBA +0.284 media, Wilcoxon p = 0.057 (was 0.221, p 0.11), borderline and in the direction of the EC deficit. | Add to L641: "(over all 1,000 paired genomes +0.28 media, p = 0.057)", or drop the nn1 all-genome clause too |
| 1069 | "and in the direction opposite to the coverage gap: the original panel, the best covered, leaves the largest share …" | The ordering is crude only. CMH is now null for both pairs, and stratum standardization reorders original and k-means: `standardised_non_ok_rate` 1.97% / 2.25% / 1.78%. The per-stratum Fisher tests also point in opposite directions (original higher in 0–19 vs random, p 0.005; lower in 50–99 vs k-means, p 0.001). The argument no longer needs the direction claim. | "These rates differ by about half a point (0.6, and 0.5 pooled), crude ordering original > k-means > random but reversed between original and k-means after stratum standardization; stratified by the number of deleted hits, neither difference is significant (CMH …)." |
| 1214 | "but score above it in EC space (EC F1 0.829 against 0.798)" | L736 and L739 say "at or above" (precision 0.791 vs 0.790, level); n = 283, no test. "above" is slightly stronger. | "but score at or above it in EC space" |
| 67 | "in EC space those answers are no worse … (…; genus accuracy 0.60)" | Genus is given without its comparator, and it is the one dimension where the recovered answers are worse (0.601 vs 0.651, L739). This is not wrong, but it reads one-sided. | "(…; genus accuracy 0.60 against 0.65)" |
| 80 | "both phenotype passes were re-run and the 1,000-genome ATP pass reproduced every ATP vector of the pre-registered pass (43,140 of 43,140 rows); in the A-C13 configuration, re-run …" | Implied but not stated: the reproduction was on the original sample. The re-run's ATP pass (43,174 rows, new sample) reproduces nothing. | "… re-run and, on the original sample, the 1,000-genome ATP pass reproduced …" |
| 529 | "By region the mean reverses only in the five regions …" (directly after the legacy sentence L528) | Ambiguous rule: the numbers 0.777/0.768/0.774 are exact-tie. The legacy pattern is the same (legacy also trails in the same five regions: −0.016, −0.018, −0.013, −0.016, −0.006), so the claim holds either way. | Move L528 after L529, or "By region, under either rule, …; under the exact-tie union prFBA's … (0.777, …)" |
| 1197, 1226 | "a +0.005 lead as shipped" | Legacy is not the shipped rule; "as shipped" here means "in the as-shipped comparison". | "a +0.005 lead in the as-shipped comparison" |
| 528 | "Under the legacy reducer, the policy the controlled comparison below selects, the mean does not reverse (population F1 0.734)." | Accurate (see section 3). Optional calibration. | add "(+0.005, 95% CI [+0.001, +0.010] over source genomes)" and, if wanted, "selected there by margins of 0.001–0.003" |

---------------------------------------------------------------------------------------------------------------

## 3. The legacy-reducer claim (L190, L527–528, L1197, L1226)

- **Values.** `region_validation_d100/picrust2_comparison/data/comparison_summary.json`,
  `by_region.ec.all.all.prfba_exclude_legacy.population.f1.value` = **0.7343** (n 67,713, coverage 0.9029);
  `picrust2` = **0.7290**; difference **+0.0053**, which prints as +0.005. Correct.
- **Bootstrap.** Reproduced make_deck.py `gap_ci_d100()` (L136–163, `GAP_CI_SEED, GAP_CI_REPS = 20261116, 10000`).
  The unit is the **source genome**, not the amplicon: a cluster bootstrap over 9,708 genomes, with per-genome
  sums of per-amplicon differences, abstention scored 0 and truth_empty dropped.
  - Legacy: +0.0054, **[+0.0006, +0.0101]**, excludes 0 (1.26% of replicates ≤ 0).
  - Union: −0.0062, [−0.0110, −0.0015].
  - For reference, an unclustered amplicon bootstrap is legacy [+0.0032, +0.0075] and union [−0.0083, −0.0040].
  - Before the repair the same code gave legacy (−0.007, +0.003) (make_deck.py comment), so the legacy lead is new
    with the repair.
- **"the policy the controlled comparison below selects".** This is consistent with how the paper describes the
  primary arm:
  - L77: "selection policy chosen by leave-one-fold-out selection, which picked legacy in every fold";
  - L595: "All five folds chose the legacy policy … so the prFBA arm below is the legacy policy";
  - the tab:sharedref caption: "the policy leave-one-fold-out selection chose (legacy in every fold)".

  The phrase is accurate. The selection margins are 0.001–0.003, well under the SESOI; I recomputed them from
  `summary_genus.json` / `summary_family.json` nested_cv: legacy beats the runner-up by 0.0014–0.0018 at genus
  and 0.0017–0.0025 at family. A reader could take "selects" as a stronger preference than the data give, which
  is why the margins are offered as an optional addition in section 2.
- **Unqualified "as-shipped mean reverses" elsewhere.**
  - L527 qualifies it ("under the shipped exact-tie union").
  - L1197 and L1226 qualify it and add the legacy figure.
  - L76, L525 and the sharedref section title concern the controlled design or the conditional margin.
  - The only rule-unqualified "below PICRUSt2" population statement is **L741** (section 2, first row).
- **Section list (L190).** It lists picrust2, discussion and limits, matching all three places the legacy as-shipped
  figure appears. `\label{sec:discussion}` was added (L1176) and resolves.

---------------------------------------------------------------------------------------------------------------

## 4. Disclosures and the provenance statement

| disclosure (L) | checked against | verdict |
|---|---|---|
| Path A (L591–592): eligibility decided 2026-10-13; 123 of 161 empty-EC exclusions were stored fetch failures; 122 pass filters, 65 of those archaea; ~115 more genomes; 17 eligible genomes changed EC sets; 23 of 30 genus-rank fold references; 130 of 63,871 test cells | `panel_eligible.json` meta generated 2026-10-14T03:45Z (2026-10-13 local); drop_counts empty_ec 161. `plan/pipelines_shared_ref.md` L39–42, L129: 122 pass, 65 archaea / 57 bacteria, ~115. RUN3 R3.1: 17 changed_eligible; R3.2: genus 23 of 30 refs. **130 recounted**: cells.csv.gz density 1.0, ec, picrust2 arm, source in the 17 changed genomes = 130 cells from 17 genomes, all 130 scored (the denominator also checks out: 63,890 rows = 63,871 scored + 19 unscoreable_min_align); this agrees with W8 source_changed 780 = 130 × 3 densities × 2 spaces | **correct**. The paper's "65 of those archaea" (of the 122) is more precise than PLAN's draft "(65 of them archaea)" of the 123 |
| D6 (L1019): 1733.9197 and 1733.9535 at exactly 25,000 contigs; "may be truncated"; 10 of 49,990 MAG queries; not re-fetched; two BV-BRC-subset (sludge) genomes capped, no real MAG hits them | `plan/pipelines_mag_deploy.md` L73–79, L284–288, L338–340 (0 hits in `mag_vs_bvbrc_sludge_skani_raw.tsv`); `two_panel/RESULTS.md` §11; 2 genomes × 5 tiers = 10 | **correct**. "may be" is right: D6's confirming count query was not run |
| D7 (L385): "on the gene-family cache of the time" | PLAN D7 | **correct** |
| Phenotype re-draw (L638–639): stratified by domain × template × ΔF1 quintile; 200 of 7,855 candidates changed; 13 of 15 strata; 442 of 1,000 and 62 of 200 kept; 106 of 19,072 ATP rows and 0 of 2,684 gapfilled rows differ | RUN4_LOCAL §2 (pheno_sample_diff.json, pheno_attrib_{atp,gf}.json) | **correct** |
| Provenance (L518): "Every result here that reads the gene-family cache was recomputed …; the exceptions … are the backend concordance check, the shared-reference panel's eligibility, and the phenotype passes and side studies run on the earlier phenotype sample" | PLAN §1 (D1, D7, D9), §6, number_diff manual_resolution, RUN5B "kept" items | **complete as scoped**; notes follow |

Notes on the completeness of the exception list:
- **The three listed exceptions are each labelled where they appear.**
  - L385: "cache of the time".
  - L591: the panel eligibility.
  - L674–679, L685–686: "on the sample of the time", "in the earlier sample", "of the earlier gapfilled sample".
- **Items the pipeline did not recompute but that I recomputed on the repaired artefacts. Both hold at print, so the
  statement stays true.**
  - **AS-15 (L530)**: rxn-vs-EC divergences, recomputed from `comparison_summary.json` (`as15.py`):
    - archaeal max |Δ| 0.0262 (V4–V5_944R, PICRUSt2);
    - PICRUSt2 distant bins 0.0137 (ident < 0.865) and 0.0112 (NSTI ≥ 0.5);
    - NSTI-unknown 9-amplicon cell +0.0489 (prFBA);
    - ident-unknown 8-amplicon cell −0.0266;
    - overall prFBA 0.8048 and PICRUSt2 0.7216, with |Δ| < 0.01 in every region and in both source strata.

    All match "0.026 / 0.014 / 0.05 / 0.03" and "0.805 / 0.722".
  - **SR-02 (L340, L595)**: the inner-CV margins 0.0014–0.0025 support "0.001–0.003". The winner is legacy in 5 of
    5 folds at both ranks.

  Recommend adding both as W6 recipes so the repo, not this report, holds the proof.
- **Not recomputed, expected unchanged, not cache-reading or effectively current.** None is a provenance exception.
  - The 0.651 / 0.800 depth-20 anchor accuracy comparators (L739) are taxonomy-only; the union never reads the gene
    provider.
  - The wobj "weighted term median 0.11%" (L310): 358 of 360 wobj rows are bit-identical (RUN4 §1).
  - Deployments: built after the repair (PLAN D4, §6).
  - GTDB arm, real MAGs, alignments and 16S front-end: no gene data.
  - `no_function_set` eligibility drops (238): the genome_function_sets key set is unchanged (RUN3 R3.1: 9,762 keys,
    set unchanged), so the D1 disclosure needs no second class.
- **Optional addition to the Path A disclosure.** PLAN's pipelines_shared_ref.md L348 notes that the kept panel
  under-represents those 65 archaea. The eligible panel has 1,113 archaea, so this is about 6%. A clause such as
  "so the archaeal arm is about 6% smaller than a re-plan would give" would state the consequence, not just the
  count.
- **Cosmetic.** In "23 of the 30 genus-rank fold references", the reader has not been told there are 30
  (5 folds × 2 domains × 3 densities). Suggest "(five folds × two domain trees × three densities)".

---------------------------------------------------------------------------------------------------------------

## 5. Wilcoxon recipe for 0.782 vs 0.780, p 0.84

**It is documented in the prFBA repo.**
- **Script.** `manuscript_numbers/derived_numbers.py`, `group_phenotype()`, keys `pheno_nn1_n_paired`,
  `pheno_nn1_ecf1_prfba`, `pheno_nn1_ecf1_nn1` and `pheno_nn1_ecf1_p`. The recipe text is embedded:
  - cells.csv.gz density 1.0, space ec, f1 of prfba_exact_tie and base_nn1 at each paired genome's (md5, fold);
  - over the prfba_exact_tie__vs__base_nn1 ATP-informative paired set;
  - scipy `wilcoxon`, default zero_method.

  Committed in **9148c20**.
- **Output.** `manuscript_numbers/derived_numbers.post.json` (tracked) holds 0.78197 / 0.78022 / p 0.83621 with
  n 759, plus sha256 pins of every input:
  - cells.csv.gz `d57709e3…`;
  - phenotype_atp_c13/per_genome.jsonl `519fd193…`, which equals the file on disk;
  - phenotype_summary.json `bf19ff02…`.
- **Pre-repair note.** `KNOWN_PRE_MISMATCH` records that the old 0.781 / p 0.64 had no surviving recipe.
- **Inputs on disk.** cells.csv.gz (56 MB) and per_genome.jsonl (92 MB) are git-ignored, as are all large
  cell tables in this repo. The summary is tracked. The pins make the recipe checkable.
- **Independent reproduction** (`nn1_wilcoxon.py`): 759 pairs, 0.7820 vs 0.7802, p = 0.836.
- **Caveat.** The p depends on zero handling, because many paired differences are exactly 0: `wilcox` (scipy
  default, the variant phenotype_fold.py's `wilcox()` uses at L1233) gives 0.84, `pratt` 0.29 and `zsplit` 0.36.
  All three are > 0.05, so "the two tie on ECs" holds.

Suggestions:
- name the variant in the text: "Wilcoxon signed-rank, zero differences dropped, p = 0.84";
- add one line to a results document for discoverability, since today only derived_numbers carries it. The
  natural home is `region_validation/shared_reference/C12_SCOPE.md` or the shared-reference RESULTS doc that
  carries the other A-C13 phenotype figures.

---------------------------------------------------------------------------------------------------------------

## 6. Build at 2387ac2

- `git archive 2387ac2 | tar -x` into `…/verify_r5/claims/build/` (images are tracked), then
  `pdflatex → bibtex → pdflatex → pdflatex`.
- **Undefined references: 0. Undefined citations: 0.** There is no "There were undefined references" line and no
  `??` in the PDF text. bibtex reports 0 warnings.
- **Overfull boxes: exactly the two pre-existing ones.** `\hbox` 23.69 pt at L133–138, and `\hbox` 16.96 pt at
  L554–571 (tab:picrust2). Other warnings are only `!h`→`!ht` float notes and one float-only page (p. 38).
  98 pages.
- **F13** is on PDF page 44; renders are `f13page-44.png` and `f13fig-44.png` in the scratch folder.
  - Title n = 67,713. Panel A: 0.82 → 0.61 against 0.73 → 0.23. Panel B coverages: 0.98, 0.99, 0.97 and 0.84 for
    the ≥ 0.995 bin. Panel C: per-region bars, prFBA ahead in all eight.
  - All of these agree with the caption and the text (L723: 0.993 / 0.841).
  - The 0.98\linewidth change clears no warning but is harmless.
  - **Cosmetic.** In panel B the orange annotation ("PICRUSt2 always answers (1.00); in the floor bin its
    precision is 0.32, where prFBA answers nothing") overprints the 0.98 / 0.99 / 0.97 bar labels. In the
    pre-change figure the bars were lower (0.97 / 0.98 / 0.96) and the labels only touched the dashed line. Fix
    in the figure producer by moving the annotation into the empty < 0.865 column or lowering it below 0.9.
- **Caption arithmetic** (pre-existing pattern, cosmetic). "548 of 6,660" EC non-answers comes from the comparison
  cells. "6,423 outright abstentions" comes from the gene-capture scoring (243 pred_empty there). In EC space
  6,660 − 6,423 = 237, not 243, so the 243 at L519 is the gene-capture count. Both round to 0.4%.

---------------------------------------------------------------------------------------------------------------

## 7. Style on changed lines

- **One sentence per line.** `style.py` checks every line of 2387ac2 not present at 0214a14. No changed line
  splits a sentence across lines. The following changed lines carry **two or three sentences** ("… . prFBA …").
  The pattern predates the rewrite but the lines were touched; fix by splitting before "prFBA":
  - L523 ("… at EC precision 0.669. prFBA, under …");
  - L533 (×2: "… are wrong. prFBA's precision …", "… is 0.231. prFBA instead …");
  - L625 ("… part of the panel. prFBA's precision …");
  - L887 ("… the cheaper marker. prFBA abstains …");
  - L1189 ("… does not have. prFBA's thresholds …").

  L834 is a `\paragraph{}` title plus one sentence, and L1073 ends in ".)"; neither is a violation. The rewrite
  did fix one old violation: "The arithmetic closes" now has its own line (L1071).
- **Section references.** Every section reference uses `\nameref`. No `Section~\ref`, bare `\ref{sec:…}`, "§" or
  hard-coded section numbers. New references (sec:recon, sec:concordance, sec:sharedref, sec:picrust2,
  sec:discussion, sec:limits, sec:magreal) all resolve.
- **Model identifiers.** None found. I searched for AI model names: Claude, Opus, Sonnet, Haiku, GPT, Anthropic,
  OpenAI, LLM. Software commit hashes in the text (ModelSEEDpy c2e31a0 and e1254ba, templates 21a8084, L676–683)
  are pre-existing commit identifiers and outside this rule.
