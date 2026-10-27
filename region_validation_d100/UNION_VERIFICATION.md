# Verification of prFBA_paper commit 8f24d9e ("Use the exact-tie union throughout")

Written 2026-10-27 by an adversarial verifier. Read-only on both repositories; this file is the only one written.

- Manuscript: `/home/freiburger/Documents/prFBA_paper/prFBA.tex` at 8f24d9e (1,138 lines), branch `depth100-standard`.
- Voiced twin: `prFBA_andrew_voice.tex`. The commit applied byte-identical change hunks to both files (the `-`/`+` lines of `git show 8f24d9e -- <file>` are identical). The two files still differ only at L121–124 and L239; L238 moved to L239 because of the new Methods sentence at L184, which is expected.
- **Every "current" substring below occurs exactly once in each file, on the same line in both, and contains no newline** (checked by script). Each correction therefore applies unchanged to both files.
- Every number was re-derived independently with `~/Documents/py_venv/bin/python`. Where a stored summary keeps 4 dp and sits on a rounding half, the value was recomputed from the cells. As a control, the main derivations were also run on the `_legacy` trees, and they reproduce the pre-commit text. That covers tab:main, McNemar, the single-genome shares, tier capture, the MAG pairing, Welch, the emptied strata, the pooled survivor rates, three-panel genus and family, and the AF-gate 31.8%.
- Not re-derived: bootstrap confidence intervals and the "no interval computed" clauses. These were read from the stored summaries (`ci95`, `ci`) or checked as absent contrast keys; they were not re-bootstrapped. The exception is tab:threepanel, which was re-run in full; see note N2.
- Repo paths are relative to `/home/freiburger/Documents/prFBA` unless absolute.

## Summary

| Severity | Count | Items |
|---|---|---|
| wrong number | 0 | every number changed or added by 8f24d9e reproduces. The added text holds about 580 numeric tokens; that count includes labels such as 16S and thresholds such as 0.987. |
| wrong claim | 2 | 1 (L634), 3 (L736) |
| inconsistency | 3 | 2 (L184), 4 (L68/L999), 5 (F14 panel C) |
| wording | 3 | 6 (L529), 7 (L737), 8 (L1105, optional) |
| pre-existing, not introduced by 8f24d9e | 1 | 9 (L363, L417 against L72; optional) |

There are also notes N1–N4, which need no manuscript change.

The two deliberate deviations from the change list are correct:

- **Family-rank gap −0.123.** `data_family/cells.csv.gz` gives prfba_legacy − picrust2 = −0.1234591, and the stored `prfba_nestedcv__vs__picrust2.ec.delta_f1_population` is −0.1235. So −0.123 is right and −0.124 was a double rounding.
- **Full-length include-self capture 99.98%.** `headline_macro_recall_pct.include.FullLength16S` = 99.979, and `regions.FullLength16S.include.macro.recall.value` = 0.9997861 gives 99.97861%. "Over 99%" and "exceeds 99%" hold, with a minimum of V4 at 99.1477%.

## Defects, most severe first

### 1. L634: the legacy reducer is not "the one rule scored identically at both depths" (wrong claim)

Current:

~~~latex
and under the legacy reducer, the one rule scored identically at both depths, conditional genus accuracy moves by less than 0.01 in every region
~~~

Corrected:

~~~latex
and under the legacy reducer, the only rule with a depth-20 \Cref{tab:main} scoring, conditional genus accuracy moves by less than 0.01 in every region
~~~

Evidence:

- **Both rules are scored at both depths on the as-shipped axis.** `region_validation/picrust2_comparison/data/comparison_summary.json` (depth 20) and `region_validation_d100/picrust2_comparison/data/comparison_summary.json` (depth 100) both carry `_meta.methods` = picrust2, prfba_exclude (exact-tie union), prfba_exclude_legacy, ….
  - `by_region.ec.all.all.prfba_exclude.conditional.f1.value`: 0.7985 (d20) and 0.8017 (d100).
  - The same key's `population.f1.value`: 0.5996 and 0.7130.
  - The union selections abstain on 15,888 and 6,457 amplicons: `picrust2_comparison/data/prfba_selection.json` at each depth, entries whose `exclude_0.987.reps` is empty. The legacy selections abstain on 16,890 and 7,344.
- **The sentence contradicts itself.** Its first clause quotes the union at both depths: "conditional F1 is 0.802 against 0.798 at depth 20". The abstract (L66), L614, L635 and L1088 also print union depth-20 figures.
- **What is unique to the legacy reducer** is the original panel's depth-20 tab:main tree. `region_validation/data/taxacc_region_summary.json` and `gene_capture_by_region.json` both have `_meta.select_mode` = "legacy". No depth-20 `default_union` tab:main tree exists for the original panel.
- **The legacy numbers themselves reproduce.** They come from `region_validation/data` against `region_validation_d100_legacy/data`.
  - Conditional genus accuracy moves by −0.0068 to +0.0091.
  - Conditional macro recall rises by +0.66 to +1.30 points.
  - Conditional F1 falls by 0.12 to 1.22 points.
  - Pooled population genus accuracy is 0.3749 at depth 20 and 0.4720 at depth 100.
- The commit message repeats the same false clause. That is out of scope for this file.

### 2. L184: the Methods exception list omits the legacy depth-20 versus depth-100 comparisons (inconsistency)

Current:

~~~latex
except the pre-registered primary prFBA arm of the controlled benchmark, whose policy was chosen by leave-one-fold-out selection (the \nameref{sec:sharedref} section).
~~~

Corrected:

~~~latex
except the pre-registered primary prFBA arm of the controlled benchmark, whose policy was chosen by leave-one-fold-out selection (the \nameref{sec:sharedref} section), and the depth-20 versus depth-100 comparisons of the \Cref{tab:main} figures (the \nameref{sec:depth} and \nameref{sec:limits} sections), which use the ordered reducer because it is the only rule with a depth-20 \Cref{tab:main} scoring.
~~~

Evidence:

- The sentence says "every result in this paper uses the exact-tie union except" the controlled-benchmark arm.
- Four passages report legacy-scored results, each correctly labelled. All of them reproduce (see item 1).
  - L634: conditional genus accuracy within 0.01, and capture +0.7 to +1.3 points.
  - L635: "population genus accuracy, compared under the legacy reducer, rises from 0.375 to 0.472".
  - L636: "the legacy-reducer depth-20 base of 0.375".
  - L1105: "under the legacy reducer, scored identically at both depths, …".
- These are the authors' deliberate exception (c), but the paper's one statement of its rule policy does not admit them.
- Both `\nameref` labels exist: `sec:depth` at L611 and `sec:limits` at L1097.

### 3. L736: under the union, "the same ordering" no longer holds in recall (wrong claim)

Current:

~~~latex
The same ordering holds in recall alone (MAG 0.816--0.824 against 0.703--0.768), so the advantage is not an artifact of precision.
~~~

Corrected:

~~~latex
The same grouping holds in recall alone (MAG 0.816--0.824 against 0.703--0.768; a margin of 0.055--0.083 against full-length 16S and the V1-anchored regions and 0.097--0.114 against the five from V3 onward), so the advantage is not an artifact of precision.
~~~

Evidence: I recomputed the 100% completeness pairing with compare_arms.py logic, using union `region_validation_d100/data/gene_capture_cells.csv` (exclude_0.987, ok) × `mag_validation/out/mag_cells_all.csv` (exclude, contamination 0, ok with non-empty precision). The 4-dp means match `paired_comparison_d100.json` `exclude|c1_x0|<region>|recall`.

| Region | MAG recall | 16S recall | Recall margin |
|---|---|---|---|
| FL | 0.821540 | 0.738997 | 0.082543 |
| V1-V2 | 0.822897 | 0.767589 | 0.055308 |
| V1-V3 | 0.822174 | 0.760151 | 0.062023 |
| V3-V4 | 0.823717 | 0.711597 | 0.112120 |
| V4 | 0.818091 | 0.711436 | 0.106655 |
| V4-V5 | 0.817440 | 0.703338 | 0.114102 |
| V4-V5_944R | 0.819804 | 0.706359 | 0.113445 |
| V7-V9 | 0.816420 | 0.719762 | 0.096658 |

- **F1 margin order:** FL < V1-V3 < V1-V2 < V7-V9 < V3-V4 < 944R < V4-V5 < V4. This is the order L735 describes, "with V4 the widest gap".
- **Recall margin order:** V1-V2 < V1-V3 < FL < V7-V9 < V4 < V3-V4 < 944R < V4-V5. Full-length 16S is no longer the smallest, and V4 is not the widest.
- **Under legacy the orders matched.** The `_legacy` recall margins were V1-V3 0.0381 ≈ V1-V2 0.0382 < FL 0.0430 < V7-V9 0.0621 < V3-V4 0.0787 < 944R 0.0795 < V4-V5 0.0836 < V4 0.0900, which is the legacy F1 order. The sentence was true before the switch.
- **Only the grouping survives:** max(FL, V1 regions) = 0.0825 < min(V3 onward) = 0.0967.
- The recall ranges in the sentence are correct: MAG 0.816420–0.823717 and 16S 0.703338–0.767589.

### 4. L68 and L999: the abstract states k-means and family bounds that the body never reports, and "at most 0.05" understates 0.052 (inconsistency, with rounding)

Abstract, current (L68):

~~~latex
and conditional family accuracy by at most 0.05, while conditional genus accuracy is higher on the $k$-means panel in every region, by up to 0.07,
~~~

Abstract, corrected:

~~~latex
and conditional family accuracy by at most 0.052, while conditional genus accuracy is higher on the $k$-means panel in every region, by up to 0.07,
~~~

Body, current (L999):

~~~latex
with intervals that exclude zero only for full-length 16S and V1--V3, and conditional gene-capture F1 is lower on the random panel by 0.013 to 0.028,
~~~

Body, corrected:

~~~latex
with intervals that exclude zero only for full-length 16S and V1--V3 (on the $k$-means panel it is higher than on the original in every region, by $+0.003$ to $+0.070$, and conditional family accuracy differs from the original's by at most 0.052 on either alternative panel), and conditional gene-capture F1 is lower on the random panel by 0.013 to 0.028,
~~~

Evidence: `three_panel/three_panel_comparison_d100.json`, cells with regime exclude_0.987 (phylum-standardized est).

- **`genus_cond_acc` `delta_kmeans_minus_original`**, all positive:

  | Region | est |
  |---|---|
  | FL | +0.0641 |
  | V1-V2 | +0.0377 |
  | V1-V3 | +0.0700 |
  | V3-V4 | +0.0027 |
  | V4 | +0.0337 |
  | V4-V5 | +0.0289 |
  | V4-V5_944R | +0.0438 |
  | V7-V9 | +0.0264 |

  So the abstract's "every region, by up to 0.07" is true, but no sentence or table in the body gives any k-means quality delta. L999 and tab:threepanel report random − original only.
- **`family_cond_acc`:** the largest |Δ| against the original is k-means − original at V1-V3, +0.0523. Random − original reaches +0.0419 at V1-V3. The body reports no family bound anywhere.
- **Rounding:** "at most 0.05" is below the true maximum of 0.0523. The `_legacy` maximum was 0.0358 (k-means − original at V3-V4), printed as 0.04.
- The abstract's other bounds hold:
  - gene-capture F1 at most 0.028: max 0.0278, random − original at V7-V9;
  - random panel "six of eight (−0.016 to +0.062)": −0.0161 at V7-V9 and +0.0622 at V1-V3.
- I re-ran `three_panel/compare_three_panels.py` with full-precision output into the scratchpad; see note N2.

### 5. F14 panel C: "p<10^{-164}" holds only at 100% and 90% completeness (inconsistency; needs a re-render, not a .tex edit)

This item is in `manuscript_figures/mag/plot_results_d100.py` at L89, whose output is copied to `images/prfba/F14_mag_completeness.pdf`.

Current:

~~~python
"%+.3f vs full-length 16S, %+.3f vs V4 (paired means)\npaired within organism, $p<10^{-164}$"
~~~

Corrected:

~~~python
"%+.3f vs full-length 16S, %+.3f vs V4 (paired means at 100%% completeness)\npaired within organism, $p<10^{-164}$"
~~~

Evidence:

- The annotation's means are the 100% values (`_c1fl`, `_c1v4` = `exclude|c1_x0|…`): +0.103 and +0.212. The panel, however, plots 30–100%.
- The full-length Wilcoxon p by completeness, from `mag_validation/out/crossover_d100.json` `FullLength16S_p` and my recomputation:

  | Completeness | FL p |
  |---|---|
  | 100% | 4.0e-165 |
  | 90% | 1.2e-166 |
  | 70% | 1.5e-160 |
  | 50% | 2.0e-164 |
  | 30% | 5.4e-159 |

- So the bound fails at 70%, 50% and 30%.
- Under legacy the annotation read p<10^-179, and it held at every level (largest full-length p 3.35e-211 at 30%).
- Alternative fix: keep the wording and use `$p<10^{-158}$`, which holds at every plotted level.
- The figure's other content matches the text: n = 1,912–2,299, "+0.05 easier" (0.053356) and "exact-tie union".

### 6. L529: "0.769 versus 0.767" is now ambiguous in tab:sharedref (wording)

Current:

~~~latex
so the 0.769 versus 0.767 conditional F1 in \Cref{tab:sharedref} is the policy difference
~~~

Corrected:

~~~latex
so the 0.769 versus 0.767 conditional F1 of the \emph{prFBA} and \emph{random abstain} rows of \Cref{tab:sharedref} is the policy difference
~~~

Evidence:

- Before 8f24d9e, only *random abstain* printed a conditional F1 of 0.767 (`summary_genus.json` `arms.base_random_abstain.ec.conditional_matched.f1` 0.7674, n 50,946).
- The commit added the row "prFBA, shipped exact-tie", which also prints 0.767 (`arms.prfba_exact_tie.ec.conditional_matched.f1` 0.7672, n 57,021).

### 7. L737: "Every comparison" also covers the recall comparisons, whose p is far weaker (wording)

Current:

~~~latex
Every comparison is paired within organism and significant at $p < 10^{-161}$.
~~~

Corrected:

~~~latex
Every F1 comparison is paired within organism and significant at $p < 10^{-161}$ (in recall, $p < 10^{-29}$).
~~~

Evidence: `paired_comparison_d100.json` at 100%.

- **F1:** the largest `exclude|c1_x0|<r>|f1.conditional.p` is 6.28e-162 (V1-V3), so the bound holds. V4, V4-V5 and 944R underflow to `null`, and V7-V9 is 4.2e-303.
- **Recall:** the `exclude|c1_x0|<r>|recall.conditional.p` values are much weaker, the largest being 5.35e-30 (V1-V2) and 3.92e-38 (V1-V3). They are reached from the preceding sentence about recall.
- This ambiguity pre-dates the commit (legacy recall p reached 2.7e-5 under "p<10^{-179}"), but the commit re-set the bound.

### 8. L1105: harmonize with item 1 (wording, optional)

Current:

~~~latex
under the legacy reducer, scored identically at both depths, conditional genus accuracy within 0.01 in every region
~~~

Corrected:

~~~latex
under the legacy reducer, the only rule with a depth-20 \Cref{tab:main} scoring, conditional genus accuracy within 0.01 in every region
~~~

Evidence:

- The clause is literally true, since the legacy reducer was scored at both depths. But it reads as the same justification as L634, and the union was also scored at both depths on the EC axis (item 1).
- The numbers reproduce:
  - conditional genus accuracy within 0.0091;
  - recall +0.66 to +1.30 points;
  - conditional gene-capture F1 −0.12 to −1.22 points.

### 9. L363 and L417 contradict the abstract (L72) on family-level metabolism (pre-existing inconsistency, not introduced by 8f24d9e; optional)

L363, current:

~~~latex
The practical reading is that family-level metabolism transfers robustly across all regions even for organisms absent from the database, while genus-level transfer is where short amplicons pay their penalty.
~~~

L363, corrected:

~~~latex
The practical reading is that family-level assignment holds up across all regions even for organisms absent from the database, although a family-tier answer recovers less of the metabolism than a genus-tier one, and genus-level transfer is where short amplicons pay their penalty.
~~~

L417 (fig:perrank caption), current:

~~~latex
so family-level metabolism transfers even for novel organisms.
~~~

L417, corrected:

~~~latex
so family-level assignment holds even for novel organisms.
~~~

Evidence:

- L72 (abstract) says that robustness "does not carry over to metabolism".
- L361 agrees with L72: family-tier answers recover 58% against 69% for genus-tier answers, from `gene_capture_cells.csv` (exclude_0.987, ok, macro recall by tier). Pooled, family is 58.218% (n 8,451) and genus 69.286% (n 51,120). By region, family runs from 55.663% (V4-V5_944R) to 60.668% (V1-V2).
- The same tension existed under legacy (61% against 72%).

## Notes (no manuscript change required)

- **N1. L181 and L183 keep ordered-reducer numbers that fall outside the three stated exceptions.** They are design rationale for the union, and all of them reproduce:
  - "mean 2.88 genomes against 1.32": mean include `n_selected` 2.8788 (union, all 67,991 cells) against 1.3178 (legacy, answered cells);
  - "V4 0.795 against 0.856": `regions.V4.include.macro.precision.value` 0.794532 against 0.856488;
  - "~19% (12,602 of 67,466)".

  With item 2 applied, L184 still reads "every result"; the authors may want "every reported result" or may consider L181/L183 covered by "remains available as an option". This is an author decision, not a correction.
- **N2. tab:threepanel half-value cells.** Five printed bounds rest on 4-dp stored values that sit exactly on a 3-dp half:

  | Row | Bound | Stored | Printed |
  |---|---|---|---|
  | FL genus coverage, random − original | CI upper | −0.1905 | −0.191 |
  | FL genus quality | CI lower | +0.0125 | +0.013 |
  | V3-V4 genus quality | CI lower | −0.0385 | −0.039 |
  | V7-V9 genus quality | CI upper | +0.0125 | +0.013 |
  | V4 gene-capture coverage | CI lower | −0.2905 | −0.291 |

  L1048's "−0.115" (stored −0.1145) is the same case. The printed values follow round-half-away-from-zero. `compare_three_panels.py` rounds `est`/`ci` to 4 dp at L252, so the stored file cannot settle these.

  I re-ran a scratch copy of `three_panel/compare_three_panels.py` with that rounding removed (`--suffix _d100`, output in the scratchpad only).
  - **Rates and point estimates** reproduce the stored file exactly: every crude and phystd rate, and all 288 `est` values.
  - **Bootstrap CI bounds do not reproduce.** Across 576 bounds the maximum |stored − re-run| is 0.0045 and the median 0.0004. For example, FL genus coverage gives [−0.2479, −0.1878] against the stored [−0.2507, −0.1905], and L1048 gives [−0.1141, −0.0999].
  - **The cause:** `boot_phystd` (L153–159) draws strata in the iteration order of `wphy`. That dict is built from the set `common_phy` of phylum-name strings, so the RNG stream is consumed in a hash-seed-dependent order, even with `SEED = 1729`.
  - **What this means for the paper:**
    - The third decimal of every tab:threepanel interval is Monte Carlo noise, so the half-value rounding above is immaterial.
    - Every printed "interval excludes zero" verdict is the same in both runs:
      - L997: the coverage gaps;
      - L999: genus quality excludes zero only for FL and V1-V3;
      - tab:threepanel: gene-capture-F1 quality.
    - Three unprinted verdicts flip between runs: random − k-means FL genus quality, k-means − original V1-V2 include gene-capture coverage, and k-means − original V7-V9 genus quality.
  - No manuscript change. However, the printed interval bounds (for example [−0.251, −0.191]) are reproducible only to about ±0.004. A re-run of the script will visibly shift the third decimals of tab:threepanel's intervals and L1048's interval. That is expected and does not mean the script is broken.
  - The script could iterate `sorted(common_phy)` (or a fixed phylum order) in `boot_phystd` to make the intervals reproducible. Any later re-run would then still change the bounds once, by up to about 0.004.
- **N3. Two "no gene families" shares are both correct.** L349's 1.6% (PGFam, 1,060/67,038) and L1082's 1.5% (EC, 1,015/67,038) are different quantities.
  - L349 comes from `gene_capture_cells.csv` exclude outcomes: pred_empty 1,060, abstained 6,407, ok 59,571.
  - L1082 comes from `comparison_cells.csv.gz`, prfba_exclude, ec: 7,422 non-answers, of which 6,407 have no reps and 1,015 have reps.

  Do not "fix" one into the other.
- **N4. Stale comments in the analysis repo, outside the manuscript.** `manuscript_figures/mag/plot_results_d100.py` L24 and `plot_crossover_d100.py` L21 read the union inputs (`paired_comparison_d100.json`, `crossover_d100.json`), but their comments still say "SCRATCH PATCH: depth-100 amplicon cells, legacy reducer". Not edited.

## Verified correct

All six re-rendered figure PDFs in `images/prfba/` are byte-identical (md5) to `manuscript_figures/out/`, and those scripts read the union trees: `rv/_common.py` has `BASE = region_validation_d100`, and `mag/*` reads `*_d100.json`.

The blocks below were each recomputed and compared with the printed text. TM = `region_validation_d100/data/taxacc_region_summary.json` `regions[r].anchor.<mode>.all.all.per_rank`; GC = `gene_capture_by_region.json`.

### tab:main: 56 cells and the caption, all correct

TM `Genus.correct` (include, excl.), `Genus.correct_population`, `Family.correct` (exclude_0.987). GC `headline_macro_recall_pct` (include, exclude_0.987) and `headline_macro_recall_pct_population`.

| Region | Genus incl. | Genus excl. | Genus pop. | Capture incl. | Capture excl. | Capture pop. | Family excl. |
|---|---|---|---|---|---|---|---|
| FL | .996484 | .726323 | .537691 | 99.97861 | 68.72283 | 63.27453 | .847912 |
| V1-V2 | .949595 | .702579 | .507398 | 99.50964 | 69.90347 | 64.31044 | .817022 |
| V1-V3 | .968789 | .706275 | .512821 | 99.64522 | 69.51789 | 64.44319 | .820248 |
| V3-V4 | .921914 | .614282 | .405550 | 99.36205 | 66.67865 | 58.36503 | .776615 |
| V4 | .890933 | .576329 | .370554 | 99.14773 | 67.30282 | 58.11545 | .759021 |
| V4-V5 | .917909 | .601488 | .390884 | 99.31124 | 66.38296 | 57.64335 | .781812 |
| V4-V5_944R | .917336 | .609185 | .400342 | 99.33083 | 66.53661 | 58.09854 | .791363 |
| V7-V9 | .938718 | .626009 | .423398 | 99.39640 | 67.29923 | 59.21378 | .796716 |

The V3-V4 population capture 58.365025 is confirmed from the cells: 8,253 cells give 58.3650254, so 58.4 is right.

### Text ranges derived from tab:main, about 30 values, all correct

The abstract (L64–65, L72), L353–356, L360–362, the four figure captions (fig:accuracy, fig:genecap, fig:leakage, fig:perrank) and the Discussion (L1076):

- **Ranges** are min and max over the table above: 0.89–1.00, over 99%, 0.58–0.73, 66–70%, 0.37–0.54, 58–64% and 0.76–0.85.
- **Precision:** V4 0.794532 → 0.79 and 0.795; FL 0.996679 → 1.00; V3-onward regions 0.7945–0.8799 → 0.79–0.88.
- **Gaps:**
  - genus 0.2470–0.2702 (FL, V1 regions) and 0.3076–0.3164 (V3 onward);
  - capture 29.6062 (V1-V2) to 32.9283 (V4-V5);
  - family − genus 11.4–18.3 points, and 0.1623–0.1827 from V3 onward ("one in six").
- **Tier capture:** family pooled 58.218 (55.663–60.668 by region); genus pooled 69.286.
- **Species (include):** 0.654326–0.877346, FL 0.991636.
- **Discussion (L1076), V3-onward regions:**
  - family conditional accuracy 0.7590–0.7967;
  - genus uncovered 0.324–0.357, "a third";
  - population capture 57.643–59.214;
  - population family accuracy 0.5359–0.5869;
  - population genus accuracy 0.3706–0.4234.
- **Discussion (L1076), full-length 16S and the V1-anchored regions:** population capture 63.275–64.443 and population genus 0.5074–0.5377.
- **Region-ordering sentences** (L354, L355 and "V4–V5 lowest 66.4") all hold.

### Denominators (L349), single-genome shares (L317) and McNemar (L357), all correct

- **Genus coverage** (TM `exclude_0.987` `Genus.coverage`): 0.642956–0.740292.
- **The further 22%:** 11,230/52,026 = 21.585%, from `taxacc_per_record.jsonl`. Here 11,230 genus-scoreable answered amplicons have an unnamed-genus anchor, and the abstentions are 5,365/52,026 = 10.3%.
- **Gene-capture zero shares** (GC `_meta.outcome_counts_total.exclude_0.987`): 7,467/67,038 = 11.14%; abstained 9.56%; pred_empty 1.58%.
- **Single-genome shares:**
  - 61,509 answered taxonomic amplicons with 41,894 single (68.110%), from `asv_region_concordance.csv`, anchor, exclude_0.987;
  - 61,534 answered gene-content amplicons with 41,916 single (68.118%), from `gene_capture_cells.csv`.
- **McNemar** (`mcnemar_region_pairs.csv`, anchor, Genus, exclude_0.987):
  - 26 of 28 pairs significant;
  - n_common 5,533–7,238, with `n_common_all8` 5,481;
  - not significant: V1-V2/V1-V3 (Holm p 0.141238) and V3-V4/V7-V9 (0.088866);
  - V4-V5/944R: b + c = 58 + 118 = 176, the fewest of any pair, with Holm p 6.3e-05;
  - FL/V1-V3 is now significant (Holm 0.006).

### Controlled benchmark: about 45 values, all correct

The exact-tie arm was added at L75–77, L79–80, L442, L515–544, tab:sharedref, L638, L1088–1091 and L1115–1116. Sources are `region_validation/shared_reference/data/summary_genus.json` and `data_family/summary_family.json`, with half values recomputed from `cells.csv.gz` (density 1.0, ec, scoreable).

| Quantity | Leave-genus-out | Leave-family-out |
|---|---|---|
| Exact-tie population F1 | 0.6849333 | 0.6234966 |
| Exact-tie forced population F1 | 0.7455385 | 0.7072431 |
| Exact-tie − PICRUSt2 | −0.0937176 [−0.109, −0.0829] | −0.1259304 [−0.1522, −0.1017] |
| Exact-tie forced − PICRUSt2 | −0.0331124 [−0.0392, −0.0291] | −0.0421839 [−0.0503, −0.0334] |
| Region range, abstaining | −0.183766 to −0.059845 | −0.247474 to −0.077339 |
| Region range, forced | −0.036436 to −0.030725 | not printed |
| Domains, abstaining | Bacteria −0.096562, Archaea −0.052791 | not printed |
| Domains, forced | Bacteria −0.034055, Archaea −0.019542 | not printed |
| Exact-tie − nn1 | −0.0601346 | −0.0834880 |
| Exact-tie − random abstain | +0.0188811 | +0.0152767 |
| Exact-tie − constant | not printed | −0.0197990 |
| Exact-tie − family donor | not printed | +0.1041749 |
| Forced − abstaining | +0.0606051 | not printed |
| Legacy (nested CV) − PICRUSt2 | not changed | −0.1234591 |

Further leave-genus-out values:

- AURC 0.7868 and cross-knob 0.7862;
- skill 0.2702 (abstaining) and 0.6922 (forced);
- conditional F1 0.7672;
- precision 0.803484 in the 0.945–0.987 band and 0.751083 in the 0.865–0.945 band;
- forced population F1 below 0.865: 0.566899;
- density slope 0.106 (legacy 0.1068, PICRUSt2 0.0301);
- a fourfold reference cut costs 0.0638 (exact-tie), 0.0643 (legacy) and 0.0181 (PICRUSt2);
- FL/V4 population F1 0.6901/0.7076 at coverage 0.8898/0.9329, the same coverage under both rules;
- the exact-tie conditional gap is −0.0288.

Further leave-family-out values: skill −0.1606, and nn1 0.7069847 (0.7070).

### MAG versus 16S: 72 + 63 table cells and about 25 text values, all correct

I recomputed with compare_arms.py and crossover.py logic: union amplicon cells × `out/mag_cells_all.csv` + `lowc_out_batches/b*/mag_cells.csv`, matching `paired_comparison_d100.json` and `crossover_d100.json`.

- **tab:magvs16s (100%):** n / MAG / 16S / Δ / median Δ / coverage MAG / coverage 16S / population MAG / population 16S.
  - FL: 1,916 / .813011 / .710014 / .102997 / .055620 / .304173 / .920699 / .247726 / .604584.
  - The other seven rows match the same way. The half cases V1-V3 16S .703525, V4-V5 Δ .192489, median Δ .157451 and population 16S .503519, and V7-V9 16S .640535 all round as printed.
- **tab:magval:** all seven rows match. The half cases are 50% V4 MAG .809462, 25% V4 MAG .808467 and 2% FL MAG .793519.
- **Text values:**
  - minimum and maximum Δ over 14 levels × 8 regions: 0.076708 (1% FL) and 0.214472 (90% V4), all positive; at 100%, 0.102997–0.212059;
  - 1% FL wins 1,129/473 at p 1.38e-79;
  - "+0.103 to +0.116" and "+0.168 to +0.212";
  - recall ranges (see item 3);
  - 16S coverage 0.863463–0.926982;
  - amplicon F1 where the MAG abstains 0.635022 (FL) and 0.544152 (V4);
  - "0.053 easier": FL paired 0.710014 − all-ok 0.656658 = 0.053356;
  - τ = 99 FL Δ 0.127441 and population 0.300823;
  - AF gate 365/1,367 = 26.70% above 0.71 (`out/af_gate_counterfactual.json`);
  - 0.945–0.987 band amplicon answer rate 0.977–0.983.
- **Table ordering:** "ordering unchanged down to 1%" holds; the full F1 order is identical at all 14 levels.
- **Figures:** the F15 annotations (+0.077, +0.103, n 1,603–2,311, 0.308 → 0.253) match.

### Panels: about 120 values, all correct

- **tab:twopanel, 16S block:** all 24 values match the crude `genecap_coverage` and `genecap_cond_f1` in `three_panel_comparison_d100.json`. The half cases are V4 coverage .863492 and .647484 and V7-V9 random coverage .547493.
- **L920–922:**
  - random below original by at most 0.029470 (V7-V9);
  - k-means within 0.022918;
  - Welch tests, which I recomputed from the cells, original against k-means: p = 2.35e-4 (FL, +0.011234), max 8.27e-3 (V1-V3); V4 2.9e-9, V4-V5 1.9e-9, V7-V9 5.9e-13; V4-V5_944R 5.9e-6 is correctly not listed;
  - coverage 0.749631 and 0.650422 against 0.920721.
- **tab:emptied:** all 24 values (shares and uncovered rates, including the 9 that changed) and the n row match `three_panel/emptied_strata_d100.json`.
- **L957 and L960–964:**
  - FL k-means − original 0.1588 standardized against 0.1711 crude;
  - survivors uncovered 236/6,879, 188/6,796 and 137/5,836;
  - per-stratum uncovered 1.47–3.66%;
  - pooled over eight regions 0.031035, 0.024324 and 0.018352 (from cells: emptied = `n_removed_self` ≥ 100);
  - fold 29.1–42.6;
  - point differences 1.08 (FL) and 1.27 (pooled);
  - CMH p 0.00407 and 0.148;
  - coverage − emptied residual at most 1.7 points;
  - pred_empty survivors 107, 76 and 62.
- **L644:** emptied 8.292% and 29.368% pooled; survivors as above.
- **tab:threepanel, 16S blocks:** all 40 values match. The half cases are FL random .479476 and V3-V4 random .397470; see N2 for the CI bounds.
- **L996–999:**
  - V4 0.6153, 0.4308 and 0.3580;
  - −0.2573 [−0.2758, −0.2393] against −0.3296 crude;
  - standardized/crude ratios 0.764–0.921;
  - random − original genus quality −0.0161 to +0.0622, positive in six of eight, intervals excluding zero for FL and V1-V3 only;
  - gene-capture F1 −0.0128 to −0.0278.
- **L1048–1049 (include):**
  - V4 0.836842 → 0.729982 (−0.1069 [−0.1145, −0.0998]);
  - V7-V9 0.906473 → 0.805086;
  - genus accuracy −0.0792 (V4) and −0.0653 (V7-V9);
  - FL −0.0159 and −0.0064.
- **L1062 and L1064, abstract L69–70:** 0.2037–0.2899; maximum gene-capture-F1 |Δ| 0.0278 (MAG arm 0.0262, `three_panel_mag_panelarm.json`); include maximum 0.1069 and FL 0.0159.

### Rule statements

These all name the correct rule:

- L63, L343, L357, L367, L726, L750, L800, L844, L929, L969, L1007;
- L314, L548, L1071, L1126;
- tab:sharedref caption, L509 (exact-tie − legacy = −0.0017).

L181's "self-exclusion precision at or above the ordered reducer's" holds in all eight regions (V4 0.5360 against 0.5343 is the narrowest). L614 and L1104 (6,457; 5,598; 86.7%) are union.
