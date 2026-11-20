# Independent number verification of the R5 manuscript re-numbering

Verifier: independent agent, 2026-11-16. Report only: no tracked file was edited, nothing was committed or checked out.

Scope: `git -C ~/Documents/prFBA_paper diff 0214a14..2387ac2 -- prFBA.tex`. The working tree equals `2387ac2`.

Line numbers are those of the **new** file (`2387ac2`) unless marked "old".

Scratch code and outputs are in
`/tmp/claude-1000/-home-freiburger-Documents-prFBA/400ec395-f172-59da-a98d-3dbe6b3a4256/scratchpad/verify_r5/numbers/`
(`code/*.py`, `values_check.json`, `tables_*.json`, `mag_unrounded.json`, `reverse.txt`, `stale*.json`).

Independence: I wrote my own line/token alignment, selector resolver, table generators and MAG pairing code. I used
`RUN5B_REPORT.md` only to find where a line came from. I did not read or run anything in `scratchpad/r5b/`.

---------------------------------------------------------------------------------------------------------------

## Verdict

- **3 blocking errors.** At 4 printed values, the tex is not the correct rounding of the unrounded source under either
  rounding rule. Each is one unit in the last digit. All are listed under B1–B3.
- **9 cosmetic or convention items** (C1–C9).
- Every other checked number is correctly rounded and traces to an artefact, with one exception. The values 122 / 65 /
  ~115 at L591–592 have only a plan-document source (C6).
- The stale sweep finds no surviving old value in the same context, except ones explicitly labelled pre-registered,
  "of the time" or pre-repair.

### Blocking (the tex is wrong at printed precision)

| # | tex line | printed | correct | unrounded source | number_diff `new_fmt` |
|---|---|---|---|---|---|
| B1 | L630 | "exact-tie $-0.083$, no interval computed" (family-rank exact-tie shortfall against nn1) | **−0.084** | 0.623665 − 0.707190 = **−0.083525**. Population EC F1 of `prfba_exact_tie` and `base_nn1`, density 1.0, computed from `region_validation/shared_reference/data_family/cells.csv.gz` (git-ignored). These reproduce `summary_family.json` `arms.{prfba_exact_tie,base_nn1}.ec.population.f1.value` = 0.6237 / 0.7072. | −0.084 (**right**). RUN5B kept −0.083 as "float noise", but the unrounded value is past the boundary. |
| B2 | L630 | "exact-tie $-0.077$ to $-0.247$" (largest region deficit, V1–V2) | **−0.248** | −0.247519, from the same cells. Stored value: `summary_family.json` `contrasts.prfba_exact_tie__vs__picrust2.ec.by_region.V1-V2.delta_f1_population` = −0.2475 (4 dp). | −0.247 (**also wrong**: float formatting of the 4-dp value) |
| B3 | L533, L577 | prFBA precision "falls only from 0.822 to 0.606" (L533), "while prFBA's falls only to 0.606" (F13 caption, L577) | **0.605** | 0.605479, from `region_validation_d100/picrust2_comparison/data/comparison_cells.csv.gz` (git-ignored): `prfba_exclude`, ec, ok, `nsti_bin` ≥0.5. Stored value: `comparison_summary.json` `by_stratum.ec.nsti_bin.>=0.5.prfba_exclude.conditional.precision.value` = 0.6055. | 0.606 (**also wrong**) |

The other end of each range is correct:
- −0.077 (−0.077345);
- the legacy range −0.076 to −0.245 (−0.075663 / −0.245253);
- −0.081 (−0.081103);
- 0.822 (0.821973).

B3 is not a stale value: the stored value moved 0.6058 → 0.6055, onto the boundary, and the print was kept. Both lines
were edited in this diff. The F13 figure itself may also show 0.606; I did not check it.

**B1–B3 are exact, not affected by per-cell rounding.**
- `comparison_cells.csv.gz` and the family `cells.csv.gz` store f1 and precision at 4 dp per cell. That is up to 5e-5
  per cell.
- So I recomputed every cell exactly from its integer `tp`, `n_truth` and `n_pred` (`code/exact.py`).
- The exact results:
  - B1 −0.0835246;
  - B2 −0.2475186;
  - B3 0.6054784.
- Each verdict stands.

### Cosmetic or convention items

1. **C1 — L78, L631: "$([-0.004, 0.040])$"** (family-rank refusal credit over a random mask, CI low bound).
   - `summary_family.json` `contrasts.prfba_nestedcv__vs__base_random_abstain.ec.ci95[0]` moved −0.0044 → **−0.0045**.
     It is a bootstrap bound stored at 4 dp, so the true value cannot be determined.
   - Half-even gives −0.004; half-up gives **−0.005**.
   - The paper's other 4-dp CI bounds follow half-up:
     - −0.0525 → −0.053 in the same L631 sentence;
     - tab:threepanel −0.0385 → −0.039 and −0.0345 → −0.035;
     - the old tex's three boundary CIs (+0.0185 → +0.019, −0.2905 → −0.291, −0.0195 → −0.020).
   - So L78 and L631 are the only exception to that convention.
   - These lines are not in the diff.
2. **C2 — L901: "(0.357 against 0.361, from 0.250 against 0.612)".** The printed operands give 0.612 − 0.250 =
   **0.362**. 0.361 is the difference of the unrounded means (0.61159 − 0.25019 = 0.36140). Either print 0.362, or
   say that the gap uses the unrounded means.
3. **C3 — tab:threepanel L1143, the new 100% panel-as-reference row.**
   - It prints Δ −0.027, but its displayed columns give 0.783 − 0.809 = −0.026.
   - The Δ is correct from the unrounded phystd difference (−0.026520).
   - The same convention applies in the unchanged rows L1126, L1130, L1138 and L1139.
   - Unlike the tab:magvs16s caption, the tab:threepanel caption does not say that Δ uses the unrounded values.
4. **C4 — producer format mismatches; the tex is fine.**
   - PH-04 p/q are printed at 2 dp (0.41 / 0.51; unrounded 0.4075 / 0.5094), but number_diff uses f3.
   - The Spearman p is printed 1.5×10⁻⁷ (1.469e-7), but number_diff uses e0 ("1e-7").
   - The PH-02 q range 0.029–0.046 is printed, but number_diff stores only the maximum (0.0461; the minimum is
     0.02887).
5. **C5 — L519, "In this comparison … (243 remain)".**
   - 243 is the gene-family-space `pred_empty` count from `gene_capture_by_region.json`.
   - In the as-shipped EC comparison, the empty-prediction cells number 6,660 − 6,423 = 237.
   - The old text took its 1,060 from the same gene-capture source, so this is wording only.
6. **C6 — provenance that is not a committed artefact.** The values are correct, but they trace only to the
   uncommitted `prFBA_rescore_2026-11` workspace (not a git repo) or to git-ignored files:

   | value (line) | source | status |
   |---|---|---|
   | 200 of 7,855; 13 of 15 strata; 442 of 1,000; 62 of 200; 106 of 19,072; 0 of 2,684 (L638–639) | `control/r4/pheno_sample_diff.json`, `pheno_attrib_{atp,gf}.json`, cited by RUN4_LOCAL_REPORT §3 | all verified ✓ |
   | 123 of 161 (L591) | `expected/rescore_work/panel_len_check.json` `original.n_old_empty` 161, `n_cur_empty` 38 | ✓ |
   | 17 genomes; 23 of 30 (L592) | `expected/changed_sets.json` (`changed_eligible` = 17); `control/r3/r32_affected_refs.json` (23 affected of 30 fold refs, plus recipe_A/B/C) | ✓ |
   | 130 of 63,871 (L592) | my own recount: tracked `fold_plan.json` × `changed_eligible` gives 24/47/24/31/4 = 130; all 130 are `scored` in the ignored `cells.csv.gz`, whose total is 63,871 | ✓ |
   | 122; 65 archaea; ~115 (L591–592) | only `plan/PLAN.md` D1; **no artefact found** | not verified |
   | 1733.9197 / 1733.9535 at 25,000 contigs (L1019) | `expected/contig_cap_d6.json` | ✓ |
   | 10 of 49,990 (L1019) | my recount from the ignored `mag_validation_rnd/manifests/*.json` | ✓ |
   | as-shipped unrounded values (0.777, 0.999, 0.935, 0.774) | ignored `comparison_cells.csv.gz`; the tracked `comparison_summary.json` holds 4 dp only | ✓ |
   | 43,174; 2,372; 107 (L680, L687) | committed `derived_numbers.post.json`, built from the untracked `per_genome.jsonl` | ✓ |
7. **C7 — a stale committed document outside the manuscript.**
   - `region_validation/shared_reference/data/objective_tuning/verification/VERIFICATION.md` L215 still says
     est_shuffled − shipped = +0.0110 [+0.0031, +0.0222].
   - The body of that file is the dated 2026-11-03 record. A 2026-11-16 delta table sits at the top, and it does not
     list this contrast.
   - So the owner can either add a row to that table or annotate the body.
   - The re-run `out_c4_c5.json` (tracked) gives 0.011091 [0.003121, 0.022253].
   - The manuscript's L319 value, 0.0111 [0.0031, 0.0223], is right.
8. **C8 — unchanged lines that cannot be verified at stored precision.**
   - L603 PICRUSt2 AURC 0.8125 is printed 0.812; L604 prFBA cross-knob AURC 0.7875 is printed 0.788.
   - Both are consistent with half-even, but half-up would give 0.813 / 0.788.
   - The raw values are unchanged from the start tag. I found no finer source.
9. **C9 — L80, the abstract.** It keeps 43,140 (now labelled "of the pre-registered pass") and states the re-run's
   0 of 8,660, but not its 43,174 rows. Those appear only at L680 and L687. This is correct, but a judgment call.

---------------------------------------------------------------------------------------------------------------

## 1. Numbers on + lines

### Counts

- **Diff size:** 187 + lines and 177 − lines.
- **My numeric-token alignment of paired lines:** 551 tokens replaced, 82 added, 14 removed, and 758 unchanged tokens
  on changed lines.
- **Of the 551 replaced tokens:**
  - 472 are explained by a number_diff value at the claim's lines (old printed → my correctly rounded new value);
  - 19 are explained by the same quantity from another inventory item: the repeats the executor swept, L402, L419,
    L463, L523, L632, L645, L740, L887, L907, L919–923 and L1163, all confirmed;
  - 60 have no inventory value: tab:magval rows 50–2%, tab:main/tab:tau per-region capture cells, and rewritten prose.
    I verified each by table regeneration (section 2) or by hand (below).
- **Of the 82 added tokens:** 10 are explained locally and 33 match an inventory value elsewhere. I checked these 33
  by hand, because several are coincidental matches (e.g. +0.005 ↔ a "CI lo"). The remaining 39 are unexplained and
  are traced below.
- **number_diff values (1,568):**
  - **My own resolver reproduces all 1,568 `new` values exactly from the live artefacts.**
  - About 190 of them are `derived:` values. I read those from the committed `derived_numbers.post.json`; I did not
    recompute them.
  - I independently re-derived these derived families:
    - MAG pairings and coverage;
    - the AF gate;
    - phenotype counts and contrasts (from the summaries);
    - surviving-uncovered and CMH;
    - the as-shipped matched set (61,044).
  - All 78 live source hashes still match `_meta.sources`.
  - This confirms collect.py's evaluation, not each spec's choice of selector. That is why sections 1–2 also go to the
    artefacts directly.
- **Located in the new tex:**
  - 997 values print exactly number_diff's `new_fmt`;
  - 26 more sit on a 4-dp rounding boundary (section 3);
  - 22 needed manual resolution: rewritten wording, historical labels, or a changed precision;
  - 378 sit on unchanged lines and are unchanged at print;
  - 106 are located by a whole-file search;
  - 30 bounds and strings all still hold, except "within 0.016", which was rewritten to 0.017 (correct: 0.01672);
  - 4 are in scientific notation (3e-83, 1.5e-4 ×2 and 1.5e-7 — correct).

### Values typed by hand, traced individually (✓ = matches the source at printed precision)

| line | printed | source (path / key) | tracked? |
|---|---|---|---|
| L79, L641, L1195, L1203, L1232 | n = 759 (750 V4); EC-F1 Δ −0.031; 0.07 more media; Wilcoxon p = q = 0.79; ρ −0.19, p 1.5×10⁻⁷ | `shared_reference/data/phenotype_atp_c13/phenotype_summary.json` `contrasts.prfba_exact_tie__vs__picrust2.{controlled,as_shipped}`: n_pairs 759, delta_atp_hamming 0.0685, p 0.7909, spearman −0.1894. The rest comes from derived `pheno_n_paired_v4` 750, `pheno_ecf1_delta_vs_picrust2` −0.03133 and `pheno_spearman_p` 1.469e-7. ✓ | yes / derived yes |
| L79, L636, L646, L1195 | four of 22 media (Ac.NO, Ac.O2, Etho.O2, For.O2); q 0.029–0.046 | same file, `mcnemar_per_medium`: 22 testable; q ≤ 0.05 for exactly these four (0.0289, 0.0396, 0.0461, 0.0461); all have c > b (favour PICRUSt2) ✓ | yes |
| L642, L1196 | 0.030 fewer media; n = 759; p = 0.41; q = 0.51; all 1,000: +0.06, p = 0.46 | `contrasts.prfba_exact_tie__vs__base_nn1`: −0.030303, 759, p 0.4075, q 0.5094, all-genome +0.06, p 0.4561 ✓ | yes |
| L1196 | EC F1 0.782 vs 0.780, p = 0.84 | derived `pheno_nn1_ecf1_{prfba,nn1,p}` 0.78197 / 0.78022 / 0.836 ✓ | derived yes |
| L651 | essentiality 0.063 vs 0.063; growth 0.020 vs 0.011 (PICRUSt2 −0.003) | `phenotype_gf_c13/phenotype_summary.json` controlled medians: 0.0632 / 0.0632; 0.0199 / 0.0108 / −0.0028 ✓ | yes |
| L80, L680, L687 | 0 of 8,660; 43,174 rows; all 2,372 grow; 107 reversals | `phenotype_gf_c13` meta `n_models` 8,660, `n_gated_out` {}; derived `pheno_atp_rows` 43,174, `pheno_gf_models_growing` 2,372, `pheno_gf_models_reversal` 107 ✓ | yes |
| L638–639 | 200 of 7,855; 13 of 15; 442 of 1,000; 62 of 200; 106 of 19,072; 0 of 2,684 | `control/r4/pheno_sample_diff.json` (delta_f1_changed 200; strata 15, 2 identical; picked_both 442 / 62) and `pheno_attrib_{atp,gf}.json` (19,072 compared / 106 changed; 2,684 / 0) ✓ | **no** (C6) |
| L591–592 | 123 of 161; 17; 23 of 30; 130 of 63,871 | see C6 ✓ | partly |
| L591–592 | 122 pass filters; 65 archaea; ~115 | PLAN.md D1 only | **no artefact** |
| L1019 | 1733.9197, 1733.9535; 25,000; 10 of 49,990 | see C6 ✓ | no |
| L527–528, L1197, L1226 | 0.723 vs 0.729 (−0.006); legacy 0.734; +0.005 | from the cells: union 0.722782, PICRUSt2 0.728971, legacy 0.734337; Δ −0.00619 / +0.00537. `comparison_summary.json` holds `prfba_exclude_legacy` ✓ | summary yes |
| L1197 | −0.120 at depth 20 | 0.6088 − 0.729 ✓ | yes |
| L516, L519 | 278 amplicons from 38 genomes; 675 more (67,713); 243 remain; 0.889 → 0.902; 0.713 → 0.723; 0.802 | cells: 278 truth_empty from 38 genomes; 67,713 − 67,038 = 675; `gene_capture_by_region.json` pred_empty 243; coverage 0.901644 ✓ | yes |
| L519 | panel capture +0.7–1.1 points; resolved ≤ 0.03 points | old (`start/`) vs new `gene_capture_by_region.json` `macro_pop` / `macro` recall: +0.672 … +1.14495 points, resolved max \|Δ\| 0.0286 ✓ | yes |
| L67, L736, L739, L1214 | 0.829 / 0.904 / 0.791 vs 0.798 / 0.834 / 0.790; 80.2%; 60.1% | `region_validation/picrust2_comparison/data/depth_probe/depth_probe_results.json` `ec_space_of_recovered` 0.8288 / 0.904 / 0.7914; anchors 195/243 = 0.8025, 134/223 = 0.6009 ✓. The answered set's 0.651 / 0.800 anchors were not recomputed (taxonomy). | yes |
| L742 | ~0.77 | recipe 0.609 + 0.221 × 0.884 × 0.829 = 0.771 ✓ | yes |
| L1068–1073 | 2.4 / 2.3 / 1.8%; at most 3.0%; none of 274; 1.8 / 1.8 / 1.3% pooled; 42–55-fold; 0.6 (0.5) points; CMH 0.117 / 0.99; 33.3 / 22.9 / 4.6%; 1 − 0.333 = 0.667 vs 0.655, etc.; within 2.3 points; 36 / 42 / 31 | `three_panel/emptied_strata_d100.json`: survivors 166/6932, 154/6798, 106/5853; k-means 50–99 rate 0.0299 = 40/1339; original 50–99 n = 274, rate 0; 1/0.0239 = 41.8 and 1/0.0181 = 55.2; CMH 0.11667 / 0.99475; non_ok − strict = 36 / 42 / 31. Max residual 0.954 − 0.931 = 0.023 ✓ | yes |
| L1095–1097 | 0.6889 (3,951) vs 0.6823 (9,702), +0.007 [−0.011, +0.024]; 0.559; 0.657 vs 0.682; 0.686 "slightly above" | derived `panels_mag_*` ✓; arithmetic 0.0066 ✓ | derived yes |
| L885–887 | 0.314 → 0.319; 0.454; 0.313; 1,391; 0.600 vs 0.807; 26.5%; 0.712; 45%; 0.743 | `mag_validation/out/af_gate_counterfactual.json`: 1,391; 0.600306; 0.806612; 369/1391 = 26.5%; level 0.7116; 1391/3121 = 44.6%; pooled (3121 × 0.8066 + 1391 × 0.6003) / 4512 = 0.7430. derived `mag_af_*` ✓ | yes |
| L831 | resolve 0.311–0.314; cond F1 0.805–0.807 (100–30%) | `mag_cells_all.csv` exclude, x0, main grid: 0.31118–0.31379; 0.80493–0.80699 ✓ | yes |
| L1028 | ≤ 0.030 (V7–V9); within 0.017; k-means 0.023 / 0.010; original ahead everywhere | `three_panel_comparison_d100.json` crude genecap_cond_f1: −0.02966 (V7–V9), k-means max 0.02314; `three_panel_mag_panelarm.json`: max 0.01672 / 0.00943; all negative ✓ | yes |
| L319 | 0.0111 [0.0031, 0.0223] | `objective_tuning/verification/out_c4_c5.json` `contrasts["est_shuffled/lin/G0/X1 - shipped"].taf_face_mid.pooled.A` = 0.011091 [0.003121, 0.022253] ✓ | yes |
| L358 | 0.0412 | `pfunc_estimators/scores.json` `composites.d100.all.contrasts_vs_E0.E4(E0).union_brier.est` = −0.0411509 ✓ | yes |
| L360 | +0.0085 [+0.0072, +0.0097] | `pfunc_estimators/verification/p4_results.json` `strata.all.{est,ci}.E3.macro_dF1_vs_union_restricted_to_union` = 0.008494 [0.007182, 0.009674] ✓ | yes |
| L529 | 0.777 / 0.768 / 0.774 vs 0.745 / 0.732 / 0.732; reversal only in the 5 regions from V3; rxn 0.805 / 0.722 | cells: 0.776505 / 0.767659 / 0.774456 vs 0.745039 / 0.731570 / 0.732227; the negative regions are exactly V3–V4, V4, V4–V5, 944R, V7–V9; rxn 0.80480 / 0.72161 ✓ | ignored cells |
| L529 | AS-15 divergences 0.026 / 0.014 / 0.05 / 0.03 (RUN5B says these were not recomputed) | recomputed on the new cells: archaeal max over regions 0.02615 (PICRUSt2, 944R); <0.865 bin 0.0137; prFBA no-NSTI cell (9) 0.0488; PICRUSt2 no-different-organism cell (8) 0.0266. **All still hold.** | ignored cells |
| L533 | 89.2% vs 93.5%; 629; 0.231; ~141 / ~410 / ~123 | cells: 629/705, 18738/20038; 0.2308; 140.9 / 409.7 / 123.1 ✓ (0.606 → see **B3**) | ignored cells |
| L577 | 8% (548 of 6,660); 147; 695 (10.8%) of 6,423; 0.234 | 548/6660 = 8.2%; 695/6423 = 10.8%; 0.233862 ✓ | yes / cells |
| L629–630 | 0.750 / 0.626 / −0.123; exact-tie 0.624, −0.126; forced 0.710, −0.040; exact-tie 0.707, −0.042; 0.7096 / 0.7072 / 0.7074; −0.081 | family cells: 0.749517 / 0.626087 / −0.12343; 0.623665 / −0.125852; 0.709579 / −0.039938; 0.707435 / −0.042082; −0.081103 ✓ (except **B1, B2**) | ignored cells / summary yes |
| L606 | exact-tie rxn −0.092 | genus cells −0.091529 ✓ (number_diff's −0.091 is wrong) | ignored cells |
| L901 | 0.661; 0.304; 0.357 vs 0.361 | see **C2** | |

---------------------------------------------------------------------------------------------------------------

## 2. Tables, regenerated with my own code

I regenerated every numeric cell of the tables that changed: **379 cells, 0 mismatches.** Two cells sit on a 4-dp CI
boundary; both follow the half-up convention.

| table | cells | source and method | result |
|---|---|---|---|
| tab:main (gene-capture columns) | 24 | `region_validation_d100/data/gene_capture_by_region.json`: `regions.<r>.include.macro.recall`, `exclude_0.987.macro.recall`, `macro_pop.recall`, unrounded, ×100 | 24/24 |
| tab:tau (capture cond / pop at 0.987 and 0.995) | 32 | `gene_capture_tau_compare.json` `regions.<r>.<mode>.recall_{resolved,population}.value`. The 0.987 values are identical to by_region's. | 32/32 |
| tab:picrust2 | 35 | `comparison_cells.csv.gz` (ec, truth_empty excluded): coverage; mean recall, precision and F1 over ok; population F1. n = 67,713 / 19,843 / 47,870. | 35/35 |
| tab:twopanel | 42 | `three_panel_comparison_d100.json` `[region, exclude_0.987, genecap_{coverage,cond_f1}].<panel>.crude`; `three_panel_mag_panelarm.json` `[exclude, comp].<panel>.{resolve_rate_NOT_COMPARABLE, cond_f1_crude}` | 42/42 |
| tab:emptied | 27 | `emptied_strata_d100.json`: shares from counts n/N; uncovered rates from round(rate × n)/n; N from the header | 27/27 |
| tab:threepanel (16S genus, 16S gene capture, MAG panel-as-reference blocks) | 84 | phystd values; Δ = unrounded random.phystd − original.phystd. I checked that it matches `est` in all 288 cells. CI bounds are 4 dp. | 84/84, 2 at a boundary (−0.0385 → −0.039, −0.0345 → −0.035, half-up) |
| tab:magvs16s | 72 | **own pairing:** `mag_cells_all.csv` (exclude, c1.0, x0) × `region_validation_d100/data/gene_capture_cells.csv` (exclude_0.987). The conditional set is both ok (MAG ok requires a precision). Population scores abstention as 0. Δ is the difference of unrounded means; the median Δ is the even-n midpoint. | 72/72 |
| tab:magval | 63 | the same pairing over `mag_cells_all.csv` plus the 20 git-ignored `lowc_out_batches/b*/mag_cells.csv`; resolve rate and cond F1 over scoreable MAGs; V4 pairing 2,006–2,448 (2,441 at 100%) | 63/63 |

The GTDB block of tab:threepanel was not regenerated. Its source is byte-identical to `start/` and its line is
unchanged; I spot-checked its CIs: 0.1066 → +0.107, …

Not regenerated, because no cell changed:
- the taxonomy columns of tab:main and tab:tau;
- tab:sharedref;
- tab:mag;
- tab:panelcomp.

---------------------------------------------------------------------------------------------------------------

## 3. Rounding: every claim whose print sits on a boundary, or where number_diff's value ≠ the tex

**Stored per-cell precision** (this decides whether a mean is truly unrounded):

| file | precision of f1 / precision | how the means below were computed |
|---|---|---|
| `mag_cells_all.csv`, `lowc_out_batches/*/mag_cells.csv` | ≥ 18 decimals | read directly |
| `region_validation_d100/data/gene_capture_cells.csv` | 8–9 decimals | read directly |
| `comparison_cells.csv.gz` and the shared-reference `cells.csv.gz` | **4 decimals** | recomputed exactly from integer tp / n_truth / n_pred (`code/exact.py`) |

- The smallest margin to a rounding boundary is V4 MAG: 0.808503 against 0.8085, 3×10⁻⁶.
- The worst-case error from the 8-dp amplicon cells is under 5×10⁻⁹, so this verdict stands too.
- Exact values for the 4-dp-cell claims:
  - FL population F1 0.7765048;
  - V1–V3 0.7744559;
  - include-self precision 0.9354923;
  - NSTI ≥0.5 precision 0.6054784;
  - FR-03 −0.0835246 / −0.2475186.
- So storing unrounded means is sufficient for the producer; the cells already carry the integers.

Settled from the unrounded source:

| claim (tex line) | tex | number_diff `new_fmt` | stored | unrounded | half-even / half-up | verdict |
|---|---|---|---|---|---|---|
| MAG V4, MAG column (L870, L918) | 0.809 | 0.808 | 0.8085 | 0.808503 | 0.809 / 0.809 | **tex right**, producer wrong |
| magvs16s V1–V3, 16S column (L868) | 0.704 | 0.705 | 0.7045 | 0.704464 | 0.704 / 0.704 | tex right, producer wrong |
| magvs16s V7–V9 median Δ (L873) | +0.125 | +0.126 | 0.1255 | 0.125462 | +0.125 / +0.125 | tex right, producer wrong |
| magval 1% FL MAG (L924) | 0.785 | 0.784 | 0.7845 | 0.784522 | 0.785 / 0.785 | tex right, producer wrong |
| as-shipped FL population F1 (L82, L529) | 0.777 | 0.776 | 0.7765 | 0.776505 | 0.777 / 0.777 | tex right, producer wrong |
| source-absent PICRUSt2 coverage (L565) | 0.999 | 1.000 | 0.9995 | 47,845/47,870 = 0.999478 | 0.999 / 0.999 | tex right, producer wrong |
| tab:emptied random 50–99 share (L1089) | 14.5 | 14.4 | 0.1445 | 1,269/8,780 = 0.144533 | 14.5 / 14.5 | tex right, producer wrong |
| panel-as-reference Δ at 100% (L1143) | −0.027 | −0.026 | −0.0265 | −0.026520 | −0.027 / −0.027 | tex right, producer wrong |
| FR-03 exact-tie vs nn1 (L630) | −0.083 | −0.084 | −0.0835 | −0.083525 | −0.084 / −0.084 | **tex wrong (B1)**; producer right |
| FR-03 exact-tie, largest region (L630) | −0.247 | −0.247 | −0.2475 | −0.247519 | −0.248 / −0.248 | **both wrong (B2)** |
| prFBA precision at NSTI ≥ 0.5 (L533, L577) | 0.606 | 0.606 | 0.6055 | 0.605479 | 0.605 / 0.605 | **both wrong (B3)** |
| SR-07 exact-tie rxn Δ (L606; line unchanged) | −0.092 | −0.091 | −0.0915 | −0.091529 | −0.092 / −0.092 | tex right, producer wrong |
| magvs16s V1–V2 MAG (L867) | 0.815 | 0.815 | 0.8145 | 0.814548 | 0.815 / 0.815 | both right |
| magvs16s V3–V4 population, 16S (L869) | 0.522 | 0.522 | 0.5225 | 0.522467 | 0.522 / 0.522 | both right |
| tab:picrust2 include-self precision (L560) | 0.935 | 0.935 | 0.9355 | 0.935492 | 0.935 / 0.935 | both right |
| as-shipped V1–V3 population F1 (L529) | 0.774 | 0.774 | 0.7745 | 0.774456 | 0.774 / 0.774 | both right |
| depth-probe family-correct (L736, L1214) | 80.2 | 80.2 | 0.8025 | 195/243 = 0.802469 | 80.2 / 80.2 | both right |
| family PICRUSt2 population F1 (L629) | 0.750 | 0.750 | 0.7495 | 0.749517 | 0.750 / 0.750 | both right |
| deployment mean selection size (L220, L1180, L1218) | 2.289 | 2.289 | 2.2895 | 8,597/3,755 = 2.289481 | 2.289 / 2.289 | both right; half-up on the 4-dp value would wrongly give 2.290 |

Indeterminate: bootstrap bounds stored only at 4 dp, with no finer source.

| claim (tex line) | tex | number_diff | stored | half-even / half-up | note |
|---|---|---|---|---|---|
| genus cond acc V3–V4 CI low (L1128) | −0.039 | −0.038 | −0.0385 | −0.038 / −0.039 | tex = half-up |
| tab:threepanel V7–V9 Δ cond F1 CI low (L1135) | −0.035 | −0.035 | −0.0345 | −0.034 / −0.035 | tex = half-up |
| FR-03 CI low vs base_constant (L631) | −0.053 | −0.052 | −0.0525 | −0.052 / −0.053 | tex = half-up |
| SR-08 PICRUSt2 vs nn1 CI low (L613) | 0.030 | 0.029 | 0.0295 | 0.030 / 0.030 | tex right under both; producer float error |
| family random-mask CI low (L78, L631) | −0.004 | −0.004 | −0.0045 | −0.004 / −0.005 | tex = half-even, **inconsistent with the paper's half-up CIs (C1)** |
| SR-05 PICRUSt2 AURC (L603) / prFBA cross-knob AURC (L604) | 0.812 / 0.788 | 0.812 / 0.787 | 0.8125 / 0.7875 | 0.812, 0.788 / 0.813, 0.788 | not a CI, but no finer source located; unchanged lines (C8) |

Format-only mismatches (producer side, not tex errors):

| claim (tex line) | tex | number_diff | note |
|---|---|---|---|
| PH-04 p / q (L642, L1196) | 0.41 / 0.51 | 0.407 / 0.509 | printed at 2 dp on purpose; set the producer's fmt to f2 |
| Spearman p (L641, L1200) | 1.5×10⁻⁷ | 1e-7 (e0) | set the producer's fmt to e1 |
| PH-02 q (L636, L1195) | 0.029–0.046 | 0.046 (maximum only) | add the minimum (0.0289) to the producer |

### Other number_diff "new" values not printed in the tex, by design

| claim | tex | why |
|---|---|---|
| 953 / 1,060 (L517) | historical (pre-repair) | the new 278 / 243 are at L516 / L519 |
| 43,140 (L80, L674) | pre-registered / "sample of the time" | the new 43,174 is at L680 / L687 |
| 8,614 (L675), 114 (L679) | labelled "of the time" | the new 8,660 / 107 are at L680 |
| 0.066 / 0.068 (L653) | pre-registered run | |

**Recommendation to the producer.** For 4-dp sources, collect.py's float formatting is unreliable at x.xxx5. Fix
`new_fmt` for the 10 "producer wrong" rows above, and recompute FR-03 and the NSTI ≥0.5 precision from the cells.
Better still, store unrounded means: in `crossover.py`, `compare_arms.py` and `comparison_summary.json` (the
by-stratum and by-region blocks), and in `summary_family.json` / `summary_genus.json`.

---------------------------------------------------------------------------------------------------------------

## 4. Stale sweep

**Method.**
- For each value whose print changed, I searched the whole new tex for the old printed string as a numeric token.
  That gave 452 occurrences, 327 of them unchanged from the old file.
- I scored each survivor by context similarity to the value's old occurrence (left 6 and right 4 content words) and
  read every hit with a score of 2 or more.
- I also grepped phrase-level stale strings: "two of 22", "q = 0.002", 0.146, 755, "near 3", "slightly conservative",
  "58--64", 67,038, 19,686, "within 0.016", 0.783, 0.854, "Ac.NO and Etho", and others.

**Survivors in the same context, all labelled historical:**
- L80 43,140 ("of the pre-registered pass"; "against 49 in the pre-registered run");
- L517 953 / 1,060 ("Before the cache was repaired");
- L519 0.889 / 0.713 (explicit before → after);
- L653 0.068 ("the pre-registered run");
- L674–675 43,140 / 8,614 ("On the sample of the time");
- L679 114 (107 under A-C12; "On that sample").

**No unlabelled survivor in the same context.** The other matches are different quantities: q = 0.002 in sec:sharedref
(controlled-benchmark BH q), 3,755 deployment ASVs, 0.768 V1–V2 / NSTI recall, and so on.

**Values whose print was kept although the new unrounded value rounds differently:** B1, B2, B3 and C1 above.

**I re-verified the executor's hand-fixed repeats:**
- L402 59–65;
- L1163 4.6 / 22.9 / 33.3;
- L523 67,713;
- L899 0.250 / 0.304 / 0.612;
- L907 2,441;
- L1072;
- L629 0.624 / 0.710;
- L527 0.722;
- L533 / L1189 0.871 / 0.768.

All are correct.

**Limitation.** The sweep can only find old values that number_diff lists. A repeat of an un-inventoried number would
not be seen.

---------------------------------------------------------------------------------------------------------------

## 5. Internal arithmetic, printed next to its operands

**Consistent:**
- L527 / L1197 / L1226: 0.723 − 0.729 = −0.006; 0.734 − 0.729 = +0.005; −0.120; the −0.092 / −0.094 assignment
  (legacy / exact-tie) agrees with L76;
- L254: 0.102 − 0.154 = −0.052; 0.928 − 0.878 = +0.050; 48% = 0.025/0.052;
- L1072: 1 − 0.333 = 0.667, 1 − 0.229 = 0.771, 1 − 0.046 = 0.954; the largest residual, 0.023, is "within 2.3 points";
- L1068: 42–55-fold = 1/0.0239 … 1/0.0181;
- L1069: 0.6 = 2.4 − 1.8; 0.5 pooled;
- L1095: +0.007 = 0.6889 − 0.6823;
- L897: +19% = 0.372/0.313;
- L1191: 9.5 + 0.4 ≈ 10%;
- L577: 548 + 147 = 695 = 10.8% of 6,423;
- L519: 67,713 − 67,038 = 675;
- tab:magvs16s and tab:magval Δ against the displayed columns: within 0.001, as the captions allow.

**Inconsistent as printed (cosmetic):**
- L901: 0.361 against 0.612 − 0.250 = 0.362 (C2);
- tab:threepanel L1143: −0.027 against 0.783 − 0.809 = −0.026 (C3). This is the table's existing convention, also seen
  at L1126, L1130, L1138 and L1139.
