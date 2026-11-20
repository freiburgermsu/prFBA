# Run 5, first half (R5.1, R5.2, R5.4, R5.5 and their commits), 2026-11-16

Executor: Run 5A agent. Plan: `plan/PLAN.md` §0-5 (Run 5), §6, §8; handoffs from RUN1, RUN2_LOCAL, RUN2_TF,
RUN3, RUN4_LOCAL and RUN4_TF reports. The manuscript (`~/Documents/prFBA_paper`) and the decks (`presentations/`)
were not touched: they are the next agent's (R5.3, R5.7).

**Outcome.**
- The cache is unchanged: `stamp_state.py --compare stamps/R4L.final.json` reads CACHE UNCHANGED, catalogue current
  (`stamps/R5A.0.json`). The d20 baseline failed on exactly the plan's 17 files before the re-stamp.
- **R5.1.** W6 gained six recipe groups covering every handoff gap; on `--state pre` 214 of 222 numbers reproduce
  the printed value, 6 have no printed value, and 2 do not reproduce (documented: the nn1 EC F1 0.781 and its
  p = 0.64 have no surviving recipe). W7 now covers every claim: 178 items, 1,568 values, 0 errors, 0 pending.
  **691 values in 116 items change at print precision; 24 items change meaning** (listed below with tex lines).
- **R5.2.** Every manuscript figure was re-rendered into `manuscript_figures/out/`. All 8 controls (pre-repair inputs
  from the start tag) are pixel-identical to the paper's PDFs. Live: F4b, F13, F14, F15 change; F1, F2, F11, F16 are
  pixel-identical. F15 had a hard-coded, now-wrong annotation: fixed (computed), re-controlled, re-rendered. The 34
  content-changed depth-20 tree figures are committed.
- **R5.4.** Results docs re-numbered or marked historical, each with a dated re-score note (list below).
- **R5.5.** `d20_baseline.sha256` re-stamped on purpose (17 lines, old -> new in the commit). `check_d100.py`: every
  check passes (checks 1-4 on TF, the rest here; see "Where each step ran").
- Commits on `hardware-aware-alignment-dispatch`, path-limited, not pushed (list below).

---------------------------------------------------------------------------------------------------------------

## Where each step ran (memory policy)

Mid-run the orchestrator capped this run's footprint on FM at about 10 GB (the user's salterns rollout needs the
headroom). Every step after that point respected it; the steps before it were small anyway.

| step | where | peak | note |
|---|---|---|---|
| stamp compare (R5A.0) | FM unit `rescore-r5a-stamp0` | 0.8 GB RSS | before the cap |
| W6 `--state pre` re-check, new groups, full pre+post | FM units `rescore-r5a-w6pre0/1b/2`, `rescore-r5a-w6full` | RSS < 1 GB; cgroup peak (page cache) 10.6-19.4 GB | before the cap; the cgroup figure is mostly reclaimable page cache from reading the cell tables |
| W7 `collect.py` (all items, several times) | FM, plain process | 1.2 GB | |
| figure renders and pixel controls | FM (subagent, then F15 by hand) | small | kept here, as the controls compare against PDFs rendered here |
| docs | FM (two subagents + me) | < 3 GB per process | agents told the cap mid-run |
| `check_d100.py` checks 5 (post, six trees), copy (three legacy trees), 6 (baseline) | FM unit `rescore-r5a-check-local`, MemoryMax 4G | RSS 40 MB | |
| `check_d100.py` checks 1-4 (align, three union trees) | **TF**, `tf_unit r5a-check-align` inside `prfba-ns` | 10.7 GB | loads the 3.5 GB hit records; inputs pushed with `tf/lists/push_R5.5_check_align.txt` (14.0 GB, 1.79 GB on the wire), VERIFY OK 22/22; TF code synced to FM HEAD bb41947 first |

| final stamp (R5A.final) | FM unit `rescore-r5a-stampfinal`, MemoryMax 4G | 0.8 GB RSS | |

No job needed more than 10 GB on FM after the cap.

**Second policy update (all Run 5a computation on TF).** It arrived after every computation above had finished:
the W6/W7 runs, the figure renders and controls, the check_d100 runs and the final stamp had all completed, and
the work was committed. After it, only git operations (the last commit, `tf_code_sync.sh`) and text edits to this
report ran on FM. Nothing was re-run on TF to replace the finished FM steps. The figures were rendered on FM
before the update, so their pixel controls compare FM renders against the paper's PDFs, which were also rendered
on FM; no TF-render control was attempted.

---------------------------------------------------------------------------------------------------------------

## R5.1 Derived numbers (W6) and the number diff (W7)

### W6 additions (`manuscript_numbers/derived_numbers.py`)

Every new recipe was run on `--state pre` (pre-repair inputs: `$W/pre_repair` for moved-aside git-ignored files,
the start-tag worktree for tracked ones) before it was used.

| group | numbers (printed) | pre control |
|---|---|---|
| `panels` (extended) | L742 emptied-record rate pooled over 8 regions 8.3% / 29.4% and FL 4.7% / 22.9% (definition of `three_panel/emptied_strata.py`) | reproduces all four |
| `identity_bins` (extended) | deck `BAND_MAG` on the depth-100 bins, as first-class `deck_BAND_MAG_d100_{0,1}` (0.081 / 0.198 pre; 0.081 / 0.200 post) beside the depth-20 ones the deck printed (0.082 / 0.202 -> 0.083 / 0.205) | d20 values reproduce; d100 values have no print |
| `phenotype` (new; port of `tools/r4_pheno_manual_numbers.py`) | 755, 747 V4, EC-F1 delta -0.028, 54 media / 22 with signal, Spearman -0.14 and p 2e-4, McNemar media at q<=0.05 (count 2, direction, names Ac.NO,Etho.O2, max q 0.002), the nn1 EC-F1 pair and its test, and the A-C13 counts (43,140 rows, 0 empty-medium ATP, 8,614 rows, 0 gated, 2,362 models / growing, 114 reversals, 200 genomes) | all reproduce except two (below) |
| `mag_depth500` (new; port of `tools/mag_depth500.py`) | 49,185 / 88 / 87 / 1 / 0 | reproduces |
| `sr_masks` (new) | legacy forced run under the random-abstention mask 0.668, like-for-like +0.019 (genus) / +0.016 (family), exact-tie vs its own mask +0.019 / +0.015, exact-tie forced gain +0.061, P(n_genomes>=2) by novelty bin 0.100 -> 0.180 and "rises monotonically" | reproduces |
| `as_selection` (new) | AS-17 603 / 248 / all below 0.987 / median 0.956; AS-18 145 / 692 / 10.8%; AS-19 43%; AS-23 median selection size 1 | reproduces |

**Not reproducible: "EC F1 0.782 against 0.781 (p = 0.64)"** (PH-04, disc.atp_vs_nn1, L1186). The prFBA mean
reproduces (0.782). No recipe tried gives nn1 0.781 or p 0.64: cells-table F1 for both arms (nn1 0.7817 -> prints
0.782), per_genome `ec_f1` for prFBA, the pre-registered `phenotype/per_genome.jsonl` and its `.bak-1015`, rxn space,
f1_accessory, jaccard, mcc; Wilcoxon (wilcox/pratt/zsplit 0.45/0.44/0.52), paired t 0.67, sign test 0.49. The recipe
keeps the cells-table F1 and the Wilcoxon signed-rank test the phenotype contrasts use, and records both as
`KNOWN_PRE_MISMATCH`. New values: 0.782 against 0.780, p = 0.84. **The manuscript agent should replace the printed p
with the recipe's (or drop it)**; "the two tie on ECs" still holds.

Full runs: `--state pre` 222 numbers (214 match, 2 documented mismatches, 6 unprinted); `--state post` written. The
40 pre-state inputs read from the live tree are the skani hits and manifests of the MAG arm, which the re-score does
not write.

### W7 (`manuscript_numbers/collect.py`, `claims_spec.json`)

- `claims_spec.json` was amended by the work-dir tool `$W/tools/r5a_spec_extend.py` (idempotent; the pre-amendment
  spec is `$W/control/r5a/claims_spec.before_r5a.json`). Every remaining manual part is now a W6 value, a JSON
  selector or expression, a hand-carried value with its source (`manual_values`: Stage O/E and calibration values the
  R4.3 agent checked by hand), or an explicit resolution (`manual_resolution`: unchanged by construction, historical,
  figure, or not recomputed, with the reason). 0 manual parts are left unresolved.
- `collect.py` adds automatic meaning flags (sign flip; a p or q crossing 0.05 or 0.01; a categorical or yes/no
  change; a CI's excludes-zero flip; a printed bound that no longer holds), item-level counts, and renders the
  manual values and resolutions. Tests: 43 pass (`tests/test_collect_numbers.py`, `tests/test_derived_numbers.py`).

**number_diff summary** (`manuscript_numbers/number_diff.{json,md}`): 178 items, 1,568 values; old resolved 1,568,
new resolved 1,568; 0 errors, 0 pending; every old value equals the print except 6 (the 2 known nn1 values in two
items each, and the known "about 45%" rounding in two items); 1,120 values changed, **691 at print precision, in
116 items**.

| stage | items | items changed at print | values changed at print | values |
|---|---|---|---|---|
| R1.3U union gene capture | 11 | 8 | 22 | 52 |
| R1.3L legacy | 6 | 5 | 7 | 20 |
| R1.4 d20 | 1 | 1 | 1 | 4 |
| R1.5 E1 0.995 | 3 | 3 | 13 | 25 |
| R1.K deployments | 7 | 0 | 0 | 38 |
| R2.1 E3 gate | 3 | 3 | 26 | 30 |
| R2.2 MAG | 4 | 3 | 13 | 29 |
| R2.3 cross-panel / MAG pairing | 42 | 32 | 407 | 584 |
| R2.4 d100 as shipped | 22 | 19 | 70 | 132 |
| R2.5 d20 as shipped | 11 | 9 | 22 | 63 |
| R3 shared reference genus | 23 | 6 | 13 | 313 |
| R3F family | 4 | 2 | 5 | 63 |
| R4P phenotype | 10 | 9 | 43 | 59 |
| R4S2 pass 2 | 4 | 3 | 17 | 37 |
| R4H post hoc | 11 | 4 | 4 | 56 |
| R5A W6 | 9 | 7 | 26 | 44 |
| D1 (Path A) | 6 | 2 | 2 | 19 (genus n_eff 6,716 -> 6,715) |
| D7 | 1 | 0 | 0 | 0 |

### Claims whose meaning changes (24 items; the manuscript agent must rewrite these, not just refresh digits)

Each is at the top of `number_diff.md` with its note and flags.

1. **Phenotype vs nn1 no longer significant** — PH-04 (L635-637), disc.atp_vs_nn1 (L1186): -0.146 media (p 0.017,
   q 0.028) -> -0.030 (p 0.41, q 0.51); all-genome p 0.28 -> 0.46 with the sign flipped. EC-F1 p = 0.64 is unsourced
   (above).
2. **McNemar: four media, not two** — PH-02 (L79, L631, L639, L1185), disc.atp_vs_picrust2 (L1185): Ac.NO q 0.029,
   Ac.O2 0.040, Etho.O2 0.046, For.O2 0.046 (all favour PICRUSt2); "q = 0.002" no longer holds.
3. **Phenotype abstract and sample** — abs-phenotype-readout (L79), PH-01 (L79, 632, 1193, 1222): 755 -> 759
   (747 -> 750 V4), 0.13 -> 0.07 media (p 0.49 -> 0.79), EC-F1 delta -0.028 -> -0.031; the stratified sample was
   re-drawn (442 of 1,000 ATP and 62 of 200 gapfilled genomes kept): disclose why. PH-03 (L634, L1190): rho -0.14
   (p 2e-4) -> -0.19 (p 1.5e-7).
4. **A-C13 counts are a re-run on a new sample** — abs-phenotype-controls (L80), AM-02 (L667-679): 43,140 -> 43,174
   rows, 0 of 8,614 -> 0 of 8,660 gated, 2,362 -> 2,372 growing models, 114 -> 107 reversals; "reproduces the
   pre-registered pass exactly (43,140 of 43,140)" describes the 10-29 sample. PH-06 (L641-646): growth log-ratio
   prFBA 0.008 -> 0.020 vs nn1 0.011, so "0.008 against 0.008" no longer reads as identical; burden 80-84 vs 82.5 ->
   78-82 vs 79.5.
5. **Depth probe now at or above the answered set in EC space** — abs-depth-ec (L67), depth_probe_quality (L728-731),
   lim.depth_probe_quality (L1204): recovered EC F1 0.783 -> 0.829 (vs 0.798), precision 0.747 -> 0.791 (vs 0.790),
   recall 0.854 -> 0.904. "score below the depth-20 answered set ... slightly lower" must be rewritten.
6. **Surviving-uncovered records** — emptied_vs_uncovered_panels (L742), surviving_uncovered (L1057-1060),
   surviving_pred_empty (L1064), tab_emptied_uncovered (L1066-1081): FL 3.4/2.8/2.3% -> 2.4/2.3/1.8%, pooled
   3.1/2.4/1.8% -> 1.8/1.8/1.3% ("near 3%" -> near 2%), stratum range 1.5-3.7% -> 0.0-3.0%, 29-43-fold -> 42-55-fold,
   spreads 1.1/1.3 -> 0.6/0.5 points; **CMH original vs random p 0.004 -> 0.117 (no longer significant)**, vs k-means
   0.15 -> 0.99; gene-less survivors 107/76/62 -> 36/42/31 (original no longer the most). L742 also: original FL
   emptied rate 4.7% -> 4.6%, pooled 8.3% -> 8.2%.
7. **"Random within 0.016 at every MAG completeness" no longer holds** — condF1_crude_bounds (L1018-1019): max gap
   0.0167.
8. **MAG multiplicity strata** — twopanel_mag_abstention_strata (L1085-1087): k-means without the >=20 stratum 0.686
   vs original 0.682 (was 0.687 vs 0.688).
9. **Cache-repair disclosure** — AS-05 (L517): must become a past-tense correction note (953 / 1,060 no longer
   describe the reported state; truth_empty 278, pred_empty 243).
10. **Gene-capture denominators** — genecap_denominators (L397, L419): score-zero 11.1% -> 9.9%, pred_empty 1.6% ->
    0.4% (0.987), 1.4% -> 0.4% (0.995); scoreable n 67,038 -> 67,713.
11. **As-shipped deficit magnitude** — disc.margin_reversal (L1184-1188, 1215-1216): -0.016 -> -0.006 (prFBA 0.713 ->
    0.723 vs PICRUSt2 0.729); the per-region ranking is unchanged (prFBA ahead only for FL, V1-V2, V1-V3), so the
    direction holds; "only narrowly" becomes narrower.
12. **Three-panel CI flips** — std_gap_share (L1095-1096), taxonomy_bootstrap_CIs (L1094-1119): 21 excludes-zero
    verdicts flipped in R2.3 (`control/r2/r23/three_panel_ci_flips.json`), all borderline. **None of the 21 cells is
    cited by any extractor** (checked against every three_panel selector); "excludes zero in all 16" still holds
    (16 of 16). Blanket sentences over cells (L997-999, L1094-1137) should be re-read.

Notable digit-only changes the agent should not miss: the abstract's as-shipped no-EC share 11.1% -> 9.8% (L74),
d20 population F1 0.600 -> 0.609, E3 gate numbers (L254-255), MAG resolve 0.308 -> 0.313 and the tab:magval /
tab:magvs16s columns, the R4.3 values in the R4_TF report (0.0282 -> 0.0283, 0.0110 [0.0031, 0.0222] -> 0.0111
[0.0031, 0.0223], 0.0411 -> 0.0412, +0.0096 -> +0.0097), the wobj values 0.0120 -> 0.0121 and 27.9 -> 27.8.

---------------------------------------------------------------------------------------------------------------

## R5.2 Figures

Details: `$W/control/r5a/figures/FIGURES_R52.md`, `figures_r52.json` (subagent), plus the F15 fix below.

**Controls** (renders from the start tag's inputs, scripts copied with every data path rewritten to `$W/start`;
strace showed no live-repo read; 150 dpi): F1, F2, F4b, F11, F13, F14, F15, F16 all **pixel-identical** to the
paper's PDFs.

**Live renders** (`manuscript_figures/out/`, git-ignored; not copied into the paper):

| paper file | out/ file | vs paper | PDF md5 |
|---|---|---|---|
| F1_accuracy_heatmap.pdf | F1_accuracy_heatmap.pdf | identical | 6ec9ee20a6cbc3741ace7007ca1b906e |
| F11_self_inclusion_exclusion_delta.pdf | same name | identical | de882a728030fbd30a21723f8bcd722c |
| F2_perrank_concordance.pdf | F2_perrank_concordance_exclude_self_ALT.pdf | identical | d7f167d524843fcc3595b3afadc0803a |
| F16_panel_space.pdf | two_panel/F_panel_space.pdf (committed copy kept) | identical | (unchanged file) |
| F4b_percent_genecapture_by_region.pdf | same name | 9,963 px differ: bars +0.7-1.1 points, every n= label | eee40cef1c3ea6bb4fd6afee16e2302c |
| F13_prfba_vs_picrust2.pdf | same name | 7,206 px: suptitle n 67,038 -> 67,713, panel B labels +0.01, panel C V3-V4 / V4-V5_944R rows swap (tie); "+0.07 to +0.10" unchanged | f7dfa4a58721f09612fb099c481dc2c6 |
| F14_mag_completeness.pdf | F_mag_completeness.pdf | 8,732 px: "+0.103 / +0.212" -> "+0.102 / +0.209", footer n 1,967-2,441 | 2f9b18cccfc2596ac1b69a7931f75ef3 |
| F15_mag_crossover.pdf | F_mag_crossover.pdf | changed: both panels, "(+0.102 at 100%)", footer n 1,645-2,448, panel-B note "0.313 -> 0.258" | 03ca35710ed5e1ea807385b42d9bf50b |

- **F15 fix.** `mag/plot_crossover_d100.py` hard-coded "0.308 -> 0.253". It is now computed from
  `crossover_d100.json` (round half up from the stored 4 dp, as the hand value was: 0.3075 -> 0.308). The control
  re-run with the edited script is pixel-identical to the paper's F15; the live figure reads "0.313 -> 0.258".
- **F14 hard-coded text still true:** "p<10^{-164}" (live p 7e-170, still below the bound), "+0.05 easier for the
  16S arm" (0.053 -> 0.055, prints +0.05), "0.56 vs 0.31", "31%".
- `fig_pairwise_hist.png`: not re-rendered (input is the 16S embedding, unchanged); its md5 equals the paper's.
- `two_panel/F_panel_space.pdf` was rewritten by the render with only its creation date changed; restored.
- **d20 tree figures** (`region_validation/figures/`, 44 tracked): all re-rendered with `--force`; 34 files changed
  pixels (F1, F2, F3, F4, F4b, F5b, F6, F7 for 7 regions, F8, F9, F11, PDF+PNG) and are committed; 5 PDFs whose only
  change was the creation date were restored (their PNGs are byte-identical). No pixel control is possible (no
  reference rendered from the June data). F5 and F12 still render as "(no data)" placeholders, as before.
- The six d100 trees' git-ignored `figures/` were re-rendered in Run 1 after their last data change.

---------------------------------------------------------------------------------------------------------------

## R5.4 Docs

Every doc below carries a dated re-score note. Two doc subagents did the MAG/two-panel group and the
PICRUSt2/16S/top-level group (each value taken from its JSON path or recomputed from the cell tables, after checking
that the same recipe on the pre-repair artefacts gives the printed value; helpers and patch scripts in
`$W/control/r5a/docs/`). I did the shared-reference docs and the historical banners.

| doc | what changed |
|---|---|
| `region_validation/shared_reference/WEIGHTED_OBJECTIVE_RESULTS.md` | re-derived from the new `recomputed_contrasts.json` / `summary.json`; every replacement keyed to its source cell (`mean [lo, hi]` triples matched to exactly one contrast x measure x stratum cell; singles checked old == printed). −0.0120 → −0.0121, 27.9 → 27.8, the n3 column, headroom 0.243 → 0.245, top-30 share 73% → 72%, vertex rows 203 → 202 / 269 → 268, n3 q 7.8e-5 → 2.8e-5, run time/memory row. Verdicts unchanged |
| `C12_SCOPE.md` | new section 7b: the A-C13 configuration re-run (43,174 ATP rows; 0 of 8,660 gated; 2,372 growing; 107 reversals; re-drawn sample; new controlled medians); top note marks the 10-27/29 body historical |
| `SPEC.md` | a dated factual note at the end of section 15 with the same counts; pre-registered and amendment text unchanged; notes that section 8's 67,038 is a pre-registration-time figure |
| `PFUNC_ESTIMATORS_RESULTS.md`, `OBJECTIVE_TUNING_RESULTS.md`, `data/pfunc_estimators/README.md` (+ the 11-03 corrections section, present), both `verification/VERIFICATION.md`, `data/calibration/README.md`, `PFUNC_E3_GATE.md` | checked only (re-numbered by the R2.1 / R4.3 TF agents); the calibration README is script-regenerated and dated 2026-11-16 in its own header |
| `mag_validation/out/RESULTS_v2.md` (~115 values), `RESULTS_tau99.md` (~39), `DEPTH_PROBE.md` (checked, unchanged) | resolve 0.3075 → 0.3134, tau-99 0.366 → 0.372, paired n/deltas, identity bins, AF sweep on one scoreable basis (0.314 ... 0.656); "9,745 of 9,746" → "9,746 of 9,746" (33953.3 repaired) |
| `mag_validation/out/RESULTS_v1_superseded.md`, `CORRECTIONS.md` | historical notes |
| `two_panel/RESULTS.md` (~84), `TWO_PANEL_RUNBOOK.md` (4) | k-means coverage edge +3.4 → +2.8 points crude; surviving-record CMH verdict no longer significant (p 0.002 → 0.71 in its d20 framing); gene-less sources 161 vs 45 → 38 vs 36 ("three times as hard" → "about equally"); **new section 11: the D6 contig-cap caveat** (random-panel genomes 1733.9197 and 1733.9535 at exactly 25,000 contigs; no random-panel note existed elsewhere) |
| `region_validation/picrust2_comparison/README.md` (~145) | n 67,038 → 67,713, exclude-self coverage 0.751 → 0.763, population 0.600 → 0.609, ceiling 2,689 → 2,719, depth probe 0.783 → 0.829 ("of ordinary quality" → "at least of ordinary quality"), MAG-arm tables; "p_holm <= 5e-14" → "<= 3e-7" (all still significant); D10 notes kept |
| `region_validation/16S_pipeline_method_paper.md` (~65) | sections 3.5-3.8, Tables 5, 5b, 8 |
| `region_validation/SELF_RECOVERY_ANALYSIS.md`, `region_validation/README.md` | d20 self-recovery failures 11,300 → 10,821 (17.1% → 16.4%), no_gene_set 974 → 262 |
| `README.md` | E3 gate block (19,158 → 19,542; Brier −0.0530 → −0.0520; ~45% → ~48%) with a dated note |
| `region_validation/data/{variant_comparison,p3_vs_old_comparison}.md` | dated "historical" notes (June studies) |
| `genecap_repaired_impact/impact_summary.md` | one dated line: the committed files now equal its "new" pass |
| historical banners: `UNION_THROUGHOUT_CHANGES.md`, `UNION_VERIFICATION.md`, `D100_VS_D20.md`, `JOBS_AUDIT_2026-11-02.md`, `EGC_C12_REPORT.md`, `C12_FIX_PLAN.md`, `PR_ModelSEEDpy_gapfill_reversals.md`, `PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md`, `MANUSCRIPT_VERIFICATION_7dc189e.md`, `impact_tables.md` | "the numbers in this document are historical", with the headline old → new moves and a pointer to `number_diff.md`; bodies untouched (UNION_THROUGHOUT_CHANGES alone is 3,100 lines of edit records) |
| historical banners (second sweep) | `PR_ModelSEEDTemplates_{directions,tierB,penalty}.md`, `ISSUE_ModelSEED_template_direction_audit.md` (the 218-model template study of 10-27/29; D9), commit `9c67cc8` |
| not changed | `METHODS.md`, `region_validation/DESIGN.md`, `RUNBOOK_d100.md` (no re-scored numbers; RUNBOOK's EMPTIED 14.3/40.4/49.9 are depth-20 alignment counts), deployment READMEs/RESULTS (D4: deployment numbers unchanged; the KBase lane updated `deployment_d100/README.md` section 11) |

Sweep after all edits, with the plan's pattern and a broader one (A-C12/A-C13 counts, 19,158, 2,689, 0.783, 0.3075,
1,916, 0.6883, 755, 0.146, ...) over every tracked .md: every remaining hit is a correct new value, a pre-registration-time figure
(SPEC section 8), inside a dated old → new note, or in a doc carrying a historical banner.

Pre-existing errors the doc agents found and did not fix (not re-score matters; flag to the user):
- method paper Table 6 prints V7-V9 Δ +0.193, but 0.819 − 0.653 = +0.167 in both old and new data;
- "~112k BV-BRC genomes with cached PGFams" (PICRUSt2 README caveat 2, method paper §3.7) is a fixed string; the
  cache holds 159,753 entries;
- `RESULTS_v2.md`'s "73x" GTDB fold disagrees with `make_deck.py`'s 14x.

---------------------------------------------------------------------------------------------------------------

## R5.5 Baseline and check_d100

- Before: `sha256sum -c` failed on exactly the plan's 17 files (`control/r5a/d20_baseline_failing.txt` ==
  `d20_baseline_expected_fail.txt`); the six hits files, the d20 union selections and taxacc,
  `space_deciles_three.json` and `three_panel_mag_gtdb.json` passed.
- Re-stamped in the same order and format (`sha256sum` of the file's own path list); old -> new in the commit
  message and `control/r5a/d20_baseline_restamp.txt`. After: 32/32 OK.
- `check_d100.py`: checks 1-4 (align) PASS for all three union trees on TF (`logs/r5a/check_align_TF.log`); copy
  checks PASS for the three legacy trees; check 5 (post) PASS for all six trees; check 6 (baseline) PASS
  (`logs/r5a/check_local.log`). This is every check `check_d100.py all` runs.

---------------------------------------------------------------------------------------------------------------

## Commits (prFBA, `hardware-aware-alignment-dispatch`, path-limited, not pushed)

| commit | stage | content |
|---|---|---|
| `9148c20` | R5.1 | W6 recipe groups, W7 meaning flags and manual rendering, amended `claims_spec.json`, `derived_numbers.{pre,post}.json`, `number_diff.{json,md}`, tests |
| `c456599` | R5.2 | `manuscript_figures/mag/plot_crossover_d100.py` (computed F15 note); 34 re-rendered d20 tree figures (live manuscript-figure md5s in the message) |
| `f021092` | R5.5 | `region_validation_d100/d20_baseline.sha256` re-stamp (17 old → new hashes and the check_d100 results in the message) |
| `9e5a08b` | R5.4 | WEIGHTED_OBJECTIVE_RESULTS.md, C12_SCOPE.md 7b, SPEC.md §15 note, ten historical banners |
| `7ca982a` | R5.1 | curated meaning-change notes (depth probe, surviving-uncovered, multiplicity, PH-06, deficit magnitude) |
| `bf79d13` | R5.4 | the two doc agents' 15 files |
| `9c67cc8` | R5.4 | historical banners on the ModelSEED template PR / audit texts (second sweep) |

TF's clone was fast-forwarded to `9c67cc8` (= FM HEAD) with `tf_code_sync.sh` at the end; its tracked tree was clean.

Messages: `logs/commit_msgs/r5a_*.txt`. Final stamp `stamps/R5A.final.json`: CACHE UNCHANGED vs R5A.0, catalogue
current, d20 baseline 32/32 OK. Untouched: `~/Documents/prFBA_paper`, `presentations/` (including the user-edits
deck), ModelSEED checkouts, codiffusion_bioreactor, EmilyKin, KBase.

New work-dir tools: `tools/r5a_spec_extend.py` (the claims_spec amendment), `tools/r5a_wobj_doc.py` (WEIGHTED doc
triples), `tools/doc_changed_leaves.py`; TF push list `tf/lists/push_R5.5_check_align.txt`.

---------------------------------------------------------------------------------------------------------------

## Handoff checklist (every "for Run 5" item of the run reports)

| item (source) | status |
|---|---|
| stamp_state --compare latest stamp | done: CACHE UNCHANGED (R5A.0, and R5A.final at the end) |
| W6 gap: L742 rates 8.3% / 29.4% have no recipe (task) | done: `panels_emptied_pooled_*`, reproduces; new 8.2% / 29.4% |
| W6 gap: deck BAND_MAG matches only d20 bins (task) | done: d20 kept (0.082/0.202 → 0.083/0.205); d100 variant added (0.081/0.200 post); the deck agent chooses |
| W6 gap: "EC F1 0.782 vs 0.781 (p=0.64)" needs nn1's per-cell F1 (RUN4_LOCAL 5.2) | recipe written; **0.781 and 0.64 not reproducible**, documented; new 0.782 vs 0.780, p 0.84 |
| W6 gap: phenotype manual-numbers recipe (`control/r4/pheno_manual_numbers_atp.json`) (RUN4_LOCAL) | done: ported into W6 `phenotype` group; reproduces every old print |
| W7 over ALL items, old / new / source / tex lines, meaning changes marked (task) | done: 178 items, 0 errors, 0 pending; 24 meaning-change items |
| collect.py `manual` parts (RUN4_TF: pfunc_calibration / stage_e still manual) | done: manual_values with source (R4.3 hand-checked values); every manual part resolved or given a reason |
| SPEC §15 and C12_SCOPE counts 43,140 → 43,174, 8,614 → 8,660, 2,362 → 2,372 (RUN4_LOCAL 5.3) | done (also 114 → 107 reversals) |
| WEIGHTED_OBJECTIVE_RESULTS.md re-derived from the new summary (RUN4_LOCAL 5.4) | done |
| PICRUSt2 READMEs quoting n = 67,038 (RUN2_LOCAL 4) | done (P2C README; no README exists under the d100 picrust2_comparison) |
| `mag_validation/out/RESULTS_*.md` not regenerated (RUN2_LOCAL 4, 6) | done |
| 44 tracked d20 tree figures (RUN1 7) | done: re-rendered; 34 content-changed files committed, 5 date-only PDFs restored, 5 PNGs byte-identical |
| June d20 variant-selection files (RUN1 7) | left as a historical variant study; their docs now say so |
| d20_baseline.sha256 re-stamp + check_d100 (RUN1 6, RUN2_LOCAL 3, R5.5) | done: failing set == plan's 17; re-stamped; all checks PASS |
| 21 three-panel CI verdict flips (`control/r2/r23/three_panel_ci_flips.json`) (RUN2_LOCAL R2.3) | done: none is cited by any extractor; listed in number_diff's meaning section for blanket sentences |
| Stage E README corrections section (RUN4_TF 5) | checked: present (re-appended by the TF agent) |
| three PFUNC/OBJECTIVE docs re-numbered by the TF agent (task) | checked, not redone |
| `genecap_repaired_impact/impact_summary.md` one line (plan §5) | done |
| D6 contig-cap caveat in random-panel notes (plan R5.4) | done: `two_panel/RESULTS.md` §11 |
| figures: F13 W4 renderer, MAG plots, F4b, F16 (task, plan R5.2) | done with pixel controls (8/8 identical); F15 hard-coded note fixed |
| `data/phenotype_summary.json` regenerated (RUN3 3) | done in Run 4 (byte copy of the gapfilled summary); nothing to do |
| TF pull hazard: `pull_R4.3_posthoc.txt` lists `weighted_objective/` (RUN4_LOCAL 6) | not touched; still applies to anyone installing that pull |
| deck constants / KBase text (R5.3) | not in scope (next agent); inputs listed below |
| disk cleanup suggestions (RUN1 9, RUN2_LOCAL 7, RUN3 7, RUN4_*) | not done (needs the user's go-ahead). New: TF holds 14 GB of pushed hit files under its shadow tree |

---------------------------------------------------------------------------------------------------------------

## What the manuscript agent (R5.3 / R5.7) must know

1. **Start from `manuscript_numbers/number_diff.md`.** Its top section lists the 24 claims whose meaning changes;
   the per-stage tables give every value old → new with source and tex lines; "(manual)" rows give resolutions.
2. **Rewrite, not refresh:** the phenotype block (nn1 significance lost; four McNemar media; re-drawn sample; the
   A-C13 counts are a re-run), the depth-probe quality sentences (recovered answers now at or above the answered
   set in EC space), the surviving-uncovered analysis (rates ~2%, CMH no longer significant), L517 (past-tense
   correction note), L397/L419 denominators, L1018-1019 ("within 0.016" → 0.017).
3. **Unsourced printed values:** "EC F1 0.782 against 0.781, p = 0.64" (L1186) has no recipe; use the recipe's
   new values (0.782 vs 0.780, Wilcoxon p = 0.84) or drop the p. "about 45%" (E3 recalibration share) was a
   rounded-input figure; new ~48%.
4. **Not recomputed** (manual_resolution "not recomputed"): AS-15's four rxn-vs-EC divergences (0.026 / 0.014 /
   0.05 / 0.03), the depth probe's 0.651 / 0.800 answered-set anchor accuracy, the ~0.96 / ~0.49 / ~0.75 probe
   projections (the P2C README now gives ~0.77 for the EC projection), the wobj "weighted term median 0.11%", SR-02's
   inner margins 0.001-0.003, the AF-gate "1,386 newly admitted incl. 19 gene-less" (the summary stores 1,391
   non-empty; share above the FL level 26.5%).
5. **Figures:** copy from `manuscript_figures/out/` (md5s above): F4b, F13, F14 (= `F_mag_completeness.pdf`), F15
   (= `F_mag_crossover.pdf`) changed; F1, F2, F11, F16 are pixel-identical (no copy needed). F13's panel-C label
   "+0.07 to +0.10" is computed and unchanged; F14's "p<10^{-164}" is still true (7e-170).
6. **Deck inputs** for R5.3 are all regenerated; the W6 deck constants are in `derived_numbers.post.json`
   (`deck_*`): AF_SWEEP 0.314/0.319/0.381/0.454/0.569/0.656, MAG_REACH 6011/20.3/45.1 (unchanged), BAND_MAG
   0.083/0.205 (d20 bins) or 0.081/0.200 (d100 bins), XPHYLUM_BENCH 6.4 (unchanged), EMPTIED from
   `emptied_strata_d100.json`.
7. **p-values:** W6 formats `pheno_spearman_p` with one digit (`e0`), so number_diff prints the new value as `1e-7`;
   the unrounded value is 1.47e-7 (the `new` field of the JSON; Run 4 printed 1.5e-7). Use the unrounded `new`
   field for any p-value.
8. **Plan errata still open** (plan §5): "460,117" → 460,116, "1 to 5,265 contigs", marker reach 6,011 (W6 still
   counts 6,011, not the plan's 6,010), F16 vs three-panel deciles.
