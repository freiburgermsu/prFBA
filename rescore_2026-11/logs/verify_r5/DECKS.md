# Verification of the R5.3 deck rebuild (commit `a6812be`), 2026-11-16

Verifier: an independent subagent. It reported only and edited, committed or checked out nothing in `~/Documents/prFBA`.
Scratch files are in `/tmp/claude-1000/-home-freiburger-Documents-prFBA/400ec395-f172-59da-a98d-3dbe6b3a4256/scratchpad/verify_r5/decks/`:
- extracted text: `txt/{pre,post,rebuild}_*.json|txt`
- rebuilds: `rebuild/`
- CI scripts: `ci/`
- deployment recount: `deploy/`
- figure regeneration: `figs/run/`
- renders: `render/`, `render/hi/`

Every job ran on this machine. The largest peaked at 1.1 GB RSS for under 20 s, under `systemd-run --user --scope -p MemoryMax=8G`.

## Verdict

**Nothing blocking.**
- The commit is path-clean.
- The committed decks are exactly what the committed generators produce from the current inputs.
- Every changed constant matches its source artefact.
- The legacy-reducer CI reproduces and is robust to seed.
- No stale value survives unlabelled.

There are four cosmetic or wording issues and two pre-existing PR-status inconsistencies (below). No slide overflows.

## 1. Commit scope

`git show --stat a6812be` changes nine paths, all intended:
- `make_deck.py` and `make_deck_condensed.py`;
- `deck_constants_d100.json`;
- `condensed/F_heatmap.png` and `condensed/figure_stats.json`;
- the four `*_2026-11.pptx` decks.

The `*_2026-10*.pptx` decks are untouched. For each of the four tracked ones, the md5 is identical at `a6812be~1`, at `HEAD` (= `a6812be`) and in the worktree. Their last commit is `817f20b`.
- `prFBA_brief_2026-10.pptx`: 960f1645…
- `prFBA_brief_panels_2026-10.pptx`: 061143d9…
- `prFBA_condensed_2026-10.pptx`: f69e5243…
- `prFBA_updates_2026-10.pptx`: 12db044f…

The untracked `prFBA_updates_2026-10.user-edits-1029.pptx` is also untouched:
- mtime 2026-10-29 11:26:52 −0500;
- ctime and birth 2026-10-29 14:34;
- md5 6e71dc53….

The worktree is clean apart from that untracked file.

## 2. Committed decks equal the committed generators' output

All four decks were rebuilt into scratch with the committed generators:
- `make_deck.py --out …`;
- `make_deck_condensed.py [--brief|--panels] --out …`.

Each took 17 s or less and about 0.5 GB. A Python audit hook logged every file opened.

**Result.** Every slide face, chart series and speaker-notes string matches the committed decks: 45/20/10/10 slides, 0 differences. Byte for byte, every zip part (slide XML, notes, charts, media) is identical. The only exception is the embedded chart `.xlsx` workbooks, which carry save timestamps.

The condensed figures were also regenerated, with a scratch copy of `make_condensed_figures.py` whose `REPO` points at the repo. All seven PNGs, `figure_stats.json` and `worked_example.json` are byte-identical to the committed ones.

**Inputs not in git** (118 files read; 102 tracked, none modified in the worktree):

| status | mtime | path |
|---|---|---|
| ignored | 11-15 22:50 | `region_validation_d100/picrust2_comparison/data/comparison_cells.csv.gz` (the bootstrap CIs) |
| ignored | 11-15 22:35 / 22:39 | `region_validation_d100/picrust2_comparison/data/prfba_selection{,_legacy}.json` |
| ignored | 11-15 22:58 | `region_validation/picrust2_comparison/data/prfba_selection.json` |
| ignored | 11-16 12:01 | `region_validation/shared_reference/data/calibration/scores.json` |
| ignored | 11-16 16:12 / 16:32 | `shared_reference/data/phenotype_{atp,gf}_c13/per_genome.jsonl` (re-run) |
| ignored | 10-18 19:00 / 19:29 | `shared_reference/data/phenotype_{atp,gf}/per_genome.jsonl` (pre-registered sample, for "kept" counts) |
| ignored | 11-15 16:49 / 11-04 18:38 | `deployment_d100/logs/kbase_upload_{bvbrc_d100,d100}.log` |
| outside repo | 10-29 | `codiffusion_bioreactor/` and `EmilyKin/bvbrc_alignment_hits_d100/asv_alignment_summary.csv` |
| outside repo | 10-29 to 11-03 | `ModelSEEDTemplates-draftloops/validation/{scan/scan_summary.json, gf_pass/compare_vs_c13.json, gf_pass/report_rows.json}` |

All of these predate the commit (21:13 on 11-16). The decks therefore cannot be rebuilt from a clone of the repo alone.

## 3. Stale-value sweep (all four decks, faces and notes)

The patterns were run on the pre-commit decks (`a6812be~1`) as well as the committed ones. A zero on the committed side means something only where the pattern hit before.
- Hits before, **0 after:**
  - `q = 0.028` (6 → 0), `2 of 22` (3 → 0), `q = 0.002` (3 → 0)
  - `0.713` (2 → 0), `0.308` (10 → 0)
  - `within resampling noise` (3 → 0), `about 45` (7 → 0)
  - `0.311`, `0.379`, `0.452`, `0.202`, `4.7%`, `33.4`, `−0.028`, `1,367`, `26.7`, `0.0096`, `−0.021`, `−0.007`, `+0.003`, `67,038`, `−0.15`, `−0.002`
  - `marks itself invalid` (1 → 0)
- **Absent from the decks even before the commit** (manuscript-only strings): `0.146`, `two media`, `near 3%`, `0.783`, `slightly conservative`, `43,140`, `within 0.016`.

The surviving matches fall into three groups.
- **Current values that only match the pattern:**
  - "45%" on slide 37's notes is the AF-gate coverage gain, 1,391/3,121 = 44.6%, computed live from `af_gate_counterfactual.json`. It is not the E3 share.
  - `0.316` (s20, PICRUSt2 precision), the `0.567`/`0.655` chart floats and `0.727` (current matched PICRUSt2 F1) are current values.
  - "755" occurs only inside "3,755".
- **Labelled historical:**
  - `6467`/`3757` appear only as the 2026-11-04 mixed-annotation GenomeSets ("keeps the earlier mixed-annotation version", "still point at the mixed-annotation k/n versions"): full deck s6 and s45, condensed s16 and s20, brief s9 and s10, panels s9 and s10.
  - `8,614` appears with "pre-repair sample" on the s29 face, and in the notes as "On the 2026-10-29 sample" (s27) or "a side study on the phenotype sample of 2026-10-29" (s29 and the brief decks' s9).
  - "row for row" is always tied to the 2026-10-29 sample.
- **Phenotype against nn1:**
  - Every deck's face reads "−0.03 vs the nearest copy · q = 0.51, n.s." in white, not the win colour.
  - The notes say "not significant, so the phenotype does not favour prFBA over that baseline either".
  - No win claim against nn1 survives.
  - The McNemar media in the notes (Ac.NO, Ac.O2, Etho.O2, For.O2; q 0.029–0.046) match the `prfba_exact_tie__vs__picrust2` entries in `phenotype_atp_c13/phenotype_summary.json`. The four media with q_bh ≤ 0.05 are exactly Ac.NO 0.0289, Ac.O2 0.0396, Etho.O2 0.0461 and For.O2 0.0461, and PICRUSt2 is ahead on all four (c > b).
- **CMH:** slide 42's notes say "do not differ significantly (CMH p = 0.117 …, 0.99 …)", which matches `emptied_strata_d100.json` (0.1167, 0.9948). No significant-CMH claim survives.

## 4. Changed constants against their sources

| constant | deck | source value | ok |
|---|---|---|---|
| AF_SWEEP | .314/.319/.381/.454/.569/.656 | `derived_numbers.post.json` deck_AF_SWEEP_0..5: .31429/.31851/.38146/.45446/.56913/.65629 | yes |
| BAND_MAG | (0.081, 0.200) | deck_BAND_MAG_d100_0/1: 0.08141, 0.20040 (the depth-20 bins now give 0.0826, 0.2045, as the comment says) | yes |
| MAG_REACH, XPHYLUM_BENCH | 6011/20.3/45.1, 6.4 | 6011, 20.285, 0.4505, 0.0637 | yes |
| EMPTIED | 4.6/22.9/33.3 | `emptied_strata_d100.json`: 336/7,268, 2,019/8,817, 2,927/8,780 (emptied_pct 4.6/22.9/33.3) | yes |
| PHENO_EC_GAP | −0.03133 | pheno_ecf1_delta_vs_picrust2 −0.031332, n = 759 | yes |
| AF2_NEW / AF2_BEAT_FL | 1,391 / 0.600 / 369 → 26.5% | `af_gate_counterfactual.json`: 1391, 0.600306, 369 at 0.7116 (369/1391 = 26.53%) | yes (read live) |
| PFE_WITHIN_UNION_F1 | all (0.0085, 0.0072, 0.0097); n≥2 (0.0249, 0.0228, 0.0265) | `p4_results.json`: 0.008494 [0.007182, 0.009674]; 0.024891 [0.022827, 0.026491] | yes (read live) |
| OTUNE_LOG_CI | (0.0011, 0.0053) | `out_c4_c5.json` log_vs_lin@est/X1: [0.001053, 0.005302] | yes (read live) |
| INCL_N_SELECTED | 2.88 / 1.32 | computed from `gene_capture_cells.csv`; legacy 1.3205 matches number_diff | yes (live) |
| GAP_CI_D100 | union −0.0110…−0.0015; legacy +0.0006…+0.0101 | reproduced, see §5 | yes |
| E3 recalibration share | 48% | (0.153604 − 0.128686) / 0.05197 = 0.4795; number_diff methods:e3_gate_recalibration_share 0.4795 | yes |
| KBASE_D100 | 261948/6468/1, 275568/3758/1, 0 failed, 2026-11-15; ref_1104 6467/3757 | `kbase_upload_bvbrc_d100.log` L267–269, L578–580; `kbase_upload_d100.log` L349, L765 (both checked at build time) | yes |
| D100_BUILT / ABSTAIN (s6) | 3,225 / 51 | `bvbrc_rebuild_summary.json` codiffusion counts 3225/51 (also the shipped build's) | yes |
| INTEGRATION_STATUS | "Manuscript: re-scored revision" | judgment (s45 face) | n/a |

Other live values checked:
- E3 gate verdict PASS, finished 2026-11-15T23:13.
- Survivor uncovered rates 2.4/2.3/1.8% (s42).

## 5. Legacy-reducer as-shipped CI

`gap_ci_d100()` was lifted verbatim into `ci/gapci.py`.
- **Population:** EC space with `truth_empty` dropped, abstention scored 0, all regions pooled.
- **Bootstrap:** 10,000 replicates resampling the 9,708 source genomes, with a ratio estimator.

The input `comparison_cells.csv.gz` has:
- 0 duplicate (amplicon, method) rows;
- 67,713 paired amplicons for both rules (`dropna` drops 0).

| seed | union gap [95% CI] | legacy gap [95% CI] |
|---|---|---|
| **20261116** (deck) | −0.00619 [−0.01098, −0.00146] | **+0.00537 [+0.00057, +0.01009]** |
| 1 | [−0.01108, −0.00144] | [+0.00048, +0.01015] |
| 2 | [−0.01092, −0.00142] | [+0.00062, +0.01010] |
| 3 | [−0.01097, −0.00143] | [+0.00060, +0.01016] |
| 12345 | [−0.01094, −0.00146] | [+0.00064, +0.01002] |
| 20261029 | [−0.01110, −0.00148] | [+0.00047, +0.01017] |

- **The deck seed reproduces the deck's numbers:** +0.0006 to +0.0101, and −0.0110 to −0.0015.
- **Seed robustness:** the legacy lower bound stays above 0 on all six seeds (0.00047–0.00064), so "(significant)" is not a seed artefact. It is marginal, though, and the notes rightly say "only just excludes zero".
- **Point gaps:** both match `comparison_summary.json` (0.7343 − 0.7290 and 0.7228 − 0.7290).
- **Pre-repair cells:** the same code on `prFBA_rescore_2026-11/pre_repair/.../comparison_cells.csv.gz` gives union [−0.0207, −0.0109] and legacy [−0.0072, +0.0027]. That matches the old hand-carried (−0.021, −0.011) and (−0.007, +0.003) at 3 dp. The generator comment's "+0.0023" is +0.0027 with this code and seed; the difference is in a code comment only.
- **The manuscript cites point values only.** `prFBA.tex` L528 gives "population F1 0.734"; L1197 and L1226 say "a +0.005 lead as shipped". No CI and no significance claim.
- **The CI is not persisted in any committed data artefact.** It exists only in:
  - the committed pptx: s31 face "(significant)", and the notes of full-deck s19, s30 and s31 and condensed s11;
  - the `a6812be` commit message;
  - the untracked `RUN5B_REPORT.md`.

  `comparison_summary.json` has no pooled-gap CI, and its `paired_tests_ec_f1` are per-region Wilcoxon tests. The tracked `D100_VS_D20.md` still carries the pre-repair "−0.002 [−0.007, +0.003] … no longer significant", but it has a "numbers are historical" banner and slide 30's notes point to it as pre-re-score. Reproducing the CI needs the git-ignored cells file plus the generator's seed (20261116) and replicate count (10,000), both recorded in `make_deck.py`.

## 6. Known inconsistencies

**PR status (pre-existing; the generator text is unchanged by `a6812be`).**

The stale wording appears in two places:
- `make_deck.py` L3615–3616 (`s06d_wobj` notes): "a fix is pushed on ModelSEEDpy branch fix-prfba-single-model (b6857a0, 09856f5), with no pull request opened yet." This appears in the full deck's slide 9 notes and the condensed deck's slide 15 notes.
- `make_deck.py` L4953 (`s06e_msprob` notes): "…pushed, with no pull request opened yet: …". This appears in the full deck's slide 10 notes.

Other text contradicts it:
- The full deck's slide 10 face says "Pull request #20 open for review".
- Slide 45 and the closing notes say "freiburgermsu/ModelSEEDpy pull request #20 (opened 2026-11-12)".

A further related wording, also pre-existing: the slide 29 face takeaway reads "Safe to adopt; branch pushed, upstream PR pending" (L5026). Slide 45 says ModelSEED/ModelSEEDTemplates "pull request #15 (opened 2026-11-12, stacked on #11)". "Pending" can be read as not yet opened.

`gh` confirms all three PRs are OPEN:
- freiburgermsu/ModelSEEDpy#20, created 2026-11-13T00:40Z (= 11-12 CST);
- ModelSEED/ModelSEEDTemplates#15;
- freiburgermsu/ModelSEEDpy#19 (slide 44's "fork PR #19 open").

**Slide 6 notes: 149/130/129.**
- **What they are:** `D100_SET_CHANGES`, from `deployment_d100/compare/compare_codif.json` `set_changes["d20_54d1150_p3on(KBase) -> d100"]` (`RESULTS.md` §3, column A). They compare per-ASV genome sets between the codiffusion 20-hit KBase build (54d1150) and the depth-100 build (the 10-29/30 mixed-annotation one) over the 3,224 ASVs both answer:
  - 149 have a different genome set;
  - 130 have a larger one;
  - 129 of those 130 are records whose 20th hit tied the best (`A_cut_records_that_grew`, of 144 such records).
- **Recomputed against the BV-BRC-only rebuild:** `compare_d20_d100.py`'s `set_changes()` was run on `codiffusion_bioreactor/synthetic_genomes_bvbrc_d100/asv_reference_selection_d100.json` (2026-11-13), with the same d20 selection and hits.
  - The result is **149 / 130 / 129, unchanged**, with 3,224 answered by both, 15 the same size, 4 shrank and 51 both abstain.
  - The same run on the shipped 10-29 selection reproduces 149/130/129.
- **Why it is unchanged:** the rebuild re-selected only two codiffusion ASVs (`717e65e6…`, `cb3f42ad…`; `bvbrc_rebuild_summary.json` selection_vs_shipped). Each is a same-size member swap inside a set that had already grown against depth 20, so the counts cannot move.
- **Conclusion:** the notes' numbers are current for the build now in KBase (261948/6468/1). Only their provenance (a comparison made before the rebuild) is undocumented. Built and abstained counts are also unchanged (3,225 / 51).

## 7. Rendering

All four decks were rendered with LibreOffice → PDF → PNG at 80 dpi. Every slide whose text changed in `a6812be` was viewed:
- full deck 5, 6, 8, 9, 11, 12, 14, 16, 18–22, 24, 25, 27, 29–31, 33, 35–37, 41–45;
- condensed 1 and 10–20;
- brief 1 and 7–10;
- brief_panels 1 and 6–10.

There is no overflow, overlap or clipped text. The only close fit is the "Manuscript: re-scored revision" chip on the condensed and brief closing slides, which wraps to two lines inside its box. Slides that changed only in notes (the s1 titles, full-deck s2, s15 and s23) have unchanged faces.

## Cosmetic and wording issues (none blocking)

1. **s31 face** (full deck) prints a bare "+0.005 (significant)" for legacy. The s19 and s30 notes say "only just excludes zero", and the lower bound is 0.0005–0.0006. A softer face qualifier would match the notes. The face also sits under the title "One rule throughout: the exact-tie union" while the table now shows legacy significantly ahead of PICRUSt2 and the union significantly behind. That tension is a data-driven judgment call, already flagged in RUN5B as Flag 1.
2. **s31 notes:** the inserted sentence "On the re-scored data the union also costs population EC F1 … legacy being the rule that leads PICRUSt2 there." now sits between "On 2026-10-27 the authors chose the union …" and "Until then the manuscript used each rule …". "Until then" has lost its antecedent (it should be "Until 2026-10-27").
3. **s27 notes** (C12 fix) say twice that the 2026-10-29 ATP pass reproduced the pre-registered pass row for row:
   - once in the rewritten s17c sentence ("On the 2026-10-29 sample, … reproduced the pre-registered pass row for row");
   - once more inside `rerun_sentence()` ("The earlier A-C13 run of 2026-10-29 reproduced … row for row; that statement is about the earlier sample …").
4. **PR wording** in the slide 9/10 notes, the condensed slide 15 notes and the slide 29 takeaway (§6). These are pre-existing, not introduced by the re-score.
5. **Code comment only:** the `make_deck.py` comment gives the pre-repair legacy CI upper bound as +0.0023. The code reproduces +0.0027, which is +0.003 at 3 dp either way.
6. **s42 notes** say "the panels do not differ significantly", scoped to the stratified CMH test (p 0.117, 0.99), which the generator used both before and after the re-score. The unstratified survivor test in the same file is significant for stratified against random: `emptied_strata_d100.json` survivors_pairwise original_vs_random Fisher p 0.023, chi² 0.027. A reader could take the phrase more broadly than the CMH scope. It is not stale.
