# Production E3 gate on the 10,000-organism panel (2026-11-12; re-scored 2026-11-15)

**Re-score (2026-11-15).** The gate was re-run on 2026-11-15 as part of the re-score of every result on
the repaired BV-BRC data. That repair (2026-11-13) replaced 1,006 shared-cache entries. The run used:
- the frozen cache (sha256 `05f03bcd...`) with a current catalogue;
- the re-scored gene-capture cells;
- code `3ffc97d`.

It ran on the second machine (TF) in its path-isolating namespace.

Controls passed before the run:
- Today's code, run on the pre-repair cache, the old catalogue and the 2026-11-12 reference cells,
  reproduced the committed `results`, `per_region_n2+`, verdict and diagnostics blocks exactly, and
  every per-amplicon row.
- The review recompute reproduced the committed `review_recompute.json`: counts exactly, floats to 1e-15.
- Two 500-amplicon smokes were identical.

The re-run reproduces the re-scored gene-capture cells in all 135,982 rows.

**The verdict is unchanged: PASS.** What moved:
- Repaired genomes now carry gene families, so 384 more n >= 2 amplicons are scored.
- k/n improves slightly, so the margin narrows: Brier -0.0530 -> -0.0520, AUROC +0.0538 -> +0.0501.

All numbers below are the re-scored ones. The 2026-11-12 values are in git (`4f7871d`, `5f7c7c1`).

**Question.** Does the production form of the Stage E estimator (E3, `pfunc_estimator.py`) beat the
shipped k/n on real panel records, in the namespace production uses (PGFam), before it becomes the
builders' default for selections of two or more genomes?

**Criterion (declared before the run; first key of `data/pfunc_e3_gate.json` `_meta`).** PASS iff,
under self-exclusion at 0.987 in the n >= 2 stratum, E3 with species dedup has a lower pooled union
Brier than k/n AND a higher pooled union AUROC, each with the 95% cluster-bootstrap CI of the
difference excluding 0 (bootstrap over source genomes, B = 2000, seed 1729). Otherwise E3 ships opt-in.

**Verdict: PASS.**

| exclude_0.987, n >= 2 (19,542 amplicons, 6,891 sources, 70.5M cells) | k/n | E3 (species) | difference [95% CI] |
|---|---|---|---|
| union Brier | 0.1536 | 0.1016 | -0.0520 [-0.0536, -0.0505] |
| union AUROC | 0.8781 | 0.9282 | +0.0501 [+0.0484, +0.0519] |
| logistic calibration slope (clip 1e-3) | 0.45 | 1.09 | |
| mean q (prevalence 0.388) | 0.532 | 0.402 | |

**Other modes and strata.**
- At 0.995 the n >= 2 result is the same in direction and size: Brier 0.1437 -> 0.1018, AUROC
  0.888 -> 0.930.
- At n = 1, k/n is identically 1 (AUROC 0.5). E3 reaches Brier 0.139 and AUROC 0.868 there.
  Stage O did not establish a flux gain at n = 1, so production keeps k/n at n = 1.

**Script and outputs.** `scripts/pfunc_e3_panel_gate.py` (streams `data/selection.json` and the hits;
re-selects each amplicon with `score_genecap.exclude_self_record` + `select_representatives` exactly as
`genecap_tau0995.py`; cells = the union of the re-selected reps' PGFams; truth = the source's PGFams;
prior members exclude the source genome and its NCBI species). `data/pfunc_e3_gate.json` (summary),
`data/pfunc_e3_gate/per_amplicon.csv.gz`, `data/pfunc_e3_gate/run_meta.json`.

**Run times.**
- Re-score run (TF, 10 workers): 944 s, cgroup peak 14.7 GB.
- The 2026-11-12 run (FM, 12 workers): 1,595 s, 21 GB.

## Caveats from the adversarial review (verified by recomputation; none changes the verdict)

**Sources of the numbers** (re-scored 2026-11-15):
- `data/pfunc_e3_gate/review_recompute.json`, the review's independent recompute, rerun with
  `scripts/pfunc_e3_gate_review_recompute.py`;
- the gate JSON;
- `data/pfunc_e3_gate/caveats_recompute.json` (`scripts/pfunc_e3_gate_caveats.py`), which holds:
  - the by-n CIs;
  - n >= 10;
  - the family and order clustering;
  - the no-lost-weight subset.

  That script uses the gate's own bootstrap draw. On the 2026-11-12 inputs it reproduces the gate's
  pooled CI and caveat 6 exactly. Its other CI bounds differ from the 2026-11-12 reviewer's ad hoc ones
  by at most 0.0003; see its docstring.

1. **About 48% of the Brier gain is recalibration.** A saturated in-sample recalibration of k/n uses the
   empirical prevalence per (k, n), which is generous to k/n.
   - It reaches Brier 0.1287 and AUROC 0.8840.
   - E3 still beats it by Brier -0.0271 and AUROC +0.0442.

   Do not present the full -0.052 as neighbour/prior information.
2. **The gain falls with n.**

   | n | Brier difference [95% CI] | AUROC difference |
   |---|---|---|
   | 2 | -0.0998 [-0.1032, -0.0965] | +0.136 |
   | 3 | -0.0617 | +0.072 |
   | >= 4 | -0.0229 [-0.0242, -0.0218] | +0.025 |
   | >= 10 | -0.0076 | |

   n = 2 is 9,265 of the 19,542 amplicons, so the pooled figure depends on the panel's n mix.
   Reweight to a deployment's n distribution before quoting an expected production effect.
3. **Genome-level re-entry.** Self-exclusion removes hits, not genomes. A genome with one 16S copy at
   >= 0.987 (removed) and another below it (kept) re-enters.
   - This affects 13.6% of n >= 2 amplicons and 1.17% of neighbour weight.
   - It works against E3: on the 16,753 clean amplicons the Brier difference is -0.0587
     [-0.0606, -0.0568] and the AUROC difference +0.0562.
4. **CIs cluster by source genome** (= species here).
   - Clustering by source family widens the Brier CI about 2.0x, to [-0.0552, -0.0488].
   - Clustering by order widens it 2.8x.

   The verdict is unaffected.
5. **Calibration.** E3 is much better calibrated than k/n (slope 0.45). It is mildly underconfident
   (slope 1.09, CI [1.07, 1.12] excludes 1) and about 0.014 high on average.
6. **Missing neighbours.** Neighbours absent from the PGFam cache (32,157 distinct genomes) were
   dropped, not fetched; the production builder fetches them. On the 13,623 amplicons that lost no
   neighbour weight to missing or empty PGFam sets, the Brier difference is -0.0496 [-0.0515, -0.0478].
7. **What the gate does not test.** Flux (that is Stage O, on the EC route), and the production
   features with no PGFam, which keep k/n: 10.8% (EmilyKin) and 7.9% (codiffusion) of the multi-genome
   genomes' features in the BV-BRC-only deployment builds (`deployment_d100/bvbrc_rebuild_summary.json`;
   about 55% in the earlier builds that mixed BV-BRC's RefSeq and PATRIC annotations).

Integration and its checks: `pfunc_estimator.py`, `build_synthetic_genomes*.py` (`--pfunc`, default e3),
`reweight_build_e3.py`, `prfba_objective.py`; README.md section "Function probabilities (E3)".
