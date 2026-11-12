# Production E3 gate on the 10,000-organism panel (2026-11-12)

**Question.** Does the production form of the Stage E estimator (E3, `pfunc_estimator.py`) beat the
shipped k/n on real panel records, in the namespace production uses (PGFam), before it becomes the
builders' default for selections of two or more genomes?

**Criterion (declared before the run; first key of `data/pfunc_e3_gate.json` `_meta`).** PASS iff,
under self-exclusion at 0.987 in the n >= 2 stratum, E3 with species dedup has a lower pooled union
Brier than k/n AND a higher pooled union AUROC, each with the 95% cluster-bootstrap CI of the
difference excluding 0 (bootstrap over source genomes, B = 2000, seed 1729). Otherwise E3 ships opt-in.

**Verdict: PASS.**

| exclude_0.987, n >= 2 (19,158 amplicons, 6,776 sources, 68.7M cells) | k/n | E3 (species) | difference [95% CI] |
|---|---|---|---|
| union Brier | 0.1547 | 0.1017 | -0.0530 [-0.0547, -0.0513] |
| union AUROC | 0.8745 | 0.9283 | +0.0538 [+0.0517, +0.0560] |
| logistic calibration slope (clip 1e-3) | 0.44 | 1.10 | |
| mean q (prevalence 0.389) | 0.528 | 0.403 | |

At 0.995 the n >= 2 result is the same in direction and size (Brier 0.1446 -> 0.1018, AUROC
0.885 -> 0.930). At n = 1 k/n is identically 1 (AUROC 0.5); E3 reaches Brier 0.139 and AUROC 0.867
there, but Stage O did not establish a flux gain at n = 1, so production keeps k/n at n = 1.

**Script and outputs.** `scripts/pfunc_e3_panel_gate.py` (streams `data/selection.json` and the hits;
re-selects each amplicon with `score_genecap.exclude_self_record` + `select_representatives` exactly as
`genecap_tau0995.py`; cells = the union of the re-selected reps' PGFams; truth = the source's PGFams;
prior members exclude the source genome and its NCBI species). `data/pfunc_e3_gate.json` (summary),
`data/pfunc_e3_gate/per_amplicon.csv.gz`, `data/pfunc_e3_gate/run_meta.json`. Run 1,595 s, cgroup peak
21 GB.

## Caveats from the adversarial review (verified by recomputation; none changes the verdict)
1. **About 45% of the Brier gain is recalibration.** A saturated in-sample recalibration of k/n (empirical
   prevalence per (k, n), generous to k/n) reaches Brier 0.1311 and AUROC 0.8806; E3 still beats it by
   Brier -0.0294 and AUROC +0.0477. Do not present the full -0.053 as neighbour/prior information.
2. **The gain falls with n.** Brier difference by n: n = 2 -0.1007 [-0.1040, -0.0976]; n = 3 -0.0623;
   n >= 4 -0.0240 [-0.0254, -0.0226]; n >= 10 -0.0083 (AUROC +0.140 / +0.076 / +0.028). n = 2 is 9,096
   of the 19,158 amplicons, so the pooled figure depends on the panel's n mix; reweight to a deployment's
   n distribution before quoting an expected production effect.
3. **Genome-level re-entry.** Self-exclusion removes hits, not genomes: a genome with one 16S copy at
   >= 0.987 (removed) and another below (kept) re-enters (13.5% of n >= 2 amplicons; 1.19% of neighbour
   weight). It works against E3: on the 16,445 clean amplicons the Brier difference is -0.0597
   [-0.0615, -0.0577] and AUROC +0.0608.
4. **CIs cluster by source genome** (= species here). Clustering by source family widens the Brier CI
   about 1.9x ([-0.0564, -0.0499]); by order 2.7x. The verdict is unaffected.
5. **Calibration.** E3 is much better calibrated than k/n (slope 0.44) but mildly underconfident (slope
   1.10, CI excludes 1) and about 0.014 high on average.
6. **Missing neighbours.** Neighbours absent from the PGFam cache (32,021 distinct genomes) were dropped,
   not fetched; the production builder fetches them. On amplicons with no neighbour weight lost to
   missing or empty PGFam sets (11,557) the Brier difference is -0.0488 [-0.0507, -0.0468].
7. **What the gate does not test.** Flux (that is Stage O, on the EC route); the ~55% of production
   features with no PGFam, which keep k/n; product-string features outside PATRIC's vocabulary.

Integration and its checks: `pfunc_estimator.py`, `build_synthetic_genomes*.py` (`--pfunc`, default e3),
`reweight_build_e3.py`, `prfba_objective.py`; README.md section "Function probabilities (E3)".
