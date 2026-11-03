# Better function probabilities, and an objective tuned to use them

Declared 2026-11-02, **post hoc**: after the weighted-objective test
(`WEIGHTED_OBJECTIVE_RESULTS.md`) and after the calibration study's exploratory top-k analysis
(`data/calibration/README.md`, item 10) had been seen. This is not a pre-registration. Every
tuned choice is made on data held apart from the data that reports it, and the configuration
this study recommends will be reported next to the shipped one, with its effect sizes on the
held-apart half.

## Why

- `p = k/n` ranks functions but is overconfident (calibration slope 0.78). Over all scored
  functions its Brier score (0.209) is worse than a constant 0.78 (0.173), because 86% of
  answered amplicons select one genome, where p = 1 for every function and presence is 0.79.
- On those single-genome selections the weighting is off by construction, yet a perfect
  signal would still move flux there (oracle at n = 1: −0.017 truth-absent flux).
- An unregistered top-k identity-weighted prevalence (top-20, τ = 0.01, genome dedup) lowered
  the universe Brier score from 0.089 to 0.063 and raised macro F1 at q ≥ 0.5 by +0.025,
  including on single-genome selections.
- In the objective, the cost `max(0.01, 1 − p)` is steep at the top and flat at the bottom.
  Gapfilled reactions get p = 0 (cost 1), but only 8.4% of them are supported by the true
  genome's own draft; gene-backed reactions at p < 1 are supported 49% of the time, at
  p = 1, 86%. The exchange term (weight 100) is 98.9% of the objective, so p only chooses
  among exchange-optimal solutions.

## Stage E: function-probability estimators (offline, no FBA)

**Population.** The controlled genus-blocked shared-reference benchmark, exactly as the
calibration study scored it: `cells.csv.gz`, arm `prfba_exact_tie`, space `ec`, outcome
`scored`, density 1.0. Truth: the source genome's EC set on the 2,235-EC universe
(`panel_ec_universe.json`). Holdout groups from `fold_plan.json`.

**Records.** Each fold's own hits, `data/prfba/f{fold}_d100/hits/`, searched against that
fold's reference only (training genera; species-dereplicated). The records hold 200 hits;
only the first 100 are used, because production records hold 100. Stage 0 exactly as
`select_references` (coverage ≥ 0.90, |SW − edlib identity| ≤ 0.02, identity ≥ 0.865), then
ranked by identity, then score.

**Estimators** (x_jf = 1 if neighbour j carries EC f; EC sets from the fold's
`_cache_ec_subset.json`):

- **E0, k/n** over the exact-tie set: the shipped rule.
- **E1, hits.** q_f = Σ_j w_j x_jf / Σ_j w_j over the top K stage-0 survivors, with
  w_j = exp(−(id_1 − id_j) / τ). Grid: K ∈ {20, 50, 100}; τ ∈ {0.0025, 0.005, 0.01, 0.02,
  0.05}; dedup by genome, or by species (the highest-identity genome of each species).
- **E2, clade prior.** π_c(f) is the share of the fold's reference genomes in clade c of the
  anchor (the first exact-tie genome) that carry f. c ∈ {genus, family}; a clade with fewer
  than 3 reference genomes backs off one rank. p = (k + α π_c) / (n + α), α ∈ {0.5, 1, 2, 4,
  8}. The anchor's genus is never the source's (genera are held out), so this is the
  nearest *available* clade's prior, and that is asserted.
- **E3, hits shrunk to the clade.** (Σ_j w_j x_jf + α π_c(f)) / (Σ_j w_j + α), K = 100,
  τ ∈ {0.005, 0.01, 0.02}, both dedups, c ∈ {genus, family}, α ∈ {0.5, 1, 2, 4}.
- **E4, recalibrated.** Out-of-fold isotonic recalibration of each family's selected variant.

**Tuning is nested.** For each outer fold, each family's variant is chosen on the other four
folds by universe log-loss (ε = 0.01) and scored on the outer fold. Reported performance is
the outer folds' predictions pooled. The hyperparameters carried into Stage O are the ones
chosen most often (ties: the smaller K, then the larger τ).

**Metrics.**
- Universe (all 2,235 ECs; an estimator's support may grow beyond the union): Brier,
  log-loss ε = 0.01, ECE over deciles.
- Union support (the union's ECs only: what the model's draft contains, the only part Stage O
  can use): Brier, WLS calibration slope, AUROC.
- Decision: macro and micro F1 of q ≥ 0.5 against the union, read against the SESOI of 0.02
  with score_folds' rule.
- Every metric overall and by selection size: n = 1, 2, 3, ≥ 4.

**Uncertainty.** Cluster bootstrap over holdout groups, score_folds' draw, 2,000 draws, seed
1729; source genome and genus as sensitivities. Paired contrasts against E0 use the same draws.

**Decision rule.** An estimator beats k/n if its nested-CV universe Brier and union-support
Brier differences from E0 both have CIs below 0.

**Assertions (all must be 0).**
- E0 reproduces `data/calibration/scores.json` (`exploratory_topk` `exact_tie_k_over_n`:
  universe Brier 0.08891, union macro F1 0.76721).
- E1 at K = 20, τ = 0.01, genome dedup reproduces its universe Brier, 0.06270.
- The source genome is never a neighbour or a clade member; no neighbour or clade member
  shares genus, species, taxon or holdout group with the source; every one is in the fold's
  reference; the anchor's genus differs from the source's.

**Caveats, stated now.**
- *Transfer.* The benchmark reference is about 7,000 species-dereplicated genomes per fold.
  Production searches about 315,000 genomes, where the top 100 hits of a well-sequenced
  species can be one species and identity weighting collapses onto the tie set. The species-
  dedup variant and, if cheap, the winner at densities 0.5 and 0.25 speak to this.
- *Support expansion* (adding relatives' functions to the draft at low p) changes gapfilling
  and is out of scope; universe metrics describe it, Stage O does not use it.

## Stage O: objective tuning (FBA)

**Sample.** The weighted-objective test's 360 amplicons (`data/weighted_objective/
population.json`), rebuilt identically: ModelSEEDpy c2e31a0, templates 21a8084, reversals
applied, complete medium, biomass ≥ 0.95 × max. Each rebuilt gapfilled model must match the
prior run's row (`gapfill.n_new_rxn`, `gapfill.growth`); `prfba_released_copy` and the
linear/ex-100 real arm's `released_check` are kept unchanged.

**Split.** A holdout group is in the *tune* half if sha256("1729|" + group_id) is even, else
in the *confirm* half. All selection happens on the tune half.

**Factors.**
- P, probability source: `kn` (the shipped k/n); `est` (the Stage E winner, at role level:
  q(role) = Σ_j w_j·1[every EC of the role ∈ EC(g_j)] / Σ_j w_j, never a product of per-EC
  values; E3's clade term is added the same way); `est_shuffled` (est permuted among the
  gene-backed reactions, seeded) for specificity.
- C, cost shape on rxn reactions: `lin` = max(0.01, 1 − p) (shipped); `log` =
  max(0.01, ln max(p, 0.01) / ln 0.01). Both give cost 1 at p = 0 and 0.01 at p = 1, so the
  gene-less draft term (80% of the internal term) is unchanged and only the shape differs.
- G, gapfilled reactions' p: 0 (shipped); the tune half's measured truth-support rate of
  gapfilled reactions; 0.2.
- X, exchange weight: 100 (shipped), 10, 1.

Every (C, G, X) cell also runs `two_tier` (gene-backed p = 1) and `oracle` (gene-backed
p = 1 if truth-supported else 0), with the same gapfilled p.

**Measures.** As in the weighted-objective test: truth-absent gene-backed flux on the optimal
face (objective ≤ opt + 1e-6), min / mid / max, divided by the per-amplicon constant D (the
two-tier, linear, ex-100 arm's gene-backed flux; fixed across all configurations). Guard
metrics: total exchange flux, gene-backed flux and flux precision against the shipped arm's;
the biomass floor met.

**Selection (tune half).** Eligible: configurations whose median total exchange flux and
gene-backed flux are within 10% of the shipped arm's. Among them, the lowest pooled TAF face
midpoint, pooled two ways and both reported:
- **A**, n ≥ 2 strata weighted by eligible source genomes, as before (primary, comparable);
- **B**, all answered amplicons, n = 1 at its population weight (86.8%). The n = 1 stratum
  has 60 rows (about 30 per half), so B is thin; it is reported, not used to select.
Ties within 0.002: the configuration with fewer departures from the shipped one. The G factor
is read only on rows where gapfilled flux varies across G's levels.

**Confirmation (confirm half).** The selected configuration against the shipped one and
against its own two-tier and shuffled arms: TAF face min / mid / max with CIs against the
SESOI of 0.01; share of the oracle's headroom captured; degeneracy by FVA on the optimal face
(finalists only, a second pass).

**Uncertainty.** Cluster bootstrap over holdout groups, 2,000 draws, as before.

## Order of work

1. This document, committed before any Stage E or Stage O number is computed.
2. Stage E (`scripts/pfunc_estimators.py` → `data/pfunc_estimators/`), and the Stage O harness
   (`scripts/objective_tuning.py`) written and smoke-tested in parallel.
3. Stage O launched with Stage E's selected estimator, as a `systemd-run --user` unit.
4. Results in `PFUNC_ESTIMATORS_RESULTS.md` and `OBJECTIVE_TUNING_RESULTS.md`, each checked by
   an independent re-derivation before any manuscript or deck text changes.
