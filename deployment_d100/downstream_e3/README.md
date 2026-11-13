# What E3 changes downstream of the depth-100 deployments

Run on 2026-11-12. The comparison is the shipped k/n builds (`*/synthetic_genomes_d100/`) against
the E3 reweighted copies (`*/synthetic_genomes_d100_e3/`, `../e3_reweight_summary.json`).

* **Descriptive only.** There is no ground truth here, only the k/n-vs-E3 difference.
* Nothing in `codiffusion_bioreactor`, `EmilyKin`, `ModelSEEDpy*`, `ModelSEEDTemplates*` or the
  shipped builds was modified. Nothing was uploaded and nothing was committed.
* The library code came from `git archive` exports of fixed commits (`common.py`):
  * ModelSEEDpy 09856f5 (`fix-prfba-single-model`);
  * MSCommunity 03bd80c;
  * ModelSEEDTemplates b6fe86b (codif-curated, part 2) and 21a8084 (c12run, part 3);
  * ModelSEEDDatabase 5389554a.

  Every results JSON records its commits, template sha256s and script hashes.

## Headline

1. **No consumer reads either d100 build today.**
   * codiffusion builds from `GDIR = "synthetic_genomes_54d1150"`
     (`build_prfba_community_models.py:41`), the 20-hit build.
   * EmilyKin's `ek_pipeline.py` reads `config.json` / `config_r13_syn.json` ->
     `EmilyKin/synthetic_genomes/genome_objects` (20-hit). `config_r13.json` points at
     `runs/r13/genome_objects_overlay`, a planned MAG overlay that does not exist yet.

   So E3 on d100 changes nothing downstream now. Parts 2 and 3 measure what would change if the
   pipelines were pointed at the E3 build.
2. **No production decision depends on p.**
   * Every member gapfill (codiffusion mediagf, legacy gapfill_models, EmilyKin mediagf) calls
     MSGapfill without `reaction_scores`, so every candidate's factor is 2 − 0 = 2.
   * Neither project ever calls `MSProbability.prFBA` (codiffusion `util.py:19` imports it and
     never uses it). The community simulations do not use p.
   * Every other production read records p in a report (EmilyKin atpcorr / gate / assimilation
     records, `ppk_max_probability`).
   * The ppk decision is presence-based: the bucket key is created for any p.
3. **p > 0 is not strictly invariant.**
   * E3 is exactly 0 for 67 features: 47 codiffusion, 20 EmilyKin
     (`results/census.json` `probability_scan.*.zero_features`). k/n is never 0.
   * The 67 sit in 3 codiffusion and 3 EmilyKin genomes. None of those genomes is in the part-2
     community or the part-3 sample, so parts 2 and 3 never meet an E3 zero (0 gene-backed
     reactions with E3 p = 0 in the 80 models).
   * **Fixed after this study (2026-11-12):** a merged feature is carried by at least one built
     source genome, so an E3 of exactly 0 is a disagreement between the genome-dict cache and the
     PGFam cache. `pfunc_estimator.reweight_asv` now keeps k/n for such a feature
     (`n_zero_kept_kn` in the diagnostics), and both `synthetic_genomes_d100_e3` builds were
     regenerated with it; p > 0 holds again everywhere. Parts 2 and 3 are unaffected (none of the
     six genomes is in them).
4. **Compartment-aware gapfill penalties** (prototype, part 2):
   * Across the 88 members of day 174's community, priced one member at a time, E3 changes the
     penalty of 8,589 candidates. That is 1.17% of the 733,964 single-member candidates (summed over
     members; e0 candidates count once per member), and 91.7% of the 9,364 member-scored
     candidates of the 41 reweighted members.
   * The median penalty change is +3.2%.
   * The cheapest member changes for 405 of the 8,524 reactions offered by at least two members.
     Only 12 of these are a change between two different unique cheapest members. 336 are E3
     breaking an exact tie.
   * Community gapfilling has no solution with stock, k/n or E3 penalties (6/6 attempts). This
     matches the demo.
5. **prFBA flux** (part 3; 80 n ≥ 2 genomes and 10 n = 1 controls; all solved):
   * **Growth:** unchanged in 80/80, because the 95% floor binds.
   * **Tuned objective:** the median share of internal flux that moves is 1.7% per genome
     (pooled 3.1%; maximum 26.7%). Exchange fluxes move in 58/80 genomes.
   * **The two optima are nearly interchangeable.** Priced at k/n costs, the E3 optimum costs a
     median 0.07% more than the k/n optimum (max 5.6%). Priced at E3 costs, the k/n optimum costs
     a median 0.14% more than the E3 optimum (max 1.9%).
   * **Released objective (ex_weight 100):** the median share of flux that moves is 0.8%.
     Exchange fluxes move in 15/80 genomes.
   * **Controls:** all 10 n = 1 controls are identical in every number.

## 1. Consumer census (`census.py` -> `results/census.json`)

* **The sweep.** It greps both repos (code `*.py` / `*.sh` and notebook code cells, outside the
  data directories) for `probability|reaction_scores|MSProbability|prFBA(`.
  * 170 hits in 35 files.
  * Each file is either a row of the table below or classified as a non-consumer, with the reason
    (`sweep_summary.files_not_consumers`): k/n producers, archived copies, design text,
    SPIEC-EASI.
  * 0 files are unclassified.
* **The table.** It is re-verified on every run: each file:line anchor must still match its
  regex (0 failures).

| consumer | file:line | path | build read | can E3 change its output? |
|---|---|---|---|---|
| draft notes writer | codif `build_prfba_community_models.py:356,361-362` | production (stage build) | 54d1150 | notes values only |
| member gapfill | codif `build_prfba_community_models.py:541,558` | production (mediagf) | 54d1150 drafts | no (no `reaction_scores`) |
| legacy gapfill | codif `gapfill_models.py:600` (`reaction_scores={}`) | legacy | bundled-template drafts | no |
| 2 − p penalties | codif `compartment_aware_gapfill.py:76,113,120` | PROTOTYPE | `models_<tag>` drafts (54d1150) | yes: part 2. A missing note defaults to p = **1.0**, not 0 (l.76) |
| demo | codif `demo_compartment_scoring.py:109,140` | demo | `models_v7` (54d1150) | penalties only |
| loader `build_model` | codif `load_genomes_to_modelseedpy.py:280,355` | old path (`modeling.py` marimo, May 2026) | pre-prFBA `genome_objects/` | notes only |
| W-L audit | codif `methane_nongrowth/atp_wl/wl_audit/audit_merge.py:194,262` | analysis | 54d1150 (`audit.py:45`) | only if re-pointed: **p ≥ 0.5** decides GAP vs `--(weak gene)` |
| dead import | codif `util.py:19` | none | - | no |
| KBase base module | codif `kbbasemodules/basemodelingmodule.py:263` | vendored, not imported | - | no |
| ppk evidence | EK `ek_pipeline.py:769-780` | production (templates) | 20-hit | recorded value only; the decision is presence |
| draft notes writer | EK `ek_pipeline.py:879,894-895` | production (build) | 20-hit / planned MAG overlay | notes only |
| atpcorr record | EK `ek_pipeline.py:2515,2555,2603` | record | 20-hit drafts | values; the `p > 0` list can change only through the 67 E3 zeros |
| assimilation / gate records | EK `ek_pipeline.py:2397,4496,5857` (order by `gene_bearing`, l.5049) | record | 20-hit | values only |
| member gapfill | EK `ek_pipeline.py:2792,3211,3232` | production | 20-hit | no (no `reaction_scores`; `curate_templates.py:433` says so) |
| other gapfill calls | codif `gapfill_syntrophic.py:210`, `regapfill_media_v2.py:112`, `modeling.py:209`; EK `ek_pipeline.py:2536` (ATP correction's internal gapfill) | production / analysis | 20-hit | no (none passes `reaction_scores`; added after the verification) |
| report | EK `report_fluxes.py:864` | report | run output | printed values |
| storage census | EK `plan_r13/maps/storage_work/gene_census.py:91,157` | analysis | 20-hit | only if re-pointed: the **p ≥ 0.5** counts |
| other planning scans | EK `vocab_scan.py:28`, `genome_scan.py:35`, `nif_probe.py:38`, `d14_verify.py:59` | analysis | 20-hit | tabulated values only |

**Probability scan** (`census.py --scan`). It covers every reweighted genome: codiffusion 1,034
genomes with 5,808,409 features, and EmilyKin 1,410 genomes with 7,646,149 features.

| | codiffusion | EmilyKin |
|---|---|---|
| replaced features | 2,758,484 | 3,423,718 |
| replaced features whose value changed | 2,708,433 | 3,358,646 |
| E3 range | [0, 1] | [0, 1] |
| features with E3 = 0 | 47 | 20 |
| minimum k/n of the replaced features | 0.0127 | 0.0130 |
| crossings of 0.5, k/n ≥ 0.5 -> E3 < 0.5 | 496,140 | 613,336 |
| crossings of 0.5, the other way | 70,377 | 73,696 |
| genomes with any 0.5 crossing | 1,033 of 1,034 | 1,405 of 1,410 |
| CDS / feature mismatches | 0 | 0 |

So the two `p ≥ 0.5` analyses would change if re-pointed.

## 2. Community gapfill penalties (`penalties.py` -> `results/penalties.json`)

**How the members were built**
* The demo's measurement (`demo_compartment_scoring.penalties_for`, imported) was reproduced on
  day 174 (Wolfe medium + H2/CO2), with the gram_neg community template as in the demo.
* Member drafts come from `build_prfba_community_models._build_one`, imported (not edited). Its
  module globals were pointed at absolute paths, and it wrote into
  `models/part2/{kn,e3}/models_v7e/`. It was run once on each build's `genome_objects`.
* Template classes come from that module's `stage_templates()` run on the d100 lineages. They are
  identical to the production `template_assignment_v7.json` for all 88 members.
* Templates: codif-curated b6fe86b.

**Member to genome object.** A member is a 20-hit ASVset. Its d100 genome object is the file named
by the same representative ASV md5.
* **All days:** 1,718 ASVsets, 0 unmapped, 646 reps reweighted. 84 reps have a d100 selection that
  differs from their 20-hit one, so a d100 draft is not the production draft. The comparison is
  d100-k/n against d100-E3.
* **Day 174:** 88 members, 0 unmapped, 41 reweighted (29.7% of the day's abundance), 6 reps with
  a changed selection. No ASVset would split under d100.

**The drafts.** 88 × 2 were built with 0 errors. They have the same reaction ids and 0 structural
differences: stoichiometry, bounds and GPR compared as a DNF set.
* 3,998 GPR strings differ only in term order. MSBuilder iterates object sets, whose order is not
  reproducible between processes.
* 39,433 reaction p values differ, all of them in the 41 reweighted members.

**Penalties.** Δfactor = factor(E3) − factor(k/n) = p_kn − p_e3; relative = penalty(E3) /
penalty(k/n) − 1.

| set | candidates | changed | factor up / down | Δfactor median [p05, p95] (range) | relative change median [p05, p95] (range) | cheapest-member change (unique / tie broken / tie created / other) |
|---|---|---|---|---|---|---|
| A: demo's top 4 (Methanobacterium.2, Methanobacteriaceae.1, Burkholderiales.1, Methanobacterium.1) | 32,075 | 416, all in Burkholderiales.1 (n = 29; 416 of its 422 member-scored) | 274 / 142 | +0.0022 [−0.018, +0.106] (−0.966 … +0.342) | +0.18% [−1.0%, +8.4%] (−49.1% … +31.1%) | 0 / 47 / 2 / 7 of 8,247 |
| B: top 4 reweighted (Burkholderiales.1 n = 29, Petrimonas.1 n = 4, Methanobacterium.3 n = 2, Acinetobacter.1 n = 37) | 30,549 | 1,050 (416 / 158 / 99 / 377) | 782 / 268 | +0.0052 [−0.011, +0.210] (−0.973 … +0.749) | +0.39% [−0.8%, +20.0%] (−49.3% … +74.9%) | 59 / 58 / 1 / 19 of 7,845 |
| C: all 88 day-174 members (single-member pricing) | 733,964 | 8,589 (41 members; median 187 per member, 75-427; 0 in the 47 n = 1 members) | 6,606 / 1,983 | +0.042 [−0.146, +0.433] (−0.973 … +0.949) | +3.2% [−10.1%, +41.7%] (−49.3% … +94.9%) | 12 / 336 / 1 / 56 of 8,524 |

* **Which candidates can change.** A member-scored candidate is the reverse direction of a
  gene-bearing draft reaction of that member. New template reactions carry no member evidence, so
  their factor stays 2.
* **Set C rests on an equivalence.** For each of the 7 members of sets A and B, its single-member
  penalty vector equals its community-compartment vector in every candidate (stock, k/n, E3).
* **Members priced differently** (demo metric), set A: stock 240, k/n 412, E3 446 reactions. The
  demo itself (20-hit v7 drafts) had 224 and 368, with 32,291 candidates.
* **Largest changes**
  * `rxn60858_c3` in Burkholderiales.1: p 0.138 -> 1.0, penalty 7,448 -> 4,000.
  * Methanobacterium.3 `rxn09333`: p 1.0 -> 0.350, penalty 2,008 -> 3,312.
* **Unique cheapest-member changes** in the 88-member community
  * `rxn01094`: Sphingopyxis.1 -> Ochrobactrum.1.
  * `rxn06078`: Burkholderiales.1 -> Acinetobacter.1.
* **Prototype caveat.** 73 candidates (rxn02028, rxn09296) have negative base penalties in
  b6fe86b, and 31 of them are among the changed ones. There, 2 − p inverts: more evidence gives a
  less negative penalty, so the candidate becomes dearer.
* **Does community gapfilling solve?** No.
  * Each of sets A and B was tried with stock, k/n and E3 penalties: 6 attempts in all.
  * Every attempt stops at "No gapfilling solution found before filtering for Wolfe activating
    bio1", after 43-53 s. The test LP is optimal with objective 0.
  * That test is independent of the penalties, so no MILP is ever priced. E3 cannot change this
    outcome.

## 3. prFBA flux (`flux.py` -> `results/flux.json`, `results/part3/{deployment}/{asv}.json`)

**Sample.** 40 n ≥ 2 genomes per deployment, 10 in each of the strata n = 2, 3-4, 5-9 and ≥ 10,
ranked by `sha256("downstream_e3|" + asv)`. Plus 5 n = 1 controls per deployment.

**Per genome**
1. Two drafts with `_build_one` (the production rule), on c12run 21a8084. They are compared
   structurally: 0 differences in 90/90, plus 7,608 order-only GPR strings.
2. One gapfill of the k/n draft, Stage O's way (inlined `phenotype_fold.gapfill_block`):
   * MSATPCorrection;
   * MSGapfill on the KBase complete medium at uptake 100;
   * reversals applied.

   All 90 grow. Median 44 new reactions for n ≥ 2; 18 reversals applied in all.
3. Two copies annotated with `prfba_objective.set_reaction_probabilities`; gapfilled reactions get
   p = 0. 88 draft reactions in 3 genomes lost their GPR: ATP correction removed them and the
   gapfiller re-added them gene-less. All 88 are in those genomes' `gapfill.new_rxn`, so their p is
   0 in both arms.
   * In the same 3 genomes many more draft reactions are absent from the gapfilled models: 593
     (EmilyKin `96e917ee…`), 490 (codiffusion `6dbf54df…`) and 89 (codiffusion `5aed6ebc…`);
     the `draft_rxn_net_removed` figures (631, 529, 95) also count the 39/43/6 reactions that came
     back gene-less.
   * In the other 87 genomes `draft_rxn_net_removed` is at most 2, and negative for most (median
     −4, min −48), because ATP-correction gapfilling adds core reactions that `gapfill.new_rxn`
     does not list; they are gene-less, so p = 0 in both arms and no result changes.
4. `run_prfba`, tuned and released, environment = the complete medium at 100.

**Runs.** 1,850.6 s wall on 8 workers. 55-241 s per genome; peak worker RSS 3.99 GB.
* The k/n tuned re-solve on a fresh copy was identical in 90/90, with a maximum flux difference of
  0.
* In all 360 solves, the reported objective equals the internal plus the exchange term recomputed
  from the fluxes, to within 1.6e-15 (relative). So the cross-pricing below uses the LP's own cost.
* The summary block of `results/flux.json` was regenerated with `flux.py --summary-only` after the
  run, which added the per-stratum and symmetric cross-pricing summaries. `_meta.code` therefore
  hashes the final `flux.py`. The per-genome computation is unchanged since the run.

**Pooled n ≥ 2 (80 genomes).** In each genome p changes on a median of 1,023 of 1,103 gene-backed
reactions. The per-genome mean Δp over the changed reactions has a median of −0.050 (mean −0.069).

| | tuned (ex_weight 1) | released (ex_weight 100) |
|---|---|---|
| growth changed | 0 / 80 | 0 / 80 |
| flux vector identical | 3 / 80 | 3 / 80 |
| re-routed share of internal flux, per genome: median [p10, p90] (max) | 0.0173 [0.0018, 0.0616] (0.2668) | 0.0076 [0.0013, 0.0268] (0.0791) |
| re-routed share, pooled | 0.0307 | 0.0116 |
| Jaccard of active reaction sets, median (min) | 0.952 (0.867) | 0.969 (0.931) |
| genomes with any exchange-flux change | 58 | 15 |
| exchange L1 change, median (max) | 1.49 (75.4) | 0 (18.9) |
| total exchange flux, sum k/n -> E3 | 21,542.3 -> 21,495.6 | 20,685.0 -> 20,685.0 |
| total uptake, sum | 10,447.5 -> 10,459.8 | 10,258.1 -> 10,257.8 |
| gene-backed flux, sum | 61,971.5 -> 61,488.7 (−0.78%) | 63,260.4 -> 63,303.6 (+0.07%) |
| share of flux through p-changed reactions, median | 0.761 -> 0.761 | 0.780 -> 0.781 |
| objective Δ, median [p10, p90]; E3 higher / lower | +23.3 [−28.8, +95.3]; 57 / 23 | +24.1 [−35.7, +100.9]; 55 / 25 |
| of which repricing alone (k/n flux at E3 prices − k/n optimum), median | +23.8 | +24.2 |
| E3 optimum's excess cost at k/n prices, relative: median (max) | 0.0007 (0.056) | 0.0000 (0.0006) |
| k/n optimum's excess cost at E3 prices, relative: median (max) | 0.0014 (0.019) | 0.0000 (0.0003) |
| flux-weighted mean p of gene-backed flux (own arm), median | 0.864 -> 0.817 | 0.837 -> 0.805 |

**Per deployment (tuned)**

| | re-routed share, median | pooled | exchange changed | gene-backed flux, sum |
|---|---|---|---|---|
| codiffusion | 0.021 | 0.034 | 29 / 40 | 26,876 -> 26,570 |
| EmilyKin | 0.017 | 0.028 | 29 / 40 | 35,095 -> 34,919 |

Under the released objective the codiffusion and EmilyKin pooled shares are 0.012 and 0.011.

**By n stratum (tuned)**

| n | re-routed share, median | pooled | mean Δp, median |
|---|---|---|---|
| 2 | 0.023 | 0.037 | −0.142 |
| 3-4 | 0.017 | 0.024 | −0.047 |
| 5-9 | 0.015 | 0.014 | −0.051 |
| ≥ 10 | 0.013 | 0.052 | −0.008 |

The ≥ 10 pooled value is driven by one genome. In codiffusion `400fe0a1…` (n = 10, gram_pos),
26.7% of flux is re-routed: gene-backed flux 1,031.8 -> 733.0, exchange 289.1 -> 330.9, growth
unchanged.

**n = 1 controls (10).** Identical drafts and identical p. Every metric under both configurations
is identical: rerouted 0, exchange 0, objective Δ 0.

**How to read the objective change.** The objective change is almost all repricing: the same flux
costs more because E3 lowers p on most changed reactions. Flux re-routing is small. LP
degeneracy means a re-route with an unchanged objective can be an alternative optimum, so the
objective numbers and the cross-pricing numbers are the clean ones. The rest is descriptive.

## Not done / caveats

* **Part 3 is not the deployments' own pipelines.**
  * It uses codiffusion's draft rule on c12run templates, with a complete-medium gapfill, for both
    deployments.
  * EmilyKin's own builder (round-12 role layer, curated r13 templates) and the day media were not
    used.
  * EmilyKin template classes use the genome lineage only: there is no MiDAS table keyed by ASV
    md5.
* **Part 2 prices candidates, as the demo does.** No community gapfill solution exists to compare.
* **Size.** `results/part2_work/single/` (52 MB of per-member penalty vectors) and `models/`
  (about 1 GB of SBML) are git-ignored. Both can be regenerated.
* **Pins.** The pinned exports live in the session scratchpad. `common.ensure_pins()` recreates
  them from the recorded commits.

## Reproduce

Each script is run under its own unit, the way the runs on 2026-11-12 were:

```
systemd-run --user --unit=ds-e3-<name>-$(date +%s) --collect -p MemoryMax=40G -p MemorySwapMax=0 \
  -p WorkingDirectory=/home/freiburger/Documents/prFBA -E PYTHONDONTWRITEBYTECODE=1 -E PYTHONHASHSEED=0 \
  /home/freiburger/Documents/py_venv/bin/python -W ignore \
  /home/freiburger/Documents/prFBA/deployment_d100/downstream_e3/<script> <args>
```

| script and arguments | what it does | time on 2026-11-12 |
|---|---|---|
| `census.py --scan` | the sweep, the table and the probability scan | about 45 s for the scan |
| `penalties.py --workers 8` | part 2 | 983 s |
| `penalties.py --post-only` | recomputes the cheapest-member block from the saved single-member vectors | seconds |
| `flux.py --workers 8` | part 3 | 1,851 s |
| `flux.py --summary-only` | re-summarises the saved per-genome files | seconds |

* The three model-building scripts are resumable: finished drafts and finished genomes are
  skipped.
* The worker pools are `multiprocessing.Pool(maxtasksperchild)` with spawn.
  `ProcessPoolExecutor(max_tasks_per_child)` deadlocked under Python 3.12 once workers retired.

## Files

| file | role |
|---|---|
| `common.py` | paths, pinned exports, provenance, canonical GPR |
| `census.py` | part 1 |
| `penalties.py` | part 2 |
| `flux.py` | part 3 |
| `results/census.json`, `census_sweep.json` | consumer table, sweep hits, probability scan |
| `results/penalties.json` | mapping, draft comparison, sets A/B/C, equivalence, cheapest member, gapfill attempts |
| `results/part2_work/` | the day-174 sets file, template assignment, draft log and comparison (per-member vectors in `single/`, ignored) |
| `results/flux.json` | sample rule, pooled / per-deployment / per-stratum summaries, per-genome table |
| `results/part3/{deployment}/{asv}.json` | per genome: drafts, gapfill, p changes, prFBA metrics, pair metrics, nonzero fluxes |
| `logs/` | run logs (git-ignored, as is every `*.log`) |

## Scope notes added after the verification
- The sweep covers the two local repositories. It skipped vendored libraries and `runs/` snapshots;
  those were checked separately: ModelSEEDpy's gapfiller reads probabilities only when
  `reaction_scores` is passed, MSCommunity has no probability reads, and `MSCommunity.gapfill`
  (which accepts `reaction_scores`) is never called.
- KBase-side readers of the uploaded k/n builds were not checked (for example
  KB-ModelSEEDReconstruction's `modelseedrecon.py:510` reads `rxn.probability`, which comes from
  annotation-ontology scores); the E3 builds are not uploaded.
