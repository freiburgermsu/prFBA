# What E3 changes downstream of the depth-100 deployments

The study compares k/n with E3 feature probabilities downstream of the two depth-100 deployments.
It has been run on two build sets. Each is defined in `common.BUILD_SETS`, and every script takes
`--builds <set>`.

| build set | k/n build | E3 build | results, models | run |
|---|---|---|---|---|
| `bvbrc` (current) | BV-BRC-only k/n rebuilds (`build_d100.py --pfunc kn`, in the session scratchpad: `bvbrc_only/kn_codif`, `bvbrc_only/kn_emilykin`) | `*/synthetic_genomes_d100_bvbrc/` | `results_bvbrc/`, `models_bvbrc/` | 2026-11-13 |
| `mixed` | the shipped `*/synthetic_genomes_d100/` | `*/synthetic_genomes_d100_e3/`, their E3 re-weightings (deleted 2026-11-13) | `results/`, `models/` | 2026-11-12 |

* **Two vocabularies.**
  * The mixed builds merged BV-BRC's PATRIC CDSs with the RefSeq CDSs it imports. That was the
    rule up to prFBA 3ce12d0.
  * The BV-BRC-only builds keep the PATRIC CDSs only (`../README.md` section 10).
  * The k/n rebuild commands are in `../bvbrc_rebuild_summary.json`
    (`deployments.*.runs.kn.command`).
* **Same genomes.** The BV-BRC-only builds reuse the shipped selection audits, so every ASV keeps
  its selected genomes, tier and n. The scripts check that the following are identical in both
  runs:
  * the part-2 member mapping, template classes and sets A and B;
  * the part-3 sample.
* **Descriptive only.** There is no ground truth here, only differences.
* **Nothing else touched.**
  * Nothing in `codiffusion_bioreactor`, `EmilyKin`, `ModelSEEDpy*`, `ModelSEEDTemplates*` or the
    builds was modified.
  * Nothing was uploaded and nothing was committed.
* **Pinned library code.** Both runs import from the same `git archive` exports of fixed commits
  (`common.py`):
  * ModelSEEDpy 09856f5 (`fix-prfba-single-model`);
  * MSCommunity 03bd80c;
  * ModelSEEDTemplates b6fe86b (codif-curated, part 2) and 21a8084 (c12run, parts 3 and 4);
  * ModelSEEDDatabase 5389554a.

  Every results JSON records its commits, template sha256s, script hashes and build set.
* **Same codiffusion code.** The four codiffusion files the scripts import have the same sha256 in
  both runs: `build_prfba_community_models.py`, `compartment_aware_gapfill.py`,
  `demo_compartment_scoring.py` and `load_genomes_to_modelseedpy.py`. This holds although
  codiffusion's HEAD moved from d4b810a to daa0505.

## Headline (BV-BRC-only builds, 2026-11-13)

1. **The k/n-vs-E3 conclusions of the mixed run hold.** The next section compares them line by
   line.
2. **Still no consumer and no decision.**
   * No consumer reads a d100 build, and no production decision depends on p.
   * The census was re-run: 19/19 table entries verified, 0 anchor failures, 0 unclassified sweep
     files.
3. **p > 0 everywhere.** No feature of either BV-BRC-only E3 build has p = 0.
4. **E3 now covers most features.**
   * It replaces 92.1% (codiffusion) and 89.1% (EmilyKin) of the features of the reweighted
     genomes. In the mixed builds it replaced 47.5% and 44.8%.
   * It replaces about as many features as before (2.76 M and 3.42 M). The RefSeq features it
     could not resolve, which kept k/n, are gone.
5. **Penalties** (part 2):
   * E3 changes 8,628 of 734,632 single-member candidates (1.17%), by a median +3.2%.
   * The unique cheapest member changes for 12 reactions.
   * Community gapfilling has no solution (0 of 6 attempts).
6. **prFBA flux** (part 3):
   * Growth is unchanged in 80/80.
   * The median re-routed share is 1.74%.
   * The two optima stay interchangeable: the median excess cost is 0.07% / 0.14%.
   * The upper tail is heavier than in the mixed run (pooled 3.95% against 3.07%). Part of it is
     alternative optima.
7. **Dropping the RefSeq features** (part 4, new; k/n drafts on the 90 part-3 genomes):
   * Features per genome halve (median ratio 0.49).
   * The drafts lose 0.76% of their reactions: 859 in 79 of 90 genomes, and none gained. They
     lose 0.72% of their gene-backed reactions.
   * 4.6% of the gene-backed reactions in both drafts lose genes from their GPR.
   * After the same complete-medium gapfill all 90 grow. Growth is unchanged at the median, and
     changes by more than 2% in 6 genomes.
   * The largest change is a fall from 2.44 to 0.90 (codiffusion `dd34ba85`, n = 4). That genome
     loses 30 draft reactions that only RefSeq product names supported. Which of them costs the
     growth was not traced.
   * The gapfill burden rises by 0.4% (4,888 -> 4,908 reactions).

## Mixed vs BV-BRC-only: do the conclusions change?

Same sample, members, templates, pins and code. Only the builds differ.

| question | mixed (`results/`) | BV-BRC-only (`results_bvbrc/`) | conclusion |
|---|---|---|---|
| does any consumer read a d100 build? | no | no | same |
| does any production decision read p? | no | no | same |
| features with p = 0 (E3 builds) | 67 (fixed afterwards, then 0) | 0 | same since the fix |
| share of features E3 replaces, codiffusion / EmilyKin | 47.5% / 44.8% | 92.1% / 89.1% | E3 now reaches most features |
| crossings of p = 0.5, k/n ≥ 0.5 -> E3 < 0.5, codiffusion / EmilyKin | 496,140 / 613,336 | 489,506 / 605,526 | same |
| part 2: candidates whose penalty E3 changes (set C) | 8,589 of 733,964 (1.17%) | 8,628 of 734,632 (1.17%) | same |
| part 2: median relative penalty change | +3.2% | +3.2% | same |
| part 2: unique cheapest-member changes | 12 of 8,524 | 12 of 8,525 (the same 12 reactions and members) | same |
| part 2: community gapfill solves | 0 / 6 | 0 / 6 | same |
| part 3: growth changed by E3 | 0 / 80 | 0 / 80 | same |
| part 3: re-routed share, tuned, median per genome | 1.73% | 1.74% | same |
| part 3: re-routed share, tuned, pooled; p90 | 3.07%; 6.2% | 3.95%; 11.2% | heavier tail (see part 3) |
| part 3: E3 optimum's excess cost at k/n prices, median (max) | 0.07% (5.6%) | 0.07% (5.6%) | same |
| part 3: k/n optimum's excess cost at E3 prices, median (max) | 0.14% (1.9%) | 0.14% (3.6%) | same median |
| part 3: re-routed share, released, median; pooled | 0.76%; 1.16% | 0.80%; 2.47% | same median, heavier tail |
| part 3: genomes with exchange changes, tuned / released | 58 / 15 | 56 / 16 | same |
| part 3: n = 1 controls identical | 10 / 10 | 10 / 10 | same |

So every conclusion about E3 carries over. E3 still changes p on almost every gene-backed reaction
of an n ≥ 2 model: a median 1,026 of 1,095, against 1,023 of 1,103 before. Its effect on penalties
and fluxes stays small.

The only shift is in the flux tail.
* Six genomes re-route at least 5 points more than in the mixed run.
* At least one of them is an alternative optimum. EmilyKin `9fc55d66` has the same objective
  (548.4) and zero excess cost at either price set, yet 11.2% of its flux moves (0.2% in the mixed
  run).
* Others reflect a different model. ATP correction pruned 490 reactions from codiffusion
  `6dbf54df`'s mixed draft but only 2 from its BV-BRC-only draft.

The vocabulary itself matters more for some models than E3 does (part 4).
* In part 3, E3 never changes growth, because the 95% floor binds.
* Dropping the RefSeq features changes the maximal growth of the gapfilled model by more than 1%
  in 9 of 90 genomes, and in one by −63%.
* At the median, though, a BV-BRC-only draft is the mixed draft minus 7 reactions.

## BV-BRC-only builds (2026-11-13, `results_bvbrc/`)

### 1. Consumer census (`census.py --builds bvbrc --scan` -> `results_bvbrc/census.json`)

* **The table** (shown in the mixed section below). All 19 entries verified, 0 anchor failures.
  codiffusion still builds from `synthetic_genomes_54d1150`, and EmilyKin's configs still point at
  the 20-hit build. The table is unchanged apart from the wording of the two rows that cited the 67
  E3 zeros.
* **The sweep.** 180 hits in 43 files, against 170 in 35 on 2026-11-12.
  * The 8 new files are codiffusion daa0505's methanogen growth-block analyses
    (`methane_nongrowth/atp_wl/growth_block/`).
  * Each reads the 20-hit build and only records or prints p. They are classified as
    non-consumers in `NOT_CONSUMERS`.
  * 0 files are unclassified.
* **Probability scan** over every reweighted genome of the two BV-BRC-only E3 builds:

| | codiffusion | EmilyKin |
|---|---|---|
| reweighted genomes | 1,034 | 1,410 |
| features (mixed E3 build in brackets) | 2,996,285 (5,808,409) | 3,840,554 (7,646,149) |
| replaced features (mixed) | 2,758,437 (2,758,484) | 3,423,698 (3,423,718) |
| replaced features whose value changed | 2,708,590 | 3,358,942 |
| E3 range | [0.0008, 1] | [0.00036, 1] |
| features with p = 0 (mixed) | 0 (47) | 0 (20) |
| minimum k/n of the replaced features | 0.0127 | 0.0130 |
| crossings of 0.5, k/n ≥ 0.5 -> E3 < 0.5 | 489,506 | 605,526 |
| crossings of 0.5, the other way | 71,890 | 75,090 |
| genomes with any 0.5 crossing | 1,033 of 1,034 | 1,405 of 1,410 |
| CDS / feature mismatches | 0 | 0 |

### 2. Community gapfill penalties (`penalties.py --builds bvbrc` -> `results_bvbrc/penalties.json`)

The method is the mixed run's (below). These parts are identical to the mixed run:
* the mapping (88 day-174 members, 0 unmapped, 41 reweighted; all days 1,718 ASVsets, 646 reps
  reweighted);
* the d100 clade ranks;
* the template classes (all equal to production);
* sets A and B.

* **Drafts.** 88 × 2 were built with 0 errors.
  * The two drafts of each member have the same reaction ids and 0 structural differences.
  * 3,679 GPR strings differ only in term order.
  * 39,692 reaction p values differ, all of them in the 41 reweighted members.
* **The BV-BRC-only k/n drafts against the mixed k/n drafts.**
  * Gene-bearing reactions fall from 74,802 to 74,187 (−0.8%).
  * 57 of the 88 members lose some (at most 43). None gains.

| set | candidates | changed | factor up / down | Δfactor median [p05, p95] (range) | relative change median [p05, p95] (range) | cheapest-member change (unique / tie broken / tie created / other) |
|---|---|---|---|---|---|---|
| A: demo's top 4 | 32,099 | 417, all in Burkholderiales.1 (417 of its 418 member-scored) | 275 / 142 | +0.0022 [−0.018, +0.106] (−0.966 … +0.342) | +0.16% [−0.9%, +8.5%] (−49.1% … +30.6%) | 0 / 47 / 1 / 7 of 8,249 |
| B: top 4 reweighted (same members) | 30,597 | 1,052 (417 / 158 / 99 / 378) | 784 / 268 | +0.0052 [−0.011, +0.207] (−0.973 … +0.749) | +0.39% [−0.8%, +20.0%] (−49.3% … +74.9%) | 59 / 59 / 1 / 18 of 7,870 |
| C: all 88 day-174 members (single-member pricing) | 734,632 | 8,628 (41 members; median 188 per member, 75-428; 0 in the 47 n = 1 members) | 6,642 / 1,986 | +0.042 [−0.146, +0.434] (−0.973 … +0.949) | +3.2% [−10.0%, +41.8%] (−49.3% … +94.9%) | 12 / 336 / 1 / 55 of 8,525 |

* **Set C.** The changed candidates are 1.17% of all single-member candidates and 92.8% of the 9,298
  member-scored candidates of the 41 reweighted members.
* **Equivalence.** For each of the 7 members of sets A and B, the single-member penalty vector equals
  its community-compartment vector in every candidate (stock, k/n, E3).
* **Members priced differently** (demo metric), set A: stock 238, k/n 408, E3 443 reactions.
* **Largest changes.**
  * `rxn60858` in Burkholderiales.1: p 0.103 -> 1.0, penalty 7,586 -> 4,000.
  * Methanobacterium.3 `rxn09333`: p 1.0 -> 0.350, penalty 2,008 -> 3,312.
* **Unique cheapest-member changes.** These are the same 12 reactions, with the same members, as
  in the mixed run. For example:
  * `rxn01094`: Sphingopyxis.1 -> Ochrobactrum.1;
  * `rxn06078`: Burkholderiales.1 -> Acinetobacter.1.
* **Prototype caveat.** It is unchanged: 73 candidates (rxn02028, rxn09296) have negative base
  penalties, and 31 of them are among the changed ones.
* **Community gapfilling.** No solution in 6 of 6 attempts (sets A and B × stock, k/n, E3). Each
  attempt stopped at "No gapfilling solution found before filtering for Wolfe activating bio1",
  after 46-83 s.
* **Run.** 1,053.8 s on 8 workers.

### 3. prFBA flux (`flux.py --builds bvbrc` -> `results_bvbrc/flux.json`, `results_bvbrc/part3/`)

**Sample.** The sample is identical to the mixed run's (`flux.py` checks it:
`_meta.sample_equals_saved_sample_of`). The pool sizes and the per-genome steps are the same.

* **Drafts.** 0 structural differences between the k/n and E3 drafts in 90/90, plus 7,737
  order-only GPR strings.
* **Gapfill.**
  * All 90 grow. The median is 44 new reactions for n ≥ 2, and 20 reversals were applied in all.
  * 60 draft reactions in 2 genomes lost their GPR (codiffusion `5aed6ebc…` 44, EmilyKin
    `96e917ee…` 16). All 60 are in those genomes' `gapfill.new_rxn`, so their p is 0 in both arms.
  * ATP correction removes many more draft reactions in the same 2 genomes. They are absent from
    the gapfilled models: 789 in `5aed6ebc…` (89 in the mixed run) and 287 in `96e917ee…` (593).
  * In codiffusion `6dbf54df…`, which lost 490 in the mixed run, it removes 2.
  * In the other 88 genomes `draft_rxn_net_removed` is at most 2 (median −5, min −48).
* **Runs.**
  * 2,037.8 s wall on 8 workers; 54-246 s per genome; peak worker RSS 1.91 GB.
  * The k/n tuned re-solve on a fresh copy was identical in 90/90.
  * The objective decomposition error is at most 1.7e-15 (relative).

**Pooled n ≥ 2 (80 genomes).**
* In each genome p changes on a median of 1,026 of 1,095 gene-backed reactions.
* The per-genome mean Δp over the changed reactions has a median of −0.050 (mean −0.070).

| | tuned (ex_weight 1) | released (ex_weight 100) |
|---|---|---|
| growth changed | 0 / 80 | 0 / 80 |
| flux vector identical | 2 / 80 | 2 / 80 |
| re-routed share of internal flux, per genome: median [p10, p90] (max) | 0.0174 [0.0006, 0.1124] (0.2645) | 0.0080 [0.0011, 0.0757] (0.2071) |
| re-routed share, pooled | 0.0395 | 0.0247 |
| Jaccard of active reaction sets, median (min) | 0.957 (0.870) | 0.973 (0.927) |
| genomes with any exchange-flux change | 56 | 16 |
| exchange L1 change, median (max) | 1.37 (75.4) | 0 (18.9) |
| total exchange flux, sum k/n -> E3 | 21,360.4 -> 21,276.2 | 20,443.8 -> 20,443.8 |
| total uptake, sum | 10,322.8 -> 10,306.1 | 10,093.2 -> 10,093.1 |
| gene-backed flux, sum | 62,916.6 -> 62,843.9 (−0.12%) | 65,559.1 -> 65,702.6 (+0.22%) |
| share of flux through p-changed reactions, median | 0.798 -> 0.798 | 0.798 -> 0.798 |
| objective Δ, median [p10, p90]; E3 higher / lower | +21.5 [−30.9, +90.9]; 59 / 21 | +24.0 [−38.3, +101.2]; 57 / 23 |
| of which repricing alone, median | +23.8 | +24.2 |
| E3 optimum's excess cost at k/n prices, relative: median (max) | 0.0007 (0.056) | 0.0000 (0.0043) |
| k/n optimum's excess cost at E3 prices, relative: median (max) | 0.0014 (0.036) | 0.0000 (0.0005) |
| flux-weighted mean p of gene-backed flux (own arm), median | 0.864 -> 0.817 | 0.837 -> 0.807 |

**Per deployment (tuned)**

| | re-routed share, median | pooled | exchange changed | gene-backed flux, sum |
|---|---|---|---|---|
| codiffusion | 0.020 | 0.046 | 28 / 40 | 26,659 -> 26,643 |
| EmilyKin | 0.017 | 0.035 | 28 / 40 | 36,258 -> 36,201 |

Under the released objective the codiffusion and EmilyKin pooled shares are 0.034 and 0.018.

**By n stratum (tuned)**

| n | re-routed share, median | pooled | mean Δp, median |
|---|---|---|---|
| 2 | 0.027 | 0.068 | −0.142 |
| 3-4 | 0.018 | 0.032 | −0.047 |
| 5-9 | 0.014 | 0.013 | −0.050 |
| ≥ 10 | 0.012 | 0.050 | −0.008 |

**The tail.** No genome re-routes 5 points less than in the mixed run, and six re-route at least 5
points more. In order (tuned, mixed -> BV-BRC-only):
* EmilyKin `b30381c0…` (n = 2): 3.7% -> 24.9%;
* codiffusion `6dbf54df…` (n = 4): 4.3% -> 21.2%;
* codiffusion `9413dd54…` (n = 2): 6.1% -> 16.3%;
* codiffusion `37e4420b…` (n = 14): 2.7% -> 13.8%;
* EmilyKin `9fc55d66…` (n = 2): 0.2% -> 11.2%;
* codiffusion `dd34ba85…` (n = 4): 0.1% -> 8.1%.

`9fc55d66…` has the same tuned objective in both arms and zero excess cost at either price set, so
its re-route is an alternative optimum. `6dbf54df…` is a different model (ATP correction, above).
codiffusion `400fe0a1…` (n = 10) is still the largest at 26.4% (26.7% before).

**n = 1 controls (10).** Identical drafts and identical p. Every metric under both configurations
is identical.

### 4. What dropping the RefSeq features does to the models (`vocab.py` -> `results_bvbrc/vocab.json`, `results_bvbrc/part4/`)

**Question.** The same genome is built from the shipped (mixed) k/n genome object and from the
BV-BRC-only k/n genome object. How do the models differ? Only k/n is compared; E3 within each
vocabulary is part 3.

**Sample.** The 90 part-3 genomes (80 n ≥ 2, 10 n = 1). Here the n = 1 genomes are not controls:
they lose their RefSeq features too.

**How it was measured.**
* Each arm's k/n draft, gapfilled model and growth come from that build set's part-3 run, with the
  same draft rule, templates (c12run 21a8084), pins and gapfill:
  * MSATPCorrection;
  * MSGapfill on the KBase complete medium at uptake 100;
  * reversals applied.
* `vocab.py` reads both arms and compares them.
* Genes are compared by function string, because feature ids are renumbered between the builds.
  A merged genome has one feature per function, so the function is the gene's identity.
* **Checks** (all pass):
  * template class equal in 90/90;
  * the BV-BRC-only function set is a subset of the mixed one in 90/90;
  * the mixed k/n draft, rebuilt here from the shipped genome object, equals the saved part-3
    draft structurally in 90/90;
  * in the 10 genomes whose features are unchanged (no RefSeq CDSs in their source genomes),
    drafts, growth and gapfill are identical.
* **Run.** 265.5 s on 8 workers.

| | mixed k/n | BV-BRC-only k/n | change |
|---|---|---|---|
| features per genome, median (n ≥ 2; n = 1) | 6,156; 2,410 | 2,875; 1,521 | ratio median 0.494 (0.31-1.0); unchanged in 10 / 90 |
| draft reactions, total | 113,576 | 112,717 | 0 gained in 90 / 90; 859 lost (0.76%) in 79 genomes; median 7 per genome, max 43 |
| gene-backed draft reactions, total | 99,589 | 98,876 | −0.72%; per-genome ratio median 0.995 (min 0.965); 32 stay in the draft gene-less |
| GPRs changed (reactions gene-backed in both) | | | 4,557 of 98,876 (4.6%); median 44 per genome (4.1%), max 183; every one loses genes only |
| draft reaction p (k/n) changed (reactions in both) | | | 319 in all; median 1 per genome; a change is a fall (median Δp −0.25) |
| drafts that grow on the complete medium | 0 / 90 | 0 / 90 | gapfilling is always needed |
| gapfilled models that grow | 90 / 90 | 90 / 90 | |
| growth (complete medium, C ≤ 100), median | 2.756 | 2.756 | identical in 56 / 90; relative change median 0 [p10 −0.26%, p90 +0.24%] (−63.3% … +8.7%); > 2% in 6 |
| MSGapfill new reactions, total (median per genome) | 4,888 (48) | 4,908 (47) | +0.4%; change per genome median 0 [p10 −1, p90 +2] (−50 … +56); more in 33, fewer in 16, equal in 41 |
| reversals applied, total | 18 | 20 | |
| all reactions gapfilling added, incl. ATP correction, total | 5,405 | 5,467 | +1.1% |
| gene-less share of the gapfilled model's rxn reactions, median | 7.53% | 7.50% | |
| lost draft reactions that the gapfill adds back gene-less | | | 20 of 859 (12 in the MSGapfill solution) |

**What is lost.**
* 178 of the 859 lost reactions are exchanges, which follow lost transporters.
* The most frequent losses are:
  * iron(II) ABC transport `rxn05555` (39 genomes);
  * glutarate-semialdehyde dehydrogenase `rxn01729` (30);
  * hydrogenobyrinate a,c-diamide `rxn11544` (28);
  * choline, arsenate / arsenite and L-fucose transporters.
* A reaction is lost when the roles that put it in the mixed draft came only from function strings
  that no PATRIC feature of the source genomes carries (RefSeq product names). The PATRIC
  annotation of the same genomes supplies no role that the template maps to these reactions.

**The growth outlier.** codiffusion `dd34ba85…` (n = 4, gram_pos) falls from 2.44 to 0.90.
* It loses 30 draft reactions. 15 of them are its glycerol-3-phosphate acyltransferase set
  (`rxn08546`-`rxn08552`, `rxn10202`-`rxn10210`).
* In the mixed draft their only genes were RefSeq's "acyl-phosphate:glycerol-3-phosphate
  O-acyltransferase PlsY" and "phosphate:acyl-ACP acyltransferase PlsX" (empty `PATRIC_id`).
* The PATRIC annotation names the same proteins "... PlsY (EC 2.3.1.n3)" and "... PlsX (EC
  2.3.1.n2)", and the c12run template does not map those names to these reactions.
* The gapfill adds 3 of the 30 back, one of them an acyltransferase (`rxn10204`). Growth on the
  complete medium drops by 63%.
* Which loss costs the growth was not traced. EmilyKin `62f901a0…` (n = 8, gram_pos) loses the
  same 15 acyltransferases and gets the same `rxn10204` back, and its growth is unchanged (2.659 ->
  2.655).
* Those two genomes account for all 30 glycerol-3-phosphate acyltransferases among the 859 lost
  reactions.

**The other changes.**
* The other genomes whose growth changes by more than 2%:
  * codiffusion `37e4420b…` −11.7%;
  * EmilyKin `bfc80b47…` +8.7%;
  * EmilyKin `96e917ee…` +7.5%;
  * codiffusion `6dbf54df…` +4.0%;
  * EmilyKin `b30381c0…` +3.2%.
* The three large gapfill swings follow ATP correction, not the lost genes directly. They are in
  the genomes where it prunes hundreds of draft reactions in one arm:
  * `6dbf54df…` 98 -> 48 new reactions;
  * `5aed6ebc…` 22 -> 78;
  * `96e917ee…` 63 -> 36.

  Their drafts differ by only 7-20 reactions (part 3, gapfill).

## Mixed-vocabulary builds (2026-11-12, `results/`)

Run on 2026-11-12. The comparison is the shipped k/n builds (`*/synthetic_genomes_d100/`) against
the E3 reweighted copies (`*/synthetic_genomes_d100_e3/`, `../e3_reweight_summary.json`). Both
sides merged BV-BRC's PATRIC CDSs and the imported RefSeq CDSs, the rule before 2026-11-13. These
numbers describe those builds, not a current one; the BV-BRC-only run above replaces them.

The scripts that made `results/` predate `--builds` (as committed at prFBA 3ce12d0; the hashes are
in each JSON's `_meta.code`). `--builds mixed` selects the same builds and directories.

### Headline (mixed)

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
     six genomes is in them). Those regenerated builds were deleted on 2026-11-13, superseded by
     the BV-BRC-only builds.
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

### 1. Consumer census (`census.py --builds mixed` -> `results/census.json`)

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

### 2. Community gapfill penalties (`penalties.py --builds mixed` -> `results/penalties.json`)

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

### 3. prFBA flux (`flux.py --builds mixed` -> `results/flux.json`, `results/part3/{deployment}/{asv}.json`)

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

* **Parts 3 and 4 are not the deployments' own pipelines.**
  * They use codiffusion's draft rule on c12run templates, with a complete-medium gapfill, for both
    deployments.
  * EmilyKin's own builder (round-12 role layer, curated r13 templates) and the day media were not
    used.
  * EmilyKin template classes use the genome lineage only: there is no MiDAS table keyed by ASV
    md5.
* **Part 2 prices candidates, as the demo does.** No community gapfill solution exists to compare.
* **Part 4 compares k/n drafts only.**
  * Each vocabulary's E3 effect is that run's part 3.
  * The growth in part 4 is the maximal growth on the complete medium. It is not the tuned prFBA
    solve, whose growth sits at the 95% floor.
  * ATP correction responds strongly to small draft differences in a few genomes (part 3), so a
    per-genome gapfill difference is not all vocabulary.
* **The mixed run cannot be re-run as is.** Its E3 builds were deleted on 2026-11-13 (superseded).
  `--builds mixed` needs them back for `census.py --scan`, `penalties.py` and `flux.py`. Regenerate
  them with `reweight_build_e3.py --build-dir <shipped synthetic_genomes_d100> --out-dir <dir>` at
  prFBA 3ce12d0; the current code refuses unstamped sources. The `--post-only` and
  `--summary-only` modes, `census.py` without `--scan`, and part 4 work from the saved files and
  the shipped k/n builds.
* **The BV-BRC-only k/n builds live in the session scratchpad** (`bvbrc_only/kn_*`). Their commands
  are in `../bvbrc_rebuild_summary.json`.
* **Size.** These are git-ignored, and all can be regenerated:
  * `results*/part2_work/single/` (52 MB of per-member penalty vectors each);
  * `models/` (about 1.3 GB of SBML);
  * `models_bvbrc/` (about 1.6 GB, with part 4's rebuilt mixed drafts).
* **Pins.** The pinned exports live in the session scratchpad. `common.ensure_pins()` recreates
  them from the recorded commits.

## Reproduce

Each script is run under its own unit, the way the runs on 2026-11-12 and 2026-11-13 were:

```
systemd-run --user --unit=ds-e3-<name>-$(date +%s) --collect -p MemoryMax=40G -p MemorySwapMax=0 \
  -p WorkingDirectory=/home/freiburger/Documents/prFBA -E PYTHONDONTWRITEBYTECODE=1 -E PYTHONHASHSEED=0 \
  /home/freiburger/Documents/py_venv/bin/python -W ignore \
  /home/freiburger/Documents/prFBA/deployment_d100/downstream_e3/<script> --builds <set> <args>
```

`--builds` is required by every script. It selects the k/n and E3 builds and the results and models
directories (`common.BUILD_SETS`). It reaches spawned workers through the command line and
`DS_E3_BUILDS`.

| script and arguments | what it does | time: mixed (11-12) / BV-BRC-only (11-13) |
|---|---|---|
| `census.py --builds <set> --scan` | the sweep, the table and the probability scan of the set's E3 builds | about 45 s / 25 s |
| `penalties.py --builds <set> --workers 8` | part 2 | 983 s / 1,054 s |
| `penalties.py --builds <set> --post-only` | recomputes the cheapest-member block from the saved single-member vectors | seconds |
| `flux.py --builds <set> --workers 8` | part 3 (checks the sample against the other sets' saved samples) | 1,851 s / 2,038 s |
| `flux.py --builds <set> --summary-only` | re-summarises the saved per-genome files with the saved sample | seconds |
| `vocab.py --builds bvbrc --against mixed --workers 8` | part 4; needs both sets' part-3 runs | - / 266 s |
| `vocab.py --builds bvbrc --against mixed --summary-only` | re-summarises part 4's per-genome files | seconds |

* The model-building scripts are resumable: finished drafts and finished genomes are skipped.
  Each build set has its own models directory, so a run never resumes from the other set's
  drafts.
* The worker pools are `multiprocessing.Pool(maxtasksperchild)` with spawn.
  `ProcessPoolExecutor(max_tasks_per_child)` deadlocked under Python 3.12 once workers retired.

## Files

| file | role |
|---|---|
| `common.py` | build sets, paths, pinned exports, provenance, canonical GPR |
| `census.py` | part 1 |
| `penalties.py` | part 2 |
| `flux.py` | part 3 |
| `vocab.py` | part 4 (mixed vs BV-BRC-only k/n models) |
| `results/` | the mixed run (2026-11-12) |
| `results_bvbrc/` | the BV-BRC-only run (2026-11-13), same layout, plus `vocab.json` and `part4/` |
| `results*/census.json`, `census_sweep.json` | consumer table, sweep hits, probability scan |
| `results*/penalties.json` | mapping, draft comparison, sets A/B/C, equivalence, cheapest member, gapfill attempts |
| `results*/part2_work/` | the day-174 sets file, template assignment, draft log and comparison, smoke output (per-member vectors in `single/`, ignored) |
| `results*/flux.json` | sample rule, pooled / per-deployment / per-stratum summaries, per-genome table |
| `results*/part3/{deployment}/{asv}.json` | per genome: drafts, gapfill, p changes, prFBA metrics, pair metrics, nonzero fluxes |
| `results_bvbrc/vocab.json` | part 4: checks, per group and deployment summaries, per-stratum and per-genome tables |
| `results_bvbrc/part4/{deployment}/{asv}.json` | per genome: features, draft differences (lost reactions, GPR changes), draft growth, gapfill burden of both arms |
| `logs/` | run logs (git-ignored, as is every `*.log`); `*_bvbrc*.log` are the 2026-11-13 runs |

## Scope notes added after the verification
- The sweep covers the two local repositories. It skipped vendored libraries and `runs/` snapshots;
  those were checked separately: ModelSEEDpy's gapfiller reads probabilities only when
  `reaction_scores` is passed, MSCommunity has no probability reads, and `MSCommunity.gapfill`
  (which accepts `reaction_scores`) is never called.
- KBase-side readers of the uploaded k/n builds were not checked (for example
  KB-ModelSEEDReconstruction's `modelseedrecon.py:510` reads `rxn.probability`, which comes from
  annotation-ontology scores); the E3 builds are not uploaded.
