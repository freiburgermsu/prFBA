# Shared-Reference k-Fold Benchmark + Phenotype Arm — Implementation Spec

**Root:** `/home/freiburger/Documents/prFBA/region_validation/shared_reference/`
(`scripts/`, `data/`; everything under `data/` is a run artefact — add
`region_validation/shared_reference/data/` to `.gitignore` and whitelist only
`data/summary_*.json`, `data/controls.json`, `data/fold_plan.json`,
`data/panel_eligible.json`, `data/phenotype_summary.json`.)

**Interpreter:** `/home/freiburger/Documents/py_venv/bin/python` — the ONLY interpreter.
PICRUSt2 binaries require `source ~/opt/picrust2/activate.sh`, which activates a
*different* Python. Never import `picrust2` from py_venv and never import `modelseedpy`
from the PICRUSt2 env. All PICRUSt2 work is done by `subprocess` through
`bash -lc 'source ~/opt/picrust2/activate.sh && <cmd>'`.

**One-sentence purpose.** Make prFBA and PICRUSt2 predict the EC content of held-out
genomes from *the same* reference genome set, using *one* annotation (BV-BRC EC
projection) as both methods' gene content and the answer key, so that neither system
authors the truth; then ask whether the EC-level difference propagates to metabolic
phenotype through a fixed ModelSEED reconstruction.

---

## 0. Decisions (auditor conflicts resolved; every later section cites a D-number)

| # | Conflict | Decision | One-line justification |
|---|---|---|---|
| D1 | Exclusion on full-16S identity (task step C) vs per-amplicon identity (A1, A3) | **Per-amplicon**, folded into group construction (D2). Full-16S identity is an *additional* edge, not the criterion. | 27.4% of genome pairs sharing a byte-identical amplicon fall below 0.987 full-16S identity; the code being replaced (`score_genecap.exclude_self_record`) already filters on amplicon identity. |
| D2 | Phylum stratification (task step C) vs GroupKFold over similarity components (A1) | **GroupKFold k=10 over connected components** of the edge graph in §3. Phylum/domain balance is a soft tiebreak in bin-packing only. | 133 of 206 phyla hold <10 genomes and 36 are singletons; stratification is arithmetically impossible, and group CV is what makes the holdout honest. |
| D3 | Prune a panel-wide tree (task step C) vs rebuild per fold (A1, A2, A3) | **Rebuild per fold, de novo, from survivors only.** No tip pruning anywhere in the fold arm. | `hmmalign --mapali` enforces an HMM/MSA checksum (place_seqs.py exits 1 on a pruned .fna with the old .hmm), and a pruned tree still carries the held-out sequence's influence in surviving branch lengths — an information channel prFBA has no counterpart to. |
| D4 | Panel scale: keep 10k (task) vs raise to ≥27k bacteria (A2) | **Keep the 10k panel; make density an explicit axis** (D5). Record ≥27k as a future extension. | Only the 10k panel has in-silico amplicons and benchmark truth generated; 27k requires re-running in-silico PCR and amplicon dedup over 17k more genomes, outside this experiment's scope. |
| D5 | Density confound: ignore vs sweep | **Nested density sweep at 25/50/100% of each fold's surviving reference, BOTH methods on every subsample**, both F1-vs-density slopes printed beside the headline. | Measured on PICRUSt2's own reference, cutting 26,828→8,600 tips costs +0.0201 F1 (p=0.011) — 29% of the gap the old headline rests on; prFBA's slope is unmeasured, so only a two-slope report is interpretable. |
| D6 | Positive control gate "F1 ~ 1" (task step E) vs density-calibrated ceiling (A2) | **Replace the 1.0 gate with the measured ceilings in §7.** The no-prune run is relabelled a *plumbing* control. | castor's max-parsimony prunes the query tip and assigns the attachment node's reconstructed state (edge cost ∝ 1/L^0.5 on a 1e-5-floored pendant edge), so PICRUSt2 can never copy the nearest tip; its include-self ceiling is 0.934 at 8,600 tips, not 1.0. |
| D7 | Conditional F1 (task step D) vs risk-coverage (A1) | **Co-primary: population EC F1 difference AND a risk-coverage curve per method, each sweeping its OWN knob.** Conditional F1 only on the matched subset both methods answered, with its denominator printed. | On the shipped data alone, letting PICRUSt2 abstain on its own worst 24.9% by NSTI shrinks prFBA's advantage from +0.070 to +0.017; bare conditional F1 compares two differently-selected populations. |
| D8 | prFBA select-mode: freeze shipped (A1) vs nested CV / legacy headline (A3) | **Nested CV** over {`legacy`, `exact_tie`, `graduated`, `exact_tie_guards`} chosen on training folds only; shipped `exact_tie` reported as a labelled sensitivity arm; inner winner logged per fold. | The shipped default was chosen by maximising gene-capture F1 on these same 67,991 amplicons, so it is fitted on every fold's labels; variants are pure re-selections over the same hits, so nested CV costs one extra `select_representatives` pass per variant. |
| D9 | PICRUSt2 domain: force from truth (A1) vs preserve native routing (A2) | **Forced-by-truth-domain is primary; min-NSTI routing is a secondary arm** with per-fold routing accuracy vs `picrust2_comparison/data/domain_routing.csv`. | prFBA has no routing step to get wrong, and the two custom trees' branch-length scales are set by independent FastTree runs, so cross-tree NSTI comparison is not commensurable. |
| D10 | Phenotype model input: `ec_rxn` reactions (A3) vs EC→template-role name match (A4) | **EC→template-role name match** (A4's `genome_from_ec`), identically for every arm. | `ec_rxn` fans out 2.8 uncurated reactions per EC and bypasses GPRs; the role route goes through curated role→complex→reaction links, and it is still one shared map applied identically to all arms, which is A3's actual requirement. |
| D11 | Phenotype readout: growth/no-growth over a defined media panel (task step F) vs core ATP-yield vector (A4) | **Primary = 54-medium core ATP-yield vector (no gapfilling).** Secondary = complete-medium growth rate + reaction essentiality on complete-medium-gapfilled models. Structural = reaction-set Jaccard + gapfill burden. | Un-gapfilled drafts give `bio1 = 0.0` with every exchange open; complete-medium-gapfilled models grow on 0 of 167 Biolog minimal media; minimal-medium gapfilling produces models that fail a thermodynamic gate. The stated growth panel returns an all-zero vector for every arm. |
| D12 | Trait table: presence/absence (forced by BV-BRC) vs counts | **Binary presence, `mp` HSP retained**, with a mechanical crash-column screen (§5.4). | Binarising cost exactly 0.0000 paired ΔF1 on 40 queries (identical recall/precision), but all-zero and some singleton columns crash castor's `mp`; dropping *only* the crashing columns preserves the default method without shrinking the predictable vocabulary wholesale. |
| D13 | Taxonomic-oracle baseline "same-genus donor" (A1) | **Oracle is defined per blocked rank:** under leave-genus-out it is a random surviving **same-family** donor plus the max-F1 surviving genome (copy-a-genome ceiling); under leave-family-out, same-order. | Under leave-genus-out no same-genus survivor exists by construction (genus is a graph edge), so A1's genus donor is undefined in the primary arm. |
| D14 | prFBA reference restriction: post-hoc hit filter vs re-alignment | **Re-align per fold against a restricted md5→seq DB** via a new `--ref-db` flag (§6.2), with `--k2 200 --gpu-cand 400`. `exclude_self_record` is used as a *canary assertion*, never as a filter. | Post-hoc filtering of 112k-DB top-20 hits reproduces the k2 truncation artifact (94.6% of self-exclusion abstentions) as a fold mechanism; re-alignment removes it by construction, and double-applying the exclusion would charge prFBA twice. |
| D15 | `preserve_branch_length` pruning helper (A2) | **Needed only in recipe-control arm A**, where the stock RAxML tree is pruned. The fold arm never prunes (D3), so `picrust2.util.prune_tree` is never called on a reference. | ete3's default pruning scales mean pairwise tip distance to 0.537x, which would fabricate the whole NSTI/novelty relationship — but only if pruning happened. |
| D16 | Template-visibility asymmetry (A4: 78.0% vs 90.9%) | **Mostly dissolved by the design; kept as a reported check, not a correction.** | Under the shared reference PICRUSt2's trait columns *are* the 2,236-EC BV-BRC universe, so both methods emit the same namespace; per-arm template-visible fraction is still logged per amplicon. |
| D17 | Marker (16S copy-number) table | **Written by the builder for set-equality bookkeeping; NOT used in the primary run.** | `hsp.py` is called directly with `--observed_trait_table` + `-n`, so NSTI comes from the EC run and no abundance/metagenome step (and no synthetic copy-number table) ever enters the comparison. |
| D18 | Phenotype ATP denominator: all held-out genomes vs only those with a non-degenerate truth vector | **Genomes whose TRUTH ATP vector is identically zero are UNSCOREABLE for ATP** (`meta.n_unscoreable_atp`, their own denominator). The contract-named primary (`arms[*].atp_hamming`, `contrasts[*].delta_atp_hamming` / `wilcoxon_p` / `n_pairs` / `q_bh`, every calibration rung, and the EC-ΔF1-vs-phenotype Spearman) is the informative subset; the all-genome value is retained beside it as a labelled sensitivity (`atp_all_genomes_sensitivity`, `delta_atp_hamming_all_genomes`). | On a silent genome `atp_hamming` has no truth signal to match: every arm, every ladder rung and the size-matched null score 0, and the only way to score above 0 is to predict any capability at all, so the metric is a pure precision penalty. This is the same cell `score_folds.py` already excludes as `unscoreable_empty_truth`. Measured on the 8-genome as-shipped panel (2 silent, no leak in that fixture): Δ 3.625 → 1.667, Wilcoxon p 0.203 → 0.625. On the leak-fixed 12-genome e2e panel (4 silent): all-genome sensitivity Δ −0.75 at p 0.281 (n=12) vs primary Δ −1.375 at p 0.1875 (n=8). (The pre-D19 e2e pair Δ 3.50 → 3.333, p 0.00195 → 0.0078, was computed on self-leaked rows and is retained only as the arithmetic of the swap, not as evidence.) Structure, growth and essentiality keep the full denominator — they are well defined on a silent genome. |
| D19 | Which prediction runs may feed a phenotype method arm | **Only runs whose reference EXCLUDES the test genomes** (`meta.mode` ∈ {`fold`, `as_shipped`}, i.e. `load_prfba`'s `PRFBA_MODES_OK`). `include_self` / `plumbing` are C2 ceiling controls and are skipped with a printed reason; two accepted runs may never claim the same amplicon md5 (hard abort), and every arm's EC F1 recomputed from the artefact is cross-checked against the cell's recorded F1. | `phenotype_fold.py` globbed `prfba/*/predictions.json` and filtered only density and `forced_answer`, so `prfba/incself_d100/` — which sorts after `prfba/f00_d100/` and covered the same 604 md5s — silently overwrote **every fold-0 prFBA prediction** with the ceiling control, whose selected rep is the source genome itself. Measured on the e2e smoke panel: 7 of 12 phenotype genomes got their own annotation (ATP Hamming 0, draft Jaccard 1.000, `base_nn1` = the truth genome), which is what produced the reported phenotype headline Δ −3.50 at p = 0.00195. With the filter, the same command gives Δ −1.375 at p = 0.1875 (n = 8 informative) and draft Jaccard 0.598 vs 0.607 — no significant difference. `load_picrust2` always had this filter. |
| D20 | A trait column dropped by the §5.4 screen with `dropped_reason == "full"` (prevalence `n`) | **Reinstated as a PICRUSt2 prediction in the scorer, resolved PER DOMAIN (per md5 via `best_domain` for a `--routing nsti` run, whose `meta.ref_name` is a list of both references).** `zero` stays predicted-0; `singleton`/`crash` (the aggressive retry) are NOT reinstated and are disclosed by C15 instead. | A trait carried by *every* reference tip reconstructs as present at every node, so castor would have predicted 1 for every query; scoring it 0 is a guaranteed false negative (measured: PICRUSt2 Archaea F1 0.6693 → 0.6997). The 'full' set is per-domain and very unequal (arc 7/10/25 vs bac 0/1/2 columns), so a run-level set would credit a bacterium-routed query with the archaeal reference's ECs; castor's behaviour on a retry-dropped column is not predictable in either direction, so repairing those would be invention. |
| D21 | An md5 the fold plan queried that is ABSENT from a launched run's prediction table | **A run defect, never a method abstention: hard control C14, and `score_folds.py` exits 2 after writing every artefact (the only control whose failure changes the exit status).** `cells.csv.gz` carries the reason in the new `abstain_reason` column. | It was previously laundered into coverage: deleting 60 of 561 md5s from one fold's bacterial run moved the primary delta −0.2610 → −0.2347 (PICRUSt2 coverage 1.000 → 0.9631) and deleting 60 from a prFBA run moved it to −0.2811, with every control still passing and the cells indistinguishable from genuine family-floor abstentions. Both shifts exceed the frozen SESOI of 0.02, and coverage is the quantity the headline rests on. |
| D22 | A cell whose required method run was never launched (e.g. `--domain bac` only, SPEC §11.1's own smoke recipe) | **`outcome = unscoreable_method_not_run` for EVERY arm**, so all denominators stay equal and the omission is visible per cell; plus a tripwire that refuses `headline.sentence` unless `F1(primary) − F1(picrust2) == delta_f1_population` to 5e-4. | Omitting those cells for PICRUSt2 only put the headline's two F1s on different populations (1,722 vs 1,922 cells): the printed values no longer subtracted to the printed paired delta (0.4345 − 0.6951 = −0.2606 against a reported −0.2443), and `skill_score`, `delta_f1_at_matched_coverage` and `accounting.per_fold` silently mixed the two. One equality catches the whole class. |
| D23 | Which reference vocabulary C10's structural recall ceiling is computed against | **Both, reported side by side: the bac+arc UNION (what prFBA's restricted `ref_db` really is) and the cell's OWN domain reference (all a forced-routed PICRUSt2 query can reach, D9), per fold and per domain.** | The union figure alone overstates PICRUSt2's archaeal ceiling by 2.6–4.9 recall points (fold 0 0.9945 vs 0.9686; fold 1 0.9923 vs 0.9547; fold 2 0.9988 vs 0.9494), i.e. in exactly the thin subpopulation §8.4 requires be reported separately; the bacterial gap is 0.0002–0.0008. Calling the union "the SHARED ceiling" was therefore prFBA's ceiling, not a shared one. |

---

## 1. Frozen configuration (recorded verbatim into every output JSON as `meta.frozen`)

```json
{
  "tau_amplicon": 0.987,
  "tau_full16s": 0.987,
  "k_folds": 10,
  "seed": 1729,
  "prfba": {"k2": 200, "gpu_cand": 400, "prefilter_k": 1000, "min_ref_len": 200,
            "force_backend": "cuda",
            "select_modes": {
              "legacy":           {"select_all": false},
              "exact_tie":        {},
              "graduated":        {"select_all_window": 0.005,
                                   "select_all_unconditional_window": 0.0},
              "exact_tie_guards": {"select_all_coherence": true,
                                   "select_all_contam": true,
                                   "select_all_cap": 6}},
            "forced_answer_knobs": {"t_family": 0.0, "family_floor": 0.0},
            "t_species": 0.987, "t_genus": 0.945, "t_family": 0.865,
            "family_floor": 0.865, "decay_tau": 0.01},
  "picrust2": {"hsp_method": "mp", "edge_exponent": 0.5, "min_align": 0.80,
               "placement_tool": "epa-ng", "max_nsti_reported": [null, 2.0],
               "calculate_nsti": true},
  "reference_qc": {"len_min": 1200, "len_max": 1800, "occupancy_min": 0.90,
                   "ec_min": 1},
  "fasttree_args": ["-nt", "-gtr", "-gamma"],
  "density_fractions": [0.25, 0.5, 1.0],
  "sesoi_delta_f1": 0.02,
  "label_prefix": "g_"
}
```

Nothing in this block may be changed after the first fold is scored. Any tuning
happens inside training folds only (D8) and is reported as a separate arm.

---

## 2. Immutable inputs (verified present; never written to)

| Path | Shape |
|---|---|
| `region_validation/data/selection_10k.json` | `{"meta": {...}, "genomes": [{genome_id, taxon_id, best_16s_md5, best_16s_len, n_full_length_16s, domain, phylum, class, order, family, genus, species}, ...]}` (10,000 records) |
| `region_validation/data/benchmark_truth.json` | `{"<genome_id>__<region>": {src_genome_id, src_taxon_id, region, domain, truth_lineage}}` (67,991 cells) |
| `region_validation/data/combined_amplicons.unique.fasta` | 60,461 records, header = bare amplicon md5 |
| `region_validation/data/amplicon_md5_fanout.json` | `{"<md5>": ["<genome_id>__<region>", ...]}` |
| `region_validation/picrust2_comparison/data/genome_function_sets.json` | `{"<genome_id>": {"ec": [...], "rxn_ec": [...], "rxn_role": [...], "n_pgfam": N}}` (49,110 genomes) |
| `region_validation/picrust2_comparison/data/function_maps.json` | `{"pgfam_ec": {...}, "pgfam_rxn": {...}, "ec_rxn": {...}, "_meta": {...}}` |
| `/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_seq.json` | `{"<md5>": "<sequence>"}` (519 MB) |
| `~/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files/bacteria/bac_ref/bac_ref.hmm` | stock bacterial 16S HMM (used only to *produce* alignments; leaks no gene content) |
| `.../archaea/arc_ref/arc_ref.hmm`, `.../bacteria/ec.txt.gz`, `.../bacteria/bac_ref/{bac_ref.fna,bac_ref.tre,bac_ref.model}` | recipe + density controls only |
| `~/Documents/ModelSEEDTemplates/templates/v7.0/{GramNegModelTemplateV7,GramPosModelTemplateV7,ArchaeaTemplateV6}.json`, `~/Documents/ModelSEEDDatabase` | phenotype arm |

Reused repo code (import, do not reimplement):
`picrust2_comparison/scripts/compare.py` → `metrics()`, `Agg`, `load_picrust2()`;
`picrust2_comparison/scripts/annotation_ceiling.py` → `acc_core()`, `load_reference_ec()`,
`load_reference_meta()`, `prf()`, `summarize()`;
`select_references.py` → `select_representatives()`, `DEFAULT_KNOBS`,
`gene_count_from_provider()`;
`region_validation/scripts/score_genecap.py` → `exclude_self_record()` (canary only, D14),
`ident_bin()`; `region_validation/scripts/_common.py` → `PATHS`, `wilson_ci`;
`picrust2_comparison/scripts/phenotype_arm.py` → every function listed in §9;
`align_hits.py` + `edlib_biopython_hits.py` → the aligner (patched per §6.2).

---

## 3. Panel eligibility and the similarity graph

### 3.1 Eligibility (applied before folding; `plan_folds.py --phase panel`)
Drop a panel genome if ANY holds:
1. absent from `genome_function_sets.json` (238 genomes);
2. `len(ec) < reference_qc.ec_min` i.e. empty EC set (161 genomes);
3. `best_16s_len` outside `[1200, 1800]` (306 are >1700 nt; these are operon
   fragments / mis-annotated features, not 16S);
4. HMM match-column occupancy < 0.90 against the stock domain HMM (computed in
   §4.1 during the panel-wide alignment);
5. the 16S sequence contains a character outside `ACGTN` after uppercasing
   (326 panel sequences contain non-ACGT). `N` is retained as-is; any other
   non-ACGTN character disqualifies the genome.
Expected survivor count ≈ 9,0xx (report exactly; the 9,601 figure in A2 predates
filters 3–5). Every dropped genome is removed from the panel, from all
references, from both methods' inputs, and from the truth — never scored.

The **same** eligibility filter is applied to a genome in its role as a reference
and in its role as a test organism, so neither method ever sees a genome the
other does not.

### 3.2 One identity definition (used by E3, E4, C6b and the novelty covariate)
**`identity` means exactly one thing everywhere in this spec:** the
`identity` field the aligner emits
(`edlib_biopython_hits._rescore` → `n_matches / aligned_len`, local
Smith-Waterman, BLASTN-like scoring `MATCH=2, MISMATCH=-3, GAP_OPEN=-5,
GAP_EXTEND=-2`), qualified by a coverage floor of `aligned_len / query_len >= 0.90`
(= prFBA's shipped `cov_min`). For E3 (full-16S vs full-16S) the coverage floor is
`aligned_len / len(shorter sequence) >= 0.90`. A hit below the coverage floor is
not an edge and is not a neighbour, at any tau.
This is precisely the quantity `select_references._tier` reads, so the exclusion
boundary is expressed in prFBA's native metric and the write-up says so (A3's
demand that one metric be named). `glocal`, `global` and MSA-derived identities are
NOT used anywhere.

### 3.2b Edges (union; `--phase plan`)
Two eligible genomes `u, v` are joined iff ANY:
- **E1 genus** — `genus(u) == genus(v)` and non-empty (from `selection_10k.json`);
- **E2 amplicon md5** — they share an amplicon md5 in any region
  (invert `amplicon_md5_fanout.json`);
- **E3 full-16S identity** ≥ `tau_full16s` (from §4.1's panel-vs-panel alignment);
- **E4 amplicon-level identity** ≥ `tau_amplicon` — some amplicon of `u` aligns to
  `v`'s 16S at `identity >= 0.987` under §3.2, or vice versa (from §4.1's
  amplicon-vs-panel alignment). This is the edge that makes D1 real.

Connected components of this graph are the **holdout groups**. A1 verified 2,986 components (1,592 singletons, largest 1,086 = 10.9% of panel) on E1∪E2
alone; E3∪E4 can only merge further, so re-report the component census.

### 3.3 Fold packing
Greedy size-descending bin packing of components into `k_folds=10` bins: sort
components by size descending, assign each to the currently smallest bin;
break ties by the bin whose domain composition and top-5-phylum composition is
furthest from the panel marginal (soft balance only, D2). **Abort** if the
largest component exceeds `1.5 * n_eligible / k_folds` — then the design cannot
be 10-folded and the orchestrator must be told, not silently re-folded.

### 3.4 Novelty ladder (the inferential axis, A1 amendment B)
- **leave-genus-out** — primary: the §3.2b graph as specified (genus is E1).
- **leave-family-out** — secondary: edges E1 replaced by `family(u)==family(v)`.
  The family graph has a 2,879-genome giant component (28.8%), so it CANNOT be
  10-folded. Run **5 holdout blocks** drawn from a stratified subsample of
  family-blocks totalling ≈2,000 test genomes; 5 reference rebuilds, not 1,235.
- **leave-order-out** — exploratory, 2 blocks, smoke-scale only.
Every endpoint is reported as a function of the rank blocked AND of the neutral
per-amplicon novelty covariate (§8.2). Claims are scoped to the rank blocked.

---

## 4. Component 1 — `scripts/plan_folds.py`

Panel QC, the similarity graph, the fold plan, and the pre-pass alignments that
feed all three. **Merged into one component** because the graph edges (E3, E4)
require the pre-pass alignment, and the fold plan requires the graph: three
scripts would pass two intermediates with no other consumer.

### 4.1 Phases

```
plan_folds.py --phase panel
  [--selection  region_validation/data/selection_10k.json]
  [--functions  region_validation/picrust2_comparison/data/genome_function_sets.json]
  [--md5-seq    /home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_seq.json]
  [--outdir     region_validation/shared_reference/data]
  [--limit-groups G --limit-genomes N]
      # SMOKE panel construction (NOT a random subsample): sample G connected
      # components of the FREE edge graph E1(genus) UNION E2(shared amplicon md5)
      # -- both are computable with no alignment -- take ALL their members, then
      # top up to N genomes with singleton components. Seeded by frozen.seed.
      # This guarantees the smoke panel CONTAINS genus-mates and md5 twins, so the
      # leakage assertions have something to detect (a flat "first N eligible"
      # subsample would make smoke-gate item 1 vacuous).
  [--label-prefix g_]          # set to "" ONLY to prove control C3 fires
  [--threads 32]
```
Writes `data/panel_eligible.json`, `data/panel_bac.fna`, `data/panel_arc.fna`
(unaligned, headers `g_<genome_id>`), `data/panel_bac.afa`, `data/panel_arc.afa`
(match-column-only aligned, §5.1), `data/panel_md5_seq.json`,
`data/panel_md5_id.json`, `data/panel_ec_universe.json`.

```
plan_folds.py --phase prepass
  [--amplicons region_validation/data/combined_amplicons.unique.fasta]
  [--outdir ...] [--k2 400] [--gpu-cand 600] [--prefilter-k 1500]
  [--inject-error 0.0] [--inject-indel-frac 0.1] [--inject-seed 1729]
```
`--inject-error R` (R in {0.005, 0.01}) additionally writes
`data/queries_err<RRR>.fna`: the SAME amplicon md5 headers with R per-base
substitutions and `R * --inject-indel-frac` indels applied to the sequence
(seeded, reference untouched). Both fold runners accept `--queries PATH` and read
this file instead of slicing the clean amplicon FASTA, so the error-injected arm
is the identical pipeline on mutated queries. Reference sequences are never
mutated.
Runs the patched aligner (§6.2) twice against `data/panel_md5_seq.json`:
(a) all panel full-16S sequences as queries (`--fasta data/panel_16s_queries.fna`,
headers = md5), (b) all amplicon md5s present in the eligible panel's cells.
Writes `data/prepass_full16s.json`, `data/prepass_amplicon.json` (raw aligner
output, `asv_top20_alignment_hits.json` schema), then distils
`data/neighbors.json`. **Asserts** for every query that the k2-th hit's identity
is `< tau` (i.e. the ≥tau neighbour set is not truncated); on violation it
re-runs that query set with `--k2` doubled, up to 3 times, then aborts.

```
plan_folds.py --phase plan [--rank genus|family|order|md5only] [--k 10] [--blocks 5]
  [--outdir ...]
```
`--rank md5only` builds the graph from E2∪E3∪E4 **without** the genus edge E1: the
**species-tier sensitivity arm**, in which surviving same-genus (and occasionally
same-species-complex) references remain, so prFBA's species/genus tiers actually
fire. Smoke scale only; written to `data/fold_plan_md5only.json`.
Writes `data/fold_plan.json` (+ `data/fold_plan_<rank>.json` for the secondary
ladders).

### 4.2 Output formats

`data/panel_eligible.json`
```json
{"meta": {"frozen": {...§1...}, "generated": "<iso8601>", "git_commit": "<sha>",
          "n_panel": 10000, "n_eligible": 0, "script": "plan_folds.py",
          "drop_counts": {"no_function_set": 0, "empty_ec": 0, "length": 0,
                          "occupancy": 0, "bad_chars": 0}},
 "genomes": {"<genome_id>": {"label": "g_<genome_id>", "taxon_id": 0,
    "md5": "<best_16s_md5>", "len": 0, "domain": "Bacteria|Archaea",
    "phylum": "", "class": "", "order": "", "family": "", "genus": "",
    "species": "", "n_ec": 0, "occupancy": 0.0, "template": "gramneg|grampos|archaea"}},
 "dropped": {"<genome_id>": "<reason>"}}
```
`template` is the pure function of true lineage from D10/§9.2 — computed here so
every downstream component reads the same value.

`data/panel_md5_seq.json` — `{"<md5>": "<uppercase sequence>"}` for eligible
genomes only. Keys are the real 16S md5s (so amplicon-vs-reference alignment is
unchanged); one md5 per eligible genome. **Assert** `len(set(md5s)) == n_eligible`
(verified 0 collisions on this panel).

`data/panel_md5_id.json` — **synthesised**, not subset from `16S_md5_ID.json`
(that file stores one header per md5 for the full BV-BRC DB, which may name a
non-panel genome). Format per key:
`"fig|<genome_id>.rna.1|   16S ribosomal RNA   [<species> | <genome_id>]"`
so `edlib_biopython_hits.parse_header` returns
`(species, genome_id, int(genome_id.split('.')[0]), 'fig|<genome_id>.rna.1|')`.
**Assert** round-trip for every entry at write time.

`data/panel_ec_universe.json`
```json
{"meta": {"n_genomes": 0, "n_ec": 0},
 "ec": ["1.1.1.1", ...],            // sorted; the FIXED trait/score vocabulary
 "prevalence": {"1.1.1.1": 0.0},    // fraction of eligible genomes carrying it
 "ec_index": {"1.1.1.1": 0}}        // column order for trait tables
```
EC strings are normalised by `annotation_ceiling.acc_core`-style stripping of an
`EC:` prefix; partial ECs containing `-` are **dropped** from the universe.

`data/neighbors.json`
```json
{"meta": {"tau_full16s": 0.987, "tau_amplicon": 0.987, "k2": 400,
          "truncation_violations": 0},
 "full16s": {"<genome_id>": [["<genome_id>", 0.9912], ...]},   // identity >= tau, self excluded
 "amplicon": {"<genome_id>": {"<region>": [["<genome_id>", 0.9931], ...]}},
 "md5_share": {"<genome_id>": [["<genome_id>", "<region>"], ...]},
 "rank1": {"<md5>": [["<genome_id>", 1.0], ...]}}   // full best-first list, for §8.2 novelty & NN1 baseline
```

`data/fold_plan.json`
```json
{"meta": {"frozen": {...}, "rank_blocked": "genus", "k": 10,
          "n_eligible": 0, "n_components": 0, "n_singletons": 0,
          "largest_component": 0, "generated": "<iso8601>", "git_commit": "<sha>"},
 "groups": {"<group_id>": ["<genome_id>", ...]},
 "genome_group": {"<genome_id>": "<group_id>"},
 "folds": [{"fold": 0,
            "test_groups": ["<group_id>", ...],
            "test_genomes": ["<genome_id>", ...],
            "reference_genomes": ["<genome_id>", ...],
            "test_cells": ["<genome_id>__<region>", ...],
            "test_md5s": ["<md5>", ...],
            "n_test_genomes": 0, "n_reference_genomes": 0,
            "n_reference_bac": 0, "n_reference_arc": 0,
            "domain_frac": {"Bacteria": 0.0, "Archaea": 0.0},
            "top5_phylum_frac": {"Bacillota": 0.0}}],
 "leakage_preflight": {"fold_0": {"md5_shared_with_reference": 0,
                                  "amplicon_identity_ge_tau": 0,
                                  "full16s_identity_ge_tau": 0,
                                  "taxon_id_shared": 0,
                                  "genus_shared": 0}}}
```
`reference_genomes` = all eligible genomes NOT in any test group of that fold.
Every `leakage_preflight` counter MUST be 0; a non-zero value aborts
`plan_folds.py` with exit 2 (the graph construction is then wrong).
`test_md5s` is the de-duplicated set of amplicon md5s of `test_cells`; because
E2 is an edge, no test md5 can appear among reference genomes' amplicons.

---

## 5. Component 2 — `scripts/build_reference.py`

Builds one PICRUSt2 custom reference directory per (fold, domain, density), and
the control references. **Control modes live here rather than in a seventh script**
because every control is "build a reference a different way and score it with the
same scorer".

```
build_reference.py
  --mode fold|recipe_control|include_self|plumbing
  --fold-plan  region_validation/shared_reference/data/fold_plan.json
  --panel      region_validation/shared_reference/data/panel_eligible.json
  --fold N                      # required for --mode fold
  --domain bac|arc|both         # default both
  --density 1.0                 # one of frozen.density_fractions
  --outdir region_validation/shared_reference/data/refs
  [--threads 32] [--fasttree-bin ~/opt/picrust2/bin/FastTree]
  [--stock-ref ~/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files]
  [--n-stock 960] [--n-holdout 500]     # recipe_control only
  [--force]
```

Reference directory name (and therefore the basename of all four files, which
`identify_ref_files` requires):
`refs/f<NN>_<dom>_d<DDD>/` e.g. `refs/f00_bac_d100/f00_bac_d100.{fna,tre,hmm,model}`
plus `traits_EC.txt.gz`, `traits_16S.txt.gz`, `manifest.json`.

### 5.1 Alignment recipe (D3, A2 amendment B)
1. Take the survivors' unaligned 16S (`panel_bac.fna` / `panel_arc.fna` subset).
2. `hmmalign --trim --dna --outformat Stockholm -o <sto> <stock domain HMM> <fna>`.
3. **Keep only match columns**: parse the `#=GC RF` line, keep positions whose RF
   character is not `.`; uppercase; convert `.` to `-`. Expected geometry:
   ~1,529 columns at ~2.4% mean gaps (bacteria), ~1,474 at ~1.5% (archaea).
   **Abort** if mean gap fraction > 0.10 or width > 2,000 — that is the 69%-gap
   failure mode A2 measured when insert columns are kept.
4. Write `NAME.fna` (aligned FASTA, `g_<genome_id>` headers).
5. `hmmbuild --dna --informat afa NAME.hmm NAME.fna` — **per fold**, from exactly
   this file (`hmmalign --mapali` enforces a checksum against it).
6. `FastTree -nt -gtr -gamma -log NAME.fasttree.log NAME.fna > NAME.tre`,
   capturing **stderr** to `NAME.fasttree.stderr`.
7. `NAME.model`: take the final `GTR_Opt` rates/freqs line from the `-log` file,
   renormalise the six rates so `gt = 1.0`, take `alpha` from the stderr line
   `Gamma(20) LogLk = ... alpha = ...` (it is NOT in the log), and write exactly
   `GTR{a/b/c/d/e/1.000000}+FU{fA/fC/fG/fT}+G4m{alpha}, noname = 1-<actual width>`.
   Never copy the stock `.model`: epa-ng reports every parameter as `(user)` and
   does not validate the range, so a copied model is a silent misspecification of
   a FastTree tree.

### 5.2 Density subsampling (D5)
For `--density f < 1.0`, subsample the fold's reference genomes to
`round(f * n)` with a seeded RNG (`seed = frozen.seed * 1000 + fold*10 + int(f*100)`),
stratified by phylum proportionally. The subsample is nested
(`d025 ⊂ d050 ⊂ d100`). The alignment, HMM, tree and model are rebuilt from the
subsample (D3 applies at every density).

### 5.3 Trait table (`traits_EC.txt.gz`)
Tab-separated, gzip. First column header literally `assembly` (matching the stock
`ec.txt.gz`); remaining headers are the **full-panel** EC universe from
`panel_ec_universe.json` in `ec_index` order — identical in every fold and every
density, so columns never drift. Row index = `g_<genome_id>` for every surviving
reference genome. Values are `0`/`1` (presence; D12). No `EC:` prefix.

### 5.4 Crash-column screen (D12, mechanical)
`hsp.py` needs a tree containing study tips, so the builder itself runs
`place_seqs.py` on 3 fixed probe queries (`refs/probe_queries.fna`: three 250 bp
V4 slices of the first three surviving reference sequences, degapped) to produce
`NAME.probe.tre`, then:
1. Drop columns with prevalence 0 or prevalence `n` among the surviving rows.
2. Run `hsp.py -t NAME.probe.tre --observed_trait_table traits_EC.txt.gz -m mp -e 0.5 -n`.
3. On the castor error `invalid 'ncol' value`, additionally drop columns with
   prevalence 1 or `n-1`, retry **once**.
4. On a second failure, **abort** the fold.
Record in `manifest.json`: `dropped_ec` (list), `dropped_reason` per EC,
`n_dropped`, `retry_used`. Every dropped EC remains in the truth and is scored as
**predicted-0** for PICRUSt2 — never as unscoreable (§8.1) — **except a
`reason == "full"` (prevalence-`n`) column, which the scorer reinstates as a
PICRUSt2 prediction per domain (D20)**; `singleton`/`crash` columns from the
step-3 retry are not reinstated and are disclosed by C15.

### 5.5 Marker table (`traits_16S.txt.gz`, D17)
Same layout; single column `16S_rRNA_Count` = `n_full_length_16s` from
`selection_10k.json` (floor 1). Written for set-equality bookkeeping only.

### 5.6 `manifest.json`
```json
{"name": "f00_bac_d100", "mode": "fold", "fold": 0, "domain": "bac",
 "density": 1.0, "n_tips": 0, "msa_width": 0, "msa_mean_gap_frac": 0.0,
 "hmm_leng": 0, "model_string": "GTR{...}, noname = 1-1529",
 "fasttree_alpha": 0.0, "fasttree_seconds": 0.0,
 "tip_labels_sha256": "<sha256 of sorted newline-joined tip labels>",
 "trait_cols": 0, "dropped_ec": [], "n_dropped": 0, "retry_used": false,
 "marker_rows": 0, "set_equality_ok": true,
 "mean_pairwise_tip_distance": 0.0,
 "genus_monophyly_frac": 0.0, "family_monophyly_frac": 0.0,
 "reference_genomes": ["<genome_id>", ...],
 "frozen": {...§1...}, "generated": "<iso8601>", "git_commit": "<sha>"}
```
`set_equality_ok` is the assertion of §7-C4 computed at build time:
`set(fna ids) == set(tree tips) == set(trait rows) == set(marker rows)`.
`genus_monophyly_frac` / `family_monophyly_frac` = fraction of panel genera /
families (with ≥2 surviving tips) recovered as a monophyletic clade — reported as
quantitative evidence the tree is not the handicap (A1's tree-quality request).

### 5.7 Control modes
- `--mode recipe_control` (**the decisive control**, A2): take `--n-stock` genomes
  from the stock bacterial reference, hold out `--n-holdout` of them, and build
  **three** references from the identical retained genomes, traits and queries:
  **arm A** = stock MSA + stock RAxML tree pruned to the retained tips with
  `preserve_branch_length=True` (the only place pruning happens, D15) + stock
  `.model`; **arm B** = stock MSA (`bac_ref.fna` is ALREADY aligned — §5.1 steps
  1–3 skipped) + §5.1 steps 5–7; **arm C** = the same retained sequences
  **degapped** and run through the FULL §5.1 recipe, steps 1–7 (hmmalign against
  the stock domain HMM, `#=GC RF` match columns, per-reference `hmmbuild`,
  FastTree, derived `.model`) — i.e. byte-for-byte the recipe every fold
  reference uses. Arm C therefore reproduces the fold geometry (measured 1,529
  columns at 2.2% mean gaps, vs the stock 1,578 at 4.5%) and its `.model` suffix
  is `1-1529`, not `1-1578`. Traits = stock `ec.txt.gz` (counts, identical for
  all three arms), truth = the same table. Queries = in-silico V4 slices of the
  held-out **stock-aligned** sequences (degapped), identical for all three arms,
  so every contrast is paired on the same queries. Output dirs `refs/recipe_A/`,
  `refs/recipe_B/`, `refs/recipe_C/`, plus `refs/recipe_control_queries.fna` and
  `refs/recipe_control_truth.json` (`{"<label>": ["<ec>", ...]}`).
  Arm A also writes `mean_pairwise_tip_distance` before and after pruning for the
  §7-C5 assertion. Each manifest records `recipe_alignment` and
  `recipe_tree_model` so the §7-C1 decomposition is machine-readable.
  **Decision (arm C added 2026-10-13):** without it both scored arms shared one
  byte-identical PICRUSt2-produced MSA (`sha256(recipe_A.fna) ==
  sha256(recipe_B.fna)`), so C1 priced only the tree+model half while §10 claimed
  it priced the alignment half too; arm C costs ~60 s at 960 tips (17 s hmmalign
  + 30 s FastTree), which is cheaper than narrowing the claim.
  Arm C re-aligns the degapped stock residues rather than the raw pre-ssu-align
  sequences because only the aligned `bac_ref.fna` is distributed; ssu-align's
  upstream curation of raw sequences is therefore outside C1's coverage.
- `--mode include_self`: build the fold reference with the fold's test genomes
  **present** (i.e. the whole eligible panel at that density) — used for the
  include-self ceiling control (§7-C2).
- `--mode plumbing`: identical to `include_self` but restricted to the smoke
  panel; used for the relabelled no-prune plumbing check (§7-C9).

**Control-mode directory names** (same `NAME/NAME.{fna,tre,hmm,model}` rule):
`refs/recipe_A/`, `refs/recipe_B/`, `refs/recipe_C/`, `refs/incself_<dom>_d<DDD>/`,
`refs/plumbing_<dom>_d<DDD>/`. `--mode recipe_control` additionally writes
`refs/recipe_control_queries.fna` and `refs/recipe_control_truth.json`
(`{"<query_label>": ["<ec>", ...]}`), and the recipe manifests carry
`mean_pairwise_tip_distance_before_prune` / `_after_prune` for C5.

---

## 6. Components 3 and 4 — the two fold runners

### 6.1 `scripts/run_picrust2_fold.py`

```
run_picrust2_fold.py
  --mode fold|recipe_control|include_self|plumbing
  --ref-dir    .../data/refs/f00_bac_d100
  --outdir     .../data/picrust2/f00_bac_d100
  [--fold-plan .../data/fold_plan.json]   # required for --mode fold|include_self
  [--panel     .../data/panel_eligible.json]
  [--fold 0] [--domain bac] [--density 1.0]
  [--amplicons region_validation/data/combined_amplicons.unique.fasta]
  [--queries PATH]   # explicit query FASTA; overrides fold-derived slicing.
                     # REQUIRED for --mode recipe_control (refs/recipe_control_queries.fna)
                     # and used by the error-injected arm (data/queries_err005.fna)
  [--truth PATH]     # explicit {"<query_id>": ["<ec>", ...]} truth; REQUIRED for
                     # --mode recipe_control (refs/recipe_control_truth.json)
  [--traits PATH]    # override <ref-dir>/traits_EC.txt.gz (recipe control uses stock ec.txt.gz)
  [--processes 30] [--chunk-size 5000] [--routing forced|nsti]
```
`--routing nsti` runs BOTH domain references for every query and keeps the
lower-NSTI prediction. It doubles the PICRUSt2 run count, so it is **restricted to
`--density 1.0`** and is a secondary arm only (D9).

Procedure (all PICRUSt2 calls through `bash -lc 'source ~/opt/picrust2/activate.sh && ...'`):
1. Query FASTA `study_seqs.fna` = the fold's `test_md5s` whose cells are in this
   `--domain` under `--routing forced` (truth domain from `benchmark_truth`), or
   all `test_md5s` under `--routing nsti` (then both domain references are run and
   the lower-NSTI prediction wins per md5). Headers = bare md5 (single field).
2. `place_seqs.py -s study_seqs.fna -r <ref-dir> -o placed.tre -p P
   --min_align 0.80 --intermediate intermediate --verbose`
3. `hsp.py -t placed.tre --observed_trait_table <ref-dir>/traits_EC.txt.gz
   -o EC_predicted.tsv.gz -m mp -e 0.5 -n -p P --seed 1729`
4. **Split `hsp.py`'s output into the two files `compare.load_picrust2` expects.**
   `hsp.py -n` appends `metadata_NSTI` as an extra column *of the EC table*; a
   naive hardlink would make `compare.load_picrust2`'s `df.to_numpy() > 0` treat
   `metadata_NSTI` as a predicted EC for every query with NSTI > 0, and its
   `r["best_domain"]` lookup (a plain `[]`, not `.get`) would raise `KeyError`
   because `best_domain` is produced only by `picrust2_pipeline.py`'s two-reference
   wrapper. Therefore write, in `--outdir`:
   - `combined_EC_predicted.tsv.gz` — `hsp.py`'s table with the `metadata_NSTI`
     column **removed**; index column header `sequence`, remaining headers are the
     EC strings of the fold's trait table (no `EC:` prefix).
   - `combined_marker_predicted_and_nsti.tsv.gz` — synthesized, tab-separated,
     columns exactly `sequence`, `metadata_NSTI`, `best_domain`,
     `closest_reference_genome`. `best_domain` = `--domain` under
     `--routing forced`; under `--routing nsti` it is the argmin-NSTI domain.
     `closest_reference_genome` is the `g_<genome_id>` label of the nearest known
     tip in `placed.tre` by patristic distance.
   With those two files present, `compare.load_picrust2(outdir, ec_index)` is
   reused verbatim.
5. **Assertions** (abort on failure, exit 3): §7-C3 (unknown-tip round trip),
   §7-C4 (set equality), §7-C7 (placement-loss accounting is recorded, not gated).

Outputs:
- `predictions.json`
```json
{"meta": {"fold": 0, "domain": "bac", "density": 1.0, "ref_name": "f00_bac_d100",
          "routing": "forced", "frozen": {...}, "n_query": 0,
          "n_placed": 0, "n_min_align_dropped": 0,
          "min_align_dropped_md5": ["<md5>", ...],
          "n_unknown_rows": 0, "unknown_rows_match_queries": true,
          "picrust2_version": "2.6.3", "wall_seconds": 0.0,
          "generated": "<iso8601>", "git_commit": "<sha>"},
 "pred": {"<md5>": {"ec": ["1.1.1.1", ...], "nsti": 0.0,
                    "best_domain": "bac", "closest_reference_genome": "g_<genome_id>"}}}
```
  `ec` is the sorted list of EC strings with value `> 0` (values are integral —
  verified — so the `>0` rule is unambiguous and needs no threshold sweep).
- `min_align_dropped.json` — `{"<md5>": {"len": 0, "aligned_len": 0, "regions": ["V4", ...]}}`
- `versions.txt`, `place_seqs.log`, `hsp.log`.

### 6.2 `scripts/run_prfba_fold.py` + the aligner override (D14)

> **BUILD ORDER.** The aligner patch below is **component 0** and must land
> first: `plan_folds.py --phase prepass` already depends on `--ref-db` /
> `--ref-id-db` / `--midas-csv`, so nothing else in this spec can run until the
> patch is in. It is a backwards-compatible change to two existing top-level
> files, not a new script.

**Patch to existing code (required, minimal, backwards-compatible):**

*`edlib_biopython_hits.py`*
- `def load_references(min_ref_len: int, db_path: str | None = None)` — use
  `db_path or DB_MD5_SEQ` instead of the module constant.
- `def enrich_and_write(..., id_db_path: str | None = None, midas_csv: str | None = None)`
  — use `id_db_path or DB_MD5_ID`, and `load_midas(midas_csv or TAXONOMY_CSV)`.
- Add CLI flags `--ref-db PATH`, `--ref-id-db PATH`, `--midas-csv PATH` to its own
  `main()` and thread them through.

*`align_hits.py`*
- Add `ap.add_argument("--ref-db", default=None, help="md5->seq JSON replacing the default BV-BRC reference DB")`
  and `ap.add_argument("--ref-id-db", default=None, help="md5->BV-BRC header JSON matching --ref-db")`.
- Immediately after `args = ap.parse_args()`:
  `if args.ref_db: P.DB_MD5_SEQ = args.ref_db` /
  `if args.ref_id_db: P.DB_MD5_ID = args.ref_id_db`.
  (Both are read at call time from the module namespace by `load_references` and
  `enrich_and_write`, so rebinding `P.DB_MD5_SEQ` is sufficient and is the
  documented mechanism; the explicit parameters above exist so the behaviour is
  not an accident of import order.)
- Record `ref_db`, `ref_id_db` in `method_selection.json["params"]`.
- **Assert** `args.gpu_cand >= args.k2` and exit 2 otherwise: the CUDA path only
  Biopython-rescores `gpu_cand` candidates, so `k2=200` with the default
  `gpu_cand=100` silently truncates at 100.

**Runner:**
```
run_prfba_fold.py
  --fold-plan .../data/fold_plan.json
  --panel     .../data/panel_eligible.json
  --fold 0 --density 1.0
  --ref-manifest .../data/refs/f00_bac_d100/manifest.json  # repeatable (bac+arc)
  --amplicons region_validation/data/combined_amplicons.unique.fasta
  --outdir .../data/prfba/f00_d100          # use .../f00_d100_forced with --forced-answer
  [--k2 200] [--gpu-cand 400] [--prefilter-k 1500] [--force-backend cuda]
  [--select-modes legacy,exact_tie,graduated,exact_tie_guards]
  [--forced-answer] [--queries PATH] [--workers 32]
```
`--forced-answer` merges `frozen.prfba.forced_answer_knobs`
(`t_family=0.0, family_floor=0.0`) into every select-mode and writes to a separate
`--outdir` (convention `..._forced`); keys under `pred` become
`<select_mode>_forced`. The two runs must never share an output directory.
`--queries PATH` overrides the fold-derived query FASTA (error-injected arm).

**Forced-answer invariant (assert per fold).** With `t_genus` and `t_species`
unchanged, lowering only `t_family` and `family_floor` to 0.0 can change a record's
outcome **only** where the shipped arm abstained (`best_identity < 0.865`): every
such cell becomes a family-tier call (`cap_family=2`, `contam_family=0.20`), and
every cell the shipped arm answered must produce a **byte-identical** EC set.
Assert this equality per fold (0 mismatches); it is what makes the post-hoc
`rank1_identity` sweep over the forced arm a clean risk-coverage curve rather than
two different methods. Note forced coverage may still be < 1.0: stage-0 reliability
(`cov_min=0.90`, `swed_max=0.02`) can empty a record regardless of the tiers.
Report that residual as **prFBA's coverage floor**.
Procedure:
1. Build `ref_db.json` = `{md5: seq}` for the **union** of the fold's surviving
   reference genomes across the `--ref-manifest`s given (so prFBA's reference is
   *exactly* the union of the two PICRUSt2 trees' tips at this density, satisfying
   §7-C4). Build `ref_id_db.json` by subsetting `panel_md5_id.json`. Write
   `midas_empty.csv` containing only the header
   `seq,Kingdom,Phylum,Class,Order,Family,Genus,Species,rel_ab`.
2. Write `queries.fna` = the fold's `test_md5s` (headers = md5).
3. `align_hits.py --fasta queries.fna --ref-db ref_db.json --ref-id-db ref_id_db.json
   --taxonomy-csv midas_empty.csv --outdir hits/ --k2 200 --gpu-cand 400
   --prefilter-k 1500 --force-backend cuda` → `hits/asv_top20_alignment_hits.json`.
4. For each `--select-modes` variant, call
   `select_references.select_representatives(record, knobs=<variant knobs merged
   with forced_answer_knobs when --forced-answer>, gene_provider=<PGFam-cache
   provider from prfba_selection.py>, gene_count_fn=SR.gene_count_from_provider(provider))`
   on every record, exactly as `prfba_selection.py` does — import its provider
   class rather than reimplementing it.
5. Predicted EC set per md5 = union over selected reps of
   `genome_function_sets[rep]["ec"]`, restricted to `panel_ec_universe`.
   **Assert** every selected rep is in the fold's reference set (0 violations).
6. **Canary assertions** (abort on failure, exit 3): §7-C6a–d, §7-C8.

Outputs — `predictions.json`
```json
{"meta": {"fold": 0, "density": 1.0, "ref_name": ["f00_bac_d100","f00_arc_d100"],
          "n_reference_md5": 0, "forced_answer": false,
          "k2": 200, "gpu_cand": 400, "backend": "cuda", "frozen": {...},
          "wall_seconds": 0.0, "generated": "<iso8601>", "git_commit": "<sha>",
          "assertions": {"exclude_self_n_removed_total": 0,
                         "midas_nonempty_records": 0,
                         "species_dedup_collisions": 0,
                         "gene_count_fn_calls": 0,
                         "reps_outside_reference": 0,
                         "k2_saturated_records": 0}},
 "pred": {"<select_mode>": {"<md5>": {
            "ec": ["1.1.1.1", ...],
            "reps": ["<genome_id>", ...], "n_reps": 0,
            "tier": "species|genus|family|below",
            "abstained": false, "abstain_reason": "",
            "best_identity": 0.0, "rank1_genome": "<genome_id>",
            "rank1_identity": 0.0, "rank1_rank": "genus|family|order|class|phylum|none"}}}}
```
`k2_saturated_records` counts records whose k2-th hit identity ≥ the family floor
(i.e. where truncation could still bite); report it, do not gate on it.
`rank1_rank` = the taxonomic rank at which `rank1_genome` and the test genome
first agree — the method-neutral novelty label (§8.2).

---

## 7. Positive controls and the numbers that falsify the setup

Every control is a mechanical assertion with an explicit falsifying value. All
results land in `data/controls.json`:
```json
{"<control_id>": {"status": "pass|fail", "measured": {...},
                  "gate": "<expression>", "falsifies": "<what a fail means>"}}
```

| ID | Control | Gate | Falsifying number |
|---|---|---|---|
| **C1** | **Recipe calibration** (`build_reference.py --mode recipe_control`, scored by `score_folds.py --mode recipe_control`): three arms over identical retained stock genomes, stock traits, stock truth and one shared set of ≥500 held-out in-silico V4 queries — **A** stock MSA + stock RAxML tree + stock model, **B** stock MSA + our FastTree/hmmbuild/derived model, **C** our full §5.1 recipe (hmmalign + RF match columns + hmmbuild + FastTree + derived model) on the same residues. Reported as a decomposition: ΔF1(C−A) = whole build-recipe dialect, ΔF1(B−A) = tree+model half with the MSA held fixed, ΔF1(C−B) = alignment half with the tree recipe held fixed (per query the three are exactly additive). | `abs(mean paired ΔF1(C−A)) < 0.01` (**gated on the full-recipe contrast**; B−A and C−B are reported, not gated. If no arm C exists the gate falls back to B−A and `controls.C1.measured.gated_contrast_covers == "tree_model_only"` says so; a built-but-unrun `refs/recipe_C` is a hard error, never a silent fallback.) | ΔF1(C−A) ≤ −0.01 ⇒ our build recipe handicaps PICRUSt2 relative to its own stock build, and every PICRUSt2 number must carry that band as an uncertainty; the B−A / C−B rows say which half is responsible. **Covers:** hmmalign + `#=GC RF` match columns vs the stock ssu-align MSA of the same residues, the per-reference `hmmbuild` that `place_seqs.py --mapali` re-aligns queries through, FastTree vs RAxML, and the derived vs stock `.model`. **Does not cover:** ssu-align's upstream curation of the raw stock sequences (only the aligned `bac_ref.fna` ships, so arm C starts from degapped stock residues), and the count-vs-binary trait axis (all three arms use the stock EC counts; the fold arm uses binary traits). Measured at n=40, 960 tips (mean F1 A 0.8221, B 0.8204, C 0.8233): B−A −0.0016 (p=0.64), C−B **+0.0029** (p=0.048), **C−A +0.0013 (p=0.62) → pass**. The alignment half is the only nominally significant term and it favours *our* recipe by 0.3 F1 points — 3.4x below the 0.01 band and 7x below the SESOI — so PICRUSt2 is not handicapped by either half of the build. |
| **C2** | **Include-self ceiling** (`--mode include_self`), reported as a *quantity* per method per density, not a pass at 1.0 (D6). Twin-free cells only (amplicon md5 maps to one genome); the twin cells are reported in `measured.prfba_twin_cells` against the analytic `source ∪ twin` ceiling. | **HARD:** prFBA source-recovery is MEASURED and `== 1.00` (an absent measurement is a FAIL, not a pass: it means the include-self run did not happen) (the source genome is among the selected reps). **SOFT (flag, never abort):** PICRUSt2 mean EC F1 `< 0.85`. | prFBA source-recovery < 1.00 ⇒ the restricted DB, the synthetic md5→header join, or the select call is wrong — abort. The PICRUSt2 threshold is deliberately SOFT and set below the 0.9343 measured at 8,600 tips because that number comes from a **different reference in a different annotation dialect** (PICRUSt2's GTDB tips, its own 3,111-EC traits); the only panel-dialect measurement is 0.8155 at 240 tips. A low value is therefore reported as the density- and dialect-conditional ceiling that caps every PICRUSt2 number, not treated as a setup failure to be tuned away. Twin cells are reported separately against the analytic `source ∪ twin` EC ceiling, never gated. |
| **C3** | **Tip-label round trip**, after every `hsp.py`. | `n_unknown_output_rows == n_query_sequences` exactly | Any mismatch ⇒ `castor_nsti.R`'s `read.table(...)$V1` typed numeric BV-BRC ids as numeric and deleted reference genomes from the NSTI known-set (3.61% of panel ids end in a trailing zero). This is why every label carries the `g_` prefix. Abort the fold. |
| **C4** | **Artefact set equality**, at build and at run. | `set(NAME.fna ids) == set(NAME.tre tips) == set(traits_EC rows) == set(traits_16S rows) == set(prFBA ref_db genomes)` | One differing id ⇒ reference genomes silently become study sequences. Abort. |
| **C5** | **Branch-length preservation**, recipe-control arm A only (D15). Arms B/C are built de novo and are recorded `gated: false` — before `==` after by construction there, so recording them instead of A made the gate unfalsifiable; every candidate is accumulated and the gate also fails when NO pruned arm was tested. | `abs(mean_pairwise_tip_distance(pruned)/mean_pairwise_tip_distance(full) − 1) < 1e-6` over a fixed 200-tip sample, on `recipe_arm == "A"`, with `n_pruned_references_tested >= 1` | ete3's default yields 0.537x and would fabricate the NSTI/novelty relationship. Abort. |
| **C6** | **Per-fold leakage**, four assertions. (a) no surviving reference shares an amplicon md5 with any test amplicon in any region; (b) `max` amplicon-level identity of any test amplicon to any surviving reference `< 0.987`; (c) no surviving reference shares `taxon_id` with a test genome; (d) `score_genecap.exclude_self_record(record, src_gid, src_tid, tau=0.987)` returns `n_removed == 0` for every record (D14 — canary, not filter). | all four counts `== 0` | Non-zero ⇒ the group graph is wrong and both methods got a free near-exact neighbour. Abort the fold. |
| **C7** | **Placement-loss accounting.** | reported, not gated | Every md5 dropped by `--min_align 0.80` is excluded from **both** methods' denominators, per region, and the count is printed. Scoring them as PICRUSt2 abstentions would charge PICRUSt2 for a shared data limitation. Such a cell's row keeps the informative `abstain_reason` sentinel but carries `outcome == unscoreable_min_align`, and C14 leg 3 therefore **does not count it**: counting it made leg 3 contradict leg 2 (which subtracts the same drops from the expected set) and voided an otherwise perfect run — one legitimate drop produced `headline.valid == false` and exit 2 with every number correct (measured: C7 total 1, all arms n=1921, picrust2 0.6955, delta −0.2611). The exemption is narrow: ONLY `unscoreable_min_align`. |
| **C8** | **prFBA leak-adjacent hooks.** | `midas_nonempty_records == 0`; `reps_outside_reference == 0`; `species_dedup_collisions` and `gene_count_fn_calls` reported | A non-empty MiDAS lineage would hand prFBA an external taxonomic prior for the held-out organism. Collisions/gene-count calls are expected to be 0 on a 10,000-distinct-species panel; non-zero is reportable, not fatal. |
| **C9** | **Plumbing (no-prune) control** — relabelled from task step E (D6). Measured over the PICRUSt2 plumbing run(s) only, and the gate now says so; prFBA's end-to-end wiring is covered substantively by C2's include-self run at coverage 1.000. | **PICRUSt2** produces a prediction for every smoke query; trait/md5 joins and domain routing resolve | Verifies wiring only. A pass licenses nothing about tree quality. |
| **C10** | **Structural recall ceiling**, against BOTH reference vocabularies (D23). | reported per fold **and per domain**: `structural_recall_ceiling` (bac+arc union) and `structural_recall_ceiling_picrust2_own_tree` (the cell's own domain reference), with `union_minus_own_domain` printed | ECs present in a held-out genome but carried by **zero** reachable reference genome. prFBA's restricted `ref_db` IS the union, so the union figure is its ceiling; a forced-routed PICRUSt2 query (D9) only reaches its own domain's tree and trait table, so the union figure overstates its archaeal ceiling by 2.6–4.9 points. Reported as two numbers, never as one "shared" number. |
| **C11** | **Phenotype identity control**, over a phenotype summary **BOUND to this run**. | two builds from the same EC set in the same batch give identical reaction sets and gapfill solution (Jaccard 1.000) **AND** `phenotype_summary.json` is bound to THIS cells table, fold plan, panel and answer key: `meta.binding` (written by `phenotype_fold.py`, contract 1) records the sha256 of `--cells`, `--fold-plan`, `--panel`, `--functions`, `--ec-universe`, of its own `per_genome.jsonl` and `identity_control.json`, the run directories it read and the genomes it scored, and `score_folds.py` re-checks every one of them (`accounting.phenotype_binding`). An unreadable, foreign or (under `--require-phenotype-binding`) unbound summary is a C11 **and** C12 **FAIL** | Measured 78/78 identical with Gurobi 13.0.1. A fail means degenerate optima are polluting every phenotype number. The binding clause closes the other half: the summary was auto-discovered by EXISTENCE beside `--cells` and trusted, so a DIFFERENT campaign's 12-genome summary dropped in there scored C11/C12 **pass**, and a truncated one made both controls **vanish** from `controls.json` (the whole block was wrapped in `except Exception: pass`), which is indistinguishable from a tree where the phenotype arm never ran. A summary written before contract 1 carries no `meta.binding`: it is reported `bound: false` with a stdout WARNING and keeps its own verdicts, and `--require-phenotype-binding` (which §11.2 passes) turns that into a failure. |
| **C12** | **Phenotype validity gate** (`phenotype_arm.sanity_gate`, `yield_check`). | **HARD** (§7.1 always listed it so; the row used to say only "report", which made it unfalsifiable): biomass and `ATPM` on water+protons both `== 0` and biomass per mmol carbon uptake below the calibrated ceiling — asserted as **(a)** the phenotype summary is BOUND to this run (C11's clause), **(b)** `meta.n_gated_out.truth == 0`, zero tolerance, and **(c)** `sum(n_gated_out) / meta.n_models <= 0.10` (`C12_MAX_GATED_FRACTION`). The per-arm gated-out counts are reported in `measured.n_gated_out` as before | The raw draft yields 1,000 ATP from water via `rxn15249_c0`+`rxn20583_c0` cycling `cpd19001/cpd19006`; `loopless_solution` does NOT remove it. The control was written as `control(True, ...)`: it could not fail for any input, in any tree — an injected `n_gated_out = {"atp_from_water": 17}` (37% of 46 models) passed, as did a phenotype arm whose identity control had already failed. Leg (b) is the silent-population rule this harness exists for: a gated TRUTH row makes every other arm of that genome `unscoreable`, so it removes the genome from every arm's denominator. Leg (c) separates "the gate is protecting the measurement" (a few degenerate drafts, reported) from "the pipeline is miscalibrated" (a third of all models gated). |
| **C13** | **Tree-quality descriptives** (not a gate). | reported | Genus/family monophyly fraction and the NSTI distribution of each custom reference vs stock, per fold. |
| **C14** | **Prediction-run completeness** (D21, D22) — the control whose failure makes `score_folds.py` exit **2** (any other HARD control failing exits **3**; see §7.1). | **HARD**, seven legs: (1) every launched run is self-consistent with its own meta (`len(pred) == n_placed == n_query − n_min_align_dropped` for PICRUSt2; for prFBA one record per query in **exactly the `meta.select_modes` mode set**, compared on base mode names because the forced runner suffixes its keys `_forced` — a missing or empty mode block is a FAIL, never a vacuous pass); (2) `set(fold_plan.test_md5s) − min_align_dropped` is covered by every launched run (domain-filtered for a forced-routing run); (3) zero cells scored with `abstain_reason ∈ {no_record, no_prediction, no_forced_record}` **and `outcome != unscoreable_min_align`** (C7); (4) `F1(primary) − F1(picrust2) == delta_f1_population` (±5e-4) and the unscoreable cells are identical across arms; (5) **plan coverage** — `{folds in fold_plan} × {densities any launched run covers}` is exactly the `(fold, density)` grid the cells table scores (`accounting.completeness.plan_coverage`), gated in `--mode fold`/`density` only (the control reference families are single-fold by construction, so there the grid is reported with `ok: null`); (6) **run inventory** — the DISCOVERED run directories equal the inventory DERIVED from the fold plan, the frozen config and this invocation's declared scope EXACTLY: no missing key, no unexpected directory, no two directories claiming the same `(fold, density, domain, routing)` / `(fold, density, forced_answer)` key (a duplicate is a hard failure naming BOTH paths in every mode, never a silent merge), and every run's own declared scope matching the plan (prFBA `meta.select_modes` present and identical on every fold; a forced-routing PICRUSt2 run naming its own domain; `meta` not contradicting the directory name) — see `accounting.completeness.inventory`, whose `derivation` block records which axis came from where; (7) **one denominator** — every `arms.<arm>.<space>.population.f1.n` equal, and no PICRUSt2 `(fold, domain)` left unlaunched (`meta.notes.picrust2_domains_not_run == {}`) unless `--expect-domains` declared the narrower scope | A truncated / resumed-stale / plan-mismatched run is otherwise read as a method with lower **coverage**: 60 deleted md5s moved the primary delta −0.2610 → −0.2347. Leg 1's mode-set clause closes a blind spot that legs 1–4 could not see: deleting the `graduated_forced` block from one forced run left C14 passing at exit 0 while `prfba_nestedcv_forced` — a headline number and D7's risk-coverage/AURC source — silently fell from n=1,922 to n=1,258 against PICRUSt2's 1,922. Leg 5 closes the coarser version: a fold the plan lists but that is absent from **both** methods' runs shrinks both denominators equally, so legs 1–4 and the D22 tripwire stay silent and the headline is quoted over 2 of 3 folds at exit 0 (at 10 folds, one OOM-killed fold = 90% of the planned population). Legs 6–7 (below) invert the question that legs 1–5 all ask of whatever runs happen to be ON DISK; a domain whose run was never *launched* is still treated differently PER CELL (`unscoreable_method_not_run` for every arm, equal denominators, D22) but in the gated modes it now FAILS leg 7 instead of being a note, because equal denominators over a 10.4%-smaller population behind a headline that never mentions Archaea is the very laundering this control exists to stop (`--expect-domains bac` declares a deliberately narrower scope). The exit status is decided on THIS invocation's own C14 (not on the failure-preserving merge), and a re-evaluation of the SAME mode may clear that mode's own recorded fail — otherwise the stale entry would keep the exit at 2 after the operator fixed the run; the cross-mode protection is unchanged (a fold-mode fail still survives an include-self pass, with `passed_in_other_modes` beside it). |
| **C15** | **Trait-column screen disclosure** (D20). | reported, not gated | Per reference: `n_dropped` by reason, `retry_used`, and which ECs were reinstated. Names any fold that used §5.4's aggressive retry, because its `singleton`/`crash` columns are not reinstatable and leave that fold's PICRUSt2 numbers carrying an undisclosed recall penalty. |

### 7.1 Exit status

`score_folds.py` writes every artefact first (the operator needs them to debug),
then returns:

| code | meaning |
|---|---|
| **0** | every HARD control of THIS invocation passed |
| **2** | **C14 failed** — the scored population is not the planned one: the numbers in `summary_*.json` are void, the headline sentence is replaced by the `RUN_INCOMPLETE` / `DENOMINATORS_DIFFER` stamp, and `headline.valid`, `meta.run_valid` and every `arms.<arm>.<space>.run_valid` are `false` |
| **3** | C14 passed but another **HARD** control failed (C1–C6, C8, C9, C11, C12), **or** the artefact is stamped invalid for a reason no control owns — an indeterminate primary claim (`primary_claim_indeterminate_*` / `primary_claim_bootstrap_below_200`, §7.2) or a control this mode always evaluates that was not emitted at all (`controls_not_emitted:<ids>`). The numbers are over the planned population, but a setup/leakage/calibration gate did not hold or the claim cannot be stated. A control is HARD unless its own `gate` text begins with `reported` (C7, C10, C13, C15 are reported-only and can never move the exit status) |
| **1** | the invocation was REFUSED before any artefact was written: `--bootstrap` below `MIN_BOOTSTRAP` (200) in a gated mode without `--allow-underpowered-bootstrap` (§7.2), or a `--cells` table that can neither be shown to match the runs on disk nor be recomputed where it sits |

Before this rule the exit status was decided by C14 alone, so a C5 fail (proven by
injection: arm A ratio 0.537) exited **0** and an orchestrator that trusted the
return code never saw it. Reading `controls.json` is still the complete record —
the exit status is now a faithful summary of it, not a substitute.

**`headline.valid` reflects EVERY hard control, not just C14.** It, `meta.run_valid`,
`meta.run_invalid_reasons` and every `arms.<arm>.<space>.run_valid` are written by ONE
function (`stamp_invalid`) from ONE reason list, AFTER every control has been
evaluated, and the headline sentence itself is prefixed `INVALID (<reasons>) -- ...`
with the provisional sentence kept behind it. Before that, validity reflected only
C14 + the inventory, so a C4 fail on a truncated `refs/*/manifest.json` exited 3 while
publishing a fully formed `"valid": true` artefact whose skill score had flipped sign
(−2.0449 → +1.5523), and a failing phenotype identity control did the same. Any
consumer reading `summary_*.json` rather than the exit code was unprotected; both
paths now agree by construction.

A control that this mode ALWAYS evaluates and that is nevertheless absent from
`controls.json` is treated as a failure, not as a pass (`controls_not_emitted`): an
absent entry read exactly like a green one, so a silently unemitted C14 was invisible.
The must-emit sets are per mode (`EXPECTED_CONTROLS`: fold/density C3, C4, C6, C7,
C10, C14; include_self adds C2; plumbing C3, C9; recipe_control C1). Controls whose
emission depends on the DATA rather than the mode (C5 needs a pruned stock reference,
C8 a prFBA run, C9 a plumbing run, C11/C12 the phenotype arm, C13/C15 a manifest
field) are listed in `CONDITIONAL_CONTROLS` and REPORTED when absent
(`accounting.controls_provenance.conditional_not_emitted`) — a gate that cried wolf on
a legitimately absent control would be worse than the hole it closes.

**Every control verdict is attributable to a run.** `controls.json` is written by a
read-modify-write merge across four modes, and its entries carried no timestamp, no
cells sha and no inventory fingerprint, so five pre-seeded `"pass"` entries made a
fresh `--mode fold` run read **15/15 pass** with four of them never evaluated by that
invocation at all. Each FRESH entry now carries
`provenance = {mode, evaluated_at, cells, cells_sha256, inventory_fp, answer_key_fp,
scorer_code_fp, git_commit}`, and `accounting.controls_provenance` classifies every
entry in the merged file **by mode first**: `fresh` (this invocation),
`from_other_modes` (C1 from recipe_control, C2 from include_self, C9 from plumbing —
legitimate; the mode, the time and `answer_key_matches` / `scorer_code_matches` are
recorded and printed), `without_provenance` (pre-seeded or pre-contract) and
`stale_same_mode` (this same mode recorded it and did NOT re-evaluate it now).
Classifying by fingerprint instead would label every legitimate cross-mode control
"different inputs/code", because `inventory_fp` contains the mode and the
mode-selected run set and so can never match across modes — a cry-wolf in every
campaign log. The stdout control list is annotated the same way
(`control C1 pass   [from --mode recipe_control @ <ts>]`, with
`-- ANSWER KEY OR CODE DIFFERS` appended when that is what the fingerprints say).
The failure-preserving cross-mode merge is unchanged; what changed is that a stale
pass is no longer indistinguishable from a fresh one.

An explicitly named `--phenotype-summary` that does not exist is a startup error
(EXIT 1), not an arm that did not run: C11/C12 would otherwise simply be absent from
`controls.json`, which reads exactly like a clean run.

### 7.2 The primary claim must be decidable, and the CI must be supportable

`exceeds_sesoi` was a boolean over three outcomes, so "no interval was computed" and
"the effect is below the smallest effect size of interest" were published in the same
word: `--bootstrap 0` produced `ci95: [None, None]`, `exceeds_sesoi: false` and a
headline reading **"SESOI 0.02 not exceeded"** at EXIT 0 with `valid: true` — the
campaign's primary conclusion inverted by one flag. Every contrast now carries
`sesoi_verdict ∈ {exceeded, not_exceeded, indeterminate_no_ci}` beside
`ci_available`, `ci_unavailable_reason` (`bootstrap_disabled` /
`fewer_than_2_holdout_groups` / `no_paired_cells`, returned by `cluster_bootstrap`
itself), `ci_underpowered` and `ci_underpowered_reason`; the sentence says
`SESOI 0.02 INDETERMINATE -- no CI was computed (<reason>), so this run cannot say
whether the SESOI is exceeded or not`; and in a gated mode an indeterminate or
underpowered primary claim is an invalidity reason (EXIT 3).

`MIN_BOOTSTRAP = 200` (the smoke recipe's own value). In `--mode fold` / `--mode
density`, `--bootstrap` below it is **REFUSED at EXIT 1 before any artefact is
written** unless `--allow-underpowered-bootstrap` is passed, in which case the run
proceeds and the claim and the run are stamped. `--bootstrap 5` otherwise produced a
perfectly quotable CI [−0.2895, −0.2224] from five resamples with no warning
anywhere.

---

## 8. Component 5 — `scripts/score_folds.py`

```
score_folds.py
  --mode fold|recipe_control|include_self|density
  --fold-plan .../data/fold_plan.json
  --panel     .../data/panel_eligible.json
  --picrust2  .../data/picrust2            # root; globbed for f*_*_d*/predictions.json
  --prfba     .../data/prfba               # root; globbed for f*_d*/predictions.json
  --truth     region_validation/data/benchmark_truth.json
  --functions region_validation/picrust2_comparison/data/genome_function_sets.json
  --maps      region_validation/picrust2_comparison/data/function_maps.json
  --neighbors .../data/neighbors.json
  --refs      .../data/refs
  --cells     .../data/cells.csv.gz
  --summary   .../data/summary_<rank>.json
  --controls  .../data/controls.json
  [--bootstrap 2000] [--rank-blocked genus] [--space ec,rxn]
  [--ec-universe ...] [--fanout ...] [--annotation-ceiling ...] [--templates ...]
  [--phenotype-summary .../phenotype_summary.json]   # default: beside --cells
  [--require-phenotype-binding]       # SPEC 11.2 passes this (C11/C12, section 7.1)
  [--allow-underpowered-bootstrap]    # --bootstrap < 200 in a gated mode (7.2)
  [--expect-folds ...] [--expect-densities ...] [--expect-domains ...]
  [--expect-nsti none|primary|all] [--expect-forced-answer none|primary|all]
  [--force]
```
Every one of `--truth`, `--functions`, `--maps`, `--ec-universe`, `--fanout`,
`--annotation-ceiling` and `--templates` is an ANSWER-KEY or vocabulary input: it is
fingerprinted, named in the artefact and re-verified after the cells pass (§8.6 /
the `<cells>.fp` paragraph above).

### 8.1 Cell table (`data/cells.csv.gz`, one row per amplicon × arm × space)
Columns, in order (superset of `compare.py:CELL_FIELDS`, so the existing
renderer keeps working):
```
amplicon, md5, region, domain, src_genome_id, fold, group_id, density,
arm, space, outcome, n_truth, n_pred, tp, recall, precision, f1, jaccard, mcc,
n_reps, tier, abstained, nsti, nsti_bin, picrust2_domain, routing,
rank1_identity, rank1_rank, ident_bin, novelty_bin,
n_truth_accessory, f1_accessory, n_truth_impossible,
template_visible_frac, min_align_dropped, unscoreable_reason,
abstain_reason, n_truth_impossible_own_domain
```
The last two are appended after the original contract (never reordered), so the
`compare.py` renderer keeps working. `score_folds.py` REFUSES a `--cells` table
that lacks them rather than treating their absence as blank: a table that cannot
be audited for completeness is exactly what C14 exists to reject.
- `arm` ∈ `prfba_<select_mode>`, `prfba_<select_mode>_forced`, `picrust2`,
  `picrust2_nsti_routed`, `base_constant`, `base_nn1`, `base_family_donor`,
  `base_clade_max`, `base_random_abstain`.
- `space` ∈ `ec`, `rxn` (reactions via `function_maps.ec_rxn`, applied identically
  to truth and to every arm). **Two different EC→reaction routes are used
  deliberately and must not be conflated:** `ec_rxn` here, for set-overlap scoring
  in the same namespace as the predecessor study; the curated EC→template-role
  route (D10) in the phenotype arm, for model *building*. Neither is used for the
  other's purpose.
- `outcome` ∈ `scored`, `abstained`, `pred_empty`, `unscoreable_empty_truth`,
  `unscoreable_min_align`, `unscoreable_method_not_run` (D22 — the cell's required
  method run was never launched; applied to EVERY arm so all denominators stay
  equal).  The three `unscoreable_*` outcomes are properties of the CELL, never of
  an arm, and `meta.unscoreable_uniform_across_arms` asserts exactly that.
- `abstain_reason` is the runner's own reason verbatim (`abstain_below_family`,
  `abstain_reliability`, `no_surviving_neighbour`, …) or one of the sentinels
  `no_record` / `no_prediction` / `no_forced_record`, which mean the md5 is absent
  from a launched run's table — a run defect that C14 fails on (D21), not an
  abstention.
- `n_truth_impossible_own_domain` is C10's second ceiling: truth ECs carried by no
  surviving genome of the cell's OWN domain reference, i.e. what a forced-routed
  PICRUSt2 query can reach (D23).  `n_truth_impossible` keeps the union semantics.
- Under `--routing forced` the two domain references' predictions are merged on
  `(md5, truth domain)`, not on md5: an amplicon md5 carried by both a bacterial
  and an archaeal genome is a query of both references, and a min-NSTI merge would
  hand such a cell the other domain's prediction — silently making the primary
  forced arm the secondary nsti-routed arm there (D9).  The min-NSTI merge is used
  for `--routing nsti` only, which is what that arm is defined to be.
- `n_truth` is computed from the held-out genome's own `ec` set restricted to
  `panel_ec_universe`. Cells with `n_truth == 0` are `unscoreable_empty_truth`
  and excluded from every mean (reusing `compare.metrics`' `None` semantics and
  `Agg`'s denominators) — identically for every arm and baseline.
- `n_truth_impossible` = C10 count for that cell.
- `f1_accessory` restricts truth and prediction to ECs with training-fold
  prevalence in [0.05, 0.95].
- `mcc` over the fixed `panel_ec_universe` (A1: MCC separates the baselines more
  than F1 does — constant 0.569 vs genus-donor 0.772 against F1 0.637 vs 0.803).

### 8.2 Neutral novelty covariate
`rank1_identity` = identity of the best surviving reference for that amplicon,
from the prFBA fold run's rank-1 hit against the restricted DB (the DB is
identical for both methods, so this is a property of the data, not of a method).
`rank1_rank` = the highest rank at which `rank1_genome` and the test genome agree.
**Both methods are binned by these two covariates only** — never by
`best_nonself_identity` or NSTI, which are each method's own confidence.
`novelty_bin` uses `score_genecap.ident_bin`'s edges
(`<0.865 / 0.865-0.945 / 0.945-0.987 / 0.987-0.995 / >=0.995`); under the primary
design the top two bins must be empty (C6b).

### 8.3 Baseline arms (A1's required controls, D13)
- `base_constant` — the top-K most prevalent ECs among that fold's **training**
  (surviving reference) genomes, `K = median(n_truth)` over training genomes.
  Measured floor on this panel: mean EC F1 0.637, MCC 0.569.
- `base_nn1` — the single surviving reference with the highest amplicon-level
  identity (`rank1_genome`), its EC set copied verbatim. **Both methods must beat
  this**; prFBA's median `n_reps = 1` means prFBA *is* this baseline plus an
  abstention policy, so the gap is prFBA's algorithmic contribution.
- `base_family_donor` — a random surviving same-family genome (seeded); under
  leave-family-out, same-order (D13).
- `base_clade_max` — the max-F1 surviving genome in the nearest surviving clade:
  the copy-a-relative ceiling.
- `base_random_abstain` — the **forced-answer** arm's EC predictions masked down
  to the shipped arm's coverage by a seeded random choice of which cells to
  abstain on (i.e. same predictions, same coverage, uninformative abstention). Measured on shipped data: 0.6006 vs 0.7985
  own-choice at coverage 0.751 — abstention is informative and is credited as a
  measured property, not as a conditional-mean boost.

### 8.4 Endpoints (D7, A1 amendment D)
- **Primary** — population EC F1 difference (`prfba_<nested-CV winner>` −
  `picrust2`) under leave-genus-out, all regions pooled, cluster-bootstrap CI over
  **holdout groups** (`--bootstrap 2000`, resample groups with replacement).
- **Co-primary** — risk-coverage curve per method, each sweeping its OWN knob from
  a single full-coverage run: prFBA on `rank1_identity` (forced-answer arm),
  PICRUSt2 on `nsti`. Summarised by **AURC** and **F1@coverage ∈ {0.5, 0.75, 0.9, 1.0}**.
  The **companion curve on the other method's knob** (prFBA on `nsti`, PICRUSt2 on
  `rank1_identity`) is reported beside it as `aurc_cross_knob` /
  `f1_at_coverage_cross_knob`, so an AURC comparison can separate "this method's
  abstention policy is better" from "this knob is a better confidence signal". The
  own-knob curves remain the D7 endpoint.
- **Skill score** beside every raw F1:
  `(F1_arm − F1_base_constant) / (F1_base_clade_max − F1_base_constant)`.
- **SESOI** `0.02` raw ΔF1 (≈12% of the 0.637–0.803 interpretable band),
  pre-registered in `frozen`.
- **Secondary** (BH-FDR across the whole secondary family): per region, per
  novelty bin, per domain (Bacteria and Archaea reported separately — never a
  pooled mean that hides a 13% subpopulation), reaction space, MCC,
  accessory-only, density 25/50/100% for both methods with both slopes.
  **Holm-over-8-regions is dropped** — it treats correlated tests on the same
  organisms as independent families.
- **Inference** — amplicon-level with clustering by holdout group: paired
  differences with a group-level cluster bootstrap (primary) and a GLMM
  `f1 ~ arm*region + rank1_identity + (1|src_genome_id) + (1|group_id)`
  (secondary, via `statsmodels` MixedLM on the paired difference if the GLMM does
  not converge). **Report ICC and n_eff next to every p-value** (measured ICC 0.66
  genome / 0.48 genus / 0.35 family; design effects 5.0 / 10.3 / 17.5;
  67,038 amplicons → n_eff ≈ 13,400 / 6,500 / 3,800). **Never** test the 10 fold
  means as replicates — overlapping training sets, no unbiased k-fold variance
  estimator. The write-up states plainly that the design is overpowered for the
  headline (SE ≈ 0.002–0.005) and that the binding risks are bias and selection,
  not variance; power statements are reserved for thin strata
  (n_eff ≈ 100 → MDE ≈ 0.04 conditional, ≈0.11 population).
- **Nested CV bookkeeping** (D8) — for each fold, the select-mode that won on the
  training folds, and how often the inner winner differs.

### 8.5 `data/summary_<rank>.json`
```json
{"meta": {"frozen": {...}, "rank_blocked": "genus", "n_folds": 0,
          "n_eligible": 0, "n_reference_genomes_median": 0,
          "n_cells_scored": 0, "n_unscoreable_empty_truth": 0,
          "n_unscoreable_min_align": 0, "generated": "<iso8601>",
          "git_commit": "<sha>"},
 "accounting": {"per_fold": [{"fold": 0, "n_cells": 0, "n_prfba_abstained": 0,
     "n_picrust2_emitted": 0, "n_truth_empty": 0, "n_matched": 0,
     "n_impossible_ec_mean": 0.0, "structural_recall_ceiling": 0.0,
     "structural_recall_ceiling_picrust2_own_tree": 0.0,
     "structural_recall_ceiling_by_domain": {"Archaea": {"...": "union + own"}},
     "n_min_align_dropped": 0, "n_method_not_run": 0, "n_no_record": {}}],
   "completeness": {"runs": {"<run>": {"self_consistent": true,
                                       "plan_complete": true, "n_missing": 0,
                                       "select_modes": [], "tables_present": [],
                                       "tables_missing": [],
                                       "tables_unexpected": []}},
                    "plan_coverage": {"gated": true, "folds_in_plan": [],
                                      "folds_scored": [],
                                      "folds_in_plan_not_scored": [],
                                      "densities_in_launched_runs": [],
                                      "densities_scored": [],
                                      "grid_missing": [], "n_grid_missing": 0,
                                      "ok": true},
                    "inventory": {"gated": true, "ok": true, "reasons": [],
                        "derivation": {"spec": "SPEC 11.2: 60 + 10 + 30 + 10",
                            "folds_source": "fold_plan.folds",
                            "folds_in_plan": [], "folds_expected": [],
                            "densities_source": "frozen.density_fractions ...",
                            "densities_expected": [],
                            "densities_in_frozen_plan": [],
                            "densities_with_a_reference_build": [],
                            "densities_in_frozen_plan_not_built": [],
                            "domains_source": "domains of the plan's own test cells",
                            "domains_expected_by_fold": {"0": ["arc", "bac"]},
                            "nsti_runs_expected": "none|primary|all",
                            "forced_answer_runs_expected": "none|primary|all",
                            "primary_density": "1.0"},
                        "n_expected": 0, "n_discovered": 0,
                        "n_run_directories": 0,
                        "counts_discovered": {"picrust2": 0, "prfba": 0},
                        "expected": ["picrust2|fold=0|d=1.0|domain=bac|forced"],
                        "discovered": {"<key>": ["<dir>"]},
                        "missing": [], "duplicates": [{"key": "", "dirs": []}],
                        "unexpected": [], "scope_mismatch": [
                            {"run": "", "dir": "", "defect": "", "declared": {},
                             "planned": {}}],
                        "select_modes_by_run": {}, 
                        "select_modes_consistent_across_runs": true,
                        "select_modes_in_frozen_plan": [],
                        "fingerprints": {"predictions": {"picrust2/<run>":
                                             {"size": 0, "sha256": "<sha>"}},
                                         "manifests": {}, "inputs": {},
                                         "answer_key": {"functions": {}, "maps": {},
                                             "ec_universe": {}, "fanout": {},
                                             "annotation_ceiling": {},
                                             "templates/<role>": {}}},
                        "inventory_fp": "<sha256 over the whole inventory>",
                        "input_paths": {"functions": "<abs path>", "maps": "...",
                                        "ec_universe": "...", "fanout": "...",
                                        "annotation_ceiling": "...",
                                        "templates": "..."},
                        "answer_key_fp": "<sha256 over the answer-key bytes>",
                        "scorer_code_fp": "<sha256 over the scorer's own code>",
                        "scorer_code_files": {"<module>.py": {"path": "", "size": 0,
                                                              "sha256": ""}}},
                    "arm_denominators": {"gated": true, "n_by_arm": {},
                                         "distinct_n_by_space": {"ec": [0]},
                                         "ok": true},
                    "n_runs_incomplete": 0, "n_missing_record_rows": 0,
                    "missing_record_rows_by_arm": {},
                    "unscoreable_by_space": {},
                    "unscoreable_uniform_across_arms": true},
   "trait_column_screen": {"<reference>": {"n_dropped": 0, "by_reason": {},
                                           "retry_used": false}},
   "full_prevalence_ec_reinstated_by_run": {"<run>": {"<reference>": 0}}},
 "arms": {"<arm>": {"<space>": {
     "population": {"recall": {"value":0.0,"n":0,"sd":0.0}, "precision": {...},
                     "f1": {...}, "jaccard": {...}, "mcc": {...}},
     "conditional_matched": {"...": "same keys", "n_matched": 0},
     "micro": {"recall": 0.0, "precision": 0.0, "f1": 0.0},
     "coverage": 0.0,
     "skill_score": 0.0,
     "aurc": 0.0,
     "f1_at_coverage": {"0.5": 0.0, "0.75": 0.0, "0.9": 0.0, "1.0": 0.0},
     "by_region": {"V4": {"...": "population block"}},
     "by_novelty_bin": {"0.945-0.987": {"...": "population block"}},
     "by_domain": {"Bacteria": {"...": "population block"}},
     "by_density": {"0.25": {"f1": 0.0}, "0.5": {"f1": 0.0}, "1.0": {"f1": 0.0}},
     "density_slope_f1_per_log10_tips": 0.0,
     "template_visible_frac_mean": 0.0,
     "run_valid": true, "run_invalid_reasons": []}}},
 "contrasts": {"<armA>__vs__<armB>": {"<space>": {
     "delta_f1_population": 0.0,
     "ci95": [0.0, 0.0], "bootstrap_n_groups": 0,
     "icc_genome": 0.0, "icc_group": 0.0, "n_eff": 0,
     "p_cluster_bootstrap": 0.0, "p_glmm": 0.0, "q_bh": 0.0,
     "exceeds_sesoi": false,
     "sesoi_verdict": "exceeded|not_exceeded|indeterminate_no_ci",
     "ci_available": true, "ci_unavailable_reason": null,
     "ci_underpowered": false, "ci_underpowered_reason": null,
     "min_bootstrap": 200,
     "delta_f1_at_matched_coverage": {"0.5": 0.0, "0.75": 0.0, "0.9": 0.0},
     "by_region": {"V4": {"delta_f1_population": 0.0, "q_bh": 0.0}},
     "by_domain": {"Archaea": {"delta_f1_population": 0.0, "q_bh": 0.0}}}}},
 "nested_cv": {"fold_0": {"inner_winner": "exact_tie",
                          "inner_f1": {"legacy": 0.0, "exact_tie": 0.0,
                                       "graduated": 0.0, "exact_tie_guards": 0.0}},
               "winner_instability": 0.0},
 "annotation_ceiling": {"recall": 0.904, "precision": 0.743,
   "source": "picrust2_comparison/data/annotation_ceiling.json"},
 "headline": {"sentence": "<rank blocked>, <n reference genomes> species-dereplicated reference genomes, coverage <c>: ...",
              "rank_blocked": "genus", "reference_density_genomes": 0,
              "coverage": 0.0, "denominators_consistent": true,
              "prediction_runs_complete": true,
              "plan_coverage_complete": true,
              "run_inventory_complete": true, "arm_denominators_equal": true,
              "run_inventory_fp": "<sha256>",
              "answer_key_fp": "<sha256>", "scorer_code_fp": "<sha256>",
              "sesoi_verdict": "exceeded", "primary_claim_decidable": true,
              "hard_controls_failed": [], "controls_evaluated_this_run": [],
              "run_invalid_reasons": [], "valid": true}}
```
`accounting` additionally carries two provenance blocks:
```
{"accounting": {
   "phenotype_binding": {"path": "", "discovered": "auto-discovered beside --cells",
       "present": true, "readable": true, "bound": true, "require_binding": false,
       "binding_contract": 1, "unbound_reason": null,
       "defects": [{"defect": "", "detail": {}}],
       "fingerprints": {"phenotype_summary": {}, "per_genome_jsonl": {},
                        "identity_control_json": {}},
       "declared": {"n_genomes": 0, "n_models": 0, "n_gated_out": {},
                    "identity_control": {}, "run_dirs": [], "binding": {}}},
   "controls_provenance": {
       "this_run": {"mode": "", "evaluated_at": "", "cells": "",
                    "cells_sha256": "", "inventory_fp": "", "answer_key_fp": "",
                    "scorer_code_fp": "", "git_commit": ""},
       "fresh": ["C14"], "expected_this_mode": ["C14"], "missing_expected": [],
       "from_other_modes": [{"id": "C1", "mode": "recipe_control",
                             "evaluated_at": "", "cells_sha256": "",
                             "inventory_fp": "", "answer_key_matches": true,
                             "scorer_code_matches": true}],
       "without_provenance": ["C2"],
       "stale_same_mode": [{"id": "", "mode": "", "inventory_fp": ""}],
       "conditional_not_emitted": {"C5": "<the condition that emits it>"},
       "n_absent_entirely": []}}}
```
`meta` additionally carries `n_unscoreable_method_not_run`,
`n_unscoreable_by_space`, `unscoreable_uniform_across_arms`,
`notes.denominators_consistent` / `notes.denominator_check` and the provenance
triple `notes.run_inventory_fp` / `notes.run_inventory_gated` /
`notes.n_run_directories`; `run_valid` / `run_invalid_reasons` carry
`duplicate_run_directories`, `run_inventory_incomplete`,
`unexpected_run_directory`, `run_scope_mismatch`, `arm_denominators_unequal` and
`picrust2_domain_never_launched` beside the older reasons.  `<cells>.fp` is written
next to every cells table (`{cells: {size, sha256}, inventory_fp, answer_key_fp,
answer_key_paths, scorer_code_fp, scorer_code_files, mode, generated, git_commit,
fingerprints}`); without `--force` the scorer reuses a cells table ONLY when that
sidecar matches the cells bytes, the inventory it just inspected (which now includes
the ANSWER KEY's bytes) **and the scorer's own `scorer_code_fp`**, otherwise it
recomputes (or refuses, if the directory is read-only).  A sidecar with no
`scorer_code_fp` (written before this contract) is never reused.  Two measured
reasons: a tampered `genome_function_sets.json` (2,822 truth ECs dropped from 20
held-out genomes) re-published the shipped headline VERBATIM through the reuse path
while the same key with `--force` moved picrust2 0.6955 → 0.6872; and a scorer edited
to halve picrust2's F1 re-published 0.6955 / delta −0.2610 from a table its own
`build_row` now scores at 0.3478 / **+0.0867** — the primary conclusion's SIGN.
`git_commit` cannot catch either (it is `rev-parse HEAD` on a tree that may be dirty).

**The answer key is inside provenance.**  `--functions` (the truth EC sets behind
every F1), `--maps`, `--ec-universe` (resolved beside `--panel`, which no reader could
otherwise reconstruct), `--fanout`, `--annotation-ceiling` and the `--templates`
files `load_template_visible` actually reads (`TEMPLATE_FILES`, never a glob) are
fingerprinted into `inventory.fingerprints.answer_key`, hashed into
`inventory.answer_key_fp` (and so into `inventory_fp`), NAMED in
`meta.notes.answer_key_paths` / `inventory.input_paths`, printed at startup, and
RE-READ and re-compared after the cells pass exactly like the run artefacts — a
mid-run rewrite lands in `inventory.inputs_changed_during_scoring` and voids the run
(EXIT 2).  Fingerprints are over CONTENT (size + sha256) and never over paths, so a
tree copied elsewhere to review it still reuses its own cells table.  The scorer's own
code identity is a FIXED file list (`SCORER_CODE_FILES`: `score_folds.py`, the three
promoted `design_audit_*` modules, `score_genecap.py`, `compare.py`,
`annotation_ceiling.py`, `_common.py`, `phenotype_arm.py`), never `sys.modules`
introspection — `load_template_visible` imports `phenotype_arm` lazily, so an
introspected list would differ between the pre-scoring reuse check and the
post-scoring sidecar and nothing would ever be reusable.  When
`headline.valid` is false the sentence is REPLACED by a `RUN_INCOMPLETE ...` or
`DENOMINATORS_DIFFER ...` stamp that carries the provisional text inside it, so a
quoted headline can never outlive a void run (D21, D22).  The arms and contrasts
tables are still written in full on a void run, so the same verdict is mirrored
where the numbers actually live: `meta.run_valid` /
`meta.run_invalid_reasons` (`denominators_differ`,
`cells_scored_without_a_record`, `launched_run_incomplete`,
`planned_fold_density_never_scored`) and, per arm and per space,
`arms.<arm>.<space>.run_valid` beside `population.f1`.  A consumer that reads an
arm number without checking `headline.valid` therefore still sees the flag.
`headline.sentence` is generated, not hand-written, and MUST name the three scope
parameters (rank blocked, reference density, coverage) — A1 amendment I.

---

## 9. Component 6 — `scripts/phenotype_fold.py`

A thin harness over the already-written
`picrust2_comparison/scripts/phenotype_arm.py` (399+ lines, imports clean in
py_venv). **Do not reimplement any of it**; import
`preload_biochemistry, load_template, template_for, genome_from_ec, build_draft,
build_gapfilled, complete_media, load_atp_media_panel, core_atp_vector,
growth_vector, essentiality, rxn_sets, jaccard, hamming, sanity_gate, yield_check,
precursor_producibility`.

```
phenotype_fold.py
  --cells    .../data/cells.csv.gz
  --fold-plan .../data/fold_plan.json
  --panel    .../data/panel_eligible.json
  --prfba    .../data/prfba
  --picrust2 .../data/picrust2
  --functions region_validation/picrust2_comparison/data/genome_function_sets.json
  --outdir   .../data/phenotype
  --n-genomes 200
  --select-mode exact_tie
  [--arms truth,prfba,picrust2,base_nn1,base_family_donor,base_constant]
  [--pipeline as_shipped,controlled]
  [--readouts atp,growth,essentiality,structure]
  [--workers 24] [--seed 1729] [--min-truth-ec 200]
```

### 9.1 Subset selection
Eligible phenotype genomes: in the fold's test set, `n_truth_ec >= --min-truth-ec`
(≥200 leaves ≈8.4k), prFBA non-abstained AND PICRUSt2 emitted (so all arms exist).
Sample `--n-genomes` **stratified by domain × template and by per-genome
Δ(prFBA, PICRUSt2) EC F1 quintile** (A1 amendment G: so the arm can actually test
whether EC disagreement propagates). Genomes with BacDive/Biolog data are
over-sampled if such a list is supplied (`--experimental-list PATH`, optional).
Abstentions and no-prediction md5s produce **no model** and are reported as
unscoreable with their own denominator — never scored as a zero phenotype.

### 9.2 Build procedure (identical for every arm; only the input EC set varies)
1. `preload_biochemistry("/home/freiburger/Documents/ModelSEEDDatabase")` once per
   worker (13 s) — **required**, else `MSATPCorrection` raises `FileNotFoundError`
   from modelseedpy's CWD-relative `config.cfg` biochem path.
2. Template is a pure function of the **held-out genome's true lineage**, read from
   `panel_eligible.json:genomes.<id>.template`:
   `Archaea → ArchaeaTemplateV6`;
   `Bacillota*|Firmicutes|Actinomycetota|Actinobacteria|Chloroflexota|Deinococcota|Thermotogota|Mycoplasmatota|Tenericutes → GramPosModelTemplateV7`;
   else `GramNegModelTemplateV7`. **Never** `auto_select_template` — it would make
   the template depend on predicted content (and downloads a classifier pickle).
3. EC → model input (D10): a template role is granted iff **every** complete EC
   number in its name is present in the input EC set; each granted role becomes one
   synthetic `MSFeature` whose id is the template role id and whose
   `ontology_terms["RAST"]` is the role name verbatim; then
   `MSBuilder(genome, template, ontology_term="RAST").build(id, "0",
   allow_all_non_grp_reactions=False, annotate_with_rast=False)`.
   **`annotate_with_rast=False` is mandatory** — the default `True` re-annotates
   the input and erases the very difference the arm measures.
4. **No gapfilling for the primary readout.** For the secondary readouts, gapfill
   **every arm identically on complete medium only** (347 compounds from
   `ModelSEEDDatabase/Media/KBaseMedia.cpd` at 100 each). Never gapfill on a
   medium that is tested, and never on a minimal medium (glucose-minimal-gapfilled
   models fail the empty-medium ATP gate and do not grow on their own target
   medium). `gapfill_model` must be explicitly `False` wherever the primary
   readout is computed.
5. Media sign convention: `ModelSEEDDatabase/Media/media_list_with_meta.txt`
   columns are `compound_id;minFlux;maxFlux;concentration` in KBase convention
   (positive = uptake); cobra exchange bounds are the **negation**
   (`lower = -maxFlux`, `upper = -minFlux`). Every other exchange defaults to
   `(0, 100)` — no uptake, free excretion.
6. Two pipeline variants (A4 fatal flaw 3), both reported:
   `as_shipped` — each model gets its own `MSATPCorrection` (what a user gets);
   `controlled` — the **truth** model's `selected_media` and `build_tests()`
   output are imposed on every arm, so the only difference is gene content. The
   gap between them is the pipeline-reaction effect and must not be attributed to
   the methods. This requires inlining `MSBuilder.build_and_gapfill`'s body
   (`msbuilder.py:977-1050`) rather than calling it — which is also how the
   per-media ATP vector it currently computes and discards (`media_eval`) is
   captured.

### 9.3 Readouts
- **PRIMARY — core ATP-yield vector** (`core_atp_vector(..., no_gapfilling=True)`):
  54 media, model restricted to Core-V6 (252 reactions), 2–3 s/model. Report the
  binary Hamming at 1e-6, the signed set difference (media lost vs gained — the
  precision-sensitive direction), Spearman of the continuous yields, and ATP on the
  `empty` medium as the validity assertion (must be 0.0).
- **SECONDARY — growth + essentiality** on complete-medium-gapfilled models:
  `bio1` ratio P/T and log-ratio, binary at 1e-5; single-reaction deletion at 1% of
  wild-type. Report disagreements over **shared** reactions and
  `|essential in T and absent from P|` as two separate quantities that are never
  merged. Essentiality is evaluated **only where the truth model grows**.
- **STRUCTURAL** — reaction-set Jaccard pre- and post-gapfill, gene-associated
  reaction Jaccard, and gapfill burden `|new reactions|` as a first-class metric.
- **Error ladders for calibration**, on the same genomes, so an EC-level ΔF1
  becomes an expected phenotype distance: (a) neighbour-genome substitution
  (nearest surviving relative, then family, then order donor's EC set as a
  stand-in prediction); (b) random dropout at 5/10/20/40% and random
  false-positive addition at 10/20/30/50%; (c) **size-matched null** — a random EC
  set of the same cardinality as each method's prediction drawn from the
  template-visible universe. (c) is what prevents a phenotype "win" that merely
  reflects predicting a larger set: measured on the 8-genome smoke panel, per unit of
  perturbation false negatives move the phenotype ~2x more than random false positives
  (dropout 1.125/2.50/2.75/4.00 at 5/10/20/40% vs false-positive 1.75/1.75/1.875/2.29
  at 10/20/30/50%), while the actual predictions are false-positive dominated. Every
  rung is averaged over the **same denominator as the primary** (D18), or the
  calibration and the arms would describe different genome sets.
- **ATP denominator (D18):** a genome whose truth ATP vector is identically zero is
  ATP-unscoreable; the contract primary is the informative subset and the all-genome
  value is a labelled sensitivity. The `nn1` donor rung and the `base_nn1` arm carry
  their donor provenance (`ladders.donor.nn1.donor_source`,
  `meta.donor_nn1_source`): in a fold run the donor is the alignment rank-1 hit
  against the fold-restricted shared DB (the intended nearest surviving relative),
  while `adapt_as_shipped.py` sets `rank1_genome = sorted(reps)[0]`, i.e. prFBA's own
  first selected rep — the fixture rung is "copy prFBA's first rep", not an
  independent neighbour.
- **Prediction source (D19):** ceiling-control runs (`--mode include_self`/`plumbing`)
  are never admissible as a method arm.
- Per-arm **template-visible fraction** of the actual predictions, with the primary
  metrics also computed on the template-visible subspace (D16).

### 9.4 Aggregation
Paired per organism over the ATP-informative genomes (D18): Wilcoxon signed-rank on
per-organism (prFBA − PICRUSt2) differences, with BH-FDR taken over those primary
p-values and the all-genome value reported beside each as a sensitivity; McNemar per medium for binary ATP/growth calls with BH-FDR and
clustering by genome. Stratified by domain, template, phylum, region and novelty
bin. **Archaea reported separately and labelled low-power**: methanogenesis is
unreachable through the EC→core-template route (Methanobacteriota score only
2–3 of 54 positive ATP media), so archaeal concordance is falsely high.

### 9.5 Outputs
`data/phenotype/per_genome.jsonl` — one JSON object per (genome, arm, pipeline):
```json
{"src_genome_id": "", "fold": 0, "arm": "prfba_exact_tie", "pipeline": "controlled",
 "template": "grampos", "domain": "Bacteria",
 "n_ec_input": 0, "template_visible_frac": 0.0,
 "n_roles": 0, "n_features": 0,
 "draft": {"n_rxn": 0, "n_metabolic_rxn": 0, "n_gene_rxn": 0, "n_genes": 0},
 "atp": {"vector": {"<medium>": 0.0}, "n_positive": 0,
         "empty_medium_atp": 0.0, "gate_ok": true},
 "gapfilled": {"n_new_rxn": 0, "new_rxn": ["rxn00001", ...], "growth": 0.0,
               "sanity_gate_ok": true, "yield_ratio": 0.0,
               "n_essential": 0, "essential": ["rxn00001", ...]},
 "vs_truth": {"atp_hamming": 0, "atp_lost": ["<medium>"], "atp_gained": ["<medium>"],
              "atp_spearman": 0.0,
              "growth_ratio": 0.0, "growth_log_ratio": 0.0, "growth_agree": true,
              "ess_shared_disagree": 0, "ess_truth_only_absent": 0,
              "rxn_jaccard_draft": 0.0, "rxn_jaccard_gapfilled": 0.0,
              "gene_rxn_jaccard": 0.0, "gapfill_burden_delta": 0},
 "wall_seconds": 0.0, "status": "ok|gated|unscoreable", "status_reason": ""}
```
`data/phenotype_summary.json`
```json
{"meta": {"frozen": {...}, "n_genomes": 0, "n_models": 0, "n_gated_out": {"<arm>": 0},
          "n_unscoreable": {"<arm>": 0}, "identity_control": {"jaccard": 1.0, "pass": true},
          "atp_primary_denominator": "informative_truth_only ...",
          "n_unscoreable_atp": {"truth_atp_allzero": 0, "truth_atp_undefined": 0,
                                "n_atp_scoreable_genomes": 0, "n_truth_rows": 0,
                                "genomes": []},
          "n_atp_informative": 0, "donor_nn1_source": "<provenance>",
          "generated": "<iso8601>", "git_commit": "<sha>"},
 "arms": {"<arm>": {"<pipeline>": {"atp_hamming": {"mean":0.0,"median":0.0,"n":0},
     "rxn_jaccard_draft": {"mean":0.0,"median":0.0,"n":0},
     "growth_log_ratio": {"mean":0.0,"median":0.0,"n":0},
     "ess_shared_disagree": {"mean":0.0,"median":0.0,"n":0},
     "gapfill_burden": {"mean":0.0,"median":0.0,"n":0},
     "atp_primary_denominator": "informative_truth_only ...",
     "atp_informative_only": {"...": "alias of the primary"},
     "atp_all_genomes_sensitivity": {"atp_hamming": {}, "atp_spearman": {},
                                     "n_genomes": 0},
     "by_domain": {"Archaea": {"...": "same keys"}}}}},
 "contrasts": {"prfba_exact_tie__vs__picrust2": {"<pipeline>": {
     "delta_atp_hamming": 0.0, "wilcoxon_p": 0.0, "q_bh": 0.0, "n_pairs": 0,
     "atp_denominator": "informative_truth_only ...",
     "delta_atp_hamming_all_genomes": 0.0, "wilcoxon_p_all_genomes": 0.0,
     "n_pairs_all_genomes": 0,
     "mcnemar_per_medium": {"<medium>": {"p": 0.0, "q_bh": 0.0}},
     "ec_f1_delta_vs_phenotype_spearman": 0.0}}},
 "ladders": {"dropout": {"0.05": {"atp_hamming_mean": 0.0, "n": 0,
                                 "atp_hamming_mean_all_genomes": 0.0,
                                 "n_all_genomes": 0}},
             "false_positive": {"0.10": {"atp_hamming_mean": 0.0}},
             "donor": {"nn1": {"atp_hamming_mean": 0.0,
                               "donor_source": "<provenance, D18>"},
                       "family": {"atp_hamming_mean": 0.0}},
             "size_matched_null": {"<arm>": {"atp_hamming_mean": 0.0}}},
 "pipeline_effect": {"<arm>": {"atp_hamming_as_shipped_minus_controlled": 0.0}}}
```

### 9.6 Compute budget (measured on this box: 64 cores, 251 GB, Gurobi 13.0.1)
Worker start-up 20 s (13 s biochemistry + 7 s templates), resident 1.1 GB; draft
build 1.5 s; core ATP vector 2–3 s; full `build_and_gapfill` 38–72 s. Set Gurobi
`Threads=1` per worker. ATP-vector arm over 1,000 genomes × 6 arms ≈ 5 core-hours;
gapfilled arm over 200 genomes × 6 arms ≈ 20 core-hours (≈40 min on 32 workers).

---

## 10. Residual confounds — eliminated or measured (nothing is left silent)

| Confound | Disposition |
|---|---|
| Reference density (shared panel ≈9k vs prFBA's 112k and PICRUSt2's 27,870) | **Measured** — D5 density sweep, both methods, both slopes printed beside every headline; every claim sentence names the density in genomes (§8.5 `headline`). |
| Build-recipe dialect (our hmmalign/FastTree vs PICRUSt2's ssu-align/RAxML) | **Measured and priced** — control C1, as a three-arm decomposition on one shared query set: alignment (C−B), tree+model (B−A) and the whole recipe (C−A, the gated number). Measured at n=40: +0.0029 / −0.0016 / **+0.0013** ΔF1 (p=0.048 / 0.64 / 0.62) — both halves are inside the band and the alignment half favours our recipe; re-run at n≥500. If |C−A| exceeds 0.01 it becomes an uncertainty band on PICRUSt2's numbers. Residual, stated not measured: ssu-align's curation of the raw sequences upstream of the distributed `bac_ref.fna`, and the fact that C1 runs stock EC *counts* while the fold arm runs *binary* traits. |
| NSTI not comparable across trees (0.42x branch-length scale under our recipe) | **Eliminated as an inferential axis** — novelty strata use the neutral `rank1_identity` / `rank1_rank` (§8.2). NSTI is used only within-arm for PICRUSt2's own risk-coverage sweep, and is never compared across arms or to the as-shipped run. PICRUSt2 reported both unfiltered and at `--max_nsti 2.0`, stating that 2.0 is uncalibrated for a custom tree. |
| Reference-sequence quality never screened | **Eliminated** — §3.1 filters 3–5 (length window 1,200–1,800 nt, HMM match occupancy ≥0.90, ACGTN only) applied identically to the PICRUSt2 reference and prFBA's DB; drop counts reported. |
| `--min_align 0.80` silently dropping long queries | **Eliminated** — same length window on queries; every drop counted per region (C7) and excluded from **both** methods. |
| Prune-after-build leaking held-out sequence into PICRUSt2's model | **Eliminated** — D3, de novo rebuild per fold; no reference pruning anywhere in the fold arm. |
| Amplicon-level identity leaks surviving a full-16S rule | **Eliminated** — D1/D2 edges E2+E4, asserted by C6a/C6b, and twin cells reported separately in C2. |
| prFBA's vanished k2=20 truncation artifact | **Measured and disclosed** — D14; `k2_saturated_records` logged, and the controlled-arm coverage is explicitly declared non-comparable to the as-shipped coverage, with the as-shipped number printed beside it. `exclude_self_record` is NOT additionally applied (no double charge). |
| Exclusion boundary == prFBA's species threshold, so no species-tier calls exist | **Measured and disclosed** — per-fold tier distribution reported; the write-up states the primary arm contains no species-tier prFBA calls, so it does not measure prFBA at its species operating point. A looser arm (exclude source genome + taxon only, no identity threshold) is run at smoke scale as a species-tier sensitivity. |
| Structural prior asymmetry (prFBA predicts a real genome's EC set; HSP predicts per-EC independently, and the truth IS one real genome) | **Measured, not removed** — it is a genuine method difference. `base_nn1` gives a copy-one-genome arm on the neutral side; PICRUSt2's binarisation is `>0` and needs no sweep (values are integral). |
| Empty / tiny EC sets (161 empty, 132 with 1–49) | **Eliminated** — §3.1 filter 2 drops them from reference, truth and both methods' inputs. |
| prFBA species-dedup and the P3 gene-count tie-break are inert on a 10,000-distinct-species panel | **Measured and disclosed** — C8 counts collisions and `gene_count_fn` calls (expected 0); the write-up states these two shipped components are untested by this design. |
| MiDAS external prior | **Eliminated** — empty MiDAS CSV passed explicitly; C8 asserts `midas_nonempty_records == 0`. The write-up notes the production pipeline has a prior the benchmark does not evaluate. |
| Domain routing / archaeal thinness (1,2xx archaeal tips vs 8,xxx bacterial) | **Eliminated for the primary** (D9 forced-by-truth-domain) and **measured** in the secondary routed arm, with per-fold routing accuracy vs `domain_routing.csv`. Bacteria and Archaea are separate primary analyses. |
| In-silico PCR queries are exact substrings of reference 16S | **Measured** — error-injected arm at smoke scale: 0.5% and 1% substitutions plus indels on **test amplicons only**, reference untouched; the degradation slope for both methods is reported. Single-16S-per-genome is stated as a scope limit (it also matches PICRUSt2's own 16S-centroid design). |
| EC is coarse; a broad prevalence gradient compresses differences | **Measured** — MCC, accessory-only (prevalence 0.05–0.95) and the skill score against `base_constant` are reported beside F1, never instead of it. |
| Truth-annotation quality (0.904 recall / 0.743 precision for a genome's own annotation) | **Measured and reported** — `summary.annotation_ceiling` is printed next to every absolute F1, and empty-truth cells are excluded identically for every arm. |
| Gapfilling converges the three networks | **Measured and disclosed** — pre-gapfill draft structure and gapfill burden reported beside every post-gapfill number; post-gapfill differences are framed as a **lower bound** on functional consequence. The primary phenotype readout needs no gapfilling at all. |
| `MSATPCorrection` makes input-dependent decisions | **Eliminated for the controlled variant, measured as `pipeline_effect`** — §9.2 step 6. |
| Template choice uses true taxonomy a user would not know | **Disclosed** — identical for every arm, so it does not bias the head-to-head; a subset is rebuilt under the other bacterial template to show the template effect is smaller than the method effect. |
| Alternative optima / degenerate solutions | **Asserted** — C11 identity control per batch; binary essentiality preferred over flux magnitudes; paired per-organism differences, never pooled means. |

---

## 11. Smoke-scale configuration (this workflow) vs full scale (orchestrator)

Both are the **same CLIs with different parameters** — no second code path.

### 11.1 Smoke (target: minutes, ≤ a few hundred genomes; the gate of A1 amendment H)
```bash
PY=~/Documents/py_venv/bin/python
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
D=$SR/data_smoke

$PY $SR/scripts/plan_folds.py --phase panel  --limit-groups 40 --limit-genomes 300 \
    --outdir $D --threads 16
$PY $SR/scripts/plan_folds.py --phase prepass --outdir $D --k2 200 --gpu-cand 300
$PY $SR/scripts/plan_folds.py --phase plan    --outdir $D --k 3

for DOM in bac arc; do
  $PY $SR/scripts/build_reference.py --mode fold --fold 0 --domain $DOM --density 1.0 \
      --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json --outdir $D/refs
done
$PY $SR/scripts/run_picrust2_fold.py --fold 0 --domain bac --density 1.0 \
    --ref-dir $D/refs/f00_bac_d100 --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json \
    --outdir $D/picrust2/f00_bac_d100 --processes 16
$PY $SR/scripts/run_prfba_fold.py --fold 0 --density 1.0 --fold-plan $D/fold_plan.json \
    --panel $D/panel_eligible.json --ref-manifest $D/refs/f00_bac_d100/manifest.json \
    --ref-manifest $D/refs/f00_arc_d100/manifest.json --outdir $D/prfba/f00_d100 \
    --select-modes legacy,exact_tie --k2 200 --gpu-cand 400

$PY $SR/scripts/build_reference.py --mode recipe_control --n-stock 960 --n-holdout 40 \
    --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json --outdir $D/refs
for ARM in A B C; do
  $PY $SR/scripts/run_picrust2_fold.py --mode recipe_control --ref-dir $D/refs/recipe_$ARM \
      --queries $D/refs/recipe_control_queries.fna --truth $D/refs/recipe_control_truth.json \
      --traits ~/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files/bacteria/ec.txt.gz \
      --outdir $D/picrust2/recipe_$ARM --processes 16
done
$PY $SR/scripts/score_folds.py --mode recipe_control \
    --picrust2 $D/picrust2 --refs $D/refs \
    --truth $D/refs/recipe_control_truth.json \
    --cells $D/cells_recipe.csv.gz --summary $D/summary_recipe.json \
    --controls $D/controls.json                                   # C1

$PY $SR/scripts/build_reference.py --mode plumbing --fold 0 --domain bac --density 1.0 \
    --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json --outdir $D/refs
$PY $SR/scripts/run_picrust2_fold.py --mode plumbing --ref-dir $D/refs/plumbing_bac_d100 \
    --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json --fold 0 --domain bac \
    --outdir $D/picrust2/plumbing_bac_d100 --processes 16          # C9
# The smoke recipe above deliberately launches ONE fold, bacteria only, and no
# --forced-answer prFBA run, so it must DECLARE that scope: the inventory gate
# (C14 leg 6) otherwise reports the other folds' and the archaeal runs as missing,
# which at full scale is exactly the failure it exists to catch.
$PY $SR/scripts/score_folds.py --mode fold --fold-plan $D/fold_plan.json \
    --panel $D/panel_eligible.json --picrust2 $D/picrust2 --prfba $D/prfba \
    --neighbors $D/neighbors.json --refs $D/refs --cells $D/cells.csv.gz \
    --summary $D/summary_genus.json --controls $D/controls.json --bootstrap 200 \
    --expect-folds 0 --expect-domains bac --expect-forced-answer none --force
$PY $SR/scripts/phenotype_fold.py --n-genomes 3 --arms truth,prfba,picrust2 \
    --pipeline controlled --readouts atp,structure --cells $D/cells.csv.gz \
    --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json \
    --prfba $D/prfba --picrust2 $D/picrust2 --outdir $D/phenotype --workers 3
```
**Smoke gate — ALL must hold before the orchestrator runs full scale:**
1. `plan_folds.py --phase plan` exits 0 with every `leakage_preflight` counter 0,
   and exits 2 when the E4 edge is deliberately disabled (the assertions must be
   shown to **fire**, then pass).
2. Every fold's reference rebuild completes: `manifest.set_equality_ok == true`,
   MSA mean gap fraction < 0.10, `.model` accepted by epa-ng.
3. C3 passes (`unknown_rows_match_queries == true`) with the `g_` prefix and
   **fails** without it (run once with `--label-prefix ""` to prove the bug is real).
4. `hsp.py -m mp` runs on the binary trait table after the §5.4 screen;
   `manifest.n_dropped` recorded.
5. C1 runs end to end at n=40 and reproduces |ΔF1| < 0.01.
6. The three baselines (`base_constant`, `base_nn1`, `base_family_donor`) compute.
7. `--ref-db` demonstrably restricts prFBA: `meta.n_reference_md5` equals the fold
   survivor count and `reps_outside_reference == 0`.
8. One phenotype triplet builds with `annotate_with_rast=False`, no gapfilling,
   `atp.gate_ok == true`, and C11 identity control Jaccard 1.000.
9. The two runner self-tests exit 0:
   `_selftest_run_picrust2_fold.py` (gate #3's own prover — it asserts the C3
   PREFLIGHT message, which is what fires on a numeric label set; asserting the
   post-hoc `C3 failed` literal made this script exit 1 on a HEALTHY harness and
   would abort an orchestrator that gates the launch on it) and
   `_selftest_score_folds_inventory.py` (C14 legs 6–7: the expected inventory
   DERIVED for the full-scale shape is exactly SPEC §11.2's 60 + 10 + 30 + 10 =
   110 runs, and each campaign failure mode — killed job with its `refs/` intact,
   never-launched domain, lost density tier, leftover retry directory, short
   `--select-modes`, dead nsti run — returns `ok: false` with the key named).

### 11.2 Full scale (orchestrator)
- `--phase panel` with no `--limit-genomes`; **`--phase prepass --k2 100
  --gpu-cand 100`** (~1–2 h); `--phase plan --k 10`.
  **MEASURED CORRECTION (believe the measurement, not this bullet).** The prepass is
  not "minutes on the 5070 Ti": 97% of aligner wall is the Biopython rescore stage,
  whose cost is `n_query × gpu_cand`, so `--k2 400` is ~36M rescored pairs ≈ 9–11 h.
  The prepass only needs hits at identity ≥ 0.987 (the τ that builds the §3.2b edge
  graph), so the campaign runs `--phase prepass --k2 100 --gpu-cand 100
  --force-backend cuda` (~1–2 h).  Verified at smoke scale: `--k2 100 --gpu-cand 100`
  on the 274-genome panel reproduces the same edge graph and the same `k=3` fold
  packing as the shipped `--k2 200 --gpu-cand 300` prepass (fold 0: 12 arc / 80 bac
  test, 182 reference = 143 bac + 39 arc; every `leakage_preflight` counter 0) and the
  same 16-arm headline to 4 decimals.  The prFBA FOLD runs keep frozen `--k2 200
  --gpu-cand 400` (§1, D14) — the correction is about the prepass only.
  The whole of stages 1–6 below is driven by `scripts/run_campaign.py`
  (`--dry-run` first; it carries this correction as its stage-2 default, records it in
  `<root>/campaign_inventory.json:deviations`, and asserts the 60 + 10 + 30 + 10 = 110
  run directories + 60 reference builds against an externally-kept
  `<root>/ledger.jsonl`).
- 10 folds × 2 domains × 3 densities = **60 reference builds**. FastTree
  `-nt -gtr -gamma` on ~8.6k × ~1.53k is the long pole (~30 s at 960 tips →
  tens of minutes at 8.6k); `hmmbuild` ~2 min at 26.8k × 1.58k, `place_seqs`
  ~1 min. FastTree is single-threaded, so run the 60 builds as **60 independent
  one-process jobs, 16–24 concurrently** (`xargs -P`); budget a few hours wall.
  A serial loop would take a day.
- 60 PICRUSt2 fold runs under `--routing forced` + 10 extra `--routing nsti` runs
  (d100 only, FIX: the routed arm is never run at d025/d050) + 30 prFBA fold runs
  (prFBA runs once per fold×density, both domains in one restricted DB) with all
  four select-modes, plus 10 `--forced-answer` runs at d100 for the risk-coverage
  curve.
- `score_folds.py --mode fold --rank-blocked genus --bootstrap 2000
  --expect-densities 0.25,0.5,1.0 --expect-domains bac,arc --expect-nsti primary
  --expect-forced-answer primary --require-phenotype-binding --force`, then
  `--rank-blocked family` over the 5 family blocks (§3.4)
  and `order` exploratory.  `--bootstrap 2000` is above `MIN_BOOTSTRAP` (200), so the
  primary claim has a supportable CI (§7.2); `--require-phenotype-binding` means the
  phenotype arm's summary must carry the `meta.binding` block that ties it to THIS
  cells table, fold plan and answer key (§7-C11), so run the fold scoring AGAIN after
  the phenotype arm finishes — the second pass reuses the cells table (its `.fp`
  sidecar still matches) and only re-derives the summary and C11/C12 — **omit
  `--force` on that second pass** (with it the 110-run cells table is rescored from
  scratch for nothing; without it the scorer verifies the sidecar against the runs,
  the answer key and its own code and reuses the table, which is what the `.fp`
  contract exists for).  Measured end to end at smoke scale: pass 1 with `--force`
  writes the table, pass 2 without it prints `reusing existing cells.csv.gz:
  fingerprints match (...)` and settles C11/C12 against the now-bound phenotype
  summary.  Those five flags are the invocation's DECLARED SCOPE
  (`--expect-domains bac,arc` makes the declared set equal the 60 runs the
  orchestrator launches unconditionally; the default — the domains of each fold's
  own test cells — would not expect an arc run for a fold that happened to draw no
  archaeal test genus, and that fold's launched arc run would then read as an extra);
  C14 leg 6 derives the expected inventory from the fold plan + `frozen` +
  that scope and asserts the discovered run directories equal it exactly — all
  **60 + 10 + 30 + 10 = 110** of them, each exactly once (`accounting.completeness.
  inventory`, with `derivation` recording which axis came from where).  A dead job
  whose `refs/` survives, a leftover `*_retry` directory, a never-launched domain, a
  fold launched with a shorter `--select-modes` and a whole missing density tier are
  then all one failure mode, not five.  `--force` is not optional: without it the
  scorer refuses to reuse a `cells.csv.gz` whose `.fp` sidecar does not fingerprint
  the runs now on disk (it recomputes where it can, refuses where it cannot).
- Every `predictions.json` and `manifest.json` the scorer reads is fingerprinted
  (size + sha256) into `accounting.completeness.inventory.fingerprints` — and
  re-read and re-compared after the cells pass, so a rewrite DURING the run lands in
  `inventory.inputs_changed_during_scoring` and voids the run
  (`run_invalid_reasons: [inputs_changed_during_scoring]`, exit 2) instead of being
  invisible; and the
  whole inventory is hashed into `meta.notes.run_inventory_fp` /
  `headline.run_inventory_fp` and into `<cells>.fp`.  The same is true of the ANSWER
  KEY and the vocabulary (`--functions`, `--maps`, `--ec-universe`, `--fanout`,
  `--annotation-ceiling`, `--templates`): fingerprinted into
  `inventory.fingerprints.answer_key`, named in `meta.notes.answer_key_paths`,
  hashed into `meta.notes.answer_key_fp` / `headline.answer_key_fp`, and re-verified
  after the cells pass.  `headline.scorer_code_fp` records the code that produced the
  numbers.  Keep those: they are the only
  way to prove after the fact WHICH BYTES produced a number (a concurrent session
  rewriting one prediction table mid-review moved picrust2 0.6955 → 0.6526 with
  every control green and nothing in the artefact able to show it).
- C1 re-run at `--n-holdout 500`, all three arms (A, B, C — arm C adds one
  hmmalign + one FastTree at 960 tips, ~1 min; `score_folds.py --mode
  recipe_control` hard-errors if `refs/recipe_C` was built but not run).
  C2 at d100 for every fold.
- Phenotype: ATP-vector arm over the whole eligible subset with `--n-genomes 1000`;
  gapfilled arm with `--n-genomes 200 --readouts atp,growth,essentiality,structure`.
  Each writes `meta.binding` (contract 1) into its `phenotype_summary.json`; pass the
  one the fold scoring should settle C11/C12 from with `--phenotype-summary` (the
  default is the file beside `--cells`, and with two phenotype arms in one tree the
  default is ambiguous), and re-run the fold scoring afterwards so C11/C12 are
  evaluated against a BOUND summary.
- Error-injection arm stays at smoke scale (0.5% / 1%).

---

## 12. Reporting discipline (A1 amendment I)

Publish **both** results side by side:
- the **controlled** k-fold result for the algorithmic claim, with the density
  sweep, both density slopes, the recipe-calibration ΔF1, the three baselines and
  the skill score;
- the **as-shipped** result for the practical claim, bracketed by
  `picrust2_comparison/data/annotation_dialect.json` (prFBA +0.050 F1 under BV-BRC
  truth, −0.204 under PICRUSt2's).

Every headline sentence names the three scope parameters that make the number
mean something: **the taxonomic rank blocked, the reference density in genomes,
and the coverage at which the comparison was made.**

**What this design may claim.** Under a single annotation vocabulary serving as
both methods' gene content and the answer key, and a reference identical between
methods consisting of ~9k species-dereplicated genomes, at genus-level novelty,
method X attains Y excess EC-F1 (and reaction-F1) over a nearest-surviving-relative
baseline, and dominates or fails to dominate the other method's risk-coverage
curve over coverage 0.5–1.0 — every number carrying the rank blocked, the
reference density, the coverage, and the ICC-corrected n_eff. The phenotype arm
may claim that model-to-model agreement through a fixed EC→role→template
bottleneck does or does not track EC-level agreement, relative to NN1-donor,
family-donor, constant and size-matched-null arms.

**What it may NOT claim.** Absolute functional-prediction accuracy (the annotation
ceiling is 0.904 recall / 0.743 precision for a genome's *own* annotation);
as-shipped superiority of either method (a database-size question, answerable only
by the uncontrolled run bracketed by the dialect test); generalisation beyond the
rank actually blocked; robustness to real sequencing data (queries are exact
substrings of reference 16S, with only the smoke-scale error-injection arm as
evidence); any NSTI-stratified cross-arm claim (NSTI is tree-specific, 0.42x under
our recipe); prFBA's species-tier behaviour or its species-dedup/P3 components;
or phenotype fidelity to real organisms absent experimental data.

---

## 13. Promote the audit reproduction scripts

A1's reproduction scripts (`/tmp/base1.py`, `/tmp/base2.py`, `/tmp/leak2.py`,
`/tmp/groups.py`, `/tmp/icc.py`, `/tmp/rc.py`) **are** the baselines, the leakage
check, the group-CV feasibility check, the ICC/n_eff computation and the
matched-coverage collapse — i.e. five of the required controls. Copy them to
`region_validation/shared_reference/scripts/design_audit_{baselines,leakage,groups,icc,riskcoverage}.py`
before `/tmp` is cleared, and have `score_folds.py` import the baseline and ICC
implementations from them rather than re-deriving the numbers. A2's persisted
evidence already lives at
`picrust2_comparison/data/kfold_audit/{audit_findings.json,recipe_control_results.json,density_results_full.json,include_self_8600.json}`
and is the source of every gate value in §7.

---

## 14. Amendment A-C12 (2026-10-27, post hoc): the phenotype arm on corrected templates

**Status.** Decided by the authors on 2026-10-27, after the pre-registered run failed
C12 (2026-10-18; `score_genus_pheno.log`, exit 3). This is a post-hoc change to the
pipeline, made with knowledge of that failure.
- The pre-amendment outputs remain the record of the pre-registered run and are reported
  as such: `data/phenotype_atp/`, `data/phenotype_gf/`, and `summary_genus.json` and
  `controls.json` as scored then (kept as `*.pre-c12fix`).
- **Diagnosis:** `EGC_C12_REPORT.md`.
- **Plan:** `C12_FIX_PLAN.md`.

**What changes (two things, and nothing else).**
1. **Templates.** Both phenotype passes read ModelSEEDTemplates commit `21a8084`, the
   five-reaction fix: catalase rxn00006, rxn10966, rxn00381, rxn00084 and rxn08167
   restricted to `>` in GramNeg/GramPos v7, Archaea v6 and Core v6.
   - **Pinned:** read from a detached checkout through `phenotype_fold.py --templates`,
     not through `phenotype_arm.MST`. MST points at a shared working tree that other
     work has since moved: since 2026-10-21 it carries rxn39860 removals.
   - **Recorded:** `meta.templates` holds the directory, the commit and the sha256 of
     every template file.
   - The pre-registered run read main `2294299`. Roles and reaction ids are identical in
     both commits; only the five reactions' direction fields differ.
2. **Gapfill reversals.** `phenotype_fold.py --integrate-reversals` applies the
   gapfiller's `reversed` set as well as its `new` set.
   - Applied as egc_rebuild.py's validated `reversed_only` variant does: open the bound in
     the direction the gapfiller chose.
   - `gapfill_block` inlines `MSBuilder.build_and_gapfill` from the ModelSEEDpy fork,
     which drops them. Upstream `MSGapfill.integrate_gapfill_solution` applies them.
   - Each gapfilled row now records `reversed_proposed` and `reversed_applied`.

**What does not change.**
- the genome selection: same seed, the same 1,000 and 200 genomes;
- the arms, pipelines, readouts, ladders and media;
- the solver settings;
- the cells table, fold plan, panel and answer key;
- the definitions and thresholds of C11 and C12: zero gated truth models, and at most 10%
  of models gated;
- the ModelSEEDpy code: the fork at `d884a50`, the commit the pre-registered run imported.
  - **Pinned** like the templates: a detached clone put first on `PYTHONPATH`.
    `phenotype_fold.py --expect-modelseedpy d884a50` stops a run that imports anything
    else, and `meta.libraries` records the path, commit and dirty flag of every library.
  - The venv's editable install is a shared working tree that has moved since 2026-10-21.
    Two of its later commits touch the build path: f395ae6 changes MSATPCorrection's
    expansion filter, and 8ac600b changes the builder's demand reactions.

The scorer's template-visible EC sets come from role names, so the re-scoring pass reads
main `2294299` (`--templates`, a pinned checkout). It thereby keeps the pre-registered
answer-key fingerprints exactly.

**Basis.** On the 218 gapfilled models of the 20 genomes holding the 49 gated rows, the
two changes together, rebuilt from template files:
- clear the gate in all 44 leaking models;
- open no new leak;
- let every model grow, with the passing models' median growth ratio 1.000.

The other 180 genomes of the gapfilled pass were not part of that validation. This run is
the test.

**Pre-commitment.** If C12 fails under this amendment, the failure is reported and the
leave-genus-out artefact stays invalid. Any further change would be a new, separately
reported amendment. No template or gapfill setting is tuned on this run's outcome.

**Runs.** Both run on the same template version:
- ATP pass: `--n-genomes 1000 --readouts atp,structure`, output `data/phenotype_atp_c12fix/`;
- gapfilled pass: `--n-genomes 200 --readouts atp,growth,essentiality,structure`, output
  `data/phenotype_gf_c12fix/`.

Then `score_folds.py` re-derives C11/C12 from the amended gapfilled summary
(`--phenotype-summary`, reusing the cells table). The genus and family summaries are
re-scored the same way.

**Execution record (2026-10-27).**
- **The first launch of both passes imported the moved fork** (c2e31a0), not `d884a50`.
  - It showed in the drafts, which were smaller than in the pre-registered run. In the 22
    drafts of genome 1005740.3, each lacked 1 to 4 demand reactions (8ac600b;
    `DM_cpd11416_c0` in all 22), and no other reaction differed.
  - Those partial outputs, 46 and 212 genomes, are kept as
    `data/phenotype_{gf,atp}_c12fix.msp-c2e31a0-partial/`. They are not used.
- **Both passes were restarted from zero on the pin.** Before the restart, a pinned
  rebuild with the pre-registered settings (main `2294299`, no reversals) reproduced every
  stored row exactly:
  - 86 of 86 rows for two ATP-pass genomes;
  - 80 of 80 rows for the gapfilled pass's two gated-truth genomes 586733.26 and
    1797949.3.
- **The amendment's basis predates both build-path commits.** The 2026-10-21 validation
  ran from 12:00 to 20:02, after the fork left `d884a50` (01:07) but before f395ae6 (22:03).
  The commits in between touch only the thermodynamics packages, a new element-closure
  package and compound formation energies, none of which the build reads.
- **`--maxtasksperchild` replaces each worker after a few genomes,** to bound memory over a
  long pass. Rows do not depend on worker history: genomes rebuilt in a fresh process match
  the pool's rows exactly, 84 of 84 (ATP pass) and 42 of 42 (gapfilled pass).

**Outcome (2026-10-28): C12 passes.**
- **ATP pass:** it reproduces the pre-registered pass exactly. All 43,140 rows and every
  statistic of its summary match.
- **Gapfilled pass:**
  - 0 of 8,614 rows are gated, against 49, 10 of them truth rows;
  - no model leaks anywhere that did not leak before, in the 20 validated genomes or the
    other 180;
  - every gapfilled model grows;
  - C11 passes.
- **Rescoring.** `score_folds.py` ran in fold mode against
  `data/phenotype_gf_c12fix/phenotype_summary.json`, with `--templates` set to the pinned
  main `2294299` and ModelSEEDpy pinned, which the scorer imports. The run is
  `scripts/c12fix_finish.sh`.
  - It reused the cells table, and every control passes, so the leave-genus-out artefact
    is valid (exit 0).
  - Every EC number is unchanged, as are the answer-key fingerprint (a9dc024d) and the
    scorer fingerprint (89a47bdf).
- **Details:** `C12_SCOPE.md`.

**The family summary is not rescored.** It was scored without the phenotype arm: it carries
no C11 or C12, and it was valid before the amendment. "The genus and family summaries are
re-scored the same way", above, holds for the genus summary only.
