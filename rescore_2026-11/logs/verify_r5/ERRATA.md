# Errata verification for the open plan errata and the unresolved values (2026-11-16)

Scope: the items handed over as open by `RUN5A_REPORT.md` (handoff 8), `RUN5B_REPORT.md` ("Flags" 3-4 and
"Not recomputed / unresolved") and `plan/PLAN.md` section 5. This was a report-only pass: no tracked file was edited,
and nothing was committed, checked out or reset.
- Manuscript: `~/Documents/prFBA_paper/prFBA.tex` at `2387ac2`. Line numbers below are that file's.
- Repo: `~/Documents/prFBA` at `a6812be` (branch `hardware-aware-alignment-dispatch`).

**Where it ran.** Everything ran on FM.
- Every computation was a light read (JSON, CSV, FASTA headers of two files) in a `systemd-run --user --scope -p MemoryMax=8G`
  scope, or a plain process under 1 GB.
- Each finished in seconds. The largest inputs were the two 16S MD5 maps (520 MB + 85 MB).
- Nothing ran on TF, and nothing wrote the cache or the repo.

**Scripts.** Copied to `~/Documents/prFBA_rescore_2026-11/logs/verify_r5/scripts/`: `md5count.py`, `md5count2.py`,
`refgenomes2.py`, `reach.py`, `reach2.py`, `table6b.py`, `table6c.py`, `as15.py`, `sr02.py`, `wobj.py`,
`anchors_cmp.py`, `anchors.py`, `depth_base.py`.
- Run each with `~/Documents/py_venv/bin/python`.
- `table6c.py` reads a scratch copy of `git show rescore-2026-11-start:region_validation/data/taxacc_region_summary.json`;
  regenerate that copy first.
- `wobj.py` is superseded by the stored `recomputed_contrasts.json` values (item 6).

Classes: **(a)** the printed value is wrong, and the correct value is given. **(b)** the printed value is right, and the
suspicion is a definitional or presentational mismatch. **(c)** the value cannot be resolved without a heavy recompute.

## Summary

| # | item | where printed | class | action |
|---|---|---|---|---|
| 1 | 460,117 unique sequences | tex L61, L135; method paper L7, L13, L23; DESIGN.md L52, L121 | (b) | keep; optional clarification |
| 2 | "1 to 5,265 contigs" | tex L822; `synth_mags.py` L11 | (a) max = **15,664**; (c) for 19.8% / 40.4% | replace max |
| 3 | marker reach 6,011 | tex L879; RESULTS_v2.md L220; `make_deck.py` L386 | (b) (6,011 has a congener; 6,010 of them are scored) | reword to 6,010 (recommended) |
| 4 | F16 vs three-panel deciles | tex L997-1000, F16 L1006, caption L1010-1011 | (b) every number correct; presentational | optional figure swap |
| 5a | method paper Table 6 V7-V9 Δ | method paper L242 (also V4 row L239) | (a) **+0.167** (and V4 family incl **0.963**) | replace |
| 5b | "~112k genomes with cached PGFams" | P2C README L266-267; method paper L278; RESULTS_v2.md L134 | (a) **~315,000** genomes (ratio roughly tenfold, not ~4x) | replace |
| 5c | RESULTS_v2 "73x" | RESULTS_v2.md L114, L246 | (a) **~14x** | replace |
| 6 | AS-15 0.026/0.014/0.05/0.03 | tex L530 | inputs changed; recomputed; **unchanged** | none |
| 6 | SR-02 inner margins 0.001-0.003 | tex L340, L595 | inputs changed; recomputed; **unchanged** | none |
| 6 | wobj weighted term median 0.11% | tex L310 | inputs changed; recomputed; **unchanged** (also 98.9%, 357/360) | none |
| 6 | depth-probe anchors 0.651 / 0.800 | tex L739 | inputs changed (85 anchors moved); recipe found; **unchanged** | none |
| 6 | depth-500 projections ~0.96 / ~0.49 | tex L742 | base re-checked; **unchanged** | none |
| + | P2C README "80.3% family-correct" | P2C README L248 | (a) **80.2%** (195/243) | replace |

---------------------------------------------------------------------------------------------------------------

## 1. "460,117" MD5-deduplicated unique sequences (tex L61, L135) vs the planner's 460,116: (b)

**Artefacts.**
- `codiffusion_bioreactor/model_inputs/16S_md5_seq.json` is the MD5 dedup map the alignment reference is built from (the
  method paper names it).
- `16S_md5_ID.json` has the same key set: 460,117 keys.
- Every key equals `md5(sequence as stored)`.
- No stored sequence contains an upper-case letter, so dedup on the upper-cased sequence gives the same 460,117.

**What the extra key is.** Exactly one key has a sequence shorter than 20 bp: `d41d8cd98f00b204e9800998ecf8427e`, the MD5 of
the **empty string**. Its header is `fig|2026780.1170.rna.42| SSU rRNA ... [Planctomycetes bacterium Dive96_bin3.198 |
2026780.1170]`.

**Why the embedding counts 460,116.**
- `manifest.json` / `index.parquet` come from `embed_16s.py`. It strips and upper-cases each record, then drops records
  shorter than `min_len` = 20 *before* dedup (`n_skipped_below_min_len` = 724).
- `index.parquet`'s shortest sequence is 60 bp, and the MD5 map holds no sequence of 1-59 bp.
- So the 724 skipped records are all empty, and they share that single empty-string digest.

**Arithmetic.**
- 460,117 = 460,116 non-empty unique sequences + 1 empty one.
- 460,117 − 816 (sequences under 200 bp, the empty one included) = **459,301**. That is the printed reference size, so
  the chain closes.

**Verdict.** Both counts are right for their definitions.
- The paper's sentence ("deduplicate by MD5 ... collapsing the set to 460,117 unique sequences") counts distinct digests
  over all 1,428,909 records, so it is literally correct.
- One of those "sequences" has no bases.
- `METHODS.md` L166, `README.md` L11 and method paper L139 (460,116) describe the embedding index and are also right.

**Proposed text** (optional; the number may stay):
- L135: "... so we deduplicate by MD5 digest of the sequence, collapsing the set to 460{,}117 distinct digests (one of
  them the empty sequence of 724 records that carry no bases)."
- Or, to print only real sequences, at L61 and L135: "460{,}116 non-empty MD5-deduplicated unique sequences". The ≥200 bp
  figure 459,301 is unaffected either way.

**Side note (wording only).** L135 says "MD5 digest of the uppercased sequence". The map's keys are MD5s of the stored,
all-lowercase sequence; `embed_16s.py` is the code that upper-cases. The count is identical, so this is optional.

## 2. "a panel spanning 1 to 5,265 contigs" (tex L822): (a) for the maximum, (c) for 19.8% / 40.4%

**Which genomes the sentence describes.** L820 defines the set: the 9,997 original-panel genomes whose deposited contigs
were fetched.
- They are `mag_validation/contigs/*.fna`, fetched by `mag_validation/fetch_contigs.py` from
  `region_validation/data/selection_10k.json`.
- The three genomes BV-BRC returned empty are excluded (`logs/fetch.log`: ok=9,997, empty=3).

**813.301 is in that set and was synthesized.**
- It is listed in `mag_validation/all_ids.txt` (line 9708).
- `mag_validation/manifests/b019.json` holds `813.301_c100_x0` with `total_src_bp` 9,275,382 and 1,556 synthetic contigs.

**Contig counts** (`grep -c '^>'` on the two files):

| genome | file | contigs | contigs ≥ 2.5 kb |
|---|---|---|---|
| 813.301 | `contigs/813.301.fna` | **15,664** | 1,462 |
| 114683.4 | `contigs/114683.4.fna` | 5,265 | 1,285 |

So 5,265 is the second-largest genome. No reading of "contigs" (all contigs, or only those above `synth_mags.py`'s 2.5 kb
`MIN_FRAG`) makes 5,265 the maximum while 813.301 is included.

**Supporting checks.**
- The planner (`plan/pipelines_mag_deploy.md` §A.6) reports "the largest original-panel contig file has 15,664
  contigs", with 5,265 second. That is evidently a ranking by contig count, not file size: 114683.4 (10.5 MB) is not
  among the largest files by size. Only the two files were re-counted here. The full scan below would also confirm
  that no third genome exceeds 15,664.
- The minimum of 1 holds: 15 of the 30 smallest files are single-contig.
- All 9,997 files were fetched in one 488 s run on 10-08 (13:03-13:11), and nothing records the 5,265 count's
  provenance. It may have been taken during that fetch: 813.301 landed at 13:10. It first appears in the
  `synth_mags.py` docstring (commit `b4f9db1`).

**Proposed text, L822:** "... would make degradation severity depend on each reference's assembly quality across a panel
spanning 1 to 15{,}664 contigs."
- The same edit applies to the `mag_validation/synth_mags.py` docstring at L11.
- `RESULTS_v1_superseded.md` L89 is historical; leave it.

**(c) The 19.8% closed and 40.4% under ten contigs** (tex L822, `synth_mags.py` L7-8, `make_deck.py` L286
`PANEL_CLOSED_PCT, PANEL_LT10_PCT`, slide text near L3443-3460).
- They cannot be checked without counting headers in all 9,997 files: 32 GB of FASTA.
- Command, if authorized:
  `systemd-run --user --scope -p MemoryMax=2G -p IOWeight=10 bash -c 'find ~/Documents/prFBA/mag_validation/contigs -name "*.fna" -print0 | xargs -0 grep -c "^>" > counts.txt'`.
- Cost: on FM this is CPU- and RAM-trivial, a few minutes of sequential reads. But it reads 32 GB, past this task's
  1 GB light-read rule.
- TF does not hold the contigs, so it would need a 32 GB push, or about 45 min of reads over its sshfs mount.
- It is (c) by the compute rule, not by cost. If the maximum was taken on an incomplete list, these two shares may be off
  too.

## 3. Marker reach "6,011" (tex L879; RESULTS_v2.md L220; deck `MAG_REACH[0]`): (b)

**What the code counts.** Producer: `manuscript_numbers/derived_numbers.py` `marker_reach()` (W6); the outputs below are
reproduced by `reach2.py`. The MAG skani hits are not re-scored, so this is not a re-score matter.

| quantity | value | definition |
|---|---|---|
| `n_congener` | **6,011** | organisms among `benchmark_truth.json`'s 9,746 scored sources that share a non-empty NCBI genus with another scored source |
| `extra.n_with_mag` | **6,010** | the same, restricted to organisms with a synthetic MAG |
| median best AF_query | 20.3% (20.285) | over the 6,010, among those with a surviving hit |
| share below the 15% AF gate | 45.1% | same population |
| no surviving hit | 15.9% | denominator 6,010 |

- The one organism dropped is **192389.10**: it is in truth but has no contigs in BV-BRC, so it has no MAG and cannot
  contribute.
- The planner's 6,010 is this `n_with_mag`. W6 prints `n_congener`, so both counts are the current code's output and they
  differ only in definition.
- The sentence quotes 6,011 as the population of statistics that are computed over 6,010.

**A literal reading of "a same-genus congener *in the reference*".**
- The reference is the 9,997 sketched genomes, which include 254 panel genomes with no amplicon truth. Their genera come
  from `selection_10k.json`, with `unclassified_*` placeholders treated as no genus; all 9,743 overlapping genera then
  agree with truth.
- Population: **6,027** organisms (+17 whose only genus-mate is one of the 254; −192389.10).
- Statistics over it: median AF 20.2%, 45.1% below the gate, 16.0% with no surviving hit.

**Proposed text (recommended; it keeps every statistic).**
- L879: "among the 6{,}010 scored panel organisms with a synthetic MAG that \emph{do} share a genus with another scored
  panel organism, the best surviving hit aligns over a median of just 20.3\% of the query, and 45.1\% fall below the
  15\% AF gate even though a genuine genus-mate is present."
- The wording names exactly W6's population: genus-mates among the 9,746 scored sources, not the 10,000-genome panel or
  the 9,997-genome reference, under which the count would be 6,027.
- L882 (15.9%) needs no change.
- In the same edit: deck `MAG_REACH = (6010, 20.3, 45.1)` (`make_deck.py` L386), and RESULTS_v2.md L220 "Among the 6,010
  scored organisms (with a synthetic MAG) that share a genus with another scored panel organism".
- W6 would then compare against `n_with_mag` (currently `n_congener`).

**Alternative.** Keep "in the reference" and print 6,027 / 20.2% / 45.1% / 16.0%. That needs a W6 recipe change.

## 4. F16 vs the three-panel deciles (tex L997-1011): (b), every number correct

**Prose values against `three_panel/space_deciles_three.json`** (deciles: 1 = sparsest). Every one reproduces:

| printed | computation | value |
|---|---|---|
| 19.2% | population[1] + population[2] | 0.0840 + 0.1079 = 0.1919 |
| 49.1% | original[1] + original[2] | 0.4914 |
| 36.6% | k-means[1] + k-means[2] | 0.3659 |
| 1.7× | excess ratio, (49.1 − 19.2) / (36.6 − 19.2) | 1.72 |
| 12.6% | population[10] | 0.1263 |
| 0.8% vs 4.1% | original[10], k-means[10] | 0.0084, 0.0411 |
| "five times" | 4.11 / 0.84 | 4.9 |
| 23.8% | random[1] + random[2] | 0.2383 |
| 12.4% | random[10] | 0.1242 |
| 7.7-12.6% (caption) | min and max population shares | 0.0772, 0.1263 |

- The original and k-means shares equal `two_panel/space_deciles.json` exactly.
- The inputs are re-score-invariant: the embedding is unchanged, and `space_deciles_three.json` passed the d20 baseline
  (R5.5).

**The mismatch is presentational only.**
- The paper's F16 is byte-identical to `two_panel/F_panel_space.pdf` (md5 `e3a53986...`). It plots the original and
  k-means panels only, with a flat 10% "population share" line.
- The random-panel shares (23.8%, 12.4%) and the per-decile population shares (19.2%, 12.6%) are not drawn in F16.
- They are drawn in the committed `three_panel/F_panel_space_three.pdf`: panels A-C are the original, k-means and random
  maps; panel D gives the shares with the population per decile. Its script asserts that the plotted shares match
  `space_deciles_three.json` within 1e-3.
- L997 already says the random panel only "can be placed on the same density deciles", so the text never claims F16
  shows it.

**Options.**
- **(A) Swap the figure.** Include `F_panel_space_three.pdf` as F16 and re-letter the caption: "(A-C) The 210{,}779
  clusterable 16S sequences ... with the original, $k$-means and random panels' picks overlaid ... (D) Share of each
  panel's on-frame picks per density decile, against the share of all clusterable sequences in that decile; the deciles
  are quantile bins of a tie-quantised distance and hold 7.7--12.6\% of sequences each."
  - Delete the sentence "The dashed 10\% line ... is only nominal".
  - Change L997 to "\Cref{fig:panelspace} shows where the three panels sit in that space."
- **(B) Keep F16** and add "(not drawn in \Cref{fig:panelspace})" after "the random panel can be placed on the same
  density deciles".

No number changes in either option.

## 5. Repo docs

### 5a. Method paper Table 6 (`region_validation/16S_pipeline_method_paper.md` L231-242): (a)

**Source.** `region_validation/data/taxacc_region_summary.json` (d20 tree, legacy reducer), anchor predictor, pooled
domain and tier, `per_rank[Family|Genus].correct`. Script: `table6b.py`.
- The V7-V9 row prints Δ **+0.193**. The correct value is 0.81930 − 0.65257 = **+0.16673 → +0.167**.
- The start tag gives the same (`table6c.py`), so this is a pre-existing typo, not a re-score change.
- Every other Δ in the table rounds the unrounded difference, as this one now does.
- A second pre-existing cell: V4 family accuracy (incl) prints **0.964**; the value is 7,122/7,392 = **0.963474 → 0.963**
  (start tag 7,070/7,338 = 0.963478 → 0.963). The 0.964 looks like a double rounding through 0.9635.
- Every other cell matches the current file at 3 dp. The prose "12-19 points" and "~0.19 for the short regions" still
  hold (range +0.110 to +0.191).

Replacement rows:
- L239: `| V4 | 0.963 | 0.790 | 0.599 | +0.191 |`
- L242: `| V7-V9 | 0.971 | 0.819 | 0.653 | +0.167 |`

### 5b. "~112k BV-BRC genomes with cached PGFams": (a)

**Where printed.**
- `region_validation/picrust2_comparison/README.md` L266-267 (caveat 2, "Database sizes differ by ~4×").
- `region_validation/16S_pipeline_method_paper.md` L278 ("differ ~4× in size (~112k ... vs 27,870)").
- The same stale count without the word "cache" is in `mag_validation/out/RESULTS_v2.md` L134 ("459,301 16S references
  from ~112,000 genomes").

**Current values** (`refgenomes2.py`, with the pipeline's own header parser `edlib_biopython_hits.parse_header`):

| quantity | value |
|---|---|
| distinct BV-BRC genomes behind the 459,301 ≥200 bp reference sequences | **314,557** |
| the same over all 460,117 digests | 315,050 (= DESIGN.md L121) |
| genomes with ≥1 full-length 16S (≥1,400 bp) | 147,286 (= DESIGN.md L52, so the parser matches the pipeline's) |
| gene-family cache keys (`stamps/R5A.final.json`) | 159,753 |
| cache keys with a non-empty PGFam list (E3 catalogue meta) | 158,790 |
| reference genomes with a cache key | 156,129 |

- "~112k" was the cache size after the June `prefetch_variant_genomes.py` run (SELF_RECOVERY_ANALYSIS.md L213 records
  "57k→112k").
- The caveat compares the databases the two methods select from. The current value of prFBA's is ~315,000 genomes, the
  manuscript's figure (tex L511, L585, L1224).
- The ratio to PICRUSt2's 27,870 is 314,557/27,870 = 11.3, "roughly tenfold", not ~4×.
- If a cache-size reading were wanted instead, it would be 159,753 (5.7×). It is not what selection draws on.

**Proposed text.**
- P2C README L266-267: "2. Database sizes differ roughly tenfold: prFBA selects from ~315,000 BV-BRC genomes (those
  behind its 459,301 unique ≥200 bp 16S sequences), PICRUSt2 from 27,870 reference genomes. Part of prFBA's precision
  advantage is database coverage, not algorithm."
- Method paper L278: "The databases differ roughly tenfold in size (~315,000 BV-BRC genomes behind the 459,301 references
  vs 27,870)."
- RESULTS_v2.md L134: "the amplicon arm 459,301 16S references from ~315,000 genomes."

**Leave as historical:** SPEC.md L40 and L1330 (pre-registered text), SELF_RECOVERY_ANALYSIS.md L213 (a true record of
June), RESULTS_v1_superseded.md L70.

### 5c. RESULTS_v2.md "73x": (a)

- **Where printed.** `mag_validation/out/RESULTS_v2.md` L114 ("replicates on a 73x larger database") and L246 ("A 73x
  denser reference (GTDB r226)").
- **Correct value.** 143,614 GTDB r226 species-representative sketches (`build_gtdb_bvbrc_map.py --expect-n 143614`;
  tex L168, L791) / 9,997 sketched panel genomes = **14.4**.
- **Where 73 came from.** 732,475 is the accession count of the bac120 + ar53 taxonomy files, and 732,475/10,000 = 73.2.
  `make_deck.py` L277-279 records this error; tex L882 already says "roughly 14-fold".
- **Proposed text.** L114 "replicates on a ~14x larger database (every CI overlapping)"; L246 "A ~14x denser reference
  (GTDB r226, 143,614 species representatives)".

### Extra finding: P2C README L248 "80.3% family-correct": (a)

- `data/depth_probe/depth_probe_results.json` stores 0.8025 at 4 dp; the underlying count is 195/243 = 0.80247 → **80.2%**.
- The tex prints 80.2% at L736 and 0.802 at L739, which is correct. Only the README rounded the stored 4-dp value half up.
- Proposed text: "anchors are 80.2% family-correct and 60.1% genus-correct."

## 6. Values RUN5B kept as printed ("Not recomputed / unresolved")

The gene-family cache repair of 2026-11-13 reached every input below except the taxonomy-only and alignment-only ones.
Each light value was recomputed, with a pre-repair control where a pre artefact exists (`$W/pre_repair/` or the start
tag).

**AS-15, "0.026 / 0.014 / 0.05 / 0.03" (tex L530).**
- Inputs changed: the d100 as-shipped comparison was re-scored in R2.4, and n went from 67,038 to 67,713.
- Recipe (`as15.py`): rxn minus EC conditional F1 for `picrust2` and `prfba_exclude` in
  `region_validation_d100/picrust2_comparison/data/comparison_summary.json`, taken over `by_region` (region × domain) and
  `by_stratum` (`ident_bin`, `nsti_bin`, `src_in_picrust2_ref`).
- On the pre-repair copy the recipe reproduces every print. Results:

| cell | pre | post | prints |
|---|---|---|---|
| Archaea max (V4-V5_944R, PICRUSt2) | −0.0264 | −0.0262 | 0.026 |
| PICRUSt2, identity bin < 0.865 | −0.0135 | −0.0137 | 0.014 |
| PICRUSt2, NSTI ≥ 0.5 (second distant bin) | −0.0115 | −0.0112 | (covered by "up to 0.014") |
| prFBA, NSTI unknown (n 9) | +0.0489 | +0.0489 | 0.05 |
| PICRUSt2, identity bin unknown (n 8) | −0.0266 | −0.0266 | 0.03 |

- Overall: rxn 0.8048 / 0.7216, which prints 0.805 / 0.722 as in the text.
- Every Bacteria, all-domain and source-in-reference cell is within 0.0087.
- **Printed values right; no edit.**

**SR-02 inner margins "0.001-0.003" (tex L340, L595).**
- Inputs changed: shared-reference re-score R3 / R4.2.
- Recipe (`sr02.py`): `nested_cv.fold_*.inner_f1` in `summary_genus.json` and `data_family/summary_family.json`.
- Legacy wins 5/5 folds at both ranks (`winner_instability` 0.0).

| | genus | family |
|---|---|---|
| margin over runner-up, post | 0.0014-0.0018 | 0.0017-0.0025 |
| margin over runner-up, pre | 0.0013-0.0018 | 0.0017-0.0025 |
| legacy minus exact-tie, post | 0.0014-0.0019 | 0.0022-0.0027 |

- **"0.001-0.003" still holds; no edit.**

**wobj "weighted term median 0.11%" (tex L310).**
- Inputs changed: the R4.3 wobj re-run.
- The value is already stored: `weighted_objective/recomputed_contrasts.json`, regenerated by `wobj_results.py` in R4.3,
  `descriptive.composition.plt1_term_share_of_objective_real_median_multi`.

| value | pre | post | prints |
|---|---|---|---|
| weighted p < 1 term, median share | 0.110% | 0.108% | **0.11%** |
| exchange share, median (multi) | 98.949% | 98.949% | 98.9% |
| ex-weight-100 arms with equal exchange term | 357/360 | 357/360 | 357/360 |

- All three print unchanged; **no edit.**
- Both the 0.11% and the 98.9% are medians over the 300 multi-genome amplicons. Over all 360 the exchange median is
  98.96%, which would print 99.0%, in both states. This is a consistent definition, not an error.

**Depth-probe answered-set anchors "0.651 / 0.800" (tex L739).**
- Inputs changed: the d20 P2C union selection `region_validation/picrust2_comparison/data/prfba_selection.json` was
  regenerated in R2.5.
- Against `$W/pre_repair`: the answered set is identical (52,103), but 177 exclude rep sets and 85 anchors differ.
- Recipe (`anchors.py`), found on the first attempt:
  - anchor = `exclude_0.987.reps[0]`, its taxon resolved by `_common.lineage_for`;
  - scored by `score_taxacc.rank_correct` against `benchmark_truth.json` `truth_lineage`;
  - pooled over the answered amplicons.
- Pre: family 0.799805 (n 32,868), genus 0.651343 (n 27,996), reproducing 0.800 / 0.651. Post: identical to 6 dp (the 85
  anchor swaps do not change correctness at family or genus).
- **Printed values right; no edit.** The recipe can go into W6 as is.
- **Flag (noted, not investigated).** The R2.5 regeneration moved 85 anchors and 177 rep sets of the d20 exact-tie
  selection while leaving the answered set identical. This is at odds with plan §6's statement that the exact-tie union
  reads no gene data.
  - `depth_probe_eval.py` passes a `gene_provider` and `gene_count_fn` to `select_representatives`. The P2C selection
    script was not checked, but it presumably does the same, so some step of the union (which one was not determined)
    consults gene data there.
  - Nothing printed changes. The provenance claim in plan §6 is too strong for this tree.

**Depth-500 projections "~0.96 / ~0.49" (tex L742).**
- RUN5B kept these as "taxonomy-based, unchanged inputs". That reason is inaccurate for the ~0.49: its base 0.375 is the
  legacy-reducer d20 population genus accuracy, whose file (`region_validation/data/taxacc_region_summary.json`) was
  rewritten in R1.4.
- Recomputed (`depth_base.py`; pooled genus `n_correct / n`, anchor, `exclude_0.987`):

| tree | pre | post | prints |
|---|---|---|---|
| d20 legacy (the base of ~0.49) | 0.3749 | 0.3749 | 0.375 |
| d100 legacy | 0.4720 | 0.4720 | 0.472 |

- With the probe's 0.601 and the recovery counts unchanged, ~0.49 stands.
- ~0.96 and the 0.766 / 0.905 coverages are alignment-only.
- **No edit.**

**The other rows of that table.**
- The AF-gate "1,386 incl. 19 gene-less" was dropped by Run 5B. Its replacement values (1,391; 26.5%; 0.712; 0.743; +45%)
  are tracked in `number_diff`.
- The L646-683 historical figures are labelled records of the pre-registered run or the earlier sample, so there is no
  current value to compute.
