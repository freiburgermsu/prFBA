# Run 5C: applying the R5 verifiers' fixes (manuscript, decks, repo docs, number tracker), 2026-11-17

Executor: Run 5C agent. Inputs: `logs/verify_r5/{NUMBERS,CLAIMS,DECKS,ERRATA}.md` and their scripts,
`RUN5A_REPORT.md`, `RUN5B_REPORT.md`, `plan/PLAN.md`.

**Outcome.**
- **Paper:** one commit on `main` of `~/Documents/prFBA_paper`, **`38706bd`** (parent `2387ac2`), not pushed.
  It touches `prFBA.tex` and `images/prfba/F13_prfba_vs_picrust2.pdf`. The untracked `k_md5.txt` and
  `*.bak-*` files were left alone.
- **Build** of the committed tree (`pdflatex prFBA && bibtex prFBA && pdflatex prFBA && pdflatex prFBA`):
  98 pages, **0 undefined references or citations**, and **only the two pre-existing overfull hboxes**
  (23.69 pt at L133–138; 16.96 pt in tab:picrust2, now L557–574).
- **prFBA:** five path-limited commits on `hardware-aware-alignment-dispatch`, not pushed:
  `68947bf` (contig recount and doc errata), `efb2243` (provenance sources), `674287f` (decks, CI JSON,
  F13), `5dcfc0f` (number tracker) and `884f43c` (tracking the legacy as-shipped figures).
- **Not touched:** no `*_2026-10*.pptx`, the untracked user-edits deck, `ModelSEEDpy` or
  `ModelSEEDTemplates`.

---------------------------------------------------------------------------------------------------------------

## Where each step ran

Everything ran on this machine. LaTeX, git, text edits and figure and deck renders stay here by rule; the
rest was light.

| step | how | size |
|---|---|---|
| contig header count, all 9,997 files (32 GB read) | `systemd-run --user` unit, `MemoryMax=2G`, `IOWeight=10`, `ionice -c3 nice -n19`. It ran as a transient service rather than `--scope`, so it could run in the background; the limits are the same | about 3 min |
| length-aware contig scan (to test the old shares' definition) | same limits | about 4 min |
| cell-table recomputations (A, C1, E, F), bootstrap and AURC re-runs | `systemd-run --user --scope -p MemoryMax=8G` | under 1 min each |
| W6 `derived_numbers.py --state pre` / `--state post` (full) | unit with `MemoryMax=8G` | about 11 min each (high load from the other project) |
| W7 `collect.py` | scope with `MemoryMax=8G` | about 20 s |
| deck builds (all four) and `soffice` render | scope with `MemoryMax=8G` | under 1 min each |

Nothing ran on TF. No step needed more than 8 GB.

---------------------------------------------------------------------------------------------------------------

## A. Blocking numbers (recomputed from integer counts before editing)

Recipe: per-cell F1 = 2·tp / (n_truth + n_pred), and per-cell precision = tp / n_pred. Abstentions score 0,
and the 4-dp per-cell scores are not used.

| where (tex at 2387ac2) | quantity | exact value | old → new |
|---|---|---|---|
| L630 | family-rank exact-tie shortfall against nn1 (`data_family/cells.csv.gz`, density 1.0, EC; 0.6236661 − 0.7071907) | −0.0835246 | `(exact-tie $-0.083$, …)` → `(exact-tie $-0.084$, …)` |
| L630 | family-rank exact-tie largest region deficit (V1–V2) | −0.2475186 | `exact-tie $-0.077$ to $-0.247$` → `exact-tie $-0.077$ to $-0.248$` |
| L533 and the F13 caption (L577) | prFBA precision at NSTI ≥ 0.5 (`comparison_cells.csv.gz`, prfba_exclude, ok, n = 629) | 0.6054784 | `falls only from 0.822 to 0.606` → `0.605`; caption `falls only to 0.606` → `0.605` |

The other range ends still hold: −0.0773 (−0.077), legacy −0.2453 / −0.0757, −0.0811 (−0.081) and
0.8220 (0.822). The F13 panel A end label reads 0.61 at 2 dp, so it agrees.

## B. Cosmetic numbers

- **C1 (L78, L631).** The most precise source was available. `score_folds.cluster_bootstrap` is seeded
  (FROZEN seed 1729, 2,000 replicates over holdout groups), so I re-ran it on the cells.
  - Every one of the 43 stored EC contrasts (`delta` and both CI bounds, genus and family) reproduces at
    4 dp.
  - The unrounded bounds decide the print:
    - random-mask refusal credit at family rank, lower bound −0.004514 → **−0.005**. At L78 and L631,
      `([-0.004, 0.040])` → `([-0.005, 0.040])`.
    - Same sentence (L631), the margin over the prevalence floor: lower bound −0.052479 → **−0.052**,
      not the −0.053 the half-up convention gave on the stored −0.0525. `([-0.053, 0.018])` →
      `([-0.052, 0.018])`.
  - **Extra, found the same way:** at L613, PICRUSt2's lead over nn1 is `[0.030, 0.040]` → **`[0.029, 0.040]`**.
    The unrounded bound is 0.029482. It was 0.029491 before the repair, so this rounding error predates the
    re-score.
  - Every other printed genus-rank and family-rank EC bound was checked against the unrounded values, and
    all are right.
  - The recipe is now W6 group `sr_ci` (see J).
- **C2 (L901).** `(0.357 against 0.361, from 0.250 against 0.612)` → `(…, from 0.250 against 0.612; both gaps
  are differences of the unrounded means)`. 0.357 = 0.66063 − 0.30400, and 0.361 = 0.61159 − 0.25019.
- **C3 (tab:threepanel caption).** "The two right-hand columns are random minus original, with 95%
  bootstrap intervals" → "…intervals; each difference is taken between the unrounded standardized rates, so
  it can differ by 0.001 from the difference of the displayed columns." `compare_three_panels.py` sets
  `est` = round(phystd_x − phystd_y, 4). The verifier also checked that it matches the table in all 288
  cells.

## C. Wording

| tex (2387ac2) | old | new |
|---|---|---|
| L67 | `genus accuracy 0.60)` | `genus accuracy 0.60 vs.\ 0.65)` |
| L80 | `both phenotype passes were re-run and the 1{,}000-genome ATP pass reproduced` | `… re-run and, on the original sample, the 1{,}000-genome ATP pass reproduced` |
| L135 | `MD5 digest of the uppercased sequence, collapsing the set to 460{,}117 unique sequences.` | `MD5 digest of the sequence as stored (lower case), collapsing the set to 460{,}117 unique sequences, one of which is the empty sequence of the 724 records that carry no bases.` (checked: 460,117 keys, every key = md5 of the stored sequence, no upper-case letter, exactly one empty sequence and none of 1–59 bp, 816 under 200 bp so 459,301 ≥ 200 bp; the 724 are the manifest's `n_skipped_below_min_len`) |
| L523, L533 (×2), L625, L887, L1189 | two- or three-sentence lines | split before "prFBA". A sentence check over every line changed since 0214a14 leaves only L838 (a `\paragraph` title plus one sentence) and L1078 (ends ".)") |
| L529 | `By region the mean reverses only in the five regions …; for full-length 16S …` | `By region, under either rule, the mean reverses only in the five regions …; under the shipped exact-tie union, for full-length 16S …`. Legacy also trails in V3–V4 to V7–V9 and leads in FL, V1–V2, V1–V3 (from `comparison_summary.json`) |
| L643 | `That matches the EC result against nn1.` | `That matches the EC result against nn1 on the cells prFBA answers.` |
| L641 | `(Wilcoxon $p = 0.79$, $q = 0.79$)` | `(Wilcoxon $p = 0.79$, $q = 0.79$; over all 1{,}000 paired genomes 0.28 more, $p = 0.057$)`. Source: `phenotype_atp_c13/phenotype_summary.json`, `delta_atp_hamming_all_genomes` 0.284, `wilcoxon_p_all_genomes` 0.0570; same sign convention as nn1's +0.06 |
| L659 | `rescored the cells table and reproduced it byte for byte` | `reused the cells table after matching its fingerprints (inventory, answer key and scorer code)` (`logs/score_genus_c13.log` L15) |
| L741 | `on the as-shipped PICRUSt2 axis prFBA's population EC F1 rises from 0.609 to 0.723, still below PICRUSt2's 0.729.` | `on the as-shipped PICRUSt2 axis, under the shipped exact-tie union, prFBA's population EC F1 rises from 0.609 to 0.723, still below PICRUSt2's 0.729 (under the legacy reducer it rises from 0.616 to 0.734, above it).` (d20 legacy 0.6164, d100 0.7343) |
| L1069 | `… and in the direction opposite to the coverage gap: the original panel, the best covered, leaves the largest share …, by a fraction of a point (…); stratified by the number of deleted hits, neither difference is significant (…)` | `…, and their order is not stable: crude, the original panel leaves the largest share of surviving full-length records uncovered (pooled over regions it ties the $k$-means panel), but standardized to the pooled mix of deleted-hit strata the $k$-means panel leaves the most (2.3\%, against 2.0\% for the original and 1.8\% for the random panel), and the crude original-versus-random difference (Fisher $p=0.023$) does not survive stratification by the number of deleted hits (CMH $p=0.117$ …, $p=0.99$ …).` (`emptied_strata_d100.json` tests) |
| L1214 | `but score above it in EC space` | `but score at or above it in EC space` |
| L1197, L1226 | `a $+0.005$ lead as shipped` | `a $+0.005$ lead in the as-shipped comparison` |
| L592 (extra, CLAIMS cosmetic) | `23 of the 30 genus-rank fold references` | `… (five folds, two domain trees, three densities)` |

## D. Errata

**(1) Contig scan.** Command: `count_contigs.sh` (committed), the awk loop the verifier's `grep -c` implies,
plus contigs ≥ 2.5 kb and total bp. The run covered all 9,997 files; none was empty.

| count | min | max | single contig | under ten |
|---|---|---|---|---|
| all deposited contigs | 1 | **15,664** (813.301); 5,265 (114683.4) is second, 4,996 third | **1,729 = 17.3%** | **3,449 = 34.5%** |
| contigs ≥ 2,500 bp (MIN_FRAG) | 1 | 2,510 | 1,763 = 17.6% | 3,658 = 36.6% |

- The printed 19.8% / 40.4% do not reproduce under either count. They also fail on the 10,000
  denominator, on `all_ids_n500.txt` and on the 2,000 GTDB subset, so they have no recoverable
  provenance. The text and docs now give the all-contig figures.
- Committed: `mag_validation/data/contig_counts.tsv` (9,997 rows, 235 KB),
  `mag_validation/data/contig_counts_summary.json`, `count_contigs.sh` and `contig_counts_summary.py`.

| text | old | new |
|---|---|---|
| L822 | `19.8\% of panel genomes are closed and 40.4\% carry fewer than ten contigs, … spanning 1 to 5{,}265 contigs` | `17.3\% of panel genomes are deposited as a single contig and 34.5\% as fewer than ten, … spanning 1 to 15{,}664 contigs` |
| `synth_mags.py` docstring | 19.8% / 40.4% / 1 to 5,265 | 17.3% / 34.5% / 1 to 15,664, with the provenance |
| deck slide 34 card and notes (`PANEL_CLOSED_PCT, PANEL_LT10_PCT`) | 19.8 / 40.4, "one closed contig", "from one contig to thousands" | 17.3 / 34.5, "one contig", "from one contig to 15,664" |

**(2) Marker reach.**
- L879: `among the 6{,}011 panel organisms that \emph{do} have a same-genus congener in the reference` →
  `among the 6{,}010 scored panel organisms with a synthetic MAG that \emph{do} share a genus with another
  scored panel organism`. The statistics are unchanged.
- Deck `MAG_REACH` 6011 → 6010, with the slide 37 notes reworded.
- RESULTS_v2.md L220: same rewording.
- W6: new key `mag_reach_n_congener_with_mag`. `deck_MAG_REACH_0` now reads it, and its pre-state mismatch
  against the deck's old 6011 is documented in `KNOWN_PRE_MISMATCH`.

**(3)** L135: see C.

**(4)** F16: the L997-area sentence now reads "The random panel (not plotted in \Cref{fig:panelspace}),
drawn without reference to the space, tracks …".

**(5) Repo docs** (commit `68947bf`):
- Method paper Table 6:
  - V4 `| 0.964 |` → `| 0.963 |` (7,122/7,392);
  - V7–V9 `+0.193` → `+0.167`;
  - the prose that depended on them: "12–19 points" → 11–19 (V1–V2 is +0.110), and "~0.19" for the short
    regions → ~0.17–0.19.
- "~112k BV-BRC genomes with cached PGFams" / "~4×" → ~315,000 genomes behind the 459,301 references,
  "roughly tenfold" (314,557 / 27,870 = 11.3; verified with the verifier's `refgenomes2.py`), in:
  - P2C README caveat 2;
  - method paper L278;
  - RESULTS_v2.md L134.
- RESULTS_v2.md "73x" → ~14x at L114 and L246 (143,614 / 9,997 = 14.4).
- P2C README 80.3% → 80.2%.
- **Not changed:**
  - `compare.py`'s caveat string "BV-BRC ~112k genomes with PGFams" (both copies). It is written into the
    two committed `comparison_summary.json`; changing the code without regenerating them would split code
    and data. Fix it on the next regeneration.
  - The docstrings of `edlib_biopython_hits.py`, `run_prfba_fold.py`, `test_ref_db_override.py` and
    `adapt_as_shipped.py`, which say "~112k-genome" for the 16S set. These are core or scorer code (the
    scorer's code fingerprint gates cell reuse), so I left them; they are comments only.
  - The `strategy_tradeoff.py` cache-size comments, and the historical SPEC / SELF_RECOVERY /
    RESULTS_v1_superseded records.

**(6) VERIFICATION.md** (`objective_tuning/verification/`):
- `+0.0110 [+0.0031, +0.0222]` → `+0.0111 [+0.0031, +0.0223]`. `out_c4_c5.json` gives 0.011091
  [0.003121, 0.022253].
- The next line's `selected − shuffled (−0.0282)` → `(−0.0283; −0.0282 on 2026-11-03)` (−0.028279).
- Both are also added as rows of the file's 2026-11-16 delta table.
- Also checked and unchanged: face minimum +0.0203, and est_shuffled − two_tier +0.0060 [−0.0011, +0.0174].

## E. Depth-20 anchors: why 85 anchors and 177 sets moved with the answered set identical

**Cause.** The cause is the **P3 gene-richest tie-break**, not the `no_gene_set` exclusion.
- `region_validation/picrust2_comparison/scripts/prfba_selection.py` calls
  `select_representatives(..., gene_provider=provider, gene_count_fn=gcf)`. `gcf` is the PGFam-cache family
  count (`_meta.p3_gene_count_fn` = "PGFam cache family counts", in both the pre-repair and current files).
- In the exact-tie union path (`select_all`), `gene_provider` is never called. The `no_gene_set` exclusion
  sits in the legacy path (stage 7) only.
- `gene_count_fn` is read by the species de-duplication tie-break (`_gc_tiebreak`: among co-optimal
  conspecifics, keep the gene-richest assembly). It is also read for the `est_genes` / `gene_set_empty`
  annotation.
- The repair changed family counts, so it changed which assembly represents a species. It could not change
  whether an amplicon is answered, its tier, or the species set.

**Test** (`scratchpad/r5c/e_anchors.py`; pre-repair file
`$W/pre_repair/region_validation/picrust2_comparison/data/prfba_selection.json` against the live one):
- answered set 52,103 = 52,103, identical;
- 177 exclude-regime rep lists changed, all the same size, and **177 / 177 keep the same NCBI-species
  multiset**;
- **171 / 177 swaps involve one of the 1,006 genomes repaired on 2026-11-13**;
- of the 85 changed anchors, 83 stay in the same NCBI species and all 85 in the same genus. The two that
  change species are the documented residual coupling: an exact inter-species `_weighted_rep` tie broken by
  `genome_id`, after P3 moved a representative.
- This is why family and genus anchor accuracy are identical to 6 dp (the verifier's `anchors.py`).
- The other 6 swaps: their genomes are not in the 11-13 repaired list. I did not determine which count
  changed for them.

**Statements to correct.**
- I grepped the manuscript and every tracked `*.md` / `*.py` for "reads no gene data", "never reads
  gene_provider", "no gene data" and similar. No repository or manuscript statement says the exact-tie
  selection reads no gene data.
- The manuscript (L189) already says the tie-break "needs per-genome gene counts, which … the as-shipped
  comparison … and the depth probe … supply", and that it "changes which assembly represents a species, not
  the selection size or the species set". That is exactly what the test shows.
- The overbroad wording is only in the work directory's `plan/PLAN.md` §6 ("The exact-tie union never reads
  `gene_provider`"). That is literally true, but it omits `gene_count_fn`, which the P2C and depth-probe
  callers pass. `score_taxacc` and the union trees' taxonomy pass neither, so C-U's byte-identity proof
  stands. PLAN.md was not edited.

## F. Provenance (commit `efb2243`, `manuscript_numbers/rescore_sources/`, README included)

**Copied:**
- `pheno_sample_diff.json`, `pheno_attrib_{atp,gf}.json` and the two Run 4 scripts (`tools/`);
- `panel_len_check.json`;
- `changed_sets.json`;
- `r32_affected_refs.json`;
- `contig_cap_d6.json`.

**Computed and committed:**

| script → output | values | match the tex |
|---|---|---|
| `path_a_counts.py` → `path_a_counts.json` | 123 of the 161 `empty_ec` drops now have a complete EC; 122 pass length [1,200, 1,800] and ACGTN (1 fails length); **65 Archaea**, 57 Bacteria; 104 share a genus with an eligible genome; occupancy drop rate 535/9,407 = 5.69%, so 122 × 0.943 = **115.06** ("an estimated 115"). Control: the pre-repair `_cache_panel_ec.json` has all 161 empty | yes |
| `count_changed_cells.py` → `changed_cells.json` | **130** cells of the 17 genomes (24/47/24/31/4 by fold) against **63,871** scored (63,890 rows, 19 unscoreable_min_align) | yes |
| `count_contig_cap_queries.py` → `contig_cap_queries.json` | **10 of 49,990** MAG queries (5 tiers × 1733.9197, 1733.9535) | yes |

No tex change followed from F. All files are small (104 KB in total).

## G. Legacy bootstrap CI

- `gap_ci_d100()` now writes **`presentations/gap_ci_d100.json`** (committed) on every build. It holds gap,
  CI, the share of replicates on each side of zero, n amplicons, n source genomes, seed 20261116, 10,000
  replicates, and the input's sha256.
  - Legacy: +0.0054 [+0.000574, +0.010085], with 1.26% of replicates ≤ 0.
  - Union: −0.0062 [−0.010975, −0.001457].
- The notes of slides 19, 30 and 31 and condensed slide 11 cite the file.
- **Slide 31 face:**
  - each population-gap cell has a second line with its CI;
  - legacy reads "+0.005 (significant, narrowly)". The rule: the bound nearer zero is within 0.001 of it.
  - The union reads "−0.006 (significant)", because its upper bound is −0.0015.
  - Title: "One rule throughout: the exact-tie union, **at a population cost**". The qualifier is computed:
    it appears only when legacy leads on that axis.
- The code comment's pre-repair CI was corrected to what the code returns: (−0.0207/−0.0109,
  −0.0072/+0.0027).

## H. Decks

**Rebuilt** with `make_deck.py` and `make_deck_condensed.py [--brief|--panels]`.

**Text diffs against the decks of `a6812be`.**

Full deck:
- slide 9 notes and 10 notes: PR #20 open;
- slides 19 and 30 notes: CI file;
- slide 27 notes: the row-for-row statement said once;
- slide 29 face: "Safe to adopt; ModelSEEDTemplates PR #15 open";
- slide 31 face and notes;
- slide 34 face and notes: contig shares;
- slide 37 notes: marker reach;
- slide 42 notes:
  - the crude and standardized order;
  - "the unstratified Fisher test separates the stratified and random panels (p = 0.023), but stratified by
    the number of deleted hits the panels do not differ significantly (Cochran-Mantel-Haenszel p = 0.117
    …, 0.99 …)";
  - **also fixed:** a hard-coded "to within a point and a half" is now computed. It reads "2.3 points"; the
    re-scored residual is 0.954 − 0.931.

Condensed deck:
- slide 11 notes: CI file;
- slide 15 notes: PR #20.

The brief and brief-panels decks rebuilt with identical text; only the embedded chart workbooks' save
timestamps differed. The committed brief decks were restored, which is why `674287f` does not touch them.

**Render.** soffice → PDF → PNG. Every slide whose face changed was viewed: 29, 31 and 34. No overflow; the
slide 31 title fits on one line. The other changed slides changed only in their notes.

**Untouched:**
- the md5s of the four `*_2026-10*` decks and the untracked user-edits deck are unchanged;
- `presentations/prFBA_updates_2026-10.user-edits-1029.pptx` is still untracked.

## I. F13 panel B

`render_comparison_d100.py`: the annotation anchor moved from axes fraction (0.97, 0.955) to (0.97, 0.99)
(top-aligned).
- A 150-dpi pixel diff against the previous render differs only inside the annotation's box (x 707–1127,
  y 105–155 of 1968 × 633). The image size is unchanged, and so are panels A and C.
- The annotation now sits clear of the 0.98 / 0.99 / 0.97 labels. Those labels still touch the dashed
  1.00 line, as they did before the re-score.
- The PDF was copied to `prFBA_paper/images/prfba/F13_prfba_vs_picrust2.pdf` (md5 3bf75ff4…).
- The prFBA commit includes the script and `manuscript_figures/out/F13_prfba_vs_picrust2.pdf`. That
  directory is git-ignored, so the file was added with `-f`.

## J. Tracker agreement (`collect.py`, `derived_numbers.py`, `claims_spec.json`)

**collect.py.**
- **Rounding.** `fmt_value` rounds half away from zero on the decimal representation
  (`Decimal(repr(v))`, `ROUND_HALF_UP`). Python's binary-float formatting turned 0.8085 into 0.808. This
  covers f/pct/pp, int and the e-formats.
- **Boundary flag.** `print_boundary` flags a value stored at one decimal more than it prints and ending in
  5. `matches_printed` keeps its "boundary" note.
- **Precise sources.** New spec keys:
  - `precise`: the same quantity from a finer source. It is used for old and new when it resolves; the
    stored value is kept as `*_stored` and must agree at the stored precision, otherwise
    `*_precise_disagrees`.
  - `precise_tol`.
- **CSV derives.** `ratio` and `f1_counts`, so means are taken over cells' integer counts.
- **Current-manuscript check.** `--tex-base 0214a14` maps the spec's line numbers onto the current tex by
  line diff. Each value's `new_fmt` is looked for at the mapped lines ±2 (`new_fmt_in_tex`). The supporting
  keys are `tex_exception` (a deliberate different print, with its reason) and `check_new_print` (a value
  first printed after the base).
- **Report.** `number_diff.md` gains a "New values against the current manuscript" section and counts in
  the summary.

**derived_numbers.py (W6).**
- `_fmt` uses the same half-up rounding.
- New `sr_ci` group:
  - the unrounded EC bootstrap bounds of every controlled contrast, at both ranks, re-running
    `score_folds.cluster_bootstrap` on the cells. All 86 bounds reproduce the stored 4-dp values, pre and
    post.
  - The two SR-05 AURCs, re-running `risk_coverage` on `score_folds`' knob tables (both reproduce):
    0.812478 and 0.787537. These show that the tex's **0.812 / 0.788 are right** (C8 is resolved; half-up
    on the stored 0.8125 would have given 0.813).
- `magval_c*` keys: the tab:magval cells, unrounded.
- `mag_reach_n_congener_with_mag`.
- `KNOWN_PRE_MISMATCH` for `deck_MAG_REACH_0` and `sr_ci_genus_picrust2__vs__base_nn1_lo`.
- **Re-run checks.** Full `--state pre` and `--state post` re-runs: 421 numbers, 222 old plus 199 new.
  - Apart from `deck_MAG_REACH_0`, no existing number changed value or formatted string in either state.
  - The pre state mismatches only the four documented keys: the two nn1 values, `deck_MAG_REACH_0` and the
    0.030 bound.

**claims_spec.json.**
- Amended by `$W/tools/r5c_spec_precise.py` (idempotent; pre-amendment copy in
  `$W/control/r5c/claims_spec.before_r5c.json`).
- **263 values** now carry a precise source:
  - tab:magvs16s → W6 pairing;
  - tab:magval → W6 cells;
  - SR EC CI bounds → `sr_ci`;
  - SR-05 → `sr_aurc`;
  - three-panel `est` → difference of the stored unrounded phystd;
  - tab:emptied shares → n/N;
  - as-shipped population F1, precision and coverage → cells by counts;
  - depth-probe and deployment shares → their counts;
  - FR-02, FR-03 and SR-07 → shared-reference cells by counts.
- **Formats:**
  - PH-04 / disc.atp_vs_nn1 p and q → f2;
  - the Spearman p → e1;
  - PH-02 and disc.atp_vs_picrust2 gain the q range's minimum;
  - "within 0.016" → a value (f3), since it no longer holds as a bound and the tex prints 0.017.
- **New printed values added:**
  - PICRUSt2 all-genome sensitivity, 0.28 and p 0.057;
  - standardized uncovered rates, 2.0 / 2.3 / 1.8%;
  - Fisher p 0.023;
  - the legacy reducer's as-shipped population EC F1 (0.734 at depth 100, 0.616 at depth 20) and its +0.005
    lead over PICRUSt2 (AS-09, disc.margin_reversal; commit `884f43c`). This closes RUN5B's Flag 1: these
    values had been printed since R5B but were untracked. All are from the cells' counts and all print as
    the tex does.

**Result** (`number_diff.{json,md}` regenerated):
- 178 items, 1,580 values; 1,580 old and 1,580 new resolved; 0 errors, 0 pending.
- 1,544 new values print exactly as the tex does at their lines.
- 33 are not checkable: strings, lists, values the tex does not print, and numbers written as words.
- **3 differ, each with its reason recorded:**
  - panels.surviving_uncovered "FL strata uncovered (low)": 0.0, printed in words ("down to none of the
    original panel's 274 records");
  - AS-12 "wins per loss (max)": 7, printed in words ("between three and seven");
  - abs-phenotype-controls: the re-run's 43,174 rows. The abstract keeps the pre-registered pass's 43,140 on
    purpose; 43,174 is printed in sec:sharedref.
- **0 unexplained.**
- The verifier's 10 producer-wrong rows and the three blocking ones all now equal the tex. Examples:
  - V4 MAG 0.809;
  - V1–V3 16S 0.704;
  - V7–V9 median Δ +0.125;
  - 1% FL MAG 0.785;
  - FL population F1 0.777;
  - source-absent coverage 0.999;
  - random 50–99 share 14.5;
  - panel-as-reference Δ −0.027;
  - SR-07 −0.092;
  - SR-08 bound 0.029;
  - −0.084, −0.248 and 0.605;
  - the p/q/Spearman/McNemar formats.
- **Still printed by convention: 3 values.** No finer source exists; the tex follows the paper's half-up
  convention. Each is a three-panel bootstrap bound that `compare_three_panels.py` stores as
  `round(x, 4)` from an unrounded float it does not write:
  - genus conditional accuracy, V3–V4, random − original, lower bound −0.0385, printed −0.039;
  - gene-capture conditional F1, V7–V9, lower bound −0.0345, printed −0.035 (it appears in two items).
- **Deciding them would need** re-running `compare_three_panels.py --suffix _d100` without the rounding.
  That reads all three depth-100 trees, and its earlier run sat inside a 24.7 GB-peak job. TF does not hold
  those trees. Not done.
- **Old-side mismatches: 9, all with recorded reasons.**
  - the two nn1 values without a surviving recipe (each in two items);
  - the "about 45%" E3 share (two items);
  - the old Spearman print at one significant digit (two items);
  - the old 0.030 bound, a pre-existing rounding error (see B).
- Since the precise sources and the new kinds of value went in, 716 values change at print
  precision, in 117 items (was 691 / 116).

**Tests.**
- `tests/test_collect_numbers.py`: 48, adding boundary formatting, `print_boundary`, the line map, and an
  end-to-end precise-source / derive / tex-check test in a throw-away paper repo.
- `tests/test_derived_numbers.py`: 14, adding half-up formatting, the bootstrap and AURC recipes.
- **62 passed.**

---------------------------------------------------------------------------------------------------------------

## Commits

| repo | commit | content |
|---|---|---|
| prFBA_paper (`main`) | `38706bd` | prFBA.tex (A, B, C, D1–D4), F13 PDF |
| prFBA | `68947bf` | contig recount (script, TSV, summary), synth_mags docstring, RESULTS_v2.md, method paper, P2C README, VERIFICATION.md |
| prFBA | `efb2243` | `manuscript_numbers/rescore_sources/` |
| prFBA | `674287f` | make_deck.py, gap_ci_d100.json, the full and condensed 2026-11 decks, render_comparison_d100.py, F13 output |
| prFBA | `5dcfc0f` | collect.py, derived_numbers.py, claims_spec.json, derived_numbers.{pre,post}.json, number_diff.{json,md}, both test files |
| prFBA | `884f43c` | claims_spec.json and number_diff.{json,md}: the legacy as-shipped figures tracked |

Nothing was pushed.

## Changes beyond the fix list (each can be reverted on its own)

**Manuscript.**
- L613: `[0.030, 0.040]` → `[0.029, 0.040]`, and L631: `[-0.053, 0.018]` → `[-0.052, 0.018]`. Both come
  from the bootstrap re-run that settled C1; the unrounded bounds are 0.029482 and −0.052479.
- L592: the 30 references are explained ("five folds, two domain trees, three densities").
- L1069: the crude Fisher p = 0.023 is now stated, beside the CMH tests.

**Repo docs.**
- Method paper: the Table 6 prose ("12–19" → 11–19; "~0.19" → ~0.17–0.19), beyond the two rows asked for.
- VERIFICATION.md: the adjacent "selected − shuffled −0.0282" → −0.0283, and two rows in its delta table.

**Decks.**
- Slide 42 notes: "to within a point and a half" (hard-coded; wrong after the re-score) → computed,
  "2.3 points".
- Slide 31: a CI line under both population-gap cells; the title qualifier "at a population cost".
- The brief and brief-panels decks were rebuilt but restored to their committed bytes, since only the
  workbook timestamps differed.

**Git.** `git add -f` was used for the F13 PDF under the git-ignored `manuscript_figures/out/`.

**Tracker.**
- W6's own `_fmt` changed to half-up; no existing W6 value changed its formatted string.
- New `check_new_print` values for numbers printed since R5B.

## Could not fix / left open

- **Three-panel bounds.** The three bootstrap bounds above are printed by convention; they need an unrounded
  `compare_three_panels.py` re-run.
- **19.8% / 40.4%.** Their origin is unknown. The text now uses the recount.
- **Code strings.** `compare.py`'s "~112k" caveat string and the "~112k-genome" docstrings in core and scorer
  code are left for the next regeneration (see D5).
- **The 6 E swaps.** They involve no genome from the 11-13 repaired list; I did not determine which count
  moved for them.
- **PLAN.md §6.** Its "never reads gene_provider" wording omits `gene_count_fn` (see E). It is a workspace
  document and was not edited.
