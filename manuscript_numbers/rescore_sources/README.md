# Sources of re-score values the manuscript prints that no other committed artefact holds

The 2026-11 re-score (work directory `~/Documents/prFBA_rescore_2026-11`, not a repository) printed
some values in `prFBA_paper/prFBA.tex` whose only source was that work directory or a git-ignored
file. This folder holds small copies or summaries of them, and the scripts that produce the ones
computed here (2026-11-17, re-score R5C). Line numbers are those of `prFBA.tex` after R5C.

| file | backs (manuscript, sec:sharedref unless stated) | how it was produced |
|---|---|---|
| `pheno_sample_diff.json` | the phenotype re-draw: 200 of the 7,855 candidate genomes change their EC-F1 difference; 13 of the 15 strata change membership (15 strata, 2 identical); 442 of 1,000 ATP and 62 of 200 gapfilled genomes kept (`n1000.picked_both`, `n200.picked_both`) | `tools/r4_pheno_sample_diff.py` (copied from `$W/tools/`), re-score Run 4; was `$W/control/r4/` |
| `pheno_attrib_atp.json`, `pheno_attrib_gf.json` | a kept genome's rows are identical unless its own inputs changed: 106 of 19,072 ATP rows and 0 of 2,684 gapfilled rows differ (`rows_compared_genomes_in_both`, `rows_changed`) | `tools/r4_pheno_attrib.py`, re-score Run 4; was `$W/control/r4/` |
| `panel_len_check.json` | Path A: 123 of the 161 genomes excluded for an empty EC set were stored fetch failures (`original.n_old_empty` 161, `n_cur_empty` 38) | the re-score planner's read-only check of the panel's PGFam lists against the repaired cache (`$W/expected/rescore_work/`); its script was not kept, `path_a_counts.py` reproduces the 123 |
| `path_a_counts.py` -> `path_a_counts.json` | Path A: 123 now have ECs; 122 of them pass the length and character filters, 65 of those archaea; "an estimated 115 more genomes" (122 x (1 - 535/9,407 occupancy drop rate) = 115.1) | `python manuscript_numbers/rescore_sources/path_a_counts.py` (reads `panel_eligible.json`, the regenerated `_cache_panel_ec.json`, `_cache_panel_md5_seq.json`, `selection_10k.json`; control: all 161 are empty in the pre-repair copy) |
| `changed_sets.json` | Path A: the 17 eligible genomes whose EC sets the repair changed (`changed_eligible`); 140 panel genomes changed (`changed_panel`) | the re-score planner (Run 3 refresh reproduced it); was `$W/expected/` |
| `r32_affected_refs.json` | Path A: 23 of the 30 genus-rank fold references carry a changed genome (`data.affected`; also family and C2) | re-score Run 3, `$W/tools/r32_checks.py`; was `$W/control/r3/` |
| `count_changed_cells.py` -> `changed_cells.json` | Path A: 130 of the 63,871 test cells (24 / 47 / 24 / 31 / 4 by fold) | `python manuscript_numbers/rescore_sources/count_changed_cells.py` (reads the git-ignored genus `cells.csv.gz`) |
| `contig_cap_d6.json` | D6 (sec:twopanel, the random panel's MAG arm): 1733.9197 and 1733.9535 hold exactly 25,000 contigs on disk against 31,169 and 36,001 at BV-BRC; two sludge genomes of the real-MAG arm are capped too | BV-BRC Content-Range count queries, 2026-11-15; was `$W/expected/` |
| `count_contig_cap_queries.py` -> `contig_cap_queries.json` | D6: 10 of the arm's 49,990 MAG queries come from those two genomes | `python manuscript_numbers/rescore_sources/count_contig_cap_queries.py` (reads the git-ignored `mag_validation_rnd/manifests/b*.json`) |

Related, elsewhere in the repository:
- `presentations/gap_ci_d100.json`: the legacy and union as-shipped population-gap bootstrap CIs the
  decks cite, written by `presentations/make_deck.py` `gap_ci_d100()`.
- `mag_validation/data/contig_counts{.tsv,_summary.json}`: the panel's deposited contig counts
  (sec:magval "1 to 15,664 contigs", 17.3% single-contig, 34.5% under ten).
- `manuscript_numbers/derived_numbers.post.json` `sr_ci_*`, `sr_aurc_*`: unrounded controlled-benchmark
  CI bounds and AURCs (the summaries store 4 dp), which decide the printed -0.005, -0.052 and 0.029.

`tools/` holds the two Run 4 scripts as they ran; they import the work directory's other tools by
absolute path and are kept as a record, not as runnable code in this repository.
