# The C12 fix run: what it involves and what it changes

2026-10-26. Background: `EGC_C12_REPORT.md` (diagnosis) and
`data/template_direction_audit/validation_summary.json` (validation).

**Outcome (2026-10-28).** The plan was run as SPEC amendment A-C12, with ModelSEEDpy pinned
to the pre-registered commit (SPEC section 14, execution record).
- C12 passes: 0 of 8,614 gapfilled rows are gated.
- The leave-genus-out artefact is valid.
- The ATP pass is unchanged.
- Results and scope are in `C12_SCOPE.md`.

This plan has been run. When it was written it was still pending, because it changes
pre-registered pipeline settings: it needed the authors' go-ahead and a SPEC amendment.

## Why C12 fails today

C12 is a hard control with three legs. It fails on leg (b), zero tolerance for gated
truth models:
- 10 truth rows (5 genomes × 2 pipelines) make ATP on water and protons;
- leg (c) passes: 49 of 8,614 rows are gated, 0.57% against a 10% ceiling.

Because C12 is hard, `summary_genus.json` marks the whole leave-genus-out artefact
invalid. That includes the EC results, which do not depend on the phenotype arm at all.

## The two changes

1. **Templates.** Run the phenotype arm on ModelSEEDTemplates `21a8084`
   (`fix-irreversible-reaction-directions`), the five-reaction fix:
   - catalase rxn00006, rxn10966, rxn00381, rxn00084 and rxn08167 restricted to their
     physiological direction `>`.

   Point `phenotype_arm.MST` at a dedicated worktree of that commit. Do not use the
   shared clone: `~/Documents/ModelSEEDTemplates` is checked out on another session's
   branch (`remove-rxn39860`).
2. **Apply the gapfiller's reversals.** `phenotype_fold.gapfill_block` inlines
   `MSBuilder.build_and_gapfill` and, like it, integrates only `res["new"]`, never
   `res["reversed"]`. Add the reversed set: open each reversed reaction in the direction
   the gapfiller chose, as `egc_rebuild.py` does for the validated `reversed_only`
   variant. The same bug is in the ModelSEEDpy fork's `MSBuilder.build_and_gapfill`
   (freiburgermsu/ModelSEEDpy d884a50), and is worth fixing there too.

Validated on the 218 gapfilled models of the 20 affected genomes, against the original
templates rebuilt with reversals applied:
- all 44 leaking models pass, and no new leak appears;
- every model grows;
- the passing models' median growth ratio is 1.000.

The other 180 gapfilled genomes were not part of that validation.

## The run

| Step | What | Size | Wall time (estimate) |
|---|---|---|---|
| 0 | SPEC amendment | — | — |
| 1 | ATP pass, re-run | 1,000 genomes, 44,348 models, no gapfilling, 1-2 min per genome per worker | about 2-3 h on 16 workers |
| 2 | gapfilled pass, re-run | 200 genomes, 8,614 models, 13-15 min per genome per worker | about 4-6 h on 8-10 workers (about 7 GB per worker) |
| 3 | C11 identity control and `score_folds` with the new phenotype summaries | — | under 1 h |
| 4 | Manuscript and deck: phenotype numbers, the C12 paragraph, and the validity statements | — | — |

**The SPEC amendment** records:
- the two changes, and that they were decided after C12 had failed;
- the diagnosis and the validation;
- that the pre-amendment phenotype results remain the record of the original
  pre-registered run.

**Why the ATP pass must be re-run too.** It does not gapfill, but the template change
also alters its un-gapfilled drafts (catalase's draft direction, for one), and both passes
must use one template version.

All runs go through `systemd-run --user` units with LimitNOFILE and MemoryMax, as for
the depth-100 campaign.

## What it means for the study

- **C12 should pass.** All 10 gated truth rows are in the 20 validated genomes, and the fix
  clears every one of them. The other 180 genomes had no gated truth model before. The
  template change could in principle create a new leak there; the run will show it, and
  C12 would then still fail.
- **The leave-genus-out artefact becomes valid** if every other hard control still passes.
  Its EC numbers do not change (population EC F1 0.687 against PICRUSt2's 0.779,
  Δ -0.092). They become quotable without the "invalid under the design's rule" caveat
  the manuscript carries now.
- **The phenotype results move.**
  - **The primary ATP readout:** all 1,000 genomes, small shifts expected through the
    drafts.
  - **The gapfilled secondary readouts** (growth, essentiality, gapfill burden) move more:
    - the 20 affected genomes return to every arm's denominator, whereas today a gated
      truth model drops its genome from all arms;
    - growth now counts the reversals the gapfiller asked for. In the 20 genomes, 23 of
      218 models had a dropped reversal, and all 10 original non-growers were such
      cases.
- **The headline phenotype conclusion is unlikely to flip, but must be re-checked.** It is
  that the readout detects no mean difference between prFBA and PICRUSt2.

The template PRs to ModelSEED are independent of this run. The run uses the same commit
the first PR proposes.
