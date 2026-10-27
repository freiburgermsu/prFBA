# `--gpu-cand 200`: the GPU prunes the shortlist without changing a record

2026-10-27. `align_hits.py` now defaults to `--gpu-cand 200`. On the CUDA path the GPU
passes only its 200 best of the 500 edlib-shortlisted candidates to the Biopython stage,
which keeps 100. The depth-100 trees were built at `--gpu-cand 500`, which rescores the
whole shortlist.

## The change

The 200 are now picked under the rescore's own total order:
1. score, descending;
2. edlib distance, ascending;
3. reference index, ascending.

The old pick was `np.argsort(-gpu_score)[:gpu_cand]`. That sort is unstable, so it split
large equal-score tie clusters arbitrarily.

The GPU score is bit-exact to Biopython's. So the new pick always contains the top 100 of
the whole shortlist under the order the rescore applies. The record is therefore the
full-shortlist record at any `--gpu-cand >= --k2`, and it still depends only on its
sequence.

## Checks

`gpu_cand_check.py` reads the run's saved `edlib_top500_shortlist.npz`, `gpu_scores.npy`
and the stored records (60,461 records, 6,046,100 hits).

| check | result |
|---|---|
| GPU score differs from the stored Biopython `align_score` | 0 hits |
| stored record differs from the canonical top 100 computed from GPU scores | 0 records |
| old argsort pick at `--gpu-cand` 100 / 150 / 200 / 250 / 300 would change | 34,524 / 3,008 / 2,180 / 1,482 / 1,004 records |

**End to end.** `align_hits.py --gpu-cand 200 --shortlist-in data/hits/raw` reran stage 3
on 40 workers from the saved shortlist and GPU scores, with the CUDA backend.
- `asv_top20_alignment_hits.json` and `asv_alignment_summary.csv` are byte-identical to
  the depth-100 run's (`cmp`).
- Stage 3 took 4,413.5 s, against 7,128.0 s for all 500 in the depth-100 run (made on a
  different day).

**Controlled timing.** `gpu_cand_timing.py` timed 150 random amplicons in one process.
- Stage 3 took 554.6 s at 500 candidates and 448.8 s at 200 (ratio 0.81). The 150
  records were identical.
- Per call, a Biopython re-score costs about 2.35 ms. An alignment traceback, with
  identity and coordinates, costs about 25.2 ms.
- So most of stage 3 is the 100 tracebacks a record needs, and pruning saves about a fifth
  of it. The larger end-to-end saving (38%) also reflects the different machine load.

## What Biopython is still for

The GPU kernel is score-only. Identity, matched bases, aligned length and reference
coordinates come from the alignment traceback, and selection runs on identity: the
0.987 / 0.945 / 0.865 bands, the coverage gate and self-exclusion.

The re-score of the 200 candidates repeats the GPU score exactly. A traceback-only stage 3
would sort by GPU score, align the top 100 and check each alignment's score against the GPU
score. That would remove the re-score, about 200 x 2.35 ms of the 2.99 s per amplicon. It is
not implemented.
