# Stage 3 on the GPU path: the GPU ranks, Biopython only aligns the record

2026-10-27. The depth-100 trees were built with `--gpu-cand 500`. That re-scored all 500
edlib-shortlisted candidates with Biopython, so the GPU stage only re-ordered them. Two
changes to `align_hits.py`, each checked against the stored depth-100 records, make the
GPU do the ranking:

1. **The rescore's own order for the GPU pick.** Candidates are ordered by
   - GPU score, descending;
   - then edlib distance, ascending;
   - then reference index, ascending.

   The old pick, `np.argsort(-gpu_score)[:gpu_cand]`, was an unstable sort. It cut large
   equal-score tie clusters arbitrarily.
2. **Traceback only.** Biopython aligns just the top `--k2` (100) candidates in that order,
   for identity, aligned length and coordinates, which the score-only kernel cannot give.
   Nothing is re-scored.
   - Each alignment's score is checked against its GPU score. On any disagreement that
     ASV falls back to a Biopython re-score of its best `--gpu-cand` (default 200)
     candidates.
   - `run_stats.n_rescore_fallback` counts these fallbacks.
   - The CPU path is unchanged: Biopython scores all 500 itself.

## Checks

`gpu_cand_check.py` reads the run's saved shortlist (`edlib_top500_shortlist.npz`), its
GPU scores (`gpu_scores.npy`) and the stored records: 60,461 records, 6,046,100 hits.

| check | result |
|---|---|
| GPU score differs from the stored Biopython `align_score` | 0 of 6,046,100 hits |
| stored record differs from the top 100 of the shortlist in the order above, from GPU scores | 0 of 60,461 records |
| the old argsort pick, at `--gpu-cand` 100 / 150 / 200 / 250 / 300, would change | 34,524 / 3,008 / 2,180 / 1,482 / 1,004 records |

**End to end.** Stage 3 was re-run on 40 workers from the saved shortlist and GPU scores,
with the CUDA backend (`--shortlist-in data/hits/raw`). The times include writing the
records:

| setting | records vs the depth-100 run | stage-3 wall |
|---|---|---|
| depth-100 run: re-score all 500, align 100 (2026-10-20) | (the reference) | 7,128 s |
| `--gpu-cand 200` re-score, align 100 | byte-identical (`cmp`, records and summary CSV) | 4,414 s |
| **traceback only, align 100 (the standard)** | **byte-identical, 0 fallbacks** | **3,425 s** |

**Per call.** `gpu_cand_timing.py` ran on 150 random amplicons, seed 7, in one process
on an otherwise idle machine:
- a Biopython re-score costs 4.3 ms per candidate, the same at ranks 1-100, 101-200 and
  201-500;
- an alignment traceback with identity and coordinates costs 16.0 ms.

| setting | per amplicon |
|---|---|
| depth-20 run (re-score 100, align 20) | 0.75 s |
| depth-100 run (re-score 500, align 100) | 3.76 s |
| `--gpu-cand 200` re-score (align 100) | 2.46 s |
| traceback only (align 100) | 1.60 s |

So the Biopython stage's growth from depth 20 to depth 100 (0.75 s to 3.76 s per
amplicon) splits as follows:
- about three fifths is re-scoring 400 more candidates (1.7 s);
- about two fifths is aligning 80 more hits (1.3 s).

Dropping the redundant re-score removes 57% of the stage in this single-process
timing, and 52% end to end.

A first, cruder timing put the two costs at 2.4 ms and 25 ms. It differenced 500-candidate
and 200-candidate runs; it is superseded by the direct measurement above.
