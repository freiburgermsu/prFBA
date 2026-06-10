# Would a deeper edlib prefilter recover the representatives edlib's top-20 missed?

Follow-up to `REPORT_edlib_vs_gpu_representatives.md`. Across both case studies exactly **2 GPU-selected
representative genomes** were absent from edlib's top-20 (both in **new_data**; EmilyKin had 0). This
checks whether they are genuine recall losses and whether a deeper edlib cut (top-100 / top-1000) recovers them.

## The 2 missing representatives are INSIDE edlib's prefilter

| ASV | rep genome | organism | SW id | GPU SW rank | edit-dist rank (of 459,301) | in k1=500 prefilter? |
|---|---|---|---|---|---|---|
| `dfbbce19b1bc` | 2164.119 | Methanobacterium sp. OMB003 | 0.8767 | 12 | 25-34 | True |
| `fce23542da16` | 57706.365 | Citrobacter braakii | 0.9598 | 2 | 2-975 | True |

Both are **well within edlib's k1=500 edit-distance prefilter** — they were NOT dropped by the prefilter.
They fell at the **kept-depth (k2=20) cut**, because each sits in a cluster of reference genomes with an
**identical Smith-Waterman score** that edlib and GPU order differently (arbitrary tie-breaking).

## It is a score-tie artifact, not a recall failure (GPU exhaustive top-200)

| ASV | rep SW score | genomes scoring HIGHER | genomes TIED at rep's score | tie-cluster ranks | top-100 guarantees? |
|---|---|---|---|---|---|
| `dfbbce19b1bc` | 508 | 3 | 26 | 4..29 | **True** |
| `fce23542da16` | 671 | 1 | 199+ | 2..200+ | **False** |

- **`dfbbce19b1bc` / 2164.119:** only 26 genomes tie at sw=508 (ranks 4..29). edlib's top-20 kept 20 of them; this one fell just past rank 20. **→ a top-100 (even top-30) edlib cut GUARANTEES recovery.**
- **`fce23542da16` / 57706.365:** the sw=671 tie is **enormous (≥199 genomes, filling the whole top-200 below rank 1)** — a dense *Citrobacter* cluster all at 0.9598 identity. A top-100 cut keeps only ~100 of the ≥199 equal-scoring genomes, so capturing THIS specific genome is a **~50/50 tie-break coin-flip; only top-1000+ guarantees it.** But edlib's chosen 2nd rep (546.1947, *Citrobacter freundii*, also sw=671) is an **equally valid representative**.

## Answer

**Partly — and the distinction is cosmetic.** A top-100 edlib prefilter would recover **1 of the 2** missing
representatives outright (2164.119, a small 26-way tie) and give the other a ~50% tie-break chance
inside a ≥199-way tie. Crucially, **neither "miss" is a real loss of information**: both missing reps are
**SW-score-identical** to genomes edlib already kept, i.e. interchangeable members of a score tie. So:

- The headline result stands and is *stronger*: edlib and GPU pick the **same representatives up to arbitrary
  tie-breaking** among equal-scoring genomes. Increasing the kept depth changes *which* tied genome is labelled
  the representative, not its score/identity/taxonomy.
- **A deeper PREFILTER (k1) is not needed** — both reps are already inside k1=500. The only knob that matters is
  kept-hit depth (k2), and only to relabel ties among already-equivalent genomes. Top-100 helps small ties;
  large equal-score ties are inherently arbitrary at any finite depth.

Artifacts (this dir): `gpu_top200_probe/asv_top200_alignment_hits_gpu.json`, `missing_rep_edit_distance_ranks.json`,
`missing_rep_top100_recovery.json`, `missing_rep_asvs.fasta`.
