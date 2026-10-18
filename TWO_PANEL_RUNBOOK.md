# Two-panel validation — state and runbook (2026-10-09)

## Part 1 — tau = 99 on the ORIGINAL panel: COMPLETE
`mag_validation/out/RESULTS_tau99.md`, `out/mag_summary_tau99.json`,
`out/paired_tau99.json`, `out/mag_cells_tau99.csv` (199,940 cells, all 20 batches).
Headline: resolve 0.307 -> 0.366, conditional F1 0.8068 -> 0.8230, genus 0.910 -> 0.895;
completeness still a null axis. tau=99 is a WEAKER exclusion than the species boundary
and is NOT matched to the amplicon arm's 0.987 -- see the caveats in that file.

## Part 2 — k-means panel: COMPLETE
`kmeans_panel/select_kmeans.py` -> `selection_kmeans.json` + `src16s_kmeans.fasta`.
k=10,000 over 210,779 full-length admissible 16S embeddings, k-means++ seed 1729,
15 minibatch Lloyd iterations, centroid-nearest genome per cluster with collisions
walked out (243). 9,506 Bacteria / 494 Archaea. Overlap with the original: 804 (8.0%).

## Part 3 — RUNNING / STAGED

### 16S arm (RUNNING, ~7 h from 13:30)
Isolated tree `region_validation_kmeans/` (BASE, RV and tmux names all redirected;
the original tree is verified untouched). Panel injected as data/selection_10k.json +
data/src16s_10k.fasta so the stock selector stage is skipped.
74,704 amplicons / 62,449 unique. Alignment in tmux `rv_align_full_km`, 48 workers,
CUDA, writes `data/hits/raw/align.done` on success.

**The driver was deliberately killed** (see incident below). After `align.done` appears:
```
cd ~/Documents/prFBA/region_validation_kmeans/scripts && bash run_validation.sh full
```
It skips the alignment (sentinel present) and runs only the POST stages.

### MAG arm (STAGED, not started)
`mag_validation_km/`: 9,996/10,000 contigs (37 GB; 4 genomes have no BV-BRC sequence),
804 hard-linked from the original panel. Reference sketch `skani_db_panel/` (4.3 GB).
Scorers now parameterized -- `score_mag.py --truth/--pgfams`, `aggregate.py --truth` --
with defaults preserving the ORIGINAL panel's behaviour exactly.
Use: `--truth region_validation_kmeans/data/benchmark_truth.json`
     `--pgfams kmeans_panel/pgfam_panel_kmeans.json`
Grid: completeness {100,90,70,50,30}%, contamination 0 only (the original's contaminated
cells are withdrawn, so 0% is the only cross-panel-comparable axis).
Persistence depth **-n 500** on BOTH panels, into sibling `hits_n500/` + `out_n500/`:
the largest ANI tie class is 304, so -n 100 truncates ~488 queries. NEVER overwrite
`mag_validation/hits/` -- that is the committed RESULTS_v2 basis.

## INCIDENT — PGFam cache poisoning (resolved)
`prefetch_pgfams.py` cached every failed fetch as `cache[gid] = []`, and BV-BRC had
begun returning HTTP 403 to the bare urllib User-Agent (verified: bare 403, any UA 200).
1,044 k-means-only entries were poisoned over 59 minutes before detection.
* ORIGINAL panel verified BIT-IDENTICAL: 10,000/10,000 present, 161 legitimate empties,
  0 lost (checked against `mag_validation/pgfam_panel.json`, built before the run).
* Repaired BY VALUE (never --force, which would refetch the 804 shared genomes at
  October vintage into the original arm's inputs). Snapshot:
  `bvbrc_cache/genome_gene_families.SNAPSHOT_before_repair.json`.
* UA added to `prefetch_gene_families.py` and `select_references.py` (`.bak-ua` backups);
  the failure path in BOTH trees' `prefetch_pgfams.py` now leaves the key ABSENT.
* Remaining 4,035 gaps filled. Cache 116,824 entries; both panels 10,000/10,000.
**The driver stays killed until the POST stage is run by hand, because
`select_references.py --gene-provider auto` would otherwise have re-poisoned the cache
unattended when the alignment ended.**

## Truth rebuild (k-means tree only)
`build_truth.py` derived the source taxid from the genome_id prefix; 126 k-means members
have a prefix absent from the prFBA taxdump. Patched to prefer the manifest taxon_id
(a provable no-op for the original panel: prefix == manifest for 10,000/10,000).
Verified: null-domain 109 -> 0, null-Phylum 124 -> 0, archaea 457 -> 476, 124 lineages
recovered (e.g. 1895630.3 had NO lineage, is Campylobacter porcelli). Prefix version kept
as `data/benchmark_truth.PREFIX.json`.

## Comparison constraints (from the 52 upheld review findings)
Report cross-panel rates as a ladder: crude -> common support -> domain- -> phylum-
standardized; only the standardized column is attributable to the selection scheme
(FullLength16S is +0.155 crude but +0.036 standardized; V4-V5_944R flips sign).
**Do NOT compare across panels:** population F1 (its denominator IS the panel);
MAG resolve rate from the panel-as-sketch arm (10,000 vs 5,048 taxids makes the k-means
sketch half as deep -- use the GTDB arm); per-phylum archaeal figures; the original's
structural 0.000 species accuracy beside any k-means number; amplicon-level counts
(join on genome_id, never best_16s_md5 -- 0/804 agree, upper vs raw md5).
Take every n from benchmark_truth.json, not placeholder counts.
