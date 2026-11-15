#!/usr/bin/env bash
# Deployment B at the depth-100 standard: the 3,950 EmilyKin ASVs (asvs.fasta, align_hits.py's
# default --fasta; MiDAS prior = edlib_biopython_hits.TAXONOMY_CSV, passed explicitly).
#   align_hits.py --force-backend cuda --k2 100 --prefilter-k 500 --gpu-cand 200 --workers 12
#   -> tieset gene cache -> prefetch_gene_families.py -> build_d100.py (parallel builder, --procs 8
#   --io-workers 8) -> deployment_stats.py.   See d100_pipeline.sh for the stages.
# --outdir is ALWAYS passed: align_hits.py's default --outdir is the 20-hit
# EmilyKin/bvbrc_alignment_hits/, which must not be overwritten. NEW outputs only:
#   /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits_d100/
#   /home/freiburger/Documents/EmilyKin/synthetic_genomes_d100/  (shipped, mixed; see below)
#   /home/freiburger/Documents/EmilyKin/synthetic_genomes_bvbrc_d100/  (BV-BRC only, SYN)
# Since 2026-11-13 every build is BV-BRC PATRIC annotation only and d100_pipeline.sh refuses to
# build into a non-empty SYN. SYN below is the BV-BRC-only build (synthetic_genomes_bvbrc_d100,
# 2026-11-13: STAGE=build re-selected from the shipped hits with the repaired tie-set cache,
# bvbrc_cache/CORRECTIONS_2026-11-13.json; audit + stats in SYN). It exists now, so build/all are
# refused until it is moved away. The SHIPPED mixed PATRIC + RefSeq build (2026-10-29/30, the
# one uploaded to KBase) stays in /home/freiburger/Documents/EmilyKin/synthetic_genomes_d100,
# with its audit and stats in OUT; it is reproducible only from git at or before 3ce12d0.
# STAGE=align|build|all (default all). Long: launch through deployment_d100/launch_dep100.sh.
export NAME=emilykin_d100
export FASTA=/home/freiburger/Documents/EmilyKin/asvs.fasta
export TAXCSV=/home/freiburger/Documents/EmilyKin/bvbrc_embedding_hits/inputs/taxonomy.csv
export OUT=/home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits_d100
export SYN=/home/freiburger/Documents/EmilyKin/synthetic_genomes_bvbrc_d100
exec /bin/bash /home/freiburger/Documents/prFBA/deployment_d100/d100_pipeline.sh
