#!/usr/bin/env bash
# Deployment B at the depth-100 standard: the 3,950 EmilyKin ASVs (asvs.fasta, align_hits.py's
# default --fasta; MiDAS prior = edlib_biopython_hits.TAXONOMY_CSV, passed explicitly).
#   align_hits.py --force-backend cuda --k2 100 --prefilter-k 500 --gpu-cand 200 --workers 12
#   -> tieset gene cache -> prefetch_gene_families.py -> build_d100.py (parallel builder, --procs 8
#   --io-workers 8) -> deployment_stats.py.   See d100_pipeline.sh for the stages.
# --outdir is ALWAYS passed: align_hits.py's default --outdir is the 20-hit
# EmilyKin/bvbrc_alignment_hits/, which must not be overwritten. NEW outputs only:
#   /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits_d100/
#   /home/freiburger/Documents/EmilyKin/synthetic_genomes_d100/
# STAGE=align|build|all (default all). Long: launch through deployment_d100/launch_dep100.sh.
export NAME=emilykin_d100
export FASTA=/home/freiburger/Documents/EmilyKin/asvs.fasta
export TAXCSV=/home/freiburger/Documents/EmilyKin/bvbrc_embedding_hits/inputs/taxonomy.csv
export OUT=/home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits_d100
export SYN=/home/freiburger/Documents/EmilyKin/synthetic_genomes_d100
exec /bin/bash /home/freiburger/Documents/prFBA/deployment_d100/d100_pipeline.sh
