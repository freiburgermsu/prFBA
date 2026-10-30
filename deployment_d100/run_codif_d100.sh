#!/usr/bin/env bash
# Deployment A at the depth-100 standard: the 3,276 codiffusion-bioreactor V4-V5 ASVs.
#   align_hits.py --force-backend cuda --k2 100 --prefilter-k 500 --gpu-cand 200 --workers 12
#   -> tieset gene cache -> prefetch_gene_families.py -> build_d100.py (parallel builder, --procs 8
#   --io-workers 8) -> deployment_stats.py.   See d100_pipeline.sh for the stages.
# NEW outputs only (the 20-hit bvbrc_alignment_hits/ and synthetic_genomes*/ are never touched;
# the v7d community models read the 20-hit selections):
#   /home/freiburger/Documents/codiffusion_bioreactor/bvbrc_alignment_hits_d100/
#   /home/freiburger/Documents/codiffusion_bioreactor/synthetic_genomes_d100/
# STAGE=align|build|all (default all). Long: launch through deployment_d100/launch_dep100.sh.
export NAME=codif_d100
export FASTA=/home/freiburger/Documents/codiffusion_bioreactor/new_data/dna-sequences_codif_all.fasta
export TAXCSV=/home/freiburger/Documents/codiffusion_bioreactor/new_data/taxonomy_codif_all.csv
export OUT=/home/freiburger/Documents/codiffusion_bioreactor/bvbrc_alignment_hits_d100
export SYN=/home/freiburger/Documents/codiffusion_bioreactor/synthetic_genomes_d100
exec /bin/bash /home/freiburger/Documents/prFBA/deployment_d100/d100_pipeline.sh
