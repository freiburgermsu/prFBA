#!/usr/bin/env bash
# Deployment A (the 3,276 codiffusion-bioreactor new_data ASVs), BUILD stage, now at the
# depth-100 standard: tie-set PGFam cache -> prefetch_gene_families.py -> build_d100.py
# (build_synthetic_genomes_parallel.py + BV-BRC User-Agent fix, --procs 8 --io-workers 8) ->
# deployment_stats.py, from codiffusion_bioreactor/bvbrc_alignment_hits_d100/ into
# codiffusion_bioreactor/synthetic_genomes_d100/ (run run_newdata_gpu_align.sh first).
# This script only forwards to deployment_d100/run_codif_d100.sh (STAGE=build).
# The old 20-hit commands (prefetch into the shared cache + build from
# asv_top20_alignment_hits_gpu.json into synthetic_genomes/, --procs 48) are in git history
# (this file at a90dd17); their outputs are kept untouched.
STAGE=build exec /bin/bash /home/freiburger/Documents/prFBA/deployment_d100/run_codif_d100.sh
