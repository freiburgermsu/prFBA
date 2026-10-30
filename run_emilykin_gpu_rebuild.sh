#!/usr/bin/env bash
# Deployment B (the 3,950 EmilyKin ASVs), now at the depth-100 standard: align_hits.py
# --force-backend cuda --k2 100 --prefilter-k 500 --gpu-cand 200 --workers 12 -> tie-set PGFam
# cache -> prefetch_gene_families.py -> build_d100.py (--procs 8 --io-workers 8) ->
# deployment_stats.py, into EmilyKin/bvbrc_alignment_hits_d100/ and EmilyKin/synthetic_genomes_d100/.
# This script only forwards to deployment_d100/run_emilykin_d100.sh (all stages; resumable).
# The old 20-hit command (build from asv_top20_alignment_hits_gpu.json into synthetic_genomes/,
# --procs 48) is in git history (this file at a90dd17); its outputs are kept untouched.
# Long run: launch the whole chain with deployment_d100/launch_dep100.sh.
exec /bin/bash /home/freiburger/Documents/prFBA/deployment_d100/run_emilykin_d100.sh
