#!/usr/bin/env bash
# Deployment A (the 3,276 codiffusion-bioreactor new_data / codif_all ASVs), ALIGNMENT stage,
# now at the depth-100 standard: align_hits.py --force-backend cuda --k2 100 --prefilter-k 500
# --gpu-cand 200 --workers 12 -> codiffusion_bioreactor/bvbrc_alignment_hits_d100/.
# This script only forwards to deployment_d100/run_codif_d100.sh (STAGE=align).
# The old 20-hit command (gpu_align.py --full --topk 20 -> bvbrc_alignment_hits/
# asv_top20_alignment_hits_gpu.json) is in git history (this file at a90dd17); its outputs are
# kept untouched. Long run: launch the whole chain with deployment_d100/launch_dep100.sh.
STAGE=align exec /bin/bash /home/freiburger/Documents/prFBA/deployment_d100/run_codif_d100.sh
