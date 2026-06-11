#!/usr/bin/env bash
# Rebuild the EmilyKin per-ASV synthetic genomes from the GPU top-20 hits (derived from the
# existing GPU top-500), replacing the edlib-based build. Uses the prFBA bvbrc_cache.
set -o pipefail
cd /home/freiburger/Documents/prFBA || exit 2
~/Documents/py_venv/bin/python build_synthetic_genomes_parallel.py \
  --hits          /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_top20_alignment_hits_gpu.json \
  --out-dir       /home/freiburger/Documents/EmilyKin/synthetic_genomes \
  --selection-out /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_reference_selection_gpu.json \
  --procs 48 --io-workers 32
echo "EMILYKIN_GPU_REBUILD_DONE exit=$?"
