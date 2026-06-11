#!/usr/bin/env bash
# Build per-ASV synthetic genomes for the new_data ASVs from the GPU hits ONLY.
# new_data is a different organism set than EmilyKin, so prefetch its candidate gene-families
# first (threaded) to avoid the lazy per-miss cache rewrite, then warm + parallel merge.
set -o pipefail
cd /home/freiburger/Documents/prFBA || exit 2
HITS=/home/freiburger/Documents/codiffusion_bioreactor/bvbrc_alignment_hits/asv_top20_alignment_hits_gpu.json
echo "=== [1/2] prefetch gene-families for new_data stage-7 candidates ==="
~/Documents/py_venv/bin/python prefetch_gene_families.py --hits "$HITS" \
  --cache bvbrc_cache/genome_gene_families.json --workers 32
echo "=== [2/2] build synthetic genomes (GPU hits only) -> codiffusion/synthetic_genomes ==="
~/Documents/py_venv/bin/python build_synthetic_genomes_parallel.py \
  --hits          "$HITS" \
  --out-dir       /home/freiburger/Documents/codiffusion_bioreactor/synthetic_genomes \
  --selection-out /home/freiburger/Documents/codiffusion_bioreactor/bvbrc_alignment_hits/asv_reference_selection_gpu.json \
  --procs 48 --io-workers 32
echo "NEWDATA_BUILD_DONE exit=$?"
