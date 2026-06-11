#!/usr/bin/env bash
# GPU-exhaustive Smith-Waterman alignment of the 3,276 new_data (codif_all) ASVs against
# the BV-BRC 16S DB -> asv_top20_alignment_hits_gpu.json (the GPU-only hits for the build).
set -o pipefail
cd /home/freiburger/Documents/prFBA || exit 2
~/Documents/py_venv/bin/python gpu_align.py --full --topk 20 --workers 56 \
  --fasta        /home/freiburger/Documents/codiffusion_bioreactor/new_data/dna-sequences_codif_all.fasta \
  --taxonomy-csv /home/freiburger/Documents/codiffusion_bioreactor/new_data/taxonomy_codif_all.csv \
  --outdir       /home/freiburger/Documents/codiffusion_bioreactor/bvbrc_alignment_hits
echo "GPU_ALIGN_DONE exit=$?"
