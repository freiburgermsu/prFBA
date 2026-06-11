#!/usr/bin/env bash
# Full per-ASV synthetic-genome build (exact selection -> warm genome cache -> parallel merge).
# Reproduces the codiffusion_bioreactor/genome_objects data, keyed per ASV. Launched in tmux.
set -o pipefail
cd /home/freiburger/Documents/prFBA || exit 2
~/Documents/py_venv/bin/python build_synthetic_genomes_parallel.py \
  --hits    /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_top20_alignment_hits.json \
  --out-dir /home/freiburger/Documents/EmilyKin/synthetic_genomes \
  --procs 56 --io-workers 32
echo "RUN_DONE exit=$?"
