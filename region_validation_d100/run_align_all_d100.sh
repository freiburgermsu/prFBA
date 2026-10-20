#!/usr/bin/env bash
# Depth-100 alignment of the three union trees, one after another (GPU + CPU exclusive).
# Idempotent/resumable: a panel whose data/hits/raw/align.done exists is skipped; a panel
# that died after the edlib stage resumes from its saved shortlist (run_validation.sh
# passes --shortlist-in automatically). Never touches the shared PGFam cache.
#   ALIGN_WORKERS (default 40) is passed through to align_hits.py --workers.
set -euo pipefail
R=/home/freiburger/Documents/prFBA
export ALIGN_WORKERS="${ALIGN_WORKERS:-40}"
for T in region_validation_d100 region_validation_kmeans_d100 region_validation_random_d100; do
  if [ -s "$R/$T/data/hits/raw/align.done" ]; then
    echo "[align-all] skip $T (align.done present)"; continue
  fi
  echo "[align-all] $T start $(date -Is) workers=${ALIGN_WORKERS}"
  STAGE=align NO_TMUX=1 bash "$R/$T/scripts/run_validation.sh" full
  echo "[align-all] $T done $(date -Is)"
done
echo "[align-all] all three panels aligned $(date -Is)"
