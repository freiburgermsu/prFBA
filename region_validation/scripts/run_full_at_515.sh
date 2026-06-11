#!/usr/bin/env bash
# Timed launcher: sleep until 5:15pm CT (machine clock is CT), then run the FULL
# 10k validation fresh. Pilot outputs (500 genomes) are archived first so the full
# run regenerates every stage at N=10000 (the shared PGFam cache in bvbrc_cache/
# persists, so prefetch resumes incrementally). Runs detached in tmux -> survives
# the controlling session closing.
set -uo pipefail
RV="/home/freiburger/Documents/prFBA/region_validation"

TARGET=$(date -d "today 17:15" +%s)
NOW=$(date +%s)
WAIT=$(( TARGET > NOW ? TARGET - NOW : 0 ))
echo "[timer] now=$(date)   target=$(date -d "@${TARGET}")   sleeping ${WAIT}s ..."
sleep "${WAIT}"
echo "[timer] firing FULL run at $(date)"

# Archive the pilot so the full run starts clean at N=10000 (idempotent).
[ -d "${RV}/data" ]    && [ ! -e "${RV}/data_pilot_500" ]    && mv "${RV}/data"    "${RV}/data_pilot_500"
[ -d "${RV}/figures" ] && [ ! -e "${RV}/figures_pilot_500" ] && mv "${RV}/figures" "${RV}/figures_pilot_500"
mkdir -p "${RV}/data" "${RV}/figures"

cd "${RV}/scripts"
bash run_validation.sh full 2>&1 | tee "${RV}/data/full_run.log"
echo "FULL_RUN_DONE=${PIPESTATUS[0]} at $(date)" | tee -a "${RV}/data/full_run.log"
