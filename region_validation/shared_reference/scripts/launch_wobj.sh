#!/bin/bash
# Full run of the weighted-objective test (WEIGHTED_OBJECTIVE_TEST.md, section 14).
# WRITTEN, NOT LAUNCHED, by its author: the orchestrator launches it on a memory-aware schedule.
#     bash scripts/launch_wobj.sh
#
# unit        wobj-run (transient systemd user unit, --collect)
# population  360 amplicons (n2=100, n3=100, n4+=100, n1=60), the harness's default
# workers     8 (one Gurobi thread each), --maxtasksperchild 2
# duration    about 1.8-2.5 h: 360 x 131.6 s / 8 = 1.65 h at the second smoke's load, x 1.42
#             under the heavier load the first smoke met; + worker restarts (180 x <= 20 s / 8)
#             + the population build (about 2 min)
# memory      estimated peak 14.7 GiB = 8 workers x 1.75 GiB + parent 0.74 GiB, where
#             1.75 GiB = the largest single-task worker peak (1.604 GiB) + the measured growth
#             of a second task in the same process (0.146 GiB); population stage alone 4.2 GiB
#             (runs first, before the pool).  MemoryMax = ceil(14.74 x 1.3) = 20G.
#             Per-process growth over 4 tasks was 1.542 -> 1.688 -> 1.906 -> 2.028 GiB, hence
#             --maxtasksperchild 2 (the harness default is 10).
set -euo pipefail
UNIT=wobj-run
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
OUT=$SR/data/weighted_objective
LOG=$SR/logs/wobj_run.log
PIN=/home/freiburger/Documents/ModelSEEDpy-c2e31a0
PY="/home/freiburger/Documents/py_venv/bin/python -W ignore"
SCRIPT=$SR/scripts/weighted_objective_test.py
MEMMAX=20G
NEED_GB=35      # MemAvailable at launch: MemoryMax 20G + 12 GiB (GTDB-Tk pauses below that) + margin

free -g
avail_gb=$(awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo)
if [ "$avail_gb" -lt "$NEED_GB" ]; then
    echo "MemAvailable ${avail_gb} GiB < ${NEED_GB} GiB: not launching" >&2
    exit 1
fi
if systemctl --user is-active --quiet "$UNIT"; then
    echo "$UNIT is already active: not launching" >&2
    exit 1
fi
if [ -e "$OUT" ]; then
    echo "$OUT exists: this launch writes a NEW output directory; refusing" >&2
    echo "(to resume a stopped run, run stage 2 below by hand with a new log name)" >&2
    exit 1
fi
if [ -e "$LOG" ]; then
    echo "$LOG exists: refusing to overwrite a log" >&2
    exit 1
fi
mkdir -p "$SR/logs"

# Stage 1 builds $OUT/population.json alone (peak 4.2 GiB, about 2 min), so the pool's parent
# does not carry the five predictions.json loads into the run.  Stage 2 runs every amplicon
# (resumable: finished keys in per_amplicon.jsonl are skipped) and writes summary.json.
systemd-run --user --unit="$UNIT" --collect -p MemoryMax="$MEMMAX" -p WorkingDirectory="$SR" \
    -E PYTHONPATH="$PIN" /bin/bash -c "\
$PY $SCRIPT --outdir $OUT --expect-modelseedpy c2e31a0 --population-only > $LOG 2>&1 && \
$PY $SCRIPT --outdir $OUT --expect-modelseedpy c2e31a0 --workers 8 --maxtasksperchild 2 >> $LOG 2>&1; \
echo EXIT=\$? >> $LOG"
echo "launched $UNIT; log $LOG; outputs $OUT/{population.json,per_amplicon.jsonl,summary.json}"
