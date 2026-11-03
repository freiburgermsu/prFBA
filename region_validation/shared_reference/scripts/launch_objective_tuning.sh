#!/bin/bash
# Full run of Stage O (PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md, section "Stage O"), to be
# launched AFTER Stage E has written data/pfunc_estimators/stage_o_weights.json.
# WRITTEN, NOT LAUNCHED, by its author:
#     bash scripts/launch_objective_tuning.sh
#
# unit        otune-run (transient systemd user unit, --collect)
# population  the weighted-objective test's 360 amplicons (data/weighted_objective/population.json),
#             90 arms each (5 P x 2 C x 3 G x 3 X), fresh model copy per arm
# workers     10 (one Gurobi thread each), --maxtasksperchild 1: every cobra Model.copy() leaks
#             its optlang solver objects upstream (~6-7 MiB per copy of a 1,500-reaction model;
#             ~92 copies per amplicon), so a worker never carries one amplicon's leak into the next
# duration    about 1.3-1.9 h: per-amplicon work = the prior run's own gapfill + truth time per
#             row (mean 69 s) + arms 2.9 s + 0.0297 s x n_rxn (fit on the 6-amplicon smoke; mean
#             n_rxn 1,339) = 112.7 s mean, + ~17 s worker start per amplicon; 360 x 129.7 s / 10
#             = 1.30 h, x 1.42 under the heavier load the weighted-objective test's first smoke met
# memory      estimated peak 25.5 GiB = 10 workers x 2.4 GiB + parent 1.5 GiB, where 2.4 GiB covers
#             the largest fresh-worker peak seen (2.23 GiB, 1,721 reactions, with the smoke's extra
#             repro and FVA copies; 1.96 GiB at 1,549 reactions without them) extrapolated to the
#             sample's largest model (1,980 reactions); parent 1.04 GiB on the smoke.  The smoke's
#             cgroup peak was 6.90 GiB for 3 workers running two tasks each (second-task peaks up
#             to 2.79 GiB, which --maxtasksperchild 1 removes).  MemoryMax = ceil(1.3 x 25.5) = 34G.
# outputs     data/objective_tuning/{per_amplicon.jsonl, meta.json}; log logs/otune_run.log
# afterwards  ~/Documents/py_venv/bin/python scripts/objective_tuning_results.py
#             (no solver; writes data/objective_tuning/summary.json and
#             data/objective_tuning/OBJECTIVE_TUNING_RESULTS.filled.md)
set -euo pipefail
UNIT=otune-run
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
OUT=$SR/data/objective_tuning
WEIGHTS=${OTUNE_WEIGHTS:-$SR/data/pfunc_estimators/stage_o_weights.json}   # Stage O carries E3: OTUNE_WEIGHTS=.../stage_o_weights_E3.json
LOG=$SR/logs/otune_run.log
PIN=/home/freiburger/Documents/ModelSEEDpy-c2e31a0
PY="/home/freiburger/Documents/py_venv/bin/python -W ignore"
SCRIPT=$SR/scripts/objective_tuning.py
MEMMAX=34G
NEED_GB=50      # MemAvailable at launch: MemoryMax 34G + 16 GiB for the machine's other jobs

free -g
if [ ! -s "$WEIGHTS" ]; then
    echo "$WEIGHTS is missing: Stage E has not written its weights; not launching" >&2
    exit 1
fi
avail_gb=$(awk '/MemAvailable/ {printf "%d", $2/1048576}' /proc/meminfo)
if [ "$avail_gb" -lt "$NEED_GB" ]; then
    echo "MemAvailable ${avail_gb} GiB < ${NEED_GB} GiB: not launching" >&2
    exit 1
fi
if systemctl --user is-active --quiet "$UNIT"; then
    echo "$UNIT is already active: not launching" >&2
    exit 1
fi
for f in "$OUT/per_amplicon.jsonl" "$OUT/meta.json"; do
    if [ -e "$f" ]; then
        echo "$f exists: this launch starts a NEW run; refusing" >&2
        echo "(a stopped run resumes by hand: the same command with a new log name; finished" >&2
        echo " keys are skipped and meta.json's spec must match)" >&2
        exit 1
    fi
done
if [ -e "$LOG" ]; then
    echo "$LOG exists: refusing to overwrite a log" >&2
    exit 1
fi
mkdir -p "$SR/logs"

# Pre-flight (parent checks only, nothing written): population, prior rows, g_meas, the weights
# file (not a stub; every key; md5 / fold / reps; no source leakage; isotonic map monotone) and
# every neighbour and clade member in its fold's EC cache.
cd "$SR"
PYTHONPATH="$PIN" $PY "$SCRIPT" --outdir "$OUT" --weights "$WEIGHTS" \
    --expect-modelseedpy c2e31a0 --validate-only

# Failure mode to watch: if the cgroup OOM killer kills a pool worker, multiprocessing.Pool
# replaces it but the task is lost and imap_unordered waits forever; the unit then stays active
# while the log gets no new "[otune] <key>" rows.  The longest amplicon should take ~3-4 min, so
# no new row for ~15 min is a stall.  Finished keys are skipped on a restart.
systemd-run --user --unit="$UNIT" --collect -p MemoryMax="$MEMMAX" -p MemorySwapMax=0 \
    -p WorkingDirectory="$SR" -E PYTHONPATH="$PIN" /bin/bash -c "\
$PY $SCRIPT --outdir $OUT --weights $WEIGHTS --expect-modelseedpy c2e31a0 --workers 10 \
--maxtasksperchild 1 > $LOG 2>&1; echo EXIT=\$? >> $LOG"
echo "launched $UNIT; log $LOG; outputs $OUT/{per_amplicon.jsonl,meta.json}"
echo "watch: systemctl --user is-active $UNIT; grep -c '^\[otune\] .* ok ' $LOG"
