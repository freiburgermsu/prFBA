#!/bin/bash
# SUPERSEDED (2026-10-29, after the pre-launch review in WEIGHTED_OBJECTIVE_TEST.md section 12):
# the full run is launched by scripts/launch_wobj.sh (unit wobj-run, log logs/wobj_run.log,
# MemoryMax 20G, --maxtasksperchild 2).  This stub refuses to run so that the old settings
# (MemoryMax 24G, --maxtasksperchild 10, log logs/weighted_objective.log) cannot be used by mistake.
echo "superseded: use scripts/launch_wobj.sh" >&2
exit 1
