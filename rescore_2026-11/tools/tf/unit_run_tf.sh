#!/bin/bash
# Runs INSIDE prfba-ns inside a tf_unit (TF). Like FM's tools/unit_run.sh: executes <cmdfile> with bash,
# samples the unit's cgroup memory every 30 s into $W/logs/<name>.mem.log, and appends one JSON line to
# $W/stamps/units_tf.jsonl with the cgroup memory.peak read from inside the cgroup, GNU time's max RSS,
# wall and rc.   usage (as the tf_unit COMMAND): bash $W/tools/unit_run_tf.sh <name> <cmdfile>
set -u
name=$1; cmdfile=$2
W=/home/freiburger/Documents/prFBA_rescore_2026-11
cg=/sys/fs/cgroup$(cut -d: -f3 /proc/self/cgroup)
memlog=$W/logs/$name.mem.log; timef=$W/logs/$name.time.txt
cp "$cmdfile" "$W/logs/$name.cmd.sh"
sample() { echo "$(date -Is) cur=$(cat "$cg/memory.current") peak=$(cat "$cg/memory.peak") avail_kB=$(awk '/MemAvailable/{print $2}' /proc/meminfo)"; }
( while :; do sample; sleep 30; done ) >> "$memlog" 2>/dev/null &
poller=$!
t0=$(date +%s); start=$(date -Is)
echo "[unit $name] start $start cgroup $cg cmdfile $cmdfile host $(hostname)"
/usr/bin/time -v -o "$timef" /bin/bash "$W/logs/$name.cmd.sh"
rc=$?
sample >> "$memlog"; kill "$poller" 2>/dev/null
peak=$(cat "$cg/memory.peak"); maxrss=$(awk -F': ' '/Maximum resident set size/{print $2}' "$timef" 2>/dev/null)
wall=$(( $(date +%s) - t0 ))
printf '{"unit": "%s", "host": "TF", "rc": %d, "start": "%s", "end": "%s", "wall_s": %d, "cgroup_peak_bytes": %s, "cgroup_peak_GB": %.2f, "maxrss_kB": %s}\n' \
  "$name" "$rc" "$start" "$(date -Is)" "$wall" "$peak" "$(echo "$peak/1073741824" | bc -l)" "${maxrss:-null}" >> $W/stamps/units_tf.jsonl
echo "[unit $name] end rc=$rc wall=${wall}s cgroup_peak=$peak maxrss_kB=${maxrss:-?}"
exit $rc
