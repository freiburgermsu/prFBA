#!/usr/bin/env bash
# Runs INSIDE a transient user unit (launched by U.sh). Executes <cmdfile> with bash, samples the
# unit's own cgroup every 30 s (memory.current / memory.peak / MemAvailable) into
# $W/logs/<name>.mem.log, and appends one JSON line to $W/stamps/units.jsonl at the end with the
# cgroup memory.peak (read from inside the cgroup before exit: systemd 255's end-of-unit
# "Memory peak" report reads an emptied cgroup and prints 256K), GNU time's max RSS, wall and rc.
#   usage: unit_run.sh <name> <cmdfile> <logfile>
set -u
name=$1; cmdfile=$2; logfile=$3
W=/home/freiburger/Documents/prFBA_rescore_2026-11
cg=/sys/fs/cgroup$(cut -d: -f3 /proc/self/cgroup)
memlog=$W/logs/$name.mem.log
timef=$W/logs/$name.time.txt
sample() { echo "$(date -Is) cur=$(cat "$cg/memory.current") peak=$(cat "$cg/memory.peak") avail_kB=$(awk '/MemAvailable/{print $2}' /proc/meminfo)"; }
( while :; do sample; sleep 30; done ) >> "$memlog" 2>/dev/null &
poller=$!
t0=$(date +%s); start=$(date -Is)
echo "[unit $name] start $start cgroup $cg cmdfile $cmdfile" >> "$logfile"
/usr/bin/time -v -o "$timef" /bin/bash "$cmdfile" >> "$logfile" 2>&1
rc=$?
sample >> "$memlog"
kill "$poller" 2>/dev/null
peak=$(cat "$cg/memory.peak")
maxrss=$(awk -F': ' '/Maximum resident set size/{print $2}' "$timef" 2>/dev/null)
wall=$(( $(date +%s) - t0 ))
printf '{"unit": "%s", "rc": %d, "start": "%s", "end": "%s", "wall_s": %d, "cgroup_peak_bytes": %s, "cgroup_peak_GB": %.2f, "maxrss_kB": %s, "cmdfile": "%s", "log": "%s"}\n' \
  "$name" "$rc" "$start" "$(date -Is)" "$wall" "$peak" "$(echo "$peak/1073741824" | bc -l)" "${maxrss:-null}" "$cmdfile" "$logfile" >> $W/stamps/units.jsonl
echo "[unit $name] end rc=$rc wall=${wall}s cgroup_peak=$peak maxrss_kB=${maxrss:-?}" >> "$logfile"
exit $rc
