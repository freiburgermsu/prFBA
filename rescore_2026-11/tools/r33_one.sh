#!/usr/bin/env bash
# one R3.3 PICRUSt2 run: r33_one.sh <root> <name> <command...>  (log in $W/logs/r3/p2_<root>_<name>.log)
root=$1; name=$2; shift 2
W=/home/freiburger/Documents/prFBA_rescore_2026-11
log=$W/logs/r3/p2_${root}_${name}.log
t0=$(date +%s)
{ echo "# $(date -Is) $root/$name"; echo "# $*"; } > $log
"$@" >> $log 2>&1; rc=$?
echo "EXIT=$rc wall=$(( $(date +%s) - t0 ))s" >> $log
echo "$(date -Is) $root/$name rc=$rc wall=$(( $(date +%s) - t0 ))s"
exit $rc
