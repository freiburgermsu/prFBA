#!/usr/bin/env bash
# The plan's U() as a detached transient user unit with a working memory record (unit_run.sh).
#   usage: U.sh <unit> <MemoryMax> <cmdfile> <logfile>
# Returns as soon as the unit starts; completion = a line for <unit> in $W/stamps/units.jsonl
# (and "[unit <unit>] end" in <logfile>).
set -euo pipefail
unit=$1; mem=$2; cmdfile=$3; logfile=$4
W=/home/freiburger/Documents/prFBA_rescore_2026-11
R=/home/freiburger/Documents/prFBA
[ -s "$cmdfile" ] || { echo "U.sh: no cmdfile $cmdfile" >&2; exit 2; }
if systemctl --user is-active --quiet "$unit.service" 2>/dev/null; then
  echo "U.sh: unit $unit is already active" >&2; exit 3
fi
cp "$cmdfile" "$W/logs/$unit.cmd.sh"          # the exact command, kept with the logs
systemd-run --user --unit="$unit" --collect -p MemoryAccounting=yes -p MemoryMax="$mem" \
  -p MemorySwapMax=0 -p LimitNOFILE=1048576 -p WorkingDirectory=$R \
  -E PYTHONHASHSEED=0 -E PYTHONUNBUFFERED=1 -E SOURCE_DATE_EPOCH=1794700800 \
  /bin/bash "$W/tools/unit_run.sh" "$unit" "$W/logs/$unit.cmd.sh" "$logfile"
