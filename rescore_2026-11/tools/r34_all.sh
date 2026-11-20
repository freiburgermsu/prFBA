#!/usr/bin/env bash
# R3.4 on FM: the 45 prFBA re-selects (genus 20, family 20, C2 5), ledger / run_c2 argv without --force,
# via r34_one.sh, N concurrent (single-threaded selection; ~10 GB RSS each at d100).
set -uo pipefail
W=/home/freiburger/Documents/prFBA_rescore_2026-11; C=$W/control/r3
N=${R34_JOBS:-5}; LIST=${R34_LIST:-$C/r34_runs.tsv}
echo "== R3.4 $(date -Is) jobs=$N list=$LIST ($(wc -l < $LIST) runs)"
xargs -P $N -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r root name cmd <<< "{}"; eval bash '$W'/tools/r34_one.sh $root $name $cmd' < $LIST
rc=$?
echo "== xargs rc=$rc $(date -Is)"
grep -h '^EXIT=' $W/logs/r3/prfba_*.log | cut -d' ' -f1 | sort | uniq -c
exit $rc
