#!/usr/bin/env bash
# R3.3: re-run run_picrust2_fold.py (ledger / run_c2 argv, NO --force) for every run whose reference traits
# changed in R3.2: genus 23 forced + 5 nsti, family 25 forced + 5 nsti, C2 10. The runner's stage stamps
# (traits_sha256) force re-placement and hsp; nothing else is touched. Each run dir is first copied
# (minus intermediate/) to pre_repair. 4 concurrent.
set -euo pipefail
LIST=${R33_LIST:-/home/freiburger/Documents/prFBA_rescore_2026-11/control/r3/r33_runs.tsv}; TAG=${R33_TAG:-all}; NJ=${R33_JOBS:-4}
W=/home/freiburger/Documents/prFBA_rescore_2026-11; R=/home/freiburger/Documents/prFBA
SR=$R/region_validation/shared_reference; C=$W/control/r3; B=$W/pre_repair/region_validation/shared_reference
cd $SR
echo "== before stamps $(date -Is)"
find data/picrust2 data_family/picrust2 data_c2/picrust2 -type f -not -path '*/intermediate/*' -print0 | sort -z | xargs -0 sha256sum > $C/r33_p2_before_$TAG.sha256
wc -l $C/r33_p2_before_$TAG.sha256
echo "== copy aside $(date -Is)"
while IFS=$'\t' read -r root name cmd; do
  d=$SR/$root/picrust2/$name; b=$B/$root/picrust2/$name
  [ -e $b/predictions.json ] && { echo "backup exists: $b" >&2; exit 2; }
  mkdir -p $b
  rsync -a --exclude='intermediate/' $d/ $b/
done < $LIST
echo "== runs (4 concurrent) $(date -Is)"
set +e
xargs -P $NJ -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r root name cmd <<< "{}"; eval bash '$W'/tools/r33_one.sh $root $name $cmd' < $LIST
rc=$?
set -e
echo "== xargs rc=$rc $(date -Is)"
find data/picrust2 data_family/picrust2 data_c2/picrust2 -type f -not -path '*/intermediate/*' -print0 | sort -z | xargs -0 sha256sum > $C/r33_p2_after_$TAG.sha256
grep -h "^EXIT=" $W/logs/r3/p2_*.log | cut -d" " -f1 | sort | uniq -c
exit $rc
