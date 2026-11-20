#!/usr/bin/env bash
# R3.2: refresh traits + manifest of the affected references (genus 23, family 25, C2 2) by re-running
# build_reference.py with the ledger argv WITHOUT --force (reuses .fna/.hmm/.tre/.probe.tre; redoes the
# crash screen, traits and manifest). Old manifest.json + traits_EC.txt.gz are MOVED to pre_repair; the
# other files the rebuild rewrites are COPIED there first.
set -euo pipefail
W=/home/freiburger/Documents/prFBA_rescore_2026-11; R=/home/freiburger/Documents/prFBA
SR=$R/region_validation/shared_reference; C=$W/control/r3; B=$W/pre_repair/region_validation/shared_reference
cd $SR
echo "== before stamps $(date -Is)"
find data/refs data_family/refs data_c2/refs -type f -not -path '*/probe_int/*' -print0 | sort -z | xargs -0 sha256sum > $C/r32_refs_before.sha256
wc -l $C/r32_refs_before.sha256
echo "== move/copy aside $(date -Is)"
while IFS=$'\t' read -r root name cmd; do
  d=$SR/$root/refs/$name; b=$B/$root/refs/$name
  [ -e $b/manifest.json ] && { echo "backup exists: $b" >&2; exit 2; }
  mkdir -p $b
  for f in traits_16S.txt.gz $name.model probe_queries.fna $name.input.fna $name.hsp_screen.log; do cp -p $d/$f $b/; done
  cp -rp $d/hsp_screen $b/
  mv $d/manifest.json $d/traits_EC.txt.gz $b/
done < $C/r32_builds.tsv
echo "== builds (12 concurrent) $(date -Is)"
set +e
awk -F'\t' '{print $1"\t"$2"\t"$3}' $C/r32_builds.tsv | \
  xargs -P 12 -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r root name cmd <<< "{}"; eval bash '$W'/tools/r32_one.sh $root $name $cmd'
rc=$?
set -e
echo "== xargs rc=$rc $(date -Is)"
echo "== after stamps $(date -Is)"
find data/refs data_family/refs data_c2/refs -type f -not -path '*/probe_int/*' -print0 | sort -z | xargs -0 sha256sum > $C/r32_refs_after.sha256
grep -h '^EXIT=' $W/logs/r3/build_*.log | sort | uniq -c
exit $rc
