#!/bin/bash
# R2.2 MAG control on TF (inside prfba-ns): score every MAG arm with the W1 --restore-old-values control
# indices (= the pre-repair indices, sha dcde696b / b597a00f / 64c2ea64) and today's code, and compare
# each batch's mag_cells.csv + mag_summary.json byte-for-byte with the committed / pre_repair old tables.
# main + tau99 are also concatenated (header of b000, then tail -n+2 of each batch in name order) and must
# give md5 a57df30b8346fcb72f14f7f33708669a / 3cb0e487aa31c2f30f1699fcbdf88c92 (plan R2.2 step 2).
set -uo pipefail
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
CI=$W/control/panel_index; O=$W/control/mag; PRE=$W/pre_repair; PY=/home/freiburger/Documents/py_venv/bin/python
cd $R; mkdir -p $O/logs
one() {  # arm batch
  local arm=$1 B=$2 a
  case $arm in
    main)  a="--hits mag_validation/hits/$B.tsv --manifest mag_validation/manifests/$B.json --pgfams $CI/pgfam_panel.pre_repair.json";;
    tau99) a="--hits mag_validation/hits/$B.tsv --manifest mag_validation/manifests/$B.json --pgfams $CI/pgfam_panel.pre_repair.json --tau 99";;
    n500)  a="--hits mag_validation/hits_n500/$B.tsv --manifest mag_validation/manifests_n500/$B.json --pgfams $CI/pgfam_panel.pre_repair.json";;
    lowc)  a="--hits mag_validation/lowc_hits/$B.tsv --manifest mag_validation/lowc_manifests/$B.json --pgfams $CI/pgfam_panel.pre_repair.json";;
    km)    a="--hits mag_validation_km/hits/$B.tsv --manifest mag_validation_km/manifests/$B.json --truth $R/region_validation_kmeans/data/benchmark_truth.json --pgfams $CI/pgfam_panel_kmeans.pre_repair.json";;
    rnd)   a="--hits mag_validation_rnd/hits/$B.tsv --manifest mag_validation_rnd/manifests/$B.json --truth $R/region_validation_random/data/benchmark_truth.json --pgfams $CI/pgfam_panel_random.pre_repair.json";;
    smoke) a="--hits mag_validation/smoke/hits.tsv --manifest mag_validation/smoke/manifest.json --pgfams $CI/pgfam_panel.pre_repair.json";;
  esac
  $PY $R/mag_validation/score_mag.py $a --outdir $O/$arm/$B > $O/logs/${arm}_$B.log 2>&1 && echo "[$arm $B] ok" || echo "[$arm $B] FAILED"
}
export -f one; export R W CI O PY
{ for arm in main tau99 n500 lowc km rnd; do
    case $arm in main|tau99) d=mag_validation/hits;; n500) d=mag_validation/hits_n500;; lowc) d=mag_validation/lowc_hits;;
                 km) d=mag_validation_km/hits;; rnd) d=mag_validation_rnd/hits;; esac
    for f in $(ls -1 $d/*.tsv | sort); do echo "$arm $(basename $f .tsv)"; done
  done; echo "smoke out"; } > $O/jobs.txt
echo "jobs: $(wc -l < $O/jobs.txt)"; t0=$(date +%s)
xargs -P ${P:-6} -L 1 bash -c 'one "$@"' _ < $O/jobs.txt | sort | uniq -c | awk '{print $NF}' | sort | uniq -c
echo "scoring wall $(( $(date +%s) - t0 )) s"
ref() { case $1 in main) echo $PRE/mag_validation/out_batches;; tau99) echo $R/mag_validation/out_batches_tau99;;
          n500) echo $R/mag_validation/out_batches_n500;; lowc) echo $PRE/mag_validation/lowc_out_batches;;
          km) echo $R/mag_validation_km/out_batches;; rnd) echo $PRE/mag_validation_rnd/out_batches;; smoke) echo $R/mag_validation/smoke;; esac; }
fail=0
while read -r arm B; do
  r=$(ref $arm)/$B; for f in mag_cells.csv mag_summary.json; do
    if ! cmp -s $O/$arm/$B/$f $r/$f; then echo "DIFFER $arm $B $f"; fail=1; fi; done
done < $O/jobs.txt
for arm in main tau99 n500 lowc km rnd smoke; do
  printf "%-6s batches=%s identical=%s\n" $arm "$(grep -c "^$arm " $O/jobs.txt)" "$(while read -r a B; do [ $a = $arm ] && cmp -s $O/$a/$B/mag_cells.csv $(ref $a)/$B/mag_cells.csv && cmp -s $O/$a/$B/mag_summary.json $(ref $a)/$B/mag_summary.json && echo x; done < $O/jobs.txt | wc -l)"
done
cat_cells() { local out=$1 d=$2; head -1 "$d/b000/mag_cells.csv" > $out; for b in $(ls -d $d/b* | sort); do tail -n+2 $b/mag_cells.csv >> $out; done; }
cat_cells $O/mag_cells_all.csv $O/main; cat_cells $O/mag_cells_tau99.csv $O/tau99
m1=$(md5sum < $O/mag_cells_all.csv | cut -c1-32); m2=$(md5sum < $O/mag_cells_tau99.csv | cut -c1-32)
echo "control mag_cells_all.csv   md5 $m1 rows $(($(wc -l < $O/mag_cells_all.csv)-1)) $([ $m1 = a57df30b8346fcb72f14f7f33708669a ] && echo MATCH || { echo MISMATCH; fail=1; })"
echo "control mag_cells_tau99.csv md5 $m2 rows $(($(wc -l < $O/mag_cells_tau99.csv)-1)) $([ $m2 = 3cb0e487aa31c2f30f1699fcbdf88c92 ] && echo MATCH || { echo MISMATCH; fail=1; })"
[ $fail = 0 ] && echo "MAG CONTROL PASS" || echo "MAG CONTROL FAIL"
exit $fail
