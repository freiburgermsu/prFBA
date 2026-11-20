#!/bin/bash
# R2.2 on FM: pull TF's MAG outputs, move FM's stale git-ignored batch tables (only copy) to $W/pre_repair,
# back up the driver logs that the pull overwrites, delete the stale tracked batch dirs (old versions in git),
# install the pulled tree, then verify FM == TF (sha256 + size + mtime) for every installed file.
set -euo pipefail
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; T=$W/tf
F=$W/from_tf/pull_R2.2_mag
[ -e $F ] && { echo "$F exists"; exit 2; }
bash $T/bin/tf_pull.sh $T/lists/pull_R2.2_mag.txt $F
for d in mag_validation/out_batches mag_validation/out_batches_tau99 mag_validation/out_batches_n500 mag_validation/lowc_out_batches mag_validation_km/out_batches mag_validation_rnd/out_batches; do
  n=$(ls -d $F/$d/b* | wc -l); c=$(ls $F/$d/b*/mag_cells.csv | wc -l); s=$(ls $F/$d/b*/mag_summary.json | wc -l)
  echo "pulled $d: $n batches, $c cells, $s summaries"; [ $n = 20 ] && [ $c = 20 ] && [ $s = 20 ] || { echo "incomplete $d"; exit 3; }
done
ls $F/mag_validation/smoke/out/mag_cells.csv $F/mag_validation/smoke/out/mag_summary.json >/dev/null
cd $R
for d in mag_validation/out_batches mag_validation/lowc_out_batches mag_validation_rnd/out_batches; do
  [ -e $W/pre_repair/$d ] && { echo "pre_repair/$d exists"; exit 4; }
  mkdir -p $W/pre_repair/$(dirname $d); mv $d $W/pre_repair/$d; echo "moved $d -> pre_repair"
done
for d in mag_validation/logs mag_validation_km/logs mag_validation_rnd/logs; do
  mkdir -p $W/pre_repair/$d; n=0
  for f in $(cd $F/$d && ls); do if [ -e $d/$f ]; then cp -p $d/$f $W/pre_repair/$d/$f; n=$((n+1)); fi; done
  echo "backed up $n old $d files that the pull overwrites"
done
rm -rf mag_validation/out_batches_tau99 mag_validation/out_batches_n500 mag_validation/smoke/out mag_validation_km/out_batches
rsync -a $F/ $R/
echo "installed $(cd $F && find . -type f | wc -l) files"
bash $T/bin/tf_verify.sh $T/lists/pull_R2.2_mag_nologs.txt
(cd $F && find mag_validation/logs mag_validation_km/logs mag_validation_rnd/logs -type f | sort) > $W/tf/lists/pulled_R2.2_logs_files.txt
bash $T/bin/tf_verify.sh $W/tf/lists/pulled_R2.2_logs_files.txt
