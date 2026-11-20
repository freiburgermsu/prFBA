#!/bin/bash
# R2.2 steps 3-5 on TF (inside prfba-ns): clear the stale batch tables, re-score every MAG arm with the
# drivers as committed (they read the canonical, now W1-rebuilt, indices), score the smoke set, then
# concatenate and run aggregate.py x4 (plan R2.2 step 5).
set -uo pipefail
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; PY=/home/freiburger/Documents/py_venv/bin/python
cd $R
echo "== indices"; sha256sum mag_validation/pgfam_panel.json kmeans_panel/pgfam_panel_kmeans.json random_panel/pgfam_panel_random.json
want="cd5c3f49b13f 22bb5115731a 7986e1087528"
got=$(sha256sum mag_validation/pgfam_panel.json kmeans_panel/pgfam_panel_kmeans.json random_panel/pgfam_panel_random.json | cut -c1-12 | tr '\n' ' ' | sed 's/ $//')
[ "$got" = "$want" ] || { echo "indices are not the W1 rebuilds: $got"; exit 2; }
echo "== clear stale batch tables (tracked: old versions in git; git-ignored: only in \$W/pre_repair)"
rm -rf mag_validation/out_batches_tau99 mag_validation/out_batches_n500 mag_validation/smoke/out mag_validation_km/out_batches
for d in mag_validation/out_batches mag_validation/lowc_out_batches; do [ -e $d ] && { echo "unexpected $d"; exit 2; }; done
[ -z "$(ls -A mag_validation_rnd/out_batches 2>/dev/null)" ] || { echo "mag_validation_rnd/out_batches not empty"; exit 2; }
mkdir -p mag_validation/logs mag_validation_km/logs mag_validation_rnd/logs
t0=$(date +%s); T=$t0
for s in "mag_validation/score_all.sh" "mag_validation/score_tau.sh 99" "mag_validation/score_n500.sh" "mag_validation/score_lowc.sh" "mag_validation_km/score_km.sh" "mag_validation_rnd/score_rnd.sh"; do
  bash $s > $W/logs/r22_driver_$(basename ${s%% *} .sh)$(echo $s | awk '{print $2}').log 2>&1
  echo "$s rc=$? $(( $(date +%s) - T )) s; ok=$(grep -c '\] ok' $W/logs/r22_driver_$(basename ${s%% *} .sh)$(echo $s | awk '{print $2}').log) failed=$(grep -c FAILED $W/logs/r22_driver_$(basename ${s%% *} .sh)$(echo $s | awk '{print $2}').log)"
  T=$(date +%s)
done
$PY mag_validation/score_mag.py --hits mag_validation/smoke/hits.tsv --manifest mag_validation/smoke/manifest.json --outdir mag_validation/smoke/out > $W/logs/r22_smoke.log 2>&1
echo "smoke rc=$? $(( $(date +%s) - T )) s; scoring total $(( $(date +%s) - t0 )) s"
nf=$(grep -l FAILED $W/logs/r22_driver_*.log | wc -l)
for d in mag_validation/out_batches mag_validation/out_batches_tau99 mag_validation/out_batches_n500 mag_validation/lowc_out_batches mag_validation_km/out_batches mag_validation_rnd/out_batches; do
  echo "$d batches=$(ls -d $d/b* | wc -l) cells_files=$(ls $d/b*/mag_cells.csv | wc -l) summaries=$(ls $d/b*/mag_summary.json | wc -l) rows=$(cat $d/b*/mag_cells.csv | grep -vc '^mag,')"; done
[ $nf = 0 ] || { echo "FAILED batches in $nf drivers"; exit 3; }
echo "== concatenate + aggregate"
cat_cells() { local out=$1 d=$2; head -1 "$d/b000/mag_cells.csv" > $out; for b in $(ls -d $d/b* | sort); do tail -n+2 $b/mag_cells.csv >> $out; done; }
cat_cells mag_validation/out/mag_cells_all.csv mag_validation/out_batches
cat_cells mag_validation/out/mag_cells_tau99.csv mag_validation/out_batches_tau99
for f in mag_validation/out/mag_cells_all.csv mag_validation/out/mag_cells_tau99.csv; do echo "$f rows=$(($(wc -l < $f)-1)) md5=$(md5sum < $f | cut -c1-32)"; done
t1=$(date +%s)
$PY mag_validation/aggregate.py > $W/logs/r22_agg1.log 2>&1; echo "aggregate full rc=$?"
$PY mag_validation/aggregate.py --batches mag_validation/out_batches_tau99 --out mag_validation/out/mag_summary_tau99.json > $W/logs/r22_agg2.log 2>&1; echo "aggregate tau99 rc=$?"
$PY mag_validation/aggregate.py --batches mag_validation/out_batches_n500 --out mag_validation/out/summary_n500.json > $W/logs/r22_agg3.log 2>&1; echo "aggregate n500 rc=$?"
$PY mag_validation/aggregate.py --batches mag_validation_km/out_batches --truth region_validation_kmeans/data/benchmark_truth.json --out mag_validation_km/out/summary_km.json > $W/logs/r22_agg4.log 2>&1; echo "aggregate km rc=$?"
echo "aggregate x4 $(( $(date +%s) - t1 )) s"
git status --porcelain | head -40
