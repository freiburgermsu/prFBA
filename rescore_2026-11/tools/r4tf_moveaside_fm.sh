#!/bin/bash
# R4.3 TF lane: move FM's old outputs of calibrate / Stage E / Stage O / FVA / wrapper / verification to
# $W/pre_repair (same repo-relative paths) before installing TF's new ones. Prints sha256 of everything moved.
# usage: r4tf_moveaside_fm.sh calibration|pfunc_estimators|objective_tuning|objective_tuning_fva
set -euo pipefail
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
SRD=region_validation/shared_reference/data; P=$W/pre_repair/$SRD
cd $R
case $1 in
  calibration) files="$SRD/calibration/README.md $SRD/calibration/assertions.json $SRD/calibration/reliability.csv $SRD/calibration/scores.json";;
  pfunc_estimators) files="$(ls $SRD/pfunc_estimators/{README.md,assertions.json,run.log,scores.json,selected.json,stage_o_weights.json,stage_o_weights_E1.json,stage_o_weights_E3.json}) $(ls $SRD/pfunc_estimators/verification/{p0_sanity.json,p1_stream.log,p2_pass.log,p2_summary.json,p3_stage_o.json,p3_stage_o.log,p4_analyze.log,p4_results.json,p5_compare.json,p6_extra.json,p8_leaky.json})";;
  objective_tuning) files="$(ls $SRD/objective_tuning/{meta.json,per_amplicon.jsonl,per_amplicon.jsonl.gz,summary.json,OBJECTIVE_TUNING_RESULTS.filled.md,wrapper_equivalence.json} $SRD/objective_tuning/verification/out_*.json)";;
  objective_tuning_fva) files="$(ls $SRD/objective_tuning_fva/{meta.json,per_amplicon.jsonl,per_amplicon.jsonl.gz})";;
  *) echo "unknown $1"; exit 2;;
esac
echo "[$(date -Is)] move aside ($1) -> $P"
for f in $files; do
  [ -e "$P/${f#$SRD/}" ] && { echo "already in pre_repair: $f"; exit 3; }
  sha256sum "$f"; mkdir -p "$(dirname "$P/${f#$SRD/}")"; mv "$f" "$P/${f#$SRD/}"
done
