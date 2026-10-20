#!/bin/bash
# C2 (include-self ceiling) for the leave-genus-out campaign, SPEC section 7 / 11.2:
# d100 only, all 5 folds x both domains.  run_campaign.py covers stages 1-6 only, so this
# script is the record of how C2 was produced.  Outputs live in a separate root
# (data_c2/) so data/'s run inventory, cells.csv.gz fingerprint and the phenotype
# binding are untouched; only data/controls.json receives the merged C2 entry.
set -u
PY=~/Documents/py_venv/bin/python
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
D=$SR/data; C=$SR/data_c2; S=$SR/scripts
mkdir -p $C/{refs,picrust2,prfba,logs}
[ -f $D/controls.json.pre-C2 ] || cp -p $D/controls.json $D/controls.json.pre-C2
echo "[c2] start $(date -Is)"

# 1. include-self references (whole eligible panel; fold-independent)
for DOM in arc bac; do
  if [ ! -f $C/refs/incself_${DOM}_d100/manifest.json ]; then
    $PY $S/build_reference.py --mode include_self --domain $DOM --density 1.0 \
      --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json --outdir $C/refs \
      --threads 16 > $C/logs/build_incself_${DOM}_d100.log 2>&1 || { echo "[c2] build $DOM failed"; exit 2; }
  fi
  $PY - "$C/refs/incself_${DOM}_d100/manifest.json" <<'EOF' || exit 2
import json, sys
m = json.load(open(sys.argv[1]))
flat = json.dumps(m)
print("[c2] manifest", sys.argv[1], {k: m.get(k) for k in ("mode", "n_tips", "set_equality_ok", "msa_mean_gap_frac") if k in m})
if m.get("mode") not in (None, "include_self") or m.get("set_equality_ok") is False:
    sys.exit("[c2] manifest check failed")
EOF
done
echo "[c2] refs done $(date -Is)"

# 2. PICRUSt2, forced routing, one per fold and domain
for F in 0 1 2 3 4; do N=$(printf %02d $F)
  for DOM in bac arc; do
    O=$C/picrust2/incself_f${N}_${DOM}_d100
    [ -f $O/predictions.json ] && continue
    $PY $S/run_picrust2_fold.py --mode include_self --routing forced --fold $F --domain $DOM \
      --density 1.0 --ref-dir $C/refs/incself_${DOM}_d100 --fold-plan $D/fold_plan.json \
      --panel $D/panel_eligible.json --outdir $O --processes 14 \
      > $C/logs/p2_incself_f${N}_${DOM}.log 2>&1 || echo "[c2] picrust2 f$N $DOM rc=$?"
  done
done
echo "[c2] picrust2 done $(date -Is)"

# 3. prFBA, include_self, the stage-6 settings from data/ledger.jsonl
for F in 0 1 2 3 4; do N=$(printf %02d $F)
  O=$C/prfba/incself_f${N}_d100
  [ -f $O/predictions.json ] && continue
  $PY $S/run_prfba_fold.py --mode include_self --fold $F --density 1.0 \
    --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json \
    --ref-manifest $C/refs/incself_bac_d100/manifest.json \
    --ref-manifest $C/refs/incself_arc_d100/manifest.json \
    --outdir $O --select-modes legacy,exact_tie,graduated,exact_tie_guards \
    --k2 200 --gpu-cand 400 --prefilter-k 1500 --force-backend cuda --workers 30 \
    > $C/logs/prfba_incself_f${N}_d100.log 2>&1 || echo "[c2] prfba f$N rc=$?"
  grep -h "\[C2\]" $C/logs/prfba_incself_f${N}_d100.log | tail -2
done
echo "[c2] prfba done $(date -Is)"

# 4. score; merges C2 into data/controls.json (never with a phenotype summary)
$PY $S/score_folds.py --mode include_self --fold-plan $D/fold_plan.json \
  --panel $D/panel_eligible.json --picrust2 $C/picrust2 --prfba $C/prfba --refs $C/refs \
  --neighbors $D/neighbors.json --cells $C/cells_incself.csv.gz \
  --summary $C/summary_incself.json --controls $D/controls.json \
  --rank-blocked genus --bootstrap 2000 --force > $C/logs/score_incself.log 2>&1
echo "[c2] score exit=$? $(date -Is)"
