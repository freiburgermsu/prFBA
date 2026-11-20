set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; CTL=/home/freiburger/Documents/prFBA_rescore_2026-11/control/region_validation_ctl; C20=/home/freiburger/Documents/prFBA_rescore_2026-11/control/d20ctl/snap20j_copy.json
CR() { $PY $W/tools/ctl_run.py --common $CTL/scripts --cache $C20 --allow-cache-writes --record "$1" -- "${@:2}"; }
st() { stat -c '[stat] %n size=%s mtime=%.6Y' "$@"; }
[ -s $CTL/data/selection.json ] || { echo "no C-L4 selection in the control tree"; exit 2; }
st $C20 $R/bvbrc_cache/genome_gene_families.json
cd $CTL/scripts
echo "== C-D20 gap close $(date -Is)"
CR $W/control/d20ctl/gap.ctl.json $R/kmeans_panel/close_exclude_rep_gap.py --tree region_validation --legacy 2>&1 | tee $CTL/data/close_exclude_rep_gap.log
grep -q "data $CTL/data" $CTL/data/close_exclude_rep_gap.log || { echo "GAP CLOSE DID NOT RUN ON THE CONTROL TREE"; exit 3; }
st $C20
echo "== C-D20 score_taxacc --legacy --force $(date -Is)"
CR $W/control/d20ctl/taxacc.ctl.json $CTL/scripts/score_taxacc.py --legacy --force
echo "== C-D20 score_genecap --legacy --force $(date -Is)"
CR $W/control/d20ctl/genecap.ctl.json $CTL/scripts/score_genecap.py --legacy --force
st $C20
echo "== C-D20 mcnemar --force $(date -Is)"
CR $W/control/d20ctl/mcnemar.ctl.json $CTL/scripts/mcnemar_regions.py --force
echo "== C-D20 diagnose_self_recovery $(date -Is)"
CR $W/control/d20ctl/diag.ctl.json $CTL/scripts/diagnose_self_recovery.py
st $C20 $R/bvbrc_cache/genome_gene_families.json
sha256sum $C20
echo "== C-D20 done $(date -Is)"
