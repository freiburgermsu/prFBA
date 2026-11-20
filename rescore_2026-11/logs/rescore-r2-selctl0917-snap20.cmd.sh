set -u
PY=/home/freiburger/Documents/py_venv/bin/python
SC=/home/freiburger/Documents/prFBA/region_validation/picrust2_comparison/scripts; D=/home/freiburger/Documents/prFBA/region_validation/picrust2_comparison/data
for M in default legacy; do
  X=""; T=d20_snap20_0917; C=$D/prfba_selection.json
  if [ $M = legacy ]; then X=--legacy; T=d20_legacy_snap20_0917; C=$D/prfba_selection_legacy.json; fi
  $PY /home/freiburger/Documents/prFBA_rescore_2026-11/tools/ctl_run.py --common /home/freiburger/Documents/prFBA/region_validation/scripts --cache /home/freiburger/Documents/prFBA_rescore_2026-11/control/r2/snap20_d20_0917.json --record /home/freiburger/Documents/prFBA_rescore_2026-11/control/r2/prfba_sel_$T.ctl.json -- $SC/prfba_selection.py --out /home/freiburger/Documents/prFBA_rescore_2026-11/control/r2/prfba_sel_$T.json --workers 24 $X
  $PY /home/freiburger/Documents/prFBA_rescore_2026-11/tools/cmp_prfba_selection.py --old $C --new /home/freiburger/Documents/prFBA_rescore_2026-11/control/r2/prfba_sel_$T.json --out /home/freiburger/Documents/prFBA_rescore_2026-11/control/r2/cmp_$T.json --cache-keys /home/freiburger/Documents/prFBA_rescore_2026-11/stamps/R1.6.cache_keys.txt --d20 && echo "[selctl0917] $T RECORDS IDENTICAL" || echo "[selctl0917] $T RECORDS DIFFER"
done
