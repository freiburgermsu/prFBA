#!/usr/bin/env bash
# W8 on the regenerated SR cell tables: r35_w8.sh <root> [<root> ...]  (data_c2 | data_family | data)
set -u
PY=/home/freiburger/Documents/py_venv/bin/python
W=/home/freiburger/Documents/prFBA_rescore_2026-11; SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
B=$W/pre_repair/region_validation/shared_reference; C=$W/control/r3
rc=0
for root in "$@"; do
  f=cells.csv.gz; [ $root = data_c2 ] && f=cells_incself.csv.gz
  echo "== W8 $root $(date -Is)"
  $PY $W/tools/attrib_sr_cells.py --root $root --old $B/$root/$f --new $SR/$root/$f --out $C/w8_cells_$root.json || rc=1
done
exit $rc
