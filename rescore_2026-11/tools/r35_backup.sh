#!/usr/bin/env bash
# R3.5 step 0: copy every scorer output that R3.5 rewrites to pre_repair (W8 input + rollback), once.
set -euo pipefail
W=/home/freiburger/Documents/prFBA_rescore_2026-11; SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
B=$W/pre_repair/region_validation/shared_reference
for f in data/cells.csv.gz data/cells.csv.gz.fp data/summary_genus.json data/controls.json data/summary_recipe.json \
         data/cells_recipe.csv.gz data_c2/cells_incself.csv.gz data_c2/cells_incself.csv.gz.fp data_c2/summary_incself.json \
         data_family/cells.csv.gz data_family/cells.csv.gz.fp data_family/summary_family.json data_family/controls.json; do
  mkdir -p $(dirname $B/$f)
  if [ -e $B/$f ]; then cmp -s $SR/$f $B/$f || { echo "backup $B/$f exists and differs" >&2; exit 2; }; continue; fi
  cp -p $SR/$f $B/$f
done
(cd $B && sha256sum data/cells.csv.gz data/cells.csv.gz.fp data/summary_genus.json data/controls.json data/summary_recipe.json data/cells_recipe.csv.gz data_c2/* data_family/cells.csv.gz* data_family/summary_family.json data_family/controls.json 2>/dev/null)
