#!/bin/bash
# R2.2 W1 verify on TF (inside prfba-ns): (1) the W1 checks (--check-only) of each installed (FM-built) new
# index against the pushed frozen cache (sha 05f03bcd = R1.6) and the old index; (2) an independent TF
# rebuild of each index (current values, and --restore-old-values) into $W/control/panel_index_tf/, which
# must be byte-identical to the installed new index and to the pre-repair (old) index respectively.
set -uo pipefail
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11; PY=/home/freiburger/Documents/py_venv/bin/python
O=$W/control/panel_index_tf; mkdir -p $O; cd $R; fail=0
for pair in mag_validation/pgfam_panel.json:pgfam_panel kmeans_panel/pgfam_panel_kmeans.json:pgfam_panel_kmeans random_panel/pgfam_panel_random.json:pgfam_panel_random; do
  p=${pair%%:*}; n=${pair##*:}; old=$W/pre_repair/$p
  $PY mag_validation/build_pgfam_panel_index.py --old-index $old --out $R/$p --check-only || { echo "CHECK FAIL $p"; fail=1; }
  $PY mag_validation/build_pgfam_panel_index.py --old-index $old --out $O/$n.new.json || fail=1
  cmp $O/$n.new.json $R/$p && echo "BYTE-IDENTICAL TF rebuild $n.new.json == installed $p" || { echo "BYTES DIFFER $n.new"; fail=1; }
  $PY mag_validation/build_pgfam_panel_index.py --old-index $old --out $O/$n.pre_repair.json --restore-old-values || fail=1
  cmp $O/$n.pre_repair.json $old && echo "BYTE-IDENTICAL TF rebuild $n.pre_repair.json == old $p" || { echo "BYTES DIFFER $n.pre_repair"; fail=1; }
done
[ $fail = 0 ] && echo "W1 VERIFY PASS" || echo "W1 VERIFY FAIL"; exit $fail
