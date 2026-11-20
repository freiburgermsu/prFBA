PY=/home/freiburger/Documents/py_venv/bin/python; R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
one() {  # tree cache legacyflag
  T=$1; C=$2; L=${3:-}
  echo "== reps old $T cache=$(basename $C) $L $(date -Is)"
  $PY $W/tools/ctl_run.py --common $R/$T/scripts --cache $C --record $W/control/w8/reps_${T}_old.ctl.json --      $W/tools/dump_genecap_reps.py --out $W/control/w8/reps_${T}_old.pkl $L --keep-tables $W/control/w8/tables_${T}_old
  echo "rc=$?"
  for f in gene_capture_cells.csv gene_capture_by_region.json gene_capture_long.csv; do
    if cmp -s $W/control/w8/tables_${T}_old/$f $R/$T/data/$f; then echo "[repro] $T $f IDENTICAL"; else echo "[repro] $T $f DIFFERS"; fi
  done
}
one region_validation_kmeans_d100 /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair.json
one region_validation_random_d100 /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair.json
one region_validation_d100_legacy /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair.json --legacy
one region_validation_kmeans_d100_legacy /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair.json --legacy
one region_validation_random_d100_legacy /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair.json --legacy
one region_validation /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair_d20.json --legacy
one region_validation_kmeans /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair.json
one region_validation_random /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair.json
stat -c '[stat] %n size=%s mtime=%.6Y' $R/bvbrc_cache/genome_gene_families.json
echo "== reps old done $(date -Is)"
