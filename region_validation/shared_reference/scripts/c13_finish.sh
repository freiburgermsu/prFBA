#!/bin/bash
# Amendment A-C13, after both passes on ModelSEEDpy c2e31a0: compare A-C13 with A-C12 (the
# library's effect) and with the pre-registered run, then rescore the genus summary against
# the A-C13 gapfilled summary. The A-C12-scored summary and controls are kept as *.a-c12.
# The scorer runs on the d884a50 pin, as in A-C12's rescoring, so the phenotype summary is
# the only input that differs.
set -e
cd /home/freiburger/Documents/prFBA/region_validation/shared_reference
D=$PWD/data
PY="$HOME/Documents/py_venv/bin/python -W ignore::SyntaxWarning"
for p in gf atp; do
  until grep -q "^EXIT=" logs/pheno_c13_$p.log 2>/dev/null; do sleep 60; done
  grep -q "^EXIT=0" logs/pheno_c13_$p.log || { echo "pass $p did not finish cleanly"; exit 1; }
done
for p in gf atp; do
  PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-d884a50 $PY scripts/c12_compare.py \
      --before $D/phenotype_${p}_c12fix/per_genome.jsonl --after $D/phenotype_${p}_c13/per_genome.jsonl \
      --out $D/c13_compare_$p.json > logs/c13_compare_$p.log 2>&1
  PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-d884a50 $PY scripts/c12_compare.py \
      --before $D/phenotype_$p/per_genome.jsonl --after $D/phenotype_${p}_c13/per_genome.jsonl \
      --out $D/c13_vs_prereg_compare_$p.json > logs/c13_vs_prereg_compare_$p.log 2>&1
done
for f in summary_genus.json controls.json; do
  [ -e $D/$f.a-c12 ] || cp -p $D/$f $D/$f.a-c12
done
PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-d884a50 $PY scripts/score_folds.py --mode fold \
    --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json \
    --picrust2 $D/picrust2 --prfba $D/prfba --refs $D/refs --neighbors $D/neighbors.json \
    --cells $D/cells.csv.gz --summary $D/summary_genus.json --controls $D/controls.json \
    --rank-blocked genus --bootstrap 2000 --expect-densities 0.25,0.5,1.0 \
    --expect-domains bac,arc --expect-nsti primary --expect-forced-answer primary \
    --require-phenotype-binding --phenotype-summary $D/phenotype_gf_c13/phenotype_summary.json \
    --templates /home/freiburger/Documents/ModelSEEDTemplates-main2294299/templates/v7.0/ \
    > logs/score_genus_c13.log 2>&1 || echo "score_folds exit $?" >> logs/score_genus_c13.log
echo FINISH_DONE
