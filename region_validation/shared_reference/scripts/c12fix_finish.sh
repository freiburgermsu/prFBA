#!/bin/bash
# Amendment A-C12, after both pinned passes: before/after comparisons, then the genus fold
# rescoring that settles C11/C12 from the amended gapfilled summary (SPEC section 14).
# The family summary is not rescored: it was scored without the phenotype arm and carries
# no C11/C12. The pre-registered summary and controls are kept as *.pre-c12fix.
set -e
cd /home/freiburger/Documents/prFBA/region_validation/shared_reference
D=$PWD/data
PY="$HOME/Documents/py_venv/bin/python -W ignore::SyntaxWarning"
export PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-d884a50    # the pin, for uniformity
for p in gf atp; do
  grep -q "^EXIT=0" logs/pheno_c12fix_$p.log || { echo "pass $p did not finish cleanly"; exit 1; }
  $PY scripts/c12_compare.py --before $D/phenotype_$p/per_genome.jsonl \
      --after $D/phenotype_${p}_c12fix/per_genome.jsonl --out $D/c12fix_compare_$p.json \
      > logs/c12fix_compare_$p.log 2>&1
done
for f in summary_genus.json controls.json; do
  [ -e $D/$f.pre-c12fix ] || cp -p $D/$f $D/$f.pre-c12fix
done
$PY scripts/score_folds.py --mode fold --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json \
    --picrust2 $D/picrust2 --prfba $D/prfba --refs $D/refs --neighbors $D/neighbors.json \
    --cells $D/cells.csv.gz --summary $D/summary_genus.json --controls $D/controls.json \
    --rank-blocked genus --bootstrap 2000 --expect-densities 0.25,0.5,1.0 \
    --expect-domains bac,arc --expect-nsti primary --expect-forced-answer primary \
    --require-phenotype-binding --phenotype-summary $D/phenotype_gf_c12fix/phenotype_summary.json \
    --templates /home/freiburger/Documents/ModelSEEDTemplates-main2294299/templates/v7.0/ \
    > logs/score_genus_c12fix.log 2>&1 || echo "score_folds exit $?" >> logs/score_genus_c12fix.log
echo FINISH_DONE
