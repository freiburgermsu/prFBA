# R4.2 genus score pass 2 (FM): ONLY the score_folds line of scripts/c13_finish.sh (A-C13 rescoring), no --force,
# pins as c13_finish.sh / Run 3 pass 1 (PYTHONPATH ModelSEEDpy-d884a50, templates ModelSEEDTemplates-main2294299),
# --require-phenotype-binding against the NEW data/phenotype_gf_c13/phenotype_summary.json (named explicitly).
# The c12_compare lines of c13_finish.sh are NOT run (c13_compare_*.json / c13_vs_prereg_*.json are A-C13 history).
set -o pipefail
cd /home/freiburger/Documents/prFBA/region_validation/shared_reference
D=$PWD/data
PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-d884a50 /home/freiburger/Documents/py_venv/bin/python -W ignore::SyntaxWarning scripts/score_folds.py --mode fold \
    --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json \
    --picrust2 $D/picrust2 --prfba $D/prfba --refs $D/refs --neighbors $D/neighbors.json \
    --cells $D/cells.csv.gz --summary $D/summary_genus.json --controls $D/controls.json \
    --rank-blocked genus --bootstrap 2000 --expect-densities 0.25,0.5,1.0 \
    --expect-domains bac,arc --expect-nsti primary --expect-forced-answer primary \
    --require-phenotype-binding --phenotype-summary $D/phenotype_gf_c13/phenotype_summary.json \
    --templates /home/freiburger/Documents/ModelSEEDTemplates-main2294299/templates/v7.0/
rc=$?; echo "EXIT=$rc"; exit $rc
