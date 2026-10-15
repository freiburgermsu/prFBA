#!/bin/bash
# SPEC 3.4 leave-FAMILY-out scoring (5 holdout blocks).  Driver's printed command is
# hardcoded to the genus names; this substitutes the family fold plan, the family
# summary name and --rank-blocked family.
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
R=$SR/data_family
exec ~/Documents/py_venv/bin/python $SR/scripts/score_folds.py --mode fold \
  --fold-plan $R/fold_plan_family.json --panel $R/panel_eligible.json \
  --picrust2 $R/picrust2 --prfba $R/prfba --refs $R/refs \
  --neighbors $R/neighbors.json --cells $R/cells.csv.gz \
  --summary $R/summary_family.json --controls $R/controls.json \
  --rank-blocked family --bootstrap 2000 --expect-densities 0.25,0.5,1.0 \
  --expect-domains bac,arc --expect-nsti primary --expect-forced-answer primary --force
