set -e
cd /home/freiburger/Documents/prFBA/region_validation_d100/scripts
/home/freiburger/Documents/py_venv/bin/python score_taxacc.py --force
/home/freiburger/Documents/py_venv/bin/python mcnemar_regions.py --force
