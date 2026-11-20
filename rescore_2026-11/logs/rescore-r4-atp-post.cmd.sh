# R4.1 ATP pass post-checks (FM): c12_compare old -> new (plan R4.1), then the row attribution control and the
# manual-number recipe (tools in $W/tools).
set -o pipefail
W=/home/freiburger/Documents/prFBA_rescore_2026-11
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
cd $SR
PY=/home/freiburger/Documents/py_venv/bin/python
$PY -W ignore::SyntaxWarning scripts/c12_compare.py --before $W/pre_repair/region_validation/shared_reference/data/phenotype_atp_c13/per_genome.jsonl \
  --after $SR/data/phenotype_atp_c13/per_genome.jsonl --out $SR/data/rescore_compare_atp.json > $SR/logs/rescore_compare_atp.log 2>&1
echo "c12_compare rc=$?"
$PY $W/tools/r4_pheno_attrib.py --pass atp --out $W/control/r4/pheno_attrib_atp.json; echo "attrib rc=$?"
$PY $W/tools/r4_pheno_manual_numbers.py $W/control/r4/pheno_manual_numbers_atp.json; echo "manual rc=$?"
