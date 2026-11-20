# R4.1 gapfilled pass post-checks (FM): c12_compare old -> new (plan R4.1), the row attribution control, and the
# sample-change analysis for both passes.
set -o pipefail
W=/home/freiburger/Documents/prFBA_rescore_2026-11
SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
cd $SR
PY=/home/freiburger/Documents/py_venv/bin/python
$PY -W ignore::SyntaxWarning scripts/c12_compare.py --before $W/pre_repair/region_validation/shared_reference/data/phenotype_gf_c13/per_genome.jsonl \
  --after $SR/data/phenotype_gf_c13/per_genome.jsonl --out $SR/data/rescore_compare_gf.json > $SR/logs/rescore_compare_gf.log 2>&1
echo "c12_compare rc=$?"
$PY $W/tools/r4_pheno_attrib.py --pass gf --out $W/control/r4/pheno_attrib_gf.json; echo "attrib rc=$?"
$PY $W/tools/r4_pheno_sample_diff.py $W/control/r4/pheno_sample_diff.json; echo "sample_diff rc=$?"
