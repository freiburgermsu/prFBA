cd /home/freiburger/Documents/prFBA
PY=/home/freiburger/Documents/py_venv/bin/python
rc=0
for t in region_validation_d100_legacy region_validation_kmeans_d100_legacy region_validation_random_d100_legacy; do echo "== copy $t"; $PY region_validation_d100/check_d100.py copy $t || rc=1; done
for t in region_validation_d100 region_validation_kmeans_d100 region_validation_random_d100 region_validation_d100_legacy region_validation_kmeans_d100_legacy region_validation_random_d100_legacy; do echo "== post $t"; $PY region_validation_d100/check_d100.py post $t || rc=1; done
echo "== baseline"; $PY region_validation_d100/check_d100.py baseline || rc=1
echo "LOCAL_CHECKS_RC=$rc"
