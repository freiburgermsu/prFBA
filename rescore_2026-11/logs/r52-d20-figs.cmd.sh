set -u
PY=/home/freiburger/Documents/py_venv/bin/python
T=/home/freiburger/Documents/prFBA_rescore_2026-11/control/r5a/figures/traces
echo "SOURCE_DATE_EPOCH=$SOURCE_DATE_EPOCH PYTHONHASHSEED=$PYTHONHASHSEED"
st() { tag=$1; shift; strace -f -qq -e trace=openat,open -o $T/d20_$tag.strace "$@"; echo "[$tag] rc=$?"; }
cd /home/freiburger/Documents/prFBA/region_validation/scripts && st acc $PY render_accuracy.py --force
cd /home/freiburger/Documents/prFBA/region_validation/scripts && st genecap $PY render_genecap.py --force
ls -la --time-style=+%F_%T /home/freiburger/Documents/prFBA/region_validation/figures
