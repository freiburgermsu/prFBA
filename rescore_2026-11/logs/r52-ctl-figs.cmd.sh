set -u
PY=/home/freiburger/Documents/py_venv/bin/python
S=/home/freiburger/Documents/prFBA_rescore_2026-11/start
C=/home/freiburger/Documents/prFBA_rescore_2026-11/control/r5a/figures/ctl
T=/home/freiburger/Documents/prFBA_rescore_2026-11/control/r5a/figures/traces
echo "SOURCE_DATE_EPOCH=$SOURCE_DATE_EPOCH PYTHONHASHSEED=$PYTHONHASHSEED"
$PY -c "import matplotlib, numpy, pandas; print('mpl', matplotlib.__version__, 'np', numpy.__version__, 'pd', pandas.__version__)"
st() { tag=$1; shift; strace -f -qq -e trace=openat,open -o $T/ctl_$tag.strace "$@"; echo "[$tag] rc=$?"; }
cd $C/rv && st F1_F11 $PY render_accuracy.py --force --only F1,F11
cd $C/rv && st F4b $PY render_genecap.py --force --only F4b
cd $C/rv && st F2alt $PY run_exclude_alt.py
cd $C/mag && st F14 $PY plot_results_d100.py
cd $C/mag && st F15 $PY plot_crossover_d100.py
cd $C/p2 && st F13 $PY render_comparison_d100.py --summary $S/region_validation_d100/picrust2_comparison/data/comparison_summary.json --outdir $C/out
cd $C/two_panel && st F16 $PY plot_space.py
ls -la $C/out
