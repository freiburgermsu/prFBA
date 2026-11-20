set -u
PY=/home/freiburger/Documents/py_venv/bin/python
R=/home/freiburger/Documents/prFBA
T=/home/freiburger/Documents/prFBA_rescore_2026-11/control/r5a/figures/traces
echo "SOURCE_DATE_EPOCH=$SOURCE_DATE_EPOCH PYTHONHASHSEED=$PYTHONHASHSEED"
st() { tag=$1; shift; strace -f -qq -e trace=openat,open -o $T/live_$tag.strace "$@"; echo "[$tag] rc=$?"; }
cd $R/manuscript_figures/rv && st F1_F11 $PY render_accuracy.py --force --only F1,F11
cd $R/manuscript_figures/rv && st F4b $PY render_genecap.py --force --only F4b
cd $R/manuscript_figures/rv && st F2alt $PY run_exclude_alt.py
cd $R/manuscript_figures/mag && st F14 $PY plot_results_d100.py
cd $R/manuscript_figures/mag && st F15 $PY plot_crossover_d100.py
cd $R/manuscript_figures/p2 && st F13 $PY render_comparison_d100.py
cd $R/two_panel && st F16 $PY plot_space.py
cp -p $R/two_panel/F_panel_space.pdf $R/two_panel/F_panel_space.png $R/manuscript_figures/out/
ls -la $R/manuscript_figures/out/
