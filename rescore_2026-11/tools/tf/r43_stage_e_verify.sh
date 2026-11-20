#!/bin/bash
# R4.3 Stage E independent verification on TF: data/pfunc_estimators/verification/p0..p8 in VERIFICATION.md's order
# (cwd = the verification dir; logs as the committed ones: p1_stream.log .. p4_analyze.log; p0/p5/p6/p7/p8 to their own logs).
set -uo pipefail
V=/home/freiburger/Documents/prFBA/region_validation/shared_reference/data/pfunc_estimators/verification
cd $V; rm -rf cache
rc_all=0
run() { local n=$1; shift; local t0=$(date +%s); echo "== $n start $(date -Is)"; { time python -u "$@"; } > $V/$n.log 2>&1; local rc=$?; echo "== $n rc=$rc $(( $(date +%s) - t0 )) s"; [ $rc = 0 ] || rc_all=$rc; }
run p0_sanity p0_sanity.py
run p1_stream p1_stream.py
run p2_pass p2_pass.py
run p3_stage_o p3_stage_o.py
run p4_analyze p4_analyze.py
run p5_compare p5_compare.py
run p6_extra p6_extra.py
run p7_k10 p7_k10.py
run p8_leaky p8_leaky.py
exit $rc_all
