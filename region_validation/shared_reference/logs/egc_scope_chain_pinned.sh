#!/bin/bash
# Draft-level scope scan and loop enumeration under the pinned ModelSEEDpy (d884a50).
set -e
cd /home/freiburger/Documents/prFBA/region_validation/shared_reference
D=$PWD/data; PY="$HOME/Documents/py_venv/bin/python -W ignore::SyntaxWarning"
TM=/home/freiburger/Documents/ModelSEEDTemplates-main2294299/templates/v7.0
TF=/home/freiburger/Documents/ModelSEEDTemplates-c12run/templates/v7.0
IN="--cells $D/cells.csv.gz --fold-plan $D/fold_plan.json --panel $D/panel_eligible.json --prfba $D/prfba --picrust2 $D/picrust2 --n-genomes 200 --expect-modelseedpy d884a50"
$PY scripts/egc_scope_scan.py $IN --workers 4 --templates $TM --templates $TF --out $D/egc_analysis/scope_scan_200.jsonl > logs/egc_scope_scan_200.log 2>&1
$PY - <<PYEOF
import json
rows = [json.loads(l) for l in open("$D/egc_analysis/scope_scan_200.jsonl")]
leak = {}
for r in rows:
    leak.setdefault((r["src_genome_id"], r["arm"]), {})["c12run" in r["templates_dir"]] = r["leak"]
fix_only = {k for k, v in leak.items() if v.get(True) and v.get(False) is False}
with open("$D/egc_analysis/scope_scan_fixonly.jsonl", "w") as fh:
    for r in rows:
        if (r["src_genome_id"], r["arm"]) in fix_only and "c12run" in r["templates_dir"]:
            fh.write(json.dumps(r, sort_keys=True) + "\n")
print("fix-only leaking drafts:", len(fix_only), "arms:", sorted({a for _, a in fix_only}))
PYEOF
ARMS=$($PY -c "import json;print(','.join(sorted({json.loads(l)['arm'] for l in open('$D/egc_analysis/scope_scan_fixonly.jsonl')})))")
$PY scripts/egc_draft_loops.py $IN --workers 4 --scan $D/egc_analysis/scope_scan_200.jsonl --arms truth --k 5 --templates $TM --out $D/egc_analysis/draft_loops_200.jsonl > logs/egc_draft_loops_200.log 2>&1
if [ -n "$ARMS" ]; then
  $PY scripts/egc_draft_loops.py $IN --workers 2 --scan $D/egc_analysis/scope_scan_fixonly.jsonl --arms $ARMS --k 5 --templates $TF --out $D/egc_analysis/draft_loops_fixonly.jsonl > logs/egc_draft_loops_fixonly.log 2>&1
fi
echo CHAIN_DONE
