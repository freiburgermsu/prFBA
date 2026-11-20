#!/usr/bin/env python
"""The .fp answer-key check: would score_folds reuse the OLD (pre-repair) cells table under the NEW
inventory without --force? Uses score_folds.cells_reuse_verdict itself, with the inventory recorded in the
new summary. Expected: old -> refused (inventory/answer key differ); new -> reusable (fingerprints match)."""
import json, sys, types
SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference"
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
B = f"{W}/pre_repair/region_validation/shared_reference"
sys.path.insert(0, f"{SR}/scripts")
import score_folds as S  # noqa: E402
out = {}
for root, cells, summ in (("data_c2", "cells_incself.csv.gz", "summary_incself.json"),
                          ("data_family", "cells.csv.gz", "summary_family.json"),
                          ("data", "cells.csv.gz", "summary_genus.json")):
    if len(sys.argv) > 1 and root not in sys.argv[1:]:
        continue
    s = json.load(open(f"{SR}/{root}/{summ}"))
    inv = s["accounting"]["completeness"]["inventory"]
    old_fp = json.load(open(f"{B}/{root}/{cells}.fp"))
    new_fp = json.load(open(f"{SR}/{root}/{cells}.fp"))
    vo = S.cells_reuse_verdict(types.SimpleNamespace(cells=f"{B}/{root}/{cells}"), inv)
    vn = S.cells_reuse_verdict(types.SimpleNamespace(cells=f"{SR}/{root}/{cells}"), inv)
    out[root] = {"old_fp_answer_key": old_fp.get("answer_key_fp"), "new_fp_answer_key": new_fp.get("answer_key_fp"),
                 "summary_answer_key": inv.get("answer_key_fp"), "old_scorer_code_fp": old_fp.get("scorer_code_fp"),
                 "new_scorer_code_fp": new_fp.get("scorer_code_fp"),
                 "old_cells_reusable": vo[0], "old_reason": vo[1], "new_cells_reusable": vn[0], "new_reason": vn[1]}
    print(root, "old reusable:", vo[0], "| new reusable:", vn[0], "| answer key", (old_fp.get("answer_key_fp") or "")[:12], "->", (inv.get("answer_key_fp") or "")[:12])
json.dump(out, open(f"{W}/control/r3/fp_check{'_' + '_'.join(sys.argv[1:]) if len(sys.argv) > 1 else ''}.json", "w"), indent=1)
ok = all((not v["old_cells_reusable"]) and v["new_cells_reusable"] and v["old_fp_answer_key"] != v["new_fp_answer_key"] for v in out.values())
sys.exit(0 if ok else 1)
