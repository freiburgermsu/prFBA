#!/usr/bin/env python
"""Rows where the selected arm's total exchange flux exceeds the shipped arm's (OBJECTIVE_TUNING_RESULTS.md §6,
'imports more on N% of rows'; the 11-03 verification counted 95/167 tune, 110/193 confirm). Ad hoc in 11-03; scripted here.
usage: otune_exchange_rise.py PER_AMPLICON_JSONL OUT_JSON [SELECTED] [SHIPPED]"""
import json, sys
p, out = sys.argv[1], sys.argv[2]
sel = sys.argv[3] if len(sys.argv) > 3 else "est/lin/G0/X1"
sh = sys.argv[4] if len(sys.argv) > 4 else "kn/lin/G0/X100"
res = {}
for line in open(p):
    if not line.strip():
        continue
    r = json.loads(line)
    if r.get("status") != "ok":
        continue
    a, b = r["arms"][sel]["flux_exchange"], r["arms"][sh]["flux_exchange"]
    h = res.setdefault(r["half"], {"n": 0, "gt_1e-9_rel": 0, "gt_0": 0, "lt_minus_1e-9_rel": 0, "max_ratio": 0.0})
    h["n"] += 1
    h["gt_0"] += a > b
    h["gt_1e-9_rel"] += a > b * (1 + 1e-9)
    h["lt_minus_1e-9_rel"] += a < b * (1 - 1e-9)
    h["max_ratio"] = max(h["max_ratio"], a / b if b else float("nan"))
for h in res.values():
    h["share_gt_1e-9_rel"] = h["gt_1e-9_rel"] / h["n"]
json.dump({"selected": sel, "shipped": sh, "by_half": res}, open(out, "w"), indent=1)
print(json.dumps(res))
