#!/usr/bin/env python
"""Old -> new headline of an as-shipped comparison_summary.json (EC and rxn, all regions)."""
import json, sys
old, new = (json.load(open(p)) for p in sys.argv[1:3])
out = {}
for space in ("ec", "rxn"):
    for m in ("picrust2", "prfba_exclude", "prfba_include", "prfba_exclude_legacy", "prfba_include_legacy"):
        def g(s):
            b = s["by_region"][space]["all"]["all"].get(m)
            if not b: return None
            c, p = b["conditional"], b["population"]
            return {"n": b["n"], "coverage": b["coverage"], "cond_f1": c["f1"]["value"] if c["f1"] else None,
                    "cond_precision": c["precision"]["value"] if c["precision"] else None,
                    "cond_recall": c["recall"]["value"] if c["recall"] else None,
                    "pop_f1": p["f1"]["value"] if p["f1"] else None, "pop_recall": p["recall"]["value"] if p["recall"] else None}
        out[f"{space}|{m}"] = {"old": g(old), "new": g(new)}
for k, v in out.items():
    o, n = v["old"], v["new"]
    if o and n:
        print(f"{k:28} " + "  ".join(f"{f}={o[f]}->{n[f]}" for f in o if o[f] != n[f] or f == "cond_f1"))
if len(sys.argv) > 3:
    json.dump(out, open(sys.argv[3], "w"), indent=1)
