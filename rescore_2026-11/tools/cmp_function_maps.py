#!/usr/bin/env python
"""Content comparison of two function_maps.json files (re-score W5 / D3 pin control).

    cmp_function_maps.py NEW.json REF.json [--json report.json]

Equal iff pgfam_ec, pgfam_rxn and ec_rxn are identical and every `_meta` key present in both is
equal except the path/provenance keys (modelseed_db, modelseed_db_git, modelseed_templates_dir,
modelseed_templates_git). Reports the number of differing keys per map, with examples.
"""
import argparse
import json
import sys

import orjson

PROVENANCE = {"modelseed_db", "modelseed_db_git", "modelseed_templates_dir", "modelseed_templates_git"}


def compare(new_path, ref_path):
    a = orjson.loads(open(new_path, "rb").read())
    b = orjson.loads(open(ref_path, "rb").read())
    rep = {"new": new_path, "ref": ref_path, "maps": {}}
    ok = True
    for m in ("pgfam_ec", "pgfam_rxn", "ec_rxn"):
        x, y = a.get(m, {}), b.get(m, {})
        diff = sorted(k for k in set(x) | set(y) if x.get(k) != y.get(k))
        rep["maps"][m] = {"n_new": len(x), "n_ref": len(y), "n_diff_keys": len(diff),
                          "examples": [{"key": k, "new": x.get(k), "ref": y.get(k)} for k in diff[:5]]}
        ok &= not diff
    ma, mb = a.get("_meta", {}), b.get("_meta", {})
    meta_diff = {k: {"new": ma.get(k), "ref": mb.get(k)} for k in (set(ma) & set(mb)) - PROVENANCE
                 if ma.get(k) != mb.get(k)}
    rep["meta_diff"] = meta_diff
    rep["meta_only_new"] = sorted(set(ma) - set(mb))
    rep["provenance_new"] = {k: ma.get(k) for k in PROVENANCE if k in ma}
    ok &= not meta_diff
    rep["equal_content"] = ok
    return ok, rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("new")
    ap.add_argument("ref")
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    ok, rep = compare(a.new, a.ref)
    if a.json:
        json.dump(rep, open(a.json, "w"), indent=1)
    print(json.dumps({k: rep[k] for k in ("equal_content", "maps", "meta_diff", "provenance_new")})[:3000])
    print("EQUAL CONTENT" if ok else "CONTENT DIFFERS")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
