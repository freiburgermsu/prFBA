#!/usr/bin/env python
"""Diff two genome_function_sets.json (and their function_maps.json) for the R2.4/R2.5 checks.

    gfs_diff.py --old-gfs O --new-gfs N --old-maps OM --new-maps NM --out report.json
                [--old-products OP --new-products NP] [--repaired ...] [--d20]
                [--cache CACHE] [--genomes-out changed_nonrepaired.json]

Maps: pgfam_ec / pgfam_rxn / ec_rxn keys added, removed and changed (an existing PGFam's
projection must never change; only PGFams whose products were newly fetched may be added), and the
BV-BRC EC vocabulary (all ECs of pgfam_ec) and the shared vocabulary (that intersected with
PICRUSt2's reference EC columns, as compare.py / annotation_ceiling.py define it) old vs new.

GFS: per genome, identical / changed / only_old / only_new. A changed genome must be a repaired
genome (its cache entry changed) or -- with --cache and --old-products -- a genome whose PGFams
include one whose product was not in the old products file (newly fetched: D10 at depth 20, or a
family of a newly selected reference). Anything else is `unexplained`.
"""
import argparse
import json
import os
import sys
from collections import Counter

import orjson

W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
P2C_SCRIPTS = "/home/freiburger/Documents/prFBA/region_validation/picrust2_comparison/scripts"


def maps_diff(om, nm):
    out = {}
    for m in ("pgfam_ec", "pgfam_rxn", "ec_rxn"):
        a, b = om[m], nm[m]
        changed = sorted(k for k in set(a) & set(b) if a[k] != b[k])
        out[m] = {"n_old": len(a), "n_new": len(b), "n_added": len(set(b) - set(a)),
                  "n_removed": len(set(a) - set(b)), "n_changed": len(changed),
                  "changed_examples": changed[:10]}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old-gfs", required=True)
    ap.add_argument("--new-gfs", required=True)
    ap.add_argument("--old-maps", required=True)
    ap.add_argument("--new-maps", required=True)
    ap.add_argument("--old-products")
    ap.add_argument("--cache")
    ap.add_argument("--repaired", default=f"{W}/expected/repaired_ids.json")
    ap.add_argument("--d20", action="store_true")
    ap.add_argument("--out", required=True)
    ap.add_argument("--genomes-out", help="JSON list of changed genomes explained by new products")
    a = ap.parse_args(argv)

    r = json.load(open(a.repaired))
    repaired = set(r["repaired_2026_11_13"]) | (set(r["poisoned_d20"]) if a.d20 else set())
    om = orjson.loads(open(a.old_maps, "rb").read())
    nm = orjson.loads(open(a.new_maps, "rb").read())
    rep = {"maps": maps_diff(om, nm)}
    sys.path.insert(0, P2C_SCRIPTS)
    from annotation_ceiling import load_reference_ec  # PICRUSt2's own EC columns
    _ref, ref_vocab = load_reference_ec()
    ov = {e for v in om["pgfam_ec"].values() for e in v}
    nv = {e for v in nm["pgfam_ec"].values() for e in v}
    rep["vocab"] = {"bvbrc_old": len(ov), "bvbrc_new": len(nv),
                    "bvbrc_added": sorted(nv - ov), "bvbrc_removed": sorted(ov - nv),
                    "shared_old": len(ov & ref_vocab), "shared_new": len(nv & ref_vocab),
                    "shared_added": sorted((nv - ov) & ref_vocab),
                    "shared_removed": sorted((ov - nv) & ref_vocab)}
    new_pgfams = set(nm["pgfam_ec"]) - set(om["pgfam_ec"])
    del om, nm

    old = orjson.loads(open(a.old_gfs, "rb").read())
    new = orjson.loads(open(a.new_gfs, "rb").read())
    rep["meta_old"], rep["meta_new"] = old.pop("_meta", None), new.pop("_meta", None)
    only_old = set(old) - set(new)
    only_new = set(new) - set(old)
    changed = {g for g in set(old) & set(new) if old[g] != new[g]}
    fields = Counter(f for g in changed for f in ("ec", "rxn_ec", "rxn_role", "n_pgfam")
                     if old[g].get(f) != new[g].get(f))
    rep["gfs"] = {"n_old": len(old), "n_new": len(new), "n_only_old": len(only_old),
                  "n_only_new": len(only_new), "n_changed": len(changed),
                  "changed_fields": dict(fields),
                  "n_identical": len(set(old) & set(new)) - len(changed),
                  "only_new_repaired": len(only_new & repaired)}
    ec_gain = Counter()
    for g in changed:
        o, n = set(old[g]["ec"]), set(new[g]["ec"])
        ec_gain["gain_only" if n > o else "loss_only" if n < o else "same_ec" if n == o else "mixed"] += 1
    rep["gfs"]["ec_direction"] = dict(ec_gain)
    ch_rep = changed & repaired
    rest = changed - repaired
    rep["gfs"]["n_changed_repaired"] = len(ch_rep)
    explained_new_products = set()
    if rest and a.cache and a.old_products:
        oprod = orjson.loads(open(a.old_products, "rb").read())
        cache = orjson.loads(open(a.cache, "rb").read())
        for g in rest:
            pg = cache.get(g) or []
            if any(p not in oprod for p in pg) and old[g]["n_pgfam"] == new[g]["n_pgfam"]:
                explained_new_products.add(g)
        del cache, oprod
    unexplained = rest - explained_new_products
    rep["gfs"]["n_changed_new_products"] = len(explained_new_products)
    rep["gfs"]["n_changed_unexplained"] = len(unexplained)
    rep["gfs"]["unexplained_examples"] = [{"genome": g, "old": {k: (v if not isinstance(v, list) else len(v)) for k, v in old[g].items()},
                                           "new": {k: (v if not isinstance(v, list) else len(v)) for k, v in new[g].items()}}
                                          for g in sorted(unexplained)[:20]]
    rep["gfs"]["n_repaired_in_both"] = len(repaired & set(old) & set(new))
    rep["gfs"]["n_repaired_unchanged"] = len((repaired & set(old) & set(new)) - changed)
    rep["n_new_pgfams_in_maps"] = len(new_pgfams)
    rep["ok"] = not unexplained
    if a.genomes_out:
        json.dump(sorted(explained_new_products), open(a.genomes_out, "w"))
    json.dump(rep, open(a.out, "w"), indent=1)
    short = {"maps": {m: {k: v for k, v in d.items() if k != "changed_examples"} for m, d in rep["maps"].items()},
             "vocab": {k: (v if not isinstance(v, list) else (len(v), v[:10])) for k, v in rep["vocab"].items()},
             "gfs": {k: v for k, v in rep["gfs"].items() if k != "unexplained_examples"}, "ok": rep["ok"]}
    print(json.dumps(short, indent=1))
    return 0 if rep["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
