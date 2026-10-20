#!/usr/bin/env python
"""Project every needed BV-BRC genome's PGFam set into EC and ModelSEED reaction space.

Needed genomes = the benchmark's source genomes (truth) + every genome prFBA selected
as a reference under either regime. Both sides of the comparison go through the SAME
projection (``function_maps.json``), so the EC-vocabulary loss of family-level products
cancels between truth and prediction; ``annotation_ceiling.py`` quantifies what the
projection itself costs against PICRUSt2's own annotation of the same genomes.

Outputs ``--out``: {genome_id: {"ec": [...], "rxn_ec": [...], "rxn_role": [...],
                                "n_pgfam": int}}
  ec       ECs parsed from the products of the genome's PGFams
  rxn_ec   ModelSEED reactions reached from those ECs (same map used for PICRUSt2)
  rxn_role ModelSEED reactions reached from PGFam products via template roles
           (prFBA-only route; reported separately, never mixed into the fair comparison)
"""
import argparse
import os
import sys
import time

import orjson

HERE = os.path.dirname(os.path.abspath(__file__))
RV = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(RV, "scripts"))
import _common as C  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maps", default=os.path.join(HERE, "..", "data", "function_maps.json"))
    ap.add_argument("--selection", default=",".join(
        [os.path.join(HERE, "..", "data", "prfba_selection.json"),
         os.path.join(HERE, "..", "data", "prfba_selection_legacy.json")]),
        help="comma-separated selection JSONs whose reference genomes to project")
    ap.add_argument("--out", default=os.path.join(HERE, "..", "data",
                                                  "genome_function_sets.json"))
    args = ap.parse_args()
    t0 = time.time()

    maps = orjson.loads(open(args.maps, "rb").read())
    pgfam_ec = maps["pgfam_ec"]
    pgfam_rxn = maps["pgfam_rxn"]
    ec_rxn = maps["ec_rxn"]
    truth = orjson.loads(open(C.PATHS.truth_json, "rb").read())
    needed = {t["src_genome_id"] for t in truth.values()}
    for path in args.selection.split(","):
        sel = orjson.loads(open(path, "rb").read())
        for key, rec in sel.items():
            if key != "_meta":
                for regime in rec.values():
                    if isinstance(regime, dict):  # regime blocks only, not scalar fields
                        needed.update(regime["reps"])
        del sel
    print(f"[genome_function_sets] {len(needed)} genomes ({time.time()-t0:.0f}s)")

    cache = orjson.loads(open(C.PATHS.pgfam_cache, "rb").read())
    out = {}
    for gid in sorted(needed):
        pgfams = cache.get(gid) or []
        ec = {e for p in pgfams for e in pgfam_ec.get(p, ())}
        rxn_role = {r for p in pgfams for r in pgfam_rxn.get(p, ())}
        rxn_ec = {r for e in ec for r in ec_rxn.get(e, ())}
        out[gid] = {"ec": sorted(ec), "rxn_ec": sorted(rxn_ec),
                    "rxn_role": sorted(rxn_role), "n_pgfam": len(pgfams)}
    del cache

    n_empty = sum(1 for v in out.values() if not v["ec"])
    counts = sorted(len(v["ec"]) for v in out.values())
    out["_meta"] = {
        "n_genomes": len(out), "n_with_empty_ec": n_empty,
        "median_n_ec": counts[len(counts) // 2] if counts else 0,
        "median_n_pgfam": sorted(v["n_pgfam"] for v in out.values()
                                 if isinstance(v, dict))[len(out) // 2],
        "maps": args.maps, "selection": args.selection,
        "pgfam_cache": C.PATHS.pgfam_cache,
    }
    with open(args.out, "wb") as fh:
        fh.write(orjson.dumps(out))
    print(f"[genome_function_sets] -> {args.out}: {orjson.dumps(out['_meta']).decode()} "
          f"({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
