#!/usr/bin/env python
"""Fetch the BV-BRC metadata the prFBA-vs-PICRUSt2 comparison needs (resumable).

``genomes``  source-genome assembly accessions (to find the 10k benchmark genomes
             inside PICRUSt2's GTDB-based reference, and for the annotation ceiling)
``pgfams``   ``family_product`` for every PGFam carried by a source genome or by any
             prFBA-selected reference genome (the function strings that get mapped
             to EC numbers and ModelSEED roles)

Both use bulk RQL POSTs with the User-Agent of ``prefetch_gene_families.py`` (BV-BRC
returns 403 to the default Python one), at most one request per ``--min-interval``
seconds, and never cache a failed batch as empty: a batch that keeps failing aborts
the run so it can be resumed.

``pgfams`` collects the PGFam ids of the sources and of every reference genome in
``--selection`` (comma-separated). The default is BOTH selection files, the same pair
``genome_function_sets.py`` projects: until 2026-11 the default was the default
selection only, so the legacy-only references' families were never fetched and the
depth-20 ``prfba_*_legacy`` arms were projected without them (re-score D10).
"""
import argparse
import json
import os
import sys
import time
import urllib.request

import orjson

HERE = os.path.dirname(os.path.abspath(__file__))
RV = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(RV, "scripts"))
import _common as C  # noqa: E402

API = "https://www.bv-brc.org/api"
# The User-Agent of prefetch_gene_families.USER_AGENT (BV-BRC 403s the bare urllib one).
UA = "prFBA/1.0 (academic; metabolic-model benchmarking)"
MIN_INTERVAL = 1.0      # seconds between two requests (rate limit; --min-interval)
_LAST = [0.0]           # monotonic time of the last request


def rql_post(collection, query, retries=5, timeout=300):
    req = urllib.request.Request(
        f"{API}/{collection}/", data=query.encode(),
        headers={"User-Agent": UA, "Accept": "application/json",
                 "Content-Type": "application/rqlquery+x-www-form-urlencoded"})
    for attempt in range(retries):
        wait = MIN_INTERVAL - (time.monotonic() - _LAST[0])
        if wait > 0:
            time.sleep(wait)
        _LAST[0] = time.monotonic()
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.load(resp)
        except Exception as exc:  # noqa: BLE001 -- retried, then re-raised
            if attempt == retries - 1:
                raise
            wait = 2 ** (attempt + 1)
            print(f"[fetch] {collection} batch failed ({exc}); retry in {wait}s", flush=True)
            time.sleep(wait)


def fetch_resumable(collection, ids, id_field, fields, out_path, batch):
    done = orjson.loads(open(out_path, "rb").read()) if os.path.exists(out_path) else {}
    todo = sorted(set(ids) - set(done))
    print(f"[fetch] {collection}: {len(ids)} ids, cached={len(ids) - len(todo)}, "
          f"to_fetch={len(todo)}", flush=True)
    t0 = time.time()
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        q = (f"in({id_field},({','.join(chunk)}))&select({','.join([id_field, *fields])})"
             f"&limit({len(chunk) + 10})")
        rows = rql_post(collection, q)
        for r in rows:
            done[str(r[id_field])] = {f: r.get(f) for f in fields}
        for missing in set(chunk) - {str(r[id_field]) for r in rows}:
            done[missing] = None  # the API answered and has no record for this id
        if (i // batch) % 10 == 0 or i + batch >= len(todo):
            tmp = out_path + ".tmp"
            with open(tmp, "wb") as fh:
                fh.write(orjson.dumps(done))
            os.replace(tmp, out_path)
            print(f"[fetch] {collection}: {min(i + batch, len(todo))}/{len(todo)} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    return done


def selected_reps(paths):
    """Every reference genome any regime of any selection file in `paths` selected."""
    reps = set()
    for path in paths:
        sel = orjson.loads(open(path, "rb").read())
        for key, rec in sel.items():
            if key != "_meta":
                for regime in rec.values():
                    if isinstance(regime, dict):  # regime blocks only, not scalar fields
                        reps.update(regime["reps"])
        del sel
    return reps


def main():
    global MIN_INTERVAL
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["genomes", "pgfams"])
    ap.add_argument("--selection", default=",".join(
        [os.path.join(HERE, "..", "data", "prfba_selection.json"),
         os.path.join(HERE, "..", "data", "prfba_selection_legacy.json")]),
        help="comma-separated selection JSONs whose reference genomes' PGFams to fetch "
             "(default: both, as genome_function_sets.py projects both)")
    ap.add_argument("--outdir", default=os.path.join(HERE, "..", "data"))
    ap.add_argument("--min-interval", type=float, default=MIN_INTERVAL,
                    help="minimum seconds between two BV-BRC requests (default %(default)s)")
    args = ap.parse_args()
    MIN_INTERVAL = args.min_interval
    truth = orjson.loads(open(C.PATHS.truth_json, "rb").read())
    sources = sorted({t["src_genome_id"] for t in truth.values()})

    if args.what == "genomes":
        fetch_resumable("genome", sources, "genome_id",
                        ["assembly_accession", "genome_name", "taxon_id", "genome_status",
                         "cds"],
                        os.path.join(args.outdir, "bvbrc_source_genomes.json"), batch=1000)
        return

    genomes = set(sources) | selected_reps(args.selection.split(","))
    cache = orjson.loads(open(C.PATHS.pgfam_cache, "rb").read())
    pgfams = set()
    for g in genomes:
        pgfams.update(cache.get(g, ()))
    del cache
    print(f"[fetch] {len(genomes)} genomes carry {len(pgfams)} distinct PGFams", flush=True)
    fetch_resumable("protein_family_ref", sorted(pgfams), "family_id", ["family_product"],
                    os.path.join(args.outdir, "pgfam_products.json"), batch=5000)


if __name__ == "__main__":
    main()
