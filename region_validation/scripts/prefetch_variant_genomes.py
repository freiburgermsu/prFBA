#!/usr/bin/env python
"""prefetch_variant_genomes.py — fetch BV-BRC PGFam sets for the genomes the
include-all VARIANT selections pick that are NOT already in the shared cache, so
gene-capture precision/recall can be scored without ``pgfam_missing`` poisoning.

The shared cache (PATHS.pgfam_cache) was built for the baseline's candidate set +
the 10k truth genomes; the union variants select many more references.  This
appends only the missing ones (existing entries untouched -> the committed
baseline gene-capture is unaffected).  Threaded, resumable, atomic save.

Interpreter: /home/freiburger/Documents/py_venv/bin/python
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, "/home/freiburger/Documents/prFBA")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import PATHS  # noqa: E402
import prefetch_gene_families as PF  # noqa: E402

SELECTIONS = [
    os.path.join(os.path.dirname(PATHS.selection_out), "selection_includeall_band.json"),
    os.path.join(os.path.dirname(PATHS.selection_out), "selection_includeall_reliable.json"),
]
WORKERS = 12


def main():
    cache = PF.load_cache(PATHS.pgfam_cache)
    have = set(cache)
    want = set()
    for p in SELECTIONS:
        sel = json.load(open(p)); sel.pop("_meta", None)
        for rec in sel.values():
            for s in rec.get("selected") or []:
                want.add(str(s["genome_id"]))
    todo = sorted(want - have)
    print(f"[prefetch-variant] cache={len(have)} want={len(want)} to_fetch={len(todo)}", flush=True)
    if not todo:
        print("[prefetch-variant] nothing to fetch")
        return
    t0 = time.perf_counter()
    done = empty = errors = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(PF.fetch_pgfams, g): g for g in todo}
        for fut in as_completed(futs):
            g = futs[fut]
            try:
                cache[g] = fut.result()
                empty += not cache[g]
            except Exception:
                cache[g] = []
                errors += 1
            done += 1
            if done % 500 == 0:
                PF.save_cache(cache, PATHS.pgfam_cache)
                el = time.perf_counter() - t0
                print(f"[prefetch-variant] {done}/{len(todo)} ({el:.0f}s, {done/max(el,1):.1f}/s) "
                      f"empty={empty} err={errors}", flush=True)
    PF.save_cache(cache, PATHS.pgfam_cache)
    print(f"[prefetch-variant] DONE fetched {len(todo)} in {time.perf_counter()-t0:.0f}s "
          f"-> cache now {len(cache)} genomes (empty={empty} err={errors})", flush=True)


if __name__ == "__main__":
    main()
