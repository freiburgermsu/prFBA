#!/usr/bin/env python
"""build_variant_selection.py — build an include-all ("union") selection over the
full 68k-amplicon hits, in the SAME JSON shape as selection.json, so the existing
scorers (score_taxacc / score_genecap) can grade it unchanged.

The union mode is a pure function of the hit record (no gene_provider needed), so
this is fast and deterministic.

    build_variant_selection.py --window 0.005 --out .../selection_includeall_band.json
    build_variant_selection.py --window 1.0   --out .../selection_includeall_reliable.json

Interpreter: /home/freiburger/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, "/home/freiburger/Documents/prFBA")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import select_references as SR  # noqa: E402
from _common import PATHS  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", default=PATHS.hits_json)
    ap.add_argument("--window", type=float, default=0.005,
                    help="identity window below best to union (0.005=in-band, >=1.0=all reliable)")
    ap.add_argument("--cap", type=int, default=0, help="0=no cap; else max genomes")
    ap.add_argument("--coherence", action="store_true", help="guarded: anchor-family cross-family firewall")
    ap.add_argument("--contam", action="store_true", help="guarded: retain per-tier off-target gene budget")
    ap.add_argument("--species-dedup", action="store_true", help="<=1 representative per species (best copy)")
    ap.add_argument("--unconditional-window", type=float, default=-1.0,
                    help="graduated: protect (keep unconditionally) members within this identity of best")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    knobs = {"select_all": True, "select_all_window": args.window, "select_all_cap": args.cap,
             "select_all_coherence": args.coherence, "select_all_contam": args.contam,
             "select_all_species_dedup": args.species_dedup,
             "select_all_unconditional_window": args.unconditional_window}
    t0 = time.perf_counter()
    hits = json.load(open(args.hits))
    hits = {k: v for k, v in hits.items() if k != "_meta"}
    audit, summary = SR.build_audit(hits, knobs=knobs)
    out = {"_meta": {
        "generated_from": args.hits,
        "scheme": "include-all union mode",
        "gene_data": "n/a (union mode does not consult gene sets; capture scored downstream)",
        "select_all_window": args.window, "select_all_cap": args.cap,
        "knobs": {**SR.DEFAULT_KNOBS, **knobs},
        "summary": summary,
    }, **audit}
    json.dump(out, open(args.out, "w"))
    print(f"[variant] window={args.window} cap={args.cap}: {summary['n_asvs']} amplicons "
          f"-> {args.out} ({time.perf_counter()-t0:.1f}s)")
    print(f"[variant]   mean_handful={summary['mean_handful']} "
          f"handful_sizes={summary['handful_sizes']} abstained={summary['abstained']}")


if __name__ == "__main__":
    main()
