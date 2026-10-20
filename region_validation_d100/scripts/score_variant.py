#!/usr/bin/env python
"""score_variant.py — score a VARIANT include-self selection with the *unmodified*
metric code in score_taxacc.py / score_genecap.py.

It does this without editing those scorers: it (1) repoints ``PATHS.selection_out``
at the variant selection JSON and suffixes every output path, and (2) optionally
monkeypatches ``select_references.select_representatives`` so the in-memory
exclude-self RE-selection also uses the variant knobs (so the generalization
scenario reflects the same selection policy as the include scenario).  The include
scenario is read straight from the pre-built variant JSON (no re-selection).

    score_variant.py --selection .../selection_includeall_band.json \
                     --suffix _ib --reselect-knobs '{"select_all":true,"select_all_window":0.005}' \
                     --which both

Interpreter: /home/freiburger/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, "/home/freiburger/Documents/prFBA")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common as C  # noqa: E402
from _common import PATHS  # noqa: E402
import select_references as SR  # noqa: E402

# output attributes to suffix (keep PATHS.selection_out / hits_json / truth_json / pgfam_cache as-is)
TAXACC_OUTS = ["taxacc_jsonl", "concordance_long_csv", "asv_region_concordance_csv",
               "self_recovery_csv", "confusion_genus_csv", "ambiguity_csv", "taxacc_summary_json"]
GENECAP_OUTS = ["genecap_cells_csv", "gene_capture_long_csv", "genecap_summary_json"]


def suffix_paths(attrs, suffix):
    for a in attrs:
        base = getattr(PATHS, a)
        root, ext = os.path.splitext(base)
        setattr(PATHS, a, root + suffix + ext)


def patch_reselect(reselect_knobs):
    """Inject variant knobs into every in-memory exclude-self re-selection."""
    orig = SR.select_representatives

    def patched(asv_record, *, knobs=None, gene_provider=None, gene_count_fn=None):
        merged = {**(knobs or {}), **reselect_knobs}
        return orig(asv_record, knobs=merged, gene_provider=gene_provider, gene_count_fn=gene_count_fn)

    SR.select_representatives = patched


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", required=True, help="variant include-self selection JSON")
    ap.add_argument("--suffix", required=True, help="appended to every output filename")
    ap.add_argument("--reselect-knobs", default="", help="JSON knobs merged into the exclude-self re-selection")
    ap.add_argument("--which", choices=["taxacc", "genecap", "both"], default="both")
    args = ap.parse_args()

    PATHS.selection_out = args.selection
    if args.reselect_knobs:
        patch_reselect(json.loads(args.reselect_knobs))

    if args.which in ("taxacc", "both"):
        suffix_paths(TAXACC_OUTS, args.suffix)
        import score_taxacc as TA
        print(f"[score_variant] taxacc on {args.selection} (suffix {args.suffix})")
        TA.run(force=True)
    if args.which in ("genecap", "both"):
        suffix_paths(GENECAP_OUTS, args.suffix)
        import score_genecap as GC
        print(f"[score_variant] genecap on {args.selection} (suffix {args.suffix})")
        print(GC.run(force=True))


if __name__ == "__main__":
    main()
