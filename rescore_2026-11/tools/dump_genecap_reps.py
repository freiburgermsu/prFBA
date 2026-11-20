#!/usr/bin/env python
"""Dump the representative genomes behind every gene-capture cell of one tree (reps provider for
attrib_diff.py --table genecap; the cell tables store only n_selected).

Run it THROUGH ctl_run.py, which points `_common.PATHS` at the tree (and at a substitute cache
and/or a moved-aside selection.json when given):

    ctl_run.py --common $R/<T>/scripts [--cache $SNAP] [--redirect PATHS.selection_out=<old sel>] \
        -- dump_genecap_reps.py --out <T>.reps.pkl [--legacy] [--keep-tables DIR]

It imports the tree's own score_genecap, wraps score_one to record `selected_rep_gids` of each
cell's selection record (the 11-13 rescore_repaired.py method), runs score_genecap.run(force=True)
with all three outputs redirected to a temporary directory (nothing in the tree is written), and
pickles {"cells": {ampliconKey: [{"self_mode": m, "_reps": (gid, ...)}, ...]}, "_meta": {...}}.
"""
import argparse
import os
import pickle
import sys
import tempfile
import time


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--keep-tables", help="also keep the scored tables in this directory")
    ap.add_argument("--legacy", action="store_true",
                    help="score with the LEGACY ordered reducer (score_genecap.py --legacy); required "
                         "for the *_legacy trees and the d20 original tree")
    a = ap.parse_args(argv)
    import _common as C
    import score_genecap as S
    S.LEGACY_SELECT = bool(a.legacy)
    PATHS = C.PATHS
    t0 = time.time()
    cells = {}
    orig = S.score_one

    def wrapped(key, self_mode, sel_record, *x, **k):
        cell = orig(key, self_mode, sel_record, *x, **k)
        cells.setdefault(key, []).append(
            {"self_mode": self_mode, "_reps": tuple(S.selected_rep_gids(sel_record or {})),
             "n_selected": cell.get("n_selected"), "outcome": cell.get("outcome")})
        return cell

    S.score_one = wrapped
    tmpdir = a.keep_tables or tempfile.mkdtemp(prefix="genecap_reps_", dir=os.path.dirname(os.path.abspath(a.out)))
    os.makedirs(tmpdir, exist_ok=True)
    saved = (PATHS.genecap_cells_csv, PATHS.gene_capture_long_csv, PATHS.genecap_summary_json)
    PATHS.genecap_cells_csv = os.path.join(tmpdir, "gene_capture_cells.csv")
    PATHS.gene_capture_long_csv = os.path.join(tmpdir, "gene_capture_long.csv")
    PATHS.genecap_summary_json = os.path.join(tmpdir, "gene_capture_by_region.json")
    try:
        msg = S.run(force=True)
    finally:
        if not a.keep_tables:
            for p in (PATHS.genecap_cells_csv, PATHS.gene_capture_long_csv, PATHS.genecap_summary_json):
                if os.path.exists(p):
                    os.remove(p)
            os.rmdir(tmpdir)
    meta = {"tree_scripts": os.path.dirname(os.path.abspath(S.__file__)), "legacy": bool(a.legacy),
            "selection": PATHS.selection_out, "hits": PATHS.hits_json, "pgfam_cache": PATHS.pgfam_cache,
            "would_have_written": list(saved), "n_amplicons": len(cells),
            "seconds": round(time.time() - t0, 1), "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    with open(a.out + ".tmp", "wb") as fh:
        pickle.dump({"cells": cells, "_meta": meta}, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(a.out + ".tmp", a.out)
    print(msg)
    print(f"[dump_genecap_reps] {len(cells)} amplicons -> {a.out} ({meta['seconds']}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
