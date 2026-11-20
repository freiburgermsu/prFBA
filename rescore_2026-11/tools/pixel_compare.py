#!/usr/bin/env python
"""Pixel comparison of two PDFs (re-score figure controls: PDF bytes carry creation dates, so
figures are compared as pixels).

    pixel_compare.py A.pdf B.pdf [--dpi 150] [--out-diff diff.png] [--json report.json]

Both PDFs are rasterised with `pdftoppm -r <dpi>` (every page). Exit 0 iff every page has the same
size and zero differing pixels. The report gives, per page, the size, the number of differing
pixels and the bounding box of the difference (in pixels and as page fractions), so a failure can
be attributed to a region of the figure (e.g. one panel title).
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image


def rasterise(pdf, dpi, outdir, tag):
    prefix = os.path.join(outdir, tag)
    subprocess.run(["pdftoppm", "-r", str(dpi), "-png", pdf, prefix], check=True)
    pages = sorted(f for f in os.listdir(outdir) if f.startswith(tag) and f.endswith(".png"))
    return [os.path.join(outdir, f) for f in pages]


def compare(a, b, dpi=150, out_diff=None):
    with tempfile.TemporaryDirectory() as td:
        pa, pb = rasterise(a, dpi, td, "a"), rasterise(b, dpi, td, "b")
        rep = {"a": os.path.abspath(a), "b": os.path.abspath(b), "dpi": dpi,
               "n_pages": [len(pa), len(pb)], "pages": []}
        ok = len(pa) == len(pb)
        for i, (x, y) in enumerate(zip(pa, pb)):
            ia = np.asarray(Image.open(x).convert("RGB"))
            ib = np.asarray(Image.open(y).convert("RGB"))
            pg = {"page": i + 1, "size_a": list(ia.shape[:2]), "size_b": list(ib.shape[:2])}
            if ia.shape != ib.shape:
                pg["same_size"] = False
                ok = False
            else:
                diff = np.any(ia != ib, axis=2)
                n = int(diff.sum())
                pg.update(same_size=True, n_diff_pixels=n, max_abs_diff=int(np.abs(
                    ia.astype(int) - ib.astype(int)).max()))
                if n:
                    ok = False
                    ys, xs = np.nonzero(diff)
                    h, w = diff.shape
                    pg["bbox_px"] = [int(xs.min()), int(ys.min()), int(xs.max()), int(ys.max())]
                    pg["bbox_frac"] = [round(xs.min() / w, 4), round(ys.min() / h, 4),
                                       round(xs.max() / w, 4), round(ys.max() / h, 4)]
                    if out_diff:
                        vis = (ia * 0.25).astype(np.uint8)
                        vis[diff] = [255, 0, 0]
                        Image.fromarray(vis).save(out_diff if len(pa) == 1 else
                                                  out_diff.replace(".png", f"_p{i + 1}.png"))
            rep["pages"].append(pg)
        rep["identical"] = ok
        return ok, rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--out-diff")
    ap.add_argument("--json")
    a = ap.parse_args(argv)
    ok, rep = compare(a.a, a.b, a.dpi, a.out_diff)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(rep, fh, indent=1)
    print(json.dumps(rep))
    print("PIXEL-IDENTICAL" if ok else "PIXELS DIFFER")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
