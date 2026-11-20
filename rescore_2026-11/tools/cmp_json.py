#!/usr/bin/env python
"""Exact JSON comparison with ignored top-level/nested keys.  cmp_json.py A B [--ignore k1,k2] [--tol 0]
Prints the number of differing leaves and the first ones; --tol gives a relative tolerance for floats
(reported separately). Exit 0 iff identical (within --tol)."""
import argparse, json, sys


def walk(x, y, p, ign, tol, out):
    if isinstance(x, dict) and isinstance(y, dict):
        for k in sorted(set(x) | set(y)):
            if k in ign:
                continue
            if k not in x or k not in y:
                out.append((p + "/" + k, "only one side")); continue
            walk(x[k], y[k], p + "/" + k, ign, tol, out)
    elif isinstance(x, list) and isinstance(y, list):
        if len(x) != len(y):
            out.append((p, "len %d vs %d" % (len(x), len(y)))); return
        for i, (u, v) in enumerate(zip(x, y)):
            walk(u, v, "%s[%d]" % (p, i), ign, tol, out)
    elif x != y:
        if tol and isinstance(x, (int, float)) and isinstance(y, (int, float)) and abs(x - y) <= tol * max(abs(x), abs(y)):
            return
        out.append((p, x, y))


ap = argparse.ArgumentParser(); ap.add_argument("a"); ap.add_argument("b")
ap.add_argument("--ignore", default=""); ap.add_argument("--tol", type=float, default=0.0)
a = ap.parse_args()
out = []
walk(json.load(open(a.a)), json.load(open(a.b)), "", set(filter(None, a.ignore.split(","))), a.tol, out)
print(json.dumps({"n_diff": len(out), "first": out[:15]}, indent=1, default=str))
print("RESULT:", "IDENTICAL" if not out else "DIFFERENT")
sys.exit(0 if not out else 1)
