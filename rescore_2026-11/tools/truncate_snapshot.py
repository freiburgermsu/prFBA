#!/usr/bin/env python
"""Write the first N entries of a PGFam cache file (file order = insertion order) to a new file:
the cache as it stood when it held N genomes (all later keys were appended after that point).
Used to emulate the cache of a committed run whose log records the cache size (re-score R2 controls).

    truncate_snapshot.py --src SNAP --n N --out OUT
"""
import argparse
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_snapshot import ENTRY  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    n, h = 0, hashlib.sha256()
    with open(a.src, "rb") as f, open(a.out + ".tmp", "wb") as out:
        assert f.read(1) == b"{"
        out.write(b"{")
        carry, done = b"", False
        while not done:
            buf = f.read(256 << 20)
            data = carry + buf
            last = 0
            for m in ENTRY.finditer(data):
                out.write((b", " if n else b"") + data[m.start():m.end()])
                h.update(m.group(1) + b"\0")
                n += 1
                last = m.end()
                if n == a.n:
                    done = True
                    break
            carry = data[last:]
            if not buf:
                break
        out.write(b"}")
    if n != a.n:
        raise SystemExit(f"only {n} entries in {a.src}")
    os.replace(a.out + ".tmp", a.out)
    json.dump({"src": a.src, "n": n, "keys_sha256": h.hexdigest(), "out": a.out},
              open(a.out + ".stamp.json", "w"), indent=1)
    print(f"[truncate] {n} entries -> {a.out}")


if __name__ == "__main__":
    main()
