#!/usr/bin/env python
"""
tieset_gene_cache.py -- write a SMALL, NEW PGFam gene cache holding exactly the genomes an
exact-tie union selection can consult, streamed (read-only) out of the shared multi-GB
cache, so the deployment build neither loads ~5 GB of JSON (~28 GB of Python lists, ~40 GB
on the provider's set path, measured by scaling a 500 MB prefix) nor rewrites the shared
file another session may be using.

WHY THIS IS EXACT. Under the shipped exact-tie union (select_references DEFAULT_KNOBS,
select_all=True, window 0.0) the gene cache enters selection ONLY through
gene_count_fn(genome_id) -- the P3 species-dedup tie-break (select_references.py
_gc_tiebreak) and est_genes / gene_set_empty of selected members. The gene_provider itself
is never called on that path (it is used only by the legacy reducer). gene_count_fn is
queried only for the genomes in each ASV's exact-tie pool after the reliability gate and
genome_id dedup -- i.e. the members a provider-free audit marks `selected` or
`species_dedup`. That pool does not depend on gene counts (P3 changes which member of a
species is the representative, never the pool). So a cache restricted to that "tie set",
with each genome's entry copied verbatim, returns the same count (or the same None for a
genome the shared cache lacks) for every query the build makes, and the selection is
identical to one made against the full shared cache. build_d100.py --tieset re-checks this
at run time: it counts gene_count_fn queries outside the tie set and aborts if any occur.

    python deployment_d100/tieset_gene_cache.py \
        --hits  <d100 hits dir>/asv_top20_alignment_hits.json \
        --source-cache bvbrc_cache/genome_gene_families.json \
        --out   <d100 hits dir>/gene_families_tieset.json
    # writes OUT and OUT.tieset.json (the tie-set genome ids + source-cache identity)

Then prefetch_gene_families.py --cache OUT fetches the selected genomes the shared cache
lacked (exactly what it would have fetched into the shared cache).

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import select_references as S  # noqa: E402

_ENTRY = re.compile(rb'"([^"]+)": \[([^\]]*)\]')


def tie_set(hits, knobs=None):
    """Genome ids gene_count_fn can be asked about under the exact-tie union."""
    audit, summary = S.build_audit(hits, knobs=knobs)
    ids = set()
    for r in audit.values():
        for d in r["decisions"]:
            if d["disposition"] == "selected" or d["stage"] == "species_dedup":
                ids.add(str(d["genome_id"]))
    return ids, summary


def stream_extract(path, want, chunk=256 << 20):
    """Copy the entries for `want` out of a single-line {gid: [..], ..} JSON, streaming."""
    out, n_seen = {}, 0
    carry = b""
    with open(path, "rb") as f:
        if f.read(1) != b"{":
            raise SystemExit(f"[abort] {path}: not a JSON object")
        while True:
            buf = f.read(chunk)
            data = carry + buf
            last = 0
            for m in _ENTRY.finditer(data):
                n_seen += 1
                gid = m.group(1).decode()
                if gid in want:
                    out.update(json.loads(b"{" + m.group(0) + b"}"))
                last = m.end()
            carry = data[last:]
            if not buf:
                break
    if carry.strip() not in (b"}", b""):
        raise SystemExit(f"[abort] {path}: {len(carry)} trailing bytes did not parse as entries")
    return out, n_seen


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--source-cache", default="/home/freiburger/Documents/prFBA/bvbrc_cache/genome_gene_families.json")
    ap.add_argument("--out", required=True)
    ap.add_argument("--force", action="store_true", help="overwrite --out if it exists")
    args = ap.parse_args()
    if os.path.exists(args.out) and not args.force:
        sys.exit(f"[abort] {args.out} exists (a prefetch may already have added to it); "
                 f"pass --force to rebuild it from the source cache")
    t0 = time.perf_counter()
    hits = json.load(open(args.hits))
    want, summary = tie_set(hits)
    del hits
    print(f"[tieset] {len(want)} genomes in the exact-tie pools of {summary['n_asvs']} ASVs "
          f"(abstain={summary['abstained']}) ({time.perf_counter()-t0:.0f}s)", flush=True)
    st = os.stat(args.source_cache)
    t1 = time.perf_counter()
    sub, n_seen = stream_extract(args.source_cache, want)
    print(f"[tieset] streamed {n_seen} entries ({st.st_size/1e9:.2f} GB) from {args.source_cache}: "
          f"{len(sub)}/{len(want)} tie-set genomes present ({time.perf_counter()-t1:.0f}s)", flush=True)
    tmp = args.out + ".tmp"
    with open(tmp, "w") as f:
        json.dump(sub, f)
    os.replace(tmp, args.out)
    side = {"hits": os.path.abspath(args.hits), "source_cache": os.path.abspath(args.source_cache),
            "source_cache_size": st.st_size, "source_cache_mtime": st.st_mtime,
            "source_cache_entries": n_seen, "n_tieset": len(want),
            "n_present_in_source": len(sub), "tieset": sorted(want)}
    with open(args.out + ".tieset.json", "w") as f:
        json.dump(side, f)
    print(f"[tieset] wrote {args.out} ({os.path.getsize(args.out)/1e6:.1f} MB) + .tieset.json "
          f"({time.perf_counter()-t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
