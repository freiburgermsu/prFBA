#!/usr/bin/env python
"""Sample self-exclusion abstentions caused by retention depth, for a deeper re-run.

Under self-exclusion, 94.6% of abstentions are amplicons whose ENTIRE retained
top-20 was deleted by the filter (every hit was the source organism, its taxon, or
>= 0.987 identity) -- the selector saw an empty record rather than a record with no
qualifying reference. The 500-deep edlib shortlist is computed during alignment but
only the top 20 are persisted (``align_hits.py --k2 20``), so the genus/family-level
relatives that would survive the filter were never written out.

This writes a region-stratified sample of those amplicons as a FASTA keyed by
amplicon md5 (the alignment stage's unit of work), to be re-aligned with --k2 500.
"""
import argparse
import json
import os
import random
import sys

import orjson

HERE = os.path.dirname(os.path.abspath(__file__))
RV = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(RV, "scripts"))
import _common as C  # noqa: E402

PATHS = C.PATHS


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selection", default=os.path.join(HERE, "..", "data",
                                                        "prfba_selection.json"))
    ap.add_argument("--outdir", default=os.path.join(HERE, "..", "data", "depth_probe"))
    ap.add_argument("--per-region", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1729)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sel = orjson.loads(open(args.selection, "rb").read())
    truth = orjson.loads(open(PATHS.truth_json, "rb").read())
    fanout = orjson.loads(open(PATHS.amplicon_fanout, "rb").read())
    md5_of = {k: md5 for md5, keys in fanout.items() for k in keys}

    by_region = {}
    for key, rec in sel.items():
        if key == "_meta":
            continue
        exc = rec["exclude_0.987"]
        if exc["n_selected"] == 0 and exc.get("best_identity") in (0.0, None):
            by_region.setdefault(truth[key]["region"], []).append(key)

    rng = random.Random(args.seed)
    picked = []
    for region in sorted(by_region):
        members = sorted(by_region[region])
        picked.extend(rng.sample(members, min(args.per_region, len(members))))
        print(f"[sample] {region:14} {len(members):6} wiped -> "
              f"{min(args.per_region, len(members))} sampled")

    # the aligner works on unique amplicon sequences, so collapse to md5s
    seqs, name = {}, None
    with open(PATHS.unique_fasta) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
            elif line:
                seqs[name] = seqs.get(name, "") + line
    md5s = sorted({md5_of[k] for k in picked})

    fasta = os.path.join(args.outdir, "wiped_sample.fna")
    with open(fasta, "w") as fh:
        for md5 in md5s:
            fh.write(f">{md5}\n{seqs[md5]}\n")
    with open(os.path.join(args.outdir, "wiped_sample.keys.json"), "w") as fh:
        json.dump({"amplicons": sorted(picked), "md5s": md5s,
                   "n_wiped_total": sum(len(v) for v in by_region.values())}, fh)
    print(f"[sample] {len(picked)} amplicons -> {len(md5s)} unique sequences -> {fasta}")


if __name__ == "__main__":
    main()
