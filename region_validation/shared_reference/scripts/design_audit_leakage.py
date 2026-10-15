#!/usr/bin/env python
"""Promoted design-audit script: do byte-identical amplicons imply >=0.987 full-16S id?

Promoted from ``/tmp/leak2.py`` before /tmp was cleared.  It is the measurement behind
SPEC.md decision D1: 27.4% of panel genome pairs that share a byte-identical in-silico
amplicon fall BELOW 0.987 full-16S overlap identity, so a full-16S exclusion rule alone
would leave an exact-amplicon twin in the reference.  That is why the holdout graph
carries the shared-md5 edge (E2) and the amplicon-level identity edge (E4), asserted by
controls C6a/C6b.

Run with no arguments to reproduce the audit printout (needs ``edlib`` and ``ijson``).
The only change from the /tmp original is structural: the body moved into ``main()`` and
``ident_overlap`` became importable.
"""
import argparse
import collections
import json
import os
import random
import statistics as st

REPO = "/home/freiburger/Documents/prFBA"
FANOUT = os.path.join(REPO, "region_validation/data/amplicon_md5_fanout.json")
SELECTION = os.path.join(REPO, "region_validation/data/selection_10k.json")
MD5_SEQ = "/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_seq.json"


def ident_overlap(a, b):
    """Best glocal (free end-gap) identity of two full-16S sequences, either direction."""
    import edlib
    best = 0.0
    for q, t in ((a, b), (b, a)):
        r = edlib.align(q, t, mode="HW", task="path")
        aln = edlib.getNiceAlignment(r, q, t)
        L = len(aln["matched_aligned"])
        m = sum(1 for x, y in zip(aln["query_aligned"], aln["target_aligned"])
                if x == y and x != "-")
        if L:
            best = max(best, m / L)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fanout", default=FANOUT)
    ap.add_argument("--selection", default=SELECTION)
    ap.add_argument("--md5-seq", default=MD5_SEQ)
    ap.add_argument("--sample", type=int, default=400)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    import ijson
    fan = json.load(open(args.fanout))
    sel = json.load(open(args.selection))["genomes"]
    meta = {g["genome_id"]: g for g in sel}
    md5of = {g["genome_id"]: g["best_16s_md5"] for g in sel}
    shared = []
    for m, lst in fan.items():
        gs = sorted(set(i.rsplit("__", 1)[0] for i in lst))
        regs = sorted(set(i.rsplit("__", 1)[1] for i in lst))
        if len(gs) >= 2:
            shared.append((m, gs, regs))
    involved = set(g for _, gs, _ in shared for g in gs)
    print("panel genomes involved in cross-genome identical amplicons:", len(involved))
    random.seed(args.seed)
    samp = random.sample(shared, min(args.sample, len(shared)))
    need = set(md5of[g] for _, gs, _ in samp for g in gs[:4] if g in md5of)
    seqs = {}
    with open(args.md5_seq, "rb") as fh:
        for k, v in ijson.kvitems(fh, ""):
            if k in need:
                seqs[k] = v if isinstance(v, str) else (v[0] if v else None)
                if len(seqs) == len(need):
                    break
    rows = []
    for m, gs, regs in samp:
        for i in range(min(len(gs), 4)):
            for j in range(i + 1, min(len(gs), 4)):
                a, b = md5of.get(gs[i]), md5of.get(gs[j])
                if not a or not b or a == b or a not in seqs or b not in seqs:
                    continue
                rows.append((ident_overlap(seqs[a], seqs[b]), regs[0],
                             meta[gs[i]]["genus"] == meta[gs[j]]["genus"],
                             meta[gs[i]]["family"] == meta[gs[j]]["family"],
                             meta[gs[i]]["phylum"] == meta[gs[j]]["phylum"],
                             abs(len(seqs[a]) - len(seqs[b]))))
    v = [r[0] for r in rows]
    print("pairs n=%d  overlap-identity: median %.4f  p10 %.4f  min %.4f"
          % (len(rows), st.median(v), sorted(v)[len(v) // 10], min(v)))
    for tau in (0.987, 0.99, 0.995):
        bl = [r for r in rows if r[0] < tau]
        print("frac below %.3f: %.3f (%d)  same-genus %.2f same-family %.2f "
              "same-phylum %.2f"
              % (tau, len(bl) / len(rows), len(bl),
                 sum(1 for r in bl if r[2]) / max(1, len(bl)),
                 sum(1 for r in bl if r[3]) / max(1, len(bl)),
                 sum(1 for r in bl if r[4]) / max(1, len(bl))))
    bl = [r for r in rows if r[0] < 0.987]
    print("regions of leaking pairs:", collections.Counter(r[1] for r in bl))
    print("min overlap-identity among leaking pairs: %.4f ; length diffs median %d"
          % (min([r[0] for r in bl]), st.median([r[5] for r in bl])))


if __name__ == "__main__":
    main()
