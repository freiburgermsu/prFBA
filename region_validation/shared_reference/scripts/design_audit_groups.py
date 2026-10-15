#!/usr/bin/env python
"""Promoted design-audit script: is the holdout-group graph 10-foldable?

Promoted from ``/tmp/groups.py`` before /tmp was cleared.  It is the feasibility
measurement behind SPEC.md decisions D2 and Section 3.4: on the free edges alone
(genus UNION shared-amplicon-md5) the 10,000-genome panel breaks into 2,986 connected
components whose largest holds 1,086 genomes (10.9% of the panel), so greedy
size-descending bin packing yields ten usable folds -- whereas the family-level graph
has a 2,879-genome (28.8%) giant component and CANNOT be 10-folded, which is why the
leave-family-out ladder uses 5 blocks instead.

The only change from the /tmp original is structural: union-find and the packer became
importable functions (``plan_folds.py`` may reuse them) and the body moved to ``main()``.
"""
import argparse
import collections
import json
import os

REPO = "/home/freiburger/Documents/prFBA"
FANOUT = os.path.join(REPO, "region_validation/data/amplicon_md5_fanout.json")
SELECTION = os.path.join(REPO, "region_validation/data/selection_10k.json")


class UnionFind:
    """Path-halving union-find over arbitrary hashable ids."""

    def __init__(self, items):
        self.par = {x: x for x in items}

    def find(self, x):
        par = self.par
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.par[ra] = rb

    def components(self):
        out = collections.defaultdict(list)
        for x in self.par:
            out[self.find(x)].append(x)
        return {k: sorted(v) for k, v in out.items()}


def greedy_pack(sizes, k):
    """Size-descending greedy bin packing -> ``{component_key: bin_index}``, bin loads.

    Deterministic: components are sorted by (-size, key) and each goes to the
    currently smallest bin, ties broken by the lowest bin index.
    """
    loads = [0] * k
    assign = {}
    for key, sz in sorted(sizes.items(), key=lambda kv: (-kv[1], str(kv[0]))):
        i = min(range(k), key=lambda j: (loads[j], j))
        assign[key] = i
        loads[i] += sz
    return assign, loads


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fanout", default=FANOUT)
    ap.add_argument("--selection", default=SELECTION)
    ap.add_argument("-k", type=int, default=10)
    args = ap.parse_args()
    fan = json.load(open(args.fanout))
    sel = json.load(open(args.selection))["genomes"]
    ids = [g["genome_id"] for g in sel]
    for rank in ("genus", "family"):
        uf = UnionFind(ids)
        bylvl = collections.defaultdict(list)
        for g in sel:
            bylvl[g[rank]].append(g["genome_id"])
        for v in bylvl.values():
            for x in v[1:]:
                uf.union(v[0], x)
        n_md5 = 0
        for _m, lst in fan.items():
            gs = [g for g in sorted(set(i.rsplit("__", 1)[0] for i in lst))
                  if g in uf.par]
            if len(gs) >= 2:
                n_md5 += 1
                for x in gs[1:]:
                    uf.union(gs[0], x)
        comp = {k: len(v) for k, v in uf.components().items()}
        sizes = sorted(comp.values(), reverse=True)
        print("%s + shared-amplicon-md5 (%d multi-genome md5s): %d components; "
              "largest 10 %s" % (rank, n_md5, len(comp), sizes[:10]))
        print("  singletons %d, <=5 %d, >50 %d; top component = %.2f%% of panel"
              % (sum(1 for s in sizes if s == 1), sum(1 for s in sizes if s <= 5),
                 sum(1 for s in sizes if s > 50), 100 * sizes[0] / len(ids)))
        _assign, loads = greedy_pack(comp, args.k)
        print("  greedy %d-fold sizes: %s" % (args.k, loads))
        limit = 1.5 * len(ids) / args.k
        print("  largest component %d vs abort limit %.0f -> %s"
              % (sizes[0], limit, "FOLDABLE" if sizes[0] <= limit else "NOT FOLDABLE"))


if __name__ == "__main__":
    main()
