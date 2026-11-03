"""Phase 6: extra probes -- per-cell (macro) AUROC on the union support, anchor-family ==
source-family share, 100-hit truncation diagnostics."""
import collections, json, os, pickle
import numpy as np
import vctx

cells = vctx.load_cells("1.0")
nC = len(cells)
nreps = cells.n_reps.astype(int).to_numpy()
with open(os.path.join(vctx.CACHE, "passB.pkl"), "rb") as fh:
    PB = pickle.load(fh)
comps, UQ, UY, UC = PB["comps"], PB["UQ"], PB["UY"], PB["UC"]
out = collections.OrderedDict()

def per_cell_auroc(q, y, c):
    order = np.lexsort((q, c))
    qs, ys, cs = q[order], y[order], c[order]
    n = len(qs)
    newcell = np.r_[True, cs[1:] != cs[:-1]]
    cell_start = np.flatnonzero(newcell)
    cstart_of = np.repeat(cell_start, np.diff(np.r_[cell_start, n]))
    newblk = newcell | np.r_[True, qs[1:] != qs[:-1]]
    bstart = np.flatnonzero(newblk)
    bend = np.r_[bstart[1:], n]
    pos_in_cell_start = bstart - cstart_of[bstart]
    pos_in_cell_end = bend - cstart_of[bstart]
    avg = (pos_in_cell_start + pos_in_cell_end + 1) / 2.0
    ranks = np.repeat(avg, bend - bstart)
    ucells, inv = np.unique(cs, return_inverse=True)
    P = np.bincount(inv, weights=ys.astype(float))
    Nn = np.bincount(inv) - P
    R = np.bincount(inv, weights=ranks * ys)
    with np.errstate(all="ignore"):
        a = (R - P * (P + 1) / 2) / (P * Nn)
    return ucells, a

res = {}
for nm in ("E0", "E1", "E3", "E4(E0)", "E4(E3)", "PREV_foldref", "CLADE_family_only"):
    c = comps.index(nm)
    uc, a = per_cell_auroc(UQ[c].astype(np.float64), UY, UC)
    ok = np.isfinite(a)
    r = {}
    for sn, m in (("all", np.ones(nC, bool)), ("n1", nreps == 1), ("multi", nreps >= 2)):
        mm = m[uc] & ok
        r[sn] = {"mean_per_cell_auroc": float(a[mm].mean()), "cells": int(mm.sum())}
    res[nm] = r
out["per_cell_auroc_union_support"] = res

# anchor family == source family (the family prior can include the source's own family)
C = vctx.Ctx(load_truth=False)
same_fam = 0; same_order = 0; nn = 0
reps0 = []
for f in range(5):
    st = pickle.load(open(os.path.join(vctx.CACHE, f"stream_f{f}.pkl"), "rb"))
    sub = cells[cells.fold == f]
    for s, m in zip(sub.src_genome_id, sub.md5):
        a = st["reps"][m]["reps"][0]
        same_fam += int(C.panel[a]["family"] == C.panel[s]["family"])
        same_order += int(C.panel[a]["order"] == C.panel[s]["order"])
        nn += 1
out["anchor_family_eq_source_family_share"] = same_fam / nn
out["anchor_order_eq_source_order_share"] = same_order / nn

# 100-hit truncation diagnostics (records)
recs = 0; sat = 0; lt100 = 0
for f in range(5):
    st = pickle.load(open(os.path.join(vctx.CACHE, f"stream_f{f}.pkl"), "rb"))
    for m, r in st["recs"].items():
        recs += 1
        lq = r["last_cov_qualified_identity_100"]
        sat += int(r["n_hits"] >= 100 and lq is not None and lq >= vctx.FLOOR)
        lt100 += int(len(r["L100"]) < 100)
out["records"] = recs
out["records_last_cov_qualified_hit_in_first100_ge_floor"] = sat
out["share"] = sat / recs
out["records_with_lt100_survivors_in_first100"] = lt100
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(vctx.OUT, "p6_extra.json"), "w"), indent=1)
