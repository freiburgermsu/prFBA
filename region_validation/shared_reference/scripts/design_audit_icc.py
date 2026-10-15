#!/usr/bin/env python
"""Promoted design-audit script: ICC, design effect and effective n of paired ΔF1.

Promoted from ``/tmp/icc.py`` before /tmp was cleared.  The audit measured that the
per-amplicon prFBA-minus-PICRUSt2 F1 difference is strongly clustered -- ICC 0.66 by
genome, 0.48 by genus, 0.35 by family, i.e. design effects 5.0 / 10.3 / 17.5 -- so
67,038 amplicons carry the information of roughly 13,400 / 6,500 / 3,800 independent
observations.  SPEC.md Section 8.4 requires ICC and n_eff to be printed next to every
p-value, and Section 13 requires score_folds.py to import this implementation rather
than re-deriving it.

The only change from the /tmp original is structural: ``icc()`` became a pure function
of a grouping and the script body moved into ``main()``.

Guarantees:
  * one-way random-effects ICC(1) from the ANOVA mean squares, with the group size
    averaged over groups of size >= 2 (singletons carry no within-group information and
    are dropped from the estimate, exactly as in the audit);
  * ``design_effect``/``n_eff`` clamp a negative ICC to 0, so a sub-Poisson grouping can
    never inflate n_eff above n.
"""
import argparse
import collections
import csv
import gzip
import json
import os
import statistics as st

REPO = "/home/freiburger/Documents/prFBA"
CELLS = os.path.join(REPO, "region_validation/picrust2_comparison/data/"
                           "comparison_cells.csv.gz")
SELECTION = os.path.join(REPO, "region_validation/data/selection_10k.json")


def icc_from_groups(groups):
    """One-way ICC(1) over ``{group_key: [values]}``.

    Returns ``{"icc", "n_clusters", "n", "mean_cluster_size", "design_effect",
    "n_eff", "mean", "sd"}``; ``icc`` is None when fewer than two usable groups.
    """
    grp = {k: v for k, v in groups.items() if len(v) >= 2}
    allv = [x for v in grp.values() for x in v]
    k, n = len(grp), len(allv)
    if k < 2 or n <= k:
        flat = [x for v in groups.values() for x in v]
        return {"icc": None, "n_clusters": len(groups), "n": len(flat),
                "mean_cluster_size": (len(flat) / len(groups)) if groups else 0.0,
                "design_effect": 1.0, "n_eff": len(flat),
                "mean": st.mean(flat) if flat else None,
                "sd": st.pstdev(flat) if len(flat) > 1 else 0.0}
    gm = st.mean(allv)
    msb = sum(len(v) * (st.mean(v) - gm) ** 2 for v in grp.values()) / (k - 1)
    msw = sum((x - st.mean(v)) ** 2 for v in grp.values() for x in v) / (n - k)
    mbar = st.mean([len(v) for v in grp.values()])
    den = msb + (mbar - 1) * msw
    val = (msb - msw) / den if den else 0.0
    deff = 1 + (mbar - 1) * max(val, 0.0)
    return {"icc": val, "n_clusters": k, "n": n, "mean_cluster_size": mbar,
            "design_effect": deff, "n_eff": n / deff if deff else n,
            "mean": gm, "sd": st.pstdev(allv)}


def icc(values, keyfn):
    """``{item: value}`` + ``item -> group key`` -> ``icc_from_groups`` block."""
    grp = collections.defaultdict(list)
    for item, v in values.items():
        g = keyfn(item)
        if g:
            grp[g].append(v)
    return icc_from_groups(grp)


def design_effect(icc_value, mean_cluster_size):
    return 1 + (mean_cluster_size - 1) * max(icc_value or 0.0, 0.0)


def n_eff(n, icc_value, mean_cluster_size):
    d = design_effect(icc_value, mean_cluster_size)
    return n / d if d else n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", default=CELLS)
    ap.add_argument("--selection", default=SELECTION)
    ap.add_argument("--arm-a", default="prfba_exclude")
    ap.add_argument("--arm-b", default="picrust2")
    ap.add_argument("--method-field", default="method")
    args = ap.parse_args()

    sel = json.load(open(args.selection))["genomes"]
    meta = {g["genome_id"]: g for g in sel}
    pr, pc = {}, {}
    with gzip.open(args.cells, "rt") as fh:
        for row in csv.DictReader(fh):
            if row["space"] != "ec":
                continue
            m = row[args.method_field]
            if m == args.arm_a:
                pr[row["amplicon"]] = row
            elif m == args.arm_b:
                pc[row["amplicon"]] = row

    def f1pop(row):
        return 0.0 if row["f1"] == "" else float(row["f1"])

    delta, dcond = {}, {}
    for a in pr:
        if a not in pc or int(pr[a]["n_truth"] or 0) == 0:
            continue
        delta[a] = f1pop(pr[a]) - f1pop(pc[a])
        if pr[a]["outcome"] in ("ok", "scored") and pc[a]["outcome"] in ("ok", "scored"):
            dcond[a] = float(pr[a]["f1"]) - float(pc[a]["f1"])

    gid = lambda a: a.rsplit("__", 1)[0]  # noqa: E731
    keys = [("genome", gid)]
    for rank in ("genus", "family", "phylum"):
        keys.append((rank, lambda a, r=rank: meta[gid(a)][r] if gid(a) in meta else None))
    for name, keyfn in keys:
        for label, d in (("population", delta), ("conditional", dcond)):
            b = icc(d, keyfn)
            if b["icc"] is None:
                continue
            print("%-11s %-11s ICC=%.3f  clusters=%5d  n=%6d  mean_m=%.2f  "
                  "design_eff=%.2f  n_eff=%.0f  mean_delta=%+.4f sd=%.3f"
                  % (name, label, b["icc"], b["n_clusters"], b["n"],
                     b["mean_cluster_size"], b["design_effect"], b["n_eff"],
                     b["mean"], b["sd"]))
    g2 = collections.defaultdict(list)
    for a in pr:
        g2[gid(a)].append(pr[a]["outcome"])
    same = collections.Counter(len(set(v)) == 1 for v in g2.values())
    print("genomes where the arm's abstain/answer is identical across all its regions:",
          dict(same))


if __name__ == "__main__":
    main()
