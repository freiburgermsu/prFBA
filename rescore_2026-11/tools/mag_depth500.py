#!/usr/bin/env python
"""Re-derive the depth-500 MAG comparison the manuscript quotes (L943-945, L1009-1013; claims
mag:twopanel_mag_depth500_original, panels:panels.mag_depth500_rerun; ad hoc in paper f8047c5):
"over 49,185 scoreable cells per regime, self-exclusion resolve 0.3071 -> 0.3089 and conditional F1
0.8066 -> 0.8061; 88 cells change outcome (87 resolve, 1 gene-less reference); no selection switches."

    mag_depth500.py --d100 <out_batches dir> --d500 <out_batches_n500 dir> [--out J]

Joins the depth-100 batch cells (restricted to the MAGs of the depth-500 run) with the depth-500 cells
on (mag, regime). scoreable = outcome != truth_empty; resolve = ok / scoreable (exclude regime);
conditional F1 = mean f1 over ok cells (exclude); changed outcome = outcome differs (both regimes);
'gene-less reference' = an ok cell whose precision is empty (the picks carry no PGFam); selection switch =
both depths resolved and (anchor, n_selected) differ.
"""
import argparse, csv, glob, json, os, collections


def load(d):
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "b*", "mag_cells.csv"))):
        for r in csv.DictReader(open(p)):
            out[(r["mag"], r["regime"])] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d100", required=True); ap.add_argument("--d500", required=True); ap.add_argument("--out")
    a = ap.parse_args()
    A, B = load(a.d100), load(a.d500)
    keys = sorted(B)
    miss = [k for k in keys if k not in A]
    res = {"n_rows_d500": len(B), "n_missing_in_d100": len(miss)}
    for reg in ("include", "exclude"):
        ks = [k for k in keys if k[1] == reg and k in A]
        sc = [k for k in ks if B[k]["outcome"] != "truth_empty"]
        sc100 = [k for k in ks if A[k]["outcome"] != "truth_empty"]
        def rate(T, S): return sum(T[k]["outcome"] == "ok" for k in S) / len(S)
        def f1(T, S):
            v = [float(T[k]["f1"]) for k in S if T[k]["outcome"] == "ok"]
            return sum(v) / len(v)
        chg = collections.Counter()
        for k in ks:
            o1, o2 = A[k]["outcome"], B[k]["outcome"]
            if o1 != o2:
                tag = f"{o1}->{o2}"
                if o2 == "ok" and B[k]["precision"] in ("", "None"):
                    tag += " (gene-less reference)"
                chg[tag] += 1
        sw = sum(1 for k in ks if A[k]["outcome"] == "ok" and B[k]["outcome"] == "ok"
                 and (A[k]["anchor"], A[k]["n_selected"]) != (B[k]["anchor"], B[k]["n_selected"]))
        anch = sum(1 for k in ks if A[k]["anchor"] and B[k]["anchor"] and A[k]["anchor"] != B[k]["anchor"])
        res[reg] = {"n_cells": len(ks), "n_scoreable_d500": len(sc), "n_scoreable_d100": len(sc100),
                    "resolve_d100": round(rate(A, sc100), 6), "resolve_d500": round(rate(B, sc), 6),
                    "cond_f1_d100": round(f1(A, sc100), 6), "cond_f1_d500": round(f1(B, sc), 6),
                    "n_outcome_changed": sum(chg.values()), "outcome_changes": dict(chg),
                    "n_selection_switches_both_resolved": sw, "n_anchor_differs_both_selected": anch}
    s = json.dumps(res, indent=1)
    print(s)
    if a.out:
        open(a.out, "w").write(s + "\n")


if __name__ == "__main__":
    main()
