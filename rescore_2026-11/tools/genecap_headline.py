#!/usr/bin/env python
"""Old -> new headline numbers of gene_capture_by_region.json files (re-score report helper).

    genecap_headline.py --pair LABEL OLD NEW [--pair ...] [--out report.json]

Per pair and self mode: outcome-count totals, the 8-region mean of macro recall over OK cells
(resolved) and over the population (abstain -> 0), and per region the population macro recall.
"""
import argparse
import json
import sys


def summarise(path):
    d = json.load(open(path))
    out = {"select_mode": d["_meta"].get("select_mode"), "outcomes": d["_meta"].get("outcome_counts_total"), "modes": {}}
    regs = d["regions"]
    modes = d["_meta"].get("self_modes") or sorted({m for r in regs.values() for m in r})
    for m in modes:
        res, pop = {}, {}
        for reg, row in regs.items():
            if m not in row:
                continue
            res[reg] = row[m]["macro"]["recall"]["value"]
            pop[reg] = row[m]["macro_pop"]["recall"]["value"]
        out["modes"][m] = {
            "mean_resolved_recall_pct": round(100 * sum(res.values()) / len(res), 3) if res else None,
            "mean_population_recall_pct": round(100 * sum(pop.values()) / len(pop), 3) if pop else None,
            "population_recall_pct": {k: round(100 * v, 3) for k, v in pop.items()},
        }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pair", nargs=3, action="append", metavar=("LABEL", "OLD", "NEW"), required=True)
    ap.add_argument("--out")
    a = ap.parse_args()
    rep = {}
    for label, old, new in a.pair:
        o, n = summarise(old), summarise(new)
        rep[label] = {"old": o, "new": n, "old_path": old, "new_path": new}
        print(f"== {label}  ({o['select_mode']} -> {n['select_mode']})")
        for m in n["modes"]:
            oo, nn = o["modes"].get(m, {}), n["modes"][m]
            print(f"  {m:14s} resolved {oo.get('mean_resolved_recall_pct')} -> {nn['mean_resolved_recall_pct']}   "
                  f"population {oo.get('mean_population_recall_pct')} -> {nn['mean_population_recall_pct']}")
            oc, nc = (o["outcomes"] or {}).get(m, {}), (n["outcomes"] or {}).get(m, {})
            print("  " + " " * 14 + " outcomes " + ", ".join(f"{k} {oc.get(k)}->{nc.get(k)}" for k in nc if oc.get(k) != nc.get(k)))
    if a.out:
        json.dump(rep, open(a.out, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
