#!/usr/bin/env python
"""Locate the completeness at which MAG performance falls to 16S performance.

Merges the original completeness grid (100/90/70/50/30%) with the low-completeness
sweep (25/20/15/10/7/5/3/2/1%), all at 0% contamination, and compares each tier
against the amplicon arm PAIRED within organism -- the 16S value is therefore a
constant computed on whichever organisms that tier pairs with, not a single global
number, so the comparison stays matched as the MAG arm's answered set shrinks.
"""
import csv, collections, json, math, os, sys
from scipy.stats import wilcoxon

MV = os.path.expanduser("~/Documents/prFBA/mag_validation")
AMP = os.path.expanduser("~/Documents/prFBA/region_validation/data/gene_capture_cells.csv")
REGIONS = ("FullLength16S", "V4")


def num(x):
    try:
        return float(x)
    except Exception:
        return None


def load_amp():
    out = collections.defaultdict(dict)
    with open(AMP) as fh:
        for r in csv.DictReader(fh):
            if r["self_mode"] != "exclude_0.987" or r["region"] not in REGIONS:
                continue
            out[r["region"]][r["src_genome_id"]] = (r["outcome"], num(r["f1"]), num(r["recall"]))
    return out


def load_mag(paths):
    """completeness -> genome_id -> (resolved, f1). pred_empty counts as unresolved."""
    out = collections.defaultdict(dict)
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh):
                if r["regime"] != "exclude" or float(r["contamination"]) != 0.0:
                    continue
                if r["outcome"] == "truth_empty":
                    continue
                ok = r["outcome"] == "ok" and r["precision"] != ""
                out[float(r["completeness"])][r["genome_id"]] = (ok, num(r["f1"]))
    return out


def main():
    amp = load_amp()
    mag = load_mag(sys.argv[1:])
    rows = []
    for comp in sorted(mag, reverse=True):
        d = mag[comp]
        n_score = len(d)
        n_ok = sum(1 for v in d.values() if v[0])
        row = dict(completeness=comp, n_scoreable=n_score, n_resolved=n_ok,
                   resolve_rate=round(n_ok / n_score, 4) if n_score else None,
                   cond_f1=round(sum(v[1] for v in d.values() if v[0] and v[1] is not None)
                                 / max(1, n_ok), 4))
        for region in REGIONS:
            A = amp[region]
            pairs = [(d[g][1], A[g][1]) for g in d
                     if d[g][0] and g in A and A[g][0] == "ok"
                     and d[g][1] is not None and A[g][1] is not None]
            if len(pairs) >= 10:
                diff = [a - b for a, b in pairs]
                ds = sorted(diff)
                mid = ds[len(ds) // 2] if len(ds) % 2 else 0.5 * (ds[len(ds) // 2 - 1] + ds[len(ds) // 2])
                try:
                    p = float(wilcoxon(diff, zero_method="zsplit").pvalue)
                except Exception:
                    p = None
                row["%s_n" % region] = len(pairs)
                row["%s_mag" % region] = round(sum(a for a, _ in pairs) / len(pairs), 4)
                row["%s_amp" % region] = round(sum(b for _, b in pairs) / len(pairs), 4)
                row["%s_delta" % region] = round(row["%s_mag" % region] - row["%s_amp" % region], 4)
                row["%s_median_delta" % region] = round(mid, 4)
                row["%s_p" % region] = p
                row["%s_mag_wins" % region] = sum(1 for x in diff if x > 1e-9)
                row["%s_amp_wins" % region] = sum(1 for x in diff if x < -1e-9)
            else:
                for k in ("n", "mag", "amp", "delta", "median_delta", "p", "mag_wins", "amp_wins"):
                    row["%s_%s" % (region, k)] = None
        rows.append(row)

    json.dump(rows, open(MV + "/out/crossover.json", "w"), indent=1)
    print("%-7s %8s %9s %9s | %-32s | %-32s" % ("compl", "resolved", "rate", "condF1",
                                                "vs FullLength16S (MAG/16S/delta)",
                                                "vs V4 (MAG/16S/delta)"))
    for r in rows:
        def blk(region):
            if r["%s_mag" % region] is None:
                return "%-32s" % "   (too few pairs)"
            return "%6.4f / %6.4f / %+7.4f  n=%-5d" % (
                r["%s_mag" % region], r["%s_amp" % region], r["%s_delta" % region], r["%s_n" % region])
        print("%6.0f%% %8d %9.4f %9.4f | %s | %s" % (
            r["completeness"] * 100, r["n_resolved"], r["resolve_rate"], r["cond_f1"],
            blk("FullLength16S"), blk("V4")))
    print("\nwrote out/crossover.json")


if __name__ == "__main__":
    main()
