#!/usr/bin/env python
"""Paired MAG-vs-amplicon comparison on the same source organisms.

The amplicon arm's per-organism gene-capture cells
(region_validation/data/gene_capture_cells.csv) and the MAG arm's cells
(mag_cells.csv) share a src_genome_id and an identical PGFam truth set (verified:
truth_size matches per genome), so every comparison below is PAIRED within
organism -- which is the whole point of reusing the same 10,000-genome panel.

Regime mapping: amplicon `include` <-> MAG `include`;
amplicon `exclude_0.987` <-> MAG `exclude` (tau rescaled from 0.987 16S identity
to 95 ANI, the species boundary in each metric's own units).

Two denominators, matching the amplicon arm:
  conditional -- organisms where BOTH arms resolved (outcome ok). Answers "when
                 both answer, which is better".
  population  -- organisms scoreable in both (ok or abstained), an abstention
                 scoring recall/F1 = 0. Answers "which recovers more metabolism".
truth_empty is excluded from both as undefined, never as a miss.
"""
import argparse, collections, csv, json, os, sys
from scipy.stats import wilcoxon

AMP = os.path.expanduser("~/Documents/prFBA/region_validation/data/gene_capture_cells.csv")
REGIME = {"include": "include", "exclude_0.987": "exclude"}


def load_amp():
    """(genome, region, regime) -> dict(f1, recall, precision, outcome)."""
    out = {}
    with open(AMP) as fh:
        for r in csv.DictReader(fh):
            rg = REGIME.get(r["self_mode"])
            if rg is None:
                continue
            def num(x):
                try: return float(x)
                except Exception: return None
            out[(r["src_genome_id"], r["region"], rg)] = dict(
                f1=num(r["f1"]), recall=num(r["recall"]), precision=num(r["precision"]),
                outcome=r["outcome"])
    return out


def load_mag(path):
    """(genome, tier, regime) -> dict(...)"""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            tier = "c%g_x%g" % (float(r["completeness"]), float(r["contamination"]))
            def num(x):
                try: return float(x)
                except Exception: return None
            out[(r["genome_id"], tier, r["regime"])] = dict(
                f1=num(r["f1"]), recall=num(r["recall"]), precision=num(r["precision"]),
                outcome=r["outcome"])
    return out


def paired(mag, amp, tier, region, regime, metric):
    """Return (cond_pairs, pop_pairs) as lists of (mag_value, amp_value)."""
    cond, pop = [], []
    for (g, t, rg), mv in mag.items():
        if t != tier or rg != regime:
            continue
        av = amp.get((g, region, regime))
        if av is None:
            continue
        if mv["outcome"] == "truth_empty" or av["outcome"] == "truth_empty":
            continue
        # amplicon side already excludes its own pred_empty (outcome != 'ok'); apply the
        # same rule to the MAG side, where pred_empty is mis-labelled 'ok' with no precision
        mag_ok = mv["outcome"] == "ok" and mv["precision"] is not None
        if mag_ok and av["outcome"] == "ok":
            a, b = mv[metric], av[metric]
            if a is not None and b is not None:
                cond.append((a, b))
        mz = mv[metric] if mag_ok else 0.0
        az = av[metric] if av["outcome"] == "ok" else 0.0
        if mz is not None and az is not None:
            pop.append((mz, az))
    return cond, pop


def stat(pairs):
    if len(pairs) < 10:
        return dict(n=len(pairs), mag=None, amp=None, median_delta=None, p=None,
                    mag_wins=None, amp_wins=None)
    d = [a - b for a, b in pairs]
    ds = sorted(d)
    mid = (ds[len(ds)//2] if len(ds) % 2 else 0.5*(ds[len(ds)//2 - 1] + ds[len(ds)//2]))
    try:
        p = float(wilcoxon(d, zero_method="zsplit").pvalue)
        if p == 0.0:
            p = None   # underflowed double precision; report as "< 1e-300", not as zero
    except Exception:
        p = None
    return dict(n=len(pairs),
                mag=round(sum(a for a, _ in pairs) / len(pairs), 4),
                amp=round(sum(b for _, b in pairs) / len(pairs), 4),
                median_delta=round(mid, 4), p=p,
                mag_wins=sum(1 for x in d if x > 1e-9),
                amp_wins=sum(1 for x in d if x < -1e-9))


def main():
    global AMP
    ap = argparse.ArgumentParser()
    ap.add_argument("--mag-cells", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--amp", default=AMP,
                    help="amplicon arm gene_capture_cells.csv (default: the depth-20 original "
                         "panel, region_validation/data/gene_capture_cells.csv)")
    a = ap.parse_args()
    AMP = a.amp
    print("amplicon arm: %s" % AMP, flush=True)
    amp, mag = load_amp(), load_mag(a.mag_cells)
    print("amplicon cells=%d  mag cells=%d" % (len(amp), len(mag)), flush=True)
    tiers = sorted({t for _, t, _ in mag}, key=lambda s: -float(s.split("_")[0][1:]))
    regions = sorted({r for _, r, _ in amp})
    res = {}
    for regime in ("include", "exclude"):
        for tier in tiers:
            for region in regions:
                for metric in ("f1", "recall"):
                    c, p = paired(mag, amp, tier, region, regime, metric)
                    res["%s|%s|%s|%s" % (regime, tier, region, metric)] = dict(
                        conditional=stat(c), population=stat(p))
    json.dump(res, open(a.out, "w"), indent=1)
    print("wrote", a.out, flush=True)

    # headline: MAG tier vs FullLength16S, population gene-capture F1
    print("\n=== population gene-capture F1: MAG tier vs FullLength16S (paired) ===")
    print("%-10s %-9s %6s %8s %8s %9s %7s %7s" % ("regime", "tier", "n", "MAG", "16S", "medianD", "MAGw", "16Sw"))
    for regime in ("include", "exclude"):
        for tier in tiers:
            s = res["%s|%s|FullLength16S|f1" % (regime, tier)]["population"]
            if s["n"] is None or not s["n"]:
                continue
            print("%-10s %-9s %6d %8s %8s %9s %7s %7s" % (
                regime, tier, s["n"], s["mag"], s["amp"], s["median_delta"],
                s["mag_wins"], s["amp_wins"]))


if __name__ == "__main__":
    main()
