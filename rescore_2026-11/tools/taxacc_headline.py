#!/usr/bin/env python
"""Old -> new taxonomy headline of taxacc_region_summary.json pairs: anchor predictor, domain all,
tier all, Genus and Species correct_population (and coverage) per self mode, 8-region mean and
min-max; plus self-recovery rate if present.  usage: taxacc_headline.py LABEL OLD NEW [...] --out F"""
import json, sys
args = sys.argv[1:]
out_path = None
if "--out" in args:
    i = args.index("--out"); out_path = args[i + 1]; args = args[:i] + args[i + 2:]
rep = {}
def get(d, reg, mode, rank, metric):
    try:
        return d["regions"][reg]["anchor"][mode]["all"]["all"]["per_rank"][rank][metric]
    except (KeyError, TypeError):
        return None
for j in range(0, len(args), 3):
    label, op, np_ = args[j:j + 3]
    o, n = json.load(open(op)), json.load(open(np_))
    regs = list(n["regions"])
    rep[label] = {}
    print(f"== {label}")
    for mode in ("include", "exclude_0.987"):
        for rank in ("Genus", "Species"):
            for metric in ("correct_population", "coverage"):
                ov = [get(o, r, mode, rank, metric) for r in regs]; nv = [get(n, r, mode, rank, metric) for r in regs]
                if any(v is None for v in ov + nv) or any(isinstance(v, dict) for v in ov + nv):
                    ov = [v.get("value") if isinstance(v, dict) else v for v in ov]; nv = [v.get("value") if isinstance(v, dict) else v for v in nv]
                if any(v is None for v in ov + nv):
                    continue
                mo, mn = sum(ov) / len(ov), sum(nv) / len(nv)
                rep[label][f"{mode}|{rank}|{metric}"] = {"old_mean": mo, "new_mean": mn, "old_min": min(ov), "old_max": max(ov),
                                                         "new_min": min(nv), "new_max": max(nv)}
                print(f"  {mode:14s} {rank:8s} {metric:19s} mean {100*mo:.3f} -> {100*mn:.3f}   range {100*min(ov):.2f}-{100*max(ov):.2f} -> {100*min(nv):.2f}-{100*max(nv):.2f}")
if out_path:
    json.dump(rep, open(out_path, "w"), indent=1)
