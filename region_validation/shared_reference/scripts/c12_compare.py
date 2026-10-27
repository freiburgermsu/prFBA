#!/usr/bin/env python
"""Before/after comparison for SPEC amendment A-C12: the phenotype rows of a pass run on the
original templates (the pre-registered run) against the same pass on templates 21a8084 with
the gapfiller's reversals applied. Rows are matched on (genome, arm, pipeline).

Reports, overall and split into the 20 genomes that held C12's 49 gated rows and the rest:
  gate      gated rows by arm, before and after; fixed (gated -> ok) and new (ok -> gated) leaks
  growth    non-growing gapfilled models; growth ratio after/before where both are ok
  reversal  rows whose gapfill proposed / applied reversed reactions
  atp       rows whose binary core-ATP vector changed; per-arm mean ATP Hamming to truth
  readouts  per-arm means of the vs-truth readouts, before and after

    python c12_compare.py --before data/phenotype_gf/per_genome.jsonl \
        --after data/phenotype_gf_c12fix/per_genome.jsonl --out data/c12fix_compare_gf.json
"""
import argparse
import csv
import json
import os
import statistics as st
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)
MAIN_ARMS = ("truth", "prfba_exact_tie", "picrust2", "base_nn1", "base_family_donor",
             "base_constant")
ATP_T = 1e-6


def rows(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                r = json.loads(line)
                out[(r["src_genome_id"], r["arm"], r["pipeline"])] = r
    return out


def gated_genomes():
    with open(os.path.join(SR, "EGC_C12_gated_models.csv")) as fh:
        return {r["genome"] for r in csv.DictReader(fh)}


def binvec(r):
    v = ((r.get("atp") or {}).get("vector")) or {}
    return {k: float(x) > ATP_T for k, x in v.items()}


def mean(xs):
    xs = [x for x in xs if x is not None]
    return (round(sum(xs) / len(xs), 4), len(xs)) if xs else (None, 0)


def compare(B, A, keys):
    res = {"n_rows": len(keys)}
    g = {"before": Counter(), "after": Counter()}
    fixed, new = Counter(), Counter()
    for k in keys:
        b, a = B[k], A[k]
        arm = k[1]
        bg, ag = b.get("status") == "gated", a.get("status") == "gated"
        g["before"][arm] += bg
        g["after"][arm] += ag
        fixed[arm] += bg and not ag
        new[arm] += ag and not bg
    res["gate"] = {"gated_before": dict(g["before"]), "gated_after": dict(g["after"]),
                   "fixed": dict(fixed), "new_leaks": dict(new),
                   "n_gated_before": sum(g["before"].values()),
                   "n_gated_after": sum(g["after"].values())}
    nb = na = 0
    ratios, n_rev_prop, n_rev_app, n_rows_gf = [], 0, 0, 0
    burden = []
    for k in keys:
        b, a = B[k].get("gapfilled") or {}, A[k].get("gapfilled") or {}
        if not a and not b:
            continue
        n_rows_gf += 1
        gb, ga = float(b.get("growth") or 0.0), float(a.get("growth") or 0.0)
        nb += gb <= 1e-5
        na += ga <= 1e-5
        if B[k].get("status") == "ok" and A[k].get("status") == "ok" and gb > 1e-5:
            ratios.append(ga / gb)
        n_rev_prop += bool(a.get("reversed_proposed"))
        n_rev_app += bool(a.get("reversed_applied"))
        if b.get("n_new_rxn") is not None and a.get("n_new_rxn") is not None:
            burden.append(a["n_new_rxn"] - b["n_new_rxn"])
    res["growth"] = {"n_gapfilled_rows": n_rows_gf, "non_growing_before": nb,
                     "non_growing_after": na,
                     "ratio_after_over_before": (
                         {"median": round(st.median(ratios), 4), "n": len(ratios),
                          "changed_gt10pct": sum(abs(r - 1) > 0.1 for r in ratios)}
                         if ratios else None),
                     "rows_with_reversals_proposed": n_rev_prop,
                     "rows_with_reversals_applied": n_rev_app,
                     "gapfill_burden_delta": ({"mean": round(st.mean(burden), 3),
                                               "median": st.median(burden),
                                               "n_changed": sum(x != 0 for x in burden)}
                                              if burden else None)}
    changed = Counter()
    seen = Counter()
    for k in keys:
        vb, va = binvec(B[k]), binvec(A[k])
        if vb or va:
            seen[k[1]] += 1
            changed[k[1]] += vb != va
    res["atp_vector_changed"] = {arm: {"changed": changed[arm], "n": seen[arm]}
                                 for arm in sorted(seen)}
    per = defaultdict(dict)
    for metric in ("atp_hamming", "growth_log_ratio", "ess_shared_disagree",
                   "rxn_jaccard_draft", "rxn_jaccard_gapfilled", "gapfill_burden_delta"):
        for arm in MAIN_ARMS:
            for pipe in ("controlled", "as_shipped"):
                ks = [k for k in keys if k[1] == arm and k[2] == pipe]
                if not ks:
                    continue
                bb = mean([(B[k].get("vs_truth") or {}).get(metric) for k in ks])
                aa = mean([(A[k].get("vs_truth") or {}).get(metric) for k in ks])
                if bb[1] or aa[1]:
                    per[metric]["%s|%s" % (arm, pipe)] = {"before": bb[0], "n_before": bb[1],
                                                           "after": aa[0], "n_after": aa[1]}
    res["readouts"] = per
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", required=True)
    ap.add_argument("--after", required=True)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    B, A = rows(a.before), rows(a.after)
    common = sorted(set(B) & set(A))
    gg = gated_genomes()
    out = {"before": os.path.abspath(a.before), "after": os.path.abspath(a.after),
           "n_before": len(B), "n_after": len(A), "n_matched": len(common),
           "only_before": len(set(B) - set(A)), "only_after": len(set(A) - set(B)),
           "all": compare(B, A, common),
           "gated_20_genomes": compare(B, A, [k for k in common if k[0] in gg]),
           "other_genomes": compare(B, A, [k for k in common if k[0] not in gg])}
    s = json.dumps(out, indent=1, sort_keys=True)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(s + "\n")
    for part in ("all", "gated_20_genomes", "other_genomes"):
        p = out[part]
        print("== %s: %d matched rows" % (part, p["n_rows"]))
        print("   gate", json.dumps(p["gate"]))
        print("   growth", json.dumps(p["growth"]))


if __name__ == "__main__":
    main()
