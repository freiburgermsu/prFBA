#!/usr/bin/env python
"""Aggregate per-batch MAG cells into the final summary.

Metric definitions match region_validation/scripts/score_taxacc.py and
score_genecap.py exactly:
  * Per-rank cells are three-state. A blank cell in mag_cells.csv is ambiguous on
    its own (it means EITHER truth-undefined OR a no-call), so scoreability is
    recovered authoritatively from benchmark_truth.json: a rank is scoreable for an
    organism iff truth_lineage[rank] is non-null. Unscoreable cells leave the
    denominator entirely and are never counted wrong.
  * `correct` = n_correct / n_called  (CONDITIONAL on the pipeline answering)
  * `correct_population` = n_correct / n_scoreable  (an abstention counts as
    not-correct) -- the anti-inflation denominator, fixed by truth, so exclusion
    abstentions cannot shrink it toward the easy cases.
  * Gene capture conditional = macro (organism-equal) mean over queries that
    resolved; population additionally scores an abstention as recall/F1/Jaccard = 0
    with precision undefined; truth_empty is excluded from both as undefined.
  * Wilson 95% CI for rank proportions; normal-approx mean CI for means of
    proportions.
"""
import argparse, collections, csv, glob, json, math, os, sys
sys.path.insert(0, os.path.expanduser("~/Documents/prFBA/region_validation/scripts"))
from _common import wilson_ci

ROOT = os.path.expanduser("~/Documents/prFBA")
MV = ROOT + "/mag_validation"
SCORED = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


def macro(vals):
    v = [x for x in vals if x is not None and not (isinstance(x, float) and math.isnan(x))]
    if not v:
        return dict(mean=None, ci=[None, None], n=0)
    m = sum(v) / len(v)
    if len(v) < 2:
        return dict(mean=round(m, 4), ci=[None, None], n=1)
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))
    h = 1.959964 * sd / math.sqrt(len(v))
    return dict(mean=round(m, 4), ci=[round(max(0, m - h), 4), round(min(1, m + h), 4)], n=len(v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batches", default=MV + "/out_batches",
                    help="directory of per-batch scorer outputs")
    ap.add_argument("--out", default=MV + "/out/mag_summary_full.json")
    ap.add_argument("--truth", default=ROOT + "/region_validation/data/benchmark_truth.json",
                    help="benchmark_truth.json for THIS panel; scoreability is read from it")
    a = ap.parse_args()
    truth = {}
    for v in json.load(open(a.truth)).values():
        truth.setdefault(v["src_genome_id"], v["truth_lineage"])
    files = sorted(glob.glob(os.path.join(a.batches, "*/mag_cells.csv")))
    print("aggregating %d batch files" % len(files), flush=True)
    cells = collections.defaultdict(list)
    n = 0
    for f in files:
        with open(f) as fh:
            for r in csv.DictReader(fh):
                cells[(float(r["completeness"]), float(r["contamination"]), r["regime"])].append(r)
                n += 1
    print("total cells: %d" % n, flush=True)

    def num(x):
        try: return float(x)
        except Exception: return None

    summary = {}
    for key in sorted(cells, key=lambda k: (k[2], -k[0], k[1])):
        comp, contam, regime = key
        rows = cells[key]
        s = dict(completeness=comp, contamination=contam, regime=regime, n=len(rows))
        for rank in SCORED:
            # scoreable comes from the TRUTH, not from the cell being blank
            sc = [r for r in rows if (truth.get(r["genome_id"]) or {}).get(rank)]
            called = [r for r in sc if r[rank] in ("0", "1")]
            corr = [r for r in called if r[rank] == "1"]
            k, nc, ns_ = len(corr), len(called), len(sc)
            s[rank] = dict(
                n_scoreable=ns_, n_called=nc, n_correct=k,
                correct=round(k / nc, 4) if nc else None,
                ci=[round(x, 4) for x in wilson_ci(k, nc)] if nc else None,
                coverage=round(nc / ns_, 4) if ns_ else None,
                correct_population=round(k / ns_, 4) if ns_ else None,
                ci_pop=[round(x, 4) for x in wilson_ci(k, ns_)] if ns_ else None)
        # score_genecap.py classifies a selection whose reps carry NO PGFams as
        # OUTCOME_PRED_EMPTY and keeps it OUT of the conditional (OK-only) macro while
        # folding it INTO the population macro. score_mag.py has no such branch, so
        # recover it here: outcome=='ok' with a blank precision IS pred_empty.
        ok = [r for r in rows if r["outcome"] == "ok" and r["precision"] != ""]
        pred_empty = [r for r in rows if r["outcome"] == "ok" and r["precision"] == ""]
        pop = [r for r in rows if r["outcome"] in ("ok", "abstained")]
        for lbl, pool in (("cond", ok), ("pop", pop)):
            for met in ("recall", "precision", "f1", "jaccard"):
                vs = []
                for r in pool:
                    if r["outcome"] == "abstained" or (r["outcome"] == "ok" and r["precision"] == ""):
                        if met != "precision":
                            vs.append(0.0)      # pred_empty recovers 0% of the metabolism
                    else:
                        vs.append(num(r[met]))
                s["genecap_%s_%s" % (lbl, met)] = macro(vs)
        s["n_ok"] = len(ok)
        s["n_pred_empty"] = len(pred_empty)
        s["n_abstained"] = sum(1 for r in rows if r["outcome"] == "abstained")
        s["n_truth_empty"] = sum(1 for r in rows if r["outcome"] == "truth_empty")
        s["coverage_answered"] = round(len(ok) / len(pop), 4) if pop else None
        s["coverage_ci"] = [round(x, 4) for x in wilson_ci(len(ok), len(pop))] if pop else None
        s["ani_classes"] = dict(collections.Counter(r["ani_class"] for r in rows))
        s["tiers"] = dict(collections.Counter(r["tier"] for r in rows))
        ns = sorted(int(r["n_selected"]) for r in rows)
        s["median_n_selected"] = ns[len(ns) // 2]
        s["mean_n_selected"] = round(sum(ns) / len(ns), 3)
        ba = sorted(x for x in (num(r["best_ani"]) for r in rows) if x is not None)
        s["median_best_ani"] = round(ba[len(ba) // 2], 2) if ba else None
        summary["c%g_x%g_%s" % (comp, contam, regime)] = s

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    out = a.out
    json.dump(dict(n_cells=n, n_batches=len(files), summary=summary), open(out, "w"), indent=1)
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()
