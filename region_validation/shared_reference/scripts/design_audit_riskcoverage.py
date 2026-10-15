#!/usr/bin/env python
"""Promoted design-audit script: risk-coverage curves, AURC, and the abstention control.

Promoted from ``/tmp/rc.py`` before /tmp was cleared.  The audit's finding is the reason
SPEC.md Section D7 makes risk-coverage co-primary: on the shipped data, letting PICRUSt2
abstain on its own worst 24.9% of queries by NSTI shrinks prFBA's advantage from +0.070
to +0.017, so a bare conditional F1 compares two differently-selected populations.  It
also measured that prFBA's abstention is informative (0.7985 own-choice vs 0.6006 at the
same coverage chosen at random), which is the ``base_random_abstain`` arm.

The only change from the /tmp original is structural: the sweep became a pure function
and the script body moved into ``main()``.  SPEC.md Section 13 requires score_folds.py
to import the implementation from here.

Guarantees:
  * each method sweeps its OWN confidence knob (prFBA: rank-1 identity, descending =
    more confident; PICRUSt2: NSTI, ascending = more confident), never the score being
    measured, so no curve is an oracle;
  * coverage is always a fraction of the full scoreable denominator, so two arms'
    curves are comparable at equal coverage;
  * a missing knob sorts last (least confident), never first.
"""
import argparse
import csv
import gzip
import os
import random
import statistics as st

REPO = "/home/freiburger/Documents/prFBA"
CELLS = os.path.join(REPO, "region_validation/picrust2_comparison/data/"
                           "comparison_cells.csv.gz")
COVERAGE_GRID = (0.5, 0.75, 0.9, 1.0)


def order_by_confidence(items, higher_is_better=True, seed=None):
    """``items`` = [(knob_or_None, payload)] -> payloads, most confident first.

    ``None`` knobs always sort last.  With ``seed`` given and every knob None the order
    is a seeded shuffle (the uninformative-confidence null), not input order.
    """
    missing = [p for k, p in items if k is None]
    known = [(k, p) for k, p in items if k is not None]
    known.sort(key=lambda kp: -kp[0] if higher_is_better else kp[0])
    if not known and seed is not None:
        rng = random.Random(seed)
        rng.shuffle(missing)
    return [p for _, p in known] + missing


def risk_coverage(scores, knobs, higher_is_better=True, grid=COVERAGE_GRID, seed=None):
    """Mean score over the most-confident ``c`` fraction, for each ``c`` in ``grid``.

    ``scores`` and ``knobs`` are parallel sequences over the FULL scoreable denominator
    (an abstention contributes its population score, normally 0.0).  Returns
    ``{"curve": [(coverage, mean_score)], "aurc": float, "at": {coverage: mean}}``.
    """
    ordered = order_by_confidence(list(zip(knobs, scores)), higher_is_better, seed)
    n = len(ordered)
    if n == 0:
        return {"curve": [], "aurc": None, "at": {c: None for c in grid}}
    cum, running = [], 0.0
    for i, v in enumerate(ordered, start=1):
        running += v
        cum.append(running / i)
    at = {}
    for c in grid:
        k = max(1, min(n, int(round(c * n))))
        at[c] = cum[k - 1]
    # AURC: trapezoid of mean-score vs coverage over the full sweep (coverage i/n).
    xs = [i / n for i in range(1, n + 1)]
    area = 0.0
    for i in range(1, n):
        area += 0.5 * (cum[i] + cum[i - 1]) * (xs[i] - xs[i - 1])
    area += 0.5 * (cum[0] + cum[0]) * xs[0]  # [0, 1/n] rectangle at the first point
    return {"curve": list(zip(xs, cum)), "aurc": area, "at": at}


def random_abstain_mask(keys, coverage, seed):
    """Deterministic set of keys kept when abstaining at random to ``coverage``."""
    rng = random.Random(seed)
    keys = sorted(keys)
    k = int(round(coverage * len(keys)))
    return set(rng.sample(keys, max(0, min(len(keys), k))))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", default=CELLS)
    ap.add_argument("--arm-a", default="prfba_exclude")
    ap.add_argument("--arm-b", default="picrust2")
    ap.add_argument("--method-field", default="method")
    ap.add_argument("--confidence-field", default="best_nonself_identity")
    args = ap.parse_args()

    pr, pc = {}, {}
    with gzip.open(args.cells, "rt") as fh:
        for row in csv.DictReader(fh):
            if row["space"] != "ec":
                continue
            if row[args.method_field] == args.arm_a:
                pr[row["amplicon"]] = row
            elif row[args.method_field] == args.arm_b:
                pc[row["amplicon"]] = row
    amps = sorted(a for a in pr if a in pc and int(pr[a]["n_truth"] or 0) > 0)
    print("scoreable", len(amps))

    def f1(row):
        return 0.0 if row["f1"] == "" else float(row["f1"])

    answered = [a for a in amps if pr[a]["outcome"] in ("ok", "scored")]
    cov = len(answered) / len(amps)
    print("%s coverage %.3f ; conditional mean F1 %.4f"
          % (args.arm_a, cov, st.mean([f1(pr[a]) for a in answered])))
    print("%s coverage %.3f ; mean F1 all %.4f"
          % (args.arm_b, sum(1 for a in amps if pc[a]["outcome"] in ("ok", "scored"))
             / len(amps), st.mean([f1(pc[a]) for a in amps])))

    b = risk_coverage([f1(pc[a]) for a in amps],
                      [float(pc[a]["nsti"]) if pc[a]["nsti"] else None for a in amps],
                      higher_is_better=False, grid=(1.0, 0.9, 0.8, cov, 0.6, 0.5))
    for c in sorted(b["at"], reverse=True):
        print("%s @coverage %.3f : mean F1 %.4f" % (args.arm_b, c, b["at"][c]))
    knob = args.confidence_field
    a_rc = risk_coverage([f1(pr[a]) for a in amps],
                         [float(pr[a][knob]) if pr[a].get(knob) else None for a in amps],
                         higher_is_better=True, grid=(cov, 0.6, 0.5, 0.4))
    for c in sorted(a_rc["at"], reverse=True):
        print("%s @coverage %.3f : mean F1 %.4f" % (args.arm_a, c, a_rc["at"][c]))
    keep = random_abstain_mask(amps, cov, 3)
    print("%s @random %.3f coverage: mean F1 %.4f (vs its own-choice %.4f)"
          % (args.arm_a, cov, st.mean([f1(pr[a]) for a in sorted(keep)]),
             st.mean([f1(pr[a]) for a in answered])))
    print("AURC: %s %.4f   %s %.4f"
          % (args.arm_a, a_rc["aurc"] or 0.0, args.arm_b, b["aurc"] or 0.0))


if __name__ == "__main__":
    main()
