#!/usr/bin/env python
"""Promoted design-audit script: the floor and ceiling baselines of EC-set prediction.

Promoted verbatim-in-substance from the audit reproduction scripts ``/tmp/base1.py``
(prevalence / donor floors) and ``/tmp/base2.py`` (accessory vocabulary, MCC) before
/tmp was cleared.  SPEC.md Section 13 requires ``score_folds.py`` to *import* the
baseline implementations from here rather than re-deriving them, so the numbers quoted
in SPEC.md Section 8.3 (constant-472 mean EC F1 0.637 / MCC 0.569, same-genus donor
0.803 / 0.772) and the numbers the harness reports come from one body of code.

The only change from the /tmp originals is structural: the analysis bodies moved into
``main()`` behind ``if __name__ == "__main__"`` and the scoring primitives became
module-level functions, because a script that runs a 1.7 GB load at import time cannot
be imported.  Run with no arguments to reproduce the audit printout.

Guarantees:
  * ``f1_set`` / ``mcc_set`` are the single definition of per-amplicon set agreement
    used by the audit and by score_folds.py (an empty prediction scores 0.0, never
    None; an empty truth is unscoreable and returns None).
  * ``constant_baseline`` and ``donor_baseline`` are fitted ONLY on the genomes handed
    to them, so passing a fold's surviving reference genomes keeps them train-only.
"""
import argparse
import collections
import json
import math
import os
import random
import statistics as st
import time

REPO = "/home/freiburger/Documents/prFBA"
FUNCTIONS = os.path.join(REPO, "region_validation/picrust2_comparison/data/"
                               "genome_function_sets.json")
SELECTION = os.path.join(REPO, "region_validation/data/selection_10k.json")


# --------------------------------------------------------------------------- #
# scoring primitives (imported by score_folds.py)
# --------------------------------------------------------------------------- #
def f1_set(truth, pred):
    """F1 of two sets; None when truth is empty (unscoreable), 0.0 when pred is."""
    if not truth:
        return None
    tp = len(truth & pred)
    if not pred or tp == 0:
        return 0.0
    r = tp / len(truth)
    p = tp / len(pred)
    return 2 * r * p / (r + p)


def mcc_set(truth, pred, universe):
    """Matthews correlation over a FIXED universe of possible functions.

    MCC credits the true negatives a 2,236-EC vocabulary is mostly made of, which is
    why it separates the baselines more sharply than F1 (SPEC Section 8.1).
    """
    n = len(universe)
    truth = truth & universe
    pred = pred & universe
    tp = len(truth & pred)
    fp = len(pred - truth)
    fn = len(truth - pred)
    tn = n - tp - fp - fn
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return (tp * tn - fp * fn) / den if den else None


def prevalence(train_ids, ec_of):
    """{function: fraction of training genomes carrying it} + the raw Counter."""
    cnt = collections.Counter()
    for gid in train_ids:
        cnt.update(ec_of[gid])
    n = len(train_ids)
    return ({e: c / n for e, c in cnt.items()} if n else {}), cnt


def constant_baseline(cnt, k):
    """The top-``k`` most prevalent training functions: the no-information floor.

    ``k`` is the median truth-set size over the training genomes (SPEC Section 8.3),
    so the floor is not handicapped by predicting an implausibly large or small set.
    """
    return set(e for e, _ in cnt.most_common(k))


def accessory_vocab(prev, lo=0.05, hi=0.95):
    """Functions whose training prevalence is in [lo, hi] -- the informative band."""
    return set(e for e, p in prev.items() if lo <= p <= hi)


def donor_baseline(rank_value, by_rank, ec_of, exclude, rng):
    """A random donor genome sharing ``rank_value`` -- copy-a-relative at that rank.

    Returns ``(donor_id, ec_set)`` or ``(None, set())`` when no donor survives, which
    is the normal case for the blocked rank itself (SPEC D13).
    """
    cands = [g for g in by_rank.get(rank_value, ()) if g not in exclude and ec_of.get(g)]
    if not cands:
        return None, set()
    donor = cands[rng.randrange(len(cands))]
    return donor, ec_of[donor]


def clade_max_baseline(truth, candidates, ec_of):
    """The best achievable single-genome copy among ``candidates``: the oracle ceiling."""
    best, best_f1 = None, -1.0
    for gid in candidates:
        v = f1_set(truth, ec_of.get(gid, set()))
        if v is not None and v > best_f1:
            best, best_f1 = gid, v
    return best, (best_f1 if best is not None else None)


# --------------------------------------------------------------------------- #
# audit reproduction (the /tmp/base1.py + /tmp/base2.py bodies)
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--functions", default=FUNCTIONS)
    ap.add_argument("--selection", default=SELECTION)
    ap.add_argument("--n-test", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1729)
    args = ap.parse_args()

    import orjson
    t = time.time()
    d = orjson.loads(open(args.functions, "rb").read())
    print("loaded", len(d), round(time.time() - t, 1), "s", flush=True)
    sel = json.load(open(args.selection))["genomes"]
    panel = [g for g in sel if g["genome_id"] in d]
    print("panel in fs", len(panel))
    ecs = {g["genome_id"]: set(d[g["genome_id"]]["ec"]) for g in panel}
    meta = {g["genome_id"]: g for g in panel}
    sizes = [len(v) for v in ecs.values()]
    print("EC set size: median", st.median(sizes), "mean", round(st.mean(sizes), 1),
          "empty", sum(1 for s in sizes if s == 0))

    rng = random.Random(args.seed)
    ids = [g["genome_id"] for g in panel if ecs[g["genome_id"]]]
    rng.shuffle(ids)
    test, train = ids[:args.n_test], ids[args.n_test:]
    prev, cnt = prevalence(train, ecs)

    median_truth = int(st.median([len(ecs[i]) for i in train]))
    for K in sorted({200, 300, 400, median_truth, 600, 800}):
        pred = constant_baseline(cnt, K)
        vals = [v for v in (f1_set(ecs[i], pred) for i in test) if v is not None]
        print("top-%d prevalence baseline: mean F1 %.3f%s"
              % (K, st.mean(vals), "  <- median truth size" if K == median_truth else ""))
    for thr in (0.5, 0.7, 0.9):
        pred = set(e for e, p in prev.items() if p >= thr)
        vals = [v for v in (f1_set(ecs[i], pred) for i in test) if v is not None]
        print("prev>=%.1f (|pred|=%d): mean F1 %.3f" % (thr, len(pred), st.mean(vals)))

    by_genus, by_family = collections.defaultdict(list), collections.defaultdict(list)
    for g in panel:
        by_genus[g["genus"]].append(g["genome_id"])
        by_family[g["family"]].append(g["genome_id"])
    rows = collections.defaultdict(list)
    for i in test:
        donor = train[rng.randrange(len(train))]
        for name, pred in (("random panel genome donor", ecs[donor]),
                           ("random SAME-GENUS donor",
                            donor_baseline(meta[i]["genus"], by_genus, ecs, {i}, rng)[1]),
                           ("random SAME-FAMILY donor",
                            donor_baseline(meta[i]["family"], by_family, ecs, {i}, rng)[1])):
            if not pred:
                continue
            v = f1_set(ecs[i], pred)
            if v is not None:
                rows[name].append(v)
    for name in sorted(rows):
        print("%s: mean F1 %.3f (n=%d)" % (name, st.mean(rows[name]), len(rows[name])))

    # ---- /tmp/base2.py: accessory vocabulary and MCC ---------------------- #
    acc = accessory_vocab(prev)
    core = set(e for e, p in prev.items() if p > 0.95)
    print("vocab sizes: total %d, core(>0.95) %d, accessory(0.05-0.95) %d, rare(<0.05) %d"
          % (len(prev), len(core), len(acc),
             sum(1 for p in prev.values() if p < 0.05)))
    V_all = set(prev)
    pred_k = constant_baseline(cnt, median_truth)
    rows = collections.defaultdict(list)
    for i in test:
        donor, dec = donor_baseline(meta[i]["genus"], by_genus, ecs, {i}, rng)
        arms = [("constant-%d" % median_truth, pred_k)]
        if donor:
            arms.append(("genus-donor", dec))
        for name, p in arms:
            for vn, V in (("full", V_all), ("accessory", acc)):
                v = f1_set(ecs[i] & V, p & V)
                if v is not None:
                    rows[(name, vn, "F1")].append(v)
                m = mcc_set(ecs[i], p, V)
                if m is not None:
                    rows[(name, vn, "MCC")].append(m)
    for k in sorted(rows):
        print(k, "mean %.3f" % st.mean(rows[k]), "n=%d" % len(rows[k]))
    have = sum(1 for i in test
               if any(x != i for x in by_genus[meta[i]["genus"]]))
    havef = sum(1 for i in test
                if any(x != i for x in by_family[meta[i]["family"]]))
    print("test genomes with >=1 other panel genome in same genus: %d/%d" % (have, len(test)))
    print("test genomes with >=1 other panel genome in same family: %d/%d" % (havef, len(test)))


if __name__ == "__main__":
    main()
