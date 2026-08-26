#!/usr/bin/env python
"""mcnemar_regions.py -- paired McNemar tests across the 16S region panel (DESIGN.md Section 6).

Cross-region accuracy differences in this benchmark are NOT independent: each of
the ~10,000 source organisms contributes up to 8 region-amplicons, so region A
and region B are scored on overlapping organisms.  Comparing their Wilson
intervals treats correlated data as independent.  The correct test is a PAIRED
one on the COMMON subset of organisms amplified in *both* regions (DESIGN Section
6 prescribes McNemar on the common-amplified subset).

For a chosen predictor (anchor|consensus), rank (default Genus), and self_mode
(include | exclude_0.987), this:

  1. reads per-record correctness from ``taxacc_per_record.jsonl`` and reduces it
     to one binary per (source organism, region): correct == (cells[rank] is True).
     A no-call / abstention (cell None) or a misassignment (False) is NOT correct
     -- the whole-population definition (matches ``correct_population``).
  2. for every region PAIR, builds the discordant counts on the pair's common
     organisms: b = A-correct & B-wrong, c = A-wrong & B-correct, and runs
     McNemar's test -- continuity-corrected chi-square when b+c is large, the
     exact binomial when b+c is small (< 25).
  3. Holm-adjusts the p-values across the 28 region pairs (per predictor x rank x
     self_mode) to control the family-wise error rate.

Outputs (per predictor x rank x self_mode):
  data/mcnemar_region_pairs.csv   region_a, region_b, n_common, acc_a, acc_b,
                                  b, c, statistic, test, p_value, p_holm, significant
  data/mcnemar_regions.json       + the common-subset marginal accuracy per region

The common-subset marginal accuracy is reported alongside so the pairwise deltas
are read on one apples-to-apples denominator (the strict all-8-regions common
subset is also reported as `n_common_all8`).

Interpreter: /home/freiburger/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import sys

from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common as C  # noqa: E402
from _common import PATHS  # noqa: E402

RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
PREDICTORS = ["anchor", "consensus"]
SELF_MODES = ["include", "exclude_0.987"]


def load_correct(jsonl_path, predictor, rank, self_mode):
    """{src_genome_id: {region: correct_bool}} for one (predictor, rank, self_mode).

    correct == cells[rank] is True; None (no-call/abstain) and False both -> not
    correct (the whole-population accuracy definition).
    """
    per_org = {}
    with open(jsonl_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            r = json.loads(line)
            org = str(r.get("src_genome_id"))
            region = r.get("region")
            sm = (r.get("self_modes") or {}).get(self_mode)
            if not region or sm is None:
                continue
            cell = ((sm.get("scores") or {}).get(predictor) or {}).get("cells", {}).get(rank)
            per_org.setdefault(org, {})[region] = (cell is True)
    return per_org


def mcnemar(b, c):
    """McNemar's paired test on discordant counts (b, c).

    Continuity-corrected chi-square (1 df) when b+c is reasonably large; the exact
    two-sided binomial(min(b,c); b+c, 0.5) when b+c < 25.  Returns (stat, p, test).
    """
    n = b + c
    if n == 0:
        return 0.0, 1.0, "none"
    if n < 25:
        p = float(stats.binomtest(min(b, c), n, 0.5, alternative="two-sided").pvalue)
        return float(min(b, c)), p, "exact_binomial"
    stat = (abs(b - c) - 1) ** 2 / n
    p = float(stats.chi2.sf(stat, 1))
    return float(stat), p, "chi2_cc"


def holm(pairs):
    """Holm-Bonferroni adjust a list of dicts carrying 'p_value'; sets 'p_holm'."""
    order = sorted(range(len(pairs)), key=lambda i: pairs[i]["p_value"])
    m = len(pairs)
    running = 0.0
    for rank_i, idx in enumerate(order):
        adj = min(1.0, (m - rank_i) * pairs[idx]["p_value"])
        running = max(running, adj)  # enforce monotonic non-decreasing
        pairs[idx]["p_holm"] = round(running, 6)
    return pairs


def analyze(per_org, alpha=0.05):
    """Pairwise McNemar over all region pairs + common-subset marginal accuracy."""
    regions = sorted({rg for d in per_org.values() for rg in d})
    # strict all-regions common subset size (reported as context)
    n_all = sum(1 for d in per_org.values() if all(rg in d for rg in regions))
    rows = []
    for a, bcol in itertools.combinations(regions, 2):
        common = [d for d in per_org.values() if a in d and bcol in d]
        n = len(common)
        n_a = sum(1 for d in common if d[a])
        n_b = sum(1 for d in common if d[bcol])
        b = sum(1 for d in common if d[a] and not d[bcol])
        c = sum(1 for d in common if not d[a] and d[bcol])
        stat, p, test = mcnemar(b, c)
        rows.append(dict(
            region_a=a, region_b=bcol, n_common=n,
            acc_a=round(n_a / n, 4) if n else None,
            acc_b=round(n_b / n, 4) if n else None,
            b=b, c=c, statistic=round(stat, 4), test=test, p_value=round(p, 6)))
    holm(rows)
    for r in rows:
        r["significant"] = bool(r["p_holm"] < alpha)
    # marginal accuracy per region on its own full set (context)
    marg = {}
    for rg in regions:
        vals = [d[rg] for d in per_org.values() if rg in d]
        marg[rg] = dict(n=len(vals), accuracy=round(sum(vals) / len(vals), 4) if vals else None)
    return dict(regions=regions, n_common_all8=n_all, marginal=marg, pairs=rows)


def run(force=False, jsonl=None):
    jsonl = jsonl or PATHS.taxacc_jsonl
    out_csv = os.path.join(os.path.dirname(PATHS.taxacc_summary_json), "mcnemar_region_pairs.csv")
    out_json = os.path.join(os.path.dirname(PATHS.taxacc_summary_json), "mcnemar_regions.json")
    if not force and C.done(out_csv) and C.done(out_json):
        return f"[mcnemar] outputs present; skip (use --force). -> {out_csv}, {out_json}"
    if not C.done(jsonl):
        raise SystemExit(f"[mcnemar] missing per-record scores: {jsonl} (run score_taxacc first)")

    result = {"_meta": dict(
        test="paired McNemar on common-amplified organisms; Holm-adjusted across 28 region pairs",
        source=jsonl, predictors=PREDICTORS, self_modes=SELF_MODES,
        correct_definition="cells[rank] is True (no-call/abstain and misassign both count as not-correct)",
    ), "by": {}}

    csv_rows = []
    for predictor in PREDICTORS:
        for rank in ("Genus", "Family"):
            for sm in SELF_MODES:
                per_org = load_correct(jsonl, predictor, rank, sm)
                res = analyze(per_org)
                key = f"{predictor}|{rank}|{sm}"
                result["by"][key] = res
                for r in res["pairs"]:
                    csv_rows.append(dict(predictor=predictor, rank=rank, self_mode=sm, **r))

    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(csv_rows[0].keys()))
        w.writeheader()
        w.writerows(csv_rows)
    with open(out_json, "w") as fh:
        json.dump(result, fh, indent=1)

    # headline print: anchor / Genus / exclude_0.987
    hk = "anchor|Genus|exclude_0.987"
    hp = result["by"][hk]
    nsig = sum(1 for r in hp["pairs"] if r["significant"])
    return (f"[mcnemar] wrote {out_csv} ({len(csv_rows)} pair-tests across "
            f"{len(result['by'])} predictor x rank x self_mode cells).\n"
            f"[mcnemar] {hk}: {nsig}/{len(hp['pairs'])} region pairs differ significantly "
            f"(Holm p<0.05); all-8 common subset n={hp['n_common_all8']}.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--jsonl", default=None, help="override per-record path")
    args = ap.parse_args()
    print(run(force=args.force, jsonl=args.jsonl))


if __name__ == "__main__":
    main()
