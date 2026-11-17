#!/usr/bin/env python
"""Path A reach of the 17 eligible genomes whose EC sets the 2026-11-13 repair changed
(prFBA_paper sec:sharedref: "... and 130 of the 63,871 test cells").

Counts the genus-rank (leave-genus-out) test cells whose source genome is one of the 17
(`changed_eligible` in rescore_sources/changed_sets.json), at density 1.0 in EC space for one arm
(picrust2; every arm is scored on the same cells), against the total of scored cells, by fold.

Inputs: region_validation/shared_reference/data/cells.csv.gz (git-ignored, 56 MB) and
        rescore_sources/changed_sets.json.
Writes manuscript_numbers/rescore_sources/changed_cells.json.
"""
import collections, csv, gzip, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.dirname(os.path.dirname(HERE))
SR = os.path.join(R, "region_validation", "shared_reference", "data")


def main():
    changed = set(json.load(open(os.path.join(HERE, "changed_sets.json")))["changed_eligible"])
    by_fold, genomes, n_scored, n_rows, outcomes = collections.Counter(), set(), 0, 0, collections.Counter()
    with gzip.open(os.path.join(SR, "cells.csv.gz"), "rt") as fh:
        for row in csv.DictReader(fh):
            if row["density"] != "1.0" or row["space"] != "ec" or row["arm"] != "picrust2":
                continue
            n_rows += 1
            outcomes[row["outcome"]] += 1
            scoreable = row["outcome"] in ("scored", "abstained", "pred_empty")
            n_scored += scoreable
            if row["src_genome_id"] in changed:
                genomes.add(row["src_genome_id"])
                if scoreable:
                    by_fold[row["fold"]] += 1
    out = {"what": "genus-rank test cells (density 1.0, EC, one arm) whose source genome's EC set the repair changed",
           "generator": "manuscript_numbers/rescore_sources/count_changed_cells.py",
           "n_changed_eligible": len(changed), "n_changed_genomes_with_cells": len(genomes),
           "n_cells_changed_truth": sum(by_fold.values()),
           "by_fold": dict(sorted(by_fold.items())),
           "n_cells_scored_total": n_scored, "n_rows_total": n_rows, "outcomes": dict(outcomes)}
    with open(os.path.join(HERE, "changed_cells.json"), "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
