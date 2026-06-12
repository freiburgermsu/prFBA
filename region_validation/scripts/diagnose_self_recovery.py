#!/usr/bin/env python
"""diagnose_self_recovery.py — attribute every self-recovery FAILURE to the exact
selection stage that dropped the source organism's hit.

A *self-recovery failure* = the source organism appears in the amplicon's top-20
hits (matched by ``genome_id`` OR by taxon = ``int(genome_id.split('.')[0])``) but
is NOT among the selected representatives of the realistic (include-self) pipeline
run.  This is the positive-control round-trip: an amplicon was extracted from a
genome that is itself in the reference DB, yet the selector did not pick that
genome (or any assembly of the same organism).

For every failure we read the audit's per-hit ``decisions`` (each carries the
``stage`` + ``disposition`` the selector assigned) and find the **terminal stage**
— the stage at which the source's *furthest-surviving* self-matching hit was
eliminated — then cross-tabulate against **taxonomic salvage**: does the kept
anchor still recover the source's true lineage at species / genus / family?  A
species-preserving ``species_dedup`` drop (a conspecific sister kept instead of
the exact source assembly) is *benign* for taxonomy; a ``family_consensus`` or
``marginal_gain`` drop that loses the genus/species is *harmful*.

reads : PATHS.selection_out   (selection.json — the include-self audit w/ decisions)
        PATHS.truth_json       (benchmark_truth.json — source true lineage)
writes: DATA/self_recovery_diagnosis.json   (aggregate cause breakdown)
        DATA/self_recovery_failures.csv      (one row per failure, fully attributed)

Interpreter: /home/freiburger/Documents/py_venv/bin/python
"""
from __future__ import annotations

import csv
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import _common as C  # noqa: E402
from _common import PATHS, _norm  # noqa: E402

# Selection stages ordered early -> late.  The "terminal" cause of a self-recovery
# failure is the LATEST stage any source-matching hit reached before exclusion
# (the source's best shot at being kept).  'complete'/'not_reached' mean the
# source survived every gate but the handful was already filled by higher-ranked
# reps (a capacity/marginal-stop effect).
STAGE_ORDER = [
    "tier_floor", "reliability_gate", "identity_floor", "identity_band",
    "family_consensus", "genome_dedup", "species_dedup",
    "marginal_gain", "contamination_budget", "tier_cap", "no_gene_set",
    "complete", "not_reached", "included",
]
STAGE_IDX = {s: i for i, s in enumerate(STAGE_ORDER)}

SCORED = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


def _taxon(gid):
    try:
        return int(str(gid).split(".")[0])
    except (ValueError, AttributeError):
        return None


def _is_self(gid, src_gid, src_taxon):
    return str(gid) == str(src_gid) or (src_taxon is not None and _taxon(gid) == src_taxon)


def salvage_level(anchor, truth_lin):
    """Deepest rank at which the kept anchor's lineage still matches source truth.

    Returns one of 'species' | 'genus' | 'family' | 'lost' (anchor present but
    diverges above family) — only called when an anchor exists.
    """
    if _norm(anchor.get("species")) and _norm(anchor.get("species")) == _norm(truth_lin.get("Species")):
        return "species"
    if _norm(anchor.get("genus")) and _norm(anchor.get("genus")) == _norm(truth_lin.get("Genus")):
        return "genus"
    if _norm(anchor.get("family")) and _norm(anchor.get("family")) == _norm(truth_lin.get("Family")):
        return "family"
    return "lost"


def run():
    print("[diag] loading truth + audit (selection.json ~700MB) ...", flush=True)
    truth = json.load(open(PATHS.truth_json))
    audit = json.load(open(PATHS.selection_out))
    audit.pop("_meta", None)

    # aggregate counters
    per_region = defaultdict(lambda: Counter())          # region -> {n, top20, selected, anchor}
    fail_stage = defaultdict(Counter)                    # region -> Counter(terminal_stage)
    fail_stage_all = Counter()
    stage_by_salvage = defaultdict(Counter)              # terminal_stage -> Counter(salvage_level)
    stage_by_salvage_region = defaultdict(lambda: defaultdict(Counter))
    tie_hist = Counter()                                 # tie-cluster size among failures
    n_failures = 0
    n_amplicons = 0

    failures = []  # CSV rows

    for key, rec in audit.items():
        src_gid = key.rpartition("__")[0]
        region = key.rpartition("__")[2]
        src_taxon = _taxon(src_gid)
        tinfo = truth.get(key, {})
        truth_lin = tinfo.get("truth_lineage") or {}
        # unscoreable (all-None truth) amplicons are excluded from the benchmark too
        if not any(truth_lin.get(r) for r in SCORED):
            continue
        n_amplicons += 1

        decisions = rec.get("decisions") or []
        selected = rec.get("selected") or []

        self_dec = [d for d in decisions if _is_self(d.get("genome_id"), src_gid, src_taxon)]
        self_in_top20 = bool(self_dec)
        self_selected = any(_is_self(s.get("genome_id"), src_gid, src_taxon) for s in selected)
        self_anchor = bool(selected) and _is_self(selected[0].get("genome_id"), src_gid, src_taxon)

        reg = per_region[region]
        reg["n"] += 1
        reg["top20"] += int(self_in_top20)
        reg["selected"] += int(self_selected)
        reg["anchor"] += int(self_anchor)
        ra = per_region["ALL"]
        ra["n"] += 1
        ra["top20"] += int(self_in_top20)
        ra["selected"] += int(self_selected)
        ra["anchor"] += int(self_anchor)

        if not (self_in_top20 and not self_selected):
            continue

        # --- this is a self-recovery FAILURE: attribute it ---
        n_failures += 1
        # terminal stage = the source-matching hit that survived FURTHEST
        terminal = max(self_dec, key=lambda d: STAGE_IDX.get(d.get("stage"), -1))
        tstage = terminal.get("stage") or "unknown"
        fail_stage[region][tstage] += 1
        fail_stage_all[tstage] += 1

        # taxonomic salvage of the kept anchor (or abstain)
        if not selected:
            salv = "abstain"
        else:
            salv = salvage_level(selected[0], truth_lin)
        stage_by_salvage[tstage][salv] += 1
        stage_by_salvage_region[region][tstage][salv] += 1

        # tie context: identities present in the decision list
        idents = [d.get("identity") or 0.0 for d in decisions]
        bmax = max(idents) if idents else 0.0
        src_ident = max((d.get("identity") or 0.0) for d in self_dec)
        tie = sum(1 for x in idents if x >= bmax - 0.002)
        tie_hist[min(tie, 20)] += 1

        anchor = selected[0] if selected else {}
        failures.append(dict(
            amplicon=key, region=region, tier=rec.get("tier"),
            src_genome_id=src_gid, src_identity=round(src_ident, 4),
            best_identity=round(bmax, 4), tie_cluster=tie,
            n_self_hits=len(self_dec), terminal_stage=tstage,
            n_selected=rec.get("n_selected", len(selected)),
            anchor_genome_id=anchor.get("genome_id"),
            anchor_species=anchor.get("species"),
            anchor_family=anchor.get("family"),
            truth_species=truth_lin.get("Species"),
            truth_genus=truth_lin.get("Genus"),
            truth_family=truth_lin.get("Family"),
            salvage=salv,
        ))

    # ---------------------------------------------------------------- outputs
    def rates(c):
        n = c["n"] or 1
        return dict(
            n=c["n"], self_in_top20=c["top20"], self_selected=c["selected"],
            self_anchor=c["anchor"],
            selected_given_top20=round(c["selected"] / (c["top20"] or 1), 4),
            anchor_given_top20=round(c["anchor"] / (c["top20"] or 1), 4),
            fail_given_top20=round(1 - c["selected"] / (c["top20"] or 1), 4),
        )

    diagnosis = dict(
        n_amplicons=n_amplicons,
        n_failures=n_failures,
        overall_fail_rate_given_top20=round(
            1 - per_region["ALL"]["selected"] / (per_region["ALL"]["top20"] or 1), 4),
        per_region={r: rates(per_region[r]) for r in sorted(per_region)},
        failure_terminal_stage_overall=dict(fail_stage_all.most_common()),
        failure_terminal_stage_by_region={
            r: dict(fail_stage[r].most_common()) for r in sorted(fail_stage)},
        terminal_stage_x_salvage=dict(
            (st, dict(stage_by_salvage[st].most_common())) for st in stage_by_salvage),
        tie_cluster_hist_among_failures=dict(sorted(tie_hist.items())),
    )
    # FullLength16S detail (the positive control)
    diagnosis["FullLength16S_terminal_stage_x_salvage"] = {
        st: dict(stage_by_salvage_region["FullLength16S"][st].most_common())
        for st in stage_by_salvage_region["FullLength16S"]}

    out_json = os.path.join(C.DATA, "self_recovery_diagnosis.json")
    with open(out_json, "w") as fh:
        json.dump(diagnosis, fh, indent=1)

    out_csv = os.path.join(C.DATA, "self_recovery_failures.csv")
    if failures:
        with open(out_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(failures[0].keys()))
            w.writeheader()
            w.writerows(failures)

    # ---- console summary ----
    print(f"\n[diag] amplicons scored: {n_amplicons}; self-recovery FAILURES "
          f"(self in top20 but not selected): {n_failures} "
          f"({100*n_failures/max(n_amplicons,1):.2f}% of all amplicons)\n")
    print("[diag] per-region fail-rate GIVEN self was in top20:")
    for r in sorted(per_region):
        if r == "ALL":
            continue
        rr = rates(per_region[r])
        print(f"  {r:14s} n={rr['n']:5d} top20={rr['self_in_top20']:5d} "
              f"selected={rr['self_selected']:5d} fail/top20={rr['fail_given_top20']:.3f}")
    rr = rates(per_region["ALL"])
    print(f"  {'ALL':14s} n={rr['n']:5d} top20={rr['self_in_top20']:5d} "
          f"selected={rr['self_selected']:5d} fail/top20={rr['fail_given_top20']:.3f}\n")
    print("[diag] FAILURE terminal-stage distribution (overall):")
    for st, n in fail_stage_all.most_common():
        print(f"  {st:22s} {n:6d}  ({100*n/max(n_failures,1):.1f}%)")
    print("\n[diag] terminal_stage x taxonomic-salvage (overall):")
    print(f"  {'stage':22s} {'species':>8s} {'genus':>7s} {'family':>7s} {'lost':>6s} {'abstain':>8s}")
    for st in sorted(stage_by_salvage, key=lambda s: -sum(stage_by_salvage[s].values())):
        c = stage_by_salvage[st]
        print(f"  {st:22s} {c['species']:8d} {c['genus']:7d} {c['family']:7d} "
              f"{c['lost']:6d} {c['abstain']:8d}")
    print(f"\n[diag] wrote {out_json}\n[diag] wrote {out_csv}")


if __name__ == "__main__":
    run()
