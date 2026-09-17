#!/usr/bin/env python
"""Compare prFBA and PICRUSt2 functional predictions on the 10k region-validation set.

Both methods are scored against the SAME truth (the source genome's BV-BRC PGFam set,
projected to EC numbers) using the SAME projection for every set involved, so no method
is credited or penalised for the projection:

  truth      source genome PGFams          -> EC -> ModelSEED reactions
  prFBA      union of selected references'  -> EC -> ModelSEED reactions
             PGFams (the synthetic genome)
  PICRUSt2   predicted EC copy numbers > 0        -> ModelSEED reactions

EC space is primary and restricted to the vocabulary both systems can express
(``shared`` universe in annotation_ceiling.json); reaction space is the secondary,
modelling-relevant view. prFBA's extra PGFam-product -> ModelSEED-role route is
reported separately (``rxn_role``), never mixed into the head-to-head.

Regimes / strata, since PICRUSt2 cannot exclude the source organism from its own
reference:
  prFBA include / exclude_0.987 (re-selected after dropping self + >=tau neighbours)
  PICRUSt2 stratified by whether the source assembly is IN its reference, by NSTI,
  and by prFBA's best-non-self identity -- an imperfect match to prFBA's regimes
  (GTDB is species-dereplicated, so "not in reference" often still has a conspecific).

Abstentions (prFBA below the family floor; PICRUSt2 sequences dropped by --min_align)
are scored both ways: conditional (cells with a prediction) and population
(abstention = recall 0 / F1 0), the denominators DESIGN.md Section 7 specifies.
"""
import argparse
import csv
import gzip
import json
import math
import os
import statistics as st
import sys
from collections import defaultdict

import orjson

HERE = os.path.dirname(os.path.abspath(__file__))
RV = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(RV, "scripts"))
sys.path.insert(0, HERE)
import _common as C  # noqa: E402
from _common import wilson_ci  # noqa: E402
from annotation_ceiling import acc_core, load_reference_ec, load_reference_meta  # noqa: E402

MODE_INCLUDE, MODE_EXCLUDE = "include", "exclude_0.987"
NSTI_BIN_EDGES = [0.0, 0.01, 0.05, 0.15, 0.5, float("inf")]
NSTI_BIN_LABELS = ["<0.01", "0.01-0.05", "0.05-0.15", "0.15-0.5", ">=0.5"]
IDENT_BIN_EDGES = [0.0, 0.865, 0.945, 0.987, 0.995, 1.0001]
IDENT_BIN_LABELS = ["<0.865", "0.865-0.945", "0.945-0.987", "0.987-0.995", ">=0.995"]
CELL_FIELDS = ["amplicon", "region", "domain", "method", "space", "outcome",
               "n_truth", "n_pred", "tp", "recall", "precision", "f1", "jaccard",
               "n_reps", "tier", "nsti", "nsti_bin", "picrust2_domain",
               "src_in_picrust2_ref", "best_nonself_identity", "ident_bin"]


def bin_of(value, edges, labels, unknown="unknown"):
    if value is None:
        return unknown
    for i in range(len(labels)):
        if edges[i] <= value < edges[i + 1]:
            return labels[i]
    return labels[-1]


def metrics(truth, pred):
    """(tp, recall, precision, f1, jaccard); None where undefined."""
    tp = len(truth & pred)
    recall = tp / len(truth) if truth else None
    precision = tp / len(pred) if pred else None
    f1 = (2 * recall * precision / (recall + precision)
          if recall and precision else (0.0 if truth and pred else None))
    jaccard = tp / len(truth | pred) if (truth or pred) else None
    return tp, recall, precision, f1, jaccard


class Agg:
    """Macro (per-amplicon mean) + micro (pooled TP/FP/FN) accumulators.

    ``n`` counts every scoreable amplicon (truth non-empty). ``called`` counts those
    with a non-empty prediction. Macro means are reported over both denominators:
    conditional (called only) and population (abstention = 0).
    """

    __slots__ = ("n", "called", "vals", "tp", "fp", "fn")

    def __init__(self):
        self.n = self.called = self.tp = self.fp = self.fn = 0
        self.vals = {"recall": [], "precision": [], "f1": [], "jaccard": []}

    def add(self, truth_n, pred_n, tp, recall, precision, f1, jaccard):
        self.n += 1
        if pred_n:
            self.called += 1
            for k, v in (("recall", recall), ("precision", precision),
                         ("f1", f1), ("jaccard", jaccard)):
                if v is not None:
                    self.vals[k].append(v)
            self.tp += tp
            self.fp += pred_n - tp
            self.fn += truth_n - tp
        else:
            self.fn += truth_n

    def block(self):
        def macro(key, population):
            vals = self.vals[key]
            if population:
                vals = vals + [0.0] * (self.n - self.called)
            if not vals:
                return None
            return {"value": round(st.mean(vals), 4), "n": len(vals),
                    "sd": round(st.pstdev(vals), 4) if len(vals) > 1 else 0.0}
        micro_r = self.tp / (self.tp + self.fn) if (self.tp + self.fn) else None
        micro_p = self.tp / (self.tp + self.fp) if (self.tp + self.fp) else None
        out = {
            "n": self.n, "n_called": self.called,
            "coverage": round(self.called / self.n, 4) if self.n else None,
            "conditional": {k: macro(k, False) for k in self.vals},
            "population": {k: macro(k, True) for k in ("recall", "f1", "jaccard")},
            "micro": {"recall": round(micro_r, 4) if micro_r is not None else None,
                      "precision": round(micro_p, 4) if micro_p is not None else None,
                      "tp": self.tp, "fp": self.fp, "fn": self.fn},
        }
        if micro_r is not None:
            lo, hi = wilson_ci(self.tp, self.tp + self.fn)
            out["micro"]["recall_ci"] = [round(lo, 4), round(hi, 4)]
        return out


def load_picrust2(outdir, ec_index):
    """{md5: (frozenset(ec ids), nsti, best_domain, closest_ref)} over every batch."""
    import numpy as np
    import pandas as pd
    preds, meta = {}, {}
    batches = sorted(d for d in os.listdir(outdir) if d.startswith("batch_"))
    if not batches:
        batches = [""]
    for b in batches:
        d = os.path.join(outdir, b)
        if not os.path.exists(os.path.join(d, "combined_EC_predicted.tsv.gz")):
            continue
        with gzip.open(os.path.join(d, "combined_marker_predicted_and_nsti.tsv.gz"),
                       "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                meta[r["sequence"]] = (float(r["metadata_NSTI"]), r["best_domain"],
                                       r.get("closest_reference_genome"))
        df = pd.read_csv(os.path.join(d, "combined_EC_predicted.tsv.gz"), sep="\t",
                         index_col=0)
        cols = np.array([ec_index.setdefault(c.replace("EC:", ""), len(ec_index))
                         for c in df.columns])
        arr = df.to_numpy() > 0
        for i, md5 in enumerate(df.index):
            preds[str(md5)] = frozenset(cols[arr[i]].tolist())
    return preds, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    D = os.path.join(HERE, "..", "data")
    ap.add_argument("--picrust2", default=os.path.join(D, "full_picrust2"))
    ap.add_argument("--selection", default=os.path.join(D, "prfba_selection.json"))
    ap.add_argument("--legacy-selection", default=os.path.join(D,
                    "prfba_selection_legacy.json"))
    ap.add_argument("--functions", default=os.path.join(D, "genome_function_sets.json"))
    ap.add_argument("--maps", default=os.path.join(D, "function_maps.json"))
    ap.add_argument("--genomes", default=os.path.join(D, "bvbrc_source_genomes.json"))
    ap.add_argument("--cells", default=os.path.join(D, "comparison_cells.csv.gz"))
    ap.add_argument("--summary", default=os.path.join(D, "comparison_summary.json"))
    ap.add_argument("--routing", default=os.path.join(D, "domain_routing.csv"))
    args = ap.parse_args()

    truth_meta = orjson.loads(open(C.PATHS.truth_json, "rb").read())
    fanout = orjson.loads(open(C.PATHS.amplicon_fanout, "rb").read())
    funcs = orjson.loads(open(args.functions, "rb").read())
    maps = orjson.loads(open(args.maps, "rb").read())
    sel = orjson.loads(open(args.selection, "rb").read())
    sel_legacy = (orjson.loads(open(args.legacy_selection, "rb").read())
                  if os.path.exists(args.legacy_selection) else {})
    bvbrc = orjson.loads(open(args.genomes, "rb").read())

    # ---- integer id spaces (EC strings and reaction ids are shared by every set) ----
    ec_index, rxn_index = {}, {}
    ec_id = lambda e: ec_index.setdefault(e, len(ec_index))  # noqa: E731
    rxn_id = lambda r: rxn_index.setdefault(r, len(rxn_index))  # noqa: E731
    ec_rxn = {ec_id(e): frozenset(rxn_id(r) for r in v) for e, v in maps["ec_rxn"].items()}
    g_ec, g_rxn_ec, g_rxn_role = {}, {}, {}
    for gid, rec in funcs.items():
        if gid == "_meta":
            continue
        g_ec[gid] = frozenset(ec_id(e) for e in rec["ec"])
        g_rxn_ec[gid] = frozenset(rxn_id(r) for r in rec["rxn_ec"])
        g_rxn_role[gid] = frozenset(rxn_id(r) for r in rec["rxn_role"])

    p2_pred, p2_meta = load_picrust2(args.picrust2, ec_index)
    print(f"[compare] PICRUSt2 predictions for {len(p2_pred)} sequences; "
          f"{len(ec_index)} ECs, {len(rxn_index)} reactions in play")

    # ---- shared EC vocabulary: expressible by BOTH annotation systems ----
    ref_ec, ref_vocab = load_reference_ec()
    ref_taxid = load_reference_meta()
    bvbrc_vocab = {e for v in maps["pgfam_ec"].values() for e in v}
    shared_ec = frozenset(ec_id(e) for e in (bvbrc_vocab & ref_vocab))
    shared_rxn = frozenset(r for e in shared_ec for r in ec_rxn.get(e, ()))
    src_in_ref = {gid: (acc_core((rec or {}).get("assembly_accession")) in ref_ec)
                  for gid, rec in bvbrc.items()}

    # ---- per-amplicon scoring -------------------------------------------------- #
    aggs = defaultdict(Agg)          # (method, space, region, domain) -> Agg
    strat = defaultdict(Agg)         # (method, space, stratum_kind, stratum) -> Agg
    paired = defaultdict(lambda: defaultdict(dict))   # region -> amplicon -> method -> f1
    routing = defaultdict(int)
    agree = defaultdict(list)        # (space, region) -> prFBA-vs-PICRUSt2 jaccard
    md5_of = {}
    for md5, keys in fanout.items():
        for k in keys:
            md5_of[k] = md5

    cells_fh = gzip.open(args.cells, "wt", newline="")
    writer = csv.DictWriter(cells_fh, fieldnames=CELL_FIELDS)
    writer.writeheader()

    for key in sorted(truth_meta):
        t = truth_meta[key]
        src, region, domain = t["src_genome_id"], t["region"], t["domain"]
        rec = sel.get(key, {})
        md5 = md5_of.get(key)
        nsti, p2_domain, _closest = p2_meta.get(md5, (None, None, None))
        in_ref = src_in_ref.get(src, False)
        bnsi = rec.get("best_nonself_identity")
        ident_bin = bin_of(bnsi, IDENT_BIN_EDGES, IDENT_BIN_LABELS)
        nsti_bin = bin_of(nsti, NSTI_BIN_EDGES, NSTI_BIN_LABELS)
        routing[(region, domain, p2_domain or "not_placed")] += 1

        preds = {}
        for regime, label in ((MODE_INCLUDE, "prfba_include"),
                              (MODE_EXCLUDE, "prfba_exclude")):
            reps = rec.get(regime, {}).get("reps", [])
            preds[label] = (reps, rec.get(regime, {}).get("tier"),
                            frozenset().union(*[g_ec[g] for g in reps]) if reps else frozenset(),
                            frozenset().union(*[g_rxn_ec[g] for g in reps]) if reps else frozenset(),
                            frozenset().union(*[g_rxn_role[g] for g in reps]) if reps else frozenset())
        for regime, label in ((MODE_INCLUDE, "prfba_include_legacy"),
                              (MODE_EXCLUDE, "prfba_exclude_legacy")):
            lrec = sel_legacy.get(key, {}).get(regime, {})
            reps = lrec.get("reps", [])
            preds[label] = (reps, lrec.get("tier"),
                            frozenset().union(*[g_ec[g] for g in reps]) if reps else frozenset(),
                            frozenset().union(*[g_rxn_ec[g] for g in reps]) if reps else frozenset(),
                            frozenset().union(*[g_rxn_role[g] for g in reps]) if reps else frozenset())
        p2_ec_pred = p2_pred.get(md5, frozenset())
        preds["picrust2"] = ([], None, p2_ec_pred,
                             frozenset(r for e in p2_ec_pred for r in ec_rxn.get(e, ())),
                             frozenset())

        truth_sets = {"ec": g_ec.get(src, frozenset()) & shared_ec,
                      "rxn": g_rxn_ec.get(src, frozenset()) & shared_rxn,
                      "rxn_role": g_rxn_role.get(src, frozenset())}
        for space in ("ec", "rxn", "rxn_role"):
            truth = truth_sets[space]
            for method, (reps, tier, pec, prxn, prole) in preds.items():
                if space == "rxn_role" and method == "picrust2":
                    continue  # prFBA-only route: no role strings on the PICRUSt2 side
                pred = {"ec": pec & shared_ec, "rxn": prxn & shared_rxn,
                        "rxn_role": prole}[space]
                if not truth:
                    outcome = "truth_empty"
                    writer.writerow({"amplicon": key, "region": region, "domain": domain,
                                     "method": method, "space": space, "outcome": outcome,
                                     "n_truth": 0, "n_pred": len(pred), "tp": "",
                                     "recall": "", "precision": "", "f1": "", "jaccard": "",
                                     "n_reps": len(reps), "tier": tier or "",
                                     "nsti": nsti if nsti is not None else "",
                                     "nsti_bin": nsti_bin, "picrust2_domain": p2_domain or "",
                                     "src_in_picrust2_ref": int(in_ref),
                                     "best_nonself_identity": bnsi if bnsi is not None else "",
                                     "ident_bin": ident_bin})
                    continue
                tp, r, p, f1, j = metrics(truth, pred)
                outcome = "abstained" if not pred else "ok"
                aggs[(method, space, region, domain)].add(len(truth), len(pred), tp, r, p, f1, j)
                aggs[(method, space, region, "all")].add(len(truth), len(pred), tp, r, p, f1, j)
                aggs[(method, space, "all", "all")].add(len(truth), len(pred), tp, r, p, f1, j)
                for kind, value in (("src_in_picrust2_ref", str(in_ref)),
                                    ("nsti_bin", nsti_bin), ("ident_bin", ident_bin)):
                    strat[(method, space, kind, value)].add(len(truth), len(pred), tp,
                                                            r, p, f1, j)
                if space == "ec" and method in ("prfba_include", "prfba_exclude", "picrust2"):
                    # (F1, made-a-call): an abstention scores 0 in the population view and
                    # is excluded from the conditional one -- which is NOT the same as
                    # "F1 > 0" (a prediction can overlap truth not at all).
                    paired[region][key][method] = (f1 if f1 is not None else 0.0,
                                                   bool(pred))
                writer.writerow({"amplicon": key, "region": region, "domain": domain,
                                 "method": method, "space": space, "outcome": outcome,
                                 "n_truth": len(truth), "n_pred": len(pred), "tp": tp,
                                 "recall": round(r, 4) if r is not None else "",
                                 "precision": round(p, 4) if p is not None else "",
                                 "f1": round(f1, 4) if f1 is not None else "",
                                 "jaccard": round(j, 4) if j is not None else "",
                                 "n_reps": len(reps), "tier": tier or "",
                                 "nsti": nsti if nsti is not None else "",
                                 "nsti_bin": nsti_bin, "picrust2_domain": p2_domain or "",
                                 "src_in_picrust2_ref": int(in_ref),
                                 "best_nonself_identity": bnsi if bnsi is not None else "",
                                 "ident_bin": ident_bin})
            if space in ("ec", "rxn") and truth:
                a = {"ec": preds["prfba_include"][2] & shared_ec,
                     "rxn": preds["prfba_include"][3] & shared_rxn}[space]
                b = {"ec": p2_ec_pred & shared_ec,
                     "rxn": preds["picrust2"][3] & shared_rxn}[space]
                if a or b:
                    agree[(space, region)].append(len(a & b) / len(a | b))
    cells_fh.close()

    # ---- paired tests: per-amplicon F1, prFBA vs PICRUSt2 ---------------------- #
    from scipy.stats import wilcoxon
    tests = {}
    for region, rows in sorted(paired.items()):
        for m in ("prfba_include", "prfba_exclude"):
            pairs = [(v[m][0], v["picrust2"][0], v[m][1]) for v in rows.values()
                     if m in v and "picrust2" in v]
            # population: a prFBA abstention already entered `paired` as F1 = 0.
            # conditional: only amplicons where prFBA actually made a prediction.
            for scope, sel_pairs in (
                    ("population", [(a, b) for a, b, _ in pairs]),
                    ("conditional", [(a, b) for a, b, called in pairs if called])):
                if len(sel_pairs) < 10:
                    continue
                x = [a for a, _ in sel_pairs]
                y = [b for _, b in sel_pairs]
                diffs = [a - b for a, b in sel_pairs]
                if sum(1 for d in diffs if d != 0) < 10:
                    continue
                stat, p = wilcoxon(x, y, zero_method="zsplit")
                tests[f"{region}|{m}_vs_picrust2|{scope}"] = {
                    "n": len(x), "n_nonzero_diff": sum(1 for d in diffs if d != 0),
                    "mean_f1_prfba": round(st.mean(x), 4),
                    "mean_f1_picrust2": round(st.mean(y), 4),
                    "median_diff": round(st.median(diffs), 4),
                    "prfba_wins": sum(1 for d in diffs if d > 0),
                    "picrust2_wins": sum(1 for d in diffs if d < 0),
                    "ties": sum(1 for d in diffs if d == 0),
                    "wilcoxon_stat": float(stat), "p_value": float(p),
                }
    # Holm-Bonferroni step-down over the family of paired tests: sort ascending,
    # multiply the k-th smallest by (n - k), then enforce monotonicity.
    n_tests = len(tests)
    running = 0.0
    for k, (name, v) in enumerate(sorted(tests.items(), key=lambda kv: kv[1]["p_value"])):
        running = max(running, min(1.0, v["p_value"] * (n_tests - k)))
        v["p_holm"] = running

    summary = {
        "_meta": {
            "question": "prFBA vs PICRUSt2 functional prediction on the same 10k "
                        "region-validation amplicons, same truth, same projection",
            "n_amplicons": len(truth_meta),
            "n_picrust2_sequences": len(p2_pred),
            "n_amplicons_without_picrust2_prediction": sum(
                1 for k in truth_meta if md5_of.get(k) not in p2_pred),
            "ec_universe_shared": len(shared_ec),
            "rxn_universe_shared": len(shared_rxn),
            "methods": sorted({k[0] for k in aggs}),
            "prfba_selection": {"default": sel.get("_meta", {}).get("select_mode"),
                                "legacy": sel_legacy.get("_meta", {}).get("select_mode")},
            "picrust2": {"version": "2.6.3", "hsp": "mp", "traits": "EC,KO",
                         "reference": "default bacteria+archaea; GTDB-based lineages, "
                                      "release not stated in the shipped metadata",
                         "n_reference_genomes": len(ref_ec),
                         "min_align_excluded_sequences": sum(
                             1 for md5 in fanout if md5 not in p2_pred)},
            "caveats": [
                "PICRUSt2 cannot exclude the source organism; its include-self analogue "
                "is src_in_picrust2_ref=True and its novel analogue src_in_picrust2_ref="
                "False, which is imperfect (GTDB is species-dereplicated).",
                "Database sizes differ: BV-BRC ~112k genomes with PGFams vs PICRUSt2's "
                f"{len(ref_ec)} reference genomes.",
                "Truth is BV-BRC/RAST annotation; see annotation_ceiling.json for the "
                "ceiling this imposes on PICRUSt2.",
                "rxn_role is prFBA's own ModelSEED-role route, reported for context only.",
            ],
        },
        "by_region": {},
        "by_stratum": {},
        "agreement_prfba_include_vs_picrust2_jaccard": {
            f"{space}|{region}": {"n": len(v), "mean": round(st.mean(v), 4),
                                  "median": round(st.median(v), 4)}
            for (space, region), v in sorted(agree.items()) if v},
        "paired_tests_ec_f1": tests,
    }
    for (method, space, region, domain), agg in sorted(aggs.items()):
        summary["by_region"].setdefault(space, {}).setdefault(region, {}) \
            .setdefault(domain, {})[method] = agg.block()
    for (method, space, kind, value), agg in sorted(strat.items()):
        summary["by_stratum"].setdefault(space, {}).setdefault(kind, {}) \
            .setdefault(value, {})[method] = agg.block()

    with open(args.summary, "w") as fh:
        json.dump(summary, fh, indent=1)
    with open(args.routing, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["region", "truth_domain", "picrust2_domain", "n_amplicons"])
        for (region, domain, p2d), n in sorted(routing.items()):
            w.writerow([region, domain, p2d, n])

    ec = summary["by_region"]["ec"]["all"]["all"]
    print(f"[compare] EC space, all regions (conditional / population recall):")
    for method in sorted(ec):
        b = ec[method]
        cond = b["conditional"]["recall"]["value"]
        pop = b["population"]["recall"]["value"]
        print(f"  {method:24} n={b['n']:6} coverage={b['coverage']:.3f} "
              f"recall={cond:.3f} / {pop:.3f} "
              f"precision={b['conditional']['precision']['value']:.3f} "
              f"F1={b['conditional']['f1']['value']:.3f}")
    print(f"[compare] -> {args.summary}, {args.cells}, {args.routing}")


if __name__ == "__main__":
    main()
