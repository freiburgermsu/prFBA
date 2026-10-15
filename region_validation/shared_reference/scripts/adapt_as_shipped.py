#!/usr/bin/env python
"""Adapter: express the AS-SHIPPED prFBA / PICRUSt2 head-to-head in the
shared-reference fold-artefact format, so `phenotype_fold.py` can run against it
unchanged.

WHY THIS EXISTS
---------------
SPEC Section 9's phenotype arm consumes the k-fold runners' artefacts
(`panel_eligible.json`, `panel_ec_universe.json`, `fold_plan.json`, `cells.csv.gz`,
`prfba/*/predictions.json`, `picrust2/*/predictions.json`).  Those do not exist yet:
the shared-reference reference builds and fold runs have not been executed.  The
phenotype arm's INPUT CONTRACT, however, is just "an EC set per (genome, arm)", and
two such sets are already on disk from the as-shipped comparison:

  truth      picrust2_comparison/data/genome_function_sets.json -> [genome]["ec"]
  prFBA      picrust2_comparison/data/prfba_selection.json -> [cell]["exclude_0.987"]
             ["reps"], unioned through genome_function_sets
  PICRUSt2   picrust2_comparison/data/full_picrust2/batch_*/combined_EC_predicted.tsv.gz
             (same loader semantics as compare.load_picrust2: value > 0 == predicted)

This script writes those into the fold-artefact shapes.  The phenotype harness
therefore accepts EITHER source with no code change; `meta.prediction_source` records
which one was used, and every output JSON produced from this fixture carries
`"prediction_source": "as_shipped"`.

WHAT IS NOT CONTROLLED HERE (read before quoting any number from this fixture)
  * The as-shipped comparison is the confounded one: BV-BRC's annotation is both the
    truth and the universe prFBA's predictions are drawn from (prFBA +0.050 F1 under
    BV-BRC truth, -0.204 under PICRUSt2's -- data/annotation_dialect.json), and the
    two methods' reference databases differ in size (~112k vs 27,870 genomes).
  * prFBA's regime is `exclude_0.987`, which carries the k2=20 truncation artifact
    (94.6% of its self-exclusion abstentions are an emptied truncated hit record).
    The `include` regime is tautological (F1 0.990) and is NOT used as a prediction.
  * There are no folds here.  `fold_plan.json` declares ONE pseudo-fold (fold 0) whose
    `test_genomes` are the fixture genomes and whose `reference_genomes` is EMPTY, so
    the reference-dependent arms (base_constant, base_nn1 via rank1, base_family_donor)
    are unavailable; run with `--arms truth,prfba,picrust2`.

EC UNIVERSE
  The 2,236-EC BV-BRC panel universe (kfold_audit/panel_ec_universe.json, the union of
  complete ECs over the 9,762 eligible panel genomes) is used verbatim, so the fixture's
  scoring vocabulary does not depend on which genomes the fixture happens to contain.
  The fraction of each method's predicted ECs falling OUTSIDE it is reported in
  `meta.out_of_universe`, because that trim is applied identically to all arms by
  phenotype_fold.build_jobs and would otherwise be invisible.
"""
import argparse
import csv
import glob
import gzip
import io
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RV = os.path.join(REPO, "region_validation")
P2D = os.path.join(RV, "picrust2_comparison", "data")

CELL_FIELDS = ["amplicon", "md5", "region", "domain", "src_genome_id", "fold",
               "group_id", "density", "arm", "space", "outcome", "n_truth", "n_pred",
               "tp", "recall", "precision", "f1", "jaccard", "mcc", "n_reps", "tier",
               "abstained", "nsti", "nsti_bin", "picrust2_domain", "routing",
               "rank1_identity", "rank1_rank", "ident_bin", "novelty_bin",
               "n_truth_accessory", "f1_accessory", "n_truth_impossible",
               "template_visible_frac", "min_align_dropped", "unscoreable_reason",
               # appended in lockstep with score_folds.CELL_FIELDS: score_folds
               # REFUSES a cells table without these (a table that cannot be audited
               # for completeness is exactly what C14 exists to reject)
               "abstain_reason", "n_truth_impossible_own_domain"]

# copied from phenotype_arm.GRAMPOS_PHYLA so the adapter never imports modelseedpy
GRAMPOS_PHYLA = {"Bacillota", "Bacillota_A", "Bacillota_B", "Bacillota_C", "Bacillota_G",
                 "Firmicutes", "Actinomycetota", "Actinobacteria", "Chloroflexota",
                 "Deinococcota", "Thermotogota", "Mycoplasmatota", "Tenericutes"}

FROZEN = {"tau_amplicon": 0.987, "k_folds": 1, "seed": 1729,
          "prediction_source": "as_shipped",
          "prfba": {"regime": "exclude_0.987", "select_mode_label": "exact_tie"},
          "picrust2": {"hsp_method": "mp", "edge_exponent": 0.5, "min_align": 0.80,
                       "placement_tool": "epa-ng"},
          "note": "as-shipped fixture; NOT the shared-reference k-fold design"}


def jload(path):
    import orjson
    with open(path, "rb") as fh:
        return orjson.loads(fh.read())


def norm_ec(values):
    out = set()
    for v in values or ():
        s = str(v).strip()
        if s.upper().startswith("EC:"):
            s = s[3:].strip()
        p = s.split(".")
        if len(p) == 4 and all(x.isdigit() for x in p):
            out.add(s)
    return out


def template_for(rec):
    if (rec.get("domain") or "").lower().startswith("arch"):
        return "archaea"
    return "grampos" if rec.get("phylum") in GRAMPOS_PHYLA else "gramneg"


def dom_key(domain):
    return "arc" if (domain or "").lower().startswith("arch") else "bac"


def load_comparison_cells(path, methods, space="ec"):
    """-> {(amplicon, method): {"outcome":.., "f1":.., "n_truth":.., "n_reps":..,
                                "tier":.., "nsti":..}}"""
    out = {}
    with gzip.open(path, "rt", newline="") as fh:
        for r in csv.DictReader(fh):
            if r.get("space") != space or r.get("method") not in methods:
                continue
            out[(r["amplicon"], r["method"])] = r
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-candidates", type=int, default=60,
                    help="genomes written into the fixture; phenotype_fold.py then "
                         "stratifies down to --n-genomes")
    ap.add_argument("--region", default="V4")
    ap.add_argument("--min-truth-ec", type=int, default=200)
    ap.add_argument("--min-archaea", type=int, default=6)
    ap.add_argument("--select-mode-label", default="exact_tie",
                    help="arm label written into cells.csv.gz / predictions.json; must "
                         "match phenotype_fold.py --select-mode")
    ap.add_argument("--selection", default=os.path.join(RV, "data", "selection_10k.json"))
    ap.add_argument("--truth", default=os.path.join(RV, "data", "benchmark_truth.json"))
    ap.add_argument("--fanout", default=os.path.join(RV, "data",
                                                     "amplicon_md5_fanout.json"))
    ap.add_argument("--prfba-selection", default=os.path.join(P2D,
                                                             "prfba_selection.json"))
    ap.add_argument("--comparison-cells", default=os.path.join(P2D,
                                                              "comparison_cells.csv.gz"))
    ap.add_argument("--picrust2", default=os.path.join(P2D, "full_picrust2"))
    ap.add_argument("--functions", default=os.path.join(P2D,
                                                        "genome_function_sets.json"))
    ap.add_argument("--ec-universe", default=os.path.join(P2D, "kfold_audit",
                                                          "panel_ec_universe.json"))
    ap.add_argument("--seed", type=int, default=1729)
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    print("[adapt] loading inputs", flush=True)
    sel = {r["genome_id"]: r for r in jload(args.selection)["genomes"]}
    truth_cells = jload(args.truth)
    fanout = jload(args.fanout)
    prfsel = jload(args.prfba_selection)
    uni_raw = jload(args.ec_universe)
    universe = norm_ec(uni_raw if isinstance(uni_raw, list) else list(uni_raw))
    print("[adapt] selection={} truth_cells={} universe={}".format(
        len(sel), len(truth_cells), len(universe)), flush=True)

    cell_md5 = {}                                   # cell -> amplicon md5
    for md5, cells in fanout.items():
        for c in cells:
            cell_md5[c] = md5

    prf_arm = "prfba_" + args.select_mode_label
    cmp_cells = load_comparison_cells(args.comparison_cells,
                                      {"prfba_exclude", "picrust2"})
    print("[adapt] comparison cells (ec space, 2 methods) = {}".format(len(cmp_cells)),
          flush=True)

    # ---- candidate cells: one region, both methods scored, prFBA answered
    cands = []
    for cell, t in truth_cells.items():
        if t.get("region") != args.region:
            continue
        gid = t["src_genome_id"]
        if gid not in sel:
            continue
        md5 = cell_md5.get(cell)
        if not md5:
            continue
        p = prfsel.get(cell)
        if not p:
            continue
        ex = p.get("exclude_0.987") or {}
        if not ex.get("reps"):
            continue
        cp = cmp_cells.get((cell, "prfba_exclude"))
        c2 = cmp_cells.get((cell, "picrust2"))
        if not cp or not c2 or cp["outcome"] != "ok" or c2["outcome"] != "ok":
            continue
        try:
            n_truth = int(float(cp["n_truth"]))
        except (TypeError, ValueError):
            continue
        if n_truth < args.min_truth_ec:
            continue
        cands.append({"cell": cell, "genome_id": gid, "md5": md5,
                      "region": args.region, "domain": t.get("domain"),
                      "reps": ex["reps"], "tier": ex.get("tier"),
                      "best_identity": ex.get("best_identity"),
                      "prf_row": cp, "p2_row": c2, "n_truth": n_truth})
    print("[adapt] candidate cells = {}".format(len(cands)), flush=True)
    if not cands:
        raise SystemExit("[adapt] no candidates -- check --region / --min-truth-ec")

    # ---- pick a spread: archaea quota first, then |dF1| deciles over bacteria so the
    # fixture contains both agreeing and disagreeing organisms (the phenotype arm's
    # whole question is whether an EC-level difference propagates)
    def f1(row):
        try:
            return float(row["f1"])
        except (TypeError, ValueError):
            return None

    for c in cands:
        a, b = f1(c["prf_row"]), f1(c["p2_row"])
        c["f1_prfba"], c["f1_picrust2"] = a, b
        c["delta_f1"] = None if (a is None or b is None) else a - b
    cands = [c for c in cands if c["delta_f1"] is not None]
    arc = sorted((c for c in cands if dom_key(c["domain"]) == "arc"),
                 key=lambda c: (-abs(c["delta_f1"]), c["genome_id"]))
    bac = sorted((c for c in cands if dom_key(c["domain"]) == "bac"),
                 key=lambda c: (c["delta_f1"], c["genome_id"]))
    chosen = arc[:args.min_archaea]
    n_bac = max(0, args.n_candidates - len(chosen))
    if n_bac and bac:
        step = max(1, len(bac) // n_bac)
        chosen += [bac[i] for i in range(0, len(bac), step)][:n_bac]
    chosen.sort(key=lambda c: c["genome_id"])
    print("[adapt] fixture genomes = {} ({} archaea)".format(
        len(chosen), sum(1 for c in chosen if dom_key(c["domain"]) == "arc")),
        flush=True)

    # ---- PICRUSt2 predicted EC per md5 (compare.load_picrust2 semantics: value > 0)
    want_md5 = {c["md5"] for c in chosen}
    p2_pred, p2_meta_rows = {}, {}
    for d in sorted(glob.glob(os.path.join(args.picrust2, "batch_*"))):
        ecf = os.path.join(d, "combined_EC_predicted.tsv.gz")
        mkf = os.path.join(d, "combined_marker_predicted_and_nsti.tsv.gz")
        if not os.path.exists(ecf):
            continue
        with gzip.open(ecf, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")[1:]
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] not in want_md5:
                    continue
                ecs = sorted({hdr[i].replace("EC:", "")
                              for i, v in enumerate(f[1:]) if float(v) > 0})
                p2_pred[f[0]] = ecs
        if os.path.exists(mkf):
            with gzip.open(mkf, "rt", newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r["sequence"] in want_md5:
                        p2_meta_rows[r["sequence"]] = r
    print("[adapt] picrust2 predictions found for {}/{} md5s".format(
        len(p2_pred), len(want_md5)), flush=True)

    # ---- EC sets for the fixture genomes AND for every prFBA rep
    want_genomes = {c["genome_id"] for c in chosen}
    for c in chosen:
        want_genomes |= set(c["reps"])
    print("[adapt] streaming genome_function_sets.json for {} genomes".format(
        len(want_genomes)), flush=True)
    import ijson
    fsub = {}
    with open(args.functions, "rb") as fh:
        for gid, rec in ijson.kvitems(fh, ""):
            if gid in want_genomes:
                fsub[gid] = {"ec": sorted(norm_ec(rec.get("ec"))),
                             "n_pgfam": rec.get("n_pgfam")}
                if len(fsub) == len(want_genomes):
                    break
    print("[adapt] got EC sets for {}/{}".format(len(fsub), len(want_genomes)),
          flush=True)

    chosen = [c for c in chosen
              if c["genome_id"] in fsub and c["md5"] in p2_pred
              and all(r in fsub for r in c["reps"])]
    print("[adapt] fixture genomes after EC/prediction join = {}".format(len(chosen)),
          flush=True)

    # ---- md5 twins.  phenotype_fold.py keys BOTH methods' predictions by amplicon md5
    # (predictions.json:pred.<md5>), so two fixture genomes carrying a byte-identical
    # amplicon write to the same key.  That is correct ONLY if their predictions agree:
    # prFBA's self-exclusion is applied relative to the SOURCE genome, so in principle
    # two twins can get different reps and the second would silently overwrite the
    # first.  Assert the agreement; drop the extra genomes of any md5 where it fails.
    twins, dropped_twins = {}, []
    for c in chosen:
        twins.setdefault(c["md5"], []).append(c)
    keep = []
    for md5, group in sorted(twins.items()):
        if len(group) == 1:
            keep.append(group[0])
            continue
        ref = (tuple(sorted(group[0]["reps"])), tuple(sorted(p2_pred[md5])))
        agree = all((tuple(sorted(g["reps"])), tuple(sorted(p2_pred[md5]))) == ref
                    for g in group)
        if agree:
            keep.extend(group)
        else:
            keep.append(group[0])
            dropped_twins.extend(g["genome_id"] for g in group[1:])
    chosen = sorted(keep, key=lambda c: c["genome_id"])
    n_twin_md5 = sum(1 for g in twins.values() if len(g) > 1)
    print("[adapt] md5 twins: {} md5s carried by >1 fixture genome, {} genomes dropped "
          "for disagreeing predictions".format(n_twin_md5, len(dropped_twins)),
          flush=True)

    # ---- artefacts
    panel = {"meta": {"frozen": FROZEN, "prediction_source": "as_shipped",
                      "n_panel": len(sel), "n_eligible": len(chosen),
                      "script": "adapt_as_shipped.py"},
             "genomes": {}, "dropped": {}}
    for c in chosen:
        r = sel[c["genome_id"]]
        panel["genomes"][c["genome_id"]] = {
            "label": "g_" + c["genome_id"], "taxon_id": r.get("taxon_id"),
            "md5": r.get("best_16s_md5"), "len": r.get("best_16s_len"),
            "n_full_length_16s": r.get("n_full_length_16s"),
            "domain": r.get("domain"), "phylum": r.get("phylum"),
            "class": r.get("class"), "order": r.get("order"),
            "family": r.get("family"), "genus": r.get("genus"),
            "species": r.get("species"),
            "n_ec": len(fsub[c["genome_id"]]["ec"]),
            "template": template_for(r)}
    with open(os.path.join(args.outdir, "panel_eligible.json"), "w") as fh:
        json.dump(panel, fh, indent=1, sort_keys=True)

    prev = {}
    for c in chosen:
        for e in set(fsub[c["genome_id"]]["ec"]) & universe:
            prev[e] = prev.get(e, 0) + 1
    ecs = sorted(universe)
    with open(os.path.join(args.outdir, "panel_ec_universe.json"), "w") as fh:
        json.dump({"meta": {"n_genomes": len(chosen), "n_ec": len(ecs),
                            "source": args.ec_universe}, "ec": ecs,
                   "prevalence": {e: round(prev.get(e, 0) / max(1, len(chosen)), 6)
                                  for e in ecs},
                   "ec_index": {e: i for i, e in enumerate(ecs)}}, fh, sort_keys=True)

    fold_plan = {"meta": {"frozen": FROZEN, "prediction_source": "as_shipped",
                          "rank_blocked": "none (as-shipped self-exclusion at 0.987)",
                          "k": 1, "n_eligible": len(chosen),
                          "script": "adapt_as_shipped.py"},
                 "groups": {c["genome_id"]: [c["genome_id"]] for c in chosen},
                 "genome_group": {c["genome_id"]: c["genome_id"] for c in chosen},
                 "folds": [{"fold": 0,
                            "test_groups": [c["genome_id"] for c in chosen],
                            "test_genomes": [c["genome_id"] for c in chosen],
                            "reference_genomes": [],
                            "test_cells": [c["cell"] for c in chosen],
                            "test_md5s": sorted({c["md5"] for c in chosen}),
                            "n_test_genomes": len(chosen),
                            "n_reference_genomes": 0}],
                 "leakage_preflight": {"fold_0": {
                     "note": "not applicable: as-shipped run, self-exclusion at 0.987 "
                             "is applied inside prFBA's own hit records"}}}
    with open(os.path.join(args.outdir, "fold_plan.json"), "w") as fh:
        json.dump(fold_plan, fh, indent=1, sort_keys=True)

    with open(os.path.join(args.outdir, "functions_subset.json"), "w") as fh:
        json.dump(fsub, fh, sort_keys=True)

    # cells.csv.gz -- the 36-column contract; only the columns phenotype_fold reads
    # are populated, the rest are blank (as score_folds.py itself writes them)
    oou = {"prfba": [0, 0], "picrust2": [0, 0]}
    raw = gzip.GzipFile(os.path.join(args.outdir, "cells.csv.gz"), "wb", mtime=0)
    with io.TextIOWrapper(raw, encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CELL_FIELDS)
        w.writeheader()
        for c in chosen:
            t = set(fsub[c["genome_id"]]["ec"]) & universe
            arms = {prf_arm: set().union(*[set(fsub[r]["ec"]) for r in c["reps"]]),
                    "picrust2": norm_ec(p2_pred[c["md5"]])}
            oou["prfba"][0] += len(arms[prf_arm] - universe)
            oou["prfba"][1] += len(arms[prf_arm])
            oou["picrust2"][0] += len(arms["picrust2"] - universe)
            oou["picrust2"][1] += len(arms["picrust2"])
            nsti = (p2_meta_rows.get(c["md5"]) or {}).get("metadata_NSTI", "")
            for arm, src in ((prf_arm, c["prf_row"]), ("picrust2", c["p2_row"])):
                p = arms[arm] & universe
                tp = len(t & p)
                row = {k: "" for k in CELL_FIELDS}
                row.update({
                    "amplicon": c["cell"], "md5": c["md5"], "region": c["region"],
                    "domain": c["domain"], "src_genome_id": c["genome_id"],
                    "fold": 0, "group_id": c["genome_id"], "density": 1.0,
                    "arm": arm, "space": "ec", "outcome": "scored",
                    "n_truth": len(t), "n_pred": len(p), "tp": tp,
                    "recall": round(tp / len(t), 6) if t else "",
                    "precision": round(tp / len(p), 6) if p else "",
                    "f1": round(2 * tp / (len(t) + len(p)), 6) if (t or p) else "",
                    "jaccard": round(tp / len(t | p), 6) if (t | p) else "",
                    "n_reps": len(c["reps"]) if arm == prf_arm else "",
                    "tier": c["tier"] if arm == prf_arm else "",
                    "abstained": 0,
                    "nsti": nsti if arm == "picrust2" else "",
                    "picrust2_domain": (p2_meta_rows.get(c["md5"]) or {}).get(
                        "best_domain", "") if arm == "picrust2" else "",
                    "routing": "as_shipped_two_reference" if arm == "picrust2" else "",
                    "rank1_identity": c["best_identity"] if arm == prf_arm else "",
                    "unscoreable_reason": ""})
                w.writerow(row)

    # prFBA predictions.json (one pseudo-fold, density 1.0, not forced)
    pdir = os.path.join(args.outdir, "prfba", "f00_d100")
    os.makedirs(pdir, exist_ok=True)
    pred = {}
    for c in chosen:
        ec = sorted(set().union(*[set(fsub[r]["ec"]) for r in c["reps"]]) & universe)
        pred[c["md5"]] = {
            "ec": ec, "reps": sorted(c["reps"]), "n_reps": len(c["reps"]),
            "tier": c["tier"], "abstained": False, "abstain_reason": "",
            "best_identity": c["best_identity"],
            "rank1_genome": sorted(c["reps"])[0],
            "rank1_identity": c["best_identity"], "rank1_rank": ""}
    with open(os.path.join(pdir, "predictions.json"), "w") as fh:
        json.dump({"meta": {"fold": 0, "density": 1.0, "mode": "as_shipped",
                            "prediction_source": "as_shipped",
                            "forced_answer": False, "frozen": FROZEN,
                            "regime": "exclude_0.987",
                            "n_reference_md5": None,
                            "source": os.path.relpath(args.prfba_selection, REPO),
                            "script": "adapt_as_shipped.py"},
                   "pred": {prf_arm.replace("prfba_", ""): pred}}, fh, sort_keys=True)

    # PICRUSt2 predictions.json, filed by the cell's TRUTH domain (what build_jobs keys on)
    by_dom = {}
    for c in chosen:
        by_dom.setdefault(dom_key(c["domain"]), {})[c["md5"]] = {
            "ec": sorted(norm_ec(p2_pred[c["md5"]]) & universe),
            "nsti": float((p2_meta_rows.get(c["md5"]) or {}).get("metadata_NSTI", 0) or 0),
            "best_domain": (p2_meta_rows.get(c["md5"]) or {}).get("best_domain", ""),
            "closest_reference_genome": (p2_meta_rows.get(c["md5"]) or {}).get(
                "closest_reference_genome", "")}
    for dom, block in by_dom.items():
        qdir = os.path.join(args.outdir, "picrust2", "f00_{}_d100".format(dom))
        os.makedirs(qdir, exist_ok=True)
        with open(os.path.join(qdir, "predictions.json"), "w") as fh:
            json.dump({"meta": {"fold": 0, "domain": dom, "density": 1.0,
                                "mode": "fold", "routing": "forced",
                                "prediction_source": "as_shipped",
                                "ref_name": "picrust2 stock default_files (as shipped)",
                                "n_query": len(block), "frozen": FROZEN,
                                "picrust2_version": "2.6.3",
                                "source": os.path.relpath(args.picrust2, REPO),
                                "script": "adapt_as_shipped.py"},
                       "pred": block}, fh, sort_keys=True)

    meta = {"prediction_source": "as_shipped", "n_genomes": len(chosen),
            "region": args.region, "prfba_arm": prf_arm,
            "prfba_regime": "exclude_0.987",
            "n_ec_universe": len(universe),
            "ec_universe_source": args.ec_universe,
            "n_twin_md5": n_twin_md5,
            "twin_md5_genomes": {m: sorted(g["genome_id"] for g in grp)
                                 for m, grp in sorted(twins.items()) if len(grp) > 1},
            "dropped_twin_genomes": sorted(dropped_twins),
            "out_of_universe": {k: {"n_outside": v[0], "n_total": v[1],
                                    "frac": round(v[0] / v[1], 6) if v[1] else None}
                                for k, v in oou.items()},
            "domains": {d: sum(1 for c in chosen if dom_key(c["domain"]) == d)
                        for d in ("bac", "arc")},
            "templates": {t: sum(1 for c in chosen
                                 if panel["genomes"][c["genome_id"]]["template"] == t)
                          for t in ("gramneg", "grampos", "archaea")},
            "delta_f1_range": [min(c["delta_f1"] for c in chosen),
                               max(c["delta_f1"] for c in chosen)],
            "genomes": [{"genome_id": c["genome_id"], "cell": c["cell"],
                         "md5": c["md5"], "domain": c["domain"],
                         "template": panel["genomes"][c["genome_id"]]["template"],
                         "n_truth_ec": len(set(fsub[c["genome_id"]]["ec"]) & universe),
                         "n_prfba_ec": len(set().union(
                             *[set(fsub[r]["ec"]) for r in c["reps"]]) & universe),
                         "n_picrust2_ec": len(norm_ec(p2_pred[c["md5"]]) & universe),
                         "tier": c["tier"], "n_reps": len(c["reps"]),
                         "f1_prfba": c["f1_prfba"], "f1_picrust2": c["f1_picrust2"],
                         "delta_f1": round(c["delta_f1"], 6)} for c in chosen]}
    with open(os.path.join(args.outdir, "fixture_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
    print("[adapt] wrote fixture to {}".format(args.outdir), flush=True)
    print("[adapt] out-of-universe predicted ECs: prfba {} picrust2 {}".format(
        meta["out_of_universe"]["prfba"], meta["out_of_universe"]["picrust2"]),
        flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
