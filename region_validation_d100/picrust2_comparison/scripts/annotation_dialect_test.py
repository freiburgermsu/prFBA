#!/usr/bin/env python
"""Whose annotation dialect is the benchmark written in?

prFBA's predictions are BV-BRC gene families and the benchmark's truth is the BV-BRC
annotation of the source genome, so prFBA and the truth speak the same annotation
dialect. PICRUSt2 predicts into its own reference's annotation of GTDB assemblies.
``annotation_ceiling.py`` shows the two dialects agree on only ~0.90 recall / ~0.74
precision for the SAME genome, which caps any PICRUSt2 score against our truth.

This runs the mirror image. For the amplicons whose source genome is in PICRUSt2's
reference, both methods are scored twice:

  truth = BV-BRC annotation of the source genome   (prFBA's home dialect)
  truth = PICRUSt2's own annotation of that genome (PICRUSt2's home dialect)

If prFBA leads under both, the advantage survives the dialect; if the ranking flips,
the margin measured in compare.py is mostly the benchmark speaking prFBA's language.
Either way the answer belongs in the paper, because it bounds what the head-to-head
can claim.
"""
import argparse
import json
import os
import statistics as st
import sys

import orjson

HERE = os.path.dirname(os.path.abspath(__file__))
RV = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(RV, "scripts"))
sys.path.insert(0, HERE)
import _common as C  # noqa: E402
from annotation_ceiling import acc_core, load_reference_ec  # noqa: E402
from compare import load_picrust2, metrics  # noqa: E402

MODE_EXCLUDE, MODE_INCLUDE = "exclude_0.987", "include"


def summarize(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    return {"n": len(vals), "mean": round(st.mean(vals), 4),
            "median": round(st.median(vals), 4)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    D = os.path.join(HERE, "..", "data")
    ap.add_argument("--picrust2", default=os.path.join(D, "full_picrust2"))
    ap.add_argument("--selection", default=os.path.join(D, "prfba_selection.json"))
    ap.add_argument("--functions", default=os.path.join(D, "genome_function_sets.json"))
    ap.add_argument("--maps", default=os.path.join(D, "function_maps.json"))
    ap.add_argument("--genomes", default=os.path.join(D, "bvbrc_source_genomes.json"))
    ap.add_argument("--out", default=os.path.join(D, "annotation_dialect.json"))
    args = ap.parse_args()

    truth_meta = orjson.loads(open(C.PATHS.truth_json, "rb").read())
    fanout = orjson.loads(open(C.PATHS.amplicon_fanout, "rb").read())
    md5_of = {k: m for m, ks in fanout.items() for k in ks}
    funcs = orjson.loads(open(args.functions, "rb").read())
    maps = orjson.loads(open(args.maps, "rb").read())
    sel = orjson.loads(open(args.selection, "rb").read())
    bvbrc = orjson.loads(open(args.genomes, "rb").read())

    ec_index = {}
    p2_pred, _p2_meta = load_picrust2(args.picrust2, ec_index)
    inv = {v: k for k, v in ec_index.items()}
    p2_pred = {m: {inv[i] for i in v} for m, v in p2_pred.items()}

    ref_ec, ref_vocab = load_reference_ec()
    bvbrc_vocab = {e for v in maps["pgfam_ec"].values() for e in v}
    shared = bvbrc_vocab & ref_vocab
    g_ec = {g: set(r["ec"]) & shared for g, r in funcs.items() if g != "_meta"}

    # scored twice: over every amplicon (abstentions = 0) and only where prFBA answered,
    # so the dialect effect is not confounded with prFBA's abstention rate
    out = {scope: {m: {t: {"recall": [], "precision": [], "f1": []}
                       for t in ("bvbrc", "picrust2")}
                   for m in ("prfba_exclude", "prfba_include", "picrust2")}
           for scope in ("population", "conditional")}
    n_amplicons = 0
    for key, t in truth_meta.items():
        src = t["src_genome_id"]
        core = acc_core((bvbrc.get(src) or {}).get("assembly_accession"))
        if core is None or core not in ref_ec:
            continue                      # only genomes PICRUSt2 also contains
        truths = {"bvbrc": g_ec.get(src, set()), "picrust2": ref_ec[core] & shared}
        if not truths["bvbrc"] or not truths["picrust2"]:
            continue
        rec = sel.get(key, {})
        preds = {}
        for mode, name in ((MODE_EXCLUDE, "prfba_exclude"), (MODE_INCLUDE, "prfba_include")):
            reps = rec.get(mode, {}).get("reps", [])
            preds[name] = set().union(*[g_ec.get(g, set()) for g in reps]) if reps else set()
        preds["picrust2"] = p2_pred.get(md5_of.get(key), set()) & shared
        n_amplicons += 1
        answered = bool(preds["prfba_exclude"])
        for method, pred in preds.items():
            for tname, T in truths.items():
                scopes = ["population"] + (["conditional"] if answered else [])
                if not pred:
                    for sc in scopes:
                        out[sc][method][tname]["recall"].append(0.0)
                        out[sc][method][tname]["f1"].append(0.0)
                    continue
                _tp, r, p, f1, _j = metrics(T, pred)
                for sc in scopes:
                    out[sc][method][tname]["recall"].append(r)
                    out[sc][method][tname]["precision"].append(p)
                    out[sc][method][tname]["f1"].append(f1 or 0.0)

    result = {
        "_meta": {
            "question": "does prFBA's measured advantage survive scoring against "
                        "PICRUSt2's own annotation instead of BV-BRC's?",
            "n_amplicons": n_amplicons,
            "note": "amplicons whose source genome is in BOTH databases; EC space, "
                    "shared vocabulary; prFBA abstentions scored as 0 (population view)",
        },
        "by_scope": {sc: {m: {t: {k: summarize(v) for k, v in d.items()}
                              for t, d in md.items()} for m, md in sd.items()}
                     for sc, sd in out.items()},
    }
    f = lambda sc, m, t: result["by_scope"][sc][m][t]["f1"]["mean"]  # noqa: E731
    result["_meta"]["prfba_minus_picrust2_f1"] = {
        f"{sc}|truth={t}": round(f(sc, "prfba_exclude", t) - f(sc, "picrust2", t), 4)
        for sc in ("conditional", "population") for t in ("bvbrc", "picrust2")}
    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=1)
    print(json.dumps(result["_meta"], indent=1))
    for sc in ("conditional", "population"):
        n = result["by_scope"][sc]["picrust2"]["bvbrc"]["f1"]["n"]
        print(f"\n=== {sc} (n={n})")
        print(f"{'method':16} {'F1 vs BV-BRC truth':>20} {'F1 vs PICRUSt2 truth':>22}")
        for m in ("picrust2", "prfba_exclude", "prfba_include"):
            print(f"{m:16} {f(sc,m,'bvbrc'):20.3f} {f(sc,m,'picrust2'):22.3f}")


if __name__ == "__main__":
    main()
