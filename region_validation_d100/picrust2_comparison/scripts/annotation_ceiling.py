#!/usr/bin/env python
"""How much of the prFBA-vs-PICRUSt2 gap is annotation disagreement, not prediction?

prFBA's functional truth is BV-BRC/RAST (PGFam products -> EC); PICRUSt2 predicts into
its own reference annotation of GTDB assemblies. For a benchmark genome that is itself
IN PICRUSt2's reference, a perfect predictor placing the amplicon exactly on that tip
still only returns PICRUSt2's row for it. Comparing that row with the BV-BRC truth set
therefore measures the annotation-system ceiling, with no prediction error involved.

Reported on two EC universes:
  full    every EC either system uses
  shared  ECs in PICRUSt2's table columns AND reachable from BV-BRC family products
          (the vocabulary both systems can express; the comparison's primary universe)

Outputs ``--out`` annotation_ceiling.json.
"""
import argparse
import csv
import gzip
import os
import statistics as st
import sys

import orjson

HERE = os.path.dirname(os.path.abspath(__file__))
RV = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(RV, "scripts"))
import _common as C  # noqa: E402

P2 = "/home/freiburger/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files"
DOMAINS = {"bacteria": "bacteria_metadata.csv", "archaea": "archaea_metadata.csv"}


def acc_core(acc):
    """GCA_003722315.1 / GCF_003722315.1 -> 003722315 (same assembly, either accession)."""
    if not acc or "_" not in acc:
        return None
    return acc.split("_", 1)[1].split(".")[0]


def load_reference_ec():
    """{accession core: set(EC)} and the full EC column vocabulary."""
    rows, vocab = {}, set()
    for dom in DOMAINS:
        with gzip.open(os.path.join(P2, dom, "ec.txt.gz"), "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")[1:]
            vocab.update(header)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                core = acc_core(f[0])
                if core is None:
                    continue
                rows[core] = {header[i] for i, v in enumerate(f[1:]) if v not in ("0", "0.0")}
    return rows, vocab


def load_reference_meta():
    """{accession core: ncbi_taxid} for every PICRUSt2 reference genome."""
    out = {}
    for dom, meta in DOMAINS.items():
        with gzip.open(os.path.join(P2, dom, meta + ".gz"), "rt") as fh:
            for r in csv.DictReader(fh):
                core = acc_core(r["accession"])
                if core:
                    out[core] = r.get("ncbi_taxid")
    return out


def prf(truth, pred):
    inter = len(truth & pred)
    return (inter / len(truth) if truth else None,
            inter / len(pred) if pred else None,
            inter / len(truth | pred) if (truth or pred) else None)


def summarize(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    return {"n": len(vals), "mean": round(st.mean(vals), 4),
            "median": round(st.median(vals), 4),
            "p10": round(sorted(vals)[len(vals) // 10], 4),
            "p90": round(sorted(vals)[9 * len(vals) // 10], 4)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--functions", default=os.path.join(HERE, "..", "data",
                                                        "genome_function_sets.json"))
    ap.add_argument("--genomes", default=os.path.join(HERE, "..", "data",
                                                      "bvbrc_source_genomes.json"))
    ap.add_argument("--maps", default=os.path.join(HERE, "..", "data", "function_maps.json"))
    ap.add_argument("--out", default=os.path.join(HERE, "..", "data",
                                                  "annotation_ceiling.json"))
    args = ap.parse_args()

    funcs = orjson.loads(open(args.functions, "rb").read())
    bvbrc = orjson.loads(open(args.genomes, "rb").read())
    maps = orjson.loads(open(args.maps, "rb").read())
    bvbrc_vocab = {e for v in maps["pgfam_ec"].values() for e in v}
    truth_meta = orjson.loads(open(C.PATHS.truth_json, "rb").read())
    sources = {t["src_genome_id"]: t for t in truth_meta.values()}

    ref_ec, ref_vocab = load_reference_ec()
    ref_taxid = load_reference_meta()
    shared_vocab = bvbrc_vocab & ref_vocab
    taxid_index = {}
    for core, tid in ref_taxid.items():
        if tid:
            taxid_index.setdefault(str(tid), []).append(core)

    stats = {"full": {"recall": [], "precision": [], "jaccard": []},
             "shared": {"recall": [], "precision": [], "jaccard": []}}
    n_in_ref = n_taxid_only = n_truth_empty = 0
    per_genome = []
    for gid in sorted(sources):
        acc = (bvbrc.get(gid) or {}).get("assembly_accession")
        core = acc_core(acc)
        in_ref = core in ref_ec if core else False
        if not in_ref:
            tid = str(sources[gid]["src_taxon_id"])
            if tid in taxid_index:
                n_taxid_only += 1
            continue
        n_in_ref += 1
        truth_ec = set(funcs.get(gid, {}).get("ec", ()))
        p2_ec = ref_ec[core]
        if not truth_ec:
            n_truth_empty += 1
            continue
        row = {"genome_id": gid, "accession": acc, "n_truth_ec": len(truth_ec),
               "n_picrust2_ec": len(p2_ec)}
        for universe, keep in (("full", None), ("shared", shared_vocab)):
            t = truth_ec if keep is None else truth_ec & keep
            p = p2_ec if keep is None else p2_ec & keep
            if not t:
                continue
            r, pr, j = prf(t, p)
            for name, v in (("recall", r), ("precision", pr), ("jaccard", j)):
                stats[universe][name].append(v)
            row[universe] = {"recall": r, "precision": pr, "jaccard": j}
        per_genome.append(row)

    out = {
        "_meta": {
            "question": "BV-BRC/RAST EC truth vs PICRUSt2's own reference EC row for the "
                        "SAME genome: annotation-system disagreement, no prediction error",
            "n_benchmark_source_genomes": len(sources),
            "n_source_assembly_in_picrust2_ref": n_in_ref,
            "n_source_not_in_ref_but_taxid_present": n_taxid_only,
            "n_skipped_truth_ec_empty": n_truth_empty,
            "n_picrust2_ref_genomes": len(ref_ec),
            "ec_vocab": {"picrust2_columns": len(ref_vocab),
                         "bvbrc_family_products": len(bvbrc_vocab),
                         "shared": len(shared_vocab)},
            "interpretation": "recall = fraction of BV-BRC truth ECs that PICRUSt2's own "
                              "annotation of the same genome also carries; this bounds any "
                              "PICRUSt2 score against BV-BRC truth",
        },
        "ceiling": {u: {m: summarize(v) for m, v in d.items()} for u, d in stats.items()},
        "per_genome": per_genome,
    }
    with open(args.out, "wb") as fh:
        fh.write(orjson.dumps(out))
    print(orjson.dumps(out["_meta"]).decode())
    print(orjson.dumps(out["ceiling"]).decode())


if __name__ == "__main__":
    main()
