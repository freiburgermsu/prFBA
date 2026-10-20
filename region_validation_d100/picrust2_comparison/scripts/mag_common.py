"""Shared paths and loaders for the PICRUSt2-on-MAGs feasibility study (``mag_*.py``).

Everything is scored the way ``compare.py`` scores the amplicon benchmark: truth is the
source genome's BV-BRC PGFam set projected to EC numbers (``genome_function_sets.json``),
and both truth and prediction are restricted to the ECs BOTH annotation systems can
express (BV-BRC family-product ECs intersected with PICRUSt2's EC columns).
"""
import csv
import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
P2C = os.path.abspath(os.path.join(HERE, ".."))
RV = os.path.abspath(os.path.join(P2C, ".."))
PRFBA = os.path.abspath(os.path.join(RV, ".."))
DATA = os.path.join(P2C, "data", "mag_arm")
MV = os.path.join(PRFBA, "mag_validation")
sys.path.insert(0, os.path.join(RV, "scripts"))
sys.path.insert(0, HERE)

from annotation_ceiling import acc_core, load_reference_ec  # noqa: E402,F401
from compare import metrics  # noqa: E402,F401

FUNCTIONS = os.path.join(P2C, "data", "genome_function_sets.json")
MAPS = os.path.join(P2C, "data", "function_maps.json")
PANEL_EC = os.path.join(DATA, "panel_ec.json")


def panel_ids():
    return sorted(f[:-4] for f in os.listdir(os.path.join(MV, "contigs")) if f.endswith(".fna"))


def load_panel_ec():
    """({genome_id: frozenset(EC)} for the 9,997 MAG-arm genomes, shared EC vocabulary).

    Cached: the full ``genome_function_sets.json`` is 1.7 GB, and every genome the MAG
    arm can score or select (truth and reference universe alike) is a panel genome.
    """
    if not os.path.exists(PANEL_EC):
        import orjson
        funcs = orjson.loads(open(FUNCTIONS, "rb").read())
        maps = orjson.loads(open(MAPS, "rb").read())
        _, ref_vocab = load_reference_ec()
        bvbrc_vocab = {e for v in maps["pgfam_ec"].values() for e in v}
        ids = panel_ids()
        out = {"_shared_ec": sorted(bvbrc_vocab & ref_vocab),
               "_missing": [g for g in ids if g not in funcs],
               "genomes": {g: sorted(funcs[g]["ec"]) for g in ids if g in funcs}}
        os.makedirs(DATA, exist_ok=True)
        tmp = PANEL_EC + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(out, fh)
        os.replace(tmp, PANEL_EC)
    d = json.load(open(PANEL_EC))
    return {g: frozenset(v) for g, v in d["genomes"].items()}, frozenset(d["_shared_ec"])


def p2_ec(outdir):
    """{study sequence id: frozenset(EC)} from a PICRUSt2 combined_EC_predicted table."""
    out = {}
    with gzip.open(os.path.join(outdir, "combined_EC_predicted.tsv.gz"), "rt") as fh:
        cols = [c.replace("EC:", "") for c in fh.readline().rstrip("\n").split("\t")[1:]]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = frozenset(c for c, v in zip(cols, f[1:]) if v not in ("0", "0.0"))
    return out


def p2_meta(outdir):
    """{study sequence id: (NSTI, best_domain, closest_reference_genome)}."""
    out = {}
    with gzip.open(os.path.join(outdir, "combined_marker_predicted_and_nsti.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["sequence"]] = (float(r["metadata_NSTI"]), r["best_domain"],
                                  r.get("closest_reference_genome"))
    return out


def fullength_md5():
    """{genome_id: md5 of its FullLength16S amplicon} -- the PICRUSt2 study id in
    ``data/full_picrust2`` whose prediction a retained 16S should reproduce."""
    fan = json.load(open(os.path.join(RV, "data", "amplicon_md5_fanout.json")))
    out = {}
    for md5, keys in fan.items():
        for k in keys:
            gid, region = k.split("__", 1)
            if region == "FullLength16S":
                out[gid] = md5
    return out
