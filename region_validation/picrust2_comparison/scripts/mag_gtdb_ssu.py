#!/usr/bin/env python
"""Real-MAG context for the PICRUSt2-on-MAGs question, from GTDB r232 metadata.

1. How often does a REAL MAG carry a 16S gene, by CheckM2 completeness, against
   isolates? (GTDB's own SSU calls: ``ssu_count``, ``ssu_length`` of the longest.)
   This is the real-world counterpart of ``mag_16s_retention.py``'s synthetic curve,
   which cannot reproduce rRNA-operon mis-assembly or coverage-biased binning.
2. How many of PICRUSt2 2.6.3's reference genomes are themselves MAGs?
3. Does a genome's 16S disagree with its genome-based taxonomy (SILVA hit of the SSU
   vs GTDB lineage) at phylum more often in MAGs than in isolates -- i.e. is the binned
   16S someone else's? Isolate-calibrated: for each GTDB phylum, a SILVA phylum label is
   "consistent" if it names >= CONSISTENT_MIN of that phylum's isolates IN THE SAME
   SOURCE (GTDB r232 uses SILVA 138.2 names, PICRUSt2's metadata older ones), which
   absorbs nomenclature drift. Phyla with < MIN_ISOLATES isolates are not scored, nor
   are SSUs whose SILVA hit names another domain: those are dominated by SILVA entries
   filed under a host (bacterial 16S labelled "Oryza sativa", "Lentinula edodes"), not
   by mis-binning, and are commoner in isolates than in MAGs in PICRUSt2's metadata.

"isolate" here means GTDB's ``ncbi_genome_category`` is empty (not "derived from
metagenome" / "single cell"), and the completeness axis is CheckM2's estimate -- not
the base-pair fraction the synthetic MAGs are cut to.

Inputs: ``data/mag_arm/gtdb/{bac120,ar53}_metadata_r232.tsv.gz`` (GTDB release files)
and PICRUSt2's ``{bacteria,archaea}_metadata.csv.gz``.
Writes ``data/mag_gtdb_ssu_summary.json``.
"""
import gzip
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mag_common as M  # noqa: E402

GTDB = os.path.join(M.DATA, "gtdb")
P2REF = os.path.expanduser("~/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files")
COLS = ["accession", "checkm2_completeness", "checkm2_contamination", "ncbi_genome_category",
        "ssu_count", "ssu_length", "gtdb_taxonomy", "ssu_silva_taxonomy"]
BINS = [50, 60, 70, 80, 90, 95, 100.01]
CONSISTENT_MIN = 0.05
MIN_ISOLATES = 20


def kind(cat):
    return {"derived from metagenome": "MAG", "none": "isolate",
            "derived from single cell": "SAG"}.get(cat, "other")


def rank(tax, i, gtdb):
    """Rank i (0 = domain, 1 = phylum) of a GTDB ('d__X;p__Y') or SILVA ('X;Y') string."""
    if not isinstance(tax, str) or tax == "none":
        return None
    parts = tax.split(";")
    if len(parts) <= i or not parts[i]:
        return None
    v = parts[i]
    if gtdb:
        return v[3:] if len(v) > 3 and v[1:3] == "__" else None
    return v


def annotate(frame):
    frame["gd"] = [rank(t, 0, True) for t in frame["gtdb_taxonomy"]]
    frame["gp"] = [rank(t, 1, True) for t in frame["gtdb_taxonomy"]]
    frame["sd"] = [rank(t, 0, False) for t in frame["ssu_silva_taxonomy"]]
    frame["sp"] = [rank(t, 1, False) for t in frame["ssu_silva_taxonomy"]]


def consistency_map(frame):
    iso = frame[(frame.kind == "isolate") & frame.gp.notna() & frame.sp.notna()]
    out = {}
    for p, d in iso.groupby("gp"):
        if len(d) >= MIN_ISOLATES:
            freq = d.sp.value_counts(normalize=True)
            out[p] = set(freq[freq >= CONSISTENT_MIN].index)
    return out


def disagreement(frame, consistent):
    has = frame[frame.sd.notna() & frame.gd.notna()]
    f = has[has.gp.isin(consistent) & (has.sd == has.gd)]
    bad = sum(sp not in consistent[gp] for gp, sp in zip(f.gp, f.sp))
    return dict(n_with_ssu=int(len(has)), phylum_n_scored=int(len(f)), phylum_mismatch=int(bad),
                phylum_mismatch_rate=round(bad / len(f), 4) if len(f) else None)


def main():
    g = pd.concat([pd.read_csv(os.path.join(GTDB, f), sep="\t", usecols=COLS, low_memory=False)
                   for f in ("bac120_metadata_r232.tsv.gz", "ar53_metadata_r232.tsv.gz")])
    g["kind"] = g["ncbi_genome_category"].map(kind)
    g["ssu_len"] = pd.to_numeric(g["ssu_length"], errors="coerce").fillna(0)
    g["has_ssu"] = g["ssu_count"] > 0
    g["has_ssu_1200"] = g["ssu_len"] >= 1200
    g["cbin"] = pd.cut(g["checkm2_completeness"], BINS, right=False,
                       labels=["50-60", "60-70", "70-80", "80-90", "90-95", "95-100"])
    out = {"source": "GTDB r232 bac120+ar53 metadata", "n_genomes": int(len(g)),
           "by_kind": {}, "by_kind_completeness": {}}
    for k, d in g.groupby("kind"):
        out["by_kind"][k] = dict(n=int(len(d)), has_ssu=round(d.has_ssu.mean(), 4),
                                 has_ssu_ge1200=round(d.has_ssu_1200.mean(), 4))
    for (k, b), d in g.groupby(["kind", "cbin"], observed=True):
        if k in ("MAG", "isolate"):
            out["by_kind_completeness"].setdefault(k, {})[str(b)] = dict(
                n=int(len(d)), has_ssu=round(d.has_ssu.mean(), 4),
                has_ssu_ge1200=round(d.has_ssu_1200.mean(), 4))

    # ---- does the binned 16S belong to the genome? (isolate-calibrated) ----
    annotate(g)
    cons = consistency_map(g)
    out["ssu_vs_genome_gtdb"] = {k: disagreement(g[g.kind == k], cons) for k in ("isolate", "MAG")}

    # ---- PICRUSt2's reference genomes ----
    g["core"] = [M.acc_core(a.split("_", 1)[1]) if isinstance(a, str) and "_" in a else None
                 for a in g["accession"]]
    cat = dict(zip(g["core"], g["kind"]))
    ref = []
    for dom, meta in (("bacteria", "bacteria_metadata.csv.gz"), ("archaea", "archaea_metadata.csv.gz")):
        r = pd.read_csv(os.path.join(P2REF, dom, meta), usecols=["accession", "checkm_completeness",
                                                                 "gtdb_taxonomy", "ssu_silva_taxonomy"])
        r["domain"] = dom
        ref.append(r)
    ref = pd.concat(ref)
    ref["core"] = [M.acc_core(a) for a in ref["accession"]]
    ref["kind"] = [cat.get(c, "not_in_r232") for c in ref["core"]]
    annotate(ref)
    rcons = consistency_map(ref)        # PICRUSt2's metadata carries older SILVA names
    out["picrust2_reference"] = dict(
        n=int(len(ref)), by_kind={k: int(v) for k, v in ref["kind"].value_counts().items()},
        ssu_vs_genome={k: disagreement(ref[ref.kind == k], rcons) for k in ("isolate", "MAG")},
        checkm_completeness_min={k: float(ref.checkm_completeness[ref.kind == k].min())
                                 for k in ("isolate", "MAG")})
    path = os.path.join(M.P2C, "data", "mag_gtdb_ssu_summary.json")
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    sys.exit(main())
