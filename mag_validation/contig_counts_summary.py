#!/usr/bin/env python
"""Summarise mag_validation/data/contig_counts.tsv (written by count_contigs.sh) into
mag_validation/data/contig_counts_summary.json: the range of deposited contig counts over the
original panel, the share of single-contig genomes and the share with fewer than ten contigs,
under both counts (all deposited contigs; contigs of at least MIN_FRAG = 2,500 bp).

These back the manuscript's sec:magval sentence on why synth_mags.py fragments first
("... a panel spanning 1 to 15,664 contigs") and synth_mags.py's docstring.
"""
import csv, json, os, statistics

HERE = os.path.dirname(os.path.abspath(__file__))
TSV = os.path.join(HERE, "data", "contig_counts.tsv")
OUT = os.path.join(HERE, "data", "contig_counts_summary.json")


def block(rows, col):
    v = [int(r[col]) for r in rows]
    n = len(v)
    top = sorted(rows, key=lambda r: -int(r[col]))[:3]
    return {"n_genomes": n, "min": min(v), "max": max(v), "median": statistics.median(v),
            "argmax": [r["genome_id"] for r in top if int(r[col]) == max(v)],
            "top3": [[r["genome_id"], int(r[col])] for r in top],
            "n_single": sum(x == 1 for x in v), "share_single": round(sum(x == 1 for x in v) / n, 6),
            "n_lt10": sum(x < 10 for x in v), "share_lt10": round(sum(x < 10 for x in v) / n, 6)}


def main():
    rows = list(csv.DictReader(open(TSV), delimiter="\t"))
    out = {"source": "mag_validation/data/contig_counts.tsv (count_contigs.sh over mag_validation/contigs/*.fna)",
           "all_contigs": block(rows, "n_contigs"),
           "contigs_ge2500": block(rows, "n_contigs_ge2500"),
           "note": ("The pre-2026-11 text printed 19.8% single-contig, 40.4% under ten and a maximum of "
                    "5,265; none of those reproduces under either count over the 9,997 fetched genomes "
                    "(5,265 is 114683.4, the second largest).")}
    json.dump(out, open(OUT, "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
