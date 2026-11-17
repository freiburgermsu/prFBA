#!/usr/bin/env python
"""D6 contig-cap caveat (prFBA_paper sec:twopanel random-panel MAG arm: "... 10 of the arm's
49,990 MAG queries"): count the synthetic MAGs of the random-panel MAG arm built from the two
genomes whose deposited contigs were fetched at exactly the 25,000-row cap (1733.9197 and
1733.9535; contig_cap_d6.json holds BV-BRC's Content-Range totals for them).

Input: mag_validation_rnd/manifests/b*.json (git-ignored). Writes contig_cap_queries.json here.
"""
import glob, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.dirname(os.path.dirname(HERE))
CAPPED = ("1733.9197", "1733.9535")


def main():
    rows = []
    for p in sorted(glob.glob(os.path.join(R, "mag_validation_rnd", "manifests", "b*.json"))):
        rows += json.load(open(p))["rows"]
    hit = [r for r in rows if r["genome_id"] in CAPPED]
    out = {"what": "random-panel MAG-arm queries built from the two contig-capped genomes",
           "generator": "manuscript_numbers/rescore_sources/count_contig_cap_queries.py",
           "n_mag_queries": len(rows), "n_genomes": len({r["genome_id"] for r in rows}),
           "capped_genomes": list(CAPPED), "n_queries_from_capped": len(hit),
           "queries_from_capped": sorted(r["mag"] for r in hit)}
    with open(os.path.join(HERE, "contig_cap_queries.json"), "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
