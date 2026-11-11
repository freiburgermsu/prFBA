#!/usr/bin/env python
"""Are the benchmark's "16S" records 16S rRNA genes? Screen them with PICRUSt2's own HMMs.

``region_validation/data/src16s_10k.fasta`` holds one representative 16S per panel
genome, taken from BV-BRC rRNA features. Found while locating 16S genes in the MAG-arm
contigs (``mag_16s_loci.py``): some records are not 16S at all -- 1191523.27's is a
1,920 bp protein-coding gene (it opens ATGCGGAAG...). In-silico PCR amplifies nothing
from such a record, so those genomes silently left the benchmark: all 179 found are
among the 254 panel genomes with no ``benchmark_truth.json`` entry, and no scored cell
comes from one.

Runs ``nhmmer --dna`` with the bacterial + archaeal reference profiles over every
record; a record with no hit at E <= 1e-10 is not a 16S, and one whose longest hit
covers < 80% of it carries a 16S plus a large non-16S remainder.

    PATH=~/opt/picrust2/bin:$PATH ~/Documents/py_venv/bin/python scripts/src16s_hmm_check.py
    # -> data/src16s_hmm_check_summary.json
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mag_common as M  # noqa: E402

E_MAX = 1e-10
COVER_MIN = 0.8


def main():
    fasta = os.path.join(M.RV, "data", "src16s_10k.fasta")
    hmm = os.path.join(M.DATA, "hmm", "ssu_bac_arc.hmm")
    tbl = os.path.join(M.DATA, "src16s_vs_hmm.tbl")
    subprocess.run(["nhmmer", "--dna", "--cpu", "32", "--noali", "--tblout", tbl, "-E", "1e-5",
                    hmm, fasta], check=True, stdout=subprocess.DEVNULL)
    seqlen, name = {}, None
    for line in open(fasta):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split("|")[0]
            seqlen[name] = 0
        else:
            seqlen[name] += len(line)
    best = {}
    for line in open(tbl):
        if line.startswith("#"):
            continue
        f = line.split()
        gid, e = f[0].split("|")[0], float(f[12])
        span = abs(int(f[9]) - int(f[8])) + 1
        if e <= E_MAX and span > best.get(gid, 0):
            best[gid] = span
    not16s = sorted(g for g in seqlen if g not in best)
    low = sorted(g for g in seqlen if g in best and best[g] < COVER_MIN * seqlen[g])
    out = dict(n_records=len(seqlen), e_max=E_MAX, cover_min=COVER_MIN,
               n_not_16s=len(not16s), n_low_cover=len(low), not_16s=not16s, low_cover=low)
    json.dump(out, open(os.path.join(M.P2C, "data", "src16s_hmm_check_summary.json"), "w"), indent=1)
    print(f"{len(seqlen)} records: {len(not16s)} not a 16S (no hit at E <= {E_MAX}); "
          f"{len(low)} with the 16S hit covering < {COVER_MIN:.0%} of the record")


if __name__ == "__main__":
    sys.exit(main())
