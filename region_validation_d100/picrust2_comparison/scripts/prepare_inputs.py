#!/usr/bin/env python
"""Write PICRUSt2 inputs for the 10k region-validation amplicons.

PICRUSt2 predicts per *sequence*, so the study FASTA is the md5-deduplicated
amplicon set (``data/combined_amplicons.unique.fasta``, 60,461 sequences; the
fanout map expands each md5 back to its ``genome__region`` amplicons). The
abundance table is a single dummy sample with count 1 per sequence: only the
per-sequence ``combined_*_predicted.tsv.gz`` tables are used downstream, never
the metagenome output.

``--pilot N`` writes a deterministic subset of up to N sequences per
(region, domain) stratum instead of the full set. ``--batches N`` splits the
sequences into N ``batch_XX/`` sub-folders; PICRUSt2 places and predicts every
study sequence independently (per-query EPA-ng placement, per-sequence domain
choice, and unknown tips do not affect castor's parsimony states elsewhere), so
batches only bound run time and checkpoint progress.
"""
import argparse
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RV = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(RV, "scripts"))
import _common as C  # noqa: E402

PATHS = C.PATHS


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pilot", type=int, default=0,
                    help="sequences per (region, domain) stratum; 0 = full set")
    ap.add_argument("--seed", type=int, default=1729)
    ap.add_argument("--batches", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    fanout = json.load(open(PATHS.amplicon_fanout))
    seqs = dict(read_fasta(PATHS.unique_fasta))
    assert set(seqs) == set(fanout), "unique FASTA and fanout map disagree"

    keep = sorted(seqs)
    if args.pilot:
        truth = json.load(open(PATHS.truth_json))
        strata = {}
        for md5 in keep:
            key = fanout[md5][0]
            t = truth[key]
            strata.setdefault((t["region"], t["domain"]), []).append(md5)
        rng = random.Random(args.seed)
        keep = []
        for stratum in sorted(strata):
            members = strata[stratum]
            keep.extend(rng.sample(members, min(args.pilot, len(members))))
        keep.sort()

    for b in range(args.batches):
        part = keep[b::args.batches]
        outdir = (os.path.join(args.outdir, f"batch_{b:02d}") if args.batches > 1
                  else args.outdir)
        os.makedirs(outdir, exist_ok=True)
        fasta = os.path.join(outdir, "study_seqs.fna")
        table = os.path.join(outdir, "study_seqs_abundance.tsv")
        with open(fasta, "w") as fh:
            for md5 in part:
                fh.write(f">{md5}\n{seqs[md5]}\n")
        with open(table, "w") as fh:
            fh.write("sequence\tall_amplicons\n")
            for md5 in part:
                fh.write(f"{md5}\t1\n")
        print(f"[prepare_inputs] {len(part)} sequences -> {fasta}, {table}")


if __name__ == "__main__":
    main()
