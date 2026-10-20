#!/usr/bin/env python
"""Locate every 16S rRNA gene in the MAG-arm source (and donor) genomes.

PICRUSt2 has no genome input: a MAG can only reach it through a 16S gene carried
on one of its contigs. This finds those genes with PICRUSt2's OWN reference
profiles (``bac_ref.hmm`` + ``arc_ref.hmm`` concatenated; domain is unknown for a
MAG, so both are searched) via ``nhmmer``, the same HMMER build PICRUSt2 uses.

Coordinates are recorded against the contig order ``synth_mags.read_fasta``
sees, so ``mag_16s_retention.py`` can replay the synthesizer and decide, for every
synthetic MAG, whether a 16S locus survived fragmentation + contig dropping.

Overlapping bacterial and archaeal hits are one locus: the envelope is the
better-scoring hit's, and both bit scores are kept (their difference is the
domain call).

    source ~/opt/picrust2/activate.sh   # nhmmer
    python mag_16s_loci.py --fasta-dir ~/Documents/prFBA/mag_validation/contigs \
        --fasta-dir ~/Documents/prFBA/mag_validation/donor_contigs \
        --hmm data/mag_arm/hmm/ssu_bac_arc.hmm --out data/mag_arm/loci.jsonl
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

E_MAX = 1e-10        # a full-length 16S scores E ~ 0; this still admits ~150 bp fragments


def contig_table(path):
    """[(name_token, length)] in file order, parsed exactly as synth_mags.read_fasta."""
    out, h, n = [], None, 0
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if h is not None:
                    out.append((h, n))
                h, n = line[1:].strip(), 0
            else:
                n += len(line.strip())
    if h is not None:
        out.append((h, n))
    # nhmmer names a target by its first whitespace-delimited header token
    return [(h.split()[0] if h else "", n) for h, n in out]


def parse_tbl(path):
    hits = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            env_a, env_b = int(f[8]), int(f[9])
            hits.append(dict(target=f[0], hmm=("arc" if f[2].startswith("archaea") else "bac"),
                             hmm_from=int(f[4]), hmm_to=int(f[5]),
                             lo=min(env_a, env_b) - 1, hi=max(env_a, env_b),   # 0-based [lo, hi)
                             strand=f[11], evalue=float(f[12]), score=float(f[13])))
    return hits


def merge(hits):
    """Collapse overlapping bac/arc hits on one contig into loci."""
    loci = []
    for h in sorted(hits, key=lambda h: -h["score"]):
        for L in loci:
            if L["target"] == h["target"] and h["lo"] < L["hi"] and L["lo"] < h["hi"]:
                L[h["hmm"]] = max(L.get(h["hmm"]) or 0.0, h["score"])
                break
        else:
            loci.append(dict(target=h["target"], lo=h["lo"], hi=h["hi"], strand=h["strand"],
                             evalue=h["evalue"], hmm_from=h["hmm_from"], hmm_to=h["hmm_to"],
                             best_hmm=h["hmm"], **{h["hmm"]: h["score"]}))
    return loci


def scan(args):
    path, hmm, tmpdir = args
    gid = os.path.basename(path)[:-4]
    table = contig_table(path)
    index = {}
    for i, (name, _) in enumerate(table):
        index.setdefault(name, []).append(i)
    with tempfile.NamedTemporaryFile(suffix=".tbl", dir=tmpdir, delete=False) as tf:
        tbl = tf.name
    try:
        # --dna: HMMER's alphabet guesser fails on a genome whose first contig opens
        # with a pyrimidine-only low-complexity run (149015.24)
        p = subprocess.run(["nhmmer", "--dna", "--cpu", "1", "--noali", "--tblout", tbl,
                            "-E", "1e-5", hmm, path], stdout=subprocess.DEVNULL,
                           stderr=subprocess.PIPE, text=True)
        if p.returncode:
            return dict(id=gid, error="nhmmer: " + p.stderr.strip()[-300:],
                        contig_lens=[n for _, n in table])
        hits = [h for h in parse_tbl(tbl) if h["evalue"] <= E_MAX]
    finally:
        os.unlink(tbl)
    loci = []
    for L in merge(hits):
        ci = index.get(L.pop("target"), [])
        if len(ci) != 1:                      # duplicate or unknown header: unplaceable
            return dict(id=gid, error="ambiguous contig name", contig_lens=[n for _, n in table])
        L["ci"] = ci[0]
        L["len"] = L["hi"] - L["lo"]
        loci.append(L)
    loci.sort(key=lambda L: (L["ci"], L["lo"]))
    return dict(id=gid, contig_lens=[n for _, n in table], loci=loci)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta-dir", action="append", required=True)
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=60)
    a = ap.parse_args()
    by_id = {}
    for d in a.fasta_dir:                     # a donor also in the panel is one genome
        for p in sorted(glob.glob(os.path.join(d, "*.fna"))):
            by_id.setdefault(os.path.basename(p)[:-4], p)
    paths = [by_id[k] for k in sorted(by_id)]
    done = set()
    if os.path.exists(a.out):
        with open(a.out) as fh:
            done = {json.loads(l)["id"] for l in fh if l.strip()}
    todo = [p for p in paths if os.path.basename(p)[:-4] not in done]
    print(f"[loci] {len(paths)} genomes, {len(done)} cached, {len(todo)} to scan", flush=True)
    tmpdir = os.path.join(os.path.dirname(os.path.abspath(a.out)), "tmp")
    os.makedirs(tmpdir, exist_ok=True)
    with open(a.out, "a") as out, Pool(a.procs) as pool:
        for i, rec in enumerate(pool.imap_unordered(scan, [(p, a.hmm, tmpdir) for p in todo], 8), 1):
            out.write(json.dumps(rec) + "\n")
            if i % 1000 == 0:
                print(f"[loci] {i}/{len(todo)}", flush=True)
    print("[loci] done", flush=True)


if __name__ == "__main__":
    sys.exit(main())
