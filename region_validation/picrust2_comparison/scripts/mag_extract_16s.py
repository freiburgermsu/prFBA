#!/usr/bin/env python
"""The one way a MAG reaches PICRUSt2: extract its 16S genes, then run the stock pipeline.

For a deterministic sample of panel genomes this regenerates synthetic MAGs with the
MAG arm's own ``synth_mags.py`` (bit-identical to what skani searched), finds 16S genes
in each MAG with ``nhmmer`` and PICRUSt2's own bac+arc profiles (``mag_16s_loci.scan``,
same thresholds as the panel scan), and writes every locus >= MIN_LEN bp as a PICRUSt2
study sequence (minus-strand hits reverse-complemented). Identical sequences are
deduplicated by md5 -- PICRUSt2 predicts per sequence -- and ``loci.json`` maps each
MAG to its loci, their md5s, and whether each came from the source genome or the
contamination donor (``synth_mags`` headers end ``|s<i>`` / ``|x<i>``).

Choosing ONE locus per MAG is left to the scorer (``mag_score.py``), so the choice
rule is explicit there and every copy's prediction stays inspectable.

``--panel`` instead writes every >= MIN_LEN bp locus of every panel source genome
(coordinates from ``mag_16s_loci.py``), so the scorer can give each of the ~140k
replayed MAGs the prediction for the exact locus it keeps, with no proxy.

nhmmer comes from the PICRUSt2 prefix, Python from py_venv (``mag_common`` needs orjson,
which the PICRUSt2 venv lacks):

    export PATH=~/opt/picrust2/bin:$PATH; PY=~/Documents/py_venv/bin/python
    $PY scripts/mag_extract_16s.py --n 400             # then:
    bash scripts/run_picrust2.sh data/mag_arm/sample/inputs data/mag_arm/sample/picrust2 1 60
    $PY scripts/mag_extract_16s.py --panel --batches 4
    TRAITS=EC bash scripts/run_picrust2.sh data/mag_arm/panel/inputs data/mag_arm/panel/picrust2 2 30
"""
import argparse
import hashlib
import json
import os
import random
import subprocess
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mag_common as M  # noqa: E402
from mag_16s_loci import scan  # noqa: E402

PY = os.path.expanduser("~/Documents/py_venv/bin/python")
GRID = [[1.0, 0], [0.5, 0], [0.3, 0], [0.1, 0], [1.0, 0.1], [0.5, 0.1]]
MIN_LEN = 250            # an amplicon-sized piece; PICRUSt2 places 250 bp V4 ASVs routinely
SEED = 1729
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    recs, h, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if h is not None:
                    recs.append((h, "".join(buf)))
                h, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if h is not None:
        recs.append((h, "".join(buf)))
    return recs


def extract(args):
    path, hmm, tmpdir = args
    rec = scan((path, hmm, tmpdir))
    contigs = read_fasta(path)
    out = []
    for L in rec.get("loci", []):
        if L["len"] < MIN_LEN:
            continue
        head, seq = contigs[L["ci"]]
        s = seq[L["lo"]:L["hi"]].upper()
        if L["strand"] == "-":
            s = s.translate(COMP)[::-1]
        out.append(dict(origin="donor" if head.rsplit("|", 1)[-1].startswith("x") else "source",
                        contig=head, lo=L["lo"], hi=L["hi"], strand=L["strand"], len=len(s),
                        best_hmm=L["best_hmm"], bac=L.get("bac"), arc=L.get("arc"),
                        md5=hashlib.md5(s.encode()).hexdigest(), seq=s))
    return os.path.basename(path)[:-4], rec.get("error"), out


def write_inputs(seqs, root, batches):
    """run_picrust2.sh layout: <root>/batch_XX/study_seqs.fna + dummy abundance table."""
    keys = sorted(seqs)
    for b in range(batches):
        d = os.path.join(root, "batch_%02d" % b)
        os.makedirs(d, exist_ok=True)
        part = keys[b::batches]
        with open(os.path.join(d, "study_seqs.fna"), "w") as fh:
            fh.writelines(">%s\n%s\n" % (k, seqs[k]) for k in part)
        with open(os.path.join(d, "study_seqs_abundance.tsv"), "w") as fh:
            fh.write("sequence\tsample\n")
            fh.writelines("%s\t1\n" % k for k in part)


def panel_loci(line):
    rec = json.loads(line)
    gid = rec["id"]
    fasta = os.path.join(M.MV, "contigs", gid + ".fna")
    if not os.path.exists(fasta):
        return gid, None, {}
    contigs = read_fasta(fasta)
    md5s, seqs = [], {}
    for L in rec.get("loci", []):
        if L["len"] < MIN_LEN:
            md5s.append(None)
            continue
        s = contigs[L["ci"]][1][L["lo"]:L["hi"]].upper()
        if L["strand"] == "-":
            s = s.translate(COMP)[::-1]
        h = hashlib.md5(s.encode()).hexdigest()
        md5s.append(h)
        seqs[h] = s
    return gid, md5s, seqs


def panel(a):
    root = os.path.join(M.DATA, "panel")
    os.makedirs(root, exist_ok=True)
    panel_set = set(M.panel_ids())
    lines = [l for l in open(os.path.join(M.DATA, "loci.jsonl"))
             if json.loads(l)["id"] in panel_set]
    md5_map, seqs = {}, {}
    with Pool(a.procs) as pool:
        for gid, md5s, s in pool.imap_unordered(panel_loci, lines, 16):
            if md5s is not None:
                md5_map[gid] = md5s
                seqs.update(s)
    json.dump(md5_map, open(os.path.join(root, "loci_md5.json"), "w"))
    write_inputs(seqs, os.path.join(root, "inputs"), a.batches)
    n = sum(1 for v in md5_map.values() for x in v if x)
    print(f"[extract] panel: {len(md5_map)} genomes, {n} loci >= {MIN_LEN} bp, "
          f"{len(seqs)} unique sequences in {a.batches} batches", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--procs", type=int, default=60)
    ap.add_argument("--outdir", default=os.path.join(M.DATA, "sample"))
    ap.add_argument("--panel", action="store_true")
    ap.add_argument("--batches", type=int, default=4)
    a = ap.parse_args()
    if a.panel:
        return panel(a)
    mags_dir = os.path.join(a.outdir, "mags")
    os.makedirs(mags_dir, exist_ok=True)
    ids = sorted(random.Random(SEED).sample(M.panel_ids(), a.n))
    with open(os.path.join(a.outdir, "ids.txt"), "w") as fh:
        fh.write("\n".join(ids) + "\n")

    manifest = os.path.join(a.outdir, "manifest.json")
    # the MAG FASTAs are disposable (deterministic, ~4.5 GB for 400 genomes): rebuild them
    # whenever any is missing
    if not os.path.exists(manifest) or not all(
            os.path.exists(os.path.join(mags_dir, r["mag"] + ".fna"))
            for r in json.load(open(manifest))["rows"]):
        # one synthesizer process per chunk, as run_batch.sh does
        chunks = [ids[i::16] for i in range(16)]
        procs = []
        for i, ch in enumerate(chunks):
            f = os.path.join(a.outdir, "ids_%02d.txt" % i)
            open(f, "w").write("\n".join(ch) + "\n")
            procs.append(subprocess.Popen(
                [PY, os.path.join(M.MV, "synth_mags.py"), "--contigs", os.path.join(M.MV, "contigs"),
                 "--donor-contigs", os.path.join(M.MV, "donor_contigs"),
                 "--donors", os.path.join(M.MV, "donor_assignment.json"), "--grid", json.dumps(GRID),
                 "--ids", f, "--outdir", mags_dir, "--manifest", f + ".manifest.json"]))
        for p in procs:
            if p.wait():
                sys.exit("synth_mags.py failed")
        rows = []
        for i in range(16):
            f = os.path.join(a.outdir, "ids_%02d.txt.manifest.json" % i)
            rows += json.load(open(f))["rows"]
            os.remove(f)
            os.remove(os.path.join(a.outdir, "ids_%02d.txt" % i))
        json.dump(dict(grid=GRID, seed=SEED, rows=rows), open(manifest, "w"))
    rows = json.load(open(manifest))["rows"]
    print(f"[extract] {len(rows)} MAGs from {len(ids)} genomes", flush=True)

    hmm = os.path.join(M.DATA, "hmm", "ssu_bac_arc.hmm")
    tmpdir = os.path.join(M.DATA, "tmp")
    os.makedirs(tmpdir, exist_ok=True)
    jobs = [(os.path.join(mags_dir, r["mag"] + ".fna"), hmm, tmpdir) for r in rows]
    loci, seqs, errors = {}, {}, {}
    with Pool(a.procs) as pool:
        for mag, err, found in pool.imap_unordered(extract, jobs, 4):
            if err:
                errors[mag] = err
            for L in found:
                seqs[L["md5"]] = L.pop("seq")
            loci[mag] = found
    meta = {r["mag"]: r for r in rows}
    json.dump(dict(min_len=MIN_LEN, errors=errors,
                   mags={m: dict(genome_id=meta[m]["genome_id"],
                                 completeness=meta[m]["completeness_target"],
                                 contamination=meta[m]["contamination_target"],
                                 src_bp=meta[m]["src_bp"], contam_bp=meta[m]["contam_bp"],
                                 loci=loci.get(m, [])) for m in sorted(meta)}),
              open(os.path.join(a.outdir, "loci.json"), "w"), indent=0)

    write_inputs(seqs, os.path.join(a.outdir, "inputs"), 1)
    with_16s = sum(1 for m in loci if loci[m])
    print(f"[extract] {with_16s}/{len(rows)} MAGs carry a >= {MIN_LEN} bp 16S locus; "
          f"{sum(map(len, loci.values()))} loci, {len(seqs)} unique sequences; "
          f"{len(errors)} nhmmer errors", flush=True)


if __name__ == "__main__":
    sys.exit(main())
