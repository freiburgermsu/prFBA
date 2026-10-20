#!/usr/bin/env python
"""Feed a synthetic MAG to PICRUSt2 2.6.3 as-is: what does it do with contigs?

PICRUSt2's study input is a FASTA of marker-gene sequences plus an abundance table;
there is no genome or MAG mode. This regenerates one MAG with the MAG arm's own
deterministic synthesizer (``mag_validation/synth_mags.py``), writes each contig as a
"study sequence" (abundance 1), and runs the stock pipeline twice:

  1. defaults -- ``place_seqs`` drops any query whose 16S-HMM-aligned length is under
     ``--min_align`` (0.8; ``--help`` prints "default: 0" because the help string
     formats 0.8 with ``%(default)d``) of its own length;
  2. ``--min_align 0`` -- every contig is forced through hmmalign, EPA-ng and castor;
  3. ``place_seqs.py`` alone on the bacterial reference, to separate the placement
     step's own time and memory from the pipeline's.

Summary (``data/mag_direct_feed_summary.json``): exit status and peak memory of each run, then
per contig: carries a 16S (nhmmer, PICRUSt2's own HMMs) or not, NSTI, routed domain,
closest reference, EC count and EC F1 against the source genome; plus the union over
contigs (what the metagenome table sums) and the same genome's FullLength16S amplicon
prediction from the benchmark run.

    python scripts/mag_direct_feed.py [--gid 3048874.6] [--procs 16]
"""
import argparse
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mag_common as M  # noqa: E402

PY = os.path.expanduser("~/Documents/py_venv/bin/python")
ACTIVATE = os.path.expanduser("~/opt/picrust2/activate.sh")
HMM = os.path.join(M.DATA, "hmm", "ssu_bac_arc.hmm")


def sh(cmd, log):
    """Run under the PICRUSt2 environment, GNU time -v appended to the log."""
    subprocess.run(["bash", "-c", f"source {ACTIVATE} && /usr/bin/time -v {cmd} > {log} 2>&1"])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gid", default="3048874.6",
                    help="panel genome; default has one 16S on a 134 kb deposited contig")
    ap.add_argument("--comp", type=float, default=1.0)
    ap.add_argument("--procs", type=int, default=16)
    a = ap.parse_args()
    wd = os.path.join(M.DATA, "direct_feed")
    os.makedirs(os.path.join(wd, "mags"), exist_ok=True)
    tag = "%s_c%d_x0" % (a.gid, round(a.comp * 100))
    mag = os.path.join(wd, "mags", tag + ".fna")
    if not os.path.exists(mag):
        open(os.path.join(wd, "ids.txt"), "w").write(a.gid + "\n")
        subprocess.run([PY, os.path.join(M.MV, "synth_mags.py"),
                        "--contigs", os.path.join(M.MV, "contigs"),
                        "--donor-contigs", os.path.join(M.MV, "donor_contigs"),
                        "--donors", os.path.join(M.MV, "donor_assignment.json"),
                        "--grid", json.dumps([[a.comp, 0]]), "--ids", os.path.join(wd, "ids.txt"),
                        "--outdir", os.path.join(wd, "mags"),
                        "--manifest", os.path.join(wd, "manifest.json")], check=True)

    # one study sequence per contig; PICRUSt2 wants single-token, newick-safe ids
    seqs, h = {}, None
    for line in open(mag):
        line = line.strip()
        if line.startswith(">"):
            h = re.sub(r"[^A-Za-z0-9_]", "_", line[1:])
            seqs[h] = []
        elif line:
            seqs[h].append(line)
    seqs = {k: "".join(v) for k, v in seqs.items()}
    fna, tsv = os.path.join(wd, "contigs.fna"), os.path.join(wd, "contigs_abundance.tsv")
    with open(fna, "w") as fh:
        fh.writelines(">%s\n%s\n" % kv for kv in seqs.items())
    with open(tsv, "w") as fh:
        fh.write("sequence\tmag\n")
        fh.writelines("%s\t1\n" % k for k in seqs)

    runs = {}
    bac_ref = os.path.expanduser("~/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files/"
                                 "bacteria/bac_ref")
    for name, extra in (("default", ""), ("minalign0", "--min_align 0"), ("place_seqs_only", None)):
        out = os.path.join(wd, "out_" + name)
        log = os.path.join(wd, name + ".log")
        if not os.path.exists(log):
            subprocess.run(["rm", "-rf", out])
            if extra is None:       # the placement step alone, to separate its cost
                os.makedirs(out)
                sh(f"place_seqs.py --study_fasta {fna} --ref_dir {bac_ref} --out_tree {out}/bac.tre "
                   f"--processes {a.procs} --intermediate {out}/int --verbose", log)
            else:
                sh(f"picrust2_pipeline.py -s {fna} -i {tsv} -o {out} -p {a.procs} "
                   f"--in_traits EC --no_pathways {extra} --verbose", log)
        txt = open(log).read()
        rc = re.search(r"Exit status: (\d+)", txt)
        rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt)
        wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", txt)
        stop = re.search(r"^(Stopping - .*)$", txt, re.M)
        runs[name] = dict(exit=int(rc.group(1)) if rc else None,
                          peak_rss_gb=round(int(rss.group(1)) / 1e6, 1) if rss else None,
                          wall=wall.group(1) if wall else None,
                          message=stop.group(1) if stop else None)

    # which contigs carry a 16S (same HMMs, same thresholds as mag_16s_loci.py)
    tbl = os.path.join(wd, "contigs_16s.tbl")
    subprocess.run(["bash", "-c", f"source {ACTIVATE} && nhmmer --dna --cpu 4 --noali "
                    f"--tblout {tbl} -E 1e-10 {HMM} {fna} > /dev/null"], check=True)
    has16s = {l.split()[0] for l in open(tbl) if not l.startswith("#")}

    ec, shared = M.load_panel_ec()
    truth = ec.get(a.gid, frozenset()) & shared
    res = dict(gid=a.gid, mag=tag, n_contigs=len(seqs),
               contig_len={"min": min(map(len, seqs.values())),
                           "max": max(map(len, seqs.values())),
                           "total": sum(map(len, seqs.values()))},
               contigs_with_16s=sorted(has16s), n_truth_ec=len(truth), runs=runs)
    out = os.path.join(wd, "out_minalign0")
    if runs["minalign0"]["exit"] == 0:
        pred, meta = M.p2_ec(out), M.p2_meta(out)
        rows = []
        for cid in sorted(pred):
            p = pred[cid] & shared
            _, r, pr, f1, _ = M.metrics(truth, p)
            nsti, dom, closest = meta[cid]
            rows.append(dict(contig=cid, length=len(seqs[cid]), has_16s=cid in has16s,
                             nsti=nsti, domain=dom, closest=closest, n_ec=len(p),
                             recall=r, precision=pr, f1=f1))
        res["per_contig"] = rows
        union = frozenset().union(*(pred[c] for c in pred)) & shared
        _, r, pr, f1, _ = M.metrics(truth, union)
        res["union_over_contigs"] = dict(n_ec=len(union), recall=r, precision=pr, f1=f1)
        res["domains"] = {d: sum(1 for x in rows if x["domain"] == d) for d in ("bac", "arc")}
        res["distinct_closest_refs"] = len({x["closest"] for x in rows})
        no16 = [x for x in rows if not x["has_16s"]]
        res["no16s_contigs"] = dict(
            n=len(no16), nsti_min=min(x["nsti"] for x in no16), nsti_max=max(x["nsti"] for x in no16),
            nsti_lt_0_15=sum(1 for x in no16 if x["nsti"] < 0.15),
            f1_mean=round(sum(x["f1"] or 0 for x in no16) / len(no16), 4))
    md5 = M.fullength_md5().get(a.gid)
    for b in sorted(os.listdir(os.path.join(M.P2C, "data", "full_picrust2"))):
        d = os.path.join(M.P2C, "data", "full_picrust2", b)
        if md5 and b.startswith("batch_") and os.path.isdir(d):
            pr = M.p2_ec(d)
            if md5 in pr:
                p = pr[md5] & shared
                _, r, prc, f1, _ = M.metrics(truth, p)
                res["fulllength16s_amplicon"] = dict(md5=md5, n_ec=len(p), recall=r,
                                                     precision=prc, f1=f1,
                                                     nsti=M.p2_meta(d)[md5][0])
                break
    json.dump(res, open(os.path.join(M.P2C, "data", "mag_direct_feed_summary.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "per_contig"}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
