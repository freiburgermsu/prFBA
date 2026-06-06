#!/usr/bin/env python
"""
true_top5_in_cosine100.py — does the cosine top-100 contain the TRUE top-5 alignment hits?

For a sample of ASVs we determine the "true" top-5 references by Biopython alignment %identity over the
WHOLE 97,624-insert reference space (edlib edit-distance prefilter -> Biopython local re-scoring), then
ask what fraction of those true top-5 fall inside the ASV's cosine top-100 hits. Reports the % of ASVs
whose cosine top-100 contains ALL 5 (and the mean fraction captured).
"""
from __future__ import annotations
import json, time, argparse
from pathlib import Path
import numpy as np, pandas as pd, edlib
from Bio.Align import PairwiseAligner

REPO = Path("/home/freiburger/Documents/codiffusion_bioreactor")
PRFBA = Path("/home/freiburger/Documents/prFBA")
STORE = PRFBA / "v4v5_store"

aln = PairwiseAligner(); aln.mode = "local"; aln.match_score = 1; aln.mismatch_score = 0
aln.open_gap_score = -1; aln.extend_gap_score = -0.5
def pid(a, b):
    c = aln.align(a, b)[0].counts(); d = c.identities + c.mismatches + c.gaps
    return 100.0 * c.identities / d if d else 0.0

def gid(fid):
    p = str(fid).split("."); return ".".join(p[:2]) if len(p) >= 2 else str(fid)

def read_fasta(p):
    out, nm, cur = {}, None, []
    for ln in Path(p).read_text().splitlines():
        if ln.startswith(">"):
            if nm is not None: out[nm] = "".join(cur)
            nm = ln[1:].strip(); cur = []
        elif ln.strip(): cur.append(ln.strip().upper())
    if nm is not None: out[nm] = "".join(cur)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=200)
    ap.add_argument("--prefilter", type=int, default=60, help="edlib top-N candidates re-scored by Biopython")
    args = ap.parse_args()
    rng = np.random.default_rng(0)

    hits = json.load(open("asv_top100_hits.json"))
    asv_seq = read_fasta(REPO / "model_inputs" / "dna-sequences.fasta")
    idx = pd.read_parquet(STORE / "index.parquet")
    mem = pd.read_parquet(STORE / "members.parquet")
    feat2row = dict(zip(mem["feature_id"], mem["row_id"]))
    print("loading reference inserts ...", flush=True)
    refs = json.load(open(PRFBA / "v4v5_refs.json"))
    feat2ins = {}
    for h, ins in refs.items():
        feat2ins[h.split("|")[1].strip() if "|" in h else h] = ins
    del refs
    # row -> insert (representative feature_id of each unique store row)
    N = len(idx)
    inserts = [None] * N
    for r, f in zip(idx["row_id"].values, idx["feature_id"].values):
        inserts[int(r)] = feat2ins.get(f)
    inserts_b = [s.encode() if s else b"" for s in inserts]      # edlib wants bytes

    asvs = [a for a in hits if a in asv_seq and asv_seq[a]]
    sample = list(rng.choice(asvs, size=min(args.sample, len(asvs)), replace=False))
    print(f"evaluating {len(sample)} ASVs vs {N:,} inserts (edlib prefilter -> Biopython top-5) ...", flush=True)

    t0 = time.time(); frac_capt, all5, captured_counts, ident_equiv5 = [], 0, [], []
    for n, a in enumerate(sample):
        s = asv_seq[a]; sb = s.encode()
        # 1. edlib edit distance to every insert (infix mode handles small length diffs)
        ed = np.fromiter((edlib.align(sb, ib, mode="HW", task="distance")["editDistance"] if ib else 10**6
                          for ib in inserts_b), dtype=np.int32, count=N)
        cand = np.argpartition(ed, args.prefilter)[:args.prefilter]      # top-N rows by smallest edit distance
        # 2. Biopython local %identity -> TRUE top-5 rows
        cand_pid = [(int(r), pid(s, inserts[int(r)])) for r in cand if inserts[int(r)]]
        ranked = sorted(cand_pid, key=lambda x: -x[1])
        true5 = {r for r, _ in ranked[:5]}
        thr5 = ranked[4][1] if len(ranked) >= 5 else (ranked[-1][1] if ranked else 0.0)  # 5th-best identity
        # 3. cosine top-100 rows
        cos_rows = [feat2row.get(h["feature_id"]) for h in hits[a]["top100"]]
        cos_rows = [r for r in cos_rows if r is not None]
        cap = len(true5 & set(cos_rows))
        captured_counts.append(cap); frac_capt.append(cap / 5.0); all5 += int(cap == 5)
        # identity-equivalent: does cosine-top-100 hold >=5 inserts as identical as the true 5th?
        n_ge = sum(1 for r in cos_rows if inserts[r] and pid(s, inserts[r]) >= thr5 - 1e-9)
        ident_equiv5.append(int(n_ge >= 5))
        if (n + 1) % 40 == 0:
            print(f"  {n+1}/{len(sample)}  all5={100*all5/(n+1):.1f}%  ({time.time()-t0:.0f}s)", flush=True)

    cc = np.array(captured_counts)
    out = {
        "n_asv": len(sample), "n_refs": int(N), "prefilter": args.prefilter,
        "pct_cosine100_contains_all5_true": round(100 * all5 / len(sample), 1),
        "pct_cosine100_contains_5_identity_equivalent": round(100 * float(np.mean(ident_equiv5)), 1),
        "pct_contains_ge4_of5": round(100 * float((cc >= 4).mean()), 1),
        "pct_contains_ge3_of5": round(100 * float((cc >= 3).mean()), 1),
        "pct_contains_true_top1": round(100 * float((cc >= 1).mean()), 1),  # at least the #1 (approx)
        "mean_fraction_of_true5_captured": round(float(np.mean(frac_capt)), 4),
        "captured_count_hist": {str(k): int((cc == k).sum()) for k in range(6)},
        "wall_seconds": round(time.time() - t0, 1),
    }
    json.dump(out, open("true_top5_in_cosine100.json", "w"), indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
