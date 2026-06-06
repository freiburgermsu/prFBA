#!/usr/bin/env python
"""
true_top5_recall_curve.py — at what cosine-top-k do we capture the TRUE (Biopython) top-5?

Extends true_top5_in_cosine100.py to a cutoff curve: for a sample of ASVs, find the true top-5
references by Biopython %identity over the whole 97,624-insert DB (edlib prefilter -> Biopython
re-score), compute the cosine ranking up to k=500 directly on the GPU, and report, for k in
{100,200,500}, the % of ASVs whose cosine top-k contains all 5 true hits (and 5 identity-equivalent).
"""
from __future__ import annotations
import json, time, sys
from pathlib import Path
import numpy as np, pandas as pd, edlib, torch
from Bio.Align import PairwiseAligner
sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import nt_embed

REPO = Path("/home/freiburger/Documents/codiffusion_bioreactor")
PRFBA = Path("/home/freiburger/Documents/prFBA")
STORE = PRFBA / "v4v5_store"
CUTOFFS = [100, 200, 500]
PREFILTER = 60

aln = PairwiseAligner(); aln.mode = "local"; aln.match_score = 1; aln.mismatch_score = 0
aln.open_gap_score = -1; aln.extend_gap_score = -0.5
def pid(a, b):
    if not b: return 0.0
    c = aln.align(a, b)[0].counts(); d = c.identities + c.mismatches + c.gaps
    return 100.0 * c.identities / d if d else 0.0

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
    rng = np.random.default_rng(0)               # same seed/sample as true_top5_in_cosine100.py
    asv_seq = read_fasta(REPO / "model_inputs" / "dna-sequences.fasta")
    idx = pd.read_parquet(STORE / "index.parquet")
    print("loading reference inserts ...", flush=True)
    refs = json.load(open(PRFBA / "v4v5_refs.json"))
    feat2ins = {(h.split("|")[1].strip() if "|" in h else h): ins for h, ins in refs.items()}
    del refs
    N = len(idx)
    inserts = [None] * N
    for r, f in zip(idx["row_id"].values, idx["feature_id"].values):
        inserts[int(r)] = feat2ins.get(f)
    inserts_b = [s.encode() if s else b"" for s in inserts]

    hit_keys = set(json.load(open("asv_top100_hits.json")))          # load once, not per ASV
    asvs = [a for a in asv_seq if a in hit_keys and asv_seq[a]]
    sample = list(rng.choice(asvs, size=min(200, len(asvs)), replace=False))

    # cosine top-500 rows per sampled ASV (GPU)
    print("embedding ASVs + cosine top-500 ...", flush=True)
    tok, model = nt_embed.load_model()
    E = torch.from_numpy(np.ascontiguousarray(np.load(STORE / "embeddings.f16.npy"))).cuda()
    cos_rows = {}
    for s in range(0, len(sample), 128):
        chunk = sample[s:s + 128]
        Q = nt_embed.embed_batch([asv_seq[a] for a in chunk], tok, model).to(E.dtype)
        top = (Q @ E.T).float().topk(max(CUTOFFS), dim=1).indices.cpu().numpy()
        for a, rowids in zip(chunk, top):
            cos_rows[a] = rowids                 # ordered store row_ids, best first

    print(f"true-top-5 (edlib->Biopython) over {N:,} inserts for {len(sample)} ASVs ...", flush=True)
    t0 = time.time()
    res = {str(k): {"all5": 0, "frac": []} for k in CUTOFFS}
    for n, a in enumerate(sample):
        s = asv_seq[a]; sb = s.encode()
        ed = np.fromiter((edlib.align(sb, ib, mode="HW", task="distance")["editDistance"] if ib else 10**6
                          for ib in inserts_b), dtype=np.int32, count=N)
        cand = np.argpartition(ed, PREFILTER)[:PREFILTER]
        ranked = sorted(((int(r), pid(s, inserts[int(r)])) for r in cand if inserts[int(r)]), key=lambda x: -x[1])
        true5 = {r for r, _ in ranked[:5]}
        rows500 = cos_rows[a]
        for k in CUTOFFS:                                  # strict capture = row membership (no alignment)
            cap = len(true5 & set(int(r) for r in rows500[:k]))
            res[str(k)]["all5"] += int(cap == 5)
            res[str(k)]["frac"].append(cap / 5.0)
        if (n + 1) % 50 == 0:
            print(f"  {n+1}/{len(sample)} ({time.time()-t0:.0f}s)  k100 all5={100*res['100']['all5']/(n+1):.0f}%  "
                  f"k500 all5={100*res['500']['all5']/(n+1):.0f}%", flush=True)

    out = {"n_asv": len(sample), "n_refs": int(N), "cutoffs": {}}
    for k in CUTOFFS:
        d = res[str(k)]
        out["cutoffs"][str(k)] = {
            "pct_contains_all5_true": round(100 * d["all5"] / len(sample), 1),
            "mean_fraction_of_true5_captured": round(float(np.mean(d["frac"])), 4),
        }
    json.dump(out, open("true_top5_recall_curve.json", "w"), indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
