#!/usr/bin/env python
"""
kmer_augmented_search.py — does adding a classical k-mer retriever (+ RRF fusion) beat NT-v2 cosine
at finding the true closest 16S hit?  NO machine learning beyond the existing embedding.

Retrievers compared, per ASV, against the TRUE ranking by Biopython %identity over the whole
97,624-insert DB (edlib prefilter -> Biopython):
  - cosine        : NT-v2 mean-pooled embedding cosine (current method)
  - kmer8 / kmer12 : TF-IDF k-mer cosine (sklearn char n-grams; classical, alignment-free, no training)
  - rrf8 / rrf12  : Reciprocal Rank Fusion of cosine + k-mer ranks, score = sum 1/(60+rank)

Metrics at k in {100,200,500}: recall of the true #1 (= the ceiling an exact reranker could reach),
fraction with all true-top-5 present, and top-1 correctness (method's #1 == true #1).
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np, pandas as pd, edlib, torch
from Bio.Align import PairwiseAligner
from sklearn.feature_extraction.text import TfidfVectorizer
import sys; sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import nt_embed

REPO = Path("/home/freiburger/Documents/codiffusion_bioreactor"); PRFBA = Path("/home/freiburger/Documents/prFBA")
STORE = PRFBA / "v4v5_store"; SAMPLE = 200; PREFILTER = 60; CUTOFFS = [100, 200, 500]; RRF_C = 60
aln = PairwiseAligner(); aln.mode = "local"; aln.match_score = 1; aln.mismatch_score = 0
aln.open_gap_score = -1; aln.extend_gap_score = -0.5
def pid(a, b):
    if not b: return -1.0
    c = aln.align(a, b)[0].counts(); d = c.identities + c.mismatches + c.gaps
    return 100.0 * c.identities / d if d else -1.0
def read_fasta(p):
    o={};nm=None;cur=[]
    for ln in Path(p).read_text().splitlines():
        if ln.startswith('>'):
            if nm: o[nm]=''.join(cur)
            nm=ln[1:].strip();cur=[]
        elif ln.strip():cur.append(ln.strip().upper())
    if nm:o[nm]=''.join(cur)
    return o
def ranks(sim_row):                          # ref -> rank position (0 = best), descending similarity
    order = np.argsort(-sim_row, kind="stable"); r = np.empty(len(order), np.int32); r[order] = np.arange(len(order)); return r

def main():
    asv_seq = read_fasta(REPO / "model_inputs" / "dna-sequences.fasta")
    idx = pd.read_parquet(STORE / "index.parquet"); N = len(idx)
    refs = json.load(open(PRFBA / "v4v5_refs.json"))
    feat2ins = {(h.split("|")[1].strip() if "|" in h else h): ins for h, ins in refs.items()}; del refs
    inserts = [feat2ins.get(f) for f in idx["feature_id"].values]
    inserts_b = [s.encode() if s else b"" for s in inserts]
    valid = np.array([bool(s) for s in inserts])

    rng = np.random.default_rng(0)
    hk = set(json.load(open("asv_top100_hits.json")))
    asvs = [a for a in asv_seq if a in hk and asv_seq[a]]
    sample = list(rng.choice(asvs, size=min(SAMPLE, len(asvs)), replace=False))
    qseq = [asv_seq[a] for a in sample]

    # ---- cosine similarity (GPU) ----
    print("embedding ASVs + cosine sim ...", flush=True)
    tok, model = nt_embed.load_model()
    E = torch.from_numpy(np.ascontiguousarray(np.load(STORE / "embeddings.f16.npy"))).cuda()
    Q = nt_embed.embed_batch(qseq, tok, model).to(E.dtype)
    cos = (Q @ E.T).float().cpu().numpy()                       # (n, N)
    del E; torch.cuda.empty_cache()

    # ---- TF-IDF k-mer cosine (classical, no ML) ----
    ref_seqs = [s if s else "N" for s in inserts]
    kmer_sim = {}
    for k in (8, 12):
        print(f"TF-IDF {k}-mer cosine ...", flush=True)
        vec = TfidfVectorizer(analyzer="char", ngram_range=(k, k), lowercase=False, norm="l2", dtype=np.float32)
        R = vec.fit_transform(ref_seqs)                         # (N, V) L2-normalized
        Qk = vec.transform(qseq)                                # (n, V)
        kmer_sim[k] = np.asarray((Qk @ R.T).todense())          # (n, N) cosine
        del R, Qk

    # ---- TRUE #1 + top-5 over whole DB (edlib -> Biopython) ----
    print(f"true ranking (edlib->Biopython) over {N:,} inserts for {len(sample)} ASVs ...", flush=True)
    t0 = time.time()
    methods = ["cosine", "kmer8", "kmer12", "rrf8", "rrf12"]
    agg = {m: {k: {"top1cap": 0, "all5": 0} for k in CUTOFFS} for m in methods}
    top1corr = {m: 0 for m in methods}
    for i, a in enumerate(sample):
        s = qseq[i]; sb = s.encode()
        ed = np.fromiter((edlib.align(sb, ib, mode="HW", task="distance")["editDistance"] if ib else 10**6
                          for ib in inserts_b), dtype=np.int32, count=N)
        cand = np.argpartition(ed, PREFILTER)[:PREFILTER]
        ranked = sorted(((int(r), pid(s, inserts[int(r)])) for r in cand if inserts[int(r)]), key=lambda x: -x[1])
        true1 = ranked[0][0]; true5 = {r for r, _ in ranked[:5]}
        rc, rk8, rk12 = ranks(cos[i]), ranks(kmer_sim[8][i]), ranks(kmer_sim[12][i])
        rrf8 = 1.0/(RRF_C+rc) + 1.0/(RRF_C+rk8); rrf12 = 1.0/(RRF_C+rc) + 1.0/(RRF_C+rk12)
        orders = {"cosine": np.argsort(-cos[i]), "kmer8": np.argsort(-kmer_sim[8][i]),
                  "kmer12": np.argsort(-kmer_sim[12][i]), "rrf8": np.argsort(-rrf8), "rrf12": np.argsort(-rrf12)}
        for m, order in orders.items():
            top1corr[m] += int(order[0] == true1)
            for k in CUTOFFS:
                topk = set(order[:k].tolist())
                agg[m][k]["top1cap"] += int(true1 in topk)
                agg[m][k]["all5"] += int(len(true5 & topk) == 5)
        if (i+1) % 50 == 0: print(f"  {i+1}/{len(sample)} ({time.time()-t0:.0f}s)", flush=True)

    n = len(sample)
    out = {"n_asv": n, "n_refs": int(N), "rrf_c": RRF_C,
           "top1_correct_pct": {m: round(100*top1corr[m]/n, 1) for m in methods},
           "recall_true1": {m: {str(k): round(100*agg[m][k]["top1cap"]/n, 1) for k in CUTOFFS} for m in methods},
           "all_true_top5_present": {m: {str(k): round(100*agg[m][k]["all5"]/n, 1) for k in CUTOFFS} for m in methods},
           "wall_seconds": round(time.time()-t0, 1)}
    json.dump(out, open("kmer_augmented_search.json", "w"), indent=2)
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
