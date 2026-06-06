#!/usr/bin/env python
"""
part2_cosine_vs_identity.py — does NT-v2 cosine actually find the highest-IDENTITY database hits?

Cosine on learned embeddings is a *proxy* for sequence similarity. This checks it empirically against
true pairwise % identity (Bio.Align) for the V4-V5 ASVs:
  A. within-candidate: for each sampled ASV, align it to its 20 cosine-top hits -> does cosine rank
     track identity rank (Spearman), and is the cosine #1 also the identity #1 among the candidates?
  B. miss check: align each sampled ASV to a large random reference pool -> is the globally
     highest-identity reference captured in the cosine top-20, and at what cosine rank?
  C. speed: brute-force cosine vs per-query alignment (the cost of the identity-based alternative).
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np, pandas as pd, torch
from Bio.Align import PairwiseAligner

REPO = Path("/home/freiburger/Documents/codiffusion_bioreactor")
PRFBA = Path("/home/freiburger/Documents/prFBA")
HITS = Path("asv_top20_hits.json")

aln = PairwiseAligner()
aln.mode = "local"; aln.match_score = 1; aln.mismatch_score = 0
aln.open_gap_score = -1; aln.extend_gap_score = -0.5

def pid(a, b):
    """% identity = identical aligned columns / alignment length of the best local alignment."""
    al = aln.align(a, b)[0]
    c = al.counts()
    denom = c.identities + c.mismatches + c.gaps
    return 100.0 * c.identities / denom if denom else 0.0


def read_fasta(p):
    names, seqs, cur, nm = [], [], [], None
    for ln in Path(p).read_text().splitlines():
        if ln.startswith(">"):
            if nm is not None: seqs.append("".join(cur))
            nm = ln[1:].strip(); names.append(nm); cur = []
        elif ln.strip(): cur.append(ln.strip().upper())
    if nm is not None: seqs.append("".join(cur))
    return dict(zip(names, seqs))


def main():
    rng = np.random.default_rng(0)
    asv_seq = read_fasta(REPO / "model_inputs" / "dna-sequences.fasta")
    hits = json.load(open(HITS))
    print("loading reference V4-V5 inserts (feature_id -> insert) ...", flush=True)
    refs = json.load(open(PRFBA / "v4v5_refs.json"))
    feat2ins = {}
    for header, ins in refs.items():
        fid = header.split("|")[1].strip() if "|" in header else header
        feat2ins[fid] = ins
    del refs

    asvs = [a for a in hits if a in asv_seq]
    sample = list(rng.choice(asvs, size=min(150, len(asvs)), replace=False))

    # ---- A. within-candidate: cosine rank vs identity rank ----
    from scipy.stats import spearmanr
    sp, top1_agree, cos1_id, idmax = [], 0, [], []
    for a in sample:
        top = hits[a]["top20"]
        cos = np.array([h["cosine"] for h in top])
        ids = np.array([pid(asv_seq[a], feat2ins.get(h["feature_id"], "")) if feat2ins.get(h["feature_id"]) else np.nan for h in top])
        ok = ~np.isnan(ids)
        if ok.sum() < 5: continue
        r = spearmanr(cos[ok], ids[ok]).correlation
        if not np.isnan(r): sp.append(r)
        top1_agree += int(np.nanargmax(ids) == 0)            # cosine #1 == identity #1 among candidates
        cos1_id.append(ids[0]); idmax.append(np.nanmax(ids))
    print(f"\nA. within top-20 candidates (n={len(sp)}):")
    print(f"   mean Spearman(cosine, %identity) = {np.mean(sp):.3f}")
    print(f"   cosine-#1 is also identity-#1 among candidates: {100*top1_agree/len(cos1_id):.1f}%")
    print(f"   %identity of cosine-#1 hit: mean={np.mean(cos1_id):.1f}  (best-in-candidates mean={np.mean(idmax):.1f})")

    # ---- B. miss check vs a random reference pool ----
    pool_feats = rng.choice(list(feat2ins), size=4000, replace=False)
    pool = [(f, feat2ins[f]) for f in pool_feats]
    miss = 0; ranks = []; sub = sample[:60]
    for a in sub:
        s = asv_seq[a]
        pid_pool = np.array([pid(s, ins) for _, ins in pool])
        gmax = pid_pool.max()
        cand_ids = np.array([pid(s, feat2ins.get(h["feature_id"], "")) if feat2ins.get(h["feature_id"]) else 0 for h in hits[a]["top20"]])
        if cand_ids.max() + 1e-6 < gmax:                     # pool had a strictly higher-identity ref
            miss += 1
    print(f"\nB. miss check (n={len(sub)} ASVs vs 4,000 random refs each):")
    print(f"   ASVs where a random ref beat the best cosine-top-20 hit on %identity: {miss}/{len(sub)} ({100*miss/len(sub):.1f}%)")

    # ---- C. speed: cosine vs alignment ----
    E = torch.from_numpy(np.ascontiguousarray(np.load(PRFBA / "v4v5_store" / "embeddings.f16.npy"))).cuda()
    q = E[:3215]                                             # stand-in query batch of same size
    torch.cuda.synchronize(); t = time.time()
    for s in range(0, q.shape[0], 512):
        (q[s:s+512] @ E.T).float().topk(20, dim=1)
    torch.cuda.synchronize(); t_cos = time.time() - t
    t = time.time(); [pid(asv_seq[sample[0]], ins) for _, ins in pool[:200]]; t_aln = (time.time() - t) / 200
    print(f"\nC. speed:")
    print(f"   brute-force cosine, 3,215 queries x 97,624 refs: {t_cos:.2f} s ({t_cos/3215*1e3:.2f} ms/query)")
    print(f"   pairwise alignment: {t_aln*1e3:.1f} ms/pair -> {t_aln*97624:.0f} s/query vs full DB (x{t_aln*97624/(t_cos/3215):.0f} slower)")


if __name__ == "__main__":
    main()
