#!/usr/bin/env python
"""
hierarchical_routing_test.py — does a coarse-to-fine "cluster -> representative -> descend" tree on the
embedding space actually reach each ASV's TRUE best-identity hit?

The proposed scheme: k-means the embeddings into ~K clusters, pick one representative per cluster, align
the query to the K representatives, descend GREEDILY into the best-scoring representative's cluster,
sub-cluster, repeat. We test whether greedy descent keeps the true #1 (Biopython %identity over the whole
DB) on the path: % of ASVs whose true #1 lies in the cluster of the best-scoring representative (1 level),
within the top-2/3 representatives (beam), and after a 2nd greedy level. Compared against flat cosine recall.
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np, pandas as pd, edlib, torch
from Bio.Align import PairwiseAligner
from sklearn.cluster import MiniBatchKMeans
import sys; sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import nt_embed

REPO = Path("/home/freiburger/Documents/codiffusion_bioreactor"); PRFBA = Path("/home/freiburger/Documents/prFBA")
STORE = PRFBA / "v4v5_store"; K = 10; SAMPLE = 150; PREFILTER = 60
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

def reps_of(rows, E, labels, centroids):
    """representative row per cluster = the member whose embedding is nearest its centroid."""
    rep = {}
    for c in range(centroids.shape[0]):
        members = rows[labels == c]                  # labels is positionally aligned to rows
        if len(members) == 0: continue
        d = (E[members] * torch.tensor(centroids[c], device=E.device)).sum(1)   # cosine to centroid
        rep[c] = int(members[int(d.argmax())])
    return rep

def main():
    asv_seq = read_fasta(REPO / "model_inputs" / "dna-sequences.fasta")
    idx = pd.read_parquet(STORE / "index.parquet"); N = len(idx)
    refs = json.load(open(PRFBA / "v4v5_refs.json"))
    feat2ins = {(h.split("|")[1].strip() if "|" in h else h): ins for h, ins in refs.items()}; del refs
    inserts = [feat2ins.get(f) for f in idx["feature_id"].values]
    inserts_b = [s.encode() if s else b"" for s in inserts]

    Efp = np.ascontiguousarray(np.load(STORE / "embeddings.f16.npy")).astype(np.float32)
    E = torch.from_numpy(Efp).cuda()
    print(f"k-means K={K} on {N:,} embeddings ...", flush=True)
    km = MiniBatchKMeans(n_clusters=K, random_state=0, n_init=3, batch_size=4096).fit(Efp)
    labels = km.labels_; centroids = km.cluster_centers_
    centroids /= (np.linalg.norm(centroids, axis=1, keepdims=True) + 1e-9)
    sizes = np.bincount(labels, minlength=K)
    all_rows = np.arange(N)
    rep = reps_of(all_rows, E, labels, centroids)

    rng = np.random.default_rng(0)
    hk = set(json.load(open("asv_top100_hits.json")))
    asvs = [a for a in asv_seq if a in hk and asv_seq[a]]
    sample = list(rng.choice(asvs, size=min(SAMPLE, len(asvs)), replace=False))
    tok, model = nt_embed.load_model()

    print(f"routing {len(sample)} ASVs (true #1 via edlib->Biopython over {N:,} inserts) ...", flush=True)
    t0 = time.time()
    g1 = b2 = b3 = g2 = flat_clustersize = 0
    for n, a in enumerate(sample):
        s = asv_seq[a]; sb = s.encode()
        # TRUE #1 over whole DB
        ed = np.fromiter((edlib.align(sb, ib, mode="HW", task="distance")["editDistance"] if ib else 10**6
                          for ib in inserts_b), dtype=np.int32, count=N)
        cand = np.argpartition(ed, PREFILTER)[:PREFILTER]
        true1 = max(((int(r), pid(s, inserts[int(r)])) for r in cand if inserts[int(r)]), key=lambda x: x[1])[0]
        true_cluster = labels[true1]
        # LEVEL 1: align query to the K representatives, rank clusters by rep %identity
        rep_score = sorted(((c, pid(s, inserts[rep[c]])) for c in rep), key=lambda x: -x[1])
        ranked_clusters = [c for c, _ in rep_score]
        g1 += int(true_cluster == ranked_clusters[0])
        b2 += int(true_cluster in ranked_clusters[:2])
        b3 += int(true_cluster in ranked_clusters[:3])
        flat_clustersize += sizes[ranked_clusters[0]]
        # LEVEL 2 (greedy): sub-cluster the best cluster, route again
        if true_cluster == ranked_clusters[0]:
            mrows = all_rows[labels == ranked_clusters[0]]
            if len(mrows) > K:
                sub = MiniBatchKMeans(n_clusters=K, random_state=0, n_init=2, batch_size=2048).fit(Efp[mrows])
                subc = sub.cluster_centers_ / (np.linalg.norm(sub.cluster_centers_, axis=1, keepdims=True) + 1e-9)
                subrep = reps_of(mrows, E, sub.labels_, subc)
                best_sub = max(subrep, key=lambda c: pid(s, inserts[subrep[c]]))
                g2 += int(sub.labels_[np.where(mrows == true1)[0][0]] == best_sub)
            else:
                g2 += 1
        if (n + 1) % 50 == 0:
            print(f"  {n+1}/{len(sample)} ({time.time()-t0:.0f}s) g1={100*g1/(n+1):.0f}% b3={100*b3/(n+1):.0f}%", flush=True)

    m = len(sample)
    out = {"n_asv": m, "n_refs": int(N), "K": K, "mean_cluster_size": int(sizes.mean()),
           "alignments_per_query_greedy": K, "alignments_flat_rerank_top500": 500,
           "greedy_1level_reaches_true1_cluster_pct": round(100 * g1 / m, 1),
           "beam2_1level_pct": round(100 * b2 / m, 1), "beam3_1level_pct": round(100 * b3 / m, 1),
           "greedy_2level_keeps_true1_pct": round(100 * g2 / m, 1),
           "note": "compare vs flat cosine recall@500 of true #1 = 88.5%; greedy descent aligns only K reps/level",
           "wall_seconds": round(time.time() - t0, 1)}
    json.dump(out, open("hierarchical_routing_test.json", "w"), indent=2)
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
