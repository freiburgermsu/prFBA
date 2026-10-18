#!/usr/bin/env python
"""Select a 10,000-genome panel by k-means over the BV-BRC 16S embedding space.

An alternative to the taxonomy-stratified selection of select_10k.py, which
deliberately sqrt-flattens availability and over-weights rare phyla and archaea.
Here the panel is an embedding-uniform sample instead: k-means with k = 10,000 over
the L2-normalised nucleotide-transformer embeddings, taking the organism whose 16S
sequence lies nearest each cluster centroid.

Design points that matter:
  * Clustering is over SEQUENCES, not genomes -- that is what "the centroid organism
    of each cluster" means, and it respects the space as embedded. A genome carrying
    several distinct operons can therefore be nearest to more than one centroid, so
    collisions are resolved by walking out from the centroid to the next-nearest
    sequence whose genome is still unclaimed.
  * Only full-length (>= 1400 bp, non-atypical) sequences are clusterable, because the
    panel must survive in-silico PCR for all eight regions. That is the same
    FULL_LEN gate select_10k.py applies, minus its 1200 bp backfill.
  * Vectors are unit-norm, so squared Euclidean distance and cosine are monotonically
    equivalent; k-means on them is spherical k-means up to a constant.
  * Deterministic: seeded k-means++ init and a fixed iteration count, SEED = 1729 to
    match the original panel.

Emits selection_kmeans.json and src16s_kmeans.fasta in exactly the schema
run_validation.sh consumes, so the rest of the 16S pipeline runs unchanged.
"""
import argparse, json, os, sys, time
import numpy as np
import pandas as pd
import torch

ROOT = os.path.expanduser("~/Documents/prFBA")
sys.path.insert(0, os.path.join(ROOT, "region_validation", "scripts"))
import _common

SEED = 1729
FULL_LEN = 1400
ALLOWED_DOMAINS = ("Archaea", "Bacteria")
MD5_SEQ = "/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_seq.json"
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
LINKEY = {"domain": "Kingdom", "phylum": "Phylum", "class": "Class", "order": "Order",
          "family": "Family", "genus": "Genus", "species": "Species"}


def log(m):
    print("[kmeans] %s" % m, flush=True)


_LIN = {}


def lineage(taxon_id):
    """select_10k.py's lineage rule, verbatim: a genome is admissible only with a
    domain in {Archaea, Bacteria} AND a named phylum; missing intermediate ranks get
    a stable `unclassified_<parentrank>:<parentname>` placeholder so genus-less clades
    form their own bin; species falls back to sp_taxon_<tid>. Returns None to drop.

    NOT _common.lineage_for(): that maps Kingdom to the modern NCBI superkingdom
    ('Pseudomonadati', 'Bacillati'), which is not the Bacteria/Archaea label that
    insilico_pcr.py reads to choose bacterial vs archaeal primers.
    """
    if taxon_id in _LIN:
        return _LIN[taxon_id]
    import taxopy
    try:
        rd = taxopy.Taxon(int(taxon_id), _common._get_taxdb()).rank_name_dictionary
    except Exception:
        _LIN[taxon_id] = None
        return None
    domain = rd.get("domain") or rd.get("superkingdom")
    phylum = rd.get("phylum")
    if domain not in ALLOWED_DOMAINS or not phylum:
        _LIN[taxon_id] = None
        return None
    lin = {"domain": domain, "phylum": phylum}
    pr, pn = "phylum", phylum
    for rank in ("class", "order", "family", "genus"):
        nm = rd.get(rank)
        if nm:
            lin[rank] = nm; pr, pn = rank, nm
        else:
            lin[rank] = "unclassified_%s:%s" % (pr, pn)
    lin["species"] = rd.get("species") or "sp_taxon_%s" % taxon_id
    _LIN[taxon_id] = lin
    return lin


def kmeans_gpu(X, k, iters, seed, batch):
    """Minibatch k-means with k-means++ seeding, on GPU, float32 accumulators."""
    g = torch.Generator(device=X.device).manual_seed(seed)
    n = X.shape[0]
    # --- k-means++ init (distance-weighted, the standard seeding) ---
    idx = [int(torch.randint(n, (1,), generator=g, device=X.device))]
    d2 = (2.0 - 2.0 * (X @ X[idx[0]]).float()).clamp_(min=0)
    for i in range(1, k):
        s = d2.sum()
        j = int(torch.multinomial(d2 / s, 1, generator=g)) if s > 0 else \
            int(torch.randint(n, (1,), generator=g, device=X.device))
        idx.append(j)
        torch.minimum(d2, (2.0 - 2.0 * (X @ X[j]).float()).clamp_(min=0), out=d2)
        if i % 2000 == 0:
            log("  kmeans++ seeded %d/%d" % (i, k))
    C = X[torch.tensor(idx, device=X.device)].clone()
    # --- minibatch Lloyd iterations ---
    for it in range(iters):
        perm = torch.randperm(n, generator=g, device=X.device)
        sums = torch.zeros_like(C, dtype=torch.float32)
        cnts = torch.zeros(k, device=X.device, dtype=torch.float32)
        for s in range(0, n, batch):
            B = X[perm[s:s + batch]]
            a = (B @ C.T).argmax(1)                  # unit-norm => max dot == min distance
            sums.index_add_(0, a, B.float())
            cnts.index_add_(0, a, torch.ones(B.shape[0], device=X.device))
        keep = cnts > 0
        C[keep] = (sums[keep] / cnts[keep, None]).to(C.dtype)
        C = torch.nn.functional.normalize(C, dim=1)   # back onto the sphere
        log("  iter %d/%d  empty clusters=%d" % (it + 1, iters, int((~keep).sum())))
    return C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=10000)
    ap.add_argument("--iters", type=int, default=15)
    ap.add_argument("--batch", type=int, default=16384)
    ap.add_argument("--outdir", default=os.path.join(ROOT, "kmeans_panel"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    df = pd.read_parquet(os.path.join(ROOT, "index.parquet"))
    mask = (df.seq_len >= FULL_LEN) & (~df.atypical_len)
    sub = df[mask].reset_index(drop=True)
    log("full-length: %d sequences over %d genomes" % (len(sub), sub.genome_id.nunique()))
    # admissibility is a property of the organism, so apply it BEFORE clustering --
    # otherwise centroids land on genomes the panel cannot use and the k is not honoured
    tids = sorted({int(x) for x in sub.taxon_id.dropna().unique()})
    log("resolving %d distinct taxa ..." % len(tids))
    ok = {t for t in tids if lineage(t) is not None}
    log("admissible taxa (domain in {Archaea,Bacteria} and named phylum): %d/%d"
        % (len(ok), len(tids)))
    sub = sub[sub.taxon_id.notna() & sub.taxon_id.astype("Int64").isin(ok)].reset_index(drop=True)
    rows = sub.row_id.to_numpy()
    log("clusterable: %d sequences over %d admissible genomes"
        % (len(sub), sub.genome_id.nunique()))

    emb = np.load(os.path.join(ROOT, "embeddings.f16.npy"), mmap_mode="r")
    log("embeddings %s %s -> taking %d rows" % (emb.shape, emb.dtype, len(rows)))
    X = torch.from_numpy(np.ascontiguousarray(emb[rows])).cuda().half()
    X = torch.nn.functional.normalize(X.float(), dim=1).half()
    log("on GPU: %s, %.1f GB" % (tuple(X.shape), X.element_size() * X.nelement() / 1e9))

    t0 = time.time()
    C = kmeans_gpu(X, a.k, a.iters, SEED, a.batch)
    log("k-means done in %.0fs" % (time.time() - t0))

    # distance from every centroid to every sequence, in centroid chunks
    log("assigning centroid-nearest sequences ...")
    order = np.empty((a.k, 64), dtype=np.int64)       # 64 nearest seqs per centroid
    for s in range(0, a.k, 512):
        sim = (C[s:s + 512] @ X.T).float()       # max dot == nearest on the unit sphere
        order[s:s + 512] = sim.topk(64, largest=True).indices.cpu().numpy()

    # one genome per cluster; on collision walk out to the next-nearest sequence
    gid = sub.genome_id.to_numpy()
    claimed, picks, unresolved = {}, [], 0
    for c in range(a.k):
        for rank in range(order.shape[1]):
            g = gid[order[c, rank]]
            if g not in claimed:
                claimed[g] = c
                picks.append((c, int(order[c, rank]), rank))
                break
        else:
            unresolved += 1
    log("resolved %d/%d clusters to distinct genomes (%d unresolved within 64 NN)"
        % (len(picks), a.k, unresolved))
    log("  centroid-nearest was already claimed for %d clusters"
        % sum(1 for _, _, r in picks if r > 0))

    # index.parquet keys on md5(sequence.upper()); 16S_md5_seq.json keys on md5 of the
    # raw (lower-case) sequence, so the two share no keys. Bridge by re-hashing.
    import hashlib
    raw = json.load(open(MD5_SEQ))
    seqs = {hashlib.md5(v.upper().encode()).hexdigest(): v for v in raw.values()}
    log("sequence map: %d raw md5 -> %d upper-md5 keys" % (len(raw), len(seqs)))
    del raw
    genomes, fasta = [], []
    for c, ri, rank in picks:
        r = sub.iloc[ri]
        tid = int(r.taxon_id)
        lin = lineage(tid)
        genomes.append(dict(
            genome_id=str(r.genome_id), taxon_id=tid,
            best_16s_md5=str(r.md5), best_16s_len=int(r.seq_len),
            n_full_length_16s=int((sub.genome_id == r.genome_id).sum()),
            **{k: lin[k] for k in RANKS},
            seed=SEED, full_len=FULL_LEN, cluster=c, centroid_rank=rank))
        fasta.append(">%s\n%s\n" % (r.genome_id, seqs[str(r.md5)]))

    sel = os.path.join(a.outdir, "selection_kmeans.json")
    json.dump(dict(meta=dict(method="kmeans", k=a.k, iters=a.iters, seed=SEED,
                             full_len=FULL_LEN, selected_n=len(genomes),
                             n_clusterable_sequences=int(len(sub)),
                             n_clusterable_genomes=int(sub.genome_id.nunique()),
                             n_unresolved_clusters=unresolved,
                             embedding=os.path.join(ROOT, "embeddings.f16.npy")),
                   genomes=genomes), open(sel, "w"))
    fa = os.path.join(a.outdir, "src16s_kmeans.fasta")
    open(fa, "w").write("".join(fasta))
    log("wrote %s (%d genomes) and %s" % (sel, len(genomes), fa))


if __name__ == "__main__":
    main()
