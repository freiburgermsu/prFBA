#!/usr/bin/env python
"""
query.py — reference new 16S sequence(s) against the embedding space by cosine similarity.

The store (embeddings.f16.npy) holds L2-normalized vectors, so cosine similarity is just a
dot product. A query sequence is embedded with the identical NT-v2 encoder (nt_embed.py),
L2-normalized, and matched against every stored vector on the GPU in one matmul.

Usage
-----
  # single sequence
  python query.py --seq ACGT...  --topk 10

  # a FASTA file (or a {header: seq} JSON), many queries at once
  python query.py --fasta new_seqs.fasta --topk 5 --out hits.csv

  # restrict to confident hits and report novelty
  python query.py --fasta new_seqs.fasta --min-sim 0.9

Interpretation
--------------
  cos >= ~0.99  near-identical 16S (often an exact/near-exact DB match; check md5)
  cos  ~0.95-0.99  same species / very close relative
  cos  ~0.85-0.95  same genus / family neighborhood
  cos  < ~0.80   no close reference -> query sits in a sparse/novel region of the space
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

STORE = Path(__file__).resolve().parent       # ../prFBA


def read_queries(args):
    """Return list of (name, sequence)."""
    if args.seq:
        return [("query", args.seq.strip().upper())]
    p = Path(args.fasta)
    text = p.read_text()
    if p.suffix == ".json" or text.lstrip().startswith("{"):
        return [(k, v.strip().upper()) for k, v in json.loads(text).items()]
    out, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf).upper()))
            name, buf = line[1:].strip(), []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf).upper()))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--seq", help="a single DNA sequence string")
    g.add_argument("--fasta", help="FASTA file or {header: seq} JSON of query sequences")
    ap.add_argument("--store", type=Path, default=STORE, help="embedding store dir (../prFBA)")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--min-sim", type=float, default=None, help="only report hits >= this cosine")
    ap.add_argument("--min-len", type=int, default=None, help="restrict references to seq_len >= this")
    ap.add_argument("--max-len", type=int, default=None, help="restrict references to seq_len <= this")
    ap.add_argument("--exclude-atypical", action="store_true", help="drop atypical-length references")
    ap.add_argument("--out", type=Path, default=None, help="write hits to CSV/JSON")
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--max-tokens", type=int, default=512)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    import nt_embed, hashlib
    import torch, pandas as pd

    idx = pd.read_parquet(args.store / "index.parquet")
    E = np.load(args.store / "embeddings.f16.npy", mmap_mode="r")
    # optional reference-set restriction by length / atypical flag
    keep = np.ones(len(idx), bool)
    if args.min_len: keep &= idx["seq_len"].values >= args.min_len
    if args.max_len: keep &= idx["seq_len"].values <= args.max_len
    if args.exclude_atypical: keep &= ~idx["atypical_len"].values.astype(bool)
    if not keep.all():
        E = np.ascontiguousarray(np.asarray(E)[keep])
        idx = idx[keep].reset_index(drop=True)
        print(f"[query] reference restricted to {len(idx):,} sequences "
              f"(len {args.min_len}-{args.max_len}, exclude_atypical={args.exclude_atypical})", file=sys.stderr)
    md5_set = dict(zip(idx["md5"], range(len(idx))))                    # md5 -> position in (filtered) idx/E
    queries = read_queries(args)
    print(f"[query] {len(queries)} sequence(s) vs {len(idx):,} reference vectors", file=sys.stderr)

    device = args.device if torch.cuda.is_available() else "cpu"
    Eg = torch.from_numpy(np.ascontiguousarray(E)).to(device)          # (N,1024) f16
    tok, model = nt_embed.load_model(device=device)

    rows = []
    for s in range(0, len(queries), args.batch_size):
        chunk = queries[s:s + args.batch_size]
        q = nt_embed.embed_batch([seq for _, seq in chunk], tok, model,
                                 device=device, max_tokens=args.max_tokens).to(Eg.dtype)
        sims = q @ Eg.T                                                 # (b, N) cosine
        k = min(args.topk, Eg.shape[0])
        top = torch.topk(sims.float(), k=k, dim=1)
        for bi, (name, seq) in enumerate(chunk):
            qmd5 = hashlib.md5(seq.encode()).hexdigest()
            exact = md5_set.get(qmd5)
            for rank, (sim, ri) in enumerate(zip(top.values[bi].tolist(), top.indices[bi].tolist()), 1):
                if args.min_sim is not None and sim < args.min_sim:
                    break
                m = idx.iloc[ri]
                rows.append({
                    "query": name, "rank": rank, "cosine": round(sim, 4),
                    "exact_md5_match": (exact == ri),
                    "row_id": int(m["row_id"]), "ref_organism": m["organism"],
                    "ref_genome_name": m["genome_name"], "ref_taxon_id": m["taxon_id"],
                    "ref_seq_len": int(m["seq_len"]), "ref_multiplicity": int(m["multiplicity"]),
                    "ref_md5": m["md5"],
                })

    res = pd.DataFrame(rows)
    if args.out:
        (res.to_csv(args.out, index=False) if args.out.suffix == ".csv"
         else res.to_json(args.out, orient="records", indent=2))
        print(f"[query] wrote {args.out}", file=sys.stderr)
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(res.to_string(index=False))


if __name__ == "__main__":
    main()
