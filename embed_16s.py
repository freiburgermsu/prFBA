#!/usr/bin/env python
"""
embed_16s.py — build the 16S rRNA embedding space.

Reads BV-BRC 16S sequences (model_inputs/BV_BRC_16S.json : {fasta_header -> sequence})
and their metadata (model_inputs/16S_metadata.json : {feature_id -> record}), embeds
every *unique* sequence with Nucleotide-Transformer-v2-500m, and writes a self-contained
embedding store to --outdir (default ../prFBA).

Design notes
------------
* Exact dedup. 16S genes are massively redundant across BV-BRC genomes; we embed each
  distinct nucleotide string once (md5 key) and carry a multiplicity + member list. This
  cuts GPU work several-fold and makes downstream density estimates meaningful (otherwise
  "density" would just measure how many genomes share an identical gene copy).
* Length-sorted batching. Unique sequences are processed shortest->longest so each padded
  batch wastes minimal compute; results are scattered back to a stable row_id.
* Resumable. Embeddings stream into a float16 memmap; progress.json records how many
  length-sorted items are finished. Re-running continues where it stopped (input + sort
  order are deterministic).

Outputs (in --outdir)
---------------------
  embeddings.f16.npy   (N_unique, 1024) L2-normalized float16  — the embedding space
  index.parquet        one row per unique sequence: row_id, md5, seq_len, multiplicity,
                       n_genomes, representative feature_id / organism / genome_id /
                       taxon_id / genome_name, atypical_len flag
  members.parquet      every input sequence: feature_id, row_id, seq_len (reverse lookup)
  sort_order.npy       length-sorted row_id order used for batching (for resume)
  manifest.json        full provenance: model, params, counts, dtype, timings, config
"""
from __future__ import annotations
import argparse, json, hashlib, re, sys, time
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parent.parent / "codiffusion_bioreactor"
DEF_FASTA = REPO / "model_inputs" / "BV_BRC_16S.json"
DEF_META = REPO / "model_inputs" / "16S_metadata.json"
DEF_OUT = Path(__file__).resolve().parent           # ../prFBA

# organism string lives in the trailing  [ ... ]  of the fasta header
_ORG_RE = re.compile(r"\[(?P<body>.*)\]\s*$")


def parse_header(header: str):
    """fig|<feature_id>|...| desc [<organism> | <genome_id>]  ->  (feature_id, organism, genome_id)."""
    parts = header.split("|")
    feature_id = parts[1].strip() if len(parts) > 1 else header.strip()
    organism, genome_id = None, None
    m = _ORG_RE.search(header)
    if m:
        body = m.group("body").strip()
        if " | " in body:
            organism, genome_id = (x.strip() for x in body.rsplit(" | ", 1))
        else:
            organism = body
    if genome_id is None and "." in feature_id:
        # feature_id like 1505.7.rna.1 -> genome_id 1505.7
        genome_id = ".".join(feature_id.split(".")[:2])
    return feature_id, organism, genome_id


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=Path, default=DEF_FASTA)
    ap.add_argument("--metadata", type=Path, default=DEF_META)
    ap.add_argument("--outdir", type=Path, default=DEF_OUT)
    ap.add_argument("--model", default=None, help="override model id")
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--max-tokens", type=int, default=512, help="token cap (~3 kb); clips only atypical mis-exports")
    ap.add_argument("--min-len", type=int, default=20, help="sequences shorter than this are recorded but not embedded")
    ap.add_argument("--atypical-len", type=int, default=2000, help="flag unique seqs longer than this as non-16S-like")
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    import nt_embed
    import torch, pandas as pd
    model_id = args.model or nt_embed.MODEL_ID
    args.outdir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    # ---- 1. load inputs -------------------------------------------------------
    log(f"loading sequences  {args.fasta}")
    fasta = json.load(open(args.fasta))
    log(f"loading metadata   {args.metadata}")
    meta = json.load(open(args.metadata))
    log(f"  {len(fasta):,} input sequences | {len(meta):,} metadata records")

    # ---- 2. dedup + join ------------------------------------------------------
    log("deduplicating by exact sequence (md5) and joining metadata ...")
    md5_to_row: dict[str, int] = {}
    rows: list[dict] = []          # per unique sequence (row_id == index)
    uniq_seqs: list[str] = []      # parallel to rows
    members = {"feature_id": [], "row_id": [], "seq_len": []}
    n_skipped = 0
    for header, raw in fasta.items():
        seq = raw.strip().upper()
        L = len(seq)
        feature_id, organism, genome_id = parse_header(header)
        if L < args.min_len:
            n_skipped += 1
            continue
        md5 = hashlib.md5(seq.encode()).hexdigest()
        row = md5_to_row.get(md5)
        if row is None:
            row = len(rows)
            md5_to_row[md5] = row
            uniq_seqs.append(seq)
            rec = meta.get(feature_id, {})
            rows.append({
                "row_id": row, "md5": md5, "seq_len": L, "multiplicity": 0,
                "genomes": set(),
                "feature_id": feature_id, "organism": organism, "genome_id": genome_id,
                "taxon_id": rec.get("taxon_id"), "genome_name": rec.get("genome_name"),
                "atypical_len": L > args.atypical_len or L < 200,
            })
        r = rows[row]
        r["multiplicity"] += 1
        if genome_id:
            r["genomes"].add(genome_id)
        members["feature_id"].append(feature_id)
        members["row_id"].append(row)
        members["seq_len"].append(L)

    N = len(rows)
    log(f"  unique sequences: {N:,}  (from {len(fasta):,}; {n_skipped:,} below min-len, "
        f"dedup factor {len(fasta)/max(N,1):.2f}x)")

    # ---- 3. length-sorted batch order (stable) --------------------------------
    order = np.array(sorted(range(N), key=lambda i: (rows[i]["seq_len"], rows[i]["md5"])), dtype=np.int64)
    np.save(args.outdir / "sort_order.npy", order)

    # ---- 4. embedding memmap + resume -----------------------------------------
    emb_path = args.outdir / "embeddings.f16.npy"
    prog_path = args.outdir / "progress.json"
    emb = np.lib.format.open_memmap(emb_path, mode="r+" if emb_path.exists() else "w+",
                                    dtype=np.float16, shape=(N, nt_embed.EMBED_DIM))
    done = 0
    if prog_path.exists():
        p = json.load(open(prog_path))
        if p.get("N") == N and p.get("model") == model_id:
            done = int(p.get("done", 0))
            log(f"resuming: {done:,}/{N:,} already embedded")

    # ---- 5. load model + embed ------------------------------------------------
    log(f"loading model {model_id}")
    tok, model = nt_embed.load_model(model_id, device=args.device)
    bs, mt = args.batch_size, args.max_tokens
    last = time.time()
    for start in range(done, N, bs):
        idx = order[start:start + bs]
        batch = [uniq_seqs[i] for i in idx]
        vecs = nt_embed.embed_batch(batch, tok, model, device=args.device, max_tokens=mt)
        emb[idx] = vecs.to(torch.float16).cpu().numpy()
        done = start + len(idx)
        if time.time() - last > 15 or done >= N:
            rate = done / max(time.time() - t_start, 1e-9)
            eta = (N - done) / max(rate, 1e-9)
            log(f"  embedded {done:,}/{N:,} ({100*done/N:.1f}%) | {rate:.0f} seq/s | ETA {eta/60:.1f} min")
            emb.flush()
            json.dump({"N": N, "done": done, "model": model_id}, open(prog_path, "w"))
            last = time.time()
    emb.flush()

    # ---- 6. write index + members + manifest ----------------------------------
    log("writing index.parquet / members.parquet / manifest.json")
    for r in rows:
        r["n_genomes"] = len(r["genomes"])
        del r["genomes"]
    idx_df = pd.DataFrame(rows)[
        ["row_id", "md5", "seq_len", "multiplicity", "n_genomes",
         "feature_id", "organism", "genome_id", "taxon_id", "genome_name", "atypical_len"]
    ]
    idx_df.to_parquet(args.outdir / "index.parquet", index=False)
    pd.DataFrame(members).to_parquet(args.outdir / "members.parquet", index=False)

    manifest = {
        "model": model_id, "embed_dim": nt_embed.EMBED_DIM,
        "pooling": "masked-mean of last hidden layer", "normalization": "L2 (unit norm)",
        "dtype_stored": "float16", "similarity": "cosine == dot product of stored vectors",
        "n_input_sequences": len(fasta), "n_unique_sequences": N,
        "n_skipped_below_min_len": n_skipped, "dedup_factor": round(len(fasta) / max(N, 1), 3),
        "batch_size": bs, "max_tokens": mt, "min_len": args.min_len, "atypical_len": args.atypical_len,
        "n_atypical_len": int(idx_df["atypical_len"].sum()),
        "seq_len_stats": {k: float(v) for k, v in idx_df["seq_len"].describe().items()},
        "multiplicity_stats": {k: float(v) for k, v in idx_df["multiplicity"].describe().items()},
        "device": args.device, "transformers": _ver("transformers"), "torch": _ver("torch"),
        "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "wall_seconds": round(time.time() - t_start, 1),
        "inputs": {"fasta": str(args.fasta), "metadata": str(args.metadata)},
    }
    json.dump(manifest, open(args.outdir / "manifest.json", "w"), indent=2)
    log(f"DONE in {(time.time()-t_start)/60:.1f} min | {N:,} unique vectors -> {emb_path}")
    print(json.dumps(manifest, indent=2))


def _ver(pkg):
    import importlib.metadata as m
    try:
        return m.version(pkg)
    except Exception:
        return None


if __name__ == "__main__":
    main()
