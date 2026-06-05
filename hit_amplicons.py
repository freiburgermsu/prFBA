#!/usr/bin/env python
"""
hit_amplicons.py — reference experimental 16S amplicons (ASVs) against the BV-BRC embedding space.

The ASVs are V4–V5 amplicons (~373 bp). We embed them with the same NT-v2 encoder and match them by
cosine similarity against a **region-matched** reference store (V4–V5 inserts excised from BV-BRC by
insilico_pcr.py + embedded by embed_16s.py) — apples-to-apples. For comparison we also score them
against the **full-length** reference store (naive), to quantify how much region-matching matters.

Outputs (to --outdir)
  asv_top20_hits.json  : { ASV : { asv_len, midas_taxonomy, best_cosine, novel, top20:[ {rank,cosine,
                          organism, genome_name, taxon_id, genome_id, feature_id, n_genomes, seq_len} ] } }
  asv_summary.csv      : one row per ASV (best hit, cosines, genus/family concordance vs MiDAS, novelty)
  findings_stats.json  : aggregate statistics for the report
"""
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import numpy as np


def read_fasta(path):
    names, seqs, cur, nm = [], [], [], None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if nm is not None: seqs.append("".join(cur))
            nm = line[1:].strip(); names.append(nm); cur = []
        elif line.strip():
            cur.append(line.strip().upper())
    if nm is not None: seqs.append("".join(cur))
    return names, seqs


def genus_of(organism):
    if not isinstance(organism, str) or not organism.strip():
        return ""
    toks = organism.replace("[", "").replace("]", "").split()
    if not toks: return ""
    g = toks[0]
    if g.lower() in ("candidatus", "uncultured", "unclassified") and len(toks) > 1:
        g = toks[1]
    return g.lower()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--amplicons", type=Path, required=True)
    ap.add_argument("--store", type=Path, required=True, help="region-matched (V4-V5) store dir")
    ap.add_argument("--naive-store", type=Path, default=None, help="full-length store dir (for comparison)")
    ap.add_argument("--asv-ids", type=Path, default=None, help="asv_IDs.csv (MiDAS taxonomy)")
    ap.add_argument("--taxonomy", type=Path, default=None, help="taxonomy.csv (structured ranks)")
    ap.add_argument("--topk", type=int, default=20)
    ap.add_argument("--novel-thr", type=float, default=0.90, help="best cosine below this = no confident match")
    ap.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    import nt_embed, torch, pandas as pd
    device = args.device if torch.cuda.is_available() else "cpu"

    # ---- load amplicons + study taxonomy -------------------------------------
    asv_ids, asv_seqs = read_fasta(args.amplicons)
    print(f"[hit] {len(asv_ids)} ASVs (len {min(map(len,asv_seqs))}-{max(map(len,asv_seqs))})")
    midas = {}        # asv -> MiDAS string ; struct -> dict of ranks
    struct = {}
    if args.asv_ids and args.asv_ids.exists():
        a = pd.read_csv(args.asv_ids)
        for _, r in a.iterrows():
            midas[str(r.iloc[0])] = str(r.get("MiDAS Taxonomy", ""))
    if args.taxonomy and args.taxonomy.exists():
        t = pd.read_csv(args.taxonomy).drop_duplicates(subset=t_seqcol(args.taxonomy))
        sc = t_seqcol(args.taxonomy)
        for _, r in t.iterrows():
            struct[str(r[sc])] = {k: str(r.get(k, "")) for k in ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]}

    # ---- load region-matched store + embed amplicons -------------------------
    idx = pd.read_parquet(args.store / "index.parquet")
    E = torch.from_numpy(np.ascontiguousarray(np.load(args.store / "embeddings.f16.npy"))).to(device)
    print(f"[hit] region-matched references: {len(idx):,}")
    tok, model = nt_embed.load_model(device=device)
    Q = embed_all(asv_seqs, tok, model, device)                      # (n_asv, 1024) unit f16

    sims = (Q @ E.T)                                                  # (n_asv, N) cosine
    top = torch.topk(sims.float(), k=min(args.topk, E.shape[0]), dim=1)
    top_vals, top_idx = top.values.cpu().numpy(), top.indices.cpu().numpy()

    # ---- naive (full-length) comparison: best cosine only --------------------
    naive_best = None
    if args.naive_store:
        En = torch.from_numpy(np.ascontiguousarray(np.load(args.naive_store / "embeddings.f16.npy"))).to(device)
        nb = []
        for s in range(0, Q.shape[0], 512):
            nb.append((Q[s:s+512] @ En.T).float().max(1).values.cpu().numpy())
        naive_best = np.concatenate(nb)
        del En; torch.cuda.empty_cache()

    # ---- assemble per-ASV records --------------------------------------------
    cols = ["organism", "genome_name", "taxon_id", "genome_id", "feature_id", "n_genomes", "seq_len", "multiplicity"]
    out, summary = {}, []
    for i, asv in enumerate(asv_ids):
        hits = []
        for rank, (ci, cv) in enumerate(zip(top_idx[i], top_vals[i]), 1):
            m = idx.iloc[int(ci)]
            hits.append({"rank": rank, "cosine": round(float(cv), 4),
                         "organism": _s(m["organism"]), "genome_name": _s(m["genome_name"]),
                         "taxon_id": _n(m["taxon_id"]), "genome_id": _s(m["genome_id"]),
                         "feature_id": _s(m["feature_id"]), "n_genomes": int(m["n_genomes"]),
                         "ref_seq_len": int(m["seq_len"])})
        best = hits[0]
        mtax = midas.get(asv, "")
        midas_genus = _midas_genus(mtax, struct.get(asv))
        midas_family = _midas_rank(struct.get(asv), "Family") or _midas_token(mtax, -2)
        bg = genus_of(best["organism"])
        gconc = bool(midas_genus and bg and midas_genus.lower() == bg)
        # family concordance: best-hit family unknown from organism; use genus->family is hard, so compare at genus only + MiDAS family vs any top hit genus chain skipped
        out[asv] = {"asv_len": len(asv_seqs[i]), "midas_taxonomy": mtax,
                    "best_cosine": best["cosine"], "novel": best["cosine"] < args.novel_thr,
                    "top%d" % args.topk: hits}
        summary.append({"asv": asv, "asv_len": len(asv_seqs[i]),
                        "midas_genus": midas_genus, "best_hit_organism": best["organism"],
                        "best_cosine": best["cosine"],
                        "naive_best_cosine": round(float(naive_best[i]), 4) if naive_best is not None else None,
                        "genus_concordant": gconc, "novel": best["cosine"] < args.novel_thr,
                        "best_taxon_id": best["taxon_id"], "best_n_genomes": best["n_genomes"]})

    args.outdir.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(args.outdir / "asv_top20_hits.json", "w"), indent=1)
    sdf = pd.DataFrame(summary)
    sdf.to_csv(args.outdir / "asv_summary.csv", index=False)

    # ---- aggregate findings ---------------------------------------------------
    bc = sdf["best_cosine"].values
    stats = {
        "n_asv": len(asv_ids), "region": "V4-V5 (515F/926R)",
        "n_region_matched_refs": int(len(idx)),
        "best_cosine": pctd(bc),
        "novel_fraction(best<%.2f)" % args.novel_thr: float((bc < args.novel_thr).mean()),
        "match_ge_0.97": float((bc >= 0.97).mean()), "match_ge_0.99": float((bc >= 0.99).mean()),
        "genus_concordance_rate(of ASVs with MiDAS genus)":
            float(sdf.loc[sdf["midas_genus"] != "", "genus_concordant"].mean()) if (sdf["midas_genus"] != "").any() else None,
        "n_asv_with_midas_genus": int((sdf["midas_genus"] != "").sum()),
    }
    if naive_best is not None:
        stats["naive_best_cosine"] = pctd(naive_best)
        stats["region_match_cosine_gain_mean"] = round(float((bc - naive_best).mean()), 4)
    json.dump(stats, open(args.outdir / "findings_stats.json", "w"), indent=2)
    print(json.dumps(stats, indent=2))


def embed_all(seqs, tok, model, device, bs=128):
    import torch
    order = np.argsort([len(s) for s in seqs])                       # length-sorted batching
    out = torch.empty((len(seqs), 1024), device=device)
    import nt_embed
    for s in range(0, len(order), bs):
        idx = order[s:s+bs]
        out[idx] = nt_embed.embed_batch([seqs[i] for i in idx], tok, model, device=device).to(out.dtype)
    return out.half()


def t_seqcol(path):
    import pandas as pd
    c = pd.read_csv(path, nrows=1).columns
    return "seq" if "seq" in c else c[0]


def _midas_genus(mtax, st):
    if st and st.get("Genus") and st["Genus"].lower() not in ("nan", ""):
        return st["Genus"].lower()
    toks = [t for t in str(mtax).split() if t]
    return toks[-1].lower() if toks else ""


def _midas_rank(st, rank):
    return st[rank].lower() if st and st.get(rank) and st[rank].lower() not in ("nan", "") else ""


def _midas_token(mtax, i):
    toks = [t for t in str(mtax).split() if t]
    return toks[i].lower() if len(toks) >= abs(i) else ""


def _s(v):
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else str(v)


def _n(v):
    try:
        return None if v is None or (isinstance(v, float) and np.isnan(v)) else int(float(v))
    except Exception:
        return None


def pctd(a):
    a = np.asarray(a, float)
    return {"mean": round(float(a.mean()), 4), "median": round(float(np.median(a)), 4),
            "p5": round(float(np.percentile(a, 5)), 4), "p25": round(float(np.percentile(a, 25)), 4),
            "p75": round(float(np.percentile(a, 75)), 4), "p95": round(float(np.percentile(a, 95)), 4),
            "min": round(float(a.min()), 4), "max": round(float(a.max()), 4)}


if __name__ == "__main__":
    main()
