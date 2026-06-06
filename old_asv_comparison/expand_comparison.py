#!/usr/bin/env python
"""
expand_comparison.py — extend the old-vs-embedding comparison.

PART 1 (fast): agreement of the embedding mapping with the prior ASV_genomeIDs.json mapping at cosine
  cutoffs k in {1,5,20,50,100} (embedding hits expanded to all genomes sharing each matched insert),
  at genome / taxon / genus level. Answers "what % of ASVs are captured by the top-k cosine hits?".

PART 2 (alignment): for each ASV, the true % sequence identity of the TOP-5 hits of
  (a) the prior all-vs-all method (its 5 most-frequent assigned genomes) and
  (b) the cosine embedding method (its 5 nearest inserts),
  to test whether the embedding's top-5 are as sequence-identical as the alignment method's top-5.

Inputs : asv_top100_hits.json, ASV_genomeIDs.json, dna-sequences.fasta, v4v5_refs.json, v4v5_store/members.parquet
Outputs: comparison_stats_topk.json, old_vs_embedding_topk.csv, identity_top5_comparison.json, identity_top5_per_asv.csv
"""
from __future__ import annotations
import json, time
from pathlib import Path
from collections import Counter
import numpy as np, pandas as pd, taxopy
from Bio.Align import PairwiseAligner

REPO = Path("/home/freiburger/Documents/codiffusion_bioreactor")
PRFBA = Path("/home/freiburger/Documents/prFBA")
STORE = PRFBA / "v4v5_store"
CUTOFFS = [1, 5, 20, 50, 100]

aln = PairwiseAligner(); aln.mode = "local"; aln.match_score = 1; aln.mismatch_score = 0
aln.open_gap_score = -1; aln.extend_gap_score = -0.5
def pid(a, b):
    if not a or not b: return np.nan
    c = aln.align(a, b)[0].counts(); d = c.identities + c.mismatches + c.gaps
    return 100.0 * c.identities / d if d else np.nan

def gid_of(fid):
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
    hits = json.load(open("asv_top100_hits.json"))
    old = json.load(open(REPO / "modeling_files" / "ASV_genomeIDs.json"))
    asv_seq = read_fasta(REPO / "model_inputs" / "dna-sequences.fasta")
    mem = pd.read_parquet(STORE / "members.parquet"); mem["genome"] = mem["feature_id"].map(gid_of)
    row2genomes = mem.groupby("row_id")["genome"].agg(set).to_dict()
    feat2row = dict(zip(mem["feature_id"], mem["row_id"]))
    nd, na, mg = PRFBA / "nodes.dmp", PRFBA / "names.dmp", PRFBA / "merged.dmp"
    db = taxopy.TaxDb(nodes_dmp=str(nd), names_dmp=str(na), merged_dmp=str(mg) if mg.exists() else None)
    gcache = {}
    def genus(t):
        if t not in gcache:
            try: gcache[t] = taxopy.Taxon(int(t), db).rank_name_dictionary.get("genus", "").lower()
            except Exception: gcache[t] = ""
        return gcache[t]

    # ---------- PART 1: top-k agreement curve ----------
    mapped = [a for a in old if a in hits]
    rows = []
    for asv in mapped:
        old_g = set(map(str, old[asv])); old_t = {g.split(".")[0] for g in old_g}
        old_gen = {genus(t) for t in old_t} - {""}
        top = hits[asv]["top100"]
        gacc = set(); r = {"asv": asv, "best_cosine": hits[asv]["best_cosine"], "n_old_genomes": len(old_g)}
        prev = 0
        for k in CUTOFFS:
            for h in top[prev:k]:
                rr = feat2row.get(h["feature_id"]); gacc |= row2genomes.get(rr, {h["genome_id"]}) if rr is not None else {h["genome_id"]}
            prev = k
            tacc = {g.split(".")[0] for g in gacc}; genacc = {genus(t) for t in tacc} - {""}
            r[f"genome_top{k}"] = bool(old_g & gacc)
            r[f"taxon_top{k}"] = bool(old_t & tacc)
            r[f"genus_top{k}"] = (bool(old_gen & genacc) if old_gen else None)
        rows.append(r)
    df = pd.DataFrame(rows); df.to_csv("old_vs_embedding_topk.csv", index=False)
    curve = {"n_compared": len(df), "cutoffs": {}}
    for k in CUTOFFS:
        curve["cutoffs"][str(k)] = {
            "genome": round(float(df[f"genome_top{k}"].mean()), 4),
            "taxon": round(float(df[f"taxon_top{k}"].mean()), 4),
            "genus": round(float(df[f"genus_top{k}"].dropna().mean()), 4),
        }
    json.dump(curve, open("comparison_stats_topk.json", "w"), indent=2)
    print("PART 1 — top-k agreement (% of", len(df), "ASVs captured):")
    print(f"{'k':>4} {'genome':>8} {'taxon':>8} {'genus':>8}")
    for k in CUTOFFS:
        c = curve["cutoffs"][str(k)]; print(f"{k:>4} {c['genome']:>8.3f} {c['taxon']:>8.3f} {c['genus']:>8.3f}")

    # ---------- PART 2: top-5 identity comparison ----------
    print("\nloading reference inserts for identity scoring ...", flush=True)
    refs = json.load(open(PRFBA / "v4v5_refs.json"))
    feat2ins = {}; genome2feats = {}
    for header, ins in refs.items():
        fid = header.split("|")[1].strip() if "|" in header else header
        feat2ins[fid] = ins; genome2feats.setdefault(gid_of(fid), []).append(fid)
    del refs
    t0 = time.time(); prows = []
    for asv in mapped:
        s = asv_seq.get(asv)
        if not s: continue
        old_top5 = [g for g, _ in Counter(map(str, old[asv])).most_common(5)]      # most-frequent old genomes
        emb_top5 = hits[asv]["top100"][:5]
        old_ids = []
        for g in old_top5:
            ins_list = [feat2ins[f] for f in genome2feats.get(g, []) if f in feat2ins][:3]
            old_ids.append(max((pid(s, x) for x in ins_list), default=np.nan))
        emb_ids = [pid(s, feat2ins.get(h["feature_id"], "")) for h in emb_top5]
        prows.append({"asv": asv,
                      "old_top5_best_pid": np.nanmax(old_ids) if np.any(~np.isnan(old_ids)) else np.nan,
                      "old_top5_mean_pid": np.nanmean(old_ids) if np.any(~np.isnan(old_ids)) else np.nan,
                      "emb_top5_best_pid": np.nanmax(emb_ids) if np.any(~np.isnan(emb_ids)) else np.nan,
                      "emb_top5_mean_pid": np.nanmean(emb_ids) if np.any(~np.isnan(emb_ids)) else np.nan})
    pdf = pd.DataFrame(prows); pdf.to_csv("identity_top5_per_asv.csv", index=False)
    both = pdf.dropna(subset=["old_top5_best_pid", "emb_top5_best_pid"])
    out = {
        "n_asv": int(len(pdf)), "n_with_both": int(len(both)),
        "old_all_vs_all_top5": {"best_mean": round(float(both.old_top5_best_pid.mean()), 2),
                                "best_median": round(float(both.old_top5_best_pid.median()), 2),
                                "mean_of_5_mean": round(float(both.old_top5_mean_pid.mean()), 2)},
        "embedding_cosine_top5": {"best_mean": round(float(both.emb_top5_best_pid.mean()), 2),
                                  "best_median": round(float(both.emb_top5_best_pid.median()), 2),
                                  "mean_of_5_mean": round(float(both.emb_top5_mean_pid.mean()), 2)},
        "emb_best_ge_old_best_pct": round(100 * float((both.emb_top5_best_pid >= both.old_top5_best_pid - 1e-9).mean()), 1),
        "emb_best_minus_old_best_mean": round(float((both.emb_top5_best_pid - both.old_top5_best_pid).mean()), 2),
        "wall_seconds": round(time.time() - t0, 1),
    }
    json.dump(out, open("identity_top5_comparison.json", "w"), indent=2)
    print("\nPART 2 — % identity of top-5 hits (n_with_both =", out["n_with_both"], "):")
    print(json.dumps({k: v for k, v in out.items() if k not in ("wall_seconds",)}, indent=2))


if __name__ == "__main__":
    main()
