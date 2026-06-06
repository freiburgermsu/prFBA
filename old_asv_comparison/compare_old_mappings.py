#!/usr/bin/env python
"""
compare_old_mappings.py — compare embedding-based ASV->genome mappings to the study's old mappings.

The codiffusion study previously assigned each ASV to a set of BV-BRC genomes (modeling_files/
ASV_genomeIDs.json; ASVs with redundant genome sets were collapsed into 1,736 representatives). Here
we map the same ASVs by NT-v2 embedding cosine against the region-matched V4-V5 BV-BRC space and ask
how well the embedding reproduces the old assignment, at three levels:
  * genome  — does any genome the embedding hit (expanded to ALL genomes sharing the matched V4-V5
              insert) appear in the old genome set?
  * taxon   — same at NCBI taxon-id level (genome_id prefix)
  * genus   — same at genus level (taxopy lineage of the taxon)
reported for the single best hit (top-1) and for the top-20.

Inputs : asv_top20_hits.json (this dir), ASV_genomeIDs.json, v4v5_store/{members,index}.parquet
Outputs: old_vs_embedding.csv (per-ASV), comparison_stats.json (aggregate)
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np, pandas as pd, taxopy

HITS = Path("asv_top20_hits.json")
OLD = Path("/home/freiburger/Documents/codiffusion_bioreactor/modeling_files/ASV_genomeIDs.json")
STORE = Path("/home/freiburger/Documents/prFBA/v4v5_store")
PRFBA = Path("/home/freiburger/Documents/prFBA")

def gid_of(feature_id):           # "1505.7.rna.1" -> genome "1505.7"
    p = str(feature_id).split(".")
    return ".".join(p[:2]) if len(p) >= 2 else str(feature_id)


def main():
    hits = json.load(open(HITS))
    old = json.load(open(OLD))
    print(f"embedding hits: {len(hits)} ASVs | old mappings: {len(old)} ASVs")

    # expand each store row -> the full set of genomes sharing that V4-V5 insert
    mem = pd.read_parquet(STORE / "members.parquet")
    mem["genome"] = mem["feature_id"].map(gid_of)
    row2genomes = mem.groupby("row_id")["genome"].agg(set).to_dict()
    feat2row = dict(zip(mem["feature_id"], mem["row_id"]))

    # genus resolution
    nd, na, mg = PRFBA / "nodes.dmp", PRFBA / "names.dmp", PRFBA / "merged.dmp"
    db = taxopy.TaxDb(nodes_dmp=str(nd), names_dmp=str(na), merged_dmp=str(mg) if mg.exists() else None) \
        if nd.exists() else taxopy.TaxDb(keep_files=True)
    genus_cache = {}
    def genus(taxid):
        if taxid not in genus_cache:
            try:
                genus_cache[taxid] = taxopy.Taxon(int(taxid), db).rank_name_dictionary.get("genus", "").lower()
            except Exception:
                genus_cache[taxid] = ""
        return genus_cache[taxid]

    rows = []
    for asv, og in old.items():
        if asv not in hits:
            continue
        old_g = set(map(str, og))
        old_t = {g.split(".")[0] for g in old_g}
        old_genus = {genus(t) for t in old_t} - {""}

        top = hits[asv]["top20"]
        def genomes_of_hit(h):
            r = feat2row.get(h["feature_id"])
            return row2genomes.get(r, {h["genome_id"]}) if r is not None else {h["genome_id"]}
        g1 = genomes_of_hit(top[0])
        gK = set().union(*(genomes_of_hit(h) for h in top))
        t1 = {g.split(".")[0] for g in g1}
        tK = {g.split(".")[0] for g in gK}
        gen1 = {genus(t) for t in t1} - {""}
        genK = {genus(t) for t in tK} - {""}

        rows.append({
            "asv": asv, "best_cosine": hits[asv]["best_cosine"],
            "n_old_genomes": len(old_g), "n_old_taxa": len(old_t),
            "genome_top1": bool(old_g & g1), "genome_top20": bool(old_g & gK),
            "taxon_top1": bool(old_t & t1), "taxon_top20": bool(old_t & tK),
            "genus_top1": bool(old_genus & gen1) if old_genus else None,
            "genus_top20": bool(old_genus & genK) if old_genus else None,
            "genome_jaccard_top20": len(old_g & gK) / max(len(old_g | gK), 1),
            "best_hit_organism": top[0]["organism"], "best_hit_genome": top[0]["genome_id"],
        })

    df = pd.DataFrame(rows)
    df.to_csv("old_vs_embedding.csv", index=False)
    def rate(col):
        s = df[col].dropna()
        return round(float(s.mean()), 4) if len(s) else None
    stats = {
        "n_compared": len(df),
        "genome_level": {"top1": rate("genome_top1"), "top20": rate("genome_top20")},
        "taxon_level": {"top1": rate("taxon_top1"), "top20": rate("taxon_top20")},
        "genus_level": {"top1": rate("genus_top1"), "top20": rate("genus_top20"),
                        "n_evaluable": int(df["genus_top1"].notna().sum())},
        "genome_jaccard_top20": {"mean": round(float(df["genome_jaccard_top20"].mean()), 4),
                                 "median": round(float(df["genome_jaccard_top20"].median()), 4)},
        "old_genomes_per_asv": {"median": float(df["n_old_genomes"].median()),
                                "mean": round(float(df["n_old_genomes"].mean()), 1)},
    }
    json.dump(stats, open("comparison_stats.json", "w"), indent=2)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
