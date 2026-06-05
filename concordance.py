#!/usr/bin/env python
"""
concordance.py — taxonomic concordance of amplicon embedding-hits vs the study's MiDAS taxonomy.

For each ASV we compare the taxonomy of its BV-BRC embedding hit(s) against the MiDAS assignment the
study made (taxonomy.csv). Because V4-V5 often cannot resolve a single species (many references tie
near cosine 1.0) and because genus names diverge between NCBI/BV-BRC and MiDAS, we report concordance
at every rank (genus→phylum), for the single best hit and for "any of the top-20", and we weight by
ASV abundance. NCBI phylum names are normalized to their classic/MiDAS synonyms before comparison.

Inputs : asv_top20_hits.json (from hit_amplicons.py), taxonomy.csv (MiDAS ranks + rel_ab)
Outputs: concordance_by_rank.json, asv_concordance.csv
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np, pandas as pd, taxopy

REPO = Path("/home/freiburger/Documents/codiffusion_bioreactor")

# NCBI (new) -> classic/MiDAS phylum names so the comparison is fair
PHYLUM_SYN = {
    "bacillota": "firmicutes", "pseudomonadota": "proteobacteria", "actinomycetota": "actinobacteria",
    "bacteroidota": "bacteroidetes", "verrucomicrobiota": "verrucomicrobia", "spirochaetota": "spirochaetes",
    "chloroflexota": "chloroflexi", "synergistota": "synergistetes", "fusobacteriota": "fusobacteria",
    "deinococcota": "deinococcus-thermus", "lentisphaerota": "lentisphaerae", "planctomycetota": "planctomycetes",
    "campylobacterota": "campylobacterota", "thermodesulfobacteriota": "thermodesulfobacteriota",
    # archaea: NCBI split Euryarchaeota; MiDAS keeps Euryarchaeota
    "methanobacteriota": "euryarchaeota", "halobacteriota": "euryarchaeota", "thermoplasmatota": "euryarchaeota",
    "methanobacteriati": "euryarchaeota",
}
def norm(s):
    s = (s or "").lower().strip()
    return PHYLUM_SYN.get(s, s)

def is_placeholder(s):
    s = (s or "").lower()
    return s in ("", "nan") or s.startswith("midas") or s.startswith("unclassified")

RANKS = ["genus", "family", "order", "class", "phylum"]


def main():
    hits = json.load(open("asv_top20_hits.json"))
    tax = pd.read_csv(REPO / "model_inputs" / "taxonomy.csv")
    seqcol = "seq" if "seq" in tax.columns else tax.columns[0]
    abund = tax.groupby(seqcol)["rel_ab"].max().to_dict() if "rel_ab" in tax.columns else {}
    midas = {}
    for _, r in tax.drop_duplicates(subset=seqcol).iterrows():
        midas[str(r[seqcol])] = {rk: str(r.get(rk.capitalize(), "")) for rk in RANKS}

    # map every top-20 taxon_id -> NCBI lineage ranks
    db = taxopy.TaxDb(keep_files=True)
    tids = {int(x["taxon_id"]) for h in hits.values() for x in h["top20"] if x["taxon_id"] is not None}
    lin = {}
    for t in tids:
        try:
            d = taxopy.Taxon(t, db).rank_name_dictionary
            lin[t] = {rk: d.get(rk, "") for rk in RANKS}
        except Exception:
            lin[t] = {rk: "" for rk in RANKS}

    rows = []
    for asv, h in hits.items():
        m = midas.get(asv, {rk: "" for rk in RANKS})
        top = h["top20"]
        best_t = top[0]["taxon_id"]
        best_lin = lin.get(int(best_t), {}) if best_t is not None else {}
        # value sets across the top-20 (for "any" concordance)
        topsets = {rk: {lin.get(int(x["taxon_id"]), {}).get(rk, "") for x in top if x["taxon_id"] is not None}
                   for rk in RANKS}
        row = {"asv": asv, "asv_len": h["asv_len"], "best_cosine": h["best_cosine"],
               "abundance": abund.get(asv, np.nan)}
        for rk in RANKS:
            mv = norm(m[rk]) if rk == "phylum" else (m[rk] or "").lower().strip()
            bv = norm(best_lin.get(rk, "")) if rk == "phylum" else best_lin.get(rk, "").lower()
            tv = {norm(v) for v in topsets[rk]} if rk == "phylum" else {v.lower() for v in topsets[rk]}
            evaluable = not is_placeholder(m[rk])
            row[f"midas_{rk}"] = m[rk]
            row[f"besthit_{rk}"] = best_lin.get(rk, "")
            row[f"{rk}_best_concordant"] = (evaluable and bool(bv) and bv == mv)
            row[f"{rk}_any20_concordant"] = (evaluable and mv in tv)
            row[f"{rk}_evaluable"] = evaluable
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv("asv_concordance.csv", index=False)

    # aggregate rates
    out = {"n_asv": len(df), "ranks": {}}
    w = df["abundance"].fillna(0).values
    for rk in RANKS:
        ev = df[f"{rk}_evaluable"].values
        nb = df.loc[ev, f"{rk}_best_concordant"]
        na = df.loc[ev, f"{rk}_any20_concordant"]
        wb = (df.loc[ev, f"{rk}_best_concordant"].values * w[ev]).sum() / max(w[ev].sum(), 1e-9)
        out["ranks"][rk] = {
            "n_evaluable": int(ev.sum()),
            "best_hit_concordance": round(float(nb.mean()), 4) if len(nb) else None,
            "any_top20_concordance": round(float(na.mean()), 4) if len(na) else None,
            "abundance_weighted_best": round(float(wb), 4),
        }
    out["placeholder_genus_fraction"] = round(float(df["midas_genus"].map(is_placeholder).mean()), 4)
    out["placeholder_family_fraction"] = round(float(df["midas_family"].map(is_placeholder).mean()), 4)
    json.dump(out, open("concordance_by_rank.json", "w"), indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
