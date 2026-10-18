#!/usr/bin/env python
"""Three-panel comparison: taxonomy-stratified vs k-means vs uniform-random.

Generalises two_panel/compare_panels.py from two panels to three, adds the taxonomic-
accuracy endpoint beside gene capture, implements the "common support" rung that the
two-panel docstring promised but never computed, and attaches genome-bootstrap CIs.

Every rate is a LADDER (only the last rung is attributable to the SELECTION SCHEME):
  crude         as measured
  common        restricted to genomes every panel could have drawn: a full-length
                (>= 1400 bp) best 16S. The original panel backfilled down to 1200 bp;
                the k-means and random panels never did, so the backfilled genomes have
                no counterpart in the other two.
  domain-std    direct standardisation to a pooled domain mix over domains in all three
  phylum-std    direct standardisation to a pooled phylum mix over phyla in all three

Endpoints, per region x regime (include | exclude_0.987), organism-macro:
  genus coverage, genus conditional accuracy        (anchor predictor, the paper's)
  family coverage, family conditional accuracy
  gene-capture coverage, gene-capture conditional F1

Deliberately NOT compared (each meaningless across panels; see TWO_PANEL_RUNBOOK.md):
  population accuracy / population F1  -- the denominator IS the panel
  species accuracy                     -- the original panel holds one genome per species,
                                          so self-excluded species accuracy is 0 by construction
  per-phylum archaeal figures          -- the panels' phyla are not nested the same way
  amplicon-level counts                -- joins are on genome_id, never on best_16s_md5

A rank whose truth is undefined for a genome is unscoreable for that genome and is
excluded from that rank's denominator (score_taxacc.py's truth_defined rule).
"""
import argparse, collections, csv, json, math, os, sys
import numpy as np
import pandas as pd

ROOT = os.path.expanduser("~/Documents/prFBA")
PANELS = {
    "original": ROOT + "/region_validation",
    "kmeans":   ROOT + "/region_validation_kmeans",
    "random":   ROOT + "/region_validation_random",
}
REGIONS = ["FullLength16S", "V1-V2", "V1-V3", "V3-V4", "V4", "V4-V5", "V4-V5_944R", "V7-V9"]
REGIMES = ["include", "exclude_0.987"]
FULL_LEN = 1400
B = 1000
SEED = 1729


def load_truth(root):
    """src_genome_id -> {domain, Phylum, Family, Genus, best_len}."""
    sel = {str(e["genome_id"]): e for e in json.load(open(root + "/data/selection_10k.json"))["genomes"]}
    out = {}
    for v in json.load(open(root + "/data/benchmark_truth.json")).values():
        g = str(v["src_genome_id"])
        if g in out:
            continue
        lin = v.get("truth_lineage") or {}
        out[g] = dict(domain=v.get("domain"), Phylum=lin.get("Phylum"),
                      Family=lin.get("Family"), Genus=lin.get("Genus"),
                      best_len=int((sel.get(g) or {}).get("best_16s_len") or 0))
    return out


def load_taxacc(root):
    d = pd.read_csv(root + "/data/asv_region_concordance.csv",
                    usecols=["src_genome_id", "region", "self_mode", "predictor",
                             "abstained", "correct_Family", "correct_Genus"],
                    dtype={"src_genome_id": str})
    return d[(d.predictor == "anchor") & d.self_mode.isin(REGIMES)]


def load_genecap(root):
    d = pd.read_csv(root + "/data/gene_capture_cells.csv",
                    usecols=["src_genome_id", "region", "self_mode", "outcome", "f1"],
                    dtype={"src_genome_id": str})
    d = d.copy()
    d["self_mode"] = d.self_mode.replace({"exclude": "exclude_0.987"})
    # score_genecap's denominator: truth_empty / pgfam_missing / no_hits are cells where
    # the question could not be asked; pred_empty is a real zero and stays.
    return d[~d.outcome.isin(["truth_empty", "pgfam_missing", "no_hits"])]


def per_genome_rows(truth, tax, gc, region, regime):
    """One record per scoreable genome: the six endpoint indicators, NaN where unscoreable."""
    rows = {}
    t = tax[(tax.region == region) & (tax.self_mode == regime)]
    for r in t.itertuples(index=False):
        tr = truth.get(r.src_genome_id)
        if not tr:
            continue
        rec = rows.setdefault(r.src_genome_id, dict(g=r.src_genome_id, **tr))
        for rank, col in (("Genus", r.correct_Genus), ("Family", r.correct_Family)):
            if tr.get(rank) in (None, "") or str(tr.get(rank)).startswith("unclassified"):
                rec[rank + "_cov"] = np.nan; rec[rank + "_acc"] = np.nan; continue
            # A call is RANK-SPECIFIC: a family-tier assignment calls the family but not the
            # genus. Given defined truth, correct_<Rank> is NaN exactly when that rank was
            # not called, which is score_taxacc's CovCorr definition. (An earlier version
            # used `abstained == 0` for every rank and overstated self-excluded genus
            # coverage by up to 0.19; caught by the scorer-reproduction check.)
            called = col == col
            rec[rank + "_cov"] = 1.0 if called else 0.0
            rec[rank + "_acc"] = (1.0 if col == 1 else 0.0) if called else np.nan
    g = gc[(gc.region == region) & (gc.self_mode == regime)]
    for r in g.itertuples(index=False):
        tr = truth.get(r.src_genome_id)
        if not tr:
            continue
        rec = rows.setdefault(r.src_genome_id, dict(g=r.src_genome_id, **tr))
        rec["gc_cov"] = 1.0 if r.outcome == "ok" else 0.0
        rec["gc_f1"] = float(r.f1) if (r.outcome == "ok" and r.f1 == r.f1) else np.nan
    return pd.DataFrame(list(rows.values()))


ENDPOINTS = [("genus_coverage", "Genus_cov"), ("genus_cond_acc", "Genus_acc"),
             ("family_coverage", "Family_cov"), ("family_cond_acc", "Family_acc"),
             ("genecap_coverage", "gc_cov"), ("genecap_cond_f1", "gc_f1")]


def std_rate(df, col, key, weights):
    """Direct standardisation over strata present in `weights`."""
    num = den = 0.0
    for s, w in weights.items():
        v = df.loc[df[key] == s, col].dropna()
        if len(v) == 0:
            continue
        num += w * v.mean(); den += w
    return num / den if den else np.nan


def ladder(df, col, wdom, wphy):
    x = df[col].dropna()
    c = df[df.best_len >= FULL_LEN]
    return dict(crude=x.mean() if len(x) else np.nan,
                common=c[col].dropna().mean() if c[col].notna().any() else np.nan,
                domstd=std_rate(c, col, "domain", wdom),
                phystd=std_rate(c, col, "Phylum", wphy),
                n=int(x.size), n_common=int(c[col].notna().sum()))


def boot_phystd(df, col, wphy, rng):
    c = df[(df.best_len >= FULL_LEN)]
    c = c[c[col].notna()]
    if c.empty:
        return np.array([])
    # stratified by phylum so every bootstrap replicate keeps the standardisation strata
    groups = {p: c.loc[c.Phylum == p, col].to_numpy() for p in wphy if (c.Phylum == p).any()}
    wsum = sum(wphy[p] for p in groups)
    out = np.empty(B)
    for b in range(B):
        acc = 0.0
        for p, v in groups.items():
            acc += wphy[p] * v[rng.integers(0, len(v), len(v))].mean()
        out[b] = acc / wsum
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=ROOT + "/three_panel/three_panel_comparison.json")
    ap.add_argument("--regions", default=",".join(REGIONS))
    a = ap.parse_args()
    rng = np.random.default_rng(SEED)

    T = {p: load_truth(r) for p, r in PANELS.items()}
    TAX = {p: load_taxacc(r) for p, r in PANELS.items()}
    GC = {p: load_genecap(r) for p, r in PANELS.items()}
    for p in PANELS:
        print("[load] %-8s truth %5d genomes | taxacc %7d cells | genecap %7d cells"
              % (p, len(T[p]), len(TAX[p]), len(GC[p])), flush=True)

    # standardisation weights: strata present in ALL THREE panels, pooled counts over the
    # common-support genomes
    def counts(p, key):
        return collections.Counter(v[key] for v in T[p].values() if v[key] and v["best_len"] >= FULL_LEN)
    dom = {p: counts(p, "domain") for p in PANELS}
    phy = {p: counts(p, "Phylum") for p in PANELS}
    common_dom = set.intersection(*(set(d) for d in dom.values()))
    common_phy = set.intersection(*(set(d) for d in phy.values()))
    wdom = {d: sum(dom[p][d] for p in PANELS) for d in common_dom}
    wphy = {f: sum(phy[p][f] for p in PANELS) for f in common_phy}
    cover = {p: sum(phy[p][f] for f in common_phy) / max(1, sum(phy[p].values())) for p in PANELS}
    print("[std] common domains %s; common phyla %d (orig %d, km %d, rnd %d); "
          "fraction of each panel's common-support genomes inside them: %s"
          % (sorted(common_dom), len(common_phy), len(phy["original"]), len(phy["kmeans"]),
             len(phy["random"]), {p: round(v, 4) for p, v in cover.items()}), flush=True)

    res = dict(meta=dict(panels=list(PANELS), regimes=REGIMES, full_len=FULL_LEN,
                         bootstrap=B, seed=SEED, predictor="anchor",
                         n_common_phyla=len(common_phy), common_domains=sorted(common_dom),
                         std_coverage=cover,
                         n_common_support={p: sum(1 for v in T[p].values() if v["best_len"] >= FULL_LEN)
                                           for p in PANELS},
                         n_truth={p: len(T[p]) for p in PANELS}),
               cells=[])
    for region in a.regions.split(","):
        for regime in REGIMES:
            dfs = {p: per_genome_rows(T[p], TAX[p], GC[p], region, regime) for p in PANELS}
            for ep, col in ENDPOINTS:
                cell = dict(region=region, regime=regime, endpoint=ep)
                boots = {}
                for p in PANELS:
                    if col not in dfs[p]:
                        continue
                    L = ladder(dfs[p], col, wdom, wphy)
                    bs = boot_phystd(dfs[p], col, wphy, rng)
                    boots[p] = bs
                    L["phystd_ci"] = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))] if bs.size else None
                    cell[p] = {k: (None if (isinstance(v, float) and math.isnan(v)) else v) for k, v in L.items()}
                for x, y in (("kmeans", "original"), ("random", "original"), ("random", "kmeans")):
                    if x in boots and y in boots and boots[x].size and boots[y].size:
                        d = boots[x] - boots[y]
                        est = cell[x]["phystd"] - cell[y]["phystd"]
                        lo, hi = np.percentile(d, [2.5, 97.5])
                        cell["delta_%s_minus_%s" % (x, y)] = dict(
                            est=round(float(est), 4), ci=[round(float(lo), 4), round(float(hi), 4)],
                            excludes_zero=bool(lo > 0 or hi < 0))
                res["cells"].append(cell)
        print("[done] %s" % region, flush=True)
    json.dump(res, open(a.out, "w"), indent=1, default=float)
    print("wrote", a.out, flush=True)


if __name__ == "__main__":
    main()
