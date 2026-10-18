#!/usr/bin/env python
"""Cross-panel comparison: taxonomy-stratified panel vs k-means-over-embeddings panel.

The two panels differ in composition as well as in selection rule, so a crude
difference confounds the two. Every rate is therefore reported as a LADDER:

  crude            -- as measured, no adjustment
  common support   -- restricted to genomes both panels could have drawn
                      (full-length >=1400bp 16S, admissible domain+phylum)
  domain-std       -- direct standardization to a common domain mix
  phylum-std       -- direct standardization to a common phylum mix

Only the phylum-standardized column is attributable to the SELECTION SCHEME. The
review measured FullLength16S at +0.155 crude but +0.036 standardized, and
V4-V5_944R flipping sign, so the crude column is reported to be discounted, not
believed.

Comparisons deliberately NOT made (each would be meaningless across panels):
  * population F1 -- its denominator IS the panel
  * MAG resolve rate from the panel-as-sketch arm -- 10,000 vs 5,048 distinct taxids
    makes the k-means sketch roughly half as deep a reference; use the GTDB arm
  * per-phylum archaeal figures -- k-means phyla are a strict subset (142 of 206)
  * the original's structural 0.000 species accuracy beside any k-means number
  * amplicon-level counts -- operons per genome differ (1.561 vs 2.298); organism-macro only

Joins are on genome_id. NEVER on best_16s_md5: the two panels hash different cases
(upper vs raw), so 0 of the 804 shared genomes agree on that key.
"""
import argparse, collections, csv, json, math, os, sys
sys.path.insert(0, os.path.expanduser("~/Documents/prFBA/region_validation/scripts"))
from _common import wilson_ci

ROOT = os.path.expanduser("~/Documents/prFBA")
RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


def load_truth(p):
    out = {}
    for v in json.load(open(p)).values():
        out.setdefault(v["src_genome_id"], v)
    return out


def load_genecap(p, region, mode="exclude_0.987"):
    """src_genome_id -> (outcome, f1, recall) for one region/regime."""
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh):
            if r["self_mode"] != mode or r["region"] != region:
                continue
            # Match score_genecap's denominator convention: truth_empty /
            # pgfam_missing / no_hits are poisoned cells (the question could not be
            # asked), excluded from both numerator and denominator. pred_empty IS
            # kept and counted as not-covered: a reference was found but predicted
            # nothing, which is a real zero. Counting truth_empty would favour
            # k-means, which has 36 of them against the original's 83.
            if r["outcome"] in ("truth_empty", "pgfam_missing", "no_hits"):
                continue
            def num(x):
                try: return float(x)
                except Exception: return None
            out[r["src_genome_id"]] = (r["outcome"], num(r["f1"]), num(r["recall"]))
    return out


def standardize(values_by_stratum, weights):
    """Direct standardization: sum_s w_s * rate_s over strata present in BOTH panels."""
    num = den = 0.0
    for s, w in weights.items():
        v = values_by_stratum.get(s)
        if v is None or v[1] == 0:
            continue
        num += w * (v[0] / v[1]); den += w
    return (num / den) if den else None


def rate_by(d, keyfn, okfn):
    agg = collections.defaultdict(lambda: [0, 0])
    for g, rec in d.items():
        k = keyfn(g)
        if k is None:
            continue
        agg[k][1] += 1
        if okfn(rec):
            agg[k][0] += 1
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orig-truth", default=ROOT + "/region_validation/data/benchmark_truth.json")
    ap.add_argument("--km-truth", default=ROOT + "/region_validation_kmeans/data/benchmark_truth.json")
    ap.add_argument("--orig-genecap", default=ROOT + "/region_validation/data/gene_capture_cells.csv")
    ap.add_argument("--km-genecap", default=ROOT + "/region_validation_kmeans/data/gene_capture_cells.csv")
    ap.add_argument("--out", default=ROOT + "/two_panel/panel_comparison.json")
    a = ap.parse_args()

    ot, kt = load_truth(a.orig_truth), load_truth(a.km_truth)
    shared = set(ot) & set(kt)
    print("sources: original %d, k-means %d, shared %d" % (len(ot), len(kt), len(shared)), flush=True)

    # a common phylum mix: strata present in both panels, weighted by the pooled count
    ophy = collections.Counter(v["truth_lineage"].get("Phylum") for v in ot.values())
    kphy = collections.Counter(v["truth_lineage"].get("Phylum") for v in kt.values())
    common_phyla = {p for p in (set(ophy) & set(kphy)) if p}
    wphy = {p: ophy[p] + kphy[p] for p in common_phyla}
    odom = collections.Counter(v.get("domain") for v in ot.values())
    kdom = collections.Counter(v.get("domain") for v in kt.values())
    common_dom = {d for d in (set(odom) & set(kdom)) if d}
    wdom = {d: odom[d] + kdom[d] for d in common_dom}
    print("common phyla %d (orig %d, km %d); common domains %s"
          % (len(common_phyla), len(ophy), len(kphy), sorted(common_dom)), flush=True)

    res = {"meta": dict(n_orig=len(ot), n_km=len(kt), n_shared=len(shared),
                        n_common_phyla=len(common_phyla))}
    regions = ["FullLength16S", "V1-V2", "V1-V3", "V3-V4", "V4", "V4-V5", "V4-V5_944R", "V7-V9"]
    rows = []
    for region in regions:
        og, kg = load_genecap(a.orig_genecap, region), load_genecap(a.km_genecap, region)
        if not kg:
            continue
        row = {"region": region}
        for lab, d, truth in (("orig", og, ot), ("km", kg, kt)):
            ok = lambda r: r[0] == "ok"
            crude = rate_by(d, lambda g: "all", ok)["all"]
            dom = rate_by(d, lambda g: (truth.get(g) or {}).get("domain"), ok)
            phy = rate_by(d, lambda g: (truth.get(g) or {}).get("truth_lineage", {}).get("Phylum"), ok)
            row[lab + "_n"] = crude[1]
            row[lab + "_crude"] = round(crude[0] / crude[1], 4) if crude[1] else None
            row[lab + "_domstd"] = standardize(dom, wdom)
            row[lab + "_phystd"] = standardize(phy, wphy)
            # conditional gene-capture F1, organism-macro, over answered organisms
            vals = [r[1] for r in d.values() if r[0] == "ok" and r[1] is not None]
            row[lab + "_condF1"] = round(sum(vals) / len(vals), 4) if vals else None
        for k in ("crude", "domstd", "phystd", "condF1"):
            a_, b_ = row.get("km_" + k), row.get("orig_" + k)
            row["delta_" + k] = round(a_ - b_, 4) if (a_ is not None and b_ is not None) else None
        rows.append(row)
        print("%-13s coverage crude %s/%s  phylum-std %s/%s  condF1 %s/%s"
              % (region, row["km_crude"], row["orig_crude"],
                 None if row["km_phystd"] is None else round(row["km_phystd"], 4),
                 None if row["orig_phystd"] is None else round(row["orig_phystd"], 4),
                 row["km_condF1"], row["orig_condF1"]), flush=True)
    res["regions"] = rows
    json.dump(res, open(a.out, "w"), indent=1)
    print("wrote", a.out, flush=True)


if __name__ == "__main__":
    main()
