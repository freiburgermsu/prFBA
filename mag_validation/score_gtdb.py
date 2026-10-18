#!/usr/bin/env python
"""Score the dense-reference (GTDB r226) arm.

Only quantities that need no accession->BV-BRC mapping: ANI, aligned fraction,
novelty class, and abstention under include / exclude(ANI >= 95). Taxonomy and
gene capture live in the panel arm, where every reference is a BV-BRC genome with
an NCBI lineage and cached PGFams.

Exclusion here is ANI-only (>= 95, the ANI species boundary). The panel arm can
also drop same-taxon references because its references carry BV-BRC taxon ids;
GTDB accessions do not, so the ANI criterion stands alone -- and since 95 ANI *is*
the species definition, it is the more principled of the two anyway.
"""
import argparse, collections, csv, json, math, os, sys
sys.path.insert(0, os.path.expanduser("~/Documents/prFBA"))
sys.path.insert(0, os.path.expanduser("~/Documents/prFBA/region_validation/scripts"))
from mag_skani_hits import ANI_KNOBS, best_hit, classify
from _common import wilson_ci

TAU = 95.0
GATE = ANI_KNOBS["cov_min"] * 100.0
FLOOR = ANI_KNOBS["family_floor"] * 100.0     # 83.0 -> below this the selector abstains


def med(v):
    v = sorted(x for x in v if x is not None)
    return v[len(v) // 2] if v else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", action="append", required=True)
    ap.add_argument("--manifest", action="append", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    meta = {}
    for mp in a.manifest:
        for r in json.load(open(mp))["rows"]:
            meta[r["mag"]] = r
    hits = collections.defaultdict(list)
    for hp in a.hits:
        with open(hp) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                q = os.path.basename(row["Query_file"])[:-4]
                hits[q].append(dict(ani=float(row["ANI"]),
                                    af_query=float(row["Align_fraction_query"]),
                                    af_ref=float(row["Align_fraction_ref"]),
                                    _dbword="GTDB"))
    print("queries=%d  synthetic MAGs=%d" % (len(hits), len(meta)), flush=True)

    w = csv.writer(open(os.path.join(a.outdir, "gtdb_cells.csv"), "w", newline=""))
    w.writerow(["mag", "genome_id", "completeness", "contamination", "regime",
                "n_hits", "n_af_qualified", "best_ani", "best_af_query", "best_af_ref",
                "ani_class", "resolves"])
    agg = collections.defaultdict(list)
    for mag, m in sorted(meta.items()):
        hh = hits.get(mag, [])
        for regime, pool in (("include", hh), ("exclude", [h for h in hh if h["ani"] < TAU])):
            b = best_hit(pool, GATE) if pool else None
            nq = sum(1 for h in pool if h["af_query"] >= GATE)
            cls, _ = classify(b or {}, len(pool), GATE)
            # the selector resolves only if the AF gate passes AND ANI clears the floor
            resolves = bool(b and b["af_pass"] and b["ani"] >= FLOOR)
            w.writerow([mag, m["genome_id"], m["completeness_target"], m["contamination_target"],
                        regime, len(pool), nq,
                        (b or {}).get("ani"), (b or {}).get("af_query"), (b or {}).get("af_ref"),
                        cls, int(resolves)])
            agg[(m["completeness_target"], m["contamination_target"], regime)].append(
                dict(nq=nq, b=b, cls=cls, resolves=resolves, n=len(pool)))

    summary = {}
    for key, rows in sorted(agg.items()):
        comp, contam, regime = key
        n = len(rows)
        k = sum(1 for r in rows if r["resolves"])
        lo, hi = wilson_ci(k, n)
        qual = [r for r in rows if r["b"] and r["b"]["af_pass"]]
        summary["c%g_x%g_%s" % (comp, contam, regime)] = dict(
            completeness=comp, contamination=contam, regime=regime, n=n,
            resolves=k, resolve_rate=k / n if n else None, ci=[lo, hi],
            median_best_ani=med([r["b"]["ani"] for r in qual]),
            median_af_query=med([r["b"]["af_query"] for r in qual]),
            median_af_ref=med([r["b"]["af_ref"] for r in qual]),
            median_n_af_qualified=med([r["nq"] for r in rows]),
            mean_n_af_qualified=round(sum(r["nq"] for r in rows) / n, 2) if n else None,
            frac_zero_af_qualified=round(sum(1 for r in rows if r["nq"] == 0) / n, 4) if n else None,
            ani_classes=dict(collections.Counter(r["cls"] for r in rows)))
    json.dump(dict(tau_ani=TAU, af_gate=GATE, ani_floor=FLOOR, summary=summary),
              open(os.path.join(a.outdir, "gtdb_summary.json"), "w"), indent=1)
    print("wrote %s/gtdb_summary.json" % a.outdir, flush=True)


if __name__ == "__main__":
    main()
