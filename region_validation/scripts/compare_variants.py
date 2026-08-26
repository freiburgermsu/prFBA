#!/usr/bin/env python
"""compare_variants.py — side-by-side comparison of the baseline selection vs the
include-all ("union") variants on the 1e4 self-recovery benchmark: self-recovery,
taxonomic concordance (anchor+consensus, include & exclude_0.987), and gene
capture (macro recall/precision/F1/Jaccard, include & exclude_0.987), per region.

Reads the summary tables produced by score_taxacc/score_genecap for each variant
(suffix '' = baseline, '_ib' = in-band union, '_ir' = all-reliable union) and
writes a markdown report + a flat CSV.

Interpreter: /home/freiburger/Documents/py_venv/bin/python
"""
from __future__ import annotations

import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import DATA  # noqa: E402

VARIANTS = [("legacy (baseline)", ""), ("exact-tie (DEFAULT)", "_et"), ("graduated", "_gr"),
            ("exact-tie+guards", "_etg"), ("in-band 0.005", "_ib")]
REGION_ORDER = ["FullLength16S", "V1-V2", "V1-V3", "V3-V4", "V4", "V4-V5", "V4-V5_944R", "V7-V9"]


def load_json(suffix, stem):
    p = os.path.join(DATA, f"{stem}{suffix}.json")
    return json.load(open(p)) if os.path.exists(p) else None


def self_recovery(suffix):
    """region -> selected_given_top20 (include, domain=all, tier=all), + n-weighted overall."""
    p = os.path.join(DATA, f"self_recovery{suffix}.csv")
    out, num, den = {}, 0.0, 0.0
    anc = {}
    for row in csv.DictReader(open(p)):
        if row["self_mode"] == "include" and row["domain"] == "all" and row["tier"] == "all":
            top20 = int(row["n_self_in_top20"])
            sel = int(row["n_self_selected"])
            out[row["region"]] = sel / top20 if top20 else None
            anc[row["region"]] = int(row["n_self_anchor"]) / top20 if top20 else None
            num += sel
            den += top20
    out["ALL"] = num / den if den else None
    return out, anc


def taxacc(suffix, predictor, self_mode, rank):
    """region -> correct concordance at `rank` for predictor/self_mode (domain=all,tier=all)."""
    s = load_json(suffix, "taxacc_region_summary")
    if not s:
        return {}
    out = {}
    for region, rb in s["regions"].items():
        try:
            out[region] = rb[predictor][self_mode]["all"]["all"]["per_rank"][rank]["correct"]
        except KeyError:
            out[region] = None
    return out


def genecap(suffix, mode, metric):
    """region -> macro `metric` value for self_mode `mode`."""
    s = load_json(suffix, "gene_capture_by_region")
    if not s:
        return {}
    out = {}
    for region, rb in s["regions"].items():
        try:
            out[region] = rb[mode]["macro"][metric]["value"]
        except (KeyError, TypeError):
            out[region] = None
    return out


def genecap_n(suffix, mode):
    s = load_json(suffix, "gene_capture_by_region")
    if not s:
        return {}
    out = {}
    for region, rb in s["regions"].items():
        try:
            out[region] = rb[mode]["macro"]["recall"]["n"]
        except (KeyError, TypeError):
            out[region] = None
    return out


def pct(x):
    return "—" if x is None else f"{100*x:.1f}"


def table(title, getter, regions=REGION_ORDER):
    rows = [f"\n### {title}\n", "| region | " + " | ".join(n for n, _ in VARIANTS) + " |",
            "|" + "---|" * (len(VARIANTS) + 1)]
    series = [getter(suf) for _, suf in VARIANTS]
    for r in regions:
        cells = " | ".join(pct(s.get(r)) for s in series)
        rows.append(f"| {r} | {cells} |")
    return "\n".join(rows)


def main():
    out = ["# Include-all (“union”) selection vs baseline — 1e4 self-recovery benchmark",
           "",
           "Three selection policies, same alignment hits, same metric code:",
           "- **baseline** = ordered reducer pipeline (family-consensus → dedup → marginal-gain → tier-cap).",
           "- **in-band union** = keep every reliable hit within 0.005 identity of the best (remove order stochasticity among co-optimal hits).",
           "- **all-reliable union** = keep every reliable hit down to the family floor.",
           ""]

    # ---- self-recovery (the headline: did we pick the source organism?) ----
    sr = {suf: self_recovery(suf) for _, suf in VARIANTS}
    out.append("## 1. Self-recovery — P(source selected | source in top-20), include scenario")
    out.append("\n| region | " + " | ".join(n for n, _ in VARIANTS) + " |")
    out.append("|" + "---|" * (len(VARIANTS) + 1))
    for r in REGION_ORDER + ["ALL"]:
        cells = " | ".join(pct(sr[suf][0].get(r)) for _, suf in VARIANTS)
        out.append(f"| {r} | {cells} |")
    out.append("\n_self-as-anchor | source in top-20 (include):_")
    out.append("\n| region | " + " | ".join(n for n, _ in VARIANTS) + " |")
    out.append("|" + "---|" * (len(VARIANTS) + 1))
    for r in REGION_ORDER:
        cells = " | ".join(pct(sr[suf][1].get(r)) for _, suf in VARIANTS)
        out.append(f"| {r} | {cells} |")

    # ---- taxonomy ----
    out.append("\n## 2. Taxonomic concordance (correct | called)")
    out.append(table("Anchor — Genus, include", lambda s: taxacc(s, "anchor", "include", "Genus")))
    out.append(table("Anchor — Species, include", lambda s: taxacc(s, "anchor", "include", "Species")))
    out.append(table("Consensus — Genus, include", lambda s: taxacc(s, "consensus", "include", "Genus")))
    out.append(table("Consensus — Genus, exclude_0.987 (honest)",
                     lambda s: taxacc(s, "consensus", "exclude_0.987", "Genus")))
    out.append(table("Anchor — Genus, exclude_0.987 (honest)",
                     lambda s: taxacc(s, "anchor", "exclude_0.987", "Genus")))

    # ---- gene capture ----
    out.append("\n## 3. Gene capture (macro, PGFam sets)")
    out.append(table("Recall — include (% gene capture)", lambda s: genecap(s, "include", "recall")))
    out.append(table("Precision — include", lambda s: genecap(s, "include", "precision")))
    out.append(table("F1 — include", lambda s: genecap(s, "include", "f1")))
    out.append(table("Jaccard — include", lambda s: genecap(s, "include", "jaccard")))
    out.append(table("Recall — exclude_0.987 (honest generalization)",
                     lambda s: genecap(s, "exclude_0.987", "recall")))
    out.append(table("Precision — exclude_0.987", lambda s: genecap(s, "exclude_0.987", "precision")))
    out.append(table("F1 — exclude_0.987", lambda s: genecap(s, "exclude_0.987", "f1")))
    out.append("\n_OK (scorable) cells per region, include:_")
    out.append(table("n OK cells — include", lambda s: genecap_n(s, "include"), REGION_ORDER))

    report = "\n".join(out) + "\n"
    rp = os.path.join(DATA, "variant_comparison.md")
    open(rp, "w").write(report)
    print(report)
    print(f"\n[compare] wrote {rp}")


if __name__ == "__main__":
    main()
