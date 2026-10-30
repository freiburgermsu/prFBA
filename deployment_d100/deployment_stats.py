#!/usr/bin/env python
"""
deployment_stats.py -- selection-size and cross-phylum census of a real-amplicon prFBA
deployment run (selection JSON [+ hits JSON] -> stats JSON).

It is the counting method behind the manuscript's deployment figures (prFBA.tex
sec:recon / Discussion / Limitations; deck DEPLOY_N1, XPHYLUM_DEPLOY). Run on the
20-hit EmilyKin build selection it reproduces them exactly:

    3,950 ASVs | 3,753 answered | 197 abstain | n_genomes=1: 2,344 (62.5% of answered)
    | 2: 599 | >=3: 810 | mean 2.034 (over ANSWERED ASVs; the selection _meta's
    mean_handful 1.933 is over all 3,950) | cross-phylum 143/451 = 31.7% (n=2) and
    462/783 = 59.0% (n>=3), among the 1,234 of 1,409 multi-genome merges with >=2
    members resolved to an NCBI lineage | 18 merges unite Bacteria with Archaea.

COUNTING RULES (each mirrors the computation that produced the published numbers;
see README.md in this directory):
  * answered = n_selected > 0; abstained = n_selected == 0.
  * n_genomes = len(selected) (= the build provenance's `genomes` list; checked by
    --provenance-dir). Equivalent (species-deduplicated) genomes are NOT counted.
  * mean n_genomes is over answered ASVs.
  * cross-phylum: for each ASV with n_genomes >= 2, map every selected genome_id to an
    NCBI taxon id, resolve its lineage with taxopy on prFBA's nodes.dmp/names.dmp
    (no merged.dmp), keep the members that resolve; a merge with < 2 resolved members
    is "unresolved" and leaves the denominator. The merge is bucketed by its FULL
    member count (n=2 vs n>=3, resolved or not). It spans phyla when its resolved
    members carry > 1 distinct NON-None phylum (an unassigned phylum is not a
    distinct value). The None-counting variant (len(set) > 1 with None allowed) is
    reported too, as `cross_phylum_none_counted` (40.6% / 66.2% on the 20-hit run).
    bacteria_archaea = merges whose resolved members include both domains
    (domain = rank 'domain' or else 'superkingdom').
  * --taxon-source index (DEFAULT; the published method): genome_id -> taxon_id from
    prFBA/index.parquet (rows with NaN taxon_id dropped; last row wins). index.parquet
    holds one representative genome per unique 16S sequence, so a selected genome that
    is not some sequence's representative stays unresolved -- this is what makes the
    denominator 1,234 of 1,409. --taxon-source selection uses the taxon_id the
    selection records for each selected genome (from the hit headers) instead; it
    resolves more merges and is NOT the published method.

    python deployment_d100/deployment_stats.py \
        --selection /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_reference_selection_gpu.json \
        --hits      /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_top20_alignment_hits_gpu.json \
        --out       deployment_d100/stats_emilykin_d20.json \
        [--provenance-dir /home/freiburger/Documents/EmilyKin/synthetic_genomes/provenance]

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import collections
import glob
import json
import os
import sys
import time

PRFBA = "/home/freiburger/Documents/prFBA"
INDEX_PARQUET = os.path.join(PRFBA, "index.parquet")
TAXDUMP_NODES = os.path.join(PRFBA, "nodes.dmp")
TAXDUMP_NAMES = os.path.join(PRFBA, "names.dmp")
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def _records(d):
    return {k: v for k, v in d.items() if not k.startswith("_")}


def _pct(a, b, nd=1):
    return round(100.0 * a / b, nd) if b else None


# ------------------------------------------------------------------ lineage resolution
class Lineages:
    def __init__(self, taxon_source, selection_recs):
        import taxopy

        self.taxdb = taxopy.TaxDb(nodes_dmp=TAXDUMP_NODES, names_dmp=TAXDUMP_NAMES)
        self.taxopy = taxopy
        self.cache = {}
        self.source = taxon_source
        if taxon_source == "index":
            import pandas as pd

            df = pd.read_parquet(INDEX_PARQUET, columns=["genome_id", "taxon_id"]).dropna()
            self.g2t = {str(a): int(b) for a, b in zip(df.genome_id, df.taxon_id)}
        elif taxon_source == "selection":
            self.g2t = {}
            for r in selection_recs.values():
                for s in r["selected"]:
                    if s.get("taxon_id") is not None:
                        self.g2t[str(s["genome_id"])] = int(s["taxon_id"])
        else:
            raise ValueError(taxon_source)

    def lineage(self, gid):
        t = self.g2t.get(str(gid))
        if t is None:
            return None
        if t in self.cache:
            return self.cache[t]
        try:
            rd = self.taxopy.Taxon(t, self.taxdb).rank_name_dictionary
        except Exception:
            self.cache[t] = None
            return None
        out = {r: ((rd.get("domain") or rd.get("superkingdom")) if r == "domain" else rd.get(r))
               for r in RANKS}
        self.cache[t] = out
        return out


def cross_phylum(merges, lin):
    """merges: {asv: [genome_id, ...]} restricted to n_genomes >= 2."""
    tie = collections.Counter()        # resolved merges per bucket (the denominators)
    cross = collections.Counter()      # >1 distinct non-None phylum (published)
    cross_none = collections.Counter() # >1 distinct phylum, None counted as a value
    none_in = collections.Counter()
    total = collections.Counter()
    unresolved = collections.Counter()
    partial = 0
    bac_arch = 0
    dom_none = 0
    for asv, gs in merges.items():
        b = "n=2" if len(gs) == 2 else "n>=3"
        total[b] += 1
        L = [y for y in (lin.lineage(g) for g in gs) if y]
        if len(L) < 2:
            unresolved[b] += 1
            continue
        if len(L) < len(gs):
            partial += 1
        tie[b] += 1
        ph = {y.get("phylum") for y in L}
        if len(ph) > 1:
            cross_none[b] += 1
        if None in ph:
            none_in[b] += 1
        if len(ph - {None}) > 1:
            cross[b] += 1
        doms = {y.get("domain") for y in L}
        if None in doms:
            dom_none += 1
        if {"Bacteria", "Archaea"} <= doms:
            bac_arch += 1
    res = {}
    for b in ("n=2", "n>=3"):
        res[b] = dict(multi_genome_merges=total[b], resolved=tie[b], unresolved=unresolved[b],
                      cross_phylum=cross[b], cross_phylum_pct=_pct(cross[b], tie[b]),
                      cross_phylum_none_counted=cross_none[b],
                      cross_phylum_none_counted_pct=_pct(cross_none[b], tie[b]),
                      merges_with_unassigned_phylum=none_in[b])
    res["all_multi"] = dict(multi_genome_merges=sum(total.values()), resolved=sum(tie.values()),
                            unresolved=sum(unresolved.values()),
                            partially_resolved=partial,
                            cross_phylum=sum(cross.values()),
                            cross_phylum_pct=_pct(sum(cross.values()), sum(tie.values())),
                            bacteria_archaea=bac_arch, any_member_domain_none=dom_none)
    return res


# ------------------------------------------------------------------ selection rule audit
def rule_audit(recs):
    """What selection rule produced this audit? union (exact-tie) vs legacy reducer, and
    whether the P3 gene-richest-conspecific tie-break was on (inferred: with P3 off the
    species representative is always the lexicographically smallest genome_id among
    co-optimal conspecifics -- same identity AND align_score)."""
    answered = [r for r in recs.values() if r["n_selected"] > 0]
    def _union_flag(r):   # select_references appends "select_all:window=<w>" in union mode
        return next((f for f in r.get("flags", []) if f.split(":")[0] == "select_all"), None)

    n_union = sum(1 for r in answered if _union_flag(r))
    union_windows = collections.Counter(_union_flag(r) for r in answered if _union_flag(r))
    stages = collections.Counter()
    for r in recs.values():
        for x in r.get("decisions", []):
            if x.get("disposition") == "excluded":
                stages[x.get("stage")] += 1
    co = smallest = 0
    windows = []
    for r in answered:
        dec = {str(x["genome_id"]): x for x in r.get("decisions", [])}
        if _union_flag(r):
            ids = [s["identity"] for s in r["selected"]]
            windows.append(round(max(ids) - min(ids), 6))
        for s in r["selected"]:
            eq = [e for e in (s.get("equivalent_genomes") or [])
                  if abs(e["identity"] - s["identity"]) < 1e-9
                  and dec.get(str(e["genome_id"]), {}).get("align_score") == s["align_score"]]
            if not eq:
                continue
            co += 1
            smallest += all(str(s["genome_id"]) < str(e["genome_id"]) for e in eq)
    legacy_stages = {"marginal_gain", "family_consensus", "tier_cap", "contamination_budget"}
    rule = ("exact-tie union" if answered and n_union == len(answered) and
            set(union_windows) == {"select_all:window=0.0"} and max(windows or [0.0]) == 0.0 else
            "legacy ordered reducer" if n_union == 0 else "mixed/other")
    return dict(
        rule=rule, answered=len(answered), answered_with_select_all_flag=n_union,
        select_all_flags=dict(union_windows),
        max_selected_identity_spread=max(windows or [0.0]),
        excluded_by_stage=dict(stages.most_common()),
        legacy_stage_exclusions=sum(stages[s] for s in legacy_stages),
        p3_cooptimal_conspecific_cases=co,
        p3_rep_is_lexicographic_smallest=smallest,
        p3_inferred=("off" if co and smallest == co else "on" if co else "undetermined"),
    )


def hits_audit(hits, recs, list_key="top20"):
    lens = [len(h.get(list_key) or []) for h in hits.values()]
    ranks_ok = all([x["rank"] for x in (h.get(list_key) or [])] ==
                   list(range(1, len(h.get(list_key) or []) + 1)) for h in hits.values())
    return dict(n_asvs=len(hits), list_key=list_key,
                hits_per_record_max=max(lens) if lens else 0,
                hits_per_record_min=min(lens) if lens else 0,
                hits_per_record_mean=round(sum(lens) / max(len(lens), 1), 3),
                records_with_zero_hits=sum(1 for x in lens if x == 0),
                ranks_contiguous_from_1=ranks_ok,
                same_asv_set_as_selection=set(hits) == set(recs))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selection", required=True, help="asv_reference_selection*.json audit")
    ap.add_argument("--hits", default=None, help="the asv_top20_alignment_hits*.json it was made from")
    ap.add_argument("--out", required=True, help="stats JSON to write (refuses to overwrite)")
    ap.add_argument("--taxon-source", default="index", choices=["index", "selection"])
    ap.add_argument("--provenance-dir", default=None,
                    help="build provenance/ dir: verify selected genomes == provenance genomes")
    ap.add_argument("--force", action="store_true", help="allow overwriting --out")
    args = ap.parse_args()
    if os.path.exists(args.out) and not args.force:
        sys.exit(f"[abort] {args.out} exists (pass --force to overwrite)")

    t0 = time.perf_counter()
    sel = json.load(open(args.selection))
    meta = sel.get("_meta", {})
    recs = _records(sel)
    n = len(recs)
    sizes = collections.Counter(r["n_selected"] for r in recs.values())
    answered = n - sizes.get(0, 0)
    tot_sel = sum(k * v for k, v in sizes.items())
    abst_flags = collections.Counter(f for r in recs.values() if r["n_selected"] == 0
                                     for f in r.get("flags", []) if f.startswith("abstain"))
    out = dict(
        selection=os.path.abspath(args.selection),
        selection_generated_from=meta.get("generated_from"),
        selection_gene_data=meta.get("gene_data"),
        hits=os.path.abspath(args.hits) if args.hits else None,
        taxon_source=args.taxon_source,
        n_asvs=n, answered=answered, abstained=sizes.get(0, 0),
        abstained_by_flag=dict(abst_flags),
        n_genomes_distribution={str(k): sizes[k] for k in sorted(sizes)},
        n_genomes_1=sizes.get(1, 0), n_genomes_2=sizes.get(2, 0),
        n_genomes_3plus=sum(v for k, v in sizes.items() if k >= 3),
        n_genomes_1_pct_of_answered=_pct(sizes.get(1, 0), answered),
        n_genomes_2_pct_of_answered=_pct(sizes.get(2, 0), answered),
        n_genomes_3plus_pct_of_answered=_pct(sum(v for k, v in sizes.items() if k >= 3), answered),
        mean_n_genomes_answered=round(tot_sel / max(answered, 1), 4),
        mean_n_genomes_all=round(tot_sel / max(n, 1), 4),
        tiers=dict(collections.Counter(r["tier"] for r in recs.values())),
        selection_rule=rule_audit(recs),
    )
    if args.hits:
        hits = json.load(open(args.hits))
        out["hits_audit"] = hits_audit(hits, recs)
        del hits
    merges = {a: [str(s["genome_id"]) for s in r["selected"]]
              for a, r in recs.items() if r["n_selected"] >= 2}
    if args.provenance_dir:
        prov = {}
        for p in glob.glob(os.path.join(args.provenance_dir, "*.json")):
            j = json.load(open(p))
            prov[j.get("asv") or os.path.basename(p)[:-5]] = [str(g) for g in (j.get("genomes") or [])]
        answered_ids = {a for a, r in recs.items() if r["n_selected"] > 0}
        mism = [a for a in answered_ids
                if a not in prov or sorted(prov[a]) != sorted(str(s["genome_id"]) for s in recs[a]["selected"])]
        out["provenance_check"] = dict(provenance_dir=os.path.abspath(args.provenance_dir),
                                       n_provenance=len(prov), answered=len(answered_ids),
                                       mismatched_or_missing=len(mism),
                                       extra_in_provenance=len(set(prov) - answered_ids),
                                       examples=sorted(mism)[:10])
    lin = Lineages(args.taxon_source, recs)
    out["cross_phylum"] = cross_phylum(merges, lin)
    out["wall_seconds"] = round(time.perf_counter() - t0, 1)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    cp = out["cross_phylum"]
    print(f"[stats] {os.path.basename(args.selection)}: n={n} answered={answered} "
          f"abstain={out['abstained']} n1={out['n_genomes_1']} ({out['n_genomes_1_pct_of_answered']}%) "
          f"n2={out['n_genomes_2']} n3+={out['n_genomes_3plus']} mean={out['mean_n_genomes_answered']} "
          f"| xphylum n=2 {cp['n=2']['cross_phylum']}/{cp['n=2']['resolved']} ({cp['n=2']['cross_phylum_pct']}%) "
          f"n>=3 {cp['n>=3']['cross_phylum']}/{cp['n>=3']['resolved']} ({cp['n>=3']['cross_phylum_pct']}%) "
          f"resolved {cp['all_multi']['resolved']}/{cp['all_multi']['multi_genome_merges']} "
          f"bac+arch={cp['all_multi']['bacteria_archaea']} | rule={out['selection_rule']['rule']} "
          f"P3={out['selection_rule']['p3_inferred']} -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
