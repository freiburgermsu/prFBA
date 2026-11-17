#!/usr/bin/env python
"""Path A disclosure counts (prFBA_paper sec:sharedref, the paragraph after "the primary analysis"):
of the 161 panel genomes the 2026-10-13 eligibility pass dropped as `empty_ec`, how many have a
complete EC on the repaired gene-family data, how many of those also pass the two remaining free
filters (16S length 1,200-1,800 bp, ACGTN only), their domains, and the number a re-plan would be
expected to admit after the last filter (HMM match-column occupancy), which was not run.

Inputs (all read-only):
  region_validation/shared_reference/data/panel_eligible.json   the kept panel and its drop reasons
  region_validation/shared_reference/data/_cache_panel_ec.json  panel EC sets, regenerated from the
        repaired genome_function_sets.json in re-score Run 3 (the eligibility pass's own copy of
        2026-10-13 is $W/pre_repair/... when it exists, used here as a control)
  region_validation/shared_reference/data/_cache_panel_md5_seq.json  panel best-16S sequences
  region_validation/data/selection_10k.json                    best_16s_len, domain, genus
The free filters are plan_folds.py phase_panel's own (FROZEN reference_qc: ec_min 1, len_min 1200,
len_max 1800, ACGTN). The expected count applies the panel's measured occupancy drop rate,
occupancy drops / genomes reaching the occupancy filter = 535 / (8,872 + 535).

Writes manuscript_numbers/rescore_sources/path_a_counts.json.
"""
import collections, json, os

R = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
W = os.environ.get("PRFBA_RESCORE_W", os.path.expanduser("~/Documents/prFBA_rescore_2026-11"))
SR = os.path.join(R, "region_validation", "shared_reference", "data")


def load(p):
    with open(p) as fh:
        return json.load(fh)


def main():
    pe = load(os.path.join(SR, "panel_eligible.json"))
    qc = pe["meta"]["frozen"]["reference_qc"]
    ec_now = load(os.path.join(SR, "_cache_panel_ec.json"))
    seq = load(os.path.join(SR, "_cache_panel_md5_seq.json"))
    sel = {g["genome_id"]: g for g in load(os.path.join(R, "region_validation", "data", "selection_10k.json"))["genomes"]}
    empty = sorted(g for g, why in pe["dropped"].items() if why == "empty_ec")
    eligible_genera = {v.get("genus") for v in pe["genomes"].values() if v.get("genus")}

    def n_ec(gid):
        rec = ec_now.get(gid)
        return None if rec is None else len(rec["ec"])

    now_ec = [g for g in empty if (n_ec(g) or 0) >= qc["ec_min"]]
    passes, fails = [], collections.Counter()
    for g in now_ec:
        s = sel[g]
        if not (qc["len_min"] <= s["best_16s_len"] <= qc["len_max"]):
            fails["length"] += 1
            continue
        q = seq.get(s["best_16s_md5"])
        if q is None:
            fails["no_sequence"] += 1
            continue
        if set(q) - set("ACGTN"):
            fails["bad_chars"] += 1
            continue
        passes.append(g)
    dom = collections.Counter(sel[g]["domain"] for g in passes)
    occ_drop, n_elig = pe["meta"]["drop_counts"]["occupancy"], pe["meta"]["n_eligible"]
    rate = occ_drop / (n_elig + occ_drop)
    out = {
        "what": "Path A disclosure counts for the shared-reference panel (prFBA_paper sec:sharedref)",
        "generator": "manuscript_numbers/rescore_sources/path_a_counts.py",
        "panel_eligibility_generated": pe["meta"]["generated"],
        "n_dropped_empty_ec": len(empty),
        "n_empty_ec_now_with_ec": len(now_ec),
        "n_empty_ec_still_empty": len(empty) - len(now_ec),
        "n_now_with_ec_passing_free_filters": len(passes),
        "free_filter_failures": dict(fails),
        "passing_by_domain": dict(dom),
        "passing_sharing_a_genus_with_an_eligible_genome": sum(1 for g in passes if sel[g].get("genus") in eligible_genera),
        "occupancy_drop_rate": round(rate, 6),
        "occupancy_drop_rate_counts": [occ_drop, n_elig + occ_drop],
        "expected_new_eligible_after_occupancy": round(len(passes) * (1 - rate), 2),
        "ids_now_with_ec_passing_free_filters": passes,
    }
    pre = os.path.join(W, "pre_repair", "region_validation", "shared_reference", "data", "_cache_panel_ec.json")
    if os.path.exists(pre):
        old = load(pre)
        out["control_pre_repair_cache_empty_ec"] = sum(1 for g in empty if len((old.get(g) or {"ec": []})["ec"]) < qc["ec_min"])
    path = os.path.join(R, "manuscript_numbers", "rescore_sources", "path_a_counts.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: v for k, v in out.items() if not k.startswith("ids_")}, indent=1))


if __name__ == "__main__":
    main()
