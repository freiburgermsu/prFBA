#!/usr/bin/env python
"""
assess_representatives.py — for the two ASV case studies (EmilyKin V4-V5, new_data/codif_all),
the FINAL representative reference organisms per ASV are the genomes selected by the
synthetic-genome pipeline run on the GPU-exhaustive Smith-Waterman alignment. This assesses
the edlib->biopython (CPU prefilter) alignment against that GPU ground truth on two questions:

  (1) COVERAGE  — does edlib's top-20 even CONTAIN each GPU-selected representative genome?
                  (if it never surfaces the genome, edlib could never select it.) Reported at
                  genome_id level (strict) and species level (lenient: same species anywhere
                  in edlib's top-20).
  (2) AGREEMENT — would running the IDENTICAL selection algorithm + gene provider on the edlib
                  top-20 have chosen the SAME representative genome(s)? Compares the edlib
                  selection audit against the GPU selection audit per ASV: identical genome set,
                  identical species set, same anchor (best representative), Jaccard, abstain.

Inputs per system: gpu selection audit, edlib selection audit, edlib top-20 hits.
Outputs (this dir): per-ASV CSV per system, a combined summary JSON, and a markdown report.

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations
import csv, json, os, statistics as st
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
EK = "/home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits"
ND = "/home/freiburger/Documents/codiffusion_bioreactor/bvbrc_alignment_hits"

SYSTEMS = {
    "EmilyKin": dict(
        gpu_sel=f"{EK}/asv_reference_selection_gpu.json",
        edlib_sel=f"{EK}/asv_reference_selection.json",       # generated_from edlib top20
        edlib_hits=f"{EK}/asv_top20_alignment_hits.json",
    ),
    "new_data": dict(
        gpu_sel=f"{ND}/asv_reference_selection_gpu.json",
        edlib_sel=f"{ND}/asv_reference_selection_edlib.json",
        edlib_hits=f"{ND}/asv_top20_alignment_hits.json",
    ),
}


def _load(p):
    d = json.load(open(p))
    d.pop("_meta", None)
    return d


def _norm_sp(s):
    """Normalize a species string for lenient matching; empty -> None (never matches)."""
    s = (s or "").strip().lower()
    return s or None


def reps(rec):
    """(ordered genome_ids, anchor genome_id, {gid: species}, set of normalized species)."""
    sel = rec.get("selected", [])
    gids = [str(s["genome_id"]) for s in sel]
    anchor = next((str(s["genome_id"]) for s in sel if s.get("role") == "anchor"), gids[0] if gids else None)
    sp_by_gid = {str(s["genome_id"]): s.get("species") for s in sel}
    species = {x for x in (_norm_sp(s.get("species")) for s in sel) if x}
    return gids, anchor, sp_by_gid, species


def edlib_topN(hits_rec):
    """(rank-by-genome_id, set of normalized species) over the edlib top-20 hits."""
    rank_by_gid, species = {}, set()
    for h in hits_rec.get("top20", []):
        gid = str(h.get("genome_id"))
        rank_by_gid.setdefault(gid, h.get("rank"))          # best (lowest) rank per genome
        sp = _norm_sp((h.get("lineage") or {}).get("Species"))
        if sp:
            species.add(sp)
    return rank_by_gid, species


def assess(name, paths):
    gpu, edl, hits = _load(paths["gpu_sel"]), _load(paths["edlib_sel"]), _load(paths["edlib_hits"])
    asvs = [a for a in gpu if a in edl]                      # both selections cover all ASVs
    rows = []
    # accumulators
    cov_gid_hit = cov_gid_tot = cov_sp_hit = cov_sp_tot = 0
    anchor_in = anchor_rank1 = 0
    anchor_rank_bins = Counter()
    both_built = ident_gid = ident_sp = same_anchor_gid = same_anchor_sp = 0
    jaccards = []
    g_abst = e_abst = abst_agree = only_gpu = only_edl = 0

    for a in asvs:
        g_gids, g_anc, g_spmap, g_sp = reps(gpu[a])
        e_gids, e_anc, e_spmap, e_sp = reps(edl[a])
        rank_by_gid, e_top_species = edlib_topN(hits.get(a, {}))
        g_built, e_built = len(g_gids) > 0, len(e_gids) > 0
        g_abst += (not g_built); e_abst += (not e_built)
        abst_agree += (g_built == e_built)
        only_gpu += (g_built and not e_built); only_edl += (e_built and not g_built)

        # --- (1) coverage of GPU reps in the edlib top-20 (only when GPU built) ---
        per_ranks = {}
        if g_built:
            for gid in g_gids:
                cov_gid_tot += 1
                r = rank_by_gid.get(gid)
                per_ranks[gid] = r
                if r is not None:
                    cov_gid_hit += 1
            for gid in g_gids:
                sp = _norm_sp(g_spmap.get(gid))
                if sp:
                    cov_sp_tot += 1
                    cov_sp_hit += (sp in e_top_species)
            ar = rank_by_gid.get(g_anc)
            anchor_in += (ar is not None)
            anchor_rank1 += (ar == 1)
            anchor_rank_bins["rank_1" if ar == 1 else
                             "rank_2_5" if (ar and ar <= 5) else
                             "rank_6_20" if ar else "absent"] += 1

        # --- (2) would edlib pick the same reps? (only when both built) ---
        gset, eset = set(g_gids), set(e_gids)
        genome_identical = species_identical = anchor_gid_match = anchor_sp_match = None
        jac = None
        if g_built and e_built:
            both_built += 1
            genome_identical = gset == eset
            species_identical = g_sp == e_sp
            anchor_gid_match = g_anc == e_anc
            anchor_sp_match = _norm_sp(g_spmap.get(g_anc)) == _norm_sp(e_spmap.get(e_anc)) \
                and _norm_sp(g_spmap.get(g_anc)) is not None
            jac = len(gset & eset) / len(gset | eset) if (gset | eset) else 1.0
            ident_gid += genome_identical
            ident_sp += species_identical
            same_anchor_gid += anchor_gid_match
            same_anchor_sp += anchor_sp_match
            jaccards.append(jac)

        rows.append(dict(
            asv=a, tier=gpu[a].get("tier"),
            gpu_n=len(g_gids), edlib_n=len(e_gids),
            gpu_reps="|".join(g_gids), edlib_reps="|".join(e_gids),
            gpu_anchor=g_anc, gpu_anchor_species=g_spmap.get(g_anc),
            edlib_anchor=e_anc, edlib_anchor_species=e_spmap.get(e_anc),
            gpu_anchor_edlib_rank=("absent" if g_built and rank_by_gid.get(g_anc) is None
                                   else rank_by_gid.get(g_anc)),
            gpu_reps_in_edlib_top20=sum(1 for r in per_ranks.values() if r is not None),
            gpu_reps_total=len(g_gids),
            gpu_rep_edlib_ranks="|".join(str(per_ranks.get(g, "absent")) for g in g_gids),
            genome_set_identical=genome_identical, species_set_identical=species_identical,
            anchor_genome_match=anchor_gid_match, anchor_species_match=anchor_sp_match,
            jaccard_genome=(round(jac, 3) if jac is not None else None),
        ))

    n = len(asvs)
    summary = dict(
        system=name, n_asvs=n,
        abstain=dict(gpu_abstain=g_abst, edlib_abstain=e_abst,
                     abstain_agreement_pct=round(100 * abst_agree / n, 2),
                     only_gpu_built=only_gpu, only_edlib_built=only_edl),
        coverage_of_gpu_reps_in_edlib_top20=dict(
            n_gpu_reps=cov_gid_tot,
            pct_reps_present_genome_level=round(100 * cov_gid_hit / max(cov_gid_tot, 1), 2),
            pct_reps_present_species_level=round(100 * cov_sp_hit / max(cov_sp_tot, 1), 2),
            pct_anchor_present=round(100 * anchor_in / max(g_built_total(rows), 1), 2),
            pct_anchor_at_edlib_rank1=round(100 * anchor_rank1 / max(g_built_total(rows), 1), 2),
            anchor_edlib_rank_distribution=dict(anchor_rank_bins),
        ),
        agreement_edlib_vs_gpu_selection=dict(
            n_both_built=both_built,
            pct_identical_genome_set=round(100 * ident_gid / max(both_built, 1), 2),
            pct_identical_species_set=round(100 * ident_sp / max(both_built, 1), 2),
            pct_same_anchor_genome=round(100 * same_anchor_gid / max(both_built, 1), 2),
            pct_same_anchor_species=round(100 * same_anchor_sp / max(both_built, 1), 2),
            mean_jaccard_genome=round(st.mean(jaccards), 3) if jaccards else None,
            median_jaccard_genome=round(st.median(jaccards), 3) if jaccards else None,
            n_asvs_differ_genome_set=both_built - ident_gid,
        ),
    )
    # write per-ASV CSV
    csv_path = os.path.join(HERE, f"representative_coverage_{name}.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    return summary, rows, csv_path


def g_built_total(rows):
    return sum(1 for r in rows if r["gpu_n"] > 0)


def md_report(summaries):
    L = ["# edlib vs GPU — final representative-organism assessment\n",
         "For each ASV the **final representative reference organism(s)** are the genome(s) the "
         "synthetic-genome pipeline selected from the **GPU-exhaustive** Smith-Waterman top-20. "
         "This compares the **edlib->biopython** alignment against that ground truth.\n",
         "Both runs use the *identical* selection algorithm + BV-BRC PGFam gene provider; only the "
         "upstream alignment (GPU exhaustive vs edlib top-500 prefilter -> biopython) differs.\n"]
    for s in summaries:
        c = s["coverage_of_gpu_reps_in_edlib_top20"]; a = s["agreement_edlib_vs_gpu_selection"]; ab = s["abstain"]
        L += [f"\n## {s['system']}  ({s['n_asvs']} ASVs)\n",
              "**(1) Coverage — are the GPU-chosen representatives even in edlib's top-20?**\n",
              f"- GPU representative genomes total: **{c['n_gpu_reps']}**",
              f"- present in edlib top-20 (genome_id): **{c['pct_reps_present_genome_level']}%**",
              f"- present in edlib top-20 (species-level): **{c['pct_reps_present_species_level']}%**",
              f"- anchor (best rep) present in edlib top-20: **{c['pct_anchor_present']}%**; "
              f"at edlib **rank 1: {c['pct_anchor_at_edlib_rank1']}%**",
              f"- anchor edlib-rank distribution: `{c['anchor_edlib_rank_distribution']}`\n",
              "**(2) Agreement — would edlib have selected the same representatives?**\n",
              f"- ASVs built by both: **{a['n_both_built']}**",
              f"- **identical representative genome set: {a['pct_identical_genome_set']}%** "
              f"({a['n_asvs_differ_genome_set']} ASVs differ)",
              f"- identical species set: **{a['pct_identical_species_set']}%**",
              f"- same anchor genome: **{a['pct_same_anchor_genome']}%**; "
              f"same anchor species: **{a['pct_same_anchor_species']}%**",
              f"- genome-set Jaccard: mean **{a['mean_jaccard_genome']}**, median **{a['median_jaccard_genome']}**",
              f"- abstain agreement: **{ab['abstain_agreement_pct']}%** "
              f"(GPU-only built {ab['only_gpu_built']}, edlib-only built {ab['only_edlib_built']})\n"]
    L += ["\n## Bottom line\n",
          "See per-system numbers above and the per-ASV CSVs (`representative_coverage_<system>.csv`) "
          "for the exact ASVs where edlib would have diverged.\n"]
    open(os.path.join(HERE, "REPORT_edlib_vs_gpu_representatives.md"), "w").write("\n".join(L))


def main():
    summaries = []
    for name, paths in SYSTEMS.items():
        s, rows, csv_path = assess(name, paths)
        summaries.append(s)
        print(f"[{name}] wrote {csv_path}")
        print(json.dumps(s, indent=2))
    json.dump(summaries, open(os.path.join(HERE, "assessment_summary.json"), "w"), indent=2)
    md_report(summaries)
    print(f"\n[done] -> {HERE}/assessment_summary.json + REPORT_edlib_vs_gpu_representatives.md")


if __name__ == "__main__":
    main()
