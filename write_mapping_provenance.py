#!/usr/bin/env python
"""
write_mapping_provenance.py — run-level provenance + abstention report for one
16S ASV -> reference-genome -> synthetic-genome mapping run.

Reads the alignment hits, the exact selection audit and the synthetic-genome build
manifest, and writes into --genome-dir:

  mapping_provenance.json   pipeline commit + script hashes, every input (md5/size/mtime),
                            per-stage method + parameters (alignment scoring, selection knobs,
                            tier thresholds, abstain rules), environment, and result summaries
  asv_mapping_summary.csv   one row per ASV (built AND abstained): tier, flags, anchor genome,
                            every selected source genome, KBase-genome stats, abstain reason
  abstained_asvs.csv        every abstained ASV: machine flags + the per-hit exclusion
  abstained_asvs.md         stage/reason evidence + a plain-English explanation

    ~/Documents/py_venv/bin/python write_mapping_provenance.py \
        --hits       .../bvbrc_alignment_hits/asv_top20_alignment_hits_gpu.json \
        --selection  .../bvbrc_alignment_hits/asv_reference_selection_gpu.54d1150.json \
        --genome-dir .../synthetic_genomes_54d1150 \
        --asv-fasta .../new_data/dna-sequences_codif_all.fasta \
        --taxonomy-csv .../new_data/taxonomy_codif_all.csv \
        [--alignment-stats .../gpu_exhaustive_full_stats.json] [--build-log LOG] [--prefetch-gaps JSON]
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import select_references as S  # noqa: E402

PIPELINE_SCRIPTS = [
    "gpu_align.py", "align_hits.py", "edlib_biopython_hits.py", "prefetch_gene_families.py",
    "select_references.py", "build_synthetic_genomes.py", "build_synthetic_genomes_parallel.py",
    "write_mapping_provenance.py", "upload_kbase_genomes.py",
]


# --------------------------------------------------------------------------- helpers
def _sh(cmd, cwd=HERE):
    try:
        return subprocess.check_output(cmd, cwd=cwd, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return None


def _hash(path, algo="md5", chunk=1 << 24):
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def _file_record(path, hash_it=True, hash_limit_mb=2048):
    if not path or not os.path.exists(path):
        return {"path": path, "exists": False}
    st = os.stat(path)
    rec = {"path": os.path.abspath(path), "exists": True, "size_bytes": st.st_size,
           "mtime": dt.datetime.fromtimestamp(st.st_mtime).astimezone().isoformat(timespec="seconds")}
    if hash_it and st.st_size <= hash_limit_mb * (1 << 20):
        rec["md5"] = _hash(path)
    elif hash_it:
        rec["md5"] = f"(skipped: > {hash_limit_mb} MB)"
    return rec


def _git_info():
    info = {"repo": HERE,
            "commit": _sh(["git", "rev-parse", "HEAD"]),
            "commit_short": _sh(["git", "rev-parse", "--short", "HEAD"]),
            "commit_date": _sh(["git", "log", "-1", "--format=%cI"]),
            "commit_subject": _sh(["git", "log", "-1", "--format=%s"]),
            "branch": _sh(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
            "remote": _sh(["git", "remote", "get-url", "origin"])}
    dirty = _sh(["git", "status", "--porcelain", "--", "*.py", "*.sh"]) or ""
    info["uncommitted_pipeline_changes"] = [l for l in dirty.splitlines() if l.strip()]
    info["clean"] = not info["uncommitted_pipeline_changes"]
    return info


def _script_hashes():
    out = {}
    for s in PIPELINE_SCRIPTS:
        p = os.path.join(HERE, s)
        if os.path.exists(p):
            out[s] = {"sha256": _hash(p, "sha256"),
                      "last_commit": _sh(["git", "log", "-1", "--format=%h %cI", "--", s])}
    return out


def _alignment_code_status(hits_path):
    """Document whether gpu_align.py's production (--full / run_full) path changed between
    the commit current when the hits were produced and HEAD."""
    mtime = dt.datetime.fromtimestamp(os.stat(hits_path).st_mtime).astimezone().isoformat()
    then = _sh(["git", "log", "-1", "--format=%h", f"--before={mtime}", "--", "gpu_align.py"])
    head = _sh(["git", "rev-parse", "--short", "HEAD"])
    rec = {"hits_mtime": mtime, "gpu_align_commit_at_hits_time": then, "head": head}
    if not then:
        return rec
    diff = _sh(["git", "diff", then, "HEAD", "--", "gpu_align.py"]) or ""
    hunks = [l for l in diff.splitlines() if l.startswith("@@")]
    rec["diff_hunks_since"] = hunks
    rec["run_full_touched"] = any("run_full" in h for h in hunks)
    # the changed lines: do any of them lie inside run_full()?
    src = open(os.path.join(HERE, "gpu_align.py")).read().splitlines()
    rf = next((i for i, l in enumerate(src, 1) if l.startswith("def run_full(")), None)
    rec["run_full_starts_at_line"] = rf
    changed_lines = []
    for h in hunks:
        m = re.match(r"@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", h)
        if m:
            a = int(m.group(1)); n = int(m.group(2) or 1)
            changed_lines.append([a, a + n - 1])
    rec["changed_line_ranges_in_head"] = changed_lines
    rec["changes_inside_run_full"] = bool(rf and any(a >= rf for a, _ in changed_lines))
    rec["verdict"] = ("gpu_align.py --full path (run_full) is byte-identical to the version that produced "
                      "the hits; re-alignment would reproduce the same file"
                      if not rec["changes_inside_run_full"] and not rec["run_full_touched"]
                      else "run_full changed since the hits were produced: RE-ALIGN before trusting these hits")
    return rec


def _sw_constants():
    src = open(os.path.join(HERE, "gpu_align.py")).read()
    out = {}
    for name in ("MATCH", "MISMATCH", "GAP_OPEN", "GAP_EXT", "MAXQ"):
        m = re.search(rf"^{name}\s*=\s*(-?\d+)", src, re.M)
        if m:
            out[name.lower()] = int(m.group(1))
    return out


def _environment():
    import importlib.metadata as md
    pk = {}
    for p in ["edlib", "biopython", "cupy-cuda13x", "cupy", "torch", "transformers", "numpy",
              "pandas", "requests", "taxopy", "modelseedpy", "cobrakbase"]:
        try:
            pk[p] = md.version(p)
        except Exception:
            pass
    gpu = _sh(["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"])
    return {"python": sys.version.split()[0], "interpreter": sys.executable,
            "platform": platform.platform(), "hostname": platform.node(),
            "cpu_count": os.cpu_count(), "gpu": gpu, "packages": pk}


def _parse_build_log(path):
    if not path or not os.path.exists(path):
        return None
    txt = open(path, errors="replace").read()
    keep = [l for l in txt.splitlines()
            if l.startswith(("[prefetch]", "  coverage:", "[select]", "[warm]", "[build]", "===", "PREFETCH", "NEWDATA"))]
    return keep[-40:]


# --------------------------------------------------------------------------- abstain reasoning
def _reason_category(stage, reason):
    if stage == "tier_floor":
        return "below_family_threshold"
    if stage == "reliability_gate" and reason and reason.startswith("coverage"):
        return "low_coverage"
    if stage == "reliability_gate" and reason and "SW-edlib" in reason:
        return "sw_edlib_discordant"
    if stage == "identity_floor":
        return "below_identity_floor"
    return stage or "unknown"


CATEGORY_TEXT = {
    "below_family_threshold": "best 16S identity is below the family threshold (t_family) — no reference "
                              "genome is close enough for a genome-level model; a model would inject false metabolism",
    "low_coverage": "alignment covers < cov_min of the ASV (partial / local alignment, not a bona fide 16S match)",
    "sw_edlib_discordant": "Smith-Waterman and edlib identities disagree by > swed_max (method-discordant — "
                           "possible chimera or low-complexity artefact)",
    "below_identity_floor": "hit identity is below the per-hit family floor",
}


def _explain(asv, rec, knobs):
    bi = rec.get("best_identity") or 0.0
    flags = rec.get("flags") or []
    decisions = rec.get("decisions") or []
    stages = Counter(_reason_category(d.get("stage"), d.get("reason")) for d in decisions)
    top = next((d for d in decisions if d.get("rank") == 1), decisions[0] if decisions else {})
    if "abstain_below_family" in flags or rec.get("tier") == "below":
        why = (f"Best 16S identity to any BV-BRC reference is {bi:.2%} (top hit: {top.get('organism')}, "
               f"genome {top.get('genome_id')}), below the family threshold of {knobs['t_family']:.1%}. "
               f"No reference genome is close enough for a genome-level model, so the pipeline abstains "
               f"rather than inject another organism's metabolism.")
        cat = "below_family_threshold"
    elif "abstain_reliability" in flags:
        parts = ", ".join(f"{n} hit(s) {CATEGORY_TEXT.get(c, c)}" for c, n in stages.most_common())
        why = (f"Best identity {bi:.2%} passes the family floor ({rec.get('tier')} tier), but every one of the "
               f"{len(decisions)} top-20 hits failed the reliability gate: {parts}. Top hit "
               f"{top.get('organism')} ({top.get('genome_id')}): {top.get('reason')}.")
        cat = stages.most_common(1)[0][0] if stages else "reliability"
    else:
        parts = ", ".join(f"{n}x {c}" for c, n in stages.most_common())
        why = f"No reference survived selection (flags={flags}); exclusion stages: {parts}."
        cat = stages.most_common(1)[0][0] if stages else "unknown"
    return cat, why, stages, top


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--asv-fasta", default=None)
    ap.add_argument("--taxonomy-csv", default=None)
    ap.add_argument("--alignment-stats", default=None)
    ap.add_argument("--build-log", default=None)
    ap.add_argument("--gene-cache", default=os.path.join(HERE, "bvbrc_cache", "genome_gene_families.json"))
    ap.add_argument("--genome-cache-dir", default=os.path.join(HERE, "bvbrc_cache", "kbase_genome_cache"))
    ap.add_argument("--prefetch-gaps", default=None, help="JSON with the candidate genome_ids whose PGFam "
                    "gene sets could not be fetched (and their classification)")
    ap.add_argument("--build-command", default=None, help="exact command line used for the build")
    args = ap.parse_args()

    gdir = os.path.abspath(args.genome_dir)
    man_path = os.path.join(gdir, "synthetic_genome_manifest.json")
    manifest = json.load(open(man_path))
    sel = json.load(open(args.selection))
    hits = json.load(open(args.hits))
    sel_meta = sel.pop("_meta", {})
    man_meta = manifest.get("_meta", {})
    knobs = dict(S.DEFAULT_KNOBS)

    asvs = [a for a in sel if not a.startswith("_")]
    flag_counts = Counter(f for a in asvs for f in (sel[a].get("flags") or []))
    tier_counts = Counter(sel[a].get("tier") for a in asvs)
    handful = Counter(sel[a].get("n_selected", 0) for a in asvs)

    # ---- abstained ASVs ----
    abst_rows = []
    for a in asvs:
        rec = sel[a]
        if rec.get("n_selected", 0) > 0:
            continue
        cat, why, stages, top = _explain(a, rec, knobs)
        abst_rows.append({
            "asv": a, "rel_ab": rec.get("rel_ab"), "asv_len": rec.get("asv_len"),
            "best_identity": rec.get("best_identity"), "tier": rec.get("tier"),
            "flags": ";".join(rec.get("flags") or []), "reason_category": cat,
            "n_hits_evaluated": len(rec.get("decisions") or []),
            "exclusion_stage_counts": ";".join(f"{c}={n}" for c, n in stages.most_common()),
            "top_hit_genome_id": top.get("genome_id"), "top_hit_organism": top.get("organism"),
            "top_hit_identity": top.get("identity"), "top_hit_coverage": top.get("coverage"),
            "top_hit_edlib_agreement": top.get("edlib_agreement"), "top_hit_family": top.get("family"),
            "top_hit_exclusion_reason": top.get("reason"),
            "midas_taxonomy": rec.get("midas_taxonomy"), "explanation": why,
        })
    abst_cat = Counter(r["reason_category"] for r in abst_rows)

    # ---- per-ASV summary (all) ----
    sum_rows = []
    for a in asvs:
        rec = sel[a]; m = manifest.get(a, {})
        selected = rec.get("selected") or []
        anchor = selected[0] if selected else {}
        row = {
            "asv": a, "status": m.get("status", "abstain" if not selected else "?"),
            "rel_ab": rec.get("rel_ab"), "asv_len": rec.get("asv_len"),
            "best_identity": rec.get("best_identity"), "tier": rec.get("tier"),
            "flags": ";".join(rec.get("flags") or []), "consensus_family": rec.get("consensus_family"),
            "n_selected": rec.get("n_selected", 0),
            "anchor_genome_id": anchor.get("genome_id"), "anchor_organism": anchor.get("organism"),
            "anchor_identity": anchor.get("identity"), "anchor_taxon_id": anchor.get("taxon_id"),
            "selected_genome_ids": "|".join(str(s.get("genome_id")) for s in selected),
            "selected_roles": "|".join(str(s.get("role")) for s in selected),
            "selected_confidence": "|".join(str(s.get("confidence")) for s in selected),
            "n_gene_set_empty": sum(1 for s in selected if s.get("gene_set_empty")),
            "n_equivalent_genomes": (rec.get("provenance") or {}).get("n_equivalent_genomes"),
            "est_union_genes": (rec.get("provenance") or {}).get("est_union_genes"),
            "scientific_name": m.get("scientific_name"), "taxonomy": m.get("taxonomy"),
            "n_features": m.get("n_features"), "n_source_genomes_built": m.get("n_genomes"),
            "build_failed_sources": json.dumps(m.get("failed")) if m.get("failed") else "",
            "genome_file": m.get("genome_file"), "provenance_file": m.get("provenance_file"),
            "midas_taxonomy": rec.get("midas_taxonomy"),
            "abstain_reason": next((r["explanation"] for r in abst_rows if r["asv"] == a), ""),
        }
        sum_rows.append(row)

    # ---- inputs ----
    inputs = {
        "asv_fasta": _file_record(args.asv_fasta),
        "taxonomy_csv": _file_record(args.taxonomy_csv),
        "alignment_hits": _file_record(args.hits),
        "selection_audit": _file_record(args.selection),
        "pgfam_gene_cache": _file_record(args.gene_cache, hash_it=False),
        "kbase_genome_cache_dir": {"path": args.genome_cache_dir,
                                   "n_cached_genomes": len(os.listdir(args.genome_cache_dir))
                                   if os.path.isdir(args.genome_cache_dir) else None},
        "build_manifest": _file_record(man_path),
    }
    if args.asv_fasta and os.path.exists(args.asv_fasta):
        inputs["asv_fasta"]["n_sequences"] = sum(1 for l in open(args.asv_fasta) if l.startswith(">"))

    astats = json.load(open(args.alignment_stats)) if args.alignment_stats and os.path.exists(args.alignment_stats) else None
    prefetch_gaps = json.load(open(args.prefetch_gaps)) if args.prefetch_gaps and os.path.exists(args.prefetch_gaps) else None
    n_hits_asvs = sum(1 for a in hits if not a.startswith("_"))

    # ---- KBase genome stats from manifest ----
    built = [a for a in asvs if manifest.get(a, {}).get("status") == "built"]
    nfeat = [manifest[a].get("n_features") for a in built if manifest[a].get("n_features") is not None]
    ngen = Counter(manifest[a].get("n_genomes") for a in built)

    prov = {
        "generated_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "generated_by": f"{os.path.basename(__file__)} (prFBA)",
        "dataset": {"name": "codiffusion_bioreactor 16S V4-V5 ASVs (codif_all)", "n_asvs": len(asvs),
                    "n_asvs_in_hits_file": n_hits_asvs, "genome_dir": gdir},
        "pipeline": {"git": _git_info(), "scripts": _script_hashes(),
                     "kbase_converter": _file_record(os.path.join(os.path.dirname(HERE), "codiffusion_bioreactor",
                                                                  "bvbrc_to_kbase_genome.py")),
                     "build_command": args.build_command},
        "inputs": inputs,
        "stages": {
            "1_alignment": {
                "method": "exhaustive GPU Smith-Waterman (gpu_align.py --full): every ASV vs every unique "
                          "BV-BRC 16S reference; top-20 by local-SW score with deterministic tie-break, "
                          "then Biopython re-alignment for identity / coverage; edlib identity recorded for the "
                          "SW-vs-edlib reliability check",
                "scoring": _sw_constants(),
                "reference_db": {"n_refs": (astats or {}).get("n_refs"),
                                 "source": "BV-BRC 16S rRNA genes (model_inputs/BV_BRC_16S.json), unique sequences"},
                "stats": astats,
                "code_status_vs_head": _alignment_code_status(args.hits),
            },
            "2_gene_family_prefetch": {
                "method": "prefetch_gene_families.py: BV-BRC PGFam gene sets for every stage-7 candidate genome "
                          "(cached in bvbrc_cache/genome_gene_families.json); genomes with an empty set are "
                          "16S-only / no assembly and are kept by the union as taxonomic reps flagged gene_set_empty",
                "unfetchable_candidates": prefetch_gaps,
            },
            "3_selection": {
                "method": "select_references.py exact-tie de-duplicated union (DEFAULT_KNOBS) with the P3 "
                          "gene-richest-conspecific species-dedup tie-break; exact gene-set novelty from the "
                          "BV-BRC PGFam provider (gene_data=bvbrc)",
                "selection_mode": "exact-tie", "p3_gene_tiebreak": True, "guarded": False,
                "gene_data": sel_meta.get("gene_data"),
                "knobs": knobs,
                "tier_thresholds": {"species": knobs["t_species"], "genus": knobs["t_genus"],
                                    "family": knobs["t_family"], "below_family_abstains": True},
                "abstain_rules": {
                    "abstain_below_family": f"best identity < t_family ({knobs['t_family']}) -> no genome-level model",
                    "abstain_reliability": f"every top-20 hit fails the reliability gate: coverage < cov_min "
                                           f"({knobs['cov_min']}) OR |SW - edlib identity| > swed_max ({knobs['swed_max']}) "
                                           f"OR identity < family_floor ({knobs['family_floor']})",
                },
                "summary": sel_meta.get("summary"),
                "tier_counts": dict(tier_counts), "handful_size_counts": {str(k): v for k, v in sorted(handful.items())},
                "flag_counts": dict(flag_counts),
                "audit_file": os.path.abspath(args.selection),
                "per_asv_audit_fields": "asv_len, best_identity, midas_taxonomy, rel_ab, tier, consensus_family, "
                                        "n_selected, flags, selected[] (genome_id, identity, role, confidence, "
                                        "est_genes, gene_set_empty, equivalent_genomes[]), provenance{}, "
                                        "decisions[] (every top-20 hit: disposition, stage, reason)",
            },
            "4_synthetic_genome_build": {
                "method": "build_synthetic_genomes_parallel.py: per-ASV union of the selected genomes' genes by "
                          "function (bvbrc_to_kbase_genome.LocalGenomeConverter.create_synthetic_genome), each "
                          "feature tagged probability = fraction of source genomes carrying the function, "
                          "selection_provenance core/accessory, source_genome_count",
                "manifest_meta": man_meta,
                "n_built": len(built), "n_abstain": man_meta.get("abstain"), "n_failed": man_meta.get("failed"),
                "n_source_genomes_per_synthetic_genome": {str(k): v for k, v in sorted(ngen.items())},
                "features_per_genome": {"min": min(nfeat) if nfeat else None, "max": max(nfeat) if nfeat else None,
                                        "median": sorted(nfeat)[len(nfeat) // 2] if nfeat else None},
                "outputs": {"genome_json": f"{gdir}/{{asv}}.json (KBaseGenomes.Genome-shaped)",
                            "per_asv_provenance": f"{gdir}/provenance/{{asv}}.json (source genomes, equivalents, "
                                                  "core/accessory functions, function prevalence)",
                            "taxonomy": f"{gdir}/taxonomies/{{asv}}.json", "manifest": man_path},
                "build_log_excerpt": _parse_build_log(args.build_log),
            },
        },
        "abstention": {"n_abstained": len(abst_rows), "by_reason_category": dict(abst_cat),
                       "category_definitions": CATEGORY_TEXT,
                       "report_csv": os.path.join(gdir, "abstained_asvs.csv"),
                       "report_md": os.path.join(gdir, "abstained_asvs.md")},
        "environment": _environment(),
    }

    with open(os.path.join(gdir, "mapping_provenance.json"), "w") as f:
        json.dump(prov, f, indent=1)

    with open(os.path.join(gdir, "asv_mapping_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sum_rows[0].keys()))
        w.writeheader(); w.writerows(sum_rows)

    with open(os.path.join(gdir, "abstained_asvs.csv"), "w", newline="") as f:
        cols = list(abst_rows[0].keys()) if abst_rows else ["asv"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader(); w.writerows(abst_rows)

    # markdown report
    g = prov["pipeline"]["git"]
    L = [f"# Abstained ASVs — {prov['dataset']['name']}", "",
         f"prFBA commit `{g.get('commit_short')}` ({g.get('commit_date')}), generated {prov['generated_at']}.", "",
         f"**{len(abst_rows)} of {len(asvs)} ASVs ({len(abst_rows)/max(len(asvs),1):.1%}) received no synthetic genome.** "
         f"An ASV is abstained when no BV-BRC 16S reference passes the selection gates; the pipeline "
         f"prefers no model over a model built from another organism's genome.", "",
         "## Rules", "",
         f"- `abstain_below_family`: best identity < t_family = {knobs['t_family']} (16S family threshold, Yarza 2014).",
         f"- `abstain_reliability`: every top-20 hit failed the reliability gate — coverage < {knobs['cov_min']} "
         f"(partial/local alignment), |SW − edlib identity| > {knobs['swed_max']} (method-discordant / possible "
         f"chimera), or identity < family_floor = {knobs['family_floor']}.", "",
         "## Counts by reason", "", "| reason category | n | definition |", "|---|---:|---|"]
    for c, n in abst_cat.most_common():
        L.append(f"| {c} | {n} | {CATEGORY_TEXT.get(c, '')} |")
    tot_ab = sum((r["rel_ab"] or 0) for r in abst_rows)
    L += ["", f"Summed relative abundance of abstained ASVs (as recorded in the hits file): {tot_ab:.4g}", "",
          "## Per-ASV", "",
          "| ASV | rel_ab | len | best id | tier | MiDAS taxonomy | top hit (BV-BRC) | reason | explanation |",
          "|---|---:|---:|---:|---|---|---|---|---|"]
    for r in sorted(abst_rows, key=lambda x: -(x["rel_ab"] or 0)):
        mt = (r["midas_taxonomy"] or "")
        if isinstance(mt, (list, dict)):
            mt = json.dumps(mt)
        L.append(f"| `{r['asv']}` | {r['rel_ab'] if r['rel_ab'] is not None else ''} | {r['asv_len']} | "
                 f"{r['best_identity']} | {r['tier']} | {str(mt)[:80]} | {r['top_hit_organism']} ({r['top_hit_genome_id']}) "
                 f"id={r['top_hit_identity']} cov={r['top_hit_coverage']} | {r['reason_category']} | {r['explanation']} |")
    with open(os.path.join(gdir, "abstained_asvs.md"), "w") as f:
        f.write("\n".join(L) + "\n")

    print(f"[provenance] {len(asvs)} ASVs | built={len(built)} abstained={len(abst_rows)} "
          f"{dict(abst_cat)} | tiers={dict(tier_counts)} | flags={dict(flag_counts)}")
    print(f"[provenance] wrote mapping_provenance.json, asv_mapping_summary.csv, abstained_asvs.csv/.md -> {gdir}")
    print(f"[provenance] alignment code status: {prov['stages']['1_alignment']['code_status_vs_head'].get('verdict')}")


if __name__ == "__main__":
    main()
