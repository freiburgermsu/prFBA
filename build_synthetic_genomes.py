#!/usr/bin/env python
"""
build_synthetic_genomes.py — turn the ASV reference-selection audit into actual per-ASV
synthetic KBase genomes (the union of the selected genomes' genes), with gene-level
core/accessory provenance, by driving codiffusion_bioreactor/bvbrc_to_kbase_genome.py's
create_synthetic_genome().

For each ASV in asv_reference_selection.json with n_selected > 0:
  1. build each selected genome_id into a source-genome dict;
  2. merge them via LocalGenomeConverter.create_synthetic_genome(asv, genomes=[...])
     (which unions features by gene FUNCTION and aggregates a consensus taxonomy);
  3. tag every merged feature core/accessory from the FUNCTION prevalence across the
     ASV's selected genomes (strategy 3) — so the downstream metabolic reconstruction can
     keep accessory genes as gap-fill candidates, not hard reactions;
  4. write the synthetic genome + a per-ASV provenance sidecar + a run manifest.
Abstained ASVs (n_selected == 0) get a taxonomy-only manifest entry, no genome.

Builders:
  --builder fast (default)  2 BV-BRC API calls/genome (product + pgfam, NO sequence
                            download); the merge unions by function, so this yields the
                            full gene union fast (~0.5 s/genome, cached).
  --builder full            the canonical BVBRCToKBaseConverter.build_kbase_genome (fetches
                            contig + protein sequences too) for downstream KBase modeling.
                            (Patched here for a BV-BRC 'go'-field that is now a list.)

EXACT BY DEFAULT: give --hits and the whole pipeline runs exact + cached — it prefetches
the candidate genomes' gene sets (bvbrc_cache/genome_gene_families.json), runs selection with exact
gene-set novelty (writing asv_reference_selection.json with gene_data:"bvbrc"), then builds.
Building from a prebuilt --selection that is gene_data:"estimated" is refused unless
--allow-estimated, so synthetic genomes are never silently built from estimate errors.

    # exact end-to-end (prefetch -> exact selection -> build), all cached:
    python build_synthetic_genomes.py \
        --hits    /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_top20_alignment_hits.json \
        --out-dir /home/freiburger/Documents/EmilyKin/synthetic_genomes \
        --limit 5          # smoke test; drop --limit for the full run

    # or build from an already-exact selection audit:
    python build_synthetic_genomes.py --selection asv_reference_selection.json --out-dir OUT

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
import urllib.request
from collections import Counter

import prefetch_gene_families as PF  # noqa: E402
import select_references as S  # noqa: E402

_MERGE_DIR = "/home/freiburger/Documents/codiffusion_bioreactor"
if _MERGE_DIR not in sys.path:
    sys.path.insert(0, _MERGE_DIR)
import bvbrc_to_kbase_genome as KB  # noqa: E402

API = "https://www.bv-brc.org/api"


# ---- disable the per-ASV matplotlib taxonomy PNG (side-effect, breaks headless batch;
#      consensus-taxonomy aggregation + the taxonomy JSON are kept) ----
def _patch_taxonomy_viz():
    KB.LocalGenomeConverter.visualize_taxonomies = lambda self, *a, **k: None


# ---- BV-BRC 'full' builder fix: the API now returns the 'go' field as a list ----
def _patch_go_field():
    orig = KB.BVBRCToKBaseConverter.fetch_genome_features

    def patched(self):
        feats = orig(self)
        for f in feats:
            v = f.get("go")
            if isinstance(v, list):
                f["go"] = ",".join(map(str, v))
        return feats

    KB.BVBRCToKBaseConverter.fetch_genome_features = patched


# --------------------------------------------------------------------------- io
def _load(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        return json.load(f)


def _api_json(core, query, timeout=60, retries=3):
    url = f"{API}/{core}/?{query}&http_accept=application/json"
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return json.load(resp)
        except Exception as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise last


# --------------------------------------------------------------------------- builders
def fast_build_genome(gid):
    """Lightweight source-genome dict reproducing the data of the prior
    codiffusion_bioreactor/genome_objects run: per-CDS functions + KBase aliases
    (PATRIC_id then FIGFAM/PGFAM/PLFAM, non-empty only) + feature_type, with NO sequences
    (dna/protein empty — exactly as that run acquired; it deliberately skipped sequence data).
    The merge keys on 'functions' and copies type/aliases through, so create_synthetic_genome
    yields per-ASV objects identical in shape to those genome_objects. Mirrors
    bvbrc_to_kbase_genome._convert_bvbrc_feature_to_kbase's field mapping."""
    core = _api_json("genome",
                     f"eq(genome_id,{gid})&select(genome_id,genome_name,gc_content,taxon_lineage_names)&limit(1)")
    core = core[0] if core else {}
    lineage = core.get("taxon_lineage_names") or []
    feats = _api_json("genome_feature",
                      f"eq(genome_id,{gid})&eq(feature_type,CDS)"
                      f"&select(patric_id,product,feature_type,figfam_id,pgfam_id,plfam_id,gene)&limit(25000)")
    features = []
    for f in feats:
        prod = (f.get("product") or "").strip()
        if not prod:
            continue
        aliases = [["PATRIC_id", f.get("patric_id", "")]]
        for fam, key in (("figfam_id", "FIGFAM"), ("pgfam_id", "PGFAM"), ("plfam_id", "PLFAM")):
            v = f.get(fam)
            if v:
                aliases.append([key, v])
        features.append({
            "type": f.get("feature_type") or "CDS",
            "functions": [prod],
            "aliases": aliases,
            "pgfam_id": f.get("pgfam_id"), "gene": f.get("gene"),  # kept for provenance (merge ignores)
            "dna_sequence": "", "protein_translation": "",
        })
    return {
        "id": str(gid), "scientific_name": core.get("genome_name", str(gid)),
        "taxonomy": "; ".join(lineage), "gc_content": float(core.get("gc_content") or 0.0) or 0.5,
        "features": features,
    }


def build_genome(gid, cache_dir, builder):
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{gid}.json")
    if os.path.exists(path):
        return _load(path)
    g = fast_build_genome(gid) if builder == "fast" else \
        KB.BVBRCToKBaseConverter(gid, base_url=API).build_kbase_genome()
    g.setdefault("id", str(gid))
    with open(path, "w") as f:
        json.dump(g, f)
    return g


# --------------------------------------------------------------------------- provenance
def function_prevalence(source_genomes):
    """How many of the ASV's source genomes carry each gene function (the merge's key)."""
    cnt = Counter()
    for g in source_genomes:
        funcs = set()
        for ft in g.get("features", []):
            funcs.update(ft.get("functions") or [])
        cnt.update(funcs)
    return cnt


def annotate_and_summarize(genome, prevalence, n_genomes):
    """Tag each merged feature by pangenome class (function prevalence across the ASV's
    selected genomes): core = in ALL, accessory = unique to one, shared = intermediate.
    Accessory genes are the gap-fill candidates for the downstream FBA reconstruction."""
    tags = Counter()
    for f in genome.get("features", []):
        funcs = f.get("functions") or []
        c = max((prevalence.get(fn, 0) for fn in funcs), default=0)
        if c == 0:
            tag = "unknown"
        elif c >= n_genomes:        # in every selected genome
            tag = "core"
        elif c == 1:                # unique to a single genome (only when n_genomes > 1)
            tag = "accessory"
        else:
            tag = "shared"
        f["selection_provenance"] = tag
        f["source_genome_count"] = int(c)
        tags[tag] += 1
    vals = list(prevalence.values())
    return dict(n_genomes=n_genomes, n_functions_union=len(vals),
                n_core=sum(1 for c in vals if c >= n_genomes),
                n_shared=sum(1 for c in vals if 1 < c < n_genomes),
                n_accessory=sum(1 for c in vals if c == 1 and c < n_genomes),
                feature_tags=dict(tags))


# --------------------------------------------------------------------------- per-ASV
def build_for_asv(asv, rec, *, genome_cache_dir, builder):
    """Runs with cwd == the output dir, so the merger's hardcoded relative side-effect
    dirs (genome_objects/, taxonomies/) land under the output dir."""
    gids = [str(s["genome_id"]) for s in rec["selected"]]
    sources, ok, failed = [], [], []
    for g in gids:
        try:
            sources.append(build_genome(g, genome_cache_dir, builder))
            ok.append(g)
        except Exception as exc:
            failed.append({"genome_id": g, "error": str(exc)[:200]})
    if not sources:
        return {"status": "failed", "tier": rec["tier"], "selected": gids,
                "failed": failed, "reason": "no source genome could be built"}

    conv = KB.LocalGenomeConverter()
    synth = conv.create_synthetic_genome(
        asv, genomes=sources, template_file=None, save_taxonomy=True, taxonomy=None,
        taxonomy_output_dir="taxonomies", genome_output_dir="genome_objects", verbose=False)

    prevalence = function_prevalence(sources)
    summary = annotate_and_summarize(synth, prevalence, len(ok))

    # P2: the 16S-indistinguishable same-species organisms each selected rep stands in for
    # (recorded by select_references; the source organism, in a self-recovery run, appears
    # here). These do NOT enter the gene union — they are honest provenance that the marker
    # cannot resolve these genomes apart.
    equivalents = [e for s in rec["selected"] for e in (s.get("equivalent_genomes") or [])]

    with open(f"{asv}.json", "w") as f:
        json.dump(synth, f)
    with open(os.path.join("provenance", f"{asv}.json"), "w") as f:
        json.dump({"asv": asv, "genomes": ok, "equivalent_genomes": equivalents,
                   "core_functions": sorted(fn for fn, c in prevalence.items()
                                            if c >= max(1, (len(ok) + 1) // 2)),
                   "accessory_functions": sorted(fn for fn, c in prevalence.items() if c == 1),
                   "function_prevalence": dict(prevalence)}, f)

    return {"status": "built", "tier": rec["tier"], "flags": rec["flags"],
            "consensus_family": rec.get("consensus_family"),
            "n_genomes": len(ok), "genomes": ok, "failed": failed,
            "n_equivalent_genomes": len(equivalents),
            "scientific_name": synth.get("scientific_name"), "taxonomy": synth.get("taxonomy"),
            "n_features": len(synth.get("features", [])), "gene_union": summary,
            "genome_file": f"{asv}.json", "provenance_file": f"provenance/{asv}.json"}


def exact_selection(hits, gene_cache, *, workers=12, prefetch=True, knobs=None):
    """Whole-pipeline exact path: prefetch the candidate genomes' gene sets (cached),
    then run selection with the exact gene-set novelty provider (no estimator)."""
    if prefetch:
        PF.populate_cache(hits, gene_cache, workers=workers)
    provider = S.bvbrc_gene_provider(gene_cache)  # exact gene-set Jaccard novelty
    # P3: gene-richest-conspecific species-dedup tie-break (fetch-free, from the loaded cache)
    audit, summary = S.build_audit(hits, gene_provider=provider,
                                   gene_count_fn=S.gene_count_from_provider(provider), knobs=knobs)
    return audit, summary


def get_selection(args, *, gene_cache, selection_out):
    """Return the {asv: record} audit to build from. Exact by default:
    - with --hits: prefetch + exact selection, written to selection_out (authoritative).
    - with --selection: build from a prebuilt audit, but REFUSE an 'estimated' one unless
      --allow-estimated (so we never silently build from estimate-error decisions)."""
    if args.hits:
        hits = json.load(open(os.path.abspath(args.hits)))
        if args.limit:  # restrict selection+prefetch to the same ASVs we will build
            hits = {a: hits[a] for a in list(hits)[:args.limit]}
        if args.estimated:
            print("[select] ESTIMATOR mode (not exact) — testing only")
            audit, summary = S.build_audit(hits)
            gene_data = "estimated"
        else:
            print("[select] exact gene-set novelty (prefetch + cache)")
            audit, summary = exact_selection(hits, gene_cache, workers=args.workers,
                                             prefetch=not args.no_prefetch)
            gene_data = "bvbrc"
        with open(selection_out, "w") as f:
            json.dump({"_meta": {"gene_data": gene_data, "generated_from": os.path.abspath(args.hits),
                                 "summary": summary}, **audit}, f, indent=1)
        print(f"[select] wrote {selection_out} (gene_data={gene_data}, {len(audit)} ASVs)")
        return audit, gene_data
    data = _load(os.path.abspath(args.selection))
    meta = data.pop("_meta", {}) or {}
    gd = meta.get("gene_data", "unknown")
    if gd == "estimated" and not args.allow_estimated:
        raise SystemExit(
            "[build] refusing to build from an ESTIMATED selection (estimate errors). "
            "Pass --hits HITS to (re)generate the exact selection, or --allow-estimated to override.")
    return data, gd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_argument_group("selection source (give --hits for the exact end-to-end pipeline)")
    src.add_argument("--hits", help="asv_top20_alignment_hits.json -> prefetch + EXACT selection + build")
    src.add_argument("--selection", help="a prebuilt asv_reference_selection.json to build from")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--gene-cache", default="bvbrc_cache/genome_gene_families.json", help="PGFam cache (exact novelty)")
    ap.add_argument("--selection-out", default=None, help="where to write the exact audit (with --hits)")
    ap.add_argument("--genome-cache-dir", default="bvbrc_cache/kbase_genome_cache")
    ap.add_argument("--builder", choices=["fast", "full"], default="fast")
    ap.add_argument("--workers", type=int, default=12, help="prefetch threads")
    ap.add_argument("--estimated", action="store_true", help="(with --hits) use the estimator — testing only")
    ap.add_argument("--no-prefetch", action="store_true", help="(with --hits) don't fetch missing gene sets")
    ap.add_argument("--allow-estimated", action="store_true", help="(with --selection) permit an estimated audit")
    ap.add_argument("--limit", type=int, default=0, help="first N ASVs (smoke test)")
    args = ap.parse_args()
    if not args.hits and not args.selection:
        ap.error("provide --hits (exact end-to-end) or --selection (prebuilt audit)")
    _patch_taxonomy_viz()
    if args.builder == "full":
        _patch_go_field()

    out_dir = os.path.abspath(args.out_dir)
    genome_cache_dir = os.path.abspath(args.genome_cache_dir)
    gene_cache = os.path.abspath(args.gene_cache)
    selection_out = os.path.abspath(args.selection_out) if args.selection_out else (
        os.path.join(os.path.dirname(os.path.abspath(args.hits)), "asv_reference_selection.json")
        if args.hits else None)

    audit, gene_data = get_selection(args, gene_cache=gene_cache, selection_out=selection_out)

    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)  # contain the merger's hardcoded relative side-effect dirs here
    for d in ("genome_objects", "taxonomies", "provenance", "ASVset_taxonomies"):
        os.makedirs(d, exist_ok=True)
    asvs = list(audit)[:args.limit] if args.limit else list(audit)
    print(f"[build] {len(asvs)} ASVs | builder={args.builder} | out={out_dir}")

    manifest, t0 = {}, time.perf_counter()
    built = abst = failed = 0
    for i, asv in enumerate(asvs, 1):
        rec = audit[asv]
        if rec["n_selected"] == 0:
            manifest[asv] = {"status": "abstain", "tier": rec["tier"], "flags": rec["flags"],
                             "midas_taxonomy": rec.get("midas_taxonomy")}
            abst += 1
            continue
        try:
            entry = build_for_asv(asv, rec, genome_cache_dir=genome_cache_dir, builder=args.builder)
        except Exception as exc:
            entry = {"status": "failed", "tier": rec["tier"], "error": str(exc)[:300]}
        manifest[asv] = entry
        built += entry["status"] == "built"
        failed += entry["status"] == "failed"
        if i % 25 == 0 or i == len(asvs):
            print(f"[build] {i}/{len(asvs)} built={built} abstain={abst} failed={failed} "
                  f"({time.perf_counter()-t0:.0f}s)", flush=True)

    sel_src = selection_out if args.hits else os.path.abspath(args.selection)
    with open("synthetic_genome_manifest.json", "w") as f:
        json.dump({"_meta": {"selection": sel_src, "gene_data": gene_data,
                             "builder": args.builder, "n_asvs": len(asvs), "built": built,
                             "abstain": abst, "failed": failed}, **manifest}, f, indent=1)
    print(f"[build] done: built={built} abstain={abst} failed={failed} "
          f"-> {os.path.join(out_dir, 'synthetic_genome_manifest.json')}")


if __name__ == "__main__":
    main()
