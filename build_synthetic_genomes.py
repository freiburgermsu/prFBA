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

    python build_synthetic_genomes.py \
        --selection /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_reference_selection.json \
        --out-dir   /home/freiburger/Documents/EmilyKin/synthetic_genomes \
        --limit 2          # smoke test; drop --limit for the full run

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
    """Lightweight source-genome dict: gc + taxonomy + CDS features (functions, pgfam).
    No sequence download. The merge dedups by 'functions', so this is the full gene union."""
    core = _api_json("genome",
                     f"eq(genome_id,{gid})&select(genome_id,genome_name,gc_content,taxon_lineage_names)&limit(1)")
    core = core[0] if core else {}
    lineage = core.get("taxon_lineage_names") or []
    feats = _api_json("genome_feature",
                      f"eq(genome_id,{gid})&eq(feature_type,CDS)&select(patric_id,product,pgfam_id,gene)&limit(25000)")
    features = []
    for f in feats:
        prod = (f.get("product") or "").strip()
        if not prod:
            continue
        features.append({
            "functions": [prod], "pgfam_id": f.get("pgfam_id"),
            "patric_id": f.get("patric_id"), "gene": f.get("gene"),
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

    with open(f"{asv}.json", "w") as f:
        json.dump(synth, f)
    with open(os.path.join("provenance", f"{asv}.json"), "w") as f:
        json.dump({"asv": asv, "genomes": ok,
                   "core_functions": sorted(fn for fn, c in prevalence.items()
                                            if c >= max(1, (len(ok) + 1) // 2)),
                   "accessory_functions": sorted(fn for fn, c in prevalence.items() if c == 1),
                   "function_prevalence": dict(prevalence)}, f)

    return {"status": "built", "tier": rec["tier"], "flags": rec["flags"],
            "consensus_family": rec.get("consensus_family"),
            "n_genomes": len(ok), "genomes": ok, "failed": failed,
            "scientific_name": synth.get("scientific_name"), "taxonomy": synth.get("taxonomy"),
            "n_features": len(synth.get("features", [])), "gene_union": summary,
            "genome_file": f"{asv}.json", "provenance_file": f"provenance/{asv}.json"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--genome-cache-dir", default="kbase_genome_cache")
    ap.add_argument("--builder", choices=["fast", "full"], default="fast")
    ap.add_argument("--limit", type=int, default=0, help="first N ASVs (smoke test)")
    args = ap.parse_args()
    _patch_taxonomy_viz()
    if args.builder == "full":
        _patch_go_field()

    out_dir = os.path.abspath(args.out_dir)
    selection = os.path.abspath(args.selection)
    genome_cache_dir = os.path.abspath(args.genome_cache_dir)
    audit = _load(selection)
    audit.pop("_meta", None)
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

    with open("synthetic_genome_manifest.json", "w") as f:
        json.dump({"_meta": {"selection": selection, "builder": args.builder,
                             "n_asvs": len(asvs), "built": built, "abstain": abst,
                             "failed": failed}, **manifest}, f, indent=1)
    print(f"[build] done: built={built} abstain={abst} failed={failed} "
          f"-> {os.path.join(out_dir, 'synthetic_genome_manifest.json')}")


if __name__ == "__main__":
    main()
