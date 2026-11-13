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
                            (Patched here for a BV-BRC 'go'-field that is now a list, and to
                            keep only annotation == 'PATRIC' feature rows.)

BV-BRC ANNOTATION ONLY (since 2026-11-13; no option to include RefSeq). BV-BRC's genome_feature
core returns two annotations per genome: annotation=PATRIC (BV-BRC's RASTtk annotation; every
CDS has a patric_id 'fig|...', ~90% a pgfam_id) and annotation=RefSeq (imported NCBI PGAP; no
patric_id, no pgfam_id). Every source genome dict is PATRIC-only: fast_build_genome queries
eq(annotation,PATRIC), the full builder keeps annotation == 'PATRIC' rows, and every dict read
from the genome cache or freshly built passes through bvbrc_only() (features with a non-empty
['PATRIC_id', ...] alias). Cache files are never rewritten; files cached before this change
still hold both annotations and are filtered on read. So the merge, k/n, core/accessory,
function_prevalence, E3's product-match route and the source digests all see PATRIC features
only. The manifest _meta records 'annotation': ANNOTATION_STAMP in every mode. Builds made
before this change (every build up to and including git 3ce12d0, e.g. the shipped
synthetic_genomes_d100 deployments and their E3 re-weightings) merged both vocabularies and
are reproducible only from git at or before 3ce12d0 (with the cache files of that time).

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

FEATURE PROBABILITY (--pfunc, default e3). The merge writes p = k/n (the share of the ASV's n
built source genomes carrying the feature's function string). `--pfunc e3` (the DEFAULT since
the 10k-panel gate passed, 2026-11-12) replaces it, for ASVs with n >= --pfunc-min-n (2)
built source genomes, by the Stage E/O estimator E3 (pfunc_estimator.py: identity-weighted
16S neighbours + NCBI family/order prior, PGFam namespace), in BOTH written copies
({asv}.json and genome_objects/{asv}.json, features and CDS); features without a PGFam and
genomes with n < min_n keep k/n. The replaced k/n values go to provenance/{asv}.json
'probability_kn', the estimator diagnostics to its 'pfunc' key and to the manifest entry,
and the run settings (with the code's sha256 and git state) to the manifest _meta 'pfunc'.
Before the build, the PGFam sets the run needs are stream-extracted from the shared cache
(never loaded whole, never rewritten) into OUT/pfunc/pgfam_store.json; missing genomes are
fetched from BV-BRC into OUT/pfunc/pgfam_fetched.json (--pfunc-no-fetch to skip; with every
needed genome cached, no network is used). The E3 clade catalogue
(bvbrc_cache/e3_clade_catalogue.json) must exist and match the PGFam cache: the run stops
BEFORE selection otherwise, printing the rebuild command (--pfunc-allow-stale-catalogue
overrides a stale one). `--pfunc kn` keeps k/n everywhere: the genome, provenance and
taxonomy files (and the manifest entries) are byte-identical to what this code writes without
E3 (pfunc=None); the manifest _meta gains a 'pfunc' key {estimator: kn, estimator_id: k/n,
gating_rule, min_n}. It does not reproduce the shipped k/n builds: they mix RefSeq features
(see BV-BRC ANNOTATION ONLY above).
reweight_build_e3.py applies E3 to an existing k/n build (copy-on-write) with the same code.

INPUT STAMPS (every --pfunc mode, for reweight_build_e3.py's drift checks; no per-ASV file
changes): the manifest _meta gains 'selection_stamp' = {sha256, size, mtime_ns} of the
selection audit as read (--selection) or written (--hits), and 'source_digests' (a record of)
OUT/pfunc/source_digests.json = {genome_id: pfunc_estimator.source_digest} for every built
source genome (the dicts as loaded, before the merge; a genome whose digest differed between
two ASVs of the run is stored as null and listed in _meta.source_digests).

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.request
from collections import Counter

import pfunc_estimator as P  # noqa: E402
import prefetch_gene_families as PF  # noqa: E402
import select_references as S  # noqa: E402

_MERGE_DIR = "/home/freiburger/Documents/codiffusion_bioreactor"
if _MERGE_DIR not in sys.path:
    sys.path.insert(0, _MERGE_DIR)
import bvbrc_to_kbase_genome as KB  # noqa: E402

API = "https://www.bv-brc.org/api"

# ---- BV-BRC annotation: PATRIC (RASTtk) only, never the imported RefSeq (NCBI PGAP) ----
# BV-BRC's genome_feature core holds two annotations per genome: annotation=PATRIC (BV-BRC's
# own RASTtk annotation: every CDS has a patric_id 'fig|...', ~90% a pgfam_id) and
# annotation=RefSeq (imported NCBI PGAP: no patric_id, no pgfam_id). Every synthetic-genome
# source dict is PATRIC-only (2026-11-13): the fast builder queries eq(annotation,PATRIC), the
# full builder drops non-PATRIC rows (_patch_go_field), and every dict read from the genome
# cache or freshly built passes through bvbrc_only(). There is no option to include RefSeq.
BVBRC_ANNOTATION = "PATRIC"
ANNOTATION_STAMP = P.ANNOTATION_STAMP   # manifest _meta['annotation'] of every build (all modes)
LAST_MIXED_COMMIT = "3ce12d0"       # builds made at or before it mix PATRIC + RefSeq features
FEATURE_COLLECTIONS = ("features", "cdss", "non_coding_features", "mrnas")


def is_bvbrc_feature(feature):
    """True for a feature of BV-BRC's own PATRIC annotation: it has a ['PATRIC_id', id] alias
    with a non-empty id. A RefSeq (PGAP) feature carries ['PATRIC_id', ''] (both builders copy
    the empty patric_id through); a feature without any PATRIC_id alias is not one either."""
    for a in feature.get("aliases") or ():
        if isinstance(a, (list, tuple)) and len(a) >= 2 and a[0] == "PATRIC_id" and a[1]:
            return True
    return False


def bvbrc_only(genome):
    """The genome dict with only BV-BRC PATRIC-annotation features (is_bvbrc_feature) in each
    of its feature collections (features; and cdss / non_coding_features / mrnas, which only
    full-builder dicts carry). Returns a new dict (the input is not modified); every other key
    is kept as is. Idempotent. Applied to every source genome a build uses: the dicts read
    from the genome cache (whose older files mix both annotations and are never rewritten)
    and every freshly built one."""
    out = dict(genome)
    for c in FEATURE_COLLECTIONS:
        if isinstance(genome.get(c), list):
            out[c] = [f for f in genome[c] if is_bvbrc_feature(f)]
    return out


# ---- disable the per-ASV matplotlib taxonomy PNG (side-effect, breaks headless batch;
#      consensus-taxonomy aggregation + the taxonomy JSON are kept) ----
def _patch_taxonomy_viz():
    KB.LocalGenomeConverter.visualize_taxonomies = lambda self, *a, **k: None


# ---- BV-BRC 'full' builder fixes: the API now returns the 'go' field as a list, and its
#      genome_feature rows carry two annotations (keep annotation == 'PATRIC' only) ----
def _patch_go_field():
    """Patch BVBRCToKBaseConverter.fetch_genome_features (in this process only; the
    codiffusion file is not edited): keep only the BV-BRC PATRIC annotation's rows
    (annotation == 'PATRIC'; the converter's query has no annotation filter, so it also
    returns the imported RefSeq rows) and join a list-valued 'go' field. Idempotent.
    fetch_source_genome applies it for every full build, so no full build sees a RefSeq row."""
    cls = KB.BVBRCToKBaseConverter
    if getattr(cls.fetch_genome_features, "_prfba_patched", False):
        return
    orig = cls.fetch_genome_features

    def patched(self):
        feats = [f for f in orig(self) if f.get("annotation") == BVBRC_ANNOTATION]
        for f in feats:
            v = f.get("go")
            if isinstance(v, list):
                f["go"] = ",".join(map(str, v))
        return feats

    patched._prfba_patched = True
    cls.fetch_genome_features = patched


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
    yields per-ASV objects identical in shape to those genome_objects (not in content: that
    run, like every build up to 3ce12d0, also merged RefSeq CDSs). Mirrors
    bvbrc_to_kbase_genome._convert_bvbrc_feature_to_kbase's field mapping.

    BV-BRC PATRIC annotation only: the CDS query carries eq(annotation,PATRIC) and the dict is
    passed through bvbrc_only() (so a row without a patric_id never gets in, whatever the API
    returns). Up to 3ce12d0 the query had no annotation filter and every dict also held the
    imported RefSeq CDSs (empty patric_id, no PGFam)."""
    core = _api_json("genome",
                     f"eq(genome_id,{gid})&select(genome_id,genome_name,gc_content,taxon_lineage_names)&limit(1)")
    core = core[0] if core else {}
    lineage = core.get("taxon_lineage_names") or []
    feats = _api_json("genome_feature",
                      f"eq(genome_id,{gid})&eq(feature_type,CDS)&eq(annotation,{BVBRC_ANNOTATION})"
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
    return bvbrc_only({
        "id": str(gid), "scientific_name": core.get("genome_name", str(gid)),
        "taxonomy": "; ".join(lineage), "gc_content": float(core.get("gc_content") or 0.0) or 0.5,
        "features": features,
    })


def load_source_genome(path):
    """A source genome dict read from the genome cache, BV-BRC PATRIC features only
    (bvbrc_only). The cache file is never rewritten: files cached up to 3ce12d0 also hold the
    RefSeq CDSs, and the filter is applied here, on every read."""
    return bvbrc_only(_load(path))


def fetch_source_genome(gid, builder):
    """A freshly built source genome dict (BV-BRC PATRIC features only), for writing to the
    genome cache: fast_build_genome, or the codiffusion converter with its feature fetch
    patched to annotation == 'PATRIC' (_patch_go_field)."""
    if builder == "fast":
        g = fast_build_genome(gid)
    else:
        _patch_go_field()
        g = KB.BVBRCToKBaseConverter(gid, base_url=API).build_kbase_genome()
    g = bvbrc_only(g)
    g.setdefault("id", str(gid))
    return g


def build_genome(gid, cache_dir, builder):
    """The source genome dict for one genome id: read from cache_dir/{gid}.json when present
    (load_source_genome), else built (fetch_source_genome) and cached. Either way it holds
    only BV-BRC PATRIC-annotation features."""
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{gid}.json")
    if os.path.exists(path):
        return load_source_genome(path)
    g = fetch_source_genome(gid, builder)
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


# --------------------------------------------------------------------------- pfunc (E3)
PFUNC_DIR = "pfunc"                      # per-run E3 inputs, under the out dir
PFUNC_STORE = "pgfam_store.json"         # slim {genome_id: sorted PGFam list} the run used
PFUNC_FETCHED = "pgfam_fetched.json"     # genomes fetched from BV-BRC (never the shared cache)
SOURCE_DIGESTS = "source_digests.json"   # {genome_id: source_digest} of every built source
#                                          genome, written in EVERY mode (kn too)


# --------------------------------------------------------------------------- input stamps
# What reweight_build_e3.py checks to tell whether a build's inputs changed since the build.
# Both are written in every --pfunc mode; neither touches a per-ASV file, so a kn build's
# genome / provenance / taxonomy files stay byte-identical to the same code's build without E3.
def file_stamp(path):
    """{sha256, size, mtime_ns} of a file (the selection audit a build read or wrote), taken
    when the build reads or writes it. None when `path` is not an existing file (the parallel
    builder's --hits --limit run writes no audit)."""
    if not path or not os.path.isfile(path):
        return None
    st = os.stat(path)
    return {"sha256": _sha256_file(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def merge_source_digests(run, new, conflicts):
    """Fold one ASV's {genome_id: source_digest} into the run's map. A genome whose digest
    differs between two ASVs of one run (its cache file changed mid-run) is recorded as None
    (never equal to a real digest, so the retro tool treats it as drift) and its id is added
    to `conflicts`."""
    for g, d in new.items():
        if g not in run:
            run[g] = d
        elif run[g] != d:
            run[g] = None
            conflicts.add(g)


def write_source_digests(out_dir, digests, conflicts=()):
    """Write OUT/pfunc/source_digests.json = {genome_id: source_digest} (sorted ids; None for
    a conflicting genome) and return the manifest _meta['source_digests'] record."""
    pdir = os.path.join(out_dir, PFUNC_DIR)
    os.makedirs(pdir, exist_ok=True)
    _dump_json({g: digests[g] for g in sorted(digests)}, os.path.join(pdir, SOURCE_DIGESTS))
    return {"file": os.path.join(PFUNC_DIR, SOURCE_DIGESTS), "n_genomes": len(digests),
            "conflicting_genome_ids": sorted(conflicts),
            "method": "pfunc_estimator.source_digest: " + P.SOURCE_DIGEST_METHOD}


class PfuncContext:
    """What build_for_asv needs to replace k/n by E3 (see the module docstring). mode 'kn'
    (or pfunc=None) leaves the build exactly as before. In 'e3' mode `store` is a
    pfunc_estimator.PGFamStore (or {gid: set}) and `catalogue` a CladeCatalogue, both
    read-only during the build (shared with fork-started workers)."""

    def __init__(self, mode="kn", min_n=P.MIN_N, dedup="species", store=None, catalogue=None):
        if mode not in ("kn", "e3"):
            raise ValueError(f"pfunc mode must be 'kn' or 'e3', not {mode!r}")
        if dedup not in ("species", "genome"):
            raise ValueError(f"pfunc dedup must be 'species' or 'genome', not {dedup!r}")
        if mode == "e3" and (store is None or catalogue is None):
            raise ValueError("pfunc mode 'e3' needs a PGFam store and a clade catalogue")
        self.mode, self.min_n, self.dedup = mode, int(min_n), dedup
        self.store, self.catalogue = store, catalogue

    @property
    def e3(self):
        return self.mode == "e3"

    def reweight(self, genome_obj, rec, sources, n_built, built_ids=None, source_digests=None):
        """pfunc_estimator.reweight_asv with this context's settings (in place). built_ids:
        the genome ids of `sources` (the clade prior is anchored on the first built rep);
        source_digests: their digests when the caller already has them (else computed)."""
        return P.reweight_asv(genome_obj, rec, sources, self.store, self.catalogue, n_built,
                              mode=self.mode, min_n=self.min_n, dedup=self.dedup,
                              built_ids=built_ids, digests=source_digests)


def add_pfunc_args(ap, default="e3"):
    """The feature-probability flags, shared by the serial and parallel builders, build_d100
    (which forwards them) and reweight_build_e3.py. Default e3 (the 10k-panel gate passed,
    2026-11-12); --pfunc kn keeps k/n everywhere (it does not reproduce builds made at or
    before 3ce12d0, which also merged RefSeq features)."""
    g = ap.add_argument_group("feature probability (--pfunc)")
    g.add_argument("--pfunc", choices=["kn", "e3"], default=default,
                   help=f"feature probability estimator (default {default}): kn = share of the "
                        "built source genomes carrying the function, everywhere; e3 = "
                        "pfunc_estimator E3 for ASVs with >= --pfunc-min-n built source genomes "
                        "(k/n elsewhere). Either way the source genomes are BV-BRC PATRIC "
                        "annotation only")
    g.add_argument("--pfunc-min-n", type=int, default=P.MIN_N,
                   help="E3 replaces k/n only when at least this many source genomes were built "
                        f"(default {P.MIN_N}; Stage O's gain is established only at n >= 2)")
    g.add_argument("--pfunc-dedup", choices=["species", "genome"], default="species",
                   help="E3 neighbour dedup (default species: one neighbour per BV-BRC species)")
    g.add_argument("--pfunc-pgfam-cache", default=P.DEFAULT_PGFAM_CACHE,
                   help="shared PGFam cache streamed (never loaded or written) for the E3 store "
                        "(default the shared bvbrc_cache/genome_gene_families.json; NOT --gene-cache, "
                        "which a d100 build points at the small tie-set cache)")
    g.add_argument("--pfunc-catalogue", default=P.DEFAULT_CATALOGUE,
                   help="E3 clade catalogue (pfunc_estimator.py build-catalogue output); must "
                        "exist and match --pfunc-pgfam-cache, else the run stops before "
                        "selection and prints the rebuild command")
    g.add_argument("--pfunc-allow-stale-catalogue", action="store_true",
                   help="proceed with a clade catalogue that does not match the PGFam cache "
                        "(size/mtime); genomes added to the cache since are then not clade "
                        "members. Recorded in the manifest _meta.pfunc.catalogue")
    g.add_argument("--pfunc-no-fetch", action="store_true",
                   help="do not fetch genomes missing from the PGFam cache from BV-BRC "
                        "(they are then left out of the neighbour list / prior and counted)")
    g.add_argument("--pfunc-fetch-workers", type=int, default=8, help="BV-BRC fetch threads")
    return g


# code whose sha256 / git state a run's _meta.pfunc records (item: estimator provenance);
# wrappers (build_d100.py) append themselves before calling a builder's main()
PFUNC_PROVENANCE_FILES = [os.path.join(P.ROOT, "pfunc_estimator.py"), os.path.abspath(__file__)]


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def code_provenance(files=None):
    """{git_head, dirty, files: {path relative to the prFBA root: {sha256, git_status}}} for
    the estimator and builder code a run used. git_status is `git status --porcelain`'s
    two-letter code ('' clean, ' M' modified, '??' untracked); dirty = any non-clean file
    (an untracked file is dirty: no commit holds it). Git failures leave None."""
    paths = []
    for p in (files or PFUNC_PROVENANCE_FILES):
        p = os.path.abspath(p)
        if p not in paths:
            paths.append(p)

    def git(*a):
        try:
            r = subprocess.run(["git", "-C", P.ROOT, *a], capture_output=True, text=True,
                               timeout=30)
            return r.stdout if r.returncode == 0 else None
        except Exception:
            return None
    head = git("rev-parse", "HEAD")
    out, dirty = {}, False
    for p in paths:
        rel = os.path.relpath(p, P.ROOT)
        st = git("status", "--porcelain", "--", rel)
        code = None if st is None else (st.splitlines()[0][:2] if st.strip() else "")
        if code:
            dirty = True
        out[rel] = {"sha256": _sha256_file(p) if os.path.exists(p) else None, "git_status": code}
    return {"git_head": head.strip() if head else None,
            "dirty": dirty if head is not None else None, "files": out}


def _catalogue_error(why, catalogue, pgfam_cache):
    return SystemExit(
        f"[pfunc] {why}\n"
        f"[pfunc] E3 (the default --pfunc) needs the clade catalogue {catalogue} to match the "
        f"PGFam cache {pgfam_cache}. Rebuild it (one streaming pass over the cache + the "
        f"taxdump; ~45 s, ~2.7 GB) with:\n    {P.catalogue_rebuild_command(pgfam_cache, catalogue)}\n"
        f"[pfunc] or pass --pfunc-allow-stale-catalogue (a stale catalogue misses the genomes "
        f"added to the cache since), or --pfunc kn for the earlier k/n build.")


def pfunc_preflight(mode, catalogue, pgfam_cache, allow_stale=False, log=print):
    """Check, BEFORE selection / cache warming / any merge, that an e3 run can prepare E3:
    the PGFam cache and the clade catalogue exist, and the catalogue is current for the cache
    (pfunc_estimator.catalogue_meta_is_current) unless allow_stale. Raises SystemExit with the
    rebuild command otherwise. Returns None in kn mode, else a snapshot dict that
    prepare_pfunc uses to tell a catalogue that was stale from the start (refused) from one
    that this run's own selection/prefetch made stale by writing new genomes into the PGFam
    cache (warned, recorded as 'staled_by_this_run')."""
    if mode != "e3":
        return None
    catalogue, pgfam_cache = os.path.abspath(catalogue), os.path.abspath(pgfam_cache)
    if not os.path.exists(pgfam_cache):
        raise SystemExit(f"[pfunc] PGFam cache {pgfam_cache} not found (--pfunc-pgfam-cache); "
                         f"E3 cannot run. Use --pfunc kn for the earlier k/n build.")
    if not os.path.exists(catalogue):
        raise _catalogue_error(f"clade catalogue {catalogue} not found", catalogue, pgfam_cache)
    meta = P.catalogue_meta(catalogue)
    current = P.catalogue_meta_is_current(meta, pgfam_cache)
    if not current and not allow_stale:
        raise _catalogue_error(
            f"clade catalogue {catalogue} is STALE: built from a cache of size "
            f"{meta.get('source_cache_size')} / mtime {meta.get('source_cache_mtime')}, the cache "
            f"is now {os.path.getsize(pgfam_cache)} bytes", catalogue, pgfam_cache)
    st = os.stat(pgfam_cache)
    if log:
        log(f"[pfunc] preflight: E3 clade catalogue {catalogue} "
            f"({'current' if current else 'STALE, allowed'}) for {pgfam_cache}")
    return {"catalogue": catalogue, "pgfam_cache": pgfam_cache, "current": current,
            "cache_stat": [st.st_size, st.st_mtime_ns]}


def anchor_candidates(rec, built_ids=None):
    """The reps whose clade members the PGFam store must hold for one ASV: the clade prior is
    anchored on the first BUILT rep (pfunc_estimator.prior_anchor). built_ids known
    (reweight_build_e3: the manifest's built genomes) -> exactly that anchor. Before the merge
    (the builders) it is unknown which reps will build -- a cached genome file can be
    unreadable, an uncached one can fail to fetch, and reps often sit in different families
    (cross-phylum merges) -- so EVERY rep is a candidate. Measured on the d100 deployments'
    n >= 2 ASVs (2026-11-12) this enlarges the store by ~13% (EmilyKin 21,003 -> 23,709
    genomes, codiffusion 20,774 -> 23,262), all of them already in the shared PGFam cache."""
    reps = [str(s["genome_id"]) for s in rec.get("selected") or []]
    if not reps:
        return []
    if built_ids is not None:
        return [P.prior_anchor(rec, built_ids)]
    return reps


def prepare_pfunc(records, out_dir, *, mode="kn", min_n=P.MIN_N, dedup="species",
                  pgfam_cache=P.DEFAULT_PGFAM_CACHE, catalogue=P.DEFAULT_CATALOGUE, fetch=True,
                  fetch_workers=8, n_by_asv=None, built_by_asv=None, allow_stale=False,
                  preflight=None, code_files=None, taxdb="auto", log=print):
    """Build the PfuncContext for a run, ONCE, before any merge (spec: one PGFam pass per run).

    records: {asv: selection-audit record} of the ASVs the run builds. An ASV is in E3 scope
    when its n >= min_n, where n = n_by_asv[asv] when given (reweight_build_e3: the built n
    from the manifest) else rec['n_selected'] (the builder, before the merge: an upper bound
    on the built n; reweight_asv re-checks the built n, so an ASV whose fetches fail to below
    min_n stays k/n). For every in-scope ASV, pfunc_estimator.required_genomes (neighbours +
    the clade members of each anchor_candidates rep: the first built rep when built_by_asv
    gives the built genome ids, else every rep) are collected; the union is filled by
    ONE PGFamStore.ensure (side file OUT/pfunc/pgfam_fetched.json, then one ijson pass over
    `pgfam_cache`, then BV-BRC fetch of the rest when `fetch`), and the store is saved to
    OUT/pfunc/pgfam_store.json. Anchors outside the catalogue are resolved with taxopy
    (taxdb='auto' loads the taxdump only if needed) and cached in the catalogue before any
    worker forks.

    Catalogue rule: missing -> SystemExit with the rebuild command; stale (not matching
    `pgfam_cache`) -> SystemExit unless allow_stale, or unless `preflight` (pfunc_preflight's
    snapshot) shows it was current when the run started and the cache has changed since
    (this run's selection/prefetch wrote to it: warned, recorded as staled_by_this_run).

    Returns (PfuncContext, meta) -- meta is the manifest _meta['pfunc'] stamp. kn: only
    {estimator, estimator_id, gating_rule, min_n}. e3: estimator, parameters, gating rule,
    code provenance (sha256 + git state of `code_files`, default PFUNC_PROVENANCE_FILES),
    catalogue/store provenance; the caller adds the counts (pfunc_counts).
    """
    meta = {**P.estimator_metadata(dedup=dedup, mode=mode, min_n=min_n), "min_n": int(min_n)}
    if mode == "kn":
        return PfuncContext("kn", min_n=min_n, dedup=dedup), meta
    meta["mode"] = mode
    t0 = time.time()
    cat_path = (os.path.abspath(catalogue) if isinstance(catalogue, str)
                else getattr(catalogue, "path", None))
    if isinstance(catalogue, str) and not os.path.exists(catalogue):
        raise _catalogue_error(f"clade catalogue {catalogue} not found", catalogue, pgfam_cache)
    cat = catalogue if isinstance(catalogue, P.CladeCatalogue) else P.load_clade_catalogue(catalogue)
    current = P.catalogue_is_current(cat, pgfam_cache)
    staled_by_run = False
    if not current:
        changed = False
        if preflight and preflight.get("current") and pgfam_cache and os.path.exists(pgfam_cache):
            st = os.stat(pgfam_cache)
            changed = [st.st_size, st.st_mtime_ns] != list(preflight.get("cache_stat") or [])
        if changed:
            staled_by_run = True
            if log:
                log(f"[pfunc] WARNING: this run's selection/prefetch wrote new genomes into "
                    f"{pgfam_cache} after the preflight, so the clade catalogue no longer "
                    f"matches it; those genomes are not clade members in this run. Rebuild the "
                    f"catalogue before the next run:\n    "
                    f"{P.catalogue_rebuild_command(pgfam_cache, cat_path or P.DEFAULT_CATALOGUE)}")
        elif allow_stale:
            if log:
                log(f"[pfunc] WARNING: the clade catalogue does not match {pgfam_cache} "
                    f"(size/mtime); --pfunc-allow-stale-catalogue given, continuing")
        else:
            raise _catalogue_error(f"clade catalogue {cat_path} is STALE for {pgfam_cache}",
                                   cat_path or P.DEFAULT_CATALOGUE, pgfam_cache)
    need, n_scope = set(), 0
    for asv, rec in records.items():
        n = (n_by_asv.get(asv, 0) if n_by_asv is not None
             else int(rec.get("n_selected") or len(rec.get("selected") or [])))
        if n < min_n or not rec.get("selected"):
            continue
        n_scope += 1
        anchors = anchor_candidates(
            rec, built_ids=(built_by_asv or {}).get(asv) if built_by_asv is not None else None)
        need |= P.required_genomes(rec, cat, dedup=dedup, taxdb=taxdb, anchors=anchors)
    pdir = os.path.join(out_dir, PFUNC_DIR)
    os.makedirs(pdir, exist_ok=True)
    store = P.PGFamStore()
    st = {}
    if need:
        st = store.ensure(need, cache_path=pgfam_cache, side_path=os.path.join(pdir, PFUNC_FETCHED),
                          fetch=fetch, workers=fetch_workers, log=log)
    store.save(os.path.join(pdir, PFUNC_STORE))
    n_taxdb = cat.n_taxdb_lookups
    P._TAXDB = None                     # the taxdump is not needed past this point
    meta.update({
        "code_provenance": code_provenance(code_files),
        "pgfam_cache": os.path.abspath(pgfam_cache) if pgfam_cache else None,
        "empty_pgfam_list": "treated as missing (16S-only assembly or stored fetch failure)",
        "catalogue": {"path": cat_path,
                      "current_vs_pgfam_cache": current, "staled_by_this_run": staled_by_run,
                      "stale_allowed": bool(allow_stale) and not current and not staled_by_run,
                      "n_genomes": len(cat.genomes),
                      "source_cache_mtime": (cat.meta or {}).get("source_cache_mtime"),
                      "generated": (cat.meta or {}).get("generated"),
                      "anchors_resolved_with_taxopy": n_taxdb},
        "scope_rule": ("ASVs with n_selected >= min_n get a PGFam store entry before the merge; "
                       "E3 is applied when the BUILT n >= min_n" if n_by_asv is None else
                       "ASVs whose built n (manifest n_genomes) >= min_n"),
        "n_asvs_in_scope": n_scope,
        "pgfam_store": {"file": os.path.join(PFUNC_DIR, PFUNC_STORE),
                        "fetched_file": os.path.join(PFUNC_DIR, PFUNC_FETCHED),
                        "genomes_needed": len(need), "genomes_available": len(store),
                        "genomes_absent": len(store.missing),
                        "absent_examples": sorted(store.missing)[:20], **st},
        "seconds_prepare": round(time.time() - t0, 1),
    })
    if log:
        log(f"[pfunc] E3 store: {len(store)}/{len(need)} genomes for {n_scope} in-scope ASVs "
            f"({len(store.missing)} absent; {st.get('fetched', 0)} fetched) in {time.time() - t0:.0f}s")
    return PfuncContext("e3", min_n=min_n, dedup=dedup, store=store, catalogue=cat), meta


_PFUNC_SUMMARY_KEYS = ("applied", "reason", "estimator", "n_built", "clade_rank", "n_neighbours",
                       "n_neighbours_used", "n_neighbours_missing", "n_neighbours_empty",
                       "weight_missing", "weight_total", "n_members", "n_members_used",
                       "feature_resolution", "n_features_estimated", "n_features_changed",
                       "mean_kn", "mean_e3")


def pfunc_summary(diag):
    """The compact per-ASV record for the manifest entry (full diagnostics: provenance)."""
    return {k: diag[k] for k in _PFUNC_SUMMARY_KEYS if k in diag}


REASON_RECORD_MISSING = "record not in selection audit"
REASON_DRIFT = "inputs drifted since build"


def reason_bucket(reason):
    """Group a not-applied reason for counting: variable texts ('n_built 1 < min_n 2',
    'error: <message>', 'source genome(s) missing ...: [...]') collapse to one key."""
    r = str(reason)
    if r.startswith("n_built"):
        return "n_built < min_n"
    if r.startswith("error:"):
        return "error"
    if r.startswith("source genome(s) missing"):
        return "source genome(s) missing from the genome cache"
    if r.startswith("n mismatch"):
        return "n mismatch"
    return r


def pfunc_counts(entries):
    """Run-level counts over manifest entries' 'pfunc' summaries (for _meta['pfunc']):
    applied, not-applied reasons (bucketed by reason_bucket), and the separately counted
    exceptional cases n_errors ('error: ...', k/n kept), n_inputs_drifted (retro tool),
    n_record_not_in_audit (retro tool)."""
    applied, reasons, ranks, res = 0, Counter(), Counter(), Counter()
    n_err = n_drift = n_norec = 0
    for e in entries:
        pf = (e or {}).get("pfunc")
        if not pf:
            continue
        if pf.get("applied"):
            applied += 1
            ranks[str(pf.get("clade_rank"))] += 1
            res.update(pf.get("feature_resolution") or {})
        else:
            r = str(pf.get("reason"))
            reasons[reason_bucket(r)] += 1
            n_err += r.startswith("error:")
            n_drift += r.startswith(REASON_DRIFT)
            n_norec += r == REASON_RECORD_MISSING
    return {"n_applied": applied, "not_applied_by_reason": dict(reasons),
            "clade_rank_when_applied": dict(ranks), "feature_resolution_when_applied": dict(res),
            "n_errors": n_err, "n_inputs_drifted": n_drift, "n_record_not_in_audit": n_norec}


def _dump_json(obj, path):
    """json.dump with the converter's default formatting, written atomically (a failed
    write leaves neither a partial target nor the temp file)."""
    tmp = f"{path}.tmp.{os.getpid()}"
    try:
        with open(tmp, "w") as f:
            json.dump(obj, f)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


def _snapshot_probabilities(genome_obj):
    """[(dict, had 'probability', value)] for every feature and CDS, to undo a partial E3."""
    return [(x, "probability" in x, x.get("probability"))
            for c in ("features", "cdss") for x in genome_obj.get(c) or ()]


def _restore_probabilities(snap):
    for x, had, v in snap:
        if had:
            x["probability"] = v
        else:
            x.pop("probability", None)


def reweighted_q(genome_obj, kn_by_fid):
    """{fid: new feature probability} for the features reweight_asv replaced."""
    feats = {f.get("id"): f for f in genome_obj.get("features") or ()}
    return {fid: feats[fid]["probability"] for fid in kn_by_fid}


# --------------------------------------------------------------------------- per-ASV
def build_for_asv(asv, rec, *, genome_cache_dir, builder, pfunc=None, digests=None):
    """Runs with cwd == the output dir, so the merger's hardcoded relative side-effect
    dirs (genome_objects/, taxonomies/) land under the output dir.

    digests: a dict to receive {genome_id: pfunc_estimator.source_digest} of the built source
    genomes (computed on the dicts as loaded, before the merge; filled only when the ASV is
    built). The builders collect it into OUT/pfunc/source_digests.json in every mode; it is
    kept out of the entry and the per-ASV files.

    The source genomes are BV-BRC PATRIC annotation only (build_genome -> bvbrc_only, on the
    cache-read and the fresh-build path), so the merge, k/n and the digests never see a RefSeq
    feature.

    pfunc: None or a kn PfuncContext -> the k/n build: the genome, provenance and taxonomy
    files and the manifest entry are byte-identical to the build without E3 (pfunc=None).
    An e3 PfuncContext ->
    after the merge (and BEFORE the core/accessory annotation, so the genome_objects copy keeps
    exactly its keys) E3 replaces k/n when >= pfunc.min_n source genomes were built (prior
    anchored on the first built rep); then the converter's genome_objects/{asv}.json is
    re-dumped, {asv}.json carries the same values, the provenance sidecar gains
    'probability_kn' and 'pfunc' (full diagnostics, incl. the per-source digests), the
    manifest entry a 'pfunc' summary. If E3 raises, the genome keeps k/n (every probability is
    restored; genome_objects/ is left as the converter wrote it) and the entry records
    pfunc {'applied': False, 'reason': 'error: ...'}: an E3 failure never loses a genome."""
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
    # what the retro tool compares against the genome cache later (the dicts as loaded)
    dig = P.source_digests(ok, sources) if digests is not None else None

    conv = KB.LocalGenomeConverter()
    synth = conv.create_synthetic_genome(
        asv, genomes=sources, template_file=None, save_taxonomy=True, taxonomy=None,
        taxonomy_output_dir="taxonomies", genome_output_dir="genome_objects", verbose=False)

    pf_diag, pf_kn = None, {}
    if pfunc is not None and pfunc.e3:
        snap = _snapshot_probabilities(synth)
        try:
            pf_diag, pf_kn = pfunc.reweight(synth, rec, sources, len(ok), built_ids=ok,
                                            source_digests=dig)
            if pf_diag.get("applied"):
                # the converter already wrote genome_objects/{asv}.json with k/n: rewrite it.
                # synth is still exactly that object (annotate_and_summarize has not run), so
                # the copy differs from a k/n build only in feature/CDS probability values.
                _dump_json(synth, os.path.join("genome_objects", f"{asv}.json"))
        except Exception as exc:
            # keep k/n: undo any probability E3 wrote (apply_to_genome validates every fid
            # before writing, but the genome_objects rewrite comes after it); the on-disk
            # genome_objects copy is the converter's k/n one (_dump_json is atomic)
            _restore_probabilities(snap)
            pf_diag = {"estimator": P.ESTIMATOR_NAME, "mode": pfunc.mode, "n_built": len(ok),
                       "min_n": pfunc.min_n, "applied": False,
                       "reason": f"error: {type(exc).__name__}: {str(exc)[:300]}"}
            pf_kn = {}

    prevalence = function_prevalence(sources)
    summary = annotate_and_summarize(synth, prevalence, len(ok))

    # P2: the 16S-indistinguishable same-species organisms each selected rep stands in for
    # (recorded by select_references; the source organism, in a self-recovery run, appears
    # here). These do NOT enter the gene union — they are honest provenance that the marker
    # cannot resolve these genomes apart.
    equivalents = [e for s in rec["selected"] for e in (s.get("equivalent_genomes") or [])]

    prov = {"asv": asv, "genomes": ok, "equivalent_genomes": equivalents,
            "core_functions": sorted(fn for fn, c in prevalence.items()
                                     if c >= max(1, (len(ok) + 1) // 2)),
            "accessory_functions": sorted(fn for fn, c in prevalence.items() if c == 1),
            "function_prevalence": dict(prevalence)}
    if pf_diag is not None and pf_diag.get("applied"):
        prov["probability_kn"] = pf_kn       # {fid: the k/n E3 replaced}
        prov["pfunc"] = pf_diag
    with open(f"{asv}.json", "w") as f:
        json.dump(synth, f)
    with open(os.path.join("provenance", f"{asv}.json"), "w") as f:
        json.dump(prov, f)

    entry = {"status": "built", "tier": rec["tier"], "flags": rec["flags"],
             "consensus_family": rec.get("consensus_family"),
             "n_genomes": len(ok), "genomes": ok, "failed": failed,
             "n_equivalent_genomes": len(equivalents),
             "scientific_name": synth.get("scientific_name"), "taxonomy": synth.get("taxonomy"),
             "n_features": len(synth.get("features", [])), "gene_union": summary,
             "genome_file": f"{asv}.json", "provenance_file": f"provenance/{asv}.json"}
    if pf_diag is not None:
        entry["pfunc"] = pfunc_summary(pf_diag)
    if digests is not None:
        digests.update(dig)
    return entry


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
    add_pfunc_args(ap)
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

    # E3 (the default) needs a current clade catalogue: fail now, not after selection
    preflight = pfunc_preflight(args.pfunc, args.pfunc_catalogue, args.pfunc_pgfam_cache,
                                allow_stale=args.pfunc_allow_stale_catalogue)
    if args.pfunc == "e3" and args.builder == "full":
        print("[pfunc] WARNING: --builder full source genomes are expected to carry no PGFAM "
              "aliases (BVBRCToKBaseConverter keeps PGFams out of the features; inferred from "
              "the converter's code, not run), so E3 would resolve no feature and keep k/n "
              "everywhere (manifest reason 'no feature resolved to a PGFam'); use --builder "
              "fast for E3", flush=True)

    audit, gene_data = get_selection(args, gene_cache=gene_cache, selection_out=selection_out)
    sel_src = selection_out if args.hits else os.path.abspath(args.selection)
    sel_stamp = file_stamp(sel_src)        # the audit as read / written now (every mode)

    os.makedirs(out_dir, exist_ok=True)
    asvs = list(audit)[:args.limit] if args.limit else list(audit)
    pfunc, pfunc_meta = prepare_pfunc(
        {a: audit[a] for a in asvs}, out_dir, mode=args.pfunc, min_n=args.pfunc_min_n,
        dedup=args.pfunc_dedup, pgfam_cache=os.path.abspath(args.pfunc_pgfam_cache),
        catalogue=os.path.abspath(args.pfunc_catalogue), fetch=not args.pfunc_no_fetch,
        fetch_workers=args.pfunc_fetch_workers,
        allow_stale=args.pfunc_allow_stale_catalogue, preflight=preflight)
    os.chdir(out_dir)  # contain the merger's hardcoded relative side-effect dirs here
    for d in ("genome_objects", "taxonomies", "provenance", "ASVset_taxonomies"):
        os.makedirs(d, exist_ok=True)
    print(f"[build] {len(asvs)} ASVs | builder={args.builder} | pfunc={args.pfunc} | out={out_dir}")

    manifest, t0 = {}, time.perf_counter()
    built = abst = failed = 0
    run_digests, dig_conflicts = {}, set()
    for i, asv in enumerate(asvs, 1):
        rec = audit[asv]
        if rec["n_selected"] == 0:
            manifest[asv] = {"status": "abstain", "tier": rec["tier"], "flags": rec["flags"],
                             "midas_taxonomy": rec.get("midas_taxonomy")}
            abst += 1
            continue
        dg = {}
        try:
            entry = build_for_asv(asv, rec, genome_cache_dir=genome_cache_dir, builder=args.builder,
                                  pfunc=pfunc, digests=dg)
        except Exception as exc:
            entry = {"status": "failed", "tier": rec["tier"], "error": str(exc)[:300]}
        if entry["status"] == "built":
            merge_source_digests(run_digests, dg, dig_conflicts)
        manifest[asv] = entry
        built += entry["status"] == "built"
        failed += entry["status"] == "failed"
        if i % 25 == 0 or i == len(asvs):
            print(f"[build] {i}/{len(asvs)} built={built} abstain={abst} failed={failed} "
                  f"({time.perf_counter()-t0:.0f}s)", flush=True)

    if pfunc.e3:
        pfunc_meta.update(pfunc_counts(manifest.values()))
    dig_meta = write_source_digests(out_dir, run_digests, dig_conflicts)
    with open("synthetic_genome_manifest.json", "w") as f:
        json.dump({"_meta": {"selection": sel_src, "selection_stamp": sel_stamp,
                             "gene_data": gene_data, "annotation": ANNOTATION_STAMP,
                             "builder": args.builder, "n_asvs": len(asvs), "built": built,
                             "abstain": abst, "failed": failed, "source_digests": dig_meta,
                             "pfunc": pfunc_meta},
                   **manifest}, f, indent=1)
    print(f"[build] done: built={built} abstain={abst} failed={failed} "
          f"-> {os.path.join(out_dir, 'synthetic_genome_manifest.json')}")


if __name__ == "__main__":
    main()
