#!/usr/bin/env python
"""
pfunc_estimator.py -- production E3 feature probabilities for prFBA synthetic genomes.

The shipped feature probability is k/n: the share of an ASV's n built source genomes that
carry the merged feature's BV-BRC product string (codiffusion_bioreactor/
bvbrc_to_kbase_genome.py:1105-1197). Stage E (region_validation/shared_reference) selected
a better estimator, E3 (`E3|K=100|tau=0.02|genome|family|alpha=2`, no isotonic map), and
Stage O showed it lowers truth-absent flux for multi-genome selections (n >= 2). This module
is the production port of E3, in the PGFam namespace, per the integration spec
(scratchpad INTEGRATION_SPEC.md, 2026-11-11).

Sources (abbreviations used in the docstrings below):
  PE = region_validation/shared_reference/scripts/pfunc_estimators.py (sha256 47d5b6c1...)
       PE:332-361  stage0_lists  (stage-0 filter, identity sort, genome dedup)
       PE:562-591  clade back-off and members (family -> order, MIN_CLADE, pi = mean)
       PE:679-745  compute_Q     (weights, fallback to reps, formula, clamping)
  OT = region_validation/shared_reference/scripts/objective_tuning.py
       OT:244-269  est_role_probabilities (q of a set-valued function from a weights record)

What E3 does, for one ASV:
  1. Neighbours: the selected reps first, then the stage-0 survivors of the (first 100) hits
     sorted by (-identity, -align_score) (stable), one per species (production default;
     the benchmark's genome dedup was a species dedup because its panel held one genome per
     species), K = 100 in total, weighted w_j = exp(-(id_1 - id_j)/tau), tau = 0.02.
  2. Clade prior: the anchor's (the first BUILT selected rep, normally selected[0]) NCBI
     family, else order, with >= 3 member
     species among the PGFam-cached genomes (one genome per species, <= 100 species);
     pi_c(f) = share of members carrying f; none admissible -> no prior.
  3. q_f = (sum_j w_j x_jf + alpha*has*pi_c(f)) / (sum_j w_j + alpha*has), alpha = 2,
     x_jf = 1[PGFams(f) & PGFams(g_j) non-empty]; q >= 1-1e-9 -> 1, q <= 1e-12 -> 0.
  4. Only for n >= 2 built source genomes. At n = 1 the feature probability stays k/n (1.0).

Memory / IO contract. The shared PGFam cache (bvbrc_cache/genome_gene_families.json, ~5 GB,
genome_id -> sorted PGFam list) is NEVER json.load-ed and NEVER rewritten here. PGFamStore
streams it ONCE per batch with ijson (call `required_genomes` for every ASV of a run, union
the ids, then one `PGFamStore.ensure(all_ids, ...)`); missing genomes can be fetched from
BV-BRC into a per-run side file (failures are never stored).

Public API (see each docstring):
  constants K, TAU, ALPHA, MIN_CLADE, MIN_N, HITS_USED, COV_MIN, SWED_MAX, FAMILY_FLOOR,
            MAX_CLADE_SPECIES, CLADE_RANKS
  species_key(species, genome_id)
  neighbours_from_decisions(rec, dedup='species', K=K, tau=TAU, reps_first=True, info=None,
                            hits_used=HITS_USED)
  neighbours_from_hits(hit_record, reps=None, dedup='genome', K=K, tau=TAU,
                       hits_used=HITS_USED, reps_first=True, info=None)
  load_taxdb(nodes_dmp=..., names_dmp=...)
  build_clade_catalogue(cache_path=..., taxdb=None, out_path=..., log=print)
  load_clade_catalogue(path=...) -> CladeCatalogue ; catalogue_is_current(cat, cache_path)
  catalogue_meta(path) ; catalogue_meta_is_current(meta, cache_path)
  catalogue_rebuild_command(cache_path, out_path)   (CLI: `pfunc_estimator.py build-catalogue`)
  clade_members(anchor_gid, catalogue, exclude_species=None, exclude_genomes=None,
                taxdb=None, min_clade=MIN_CLADE, max_species=MAX_CLADE_SPECIES)
  PGFamStore (ensure / fetch_missing / save / load / get)
  prior_anchor(rec, built_ids=None)
  required_genomes(rec, catalogue, dedup='species', K=K, tau=TAU, taxdb=None,
                   hits_used=HITS_USED, anchors=None)
  e3_probabilities(feature_pgfams, neighbours, members, store, alpha=ALPHA, ...)
  resolve_feature_pgfams(genome_obj, sources)
  source_digest(src) ; source_digests(genome_ids, sources)
  apply_to_genome(genome_obj, q_by_fid)
  reweight_asv(genome_obj, rec, sources, store, catalogue, n_built, mode='e3', ...,
               built_ids=None)
  estimator_metadata(dedup='species', mode='e3', min_n=MIN_N)
  probability_description(pfunc_meta=None, entry_pfunc=None)   (metadata text)
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import time
from collections import Counter

import select_references as S

ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PGFAM_CACHE = os.path.join(ROOT, "bvbrc_cache", "genome_gene_families.json")
DEFAULT_CATALOGUE = os.path.join(ROOT, "bvbrc_cache", "e3_clade_catalogue.json")
TAXDUMP_NODES = os.path.join(ROOT, "nodes.dmp")
TAXDUMP_NAMES = os.path.join(ROOT, "names.dmp")

# --------------------------------------------------------------------------- constants
K = 100                 # neighbours kept (PE:136 K3 = 100; selected.json families.E3)
TAU = 0.02              # identity-decay length scale (PE:135 TAU3; selected tau)
ALPHA = 2.0             # clade-prior pseudo-count, in units of neighbour weight (selected alpha)
MIN_CLADE = 3           # admissible clade: >= 3 member species (PE:144)
MIN_N = 2               # E3 replaces k/n only for >= 2 built source genomes (Stage O scope)
HITS_USED = 100         # first 100 hits (PE:129 HITS_USED); d100 records hold exactly 100
MAX_CLADE_SPECIES = 100  # production prior: at most 100 member species (sha1-ordered cap)
CLADE_RANKS = ("family", "order")   # Stage E start = family, back off to order (PE:569-583)
# stage-0 knobs: calibrate_pfunc.STAGE0 == select_references.DEFAULT_KNOBS (PE:128)
COV_MIN = S.DEFAULT_KNOBS["cov_min"]            # 0.90
SWED_MAX = S.DEFAULT_KNOBS["swed_max"]          # 0.02
FAMILY_FLOOR = S.DEFAULT_KNOBS["family_floor"]  # 0.865
Q_ONE = 1 - 1e-9        # PE:743  q >= 1-1e-9 -> 1.0
Q_ZERO = 1e-12          # PE:744  q <= 1e-12  -> 0.0

ESTIMATOR_NAME = "E3"
ESTIMATOR_ID = "E3|K=100|tau=0.02|family->order|alpha=2|PGFam"


# =========================================================================== neighbours
def species_key(species, genome_id) -> str:
    """Production species-dedup key: select_references._norm(species), or 'genome:<id>' when
    the BV-BRC species is blank (decisions carry no taxon_id, so PE:356-357's 'taxon:<id>'
    fallback is not available; spec 'Amendments')."""
    return S._norm(species) or f"genome:{genome_id}"


def _weights(ids, tau):
    """w_j = exp(-(id_1 - id_j)/tau), id_1 = the best identity in the list (PE:706-711;
    PE uses ids[0], which is the maximum because the list is identity-sorted; with reps
    first the maximum is taken explicitly, and the exact-tie reps get w = 1)."""
    if not ids:
        return []
    id1 = max(ids)
    return [math.exp(-(id1 - x) / tau) for x in ids]


def _assemble(reps, survivors, dedup, K, tau, reps_first, info):
    """reps: [(gid, identity|None, species)], survivors: [(gid, identity, species)] already
    sorted. Returns [(gid, w)]. Shared by both neighbour builders."""
    if dedup not in ("species", "genome"):
        raise ValueError(f"dedup must be 'species' or 'genome', not {dedup!r}")

    def key(g, sp):
        return species_key(sp, g) if dedup == "species" else g

    out, seen, seen_g = [], set(), set()
    if reps_first:
        for g, x, sp in reps:
            if g in seen_g:
                continue
            seen_g.add(g)
            seen.add(key(g, sp))
            out.append([g, x])
    for g, x, sp in survivors:
        if len(out) >= K:
            break
        if g in seen_g:          # a rep, or (genome mode) an already-kept genome
            continue
        k = key(g, sp)
        if k in seen:            # species (or genome) already represented
            continue
        seen.add(k)
        seen_g.add(g)
        out.append([g, x])
    out = out[:K]
    if info is not None:
        info["n_reps"] = len(reps)
        info["n_stage0"] = len(survivors)
        info["n_neighbours"] = len(out)
        info["fallback_reps"] = False
    if not out:
        # PE:702-704: no stage-0 survivor in the first 100 hits -> the reps, each w = 1
        fb, seen_fb = [], set()
        for g, _x, _sp in reps:
            if g not in seen_fb:
                seen_fb.add(g)
                fb.append((g, 1.0))
        if info is not None:
            info["fallback_reps"] = bool(fb)
            info["n_neighbours"] = len(fb[:K])
        return fb[:K]
    known = [x for _, x in out if x is not None]
    if not known:
        return [(g, 1.0) for g, _ in out]
    id1 = max(known)
    # a rep whose identity is unknown is treated as an exact tie (w = 1)
    ids = [id1 if x is None else x for _, x in out]
    return list(zip([g for g, _ in out], _weights(ids, tau)))


def neighbours_from_decisions(rec, dedup="species", K=K, tau=TAU, reps_first=True, info=None,
                              hits_used=HITS_USED):
    """E3 neighbour list from a select_references audit record (what build_for_asv has).

    rec['selected'] : reps, anchor first (genome_id, identity, species, ...)
    rec['decisions']: every hit in rank order with identity (4 dp), align_score, coverage
                      (round(min(aligned_len/asv_len, 1), 4), the stage-0 formula),
                      edlib_agreement (round(|SW-edlib|, 4) from raw identities, None passes),
                      genome_id, species, stage.
    Only the first `hits_used` (100) decisions by rank enter, as PE's top[:100] (PE:409) and
    neighbours_from_hits; every d100 audit holds exactly 100, so the cap matters only for an
    audit made at a deeper record depth.
    Stage 0 (PE:332-347, select_references DEFAULT_KNOBS): coverage >= 0.90, agreement
    <= 0.02 (None passes), identity >= 0.865. Decisions store identity rounded to 4 dp, so a
    raw identity in [0.86495, 0.865) would pass on the field alone; such a hit carries
    stage == 'identity_floor' (select_references l.311-313 tests the raw identity), and is
    rejected here too, so the survivor set equals select_references' stage-0 survivors.
    Survivors are sorted by (-identity, -align_score) (stable over rank order, PE:349) and
    deduplicated (dedup='species': species_key; dedup='genome': genome_id, PE:350-355),
    after the reps (reps_first=True: production rule, reps' species/genomes are skipped),
    K in total. Weights: exp(-(id_1 - id_j)/tau). Fallback (no neighbour at all): the reps
    at w = 1 (PE:702-704). Identities are the audit's 4-dp values (raw identities are not in
    the audit); the hit identities are themselves 4-dp in both the benchmark and the d100
    deployments, so this equals neighbours_from_hits(reps first) there (fidelity review).

    Returns [(genome_id, w)]. `info` (dict, optional) receives n_reps, n_stage0,
    n_neighbours, fallback_reps.
    """
    reps = []
    for s in rec.get("selected") or []:
        reps.append((str(s["genome_id"]), s.get("identity"), s.get("species")))
    decs = list(rec.get("decisions") or [])
    # first `hits_used` by rank (stable: decisions are stored in rank order; a missing rank
    # keeps its stored position)
    order = sorted(range(len(decs)), key=lambda i: (decs[i].get("rank")
                                                    if decs[i].get("rank") is not None else i, i))
    decs = [decs[i] for i in order[:hits_used]]
    surv = []
    for d in decs:
        cov = d.get("coverage")
        agree = d.get("edlib_agreement")
        ident = d.get("identity")
        if cov is None or ident is None or cov < COV_MIN:
            continue
        if agree is not None and agree > SWED_MAX:
            continue
        if ident < FAMILY_FLOOR or d.get("stage") == "identity_floor":
            continue
        surv.append(d)
    surv.sort(key=lambda d: (-d["identity"], -d["align_score"]))
    survivors = [(str(d["genome_id"]), float(d["identity"]), d.get("species")) for d in surv]
    return _assemble(reps, survivors, dedup, K, tau, reps_first, info)


def stage0_hits(top, asv_len):
    """select_references stage 0 applied to raw hits exactly as PE:337-348:
    coverage = round(min(aligned_len/L, 1), 4) >= cov_min; agree = round(|identity -
    edlib_identity|, 4) <= swed_max when edlib_identity is not None; identity >=
    family_floor. Returns the survivors sorted by (-identity, -align_score) (stable)."""
    L = asv_len or 1
    surv = []
    for h in top:
        ident = h["identity"]
        cov = round(min(h["aligned_len"] / L, 1.0), 4)
        ei = h.get("edlib_identity")
        agree = None if ei is None else round(abs(ident - ei), 4)
        if cov < COV_MIN:
            continue
        if agree is not None and agree > SWED_MAX:
            continue
        if ident < FAMILY_FLOOR:
            continue
        surv.append(h)
    surv.sort(key=lambda h: (-h["identity"], -h["align_score"]))
    return surv


def neighbours_from_hits(hit_record, reps=None, dedup="genome", K=K, tau=TAU,
                         hits_used=HITS_USED, reps_first=True, info=None):
    """E3 neighbour list from a raw alignment-hit record (align_hits / fold hits format:
    asv_len, top20 = [{rank, genome_id, identity, edlib_identity, aligned_len, align_score,
    lineage{Species}}...]). Used by the panel gate, where reps come from a re-selection.

    Stage 0 exactly as PE:332-361 on top20[:hits_used] (stage0_hits), sort (-identity,
    -align_score) stable, dedup by genome_id (dedup='genome', PE:350-355, the Stage E/O
    configuration) or by species_key(lineage.Species, genome_id) (dedup='species', the
    production rule), K total, w = exp(-(id_1 - id_j)/tau).

    reps: optional list of genome ids or dicts with 'genome_id' (and optionally 'species',
    'identity'). With reps_first=True they lead the list (production rule); with
    reps_first=False the list is PE's pure identity order and the reps serve only as the
    fallback when nothing survives stage 0 (PE:702-704, w = 1). A rep's identity is its best
    identity among this record's hits (else the dict's identity, else treated as a tie);
    its species is the dict's 'species', else its hit lineage Species.
    """
    top = (hit_record.get("top20") or [])
    surv = stage0_hits(top[:hits_used], hit_record.get("asv_len"))
    # a rep's identity: its first (= best) stage-0 survivor copy, else its best hit copy
    best_id, best_sp = {}, {}
    for h in surv:
        g = str(h["genome_id"])
        if g not in best_id:
            best_id[g] = float(h["identity"])
            best_sp[g] = (h.get("lineage") or {}).get("Species")
    for h in top:
        g = str(h["genome_id"])
        if g in best_id and g in best_sp:
            continue
        x = float(h["identity"])
        if g not in best_id or x > best_id[g]:
            best_id[g] = x
        best_sp.setdefault(g, (h.get("lineage") or {}).get("Species"))
    rep_list = []
    for r in reps or []:
        if isinstance(r, dict):
            g = str(r["genome_id"])
            x = best_id.get(g, r.get("identity"))
            sp = r.get("species") if r.get("species") is not None else best_sp.get(g)
        else:
            g = str(r)
            x = best_id.get(g)
            sp = best_sp.get(g)
        rep_list.append((g, None if x is None else float(x), sp))
    survivors = [(str(h["genome_id"]), float(h["identity"]),
                  (h.get("lineage") or {}).get("Species")) for h in surv]
    return _assemble(rep_list, survivors, dedup, K, tau, reps_first, info)


# =========================================================================== taxonomy
def load_taxdb(nodes_dmp=TAXDUMP_NODES, names_dmp=TAXDUMP_NAMES):
    """taxopy TaxDb from the prFBA-root NCBI taxdump, without merged.dmp -- the same source
    and settings the hit lineage uses (edlib_biopython_hits.lineage_for;
    region_validation_d100/scripts/_common._get_taxdb). keep_files=True: never delete the
    dump files. Several GB of RAM and about a minute to load."""
    import taxopy
    return taxopy.TaxDb(nodes_dmp=nodes_dmp, names_dmp=names_dmp, keep_files=True)


def _taxid_of(genome_id):
    try:
        return int(str(genome_id).split(".")[0])
    except (TypeError, ValueError):
        return None


def resolve_taxid(taxid, taxdb):
    """(species_taxid, family_name, order_name) for an NCBI taxid via taxopy; missing ranks
    are None (the production analogue of the panel's 'unclassified_*' placeholders, spec
    'Amendments'); an unknown taxid (or any taxopy failure) -> (None, None, None)."""
    import taxopy
    if taxid is None:
        return (None, None, None)
    try:
        t = taxopy.Taxon(int(taxid), taxdb)
        rn = t.rank_name_dictionary
        rt = t.rank_taxid_dictionary
        sp = rt.get("species")
        return (int(sp) if sp is not None else None, rn.get("family") or None,
                rn.get("order") or None)
    except Exception:
        return (None, None, None)


def build_clade_catalogue(cache_path=DEFAULT_PGFAM_CACHE, taxdb=None, out_path=DEFAULT_CATALOGUE,
                          log=print):
    """Precompute genome -> (species taxid, family, order) for every genome in the shared
    PGFam cache, the universe of clade members for the E3 prior (spec decision 3).

    Streams the cache with ijson.kvitems (one pass; values are only measured for emptiness,
    never kept) -- the cache is never loaded or written. Each genome's taxid is its
    genome_id prefix (as edlib_biopython_hits.py:144); species/family/order come from taxopy
    on the prFBA-root nodes.dmp/names.dmp without merged.dmp (as the hit lineage).
    Genomes whose cached PGFam list is empty (16S-only assemblies or stored failures) are
    excluded from 'genomes' (they would count as members lacking every function) and
    counted in _meta. Genomes with an unknown taxid or no species rank are kept with None
    fields; clade_members never selects them.

    Output JSON (compact): {"_meta": {...}, "genomes": {gid: [species_taxid|null,
    family|null, order|null]}, "taxa": {taxid: [species_taxid|null, family|null,
    order|null]}}; 'taxa' covers every prefix taxid seen, for anchors outside the cache.
    Written atomically (tmp + os.replace). Returns the _meta dict.
    """
    import ijson
    t0 = time.time()
    st = os.stat(cache_path)
    gids, empty = [], []
    with open(cache_path, "rb") as fh:
        for gid, fams in ijson.kvitems(fh, "", use_float=True):
            if fams:
                gids.append(gid)
            else:
                empty.append(gid)
            if log and (len(gids) + len(empty)) % 20000 == 0:
                log(f"[catalogue] streamed {len(gids) + len(empty)} keys "
                    f"({time.time() - t0:.0f}s)")
    t_stream = time.time() - t0
    if log:
        log(f"[catalogue] {len(gids)} non-empty + {len(empty)} empty keys in {t_stream:.0f}s")
    if taxdb is None:
        if log:
            log("[catalogue] loading taxdump ...")
        taxdb = load_taxdb()
    taxa = {}
    for g in list(gids) + list(empty):
        t = _taxid_of(g)
        if t is not None and t not in taxa:
            taxa[t] = resolve_taxid(t, taxdb)
    genomes = {}
    n_bad_id = n_unres = n_nosp = n_nofam = n_noord = 0
    for g in gids:
        t = _taxid_of(g)
        if t is None:
            n_bad_id += 1
            genomes[g] = [None, None, None]
            continue
        sp, fam, order = taxa[t]
        if sp is None and fam is None and order is None:
            n_unres += 1
        if sp is None:
            n_nosp += 1
        if fam is None:
            n_nofam += 1
        if order is None:
            n_noord += 1
        genomes[g] = [sp, fam, order]
    meta = {
        "what": "E3 clade catalogue: PGFam-cached genome -> NCBI species taxid, family, order",
        "source_cache": os.path.abspath(cache_path),
        "source_cache_size": st.st_size,
        "source_cache_mtime": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(st.st_mtime)),
        "taxdump": {"nodes": TAXDUMP_NODES, "names": TAXDUMP_NAMES, "merged": None},
        "n_cache_keys": len(gids) + len(empty),
        "n_genomes": len(genomes),
        "n_empty_pgfam_lists_excluded": len(empty),
        "empty_pgfam_genomes": sorted(empty),
        "n_genome_id_without_numeric_taxid": n_bad_id,
        "n_unresolved_taxid": n_unres,
        "n_without_species_rank": n_nosp,
        "n_without_family": n_nofam,
        "n_without_order": n_noord,
        "n_taxa": len(taxa),
        "seconds_stream": round(t_stream, 1),
        "seconds_total": round(time.time() - t0, 1),
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "generator": "pfunc_estimator.build_clade_catalogue",
    }
    obj = {"_meta": meta, "genomes": genomes, "taxa": {str(k): list(v) for k, v in taxa.items()}}
    tmp = out_path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, separators=(",", ":"))
    os.replace(tmp, out_path)
    if log:
        log(f"[catalogue] wrote {out_path}: {len(genomes)} genomes, {len(taxa)} taxa "
            f"({time.time() - t0:.0f}s)")
    return meta


class CladeCatalogue:
    """In-memory clade catalogue: genomes {gid: (species_taxid, family, order)} (the clade
    member universe = non-empty PGFam-cached genomes) and taxa {taxid: (...)} (anchor
    lookup for genomes outside the cache). Build from build_clade_catalogue's JSON
    (load_clade_catalogue) or directly from dicts (tests). The per-rank species index is
    built lazily on first use."""

    def __init__(self, genomes, taxa=None, meta=None):
        self.genomes = {str(g): tuple(v) for g, v in genomes.items()}
        self.taxa = {int(t): tuple(v) for t, v in (taxa or {}).items()}
        self.meta = meta or {}
        self.path = None
        self._index = None
        self.n_taxdb_lookups = 0

    def _build_index(self):
        # rank -> label -> species_taxid -> sorted genome ids (string order)
        idx = {r: {} for r in CLADE_RANKS}
        for g, (sp, fam, order) in self.genomes.items():
            if sp is None:
                continue
            for r, lab in (("family", fam), ("order", order)):
                if lab is None:
                    continue
                idx[r].setdefault(lab, {}).setdefault(sp, []).append(g)
        self._index = {}
        for r, labs in idx.items():
            self._index[r] = {}
            for lab, sps in labs.items():
                # species in the deterministic cap order: smallest sha1(str(species_taxid))
                order_sp = sorted(sps, key=lambda s: hashlib.sha1(str(s).encode()).hexdigest())
                self._index[r][lab] = [(s, sorted(sps[s])) for s in order_sp]

    def ranks_for(self, genome_id, taxdb=None):
        """(species_taxid, family, order) of a genome: its catalogue entry, else its prefix
        taxid in 'taxa', else a taxopy lookup when `taxdb` is given (a TaxDb, or 'auto' to
        load the prFBA-root taxdump once), else (None, None, None)."""
        g = str(genome_id)
        v = self.genomes.get(g)
        if v is not None:
            return v
        t = _taxid_of(g)
        if t is not None and t in self.taxa:
            return self.taxa[t]
        if taxdb is not None and t is not None:
            if isinstance(taxdb, str) and taxdb == "auto":
                taxdb = _auto_taxdb()
            v = resolve_taxid(t, taxdb)
            self.taxa[t] = v
            self.n_taxdb_lookups += 1
            return v
        return (None, None, None)


_TAXDB = None


def _auto_taxdb():
    global _TAXDB
    if _TAXDB is None:
        _TAXDB = load_taxdb()
    return _TAXDB


def load_clade_catalogue(path=DEFAULT_CATALOGUE) -> CladeCatalogue:
    """Load build_clade_catalogue's JSON (tens of MB) into a CladeCatalogue (its `path`
    attribute records where it came from)."""
    with open(path) as f:
        obj = json.load(f)
    cat = CladeCatalogue(obj["genomes"], obj.get("taxa"), obj.get("_meta"))
    cat.path = os.path.abspath(path)
    return cat


def catalogue_meta(path=DEFAULT_CATALOGUE):
    """The catalogue's _meta only, streamed with ijson (build_clade_catalogue writes _meta
    first, so this reads a few hundred bytes, not the whole file). {} when absent."""
    import ijson
    with open(path, "rb") as fh:
        for k, v in ijson.kvitems(fh, "", use_float=True):
            if k == "_meta":
                return v or {}
    return {}


def catalogue_meta_is_current(meta, cache_path=DEFAULT_PGFAM_CACHE):
    """catalogue_is_current on a catalogue's _meta dict: the size and mtime (to the second,
    local time, as build_clade_catalogue records them) of `cache_path` now equal the ones the
    catalogue was built from."""
    m = meta or {}
    try:
        st = os.stat(cache_path)
    except OSError:
        return False
    mt = time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(st.st_mtime))
    return m.get("source_cache_size") == st.st_size and m.get("source_cache_mtime") == mt


def catalogue_is_current(catalogue, cache_path=DEFAULT_PGFAM_CACHE):
    """True when the catalogue was built from the shared cache as it is now (same size and
    mtime as recorded in _meta). In a stale catalogue, genomes added to the cache since (e.g.
    newly prefetched reps) are not clade members, and anchors among them resolve through
    'taxa' / taxopy instead; the builders therefore refuse a stale catalogue unless
    --pfunc-allow-stale-catalogue. Rebuild with catalogue_rebuild_command()."""
    return catalogue_meta_is_current(_as_catalogue(catalogue).meta, cache_path)


def catalogue_rebuild_command(cache_path=DEFAULT_PGFAM_CACHE, out_path=DEFAULT_CATALOGUE):
    """The shell command that rebuilds the clade catalogue (this module's CLI, wrapped in
    systemd-run as the shared machine requires; measured 2026-11-11: 45 s, 2.7 GB peak RSS)."""
    py = os.path.expanduser("~/Documents/py_venv/bin/python")
    return (f"systemd-run --user --unit=e3-catalogue-$(date +%s) --collect -p MemoryMax=8G "
            f"-p MemorySwapMax=0 -p WorkingDirectory={ROOT} --wait --pipe {py} "
            f"{os.path.join(ROOT, 'pfunc_estimator.py')} "
            f"build-catalogue --cache {os.path.abspath(cache_path)} "
            f"--out {os.path.abspath(out_path)}")


def _as_catalogue(catalogue):
    if isinstance(catalogue, CladeCatalogue):
        return catalogue
    if isinstance(catalogue, dict) and "genomes" in catalogue:
        return CladeCatalogue(catalogue["genomes"], catalogue.get("taxa"), catalogue.get("_meta"))
    return CladeCatalogue(catalogue)


def clade_members(anchor_gid, catalogue, exclude_species=None, exclude_genomes=None,
                  taxdb=None, min_clade=MIN_CLADE, max_species=MAX_CLADE_SPECIES):
    """The anchor's E3 clade: (rank | None, [member genome ids]).

    Back-off (PE:569-583, start at family): family, then order. A rank is skipped when the
    anchor's label is None (production analogue of the panel's placeholders, PE:262-263) or
    fewer than `min_clade` member species remain. Members (spec decision 3, production
    substitute for the ~7k species-dereplicated fold reference): catalogue genomes with the
    anchor's label and a species taxid, minus `exclude_genomes` and every genome of a species
    in `exclude_species` (panel gate: the source genome and its species); then ONE genome per
    species, the smallest genome_id (string order); then at most `max_species` species, those
    with the smallest sha1(str(species_taxid)). Unlike PE (where the members were exactly
    the fold reference), the anchor and the reps are members when their clade contains them.
    No admissible rank -> (None, []), i.e. no prior (alpha dropped).
    """
    cat = _as_catalogue(catalogue)
    if cat._index is None:
        cat._build_index()
    _sp, fam, order = cat.ranks_for(anchor_gid, taxdb=taxdb)
    exs = {int(s) for s in (exclude_species or ()) if s is not None}
    exg = {str(g) for g in (exclude_genomes or ())}
    for rank, lab in (("family", fam), ("order", order)):
        if lab is None:
            continue
        members = []
        for sp, gl in cat._index[rank].get(lab, ()):
            if sp in exs:
                continue
            g = next((x for x in gl if x not in exg), None)
            if g is None:
                continue
            members.append(g)
            if len(members) >= max_species:
                break
        if len(members) >= min_clade:
            return rank, members
    return None, []


# =========================================================================== PGFam store
class PGFamStore:
    """Slim in-memory PGFam sets for the genomes one run needs: {genome_id: frozenset of
    interned PGFam strings} (pattern: region_validation_d100/scripts/genecap_tau0995.py
    Provider / fill_provider). Fill it ONCE per batch with ensure(); do not call ensure()
    per ASV (each call that needs the shared cache streams all ~5 GB of it).

    Absent genomes: get() returns None. Genomes fetched from BV-BRC go to a per-run side
    JSON (never the shared cache; a failed fetch is never stored). save()/load() persist the
    whole store as a slim per-run cache file.
    """

    def __init__(self, sets=None):
        self._sets = {}
        self._intern = {}
        self.missing = set()        # genomes asked for and found nowhere
        self.stats = Counter()
        for g, fams in (sets or {}).items():
            self.add(g, fams)

    # -- core
    def add(self, gid, fams):
        it = self._intern
        self._sets[str(gid)] = frozenset(it.setdefault(f, f) for f in fams)
        self.missing.discard(str(gid))

    def get(self, gid):
        return self._sets.get(str(gid))

    def __contains__(self, gid):
        return str(gid) in self._sets

    def __len__(self):
        return len(self._sets)

    def genome_ids(self):
        return set(self._sets)

    # -- filling
    def ensure(self, gids, cache_path=DEFAULT_PGFAM_CACHE, side_path=None, fetch=False,
               workers=8, log=print):
        """Make every genome in `gids` available if possible: side file first, then ONE
        streaming pass over the shared cache (ijson.kvitems, stops early once every wanted
        genome is found), then (fetch=True) BV-BRC fetch of the rest into `side_path`.
        Returns a stats dict; ids still absent are added to self.missing."""
        want = {str(g) for g in gids} - set(self._sets)
        st = Counter(wanted=len(want))
        if want and side_path and os.path.exists(side_path):
            with open(side_path) as f:
                side = json.load(f)
            for g in list(want):
                if g in side:
                    self.add(g, side[g])
                    want.discard(g)
                    st["from_side_file"] += 1
        if want and cache_path:
            import ijson
            t0 = time.time()
            n = 0
            with open(cache_path, "rb") as fh:
                for g, fams in ijson.kvitems(fh, "", use_float=True):
                    n += 1
                    if g in want:
                        self.add(g, fams)
                        want.discard(g)
                        st["from_shared_cache"] += 1
                        if not want:
                            break
            st["cache_keys_streamed"] = n
            if log:
                log(f"[pgfam] streamed {n} cache keys in {time.time() - t0:.0f}s, "
                    f"{st['from_shared_cache']} wanted genomes found, {len(want)} absent")
        if want and fetch:
            if not side_path:
                raise ValueError("fetch=True needs side_path (never the shared cache)")
            got, err = self.fetch_missing(want, side_path, workers=workers, log=log)
            st["fetched"] = len(got)
            st["fetch_failed"] = err
            want -= set(got)
        self.missing |= want
        st["absent"] = len(want)
        self.stats.update(st)
        return dict(st)

    def fetch_missing(self, gids, side_path, workers=8, log=print):
        """Fetch PGFam sets with prefetch_gene_families.fetch_pgfams (BV-BRC API) into the
        store and merge the successes into `side_path` atomically. An exception is a failure
        and is never stored; an empty list from a successful call is stored (as the shared
        cache does) but counts as unusable in e3_probabilities. Returns (got, n_failed)."""
        from concurrent.futures import ThreadPoolExecutor
        import prefetch_gene_families as PF

        def one(g):
            try:
                return g, PF.fetch_pgfams(g)
            except Exception:
                return g, None
        got, err = {}, 0
        with ThreadPoolExecutor(max(1, int(workers))) as ex:
            for i, (g, fams) in enumerate(ex.map(one, sorted({str(x) for x in gids})), 1):
                if fams is None:
                    err += 1
                else:
                    got[g] = fams
                if log and i % 200 == 0:
                    log(f"[pgfam] fetched {i}/{len(gids)} (failed {err})")
        for g, fams in got.items():
            self.add(g, fams)
        if got:
            side = {}
            if os.path.exists(side_path):
                with open(side_path) as f:
                    side = json.load(f)
            side.update(got)
            _atomic_json(side, side_path)
        return got, err

    # -- persistence
    def save(self, path):
        """Write the slim per-run cache {gid: sorted PGFam list} atomically."""
        _atomic_json({g: sorted(s) for g, s in self._sets.items()}, path)

    @classmethod
    def load(cls, path):
        with open(path) as f:
            return cls(json.load(f))


def _atomic_json(obj, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, separators=(",", ":"))
    os.replace(tmp, path)


def _store_get(store, gid):
    if store is None:
        return None
    return store.get(str(gid))


# =========================================================================== estimator
def e3_probabilities(feature_pgfams, neighbours, members, store, alpha=ALPHA,
                     min_clade=MIN_CLADE, empty_as_missing=True):
    """E3 q for each feature (PE:742-745 semantics; set-valued features as OT:244-269 with
    'any of the feature's PGFams' in place of 'all of the role's ECs', spec amendment):

        q_f = (sum_j w_j 1[S_f & PGF(g_j)] + alpha*has*pi(S_f)) / (sum_j w_j + alpha*has)
        pi(S_f) = share of members m with S_f & PGF(m) non-empty
        q >= 1-1e-9 -> 1.0 ; q <= 1e-12 -> 0.0

    feature_pgfams: {fid: set of PGFam ids}; features with an empty set get no q.
    neighbours: [(gid, w)]; members: [gid] (clade_members output; [] = no clade).
    store: PGFamStore or a plain {gid: set} mapping. Neighbours/members absent from the
    store (or with an empty PGFam set when empty_as_missing, the default: an empty set is a
    16S-only assembly or a stored fetch failure, and counting it would assert that the
    genome lacks every function; 1,802 of the shared cache's 159,753 keys were empty on
    2026-11-11) are excluded and counted. The
    prior needs >= min_clade usable members (the admission rule, re-applied after the
    exclusions), else has = 0 and alpha leaves numerator and denominator (E3 = E1, as PE
    and the Stage O file's alpha = 0 for null clades). If no neighbour is usable and there
    is no prior, no q can be formed: returns {} (the caller keeps k/n) with
    diagnostics['no_estimate'] = True.

    Returns (q_by_fid, diagnostics).
    """
    nb, nb_missing, nb_empty, w_missing = [], 0, 0, 0.0
    for g, w in neighbours:
        s = _store_get(store, g)
        if s is None:
            nb_missing += 1
            w_missing += float(w)
        elif not s and empty_as_missing:
            nb_empty += 1
            w_missing += float(w)
        else:
            nb.append((s, float(w)))
    mem, m_missing, m_empty = [], 0, 0
    for g in members or ():
        s = _store_get(store, g)
        if s is None:
            m_missing += 1
        elif not s and empty_as_missing:
            m_empty += 1
        else:
            mem.append(s)
    has = 1.0 if (mem and len(mem) >= min_clade) else 0.0
    wsum = sum(w for _, w in nb)
    a = float(alpha) * has
    den = wsum + a
    diag = {
        "n_neighbours": len(neighbours), "n_neighbours_used": len(nb),
        "n_neighbours_missing": nb_missing, "n_neighbours_empty": nb_empty,
        "weight_total": round(wsum + w_missing, 6), "weight_missing": round(w_missing, 6),
        "n_members": len(members or ()), "n_members_used": len(mem) if has else 0,
        "n_members_missing": m_missing, "n_members_empty": m_empty,
        "clade_used": bool(has), "clade_dropped_after_missing": bool(members) and not has,
        "alpha_effective": a, "no_estimate": False,
    }
    if den <= 0:
        diag["no_estimate"] = True
        diag["n_features_estimated"] = 0
        return {}, diag
    # carriers: PGFam -> neighbour indices / member indices, restricted to the features' union
    U = set()
    for S_f in feature_pgfams.values():
        if S_f:
            U.update(S_f)
    cn, cm = {}, {}
    for j, (s, _w) in enumerate(nb):
        for p in (U & s if len(U) < len(s) else s & U):
            cn.setdefault(p, []).append(j)
    if has:
        for i, s in enumerate(mem):
            for p in (U & s if len(U) < len(s) else s & U):
                cm.setdefault(p, []).append(i)
    wv = [w for _, w in nb]
    nm = len(mem)
    q = {}
    for fid, S_f in feature_pgfams.items():
        if not S_f:
            continue
        if len(S_f) == 1:
            p = next(iter(S_f))
            num = sum(wv[j] for j in cn.get(p, ()))
            km = len(cm.get(p, ())) if has else 0
        else:
            js, ms = set(), set()
            for p in S_f:
                js.update(cn.get(p, ()))
                if has:
                    ms.update(cm.get(p, ()))
            num = sum(wv[j] for j in js)
            km = len(ms)
        if has:
            num += a * (km / nm)
        x = num / den
        if x >= Q_ONE:
            x = 1.0
        elif x <= Q_ZERO:
            x = 0.0
        q[fid] = x
    diag["n_features_estimated"] = len(q)
    return q, diag


# =========================================================================== genome objects
def _pgfam_aliases(feature):
    out = set()
    for a in feature.get("aliases") or ():
        if isinstance(a, (list, tuple)) and len(a) >= 2 and a[0] == "PGFAM" and a[1]:
            out.add(str(a[1]))
    v = feature.get("pgfam_id")
    if v:
        out.add(str(v))
    return out


def _functions(feature):
    fs = feature.get("functions")
    if fs:
        return [f for f in fs if f]
    f = feature.get("function")
    return [f] if f else []


def resolve_feature_pgfams(genome_obj, sources):
    """Feature -> PGFam set for a merged synthetic genome (spec 'Amendments'):
      (1) own_alias:      the feature's own ['PGFAM', id] alias(es) (copied by the merge from
                          the CDS that introduced the product string);
      (2) product_match:  else the pgfam ids (top-level 'pgfam_id' and PGFAM aliases) of every
                          source-genome CDS whose function string equals the feature's
                          function string (the merge key, functions[0]);
      (3) unresolved:     else no set -- the feature keeps k/n.
    sources: the ASV's built source genome dicts (build_for_asv's `sources`).
    Returns ({fid: set}, counts{own_alias, product_match, unresolved, features}).
    """
    by_product = {}
    for src in sources or ():
        for f in src.get("features") or ():
            pg = _pgfam_aliases(f)
            if not pg:
                continue
            for fn in _functions(f):
                by_product.setdefault(fn, set()).update(pg)
    out, counts = {}, Counter()
    for f in genome_obj.get("features") or ():
        fid = f.get("id")
        counts["features"] += 1
        own = _pgfam_aliases(f)
        if own:
            out[fid] = own
            counts["own_alias"] += 1
            continue
        pm = set()
        for fn in _functions(f):
            pm |= by_product.get(fn, set())
        if pm:
            out[fid] = pm
            counts["product_match"] += 1
        else:
            counts["unresolved"] += 1
    for k in ("own_alias", "product_match", "unresolved", "features"):
        counts.setdefault(k, 0)
    return out, dict(counts)


SOURCE_DIGEST_METHOD = ("sha256 of the JSON list (compact, ensure_ascii=False) of sorted "
                        "[function string, PGFam id or ''] pairs over the genome's features "
                        "(one pair per function x PGFam: the feature's PGFAM alias(es) and "
                        "top-level pgfam_id)")


def source_digest(src):
    """Digest of what E3 reads from one source genome dict: every (function, PGFam) pair
    (SOURCE_DIGEST_METHOD). Feature order does not matter; a changed product string or PGFAM
    alias does. Every build writes it per built source genome to OUT/pfunc/source_digests.json;
    reweight_build_e3.py recomputes it from the genome cache and compares, to catch a cache
    file re-fetched since the build (a changed alias alone leaves function_prevalence
    unchanged)."""
    pairs = []
    for f in src.get("features") or ():
        fams = sorted(_pgfam_aliases(f)) or [""]
        for fn in _functions(f):
            for pg in fams:
                pairs.append([fn, pg])
    pairs.sort()
    blob = json.dumps(pairs, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def source_digests(genome_ids, sources):
    """{genome_id: source_digest} for the built source genomes, in build order."""
    return {str(g): source_digest(s) for g, s in zip(genome_ids, sources)}


def apply_to_genome(genome_obj, q_by_fid):
    """Write q into features[i]['probability'] and the matching cdss '{fid}.CDS'
    'probability' (both are written by the merge, bvbrc_to_kbase_genome.py:1185-1197; cdss
    entries are copies, so both are set). No feature or CDS is added or removed; features
    not in q_by_fid are untouched. Returns the replaced values {fid: old feature
    probability} (k/n, for the provenance sidecar). Raises KeyError for an unknown fid."""
    feats = {f.get("id"): f for f in genome_obj.get("features") or ()}
    cdss = {c.get("id"): c for c in genome_obj.get("cdss") or ()}
    missing = [fid for fid in q_by_fid if fid not in feats]
    if missing:
        raise KeyError(f"{len(missing)} feature ids not in the genome, e.g. {missing[:3]}")
    old = {}
    for fid, q in q_by_fid.items():
        f = feats[fid]
        old[fid] = f.get("probability")
        f["probability"] = float(q)
        c = cdss.get(f"{fid}.CDS")
        if c is not None:
            c["probability"] = float(q)
    return old


# =========================================================================== per ASV
def prior_anchor(rec, built_ids=None):
    """The genome that anchors the clade prior: the first selected rep that was BUILT
    (selected order; selected[0] when it built, which is the usual case). built_ids None
    (unknown) -> selected[0]. A rep that failed to build stays a neighbour (its PGFam set is
    real 16S-hit evidence) but does not anchor the prior. None for a record without reps."""
    reps = [str(s["genome_id"]) for s in (rec or {}).get("selected") or []]
    if not reps:
        return None
    if built_ids is not None:
        b = {str(g) for g in built_ids}
        for g in reps:
            if g in b:
                return g
    return reps[0]


def required_genomes(rec, catalogue, dedup="species", K=K, tau=TAU, taxdb=None,
                     hits_used=HITS_USED, anchors=None):
    """Genome ids reweight_asv will look up for this ASV (neighbours + clade members), for
    one batch PGFamStore.ensure() over a whole run. [] for a record without reps.

    anchors: the genome ids that may anchor the clade prior (prior_anchor); default
    [selected[0]]. A builder that cannot know before the merge which reps will build passes
    every rep (build_synthetic_genomes.anchor_candidates) so the store holds the clade members
    of whichever rep turns out to be the first built one."""
    if not rec or not rec.get("selected"):
        return set()
    nb = neighbours_from_decisions(rec, dedup=dedup, K=K, tau=tau, hits_used=hits_used)
    out = {g for g, _ in nb}
    for a in (anchors if anchors else [str(rec["selected"][0]["genome_id"])]):
        _rank, mem = clade_members(str(a), catalogue, taxdb=taxdb)
        out |= set(mem)
    return out


def reweight_asv(genome_obj, rec, sources, store, catalogue, n_built, mode="e3", min_n=MIN_N,
                 dedup="species", K=K, tau=TAU, alpha=ALPHA, taxdb=None, built_ids=None,
                 hits_used=HITS_USED, digests=None):
    """Replace k/n by E3 on one built synthetic genome, in place (spec decisions 1-4).

    No-op (returns ({'applied': False, 'reason': ...}, {})) when mode == 'kn', when
    n_built < min_n (n = successfully built source genomes, the n of k/n; Stage O's gain is
    established only at n >= 2), when the record has no reps, when no feature resolves to a
    PGFam ('no feature resolved to a PGFam'), or when no estimate can be formed ('no
    estimate (no usable neighbour, no clade)'). Otherwise: neighbours_from_decisions(rec,
    dedup) -> clade_members(anchor = prior_anchor(rec, built_ids)) ->
    resolve_feature_pgfams(genome_obj, sources) -> e3_probabilities -> apply_to_genome.
    Features without a PGFam keep k/n.

    built_ids: the genome ids of `sources` (the reps that built, in build order). The reps
    are neighbours from the PGFam store whether or not they built (a rep that failed to
    build is still a 16S hit with a real PGFam set); the clade prior is anchored on the
    first BUILT rep. Diagnostics record 'reps_not_built' (selected reps absent from
    built_ids) and 'anchor_built' (selected[0] built; None when built_ids is not given), and
    'source_digests' ({gid: source_digest} of the built sources; the `digests` argument
    when the caller already computed them on the dicts as loaded, else computed here). The
    builders also write them, in every mode, to OUT/pfunc/source_digests.json, which is what
    reweight_build_e3.py checks for genome-cache drift.

    Returns (diagnostics, kn_by_fid): kn_by_fid = {fid: replaced k/n} for the provenance
    sidecar (key 'probability_kn'); diagnostics records the estimator, the parameters, the
    resolution counts, the neighbour/member coverage and the clade rank.
    """
    if mode not in ("e3", "kn"):
        raise ValueError(f"mode must be 'e3' or 'kn', not {mode!r}")
    base = {"estimator": ESTIMATOR_NAME if mode == "e3" else "kn", "mode": mode,
            "n_built": int(n_built), "min_n": int(min_n)}
    if mode == "kn":
        return {**base, "applied": False, "reason": "mode kn"}, {}
    if n_built < min_n:
        return {**base, "applied": False, "reason": f"n_built {n_built} < min_n {min_n}"}, {}
    if not rec or not rec.get("selected"):
        return {**base, "applied": False, "reason": "no selected reps"}, {}
    ninfo = {}
    nb = neighbours_from_decisions(rec, dedup=dedup, K=K, tau=tau, info=ninfo,
                                   hits_used=hits_used)
    reps = [str(s["genome_id"]) for s in rec["selected"]]
    if built_ids is not None:
        built = [str(g) for g in built_ids]
        bset = set(built)
        reps_not_built = [g for g in reps if g not in bset]
        anchor_built = reps[0] in bset
    else:
        built = [str((s or {}).get("id")) for s in sources or ()]
        reps_not_built, anchor_built = None, None
    anchor = prior_anchor(rec, built_ids)
    cat = _as_catalogue(catalogue)
    rank, mem = clade_members(anchor, cat, taxdb=taxdb)
    a_sp, a_fam, a_ord = cat.ranks_for(anchor)
    a_tid = _taxid_of(anchor)
    fpg, counts = resolve_feature_pgfams(genome_obj, sources)
    q, ediag = e3_probabilities(fpg, nb, mem, store, alpha=alpha)
    diag = {**base, "dedup": dedup, "K": K, "tau": tau, "alpha": alpha,
            "anchor": anchor, "anchor_selected": reps[0], "anchor_built": anchor_built,
            "reps_not_built": reps_not_built,
            "anchor_in_catalogue": anchor in cat.genomes,
            "anchor_taxid_in_catalogue_taxa": a_tid in cat.taxa,
            "anchor_family": a_fam, "anchor_order": a_ord,
            "clade_rank": rank if ediag["clade_used"] else None,
            "clade_rank_before_store": rank, "neighbour_build": ninfo,
            "feature_resolution": counts, **ediag,
            "source_digests": (dict(digests) if digests is not None
                               else source_digests(built, sources or ()))}
    if not fpg:
        return {**diag, "applied": False, "reason": "no feature resolved to a PGFam"}, {}
    if not q:
        return {**diag, "applied": False, "reason": "no estimate (no usable neighbour, no clade)"}, {}
    # A merged feature is carried by at least one BUILT source genome, so its probability must
    # stay > 0. An E3 of exactly 0 means no neighbour and no clade member in the PGFam store
    # carries its family, i.e. the genome-dict cache (which resolved the feature's PGFam) and the
    # PGFam cache disagree about the carrier (e.g. an empty or differing PGFam list). Such a
    # feature keeps k/n (found on 67 features of 6 genomes of the two d100 deployments, 2026-11-12).
    zero = sorted(f for f, v in q.items() if v <= 0.0)
    if zero:
        q = {f: v for f, v in q.items() if v > 0.0}
    diag["n_zero_kept_kn"] = len(zero)
    if not q:
        return {**diag, "applied": False, "reason": "every estimate was 0 (kept k/n)"}, {}
    kn = apply_to_genome(genome_obj, q)
    changed = sum(1 for f, v in kn.items() if v is None or abs(float(v) - q[f]) > 1e-12)
    diag.update(applied=True, reason=None, n_features_changed=changed,
                mean_kn=round(sum(float(v or 0) for v in kn.values()) / len(kn), 6),
                mean_e3=round(sum(q.values()) / len(q), 6))
    return diag, kn


def estimator_metadata(dedup="species", mode="e3", min_n=MIN_N):
    """The dict a build stamps into its manifest _meta (spec decision 5). mode 'kn' gives
    only {estimator, estimator_id, gating_rule}: a k/n build carries no E3 parameters."""
    if mode != "e3":
        return {"estimator": "kn", "estimator_id": "k/n",
                "gating_rule": "k/n everywhere (--pfunc kn)"}
    return {
        "estimator": ESTIMATOR_NAME,
        "estimator_id": ESTIMATOR_ID,
        "gating_rule": (f"E3 replaces k/n only when >= {min_n} source genomes were built; "
                        f"otherwise (and for features without a PGFam, or whose estimate is 0 although a "
                        f"built source genome carries them) probability = k/n"),
        "K": K, "tau": TAU, "alpha": ALPHA, "min_clade_species": MIN_CLADE,
        "clade_ranks": list(CLADE_RANKS), "max_clade_species": MAX_CLADE_SPECIES,
        "hits_used": HITS_USED,
        "neighbour_dedup": dedup,
        "prior_anchor": "first BUILT selected rep (selected[0] unless it failed to build); "
                        "reps that failed to build stay neighbours",
        "stage0": {"cov_min": COV_MIN, "swed_max": SWED_MAX, "family_floor": FAMILY_FLOOR},
        "function_namespace": "PGFam (feature PGFAM alias, else source-CDS product match)",
        "clamp": {"to_one_at": Q_ONE, "to_zero_at": Q_ZERO},
        "source_digest": SOURCE_DIGEST_METHOD,
        "source": "region_validation/shared_reference/scripts/pfunc_estimators.py "
                  "(Stage E, selected E3) ; objective_tuning.py (Stage O)",
    }


# =========================================================================== descriptions
# (added by the build integration: one source for the text the metadata writers print)
KN_DESCRIPTION = "fraction of source genomes carrying the function"

# Which BV-BRC annotation a build's feature functions come from. A build stamps its manifest
# _meta['annotation'] with ANNOTATION_STAMP (build_synthetic_genomes.ANNOTATION_STAMP is this
# string): since 2026-11-13 every source genome is BV-BRC's own PATRIC (RASTtk) annotation only.
# A manifest without it is a build made at or before git 3ce12d0, which also merged the RefSeq
# (NCBI PGAP) CDSs that BV-BRC imports.
ANNOTATION_STAMP = "BV-BRC PATRIC only (eq(annotation,PATRIC); cached dicts filtered on read)"
FUNCTIONS_BVBRC = ("BV-BRC (PATRIC) annotations: the RASTtk function strings of BV-BRC's own "
                   "annotation of each source genome (BV-BRC's imported RefSeq/PGAP features "
                   "excluded)")
FUNCTIONS_MIXED = ("BV-BRC genome_feature CDS products of both annotations BV-BRC serves, its own "
                   "PATRIC (RASTtk) annotation and the imported RefSeq (PGAP) one (a build made at "
                   "or before prFBA 3ce12d0, before the BV-BRC-only rule)")


def is_bvbrc_only_build(man_meta):
    """True when a build manifest's _meta records the BV-BRC PATRIC-only annotation stamp."""
    return (man_meta or {}).get("annotation") == ANNOTATION_STAMP


def functions_description(man_meta=None):
    """What a build's feature functions are, for metadata text (upload_kbase_genomes,
    write_mapping_provenance): FUNCTIONS_BVBRC for a build stamped BV-BRC PATRIC-only, else
    FUNCTIONS_MIXED (every build made before the stamp merged both annotations)."""
    return FUNCTIONS_BVBRC if is_bvbrc_only_build(man_meta) else FUNCTIONS_MIXED


def probability_description(pfunc_meta=None, entry_pfunc=None):
    """What Feature.probability means in a build, for metadata text (upload_kbase_genomes,
    write_mapping_provenance). `pfunc_meta` is the manifest `_meta['pfunc']` stamp; None (a
    build made before E3) or estimator 'kn' gives KN_DESCRIPTION, the text those writers
    printed before E3, unchanged. `entry_pfunc` (optional) is one ASV's manifest 'pfunc'
    summary, to describe that genome specifically (applied or not, and why)."""
    if not pfunc_meta or pfunc_meta.get("estimator") != ESTIMATOR_NAME:
        return KN_DESCRIPTION
    min_n = pfunc_meta.get("min_n", MIN_N)
    e3 = (f"E3 estimate ({pfunc_meta.get('estimator_id', ESTIMATOR_ID)}, neighbour dedup "
          f"{pfunc_meta.get('neighbour_dedup')}): the identity-weighted share of the ASV's 16S "
          f"neighbours (selected references, then stage-0 hits; w = exp(-(id1 - id)/"
          f"{pfunc_meta.get('tau', TAU)})) whose BV-BRC PGFam set contains the feature's PGFam, "
          f"shrunk toward the anchor's NCBI family (else order) PGFam prevalence "
          f"(alpha = {pfunc_meta.get('alpha', ALPHA)})")
    kn = f"{KN_DESCRIPTION} (k/n)"
    if entry_pfunc is None:
        return (f"{e3}, for genomes built from >= {min_n} source genomes; features without a "
                f"PGFam, and genomes built from fewer source genomes, keep the {kn}; replaced "
                f"k/n values are kept in provenance/{{asv}}.json 'probability_kn'")
    if entry_pfunc.get("applied"):
        return (f"{e3}, for the {entry_pfunc.get('n_features_estimated')} features with a PGFam; "
                f"the other features keep the {kn}; replaced k/n values are kept in "
                f"provenance 'probability_kn'")
    return f"{kn}; E3 not applied ({entry_pfunc.get('reason')})"


# =========================================================================== CLI
def main(argv=None):
    """`python pfunc_estimator.py build-catalogue [--cache P] [--out P]`: (re)build the E3
    clade catalogue from the shared PGFam cache (streamed, never loaded or written) and the
    prFBA-root taxdump. ~45 s and ~2.7 GB RSS (2026-11-11); run it under systemd-run on the
    shared machine (catalogue_rebuild_command() prints the full command)."""
    import argparse
    ap = argparse.ArgumentParser(description=main.__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build-catalogue", help="rebuild the E3 clade catalogue")
    b.add_argument("--cache", default=DEFAULT_PGFAM_CACHE, help="shared PGFam cache (read only)")
    b.add_argument("--out", default=DEFAULT_CATALOGUE, help="catalogue JSON to write")
    sub.add_parser("rebuild-command", help="print the systemd-run rebuild command")
    args = ap.parse_args(argv)
    if args.cmd == "rebuild-command":
        print(catalogue_rebuild_command())
        return 0
    if os.path.abspath(args.out) == os.path.abspath(args.cache):
        raise SystemExit("--out must not be the PGFam cache")
    meta = build_clade_catalogue(cache_path=args.cache, out_path=args.out)
    print(json.dumps({k: v for k, v in meta.items() if k != "empty_pgfam_genomes"}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
