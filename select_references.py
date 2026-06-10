#!/usr/bin/env python
"""
select_references.py — choose a representative HANDFUL of BV-BRC reference genomes per
16S ASV from its top-20 Smith-Waterman hits, and emit a JSON audit that justifies every
inclusion/exclusion decision.

This is the stage between the alignment hits (`asv_top20_alignment_hits.json`, produced by
align_hits.py / edlib_biopython_hits.py / gpusw) and the synthetic-genome merge
(`bvbrc_to_kbase_genome.py create_synthetic_genome --genomes ...`): the selected
genome_ids are exactly what feed `--genomes`.

Selection is a function of alignment identity/score AND contamination control, via an
ordered pipeline whose ORDER is the load-bearing design choice:

  0  reliability gate     coverage + |SW-edlib| agreement   (drop partial/chimeric)
  1  family floor + tier  best identity vs 16S rank thresholds (abstain below family)
  2  identity band        relative tie-cloud, on IDENTITY not score
  3  family-consensus     drop taxonomically-inconsistent hits (incl. blank-family rule)  <-- before dedup
  4  genome_id dedup      collapse multi-copy 16S operons (best copy)
  5  species dedup        <=1 representative per species
  6  identity-decay rank  order remaining candidates by a decayed representativeness score
  7  marginal-gain stop   greedily add a genome only while it adds enough NOVEL genes and
                          stays within a contamination budget (the data-driven cap)
  8  provenance tagging   anchor vs congener, core vs accessory, confidence per genome

Strategies implemented (1-4 from the design discussion):
  (1) marginal pangenome-gain stopping rule           -> stage 7
  (2) identity-decay weighting                         -> stage 6/7 representativeness
  (3) core/accessory provenance tagging                -> stage 8 + per-gene map helper
  (4) tighter blank-family ("uncultured" lineage) rule -> stage 3

Gene-level strategies (1 & 3) need each genome's gene set. Provide a `gene_provider`
(genome_id -> set of orthologous gene ids, e.g. BV-BRC PGFam/PLFam) for exact novelty;
without one, a taxonomy+identity NOVELTY ESTIMATOR is used and every record is tagged
`gene_data: "estimated"` so the provenance stays honest. See `bvbrc_gene_provider`.

    python select_references.py \
        --hits  /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_top20_alignment_hits.json \
        --out   /home/freiburger/Documents/EmilyKin/bvbrc_alignment_hits/asv_reference_selection.json

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import json
import os
import re
import time
from collections import Counter

# --------------------------------------------------------------------------- knobs
DEFAULT_KNOBS = dict(
    # stage 0 reliability
    cov_min=0.90,            # min aligned_len/asv_len
    swed_max=0.02,           # max |SW identity - edlib identity|
    # stage 1 tiers (16S rank thresholds, Yarza 2014 / Stackebrandt)
    t_species=0.987, t_genus=0.945, t_family=0.865,
    family_floor=0.865,      # below this -> abstain (no genome-level model)
    # stage 2 band
    band=0.005,              # keep hits within this identity of the best
    # stage 3 consensus
    consensus_rank="Family", # genus is blank for ~59% of hits; family is stable
    override_margin=0.03,    # an off-consensus hit may override only if this much better AND in MiDAS
    blank_family_tol=0.001,  # keep a blank-family hit only if within this id of a named consensus hit
    # stage 6 representativeness (rep score weights + identity-decay)
    w_id=0.55, w_score=0.30, w_cov=0.15,
    uncult_factor=0.98,      # gentle sub-tie nudge away from uncultured (NOT a gate)
    decay_tau=0.01,          # identity-decay length scale (1% identity) for weighting
    # stage 7 marginal-gain stopping rule + contamination budget (the data-driven cap)
    novel_min=0.15,          # add a genome only if it contributes >= this fraction of NEW genes
    max_low_conf=1,          # at most this many uncultured/MAG genomes in the handful
    hard_cap=6,              # absolute safety ceiling across all tiers
    # confidence-aware ceilings: the marginal-gain rule is primary, but the LESS certain the
    # ASV's identity, the tighter the speculative-breadth bound (a 90% family-level ASV must
    # not accrue 6 genomes' worth of guessed metabolism).
    cap_species=2, cap_genus=5, cap_family=2,
    contam_species=0.25, contam_genus=0.40, contam_family=0.20,  # off-target gene budget / anchor, by tier
    # novelty ESTIMATOR (used only when no gene_provider): fraction of genes SHARED by rank
    ov_species=0.90, ov_genus=0.75, ov_family=0.55, ov_other=0.45,
    accum_decay=0.12,        # each prior inclusion lowers a new genome's estimated novelty
    default_gene_count=3000, # assumed CDS count when unknown (BV-BRC num_cdss can fill this)
)

# what counts as an unculturable / unusable "organism" name (lineage often blank)
_UNCULT = re.compile(r"uncultured|metagenome|environmental sample|unidentified|unclassified|^\s*$", re.I)
# taxon-name prefixes to strip before any lineage/MiDAS comparison
_PREFIX = re.compile(r"^(Candidatus\s+|Ca[._]\s*|midas_[a-z]_)\s*", re.I)
# a real binomial "Genus species" (not 'sp.' / 'bacterium')
_BINOMIAL = re.compile(r"^[A-Z][a-z]+ [a-z]{3,}$")
# MAG / metagenome-bin organism signals
_MAG = re.compile(r"\b(MAG|metagenome|bin[_.\s-]?\d|UBA\d|SRR\d|sludge|enrichment)\b", re.I)


# --------------------------------------------------------------------------- helpers
def _norm(s) -> str:
    if not s:
        return ""
    return _PREFIX.sub("", str(s).strip()).casefold()


def _is_uncult(org) -> bool:
    return bool(_UNCULT.search(org or ""))


def _is_mag(genome_id, org) -> bool:
    # numeric "taxid.assembly" => an isolate/RefSeq assembly; non-numeric prefix => bin/MAG
    gid = str(genome_id)
    if not re.match(r"^\d+\.\d+$", gid):
        return True
    return bool(_MAG.search(org or ""))


def _midas_set(m) -> set:
    return {_norm(t) for t in (m or "").split()}


def _tier(bi, k) -> str:
    if bi >= k["t_species"]:
        return "species"
    if bi >= k["t_genus"]:
        return "genus"
    if bi >= k["t_family"]:
        return "family"
    return "below"


def _rep_score(h, asv_len, k) -> float:
    cov = min(1.0, h["aligned_len"] / asv_len)
    nscore = h["align_score"] / (2 * asv_len)
    s = k["w_id"] * h["identity"] + k["w_score"] * nscore + k["w_cov"] * cov
    if _is_uncult(h.get("organism")):
        s *= k["uncult_factor"]
    return s


def _weighted_rep(h, best_id, asv_len, k) -> float:
    """rep score with identity-decay relative to the ASV's best hit (strategy 2)."""
    import math

    decay = math.exp(-(best_id - h["identity"]) / max(k["decay_tau"], 1e-9))
    return _rep_score(h, asv_len, k) * decay


def _confidence(rec, cons, midas) -> str:
    if rec["uncultured"] or not _norm(rec["family"]):
        return "low"
    fam_ok = (cons is None) or (_norm(rec["family"]) == cons) or (_norm(rec["family"]) in midas)
    sp = (rec["species"] or "").strip()
    named = bool(_BINOMIAL.match(sp)) and "bacterium" not in sp.lower()
    if named and fam_ok:
        return "high"
    if fam_ok:
        return "medium"
    return "low"


def _unnamed(sp) -> bool:
    """A strain with no real species epithet ('Genus sp. XYZ' / '... bacterium' / blank)."""
    s = (sp or "").lower()
    return (not s) or ("sp." in s) or ("bacterium" in s)


def _estimate_novelty(rec, selected, k) -> float:
    """Estimated fraction of NEW genes a candidate adds vs the already-selected set,
    from the closest taxonomic relationship + an accumulation decay. Placeholder for an
    exact gene-set Jaccard; only used when no gene_provider is supplied."""
    cf, cg, cs = _norm(rec["family"]), _norm(rec["genus"]), _norm(rec["species"])
    cun = _unnamed(rec["species"])
    best_overlap = 0.0
    for r in selected:
        if cs and cs == _norm(r["species"]):
            ov = k["ov_species"]
        elif cg and cg == _norm(r["genus"]):
            # two UNNAMED strains of the same genus are likely conspecific (near-identical
            # gene content) -> treat as species-level overlap so they don't both get in
            ov = k["ov_species"] if (cun and _unnamed(r["species"])) else k["ov_genus"]
        elif cf and cf == _norm(r["family"]):
            ov = k["ov_family"]
        else:
            ov = k["ov_other"]
        best_overlap = max(best_overlap, ov)
    return (1.0 - best_overlap) * max(0.0, 1.0 - k["accum_decay"] * len(selected))


# --------------------------------------------------------------------------- core
def select_representatives(asv_record, *, knobs=None, gene_provider=None, gene_count_fn=None):
    """Reduce one ASV's top-20 hits to a representative handful, tracking every decision.

    Parameters
    ----------
    asv_record : dict with keys asv_len, best_identity, midas_taxonomy, rel_ab, top20[...]
    knobs      : overrides for DEFAULT_KNOBS (every threshold is tunable)
    gene_provider : optional genome_id -> set(orthologous gene ids). If given, novelty and
                    the gene union are computed exactly (strategies 1 & 3). Else estimated.
    gene_count_fn : optional genome_id -> int CDS count (e.g. BV-BRC num_cdss); else default.

    Returns a dict: handful + per-hit decisions + provenance + flags (see module docstring).
    """
    k = {**DEFAULT_KNOBS, **(knobs or {})}
    L = asv_record["asv_len"]
    bi = asv_record["best_identity"]
    midas = _midas_set(asv_record.get("midas_taxonomy"))
    gene_data = "bvbrc" if gene_provider else "estimated"

    # base per-hit records (one per top-20 hit; disposition filled as it flows)
    recs = []
    score1_gid = None
    score1 = -1
    for h in asv_record["top20"]:
        lin = h.get("lineage") or {}
        ei = h.get("edlib_identity")
        agree = None if ei is None else round(abs(h["identity"] - ei), 4)
        rec = dict(
            rank=h["rank"], genome_id=h["genome_id"], taxon_id=h.get("taxon_id"),
            organism=h.get("organism"), identity=round(h["identity"], 4),
            align_score=h["align_score"], coverage=round(min(h["aligned_len"] / L, 1.0), 4),
            edlib_agreement=agree, family=lin.get("Family"), genus=lin.get("Genus"),
            species=lin.get("Species"), uncultured=_is_uncult(h.get("organism")),
            disposition=None, stage=None, reason=None, _h=h,
        )
        recs.append(rec)
        if h["align_score"] > score1:
            score1, score1_gid = h["align_score"], h["genome_id"]

    def excl(rec, stage, reason):
        rec["disposition"], rec["stage"], rec["reason"] = "excluded", stage, reason

    flags = []

    def finish(handful, tier, cons):
        # mark any rec never given a disposition as excluded (shouldn't happen)
        for r in recs:
            if r["disposition"] is None:
                excl(r, "not_reached", "not evaluated (handful already complete)")
        # review flag: was the score-#1 hit excluded by a taxonomic/reliability gate?
        s1 = next((r for r in recs if r["genome_id"] == score1_gid), None)
        if s1 and s1["disposition"] == "excluded" and s1["stage"] in ("family_consensus", "reliability_gate"):
            flags.append(f"score_top_dropped:{s1['stage']}")
        if handful and all(r["uncultured"] for r in handful):
            flags.append("uncultured_only")
        if tier == "family" and handful:
            flags.append("family_level_draft")
        prov = _provenance(handful, cons, gene_data)
        return dict(
            asv_len=L, best_identity=round(bi, 4),
            midas_taxonomy=asv_record.get("midas_taxonomy"),
            rel_ab=asv_record.get("rel_ab"),
            tier=tier, consensus_family=cons, n_selected=len(handful), flags=flags,
            selected=[_clean_selected(r) for r in handful],
            provenance=prov,
            decisions=[_clean_decision(r) for r in recs],
        )

    tier = _tier(bi, k)

    # ---- stage 1: tier floor (abstain below family) ----
    if tier == "below":
        for r in recs:
            excl(r, "tier_floor",
                 f"ASV best identity {bi:.4f} < family threshold {k['t_family']}: too divergent "
                 f"for a genome-level model (would inject false metabolism)")
        flags.append("abstain_below_family")
        return finish([], "below", None)

    # ---- stage 0: reliability gate + per-hit family floor ----
    survivors = []
    for r in recs:
        if r["coverage"] < k["cov_min"]:
            excl(r, "reliability_gate",
                 f"coverage {r['coverage']:.3f} < {k['cov_min']} (partial/local alignment)")
        elif r["edlib_agreement"] is not None and r["edlib_agreement"] > k["swed_max"]:
            excl(r, "reliability_gate",
                 f"|SW-edlib| {r['edlib_agreement']:.3f} > {k['swed_max']} (method-discordant / possible chimera)")
        elif r["_h"]["identity"] < k["family_floor"]:
            excl(r, "identity_floor",
                 f"identity {r['_h']['identity']:.4f} < family floor {k['family_floor']}")
        else:
            survivors.append(r)
    if not survivors:
        flags.append("abstain_reliability")
        return finish([], tier, None)

    # ---- stage 2: relative IDENTITY band ----
    bmax = max(r["_h"]["identity"] for r in survivors)
    in_band = []
    for r in survivors:
        if bmax - r["_h"]["identity"] <= k["band"]:
            in_band.append(r)
        else:
            excl(r, "identity_band",
                 f"identity {r['_h']['identity']:.4f} is >{k['band']} below the best {bmax:.4f} "
                 f"(outside the tie cloud)")

    # ---- stage 3: family-consensus + blank-family rule (strategy 4) ----
    named_fams = [_norm(r["family"]) for r in in_band if _norm(r["family"])]
    cons = Counter(named_fams).most_common(1)[0][0] if named_fams else None
    named_consistent_ids = [
        r["_h"]["identity"] for r in in_band
        if _norm(r["family"]) and (_norm(r["family"]) == cons or _norm(r["family"]) in midas)
    ]
    best_consistent = max(named_consistent_ids) if named_consistent_ids else None
    consistent = []
    for r in in_band:
        f = _norm(r["family"])
        if f == "":
            # blank/unclassified family: keep only if sole in-band, or essentially tied
            # with a named consensus-consistent hit (likely the same organism, unannotated)
            if len(in_band) == 1:
                consistent.append(r)
            elif best_consistent is not None and abs(r["_h"]["identity"] - best_consistent) <= k["blank_family_tol"]:
                consistent.append(r)
            else:
                excl(r, "family_consensus",
                     f"unclassified family and not within {k['blank_family_tol']} identity of a named "
                     f"'{cons}' hit (probable silent off-target)")
        elif cons is None or f == cons or f in midas:
            consistent.append(r)
        elif (r["_h"]["identity"] - bmax) >= k["override_margin"] and f in midas:
            consistent.append(r)
        else:
            excl(r, "family_consensus",
                 f"family '{r['family']}' != band consensus '{cons}' and not in MiDAS prior "
                 f"(taxonomically inconsistent off-target)")
    if not consistent:
        consistent = in_band

    # ---- stage 4: genome_id dedup (multi-copy 16S) ----
    by_gid = {}
    for r in sorted(consistent, key=lambda x: (-x["_h"]["identity"], -x["_h"]["align_score"])):
        g = r["genome_id"]
        if g in by_gid:
            excl(r, "genome_dedup",
                 f"multi-copy 16S: same genome {g} as rank {by_gid[g]['rank']} (best copy kept)")
        else:
            by_gid[g] = r

    # ---- stage 5: species dedup ----
    by_sp = {}
    for r in sorted(by_gid.values(), key=lambda x: -_rep_score(x["_h"], L, k)):
        sp = _norm(r["species"]) or r["genome_id"]
        if sp in by_sp:
            excl(r, "species_dedup",
                 f"redundant species '{r['species']}' (representative kept: {by_sp[sp]['genome_id']})")
        else:
            by_sp[sp] = r

    # ---- stage 6: order candidates by identity-decayed representativeness ----
    candidates = sorted(by_sp.values(), key=lambda x: -_weighted_rep(x["_h"], bi, L, k))

    # ---- stage 7: marginal-gain greedy inclusion + contamination budget (strategies 1 & 2) ----
    tier_cap = min(k["hard_cap"], {"species": k["cap_species"], "genus": k["cap_genus"],
                                   "family": k["cap_family"]}[tier])
    tier_budget = {"species": k["contam_species"], "genus": k["contam_genus"],
                   "family": k["contam_family"]}[tier]
    handful, union, anchor_genes, contam_genes, low_conf = [], set(), 0, 0.0, 0
    for r in candidates:
        h = r["_h"]
        if gene_provider:
            gset = gene_provider(r["genome_id"])
            if not gset:  # 16S-only entry / no assembly in BV-BRC -> unusable for a gene union
                excl(r, "no_gene_set",
                     "no usable CDS set in BV-BRC for this genome_id (16S-only entry / no assembly)")
                continue
            gc = len(gset)
        else:
            gset = None
            gc = (gene_count_fn(r["genome_id"]) if gene_count_fn else None) or k["default_gene_count"]
        is_low = r["uncultured"] or _is_mag(r["genome_id"], r["organism"])
        if not handful:
            novel_frac, novel_genes = 1.0, gc
            r["reason"] = "anchor: best representative of the ASV"
            role = "anchor"
        else:
            if gene_provider:
                novel = gset - union
                novel_frac = len(novel) / max(len(gset), 1)
                novel_genes = len(novel)
            else:
                novel_frac = _estimate_novelty(r, handful, k)
                novel_genes = int(round(novel_frac * gc))
            cand_contam = (1.0 - h["identity"]) * novel_genes
            if novel_frac < k["novel_min"]:
                excl(r, "marginal_gain",
                     f"adds only ~{novel_frac*100:.0f}% new genes (< {k['novel_min']*100:.0f}%): "
                     f"redundant with the existing union")
                continue
            if low_conf + (1 if is_low else 0) > k["max_low_conf"]:
                excl(r, "contamination_budget",
                     f"would exceed max {k['max_low_conf']} low-confidence (uncultured/MAG) genome(s) in the union")
                continue
            if anchor_genes and (contam_genes + cand_contam) > tier_budget * anchor_genes:
                excl(r, "contamination_budget",
                     f"estimated off-target genes (~{int(contam_genes + cand_contam)}) would exceed the "
                     f"{tier} budget {tier_budget*100:.0f}% of the anchor genome (~{int(tier_budget*anchor_genes)})")
                continue
            if len(handful) >= tier_cap:
                excl(r, "tier_cap", f"{tier}-tier ceiling {tier_cap} reached (confidence-scaled)")
                continue
            role = "congener"
            r["reason"] = f"adds ~{novel_frac*100:.0f}% new genes within the contamination budget (pangenome breadth)"
        # include
        r["disposition"], r["stage"] = "selected", "included"
        r["role"] = role
        r["confidence"] = _confidence(r, cons, midas)
        r["est_genes"] = gc
        r["est_novel_genes"] = novel_genes
        r["est_novel_fraction"] = round(novel_frac, 3)
        if gene_provider:
            union |= gset
        if role == "anchor":
            anchor_genes = gc
        else:
            contam_genes += (1.0 - h["identity"]) * novel_genes
        low_conf += 1 if is_low else 0
        handful.append(r)

    # candidates that survived to stage 7 but were never visited (handful filled) -> mark
    for r in candidates:
        if r["disposition"] is None:
            excl(r, "complete", "handful already satisfied by higher-ranked representatives")

    if not handful:
        flags.append("abstain_reliability")
    return finish(handful, tier, cons)


# --------------------------------------------------------------------------- provenance
def _provenance(handful, cons, gene_data):
    if not handful:
        return dict(gene_data=gene_data, anchor_genome=None, est_union_genes=0,
                    est_core_genes=0, est_accessory_genes=0, n_low_confidence=0,
                    consensus_family=cons)
    anchor = handful[0]
    core = anchor.get("est_genes", 0)
    union = core + sum(r.get("est_novel_genes", 0) for r in handful[1:])
    nlc = sum(1 for r in handful if r["confidence"] == "low")
    return dict(
        gene_data=gene_data, anchor_genome=anchor["genome_id"],
        est_union_genes=union, est_core_genes=core, est_accessory_genes=union - core,
        n_low_confidence=nlc, consensus_family=cons,
    )


def _clean_selected(r):
    return dict(
        rank=r["rank"], genome_id=r["genome_id"], taxon_id=r["taxon_id"],
        organism=r["organism"], identity=r["identity"], align_score=r["align_score"],
        family=r["family"], genus=r["genus"], species=r["species"],
        role=r.get("role"), confidence=r.get("confidence"),
        est_novel_fraction=r.get("est_novel_fraction"), est_genes=r.get("est_genes"),
    )


def _clean_decision(r):
    out = dict(
        rank=r["rank"], genome_id=r["genome_id"], organism=r["organism"],
        identity=r["identity"], align_score=r["align_score"], coverage=r["coverage"],
        edlib_agreement=r["edlib_agreement"], family=r["family"], genus=r["genus"],
        species=r["species"], disposition=r["disposition"], stage=r["stage"],
        reason=r["reason"],
    )
    if r["disposition"] == "selected":
        out["role"] = r.get("role")
        out["confidence"] = r.get("confidence")
    return out


# --------------------------------------------------------------------------- BV-BRC genes
def bvbrc_gene_provider(cache_path="bvbrc_cache/genome_gene_families.json", family="pgfam_id", timeout=60):
    """Return a cached genome_id -> set(orthologous-family ids) function backed by the
    BV-BRC API. Use as gene_provider= for EXACT novelty / core-accessory tagging.

    PGFam (global protein families) make 'same gene' comparable across genomes. Results
    are cached to disk so the ~21.7k-genome corpus is fetched at most once.
    """
    import os
    import urllib.parse
    import urllib.request

    cache = {}
    if os.path.exists(cache_path):
        cache = {gid: set(v) for gid, v in json.load(open(cache_path)).items()}

    def _fetch(gid):
        url = "https://www.bv-brc.org/api/genome_feature/"
        q = (f"eq(genome_id,{gid})&eq(feature_type,CDS)&select({family})"
             f"&limit(50000)&http_accept=application/json")
        req = urllib.request.Request(url, data=q.encode(), headers={
            "Content-Type": "application/x-www-form-urlencoded", "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            rows = json.load(resp)
        return {r[family] for r in rows if r.get(family)}

    def provider(gid):
        gid = str(gid)
        if gid not in cache:
            try:
                cache[gid] = _fetch(gid)
            except Exception:
                cache[gid] = set()
            json.dump({g: sorted(v) for g, v in cache.items()}, open(cache_path, "w"))
        return cache[gid]

    return provider


# --------------------------------------------------------------------------- driver
def build_audit(hits, *, knobs=None, gene_provider=None, gene_count_fn=None, limit=0):
    """Run selection over every ASV; return (audit dict, summary dict)."""
    asvs = list(hits)
    if limit:
        asvs = asvs[:limit]
    audit = {}
    tiers, sizes, flagc = Counter(), Counter(), Counter()
    n_excl_stage = Counter()
    for asv in asvs:
        res = select_representatives(hits[asv], knobs=knobs,
                                     gene_provider=gene_provider, gene_count_fn=gene_count_fn)
        audit[asv] = res
        tiers[res["tier"]] += 1
        sizes[res["n_selected"]] += 1
        for fl in res["flags"]:
            flagc[fl.split(":")[0]] += 1
        for d in res["decisions"]:
            if d["disposition"] == "excluded":
                n_excl_stage[d["stage"]] += 1
    n = len(asvs)
    summary = dict(
        n_asvs=n,
        tiers={t: tiers[t] for t in ("species", "genus", "family", "below")},
        handful_sizes=dict(sorted(sizes.items())),
        mean_handful=round(sum(s * c for s, c in sizes.items()) / max(n, 1), 3),
        abstained=sizes.get(0, 0),
        flags=dict(flagc),
        exclusions_by_stage=dict(n_excl_stage.most_common()),
    )
    return audit, summary


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True, help="asv_top20_alignment_hits.json (SW hits)")
    ap.add_argument("--out", required=True, help="output audit JSON path")
    ap.add_argument("--gene-provider", choices=["auto", "estimate", "bvbrc"], default="auto",
                    help="auto = exact (bvbrc) if the gene cache exists else estimator; "
                         "bvbrc = exact (lazily fetches misses); estimate = offline taxonomy estimator")
    ap.add_argument("--gene-cache", default="bvbrc_cache/genome_gene_families.json")
    ap.add_argument("--limit", type=int, default=0, help="first N ASVs (smoke test)")
    args = ap.parse_args()

    t0 = time.perf_counter()
    hits = json.load(open(args.hits))
    use_exact = args.gene_provider == "bvbrc" or (
        args.gene_provider == "auto" and os.path.exists(args.gene_cache))
    gp = bvbrc_gene_provider(args.gene_cache) if use_exact else None
    if args.gene_provider == "auto":
        print(f"[selection] gene_data={'bvbrc (exact)' if gp else 'estimated'} "
              f"(cache {'found' if gp else 'absent'}: {args.gene_cache})")
    audit, summary = build_audit(hits, gene_provider=gp, limit=args.limit)

    out = {"_meta": {
        "generated_from": args.hits,
        "scheme": "local SW match+2/mismatch-3/gap_open-5/gap_extend-2",
        "gene_data": "bvbrc" if gp else "estimated",
        "knobs": DEFAULT_KNOBS,
        "summary": summary,
    }, **audit}
    json.dump(out, open(args.out, "w"), indent=1)
    print(f"[selection] {summary['n_asvs']} ASVs -> {args.out} ({time.perf_counter()-t0:.1f}s)")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
