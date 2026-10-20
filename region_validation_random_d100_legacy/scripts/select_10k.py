#!/usr/bin/env python
"""S1 — deterministic selection of N taxonomically diverse, full-length-16S BV-BRC genomes.

Implements DESIGN.md §2 verbatim (the "final algorithm — 10,000 diverse" block),
with the USER DECISIONS banner honored:

  * **Archaea are KEPT** (stratified by domain, never dropped) — the √-flattening
    over-weights the ~4,964-archaea pool on purpose so that V4 / V4-V5 archaeal
    performance is measurable downstream.
  * Every intermediate is written under ``region_validation/data/`` via ``PATHS``.
  * The stage is **resumable**: if both outputs already exist (and are non-empty)
    it short-circuits unless ``--force`` is given.

Determinism rule (DESIGN §2): iterate only over ``sorted(...)``; a single
``random.Random(SEED)`` is consumed in sorted order ⇒ byte-identical output across
runs and machines.  ``--n`` controls the budget (default 10000; the pilot uses
``--n 500``); ``SEED`` / ``α`` / ``CAP_FRAC`` are pinned constants.

Sole interpreter: ``/home/freiburger/Documents/py_venv/bin/python``.

Inputs (never hardcoded elsewhere — see ``_common`` / the module constants below):
  * ``/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_ID.json``
      ``{operon_md5: representative_header}``  (header carries genome_id + organism)
  * ``/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_seq.json``
      ``{operon_md5: sequence}``
  * prFBA taxdump ``/home/freiburger/Documents/prFBA/{nodes,names}.dmp`` (taxopy 0.14.0)

Outputs:
  * ``PATHS.selection_json``  — the per-genome manifest (one entry per selected
    genome: genome_id, taxon_id, best_16s_md5, best_16s_len, n_full_length_16s,
    domain..species, and the seed/alpha/cap_frac provenance block).
  * ``PATHS.src16s_fasta``  — one FASTA record per full-length 16S operon of each
    selected genome; header ``>{genome_id}|{operon_md5}`` (genome_id-bearing, so
    downstream in-silico PCR keeps provenance).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import re
import sys

# Make `import _common` work regardless of cwd (agent threads reset cwd).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common
from _common import PATHS, RANKS, done

# --------------------------------------------------------------------------- #
# Inputs (the only place these source paths are named)
# --------------------------------------------------------------------------- #
MODEL_INPUTS = "/home/freiburger/Documents/codiffusion_bioreactor/model_inputs"
MD5_ID_JSON = os.path.join(MODEL_INPUTS, "16S_md5_ID.json")
MD5_SEQ_JSON = os.path.join(MODEL_INPUTS, "16S_md5_seq.json")

# --------------------------------------------------------------------------- #
# Pinned constants (DESIGN §2)
# --------------------------------------------------------------------------- #
SEED = 1729
FULL_LEN = 1400
FALLBACK_LEN = 1200
ALPHA = 0.5          # √-flattening exponent on availability
CAP_FRAC = 0.35      # no child may take more than ceil(CAP_FRAC * budget)

# Allowed domains (USER DECISION: keep Archaea, stratify; drop everything else).
ALLOWED_DOMAINS = ("Archaea", "Bacteria")

# genome_id shape: <taxon_id>.<assembly_index>
_GID_RE = re.compile(r"^\d+\.\d+$")


# --------------------------------------------------------------------------- #
# Step 1 — candidate table (one pass over 16S_md5_ID / 16S_md5_seq)
# --------------------------------------------------------------------------- #
def parse_genome_id(header: str):
    """Extract genome_id from a representative header (DESIGN §2 step 1).

    ``header[rfind('[')+1 : rfind(']')]`` isolates the bracketed organism block;
    ``rsplit('|', 1)[-1]`` takes the trailing ``... | <genome_id>`` field;
    validate against ``^\\d+\\.\\d+$``.  Returns ``None`` if malformed.
    """
    lb = header.rfind("[")
    rb = header.rfind("]")
    if lb < 0 or rb < 0 or rb <= lb:
        return None
    gid = header[lb + 1 : rb].rsplit("|", 1)[-1].strip()
    return gid if _GID_RE.match(gid) else None


def build_candidate_table(md5_id, md5_seq):
    """Collapse per-operon md5 records into a per-genome candidate table.

    Returns ``{genome_id: {best_len, best_md5, n_16s, n_full, full_md5s}}`` where
    ``full_md5s`` is the set of this genome's full-length (>= FULL_LEN) operon
    md5s — needed both for the n_full count and for the src16s FASTA emission.

    Tie rule for ``best_md5`` (DESIGN §2 step 1): longest seq; tie -> the
    lexicographically-smallest md5 (deterministic).
    """
    table: dict[str, dict] = {}
    missing_seq = 0
    bad_header = 0

    for md5, header in md5_id.items():
        gid = parse_genome_id(header)
        if gid is None:
            bad_header += 1
            continue
        seq = md5_seq.get(md5)
        if seq is None:
            missing_seq += 1
            continue
        ln = len(seq)

        rec = table.get(gid)
        if rec is None:
            rec = {
                "best_len": -1,
                "best_md5": None,
                "n_16s": 0,
                "n_full": 0,
                "full_md5s": set(),
            }
            table[gid] = rec

        rec["n_16s"] += 1
        if ln >= FULL_LEN:
            rec["n_full"] += 1
            rec["full_md5s"].add(md5)

        # best = longest; tie -> lexicographically-smallest md5.
        if ln > rec["best_len"] or (ln == rec["best_len"] and md5 < rec["best_md5"]):
            rec["best_len"] = ln
            rec["best_md5"] = md5

    return table, {"bad_header": bad_header, "missing_seq": missing_seq}


# --------------------------------------------------------------------------- #
# Step 2 — lineage with rank fallback (taxopy)
# --------------------------------------------------------------------------- #
def lineage_with_fallback(taxon_id, taxdb, cache):
    """Resolve a taxon to the 7 lower-case ranks with stable placeholders.

    DESIGN §2 step 2:
      * ``domain`` read with ``superkingdom`` fallback; require it in
        {Bacteria, Archaea} and a non-null ``phylum`` else return ``None`` (drop).
      * each missing *intermediate* rank (class..genus) gets a stable placeholder
        ``unclassified_<parentrank>:<parentname>`` so genus-less / family-less
        environmental clades form their OWN bin instead of collapsing together.
      * ``species`` falls back to ``sp_taxon_<tid>``.

    Returns an ordered dict over ``RANKS`` (domain..species) or ``None`` to drop.
    Cached per taxon_id.
    """
    if taxon_id in cache:
        return cache[taxon_id]

    import taxopy

    try:
        rd = taxopy.Taxon(int(taxon_id), taxdb).rank_name_dictionary
    except Exception:
        cache[taxon_id] = None
        return None

    domain = rd.get("domain") or rd.get("superkingdom")
    phylum = rd.get("phylum")
    if domain not in ALLOWED_DOMAINS or not phylum:
        cache[taxon_id] = None
        return None

    lin = {"domain": domain, "phylum": phylum}
    # Intermediate ranks: class, order, family, genus.  Fill placeholders off the
    # deepest *named* ancestor so each gap forms its own stable bin.
    parent_rank, parent_name = "phylum", phylum
    for rank in ("class", "order", "family", "genus"):
        name = rd.get(rank)
        if name:
            lin[rank] = name
            parent_rank, parent_name = rank, name
        else:
            lin[rank] = f"unclassified_{parent_rank}:{parent_name}"
    # species fallback keyed on the concrete taxon id (stable, unique).
    lin["species"] = rd.get("species") or f"sp_taxon_{taxon_id}"

    cache[taxon_id] = lin
    return lin


# --------------------------------------------------------------------------- #
# Step 3 — species dedup -> one representative genome per (domain,phylum,species)
# --------------------------------------------------------------------------- #
def rep_sort_key(entry):
    """Deterministic ladder to pick the species representative (DESIGN §2 step 3).

    ``(-n_full, -best_len, -int(assembly_index), genome_id)`` — most full-length
    operons first, then longest best operon, then newest assembly, then
    genome_id as the final stable tiebreak.
    """
    gid = entry["genome_id"]
    assembly_index = int(gid.split(".", 1)[1])
    return (-entry["n_full"], -entry["best_len"], -assembly_index, gid)


# --------------------------------------------------------------------------- #
# Step 4 — stratified, capped, deterministic allocation
# --------------------------------------------------------------------------- #
def _largest_remainder(children, budget, avail, rng):
    """Hamilton largest-remainder allocation of ``budget`` across ``children``.

    ``weight(child) = avail(child)**ALPHA`` (√-flattening); ``share = budget*w/Σw``;
    floor each share, distribute the leftover to the largest fractional remainders
    (ties broken on the sorted child key, then on the integer floor for stability).
    Each child is hard-capped at ``min(avail, ceil(CAP_FRAC*budget))``; a cap is
    relaxed only when *every* uncapped child is exhausted and budget remains.

    ``children`` MUST already be sorted (callers pass ``sorted(...)``).  Returns
    ``{child: alloc}`` summing to ``min(budget, Σ avail)``.
    """
    if budget <= 0 or not children:
        return {c: 0 for c in children}

    total_avail = sum(avail[c] for c in children)
    if total_avail <= budget:
        # Everything fits — take all of every child.
        return {c: avail[c] for c in children}

    cap = math.ceil(CAP_FRAC * budget)
    caps = {c: min(avail[c], cap) for c in children}

    weights = {c: (avail[c] ** ALPHA) for c in children}
    wsum = sum(weights.values())

    # Initial floored shares, clamped to caps.
    alloc = {}
    remainders = {}
    for c in children:
        raw = budget * weights[c] / wsum if wsum > 0 else 0.0
        fl = int(math.floor(raw))
        fl = min(fl, caps[c])
        alloc[c] = fl
        remainders[c] = raw - math.floor(raw)

    assigned = sum(alloc.values())
    leftover = budget - assigned

    # Distribute leftover by largest fractional remainder, respecting caps;
    # relax caps (round-robin in sorted order) only once all caps are saturated.
    if leftover > 0:
        # Sort children by (-remainder, sorted-key) for a deterministic order.
        ordered = sorted(children, key=lambda c: (-remainders[c], c))
        # Pass 1: capped largest-remainder.
        for c in ordered:
            if leftover <= 0:
                break
            room = caps[c] - alloc[c]
            if room > 0:
                give = min(room, leftover)
                alloc[c] += give
                leftover -= give
        # Pass 2: cap relaxation up to true availability (clades genuinely
        # exhausted otherwise).  Round-robin in sorted key order for determinism.
        while leftover > 0:
            progressed = False
            for c in sorted(children):
                if leftover <= 0:
                    break
                if alloc[c] < avail[c]:
                    alloc[c] += 1
                    leftover -= 1
                    progressed = True
            if not progressed:
                break  # truly exhausted (shouldn't happen: total_avail>budget)

    return alloc


def allocate(tree, budget, rng):
    """Recursively allocate ``budget`` genomes across the domain->...->species tree.

    ``tree`` is the nested dict built by :func:`build_tree`:
      * internal node -> ``{"children": {key: subtree}, "avail": int}``
      * leaf (species) -> ``{"species_genomes": [candidate, ...], "avail": int}``

    Returns the flat list of chosen candidate dicts.  Iterates only over
    ``sorted(children)`` and consumes ``rng`` in sorted order ⇒ deterministic.
    """
    chosen = []

    children = tree.get("children")
    if children is not None:
        keys = sorted(children)
        avail = {k: children[k]["avail"] for k in keys}
        alloc = _largest_remainder(keys, budget, avail, rng)
        for k in keys:
            sub_budget = alloc[k]
            if sub_budget > 0:
                chosen.extend(allocate(children[k], sub_budget, rng))
        return chosen

    # Leaf: a species bin holding one representative genome.
    genomes = tree["species_genomes"]
    if budget >= len(genomes):
        chosen.extend(genomes)
    else:
        shuffled = list(sorted(genomes, key=lambda g: g["genome_id"]))
        rng.shuffle(shuffled)
        prefix = shuffled[:budget]
        # re-sort the chosen prefix for a stable manifest order.
        chosen.extend(sorted(prefix, key=lambda g: g["genome_id"]))
    return chosen


def build_tree(reps):
    """Build the nested domain->phylum->class->order->family->genus->species tree.

    ``reps`` is the list of per-species representative candidate dicts (each with a
    ``lineage`` over RANKS).  Each leaf is one species holding exactly one
    representative genome (species dedup already collapsed it).  ``avail`` is the
    number of distinct species reachable below a node (the allocation weight base).
    """
    root = {"children": {}, "avail": 0}
    for rep in reps:
        lin = rep["lineage"]
        node = root
        for rank in RANKS[:-1]:  # domain..genus internal levels
            key = lin[rank]
            child = node["children"].get(key)
            if child is None:
                child = {"children": {}, "avail": 0}
                node["children"][key] = child
            node = child
        # species leaf
        sp_key = lin["species"]
        leaf = node["children"].get(sp_key)
        if leaf is None:
            leaf = {"species_genomes": [], "avail": 0}
            node["children"][sp_key] = leaf
        leaf["species_genomes"].append(rep)

    # Compute avail bottom-up (= number of distinct species below each node).
    def _fill_avail(node):
        if "species_genomes" in node:
            node["avail"] = 1  # one species
            return 1
        total = 0
        for child in node["children"].values():
            total += _fill_avail(child)
        node["avail"] = total
        return total

    _fill_avail(root)
    return root


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def select(n, force):
    """Run the full §2 selection for budget ``n``; write manifest + src16s FASTA."""
    out_json = PATHS.selection_json
    out_fasta = PATHS.src16s_fasta

    if not force and done(out_json) and done(out_fasta):
        return (
            f"[skip] selection already present for budget; "
            f"{out_json} and {out_fasta} both exist (use --force to rebuild)."
        )

    print(f"[load] {MD5_ID_JSON}", file=sys.stderr)
    with open(MD5_ID_JSON) as fh:
        md5_id = json.load(fh)
    print(f"[load] {MD5_SEQ_JSON}", file=sys.stderr)
    with open(MD5_SEQ_JSON) as fh:
        md5_seq = json.load(fh)

    print("[step1] building candidate table", file=sys.stderr)
    table, diag = build_candidate_table(md5_id, md5_seq)
    print(
        f"[step1] {len(table)} genomes "
        f"(bad_header={diag['bad_header']}, missing_seq={diag['missing_seq']})",
        file=sys.stderr,
    )

    print("[step2/3] lineage + species dedup", file=sys.stderr)
    import taxopy

    taxdb = taxopy.TaxDb(
        nodes_dmp=_common.TAXDUMP_NODES, names_dmp=_common.TAXDUMP_NAMES
    )
    lin_cache: dict = {}

    # Per-species best representative, with full-length enforced BEFORE allocation.
    species_best: dict[tuple, dict] = {}
    n_dropped_short = n_dropped_lineage = 0
    fallback_pool: dict[tuple, dict] = {}

    for gid in sorted(table):
        rec = table[gid]
        taxon_id = int(gid.split(".", 1)[0])
        lin = lineage_with_fallback(taxon_id, taxdb, lin_cache)
        if lin is None:
            n_dropped_lineage += 1
            continue

        cand = {
            "genome_id": gid,
            "taxon_id": taxon_id,
            "best_md5": rec["best_md5"],
            "best_len": rec["best_len"],
            "n_full": rec["n_full"],
            "full_md5s": rec["full_md5s"],
            "lineage": lin,
        }
        sp_key = (lin["domain"], lin["phylum"], lin["species"])

        if rec["best_len"] >= FULL_LEN:
            prev = species_best.get(sp_key)
            if prev is None or rep_sort_key(cand) < rep_sort_key(prev):
                species_best[sp_key] = cand
        elif rec["best_len"] >= FALLBACK_LEN:
            # candidate fallback rep, used only to backfill species with no
            # >=FULL_LEN member.
            prev = fallback_pool.get(sp_key)
            if prev is None or rep_sort_key(cand) < rep_sort_key(prev):
                fallback_pool[sp_key] = cand
        else:
            n_dropped_short += 1

    # Optional second pass: admit a 1200<=len<1400 representative ONLY for species
    # with no full-length member (16S-poor-clade backfill, DESIGN §2 step 5).
    n_backfill = 0
    for sp_key, cand in fallback_pool.items():
        if sp_key not in species_best:
            species_best[sp_key] = cand
            n_backfill += 1

    reps = list(species_best.values())
    print(
        f"[step3] {len(reps)} species reps "
        f"(full-length {len(reps) - n_backfill}, fallback-backfill {n_backfill}; "
        f"dropped short={n_dropped_short}, lineage={n_dropped_lineage})",
        file=sys.stderr,
    )

    n_domains = len({r["lineage"]["domain"] for r in reps})
    n_phyla = len({(r["lineage"]["domain"], r["lineage"]["phylum"]) for r in reps})
    print(f"[step3] {n_domains} domains / {n_phyla} phyla in rep pool", file=sys.stderr)

    print(f"[step4] stratified allocation to budget={n}", file=sys.stderr)
    rng = random.Random(SEED)
    tree = build_tree(reps)
    target = min(n, len(reps))
    if target < n:
        print(
            f"[warn] only {len(reps)} species reps available (< requested {n}); "
            f"selecting all of them.",
            file=sys.stderr,
        )
    chosen = allocate(tree, target, rng)

    # Largest-remainder makes the root sum exact (DESIGN §2 step 5 guarantee).
    assert len(chosen) == target, f"alloc sum {len(chosen)} != target {target}"

    # Stable manifest order.
    chosen.sort(key=lambda c: c["genome_id"])

    # ---- write manifest (PATHS.selection_json) -------------------------------- #
    provenance = {
        "seed": SEED,
        "alpha": ALPHA,
        "cap_frac": CAP_FRAC,
        "full_len": FULL_LEN,
        "fallback_len": FALLBACK_LEN,
    }
    manifest = {
        "meta": {
            "target_n": n,
            "selected_n": len(chosen),
            "n_species_reps": len(reps),
            "n_backfill_reps": n_backfill,
            "n_domains": n_domains,
            "n_phyla": n_phyla,
            **provenance,
            "inputs": {"md5_id": MD5_ID_JSON, "md5_seq": MD5_SEQ_JSON},
        },
        "genomes": [],
    }
    for c in chosen:
        lin = c["lineage"]
        manifest["genomes"].append(
            {
                "genome_id": c["genome_id"],
                "taxon_id": c["taxon_id"],
                "best_16s_md5": c["best_md5"],
                "best_16s_len": c["best_len"],
                "n_full_length_16s": c["n_full"],
                "domain": lin["domain"],
                "phylum": lin["phylum"],
                "class": lin["class"],
                "order": lin["order"],
                "family": lin["family"],
                "genus": lin["genus"],
                "species": lin["species"],
                **provenance,
            }
        )

    tmp_json = out_json + ".tmp"
    with open(tmp_json, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=False)
    os.replace(tmp_json, out_json)

    # ---- write src16s FASTA (PATHS.src16s_fasta) ------------------------------ #
    # ONE record per selected genome: its representative (best) full-length 16S
    # operon, header '>{genome_id}'. One organism = one query per region downstream;
    # intra-genomic operons are ~identical across the V-regions, and a single id per
    # genome keeps the amplicon ids '{genome_id}__{region}' unique (multi-operon
    # emission collides). Emitted in genome_id-sorted order ⇒ deterministic.
    n_records = 0
    tmp_fasta = out_fasta + ".tmp"
    with open(tmp_fasta, "w") as fh:
        for c in chosen:
            seq = md5_seq[c["best_md5"]]
            fh.write(f">{c['genome_id']}\n{seq}\n")
            n_records += 1
    os.replace(tmp_fasta, out_fasta)

    return (
        f"Selected {len(chosen)}/{n} genomes "
        f"({n_domains} domains, {n_phyla} phyla, seed={SEED}); "
        f"wrote {n_records} representative 16S genes to FASTA."
    )


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="S1 — deterministic 10k diverse-genome selection (DESIGN §2)."
    )
    ap.add_argument(
        "--n",
        type=int,
        default=10000,
        help="selection budget (default 10000; the pilot uses --n 500).",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="rebuild even if the output files already exist.",
    )
    args = ap.parse_args(argv)

    summary = select(args.n, args.force)
    print(summary)


if __name__ == "__main__":
    main()
