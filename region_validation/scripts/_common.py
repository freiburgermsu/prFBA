#!/usr/bin/env python
"""Shared contract module for the 16S hypervariable-region validation pipeline.

Every validation script (S1-S12, render_region_*) imports this module for the
single source of truth on:

  * directory layout + an absolute path for *every* intermediate (``PATHS``),
  * the canonical rank ordering (``RANKS`` / ``SCORED``),
  * the corrected 8-region primer panel with runtime-compiled regexes
    (``load_panel`` / ``revcomp`` / ``iupac_to_regex``),
  * taxon-name normalization (``_norm``, mirroring ``select_references._norm``),
  * lineage resolution through the prFBA taxdump (``lineage_for``, mirroring
    ``edlib_biopython_hits.lineage_for``),
  * Wilson confidence intervals (``wilson_ci``),
  * resumability helpers (``done`` / ``ensure_dirs``).

Sole interpreter: ``/home/freiburger/Documents/py_venv/bin/python``.

Nothing here triggers a slow import (taxopy ``TaxDb`` and the region-panel JSON
are loaded lazily on first use), so ``import _common`` is cheap and side-effect
free apart from idempotent ``mkdir`` of the output directories.
"""
from __future__ import annotations

import json
import os
import re

# --------------------------------------------------------------------------- #
# Directory constants
# --------------------------------------------------------------------------- #
BASE = "/home/freiburger/Documents/prFBA/region_validation"
DATA = os.path.join(BASE, "data")
SCRIPTS = os.path.join(BASE, "scripts")
FIGS = os.path.join(BASE, "figures")

# The corrected, collapsed 8-region panel lives here (V6-V8/V7-V9 mislabel
# fixed, 1100F/1492R collapsed to a single V7-V9 row).
REGION_PANEL_JSON = os.path.join(BASE, "region_panel_final.json")

# prFBA taxdump (taxopy).  Same files used by edlib_biopython_hits.lineage_for.
TAXDUMP_NODES = "/home/freiburger/Documents/prFBA/nodes.dmp"
TAXDUMP_NAMES = "/home/freiburger/Documents/prFBA/names.dmp"


class PATHS:
    """Absolute path for every pipeline intermediate.

    All intermediates live under ``DATA`` except the shared PGFam cache, which
    is deliberately the prFBA-wide cache (one namespace for truth and predicted
    gene sets => apples-to-apples recall/precision).  Attributes are plain class
    attributes (no instantiation needed): ``PATHS.selection_json`` etc.
    """

    # --- S1: deterministic 10k selection ---------------------------------- #
    selection_json = os.path.join(DATA, "selection_10k.json")
    selection_manifest = os.path.join(DATA, "selection_10k.manifest.json")
    src16s_fasta = os.path.join(DATA, "src16s_10k.fasta")

    # --- S3/S4: in-silico PCR -> combined amplicons + dedup --------------- #
    combined_fasta = os.path.join(DATA, "combined_amplicons.fasta")
    combined_stats = os.path.join(DATA, "combined_amplicons.stats.json")
    unique_fasta = os.path.join(DATA, "combined_amplicons.unique.fasta")
    amplicon_fanout = os.path.join(DATA, "amplicon_md5_fanout.json")

    # --- S5: GPU alignment pass ------------------------------------------- #
    hits_dir = os.path.join(DATA, "hits")
    hits_json = os.path.join(hits_dir, "asv_top20_alignment_hits.json")

    # --- S6: selection (include-self headline + exclude-self cascade) ----- #
    selection_out = os.path.join(DATA, "selection.json")

    # --- S7: benchmark truth ---------------------------------------------- #
    truth_json = os.path.join(DATA, "benchmark_truth.json")

    # --- S8: PGFam gene-family cache (prFBA-wide; NOT under DATA) ---------- #
    pgfam_cache = "/home/freiburger/Documents/prFBA/bvbrc_cache/genome_gene_families.json"

    # --- S9: taxonomic-accuracy tables ------------------------------------ #
    taxacc_jsonl = os.path.join(DATA, "taxacc_per_record.jsonl")
    concordance_long_csv = os.path.join(DATA, "concordance_long.csv")
    asv_region_concordance_csv = os.path.join(DATA, "asv_region_concordance.csv")
    self_recovery_csv = os.path.join(DATA, "self_recovery.csv")
    confusion_genus_csv = os.path.join(DATA, "confusion_genus.csv")
    ambiguity_csv = os.path.join(DATA, "ambiguity.csv")
    taxacc_summary_json = os.path.join(DATA, "taxacc_region_summary.json")

    # --- S10: gene-capture tables ----------------------------------------- #
    genecap_cells_csv = os.path.join(DATA, "gene_capture_cells.csv")
    gene_capture_long_csv = os.path.join(DATA, "gene_capture_long.csv")
    genecap_summary_json = os.path.join(DATA, "gene_capture_by_region.json")

    # --- S11: region metadata + amplification ----------------------------- #
    region_meta_csv = os.path.join(DATA, "region_meta.csv")
    amplification_by_phylum_csv = os.path.join(DATA, "amplification_by_phylum.csv")


# --------------------------------------------------------------------------- #
# Rank ordering
# --------------------------------------------------------------------------- #
# Selection / lineage-key convention is lower-case 7-rank (matches DESIGN §2).
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# Scored ranks for accuracy: Phylum..Species.  Kingdom/domain is a sanity gate,
# not scored (DESIGN §6).  Capitalized to match lineage_for() dict keys.
SCORED = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


# --------------------------------------------------------------------------- #
# Figure constants — consistent region ordering + palette + depth ramp
# --------------------------------------------------------------------------- #
# DESIGN §8: every figure shares one REGION_ORDER (short -> long amplicon), one
# REGION_PALETTE (a fixed color per region used *everywhere*), and a DEPTH_RAMP
# (a sequential ramp over the deepest-correct-rank stacked bar).  These are
# derived deterministically from region_panel_final.json (amp_bp) so the order
# is "by amplicon length" without hard-coding a region list.  Loaded lazily.
_FIG_CONST_CACHE: dict = {}


def _build_fig_consts() -> dict:
    """Build (and cache) REGION_ORDER / REGION_PALETTE / DEPTH_RAMP.

    * ``REGION_ORDER`` — region names sorted by published amplicon length
      (``amp_bp``) ascending, tie-broken on name, so the y-axis of every
      heatmap runs short -> long exactly once across the whole figure suite.
    * ``REGION_PALETTE`` — one stable color per region, sampled evenly along
      ``viridis`` in REGION_ORDER, so a region keeps its color in every figure.
    * ``DEPTH_RAMP`` — ordered ``{label: color}`` for the deepest-correct-rank
      stacked bar: a light->dark sequential ramp over Phylum..Species plus a
      neutral grey for the "none / below-phylum" floor.
    """
    if _FIG_CONST_CACHE:
        return _FIG_CONST_CACHE

    import matplotlib.cm as _cm
    import numpy as _np

    with open(REGION_PANEL_JSON) as fh:
        _panel = json.load(fh)
    _regions = _panel["regions"]
    # by amplicon length (short -> long), name as deterministic tiebreak.
    order = [r["name"] for r in sorted(_regions, key=lambda r: (r.get("amp_bp", 0), r["name"]))]
    amp_bp = {r["name"]: r.get("amp_bp") for r in _regions}

    _vir = _cm.get_cmap("viridis")
    n = max(len(order), 1)
    palette = {
        name: _vir(0.0 if n == 1 else i / (n - 1))
        for i, name in enumerate(order)
    }

    # Deepest-correct-rank ramp: "none" (no resolution) -> Species (deepest).
    depth_labels = ["none"] + list(SCORED)
    _ylgnbu = _cm.get_cmap("YlGnBu")
    ramp = {"none": (0.85, 0.85, 0.85, 1.0)}  # neutral grey floor
    for i, lab in enumerate(SCORED):
        ramp[lab] = _ylgnbu(0.25 + 0.7 * (i / max(len(SCORED) - 1, 1)))

    _FIG_CONST_CACHE.update(
        REGION_ORDER=order,
        REGION_AMP_BP=amp_bp,
        REGION_PALETTE=palette,
        DEPTH_LABELS=depth_labels,
        DEPTH_RAMP=ramp,
    )
    return _FIG_CONST_CACHE


def __getattr__(name):
    """Lazily expose REGION_ORDER / REGION_PALETTE / DEPTH_RAMP / DEPTH_LABELS /
    REGION_AMP_BP as module attributes (built from the panel on first access)."""
    if name in ("REGION_ORDER", "REGION_PALETTE", "DEPTH_RAMP", "DEPTH_LABELS", "REGION_AMP_BP"):
        return _build_fig_consts()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# --------------------------------------------------------------------------- #
# Directory creation (idempotent)
# --------------------------------------------------------------------------- #
def ensure_dirs() -> None:
    """Create DATA, FIGS and the hits/ subdir; idempotent."""
    for d in (DATA, FIGS, PATHS.hits_dir):
        os.makedirs(d, exist_ok=True)


# Make DATA / FIGS / hits_dir creation idempotent at import time so any importer
# can immediately write its output without a manual mkdir.
ensure_dirs()


# --------------------------------------------------------------------------- #
# Resumability
# --------------------------------------------------------------------------- #
def done(path: str) -> bool:
    """True iff ``path`` exists and is non-empty (a usable checkpoint)."""
    try:
        return os.path.getsize(path) > 0
    except OSError:
        return False


# --------------------------------------------------------------------------- #
# Taxon-name normalization (mirrors select_references._norm)
# --------------------------------------------------------------------------- #
# Strip leading "Candidatus ", "Ca." / "Ca_" and "midas_<letter>_" prefixes,
# then casefold.  Kept byte-identical to select_references._PREFIX so that the
# truth and predicted sides compare under the exact same normalization.
_PREFIX = re.compile(r"^(Candidatus\s+|Ca[._]\s*|midas_[a-z]_)\s*", re.I)


def _norm(s) -> str:
    """Normalize a taxon name for ``==`` comparison (mirror of select_references._norm)."""
    if not s:
        return ""
    return _PREFIX.sub("", str(s).strip()).casefold()


# --------------------------------------------------------------------------- #
# IUPAC primer handling + region panel
# --------------------------------------------------------------------------- #
# Canonical IUPAC degeneracy map.  Mirrors the "iupac" block in
# region_panel_final.json but is defined inline so primer helpers work even if
# the JSON is loaded lazily / out of order.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

# Watson-Crick complement including IUPAC degeneracy codes (for revcomp()).
_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}


def revcomp(seq: str) -> str:
    """Reverse complement of an IUPAC nucleotide string (uppercased)."""
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def iupac_to_regex(primer: str) -> str:
    """Convert an IUPAC primer string to a regex character-class pattern.

    Each base becomes a bracketed character class of its allowed bases, so the
    pattern matches any concrete sequence the degenerate primer can prime.
    """
    out = []
    for b in primer.upper():
        bases = IUPAC.get(b, b)
        out.append(f"[{bases}]" if len(bases) > 1 else bases)
    return "".join(out)


# Lazy singleton for the parsed/compiled panel.
_PANEL_CACHE = None


def load_panel():
    """Return the 8 regions from ``region_panel_final.json``, regex-compiled.

    Each returned dict carries the original panel fields plus:

      * ``fwd_re``  -- compiled regex for the forward primer on the sense strand.
      * ``rev_revcomp`` -- the reverse-complement of the reverse primer (the
        pattern that appears on the sense strand downstream of the insert).
      * ``rev_re``  -- compiled regex for ``rev_revcomp`` on the sense strand.

    The ``revc`` / runtime patterns are computed here via ``revcomp`` +
    ``iupac_to_regex`` (the JSON's published ``rev`` is 5'->3'; the extractor
    matches its reverse-complement on the sense strand).  Result is cached.
    """
    global _PANEL_CACHE
    if _PANEL_CACHE is not None:
        return _PANEL_CACHE

    with open(REGION_PANEL_JSON) as fh:
        panel = json.load(fh)

    regions = []
    for r in panel["regions"]:
        reg = dict(r)
        rev_rc = revcomp(r["rev"])
        reg["rev_revcomp"] = rev_rc
        reg["fwd_re"] = re.compile(iupac_to_regex(r["fwd"]))
        reg["rev_re"] = re.compile(iupac_to_regex(rev_rc))
        regions.append(reg)

    _PANEL_CACHE = regions
    return _PANEL_CACHE


# --------------------------------------------------------------------------- #
# Lineage resolution (mirrors edlib_biopython_hits.lineage_for)
# --------------------------------------------------------------------------- #
# Lazy TaxDb singleton + per-taxon lineage cache (taxopy is slow to construct).
_TAXDB = None
_LIN_CACHE: dict = {}
# Keys returned by lineage_for() — capitalized, Kingdom-first (7 ranks).
_LINEAGE_KEYS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def _get_taxdb():
    """Build the taxopy TaxDb once (lazy singleton) from the prFBA taxdump."""
    global _TAXDB
    if _TAXDB is None:
        import taxopy

        _TAXDB = taxopy.TaxDb(nodes_dmp=TAXDUMP_NODES, names_dmp=TAXDUMP_NAMES)
    return _TAXDB


def lineage_for(taxon_id) -> dict:
    """Resolve a taxon id to ``{Kingdom,Phylum,Class,Order,Family,Genus,Species}``.

    Mirrors ``edlib_biopython_hits.lineage_for``:
      * ``Kingdom = superkingdom or kingdom or domain`` (modern NCBI yields
        e.g. 'Pseudomonadati' for bacteria — fine, both truth and predicted
        sides go through this same function).
      * ``None`` taxon id -> all-None lineage.
      * any taxopy failure -> all-None lineage (logged-as-None, not raised).

    Module-level cache keyed by ``taxon_id``.
    """
    if taxon_id is None:
        return {k: None for k in _LINEAGE_KEYS}
    if taxon_id in _LIN_CACHE:
        return _LIN_CACHE[taxon_id]

    import taxopy

    try:
        rd = taxopy.Taxon(int(taxon_id), _get_taxdb()).rank_name_dictionary
        lin = {
            "Kingdom": rd.get("superkingdom") or rd.get("kingdom") or rd.get("domain"),
            "Phylum": rd.get("phylum"),
            "Class": rd.get("class"),
            "Order": rd.get("order"),
            "Family": rd.get("family"),
            "Genus": rd.get("genus"),
            "Species": rd.get("species"),
        }
    except Exception:
        lin = {k: None for k in _LINEAGE_KEYS}

    _LIN_CACHE[taxon_id] = lin
    return lin


# --------------------------------------------------------------------------- #
# Uncertainty: Wilson confidence interval
# --------------------------------------------------------------------------- #
def wilson_ci(k: int, n: int, conf: float = 0.95):
    """Wilson 95% CI for a binomial proportion ``k/n`` via statsmodels.

    Returns ``(lo, hi)``.  For ``n == 0`` returns ``(0.0, 1.0)`` (no information)
    so callers never divide by zero.
    """
    if n <= 0:
        return (0.0, 1.0)
    from statsmodels.stats.proportion import proportion_confint

    lo, hi = proportion_confint(k, n, alpha=1.0 - conf, method="wilson")
    return (float(lo), float(hi))


__all__ = [
    "BASE", "DATA", "SCRIPTS", "FIGS", "REGION_PANEL_JSON",
    "TAXDUMP_NODES", "TAXDUMP_NAMES", "PATHS",
    "RANKS", "SCORED",
    "ensure_dirs", "done",
    "_norm", "IUPAC", "revcomp", "iupac_to_regex", "load_panel",
    "lineage_for", "wilson_ci",
    "REGION_ORDER", "REGION_AMP_BP", "REGION_PALETTE", "DEPTH_LABELS", "DEPTH_RAMP",
]
