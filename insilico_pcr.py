#!/usr/bin/env python
"""
insilico_pcr.py — extract 16S sub-region amplicons from reference sequences by in-silico PCR.

The experimental ASVs in this study are V4–V5 amplicons (E. coli ~533–906, ~373 bp insert)
generated with the universal 515F / 926R primer pair (covers bacteria *and* archaea). To compare
amplicons against the BV-BRC reference space *like-with-like* (region-matched), we excise the same
insert from every reference 16S gene: locate the forward primer and the reverse-complement of the
reverse primer (IUPAC-degenerate, fuzzy ≤emax-mismatch fallback for divergent taxa) and keep the
sequence between them.

This module supports two modes:

  (1) Single-region legacy mode (`--fasta/--out`): one region, {header: insert} JSON out
      (back-compat; the original V4–V5 behaviour, unchanged defaults).

  (2) Multi-region panel mode (`--fasta/--combined-fasta [--panel ...]`): scan every genome
      against the whole REGION_PANEL and emit ONE combined FASTA where each record id encodes
      `{genome_id}__{region}` and the sequence is the amplicon (primers included by default).
      This is the front-loaded combined output for the hypervariable-region validation study.

Default V4–V5 primers (Parada/Quince, EMP):
  515F  GTGYCAGCMGCCGCGGTAA
  926R  CCGYCAATTYMTTTRAGTTT      (its rev-comp is matched on the sense strand)

Input  : {header: full_16S_sequence} JSON (e.g. model_inputs/BV_BRC_16S.json), or a FASTA file.
Output : combined FASTA (panel mode) or {header: region_insert} JSON (legacy mode), + a stats JSON.
"""
from __future__ import annotations
import argparse, json, os, re, sys, time
from pathlib import Path
import regex

# The validation pipeline's single source of truth for paths/constants lives in
# region_validation/scripts/_common.py.  Import it for PATHS (so the combined
# FASTA and its stats sibling are never hard-coded here) and the authoritative
# region_panel_final.json location.  It is optional: when insilico_pcr.py is run
# standalone (outside the validation tree) the module degrades to the
# byte-identical hard-coded fallbacks below, so nothing here requires _common.
_COMMON = None
try:
    sys.path.insert(0, "/home/freiburger/Documents/prFBA/region_validation/scripts")
    import _common as _COMMON  # type: ignore
except Exception:
    _COMMON = None

# ────────────────────────────────────────────────────────────────────────────
# IUPAC handling
# ────────────────────────────────────────────────────────────────────────────
# Forward (degenerate) IUPAC code -> regex character class.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
    "N": "[ACGT]",
}
# Complement of each IUPAC code (for reverse-complementing reverse primers).
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def iupac_to_regex(primer: str) -> str:
    """Translate a degenerate IUPAC primer (5'->3') into a regex pattern.

    e.g. GTGYCAGCMGCCGCGGTAA -> GTG[CT]CAGC[AC]GCCGCGGTAA
    """
    primer = primer.strip().upper()
    try:
        return "".join(IUPAC[b] for b in primer)
    except KeyError as e:  # surface the offending base, don't silently mistranslate
        raise ValueError(f"non-IUPAC base {e} in primer {primer!r}") from None


def revcomp(seq: str) -> str:
    """Reverse-complement an IUPAC sequence (degeneracy-preserving)."""
    return seq.strip().upper().translate(COMP)[::-1]


# ────────────────────────────────────────────────────────────────────────────
# Region panel — the 8 corrected regions (authoritative: region_panel_final.json)
# ────────────────────────────────────────────────────────────────────────────
# Each region: forward primer matched on the sense strand, reverse primer matched as its
# reverse-complement on the sense strand. `lo`/`hi` bound the INSERT (between primers, primers
# excluded) for amplification sanity; the emitted amplicon INCLUDES the primers by default.
#
# Insert bounds derive from E. coli coords (rev_primer_start − fwd_primer_end) with generous
# slack for indel/length variation across taxa. They are deliberately wide enough to admit real
# biological length variation but tight enough to reject spurious primer matches.
#
# CORRECTIONS (verifier PRIMERS, DESIGN §3): the old "V6-V8" entry was a mislabel — 1100F's
# 3' end (~1114) sits downstream of V6 (ends 1043), so 1100F/1492R amplifies V7–V9. The two
# byte-identical 1100F/1492R rows are collapsed to ONE core "V7-V9" row. The corrected primers
# (926R sense-revcomp uses K; 1492R single-R) are the published 5'->3' strings below — revcomp()
# computes the sense-strand pattern at runtime, so no manual revcomp typo can creep in.
#
# This panel is kept byte-identical to region_panel_final.json (the single source of truth);
# if that JSON is present it is loaded at import time so the two never drift.
REGION_PANEL_JSON = Path(
    _COMMON.REGION_PANEL_JSON if _COMMON is not None
    else "/home/freiburger/Documents/prFBA/region_validation/region_panel_final.json"
)

# Hard-coded fallback (mirrors region_panel_final.json exactly: 8 corrected regions, the
# 1100F/1492R duplicate collapsed to a single core "V7-V9").  `lo`/`hi` == insert_lo/insert_hi.
_REGION_PANEL_FALLBACK = [
    {
        "name": "FullLength16S",        # 27F / 1492R — positive control
        "fwd_name": "27F", "fwd": "AGAGTTTGATCMTGGCTCAG",
        "rev_name": "1492R", "rev": "TACGGYTACCTTGTTACGACTT",
        "lo": 1380, "hi": 1520, "priority": "core", "domain": "Bacteria",
        "ecoli": "8-1507",
    },
    {
        "name": "V1-V2",                # 27F / 338R
        "fwd_name": "27F", "fwd": "AGAGTTTGATCMTGGCTCAG",
        "rev_name": "338R", "rev": "TGCTGCCTCCCGTAGGAGT",
        "lo": 270, "hi": 340, "priority": "core", "domain": "Bacteria",
        "ecoli": "8-355",
    },
    {
        "name": "V1-V3",                # 27F / 534R
        "fwd_name": "27F", "fwd": "AGAGTTTGATCMTGGCTCAG",
        "rev_name": "534R", "rev": "ATTACCGCGGCTGCTGG",
        "lo": 420, "hi": 520, "priority": "core", "domain": "Bacteria",
        "ecoli": "8-534",
    },
    {
        "name": "V3-V4",                # 341F / 805R (Klindworth / Illumina)
        "fwd_name": "341F", "fwd": "CCTACGGGNGGCWGCAG",
        "rev_name": "805R", "rev": "GACTACHVGGGTATCTAATCC",
        "lo": 390, "hi": 480, "priority": "core", "domain": "Bacteria",
        "ecoli": "341-805",
    },
    {
        "name": "V4",                   # 515F-Parada / 806R-Apprill (EMP; bact+arch)
        "fwd_name": "515F-Parada", "fwd": "GTGYCAGCMGCCGCGGTAA",
        "rev_name": "806R-Apprill", "rev": "GGACTACNVGGGTWTCTAAT",
        "lo": 240, "hi": 320, "priority": "core", "domain": "Bact+Arch",
        "ecoli": "515-806",
    },
    {
        "name": "V4-V5",                # 515F-Parada / 926R-Quince (study anchor; bact+arch)
        "fwd_name": "515F-Parada", "fwd": "GTGYCAGCMGCCGCGGTAA",
        "rev_name": "926R-Quince", "rev": "CCGYCAATTYMTTTRAGTTT",
        "lo": 300, "hi": 450, "priority": "core", "domain": "Bact+Arch",
        "ecoli": "515-926",
    },
    {
        "name": "V4-V5_944R",           # 515F-Parada / 944R — sensitivity variant (bact+arch)
        "fwd_name": "515F-Parada", "fwd": "GTGYCAGCMGCCGCGGTAA",
        "rev_name": "944R", "rev": "GAATTAAACCACATGCTC",
        "lo": 320, "hi": 470, "priority": "variant", "domain": "Bact+Arch",
        "ecoli": "515-944",
    },
    {
        "name": "V7-V9",                # 1100F / 1492R (was the mislabelled "V6-V8" + "V7-V9")
        "fwd_name": "1100F", "fwd": "YAACGAGCGCAACCC",
        "rev_name": "1492R", "rev": "TACGGYTACCTTGTTACGACTT",
        "lo": 340, "hi": 440, "priority": "core", "domain": "Bacteria",
        "ecoli": "1100-1507",
    },
]


def _load_region_panel():
    """Return the 8 corrected regions, preferring region_panel_final.json (single source).

    The JSON carries `insert_lo`/`insert_hi`; we map them to the `lo`/`hi` the extractor
    expects and default a missing `domain` to 'Bacteria'.  If the JSON is absent or
    unreadable we fall back to the byte-identical hard-coded panel above so the module
    still runs standalone.
    """
    try:
        panel = json.loads(REGION_PANEL_JSON.read_text())
        out = []
        for r in panel["regions"]:
            reg = dict(r)
            reg["lo"] = r.get("insert_lo", r.get("lo"))
            reg["hi"] = r.get("insert_hi", r.get("hi"))
            reg.setdefault("domain", "Bacteria")
            out.append(reg)
        # Sanity: the collapsed panel must be exactly 8 regions with no duplicate name.
        names = [r["name"] for r in out]
        if len(out) == 8 and len(set(names)) == 8 and "V6-V8" not in names:
            return out
    except Exception:
        pass
    return [dict(r) for r in _REGION_PANEL_FALLBACK]


REGION_PANEL = _load_region_panel()

# Back-compat exports (original hard-coded V4–V5 expansions; identical to panel V4-V5).
FWD_515F = iupac_to_regex("GTGYCAGCMGCCGCGGTAA")            # GTG[CT]CAGC[AC]GCCGCGGTAA
REVC_926R = iupac_to_regex(revcomp("CCGYCAATTYMTTTRAGTTT"))  # AAACT[CT]AAA[GT][AG]AATTG[AG]CGG


# ────────────────────────────────────────────────────────────────────────────
# Extractor
# ────────────────────────────────────────────────────────────────────────────
def _edits(m) -> int:
    """Total primer mismatches for a `regex` match (0 for an exact-degenerate hit).

    Fuzzy matches expose `.fuzzy_counts == (substitutions, insertions, deletions)`; their
    sum is the edit distance.  Exact (non-fuzzy) patterns have no `fuzzy_counts` → 0 edits.
    """
    fc = getattr(m, "fuzzy_counts", None)
    return int(sum(fc)) if fc else 0


def make_extractor(fwd=FWD_515F, revc=REVC_926R, emax=3, lo=300, hi=450, include_primers=True):
    """Build an in-silico-PCR extractor for one (forward-regex, revcomp-reverse-regex) pair.

    `fwd`/`revc` are regex patterns (already IUPAC-expanded). Forward primer is matched on the
    sense strand; `revc` is the reverse-complement of the reverse primer, also matched on the
    sense strand (downstream of the forward primer).

    Returns a function ``seq -> (amplicon, fwd_edits, rev_edits) | None``.  Matching is
    exact-degenerate first, then a fuzzy ``{e<=emax}`` fallback for divergent taxa; the per-primer
    edit distances are reported so a 3+3-mismatch "amplification" can be told apart from a perfect
    match downstream.  ``include_primers=True`` returns the full primer-to-primer amplicon; False
    returns just the inter-primer insert (legacy behaviour).
    """
    fx, rx = regex.compile(fwd), regex.compile(revc)
    ff, rf = regex.compile(f'(?:{fwd}){{e<={emax}}}'), regex.compile(f'(?:{revc}){{e<={emax}}}')

    def extract(seq: str):
        f = fx.search(seq) or ff.search(seq)            # exact-degenerate first, fuzzy fallback
        if not f:
            return None
        r = rx.search(seq, f.end()) or rf.search(seq, f.end())
        if not r:
            return None
        insert_len = r.start() - f.end()                # primers excluded → amplification sanity
        if not (lo <= insert_len <= hi):
            return None
        amp = seq[f.start():r.end()] if include_primers else seq[f.end():r.start()]
        return amp, _edits(f), _edits(r)
    return extract


def build_panel_extractors(panel=REGION_PANEL, emax=3, include_primers=True):
    """Compile one extractor per panel region. Returns [(region_dict, extract_fn), ...]."""
    out = []
    for reg in panel:
        fwd_rx = iupac_to_regex(reg["fwd"])
        revc_rx = iupac_to_regex(revcomp(reg["rev"]))
        fn = make_extractor(fwd_rx, revc_rx, emax=emax, lo=reg["lo"], hi=reg["hi"],
                            include_primers=include_primers)
        out.append((reg, fn))
    return out


# ────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ────────────────────────────────────────────────────────────────────────────
def load_sequences(path: Path) -> dict:
    """Load {header: sequence}. Accepts a JSON dict or a FASTA file (.fa/.fasta/.fna)."""
    if path.suffix.lower() in (".fa", ".fasta", ".fna", ".ffn"):
        seqs, hdr, buf = {}, None, []
        for line in open(path):
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    seqs[hdr] = "".join(buf)
                hdr, buf = line[1:], []
            elif hdr is not None:
                buf.append(line)
        if hdr is not None:
            seqs[hdr] = "".join(buf)
        return seqs
    return json.load(open(path))


# `fig|1055102.5.rna.10|QD227_19370| SSU rRNA ...` → `fig|1055102.5.rna.10`
# Falls back to the first whitespace-delimited token for non-BV-BRC headers.
_FIG = re.compile(r'^(fig\|[^|]+)')


def genome_id_from_header(header: str) -> str:
    m = _FIG.match(header)
    if m:
        return m.group(1)
    tok = header.split()[0] if header.split() else header
    # BV-BRC validation headers are ">{genome_id}|{md5}" — keep only the genome_id.
    return tok.split("|", 1)[0]


# Optional {genome_id: domain} map (e.g. from the selection manifest); when absent, the
# domain is parsed from the header text or recorded as 'Unknown'.  Populated by run_panel
# from --domain-map; the authoritative domain stratification still happens later via
# lineage_for(taxon_id) in the scoring scripts (DESIGN §6).
_DOMAIN_MAP: dict = {}
_DOMAIN_TOKEN = re.compile(r"\b(Bacteria|Archaea)\b", re.I)


def domain_from_header(header: str) -> str:
    """Best-effort source domain for a 16S header: map lookup, then header token, else 'Unknown'."""
    gid = genome_id_from_header(header)
    if gid in _DOMAIN_MAP:
        return _DOMAIN_MAP[gid]
    m = _DOMAIN_TOKEN.search(header)
    if m:
        return m.group(1).capitalize()
    return "Unknown"


def _load_domain_map(path) -> dict:
    """Build {genome_id: domain} from a selection-manifest JSON ({'genomes':[{genome_id,domain}]}).

    Accepts either that manifest shape or a flat {genome_id: domain} dict.  Returns {} on any
    failure (the per-domain split then degrades gracefully to header-parsed / 'Unknown').
    """
    try:
        d = json.loads(Path(path).read_text())
    except Exception:
        return {}
    out = {}
    if isinstance(d, dict) and "genomes" in d:
        for g in d["genomes"]:
            gid, dom = g.get("genome_id"), g.get("domain")
            if gid and dom:
                out[gid] = dom
    elif isinstance(d, dict):
        out = {k: v for k, v in d.items() if isinstance(v, str)}
    return out


def _nonempty(path) -> bool:
    """True iff `path` exists and is non-empty (a usable resumability checkpoint)."""
    try:
        return os.path.getsize(path) > 0
    except OSError:
        return False


def fasta_record(rec_id: str, seq: str, width: int = 0) -> str:
    if width and width > 0:
        body = "\n".join(seq[i:i + width] for i in range(0, len(seq), width))
    else:
        body = seq
    return f">{rec_id}\n{body}\n"


# ────────────────────────────────────────────────────────────────────────────
# Modes
# ────────────────────────────────────────────────────────────────────────────
def run_panel(args):
    """Multi-region panel mode → ONE combined FASTA of {genome_id}__{region} amplicons."""
    import numpy as np
    t0 = time.time()

    # Resumability: if the combined FASTA (and its stats sidecar) already exist
    # and are non-empty, short-circuit unless --force.  Makes the extraction
    # stage re-runnable across sessions without recomputation (DESIGN §1.5).
    stats_path = _combined_stats_path(args.combined_fasta)
    if not args.force and _nonempty(args.combined_fasta) and _nonempty(stats_path):
        print(f"[pcr] {args.combined_fasta} already exists (non-empty); skipping "
              f"(use --force to regenerate). stats: {stats_path}", flush=True)
        return

    panel = REGION_PANEL
    if args.regions:
        want = {r.strip() for r in args.regions.split(",")}
        panel = [r for r in panel if r["name"] in want]
        missing = want - {r["name"] for r in panel}
        if missing:
            raise SystemExit(f"unknown region(s): {sorted(missing)}; "
                            f"available: {[r['name'] for r in REGION_PANEL]}")
    if args.core_only:
        panel = [r for r in panel if r["priority"] == "core"]

    extractors = build_panel_extractors(panel, emax=args.emax, include_primers=not args.insert_only)
    print(f"[pcr] loading {args.fasta}", flush=True)
    seqs = load_sequences(args.fasta)

    # Per-region accounting holds an overall length list plus a per-domain split (the source
    # domain of each genome's 16S — Bacteria/Archaea/Unknown), so extract_rate and the
    # amplicon-length quantiles can be reported per region AND per domain (DESIGN §4).
    per_region = {
        reg["name"]: {"hit": 0, "lens": [], "by_domain": {}}
        for reg, _ in extractors
    }
    # Count input records per source domain for the per-domain extract_rate denominator.
    n_in_by_domain: dict = {}
    n_in = 0
    n_seen_genomes = set()
    n_revcomp_rescued = 0
    out_fh = open(args.combined_fasta, "w")
    n_records = 0
    for header, raw in seqs.items():
        n_in += 1
        seq = raw.strip().upper()
        gid = genome_id_from_header(header)
        n_seen_genomes.add(gid)
        dom = domain_from_header(header)
        n_in_by_domain[dom] = n_in_by_domain.get(dom, 0) + 1
        # Both-strand pass: scan the sense strand; only if NOTHING amplifies do we retry the
        # reverse complement (catches minus-strand-deposited 16S genes — otherwise a silent
        # selection bias).  rc(seq) is computed at most once and reused across all regions.
        rc = None
        for reg, fn in extractors:
            res = fn(seq)
            if res is None:
                if rc is None:
                    rc = revcomp(seq)
                res = fn(rc)
                if res is not None:
                    n_revcomp_rescued += 1
            if res is None:
                continue  # non-amplification on either strand → skip this (genome, region)
            amp, fwd_edits, rev_edits = res
            # Provenance-native id; carry the per-primer edit counts in the FASTA description so a
            # 3+3-mismatch "amplification" is distinguishable downstream without a sidecar file.
            out_fh.write(fasta_record(
                f"{gid}__{reg['name']} fwd_edits={fwd_edits} rev_edits={rev_edits}",
                amp, args.wrap))
            n_records += 1
            rstat = per_region[reg["name"]]
            rstat["hit"] += 1
            rstat["lens"].append(len(amp))
            rstat["by_domain"].setdefault(dom, []).append(len(amp))
        if n_in % 200000 == 0:
            print(f"[pcr] {n_in:,} scanned, {n_records:,} amplicons "
                  f"({time.time()-t0:.0f}s)", flush=True)
    out_fh.close()

    def lstats(L):
        if not L:
            return None
        a = np.array(L)
        return {"n": int(a.size), "median": int(np.median(a)), "mean": round(float(a.mean()), 1),
                "p5": int(np.percentile(a, 5)), "p95": int(np.percentile(a, 95)),
                "min": int(a.min()), "max": int(a.max())}

    def domain_block(rstat):
        """Per-domain {extract_rate, amplicon_len p5/median/p95} for one region."""
        out = {}
        for dom, n_dom_in in sorted(n_in_by_domain.items()):
            lens = rstat["by_domain"].get(dom, [])
            out[dom] = {
                "n_input": n_dom_in,
                "n_amplicons": len(lens),
                "extract_rate": round(len(lens) / max(n_dom_in, 1), 4),
                "amplicon_len": lstats(lens),
            }
        return out

    stats = {
        "mode": "panel",
        "amplicon_includes_primers": not args.insert_only,
        "emax": args.emax,
        "n_input_records": n_in,
        "n_unique_genomes": len(n_seen_genomes),
        "n_amplicons_total": n_records,
        "n_revcomp_strand_rescued": n_revcomp_rescued,
        "n_input_by_domain": dict(sorted(n_in_by_domain.items())),
        "regions": {
            reg["name"]: {
                "fwd": f"{reg['fwd_name']} {reg['fwd']}",
                "rev": f"{reg['rev_name']} {reg['rev']}",
                "rev_revcomp": revcomp(reg["rev"]),
                "ecoli_coords": reg["ecoli"],
                "insert_bounds": [reg["lo"], reg["hi"]],
                "priority": reg["priority"],
                "domain": reg.get("domain", "Bacteria"),
                "n_amplicons": per_region[reg["name"]]["hit"],
                "extract_rate": round(per_region[reg["name"]]["hit"] / max(n_in, 1), 4),
                "amplicon_len": lstats(per_region[reg["name"]]["lens"]),
                "by_domain": domain_block(per_region[reg["name"]]),
            }
            for reg, _ in extractors
        },
        "wall_seconds": round(time.time() - t0, 1),
    }
    stats_path = _combined_stats_path(args.combined_fasta)
    json.dump(stats, open(stats_path, "w"), indent=2)
    print(f"[pcr] wrote {n_records:,} amplicons -> {args.combined_fasta}", flush=True)
    print(f"[pcr] wrote stats -> {stats_path}", flush=True)
    print(json.dumps(stats, indent=2))


def _combined_stats_path(combined_fasta: Path) -> Path:
    """Resolve the stats sidecar for a combined FASTA.

    The validation I/O contract names it PATHS.combined_stats
    (``combined_amplicons.stats.json``).  When the caller writes to the canonical
    PATHS.combined_fasta we honour that exact path; for any other --combined-fasta
    (e.g. the EmilyKin staging tree) we write a sibling ``<stem>.stats.json`` so
    the stats always live next to their FASTA.
    """
    if _COMMON is not None and os.path.abspath(str(combined_fasta)) == os.path.abspath(
        _COMMON.PATHS.combined_fasta
    ):
        return Path(_COMMON.PATHS.combined_stats)
    return combined_fasta.with_suffix(".stats.json")


def run_legacy(args):
    """Single-region legacy mode → {header: insert} JSON (original V4–V5 behaviour)."""
    import numpy as np
    t0 = time.time()
    if not args.force and _nonempty(args.out) and _nonempty(args.out.with_suffix(".stats.json")):
        print(f"[pcr] {args.out} already exists (non-empty); skipping "
              f"(use --force to regenerate).", flush=True)
        return
    print(f"[pcr] loading {args.fasta}", flush=True)
    fasta = load_sequences(args.fasta)
    extract = make_extractor(emax=args.emax, lo=args.lo, hi=args.hi, include_primers=False)
    out, n_in, n_full, n_full_hit, lens = {}, 0, 0, 0, []
    for header, raw in fasta.items():
        n_in += 1
        seq = raw.strip().upper()
        is_full = 1400 <= len(seq) <= 1600
        n_full += is_full
        ins = extract(seq)
        if ins:
            out[header] = ins
            lens.append(len(ins))
            n_full_hit += is_full
        if n_in % 200000 == 0:
            print(f"[pcr] {n_in:,} scanned, {len(out):,} extracted ({time.time()-t0:.0f}s)", flush=True)

    json.dump(out, open(args.out, "w"))
    L = np.array(lens)
    stats = {
        "primers": {"515F": "GTGYCAGCMGCCGCGGTAA", "926R": "CCGYCAATTYMTTTRAGTTT"},
        "region": "V4-V5 (E. coli ~533-906)",
        "n_input": n_in, "n_extracted": len(out), "extract_rate": round(len(out) / max(n_in, 1), 4),
        "n_full_length_input": n_full, "n_full_length_extracted": n_full_hit,
        "full_length_extract_rate": round(n_full_hit / max(n_full, 1), 4),
        "insert_len": {"median": int(np.median(L)), "mean": round(float(L.mean()), 1),
                       "p5": int(np.percentile(L, 5)), "p95": int(np.percentile(L, 95)),
                       "min": int(L.min()), "max": int(L.max())},
        "wall_seconds": round(time.time() - t0, 1),
    }
    json.dump(stats, open(args.out.with_suffix(".stats.json"), "w"), indent=2)
    print(json.dumps(stats, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=Path, required=True,
                    help="{header: full 16S sequence} JSON, or a FASTA file")
    ap.add_argument("--emax", type=int, default=3, help="max primer mismatches (fuzzy fallback)")
    # Legacy single-region mode:
    ap.add_argument("--out", type=Path, help="(legacy mode) {header: region insert} JSON to write")
    ap.add_argument("--lo", type=int, default=300, help="(legacy mode) min V4-V5 insert bp")
    ap.add_argument("--hi", type=int, default=450, help="(legacy mode) max V4-V5 insert bp")
    # Panel mode:
    ap.add_argument("--combined-fasta", type=Path,
                    help="(panel mode) ONE FASTA of >{genome_id}__{region} amplicons to write")
    ap.add_argument("--regions", type=str, default="",
                    help="(panel mode) comma-separated subset of region names (default: all)")
    ap.add_argument("--core-only", action="store_true",
                    help="(panel mode) drop priority='variant' regions (944R, V7-V9)")
    ap.add_argument("--insert-only", action="store_true",
                    help="(panel mode) emit inter-primer insert instead of primer-inclusive amplicon")
    ap.add_argument("--wrap", type=int, default=0,
                    help="(panel mode) FASTA line width (0 = one-line sequences)")
    ap.add_argument("--domain-map", type=Path, default=None,
                    help="(panel mode) selection-manifest JSON ({'genomes':[{genome_id,domain}]}) "
                         "or flat {genome_id: domain}; populates the per-domain stats split. "
                         "Defaults to PATHS.selection_json when present.")
    ap.add_argument("--force", action="store_true",
                    help="regenerate even if the output (and its stats sidecar) already exist")
    args = ap.parse_args()

    if args.combined_fasta:
        # Wire the optional {genome_id: domain} map for the per-domain stats split.
        # Authoritative domain stratification still happens later via lineage_for()
        # in scoring; this is a best-effort split so combined_stats carries Bact/Arch.
        dm_path = args.domain_map
        if dm_path is None and _COMMON is not None and _nonempty(_COMMON.PATHS.selection_json):
            dm_path = Path(_COMMON.PATHS.selection_json)
        if dm_path is not None:
            _DOMAIN_MAP.update(_load_domain_map(dm_path))
            if _DOMAIN_MAP:
                print(f"[pcr] domain map: {len(_DOMAIN_MAP):,} genome_id->domain "
                      f"from {dm_path}", flush=True)
        run_panel(args)
    elif args.out:
        run_legacy(args)
    else:
        ap.error("supply --combined-fasta (panel mode) or --out (legacy single-region mode)")


if __name__ == "__main__":
    main()
