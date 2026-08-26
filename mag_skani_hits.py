#!/usr/bin/env python
"""
mag_skani_hits.py — whole-genome (MAG) reference hits for prFBA via skani ANI,
against one or more sketched reference databases (GTDB in particular).

This is the MAG analog of the 16S ASV hits stage (`align_hits.py` /
`hit_amplicons.py`): where an ASV is referenced into BV-BRC by 16S identity, a
MAG is referenced into a genome database by skani average nucleotide identity.
The output JSON is schema-compatible with `asv_top20_alignment_hits.json`, so
`select_references.py` can consume it unchanged (see ANI_KNOBS below for the
ANI-recalibrated thresholds — the 16S Yarza tiers do not apply to genome ANI).

Databases
---------
* GTDB       : the pre-sketched skani database from the skani authors
               (e.g. skani_db/gtdb_r226/), taxonomy joined from the GTDB
               bac120+ar53 taxonomy TSVs (accession -> d__..;s__ lineage).
* BV-BRC     : any directory of BV-BRC `<genome_id>.fna` genomes sketched with
               the `sketch` subcommand (e.g. ~/Documents/sludge_genomes),
               taxonomy joined from a JSON map (genome_id -> rank dict, e.g.
               ~/Documents/sludge/gtdb_taxonomy.json). BV-BRC hits carry real
               BV-BRC genome_ids, so the downstream synthetic-genome merge
               (`select_references.py` -> `build_synthetic_genomes.py`) works
               end-to-end with the PGFam gene provider.

Best-hit rule (AF-aware, from the EmilyKin server handoff)
----------------------------------------------------------
skani ranks hits by ANI, which lets a spurious high-ANI / near-zero-AF hit
(one reference did this to ~24 EmilyKin MAGs) outrank the true genome-wide
match. The "best hit" is therefore the highest-ANI hit with
Align_fraction_query >= --min-af-select (default 15, skani's own --min-af
default); if no hit clears it, the MAG is `low_AF_unresolved`.

ANI novelty classes (collaborator wording from Table S6):
  >=95            same species as closest reference
  85-95           putatively species-distinct from closest reference
  <85             stronger divergence, no close ANI-level match
  low_AF_unresolved  ANI unreliable because too little genome aligned
  no_hit          below skani's ~80% ANI detection floor

Usage
-----
  # 1. sketch a FASTA-directory database (once per DB; GTDB comes pre-sketched)
  python mag_skani_hits.py sketch --fasta-dir ~/Documents/sludge_genomes \
      --out skani_db/bvbrc_sludge

  # 2. search MAGs against a database, join taxonomy, emit hits + summary
  python mag_skani_hits.py search \
      --mags ~/Documents/EmilyKin/dereplicated_genomes \
      --db skani_db/skani_gtdb_r226-v0.3 --db-name gtdb_r226 \
      --taxonomy skani_db/gtdb_r226/bac120_taxonomy_r226.tsv.gz \
      --taxonomy skani_db/gtdb_r226/ar53_taxonomy_r226.tsv.gz \
      --gtdbtk ~/Documents/EmilyKin/gtdbtk.bac120.summary.tsv \
      --outdir ~/Documents/EmilyKin/mag_skani_hits --topk 20

  # 3. per-rank concordance of the AF-aware best hit vs the GTDB-Tk call
  python mag_skani_hits.py concordance \
      --hits ~/Documents/EmilyKin/mag_skani_hits/mag_top20_skani_hits.gtdb_r226.json \
      --gtdbtk ~/Documents/EmilyKin/gtdbtk.bac120.summary.tsv \
      --outdir ~/Documents/EmilyKin/mag_skani_hits

Feeding select_references.py (prFBA integration)
------------------------------------------------
The hits JSON doubles as the selector input:  `identity` = ANI/100,
`align_score` = ANI, `aligned_len` = AF_query/100 * genome length (so the
selector's coverage = AF_query/100), `edlib_identity` = None (agreement gate
passes). Recalibrate the 16S knobs for genome ANI:

    from mag_skani_hits import ANI_KNOBS
    select_representatives(record, knobs=ANI_KNOBS, gene_provider=...)

    python select_references.py --hits mag_top20_skani_hits.bvbrc_sludge.json \
        --out mag_reference_selection.json --knobs-json ani_knobs.json
"""

import argparse
import csv
import glob
import gzip
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SKANI_DEFAULT = os.path.join(HERE, "bin", "skani")

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
GTDB_PREFIX = dict(zip(["d__", "p__", "c__", "o__", "f__", "g__", "s__"], RANKS))
# lineage keys used by the ASV hits schema (select_references reads Family/Genus/Species)
LINEAGE_KEYS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")

# select_references.py DEFAULT_KNOBS are 16S-calibrated (Yarza tiers). For genome
# ANI: species boundary 95% ANI; below ~83% skani ANI is at its detection floor and
# rank inference is unreliable, so the floor abstains there. Coverage for a MAG is
# AF_query/100, so cov_min is the AF gate (0.15), not the 16S 0.90.
ANI_KNOBS = dict(
    cov_min=0.15,            # reliability gate = AF_query >= 15%
    t_species=0.95,          # ANI species boundary
    t_genus=0.85,            # "putatively species-distinct" band lower edge
    t_family=0.83,           # skani reliable-detection floor
    family_floor=0.83,       # abstain below the detection floor
    band=0.005,              # 0.5 ANI points = an ANI tie
    ov_species=0.90, ov_genus=0.60, ov_family=0.45,  # gene-overlap estimator, genome scale
)


# ----------------------------------------------------------------------- utils
def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def log(msg):
    print(f"[mag_skani] {msg}", file=sys.stderr, flush=True)


def find_skani(arg):
    for cand in ([arg] if arg else []) + [SKANI_DEFAULT, "skani"]:
        try:
            subprocess.run([cand, "--version"], capture_output=True, check=True)
            return cand
        except (OSError, subprocess.CalledProcessError):
            continue
    sys.exit("skani binary not found (tried --skani, prFBA/bin/skani, PATH)")


def mag_files(mags):
    if os.path.isdir(mags):
        files = sorted(sum((glob.glob(os.path.join(mags, e))
                            for e in ("*.fa", "*.fna", "*.fasta")), []))
    else:
        files = sorted(glob.glob(mags))
    if not files:
        sys.exit(f"no FASTA files found under {mags}")
    return files


def stem(path):
    return re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", os.path.basename(path))


def fasta_stats(path):
    """(total_len, n_contigs) — buffered byte scan, no per-record parsing."""
    n, total = 0, 0
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
            else:
                total += len(line.strip())
    return total, n


# ------------------------------------------------------------------- taxonomy
def load_taxonomy_tsv(paths):
    """GTDB taxonomy TSV(s): 'GB_GCA_...|RS_GCF_...  d__..;s__..' -> bare accession -> lineage."""
    tax = {}
    for path in paths:
        with _open(path) as fh:
            for line in fh:
                key, lineage = line.rstrip("\n").split("\t", 1)
                tax[re.sub(r"^(GB|RS)_", "", key)] = lineage
    return tax


def load_taxonomy_json(path):
    """JSON map genome_id (or 'genome_id__suffix') -> rank dict -> genome_id -> lineage string."""
    raw = json.load(open(path))
    tax = {}
    for key, ranks in raw.items():
        gid = key.split("__", 1)[0]
        parts = []
        for pref, rank in GTDB_PREFIX.items():
            val = ranks.get(rank.capitalize()) or ranks.get(rank) or ""
            if rank == "domain" and not val:
                val = ranks.get("Domain") or ranks.get("Kingdom") or ""
            parts.append(pref + str(val))
        tax[gid] = ";".join(parts)
    return tax


def split_lineage(lineage):
    """GTDB lineage string -> {'domain': 'Bacteria', ...} (values without prefixes)."""
    out = {r: "" for r in RANKS}
    for token in (lineage or "").split(";"):
        token = token.strip()
        for pref, rank in GTDB_PREFIX.items():
            if token.startswith(pref):
                out[rank] = token[len(pref):]
    return out


def compat_lineage(ranks):
    """rank dict -> ASV-schema lineage {'Kingdom':...,'Species':...} for select_references."""
    return dict(zip(LINEAGE_KEYS,
                    (ranks[r] for r in RANKS)))


def ref_id(ref_file, db_kind):
    """gtdb: GCA_020442275.1 from .../GCA/020/442/275/GCA_020442275.1_genomic.fna.gz
       bvbrc: 1002339.3 from .../1002339.3.fna"""
    base = os.path.basename(ref_file)
    if db_kind == "gtdb":
        m = ACC_RE.search(base)
        if m:
            return m.group(1)
    return re.sub(r"(_genomic)?\.(fa|fna|fasta)(\.gz)?$", "", base)


def load_gtdbtk(path):
    """gtdbtk.bac120.summary.tsv -> mag -> {'classification':..., 'closest': accession}."""
    out = {}
    with _open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            closest = row.get("closest_genome_reference") or ""
            out[row["user_genome"]] = dict(
                classification=row.get("classification") or "",
                closest="" if closest in ("N/A", "") else re.sub(r"^(GB|RS)_", "", closest))
    return out


# ------------------------------------------------------------------ sketch cmd
def cmd_sketch(a):
    skani = find_skani(a.skani)
    files = mag_files(a.fasta_dir)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    lst = a.out.rstrip("/") + ".files.txt"
    with open(lst, "w") as fh:
        fh.write("\n".join(files) + "\n")
    log(f"sketching {len(files)} genomes -> {a.out}")
    t0 = time.time()
    subprocess.run([skani, "sketch", "-l", lst, "-o", a.out, "-t", str(a.threads)],
                   check=True)
    log(f"sketch done in {time.time()-t0:.0f}s")


# ------------------------------------------------------------------ search cmd
def classify(best, n_hits, af_min):
    if n_hits == 0:
        return "no_hit", "no skani hit (below the ~80% ANI detection floor)"
    if not best["af_pass"]:
        return "low_AF_unresolved", "ANI unreliable because too little genome aligned"
    ani = best["ani"]
    if ani >= 95:
        return ">=95", "same species as closest %s representative" % best["_dbword"]
    if ani >= 85:
        return "85-95", "putatively species-distinct from closest %s representative" % best["_dbword"]
    return "<85", "stronger divergence; no close ANI-level match"


def best_hit(hits, af_min):
    """Highest-ANI hit with AF_query >= af_min; else highest-ANI overall (flagged)."""
    good = [h for h in hits if h["af_query"] >= af_min]
    pool = good if good else hits
    b = dict(max(pool, key=lambda h: h["ani"]))
    b["af_pass"] = bool(good)
    return b


def cmd_search(a):
    skani = find_skani(a.skani)
    mags = mag_files(a.mags)
    os.makedirs(a.outdir, exist_ok=True)
    tag = a.db_name
    raw_tsv = os.path.join(a.outdir, f"mag_vs_{tag}_skani_raw.tsv")
    err_log = os.path.join(a.outdir, f"skani_search.{tag}.log")

    # taxonomy
    tax = {}
    for tpath in a.taxonomy or []:
        loader = load_taxonomy_json if tpath.endswith(".json") else lambda p: load_taxonomy_tsv([p])
        tax.update(loader(tpath))
    log(f"taxonomy records: {len(tax)}")
    gtdbtk = load_gtdbtk(a.gtdbtk) if a.gtdbtk else {}

    # skani search (skip if raw tsv already present and --resume)
    t0 = time.time()
    if not (a.resume and os.path.exists(raw_tsv)):
        lst = os.path.join(a.outdir, f"query_mags.{tag}.txt")
        with open(lst, "w") as fh:
            fh.write("\n".join(mags) + "\n")
        cmd = [skani, "search", "-d", a.db, "--ql", lst, "-o", raw_tsv,
               "-t", str(a.threads), "-n", str(a.topk), "--min-af", str(a.min_af_search)]
        log("running: " + " ".join(cmd))
        with open(err_log, "w") as eh:
            subprocess.run(cmd, check=True, stderr=eh)
    search_s = round(time.time() - t0, 1)
    stderr_txt = open(err_log).read() if os.path.exists(err_log) else ""
    learned = "Learned ANI mode" in stderr_txt

    # parse + join
    hits = {}
    n_rows = 0
    with open(raw_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n_rows += 1
            rid = ref_id(row["Ref_file"], a.db_kind)
            lineage = tax.get(rid, "")
            ranks = split_lineage(lineage)
            hits.setdefault(stem(row["Query_file"]), []).append(dict(
                ani=float(row["ANI"]),
                af_ref=float(row["Align_fraction_ref"]),
                af_query=float(row["Align_fraction_query"]),
                ref_id=rid, ref_file=row["Ref_file"],
                ref_name=row.get("Ref_name", ""),
                lineage_gtdb=lineage, ranks=ranks, _dbword=a.db_kind.upper()
                if a.db_kind == "gtdb" else a.db_kind))
    log(f"{n_rows} hit rows over {len(hits)}/{len(mags)} MAGs "
        f"({search_s}s{', Learned-ANI mode' if learned else ''})")

    # per-MAG records (ASV-schema compatible) + summary rows
    records, summary = {}, []
    for path in mags:
        mag = stem(path)
        mlen, ncontig = fasta_stats(path)
        mh = sorted(hits.get(mag, []), key=lambda h: -h["ani"])
        b = best_hit(mh, a.min_af_select) if mh else None
        ani_class, interp = classify(b or {}, len(mh), a.min_af_select)
        gt = gtdbtk.get(mag, {})
        gt_ranks = split_lineage(gt.get("classification", ""))

        top = []
        for i, h in enumerate(mh, 1):
            org = (h["ranks"]["species"] or h["ranks"]["genus"]
                   or h["ref_name"].split(",")[0])
            top.append(dict(
                rank=i,
                # --- select_references.py compatible fields ---
                identity=round(h["ani"] / 100.0, 4),
                align_score=h["ani"],
                aligned_len=int(round(h["af_query"] / 100.0 * mlen)),
                edlib_identity=None, edlib_distance=None,
                n_matches=None, md5=None, feature_id=None, taxon_id=None,
                genome_id=h["ref_id"], organism=org,
                lineage=compat_lineage(h["ranks"]),
                # --- skani-native fields ---
                ani=h["ani"], af_ref=h["af_ref"], af_query=h["af_query"],
                af_pass=h["af_query"] >= a.min_af_select,
                ref_file=h["ref_file"], ref_name=h["ref_name"],
                lineage_gtdb=h["lineage_gtdb"],
            ))

        best_ident = round(b["ani"] / 100.0, 4) if (b and b["af_pass"]) else None
        records[mag] = dict(
            asv_len=mlen,                       # selector-required key (= genome bp)
            mag_len=mlen, n_contigs=ncontig,
            best_identity=best_ident,
            best_align_score=b["ani"] if b else None,
            best_ani=b["ani"] if b else None,
            best_af_query=b["af_query"] if b else None,
            best_af_ref=b["af_ref"] if b else None,
            best_ref=b["ref_id"] if b else None,
            af_pass=b["af_pass"] if b else False,
            ani_class=ani_class, novelty_interpretation=interp,
            midas_taxonomy=" ".join(v for v in gt_ranks.values() if v) or None,
            gtdbtk_taxonomy=gt.get("classification") or None,
            gtdbtk_closest=gt.get("closest") or None,
            db=tag, n_hits=len(mh),
            top20=top,
        )
        summary.append(dict(
            mag=mag, mag_len=mlen, n_contigs=ncontig, n_hits=len(mh),
            best_ref=b["ref_id"] if b else "",
            best_ref_species=b["ranks"]["species"] if b else "",
            best_ani=b["ani"] if b else "",
            best_af_query=b["af_query"] if b else "",
            best_af_ref=b["af_ref"] if b else "",
            af_pass=int(b["af_pass"]) if b else 0,
            ani_class=ani_class, novelty_interpretation=interp,
            gtdbtk_species=gt_ranks["species"], gtdbtk_genus=gt_ranks["genus"],
            gtdbtk_closest=gt.get("closest", ""),
            same_ref_as_gtdbtk=int(bool(b and gt.get("closest")
                                        and b["ref_id"] == gt["closest"])),
        ))

    hits_json = os.path.join(a.outdir, f"mag_top{a.topk}_skani_hits.{tag}.json")
    with open(hits_json, "w") as fh:
        json.dump(records, fh, indent=1)
    sum_csv = os.path.join(a.outdir, f"mag_skani_summary.{tag}.csv")
    with open(sum_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    classes = {}
    for s in summary:
        classes[s["ani_class"]] = classes.get(s["ani_class"], 0) + 1
    stats = dict(db=tag, db_path=a.db, db_kind=a.db_kind, n_mags=len(mags),
                 n_hit_rows=n_rows, n_mags_with_hits=len(hits), topk=a.topk,
                 min_af_search=a.min_af_search, min_af_select=a.min_af_select,
                 learned_ani_mode=learned, search_seconds=search_s,
                 ani_class_counts=classes, taxonomy_records=len(tax),
                 skani=skani, outputs=[hits_json, sum_csv, raw_tsv])
    with open(os.path.join(a.outdir, f"run_stats.{tag}.json"), "w") as fh:
        json.dump(stats, fh, indent=1)
    log(f"ani_class: {classes}")
    log(f"wrote {hits_json}\nwrote {sum_csv}")
    if a.gtdbtk:
        do_concordance(hits_json, a.gtdbtk, a.outdir, tag)


# ------------------------------------------------------------- concordance cmd
def do_concordance(hits_json, gtdbtk_path, outdir, tag):
    """Per-rank agreement of the AF-aware best hit's lineage vs the GTDB-Tk call."""
    records = json.load(open(hits_json))
    gtdbtk = load_gtdbtk(gtdbtk_path)
    rows, agg = [], {r: dict(compared=0, agree=0, mag_blank=0, ref_blank=0)
                     for r in RANKS}
    n_af_pass = 0
    for mag, rec in records.items():
        gt = split_lineage(gtdbtk.get(mag, {}).get("classification", ""))
        best = next((h for h in rec["top20"] if h["af_pass"]), None) \
            if rec.get("af_pass") else None
        row = dict(mag=mag, ani_class=rec["ani_class"],
                   best_ref=rec.get("best_ref") or "",
                   best_ani=rec.get("best_ani") or "")
        if best:
            n_af_pass += 1
        ref = split_lineage(best["lineage_gtdb"]) if best else {r: "" for r in RANKS}
        for r in RANKS:
            g, f = gt[r], ref[r]
            if not best:
                row[r] = ""
                continue
            if not g:
                agg[r]["mag_blank"] += 1
                row[r] = "gtdbtk_blank"
            elif not f:
                agg[r]["ref_blank"] += 1
                row[r] = "ref_blank"
            else:
                agg[r]["compared"] += 1
                ok = g == f
                agg[r]["agree"] += ok
                row[r] = "agree" if ok else f"{g}|{f}"
        rows.append(row)

    for r in RANKS:
        c = agg[r]
        c["fraction"] = round(c["agree"] / c["compared"], 4) if c["compared"] else None
    out_csv = os.path.join(outdir, f"mag_skani_concordance.{tag}.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["mag", "ani_class", "best_ref", "best_ani"] + RANKS)
        w.writeheader()
        w.writerows(rows)
    out_json = os.path.join(outdir, f"concordance_by_rank.{tag}.json")
    with open(out_json, "w") as fh:
        json.dump(dict(n_mags=len(records), n_af_pass_best=n_af_pass,
                       by_rank=agg), fh, indent=1)
    log("concordance (agree/compared): " +
        ", ".join(f"{r}={agg[r]['agree']}/{agg[r]['compared']}" for r in RANKS))
    log(f"wrote {out_csv}\nwrote {out_json}")


def cmd_concordance(a):
    tag = re.sub(r"^mag_top\d+_skani_hits\.|\.json$", "",
                 os.path.basename(a.hits))
    do_concordance(a.hits, a.gtdbtk, a.outdir, tag)


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skani", help="skani binary (default prFBA/bin/skani, then PATH)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("sketch", help="sketch a FASTA directory into a skani DB")
    s.add_argument("--fasta-dir", required=True)
    s.add_argument("--out", required=True, help="output sketch directory")
    s.add_argument("--threads", type=int, default=os.cpu_count())
    s.set_defaults(fn=cmd_sketch)

    s = sub.add_parser("search", help="search MAGs against a sketched DB")
    s.add_argument("--mags", required=True, help="MAG FASTA dir or glob")
    s.add_argument("--db", required=True, help="skani sketch directory")
    s.add_argument("--db-name", required=True, help="tag for output filenames, e.g. gtdb_r226")
    s.add_argument("--db-kind", choices=["gtdb", "bvbrc"], default="gtdb",
                   help="how to parse reference ids from paths")
    s.add_argument("--taxonomy", action="append",
                   help="GTDB taxonomy tsv(.gz) (repeatable) or JSON rank-map")
    s.add_argument("--gtdbtk", help="gtdbtk summary tsv (adds context + auto-concordance)")
    s.add_argument("--outdir", required=True)
    s.add_argument("--topk", type=int, default=20)
    s.add_argument("--min-af-search", type=float, default=0.0,
                   help="skani --min-af during search (keep 0; filter at selection)")
    s.add_argument("--min-af-select", type=float, default=15.0,
                   help="AF_query %% floor for the AF-aware best hit")
    s.add_argument("--threads", type=int, default=os.cpu_count())
    s.add_argument("--resume", action="store_true",
                   help="reuse an existing raw tsv instead of re-searching")
    s.set_defaults(fn=cmd_search)

    s = sub.add_parser("concordance", help="per-rank best-hit vs GTDB-Tk agreement")
    s.add_argument("--hits", required=True, help="mag_top*_skani_hits.*.json")
    s.add_argument("--gtdbtk", required=True)
    s.add_argument("--outdir", required=True)
    s.set_defaults(fn=cmd_concordance)

    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
