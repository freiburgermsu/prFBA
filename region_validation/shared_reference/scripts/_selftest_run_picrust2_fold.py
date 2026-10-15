#!/usr/bin/env python
"""Self-test for ``run_picrust2_fold.py`` -- NOT a spec component.

``run_picrust2_fold.py`` consumes ``build_reference.py``'s reference directories and
``plan_folds.py``'s fold plan, neither of which exists yet when it is written.  This
script synthesises the smallest honest stand-ins (a 40-tip bacterial reference, a
25-tip archaeal one, a 160-tip stock-aligned recipe-control reference, a 5-genome
fold plan) and then runs every code path, including the ones that must FAIL:

  1  --mode fold                 happy path, one query poisoned so --min_align drops it
  2  resume                      a second run recomputes nothing
  3  C3 fires (exit 3)           same reference with the g_ prefix stripped; the
                                 trailing-zero genome ids are exactly the 3.61% that
                                 castor_nsti.R's read.table coerces away
  4  C4 fires (exit 3)           one trait row renamed to a non-tip
  5  --mode recipe_control       stock ec.txt.gz: 26,868 rows -> 160 tips (castor's
                                 num_known/nrow >= 0.1 guard) + the §5.4 column screen
  6  --routing nsti              both domains, sibling auto-derivation, argmin-NSTI
                                 routing, union columns with 0-fill
  7  include_self / plumbing     mode labels resolve
  8  usage guards (exit 2)
  9  compare.load_picrust2       reads the two combined files verbatim and agrees with
                                 predictions.json on every EC set, NSTI and label

Everything is written under ``--workdir`` (default /tmp/sr_fix_selftest); nothing in
the repo is touched.  Run:

  ~/Documents/py_venv/bin/python _selftest_run_picrust2_fold.py
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import random
import subprocess
import sys

import orjson

HERE = os.path.dirname(os.path.abspath(__file__))
RUNNER = os.path.join(HERE, "run_picrust2_fold.py")
RV = os.path.abspath(os.path.join(HERE, "..", ".."))
PY = os.path.expanduser("~/Documents/py_venv/bin/python")
MD5SEQ = ("/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/"
          "16S_md5_seq.json")
STOCK = os.path.expanduser(
    "~/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files")
SEED = 1729
FAILURES = []


def p2(cmd: str) -> None:
    subprocess.run(["bash", "-lc", f"source ~/opt/picrust2/activate.sh && {cmd}"],
                   check=True)


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: got {got!r}, want {want!r}")
    if not ok:
        FAILURES.append(name)


def run(args, expect: int) -> str:
    out = subprocess.run([PY, RUNNER] + args, capture_output=True, text=True)
    tail = (out.stdout + out.stderr).strip().splitlines()
    print(f"  $ run_picrust2_fold.py {' '.join(args[:4])} ... -> exit {out.returncode}")
    if tail:
        print(f"    {tail[-1][:160]}")
    check(f"exit code for {args[1]}/{args[args.index('--routing') + 1] if '--routing' in args else 'forced'}",
          out.returncode, expect)
    return "\n".join(tail)


# --------------------------------------------------------------------------- #
# fixture construction
# --------------------------------------------------------------------------- #
def match_column_align(unaligned: str, hmm: str, labels, sto: str):
    """hmmalign -> keep only #=GC RF match columns (SPEC §5.1 steps 2-3)."""
    p2(f"hmmalign --trim --dna --outformat Stockholm -o {sto} {hmm} {unaligned}")
    rows, rf = {}, ""
    for line in open(sto):
        line = line.rstrip()
        if line.startswith("#=GC RF"):
            rf += line.split()[2]
            continue
        if not line or line[0] == "#" or line[:2] == "//":
            continue
        k, v = line.split(None, 1)
        rows[k] = rows.get(k, "") + v.strip()
    keep = [i for i, c in enumerate(rf) if c != "."]
    out, gaps = {}, []
    for lab in labels:
        s = "".join(rows[lab][i] for i in keep).upper().replace(".", "-")
        gaps.append(s.count("-") / len(s))
        out[lab] = s
    mean_gap = sum(gaps) / len(gaps)
    print(f"    MSA width {len(keep)} mean gap {mean_gap:.4f}")
    assert len(keep) < 2000 and mean_gap < 0.10, "fixture MSA geometry is wrong"
    return out, len(keep)


def write_reference(rd, name, aligned, width, stock_model, trait_cols, ec_of, domain):
    os.makedirs(rd, exist_ok=True)
    with open(f"{rd}/{name}.fna", "w") as fh:
        for lab, s in aligned.items():
            fh.write(f">{lab}\n{s}\n")
    p2(f"hmmbuild --dna --informat afa {rd}/{name}.hmm {rd}/{name}.fna > /dev/null")
    p2(f"FastTree -nt -gtr -gamma {rd}/{name}.fna > {rd}/{name}.tre "
       f"2> {rd}/{name}.fasttree.stderr")
    # fixture .model: stock GTR string, correct width only. Deriving the rates from
    # the FastTree log is build_reference.py's §5.1 step 7, not this component's.
    open(f"{rd}/{name}.model", "w").write(
        stock_model.split(", noname")[0] + f", noname = 1-{width}\n")
    labs = list(aligned)
    # the §5.4 crash-column screen (build_reference.py owns the real one): castor's
    # mp dies with `invalid 'ncol' value` on columns that are constant across tips
    prev = {e: sum(1 for l in labs if e in ec_of[l]) for e in trait_cols}
    cols = [e for e in trait_cols if 0 < prev[e] < len(labs)]
    with gzip.open(f"{rd}/traits_EC.txt.gz", "wt", newline="") as fh:
        fh.write("assembly\t" + "\t".join(cols) + "\n")
        for l in labs:
            fh.write(f"{l}\t" + "\t".join("1" if e in ec_of[l] else "0"
                                          for e in cols) + "\n")
    json.dump({"name": name, "mode": "fold", "fold": 0, "domain": domain,
               "density": 1.0, "n_tips": len(labs), "msa_width": width,
               "trait_cols": len(cols), "set_equality_ok": True,
               "frozen": {"seed": SEED, "label_prefix": "g_",
                          "picrust2": {"hsp_method": "mp", "edge_exponent": 0.5,
                                       "min_align": 0.80,
                                       "placement_tool": "epa-ng",
                                       "calculate_nsti": True}}},
              open(f"{rd}/manifest.json", "w"), indent=1)
    print(f"    {name}: {len(labs)} tips, {len(cols)}/{len(trait_cols)} trait columns")
    return cols


def build_fixture(W):
    rng = random.Random(SEED)
    sel = orjson.loads(open(f"{RV}/data/selection_10k.json", "rb").read())["genomes"]
    funcs = orjson.loads(open(
        f"{RV}/picrust2_comparison/data/genome_function_sets.json", "rb").read())

    def pool_for(domain, n):
        cand = [g for g in sel if g["domain"] == domain
                and 1200 <= g["best_16s_len"] <= 1800
                and g["genome_id"] in funcs and funcs[g["genome_id"]]["ec"]]
        seen, pool = set(), []
        for g in sorted(cand, key=lambda r: r["genome_id"]):
            if g.get("genus") and g["genus"] in seen:
                continue
            seen.add(g.get("genus"))
            pool.append(g)
        return pool

    bac_pool = pool_for("Bacteria", 45)
    # C3 needs a reference id R's read.table would coerce: one ending in a zero
    tz = next(g for g in bac_pool
              if g["genome_id"].replace(".", "").isdigit()
              and g["genome_id"].split(".")[1].endswith("0"))
    chosen = [tz] + [g for g in rng.sample(bac_pool, 45) if g is not tz]
    ref_g, test_g = chosen[:40], chosen[40:45]
    arc_g = random.Random(SEED).sample(pool_for("Archaea", 25), 25)

    db = orjson.loads(open(MD5SEQ, "rb").read())
    seq = {g["best_16s_md5"]: db[g["best_16s_md5"]].upper()
           for g in chosen + arc_g}
    del db
    ec_of = {f"g_{g['genome_id']}": {e for e in funcs[g["genome_id"]]["ec"]
                                     if "-" not in e} for g in chosen + arc_g}

    # ---- bacterial fold reference ---------------------------------------- #
    print("  building the bacterial fold reference")
    un = f"{W}/panel_bac.fna"
    open(un, "w").write("".join(f">g_{g['genome_id']}\n{seq[g['best_16s_md5']]}\n"
                                for g in ref_g))
    labs = [f"g_{g['genome_id']}" for g in ref_g]
    aligned, width = match_column_align(
        un, f"{STOCK}/bacteria/bac_ref/bac_ref.hmm", labs, f"{W}/panel_bac.sto")
    prev = {}
    for l in labs:
        for e in ec_of[l]:
            prev[e] = prev.get(e, 0) + 1
    usable = sorted(e for e, c in prev.items() if 2 <= c <= len(labs) - 2)
    universe = sorted(rng.sample(usable, min(300, len(usable))))
    stock_bac_model = open(
        f"{STOCK}/bacteria/bac_ref/bac_ref.model").read().strip()
    write_reference(f"{W}/refs/f00_bac_d100", "f00_bac_d100", aligned, width,
                    stock_bac_model, universe, ec_of, "bac")

    # ---- archaeal fold reference (same EC universe: SPEC §5.3) ------------ #
    print("  building the archaeal fold reference")
    un = f"{W}/panel_arc.fna"
    open(un, "w").write("".join(f">g_{g['genome_id']}\n{seq[g['best_16s_md5']]}\n"
                                for g in arc_g))
    alabs = [f"g_{g['genome_id']}" for g in arc_g]
    aaln, awidth = match_column_align(
        un, f"{STOCK}/archaea/arc_ref/arc_ref.hmm", alabs, f"{W}/panel_arc.sto")
    write_reference(f"{W}/refs/f00_arc_d100", "f00_arc_d100", aaln, awidth,
                    open(f"{STOCK}/archaea/arc_ref/arc_ref.model").read().strip(),
                    universe, ec_of, "arc")

    # ---- fold plan + panel ------------------------------------------------ #
    truth = orjson.loads(open(f"{RV}/data/benchmark_truth.json", "rb").read())
    fanout = orjson.loads(open(f"{RV}/data/amplicon_md5_fanout.json", "rb").read())
    cell_md5 = {c: m for m, cs in fanout.items() for c in cs}
    cells, md5s = [], []
    for g in test_g:
        for reg in ("V4", "V3-V4", "V1-V2"):
            c = f"{g['genome_id']}__{reg}"
            if c in truth and c in cell_md5:
                cells.append(c)
                md5s.append(cell_md5[c])
    md5s = sorted(set(md5s))
    frozen = {"seed": SEED, "label_prefix": "g_", "k_folds": 3,
              "tau_amplicon": 0.987,
              "picrust2": {"hsp_method": "mp", "edge_exponent": 0.5,
                           "min_align": 0.80, "placement_tool": "epa-ng",
                           "calculate_nsti": True}}
    json.dump({"meta": {"frozen": frozen, "rank_blocked": "genus", "k": 3,
                        "n_eligible": len(chosen),
                        "note": "SYNTHETIC SELF-TEST FIXTURE, not a real fold plan"},
               "folds": [{"fold": 0,
                          "test_genomes": [g["genome_id"] for g in test_g],
                          "test_cells": cells, "test_md5s": md5s,
                          "reference_genomes": [g["genome_id"] for g in ref_g],
                          "n_test_genomes": len(test_g),
                          "n_reference_genomes": len(ref_g)}]},
              open(f"{W}/fold_plan.json", "w"), indent=1)
    json.dump({"meta": {"n_eligible": len(chosen) + len(arc_g)},
               "genomes": {g["genome_id"]: {"label": f"g_{g['genome_id']}",
                                            "md5": g["best_16s_md5"],
                                            "domain": g["domain"],
                                            "template": "gramneg"}
                           for g in chosen + arc_g}},
              open(f"{W}/panel_eligible.json", "w"), indent=1)

    # ---- query source, with ONE md5 poisoned so --min_align drops it ------ #
    amp, name = {}, None
    for line in open(f"{RV}/data/combined_amplicons.unique.fasta"):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            amp[name] = ""
        elif name in amp:
            amp[name] += line
    poison = md5s[0]
    with open(f"{W}/queries_poisoned.fna", "w") as fh:
        for m in md5s:
            fh.write(f">{m}\n{'GGCCGGCCATAT' * 25 if m == poison else amp[m]}\n")
    print(f"    fold plan: {len(test_g)} test genomes, {len(cells)} cells, "
          f"{len(md5s)} md5s (1 poisoned: {poison[:8]})")

    # ---- bare-label copy: C3 must fire ----------------------------------- #
    nb = f"{W}/refs/nopfx_bac_d100"
    os.makedirs(nb, exist_ok=True)
    src = f"{W}/refs/f00_bac_d100/f00_bac_d100"
    open(f"{nb}/nopfx_bac_d100.fna", "w").write(
        open(src + ".fna").read().replace(">g_", ">"))
    open(f"{nb}/nopfx_bac_d100.tre", "w").write(
        open(src + ".tre").read().replace("g_", ""))
    open(f"{nb}/nopfx_bac_d100.model", "w").write(open(src + ".model").read())
    with gzip.open(src.rsplit("/", 1)[0] + "/traits_EC.txt.gz", "rt",
                   newline="") as f, \
            gzip.open(f"{nb}/traits_EC.txt.gz", "wt", newline="") as o:
        rd, wr = csv.reader(f, delimiter="\t"), csv.writer(
            o, delimiter="\t", lineterminator="\n")
        wr.writerow(next(rd))
        for r in rd:
            r[0] = r[0][2:] if r[0].startswith("g_") else r[0]
            wr.writerow(r)
    p2(f"hmmbuild --dna --informat afa {nb}/nopfx_bac_d100.hmm "
       f"{nb}/nopfx_bac_d100.fna > /dev/null")

    # ---- corrupted trait row: C4 must fire -------------------------------- #
    bc = f"{W}/refs/badc4_bac_d100"
    os.makedirs(bc, exist_ok=True)
    for ext in ("fna", "tre", "hmm", "model"):
        open(f"{bc}/badc4_bac_d100.{ext}", "w").write(open(f"{src}.{ext}").read())
    with gzip.open(src.rsplit("/", 1)[0] + "/traits_EC.txt.gz", "rt",
                   newline="") as f, \
            gzip.open(f"{bc}/traits_EC.txt.gz", "wt", newline="") as o:
        rd, wr = csv.reader(f, delimiter="\t"), csv.writer(
            o, delimiter="\t", lineterminator="\n")
        wr.writerow(next(rd))
        rows = list(rd)
        rows[0][0] = "g_999999.9"
        for r in rows:
            wr.writerow(r)

    # ---- recipe control: stock ALIGNED seqs, 160 retained / 20 held out --- #
    print("  building the recipe-control reference (arm B recipe)")
    aln, name = {}, None
    for line in open(f"{STOCK}/bacteria/bac_ref/bac_ref.fna"):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            aln[name] = []
        else:
            aln[name].append(line)
    aln = {k: "".join(v) for k, v in aln.items()}
    with gzip.open(f"{STOCK}/bacteria/ec.txt.gz", "rt") as fh:
        rd = csv.reader(fh, delimiter="\t")
        ec_head = next(rd)
        stock_ec = {r[0]: r for r in rd}
    ids = sorted(set(aln) & set(stock_ec))
    pick = random.Random(SEED).sample(ids, 180)
    retain, holdout = pick[:160], pick[160:]
    swidth = len(aln[retain[0]])
    rdir = f"{W}/refs/recipe_B"
    os.makedirs(rdir, exist_ok=True)
    with open(f"{rdir}/recipe_B.fna", "w") as fh:
        for i in retain:
            fh.write(f">{i}\n{aln[i]}\n")
    p2(f"hmmbuild --dna --informat afa {rdir}/recipe_B.hmm {rdir}/recipe_B.fna "
       f"> /dev/null")
    p2(f"FastTree -nt -gtr -gamma {rdir}/recipe_B.fna > {rdir}/recipe_B.tre "
       f"2> {rdir}/recipe_B.ft.stderr")
    open(f"{rdir}/recipe_B.model", "w").write(
        stock_bac_model.split(", noname")[0] + f", noname = 1-{swidth}\n")
    json.dump({"name": "recipe_B", "mode": "recipe_control", "domain": "bac",
               "n_tips": len(retain), "msa_width": swidth},
              open(f"{rdir}/manifest.json", "w"), indent=1)
    with open(f"{W}/refs/recipe_control_queries.fna", "w") as fh:
        for i in holdout:
            s = aln[i].replace("-", "").replace(".", "").upper()
            fh.write(f">{i}\n{s[len(s) // 3: len(s) // 3 + 290]}\n")
    json.dump({i: [ec_head[j + 1] for j, v in enumerate(stock_ec[i][1:])
                   if v not in ("0", "0.0")] for i in holdout},
              open(f"{W}/refs/recipe_control_truth.json", "w"))
    print(f"    recipe_B: {len(retain)} tips, {len(holdout)} held-out queries, "
          f"stock trait table {len(stock_ec)} rows x {len(ec_head) - 1} cols")
    return poison


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", default="/tmp/sr_fix_selftest")
    ap.add_argument("--processes", type=int, default=8)
    ap.add_argument("--keep-fixture", action="store_true",
                    help="reuse an existing fixture instead of rebuilding it")
    a = ap.parse_args()
    W = os.path.abspath(a.workdir)
    os.makedirs(W, exist_ok=True)
    sys.path.insert(0, os.path.join(RV, "picrust2_comparison", "scripts"))

    print("== 0. fixture ==")
    if a.keep_fixture and os.path.exists(f"{W}/fold_plan.json"):
        poison = json.load(open(f"{W}/fixture.json"))["poison"]
        print("  reusing the existing fixture")
    else:
        poison = build_fixture(W)
        json.dump({"poison": poison}, open(f"{W}/fixture.json", "w"))

    base = ["--fold-plan", f"{W}/fold_plan.json", "--panel",
            f"{W}/panel_eligible.json", "--queries", f"{W}/queries_poisoned.fna",
            "--processes", str(a.processes)]
    O = f"{W}/picrust2"
    subprocess.run(["rm", "-rf", O], check=True)

    print("== 1. --mode fold (happy path, 1 poisoned query) ==")
    run(["--mode", "fold", "--fold", "0", "--domain", "bac", "--density", "1.0",
         "--ref-dir", f"{W}/refs/f00_bac_d100", "--outdir", f"{O}/f00_bac_d100"]
        + base, 0)
    d = json.load(open(f"{O}/f00_bac_d100/predictions.json"))
    m = d["meta"]
    check("n_min_align_dropped", m["n_min_align_dropped"], 1)
    check("poisoned md5 is the dropped one", m["min_align_dropped_md5"], [poison])
    check("n_placed == n_query - dropped", m["n_placed"],
          m["n_query"] - m["n_min_align_dropped"])
    check("C3 unknown_rows_match_queries", m["unknown_rows_match_queries"], True)
    check("C3 status", m["controls"]["C3"]["bac"]["status"], "pass")
    check("C4 status", m["controls"]["C4"]["f00_bac_d100"]["status"], "pass")
    check("every prediction has a nearest-tip label",
          sum(1 for v in d["pred"].values()
              if not v["closest_reference_genome"].startswith("g_")), 0)
    check("no empty prediction", sum(1 for v in d["pred"].values() if not v["ec"]), 0)
    drops = json.load(open(f"{O}/f00_bac_d100/min_align_dropped.json"))
    check("min_align_dropped.json records len/aligned_len/regions",
          sorted(drops[poison]), ["aligned_len", "len", "regions"])

    print("== 2. resume: nothing is recomputed ==")
    out = run(["--mode", "fold", "--fold", "0", "--domain", "bac", "--ref-dir",
               f"{W}/refs/f00_bac_d100", "--outdir", f"{O}/f00_bac_d100"] + base, 0)
    check("placement reused", "reusing existing placement" in out, True)
    check("hsp reused", "reusing existing" in out, True)

    print("== 3. C3 must FIRE without the g_ prefix (exit 3) ==")
    out = run(["--mode", "fold", "--fold", "0", "--domain", "bac", "--ref-dir",
               f"{W}/refs/nopfx_bac_d100", "--outdir", f"{O}/nopfx"] + base, 3)
    # The C3 PREFLIGHT catches the numeric label set before castor ever runs, so
    # the message is run_picrust2_fold.py's "C3 preflight failed for domain bac"
    # (the post-hoc "C3 failed for domain ..." only fires when the preflight
    # passed and the round trip still lost tips).  Asserting the stale literal
    # made SPEC 11.1 smoke gate #3 exit 1 on a HEALTHY harness, which would abort
    # an orchestrator that gates the launch on this script's status.
    check("C3 fired before castor ran",
          ("C3 preflight failed" in out) or ("C3 failed" in out), True)
    check("C3 named the trailing-zero tips", "trailing zero" in out, True)

    print("== 4. C4 must FIRE on a non-tip trait row (exit 3) ==")
    out = run(["--mode", "fold", "--fold", "0", "--domain", "bac", "--ref-dir",
               f"{W}/refs/badc4_bac_d100", "--outdir", f"{O}/badc4"] + base, 3)
    check("C4 named the stray row", "g_999999.9" in out, True)

    print("== 5. --mode recipe_control (stock 26,868-row trait table) ==")
    run(["--mode", "recipe_control", "--ref-dir", f"{W}/refs/recipe_B",
         "--queries", f"{W}/refs/recipe_control_queries.fna",
         "--truth", f"{W}/refs/recipe_control_truth.json",
         "--traits", f"{STOCK}/bacteria/ec.txt.gz",
         "--outdir", f"{O}/recipe_B", "--processes", str(a.processes)], 0)
    rm = json.load(open(f"{O}/recipe_B/predictions.json"))["meta"]
    check("trait rows subset to the tree's tips",
          rm["trait_table"]["bac"]["n_rows_used"], 160)
    check("the §5.4 column screen ran", rm["trait_table"]["bac"]["screen"]["applied"],
          True)
    check("recipe_control emitted every query", rm["n_placed"], rm["n_query"])

    print("== 6. --routing nsti (both domains, sibling auto-derivation) ==")
    run(["--mode", "fold", "--fold", "0", "--routing", "nsti", "--density", "1.0",
         "--ref-dir", f"{W}/refs/f00_bac_d100", "--outdir", f"{O}/f00_nsti"] + base, 0)
    nd = json.load(open(f"{O}/f00_nsti/predictions.json"))
    check("both domain references ran", sorted(nd["meta"]["ref_name"]),
          ["f00_arc_d100", "f00_bac_d100"])

    def nsti_of(path):
        with gzip.open(path, "rt") as fh:
            rd = csv.reader(fh, delimiter="\t")
            h = next(rd)
            i = h.index("metadata_NSTI")
            return {r[0]: float(r[i]) for r in rd}, set(h[1:]) - {
                "metadata_NSTI", "closest_reference_genome"}

    nb, cb = nsti_of(f"{O}/f00_nsti/bac/EC_predicted.tsv.gz")
    na, ca = nsti_of(f"{O}/f00_nsti/arc/EC_predicted.tsv.gz")
    check("argmin-NSTI routing violations",
          sum(1 for k, v in nd["pred"].items()
              if v["best_domain"] != ("bac" if nb[k] <= na[k] else "arc")
              or abs(v["nsti"] - min(nb[k], na[k])) > 1e-12), 0)
    with gzip.open(f"{O}/f00_nsti/combined_EC_predicted.tsv.gz", "rt") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    union = set(rows[0][1:])
    check(f"union columns == bac | arc ({len(cb)} | {len(ca)})",
          (len(union), union == cb | ca), (len(cb | ca), True))
    absent = [i for i, c in enumerate(rows[0][1:]) if c not in ca]
    arc_won = {k for k, v in nd["pred"].items() if v["best_domain"] == "arc"}
    check("arc-won rows are 0-filled where the arc table lacks a column",
          all(r[1:][i] == "0" for r in rows[1:] if r[0] in arc_won for i in absent),
          True)

    print("== 7. include_self / plumbing ==")
    for mode in ("include_self", "plumbing"):
        run(["--mode", mode, "--fold", "0", "--domain", "bac", "--ref-dir",
             f"{W}/refs/f00_bac_d100", "--outdir", f"{O}/{mode}"] + base, 0)

    print("== 8. usage guards (exit 2) ==")
    run(["--mode", "recipe_control", "--ref-dir", f"{W}/refs/recipe_B",
         "--outdir", f"{W}/x1"], 2)
    run(["--mode", "fold", "--ref-dir", f"{W}/refs/f00_bac_d100",
         "--outdir", f"{W}/x2"], 2)
    run(["--mode", "fold", "--routing", "nsti", "--density", "0.5", "--ref-dir",
         f"{W}/refs/f00_bac_d100", "--fold-plan", f"{W}/fold_plan.json",
         "--outdir", f"{W}/x3"], 2)

    print("== 9. compare.load_picrust2 reads the outputs verbatim ==")
    import compare
    for sub in ("f00_bac_d100", "recipe_B", "f00_nsti"):
        ec_index = {}
        preds, meta = compare.load_picrust2(f"{O}/{sub}", ec_index)
        d = json.load(open(f"{O}/{sub}/predictions.json"))
        inv = {v: k for k, v in ec_index.items()}
        check(f"{sub}: n preds", len(preds), len(d["pred"]))
        check(f"{sub}: EC sets identical to predictions.json",
              all(sorted(inv[i] for i in preds[k]) == d["pred"][k]["ec"]
                  for k in preds), True)
        check(f"{sub}: NSTI + best_domain + label identical",
              all(abs(meta[k][0] - d["pred"][k]["nsti"]) < 1e-12
                  and meta[k][1] == d["pred"][k]["best_domain"]
                  and meta[k][2] == d["pred"][k]["closest_reference_genome"]
                  for k in preds), True)

    print(f"\n== {'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + ', '.join(FAILURES)} ==")
    sys.exit(1 if FAILURES else 0)


if __name__ == "__main__":
    main()
