#!/usr/bin/env python
"""Pull one real benchmark amplicon through every stage, for the condensed deck's workflow
figure: its 100-hit record, its selection with and without the source organism, and the EC
sets of the selected genomes and of the source (the truth).

The amplicon is chosen by a fixed rule, not by hand: a V4 amplicon from the original panel
whose source is its own single selection under self-inclusion and which, with the source
removed at tau = 0.987, falls to a genus-tier tie of four genomes after 2-10 self hits are
stripped; the first such amplicon in sorted order. Depth-100 records, exact-tie union.

    python presentations/condensed/extract_worked_example.py   # -> worked_example.json
"""
import json
import os

import ijson

R = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
D = R + "/region_validation_d100"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "worked_example.json")

sel = json.load(open(D + "/picrust2_comparison/data/prfba_selection.json"))
pick = None
for a in sorted(k for k in sel if k != "_meta"):
    v = sel[a]
    inc, exc = v["include"], v["exclude_0.987"]
    src = a.split("__")[0]
    if (a.endswith("__V4") and inc["n_selected"] == 1 and inc["reps"] == [src]
            and exc["tier"] == "genus" and exc["n_selected"] == 4
            and 2 <= exc["n_removed_self"] <= 10):
        pick = a
        break
assert pick, "no amplicon meets the rule"
src = pick.split("__")[0]
fan = json.load(open(D + "/data/amplicon_md5_fanout.json"))
md5 = next(m for m, ams in fan.items() if pick in ams)
print("[worked] amplicon", pick, "md5", md5)

rec = None
with open(D + "/data/hits/raw/asv_top20_alignment_hits.json", "rb") as fh:
    for k, v in ijson.kvitems(fh, "", use_float=True):
        if k == md5:
            rec = v
            break
assert rec, "record not found"
hits = [{f: h.get(f) for f in ("rank", "identity", "align_score", "genome_id", "taxon_id",
                               "organism")} | {"genus": (h.get("lineage") or {}).get("Genus"),
                                               "phylum": (h.get("lineage") or {}).get("Phylum"),
                                               "species": (h.get("lineage") or {}).get("Species")}
        for h in rec["top20"]]
want = set(sel[pick]["exclude_0.987"]["reps"]) | {src}
ec = {}
with open(D + "/picrust2_comparison/data/genome_function_sets.json", "rb") as fh:
    for k, v in ijson.kvitems(fh, ""):
        if k in want:
            ec[k] = sorted(v["ec"])
            if len(ec) == len(want):
                break
row = None
with open(D + "/data/taxacc_per_record.jsonl") as fh:
    for line in fh:
        if line.startswith('{"amplicon": "%s"' % pick):
            row = json.loads(line)
            break
json.dump({"amplicon": pick, "md5": md5, "source": src, "n_hits": len(hits),
           "n_edlib_candidates": rec.get("n_edlib_candidates"), "asv_len": rec.get("asv_len"),
           "selection": sel[pick], "hits": hits, "ec": ec,
           "truth_lineage": row["truth_lineage"] if row else None,
           "exclude_anchor_scores": (row["self_modes"]["exclude_0.987"]["scores"]["anchor"]
                                     if row else None),
           "rule": __doc__.split("\n\n")[1]}, open(OUT, "w"), indent=1)
print("[worked] wrote", OUT, "| hits", len(hits), "| EC sets", {k: len(v) for k, v in ec.items()})
