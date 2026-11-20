#!/usr/bin/env python
"""Detail of the substantive C-L4 differences (records whose selection or decisions changed):
for each, the decisions entries whose (disposition, stage) moved, with genome ids; and whether
those genomes are among the 4 SNAP20J extra resets or the 20 poisoned keys."""
import json, sys, collections
import ijson
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
old_p, new_p, out_p = sys.argv[1], sys.argv[2], sys.argv[3]
rep = json.load(open(f"{W}/control/d20_code_delta.json"))
keys = json.load(open(f"{W}/control/d20_code_delta.keys.json"))["changed"]
j24 = set(json.load(open(f"{W}/expected/d20_june_reset_24.json")) if isinstance(json.load(open(f"{W}/expected/d20_june_reset_24.json")), list) else json.load(open(f"{W}/expected/d20_june_reset_24.json")).get("keys", []))
p20 = set(json.load(open(f"{W}/expected/repaired_ids.json")).get("poisoned_d20", []))
def load(p, want):
    out = {}
    with open(p, "rb") as fh:
        for k, v in ijson.kvitems(fh, "", use_float=True):
            if k in want: out[k] = v
    return out
want = set(keys)
o, n = load(old_p, want), load(new_p, want)
subst = []
moves = collections.Counter(); genomes = collections.Counter()
for k in keys:
    a, b = o[k], n[k]
    da = {(d.get("genome_id"), d.get("rank")): (d.get("disposition"), d.get("stage")) for d in a.get("decisions", [])}
    db = {(d.get("genome_id"), d.get("rank")): (d.get("disposition"), d.get("stage")) for d in b.get("decisions", [])}
    sa = [s["genome_id"] for s in a.get("selected", [])]; sb = [s["genome_id"] for s in b.get("selected", [])]
    mv = [(g, da.get(g), db.get(g)) for g in set(da) | set(db) if da.get(g) != db.get(g)]
    if not mv and sa == sb and a.get("flags") == b.get("flags"):
        continue
    for g, x, y in mv:
        moves[f"{x} -> {y}"] += 1; genomes[g[0]] += 1
    subst.append({"key": k, "selected_old": sa, "selected_new": sb, "flags_old": a.get("flags"), "flags_new": b.get("flags"),
                  "moves": [[g[0], g[1], x, y] for g, x, y in mv]})
res = {"n_substantive": len(subst), "moves": moves.most_common(), "genomes_moved": genomes.most_common(),
       "genomes_in_snap20j_extra_4": sorted(set(genomes) & j24), "genomes_in_poisoned_20": sorted(set(genomes) & p20),
       "j24": sorted(j24), "records": subst}
json.dump(res, open(out_p, "w"), indent=1)
print(json.dumps({k: res[k] for k in ("n_substantive", "moves", "genomes_moved", "genomes_in_snap20j_extra_4", "genomes_in_poisoned_20")}, indent=1))
