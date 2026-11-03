"""Phase 3: independent recomputation of the Stage O weight files for the 360 amplicons of
data/weighted_objective/population.json (neighbours, weights, clade rank + members, alpha,
reps), and of the consumer q, for stage_o_weights.json, stage_o_weights_E1.json and
stage_o_weights_E3.json."""
import collections
import json
import math
import os
import random

import numpy as np

import p2_pass as P
import vctx

P.setup()
G = P.G
C = G["C"]
cells = G["cells"]
cell_of = {a: i for i, a in enumerate(cells.amplicon)}
pop = vctx.jl(os.path.join(vctx.DATA, "weighted_objective", "population.json"))["jobs"]
D = os.path.join(vctx.DATA, "pfunc_estimators")
files = {k: vctx.jl(os.path.join(D, k)) for k in
         ("stage_o_weights.json", "stage_o_weights_E1.json", "stage_o_weights_E3.json")}
out = collections.OrderedDict()
A = collections.Counter()
for job in pop:
    i = cell_of.get(job["key"])
    if i is None:
        A["pop_key_not_a_scored_cell"] += 1
        continue
    A["pop_md5_mismatch"] += int(cells.md5.iat[i] != job["md5"])
    A["pop_fold_mismatch"] += int(int(cells.fold.iat[i]) != int(job["fold"]))
    A["pop_src_mismatch"] += int(cells.src_genome_id.iat[i] != job["src_genome_id"])
    A["pop_group_mismatch"] += int(cells.group_id.iat[i] != job["group_id"])
    A["pop_reps_mismatch"] += int([C.gids[x] for x in G["reps"][i]] != list(job["reps"]))
out["population_join"] = dict(A)


def expected(i, K, tau, start, alpha):
    g, idv = G["l100"][i]
    g, idv = g[:K], idv[:K]
    w = np.exp(-(idv[0] - idv) / tau)
    nb = [(C.gids[x], float(wx)) for x, wx in zip(g, w)]
    clade, a = None, 0.0
    if start is not None:
        j = G["cl_of"][start][i]
        if j >= 0:
            rk = P.CHAIN[G["cl_rank"][start][i]]
            clade = {"rank": rk, "members": sorted(C.gids[x] for x in G["members"][j])}
            a = float(alpha)
    return nb, clade, a


def compare(fname, K, tau, start, alpha, keys):
    obj = files[fname]
    B = collections.Counter()
    maxw = 0.0
    for key in keys:
        i = cell_of[key]
        a = obj["amplicons"].get(key)
        if a is None:
            B["missing_in_file"] += 1
            continue
        nb, clade, al = expected(i, K, tau, start, alpha)
        fn = [x[0] for x in a["neighbours"]]
        fw = [x[1] for x in a["neighbours"]]
        if fn != [x[0] for x in nb]:
            B["neighbour_ids_or_order_differ"] += 1
        else:
            d = max(abs(u - v[1]) for u, v in zip(fw, nb))
            maxw = max(maxw, d)
            B["weights_differ_gt_1e-9"] += int(d > 1e-9)
        fc = a["clade"]
        if (fc is None) != (clade is None):
            B["clade_presence_differs"] += 1
        elif fc is not None:
            B["clade_rank_differs"] += int(fc["rank"] != clade["rank"])
            B["clade_members_differ"] += int(sorted(fc["members"]) != clade["members"])
        B["alpha_differs"] += int(abs(a["alpha"] - al) > 0)
        B["reps_differ"] += int(a["reps"] != [C.gids[x] for x in G["reps"][i]])
        B["md5_differs"] += int(a["md5"] != cells.md5.iat[i])
        B["fold_differs"] += int(int(a["fold"]) != int(cells.fold.iat[i]))
        B["checked"] += 1
    B["max_abs_weight_diff"] = maxw
    return dict(B)


keys = [j["key"] for j in pop]
random.seed(20261103)
sample40 = random.sample(keys, 40)
out["E4_main_file_all360"] = compare("stage_o_weights.json", 100, 0.02, "family", 2.0, keys)
out["E1_file_all360"] = compare("stage_o_weights_E1.json", 50, 0.05, None, None, keys)
out["E3_file_all360"] = compare("stage_o_weights_E3.json", 100, 0.02, "family", 2.0, keys)
out["random40_keys"] = sample40
out["E4_main_file_random40"] = compare("stage_o_weights.json", 100, 0.02, "family", 2.0, sample40)

# E3 file == main file minus the isotonic map?
m, e3 = files["stage_o_weights.json"], files["stage_o_weights_E3.json"]
out["E3_file_amplicons_identical_to_main"] = (m["amplicons"] == e3["amplicons"])
out["E3_file_meta_isotonic"] = e3["meta"].get("isotonic")
out["meta_main"] = {k: v for k, v in m["meta"].items() if k != "isotonic"}
out["meta_E3"] = {k: v for k, v in e3["meta"].items() if k != "isotonic"}
iso = m["meta"]["isotonic"]
out["main_isotonic"] = {"knots": len(iso["x"]), "r0": iso["y"][0], "r_at_0.5":
                        float(np.interp(0.5, iso["x"], iso["y"])), "r1": iso["y"][-1],
                        "max_y": max(iso["y"])}

# file-derived counts in the draft
cl = collections.Counter((a["clade"]["rank"] if a["clade"] else "none")
                         for a in m["amplicons"].values())
out["clade_rank_counts_main_file"] = dict(cl)
out["median_neighbours"] = float(np.median([len(a["neighbours"]) for a in m["amplicons"].values()]))
out["median_clade_members"] = float(np.median([len(a["clade"]["members"])
                                               for a in m["amplicons"].values() if a["clade"]]))
out["median_clade_members_counting_none_as_0"] = float(np.median(
    [len(a["clade"]["members"]) if a["clade"] else 0 for a in m["amplicons"].values()]))
out["alpha_zero_iff_clade_null"] = all((a["alpha"] == 0.0) == (a["clade"] is None)
                                       for a in m["amplicons"].values())

# consumer q from the E3 file vs my E3 q (K=100, tau=0.02, family, alpha=2), all 360
V = [dict(name="E3", fam="E3", tau=0.02, K=100, start="family", alpha=2.0, grid=True)]
mx = 0.0
for key in keys:
    i = cell_of[key]
    a = e3["amplicons"][key]
    num = np.zeros(C.NU)
    den = 0.0
    for g, w in a["neighbours"]:
        num += w * C.E[C.gi[g]]
        den += w
    if a["clade"] is not None:
        mem = [C.gi[g] for g in a["clade"]["members"]]
        num += a["alpha"] * C.E[mem].mean(0)
        den += a["alpha"]
    q = num / den
    Q, y, u = P.cell_Q(i, V)
    mx = max(mx, float(np.abs(q - Q[0]).max()))
out["consumer_q_from_E3_file_vs_my_E3_max_abs_diff"] = mx

# leakage on the 360 file entries (neighbours + members)
L = collections.Counter()
for key in keys:
    i = cell_of[key]
    s = G["src"][i]
    f = int(cells.fold.iat[i])
    a = m["amplicons"][key]
    gl = [C.gi[g] for g, _ in a["neighbours"]] + ([C.gi[g] for g in a["clade"]["members"]]
                                                 if a["clade"] else [])
    gl = np.array(gl, np.int64)
    L["source_in_file_lists"] += int((gl == s).any())
    for c_ in ("genus", "species", "taxon", "group"):
        L[f"shares_{c_}"] += int((C.code[c_][gl] == C.code[c_][s]).any())
    L["outside_fold_reference"] += int((~G["inref"][f][gl]).sum())
out["leakage_on_file_entries"] = dict(L)
print(json.dumps(out, indent=1, default=float))
with open(os.path.join(vctx.OUT, "p3_stage_o.json"), "w") as fh:
    json.dump(out, fh, indent=1, default=float)
