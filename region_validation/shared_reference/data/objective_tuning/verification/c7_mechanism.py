"""Where does the truth-absent flux go?  Per confirm row, at the returned vertex, the change
(arm a minus arm b) in each flux class, in units of D; pooled A (and B) with the cluster
bootstrap.  Classes: gene-backed truth-absent (= flux_absent), gene-backed truth-supported
(= flux_gene_backed - flux_absent), gapfilled (flux_gapfilled), gene-less draft (= flux_rxn
- flux_gene_backed - flux_gapfilled; equals internal_term_split.geneless_draft to 3e-13),
exchange (flux_exchange).  Also the face width, and row-level sign patterns."""
import json
import os

import numpy as np

import vlib as L

rows = L.load_rows()
WA, WB, _ = L.weights()
POOLS = {"A": WA, "B": WB}
conf = [r for r in rows if r["half"] == "confirm"]
tune = [r for r in rows if r["half"] == "tune"]


def cls(r, arm):
    a = r["arms"][arm]
    D = r["model"]["D_flux_gb_two_tier"]
    gb, ab, gf, rx, ex = (a["flux_gene_backed"], a["flux_absent"], a["flux_gapfilled"],
                          a["flux_rxn"], a["flux_exchange"])
    return {"absent": ab / D, "supported_gb": (gb - ab) / D, "gene_backed": gb / D,
            "gapfilled": gf / D, "geneless_draft": (rx - gb - gf) / D,
            "all_internal": rx / D, "exchange": ex / D,
            "face_width": a["taf_face_width"], "face_mid": a["taf_face_mid"],
            "face_min": a["taf_face_min"], "face_max": a["taf_face_max"],
            "vertex_minus_facemin": (ab / D) - a["taf_face_min"],
            "n_active_gb": a["n_active_gb"], "n_active_absent": a["n_active_absent"]}


KEYS = ("absent", "supported_gb", "gene_backed", "gapfilled", "geneless_draft",
        "all_internal", "exchange", "face_width", "vertex_minus_facemin", "n_active_gb",
        "n_active_absent")
COMPS = {"selected-shipped": (L.SELECTED, L.SHIPPED),
         "est@X100 - shipped": ("est/lin/G0/X100", L.SHIPPED),
         "kn X1 - kn X100": ("kn/lin/G0/X1", L.SHIPPED),
         "est X1 - est X100": (L.SELECTED, "est/lin/G0/X100"),
         "oracle@X1 - two_tier@X1": ("oracle/lin/G0/X1", "two_tier/lin/G0/X1"),
         "selected - two_tier@X1": (L.SELECTED, "two_tier/lin/G0/X1"),
         "selected - est_shuffled@X1": (L.SELECTED, "est_shuffled/lin/G0/X1")}
out = {}
for half, rs in (("confirm", conf), ("tune", tune)):
    out[half] = {}
    for name, (a, b) in COMPS.items():
        rec = {}
        for k in KEYS:
            c = L.paired(rs, lambda r: cls(r, a)[k], lambda r: cls(r, b)[k], pools=POOLS)
            rec[k] = {"A": c["pooled"]["A"], "B": c["pooled"]["B"],
                      "per": {s: {"mean": p["mean"], "ci": p["ci"]} for s, p in c["per"].items()}}
        out[half][name] = rec

# row-level: on n>=2 confirm rows, how often does absent fall while gene-less rises, etc.
pat = {}
for name, (a, b) in COMPS.items():
    n = 0
    cnt = {"absent_down": 0, "absent_down_supported_up": 0, "absent_down_supported_down": 0,
           "absent_down_geneless_up": 0, "absent_down_gapf_up": 0, "absent_down_exchange_up": 0,
           "absent_down_gb_down_by_more_than_absent": 0}
    big = []
    for r in conf:
        if r["stratum"] not in L.MULTI:
            continue
        n += 1
        x, y = cls(r, a), cls(r, b)
        d = {k: x[k] - y[k] for k in ("absent", "supported_gb", "geneless_draft", "gapfilled",
                                      "exchange", "gene_backed")}
        if d["absent"] < -1e-9:
            cnt["absent_down"] += 1
            cnt["absent_down_supported_up"] += d["supported_gb"] > 1e-9
            cnt["absent_down_supported_down"] += d["supported_gb"] < -1e-9
            cnt["absent_down_geneless_up"] += d["geneless_draft"] > 1e-9
            cnt["absent_down_gapf_up"] += d["gapfilled"] > 1e-9
            cnt["absent_down_exchange_up"] += d["exchange"] > 1e-9
            cnt["absent_down_gb_down_by_more_than_absent"] += d["gene_backed"] < d["absent"] - 1e-9
        big.append((d["gene_backed"], r["key"], r["group_id"], r["stratum"], d))
    big.sort(key=lambda t: t[0])
    pat[name] = {"n_rows_n>=2": n, **{k: int(v) for k, v in cnt.items()},
                 "largest_gene_backed_drops": [(round(t[0], 4), t[1], t[2], t[3],
                                                {k: round(v, 4) for k, v in t[4].items()})
                                               for t in big[:6]]}
out["row_patterns_confirm_n>=2"] = pat


# OLS-free accounting: share of the pooled-A drop in absent flux that is matched by a rise in
# supported gene-backed flux vs by a rise in gene-less / gapfilled / exchange
def A(half, name, k):
    return out[half][name][k]["A"]["mean"]


acc = {}
for half in ("confirm", "tune"):
    acc[half] = {}
    for name in COMPS:
        acc[half][name] = {k: A(half, name, k) for k in ("absent", "supported_gb", "gene_backed",
                                                          "gapfilled", "geneless_draft",
                                                          "all_internal", "exchange",
                                                          "face_width")}
out["pooled_A_accounting"] = acc
json.dump(out, open(os.path.join(L.HERE, "out_c7.json"), "w"), indent=1, default=str)
for half in ("confirm", "tune"):
    print("=== %s, pooled A (B), units of D, returned vertex" % half)
    for name in COMPS:
        r = out[half][name]
        print("  %-28s" % name + " ".join("%s %s" % (k, L.fmt(r[k]["A"])) for k in
                                         ("absent", "supported_gb", "gapfilled", "geneless_draft",
                                          "exchange", "face_width")))
        print("  %-28s" % "  B:" + " ".join("%s %s" % (k, L.fmt(r[k]["B"])) for k in
                                           ("absent", "supported_gb", "gapfilled", "geneless_draft",
                                            "exchange")))
print(json.dumps(pat, indent=1, default=str)[:6000])
