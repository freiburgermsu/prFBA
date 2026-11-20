#!/usr/bin/env python
"""Re-score R4.3: the numbers of PFUNC_ESTIMATORS_RESULTS.md, each read from its artefact by an explicit path.

    pfe_doc_numbers.py --doc DOC --old-dir OLD_PFE_DIR --new-dir NEW_PFE_DIR --out-doc OUT --report REPORT.json

OLD_/NEW_PFE_DIR: a data/pfunc_estimators directory (scores.json, selected.json, verification/{p2_summary,
p3_stage_o,p4_results}.json, verification/p7_k10.log or the p7 numbers in --p7-old/--p7-new). For every
(line, token index) in SPEC the getter is evaluated on OLD and printed at the token's precision: it must equal the
printed token (CONTROL). The same getter on NEW gives the replacement (style kept: unicode minus, '+', commas).
A control mismatch is reported and that token is left unchanged.
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from doc_remap import NUM, norm, restyle  # noqa: E402

EST5 = ("E0", "E1", "E2", "E3", "E4(E3)")
EST4 = ("E1", "E2", "E3", "E4(E3)")
EST3 = ("E1", "E3", "E4(E3)")


class Src:
    def __init__(self, d, p7):
        j = lambda f: json.load(open(os.path.join(d, f)))
        self.s = j("scores.json")
        self.sel = j("selected.json")
        self.p2 = j("verification/p2_summary.json")
        self.p3 = j("verification/p3_stage_o.json")
        self.p4 = j("verification/p4_results.json")
        self.p7 = p7

    def C(self, st, d="d100"):
        return self.s["composites"][d][st]


def SC(st, est, met, part="est", d="d100"):
    return lambda S: S.C(st, d)[est][met][part]


def CV(st, est, met, part="est", d="d100"):
    return lambda S: S.C(st, d)["contrasts_vs_E0"][est][met][part]


def PW(st, pair, met, part="est", d="d100"):
    return lambda S: S.C(st, d)["contrasts_pairwise"][pair][met][part]


def P4E(st, name, met):
    return lambda S: S.p4["strata"][st]["est"][name][met]


def P4C(st, name, met, i):
    return lambda S: S.p4["strata"][st]["ci"][name][met][i]


def P4K(st, pair, met, i):
    return lambda S: S.p4["strata"][st]["contrasts"][pair][met][i]


def R2(g, first):
    """the doc's double rounding: round to `first` decimals (the README print) before the token's precision."""
    from decimal import Decimal, ROUND_HALF_UP
    def h(S):
        x = Decimal(repr(round(float(g(S)), first)))
        return x
    h.double = True
    return h


def tri(f, *a):
    return [f(*a, p) for p in ("est", "lo", "hi")]


def p4tri(st, pair, met):
    return [P4K(st, pair, met, i) for i in range(3)]


def p4eci(st, name, met):
    return [P4E(st, name, met), P4C(st, name, met, 0), P4C(st, name, met, 1)]


def probe_range(fam, which):
    """min / max over outer folds of (best probe train log-loss - best grid train log-loss)."""
    def g(S):
        t = S.p4["train_ll_by_variant_and_outer_fold"]
        grid = [k for k in t if k.startswith(fam + "|")]
        probe = [k for k in t if k.startswith(fam + "x|") and (fam != "E3" or "|family|" in k)]
        diffs = []
        for f in t[grid[0]]:
            best = min(t[k][f] for k in grid)
            diffs.append(min(t[k][f] for k in probe) - best)
        return (min if which == "min" else max)(diffs)
    return g


def oof_r1(which):
    def g(S):
        m = S.sel["E4_isotonic_maps_out_of_fold"]["d100"]["E4(E3)"]
        r = [max(v["y"]) for v in m.values()]
        return (min if which == "min" else max)(r)
    return g


N = None
SPEC = {
    11: {1: SC("all", "E4(E0)", "union_auroc"), 2: SC("all", "E4(E0)", "d_macro_F1_vs_union"),
         8: SC("all", "E0", "universe_brier"), 9: SC("all", "E1", "universe_brier"),
         10: SC("all", "E0", "union_brier"), 11: SC("all", "E1", "union_brier"),
         12: SC("all", "E4(E0)", "union_brier"), 13: CV("all", "E4(E0)", "union_brier"),
         **dict(zip(range(14, 23), p4tri("all", "E3 - E4(E0)", "universe_brier")
                    + p4tri("all", "E3 - E4(E0)", "union_brier") + p4tri("all", "E3 - E4(E0)", "union_auroc")))},
    15: dict(enumerate(tri(PW, "all", "E3 - E1", "universe_brier"))),
    16: dict(enumerate(tri(PW, "all", "E3 - E1", "union_brier"))),
    17: {0: SC("all", "E1", "union_slope_wls"), 1: SC("all", "E3", "union_slope_wls"),
         **dict(zip(range(5, 11), p4tri("all", "E3 - E1|K=100|tau=0.02", "universe_brier")
                    + p4tri("all", "E3 - E1|K=100|tau=0.02", "union_brier"))),
         11: P4E("all", "E1|K=100|tau=0.02", "union_wls_slope"), 12: P4E("all", "E3", "union_wls_slope"),
         **dict(zip(range(13, 16), p4tri("all", "E3 - E1|K=100|tau=0.02", "macro_dF1_vs_union")))},
    19: {0: SC("all", "E2", "union_brier"), 1: SC("all", "E1", "union_brier")},
    23: dict(enumerate(tri(PW, "all", "E4(E3) - E3", "universe_brier"))),
    24: dict(enumerate(tri(PW, "all", "E4(E3) - E3", "union_brier"))),
    26: {**dict(enumerate(tri(PW, "all", "E4(E3) - E3", "universe_logloss_eps0.01"))),
         3: SC("all", "E4(E3)", "union_slope_wls"), 4: SC("all", "E3", "union_slope_wls"),
         5: SC("all", "E4(E3)", "universe_slope_ols"), 6: SC("all", "E3", "universe_slope_ols"),
         **dict(zip(range(8, 11), tri(PW, "all", "E4(E3) - E3", "macro_F1_q_ge_0.5")))},
    28: {3: SC("n1", "E3", "union_auroc"), 4: SC("n1", "E0", "union_brier"), 5: R2(SC("n1", "E3", "union_brier"), 4),
         **dict(zip((6, 7, 8), p4eci("n1", "PREV_foldref", "union_auroc"))),
         **dict(zip((9, 10, 11), p4tri("n1", "E3 - PREV_foldref", "union_auroc"))),
         12: SC("n1", "E4(E0)", "union_mean_q"), 13: SC("n1", "E4(E0)", "union_brier"),
         **dict(zip((15, 16, 17), p4tri("n1", "E3 - E4(E0)", "union_brier"))),
         18: CV("n1", "E3", "union_brier"),
         **dict(zip((20, 21, 22), tri(CV, "multi", "E3", "union_brier"))),
         23: SC("multi", "E0", "union_slope_wls"), 24: R2(SC("multi", "E3", "union_slope_wls"), 3)},
    32: {**dict(zip((0, 1, 2), tri(SC, "all", "E3", "d_macro_F1_vs_union"))),
         4: P4E("all", "E3", "pred_q_ge_0.5_outside_union_per_cell"), 5: P4E("all", "E3", "tp_outside_union_per_cell"),
         **dict(zip((6, 7, 8), p4eci("all", "E3", "macro_dF1_vs_union_restricted_to_union"))),
         9: P4E("n1", "E3", "macro_dF1_vs_union_restricted_to_union"),
         **dict(zip((11, 12, 13), p4eci("multi", "E3", "macro_dF1_vs_union_restricted_to_union"))),
         15: P4E("all", "E1", "macro_dF1_vs_union_restricted_to_union")},
    33: dict(enumerate(tri(SC, "all", "E4(E3)", "d_macro_F1_vs_union"))),
    34: {0: R2(SC("all", "E1", "d_macro_F1_vs_union"), 4), 1: SC("all", "E1", "d_macro_F1_vs_union", "lo"),
         2: SC("all", "E1", "d_macro_F1_vs_union", "hi")},
    35: dict(enumerate(tri(SC, "all", "E2", "d_macro_F1_vs_union"))),
    37: {1: SC("all", "E0", "d_macro_F1_vs_union")},
    39: {0: CV("all", "E3*", "universe_brier", d="d050"), 1: CV("all", "E3*", "universe_brier", d="d025"),
         2: CV("all", "E3*", "union_brier", d="d050"), 3: CV("all", "E3*", "union_brier", d="d025"),
         **dict(zip((4, 5, 6), [PW("all", "E4(E3*) - E3*", "union_brier", p, d="d025") for p in ("est", "lo", "hi")]))},
    61: dict(enumerate(sum([tri(SC, "all", e, "universe_brier") for e in EST5], []))),
    62: dict(enumerate(sum([tri(CV, "all", e, "universe_brier") for e in EST4], []))),
    63: dict(enumerate(sum([tri(SC, "all", e, "universe_logloss_eps0.01") for e in EST5], []))),
    64: dict(enumerate(sum([tri(CV, "all", e, "universe_logloss_eps0.01") for e in EST4], []))),
    65: dict(enumerate([SC("all", e, "universe_ece_decile") for e in EST5])),
    66: dict(enumerate(sum([tri(SC, "all", e, "union_brier") for e in EST5], []))),
    67: dict(enumerate(sum([tri(CV, "all", e, "union_brier") for e in EST4], []))),
    68: dict(enumerate(sum([tri(SC, "all", e, "union_slope_wls") for e in EST5], []))),
    69: dict(enumerate(sum([tri(SC, "all", e, "union_auroc") for e in EST5], []))),
    70: dict(zip(range(1, 16), sum([tri(SC, "all", e, "d_macro_F1_vs_union") for e in EST5], []))),
    71: {0: SC("all", "E0", "d_micro_F1_vs_union"),
         **dict(zip(range(1, 13), sum([tri(SC, "all", e, "d_micro_F1_vs_union") for e in EST4], [])))},
    72: dict(zip(range(1, 5), [SC("all", e, "support_mean") for e in ("E0", "E1", "E2", "E3")])),
    74: {0: SC("all", "E4(E3)", "universe_ece_decile", "lo"), 1: SC("all", "E4(E3)", "universe_ece_decile", "hi"),
         2: CV("all", "E4(E0)", "universe_brier"), 3: CV("all", "E4(E0)", "union_brier")},
    80: dict(enumerate(sum([tri(CV, "n1", e, "universe_brier") for e in EST3], []))),
    81: dict(enumerate([CV("n1", "E1", "universe_logloss_eps0.01")]
                       + sum([tri(CV, "n1", e, "universe_logloss_eps0.01") for e in ("E3", "E4(E3)")], []))),
    82: dict(enumerate([SC("n1", "E0", "union_brier")] + sum([tri(CV, "n1", e, "union_brier") for e in EST3], []))),
    83: dict(enumerate([SC("n1", "E0", "union_auroc"), SC("n1", "E1", "union_auroc")]
                       + sum([tri(SC, "n1", e, "union_auroc") for e in ("E3", "E4(E3)")], []))),
    84: dict(enumerate([SC("n1", "E1", "union_slope_wls")]
                       + sum([tri(SC, "n1", e, "union_slope_wls") for e in ("E3", "E4(E3)")], []))),
    85: dict(enumerate(sum([tri(SC, "n1", e, "d_macro_F1_vs_union") for e in EST3], []))),
    91: dict(enumerate([SC("multi", "E0", "universe_brier")] + sum([tri(CV, "multi", e, "universe_brier") for e in EST3], []))),
    92: dict(enumerate([CV("multi", "E1", "universe_logloss_eps0.01")]
                       + sum([tri(CV, "multi", e, "universe_logloss_eps0.01") for e in ("E3", "E4(E3)")], []))),
    93: dict(enumerate([SC("multi", "E0", "union_brier")] + sum([tri(CV, "multi", e, "union_brier") for e in EST3], []))),
    94: dict(enumerate([SC("multi", "E0", "union_slope_wls"), SC("multi", "E1", "union_slope_wls")]
                       + sum([tri(SC, "multi", e, "union_slope_wls") for e in ("E3", "E4(E3)")], []))),
    95: dict(enumerate([SC("multi", "E0", "union_auroc"), SC("multi", "E1", "union_auroc")]
                       + sum([tri(SC, "multi", e, "union_auroc") for e in ("E3", "E4(E3)")], []))),
    96: dict(enumerate([SC("multi", "E0", "d_macro_F1_vs_union")]
                       + sum([tri(SC, "multi", e, "d_macro_F1_vs_union") for e in EST3], []))),
    98: {0: CV("n2", "E3", "universe_brier"), 2: CV("n3", "E3", "universe_brier"), 4: CV("n4+", "E3", "universe_brier"),
         6: CV("n2", "E3", "union_brier"), 7: CV("n3", "E3", "union_brier"), 8: CV("n4+", "E3", "union_brier"),
         11: PW("n1", "E3 - E1", "union_brier"), **dict(zip((13, 14, 15), tri(PW, "multi", "E3 - E1", "union_brier")))},
    104: {0: SC("all", "E0", "universe_brier", d="d050"), 1: SC("all", "E0", "union_brier", d="d050"),
          2: SC("all", "E0", "universe_brier", d="d025"), 3: SC("all", "E0", "union_brier", d="d025")},
    105: dict(enumerate(tri(CV, "all", "E3*", "universe_brier", d="d050") if False else
                        [CV("all", "E3*", "universe_brier", p, d) for d in ("d050", "d025") for p in ("est", "lo", "hi")])),
    106: dict(enumerate([CV("all", "E3*", "union_brier", p, d) for d in ("d050", "d025") for p in ("est", "lo", "hi")])),
    107: {0: PW("all", "E3* - E1*", "union_brier", d="d050"), 1: PW("all", "E3* - E1*", "union_brier", d="d025")},
    108: dict(enumerate([PW("all", "E4(E3*) - E3*", "union_brier", p, d) for d in ("d050", "d025") for p in ("est", "lo", "hi")])),
    119: {6: probe_range("E1", "min"), 7: probe_range("E1", "max"), 13: probe_range("E3", "min"), 14: probe_range("E3", "max"),
          **dict(zip(range(17, 23), p4tri("all", "E3x|tau=0.05|family|alpha=2.0 - E3", "universe_brier")
                     + p4tri("all", "E3x|tau=0.05|family|alpha=2.0 - E3", "union_brier")))},
    120: {4: lambda S: S.p3["main_isotonic"]["knots"], 6: lambda S: S.p3["main_isotonic"]["r0"],
          8: lambda S: S.p3["main_isotonic"]["r_at_0.5"], 10: lambda S: S.p3["main_isotonic"]["r1"],
          12: oof_r1("min"), 13: oof_r1("max"), 14: lambda S: S.p3["main_isotonic"]["max_y"]},
    121: {0: lambda S: S.p3["clade_rank_counts_main_file"]["family"], 1: lambda S: S.p3["clade_rank_counts_main_file"]["order"],
          2: lambda S: S.p3["clade_rank_counts_main_file"]["none"]},
    122: {0: lambda S: S.p3["median_neighbours"], 2: lambda S: S.p3["median_clade_members"]},
    128: {4: lambda S: S.s["reproduction"]["E1_200hit_K10_tau0.01_genome_q_ge_0.5_macro_F1"]["value"],
          5: lambda S: S.s["reproduction"]["E1_200hit_K10_tau0.01_genome_q_ge_0.5_macro_F1"]["target"],
          8: lambda S: S.p7["n_near"]},
    129: {0: lambda S: S.s["reproduction"]["E0_universe_brier"]["value"],
          1: lambda S: S.s["reproduction"]["E0_union_macro_F1"]["value"]},
    131: {8: lambda S: S.s["reproduction"]["E1_200hit_K20_tau0.01_genome_universe_brier"]["value"]},
    134: {2: lambda S: S.s["reproduction"]["E1_100hit_K20_tau0.01_genome_universe_brier"]["value"]},
}


def fmt_like(x, tok):
    dec = len(tok.split(".")[1]) if "." in tok else 0
    from decimal import Decimal, ROUND_HALF_UP
    if isinstance(x, Decimal):
        s = str(x.quantize(Decimal(1).scaleb(-dec), rounding=ROUND_HALF_UP))
        return s.lstrip("-") if float(s) == 0 else s
    s = f"{float(x):.{dec}f}"
    return s.lstrip("-") if float(s) == 0 else s


def p7_numbers(path):
    t = open(path).read()
    m = re.search(r"dot-product:\s*([0-9.]+)\s+cumulative-sum:\s*([0-9.]+).*?:\s*(\d+)", t)
    return {"dot": float(m.group(1)), "cum": float(m.group(2)), "n_near": int(m.group(3))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--doc", required=True)
    ap.add_argument("--old-dir", required=True)
    ap.add_argument("--new-dir", required=True)
    ap.add_argument("--p7-old", default=None, help="n_near etc. for the old run (int), if no old p7 log")
    ap.add_argument("--out-doc", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()
    p7o = {"n_near": int(a.p7_old)} if a.p7_old else p7_numbers(os.path.join(a.old_dir, "verification/p7_k10.log"))
    p7n = p7_numbers(os.path.join(a.new_dir, "verification/p7_k10.log"))
    So, Sn = Src(a.old_dir, p7o), Src(a.new_dir, p7n)
    doc = open(a.doc).read().splitlines(keepends=True)
    rep = {"control_fail": [], "changed": [], "unchanged": 0}
    for ln, spec in SPEC.items():
        line = doc[ln - 1]
        ms = list(NUM.finditer(line))
        out, last = [], 0
        for i, m in enumerate(ms):
            g = spec.get(i)
            if g is None:
                continue
            sign, intp, frac = m.groups()
            tok = norm(sign, intp, frac)
            try:
                old = fmt_like(g(So), tok)
            except Exception as e:  # noqa: BLE001
                rep["control_fail"].append({"line": ln, "i": i, "tok": m.group(0), "error": repr(e)})
                continue
            if float(old) != float(tok):
                rep["control_fail"].append({"line": ln, "i": i, "tok": m.group(0), "old_value": old})
                continue
            new = fmt_like(g(Sn), tok)
            if float(new) == float(tok) and new.lstrip("-") == tok.lstrip("-"):
                rep["unchanged"] += 1
                continue
            rep["changed"].append({"line": ln, "i": i, "old": m.group(0), "new": restyle(new, sign, intp)})
            out.append(line[last:m.start()] + restyle(new, sign, intp))
            last = m.end()
        out.append(line[last:])
        doc[ln - 1] = "".join(out)
    open(a.out_doc, "w").write("".join(doc))
    json.dump(rep, open(a.report, "w"), indent=1, ensure_ascii=False)
    print(f"control failures {len(rep['control_fail'])}; changed {len(rep['changed'])}; unchanged {rep['unchanged']}")
    for c in rep["control_fail"]:
        print("CONTROL", c)


if __name__ == "__main__":
    main()
