#!/usr/bin/env python
"""Re-score R4.3: update OBJECTIVE_TUNING_RESULTS.md from re-run Stage O outputs.

    otune_doc_numbers.py --doc DOC --old-gen OLD_FILLED --new-gen NEW_FILLED --old OLD_DIR --new NEW_DIR
                         --out-doc OUT --report REPORT.json

OLD_/NEW_DIR: a data/objective_tuning directory (summary.json, verification/out_c*.json).
OLD_/NEW_GEN: objective_tuning_results.py's generated filled.md for the old rows (control run) and the new rows.
1. A doc line equal to a line of OLD_GEN is replaced by NEW_GEN's line at the same index (the two generated files
   must share their non-digit skeleton line by line; checked).
2. The hand-written lines (the answer section and the 2026-11-03 corrections) are mapped token by token through
   SPEC: each getter, evaluated on OLD, must print as the doc's token (CONTROL); evaluated on NEW it gives the
   replacement. Text tokens (row keys) in SPEC_TEXT are checked the same way.
Every doc line with digits that neither step covers is listed in the report (uncovered), so nothing is silently left.
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from doc_remap import NUM, norm, restyle  # noqa: E402

SKEL = re.compile(r"[0-9]")


class Src:
    def __init__(self, d, rise=None):
        self.rise = json.load(open(rise)) if rise else None
        self.s = json.load(open(os.path.join(d, "summary.json")))
        self.c = {}
        for f in os.listdir(os.path.join(d, "verification")):
            if f.startswith("out_") and f.endswith(".json"):
                self.c[f[4:-5]] = json.load(open(os.path.join(d, "verification", f)))


def g(*path, src="s"):
    def h(S):
        o = S.s if src == "s" else S.c[src]
        for k in path:
            o = o[k]
        return o
    return h


def CON(name, met, pool="A", part="mean", i=None):
    p = ("confirm", "contrasts", name, "metrics", met, "pooled", pool, part) + ((i,) if i is not None else ())
    return g(*p)


def STR(name, met, st, part="mean", i=None):
    p = ("confirm", "contrasts", name, "metrics", met, "per_stratum", st, part) + ((i,) if i is not None else ())
    return g(*p)


def C45(name, met, part="mean", i=None, pool="A"):
    p = ("contrasts", name, met, "pooled", pool, part) + ((i,) if i is not None else ())
    return g(*p, src="c4_c5")


def tri(f, *a, **k):
    return [f(*a, part="mean", **k), f(*a, part="ci", i=0, **k), f(*a, part="ci", i=1, **k)]


def x100(h):
    return lambda S: 100 * h(S)


def frac_pct(h):
    return lambda S: 100 * h(S)


EST100 = "est_alone@X100 (est/lin/G0/X100 - kn/lin/G0/X100)"
X1EST = "X1_given_est (est/lin/G0/X1 - est/lin/G0/X100)"
LOGLIN = "log_vs_lin@est/X1 (est/log/G0/X1 - est/lin/G0/X1)"
GMEAS = "Gmeas_vs_G0@est/lin/X1"
SS, ST, SH = "selected-shipped", "selected-two_tier", "selected-shuffled"
DEG = "faces.opt.n_var_1e-3"
ACC = ("pooled_A_accounting", "confirm")
RR = ("selected_row_ratios",)


def acc(name, k):
    return g(*ACC, name, k, src="c7")


def rr(half, met, k, i=None):
    return g(*RR, half, met, k, *(() if i is None else (i,)), src="c2_c3")


def grp_share(which):
    def h(S):
        v = [x["share_rows_grp00000"] for k, x in S.c["c8"]["grp00000_bootstrap_presence"].items() if k != "n1"]
        return 100 * (min(v) if which == "min" else max(v))
    return h


def logo(which):
    return lambda S: S.c["c8"]["logo_confirm"][SS][which][0]


def excl_n(nkey):
    return lambda S: S.c["c8"]["exclusions"][nkey]


SPEC = {
    13: {0: g("meta", "n_by_half", "tune")},
    14: {2: g("meta", "n_by_half", "confirm")},
    18: dict(enumerate(tri(CON, SS, "taf_face_mid"))),
    19: {1: CON(SS, "taf_face_min"), 2: CON(SS, "taf_face_max")},
    20: {0: g("tune", "contrast_vs_shipped_tune", "taf_face_mid", "A", "mean")},
    22: {**dict(enumerate(tri(CON, SH, "taf_face_mid"))), 3: CON(ST, "taf_face_mid")},
    23: {0: g("confirm", "headroom_confirm", "A", "selected_own_cell"), 1: g("confirm", "headroom_confirm", "A", "shipped")},
    26: {**dict(enumerate(tri(C45, EST100, "taf_face_mid"))), 4: C45(X1EST, "taf_face_mid")},
    27: {0: C45(LOGLIN, "taf_face_mid")},
    28: {0: g("meta", "g_meas"), 3: C45(GMEAS, "taf_face_mid")},
    30: {0: g("confirm", "degeneracy", SS, DEG, "pooled", "A", "arm_mean")},
    31: {0: g("confirm", "degeneracy", SS, DEG, "pooled", "A", "comparator_mean"),
         1: g("confirm", "degeneracy", SH, DEG, "pooled", "A", "comparator_mean")},
    33: {0: x100(g("meta", "weights_B", "n1"))},
    34: {0: STR(SS, "taf_face_mid", "n1"), 1: STR(SS, "taf_face_mid", "n1", "ci", 0), 2: STR(SS, "taf_face_mid", "n1", "ci", 1),
         3: g("confirm", "n_by_stratum", "n1"), 4: g("tune_selected_minus_shipped", "per", "n1", "mean", src="c2_c3")},
    38: {0: g("grp00000", "n_rows", src="c8"), 1: g("meta", "n_by_half", "confirm"),
         2: g("grp00000", "contrasts", SS, "taf_face_mid", "without_grp00000_A", "mean", src="c8")},
    39: {0: g("grp00000", "contrasts", SS, "taf_face_mid", "without_grp00000_A", "ci", 0, src="c8"),
         1: g("grp00000", "contrasts", SS, "taf_face_mid", "without_grp00000_A", "ci", 1, src="c8")},
    54: {0: g("meta", "n_by_half", "tune") , },  # placeholder, replaced below
    92: {7: x100(g("eligibility", "max_abs_ratio_dev", "flux_exchange", src="c2_c3")),
         8: x100(g("eligibility", "max_abs_ratio_dev", "flux_gene_backed", src="c2_c3")),
         9: rr("tune", "flux_gene_backed", "quantiles_0_1_5_25_50_75_95_99_100", 0),
         10: rr("tune", "flux_gene_backed", "quantiles_0_1_5_25_50_75_95_99_100", 8),
         11: rr("tune", "flux_gene_backed", "n_outside_10pct"), 12: rr("tune", "flux_gene_backed", "n"),
         14: rr("confirm", "flux_gene_backed", "quantiles_0_1_5_25_50_75_95_99_100", 0),
         15: rr("confirm", "flux_gene_backed", "quantiles_0_1_5_25_50_75_95_99_100", 8),
         16: rr("confirm", "flux_gene_backed", "n_outside_10pct"), 17: rr("confirm", "flux_gene_backed", "n"),
         18: rr("tune", "flux_exchange", "quantiles_0_1_5_25_50_75_95_99_100", 0),
         19: rr("tune", "flux_exchange", "quantiles_0_1_5_25_50_75_95_99_100", 8),
         20: rr("tune", "flux_exchange", "n_outside_10pct"), 21: rr("tune", "flux_exchange", "n"),
         22: rr("confirm", "flux_exchange", "quantiles_0_1_5_25_50_75_95_99_100", 0),
         23: rr("confirm", "flux_exchange", "quantiles_0_1_5_25_50_75_95_99_100", 8),
         24: rr("confirm", "flux_exchange", "n_outside_10pct"), 25: rr("confirm", "flux_exchange", "n")},
    201: {2: x100(g("meta", "weights_B", "n1")), 3: STR(SS, "taf_face_mid", "n1"), 4: STR(SS, "taf_face_mid", "n1", "ci", 0),
          5: STR(SS, "taf_face_mid", "n1", "ci", 1), 6: g("confirm", "n_by_stratum", "n1"),
          7: g("confirm", "n_groups_by_stratum", "n1"),
          8: g("tune_selected_minus_shipped", "per", "n1", "mean", src="c2_c3"),
          9: g("tune_selected_minus_shipped", "per", "n1", "ci", 0, src="c2_c3"),
          10: g("tune_selected_minus_shipped", "per", "n1", "ci", 1, src="c2_c3"),
          **dict(zip((11, 12, 13), tri(CON, SS, "taf_face_mid", pool="B"))),
          15: CON(SS, "taf_face_min"), **dict(zip((16, 17, 18), tri(CON, SS, "taf_point_D")))},
    203: {**dict(enumerate(tri(C45, EST100, "taf_face_mid"))),
          3: acc("est@X100 - shipped", "supported_gb"), 4: lambda S: -acc("est@X100 - shipped", "face_width")(S),
          **dict(zip((5, 6, 7), tri(C45, X1EST, "taf_face_mid"))),
          8: acc("est X1 - est X100", "exchange"), 9: acc("est X1 - est X100", "supported_gb"),
          10: acc("est X1 - est X100", "geneless_draft"), 11: acc("est X1 - est X100", "gapfilled"),
          12: g("precision_decomposition_confirm", "est X1-X100", "precision_flux", "A", "mean", src="c8"),
          13: lambda S: C45(EST100, "taf_face_mid")(S) + C45("X1_alone@kn (kn/lin/G0/X1 - kn/lin/G0/X100)", "taf_face_mid")(S),
          14: CON(SS, "taf_face_mid"),
          15: g(SS, "volume_A", "mean", src="c13"), 16: g(SS, "volume_A", "ci", 0, src="c13"),
          17: g(SS, "volume_A", "ci", 1, src="c13")},
    205: {0: g("confirm", "headroom_confirm", "A", "selected_own_cell"), 1: g("confirm", "headroom_confirm", "A", "shipped"),
          2: g("headroom", "A", "selected_in_shipped_frame", src="c4_c5"),
          3: g("tune", "headroom_tune", "A", "captured_shipped_cell_frame")},
    207: {1: g("confirm", "degeneracy", SS, DEG, "pooled", "A", "arm_mean"),
          2: g("confirm", "degeneracy", SS, DEG, "pooled", "A", "comparator_mean"),
          3: g("confirm", "degeneracy", ST, DEG, "pooled", "A", "comparator_mean"),
          4: g("confirm", "degeneracy", SH, DEG, "pooled", "A", "comparator_mean")},
    227: {0: g("meta", "n_by_half", "confirm"),
          2: g("prior_ex1_vs_stageO_X1", "real_ex1==kn/lin/G0/X1 (face mid, |d|<=1e-9)", src="c10"),
          3: g("prior_ex1_vs_stageO_X1", "n", src="c10"),
          4: g("prior_confirm_real_ex1_minus_two_tier_ex1", "A", "mean", src="c10"),
          5: g("prior_confirm_real_ex1_minus_two_tier_ex1", "A", "ci", 0, src="c10"),
          6: g("prior_confirm_real_ex1_minus_two_tier_ex1", "A", "ci", 1, src="c10")},
    228: {1: g("grp00000", "n_rows", src="c8"), 2: g("grp00000", "n_rows", src="c8"), 3: g("meta", "n_by_half", "confirm"),
          4: grp_share("min"), 5: grp_share("max"),
          6: g("grp00000", "contrasts", SS, "taf_face_mid", "without_grp00000_A", "mean", src="c8"),
          7: g("grp00000", "contrasts", SS, "taf_face_mid", "without_grp00000_A", "ci", 0, src="c8"),
          8: g("grp00000", "contrasts", SS, "taf_face_mid", "without_grp00000_A", "ci", 1, src="c8"),
          9: g("grp00000", "contrasts", SS, "taf_face_mid", "grp00000_only_A", src="c8"),
          10: logo("min_A"), 11: logo("max_A")},
    229: {0: lambda S: 100 * S.rise["by_half"]["confirm"]["share_gt_1e-9_rel"],
          1: lambda S: 100 * (S.c["c2_c3"]["selected_row_ratios"]["confirm"]["flux_exchange"]["quantiles_0_1_5_25_50_75_95_99_100"][8] - 1)},
}
SPEC.pop(54)   # line 54 is handled with line 70 (prior-LP rows), see MANUAL


def fmt_like(x, tok):
    dec = len(tok.split(".")[1]) if "." in tok else 0
    s = f"{float(x):.{dec}f}"
    return s.lstrip("-") if float(s) == 0 else s


def main():
    ap = argparse.ArgumentParser()
    for k in ("doc", "old-gen", "new-gen", "old", "new", "out-doc", "report", "old-rise", "new-rise"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    So, Sn = Src(a.old, a.old_rise), Src(a.new, a.new_rise)
    doc = open(a.doc).read().splitlines(keepends=True)
    og = open(a.old_gen).read().splitlines(keepends=True)
    ng = open(a.new_gen).read().splitlines(keepends=True)
    rep = {"skeleton_mismatch": [], "gen_lines_replaced": 0, "gen_lines_changed": [], "control_fail": [],
           "changed": [], "uncovered": []}
    if len(og) != len(ng):
        rep["skeleton_mismatch"].append(f"line counts {len(og)} vs {len(ng)}")
    for i, (x, y) in enumerate(zip(og, ng)):
        if SKEL.sub("", x) != SKEL.sub("", y):
            rep["skeleton_mismatch"].append({"gen_line": i + 1, "old": x.strip()[:160], "new": y.strip()[:160]})
    start = next(i for i, l in enumerate(doc) if l.startswith("## 1. What was run"))
    gstart = next(i for i, l in enumerate(og) if l.startswith("## 1. What was run"))
    covered = set(SPEC)
    for i in range(start, len(doc)):
        ln = i + 1
        if ln in SPEC:
            continue
        j = gstart + (i - start)
        # find the identical old generated line (same offset first, else anywhere unique)
        cand = [j] if j < len(og) and og[j] == doc[i] else [k for k, l in enumerate(og) if l == doc[i]]
        if len(cand) >= 1 and NUM.search(doc[i]):
            k = cand[0]
            if k < len(ng) and SKEL.sub("", og[k]) == SKEL.sub("", ng[k]):
                rep["gen_lines_replaced"] += 1
                if ng[k] != doc[i]:
                    rep["gen_lines_changed"].append({"line": ln, "old": doc[i].strip()[:200], "new": ng[k].strip()[:200]})
                doc[i] = ng[k]
                covered.add(ln)
    for ln, spec in SPEC.items():
        line = doc[ln - 1]
        out, last = [], 0
        for i, m in enumerate(NUM.finditer(line)):
            gtr = spec.get(i)
            if gtr is None:
                continue
            sign, intp, frac = m.groups()
            tok = norm(sign, intp, frac)
            try:
                old = fmt_like(gtr(So), tok)
            except Exception as e:  # noqa: BLE001
                rep["control_fail"].append({"line": ln, "i": i, "tok": m.group(0), "error": repr(e)})
                continue
            if float(old) != float(tok):
                rep["control_fail"].append({"line": ln, "i": i, "tok": m.group(0), "old_value": old})
                continue
            try:
                new = fmt_like(gtr(Sn), tok)
            except Exception as e:  # noqa: BLE001
                rep["control_fail"].append({"line": ln, "i": i, "tok": m.group(0), "new_error": repr(e)})
                continue
            if float(new) != float(tok) or new.lstrip("-") != tok.lstrip("-"):
                rep["changed"].append({"line": ln, "i": i, "old": m.group(0), "new": restyle(new, sign, intp)})
                out.append(line[last:m.start()] + restyle(new, sign, intp))
                last = m.end()
        out.append(line[last:])
        doc[ln - 1] = "".join(out)
    for i, l in enumerate(doc):
        if i + 1 not in covered and i >= 10 and NUM.search(l):
            rep["uncovered"].append({"line": i + 1, "text": l.strip()[:200]})
    open(a.out_doc, "w").write("".join(doc))
    json.dump(rep, open(a.report, "w"), indent=1, ensure_ascii=False)
    print(f"skeleton mismatches {len(rep['skeleton_mismatch'])}; generated lines replaced {rep['gen_lines_replaced']} "
          f"(changed {len(rep['gen_lines_changed'])}); control failures {len(rep['control_fail'])}; "
          f"hand tokens changed {len(rep['changed'])}; uncovered lines {len(rep['uncovered'])}")


if __name__ == "__main__":
    main()
