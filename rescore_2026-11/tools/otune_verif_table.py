#!/usr/bin/env python
"""Key results of the Stage O independent verification (data/objective_tuning/verification/out_*.json), old -> new,
as a markdown table for the dated re-score section of VERIFICATION.md.
usage: otune_verif_table.py OLD_VERIF_DIR NEW_VERIF_DIR [OUT_MD]"""
import json
import os
import re
import sys


def load(d):
    return {f[4:-5]: json.load(open(os.path.join(d, f))) for f in os.listdir(d) if f.startswith("out_") and f.endswith(".json")}


def ci(x, k=4):
    m, c = (x["mean"], x["ci"]) if isinstance(x, dict) else (x[0], x[1])
    return f"{m:+.{k}f} [{c[0]:+.{k}f}, {c[1]:+.{k}f}]".replace("-", "−")


ROWS = [
    ("confirm-half numbers compared with summary.json; max abs diff",
     lambda c: f"{c['c4_c5']['cmp_vs_summary']['n_numbers_compared']}; {c['c4_c5']['cmp_vs_summary']['max_abs_diff']}"),
    ("rows tune / confirm; groups tune / confirm; groups spanning halves",
     lambda c: "{} / {}; {} / {}; {}".format(c["c1_c6"]["claim1"]["rows_by_half_recomputed"]["tune"],
                                             c["c1_c6"]["claim1"]["rows_by_half_recomputed"]["confirm"],
                                             c["c1_c6"]["claim1"]["groups_by_half"]["tune"],
                                             c["c1_c6"]["claim1"]["groups_by_half"]["confirm"],
                                             len(c["c1_c6"]["claim1"]["groups_spanning_halves"]))),
    ("selected (tune half, pooled A) and naive joint selection",
     lambda c: "{} ({:.4f}); naive {}".format(c["c2_c3"]["selection"]["selected"], c["c2_c3"]["selection"]["pooled_A"],
                                            c["c2_c3"]["naive_selection"]["selected"])),
    ("tune-half selection bootstrap: share selecting the selected arm",
     lambda c: "{:.3f}".format(c["c12"]["per_stratum_seed1729"]["winner_freq"].get("est/lin/G0/X1", 0) / 2000)),
    ("selected − shipped, confirm, pooled A (face mid)", lambda c: c["c12"]["confirm_selected_minus_shipped"]["A"]),
    ("same, tune half", lambda c: c["c12"]["tune_selected_minus_shipped"]["A"]),
    ("n = 1, confirm / tune", lambda c: "{:+.4f} / {:+.4f}".format(c["c12"]["confirm_selected_minus_shipped"]["per"]["n1"][0],
                                                                     c["c12"]["tune_selected_minus_shipped"]["per"]["n1"][0])),
    ("leave one group out: pooled-A range", lambda c: "[{:+.4f}, {:+.4f}]".format(c["c8"]["logo_confirm"]["selected-shipped"]["min_A"][0],
                                                                                   c["c8"]["logo_confirm"]["selected-shipped"]["max_A"][0])),
    ("grp00000 confirm rows; without it (A)",
     lambda c: "{}; {}".format(c["c8"]["grp00000"]["n_rows"],
                               ci(c["c8"]["grp00000"]["contrasts"]["selected-shipped"]["taf_face_mid"]["without_grp00000_A"]))),
    ("excluding the 2 prior-match failures / the rows whose D is not bitwise the prior's (A)",
     lambda c: "{} / {}".format(ci(c["c8"]["exclusions"]["excl_2_prior_match_false"]["selected-shipped"]["A"]),
                                ci(c["c8"]["exclusions"]["excl_3_D_not_bitwise"]["selected-shipped"]["A"]))),
    ("rows whose rebuilt LP differs from the prior run (D not bitwise)",
     lambda c: ", ".join(d[0] for d in c["c1_c6"]["claim6"]["D_not_bitwise_prior"])),
    ("g_meas from tune rows (num / den)", lambda c: "{} / {} = {:.6f}".format(*c["c10"]["g_meas_recomputed"]["tune"])),
    ("FVA pass: rows; metric mismatches vs the main run", lambda c: "{}; {}".format(c["c10"]["fva_vs_main"]["n_rows"],
                                                                                   c["c10"]["fva_vs_main"]["n_metric_mismatches_vs_main_run"])),
    ("volume term of selected − shipped at the returned vertex (A)", lambda c: ci(c["c13"]["selected-shipped"]["volume_A"])),
    ("headroom captured, confirm A: selected own cell / shipped / selected in shipped frame",
     lambda c: "{:.3f} / {:.3f} / {:.3f}".format(c["c4_c5"]["headroom"]["A"]["selected_own_cell"],
                                                 c["c4_c5"]["headroom"]["A"]["shipped"],
                                                 c["c4_c5"]["headroom"]["A"]["selected_in_shipped_frame"])),
]


def main():
    o, n = load(sys.argv[1]), load(sys.argv[2])
    lines = ["| check | 2026-11-03 | 2026-11-16 |", "|---|---|---|"]
    for name, f in ROWS:
        def val(c):
            try:
                return re.sub(r"(?<![A-Za-z0-9_])-(?=\d)", "−", str(f(c)))
            except Exception as e:  # noqa: BLE001
                return f"n/a ({type(e).__name__}: {e})"
        a, b = val(o), val(n)
        lines.append(f"| {name} | {a} | {b if b != a else 'unchanged'} |")
    out = "\n".join(lines) + "\n"
    if len(sys.argv) > 3:
        open(sys.argv[3], "w").write(out)
    print(out)


if __name__ == "__main__":
    main()
