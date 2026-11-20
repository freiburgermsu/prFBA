#!/usr/bin/env python
"""Cell-level attribution check for the 2026-11 re-score (plan W8).

    attrib_diff.py --table {genecap,mag,comparison,sr} --old OLD --new NEW
                   [--repaired $W/expected/repaired_ids.json] [--d20]
                   [--reps-pickle P ...] [--selection S ...] [--legacy-selection S ...]
                   [--sr-prfba-dir DIR ...] [--affected FILE] [--out report.json]

Joins the old and new cell tables on their key columns and asserts that every row whose values
changed can be attributed to the repaired BV-BRC data:

  source_repaired     the row's source genome is a repaired genome (its truth changed)
  anchor_repaired     the stored anchor (MAG cells) is repaired
  rep_repaired        a selected representative (old or new, from a reps provider) is repaired
  reference_changed   SR arms that read the fold reference (PICRUSt2, baselines) in a fold/density
                      listed in --affected (its reference traits changed)
  multi_pick_unstored n_selected (or n_reps) >= 2 in old or new and no reps are known: a repaired
                      non-anchor pick cannot be excluded

Anything else is `unattributed`. Rows present on one side only, and rows whose identity columns
(those no gene-data change can move) differ, are failures too. The repaired set is
`repaired_2026_11_13` of --repaired (the 1,006 shared-cache entries replaced on 2026-11-13) plus,
with --d20, `poisoned_d20` (the 20 poisoned keys the 09-17 depth-20 chain saw as empty).
Exit 0 iff there is no unattributed, one-sided or identity-changed row.

Reps providers (all optional; old and new are unioned):
  --reps-pickle      the 11-13 re-score pickles (cells_old.pkl / cells_new.pkl: '_reps' per cell)
  --selection        prfba_selection.json (prfba_include / prfba_exclude methods)
  --legacy-selection prfba_selection_legacy.json (the *_legacy methods)
  --sr-prfba-dir     a shared-reference data/prfba dir (per fold/density/mode predictions.json reps)
"""
import argparse
import gzip
import json
import os
import pickle
import sys
import time

import pandas as pd

W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
RANKS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

PROFILES = {
    "genecap": {
        "keys": ["ampliconKey", "self_mode"],
        "source": "src_genome_id", "anchor": None, "nsel": ["n_selected", "n_reps"],
        "identity": ["src_genome_id", "src_taxon_id", "region", "domain", "best_identity",
                     "ident_bin", "best_nonself_identity"],
    },
    "mag": {
        "keys": ["mag", "regime"],
        "source": "genome_id", "anchor": "anchor", "nsel": ["n_selected"],
        "identity": ["genome_id", "completeness", "contamination", "ani_class", "tier", "best_ani",
                     "n_hits", "n_selected", "anchor"] + RANKS,
    },
    "comparison": {
        "keys": ["amplicon", "method", "space"],
        "source": None, "anchor": None, "nsel": ["n_reps"],
        "identity": ["region", "domain", "nsti", "nsti_bin", "picrust2_domain", "src_in_picrust2_ref",
                     "best_nonself_identity", "ident_bin"],
    },
    "sr": {
        "keys": ["amplicon", "fold", "density", "arm", "space"],
        "source": "src_genome_id", "anchor": None, "nsel": ["n_reps"],
        "identity": ["md5", "region", "domain", "src_genome_id", "group_id", "nsti", "nsti_bin",
                     "picrust2_domain", "routing", "rank1_identity", "rank1_rank", "ident_bin",
                     "novelty_bin"],
    },
}


def read_table(path):
    comp = "gzip" if path.endswith(".gz") else None
    return pd.read_csv(path, dtype=str, keep_default_na=False, compression=comp)


def load_repaired(path, d20):
    d = json.load(open(path))
    s = set(d["repaired_2026_11_13"])
    if d20:
        s |= set(d.get("poisoned_d20", []))
    return s


# ------------------------------------------------------------------ reps providers
def reps_from_pickles(paths):
    """{(ampliconKey, self_mode): set(reps)} from the 11-13 re-score pickles."""
    out = {}
    for p in paths:
        d = pickle.load(open(p, "rb"))
        cells = d["cells"] if isinstance(d, dict) and "cells" in d else d
        for key, cs in cells.items():
            for c in (cs if isinstance(cs, (list, tuple)) else [cs]):
                r = c.get("_reps")
                if r is not None:
                    out.setdefault((key, c["self_mode"]), set()).update(r)
    return out


def reps_from_selections(default_paths, legacy_paths):
    """{(amplicon, method): set(reps)} from prfba_selection{,_legacy}.json."""
    import orjson
    out = {}
    for paths, suffix in ((default_paths, ""), (legacy_paths, "_legacy")):
        for p in paths:
            sel = orjson.loads(open(p, "rb").read())
            for amp, rec in sel.items():
                if not isinstance(rec, dict):
                    continue
                for mode, method in (("include", "prfba_include"), ("exclude_0.987", "prfba_exclude")):
                    r = (rec.get(mode) or {}).get("reps")
                    if r:
                        out.setdefault((amp, method + suffix), set()).update(r)
    return out


def reps_from_sr(prfba_dirs):
    """{(md5, fold, density, arm): set(reps)} from shared-reference predictions.json files."""
    import glob
    import orjson
    out = {}
    for d in prfba_dirs:
        for pj in sorted(glob.glob(os.path.join(d, "f*_d*", "predictions.json"))):
            run = os.path.basename(os.path.dirname(pj))
            forced = run.endswith("_forced")
            p = orjson.loads(open(pj, "rb").read())
            fold = str(p["meta"]["fold"])
            dens = p["meta"]["density"]
            for mode, preds in p["pred"].items():
                arm = f"prfba_{mode}" + ("_forced" if forced else "")
                for md5, rec in preds.items():
                    r = rec.get("reps")
                    if r:
                        out.setdefault((md5, fold, _dens_key(dens), arm), set()).update(r)
    return out


def _dens_key(x):
    return f"{float(x):.4f}"


# ------------------------------------------------------------------ core
def attribute(old, new, prof, repaired, reps=None, rep_key=None, affected=None, max_examples=20,
              extra_genomes=None, extra_keys=None):
    keys = prof["keys"]
    for name, df in (("old", old), ("new", new)):
        dup = int(df.duplicated(keys).sum())
        if dup:
            raise SystemExit(f"[attrib] {dup} duplicate keys in the {name} table on {keys}")
    cols_old, cols_new = list(old.columns), list(new.columns)
    rep = {"n_old": len(old), "n_new": len(new), "columns_equal": cols_old == cols_new}
    common_cols = [c for c in cols_old if c in cols_new]
    m = old[common_cols].merge(new[common_cols], on=keys, how="outer", suffixes=("_o", "_n"),
                               indicator=True)
    only_old = m[m["_merge"] == "left_only"]
    only_new = m[m["_merge"] == "right_only"]
    both = m[m["_merge"] == "both"]
    rep.update(n_common=len(both), n_only_old=len(only_old), n_only_new=len(only_new),
               only_old_examples=only_old[keys].head(max_examples).to_dict("records"),
               only_new_examples=only_new[keys].head(max_examples).to_dict("records"))
    ident = [c for c in prof["identity"] if c in common_cols and c not in keys]
    values = [c for c in common_cols if c not in keys and c not in ident]
    id_mask = pd.Series(False, index=both.index)
    for c in ident:
        id_mask |= both[c + "_o"] != both[c + "_n"]
    rep["identity_columns"] = ident
    rep["n_identity_changed"] = int(id_mask.sum())
    rep["identity_changed_examples"] = both.loc[id_mask, keys].head(max_examples).to_dict("records")
    ch_mask = pd.Series(False, index=both.index)
    col_counts = {}
    for c in values:
        d = both[c + "_o"] != both[c + "_n"]
        col_counts[c] = int(d.sum())
        ch_mask |= d
    rep["value_columns"] = values
    rep["changed_by_column"] = {k: v for k, v in col_counts.items() if v}
    changed = both[ch_mask]
    rep["n_changed"] = len(changed)

    first, anyc = {}, {}
    unattributed = []
    src_col, anc_col = prof["source"], prof["anchor"]
    for _, r in changed.iterrows():
        reasons = []
        src = r[src_col + "_o"] if src_col and src_col not in keys else (r[src_col] if src_col else None)
        if src is None:                                   # comparison: the source is in the key
            src = r["amplicon"].split("__")[0]
        if src in repaired:
            reasons.append("source_repaired")
        if anc_col:
            a = {r[anc_col + "_o"], r[anc_col + "_n"]} if anc_col not in keys else {r[anc_col]}
            if a & repaired:
                reasons.append("anchor_repaired")
        known = None
        if reps is not None and rep_key is not None:
            known = reps.get(rep_key(r))
            if known and known & repaired:
                reasons.append("rep_repaired")
        # documented non-repair reasons (after the repair reasons, so they only count when no
        # repair reason applies first): extra genome sets (source / rep) and extra key sets
        for name, gs in (extra_genomes or {}).items():
            if src in gs:
                reasons.append(f"source_{name}")
            if known and known & gs:
                reasons.append(f"rep_{name}")
        for name, ks in (extra_keys or {}).items():
            if r[keys[0]] in ks or tuple(r[k] for k in keys[:2]) in ks:
                reasons.append(name)
        if affected is not None and prof is PROFILES["sr"] and not r["arm"].startswith("prfba_"):
            if f"{r['fold']}|{_dens_key(r['density'])}" in affected:
                reasons.append("reference_changed")
        nsel = 0
        for c in prof["nsel"]:
            for suf in ("_o", "_n"):
                v = r.get(c + suf, "")
                try:
                    nsel = max(nsel, int(float(v)))
                except (TypeError, ValueError):
                    pass
        if not reasons and nsel >= 2 and not known:
            reasons.append("multi_pick_unstored")
        if not reasons:
            reasons.append("unattributed" if known is not None or nsel == 0
                           else "unattributed_reps_unknown")
        first[reasons[0]] = first.get(reasons[0], 0) + 1
        for x in reasons:
            anyc[x] = anyc.get(x, 0) + 1
        if reasons[0].startswith("unattributed") and len(unattributed) < max_examples:
            ex = {k: r[k] for k in keys}
            ex.update({c: [r[c + "_o"], r[c + "_n"]] for c in values if r[c + "_o"] != r[c + "_n"]})
            ex["n_selected_max"] = nsel
            unattributed.append(ex)
    rep["by_reason_first"] = dict(sorted(first.items()))
    rep["by_reason_any"] = dict(sorted(anyc.items()))
    rep["n_unattributed"] = sum(v for k, v in first.items() if k.startswith("unattributed"))
    rep["unattributed_examples"] = unattributed
    rep["ok"] = (rep["n_unattributed"] == 0 and rep["n_identity_changed"] == 0
                 and rep["n_only_old"] == 0 and rep["n_only_new"] == 0)
    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, choices=sorted(PROFILES))
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--repaired", default=os.path.join(W, "expected", "repaired_ids.json"))
    ap.add_argument("--d20", action="store_true", help="also count the 20 d20 poisoned keys as repaired")
    ap.add_argument("--reps-pickle", action="append", default=[])
    ap.add_argument("--selection", action="append", default=[])
    ap.add_argument("--legacy-selection", action="append", default=[])
    ap.add_argument("--sr-prfba-dir", action="append", default=[])
    ap.add_argument("--affected", help="SR: JSON list of 'fold|density' whose reference traits changed")
    ap.add_argument("--extra-genomes", action="append", default=[], metavar="NAME=FILE",
                    help="documented non-repair genome set (JSON list): reasons source_NAME / rep_NAME")
    ap.add_argument("--extra-keys", action="append", default=[], metavar="NAME=FILE",
                    help="documented non-repair key set (JSON list of first-key values or "
                         "[key1, key2] pairs): reason NAME")
    ap.add_argument("--out")
    ap.add_argument("--max-examples", type=int, default=20)
    a = ap.parse_args(argv)
    t0 = time.time()
    prof = PROFILES[a.table]
    repaired = load_repaired(a.repaired, a.d20)
    old, new = read_table(a.old), read_table(a.new)
    reps = rep_key = None
    if a.table == "genecap" and a.reps_pickle:
        reps = reps_from_pickles(a.reps_pickle)
        rep_key = lambda r: (r["ampliconKey"], r["self_mode"])  # noqa: E731
    elif a.table == "comparison" and (a.selection or a.legacy_selection):
        reps = reps_from_selections(a.selection, a.legacy_selection)
        rep_key = lambda r: (r["amplicon"], r["method"])  # noqa: E731
    elif a.table == "sr" and a.sr_prfba_dir:
        reps = reps_from_sr(a.sr_prfba_dir)
        rep_key = lambda r: (r["md5_o"], r["fold"], _dens_key(r["density"]), r["arm"])  # noqa: E731
    affected = None
    if a.affected:
        affected = {f"{x.split('|')[0]}|{_dens_key(x.split('|')[1])}" for x in json.load(open(a.affected))}
    extra_genomes = {x.split("=", 1)[0]: set(json.load(open(x.split("=", 1)[1]))) for x in a.extra_genomes}
    extra_keys = {}
    for x in a.extra_keys:
        name, path = x.split("=", 1)
        extra_keys[name] = {tuple(v) if isinstance(v, list) else v for v in json.load(open(path))}
    rep = attribute(old, new, prof, repaired, reps, rep_key, affected, a.max_examples,
                    extra_genomes, extra_keys)
    rep["extra_genomes"] = {k: len(v) for k, v in extra_genomes.items()}
    rep["extra_keys"] = {k: len(v) for k, v in extra_keys.items()}
    rep = dict({"table": a.table, "old": os.path.abspath(a.old), "new": os.path.abspath(a.new),
                "repaired": os.path.abspath(a.repaired), "n_repaired": len(repaired), "d20": a.d20,
                "reps_provider": bool(reps), "n_reps_keys": len(reps) if reps else 0,
                "seconds": None}, **rep)
    rep["seconds"] = round(time.time() - t0, 1)
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(rep, fh, indent=1, default=str)
    brief = {k: rep[k] for k in ("table", "n_old", "n_new", "n_only_old", "n_only_new",
                                 "n_identity_changed", "n_changed", "by_reason_first", "n_unattributed", "ok")}
    print(json.dumps(brief))
    print("ATTRIBUTION PASS" if rep["ok"] else "ATTRIBUTION FAIL")
    return 0 if rep["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
