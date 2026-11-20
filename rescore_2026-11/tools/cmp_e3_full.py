#!/usr/bin/env python
"""Compare two E3-gate runs (re-score R2.1: C-E3 control vs committed; smoke1 vs smoke2).

    cmp_e3_full.py --a-json A.json --b-json B.json [--a-csv A.csv.gz --b-csv B.csv.gz] [--out report.json]

JSON: every top-level block except _meta compared exactly (diagnostics with timing keys
'*t_e3_ms*'/'*t_select*' masked); _meta: the content keys (store, reproduction, neighbours, keys,
streamed, provider misses, catalogue.is_current). CSV: same key set ((key, mode) rows), same columns,
every column except t_select_ms/t_e3_ms equal. Exit 0 iff everything compared is identical.
"""
import argparse, csv, gzip, json, sys

TIMING_COLS = {"t_select_ms", "t_e3_ms"}


def mask(x):
    if isinstance(x, dict):
        return {k: ("<timing>" if ("t_e3_ms" in k or "t_select" in k) else mask(v)) for k, v in x.items()}
    if isinstance(x, list):
        return [mask(v) for v in x]
    return x


def diff(x, y, p=""):
    if type(x) != type(y):
        return [(p, x, y)]
    if isinstance(x, dict):
        d = [(p + "/" + k, "only one side", None) for k in sorted(set(x) ^ set(y))]
        for k in sorted(set(x) & set(y)):
            d += diff(x[k], y[k], p + "/" + k)
        return d
    if isinstance(x, list):
        if len(x) != len(y):
            return [(p, "len %d" % len(x), "len %d" % len(y))]
        d = []
        for i, (u, v) in enumerate(zip(x, y)):
            d += diff(u, v, "%s[%d]" % (p, i))
        return d
    return [] if x == y else [(p, x, y)]


def rows(p):
    with gzip.open(p, "rt", newline="") as fh:
        r = csv.DictReader(fh)
        return r.fieldnames, {(x["key"], x["mode"]): x for x in r}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a-json", required=True); ap.add_argument("--b-json", required=True)
    ap.add_argument("--a-csv"); ap.add_argument("--b-csv"); ap.add_argument("--out")
    a = ap.parse_args()
    A, B = json.load(open(a.a_json)), json.load(open(a.b_json))
    out, ok = {"a": a.a_json, "b": a.b_json, "blocks": {}}, True
    for k in sorted((set(A) | set(B)) - {"_meta"}):
        d = diff(mask(A.get(k)), mask(B.get(k)), k)
        out["blocks"][k] = {"identical": not d, "n_diff": len(d), "first": d[:8]}
        ok &= not d
    ma, mb = A["_meta"], B["_meta"]
    for k in ("store", "reproduction_vs_gene_capture_cells", "neighbours_absent_from_cache", "n_selection_keys",
              "n_amplicons_streamed", "selection_keys_without_hit_record", "provider_misses", "pass_criterion"):
        same = ma.get(k) == mb.get(k)
        out["blocks"]["_meta." + k] = {"identical": same, "a": ma.get(k) if not same else "==", "b": mb.get(k) if not same else "=="}
        ok &= same
    out["a_catalogue_is_current"] = ma.get("catalogue", {}).get("is_current")
    out["b_catalogue_is_current"] = mb.get("catalogue", {}).get("is_current")
    out["a_reproduction"] = ma.get("reproduction_vs_gene_capture_cells", {}).get("counts")
    if a.a_csv:
        fa, ra = rows(a.a_csv); fb, rb = rows(a.b_csv)
        cols = [c for c in fa if c not in TIMING_COLS]
        nd = [(k, c, ra[k][c], rb[k][c]) for k in sorted(set(ra) & set(rb)) for c in cols if ra[k][c] != rb.get(k, {}).get(c)]
        out["csv"] = {"n_a": len(ra), "n_b": len(rb), "same_keys": set(ra) == set(rb), "same_columns": fa == fb,
                      "columns_compared": len(cols), "n_cell_diffs": len(nd),
                      "n_rows_differ": len({x[0] for x in nd}), "first": nd[:10]}
        ok &= (set(ra) == set(rb)) and fa == fb and not nd
    out["IDENTICAL"] = ok
    s = json.dumps(out, indent=1, default=str)
    if a.out:
        open(a.out, "w").write(s)
    print(s[:6000])
    print("RESULT:", "IDENTICAL" if ok else "DIFFERENT")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
