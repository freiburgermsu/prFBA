#!/usr/bin/env python
"""Audit and correct ModelSEED template reaction directions against the ModelSEED
database's irreversible directions.

audit:  for every v7.0 template reaction the database marks irreversible ('>' or '<'),
        check that the template writes the reaction in the database's orientation
        (coefficient signs of the shared compounds, protons ignored; a reaction written
        reversed has its target mirrored, a mixed or empty comparison is skipped), then
        classify each of the two direction fields:
          match                 the template already allows only the database direction
          template_opposite     the template allows only the forbidden direction (tier B)
          template_reversible   the template allows both directions (tier C)
        Also flags non-default bounds and checks whether negative reverse penalties sit
        on reactions whose forward direction is thermodynamically favoured.

apply:  writes patched copies of the template files for a tier set / field set /
        template subset, editing only the direction fields of the selected reactions
        (byte-minimal, so a word diff stays readable), and asserts after a full parse
        that nothing else changed.

table:  per-reaction table of one tier (what it changes per field and template group).

    python template_direction_audit.py audit --templates DIR --out OUTDIR
    python template_direction_audit.py apply --templates DIR --audit OUTDIR/direction_audit.tsv \
        --tiers B,C --fields direction,GapfillDirection [--only Core-V6.json] \
        [--reactions rxn00206,rxn00380] --dest DIR2
    python template_direction_audit.py table --templates DIR --audit OUTDIR/direction_audit.tsv \
        --tier B --out OUTDIR/tierB_reactions.tsv
"""
import argparse
import csv
import glob
import json
import os
import re
import shutil
import statistics

MSD = "/home/freiburger/Documents/ModelSEEDDatabase/Biochemistry"
COMP = {"0": "c", "1": "e", "2": "p"}
H = "cpd00067"
DEFAULT_BOUNDS = {-1000, 0, 1000}
VALS = {">": {"direction": '">"', "GapfillDirection": '">"', "maxforflux": "1000",
              "maxrevflux": "0.0", "lower_bound": "0", "upper_bound": "1000"},
        "<": {"direction": '"<"', "GapfillDirection": '"<"', "maxforflux": "0",
              "maxrevflux": "1000.0", "lower_bound": "-1000", "upper_bound": "0"}}
BOUND_FIELDS = ("maxforflux", "maxrevflux", "lower_bound", "upper_bound")
MIRROR = {">": "<", "<": ">"}
START = re.compile(r'\{\s*"GapfillDirection"')


def load_db(msd=None):
    msd = msd or MSD
    db = {}
    for f in sorted(glob.glob(os.path.join(msd, "reaction_*.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                db[r["id"]] = r
    grades = {}
    gf = os.path.join(msd, "Thermodynamics/SourceGrading/results/thermo_grades/reaction_grades.tsv")
    with open(gf) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            grades[r["rxn"]] = r
    return db, grades


def db_stoich(s):
    out = {}
    for part in (s or "").split(";"):
        bits = part.split(":")
        if len(bits) < 3:
            continue
        try:
            out[(bits[1], COMP.get(bits[2], bits[2]))] = float(bits[0])
        except ValueError:
            continue
    return out


def tmpl_stoich(r):
    out = {}
    for x in r.get("templateReactionReagents", []):
        cid = x["templatecompcompound_ref"].rsplit("/", 1)[-1]
        cpd, comp = cid.rsplit("_", 1)
        out[(cpd, comp)] = float(x["coefficient"])
    return out


def orientation(tr, dr):
    shared = [k for k in tr if k in dr and k[0] != H]
    if not shared:
        return 0, 0
    agree = sum(1 for k in shared if (tr[k] > 0) == (dr[k] > 0))
    if agree == len(shared):
        return 1, len(shared)
    if agree == 0:
        return -1, len(shared)
    return 0, len(shared)


def best_dg(g):
    if not g:
        return None, ""
    src = (g.get("best_source") or "").strip()
    col = {"EQ": "dg_EQ", "GC": "dg_GC", "DGP": "dg_DGP", "openTECR": "dg_openTECR"}.get(src)
    for c in ([col] if col else []) + ["dg_EQ", "dg_openTECR", "dg_DGP", "dg_GC"]:
        v = (g.get(c) or "").strip()
        if v:
            try:
                return float(v), g.get("best_grade") or ""
            except ValueError:
                pass
    return None, g.get("best_grade") or ""


def kind(cur, target):
    if cur == target:
        return "match"
    if cur == "=":
        return "template_reversible"
    if cur in (">", "<"):
        return "template_opposite"
    return "other"


def cmd_audit(a):
    db, grades = load_db(a.db)
    rows, pen = [], []
    for f in sorted(glob.glob(os.path.join(a.templates, "*.json"))):
        t = json.load(open(f))
        for r in t["reactions"]:
            base = r["id"].rsplit("_", 1)[0]
            d = db.get(base)
            if not d:
                continue
            o, nshared = orientation(tmpl_stoich(r), db_stoich(d["stoichiometry"]))
            dg, grade = best_dg(grades.get(base))
            dg_t = None if dg is None or o == 0 else dg * o     # dG in the template's orientation
            if (r.get("reverse_penalty") or 0) < 0:
                # ModelSEEDpy builds a reverse gapfill candidate (and so uses the reverse
                # penalty) only when the template's lower bound, which follows `direction`,
                # is negative and GapfillDirection is not '>'
                lb = r.get("lower_bound", -1000 if r.get("direction") in ("=", "<") else 0)
                reach = lb < 0 and r.get("GapfillDirection") != ">"
                cheaper = abs(r.get("base_cost", 0) + r["reverse_penalty"]) < abs(
                    r.get("base_cost", 0) + (r.get("forward_penalty") or 0))
                pen.append((r["reverse_penalty"], dg_t, grade, reach, cheaper, os.path.basename(f)))
            if d["reversibility"] not in (">", "<"):
                continue
            target = None if o == 0 else (d["reversibility"] if o == 1 else MIRROR[d["reversibility"]])
            nondef = [k for k in ("lower_bound", "upper_bound")
                      if r.get(k) is not None and r[k] not in DEFAULT_BOUNDS]
            for fld in ("direction", "GapfillDirection"):
                cur = r.get(fld)
                k = "orientation_unresolved" if target is None else kind(cur, target)
                if k == "match":
                    continue
                rows.append({"template": os.path.basename(f), "reaction": r["id"], "field": fld,
                             "current": cur, "db_reversibility": d["reversibility"],
                             "orientation": o, "n_shared_compounds": nshared,
                             "target": target or "", "kind": k,
                             "tier": {"template_opposite": "B", "template_reversible": "C"}.get(k, "-"),
                             "nondefault_bounds": "|".join(nondef),
                             "dG_template_orientation": "" if dg_t is None else round(dg_t, 3),
                             "thermo_grade": grade, "name": r.get("name", "")})
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "direction_audit.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    # negative reverse penalties vs the thermodynamically favoured direction, over all entries
    # and over the entries the gapfiller can actually use (a reverse candidate exists)
    def pen_stats(sub):
        dg = [(p, g) for p, g, gr, *_ in sub if g is not None and gr]
        fwd = sum(1 for p, g in dg if g < 0)
        return {"n_entries": len(sub), "n_with_dG": len(dg),
                "n_reverse_cheaper": sum(1 for x in sub if x[4]),
                "n_forward_favoured": fwd, "n_reverse_favoured": len(dg) - fwd,
                "share_forward_favoured": round(fwd / len(dg), 4) if dg else None,
                "spearman_penalty_vs_dG": spearman([p for p, _ in dg], [g for _, g in dg])}
    by_t = {}
    for x in pen:
        by_t.setdefault(x[5], []).append(x)
    summary = {"rows": len(rows), "by_template_field_kind": {}, "negative_reverse_penalty": {
        "all": pen_stats(pen), "reachable": pen_stats([x for x in pen if x[3]]),
        "by_template": {t: {"all": pen_stats(v), "reachable": pen_stats([x for x in v if x[3]])}
                        for t, v in sorted(by_t.items())}}}
    for r in rows:
        key = "%s|%s|%s" % (r["template"], r["field"], r["kind"])
        summary["by_template_field_kind"][key] = summary["by_template_field_kind"].get(key, 0) + 1
    json.dump(summary, open(os.path.join(a.out, "direction_audit_summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


def spearman(x, y):
    if len(x) < 3:
        return None
    rank = lambda v: {i: r for r, i in enumerate(sorted(range(len(v)), key=v.__getitem__))}  # noqa: E731
    rx, ry = rank(x), rank(y)
    mx, my = statistics.mean(rx.values()), statistics.mean(ry.values())
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(len(x)))
    den = (sum((rx[i] - mx) ** 2 for i in range(len(x))) * sum((ry[i] - my) ** 2 for i in range(len(x)))) ** 0.5
    return round(num / den, 4) if den else None


def reaction_spans(s):
    """{reaction id: (start, end)} of every top-level reaction object in a template's text."""
    spans, pos = {}, 0
    for m in START.finditer(s):
        if m.start() < pos:
            continue
        depth, j = 0, m.start()
        while True:
            c = s[j]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        spans[json.loads(s[m.start():j + 1])["id"]] = (m.start(), j + 1)
        pos = j + 1
    return spans


def penalties_to_negate(t, mode, db, grades):
    """Ids of reactions whose negative reverse_penalty is negated: every one ('all'), or only
    those whose forward direction the database's dG favours ('forward_favoured')."""
    out = set()
    for r in t["reactions"]:
        if (r.get("reverse_penalty") or 0) >= 0:
            continue
        if mode == "forward_favoured":
            base = r["id"].rsplit("_", 1)[0]
            d = db.get(base)
            if not d:
                continue
            o, _ = orientation(tmpl_stoich(r), db_stoich(d["stoichiometry"]))
            dg, grade = best_dg(grades.get(base))
            if dg is None or o == 0 or not grade or dg * o >= 0:
                continue
        out.add(r["id"])
    return out


def cmd_apply(a):
    tiers = set(a.tiers.split(",")) if a.tiers else set()
    fields = set(a.fields.split(","))
    pen_mode = "all" if a.negate_negative_reverse_penalty else a.negate_reverse_penalty_where
    db = grades = None
    if pen_mode == "forward_favoured":
        db, grades = load_db()
    want = {}
    with open(a.audit) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["tier"] in tiers and r["field"] in fields and not r["nondefault_bounds"]:
                if a.only and r["template"] not in a.only.split(","):
                    continue
                if a.reactions and r["reaction"].rsplit("_", 1)[0] not in a.reactions.split(","):
                    continue
                want.setdefault(r["template"], {}).setdefault(r["reaction"], {})[r["field"]] = r["target"]
    # explicit directions (--set rxnA=>,rxnB=<): both fields, in every template carrying it
    explicit = dict(x.split("=", 1) for x in a.set.split(",")) if a.set else {}
    if explicit:
        assert all(v in (">", "<") for v in explicit.values()), explicit
        for f in sorted(glob.glob(os.path.join(a.templates, "*.json"))):
            ids = {r["id"] for r in json.load(open(f))["reactions"]}
            for base, d in explicit.items():
                if base + "_c" in ids:
                    want.setdefault(os.path.basename(f), {})[base + "_c"] = {
                        "direction": d, "GapfillDirection": d}
    os.makedirs(a.dest, exist_ok=True)
    report = {}
    for f in sorted(glob.glob(os.path.join(a.templates, "*.json"))):
        name = os.path.basename(f)
        dst = os.path.join(a.dest, name)
        if name not in want and not pen_mode:
            if os.path.abspath(f) != os.path.abspath(dst):
                shutil.copyfile(f, dst)
            continue
        want.setdefault(name, {})
        s = open(f).read()
        before = json.loads(s)
        negate = penalties_to_negate(before, pen_mode, db, grades) if pen_mode else set()
        spans = reaction_spans(s)
        assert set(spans) == {r["id"] for r in before["reactions"]}, name
        # edit each selected reaction object in place, last first so earlier offsets hold
        for rid in sorted(set(want[name]) | negate, key=lambda x: spans[x][0], reverse=True):
            start, end = spans[rid]
            new = s[start:end]
            if rid in negate:
                # reverse_penalty values track the FORWARD dG (negative where forward is favoured);
                # negating them makes the thermodynamically disfavoured reverse direction cost more
                new, n = re.subn(r'("reverse_penalty": )-([0-9.]+(?:[eE]-?[0-9]+)?)',
                                 lambda m: m.group(1) + m.group(2), new, count=1)
                assert n == 1, (name, rid)
            for fld, target in want[name].get(rid, {}).items():
                edits = {fld: VALS[target][fld]}
                if fld == "direction":            # the draft bounds follow `direction`
                    edits.update({k: VALS[target][k] for k in BOUND_FIELDS})
                for k, v in edits.items():
                    new, n = re.subn(r'("%s": )("[^"]*"|-?[0-9.]+)' % k, lambda m: m.group(1) + v,
                                     new, count=1)
                    assert n == 1, (name, rid, k)
            s = s[:start] + new + s[end:]
        after = json.loads(s)
        allowed = {"direction", "GapfillDirection", *BOUND_FIELDS}
        if negate:
            allowed.add("reverse_penalty")
        rb = {r["id"]: r for r in before["reactions"]}
        ra = {r["id"]: r for r in after["reactions"]}
        assert set(rb) == set(ra)
        nchanged = 0
        for rid in rb:
            diff = {k for k in set(rb[rid]) | set(ra[rid]) if rb[rid].get(k) != ra[rid].get(k)}
            if diff:
                assert diff <= allowed, (name, rid, diff)
                assert rid in want[name] or diff == {"reverse_penalty"}, (name, rid, diff)
                if "reverse_penalty" in diff:
                    assert rid in negate and rb[rid]["reverse_penalty"] < 0
                    assert ra[rid]["reverse_penalty"] == -rb[rid]["reverse_penalty"]
                nchanged += 1
                for fld, target in want[name].get(rid, {}).items():
                    assert ra[rid][fld] == target, (name, rid, fld)
        assert all(ra[rid]["reverse_penalty"] == -rb[rid]["reverse_penalty"] for rid in negate)
        for k in before:
            if k != "reactions":
                assert before[k] == after[k], (name, k)
        open(dst, "w").write(s)
        report[name] = {"reactions_changed": nchanged, "penalties_negated": len(negate)}
    print(json.dumps({"dest": a.dest, "reactions_changed": report}, indent=1))

GROUPS = (("genome-scale", ("GramNegModelTemplateV7.json", "GramPosModelTemplateV7.json",
                            "ArchaeaTemplateV6.json")), ("Core", ("Core-V6.json",)))


def cmd_table(a):
    """Per-reaction table of one tier: for each field and template group, what the tier
    changes (current -> target), what it leaves for another tier, or that it already
    matches. Genome-scale templates that disagree are listed separately."""
    rows = list(csv.DictReader(open(a.audit), delimiter="\t"))
    rxns = sorted({r["reaction"] for r in rows if r["tier"] == a.tier})
    present = {}
    for f in glob.glob(os.path.join(a.templates, "*.json")):
        present[os.path.basename(f)] = {r["id"] for r in json.load(open(f))["reactions"]}
    by = {(r["template"], r["reaction"], r["field"]): r for r in rows}
    out = []
    for rx in rxns:
        info = next(r for r in rows if r["reaction"] == rx and r["tier"] == a.tier)
        o = int(info["orientation"])
        dg = info["dG_template_orientation"]
        rec = {"reaction": rx.rsplit("_", 1)[0], "name": info["name"],
               "database": info["db_reversibility"],
               "dG_forward_kcal": "" if dg == "" else round(float(dg) * o, 2),
               "grade": info["thermo_grade"]}
        for gname, tmpls in GROUPS:
            for fld in ("direction", "GapfillDirection"):
                cells = {}
                for t in tmpls:
                    if rx not in present.get(t, ()):
                        continue
                    r = by.get((t, rx, fld))
                    if r is None:
                        c = "ok"
                    elif r["tier"] == a.tier:
                        c = "%s -> %s" % (r["current"], r["target"])
                    else:
                        c = "%s (tier %s)" % (r["current"], r["tier"])
                    cells.setdefault(c, []).append(t.split("Model")[0].split("Template")[0])
                rec["%s %s" % (gname, fld)] = ("-" if not cells else next(iter(cells)) if len(cells) == 1
                                               else "; ".join("%s: %s" % ("/".join(v), k) for k, v in cells.items()))
        out.append(rec)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print("%d reactions -> %s" % (len(out), a.out))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("audit")
    p.add_argument("--templates", required=True)
    p.add_argument("--db", default=None, help="ModelSEEDDatabase Biochemistry directory (default %s)" % MSD)
    p.add_argument("--out", required=True)
    p = sub.add_parser("apply")
    p.add_argument("--templates", required=True)
    p.add_argument("--audit", required=True)
    p.add_argument("--tiers", default="")
    p.add_argument("--negate-negative-reverse-penalty", action="store_true",
                   help="negate every negative reverse_penalty")
    p.add_argument("--negate-reverse-penalty-where", choices=["forward_favoured"], default=None,
                   help="negate a negative reverse_penalty only where the database dG favours "
                        "the forward direction")
    p.add_argument("--fields", default="direction,GapfillDirection")
    p.add_argument("--only", default=None)
    p.add_argument("--set", default=None,
                   help="explicit directions, e.g. rxn00006=>,rxn08167=> (both fields, all templates)")
    p.add_argument("--reactions", default=None,
                   help="comma-separated database ids (rxnNNNNN): restrict the tier set to these")
    p.add_argument("--dest", required=True)
    p = sub.add_parser("table")
    p.add_argument("--templates", required=True)
    p.add_argument("--audit", required=True)
    p.add_argument("--tier", default="B")
    p.add_argument("--out", required=True)
    a = ap.parse_args()
    {"audit": cmd_audit, "apply": cmd_apply, "table": cmd_table}[a.cmd](a)


if __name__ == "__main__":
    main()
