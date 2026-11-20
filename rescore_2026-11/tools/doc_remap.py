#!/usr/bin/env python
"""Propose old -> new replacements for the numbers of a hand-written results doc (re-score R4.3).

    doc_remap.py --doc DOC [--text OLD NEW ...] [--json OLD NEW ...] [--lines A-B] --out-json REPORT [--out-doc NEWDOC]

For every number in DOC (skipping digits glued to letters, e.g. E4, d050, f00):
  text sources: the token is looked up in OLD (a script-generated text whose NEW version has the same line and
    token structure); each occurrence maps to the NEW token at the same line/token position;
  json sources: leaves of OLD whose value, printed at the token's precision, equals the token map to the NEW
    leaf at the same path.
A token is 'auto' if every text occurrence (else every json match) maps to one new string; 'ambiguous' if they
disagree; 'unmatched' if none is found. --out-doc writes DOC with the auto replacements applied (style kept:
unicode minus, explicit '+', thousands separators). Nothing is written over DOC itself.
"""
import argparse
import json
import re

NUM = re.compile(r"(?<![A-Za-z0-9_.])([−+\-]?)(\d{1,3}(?:,\d{3})+|\d+)(\.\d+)?(?![A-Za-z0-9_])")


def norm(sign, intp, frac):
    s = "-" if sign in ("−", "-") else ""
    return s + intp.replace(",", "") + (frac or "")


def text_tokens(text):
    out = []
    for ln, line in enumerate(text.splitlines()):
        for m in NUM.finditer(line):
            out.append((ln, m.start(), norm(*m.groups())))
    return out


def text_map(old, new):
    """{normalised old token: [new tokens at the same position]} for lines with equal token counts."""
    ol, nl = old.splitlines(), new.splitlines()
    mp = {}
    note = ""
    if len(ol) != len(nl):
        note = f"line counts differ {len(ol)} vs {len(nl)}: aligned the first {min(len(ol), len(nl))}; "
    bad = 0
    for a, b in zip(ol, nl):
        ta = [norm(*m.groups()) for m in NUM.finditer(a)]
        tb = [norm(*m.groups()) for m in NUM.finditer(b)]
        if len(ta) != len(tb):
            bad += 1
            continue
        for x, y in zip(ta, tb):
            mp.setdefault(x, []).append(y)
    return mp, note + f"{bad} lines with differing token counts skipped"


def leaves(o, p=()):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from leaves(v, p + (str(k),))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from leaves(v, p + (i,))
    elif isinstance(o, (int, float)) and not isinstance(o, bool):
        yield p, o


def get(o, p):
    for k in p:
        try:
            o = o[k] if not isinstance(o, list) else o[int(k)]
        except (KeyError, IndexError, TypeError, ValueError):
            return None
    return o


def fmt_like(x, tok):
    dec = len(tok.split(".")[1]) if "." in tok else 0
    return f"{x:.{dec}f}"


def restyle(new_norm, old_raw_sign, old_intp):
    neg = new_norm.startswith("-")
    body = new_norm.lstrip("-")
    if "," in old_intp:
        ip, _, fp = body.partition(".")
        body = f"{int(ip):,}" + ("." + fp if fp else "")
    if neg:
        return ("−" if old_raw_sign == "−" or old_raw_sign == "" else "-") + body
    return ("+" if old_raw_sign == "+" else "") + body


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--doc", required=True)
    ap.add_argument("--text", nargs=2, action="append", default=[])
    ap.add_argument("--json", nargs=2, action="append", default=[])
    ap.add_argument("--lines", default=None)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-doc", default=None)
    a = ap.parse_args()
    tmaps, notes = [], []
    for o, n in a.text:
        mp, note = text_map(open(o).read(), open(n).read())
        tmaps.append(mp)
        notes.append(f"{o}: {note}")
    jsrc = []
    for o, n in a.json:
        jo, jn = json.load(open(o)), json.load(open(n))
        jsrc.append((o, list(leaves(jo)), jn))
    lo, hi = (0, 10 ** 9)
    if a.lines:
        lo, hi = (int(x) for x in a.lines.split("-"))
    doc = open(a.doc).read().splitlines(keepends=True)
    rep, out = [], []
    for ln, line in enumerate(doc, 1):
        if not (lo <= ln <= hi):
            out.append(line)
            continue
        pieces, last = [], 0
        for m in NUM.finditer(line):
            sign, intp, frac = m.groups()
            tok = norm(sign, intp, frac)
            cands = []
            for mp in tmaps:
                cands += mp.get(tok, [])
            src = "text"
            if not cands:
                src = "json"
                try:
                    val = float(tok)
                except ValueError:
                    val = None
                if val is not None and (frac or abs(val) >= 1):
                    for o, lv, jn in jsrc:
                        for p, x in lv:
                            if fmt_like(x, tok) == tok or (tok.startswith("-") is False and sign == "+" and fmt_like(x, tok) == tok):
                                y = get(jn, p)
                                if isinstance(y, (int, float)):
                                    cands.append(fmt_like(y, tok))
            uniq = sorted(set(cands))
            status = "unmatched" if not uniq else ("auto" if len(uniq) == 1 else "ambiguous")
            new = uniq[0] if status == "auto" else tok
            if status != "unmatched" and (uniq != [tok]):
                rep.append({"line": ln, "old": m.group(0), "status": status, "src": src,
                            "candidates": dict((u, cands.count(u)) for u in uniq),
                            "context": line.strip()[max(0, m.start() - 60):m.end() + 30]})
            if status == "auto" and new != tok:
                pieces.append(line[last:m.start()] + restyle(new, sign, intp))
                last = m.end()
        pieces.append(line[last:])
        out.append("".join(pieces))
    json.dump({"notes": notes, "changes": rep}, open(a.out_json, "w"), indent=1, ensure_ascii=False)
    if a.out_doc:
        open(a.out_doc, "w").write("".join(out))
    n_auto = sum(1 for r in rep if r["status"] == "auto")
    print(f"{len(rep)} tokens with a change or ambiguity; auto {n_auto}; ambiguous {len(rep) - n_auto}; notes {notes}")


if __name__ == "__main__":
    main()
