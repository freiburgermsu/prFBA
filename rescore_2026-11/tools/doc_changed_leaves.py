#!/usr/bin/env python
"""R5.4 helper: which numbers of a doc can have moved. For each number token in DOC, list the JSON leaves
(old/new file pairs) whose OLD value prints as the token at the token's precision AND whose value changed
between old and new. A token no changed leaf matches is reported as 'unchanged' (no changed source can
print it); a token some changed leaf matches is a candidate to review (with the matching paths, changed and
unchanged). usage: doc_changed_leaves.py DOC OLD1 NEW1 [OLD2 NEW2 ...] > report"""
import json, re, sys
NUM = re.compile(r"(?<![A-Za-z0-9_.])([−+\-]?)(\d{1,3}(?:,\d{3})+|\d+)(\.\d+)?(?![A-Za-z0-9_])")


def leaves(o, p=()):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from leaves(v, p + (str(k),))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from leaves(v, p + (i,))
    elif isinstance(o, (int, float)) and not isinstance(o, bool):
        yield p, float(o)


def get(o, p):
    for k in p:
        try:
            o = o[k] if not isinstance(o, list) else o[int(k)]
        except (KeyError, IndexError, TypeError, ValueError):
            return None
    return o


def fmts(v, nd, pct):
    x = v * 100 if pct else v
    return ["%.*f" % (nd, x)]


doc = open(sys.argv[1]).read().splitlines()
pairs = [(json.load(open(a)), json.load(open(b)), a) for a, b in zip(sys.argv[2::2], sys.argv[3::2])]
L = []
for o, n, name in pairs:
    for p, v in leaves(o):
        nv = get(n, p)
        L.append((name, p, v, nv if isinstance(nv, (int, float)) and not isinstance(nv, bool) else None))
for ln, line in enumerate(doc, 1):
    for m in NUM.finditer(line):
        sign, ip, fr = m.groups()
        tok = ("-" if sign in "−-" and sign else "") + ip.replace(",", "") + (fr or "")
        nd = len(fr) - 1 if fr else 0
        pct = line[m.end():m.end() + 1] == "%"
        t = float(tok)
        hits_ch, hits_un = [], []
        for name, p, v, nv in L:
            for x, scale in ((v, 1), (abs(v), 1)):
                if "%.*f" % (nd, x * (100 if pct else 1)) in (tok, tok.lstrip("-")):
                    if nv is None:
                        continue
                    newx = nv if x == v else abs(nv)
                    newf = "%.*f" % (nd, newx * (100 if pct else 1))
                    (hits_ch if newf != "%.*f" % (nd, x * (100 if pct else 1)) else hits_un).append(
                        (name.split("/")[-1], "/".join(map(str, p)), newf))
                    break
        if hits_ch:
            print("L%d %s%s | ctx: %s" % (ln, m.group(0), "%" if pct else "", line[max(0, m.start() - 60):m.end() + 20]))
            for h in hits_ch[:6]:
                print("    CHANGED-LEAF %s %s -> %s" % (h[0], h[1], h[2]))
            for h in hits_un[:3]:
                print("    unchanged-leaf %s %s" % (h[0], h[1]))
            if len(hits_un) > 3:
                print("    ... %d more unchanged leaves" % (len(hits_un) - 3))
