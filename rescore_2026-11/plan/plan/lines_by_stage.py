import json, re
d = json.load(open("claims_by_stage.json"))
def strip_paren(s):
    out, depth = [], 0
    for ch in s:
        if ch == "(": depth += 1; continue
        if ch == ")": depth = max(0, depth - 1); continue
        if depth == 0: out.append(ch)
    s = "".join(out)
    return re.sub(r"'[^']*'|\"[^\"]*\"", "", s)
def lines(s):
    s = strip_paren(s)
    ls = set()
    for a, b in re.findall(r"(?<![\d.])(\d{2,4})(?:\s*-\s*(\d{2,4}))?(?![\d.%])", s):
        a = int(a); b = int(b) if b else a
        if 50 <= a <= 1241 and a <= b <= 1241 and b - a < 80:
            ls.update(range(a, b + 1))
    return ls
def ranges(ls):
    ls = sorted(ls); out = []; i = 0
    while i < len(ls):
        j = i
        while j + 1 < len(ls) and ls[j + 1] == ls[j] + 1: j += 1
        out.append(f"{ls[i]}" if i == j else f"{ls[i]}-{ls[j]}"); i = j + 1
    return ", ".join(out)
res = {}
for st, items in d.items():
    L = set()
    for it in items: L |= lines(it["tex_lines"])
    res[st] = ranges(L)
json.dump(res, open("lines_by_stage.json", "w"), indent=1)
for st in ["R1.3U","R1.3L","R1.4","R1.5","R1.K","R2.1","R2.2","R2.3","R2.4","R2.5","R3","R3F","R4P","R4S2","R4H","R5A","D1","D7"]:
    print(st, "|", len(d[st]), "|", res[st])
