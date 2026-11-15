#!/usr/bin/env python
"""Number differ for the 2026-11 re-score (plan item W7).

For every manuscript claim that depends on the repaired BV-BRC data (the 178 dependent items of
the re-score plan's claims_by_stage.json), read each cited value from its artefact twice:

  old  = the artefact at the start tag (rescore-2026-11-start): tracked files from the tag;
         git-ignored files from $W/pre_repair/<repo-relative path> when the re-score moved them
         aside, otherwise from the live tree (flagged ``old_source: live-ignored``).
  new  = the artefact in the live tree (HEAD + whatever the re-score regenerated).

and print old -> new with the tex lines, so the manuscript update (R5.7) is mechanical.

The extractors live in manuscript_numbers/claims_spec.json, one entry per item key
("<inventory>:<key>"):

  {"stage", "inv", "key", "tex_lines", "claim",
   "values": [ VALUE, ... ],
   "manual": "<reason>"            # optional: the parts no extractor covers (or the whole item)
  }

  VALUE = {"label": "...", "printed": "<number as printed in the manuscript>", "fmt": FMT,
           "artifact": "<repo-relative path | absolute path | derived:<key> | derived>",
           "selector": [STEP, ...],                       # JSON artefacts
           "reduce": "min|max|minabs|maxabs|sum|mean|len",  # after an "each" step
           "recipe": "<name>", "args": {...},              # CSV/CSV.gz/text artefacts
           "expr": "a - b", "vars": {"a": VALUE, "b": VALUE}, # arithmetic over sub-values
           "note": "..."}

  STEP = "key" | 0 (list index) | {"where": {"field.sub": value, ...}} (unique list element)
         | {"each": ["k1", "k2"] | "*"} (map the rest over keys/elements; returns a list)
         | {"len": true}

  FMT  = f3 | +f3 | absf3 | pct1 | +pct1 | abspct1 | pp1 (already in percent) | int (thousands
         separators) | e1 (scientific, 1 decimal mantissa) | lt / le / gt (printed is a bound)
         | str

"derived:<key>" reads manuscript_numbers/derived_numbers.pre.json (old) and
derived_numbers.post.json (new), written by derived_numbers.py (W6) for the numbers no artefact
stores; a missing file or key leaves the value unresolved ("pending"), never an error.

Control (before anything is re-run): old == new everywhere, every non-manual selector resolves,
and the old value formatted with FMT equals the printed number. Mismatches are reported with a
reason (rounding / differs), never forced.

Usage:
  ~/Documents/py_venv/bin/python manuscript_numbers/collect.py
      [--spec manuscript_numbers/claims_spec.json] [--tag rescore-2026-11-start]
      [--w /home/freiburger/Documents/prFBA_rescore_2026-11] [--only KEY ...]
      [--out-json manuscript_numbers/number_diff.json] [--out-md manuscript_numbers/number_diff.md]
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import statistics
import subprocess
import sys
import time

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.join(R, "manuscript_numbers")
DEFAULT_TAG = "rescore-2026-11-start"
DEFAULT_W = "/home/freiburger/Documents/prFBA_rescore_2026-11"


class SelectorError(Exception):
    pass


class Pending(Exception):
    """A value whose source does not exist yet (derived numbers before W6 has run)."""


# ------------------------------------------------------------------ selector evaluation
def _get(obj, key):
    if isinstance(obj, list):
        try:
            return obj[int(key)]
        except (ValueError, IndexError) as e:
            raise SelectorError(f"list index {key!r}: {e}")
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        if not isinstance(key, str) and str(key) in obj:
            return obj[str(key)]
        raise SelectorError(f"missing key {key!r} (have {list(obj)[:8]}...)")
    raise SelectorError(f"cannot index {type(obj).__name__} with {key!r}")


def _dotted(obj, path):
    for part in str(path).split("."):
        obj = _get(obj, part)
    return obj


def _match(elem, cond):
    for k, v in cond.items():
        try:
            got = _dotted(elem, k)
        except SelectorError:
            return False
        if got != v and str(got) != str(v):
            return False
    return True


def select(obj, steps):
    """Apply selector steps; an {"each": ...} step maps the remaining steps and returns a list."""
    for i, st in enumerate(steps):
        if isinstance(st, dict):
            if "where" in st:
                if not isinstance(obj, list):
                    raise SelectorError("'where' on a non-list")
                hits = [e for e in obj if isinstance(e, dict) and _match(e, st["where"])]
                if st.get("pick") in ("first", "last") and hits:
                    hits = [hits[0] if st["pick"] == "first" else hits[-1]]
                if len(hits) != 1:
                    raise SelectorError(f"'where' {st['where']} matched {len(hits)} elements")
                obj = hits[0]
            elif "each" in st:
                keys = st["each"]
                if keys == "*":
                    keys = list(obj) if isinstance(obj, dict) else list(range(len(obj)))
                return [select(_get(obj, k), steps[i + 1:]) for k in keys]
            elif "len" in st:
                obj = len(obj)
            else:
                raise SelectorError(f"unknown step {st}")
        else:
            obj = _get(obj, st)
    return obj


def reduce_values(vals, how, skip_null=False):
    if how is None:
        return vals
    if how in ("maxkey", "minkey"):
        if not isinstance(vals, dict):
            raise SelectorError(f"reduce {how!r} needs a dict")
        keys = [float(k) for k in vals]
        return max(keys) if how == "maxkey" else min(keys)
    if not isinstance(vals, list):
        raise SelectorError(f"reduce {how!r} needs an 'each' step")
    if how == "len":
        return len(vals)
    if skip_null:
        vals = [v for v in vals if v is not None]
        if not vals:
            raise SelectorError("all values null")
    vals = [float(v) for v in vals]
    if how == "min":
        return min(vals)
    if how == "max":
        return max(vals)
    if how == "minabs":
        return min(vals, key=abs)
    if how == "maxabs":
        return max(vals, key=abs)
    if how == "sum":
        return sum(vals)
    if how == "mean":
        return sum(vals) / len(vals)
    if how == "len":
        return len(vals)
    raise SelectorError(f"unknown reduce {how!r}")


# ------------------------------------------------------------------ formatting
_SUP = str.maketrans({"−": "-", "–": "-", " ": ""})


def norm_printed(s):
    s = str(s).translate(_SUP)
    for a, b in (("{,}", ","), ("\\,", ""), ("\\%", "%"), ("$", ""), ("\\times", "x"),
                 ("{", ""), ("}", ""), (" ", "")):
        s = s.replace(a, b)
    return s


def _num(s):
    s = norm_printed(s).replace(",", "").replace("%", "")
    return float(s)


def _sci(s):
    """'2.3e-4', '2.3x10^-4', '1e-79' -> (mantissa, exponent)."""
    s = norm_printed(s).lower().replace("x10^", "e").replace("*10^", "e").replace("^", "")
    m = re.fullmatch(r"([+-]?\d+(?:\.\d+)?)e([+-]?\d+)", s)
    if not m:
        raise ValueError(f"not scientific: {s!r}")
    return float(m.group(1)), int(m.group(2))


def fmt_value(v, fmt):
    """Render v the way the manuscript prints it under FMT (bound formats return a repr)."""
    if v is None:
        return None
    if fmt == "str":
        return str(v)
    v = float(v)
    m = re.fullmatch(r"(\+|abs)?(f|pct|pp)(\d)", fmt)
    if m:
        sign, kind, nd = m.group(1), m.group(2), int(m.group(3))
        x = abs(v) if sign == "abs" else v
        if kind == "pct":
            x *= 100
        s = f"{x:+.{nd}f}" if sign == "+" else f"{x:.{nd}f}"
        if s in ("-0." + "0" * nd, "+0." + "0" * nd) and sign != "+":
            s = "0." + "0" * nd
        return s
    if fmt in ("int", "absint"):
        x = abs(v) if fmt == "absint" else v
        return f"{int(round(x)):,}"
    m = re.fullmatch(r"e(\d)", fmt)
    if m:
        mant, exp = f"{v:.{int(m.group(1))}e}".split("e")
        return f"{mant}e{int(exp)}"
    if fmt in ("lt", "le", "gt", "ge"):
        return repr(v)
    raise ValueError(f"unknown fmt {fmt!r}")


def _scaled(v, fmt):
    """The number fmt_value rounds (percent formats multiply by 100; abs formats drop the sign)."""
    x = float(v)
    if "abs" in fmt:
        x = abs(x)
    if "pct" in fmt:
        x *= 100
    return x


def matches_printed(v, printed, fmt):
    """(match: bool|None, reason: str|None)."""
    if v is None or printed is None:
        return None, None
    try:
        if fmt in ("lt", "le", "gt", "ge"):
            b = _sci(printed) if re.search(r"[eE^]", norm_printed(printed)) else None
            bound = b[0] * 10 ** b[1] if b else _num(printed)
            ok = {"lt": v < bound, "le": v <= bound, "gt": v > bound, "ge": v >= bound}[fmt]
            return ok, None if ok else "bound violated"
        if fmt == "str":
            ok = str(v) == str(printed)
            return ok, None if ok else "differs"
        if fmt.startswith("e"):
            mp, ep = _sci(printed)
            mv, ev = _sci(fmt_value(v, fmt))
            ok = (mp, ep) == (mv, ev)
            return ok, None if ok else ("rounding" if ep == ev and abs(mp - mv) <= 0.11 else "differs")
        got = norm_printed(fmt_value(v, fmt)).replace(",", "")
        want = norm_printed(printed).replace(",", "").replace("%", "")
        if not fmt.startswith("+"):
            want = want.lstrip("+")
        if got == want:
            return True, None
        try:
            nd = len(want.split(".")[1]) if "." in want else 0
            x = _scaled(v, fmt)
            # the stored value sits on the rounding boundary of the printed one (e.g. a 4-dp 0.3075
            # printed 0.308 from the unrounded 0.30751): consistent, flagged "boundary"
            if abs(x - float(want)) <= 0.5 * 10 ** (-nd) + 1e-9 * max(1.0, abs(x)):
                return True, "boundary"
            # near miss: within one unit of the last printed decimal -> rounding, else differs
            if abs(float(got) - float(want)) <= 1.0001 * 10 ** (-nd):
                return False, "rounding"
        except ValueError:
            pass
        return False, "differs"
    except (ValueError, TypeError) as e:
        return False, f"unparseable printed value ({e})"


# ------------------------------------------------------------------ sources (old / new)
class Sources:
    def __init__(self, tag, w, root=R, derived_dir=None):
        self.tag, self.w, self.root = tag, w, root
        self.derived_dir = derived_dir or os.path.join(root, "manuscript_numbers")
        self.tag_commit = self._git("rev-parse", f"{tag}^{{commit}}").strip()
        self.tracked = set(self._git("ls-tree", "-r", "--name-only", self.tag_commit).splitlines())
        self.start = os.path.join(w, "start")
        self.start_ok = False
        try:
            head = subprocess.run(["git", "-C", self.start, "rev-parse", "HEAD"], capture_output=True,
                                  text=True, check=True, env=_git_env()).stdout.strip()
            self.start_ok = head == self.tag_commit
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass
        self._cache = {}
        self.used = {}           # (version, path) -> {"source", "path", "sha256"}

    def _git(self, *a):
        return subprocess.run(["git", "-C", self.root, *a], capture_output=True, text=True, check=True,
                              env=_git_env()).stdout

    def locate(self, rel, version):
        """-> (filesystem path or None, source label, bytes or None)."""
        if os.path.isabs(rel):
            return (rel if os.path.exists(rel) else None), "external", None
        if version == "new":
            p = os.path.join(self.root, rel)
            return (p if os.path.exists(p) else None), "live", None
        if rel in self.tracked:
            if self.start_ok:
                return os.path.join(self.start, rel), "tag-worktree", None
            data = subprocess.run(["git", "-C", self.root, "show", f"{self.tag_commit}:{rel}"],
                                  capture_output=True, check=True, env=_git_env()).stdout
            return None, "tag-git-show", data
        p = os.path.join(self.w, "pre_repair", rel)
        if os.path.exists(p):
            return p, "pre_repair", None
        p = os.path.join(self.root, rel)
        return (p if os.path.exists(p) else None), "live-ignored", None

    def raw(self, rel, version):
        key = (rel, version)
        if key not in self._cache:
            path, src, data = self.locate(rel, version)
            if data is None:
                if path is None:
                    raise SelectorError(f"{version}: artefact not found: {rel}")
                with open(path, "rb") as fh:
                    data = fh.read()
            if rel.endswith(".gz"):
                data = gzip.decompress(data)
            self.used[key] = {"source": src, "path": path, "sha256": hashlib.sha256(data).hexdigest()}
            self._cache[key] = data
        return self._cache[key]

    def json(self, rel, version):
        k = ("json", rel, version)
        if k not in self._cache:
            self._cache[k] = json.loads(self.raw(rel, version))
        return self._cache[k]

    def text(self, rel, version):
        return self.raw(rel, version).decode("utf-8", errors="replace")

    def table(self, rel, version, usecols=None):
        """pandas DataFrame of a CSV(.gz); columns as str unless numeric-convertible."""
        k = ("csv", rel, version, tuple(usecols) if usecols else None)
        if k not in self._cache:
            import pandas as pd
            self._cache[k] = pd.read_csv(io.BytesIO(self.raw(rel, version)), usecols=usecols,
                                         low_memory=False, keep_default_na=True,
                                         dtype={"genome_id": str, "src_genome_id": str, "anchor": str})
        return self._cache[k]

    def derived(self, key, version, field="value", printed=None, tex_lines=None):
        name = "derived_numbers.pre.json" if version == "old" else "derived_numbers.post.json"
        p = os.path.join(self.derived_dir, name)
        if not os.path.exists(p):
            raise Pending(f"{name} not written yet (W6)")
        k = ("derived", version)
        if k not in self._cache:
            with open(p, "rb") as fh:
                data = fh.read()
            self.used[(f"manuscript_numbers/{name}", version)] = {
                "source": "derived", "path": p, "sha256": hashlib.sha256(data).hexdigest()}
            self._cache[k] = _index_derived(json.loads(data))
        entry = self._cache[k].get(key)
        if entry is None and printed is not None:
            entry = _derived_by_printed(self._cache[("derived", "old")] if ("derived", "old") in self._cache else self._cache[k],
                                        self._cache[k], printed, tex_lines)
        if entry is None:
            raise Pending(f"{name} has no key {key!r}")
        if not isinstance(entry, dict):
            return entry
        if field not in entry:
            raise Pending(f"{name}[{key!r}] has no field {field!r}")
        return entry[field]


def _derived_by_printed(old_index, index, printed, tex_lines):
    """Fallback when W6 named an entry differently: the unique W6 entry (located in the pre-state
    index, read in `index`) whose printed string equals ours, preferring one whose tex lines overlap."""
    want = norm_printed(printed).replace(",", "")
    hits = [k for k, e in old_index.items()
            if isinstance(e, dict) and norm_printed(str(e.get("printed", ""))).replace(",", "") == want]
    if len(hits) > 1 and tex_lines:
        lines = set(re.findall(r"\d+", str(tex_lines)))
        hits = [k for k in hits if lines & set(re.findall(r"\d+", str(old_index[k].get("tex_lines", ""))))] or hits
    if len(hits) != 1:
        return None
    e = index.get(hits[0])
    if isinstance(e, dict):
        e = dict(e, _matched_by="printed", _w6_key=hits[0])
    return e


def _index_derived(d):
    """Tolerant index of derived_numbers.*.json: {key: entry} whatever the container shape."""
    if isinstance(d, dict):
        for c in ("numbers", "items", "values", "results"):
            if c in d:
                return _index_derived(d[c])
        return {k: v for k, v in d.items() if not k.startswith("_")}
    if isinstance(d, list):
        return {e["key"]: e for e in d if isinstance(e, dict) and "key" in e}
    return {}


def _git_env():
    e = dict(os.environ)
    e["GIT_OPTIONAL_LOCKS"] = "0"
    return e


# ------------------------------------------------------------------ CSV / text recipes
def _filter(df, filt):
    """filt: {col: value | [op, value] | {"in": [...]}}; ops ==, !=, >, >=, <, <=, notna, isna."""
    import pandas as pd
    mask = pd.Series(True, index=df.index)
    for col, cond in (filt or {}).items():
        s = df[col]
        if isinstance(cond, dict) and "in" in cond:
            vals = cond["in"]
            mask &= s.isin(vals) | s.astype(str).isin([str(x) for x in vals])
        elif isinstance(cond, list):
            op, val = cond[0], (cond[1] if len(cond) > 1 else None)
            if op == "eqcol":
                mask &= s.astype(str) == df[val].astype(str)
            elif op == "notna":
                mask &= s.notna()
            elif op == "isna":
                mask &= s.isna()
            else:
                num = pd.to_numeric(s, errors="coerce")
                mask &= {"==": num == val, "!=": num != val, ">": num > val, ">=": num >= val,
                         "<": num < val, "<=": num <= val}[op]
        else:
            mask &= (s == cond) | (s.astype(str) == str(cond))
    return df[mask]


def recipe(src, spec, version):
    name, a = spec["recipe"], spec.get("args", {})
    art = spec["artifact"]
    if name == "regex":
        m = re.search(a["pattern"], src.text(art, version), re.S | re.M)
        if not m:
            raise SelectorError(f"regex not found: {a['pattern'][:60]}")
        return float(m.group(a.get("group", 1)).replace(",", ""))
    import pandas as pd
    cols = a.get("usecols")
    df = src.table(art, version, cols)
    for newcol, how in (a.get("derive") or {}).items():
        if how[0] == "split":
            df = df.assign(**{newcol: df[how[1]].astype(str).str.split(how[2]).str[int(how[3])]})
        else:
            raise SelectorError(f"unknown derive {how!r}")
    if a.get("dedup"):
        df = df.drop_duplicates(subset=a["dedup"])
    d = _filter(df, a.get("filter"))
    if a.get("within"):
        w = a["within"]
        keys = set(_filter(df, w["filter"])[w["on"]].astype(str))
        d = d[d[w["on"]].astype(str).isin(keys)]
    if name == "csv_count":
        return int(len(d))
    if name == "csv_nunique":
        return int(d[a["column"]].nunique())
    col = a.get("column")
    if name in ("csv_mean", "csv_median", "csv_sum"):
        x = pd.to_numeric(d[col], errors="coerce")
        if a.get("fillna") is not None:
            x = x.fillna(a["fillna"])
        else:
            x = x.dropna()
        if len(x) == 0:
            raise SelectorError("empty selection")
        return float({"csv_mean": x.mean, "csv_median": x.median, "csv_sum": x.sum}[name]())
    if name == "csv_share":
        if len(d) == 0:
            raise SelectorError("empty selection")
        return float(len(_filter(d, a["where"])) / len(d))
    raise SelectorError(f"unknown recipe {name!r}")


# ------------------------------------------------------------------ evaluation
_SAFE = {"abs": abs, "min": min, "max": max, "round": round, "sqrt": math.sqrt,
         "median": statistics.median, "len": len, "sum": sum}


def evaluate(src, spec, version):
    if "expr" in spec:
        env = {k: evaluate(src, v, version) for k, v in spec["vars"].items()}
        return eval(spec["expr"], {"__builtins__": {}}, {**_SAFE, **env})  # noqa: S307 (spec is ours)
    art = spec["artifact"]
    if version == "new" and spec.get("artifact_new") and os.path.exists(os.path.join(src.root, spec["artifact_new"])):
        spec = dict(spec, artifact=spec["artifact_new"])
        art = spec["artifact"]
    if art == "derived" or art.startswith("derived:"):
        key = spec.get("key") or art.split(":", 1)[1]
        val = src.derived(key, version, spec.get("field", "value"), printed=spec.get("printed"),
                          tex_lines=spec.get("_tex_lines"))
        if spec.get("selector"):
            val = select(val, spec["selector"])
        return reduce_values(val, spec.get("reduce"), spec.get("skip_null", False))
    if "recipe" in spec:
        return recipe(src, spec, version)
    val = select(src.json(art, version), spec.get("selector", []))
    return reduce_values(val, spec.get("reduce"), spec.get("skip_null", False))


def _artifacts(spec):
    if "expr" in spec:
        out = []
        for v in spec["vars"].values():
            out += _artifacts(v)
        return out
    return [spec.get("artifact"), spec.get("artifact_new")]


def _same(a, b):
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return a == b or (math.isfinite(a) and math.isfinite(b) and abs(a - b) <= 1e-12 * max(1, abs(a)))
    return a == b


def tex_line_set(tex_lines, pad=2):
    """'420-442; 1175 (57.6-59.2%)' -> {418..444, 1173..1177}; numbers in parentheses are ignored."""
    txt = re.sub(r"\([^)]*\)", "", str(tex_lines))
    txt = re.sub(r"'[^']*'", "", txt)
    out = set()
    for a, b in re.findall(r"(\d{1,4})(?:\s*[-\u2013]\s*(\d{1,4}))?", txt):
        lo, hi = int(a), int(b or a)
        if hi < lo or hi - lo > 200:
            hi = lo
        out.update(range(max(1, lo - pad), hi + pad + 1))
    return out


def printed_in_tex(tex, printed, fmt, lines, tex_text=None):
    """Is the printed number in the manuscript at (or near) the cited lines? None for non-numeric checks.
    tex_text: a literal to look for instead (numbers the text spells out: 'eight amplicons', 'a quarter')."""
    if tex is None or printed is None:
        return None
    if tex_text:
        pool = [tex[i - 1] for i in sorted(lines) if 0 < i <= len(tex)] if lines else tex
        return any(tex_text in ln for ln in pool)
    if fmt == "str":
        return None
    want = norm_printed(printed).lstrip("+").replace("%", "")
    if not re.search(r"\d", want):
        return None
    pool = [tex[i - 1] for i in sorted(lines) if 0 < i <= len(tex)] if lines else tex
    for ln in pool:
        t = ln.translate(_SUP)
        t = re.sub(r"(\d+(?:\.\d+)?)\s*\\times\s*10\^\{?\s*(-?\d+)\s*\}?", r"\1e\2", t)
        t = re.sub(r"(?<![\d.])10\^\{?\s*(-?\d+)\s*\}?", r"1e\1", t)
        for x, y in (("{,}", ","), ("\\,", ""), ("\\%", ""), ("$", ""), ("%", ""), ("{", " "), ("}", " ")):
            t = t.replace(x, y)
        for m in re.finditer(re.escape(want), t):
            before = t[m.start() - 1] if m.start() > 0 else ""
            after = t[m.end()] if m.end() < len(t) else ""
            if not (before.isdigit() or (before == "." and want[0].isdigit())) and not after.isdigit():
                return True
    return False


def run_item(src, item, tex=None):
    rows = []
    for v in item.get("values", []):
        v = dict(v, _tex_lines=item.get("tex_lines"))
        fmt = v.get("fmt", "f3")
        row = {"label": v.get("label"), "artifacts": sorted(set(filter(None, _artifacts(v)))),
               "selector": v.get("selector") or v.get("recipe") or v.get("expr"),
               "printed": v.get("printed"), "fmt": fmt, "note": v.get("note")}
        for ver in ("old", "new"):
            try:
                val = evaluate(src, v, ver)
                row[ver] = val
                row[f"{ver}_fmt"] = fmt_value(val, fmt) if not isinstance(val, list) else None
                row[f"{ver}_status"] = "ok"
            except Pending as e:
                row[ver], row[f"{ver}_fmt"], row[f"{ver}_status"] = None, None, f"pending: {e}"
            except (SelectorError, KeyError, ValueError, TypeError, ZeroDivisionError,
                    subprocess.CalledProcessError) as e:
                row[ver], row[f"{ver}_fmt"], row[f"{ver}_status"] = None, None, f"error: {e}"
        ok_old = row["old_status"] == "ok"
        ok_new = row["new_status"] == "ok"
        row["changed"] = (not _same(row["old"], row["new"])) if ok_old and ok_new else None
        row["changed_printed"] = (row["old_fmt"] != row["new_fmt"]) if ok_old and ok_new else None
        m, why = matches_printed(row["old"], row.get("printed"), fmt) if ok_old else (None, None)
        row["old_matches_printed"], row["mismatch"] = m, (why if m is False else None)
        row["match_note"] = why if m is True else None
        row["printed_in_tex"] = printed_in_tex(tex, v.get("printed"), fmt, tex_line_set(item.get("tex_lines")), v.get("tex_text"))
        if m is False and v.get("known_mismatch"):
            row["mismatch"] = f"{why}; known: {v['known_mismatch']}"
        rows.append(row)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--spec", default=os.path.join(HERE, "claims_spec.json"))
    ap.add_argument("--tag", default=DEFAULT_TAG)
    ap.add_argument("--w", default=DEFAULT_W)
    ap.add_argument("--only", nargs="*", help="item ids (<inv>:<key>) or bare keys")
    ap.add_argument("--out-json", default=os.path.join(HERE, "number_diff.json"))
    ap.add_argument("--out-md", default=os.path.join(HERE, "number_diff.md"))
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--root", default=R, help=argparse.SUPPRESS)
    ap.add_argument("--tex", default=os.path.expanduser("~/Documents/prFBA_paper/prFBA.tex"),
                    help="manuscript, read-only: checks that each printed number is at the cited lines")
    a = ap.parse_args(argv)

    t0 = time.time()
    spec = json.load(open(a.spec))
    items = spec["items"]
    if a.only:
        want = set(a.only)
        items = {k: v for k, v in items.items() if k in want or v["key"] in want}
    src = Sources(a.tag, a.w, root=a.root)
    tex = open(a.tex, encoding="utf-8").read().splitlines() if a.tex and os.path.exists(a.tex) else None
    out = []
    for iid, it in items.items():
        rows = run_item(src, it, tex)
        out.append({"id": iid, "stage": it["stage"], "inv": it["inv"], "key": it["key"],
                    "tex_lines": it["tex_lines"], "claim": it["claim"],
                    "extractor": ("manual" if not rows else "partial" if it.get("manual") else "spec"),
                    "manual": it.get("manual"), "values": rows})
        if not a.quiet:
            for r in rows:
                flag = ("CHANGED " if r["changed"] else "") + (
                    "" if r["old_matches_printed"] in (True, None) else f"[{r['mismatch']}] ")
                print(f"{it['stage']:6} {iid:55.55} L{str(it['tex_lines'])[:14]:14} "
                      f"{str(r['label'])[:28]:28} printed {str(r['printed']):>10} "
                      f"old {str(r['old_fmt']):>10} -> new {str(r['new_fmt']):>10} {flag}")
            if not rows:
                print(f"{it['stage']:6} {iid:55.55} MANUAL: {it.get('manual')}")

    vals = [r for o in out for r in o["values"]]
    summary = {
        "n_items": len(out),
        "n_items_spec": sum(o["extractor"] == "spec" for o in out),
        "n_items_partial": sum(o["extractor"] == "partial" for o in out),
        "n_items_manual": sum(o["extractor"] == "manual" for o in out),
        "n_values": len(vals),
        "n_old_resolved": sum(r["old_status"] == "ok" for r in vals),
        "n_new_resolved": sum(r["new_status"] == "ok" for r in vals),
        "n_pending": sum(r["old_status"].startswith("pending") or r["new_status"].startswith("pending")
                         for r in vals),
        "n_errors": sum(r["old_status"].startswith("error") or r["new_status"].startswith("error")
                        for r in vals),
        "n_changed": sum(bool(r["changed"]) for r in vals),
        "n_changed_at_print_precision": sum(bool(r["changed_printed"]) for r in vals),
        "n_old_matches_printed": sum(r["old_matches_printed"] is True for r in vals),
        "n_old_matches_printed_at_rounding_boundary": sum(r.get("match_note") == "boundary" for r in vals),
        "n_old_mismatch_printed": sum(r["old_matches_printed"] is False for r in vals),
        "mismatch_reasons": {},
        "n_printed_found_in_tex": sum(r.get("printed_in_tex") is True for r in vals),
        "n_printed_not_found_in_tex": sum(r.get("printed_in_tex") is False for r in vals),
    }
    for r in vals:
        if r["old_matches_printed"] is False:
            k = (r["mismatch"] or "").split(";")[0]
            summary["mismatch_reasons"][k] = summary["mismatch_reasons"].get(k, 0) + 1
    meta = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "tag": a.tag, "tag_commit": src.tag_commit,
            "head": subprocess.run(["git", "-C", a.root, "rev-parse", "HEAD"], capture_output=True, text=True,
                                   env=_git_env()).stdout.strip(),
            "w": a.w, "start_worktree_ok": src.start_ok, "spec": os.path.relpath(a.spec, a.root),
            "tex": a.tex, "tex_sha256": (hashlib.sha256(open(a.tex, "rb").read()).hexdigest() if tex is not None else None),
            "spec_sha256": hashlib.sha256(open(a.spec, "rb").read()).hexdigest(),
            "seconds": round(time.time() - t0, 1),
            "sources": {f"{ver}:{p}": u for (p, ver), u in sorted(src.used.items())},
            "ignored_read_live_as_old": sorted(p for (p, ver), u in src.used.items()
                                               if ver == "old" and u["source"] == "live-ignored")}
    with open(a.out_json + ".tmp", "w") as fh:
        json.dump({"_meta": meta, "summary": summary, "items": out}, fh, indent=1, default=str)
    os.replace(a.out_json + ".tmp", a.out_json)
    with open(a.out_md + ".tmp", "w") as fh:
        fh.write(render_md(meta, summary, out))
    os.replace(a.out_md + ".tmp", a.out_md)
    print(json.dumps(summary, indent=1))
    if meta["ignored_read_live_as_old"]:
        print("WARNING: git-ignored artefacts read from the live tree as 'old' (no $W/pre_repair copy): "
              + ", ".join(meta["ignored_read_live_as_old"]))
    return 0


def render_md(meta, summary, out):
    L = ["# Manuscript number diff (W7)", "",
         f"Generated {meta['generated']} by `manuscript_numbers/collect.py`; old = tag `{meta['tag']}` "
         f"(`{meta['tag_commit'][:7]}`; git-ignored files from `$W/pre_repair` or, failing that, the live "
         f"tree), new = live tree at `{meta['head'][:7]}`.", "",
         "| count | value |", "|---|---|"]
    for k, v in summary.items():
        if k != "mismatch_reasons":
            L.append(f"| {k} | {v} |")
    L.append(f"| mismatch_reasons | {summary['mismatch_reasons']} |")
    if meta["ignored_read_live_as_old"]:
        L += ["", "Git-ignored artefacts read live as 'old': " +
              ", ".join(f"`{p}`" for p in meta["ignored_read_live_as_old"])]
    stages = []
    for o in out:
        if o["stage"] not in stages:
            stages.append(o["stage"])
    for st in stages:
        L += ["", f"## {st}", "",
              "| item | tex lines | value | printed | old | new | flag |", "|---|---|---|---|---|---|---|"]
        rows = []
        for o in (x for x in out if x["stage"] == st):
            for r in o["values"]:
                flag = []
                if r["changed_printed"]:
                    flag.append("**CHANGED**")
                elif r["changed"]:
                    flag.append("changed below print precision")
                if r["old_matches_printed"] is False:
                    flag.append(f"old≠printed ({r['mismatch']})")
                if r.get("printed_in_tex") is False:
                    flag.append("printed not found at the cited tex lines")
                for ver in ("old", "new"):
                    if r[f"{ver}_status"] != "ok":
                        flag.append(f"{ver}: {r[f'{ver}_status'][:60]}")
                rows.append((0 if r["changed_printed"] else 1 if r["changed"] else 2,
                             f"| {o['key']} | {o['tex_lines']} | {r['label']} | {r['printed']} | "
                             f"{r['old_fmt']} | {r['new_fmt']} | {'; '.join(flag)} |"))
            if o["extractor"] != "spec":
                rows.append((3, f"| {o['key']} | {o['tex_lines']} | (manual) | | | | {o['manual']} |"))
        L += [r for _, r in sorted(rows, key=lambda x: x[0])]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    sys.exit(main())
