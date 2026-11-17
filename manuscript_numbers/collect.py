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
           "precise": VALUE,          # optional: the same quantity from a finer source (cells,
                                      # integer counts, a W6 recomputation); used for old/new when
                                      # it resolves, the stored value is kept as old/new_stored and
                                      # must agree with it at the stored precision
           "precise_tol": 1e-4,       # optional: allowed |precise - stored| (default half a unit
                                      # of the stored last decimal, at least 5e-5)
           "tex_exception": "<why>",  # optional: the current tex deliberately prints something
                                      # else here (a historical value, a number written as words)
           "check_new_print": true,   # optional, with "printed": null: a value first printed after
                                      # the base commit, checked against the current tex only
           "note": "..."}

  STEP = "key" | 0 (list index) | {"where": {"field.sub": value, ...}} (unique list element)
         | {"each": ["k1", "k2"] | "*"} (map the rest over keys/elements; returns a list)
         | {"len": true}

  FMT  = f3 | +f3 | absf3 | pct1 | +pct1 | abspct1 | pp1 (already in percent) | int (thousands
         separators) | e1 (scientific, 1 decimal mantissa) | lt / le / gt (printed is a bound)
         | str

  Rounding is the manuscript's: half away from zero on the decimal representation of the value
  (Decimal(repr(v)), ROUND_HALF_UP), not Python's binary-float formatting, which rounds a stored
  0.8085 down. A value stored at one decimal more than it prints and ending in 5 sits on the
  boundary: its print is decided only by a finer source (`precise`); without one, half-up is the
  paper's convention (bootstrap bounds stored at 4 dp), and the row is flagged
  `new_print_indeterminate`.

CSV recipes may derive columns before filtering: ["split", col, sep, i], ["ratio", num, den]
(NaN where den is 0) and ["f1_counts", tp, n_truth, n_pred] (2 tp / (n_truth + n_pred), 0 where
the denominator is 0), so a mean over cells is computed from their integer counts rather than
from the 4-dp per-cell scores.

Tex agreement (new, R5C): with --tex-base REV (the paper commit whose line numbers the spec's
tex_lines use, 0214a14), the cited lines are mapped onto the current manuscript (--tex) and each
value's new_fmt is looked for there (`new_fmt_in_tex`); a miss without a `tex_exception` is
listed in number_diff.md. The old printed value is looked for in the base manuscript.

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
import difflib
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
from decimal import Decimal, ROUND_HALF_UP

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


def _dec(v):
    """The value as the decimal it is written as (repr), so 0.8085 is 0.8085, not 0.80849999..."""
    return Decimal(repr(float(v)))


def _q(d, nd):
    """Round a Decimal half away from zero to nd decimals (the manuscript's rule)."""
    return d.quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP)


def fmt_value(v, fmt):
    """Render v the way the manuscript prints it under FMT (bound formats return a repr).
    Rounding is half away from zero on the decimal representation (see the module docstring)."""
    if v is None:
        return None
    if fmt == "str":
        return str(v)
    v = float(v)
    m = re.fullmatch(r"(\+|abs)?(f|pct|pp)(\d)", fmt)
    if m:
        sign, kind, nd = m.group(1), m.group(2), int(m.group(3))
        x = _dec(v)
        if sign == "abs":
            x = abs(x)
        if kind == "pct":
            x *= 100
        x = _q(x, nd)
        s = f"{x:+.{nd}f}" if sign == "+" else f"{x:.{nd}f}"
        if s in ("-0." + "0" * nd, "+0." + "0" * nd, "-0", "+0") and sign != "+":
            s = "0." + "0" * nd if nd else "0"
        return s
    if fmt in ("int", "absint"):
        x = _dec(v)
        if fmt == "absint":
            x = abs(x)
        return f"{int(_q(x, 0)):,}"
    m = re.fullmatch(r"e(\d)", fmt)
    if m:
        nd = int(m.group(1))
        if v == 0:
            return f"{0:.{nd}f}e0"
        d = _dec(v)
        exp = d.adjusted()
        mant = _q(d.scaleb(-exp), nd)
        if abs(mant) >= 10:                       # 9.96 -> 10.0: renormalise
            exp += 1
            mant = _q(d.scaleb(-exp), nd)
        return f"{mant:.{nd}f}e{exp}"
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


def print_boundary(v, fmt):
    """True when v, written at its own decimal precision, has exactly one decimal more than FMT
    prints and that decimal is 5: the print then depends on digits the stored value does not
    carry (a 4-dp 0.8085 printed at 3 dp)."""
    m = re.fullmatch(r"(\+|abs)?(f|pct|pp)(\d)", str(fmt))
    if v is None or isinstance(v, (bool, list, str)) or not m:
        return False
    try:
        d = _dec(v)
    except (TypeError, ValueError):
        return False
    if m.group(2) == "pct":
        d *= 100
    nd = int(m.group(3))
    t = d.normalize().as_tuple()
    return -t.exponent == nd + 1 and t.digits[-1] == 5


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
            # an exact match can still sit on a rounding boundary (a 4-dp 0.3075 printed 0.308):
            # the print is then the half-up convention's, not determined by the stored value
            return True, ("boundary" if print_boundary(v, fmt) else None)
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
        elif how[0] == "ratio":
            num, den = (pd.to_numeric(df[c], errors="coerce") for c in how[1:3])
            df = df.assign(**{newcol: num / den.where(den != 0)})
        elif how[0] == "f1_counts":
            tp, nt, npred = (pd.to_numeric(df[c], errors="coerce").fillna(0) for c in how[1:4])
            den = nt + npred
            df = df.assign(**{newcol: (2 * tp / den.where(den != 0)).fillna(0.0)})
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


def tex_line_map(base, cur):
    """{base line number: [current line numbers]} (1-based) from a line diff of two versions of
    the manuscript: an unchanged line maps to itself, a rewritten one to the block replacing it."""
    sm = difflib.SequenceMatcher(None, base, cur, autojunk=False)
    out = {}
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op == "equal":
            for k in range(i2 - i1):
                out[i1 + k + 1] = [j1 + k + 1]
        elif op in ("replace", "delete"):
            block = list(range(j1 + 1, j2 + 1)) or [min(j1 + 1, len(cur))]
            for k in range(i1, i2):
                out[k + 1] = block
    return out


def mapped_lines(tex_lines, lmap, pad=2):
    """The current-manuscript lines a claim's (base-numbered) tex_lines map to, padded."""
    if lmap is None:
        return tex_line_set(tex_lines, pad)
    out = set()
    for ln in tex_line_set(tex_lines, pad):
        for m in lmap.get(ln, []):
            out.update(range(max(1, m - pad), m + pad + 1))
    return out


def _stored_decimals(x):
    try:
        t = _dec(x).normalize().as_tuple()
    except (TypeError, ValueError):
        return None
    return max(0, -t.exponent)


def _eval_errors():
    return (SelectorError, KeyError, ValueError, TypeError, ZeroDivisionError,
            subprocess.CalledProcessError, IndexError)


def new_print_check(tex_cur, lmap, item, v, row, fmt):
    """Is the new value, as number_diff formats it, what the current manuscript prints at the
    claim's (mapped) lines?  -> {"new_fmt_in_tex": True | False | None, "new_fmt_tex_where": ...}.
    None: not checkable (no manuscript, a list, a string, or the new value did not resolve).
    Bound formats: the bound must still hold and the printed bound must be at the lines."""
    out = {"new_fmt_in_tex": None, "new_fmt_tex_where": None}
    if tex_cur is None or row.get("new_status") != "ok" or isinstance(row.get("new"), list):
        return out
    if fmt == "str" or v.get("tex_text"):
        return out
    if v.get("printed") is None and not v.get("check_new_print"):
        return out                     # a value the manuscript does not print
    lines = mapped_lines(item.get("tex_lines"), lmap)
    if fmt in ("lt", "le", "gt", "ge"):
        holds = matches_printed(row["new"], v.get("printed"), fmt)[0]
        found = printed_in_tex(tex_cur, v.get("printed"), fmt, lines)
        out["new_fmt_in_tex"] = bool(holds and found)
        out["new_fmt_tex_where"] = ("lines" if out["new_fmt_in_tex"] else
                                    "bound no longer holds" if not holds else "bound not at lines")
    else:
        nf = row.get("new_fmt")
        if nf is None:
            return out
        if printed_in_tex(tex_cur, nf, fmt, lines):
            out["new_fmt_in_tex"], out["new_fmt_tex_where"] = True, "lines"
        else:
            out["new_fmt_in_tex"] = False
            out["new_fmt_tex_where"] = ("elsewhere in the file" if printed_in_tex(tex_cur, nf, fmt, None)
                                        else "absent")
    if v.get("tex_exception"):
        out["tex_exception"] = v["tex_exception"]
    return out


def run_item(src, item, tex=None, tex_cur=None, lmap=None):
    rows = []
    for v in item.get("values", []):
        v = dict(v, _tex_lines=item.get("tex_lines"))
        fmt = v.get("fmt", "f3")
        row = {"label": v.get("label"), "artifacts": sorted(set(filter(None, _artifacts(v)))),
               "selector": v.get("selector") or v.get("recipe") or v.get("expr"),
               "printed": v.get("printed"), "fmt": fmt, "note": v.get("note")}
        if v.get("precise"):
            row["precise_artifacts"] = sorted(set(filter(None, _artifacts(v["precise"]))))
        for ver in ("old", "new"):
            try:
                val = evaluate(src, v, ver)
                row[f"{ver}_status"] = "ok"
                if v.get("precise") and not isinstance(val, list):
                    # the same quantity from a finer source decides the print; the stored value
                    # must agree with it at the stored precision
                    try:
                        pv = evaluate(src, dict(v["precise"], _tex_lines=item.get("tex_lines")), ver)
                        row[f"{ver}_stored"] = val
                        row[f"{ver}_source"] = "precise"
                        k = _stored_decimals(val)
                        tol = max(0.5 * 10 ** -k, 5e-5) if k is not None else 5e-5
                        tol = v.get("precise_tol", tol)   # e.g. a stored difference of two rounded values
                        if isinstance(val, (int, float)) and abs(float(pv) - float(val)) > tol + 1e-12:
                            row[f"{ver}_precise_disagrees"] = f"{pv!r} vs stored {val!r}"
                        val = pv
                    except Pending as e:
                        row[f"{ver}_source"] = f"stored (precise pending: {e})"
                    except _eval_errors() as e:
                        row[f"{ver}_source"] = f"stored (precise unavailable: {e})"
                row[ver] = val
                row[f"{ver}_fmt"] = fmt_value(val, fmt) if not isinstance(val, list) else None
                row[f"{ver}_print_indeterminate"] = print_boundary(val, fmt)
            except Pending as e:
                row[ver], row[f"{ver}_fmt"], row[f"{ver}_status"] = None, None, f"pending: {e}"
            except _eval_errors() as e:
                row[ver], row[f"{ver}_fmt"], row[f"{ver}_status"] = None, None, f"error: {e}"
        ok_old = row["old_status"] == "ok"
        ok_new = row["new_status"] == "ok"
        row["changed"] = (not _same(row["old"], row["new"])) if ok_old and ok_new else None
        row["changed_printed"] = (row["old_fmt"] != row["new_fmt"]) if ok_old and ok_new else None
        m, why = matches_printed(row["old"], row.get("printed"), fmt) if ok_old else (None, None)
        row["old_matches_printed"], row["mismatch"] = m, (why if m is False else None)
        row["match_note"] = why if m is True else None
        row["printed_in_tex"] = printed_in_tex(tex, v.get("printed"), fmt, tex_line_set(item.get("tex_lines")), v.get("tex_text"))
        row.update(new_print_check(tex_cur, lmap, item, v, row, fmt))
        if m is False and v.get("known_mismatch"):
            row["mismatch"] = f"{why}; known: {v['known_mismatch']}"
        row["meaning_flags"] = meaning_flags(row) if ok_old and ok_new else []
        rows.append(row)
    _ci_flags(rows)
    return rows


_PQ = re.compile(r"(^|[^a-z])(p|q|p-value|q_bh|wilcoxon p|mcnemar p)([^a-z]|$)", re.I)


def meaning_flags(row):
    """Automatic flags for a change that alters what a sentence says, not just its digits."""
    old, new, fmt, lab = row["old"], row["new"], row["fmt"], str(row.get("label") or "")
    out = []
    if isinstance(old, list) or isinstance(new, list):
        return out
    if fmt == "str" or isinstance(old, bool) or isinstance(new, bool):
        if str(old) != str(new):
            out.append(f"categorical change {old} -> {new}")
        return out
    try:
        o, n = float(old), float(new)
    except (TypeError, ValueError):
        return out
    if fmt in ("lt", "le", "gt", "ge"):
        po = matches_printed(o, row.get("printed"), fmt)[0]
        pn = matches_printed(n, row.get("printed"), fmt)[0]
        if po and pn is False:
            out.append(f"printed bound {fmt} {row.get('printed')} no longer holds")
        return out
    if o * n < 0 and "abs" not in fmt:
        out.append("sign flip")
    if _PQ.search(lab) and "spearman rho" not in lab.lower():
        for thr in (0.05, 0.01):
            if (o <= thr) != (n <= thr):
                out.append(f"significance crosses {thr} ({o:.3g} -> {n:.3g})")
                break
    return out


def _ci_flags(rows):
    """Pair consecutive 'CI ...' rows of one item (lower, upper) and flag a flip of 'excludes zero'."""
    ci = [r for r in rows if "ci" in str(r.get("label") or "").lower().split()
          or str(r.get("label") or "").strip().lower().startswith("ci")]
    for a, b in zip(ci[0::2], ci[1::2]):
        try:
            lo_o, hi_o, lo_n, hi_n = (float(a["old"]), float(b["old"]), float(a["new"]), float(b["new"]))
        except (TypeError, ValueError):
            continue
        if "abs" in (a["fmt"] + b["fmt"]):
            continue          # absolute-value bounds: the sign is printed elsewhere
        ex_o = lo_o > 0 or hi_o < 0 if lo_o <= hi_o else lo_o < 0 or hi_o > 0
        ex_n = lo_n > 0 or hi_n < 0 if lo_n <= hi_n else lo_n < 0 or hi_n > 0
        if ex_o != ex_n:
            msg = "CI excludes zero: %s -> %s" % (ex_o, ex_n)
            for r in (a, b):
                r["meaning_flags"] = r.get("meaning_flags", []) + [msg]


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
                    help="manuscript, read-only: checks that each new value prints there")
    ap.add_argument("--tex-base", default="0214a14",
                    help="paper commit the spec's tex_lines and printed values refer to (read with git "
                         "show from the --tex file's repository); '' = use --tex for both checks")
    a = ap.parse_args(argv)

    t0 = time.time()
    spec = json.load(open(a.spec))
    items = spec["items"]
    if a.only:
        want = set(a.only)
        items = {k: v for k, v in items.items() if k in want or v["key"] in want}
    src = Sources(a.tag, a.w, root=a.root)
    tex_cur = open(a.tex, encoding="utf-8").read().splitlines() if a.tex and os.path.exists(a.tex) else None
    tex, lmap, base_sha = tex_cur, None, None
    if tex_cur is not None and a.tex_base:
        try:
            blob = subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(a.tex)), "show",
                                   f"{a.tex_base}:{os.path.basename(a.tex)}"], capture_output=True,
                                  check=True, env=_git_env()).stdout
            tex = blob.decode("utf-8").splitlines()
            base_sha = hashlib.sha256(blob).hexdigest()
            lmap = tex_line_map(tex, tex_cur)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"WARNING: cannot read {a.tex_base}:{os.path.basename(a.tex)}; the old-print and new-print "
                  "checks both use --tex")
    out = []
    for iid, it in items.items():
        rows = run_item(src, it, tex, tex_cur, lmap)
        mvals = it.get("manual_values") or []
        auto = sorted({f for r in rows for f in r.get("meaning_flags", [])})
        out.append({"id": iid, "stage": it["stage"], "inv": it["inv"], "key": it["key"],
                    "tex_lines": it["tex_lines"], "claim": it["claim"],
                    "extractor": ("manual" if not rows and not mvals else
                                  "partial" if (it.get("manual") or mvals) else "spec"),
                    "manual": it.get("manual"), "manual_resolution": it.get("manual_resolution"),
                    "manual_values": mvals, "meaning_change": it.get("meaning_change"),
                    "meaning_flags_auto": auto,
                    "changed_at_print": (any(bool(r["changed_printed"]) for r in rows)
                                         or any(str(m.get("new")) not in (str(m.get("printed")), "unchanged at print",
                                                                          "as printed")
                                                and not str(m.get("new")).startswith("see ") for m in mvals)),
                    "values": rows})
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
        "n_items_changed_at_print": sum(bool(o["changed_at_print"]) for o in out),
        "n_items_meaning_change_curated": sum(bool(o["meaning_change"]) for o in out),
        "n_items_meaning_flag_auto": sum(bool(o["meaning_flags_auto"]) for o in out),
        "n_items_meaning_change_any": sum(bool(o["meaning_change"] or o["meaning_flags_auto"]) for o in out),
        "n_items_with_manual_values": sum(bool(o["manual_values"]) for o in out),
        "n_items_manual_unresolved": sum(bool(o["manual"]) and not o["manual_resolution"] for o in out),
        "manual_resolution_status": {},
        "n_new_from_precise_source": sum(r.get("new_source") == "precise" for r in vals),
        "n_new_precise_disagrees_with_stored": sum(bool(r.get("new_precise_disagrees")) for r in vals),
        "n_new_print_indeterminate": sum(bool(r.get("new_print_indeterminate")) for r in vals),
        "n_new_fmt_in_tex": sum(r.get("new_fmt_in_tex") is True for r in vals),
        "n_new_fmt_not_in_tex": sum(r.get("new_fmt_in_tex") is False for r in vals),
        "n_new_fmt_not_in_tex_with_exception": sum(r.get("new_fmt_in_tex") is False and bool(r.get("tex_exception"))
                                                   for r in vals),
        "n_new_fmt_not_in_tex_unexplained": sum(r.get("new_fmt_in_tex") is False and not r.get("tex_exception")
                                                for r in vals),
    }
    for o in out:
        st = (o["manual_resolution"] or {}).get("status")
        if st:
            summary["manual_resolution_status"][st] = summary["manual_resolution_status"].get(st, 0) + 1
    for r in vals:
        if r["old_matches_printed"] is False:
            k = (r["mismatch"] or "").split(";")[0]
            summary["mismatch_reasons"][k] = summary["mismatch_reasons"].get(k, 0) + 1
    meta = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "tag": a.tag, "tag_commit": src.tag_commit,
            "head": subprocess.run(["git", "-C", a.root, "rev-parse", "HEAD"], capture_output=True, text=True,
                                   env=_git_env()).stdout.strip(),
            "w": a.w, "start_worktree_ok": src.start_ok, "spec": os.path.relpath(a.spec, a.root),
            "tex": a.tex, "tex_sha256": (hashlib.sha256(open(a.tex, "rb").read()).hexdigest() if tex_cur is not None else None),
            "tex_base": a.tex_base or None, "tex_base_sha256": base_sha,
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
    mc = [o for o in out if o["meaning_change"] or o["meaning_flags_auto"]]
    L += ["", "## Claims whose meaning changes (rewrite, not just refresh)", "",
          "Curated notes (`meaning_change` in claims_spec.json) and automatic flags (sign flip, a p or q crossing "
          "0.05 / 0.01, a categorical or yes/no value changing, a CI's 'excludes zero' flipping, a printed bound no "
          "longer holding). Automatic flags are candidates: read the row before rewriting.", "",
          "| item | tex lines | curated note | automatic flags (value: flag) |", "|---|---|---|---|"]
    for o in mc:
        auto = "; ".join(f"{r['label']}: {', '.join(r['meaning_flags'])}" for r in o["values"] if r.get("meaning_flags"))
        L.append(f"| {o['id']} | {o['tex_lines']} | {o['meaning_change'] or ''} | {auto} |")
    miss = [(o, r) for o in out for r in o["values"] if r.get("new_fmt_in_tex") is False]
    L += ["", "## New values against the current manuscript", "",
          f"Each new value, formatted as the manuscript prints it, looked for at the claim's lines mapped from "
          f"`{meta.get('tex_base')}` onto the current tex (\u00b12 lines). "
          f"{summary.get('n_new_fmt_in_tex', 0)} found; {summary.get('n_new_fmt_not_in_tex', 0)} not found, "
          f"{summary.get('n_new_fmt_not_in_tex_with_exception', 0)} of them with a recorded reason "
          f"(`tex_exception`). Values from a finer source than the stored one: "
          f"{summary.get('n_new_from_precise_source', 0)}; still on a rounding boundary (printed by the "
          f"half-up convention): {summary.get('n_new_print_indeterminate', 0)}.", "",
          "| item | tex lines | value | new | where | reason |", "|---|---|---|---|---|---|"]
    for o, r in miss:
        L.append(f"| {o['key']} | {o['tex_lines']} | {r['label']} | {r['new_fmt']} | {r.get('new_fmt_tex_where')} | "
                 f"{r.get('tex_exception') or '**unexplained**'} |")
    ind = [(o, r) for o in out for r in o["values"] if r.get("new_print_indeterminate")]
    if ind:
        L += ["", "Values still on a rounding boundary (no finer source; printed half-up):", "",
              "| item | value | new | printed as |", "|---|---|---|---|"]
        for o, r in ind:
            L.append(f"| {o['key']} | {r['label']} | {r['new']} | {r['new_fmt']} |")
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
                if r.get("meaning_flags"):
                    flag.append("MEANING: " + ", ".join(r["meaning_flags"]))
                for ver in ("old", "new"):
                    if r[f"{ver}_status"] != "ok":
                        flag.append(f"{ver}: {r[f'{ver}_status'][:60]}")
                rows.append((0 if r["changed_printed"] else 1 if r["changed"] else 2,
                             f"| {o['key']} | {o['tex_lines']} | {r['label']} | {r['printed']} | "
                             f"{r['old_fmt']} | {r['new_fmt']} | {'; '.join(flag)} |"))
            for m in o["manual_values"]:
                rows.append((1, f"| {o['key']} | {o['tex_lines']} | {m['label']} | {m['printed']} | "
                                f"{m['printed']} | {m['new']} | manual value; source: {m['source']} |"))
            if o["manual"] or o["manual_resolution"]:
                res = o["manual_resolution"] or {}
                txt = (f"**{res.get('status')}**" + (f" ({res['parts']})" if res.get("parts") else "")
                       + f": {res.get('why')}") if res else f"UNRESOLVED: {o['manual']}"
                rows.append((3, f"| {o['key']} | {o['tex_lines']} | (manual) | | | | {txt} |"))
        L += [r for _, r in sorted(rows, key=lambda x: x[0])]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    sys.exit(main())
