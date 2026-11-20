#!/usr/bin/env python
"""Run an unmodified prFBA script against a substitute PGFam cache (and E3 catalogue, and output
paths) without editing the repo: the control harness of the 2026-11 re-score (plan W2).

    ctl_run.py --common <scripts dir> --cache <snapshot> [--catalogue <json>]
               [--redirect TARGET=PATH ...] [--entry main] [--allow-cache-writes]
               -- <script.py> [script args ...]

What it does, in order:
  1. puts <scripts dir> first on sys.path (then the script's own directory), imports `_common`
     and sets `_common.PATHS.pgfam_cache = <snapshot>`: every tree script looks the cache up there
     at run time (score_genecap, prfba_selection, genome_function_sets, depth_probe_eval,
     fetch_bvbrc_annotations, pfunc_e3_panel_gate).
  2. rebinds the real cache path (absolute, and the relative `bvbrc_cache/genome_gene_families.json`
     that select_references uses as a default) and, with --catalogue, the real
     `bvbrc_cache/e3_clade_catalogue.json` to the substitutes wherever a repo module holds them as
     a module constant, a class attribute or a function default (pfunc_estimator binds
     DEFAULT_PGFAM_CACHE / DEFAULT_CATALOGUE as defaults at import, so changing the module
     attribute alone would not reach `load_clade_catalogue()` or `catalogue_is_current(cat)`).
     Modules imported later are rebound by an import hook the moment they load.
  3. installs an audit hook that makes any `open()` of the real cache (or of the real catalogue
     when --catalogue is given) raise PermissionError, in this process and in forked workers, so a
     control can never read the live data by a path this harness did not see. Unless
     --allow-cache-writes, opening the substitute cache (or a sibling temp file) for writing also
     raises: a lazy-miss write would change the snapshot that other controls share.
  4. --redirect TARGET=PATH, repeatable. TARGET is `PATHS.<attr>` (the `_common.PATHS` class),
     `<module>:<attr>` (an importable module's global; function defaults equal to the old value
     are rebound too), or `main:<NAME>` (a global of the script itself; needs --entry, because a
     script's module-level constants exist only once it has executed; a value starting with { or
     [ is parsed as JSON, e.g. main:REF_CELLS='{"exclude_0.987": "...", "exclude_0.995": "..."}').
     The E3 gate reads OUT_JSON, OUT_DIR and REF_CELLS as globals inside main(), so C-E3 is
     `--entry main --redirect main:OUT_JSON=... --redirect main:OUT_DIR=... --redirect main:REF_CELLS=...`.
  5. runs the script: by default `runpy.run_path(script, run_name="__main__")`. With --entry FUNC
     the script is executed as module `__ctl_main__` (kept in sys.modules, so multiprocessing can
     pickle its functions), the `main:` redirects are applied to its globals, and FUNC() is
     called with sys.argv = [script, args...].
  6. stamps the real cache (size, mtime) before and after and fails if it changed.

A JSON record of every rebinding is written to --record (default: stderr summary only).
"""
import argparse
import importlib.abc
import importlib.machinery
import json
import os
import runpy
import sys
import time
import types

REPO = "/home/freiburger/Documents/prFBA"
REAL_CACHE = os.path.join(REPO, "bvbrc_cache", "genome_gene_families.json")
REAL_CACHE_REL = "bvbrc_cache/genome_gene_families.json"
REAL_CATALOGUE = os.path.join(REPO, "bvbrc_cache", "e3_clade_catalogue.json")

RECORD = {"rebound": [], "redirects": [], "forbidden_opens": []}
_MAP = {}            # old str -> new str (exact match)
_ROOTS = [REPO]      # modules whose file is under one of these get rebound


def _say(msg):
    print(f"[ctl_run] {msg}", file=sys.stderr, flush=True)


def _rebind_callable(fn, where):
    n = 0
    f = getattr(fn, "__func__", fn)
    d = getattr(f, "__defaults__", None)
    if d and any(isinstance(x, str) and x in _MAP for x in d):
        f.__defaults__ = tuple(_MAP.get(x, x) if isinstance(x, str) else x for x in d)
        RECORD["rebound"].append({"where": where, "kind": "defaults"})
        n += 1
    kd = getattr(f, "__kwdefaults__", None)
    if kd and any(isinstance(v, str) and v in _MAP for v in kd.values()):
        f.__kwdefaults__ = {k: (_MAP.get(v, v) if isinstance(v, str) else v) for k, v in kd.items()}
        RECORD["rebound"].append({"where": where, "kind": "kwdefaults"})
        n += 1
    return n


def rebind_namespace(ns, modname):
    """Replace exact-match path strings in a module namespace: globals, function defaults, and
    attributes / method defaults of classes defined in that module."""
    if not _MAP:
        return 0
    n = 0
    for name, val in list(ns.items()):
        if isinstance(val, str) and val in _MAP:
            ns[name] = _MAP[val]
            RECORD["rebound"].append({"where": f"{modname}.{name}", "kind": "global",
                                      "old": val, "new": _MAP[val]})
            n += 1
        elif isinstance(val, types.FunctionType) and val.__module__ == modname:
            n += _rebind_callable(val, f"{modname}.{name}")
        elif isinstance(val, type) and getattr(val, "__module__", None) == modname:
            for an, av in list(vars(val).items()):
                if isinstance(av, str) and av in _MAP:
                    setattr(val, an, _MAP[av])
                    RECORD["rebound"].append({"where": f"{modname}.{name}.{an}", "kind": "class attr",
                                              "old": av, "new": _MAP[av]})
                    n += 1
                elif isinstance(av, (types.FunctionType, staticmethod, classmethod)):
                    n += _rebind_callable(av, f"{modname}.{name}.{an}")
    return n


def _under_roots(path):
    if not path:
        return False
    p = os.path.realpath(path)
    return any(p == r or p.startswith(r.rstrip("/") + "/") for r in _ROOTS)


class _RebindLoader(importlib.abc.Loader):
    def __init__(self, inner):
        self.inner = inner

    def create_module(self, spec):
        return self.inner.create_module(spec)

    def exec_module(self, module):
        self.inner.exec_module(module)
        rebind_namespace(vars(module), module.__name__)


class _RebindFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        spec = importlib.machinery.PathFinder.find_spec(name, path, target)
        if spec is not None and spec.loader is not None and _under_roots(spec.origin) \
                and hasattr(spec.loader, "exec_module"):
            spec.loader = _RebindLoader(spec.loader)
        return spec


def install_audit(forbid_read, forbid_write_prefixes):
    forbid_read = {os.path.realpath(p) for p in forbid_read}
    prefixes = [os.path.realpath(p) for p in forbid_write_prefixes]

    def _real(p):
        if not isinstance(p, (str, bytes, os.PathLike)):
            return None
        try:
            return os.path.realpath(os.fsdecode(p))
        except Exception:
            return None

    def _guarded(rp):
        return rp is not None and any(rp == q or rp.startswith(q + ".") or rp.startswith(q + "~")
                                      for q in prefixes)

    def hook(event, args):
        # a rename / replace / remove / truncate that would replace or alter the substitute cache
        # (the usual atomic-write pattern opens a temp file, then os.replace()s it over the target)
        if event in ("os.rename", "os.remove", "os.truncate", "os.link", "os.symlink") and args:
            targets = [args[1]] if event in ("os.rename", "os.link", "os.symlink") and len(args) > 1 else [args[0]]
            if event == "os.rename":
                targets.append(args[0])                 # moving the snapshot away is a write too
            for t in targets:
                rt = _real(t)
                if rt in forbid_read or (prefixes and rt in {os.path.realpath(q) for q in prefixes}):
                    RECORD["forbidden_opens"].append(f"{rt} ({event})")
                    raise PermissionError(f"ctl_run: {event} on {rt} is forbidden in a control run "
                                          f"(pass --allow-cache-writes to permit lazy-miss writes)")
            return
        if event != "open" or not args:
            return
        p = args[0]
        if not isinstance(p, (str, bytes, os.PathLike)):
            return
        try:
            rp = os.path.realpath(os.fsdecode(p))
        except Exception:
            return
        if rp in forbid_read:
            RECORD["forbidden_opens"].append(rp)
            raise PermissionError(f"ctl_run: the live file {rp} is forbidden in a control run")
        mode = args[1] if len(args) > 1 else "r"
        flags = args[2] if len(args) > 2 else 0
        writing = (isinstance(mode, str) and any(c in mode for c in "wax+")) or \
                  (isinstance(flags, int) and flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT))
        if writing and any(rp == q or rp.startswith(q + ".") or rp.startswith(q + "~") for q in prefixes):
            RECORD["forbidden_opens"].append(rp + " (write)")
            raise PermissionError(f"ctl_run: writing the substitute cache {rp} is forbidden "
                                  f"(pass --allow-cache-writes to permit lazy-miss writes)")
    sys.addaudithook(hook)


def _stat(path):
    st = os.stat(path)
    return {"size": st.st_size, "mtime": st.st_mtime}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if "--" not in argv:
        raise SystemExit("usage: ctl_run.py [options] -- script.py [args]")
    i = argv.index("--")
    own, rest = argv[:i], argv[i + 1:]
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--common", help="scripts dir whose _common.PATHS is patched")
    ap.add_argument("--cache", help="substitute PGFam cache (e.g. $SNAP)")
    ap.add_argument("--catalogue", help="substitute E3 clade catalogue (plain JSON)")
    ap.add_argument("--redirect", action="append", default=[], metavar="TARGET=PATH")
    ap.add_argument("--entry", help="call this function of the script instead of running it as __main__")
    ap.add_argument("--allow-cache-writes", action="store_true")
    ap.add_argument("--record", help="write a JSON record of the rebinding here")
    a = ap.parse_args(own)
    if not rest:
        ap.error("no script given after --")
    script = os.path.abspath(rest[0])
    sargs = rest[1:]
    t0 = time.time()

    if a.cache:
        a.cache = os.path.abspath(a.cache)
        if not os.path.exists(a.cache):
            raise SystemExit(f"substitute cache not found: {a.cache}")
        _MAP[REAL_CACHE] = a.cache
        _MAP[REAL_CACHE_REL] = a.cache
    if a.catalogue:
        a.catalogue = os.path.abspath(a.catalogue)
        if a.catalogue.endswith(".gz"):
            raise SystemExit("--catalogue must be plain JSON (gunzip it first)")
        if not os.path.exists(a.catalogue):
            raise SystemExit(f"substitute catalogue not found: {a.catalogue}")
        _MAP[REAL_CATALOGUE] = a.catalogue
    if a.common:
        _ROOTS.append(os.path.realpath(a.common))
    _ROOTS.append(os.path.realpath(os.path.dirname(script)))

    real_cache = str(REAL_CACHE)          # fixed now: the rebinding below must never reach it
    real_before = _stat(real_cache)
    forbid = ([REAL_CACHE] if a.cache else []) + ([REAL_CATALOGUE] if a.catalogue else [])
    install_audit(forbid, [] if (a.allow_cache_writes or not a.cache) else [a.cache])

    sys.meta_path.insert(0, _RebindFinder())
    # sys.path: the common dir first, then the script's dir (as `python script.py` would)
    sdir = os.path.dirname(script)
    for d in [sdir] + ([os.path.abspath(a.common)] if a.common else []):
        while d in sys.path:
            sys.path.remove(d)
        sys.path.insert(0, d)
    # modules already loaded from the repo (none normally) get rebound now
    me = os.path.realpath(__file__)
    for name, mod in list(sys.modules.items()):
        f = getattr(mod, "__file__", None)
        if _under_roots(f) and os.path.realpath(f) != me:      # never rebind this harness
            rebind_namespace(vars(mod), name)

    C = None
    if a.common:
        import importlib
        C = importlib.import_module("_common")
        if os.path.realpath(os.path.dirname(C.__file__)) != os.path.realpath(a.common):
            raise SystemExit(f"_common resolved to {C.__file__}, not under {a.common}")
        if a.cache:
            old = C.PATHS.pgfam_cache
            C.PATHS.pgfam_cache = a.cache
            RECORD["redirects"].append({"target": "PATHS.pgfam_cache", "old": old, "new": a.cache})

    main_redirects = []
    for r in a.redirect:
        tgt, _, path = r.partition("=")
        if not path:
            raise SystemExit(f"--redirect needs TARGET=PATH: {r}")
        if tgt.startswith("PATHS."):
            if C is None:
                raise SystemExit("PATHS.* redirects need --common")
            attr = tgt.split(".", 1)[1]
            if not hasattr(C.PATHS, attr):
                raise SystemExit(f"_common.PATHS has no attribute {attr}")
            old = getattr(C.PATHS, attr)
            setattr(C.PATHS, attr, path)
            RECORD["redirects"].append({"target": tgt, "old": old, "new": path})
        elif tgt.startswith("main:"):
            if not a.entry:
                raise SystemExit("main:NAME redirects need --entry")
            if path.lstrip().startswith(("{", "[")):          # a JSON value, e.g. a dict of paths
                path = json.loads(path)
            main_redirects.append((tgt[5:], path))
        elif ":" in tgt:
            modname, attr = tgt.split(":", 1)
            import importlib
            mod = importlib.import_module(modname)
            if not hasattr(mod, attr):
                raise SystemExit(f"{modname} has no attribute {attr}")
            old = getattr(mod, attr)
            setattr(mod, attr, path)
            if isinstance(old, str):
                _MAP[old] = path
                rebind_namespace(vars(mod), mod.__name__)
            RECORD["redirects"].append({"target": tgt, "old": old, "new": path})
        else:
            raise SystemExit(f"unknown redirect target {tgt!r}")

    _say(f"script {script} | cache {a.cache} | catalogue {a.catalogue} | "
         f"{len(RECORD['redirects'])} redirects | rebinding {len(_MAP)} path strings")
    sys.argv = [script] + sargs
    rc = 0
    try:
        if a.entry:
            code = compile(open(script, "rb").read(), script, "exec")
            mod = types.ModuleType("__ctl_main__")
            mod.__file__ = script
            sys.modules["__ctl_main__"] = mod
            exec(code, vars(mod))
            rebind_namespace(vars(mod), "__ctl_main__")
            for name, path in main_redirects:
                if name not in vars(mod):
                    raise SystemExit(f"the script has no global {name}")
                old = vars(mod)[name]
                vars(mod)[name] = path
                RECORD["redirects"].append({"target": f"main:{name}", "old": old, "new": path})
            fn = vars(mod).get(a.entry)
            if not callable(fn):
                raise SystemExit(f"the script has no callable {a.entry}")
            out = fn()
            rc = out if isinstance(out, int) else 0
        else:
            runpy.run_path(script, run_name="__main__")
    except SystemExit as e:
        rc = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
        if e.code not in (None, 0) and not isinstance(e.code, int):
            print(e.code, file=sys.stderr)
    finally:
        real_after = _stat(real_cache)
        RECORD.update(script=script, args=sargs, cache=a.cache, catalogue=a.catalogue,
                      real_cache_before=real_before, real_cache_after=real_after,
                      seconds=round(time.time() - t0, 1), rc=rc)
        if a.record:
            with open(a.record, "w") as fh:
                json.dump(RECORD, fh, indent=1, default=str)
        _say(f"done rc={rc} in {RECORD['seconds']}s; {len(RECORD['rebound'])} rebinds; "
             f"forbidden opens {len(RECORD['forbidden_opens'])}")
        if real_after != real_before:
            _say("ERROR: the real cache changed during the control run")
            rc = rc or 4
    return rc


if __name__ == "__main__":
    sys.exit(main())
