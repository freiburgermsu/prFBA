#!/usr/bin/env python
"""repro.py: the failure points of ModelSEEDpy MSProbability.prFBA, reproduced.

Run it once per ModelSEEDpy tree, from this folder, with the only interpreter allowed here:

  PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-c2e31a0 \
      ~/Documents/py_venv/bin/python repro.py > repro_c2e31a0.log 2>&1
  PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-prfba \
      ~/Documents/py_venv/bin/python repro.py > repro_b6857a0.log 2>&1

The header of each log proves which tree was imported (modelseedpy.__file__, the source file of
MSProbability.prFBA, its md5, and the tree's HEAD). Case ids match INVESTIGATION.md.

Fixtures (this folder)
  fixture_g100053.8_gf.json.gz     A gapfilled single-genome ModelSEED model: genome 100053.8,
                                   gramneg template, 515 ECs, 1261 reactions, 61 of them added by
                                   gapfilling. Built by phenotype_fold.gapfill_block (as_shipped
                                   pipeline, complete medium, ModelSEEDpy c2e31a0) for the weighted-
                                   objective probe (WEIGHTED_OBJECTIVE_TEST.md section 2). As built
                                   it carries ATPM_c0 as its objective (MSATPCorrection.build_tests).
  fixture_g100053.8_gf.meta.json   The gapfill block: the 61 gapfilled reaction ids (new_rxn).
  fixture_kbase_complete_media.txt The 347 compound ids of ModelSEEDDatabase Media/KBaseMedia.cpd;
                                   the environment is {cpd: 100.0}, as in the probe.
  A two-member community and e_coli_core (cobra's bundled "textbook", BiGG ids) are built in code.

Probabilities follow the production rule on a single genome: p = 1.0 for a gene-associated "rxn"
reaction, 0.0 for a gene-less one, stored as str in rxn.notes["probability"].

Parts
  A  The released call forms (both trees): bare Model, one-element list, environment=None, two models.
  B  c2e31a0 only. Defects hidden behind an earlier crash, reached by line-count-preserving edits
     of the released module source (F1: unwrap one model on l.165; F2: drop the duplicate
     ElementUptakePkg build on l.173-174). The edited module is compiled under a pseudo file name,
     so every traceback line number is the released file's own. The edits exist only to get past
     earlier defects; nothing else changes. On b6857a0 the edit targets are absent and Part B
     reports n/a.
  C  b6857a0's keywords and defaults (both trees; on c2e31a0 they fail at the first defect or with
     "unexpected keyword argument").
  D  Further problems (both trees; c2e31a0 through the F1+F2 variant, b6857a0 through the real
     function; the numbers must agree).

Nothing here writes outside this folder except prFBA.lp in a temporary directory (case C5). No
ModelSEEDpy tree is modified; the only git commands are `git rev-parse` and `git show`.
"""
import contextlib
import gzip
import hashlib
import inspect
import io
import json
import linecache
import math
import os
import subprocess
import sys
import tempfile
import time
import traceback
import types

sys.stdout.reconfigure(line_buffering=True)        # keep library warnings (stderr) in order in the log

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE = os.path.join(HERE, "fixture_g100053.8_gf.json.gz")
FIXTURE_META = os.path.join(HERE, "fixture_g100053.8_gf.meta.json")
MEDIA = os.path.join(HERE, "fixture_kbase_complete_media.txt")
KNOWN_MD5 = {"8d40e57a9978b25f816bf7193448a88a": "c2e31a0 (released; byte-identical at e1254ba)",
             "1014bf2e33472be88efdc93b60867677": "b6857a0 (fix-prfba-single-model)"}
SUMMARY = []          # (case id, one-line outcome) for the table at the end


# =========================================================================== environment
def git(tree, *args):
    try:
        return subprocess.run(["git", "-C", tree, *args], capture_output=True, text=True,
                              timeout=60).stdout
    except Exception as exc:                       # noqa: BLE001
        return "<git failed: {}>".format(exc)


def md5_bytes(b):
    return hashlib.md5(b).hexdigest()


import modelseedpy                                   # noqa: E402  (the PYTHONPATH tree)
import cobra                                         # noqa: E402
import cobra.io                                      # noqa: E402
import optlang                                       # noqa: E402

with contextlib.redirect_stdout(io.StringIO()) as _banner_mscomm:
    import mscommunity                               # noqa: E402  (prints package versions)
    import mscommunity.commhelper                    # noqa: E402,F401
    import mscommunity.mscommsim as MCS              # noqa: E402
    from modelseedpy.core.msprobability import MSProbability   # noqa: E402
from modelseedpy.core.msmodelutl import MSModelUtil              # noqa: E402
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager  # noqa: E402

MSPROB_FILE = inspect.getsourcefile(MSProbability.prFBA)
TREE = os.path.dirname(os.path.dirname(os.path.abspath(modelseedpy.__file__)))
MSPROB_SRC = open(MSPROB_FILE).read()
MSPROB_MD5 = md5_bytes(MSPROB_SRC.encode())
HEAD = git(TREE, "rev-parse", "HEAD").strip()
HEAD_SRC = git(TREE, "show", "HEAD:modelseedpy/core/msprobability.py")
IS_FIXED = "missing_probability" in inspect.signature(MSProbability.prFBA).parameters


def header():
    print("=" * 100)
    print("repro.py: MSProbability.prFBA failure points    run at", time.strftime("%Y-%m-%d %H:%M:%S %z"))
    print("=" * 100)
    print("python                :", sys.executable, sys.version.split()[0])
    print("PYTHONPATH            :", os.environ.get("PYTHONPATH"))
    print("modelseedpy.__file__  :", modelseedpy.__file__)
    print("prFBA source file     :", MSPROB_FILE)
    print("msprobability.py md5  :", MSPROB_MD5, "->", KNOWN_MD5.get(MSPROB_MD5, "UNKNOWN"))
    print("tree HEAD             :", HEAD, "| working file == HEAD blob:", md5_bytes(HEAD_SRC.encode()) == MSPROB_MD5)
    print("prFBA signature       :", inspect.signature(MSProbability.prFBA))
    print("prFBA def block starts:", "line", inspect.getsourcelines(MSProbability.prFBA)[1], "(its @staticmethod line)")
    mc_root = os.path.dirname(os.path.dirname(mscommunity.__file__))
    print("mscommunity           :", mscommunity.__file__, "HEAD", git(mc_root, "rev-parse", "--short", "HEAD").strip())
    print("cobra / optlang       :", cobra.__version__, "/", optlang.__version__)
    print("solver                :", cobra.Configuration().solver.__name__)
    print("stdout of importing mscommunity + msprobability (not prFBA's):", repr(_banner_mscomm.getvalue()[:120]))
    # the solver prints a license banner on first use: absorb it before any stdout measurement
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        cobra.io.load_model("textbook").slim_optimize()
    print("solver warm-up stdout (absorbed, not prFBA's):", len(buf.getvalue()), "chars")
    print()


# =========================================================================== fixtures
ENV = {c.strip(): 100.0 for c in open(MEDIA) if c.strip()}
META = json.load(open(FIXTURE_META))
GAPFILLED = set(META["block"]["new_rxn"])          # ids without the compartment index, e.g. rxn00068_c


def load_probe(objective="bio1", strip_gapfilled_notes=False, graded=False):
    """The gapfilled probe model with production-rule probabilities on its rxn reactions.
    objective: "bio1", or None to keep the objective it was built with (ATPM_c0)."""
    with gzip.open(FIXTURE, "rt") as fh:
        m = cobra.io.load_json_model(fh)
    for r in m.reactions:
        if r.id.startswith("rxn"):
            if graded and r.genes:   # deterministic fractional p for N3 only
                p = (0.25, 0.5, 0.75, 1.0)[int(hashlib.md5(r.id.encode()).hexdigest(), 16) % 4]
            else:
                p = 1.0 if r.genes else 0.0
            r.notes["probability"] = str(p)
    if strip_gapfilled_notes:
        for r in m.reactions:
            if r.id[:-1] in GAPFILLED:
                r.notes.pop("probability", None)
    if objective:
        m.objective = objective
    return m


def member(mid, carbon, uptake):
    """Three-reaction member (uptake, rxn transport, biomass) with formulas for the carbon cap:
    the two-member fixture of MSCommunity tests/test_f_biomass_drains.py, as used by b6857a0's
    tests/core/test_msprobability_prfba.py."""
    form = {carbon: "C6H12O6", "cpd00002": "C10H13N5O13P3", "cpd00001": "H2O",
            "cpd00008": "C10H13N5O10P2", "cpd00009": "HO4P", "cpd00067": "H", "cpd11416": "R"}
    m = cobra.Model(mid)
    mets = {k: cobra.Metabolite(k, name=k, compartment=k[-2:], formula=form[k[:-3]])
            for k in (f"{carbon}_e0", f"{carbon}_c0", "cpd00002_c0", "cpd00001_c0",
                      "cpd00008_c0", "cpd00009_c0", "cpd00067_c0", "cpd11416_c0")}
    m.add_metabolites(list(mets.values()))
    ex = cobra.Reaction(f"EX_{carbon}_e0", lower_bound=-float(uptake), upper_bound=1000)
    ex.add_metabolites({mets[f"{carbon}_e0"]: -1})
    tr = cobra.Reaction("rxn10000_c0", lower_bound=0, upper_bound=1000)
    tr.add_metabolites({mets[f"{carbon}_e0"]: -1, mets[f"{carbon}_c0"]: 1})
    tr.notes["probability"] = "0.5"
    bio = cobra.Reaction("bio1", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({mets[f"{carbon}_c0"]: -1, mets["cpd11416_c0"]: 1})
    sink = cobra.Reaction("SK_cpd11416_c0", lower_bound=0, upper_bound=1000)
    sink.add_metabolites({mets["cpd11416_c0"]: -1})
    m.add_reactions([ex, tr, bio, sink])
    m.objective = "bio1"
    return m


def community():
    return [member("memA", "cpd00027", 10), member("memB", "cpd00082", 4)]


COMM_ENV = {"cpd00027": 10, "cpd00082": 4}


def purge():
    """MSModelUtil.mdlutls and MSPackageManager.pkgmgrs are class-level dicts keyed by model:
    every model passed stays alive (and registered) for the life of the process. Cleared between
    cases so each case starts clean (except N7, which measures the growth)."""
    MSModelUtil.mdlutls.clear()
    MSPackageManager.pkgmgrs.clear()


# =========================================================================== reporting
def short_path(p):
    for root, tag in ((TREE, "<tree>"), (HERE, "<here>")):
        if p.startswith(root):
            return tag + p[len(root):]
    return p


def flux(sol, rid):
    try:
        return float(sol.fluxes[rid])
    except Exception:                               # noqa: BLE001
        return None


def sol_line(sol, extra=("bio1", "ATPM_c0")):
    parts = ["status={}".format(sol.status), "objective={:.6f}".format(sol.objective_value)
             if sol.objective_value is not None and not math.isnan(sol.objective_value)
             else "objective={}".format(sol.objective_value)]
    for rid in extra:
        v = flux(sol, rid)
        if v is not None:
            parts.append("{}={:.6f}".format(rid, v))
    return " ".join(parts)


def case(cid, title, fn):
    """Run fn(note) with stdout captured; print the outcome, the full traceback of any exception,
    what the call printed, and the lines fn put in `note` (state read after the call, printed
    whether or not it raised). fn returns a list of lines to report (or a Solution)."""
    print("-" * 100)
    print("[{}] {}".format(cid, title))
    if not cid.startswith("N7"):
        purge()
    buf = io.StringIO()
    note = []
    t0 = time.time()
    try:
        with contextlib.redirect_stdout(buf):
            res = fn(note)
        out = buf.getvalue()
        lines = [sol_line(res)] if hasattr(res, "fluxes") else list(res or [])
        print("  outcome : returned ({:.1f}s)".format(time.time() - t0))
        for ln in lines:
            print("    " + ln)
        SUMMARY.append((cid, "returned: " + (lines[0] if lines else "")))
    except Exception as exc:                        # noqa: BLE001
        out = buf.getvalue()
        name = type(exc).__module__ + "." + type(exc).__name__
        msg = str(exc)
        print("  outcome : RAISED {}: {}".format(name, msg[:600]))
        print("  traceback (most recent call last):")
        for fr in traceback.extract_tb(exc.__traceback__):
            print("    {}:{} in {}".format(short_path(fr.filename), fr.lineno, fr.name))
            if fr.line:
                print("        " + fr.line.strip()[:200])
        SUMMARY.append((cid, "RAISED {}: {}".format(type(exc).__name__, msg[:110])))
    for ln in note:
        print("  after   : " + ln)
    print("  stdout during the call: {} chars{}".format(
        len(out), "" if not out else " | first 160: " + repr(out[:160])))


def released_variant(tags):
    """The released module source with line-count-preserving edits, compiled under a pseudo file
    name (registered in linecache, so tracebacks show the edited text at the released line
    numbers). Returns prFBA, or None when an edit target is absent (the tree is not c2e31a0)."""
    edits = {
        "F1": ("        mdlUtil = MSModelUtil(model_s_ if len(model_s_) == 1 else",
               "        model_s_ = model_s_ if isinstance(model_s_, (list, tuple)) else [model_s_]; "
               "mdlUtil = MSModelUtil(model_s_[0] if len(model_s_) == 1 else"),
        "F2": ("        # constrain carbon consumption\n"
               "        elepkg = ElementUptakePkg(mdlUtil.model)  ;  elepkg.build_package({\"C\": 100})\n",
               "        # [F2: the duplicate ElementUptakePkg build of l.173-174 is removed]\n"
               "        pass\n"),
    }
    src = MSPROB_SRC
    for t in tags:
        old, new = edits[t]
        if src.count(old) != 1:
            return None
        src = src.replace(old, new)
    assert src.count("\n") == MSPROB_SRC.count("\n"), "edits must preserve line numbers"
    fname = "msprobability.py[c2e31a0+{}]".format("+".join(tags))
    linecache.cache[fname] = (len(src), None, src.splitlines(True), fname)
    mod = types.ModuleType("msprobability_variant_" + "_".join(tags))
    mod.__file__ = fname
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(src, fname, "exec"), mod.__dict__)
    return mod.MSProbability.prFBA


def constraint_names(m):
    return {c.name for c in m.constraints}


def floor_vars(m):
    if "min_value" not in m.constraints:
        return None
    return sorted(v.name.split("_reverse_")[0] + ("(rev)" if "_reverse_" in v.name else "")
                  for v in m.constraints["min_value"].variables)


# =========================================================================== Part A
def part_a():
    print("\n" + "#" * 100 + "\n# PART A: the released call forms (both trees)\n" + "#" * 100)
    P = MSProbability.prFBA

    def fp1a(note):
        m = load_probe("bio1")
        return [sol_line(P(m, environment=ENV))]
    case("FP1a", "prFBA(model, environment=ENV): a bare cobra.Model (bio1 set as objective)", fp1a)

    def fp1b(note):
        m = load_probe("bio1")
        return [sol_line(P([m], environment=ENV))]
    case("FP1b", "prFBA([model], environment=ENV): a one-element list", fp1b)

    def fp1c(note):
        m = load_probe("bio1")
        return [sol_line(P([m]))]
    case("FP1c", "prFBA([model]): one-element list, environment=None (the default)", fp1c)

    def fp2b(note):
        m = load_probe("bio1")
        return ["add_medium(None) returned {}".format(type(MSModelUtil(m).add_medium(None)).__name__)]
    case("FP2b", "MSModelUtil(model).add_medium(None) by itself (the call prFBA l.166 makes when environment=None)", fp2b)

    def fp3b(note):
        ms = community()
        return [sol_line(P(ms, environment=COMM_ENV), extra=("bio1",))]
    case("FP3b", "prFBA([memA, memB], environment=...): the community path (two members)", fp3b)


# =========================================================================== Part B
def part_b():
    print("\n" + "#" * 100 + "\n# PART B: c2e31a0 only. Defects hidden behind FP1/FP3, reached by "
          "line-preserving source edits\n" + "#" * 100)
    v1, v12, v2 = released_variant(["F1"]), released_variant(["F1", "F2"]), released_variant(["F2"])
    if v1 is None:
        print("n/a on this tree: the released text of l.165 / l.173-174 is absent (the defects are "
              "fixed in source; see Parts A and C).")
        for cid in ("FP3a", "FP2", "FP4", "FP5", "FP6", "FP3c"):
            SUMMARY.append((cid, "n/a (edit target absent: fixed in source)"))
        return

    def fp3a(note):
        m = load_probe("bio1")
        return [sol_line(v1(m, environment=ENV))]
    case("FP3a", "[F1] bare model, environment=ENV: the second ElementUptakePkg on the single-model path", fp3a)

    def fp2(note):
        m = load_probe("bio1")
        return [sol_line(v1(m))]
    case("FP2", "[F1] bare model, environment=None (the default)", fp2)

    def fp4(note):
        m = load_probe("bio1")
        note.append("the stdout below is prFBA's own print of every constraint name (l.194)")
        return [sol_line(v12(m, environment=ENV))]
    case("FP4", "[F1+F2] bare model, environment=ENV, bio1 objective: runs; what does it print?", fp4)

    def fp5(note):
        m = load_probe("bio1", strip_gapfilled_notes=True)
        n_unnoted = sum(1 for r in m.reactions if r.id.startswith("rxn") and "probability" not in r.notes)
        before = (constraint_names(m), str(m.objective.expression), dict(m.medium))
        note.append("rxn reactions without notes['probability'] before the call: {} (the gapfilled ones)".format(n_unnoted))
        try:
            return [sol_line(v12(m, environment=ENV))]
        finally:
            added = sorted(constraint_names(m) - before[0])
            note.append("the caller's model was changed before the error: constraints added {}, objective "
                        "changed {}, medium changed {} ({} -> {} entries)".format(
                            added, str(m.objective.expression) != before[1], dict(m.medium) != before[2],
                            len(before[2]), len(m.medium)))
    case("FP5", "[F1+F2] gapfilled reactions carry no notes['probability'] (production order)", fp5)

    def fp6(note):
        m = load_probe(None)
        obj = str(m.objective.expression)
        sol = v12(m, environment=ENV)
        return [sol_line(sol), "objective on entry: " + obj,
                "floor (min_value) is on: {} lb={:.6f}".format(floor_vars(m), m.constraints["min_value"].lb)]
    case("FP6", "[F1+F2] the model as built (objective ATPM_c0): which objective gets the 95% floor?", fp6)

    def fp3c(note):
        ms = community()
        return [sol_line(v2(ms, environment=COMM_ENV), extra=("bio1",))]
    case("FP3c", "[F2] community path once the duplicate build is removed", fp3c)


# =========================================================================== Part C
def part_c():
    print("\n" + "#" * 100 + "\n# PART C: b6857a0 keywords and defaults (both trees)\n" + "#" * 100)
    P = MSProbability.prFBA

    def c1(note):
        m = load_probe("bio1", strip_gapfilled_notes=True)
        before = (constraint_names(m), str(m.objective.expression), dict(m.medium))
        try:
            return [sol_line(P(m, environment=ENV))]
        finally:
            note.append("the caller's model is unchanged: constraints {}, objective {}, medium {}".format(
                constraint_names(m) == before[0], str(m.objective.expression) == before[1],
                dict(m.medium) == before[2]))
    case("C1", "default missing_probability on un-noted gapfilled reactions", c1)

    def c2(note):
        m = load_probe("bio1", strip_gapfilled_notes=True)
        sol = P(m, environment=ENV, missing_probability=0.0)
        m0 = load_probe("bio1")      # the same reactions annotated "0.0"
        sol0 = P(m0, environment=ENV)
        rel = abs(sol.objective_value - sol0.objective_value) / max(1.0, abs(sol0.objective_value))
        return [sol_line(sol), "annotated-with-0.0 run: " + sol_line(sol0), "relative difference: {:.3g}".format(rel)]
    case("C2", "missing_probability=0.0 on un-noted gapfilled reactions", c2)

    def c3(note):
        m = load_probe(None)
        sol = P(m, environment=ENV, biomass_id="bio1")
        return [sol_line(sol), "floor on: {} lb={:.6f}".format(floor_vars(m), m.constraints["min_value"].lb)]
    case("C3", "biomass_id='bio1' on the model as built (objective ATPM_c0)", c3)

    def c4(note):
        m = load_probe(None)
        sol = P(m, environment=ENV)
        return [sol_line(sol), "floor on: {}".format(floor_vars(m))]
    case("C4", "default biomass_id=None on the model as built: the floor stays on ATPM_c0", c4)

    def c5(note):
        m_quiet, m_loud = load_probe("bio1"), load_probe("bio1")
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as td:
            os.chdir(td)
            try:
                b1 = io.StringIO()
                with contextlib.redirect_stdout(b1):
                    P(m_quiet, environment=ENV)
                lp1 = os.path.exists("prFBA.lp")
                b2 = io.StringIO()
                with contextlib.redirect_stdout(b2):
                    P(m_loud, environment=ENV, printLP=True)
                lp2 = os.path.exists("prFBA.lp")
            finally:
                os.chdir(cwd)
        return ["default: printed {} chars, prFBA.lp written: {}".format(len(b1.getvalue()), lp1),
                "printLP=True: printed {} chars ('min_value' in it: {}), prFBA.lp written: {}".format(
                    len(b2.getvalue()), "min_value" in b2.getvalue(), lp2)]
    case("C5", "printing: default vs printLP=True (run in a temporary directory)", c5)

    def c6(note):
        m = load_probe("bio1")
        own = m.slim_optimize()      # bio1 max under the fixture's own exchange bounds (no change to m)
        sol = P(m)
        return [sol_line(sol), "own exchange bounds: bio1 max = {:.6f} (the fixture as built)".format(own),
                "floor lb = {}".format(m.constraints["min_value"].lb),
                "nonzero fluxes in the returned solution: {}".format(int((sol.fluxes.abs() > 1e-9).sum()))]
    case("C6", "environment=None (default) on the gapfilled model whose own bounds admit no growth", c6)

    def c7(note):
        m = load_probe("bio1")
        P(m, environment=ENV)
        note.append("first call returned; calling again on the same object")
        P(m, environment=ENV)
        return ["second call returned"]
    case("C7", "a second call on the same model object", c7)

    def c8(note):
        return [sol_line(P([], environment=ENV))]
    case("C8", "an empty list", c8)


# =========================================================================== Part D
def lp_state(m):
    """|v| per reaction (forward + reverse primal) after the last solve."""
    return {r.id: max(0.0, r.forward_variable.primal) + max(0.0, r.reverse_variable.primal)
            for r in m.reactions}


def resolve_with(m, coef):
    m.objective.set_linear_coefficients(coef)
    m.solver.update()
    m.slim_optimize(error_value=float("nan"))
    return m.solver.objective.value, lp_state(m), {r.id: r.flux for r in m.reactions}


def face_max(m, opt, variables, slack=1e-6):
    """Max of sum(variables) over {current objective <= opt + slack}: how much flux the optimal
    face allows through these variables."""
    expr = m.solver.objective.expression
    with m:
        cap = m.problem.Constraint(expr, ub=opt + slack, name="face_cap")
        m.add_cons_vars(cap)
        m.objective = m.problem.Objective(optlang.symbolics.Zero, direction="max")
        m.objective.set_linear_coefficients({v: 1.0 for v in variables})
        val = m.slim_optimize(error_value=float("nan"))
    return val


def part_d():
    print("\n" + "#" * 100 + "\n# PART D: further problems (both trees)\n" + "#" * 100)
    v12 = released_variant(["F1", "F2"])
    P = v12 if v12 is not None else MSProbability.prFBA
    print("runner: {}".format("c2e31a0 source with F1+F2 (Part B)" if v12 is not None
                              else "MSProbability.prFBA of this tree"))
    min_prob, prob_exp, ex_weight = 0.01, 1, 100

    # ---- N1 silent zero: the floor is 0 when the medium admits no growth
    def n1(note):
        m = load_probe("bio1")
        sol = P(m, environment={})
        return [sol_line(sol), "floor lb = {}; every flux zero: {}".format(
            m.constraints["min_value"].lb, bool((sol.fluxes.abs() < 1e-9).all()))]
    case("N1", "environment={} (every uptake closed): the floor becomes 0 and the call 'succeeds'", n1)

    # ---- N2 NaN floor: slim_optimize() of an infeasible model
    def n2(note):
        m = cobra.io.load_model("textbook")      # ATPM lower bound 8.39: infeasible without uptake
        for r in m.reactions:
            r.notes["probability"] = "1.0"
        sol = P(m, environment={})
        lb = m.constraints["min_value"].lb if "min_value" in m.constraints else None
        return [sol_line(sol, extra=("Biomass_Ecoli_core",)), "floor lb = {}".format(lb)]
    case("N2", "e_coli_core, environment={} (infeasible: ATPM lb 8.39 without uptake): slim_optimize -> NaN", n2)

    # ---- N3 comment vs code: the three readings of the objective
    def n3(note):
        out = []
        for graded in (False, True):
            m = load_probe("bio1", graded=graded)
            sol = P(m, environment=ENV)
            opt_code = sol.objective_value
            v_code = lp_state(m)
            rx = [r for r in m.reactions if r.id.startswith("rxn")]
            p = {r.id: float(r.notes["probability"]) for r in rx}
            certain = [r for r in rx if p[r.id] >= 1.0]
            fwdrev = lambda rs: [v for r in rs for v in (r.forward_variable, r.reverse_variable)]
            tag = "graded p in {0.25,0.5,0.75,1} on gene-backed" if graded else "two-tier p (single genome)"
            out.append("--- {}: {} rxn reactions, {} with p = 1".format(tag, len(rx), len(certain)))

            def stats(name, opt, v, extra=""):
                act = sum(1 for r in rx if v[r.id] > 1e-6)
                l1 = sum(abs(v[k] - v_code[k]) for k in v)
                out.append("  {:<34s} objective={:<14.6f} active rxn={:<4d} sum|v| rxn={:<11.4f} "
                           "sum|v| p=1 rxn={:<11.4f} L1 to code={:.4f}{}".format(
                               name, opt, act, sum(v[r.id] for r in rx), sum(v[r.id] for r in certain), l1, extra))
            stats("code: max(min_prob, 1-p^e)", opt_code, v_code)
            # additive reading: coefficient (1 - p^e) + min_prob
            o_add, v_add, _ = resolve_with(m, {var: (1 - p[r.id] ** prob_exp) + min_prob
                                               for r in rx for var in (r.forward_variable, r.reverse_variable)})
            stats("additive: (1-p^e) + min_prob", o_add, v_add)
            # literal comment: sum((1-p^e) v + min_prob) -> the + min_prob is a constant
            o_lit, v_lit, _ = resolve_with(m, {var: (1 - p[r.id] ** prob_exp)
                                               for r in rx for var in (r.forward_variable, r.reverse_variable)})
            fm = face_max(m, o_lit, fwdrev(certain))
            stats("literal comment: (1-p^e)v + const", o_lit, v_lit,
                  " | face max sum|v| p=1 rxn={:.1f}".format(fm))
            # restore the code coefficients (and confirm the face max under them)
            o_back, _, _ = resolve_with(m, {var: max(min_prob, 1 - p[r.id] ** prob_exp)
                                            for r in rx for var in (r.forward_variable, r.reverse_variable)})
            fm_code = face_max(m, o_back, fwdrev(certain))
            out.append("  code coefficients restored: objective={:.6f}; face max sum|v| p=1 rxn={:.4f}".format(o_back, fm_code))
        return out
    case("N3", "the comment's objective vs the code's: additive and literal readings on the same constraints", n3)

    # ---- N4 unpriced reaction classes
    def n4(note):
        out = []
        m = load_probe("bio1")
        sol = P(m, environment=ENV)
        opt = sol.objective_value
        coefs = m.objective.get_linear_coefficients([v for r in m.reactions for v in (r.forward_variable, r.reverse_variable)])
        unpriced = [r for r in m.reactions if coefs[r.forward_variable] == 0 and coefs[r.reverse_variable] == 0]
        out.append("probe model: {} reactions; unpriced (coefficient 0 both directions): {}".format(
            len(m.reactions), len(unpriced)))
        for r in unpriced:
            out.append("  {:<16s} bounds={} flux={:.6f} | {}".format(r.id, r.bounds, sol.fluxes[r.id],
                                                                   r.build_reaction_string()[:70]))
        cel = {v.name for v in m.constraints["C_elements"].variables}
        out.append("C_elements (carbon cap) counts any SK_ reaction: {}; counts {} EX_ variables".format(
            any(n.startswith("SK_") for n in cel), sum(1 for n in cel if n.startswith("EX_"))))
        sk = [m.reactions.get_by_id(i) for i in ("SK_cpd02701_c0", "SK_cpd03091_c0")]
        o2, _, f2 = resolve_with(m, {r.forward_variable: ex_weight for r in sk})
        out.append("price SK_cpd02701_c0 / SK_cpd03091_c0 like an exchange (100): objective {:.6f} -> {:.6f}; "
                   "their flux {} -> {}".format(opt, o2, [round(sol.fluxes[r.id], 6) for r in sk],
                                                [round(f2[r.id], 6) for r in sk]))
        # a model whose reaction ids are not ModelSEED's: nothing internal is priced
        objs, fluxes = {}, {}
        for pval in ("1.0", "0.0"):
            t = cobra.io.load_model("textbook")
            for r in t.reactions:
                r.notes["probability"] = pval
            s = P(t, environment=dict(t.medium))
            objs[pval] = s.objective_value
            fluxes[pval] = s.fluxes
            if pval == "1.0":
                internal = [r for r in t.reactions if not r.id.startswith("EX_")]
                c = t.objective.get_linear_coefficients([r.forward_variable for r in internal])
                n_priced = sum(1 for x in c.values() if x != 0)
                at_vertex = sum(abs(s.fluxes[r.id]) for r in internal if r.id != "Biomass_Ecoli_core")
                fmx = face_max(t, s.objective_value, [v for r in internal if r.id != "Biomass_Ecoli_core"
                                                      for v in (r.forward_variable, r.reverse_variable)])
        out.append("e_coli_core (BiGG ids, every reaction annotated): internal reactions priced: {} of {}".format(
            n_priced, len(internal)))
        out.append("  objective with every p = 1: {:.6f}; with every p = 0: {:.6f}; identical flux vectors: {}".format(
            objs["1.0"], objs["0.0"], bool((fluxes["1.0"] - fluxes["0.0"]).abs().max() < 1e-9)))
        out.append("  internal sum|v| at the returned vertex {:.2f}; max on the optimal face {:.2f}".format(at_vertex, fmx))
        return out
    case("N4", "reactions left at zero cost (not 'rxn', not 'EX_'), and a non-ModelSEED model", n4)

    # ---- N5 exchanges priced in both directions
    def n5(note):
        out = []
        m = load_probe("bio1")
        sol = P(m, environment=ENV)
        opt = sol.objective_value
        ex = [r for r in m.reactions if r.id.startswith("EX_")]
        rx = [r for r in m.reactions if r.id.startswith("rxn")]
        v = lp_state(m)
        coefs = m.objective.get_linear_coefficients([x for r in rx for x in (r.forward_variable, r.reverse_variable)])
        ex_term = sum(ex_weight * v[r.id] for r in ex)
        int_term = sum(coefs[r.forward_variable] * max(0, r.forward_variable.primal)
                       + coefs[r.reverse_variable] * max(0, r.reverse_variable.primal) for r in rx)
        upt = sum(max(0.0, r.reverse_variable.primal) for r in ex)
        sec = sum(max(0.0, r.forward_variable.primal) for r in ex)
        out.append("objective {:.6f} = exchange term {:.6f} ({:.2%}) + internal term {:.6f}".format(
            opt, ex_term, ex_term / opt, int_term))
        out.append("exchange flux: uptake {:.4f}, secretion {:.4f}; secretion share of the exchange term {:.2%}".format(
            upt, sec, sec / (upt + sec)))
        top = sorted(ex, key=lambda r: -v[r.id])[:8]
        out.append("largest exchanges (net flux; negative = uptake): " + ", ".join(
            "{} {:+.3f}".format(r.id, sol.fluxes[r.id]) for r in top))
        hw = sum(ex_weight * v[i] for i in ("EX_cpd00067_e0", "EX_cpd00001_e0") if i in v)
        out.append("H+ (cpd00067) + H2O (cpd00001) exchanges: {:.2%} of the exchange term".format(hw / ex_term))
        sk_before = {i: round(sol.fluxes[i], 6) for i in ("SK_cpd02701_c0", "SK_cpd03091_c0")}
        # counterfactual: price uptake only (secretion free), same constraints
        o_u, v_u, f_u = resolve_with(m, {r.forward_variable: 0.0 for r in ex})
        upt_u = sum(max(0.0, r.reverse_variable.primal) for r in ex)
        sec_u = sum(max(0.0, r.forward_variable.primal) for r in ex)
        out.append("exchange flux, uptake-only pricing: uptake {:.4f} (code {:.4f}), secretion {:.4f} (code {:.4f}); "
                   "bio1 {:.6f}".format(upt_u, upt, sec_u, sec, f_u["bio1"]))
        act = lambda vv: sum(1 for r in rx if vv[r.id] > 1e-6)
        out.append("uptake-only pricing: objective {:.6f}; L1 change of the flux vector {:.4f}; active rxn {} -> {}; "
                   "internal sum|v| {:.4f} -> {:.4f}".format(o_u, sum(abs(f_u[k] - sol.fluxes[k]) for k in f_u),
                                                         act(v), act(v_u), sum(v[r.id] for r in rx), sum(v_u[r.id] for r in rx)))
        out.append("  secreted (code): " + ", ".join("{} {:.3f}".format(r.id, sol.fluxes[r.id]) for r in ex if sol.fluxes[r.id] > 1e-6))
        out.append("  secreted (uptake-only): " + ", ".join("{} {:.3f}".format(r.id, f_u[r.id]) for r in ex if f_u[r.id] > 1e-6))
        out.append("  SK_ sink flux: code {} | uptake-only {}".format(
            sk_before, {i: round(f_u[i], 6) for i in sk_before}))
        return out
    case("N5", "exchanges priced ex_weight in both directions: what the exchange term prices", n5)

    # ---- N6 community path: where the kinetic rows go; commkinetics; abundances
    def n6(note):
        out = []
        seen = []
        Orig = MCS.MSCommunity

        class Spy(Orig):
            def __init__(self, model=None, *a, **k):
                super().__init__(model, *a, **k)
                seen.append({"arg": model, "util": self.util.model, "kw": dict(k)})
        kin = lambda mm: sorted(c.name for c in mm.constraints if "commkin" in c.name.lower())
        MCS.MSCommunity = Spy
        try:
            res = {}
            for label, kw in (("kinetics_coef=1000 (default)", {}), ("kinetics_coef=1", {"kinetics_coef": 1}),
                              ("commkinetics=True", {"commkinetics": True}),
                              ("abundances memA 0.9 / memB 0.1", {"abundances": {"memA": 0.9, "memB": 0.1}})):
                purge()
                del seen[:]
                ms = community()
                captured = {}
                orig_bfsm = mscommunity.commhelper.build_from_species_models

                def bfsm(*a, **k):
                    captured["model"] = orig_bfsm(*a, **k)
                    return captured["model"]
                mscommunity.commhelper.build_from_species_models = bfsm
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        sol = P(ms, environment=COMM_ENV, **kw)
                finally:
                    mscommunity.commhelper.build_from_species_models = orig_bfsm
                solved = captured["model"]
                res[label] = sol.objective_value
                line = "{:<34s} {} | MSCommunity built: {}".format(label, sol_line(sol, extra=("bio1",)), len(seen))
                if seen:
                    s0 = seen[0]
                    line += (" | MSCommunity was handed the solved model: {}; it constrains that object: {}"
                             " | commKin rows on the solved model: {} ({} constraints), on MSCommunity's "
                             "copy: {} ({} constraints)").format(
                        s0["arg"] is solved, s0["util"] is solved, kin(solved), len(solved.constraints),
                        kin(s0["util"]), len(s0["util"].constraints))
                else:
                    line += " | commKin rows on solved model: {}".format(kin(solved))
                out.append(line)
                if label.startswith("abundances"):
                    bio_s = {mt.id: c for mt, c in solved.reactions.bio1.metabolites.items()}
                    bio_c = {mt.id: c for mt, c in seen[0]["util"].reactions.bio1.metabolites.items()}
                    out.append("  community bio1 on the solved model: {} | on MSCommunity's copy: {}".format(bio_s, bio_c))
                if label.startswith("kinetics_coef=1000"):
                    unnoted = [r.id for r in solved.reactions if r.id.startswith("rxn") and "probability" not in r.notes]
                    unpriced = [r.id for r in solved.reactions if not r.id.startswith(("rxn", "EX_"))]
                    out.append("  merged (solved) model: {} reactions; rxn reactions without a note: {}; unpriced: {}".format(
                        len(solved.reactions), unnoted, unpriced))
                    copy_only = sorted({r.id for r in seen[0]["util"].reactions} - {r.id for r in solved.reactions})
                    out.append("  reactions only on MSCommunity's copy: {} (notes: {})".format(
                        copy_only, {i: seen[0]["util"].reactions.get_by_id(i).notes.get("probability") for i in copy_only}))
            out.append("objective with kinetics_coef=1000 vs 1: {:.6f} vs {:.6f} (identical: {})".format(
                res["kinetics_coef=1000 (default)"], res["kinetics_coef=1"],
                abs(res["kinetics_coef=1000 (default)"] - res["kinetics_coef=1"]) < 1e-9))
            src = inspect.getsource(mscommunity.commhelper.build_from_species_models)
            out.append("build_from_species_models mentions 'commkinetics' {} time(s) (its signature only)".format(
                src.count("commkinetics")))
        finally:
            MCS.MSCommunity = Orig
        return out
    case("N6", "community path: the model MSCommunity constrains vs the model prFBA solves", n6)

    # ---- N7 in-place mutation and process-lifetime caches
    def n7(note):
        purge()
        n0 = (len(MSPackageManager.pkgmgrs), len(MSModelUtil.mdlutls))
        for _ in range(3):
            P(load_probe("bio1"), environment=ENV)
        n1 = (len(MSPackageManager.pkgmgrs), len(MSModelUtil.mdlutls))
        out = ["MSPackageManager.pkgmgrs / MSModelUtil.mdlutls entries: before {} after 3 calls {}".format(n0, n1)]
        m = load_probe("bio1")
        P(m, environment=ENV)
        out.append("the caller's model after one call: objective direction={}, constraints added={}, "
                   "medium entries={}".format(m.objective.direction,
                                              sorted(n for n in constraint_names(m) if n in ("min_value", "C_elements")),
                                              len(m.medium)))
        try:
            P(m, environment=ENV)
            out.append("second call on the same object: returned")
        except Exception as exc:                    # noqa: BLE001
            out.append("second call on the same object: {}: {}".format(type(exc).__name__, exc))
        purge()
        return out
    case("N7", "in-place mutation and class-level caches", n7)

    # ---- N8 probability parsing
    def n8(note):
        m = load_probe("bio1")
        picks = sorted(r.id for r in m.reactions if r.id.startswith("rxn") and r.genes)[:4]
        vals = ("nan", "1.5", "-0.5", "0.999")
        for rid, val in zip(picks, vals):
            m.reactions.get_by_id(rid).notes["probability"] = val
        P(m, environment=ENV)
        c = m.objective.get_linear_coefficients([m.reactions.get_by_id(r).forward_variable for r in picks])
        return ["notes['probability']={!r:<7s} -> coefficient {:.4f}".format(val, c[m.reactions.get_by_id(rid).forward_variable])
                for rid, val in zip(picks, vals)]
    case("N8", "unvalidated probability strings (coefficient = max(min_prob, 1 - float(p)**prob_exp))", n8)


# =========================================================================== main
def main():
    header()
    part_a()
    part_b()
    part_c()
    part_d()
    print("\n" + "=" * 100 + "\nSUMMARY ({}; msprobability.py md5 {})".format(KNOWN_MD5.get(MSPROB_MD5, "?"), MSPROB_MD5))
    print("=" * 100)
    for cid, line in SUMMARY:
        print("{:<6s} {}".format(cid, line[:150]))


if __name__ == "__main__":
    main()
