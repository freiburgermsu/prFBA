"""Shared paths, pinned library sources and provenance for the downstream-E3 measurement.

Everything here only READS the deployments (codiffusion_bioreactor, EmilyKin) and the library
trees; all outputs go under OUT (this directory) or SCRATCH.

Pinned sources. The shared library trees move (other sessions switch branches and edit them),
so each run imports from `git archive` exports of fixed commits, made once into PINS:
  ModelSEEDpy        09856f5  branch fix-prfba-single-model (c2e31a0 + the MSProbability.prFBA
                              fix prfba_objective.run_prfba requires)
  MSCommunity        03bd80c  (build_from_species_models; ~/Documents/MSCommunity HEAD 0532529
                              carries uncommitted edits to commhelper.py)
  ModelSEEDTemplates b6fe86b  codif-curated: codiffusion's production templates
                              (rebuild_pipeline.sh, TAG v7e)                -> part 2
  ModelSEEDTemplates 21a8084  c12run: the templates of Stage O, where the tuned prFBA objective
                              was selected                                  -> part 3
  ModelSEEDDatabase  5389554a Biochemistry/ + Media/KBaseMedia.cpd (biochemistry preload for
                              MSATPCorrection; the complete medium of phenotype_arm.complete_media)
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time

DOC = "/home/freiburger/Documents"
PRFBA = f"{DOC}/prFBA"
OUT = f"{PRFBA}/deployment_d100/downstream_e3"
SCRATCH = ("/tmp/claude-1000/-home-freiburger-Documents-prFBA/"
           "400ec395-f172-59da-a98d-3dbe6b3a4256/scratchpad/downstream")
PINS = os.environ.get("DS_E3_PINS", f"{SCRATCH}/pins")
PY = f"{DOC}/py_venv/bin/python"
CODIF = f"{DOC}/codiffusion_bioreactor"
EK = f"{DOC}/EmilyKin"

BVBRC_KN = ("/tmp/claude-1000/-home-freiburger-Documents-prFBA/"
            "400ec395-f172-59da-a98d-3dbe6b3a4256/scratchpad/bvbrc_only")

# Build sets: which k/n and E3 build each deployment is measured on, and where that measurement's
# results and models go. Every script takes --builds NAME (required).
#   mixed  the shipped k/n builds and their E3 re-weightings (2026-11-12). Both merged BV-BRC's
#          PATRIC CDSs and the imported RefSeq CDSs, the rule up to prFBA 3ce12d0. The two _e3
#          directories were deleted on 2026-11-13 (superseded); results/ and models/ are their
#          record. Regenerable from prFBA 3ce12d0 with reweight_build_e3.py.
#   bvbrc  the BV-BRC (PATRIC) annotation only rebuilds of 2026-11-13: build_d100.py --selection
#          <shipped audit>, once with --pfunc kn (BVBRC_KN, the session scratchpad) and once with
#          E3 (synthetic_genomes_d100_bvbrc). Selections are the shipped ones.
BUILD_SETS = {
    "mixed": {
        "label": "mixed PATRIC + RefSeq vocabulary (shipped k/n builds vs their E3 re-weighting)",
        "builds": {"codiffusion": {"kn": f"{CODIF}/synthetic_genomes_d100",
                                   "e3": f"{CODIF}/synthetic_genomes_d100_e3"},
                   "emilykin": {"kn": f"{EK}/synthetic_genomes_d100",
                                "e3": f"{EK}/synthetic_genomes_d100_e3"}},
        "results": f"{OUT}/results", "models": f"{OUT}/models"},
    "bvbrc": {
        "label": "BV-BRC (PATRIC) annotation only (k/n rebuild vs E3 rebuild)",
        "builds": {"codiffusion": {"kn": f"{BVBRC_KN}/kn_codif",
                                   "e3": f"{CODIF}/synthetic_genomes_d100_bvbrc"},
                   "emilykin": {"kn": f"{BVBRC_KN}/kn_emilykin",
                                "e3": f"{EK}/synthetic_genomes_d100_bvbrc"}},
        "results": f"{OUT}/results_bvbrc", "models": f"{OUT}/models_bvbrc"},
}
BUILD_SET_ENV = "DS_E3_BUILDS"


def _builds_from_argv(argv):
    for i, x in enumerate(argv):
        if x == "--builds" and i + 1 < len(argv):
            return argv[i + 1]
        if x.startswith("--builds="):
            return x.split("=", 1)[1]
    return None


def use_build_set(name):
    """Select a build set for this process and every worker it spawns (through the environment).

    It runs once at import, from --builds on the command line or else DS_E3_BUILDS, so the
    scripts' module-level paths are right in the parent and in spawned workers alike (a spawned
    worker re-imports the main module with the parent's sys.argv and environment)."""
    global BUILD_SET, BUILDS, RESULTS, MODELS_ROOT
    if name is not None and name not in BUILD_SETS:
        raise SystemExit(f"--builds must be one of {sorted(BUILD_SETS)}, not {name!r}")
    BUILD_SET = name
    spec = BUILD_SETS.get(name) or {}
    BUILDS = spec.get("builds")
    RESULTS = spec.get("results")
    MODELS_ROOT = spec.get("models")
    if name is not None:
        os.environ[BUILD_SET_ENV] = name


def require_builds(arms=("kn", "e3")):
    """Fail early, with the reason, if a build directory of the selected set is missing."""
    missing = [BUILDS[d][a] for d in BUILDS for a in arms if not os.path.isdir(BUILDS[d][a])]
    if missing:
        raise SystemExit(f"build set {BUILD_SET!r}: missing build directories {missing}"
                         + (" (the mixed E3 builds were deleted on 2026-11-13; regenerate them "
                            "with reweight_build_e3.py at prFBA 3ce12d0)" if BUILD_SET == "mixed" else ""))


use_build_set(_builds_from_argv(sys.argv) or os.environ.get(BUILD_SET_ENV))

PIN_SPECS = {
    # name: (repository, commit, paths inside the repo to export)
    "ModelSEEDpy-09856f5": (f"{DOC}/ModelSEEDpy-prfba", "09856f5", ["modelseedpy"]),
    "MSCommunity-03bd80c": (f"{DOC}/MSCommunity", "03bd80c", ["mscommunity"]),
    "ModelSEEDTemplates-b6fe86b": (f"{DOC}/ModelSEEDTemplates-codif-curated", "b6fe86b",
                                   ["templates/v7.0/GramNegModelTemplateV7.json",
                                    "templates/v7.0/GramPosModelTemplateV7.json",
                                    "templates/v7.0/ArchaeaTemplateV6.json",
                                    "templates/v7.0/Core-V6.json"]),
    "ModelSEEDTemplates-21a8084": (f"{DOC}/ModelSEEDTemplates-c12run", "21a8084",
                                   ["templates/v7.0/GramNegModelTemplateV7.json",
                                    "templates/v7.0/GramPosModelTemplateV7.json",
                                    "templates/v7.0/ArchaeaTemplateV6.json",
                                    "templates/v7.0/Core-V6.json"]),
    "ModelSEEDDatabase-5389554a": (f"{DOC}/ModelSEEDDatabase", "5389554a",
                                   ["Biochemistry", "Media/KBaseMedia.cpd"]),
}

TPL_ROOT = {"codif_curated": f"{PINS}/ModelSEEDTemplates-b6fe86b/templates/v7.0",
            "c12run": f"{PINS}/ModelSEEDTemplates-21a8084/templates/v7.0"}
MSD = f"{PINS}/ModelSEEDDatabase-5389554a"
PYTHONPATH = f"{PINS}/ModelSEEDpy-09856f5:{PINS}/MSCommunity-03bd80c"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sha256_file(path, n=None):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n] if n else h.hexdigest()


def git(repo, *args):
    env = dict(os.environ, GIT_OPTIONAL_LOCKS="0")
    return subprocess.run(["git", "-C", repo, *args], capture_output=True, text=True,
                          env=env, timeout=600).stdout.strip()


def ensure_pins():
    """Create the git-archive exports (read-only on the source repositories)."""
    os.makedirs(PINS, exist_ok=True)
    env = dict(os.environ, GIT_OPTIONAL_LOCKS="0")
    for name, (repo, commit, paths) in PIN_SPECS.items():
        dest = f"{PINS}/{name}"
        stamp = f"{dest}/.export_ok"
        if os.path.exists(stamp):
            continue
        os.makedirs(dest, exist_ok=True)
        full = git(repo, "rev-parse", commit)
        arch = subprocess.Popen(["git", "-C", repo, "archive", "--format=tar", full, *paths],
                                stdout=subprocess.PIPE, env=env)
        subprocess.run(["tar", "-x", "-C", dest], stdin=arch.stdout, check=True)
        if arch.wait() != 0:
            raise RuntimeError(f"git archive failed for {name}")
        with open(stamp, "w") as fh:
            json.dump({"repo": repo, "commit": full, "paths": paths,
                       "exported": time.strftime("%Y-%m-%dT%H:%M:%S")}, fh)
        log(f"pinned {name} <- {repo}@{full[:10]}")


def pin_record():
    out = {}
    for name in PIN_SPECS:
        stamp = f"{PINS}/{name}/.export_ok"
        out[name] = json.load(open(stamp)) if os.path.exists(stamp) else None
    for key, root in TPL_ROOT.items():
        out[f"templates_{key}_sha256"] = {f: sha256_file(f"{root}/{f}", 16)
                                          for f in sorted(os.listdir(root)) if f.endswith(".json")}
    return out


def check_imports():
    """Raise unless modelseedpy and mscommunity come from the pinned exports."""
    import modelseedpy
    import mscommunity
    bad = []
    if not modelseedpy.__file__.startswith(f"{PINS}/ModelSEEDpy-09856f5/"):
        bad.append(f"modelseedpy from {modelseedpy.__file__}")
    if not mscommunity.__file__.startswith(f"{PINS}/MSCommunity-03bd80c/"):
        bad.append(f"mscommunity from {mscommunity.__file__}")
    if bad:
        raise RuntimeError("unpinned import: " + "; ".join(bad) + f" (PYTHONPATH={PYTHONPATH})")
    return {"modelseedpy": modelseedpy.__file__, "mscommunity": mscommunity.__file__}


def code_record(scripts):
    rec = {"prfba_head": git(PRFBA, "rev-parse", "--short", "HEAD"),
           "prfba_dirty_files": len(git(PRFBA, "status", "--porcelain").splitlines()),
           "codiffusion_head": git(CODIF, "rev-parse", "--short", "HEAD"),
           "emilykin_head": git(EK, "rev-parse", "--short", "HEAD"),
           "python": sys.version.split()[0],
           "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED"),
           "build_set": BUILD_SET, "builds": BUILDS}
    files = list(scripts) + [f"{OUT}/common.py", f"{PRFBA}/prfba_objective.py",
                             f"{CODIF}/build_prfba_community_models.py",
                             f"{CODIF}/compartment_aware_gapfill.py",
                             f"{CODIF}/demo_compartment_scoring.py",
                             f"{CODIF}/load_genomes_to_modelseedpy.py"]
    rec["sha256"] = {os.path.relpath(f, DOC): sha256_file(f, 16) for f in files}
    return rec


def sample_rank(asv, salt="downstream_e3"):
    return hashlib.sha256(f"{salt}|{asv}".encode()).hexdigest()


def n_stratum(n):
    return "1" if n == 1 else "2" if n == 2 else "3-4" if n <= 4 else "5-9" if n <= 9 else ">=10"


def dump(obj, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=False, default=str)
    os.replace(tmp, path)


def gpr_canonical(rule):
    """A GPR string as a DNF set of gene sets, so that 'a or b' == 'b or a'. MSBuilder writes
    GPRs by iterating sets of cobra objects, whose order follows object addresses and so
    differs between processes."""
    import ast
    rule = (rule or "").strip()
    if not rule:
        return frozenset()
    from cobra.core.gene import GPR
    try:
        tree = GPR.from_string(rule)
    except Exception:
        return ("unparsed", rule)

    def dnf(node):
        if isinstance(node, (ast.Module, ast.Expression)):
            body = node.body
            return dnf(body[0].value if isinstance(body, list) else body)
        if isinstance(node, ast.Expr):
            return dnf(node.value)
        if isinstance(node, ast.Name):
            return {frozenset([node.id])}
        if isinstance(node, ast.BoolOp):
            parts = [dnf(v) for v in node.values]
            if isinstance(node.op, ast.Or):
                return set().union(*parts)
            out = {frozenset()}
            for p in parts:
                out = {a | b for a in out for b in p}
            return out
        raise ValueError(f"unexpected GPR node {type(node).__name__}")
    try:
        return frozenset(dnf(tree))
    except Exception:
        return ("unparsed", rule)
