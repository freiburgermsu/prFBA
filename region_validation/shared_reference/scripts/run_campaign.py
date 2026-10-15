#!/usr/bin/env python
"""Full-scale campaign driver for SPEC 11.2 (stages 1-6) with EXTERNAL accounting.

    ~/Documents/py_venv/bin/python scripts/run_campaign.py --root <campaign root> [--dry-run]

WHY PYTHON AND NOT BASH.  The four things the reviewers asked for are all things bash
cannot do without inventing a parser: (a) the resumability predicate is "predictions.json
PARSES and its ``meta`` agrees with the job I was about to launch" -- ``run_picrust2_fold``
re-runs and REWRITES predictions.json on a bare re-invocation, so file existence alone is
not "done" and a stale run built for a different ``--select-modes`` must be detected from
inside the JSON; (b) the ledger needs sha256+size per output in the scorer's own
``file_fp`` shape; (c) the inventory assertion is a set-diff over derived keys printed in
``score_folds._ikey`` format; (d) per-stage concurrency differs per stage and the stage
must be DRAINED (not aborted) on a failure.  ``xargs -P`` + ``jq`` would re-implement all
of it in a shell that silently loses exit statuses through pipes.  This file is pure
stdlib and runs under the only interpreter (py_venv); it never imports ``picrust2`` or
``modelseedpy`` and never imports the harness it is accounting for.

STAGES (SPEC 11.2; --stages selects a subset, default 1-6)

  1  plan_folds --phase panel                       (no --limit-genomes at full scale)
  2  plan_folds --phase prepass --k2 100 --gpu-cand 100 --force-backend cuda
  3  plan_folds --phase plan --k 10 --rank genus
  4  60 build_reference jobs   (10 folds x {bac,arc} x {0.25,0.5,1.0}), 18 concurrent
  5  60 run_picrust2_fold --routing forced  +  10 --routing nsti (d100 only)
  6  30 run_prfba_fold  +  10 --forced-answer (d100, --hits-from the base fold run)

MEASURED CORRECTION TO SPEC 11.2 (stage 2).  The spec calls the prepass "minutes on the
5070 Ti" with ``--k2 400``.  97% of aligner wall is the Biopython rescore stage and its
cost is ``n_query x gpu_cand``, so ``--k2 400`` is ~36M rescored pairs ~= 9-11 h.  The
prepass only needs hits at identity >= 0.987 (the tau that builds the edge graph), so this
driver runs ``--k2 100 --gpu-cand 100`` (~1-2 h) and records the deviation in
``<root>/campaign_inventory.json:deviations``.  The prFBA FOLD runs keep frozen
``--k2 200 --gpu-cand 400`` (SPEC 1 / D14) -- the correction is about the prepass only.

CONTRACTS (verified against the runners, not assumed)

* ``build_reference.py`` composes ``<outdir>/f{NN}_{bac|arc}_d{ddd}`` from
  --fold/--domain/--density, so --outdir is ``<root>/refs`` (NOT the campaign root) --
  ``refdir = os.path.join(args.outdir, ref_name(...))``, default ``SR/data/refs``.
* ``run_picrust2_fold.py --outdir`` IS the run directory.  A ``--routing nsti`` run names
  no domain (``meta.domain == "both"``), and ``score_folds._P2_NAME`` (``^f(\\d+)_(bac|arc)_d(\\d+)$``)
  raises ``declared_scope_contradicts_directory_name`` for any forced-looking name whose
  meta disagrees -- so the nsti run directory is ``f{NN}_nsti``, which that regex must not
  match, and both --ref-dir are passed explicitly rather than relying on sibling derivation.
* ``run_prfba_fold.py --outdir`` is the run directory and its basename MUST end in
  ``_forced`` when ``--forced-answer`` is passed (run_prfba_fold.py:336).  ``--hits-from``
  adopts the base fold run's alignment only when ``align_stamp.json`` is byte-identical,
  so every stamp knob (--k2/--gpu-cand/--prefilter-k/--min-ref-len/--force-backend and the
  reference DB) is passed IDENTICALLY to the base and forced runs; a stamp mismatch is not
  an error, it silently re-aligns (7x), so the ledger records ``hits_adopted`` read back
  out of the job log.
* ``--select-modes`` is passed explicitly to every prFBA run, base and forced: the scorer
  fails the inventory with ``select_modes_differ_across_runs`` if one run's set differs.
* Stage 4 is serialised for ONE warm-up build (the cheapest: lowest density, arc, first
  fold) before the pool starts, because ``build_reference`` copies its probe FASTA to the
  shared ``<refs>/probe_queries.fna`` under a non-atomic ``exists() -> copyfile`` guard;
  60 concurrent builds race there.  Nothing reads that file, but a shredded artefact in
  refs/ is not something to leave for a reviewer to find.  The warm-up also fails fast.

ACCOUNTING

* ``<root>/ledger.jsonl`` -- one JSON object appended (fsync'd, lock-serialised) as each
  job finishes, whether it RAN, was SKIPPED or FAILED: job key, resolved argv, exit
  status, wall seconds, and {size, sha256} of its primary outputs in the scorer's own
  ``file_fp`` shape, so ``accounting.completeness.inventory.fingerprints`` can be diffed
  against it key by key.
* ``<root>/campaign_inventory.json`` -- the end-of-campaign assertion: the expected run
  set DERIVED from the scale axes (never from disk), the discovered set read off disk
  independently of the ledger, the diff in both directions, duplicate/incomplete run
  directories, the 60 reference builds, and a LEDGER-vs-DISK re-fingerprint that catches
  an output rewritten after its ledger row was written.  Mismatch -> exit 1.

Exit: 0 ok, 1 inventory mismatch, 2 usage, 3 a stage lost jobs, 130 interrupted.
"""
from __future__ import annotations

import argparse
import collections
import datetime as _dt
import fnmatch
import hashlib
import json
import os
import re
import shlex
import signal
import socket
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)                      # region_validation/shared_reference
RV = os.path.dirname(SR)
REPO = os.path.dirname(RV)
PY = sys.executable
PLAN_FOLDS = os.path.join(HERE, "plan_folds.py")
BUILD_REF = os.path.join(HERE, "build_reference.py")
RUN_P2 = os.path.join(HERE, "run_picrust2_fold.py")
RUN_PF = os.path.join(HERE, "run_prfba_fold.py")
EXPECTED_PY = os.path.expanduser("~/Documents/py_venv/bin/python")

PRIMARY_DENSITY = 1.0
CONTROL_FAMILY = re.compile(r"^(incself|plumbing|recipe)")
_P2_NAME = re.compile(r"^f(\d+)_(bac|arc)_d(\d+)$")   # score_folds.py:766, verbatim

EXIT_OK, EXIT_INVENTORY, EXIT_USAGE, EXIT_STAGE, EXIT_INTR = 0, 1, 2, 3, 130


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def now_iso() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


def log(msg: str) -> None:
    print(f"[campaign {_dt.datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def dkey(d) -> str:
    """score_folds.dkey: 1.0 -> '1.0', 0.5 -> '0.5', 0.25 -> '0.25'."""
    return str(float(d))


def ddd(d) -> str:
    return "d%03d" % int(round(float(d) * 100))


def ikey(arm, fold, dk, domain="", extra="") -> str:
    """score_folds._ikey, verbatim, so the two inventories diff literally."""
    return "%s|fold=%s|d=%s|domain=%s|%s" % (arm, fold, dk, domain or "-", extra)


def file_fp(path):
    """score_folds.file_fp: {size, sha256} or None."""
    try:
        h, n = hashlib.sha256(), 0
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
                n += len(chunk)
        return {"size": n, "sha256": h.hexdigest()}
    except OSError:
        return None


def jload(path):
    with open(path) as fh:
        return json.load(fh)


def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"],
                              capture_output=True, text=True,
                              timeout=20).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def rel(path, root) -> str:
    try:
        return os.path.relpath(path, root)
    except ValueError:
        return path


# --------------------------------------------------------------------------- #
# the job
# --------------------------------------------------------------------------- #
class Job:
    """One process.  ``done()`` is OUR predicate over the job's OWN outputs."""

    def __init__(self, stage, key, argv, outputs, log_name, *, force_flag="--force",
                 done=None, inventory_key=None, run_dir=None, note=""):
        self.stage = stage
        self.key = key
        self.argv = list(argv)
        self.outputs = list(outputs)          # absolute paths; [0] is the primary
        self.log_name = log_name
        self.force_flag = force_flag
        self._done = done or (lambda: all(os.path.exists(p) for p in self.outputs))
        self.inventory_key = inventory_key
        self.run_dir = run_dir
        self.note = note

    def done(self) -> bool:
        try:
            return bool(self._done())
        except Exception:
            return False


def meta_of(path):
    """``predictions.json`` -> (meta, n_pred) or (None, 0) if absent/corrupt."""
    try:
        blob = jload(path)
    except Exception:
        return None, 0
    if not isinstance(blob, dict):
        return None, 0
    return blob.get("meta") or {}, len(blob.get("pred") or {})


def p2_done(pred_path, *, fold, domain, density, routing, mode="fold"):
    def _f():
        meta, n = meta_of(pred_path)
        if meta is None or n == 0:
            return False
        return (meta.get("mode") == mode
                and meta.get("fold") == fold
                and str(meta.get("domain")) == domain
                and abs(float(meta.get("density", -1)) - density) < 1e-9
                and str(meta.get("routing") or "forced") == routing)
    return _f


def pf_done(pred_path, *, fold, density, forced, modes, k2, gpu_cand):
    def _f():
        meta, n = meta_of(pred_path)
        if meta is None or n == 0:
            return False
        return (meta.get("mode") == "fold"
                and meta.get("fold") == fold
                and abs(float(meta.get("density", -1)) - density) < 1e-9
                and bool(meta.get("forced_answer")) == bool(forced)
                and [str(m) for m in (meta.get("select_modes") or [])] == list(modes)
                and int(meta.get("k2", -1)) == int(k2)
                and int(meta.get("gpu_cand", -1)) == int(gpu_cand))
    return _f


def json_parses(path):
    def _f():
        try:
            jload(path)
            return True
        except Exception:
            return False
    return _f


# --------------------------------------------------------------------------- #
# the ledger
# --------------------------------------------------------------------------- #
class Ledger:
    """Append-only JSONL, fsync'd per row, serialised across the worker threads."""

    def __init__(self, path, root, enabled=True):
        self.path = path
        self.root = root
        self.enabled = enabled
        self._lock = threading.Lock()
        self.rows = []

    def append(self, row):
        row = dict(row)
        row.setdefault("ts", now_iso())
        with self._lock:
            self.rows.append(row)
            if not self.enabled:
                return row
            with open(self.path, "a") as fh:
                fh.write(json.dumps(row, sort_keys=True) + "\n")
                fh.flush()
                os.fsync(fh.fileno())
        return row

    def last_by_key(self):
        out = {}
        for r in self.rows:
            out[r["key"]] = r
        return out

    def load_existing(self):
        """Rows from previous invocations, so --dry-run and the inventory can see them."""
        if not os.path.exists(self.path):
            return []
        out = []
        with open(self.path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    out.append(json.loads(line))
                except Exception:
                    continue
        return out


# --------------------------------------------------------------------------- #
# process management
# --------------------------------------------------------------------------- #
_LIVE = {}
_LIVE_LOCK = threading.Lock()
_STOP = threading.Event()


def _on_signal(signum, _frame):
    _STOP.set()
    log(f"signal {signum}: terminating {len(_LIVE)} live child process(es); "
        f"the ledger will record them as killed")
    with _LIVE_LOCK:
        for p in list(_LIVE.values()):
            try:
                p.terminate()
            except Exception:
                pass


def execute(job, ledger, root, logdir, forced):
    """Run (or skip) one job and append exactly one ledger row.  Returns the row."""
    log_path = os.path.join(logdir, job.log_name)
    common = {"stage": job.stage, "key": job.key, "argv": job.argv,
              "log": rel(log_path, root), "inventory_key": job.inventory_key,
              "run_dir": rel(job.run_dir, root) if job.run_dir else None,
              "forced": bool(forced)}
    if not forced and job.done():
        return ledger.append(dict(common, status="skipped", rc=0, wall_s=0.0,
                                  outputs=output_fps(job, root),
                                  reason="output already satisfies the completion "
                                         "predicate; --force to rerun"))
    if _STOP.is_set():
        return ledger.append(dict(common, status="not_started", rc=None, wall_s=0.0,
                                  outputs=output_fps(job, root),
                                  reason="campaign was interrupted before this job "
                                         "started"))
    argv = list(job.argv)
    if forced and job.force_flag and job.force_flag not in argv:
        argv.append(job.force_flag)
    common["argv"] = argv
    os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)
    t0 = time.time()
    rc, killed = None, False
    with open(log_path, "wb") as fh:
        fh.write((f"# {now_iso()} {job.key}\n# cwd {os.getcwd()}\n# "
                  + " ".join(shlex.quote(a) for a in argv) + "\n").encode())
        fh.flush()
        try:
            p = subprocess.Popen(argv, stdout=fh, stderr=subprocess.STDOUT,
                                 cwd=REPO, close_fds=True)
        except OSError as exc:
            return ledger.append(dict(common, status="failed", rc=None,
                                      wall_s=round(time.time() - t0, 1),
                                      outputs=output_fps(job, root),
                                      reason=f"spawn failed: {exc}"))
        with _LIVE_LOCK:
            _LIVE[job.key] = p
        try:
            rc = p.wait()
        finally:
            with _LIVE_LOCK:
                _LIVE.pop(job.key, None)
    wall = round(time.time() - t0, 1)
    if rc is not None and rc < 0:
        killed = True
    status = "ok" if rc == 0 else ("killed" if killed else "failed")
    row = dict(common, status=status, rc=rc, wall_s=wall,
               outputs=output_fps(job, root),
               done_predicate_ok=job.done(),
               log_tail=tail(log_path, 12) if rc != 0 else None,
               hits_adopted=hits_adopted(log_path) if job.key.startswith("prfba:")
               else None)
    if rc == 0 and not row["done_predicate_ok"]:
        row["status"] = "failed"
        row["reason"] = ("exit 0 but the completion predicate still fails: the run "
                         "did not write the output this job is accountable for")
    return ledger.append(row)


def output_fps(job, root):
    return {rel(p, root): file_fp(p) for p in job.outputs}


def tail(path, n):
    try:
        with open(path, errors="replace") as fh:
            return [ln.rstrip("\n") for ln in fh.readlines()[-n:]]
    except OSError:
        return []


def hits_adopted(log_path):
    """True when --hits-from really adopted the donor alignment (a silent 7x otherwise)."""
    try:
        with open(log_path, errors="replace") as fh:
            txt = fh.read()
    except OSError:
        return None
    if "[hits-from]" not in txt:
        return None
    return "[hits-from] adopted the alignment" in txt


# --------------------------------------------------------------------------- #
# job construction
# --------------------------------------------------------------------------- #
def nsti_folds(a):
    """The folds that get a --routing nsti run.

    The routed arm compares both domain references by NSTI, so it is defined only when
    BOTH domains are in scope; a single-domain campaign (SPEC 11.1's smoke shape) has no
    routing decision to make.  build_jobs and expected_inventory both call this, so the
    job list and the assertion can never disagree about it.
    """
    folds = list(range(a.folds)) if a.fold_list is None else a.fold_list
    return folds if set(a.domains) == {"bac", "arc"} else []


def build_jobs(a):
    """-> {stage_number: [Job, ...]} for the declared scale axes (never from disk)."""
    root, refs = a.root, os.path.join(a.root, "refs")
    p2root, pfroot = os.path.join(root, "picrust2"), os.path.join(root, "prfba")
    folds = list(range(a.folds)) if a.fold_list is None else a.fold_list
    dens = a.densities
    doms = a.domains
    modes = a.select_modes.split(",")
    J = collections.OrderedDict((s, []) for s in range(1, 7))

    # ---- stage 1-3: plan_folds ------------------------------------------- #
    panel_out = os.path.join(root, "panel_eligible.json")
    argv = [PY, PLAN_FOLDS, "--phase", "panel", "--outdir", root,
            "--threads", str(a.panel_threads)]
    if a.limit_groups:
        argv += ["--limit-groups", str(a.limit_groups)]
    if a.limit_genomes:
        argv += ["--limit-genomes", str(a.limit_genomes)]
    J[1].append(Job(1, "plan:panel", argv,
                    [panel_out, os.path.join(root, "panel_bac.fna"),
                     os.path.join(root, "panel_arc.fna"),
                     os.path.join(root, "panel_md5_seq.json"),
                     os.path.join(root, "panel_md5_id.json"),
                     os.path.join(root, "panel_ec_universe.json"),
                     os.path.join(root, "panel_16s_queries.fna")],
                    "stage1_plan_panel.log", done=json_parses(panel_out)))

    nb = os.path.join(root, "neighbors.json")
    J[2].append(Job(2, "plan:prepass",
                    [PY, PLAN_FOLDS, "--phase", "prepass", "--outdir", root,
                     "--k2", str(a.prepass_k2), "--gpu-cand", str(a.prepass_gpu_cand),
                     "--force-backend", a.force_backend,
                     "--threads", str(a.panel_threads)],
                    [nb, os.path.join(root, "prepass_amplicon.json"),
                     os.path.join(root, "prepass_full16s.json"),
                     os.path.join(root, "panel_amplicon_queries.fna")],
                    "stage2_plan_prepass.log", done=json_parses(nb),
                    note="MEASURED CORRECTION to SPEC 11.2: k2/gpu-cand 100 not 400"))

    plan_name = ("fold_plan.json" if a.rank == "genus" else f"fold_plan_{a.rank}.json")
    plan_path = os.path.join(root, plan_name)
    J[3].append(Job(3, "plan:plan",
                    [PY, PLAN_FOLDS, "--phase", "plan", "--outdir", root,
                     "--k", str(a.k), "--rank", a.rank],
                    [plan_path], "stage3_plan_plan.log",
                    done=json_parses(plan_path)))

    # ---- stage 4: 60 reference builds ------------------------------------ #
    for f in folds:
        for dom in doms:
            for d in dens:
                name = f"f{f:02d}_{dom}_{ddd(d)}"
                man = os.path.join(refs, name, "manifest.json")
                J[4].append(Job(
                    4, f"build:{name}",
                    [PY, BUILD_REF, "--mode", "fold", "--fold", str(f),
                     "--domain", dom, "--density", repr(float(d)),
                     "--fold-plan", plan_path, "--panel", panel_out,
                     "--outdir", refs, "--threads", str(a.build_threads)],
                    [man], f"stage4_build_{name}.log", done=json_parses(man),
                    run_dir=os.path.join(refs, name)))

    # ---- stage 5: 60 forced PICRUSt2 + 10 nsti --------------------------- #
    for f in folds:
        for dom in doms:
            for d in dens:
                name = f"f{f:02d}_{dom}_{ddd(d)}"
                rundir = os.path.join(p2root, name)
                pred = os.path.join(rundir, "predictions.json")
                J[5].append(Job(
                    5, f"picrust2:{name}",
                    [PY, RUN_P2, "--mode", "fold", "--routing", "forced",
                     "--fold", str(f), "--domain", dom, "--density", repr(float(d)),
                     "--ref-dir", os.path.join(refs, name),
                     "--fold-plan", plan_path, "--panel", panel_out,
                     "--outdir", rundir, "--processes", str(a.p2_processes)],
                    [pred], f"stage5_picrust2_{name}.log",
                    done=p2_done(pred, fold=f, domain=dom, density=float(d),
                                 routing="forced"),
                    inventory_key=ikey("picrust2", f, dkey(d), dom, "forced"),
                    run_dir=rundir))
    for f in nsti_folds(a):
        name = f"f{f:02d}_nsti"
        rundir = os.path.join(p2root, name)
        pred = os.path.join(rundir, "predictions.json")
        refargs = []
        for dom in ("bac", "arc"):
            refargs += ["--ref-dir",
                        os.path.join(refs, f"f{f:02d}_{dom}_{ddd(PRIMARY_DENSITY)}")]
        J[5].append(Job(
            5, f"picrust2:{name}",
            [PY, RUN_P2, "--mode", "fold", "--routing", "nsti", "--fold", str(f),
             "--density", repr(PRIMARY_DENSITY)] + refargs
            + ["--fold-plan", plan_path, "--panel", panel_out, "--outdir", rundir,
               "--processes", str(a.p2_processes)],
            [pred], f"stage5_picrust2_{name}.log",
            done=p2_done(pred, fold=f, domain="both", density=PRIMARY_DENSITY,
                         routing="nsti"),
            inventory_key=ikey("picrust2", f, dkey(PRIMARY_DENSITY), "", "nsti"),
            run_dir=rundir,
            note="D9 secondary arm: d100 only; directory name must not match "
                 "score_folds._P2_NAME"))

    # ---- stage 6: 30 prFBA fold runs + 10 forced-answer ------------------ #
    def pf_argv(f, d, rundir, forced=False, hits_from=None):
        argv = [PY, RUN_PF, "--mode", "fold", "--fold", str(f),
                "--density", repr(float(d)),
                "--fold-plan", plan_path, "--panel", panel_out]
        for dom in a.domains:
            argv += ["--ref-manifest",
                     os.path.join(refs, f"f{f:02d}_{dom}_{ddd(d)}", "manifest.json")]
        argv += ["--outdir", rundir,
                 "--select-modes", a.select_modes,
                 "--k2", str(a.prfba_k2), "--gpu-cand", str(a.prfba_gpu_cand),
                 "--prefilter-k", str(a.prfba_prefilter_k),
                 "--force-backend", a.force_backend,
                 "--workers", str(a.prfba_workers)]
        if forced:
            argv += ["--forced-answer"]
        if hits_from:
            argv += ["--hits-from", hits_from]
        return argv

    for f in folds:
        for d in dens:
            name = f"f{f:02d}_{ddd(d)}"
            rundir = os.path.join(pfroot, name)
            pred = os.path.join(rundir, "predictions.json")
            J[6].append(Job(
                6, f"prfba:{name}", pf_argv(f, d, rundir),
                [pred], f"stage6_prfba_{name}.log",
                done=pf_done(pred, fold=f, density=float(d), forced=False,
                             modes=modes, k2=a.prfba_k2, gpu_cand=a.prfba_gpu_cand),
                inventory_key=ikey("prfba", f, dkey(d), "", "fold"), run_dir=rundir))
    for f in folds:
        base = os.path.join(pfroot, f"f{f:02d}_{ddd(PRIMARY_DENSITY)}")
        name = f"f{f:02d}_{ddd(PRIMARY_DENSITY)}_forced"
        rundir = os.path.join(pfroot, name)
        pred = os.path.join(rundir, "predictions.json")
        J[6].append(Job(
            6, f"prfba:{name}",
            pf_argv(f, PRIMARY_DENSITY, rundir, forced=True, hits_from=base),
            [pred], f"stage6_prfba_{name}.log",
            done=pf_done(pred, fold=f, density=PRIMARY_DENSITY, forced=True,
                         modes=modes, k2=a.prfba_k2, gpu_cand=a.prfba_gpu_cand),
            inventory_key=ikey("prfba", f, dkey(PRIMARY_DENSITY), "",
                               "forced_answer"),
            run_dir=rundir,
            note="basename must end in _forced (run_prfba_fold.py:336); reuses the "
                 "base run's alignment via --hits-from"))
    return J


# --------------------------------------------------------------------------- #
# expected inventory (DERIVED from the axes; never from disk)
# --------------------------------------------------------------------------- #
def expected_inventory(a):
    folds = list(range(a.folds)) if a.fold_list is None else a.fold_list
    runs, builds = {}, {}
    for f in folds:
        for d in a.densities:
            for dom in a.domains:
                nm = f"f{f:02d}_{dom}_{ddd(d)}"
                runs[ikey("picrust2", f, dkey(d), dom, "forced")] = f"picrust2/{nm}"
                builds[nm] = f"refs/{nm}"
            runs[ikey("prfba", f, dkey(d), "", "fold")] = f"prfba/f{f:02d}_{ddd(d)}"
        runs[ikey("prfba", f, dkey(PRIMARY_DENSITY), "", "forced_answer")] = \
            f"prfba/f{f:02d}_{ddd(PRIMARY_DENSITY)}_forced"
    for f in nsti_folds(a):
        runs[ikey("picrust2", f, dkey(PRIMARY_DENSITY), "", "nsti")] = \
            f"picrust2/f{f:02d}_nsti"
    return runs, builds


def discover_runs(root):
    """Read the run directories off DISK, independently of the ledger and the scorer.

    Mirrors score_folds.discover_picrust2/discover_prfba key derivation (from meta, not
    from the directory name) and run_kind's reference-family classification, so the two
    inventories are literally diffable -- but implements it here so this accounting does
    not depend on the component it is auditing.
    """
    found = collections.defaultdict(list)     # ikey -> [relative dir]
    incomplete, control, unreadable = [], [], []
    for arm, sub in (("picrust2", "picrust2"), ("prfba", "prfba")):
        base = os.path.join(root, sub)
        if not os.path.isdir(base):
            continue
        for name in sorted(os.listdir(base)):
            d = os.path.join(base, name)
            if not os.path.isdir(d):
                continue
            pred = os.path.join(d, "predictions.json")
            relpath = f"{sub}/{name}"
            if not os.path.exists(pred):
                incomplete.append({"dir": relpath,
                                   "why": "no predictions.json (a dead job; the "
                                          "scorer's discover_* skips it silently, so "
                                          "it reads as a MISSING run, not a broken one)"})
                continue
            meta, npred = meta_of(pred)
            if meta is None:
                unreadable.append({"dir": relpath,
                                   "why": "predictions.json does not parse"})
                continue
            rn = meta.get("ref_name") or name
            ref_names = [str(x) for x in (rn if isinstance(rn, (list, tuple)) else [rn])
                         if x]
            if any(CONTROL_FAMILY.match(x) for x in ref_names) or \
                    CONTROL_FAMILY.match(name):
                control.append({"dir": relpath, "ref_names": ref_names,
                                "mode": meta.get("mode")})
                continue
            fold = meta.get("fold")
            dk = dkey(meta.get("density") if meta.get("density") is not None
                      else PRIMARY_DENSITY)
            if arm == "picrust2":
                routing = str(meta.get("routing") or "forced")
                dom = {"bac": "bac", "arc": "arc", "Bacteria": "bac",
                       "Archaea": "arc"}.get(str(meta.get("domain")), "")
                k = ikey("picrust2", fold, dk, dom if routing == "forced" else "",
                         routing)
            else:
                k = ikey("prfba", fold, dk, "",
                         "forced_answer" if meta.get("forced_answer") else "fold")
            found[k].append(relpath)
            if npred == 0:
                incomplete.append({"dir": relpath,
                                   "why": "predictions.json holds 0 predictions"})
    return found, incomplete, control, unreadable


def assert_inventory(a, ledger_rows, write=True):
    root = a.root
    exp_runs, exp_builds = expected_inventory(a)
    found, incomplete, control, unreadable = discover_runs(root)

    missing = sorted(k for k in exp_runs if k not in found)
    unexpected = sorted(k for k in found if k not in exp_runs)
    duplicates = sorted(k for k, v in found.items() if len(v) > 1)
    # the directory a key was found in must be the one the campaign addressed
    misplaced = [{"key": k, "expected_dir": exp_runs[k], "found_in": found[k]}
                 for k in sorted(exp_runs)
                 if k in found and found[k] != [exp_runs[k]]]

    built, build_missing = {}, []
    for nm, reldir in sorted(exp_builds.items()):
        man = os.path.join(root, reldir, "manifest.json")
        fp = file_fp(man)
        if fp is None:
            build_missing.append(reldir + ("/manifest.json (directory exists but the "
                                           "build did not finish)"
                                           if os.path.isdir(os.path.join(root, reldir))
                                           else "/manifest.json (absent)"))
        else:
            built[nm] = fp
    extra_builds = []
    refs_dir = os.path.join(root, "refs")
    if os.path.isdir(refs_dir):
        for nm in sorted(os.listdir(refs_dir)):
            p = os.path.join(refs_dir, nm)
            if not os.path.isdir(p) or nm in exp_builds or CONTROL_FAMILY.match(nm):
                continue
            extra_builds.append(f"refs/{nm}")

    # ---- ledger <-> disk: an output rewritten after its row was written --- #
    last = {}
    for r in ledger_rows:
        last[r.get("key")] = r
    drift, unledgered = [], []
    for key, reldir in sorted(list(exp_runs.items())):
        jk = ("picrust2:" if key.startswith("picrust2") else "prfba:") + \
            os.path.basename(reldir)
        row = last.get(jk)
        if row is None:
            unledgered.append(jk)
            continue
        for relp, fp in (row.get("outputs") or {}).items():
            now = file_fp(os.path.join(root, relp))
            if fp is None and now is None:
                continue
            if fp != now:
                drift.append({"key": jk, "output": relp, "ledger": fp, "on_disk": now})
    for nm, reldir in sorted(exp_builds.items()):
        row = last.get(f"build:{nm}")
        if row is None:
            unledgered.append(f"build:{nm}")
            continue
        for relp, fp in (row.get("outputs") or {}).items():
            now = file_fp(os.path.join(root, relp))
            if fp != now and not (fp is None and now is None):
                drift.append({"key": f"build:{nm}", "output": relp, "ledger": fp,
                              "on_disk": now})

    reasons = []
    if missing:
        reasons.append("run_inventory_incomplete")
    if unexpected:
        reasons.append("unexpected_run_directory")
    if duplicates:
        reasons.append("duplicate_run_directories")
    if misplaced:
        reasons.append("run_in_unexpected_directory")
    if incomplete or unreadable:
        reasons.append("incomplete_run_directory")
    if build_missing:
        reasons.append("reference_builds_incomplete")
    if extra_builds:
        reasons.append("unexpected_reference_build")
    if drift:
        reasons.append("output_changed_after_its_ledger_row")
    if unledgered:
        reasons.append("output_without_a_ledger_row")

    inv = {
        "generated": now_iso(),
        "root": root,
        "git_commit": git_commit(),
        "ok": not reasons,
        "reasons": reasons,
        "scale": {"folds": (list(range(a.folds)) if a.fold_list is None
                            else a.fold_list),
                  "domains": a.domains, "densities": [dkey(d) for d in a.densities],
                  "select_modes": a.select_modes.split(","),
                  "primary_density": dkey(PRIMARY_DENSITY)},
        "derivation": "expected = folds x domains x densities (forced PICRUSt2) + "
                      "folds (nsti, d100) + folds x densities (prFBA) + folds "
                      "(prFBA forced-answer, d100); keys in score_folds._ikey format "
                      "so accounting.completeness.inventory.expected diffs literally. "
                      "At SPEC 11.2 scale that is 60 + 10 + 30 + 10 = 110 runs and 60 "
                      "reference builds.",
        "n_expected_runs": len(exp_runs),
        "n_discovered_runs": len(found),
        "n_run_directories": sum(len(v) for v in found.values()),
        "n_expected_builds": len(exp_builds),
        "n_builds_present": len(built),
        "expected_runs": exp_runs,
        "discovered_runs": {k: sorted(v) for k, v in sorted(found.items())},
        "missing_runs": missing,
        "unexpected_runs": [{"key": k, "dirs": found[k]} for k in unexpected],
        "duplicate_runs": [{"key": k, "dirs": sorted(found[k])} for k in duplicates],
        "misplaced_runs": misplaced,
        "incomplete_run_directories": incomplete,
        "unreadable_run_directories": unreadable,
        "out_of_scope_control_runs": control,
        "missing_reference_builds": build_missing,
        "unexpected_reference_builds": extra_builds,
        "build_fingerprints": built,
        "ledger_vs_disk_drift": drift,
        "outputs_without_a_ledger_row": sorted(unledgered),
        "ledger": rel(os.path.join(root, "ledger.jsonl"), root),
        "n_ledger_rows": len(ledger_rows),
        "deviations": DEVIATIONS,
    }
    if write:
        inv.setdefault("gates_exit_status", None)
        p = os.path.join(root, "campaign_inventory.json")
        with open(p, "w") as fh:
            json.dump(inv, fh, indent=1, sort_keys=True)
            fh.write("\n")
        inv["written_to"] = p
    return inv


def print_inventory(inv, gating=True):
    print("\n==================== INVENTORY ASSERTION ====================")
    if not gating:
        print("INFORMATIONAL ONLY (stage 6 was not selected, so the campaign is "
              "expected to be incomplete; this verdict does NOT set the exit code)")
    print(f"expected runs   : {inv['n_expected_runs']}  "
          f"(discovered keys {inv['n_discovered_runs']}, "
          f"run directories {inv['n_run_directories']})")
    print(f"expected builds : {inv['n_expected_builds']}  "
          f"(present {inv['n_builds_present']})")
    for label, items in (("MISSING runs", inv["missing_runs"]),
                         ("UNEXPECTED runs", inv["unexpected_runs"]),
                         ("DUPLICATE runs", inv["duplicate_runs"]),
                         ("MISPLACED runs", inv["misplaced_runs"]),
                         ("INCOMPLETE run dirs", inv["incomplete_run_directories"]),
                         ("UNREADABLE run dirs", inv["unreadable_run_directories"]),
                         ("MISSING builds", inv["missing_reference_builds"]),
                         ("UNEXPECTED builds", inv["unexpected_reference_builds"]),
                         ("LEDGER vs DISK drift", inv["ledger_vs_disk_drift"]),
                         ("outputs with no ledger row",
                          inv["outputs_without_a_ledger_row"])):
        if items:
            print(f"  {label} ({len(items)}):")
            for it in items[:20]:
                print(f"    - {it if isinstance(it, str) else json.dumps(it, sort_keys=True)}")
            if len(items) > 20:
                print(f"    ... {len(items) - 20} more")
    if inv["out_of_scope_control_runs"]:
        print(f"  control-family run dirs present (not this driver's scope, "
              f"not a defect): "
              f"{[c['dir'] for c in inv['out_of_scope_control_runs']]}")
    print(f"VERDICT: {'OK' if inv['ok'] else 'MISMATCH ' + ','.join(inv['reasons'])}"
          f"{'' if gating else '  [informational]'}")
    if inv.get("written_to"):
        print(f"written: {inv['written_to']}")
    print("=============================================================\n")


# --------------------------------------------------------------------------- #
# preflight: a (fold, domain) with no test genome of that domain
# --------------------------------------------------------------------------- #
def preflight_domains(a):
    plan_name = ("fold_plan.json" if a.rank == "genus"
                 else f"fold_plan_{a.rank}.json")
    plan_p = os.path.join(a.root, plan_name)
    panel_p = os.path.join(a.root, "panel_eligible.json")
    if not (os.path.exists(plan_p) and os.path.exists(panel_p)):
        return {"status": "deferred",
                "why": "no fold plan / panel on disk yet (stages 1-3 have not run); "
                       "the check reruns automatically before stage 5"}
    try:
        plan, panel = jload(plan_p), jload(panel_p)
    except Exception as exc:
        return {"status": "unreadable", "why": str(exc)}
    gen = panel.get("genomes") or {}
    short = {"Bacteria": "bac", "Archaea": "arc"}
    per = {}
    empty = []
    folds = list(range(a.folds)) if a.fold_list is None else a.fold_list
    for rec in plan.get("folds") or []:
        f = int(rec.get("fold"))
        if f not in folds:
            continue
        c = collections.Counter(
            short.get((gen.get(g) or {}).get("domain"), "other")
            for g in (rec.get("test_genomes") or []))
        per[str(f)] = dict(c)
        for dom in a.domains:
            if not c.get(dom):
                empty.append(f"fold {f} domain {dom}")
    return {"status": "ok" if not empty else "empty_domain_folds",
            "test_genomes_by_domain": per,
            "folds_with_no_test_genome_of_a_declared_domain": empty,
            "consequence": "run_picrust2_fold exits 2 ('no fold-N test md5 has truth "
                           "domain ...') for such a (fold, domain); SPEC 11.2 launches "
                           "all 60 unconditionally and declares --expect-domains "
                           "bac,arc, so this must be resolved before stage 5, not "
                           "silently skipped" if empty else "none"}


# --------------------------------------------------------------------------- #
DEVIATIONS = [
    "SPEC 11.2 prepass '--k2 400 ... minutes on the 5070 Ti' is wrong by measurement: "
    "97% of aligner wall is the Biopython rescore and its cost is n_query x gpu_cand, so "
    "k2 400 ~= 36M pairs ~= 9-11 h.  Stage 2 runs --k2 100 --gpu-cand 100 (~1-2 h); the "
    "prepass only needs hits at identity >= 0.987.  The prFBA FOLD runs keep frozen "
    "k2 200 / gpu_cand 400.",
    "build_reference.py --outdir is the REFS parent (refdir = outdir/f{NN}_{dom}_d{ddd}, "
    "default SR/data/refs), so the driver passes <root>/refs, not the campaign root.",
    "The --routing nsti run directory is <root>/picrust2/f{NN}_nsti: score_folds._P2_NAME "
    "(^f(\\d+)_(bac|arc)_d(\\d+)$) must NOT match it, or the run's domain-less meta reads "
    "as declared_scope_contradicts_directory_name.",
    "Stage 4 serialises ONE warm-up build before the 18-way pool: build_reference copies "
    "its probe FASTA to the shared <refs>/probe_queries.fna under a non-atomic "
    "exists()->copyfile guard, which 60 concurrent builds race on (nothing reads that "
    "file; the serialised warm-up also fails fast before 18 trees are spawned).",
    "Stages 7+ (score_folds, C1/C2 recipe and plumbing controls, phenotype arms) are NOT "
    "in this driver: SPEC 11.2 scopes them to the orchestrator after review.  The exact "
    "score_folds invocation is printed at the end.",
]


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="SPEC 11.2 campaign driver (stages 1-6) with external accounting",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--root", required=True,
                    help="campaign root; plan artefacts land here and refs/, picrust2/, "
                         "prfba/, logs/, ledger.jsonl are created under it.  REQUIRED "
                         "(no default) so nobody launches full scale by accident")
    ap.add_argument("--stages", default="1-6",
                    help="stage subset, e.g. '4-6' or '1,2,3'")
    ap.add_argument("--jobs", default="",
                    help="comma-separated job-key globs to restrict to "
                         "(e.g. 'picrust2:f01_arc_d100,build:f0*')")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the resolved job list and the expected inventory; "
                         "run nothing, write nothing")
    ap.add_argument("--force", action="store_true",
                    help="rerun the selected jobs even if complete (adds each runner's "
                         "own --force)")
    ap.add_argument("--keep-going", action="store_true",
                    help="do not stop at the first stage that lost jobs (default is "
                         "fail-fast AFTER draining the stage)")
    # scale axes
    ap.add_argument("--folds", type=int, default=10)
    ap.add_argument("--fold-list", default="",
                    help="explicit fold ids instead of range(--folds)")
    ap.add_argument("--densities", default="0.25,0.5,1.0")
    ap.add_argument("--domains", default="bac,arc")
    ap.add_argument("--k", type=int, default=10, help="plan_folds --phase plan --k")
    ap.add_argument("--rank", default="genus",
                    choices=["genus", "family", "order", "md5only"])
    ap.add_argument("--select-modes",
                    default="legacy,exact_tie,graduated,exact_tie_guards")
    # stage 1-2 knobs
    ap.add_argument("--limit-groups", type=int, default=0)
    ap.add_argument("--limit-genomes", type=int, default=0)
    ap.add_argument("--panel-threads", type=int, default=32)
    ap.add_argument("--prepass-k2", type=int, default=100)
    ap.add_argument("--prepass-gpu-cand", type=int, default=100)
    ap.add_argument("--force-backend", default="cuda",
                    choices=["auto", "cuda", "metal", "cpu"])
    # prFBA frozen aligner knobs (SPEC 1 / D14)
    ap.add_argument("--prfba-k2", type=int, default=200)
    ap.add_argument("--prfba-gpu-cand", type=int, default=400)
    ap.add_argument("--prfba-prefilter-k", type=int, default=1500)
    # concurrency, tuned to the measured profile
    ap.add_argument("--build-jobs", type=int, default=18,
                    help="FastTree is single-threaded -> many one-process jobs")
    ap.add_argument("--build-threads", type=int, default=3,
                    help="only hmmbuild/place_seqs use it; 18 x 3 ~= 54 of 64 cores")
    ap.add_argument("--p2-jobs", type=int, default=4)
    ap.add_argument("--p2-processes", type=int, default=14)
    ap.add_argument("--prfba-jobs", type=int, default=2,
                    help="CPU-bound on the Biopython rescore (and one GPU): few jobs, "
                         "each saturating cores")
    ap.add_argument("--prfba-workers", type=int, default=30)
    ap.add_argument("--prfba-forced-jobs", type=int, default=4,
                    help="the forced-answer runs adopt --hits-from, so they are "
                         "select-only and cheap")
    a = ap.parse_args(argv)

    a.root = os.path.abspath(os.path.expanduser(a.root))
    a.domains = [d for d in a.domains.replace(" ", "").split(",") if d]
    for d in a.domains:
        if d not in ("bac", "arc"):
            ap.error(f"--domains takes bac,arc; got {d}")
    try:
        a.densities = [float(x) for x in a.densities.replace(" ", "").split(",") if x]
    except ValueError:
        ap.error("--densities must be floats")
    if not any(abs(d - PRIMARY_DENSITY) < 1e-9 for d in a.densities):
        ap.error("--densities must include the primary density 1.0: the nsti and "
                 "forced-answer arms are d100-only (SPEC D9 / 11.2)")
    a.fold_list = ([int(x) for x in a.fold_list.replace(" ", "").split(",") if x]
                   or None) if a.fold_list else None
    stages = set()
    for part in a.stages.replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-", 1)
            stages |= set(range(int(lo), int(hi) + 1))
        else:
            stages.add(int(part))
    bad = sorted(s for s in stages if s not in range(1, 7))
    if bad:
        ap.error(f"--stages names stage(s) outside 1-6: {bad}")
    a.stage_set = stages
    a.job_globs = [g for g in a.jobs.replace(" ", "").split(",") if g]
    return a


def selected(job, a):
    if job.stage not in a.stage_set:
        return False
    if not a.job_globs:
        return True
    return any(fnmatch.fnmatch(job.key, g) for g in a.job_globs)


STAGE_CONC = {1: 1, 2: 1, 3: 1}
STAGE_LABEL = {
    1: "plan_folds --phase panel",
    2: "plan_folds --phase prepass (MEASURED CORRECTION: k2/gpu-cand 100)",
    3: "plan_folds --phase plan",
    4: "build_reference (10 folds x domains x densities)",
    5: "run_picrust2_fold (forced + nsti)",
    6: "run_prfba_fold (fold + forced-answer)",
}


def stage_groups(stage, jobs, a):
    """[(label, concurrency, [Job...])] -- the ordered sub-pools of one stage."""
    if stage in (1, 2, 3):
        return [(STAGE_LABEL[stage], 1, jobs)]
    if stage == 4:
        pool = sorted(jobs, key=lambda j: j.key)
        # The shared <refs>/probe_queries.fna is written under a non-atomic
        # exists()->copyfile guard, so exactly one build must finish before the pool
        # starts.  Once that file is there the guard never fires again and the whole
        # stage can go wide.
        warm = []
        if pool and not os.path.exists(os.path.join(a.root, "refs",
                                                    "probe_queries.fna")):
            warm = sorted(pool, key=lambda j: (            # cheapest build first
                float(j.argv[j.argv.index("--density") + 1]),
                0 if j.argv[j.argv.index("--domain") + 1] == "arc" else 1,
                j.key))[:1]
        rest = [j for j in pool if j not in warm]
        out = []
        if warm:
            out.append(("warm-up build (serialised: shared probe_queries.fna race)",
                        1, warm))
        if rest:
            out.append(("reference builds", max(1, a.build_jobs), rest))
        return out
    if stage == 5:
        forced = [j for j in jobs if not j.key.endswith("_nsti")]
        nsti = [j for j in jobs if j.key.endswith("_nsti")]
        out = []
        if forced:
            out.append(("PICRUSt2 --routing forced", max(1, a.p2_jobs), forced))
        if nsti:
            out.append(("PICRUSt2 --routing nsti (d100)", max(1, a.p2_jobs), nsti))
        return out
    if stage == 6:
        base = [j for j in jobs if not j.key.endswith("_forced")]
        fa = [j for j in jobs if j.key.endswith("_forced")]
        out = []
        if base:
            out.append(("prFBA fold runs", max(1, a.prfba_jobs), base))
        if fa:
            out.append(("prFBA --forced-answer (--hits-from the base run)",
                        max(1, a.prfba_forced_jobs), fa))
        return out
    return [(STAGE_LABEL[stage], 1, jobs)]


def stale_dependents(a, J):
    """Jobs that --force on an UPSTREAM job invalidates but does not rerun.

    A rebuilt reference or a re-aligned base prFBA run leaves its consumers on disk
    complete-looking and silently stale -- the exact shape of defect the inventory gate
    exists for, except that here the directory IS present and its meta DOES match.  The
    driver cannot know whether the orchestrator wants them rerun, so it names them.
    """
    if not a.force:
        return []
    sel = {j.key for stage in J for j in J[stage] if selected(j, a)}
    out = {}
    for key in sorted(sel):
        if key.startswith("build:"):
            f, dom, dd = key.split(":")[1].split("_")
            deps = [f"picrust2:{f}_{dom}_{dd}", f"prfba:{f}_{dd}"]
            if dd == ddd(PRIMARY_DENSITY):
                deps += [f"picrust2:{f}_nsti", f"prfba:{f}_{dd}_forced"]
        elif key.startswith("prfba:") and not key.endswith("_forced"):
            f, dd = key.split(":")[1].split("_")
            deps = [f"prfba:{f}_{dd}_forced"] if dd == ddd(PRIMARY_DENSITY) else []
        else:
            continue
        for d in deps:
            if d not in sel:
                out.setdefault(d, []).append(key)
    return [{"stale": k, "because_forced": v} for k, v in sorted(out.items())]


def warn_stale(a, J, ledger=None):
    stale = stale_dependents(a, J)
    if not stale:
        return stale
    print("[campaign] WARNING: --force reruns the selected jobs only.  These "
          "DOWNSTREAM runs stay on disk, keep passing the completion predicate and "
          "the inventory gate, and are now stale:", flush=True)
    for rec in stale:
        print(f"    {rec['stale']}  (rebuilt upstream: {', '.join(rec['because_forced'])})",
              flush=True)
    print("[campaign] rerun them too (--jobs '<key>,...' --force) or accept that the "
          "campaign mixes pre- and post-force inputs.", flush=True)
    if ledger is not None:
        ledger.append({"stage": 0, "key": "campaign:stale_dependents_warning",
                       "status": "ok", "rc": 0, "wall_s": 0.0, "argv": [],
                       "outputs": {}, "stale_dependents": stale})
    return stale


def dry_run(a, J):
    total_runs = total_builds = 0
    print("=" * 78)
    print(f"DRY RUN -- campaign root {a.root}")
    print(f"interpreter {PY}")
    folds = list(range(a.folds)) if a.fold_list is None else a.fold_list
    print(f"axes: folds={folds} domains={a.domains} "
          f"densities={[dkey(d) for d in a.densities]} "
          f"select_modes={a.select_modes}")
    print("=" * 78)
    for stage in sorted(a.stage_set):
        jobs = [j for j in J[stage] if selected(j, a)]
        if not jobs:
            continue
        print(f"\n### STAGE {stage}: {STAGE_LABEL[stage]}  "
              f"({len(jobs)} job(s))")
        for label, conc, group in stage_groups(stage, jobs, a):
            print(f"  -- {label}: {len(group)} job(s), {conc} concurrent")
            for j in group:
                state = "DONE(skip)" if j.done() else "run"
                print(f"   [{state:10s}] {j.key}")
                print(f"       argv: {' '.join(shlex.quote(x) for x in j.argv)}")
                print(f"       outputs: {', '.join(rel(p, a.root) for p in j.outputs)}")
                if j.inventory_key:
                    print(f"       inventory_key: {j.inventory_key}")
                if j.note:
                    print(f"       note: {j.note}")
                total_runs += 1 if j.inventory_key else 0
                total_builds += 1 if j.key.startswith("build:") else 0
    exp_runs, exp_builds = expected_inventory(a)
    print("\n### EXPECTED INVENTORY (derived from the axes, not from disk)")
    print(f"  run directories: {len(exp_runs)}"
          f"   (forced "
          f"{sum(1 for k in exp_runs if k.endswith('|forced'))}"
          f" + nsti {sum(1 for k in exp_runs if k.endswith('|nsti'))}"
          f" + prfba fold {sum(1 for k in exp_runs if k.endswith('|fold'))}"
          f" + prfba forced-answer "
          f"{sum(1 for k in exp_runs if k.endswith('|forced_answer'))})")
    print(f"  reference builds: {len(exp_builds)}")
    for k in sorted(exp_runs):
        print(f"    {k}  ->  {exp_runs[k]}")
    print(f"  jobs in the selected stages that produce a run directory: {total_runs}")
    print(f"  jobs in the selected stages that produce a reference build: "
          f"{total_builds}")
    print("\n### PREFLIGHT (fold x domain test-genome coverage)")
    print(json.dumps(preflight_domains(a), indent=1, sort_keys=True))
    warn_stale(a, J)
    print("\n### DEVIATIONS RECORDED")
    for d in DEVIATIONS:
        print(f"  - {d}")
    print("\nnothing was executed (--dry-run).")


def main(argv=None) -> int:
    a = parse_args(argv)
    if os.path.realpath(PY) != os.path.realpath(EXPECTED_PY):
        print(f"[campaign] refusing to run under {PY}: the only interpreter is "
              f"{EXPECTED_PY}", file=sys.stderr)
        return EXIT_USAGE
    J = build_jobs(a)

    if a.dry_run:
        dry_run(a, J)
        return EXIT_OK

    os.makedirs(a.root, exist_ok=True)
    logdir = os.path.join(a.root, "logs")
    os.makedirs(logdir, exist_ok=True)
    ledger = Ledger(os.path.join(a.root, "ledger.jsonl"), a.root)
    previous = ledger.load_existing()
    signal.signal(signal.SIGINT, _on_signal)
    signal.signal(signal.SIGTERM, _on_signal)

    ledger.append({"stage": 0, "key": "campaign:start", "status": "ok", "rc": 0,
                   "wall_s": 0.0, "argv": [PY, os.path.abspath(__file__)] + list(
                       argv if argv is not None else sys.argv[1:]),
                   "outputs": {}, "host": socket.gethostname(), "pid": os.getpid(),
                   "git_commit": git_commit(),
                   "scale": {"folds": (list(range(a.folds)) if a.fold_list is None
                                       else a.fold_list),
                             "domains": a.domains,
                             "densities": [dkey(d) for d in a.densities],
                             "stages": sorted(a.stage_set),
                             "select_modes": a.select_modes.split(",")},
                   "n_previous_ledger_rows": len(previous),
                   "deviations": DEVIATIONS})

    warn_stale(a, J, ledger)
    rc_final = EXIT_OK
    dead_total = []
    for stage in sorted(a.stage_set):
        jobs = [j for j in J[stage] if selected(j, a)]
        if not jobs:
            log(f"stage {stage}: no selected jobs, skipping")
            continue
        if stage == 5:
            pf = preflight_domains(a)
            if pf["status"] == "empty_domain_folds":
                log("stage 5 PREFLIGHT: " + json.dumps(pf, sort_keys=True))
                log("a declared (fold, domain) has no test genome of that domain; "
                    "run_picrust2_fold will exit 2 for it.  Resolve the scope before "
                    "stage 5 (SPEC 11.2 launches all 60 unconditionally).")
                ledger.append({"stage": 5, "key": "campaign:preflight_domains",
                               "status": "failed", "rc": None, "wall_s": 0.0,
                               "argv": [], "outputs": {}, "preflight": pf})
                dead_total.append("campaign:preflight_domains")
                rc_final = EXIT_STAGE
                break
            ledger.append({"stage": 5, "key": "campaign:preflight_domains",
                           "status": "ok", "rc": 0, "wall_s": 0.0, "argv": [],
                           "outputs": {}, "preflight": pf})
        log(f"===== STAGE {stage}: {STAGE_LABEL[stage]} ({len(jobs)} job(s)) =====")
        dead, t_stage = [], time.time()
        for label, conc, group in stage_groups(stage, jobs, a):
            log(f"-- {label}: {len(group)} job(s), {conc} concurrent")
            done_n = [0]
            n = len(group)
            count_lock = threading.Lock()

            def one(job):
                row = execute(job, ledger, a.root, logdir, a.force)
                with count_lock:
                    done_n[0] += 1
                extra = ""
                if row.get("hits_adopted") is False:
                    extra = "  [hits-from did NOT adopt: it re-aligned]"
                log(f"   ({done_n[0]}/{n}) {row['status']:9s} rc={row['rc']} "
                    f"{row['wall_s']:>8.1f}s  {job.key}{extra}")
                if row["status"] not in ("ok", "skipped"):
                    for ln in (row.get("log_tail") or []):
                        log(f"        | {ln}")
                return row

            if conc <= 1:
                rows = [one(j) for j in group]
            else:
                with ThreadPoolExecutor(max_workers=conc) as ex:
                    rows = list(ex.map(one, group))
            bad = [r for r in rows if r["status"] not in ("ok", "skipped")]
            dead += [(r["key"], r.get("log")) for r in bad]
            if bad and label.startswith("warm-up"):
                log("the serialised warm-up build failed: aborting stage 4 before "
                    "spawning the pool")
                break
        log(f"===== STAGE {stage} finished in {time.time() - t_stage:.0f}s: "
            f"{len(jobs) - len(dead)}/{len(jobs)} ok/skipped, {len(dead)} lost =====")
        if dead:
            dead_total += [k for k, _ in dead]
            print(f"[campaign] STAGE {stage} LOST {len(dead)} JOB(S):", flush=True)
            for k, lg in dead:
                print(f"  - {k}   log: {os.path.join(a.root, lg) if lg else '(none)'}",
                      flush=True)
            rc_final = EXIT_STAGE
            if not a.keep_going:
                print("[campaign] fail-fast: not starting the next stage.  The ledger "
                      "is complete for every job that ran.", flush=True)
                break
        if _STOP.is_set():
            rc_final = EXIT_INTR
            break

    rows = previous + ledger.rows
    # The assertion is over the WHOLE campaign, so it can only gate the exit code when
    # the last stage was in scope.  `--stages 1-3` (the orchestrator's natural first
    # move: inspect the fold plan before spending hours on stage 4) would otherwise
    # exit 1 with 110 missing runs, which reads as a failure of the thing that just
    # succeeded.
    gating = 6 in a.stage_set and not a.job_globs
    inv = assert_inventory(a, rows)
    inv["gates_exit_status"] = gating
    inv["gates_exit_status_why"] = (
        "stage 6 selected and no --jobs filter" if gating else
        "stage 6 not selected or --jobs restricted the campaign: the inventory is "
        "reported but does not set the exit status")
    print_inventory(inv, gating)
    if dead_total:
        print(f"[campaign] jobs lost this invocation: {dead_total}")
    print("Next, outside this driver (SPEC 11.2, orchestrator scope):")
    print(f"  {EXPECTED_PY} {os.path.join(HERE, 'score_folds.py')} --mode fold \\\n"
          f"      --fold-plan {os.path.join(a.root, 'fold_plan.json')} "
          f"--panel {os.path.join(a.root, 'panel_eligible.json')} \\\n"
          f"      --picrust2 {os.path.join(a.root, 'picrust2')} "
          f"--prfba {os.path.join(a.root, 'prfba')} "
          f"--refs {os.path.join(a.root, 'refs')} \\\n"
          f"      --neighbors {os.path.join(a.root, 'neighbors.json')} "
          f"--cells {os.path.join(a.root, 'cells.csv.gz')} \\\n"
          f"      --summary {os.path.join(a.root, 'summary_genus.json')} "
          f"--controls {os.path.join(a.root, 'controls.json')} \\\n"
          f"      --rank-blocked genus --bootstrap 2000 "
          f"--expect-densities {','.join(dkey(d) for d in a.densities)} \\\n"
          f"      --expect-domains {','.join(a.domains)} --expect-nsti primary "
          f"--expect-forced-answer primary --force")
    if _STOP.is_set():
        return EXIT_INTR
    if rc_final != EXIT_OK:
        return rc_final
    if not gating:
        return EXIT_OK
    return EXIT_OK if inv["ok"] else EXIT_INVENTORY


if __name__ == "__main__":
    sys.exit(main())
