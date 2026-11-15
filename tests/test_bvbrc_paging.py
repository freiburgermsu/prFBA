"""Offline tests for the paged BV-BRC PATRIC CDS fetch (2026-11-13), the genome-row request's
User-Agent, and the uploader's two refusals: builds without the BV-BRC PATRIC-only annotation
stamp, and builds holding a feature whose PATRIC_id is not 'fig|<g>.' of one of its ASV's source
genomes (the feature gate).

Up to 2026-11-13 fast_build_genome and prefetch_gene_families.fetch_pgfams sent ONE request with
limit(25000) and no paging, so every genome over 25,000 CDS rows was cached silently truncated.
prefetch_gene_families.fetch_patric_cds_rows now pages limit(PAGE_SIZE,start) sorted by
feature_id and raises unless the rows provably equal the server's Content-Range total. No test
touches the network: the transport (prefetch_gene_families._get, or urlopen beneath it) is
replaced by a fake BV-BRC that serves pages from an in-memory row list. Run:
    ~/Documents/py_venv/bin/python -m pytest -q tests/test_bvbrc_paging.py
"""
import copy
import glob
import importlib.util
import json
import os
import re
import subprocess
import sys
import types
import urllib.error

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import build_synthetic_genomes as BSG  # noqa: E402
import pfunc_estimator as P  # noqa: E402
import prefetch_gene_families as PF  # noqa: E402

LIMIT_RE = re.compile(r"&limit\((\d+),(\d+)\)")
CAPAWARE = sorted(glob.glob(os.path.join(ROOT, "region_validation*", "scripts", "prefetch_pgfams.py")))


def cds_rows(n, gid="7.1", pg_every=3, extra_pg_at=None):
    """n BV-BRC PATRIC CDS rows sorted by feature_id (zero-padded so string order = index
    order), every pg_every-th with a PGFam PGF_<i % 50>; extra_pg_at: {index: pgfam}."""
    rows = []
    for i in range(n):
        r = {"feature_id": f"PATRIC.{gid}.CDS.{i:07d}", "patric_id": f"fig|{gid}.peg.{i + 1}",
             "product": f"Prod {i % 997}", "feature_type": "CDS"}
        if i % pg_every == 0:
            r["pgfam_id"] = f"PGF_{i % 50:08d}"
        if extra_pg_at and i in extra_pg_at:
            r["pgfam_id"] = extra_pg_at[i]
        rows.append(r)
    return rows


class FakeBVBRC:
    """Serves genome_feature pages from `rows` for limit(n,start) like BV-BRC: (page, total).
    Faults: max_page (the server silently returns at most this many rows per page),
    fail_at (raise URLError for that start offset), no_total (no Content-Range),
    total_shift (total reported +1 from the second page on), restart (pages after the first
    restart at 0: overlapping pages AND too many rows), unstable (pages after the first start one
    row early but still hold the right number of rows: one row repeated, one missing, as an
    unsorted query across replicas can do)."""

    def __init__(self, rows, *, max_page=None, fail_at=None, no_total=False, total_shift=False,
                 restart=False, unstable=False):
        self.rows, self.urls = rows, []
        self.max_page, self.fail_at, self.no_total = max_page, fail_at, no_total
        self.total_shift, self.restart, self.unstable = total_shift, restart, unstable

    def get(self, url, timeout=60, retries=3):
        self.urls.append(url)
        m = LIMIT_RE.search(url)
        assert m, f"unpaged query: {url}"
        n, start = int(m.group(1)), int(m.group(2))
        if self.fail_at is not None and start == self.fail_at:
            raise urllib.error.URLError(f"connection reset at start={start}")
        if self.max_page:
            n = min(n, self.max_page)
        total = len(self.rows) + (1 if self.total_shift and len(self.urls) > 1 else 0)
        if self.restart and start > 0:
            start = 0
        if self.unstable and start > 0:
            n, start = min(n, len(self.rows) - start), start - 1
        return copy.deepcopy(self.rows[start:start + n]), (None if self.no_total else total)


@pytest.fixture
def fake(monkeypatch):
    """Install a FakeBVBRC as prefetch_gene_families' transport; returns a factory."""
    def make(rows, **kw):
        f = FakeBVBRC(rows, **kw)
        monkeypatch.setattr(PF, "_get", f.get)
        return f
    return make


# --------------------------------------------------------------------------- paging
def test_pages_until_a_short_page_two_and_a_half_pages(fake):
    rows = cds_rows(62_500)                          # 2.5 pages of the real PAGE_SIZE (25,000)
    f = fake(rows)
    got = PF.fetch_patric_cds_rows("7.1", ["patric_id", "pgfam_id"])
    assert got == rows
    assert [LIMIT_RE.search(u).groups() for u in f.urls] == [
        ("25000", "0"), ("25000", "25000"), ("25000", "50000")]
    for u in f.urls:                                  # every page: PATRIC CDS only, stable order
        assert u.startswith(PF.API + "?eq(genome_id,7.1)&eq(feature_type,CDS)&eq(annotation,PATRIC)")
        assert "&select(feature_id,patric_id,pgfam_id)&sort(+feature_id)&limit(" in u
        assert u.endswith("&http_accept=application/json")


def test_page_size_is_a_module_constant_and_exact_multiples_and_empty_genomes(fake, monkeypatch):
    monkeypatch.setattr(PF, "PAGE_SIZE", 4)
    f = fake(cds_rows(10))                            # 2.5 pages of 4
    assert len(PF.fetch_patric_cds_rows("7.1", ["pgfam_id"])) == 10
    assert [LIMIT_RE.search(u).groups() for u in f.urls] == [("4", "0"), ("4", "4"), ("4", "8")]
    f = fake(cds_rows(8))                             # exactly 2 pages: stops at the total
    assert len(PF.fetch_patric_cds_rows("7.1", ["pgfam_id"])) == 8 and len(f.urls) == 2
    f = fake([])                                      # a genome with no PATRIC CDS
    assert PF.fetch_patric_cds_rows("7.1", ["pgfam_id"]) == [] and len(f.urls) == 1


def test_fetch_pgfams_covers_every_page(fake):
    # a PGFam that only the third (half) page carries, and one only the second page carries
    rows = cds_rows(62_500, extra_pg_at={30_001: "PGF_ONLY_PAGE2", 61_000: "PGF_ONLY_PAGE3"})
    fake(rows)
    got = PF.fetch_pgfams("7.1")
    assert got == sorted({r["pgfam_id"] for r in rows if r.get("pgfam_id")})
    assert "PGF_ONLY_PAGE2" in got and "PGF_ONLY_PAGE3" in got
    assert len(got) == 52
    # what the old single limit(25000) window would have returned: the first page only
    assert set(got) - {r["pgfam_id"] for r in rows[:25_000] if r.get("pgfam_id")} == \
        {"PGF_ONLY_PAGE2", "PGF_ONLY_PAGE3"}


def test_fast_build_genome_keeps_every_page(fake, monkeypatch, tmp_path):
    monkeypatch.setattr(PF, "PAGE_SIZE", 4)
    monkeypatch.setattr(BSG, "_api_json", lambda core, q, **kw: [
        {"genome_name": "Org 7.1", "gc_content": 50.0, "taxon_lineage_names": ["Bacteria"]}])
    rows = cds_rows(10)
    rows[9]["product"] = "Prod only on page 3"
    f = fake(rows)
    g = BSG.fast_build_genome("7.1")
    assert len(f.urls) == 3
    assert [ft["aliases"][0][1] for ft in g["features"]] == [r["patric_id"] for r in rows]
    assert g["features"][-1]["functions"] == ["Prod only on page 3"]
    assert all("feature_id" not in ft for ft in g["features"])   # the cached shape is unchanged
    assert g["features"][0]["aliases"] == [["PATRIC_id", "fig|7.1.peg.1"], ["PGFAM", "PGF_00000000"]]


# --------------------------------------------------------------------------- guards
def test_a_page_failing_mid_way_raises_and_nothing_is_returned(fake):
    fake(cds_rows(62_500), fail_at=25_000)            # page 2 fails after its retries
    with pytest.raises(urllib.error.URLError):
        PF.fetch_patric_cds_rows("7.1", ["pgfam_id"])
    with pytest.raises(urllib.error.URLError):
        PF.fetch_pgfams("7.1")


def test_a_failed_page_is_never_cached(fake, monkeypatch, tmp_path):
    """populate_cache: the genome whose second page fails stays ABSENT from the cache (retried
    next run), the complete one is cached; build_genome writes no cache file for it."""
    good, bad = cds_rows(30, gid="1.1"), cds_rows(30, gid="2.1")
    monkeypatch.setattr(PF, "PAGE_SIZE", 10)

    def get(url, timeout=60, retries=3):
        rows = good if "eq(genome_id,1.1)" in url else bad
        n, start = map(int, LIMIT_RE.search(url).groups())
        if rows is bad and start == 10:
            raise urllib.error.URLError("reset")
        return copy.deepcopy(rows[start:start + n]), len(rows)
    monkeypatch.setattr(PF, "_get", get)
    monkeypatch.setattr(PF, "candidate_genome_ids", lambda hits, knobs=None: {"1.1", "2.1"})
    cache_path = str(tmp_path / "gene_families.json")
    cache, cand = PF.populate_cache({}, cache_path, workers=2, verbose=False)
    on_disk = json.load(open(cache_path))
    assert set(on_disk) == {"1.1"} == set(cache)
    assert on_disk["1.1"] == sorted({r["pgfam_id"] for r in good if r.get("pgfam_id")})
    monkeypatch.setattr(BSG, "_api_json", lambda core, q, **kw: [{"genome_name": "Org"}])
    with pytest.raises(urllib.error.URLError):
        BSG.build_genome("2.1", str(tmp_path / "kbase_cache"), "fast")
    assert os.listdir(tmp_path / "kbase_cache") == []


@pytest.mark.parametrize("fault, message", [
    ({"max_page": 10_000}, "fetched 10000 rows but the server reports 62500"),
    ({"no_total": True}, "no Content-Range total"),
    ({"total_shift": True}, "total changed between pages"),
    ({"restart": True}, "fetched 75000 rows but the server reports 62500"),
    ({"unstable": True}, "1 repeated or missing feature_id over 62500 rows"),
])
def test_incomplete_fetches_raise(fake, fault, message):
    """A server that silently caps its pages, sends no total, changes its total, or serves
    overlapping pages: IncompleteFetch, never a partial set."""
    fake(cds_rows(62_500), **fault)
    with pytest.raises(PF.IncompleteFetch, match=re.escape(message)):
        PF.fetch_patric_cds_rows("7.1", ["pgfam_id"])
    with pytest.raises(PF.IncompleteFetch):
        PF.fetch_pgfams("7.1")


def test_missing_feature_id_and_oversized_page_raise(fake, monkeypatch):
    rows = cds_rows(5)
    del rows[2]["feature_id"]
    fake(rows)
    with pytest.raises(PF.IncompleteFetch, match="repeated or missing feature_id"):
        PF.fetch_patric_cds_rows("7.1", ["pgfam_id"])
    monkeypatch.setattr(PF, "_get", lambda url, **kw: (cds_rows(5), 5))
    with pytest.raises(PF.IncompleteFetch, match="page of 5 rows > limit 4"):
        PF.fetch_patric_cds_rows("7.1", ["pgfam_id"], page_size=4)


def test_transport_sends_the_user_agent_reads_content_range_and_retries(monkeypatch):
    calls = []

    class Resp:
        def __init__(self, body, cr):
            self.body, self.headers = body, ({"Content-Range": cr} if cr else {})

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self, *a):
            return json.dumps(self.body).encode()

    plan = [urllib.error.URLError("transient"), Resp([{"feature_id": "a"}], "items 0-0/1")]

    def urlopen(req, timeout=None):
        calls.append(req.get_header("User-agent"))
        x = plan.pop(0)
        if isinstance(x, Exception):
            raise x
        return x
    monkeypatch.setattr(PF.time, "sleep", lambda s: None)
    monkeypatch.setattr(PF.urllib.request, "urlopen", urlopen)
    url = PF.API + "?eq(genome_id,7.1)"
    assert PF._get(url) == ([{"feature_id": "a"}], 1)
    assert calls == [PF.USER_AGENT, PF.USER_AGENT]
    plan[:] = [urllib.error.URLError("down")] * 3
    with pytest.raises(urllib.error.URLError, match="down"):
        PF._get(url, retries=3)
    assert plan == [] and len(calls) == 5
    assert [PF._content_range_total(v) for v in
            ("items 0-24999/39609", "items 0-0/0", "items 40000-40000/39609", None, "", "items */*")] \
        == [39609, 0, 39609, None, None, None]


def test_api_json_sends_the_prefetch_user_agent_and_retries(monkeypatch):
    """build_synthetic_genomes._api_json (fast_build_genome's genome-row request) sends
    prefetch_gene_families.USER_AGENT on every attempt (up to 2026-11-13 it sent none)."""
    calls = []

    class Resp:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self, *a):
            return json.dumps([{"genome_name": "Org 7.1"}]).encode()

    plan = [urllib.error.URLError("transient"), Resp()]

    def urlopen(req, timeout=None):
        calls.append((req.get_header("User-agent"), req.full_url))
        x = plan.pop(0)
        if isinstance(x, Exception):
            raise x
        return x
    monkeypatch.setattr(BSG.time, "sleep", lambda s: None)
    monkeypatch.setattr(BSG.urllib.request, "urlopen", urlopen)
    assert BSG._api_json("genome", "eq(genome_id,7.1)&limit(1)") == [{"genome_name": "Org 7.1"}]
    url = f"{BSG.API}/genome/?eq(genome_id,7.1)&limit(1)&http_accept=application/json"
    assert calls == [(PF.USER_AGENT, url), (PF.USER_AGENT, url)]
    assert PF.USER_AGENT.startswith("prFBA/")


# --------------------------------------------------------------------------- region_validation copies
def _load_capaware(path, monkeypatch):
    """Import one region_validation*/scripts/prefetch_pgfams.py with a stub _common (the real one
    mkdirs its data dirs on import); sys.path changes are undone at teardown."""
    monkeypatch.setattr(sys, "path", list(sys.path))
    stub = types.ModuleType("_common")
    stub.PATHS = types.SimpleNamespace(selection_json="", selection_out="", pgfam_cache="")
    stub.done = lambda p: False
    monkeypatch.setitem(sys.modules, "_common", stub)
    name = "capaware_" + os.path.relpath(path, ROOT).split(os.sep)[0]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.PGF is PF
    return mod


def test_there_are_nine_capaware_copies_differing_only_in_save_every():
    assert len(CAPAWARE) == 9
    norm = {re.sub(r"^SAVE_EVERY = .*$", "", open(p).read(), flags=re.M) for p in CAPAWARE}
    assert len(norm) == 1


@pytest.mark.parametrize("path", CAPAWARE, ids=lambda p: os.path.relpath(p, ROOT).split(os.sep)[0])
def test_capaware_copies_page_patric_only_and_flag_genomes_over_one_page(path, fake, monkeypatch, tmp_path):
    mod = _load_capaware(path, monkeypatch)
    rows = cds_rows(62_500, extra_pg_at={61_000: "PGF_ONLY_PAGE3"})
    f = fake(rows)
    pgfams, n_rows, capped = mod.fetch_pgfams_capaware("7.1")
    assert pgfams == sorted({r["pgfam_id"] for r in rows if r.get("pgfam_id")}) and \
        "PGF_ONLY_PAGE3" in pgfams
    assert (n_rows, capped) == (62_500, True) and len(f.urls) == 3
    assert all("eq(annotation,PATRIC)" in u and "sort(+feature_id)" in u for u in f.urls)
    fake(cds_rows(25_000))                            # exactly one full page: complete, not flagged
    assert mod.fetch_pgfams_capaware("7.1")[1:] == (25_000, False)
    fake(cds_rows(100))
    assert mod.fetch_pgfams_capaware("7.1")[1:] == (100, False)
    # a page failing mid-way: run_pass caches nothing for that genome
    fake(cds_rows(62_500), fail_at=25_000)
    with pytest.raises(urllib.error.URLError):
        mod.fetch_pgfams_capaware("7.1")
    cache, cache_path = {}, str(tmp_path / "c.json")
    n, capped_ids = mod.run_pass("t", {"7.1"}, cache, cache_path)
    assert (n, capped_ids, cache) == (1, [], {}) and json.load(open(cache_path)) == {}


# --------------------------------------------------------------------------- uploader refusal
def feat(fid, *pids, coll="features"):
    """A built (merged) feature / CDS: the converter copies the source feature's aliases."""
    al = [["PATRIC_id", p] for p in pids] + [["PGFAM", "PGF_00000001"]]
    f = {"id": fid, "type": "CDS" if coll == "cdss" else "gene", "functions": [f"Prod {fid}"],
         "aliases": al}
    return f


def genome_obj(asv, own_pids):
    """A built genome whose every feature and CDS carries one of own_pids."""
    feats = [feat(f"{asv}_{i + 1}", p) for i, p in enumerate(own_pids)]
    cdss = [feat(f"{asv}_{i + 1}.CDS", p, coll="cdss") for i, p in enumerate(own_pids)]
    return {"id": asv, "features": feats, "cdss": cdss, "non_coding_features": [], "mrnas": []}


GOOD = {"a1": (["7.1"], ["fig|7.1.peg.1", "fig|7.1.rna.2"]),
        "a2": (["294.202", "294.71"], ["fig|294.202.peg.1", "fig|294.71.peg.5"])}


def _build_dir(tmp_path, meta_extra, asvs=None):
    """A build directory as build_for_asv writes it: manifest (entry 'genomes'), {asv}.json and
    provenance/{asv}.json (whose 'genomes' the feature gate reads). asvs: {asv: (source genome
    ids, [PATRIC_id per feature]) or (ids, genome dict)}; default GOOD."""
    asvs = GOOD if asvs is None else asvs
    gdir = tmp_path / "build"
    gdir.mkdir()
    (gdir / "provenance").mkdir()
    sel = tmp_path / "sel.json"
    sel.write_text(json.dumps({"_meta": {"generated_from": "h"},
                               **{a: {"selected": []} for a in asvs}}))
    man = {"_meta": {"selection": str(sel), "n_asvs": len(asvs), **meta_extra}}
    for a, (srcs, body) in asvs.items():
        g = body if isinstance(body, dict) else genome_obj(a, body)
        man[a] = {"status": "built", "n_features": len(g.get("features") or []), "genomes": srcs}
        (gdir / f"{a}.json").write_text(json.dumps(g))
        (gdir / "provenance" / f"{a}.json").write_text(json.dumps({"asv": a, "genomes": srcs}))
    (gdir / "synthetic_genome_manifest.json").write_text(json.dumps(man))
    return gdir


class _Reached(Exception):
    pass


@pytest.fixture
def uploader(monkeypatch, tmp_path):
    """upload_kbase_genomes with every way out to KBase replaced: the token file does not
    exist, and constructing the workspace client or posting raises _Reached."""
    U = pytest.importorskip("upload_kbase_genomes")
    monkeypatch.setattr(U, "TOKEN_FILE", str(tmp_path / "no_such_token"))

    def no_ws(*a, **k):
        raise _Reached("workspace client constructed")
    monkeypatch.setattr(U, "WS", no_ws)
    monkeypatch.setattr(U.requests, "post", no_ws)
    return U


@pytest.mark.parametrize("meta", [{}, {"annotation": None},
                                  {"annotation": "BV-BRC genome_feature CDS (mixed)"}],
                         ids=["no_stamp", "null_stamp", "other_stamp"])
@pytest.mark.parametrize("dry", [False, True], ids=["upload", "dry_run"])
def test_uploader_refuses_builds_without_the_bvbrc_only_stamp(uploader, monkeypatch, tmp_path, meta, dry):
    gdir = _build_dir(tmp_path, meta)
    before = sorted(os.listdir(gdir))
    monkeypatch.setattr(sys, "argv", ["upload_kbase_genomes.py", "--genome-dir", str(gdir),
                                      "--ws-id", "1"] + (["--dry-run"] if dry else []))
    with pytest.raises(SystemExit) as ei:
        uploader.main()
    msg = ei.value.code
    assert isinstance(msg, str) and msg.startswith("[upload] REFUSED")   # exit status 1
    assert "not a BV-BRC PATRIC-only build" in msg and "Nothing was uploaded" in msg
    assert sorted(os.listdir(gdir)) == before                          # no kbase_upload_log.json


def test_uploader_lets_a_stamped_build_through_the_gate(uploader, monkeypatch, tmp_path):
    gdir = _build_dir(tmp_path, {"annotation": P.ANNOTATION_STAMP})
    monkeypatch.setattr(sys, "argv", ["upload_kbase_genomes.py", "--genome-dir", str(gdir),
                                      "--ws-id", "1", "--dry-run"])
    with pytest.raises(FileNotFoundError, match="no_such_token"):      # next step: the token
        uploader.main()
    uploader.require_bvbrc_only_build({"annotation": P.ANNOTATION_STAMP}, str(gdir))  # no exit


def test_uploader_refusal_exit_status_as_a_script(tmp_path):
    """The real process: exit status 1 and the refusal on stderr. HOME points at tmp_path, so
    even a broken gate could not read the user's KBase token."""
    gdir = _build_dir(tmp_path, {})
    env = {**os.environ, "HOME": str(tmp_path)}
    r = subprocess.run([sys.executable, os.path.join(ROOT, "upload_kbase_genomes.py"),
                        "--genome-dir", str(gdir), "--ws-id", "1", "--dry-run"],
                       capture_output=True, text=True, env=env, cwd=str(tmp_path))
    assert r.returncode == 1, (r.stdout, r.stderr)
    assert "[upload] REFUSED" in r.stderr and "[ws]" not in r.stdout
    assert not os.path.exists(gdir / "kbase_upload_log.json")


# --------------------------------------------------------------------------- uploader feature gate
STAMP = {"annotation": P.ANNOTATION_STAMP}


def _run_main(uploader, monkeypatch, gdir, extra=()):
    monkeypatch.setattr(sys, "argv", ["upload_kbase_genomes.py", "--genome-dir", str(gdir),
                                      "--ws-id", "1", *extra])
    return uploader.main()


def _bad_genome(kind):
    """a2's genome (sources 294.202, 294.71) with ONE offending CDS of the given kind."""
    g = genome_obj("a2", ["fig|294.202.peg.1", "fig|294.71.peg.5"])
    bad = {"foreign": feat("a2_9.CDS", "fig|294.7.peg.2281", coll="cdss"),
           "prefix_collision": feat("a2_9.CDS", "fig|294.2020.peg.1", coll="cdss"),
           "empty": feat("a2_9.CDS", "", coll="cdss"),
           "no_alias": {"id": "a2_9.CDS", "functions": ["Prod x"], "aliases": [["PGFAM", "PGF_1"]]},
           "own_plus_foreign": feat("a2_9.CDS", "fig|294.202.peg.3", "fig|294.7.peg.3",
                                    coll="cdss")}[kind]
    g["cdss"].append(bad)
    return g


@pytest.mark.parametrize("kind", ["foreign", "prefix_collision", "empty", "no_alias",
                                  "own_plus_foreign"])
@pytest.mark.parametrize("dry", [False, True], ids=["upload", "dry_run"])
def test_feature_gate_refuses_a_stamped_build_with_a_non_own_patric_id(uploader, monkeypatch,
                                                                         tmp_path, kind, dry):
    """A STAMPED build (it passes the annotation gate) with one CDS whose PATRIC_id is not
    'fig|<g>.' of a source genome of its ASV (another genome's id, a prefix-colliding one
    (294.2020 is not 294.202), an empty one, none, or an own id next to a foreign one) is
    refused, exit status 1, before the token is read; nothing is written. --limit does not
    narrow the gate: the offending ASV is the second one."""
    gdir = _build_dir(tmp_path, STAMP, {"a1": GOOD["a1"],
                                        "a2": (GOOD["a2"][0], _bad_genome(kind))})
    before = sorted(os.listdir(gdir))
    with pytest.raises(SystemExit) as ei:
        _run_main(uploader, monkeypatch, gdir, ["--limit", "1"] + (["--dry-run"] if dry else []))
    msg = ei.value.code
    assert isinstance(msg, str) and msg.startswith("[upload] REFUSED: 1 of 2 built genomes")
    assert "\n  a2: 1 of 5 features/CDSs (" in msg and "  a1:" not in msg
    assert "Nothing was uploaded" in msg and "cdss/a2_9.CDS" in msg
    want = {"foreign": "another genome's id fig|294.7. x1",
            "prefix_collision": "another genome's id fig|294.2020. x1",
            "empty": "empty id x1", "no_alias": "no PATRIC_id x1",
            "own_plus_foreign": "another genome's id fig|294.7. x1"}[kind]
    assert want in msg
    assert sorted(os.listdir(gdir)) == before                          # no kbase_upload_log.json


def test_feature_gate_passes_own_ids_of_every_source_genome(uploader, monkeypatch, tmp_path):
    """Own ids of any of the ASV's source genomes pass, in every collection (rna ids too); the
    next step is the token (FileNotFoundError), so the gate let the build through."""
    gdir = _build_dir(tmp_path, STAMP)
    g = json.load(open(gdir / "a1.json"))
    g["non_coding_features"] = [feat("a1_r", "fig|7.1.rna.9", coll="non_coding_features")]
    (gdir / "a1.json").write_text(json.dumps(g))
    with pytest.raises(FileNotFoundError, match="no_such_token"):
        _run_main(uploader, monkeypatch, gdir, ["--dry-run"])
    logged = []
    man = json.load(open(gdir / "synthetic_genome_manifest.json"))
    assert uploader.require_own_patric_ids(str(gdir), man, workers=1, log=logged.append) == \
        {"n_genomes": 2, "n_checked": 9}
    assert logged and "every PATRIC_id is a source genome's own" in logged[0]


def test_feature_gate_a_genome_of_another_asvs_sources_is_foreign(uploader, tmp_path):
    """The prefixes are per ASV: a1 holding a2's source genome's id is refused."""
    gdir = _build_dir(tmp_path, STAMP, {"a1": (["7.1"], ["fig|7.1.peg.1", "fig|294.202.peg.1"]),
                                        "a2": GOOD["a2"]})
    man = json.load(open(gdir / "synthetic_genome_manifest.json"))
    with pytest.raises(SystemExit, match=r"a1: 2 of 4 features/CDSs \(another genome's id "
                                         r"fig\|294\.202\. x2\)"):
        uploader.require_own_patric_ids(str(gdir), man, workers=1)


@pytest.mark.parametrize("damage, message", [
    ("no_provenance", "provenance/a2.json unreadable"),
    ("no_genomes", "provenance/a2.json has no source-genome list"),
    ("manifest_differs", "provenance/a2.json 'genomes' differs from the manifest entry's"),
    ("no_genome_file", "a2.json unreadable"),
    ("bad_json", "a2.json unreadable"),
])
def test_feature_gate_refuses_what_it_cannot_check(uploader, tmp_path, damage, message):
    gdir = _build_dir(tmp_path, STAMP)
    if damage == "no_provenance":
        os.remove(gdir / "provenance" / "a2.json")
    elif damage == "no_genomes":
        (gdir / "provenance" / "a2.json").write_text(json.dumps({"asv": "a2", "genomes": []}))
    elif damage == "manifest_differs":
        (gdir / "provenance" / "a2.json").write_text(json.dumps({"asv": "a2", "genomes": ["294.202"]}))
    elif damage == "no_genome_file":
        os.remove(gdir / "a2.json")
    else:
        (gdir / "a2.json").write_text('{"id": "a2", "features": [')
    man = json.load(open(gdir / "synthetic_genome_manifest.json"))
    with pytest.raises(SystemExit) as ei:
        uploader.require_own_patric_ids(str(gdir), man, workers=1)
    assert f"a2: cannot be checked: {message}" in ei.value.code
    assert ei.value.code.startswith("[upload] REFUSED: 1 of 2 built genomes")


def test_feature_gate_process_pool_matches_serial(uploader, tmp_path):
    """The pooled pass (>= 4 x workers genomes) gives the serial result, for a passing and a
    refused build; abstained / failed entries are not genomes and are skipped."""
    asvs = {f"g{i:02d}": (["7.1"], [f"fig|7.1.peg.{i}"]) for i in range(12)}
    asvs["g05"] = (["7.1"], ["fig|294.71.peg.2281"])
    gdir = _build_dir(tmp_path, STAMP, asvs)
    man = json.load(open(gdir / "synthetic_genome_manifest.json"))
    man["abst"] = {"status": "abstained"}
    man["fail"] = {"status": "failed"}
    msgs = []
    for w in (1, 2):
        with pytest.raises(SystemExit) as ei:
            uploader.require_own_patric_ids(str(gdir), man, workers=w)
        msgs.append(ei.value.code)
    assert msgs[0] == msgs[1] and msgs[0].startswith("[upload] REFUSED: 1 of 12 built genomes")
    assert "g05: 2 of 2 features/CDSs (another genome's id fig|294.71. x2)" in msgs[0]
    del man["g05"]
    assert uploader.require_own_patric_ids(str(gdir), man, workers=1, log=lambda m: None) == \
        uploader.require_own_patric_ids(str(gdir), man, workers=2, log=lambda m: None) == \
        {"n_genomes": 11, "n_checked": 22}


def test_feature_gate_exit_status_as_a_script(tmp_path):
    """The real process: exit status 1 and the feature-gate refusal on stderr, for a stamped
    build; HOME points at tmp_path, so even a broken gate could not read the user's token."""
    gdir = _build_dir(tmp_path, STAMP, {"a1": GOOD["a1"], "a2": (GOOD["a2"][0], _bad_genome("foreign"))})
    env = {**os.environ, "HOME": str(tmp_path)}
    r = subprocess.run([sys.executable, os.path.join(ROOT, "upload_kbase_genomes.py"),
                        "--genome-dir", str(gdir), "--ws-id", "1", "--dry-run"],
                       capture_output=True, text=True, env=env, cwd=str(tmp_path))
    assert r.returncode == 1, (r.stdout, r.stderr)
    assert "[upload] REFUSED: 1 of 2 built genomes" in r.stderr and "[ws]" not in r.stdout
    assert not os.path.exists(gdir / "kbase_upload_log.json")
