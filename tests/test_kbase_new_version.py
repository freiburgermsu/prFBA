"""Offline tests for upload_kbase_genomes.py --new-version (2026-11-15): a rebuild of ASVs whose
genomes are already in a KBase workspace is saved as a NEW VERSION of each same-named object, the
new GenomeSet references those exact versions (ws/obj/ver), the earlier objects' versions and the
earlier GenomeSet are left as they are, and an interrupted run resumes by the build stamp
(sha256 of {asv}.json in the latest version's metadata), not by name. No test touches the
network: upload_kbase_genomes.WS is replaced by FakeWS, an in-memory workspace that versions
objects by name as the Workspace service does. Run:
    ~/Documents/py_venv/bin/python -m pytest -q tests/test_kbase_new_version.py
"""
import hashlib
import json
import os
import re
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import pfunc_estimator as P  # noqa: E402
import upload_kbase_genomes as U  # noqa: E402

WSID = 7
VERSIONED = re.compile(r"^\d+/\d+/\d+$")
SOURCES = {"a1": ["7.1"], "a2": ["294.202", "294.71"], "a3": ["9.9"]}


# --------------------------------------------------------------------------- fake workspace
class FakeWS(U.WS):
    """In-memory Workspace: objects[name] = list of versions, each {'info': object_info,
    'data': ...}; saving to an existing name appends version n+1 with the same object id.
    `calls` records every (method, params); `fail_names` makes save_objects raise for a batch
    holding any of those names (as the service rejects the whole call)."""

    def __init__(self, token=None, url=None):
        self.objects, self.next_id, self.calls, self.fail_names = {}, 1, [], set()

    # helpers for the tests
    def put(self, name, typ, data, meta):
        return self._save_one({"name": name, "type": typ, "data": data, "meta": meta})

    def versions(self, name):
        return [v["info"][4] for v in self.objects.get(name, [])]

    def latest(self, name):
        return self.objects[name][-1]

    def saved_names(self):
        return [o["name"] for m, p in self.calls if m == "Workspace.save_objects" for o in p[0]["objects"]]

    def _save_one(self, o):
        vers = self.objects.setdefault(o["name"], [])
        if vers:
            objid, ver = vers[0]["info"][0], len(vers) + 1
        else:
            objid, ver = self.next_id, 1
            self.next_id += 1
        typ = o["type"] + ("-17.2" if o["type"] == U.GENOME_TYPE else "-2.1" if o["type"] == U.GENOMESET_TYPE else "-1.0")
        info = [objid, o["name"], typ, "2026-11-15T00:00:00+0000", ver, "me", WSID, "me:narrative_1",
                "chk", len(json.dumps(o.get("data"))), dict(o.get("meta") or {})]
        vers.append({"info": info, "data": json.loads(json.dumps(o.get("data"))),
                     "provenance": o.get("provenance")})
        return info

    def call(self, method, params, timeout=3600, retries=3):
        self.calls.append((method, params))
        p = params[0]
        if method == "Workspace.get_workspace_info":
            return [WSID, "me:narrative_1", "me", "", 0, "a", "n", "unlocked", {"narrative_nice_name": "test"}]
        if method == "Workspace.list_objects":
            rows = sorted((v[-1]["info"] for v in self.objects.values()), key=lambda i: i[0])
            if p.get("type"):
                rows = [i for i in rows if i[2].startswith(p["type"])]
            rows = [i for i in rows if i[0] >= p.get("minObjectID", 1)][:p.get("limit", 10000)]
            return [list(i[:10]) + ([dict(i[10])] if p.get("includeMetadata") else [None]) for i in rows]
        if method == "Workspace.get_object_info3":
            out = []
            for spec in p["objects"]:
                v = self.objects.get(spec["name"])
                out.append(list(v[-1]["info"]) if v else None)
            return {"infos": out, "paths": [None] * len(out)}
        if method == "Workspace.get_objects2":
            data = []
            for spec in p["objects"]:
                ws_, oid, ver = (int(x) for x in spec["ref"].split("/"))
                v = next(v for vs in self.objects.values() for v in vs
                         if v["info"][0] == oid and v["info"][4] == ver)
                d = v["data"]
                if spec.get("included") == ["/elements/*/ref"]:
                    d = {"elements": {k: {"ref": e["ref"]} for k, e in d["elements"].items()}}
                data.append({"info": v["info"], "data": d})
            return {"data": data}
        if method == "Workspace.save_objects":
            if any(o["name"] in self.fail_names for o in p["objects"]):
                raise RuntimeError("simulated save failure")
            return [self._save_one(o) for o in p["objects"]]
        raise AssertionError(f"unexpected call {method}")


# --------------------------------------------------------------------------- fixtures
def genome(asv, srcs):
    feats = [{"id": f"{asv}_{i}", "type": "gene", "functions": [f"Prod {i}"],
              "aliases": [["PATRIC_id", f"fig|{g}.peg.{i}"], ["PGFAM", "PGF_1"]], "probability": 0.5}
             for i, g in enumerate(srcs, 1)]
    return {"id": asv, "features": feats, "cdss": [], "non_coding_features": [], "mrnas": []}


def make_build(tmp_path, name="build_new", asvs=("a1", "a2")):
    """A stamped (BV-BRC PATRIC-only) build dir as the builder writes it."""
    gdir = tmp_path / name
    (gdir / "provenance").mkdir(parents=True)
    sel = tmp_path / f"{name}_sel.json"
    sel.write_text(json.dumps({"_meta": {"generated_from": "hits.json"},
                               **{a: {"selected": [{"genome_id": SOURCES[a][0], "organism": "X"}],
                                      "tier": 1, "rel_ab": 0.1} for a in asvs}}))
    man = {"_meta": {"selection": str(sel), "n_asvs": len(asvs) + 1, "abstain": 1,
                     "annotation": P.ANNOTATION_STAMP,
                     "selection_stamp": {"sha256": "ab" * 32}}}
    for a in asvs:
        (gdir / f"{a}.json").write_text(json.dumps(genome(a, SOURCES[a])))
        (gdir / "provenance" / f"{a}.json").write_text(json.dumps({"asv": a, "genomes": SOURCES[a]}))
        man[a] = {"status": "built", "n_features": len(SOURCES[a]), "genomes": SOURCES[a]}
    man["zz"] = {"status": "abstain"}
    (gdir / "synthetic_genome_manifest.json").write_text(json.dumps(man))
    return gdir


def sha(gdir, asv):
    return hashlib.sha256((gdir / f"{asv}.json").read_bytes()).hexdigest()


@pytest.fixture
def ws(monkeypatch, tmp_path):
    """FakeWS holding the EARLIER upload: Genomes a1, a2 at version 1 (metadata 'asv', no
    stamp), an unrelated 20-hit-style Genome under another name, and the earlier GenomeSet
    referencing a1/a2 version 1; plus that upload's kbase_upload_log.json."""
    fake = FakeWS()
    old = {a: fake.put(a, U.GENOME_TYPE, {"id": a, "old": True}, {"asv": a}) for a in ("a1", "a2")}
    fake.put("Methanobacterium.1", U.GENOME_TYPE, {"id": "m"}, {"asv": "a3"})
    refs = {a: U.upa_of(i) for a, i in old.items()}
    fake.put("old_set", U.GENOMESET_TYPE, {"description": "old", "elements": {a: {"ref": r} for a, r in refs.items()}},
             {"genome_dir": "/old/build"})
    prior = tmp_path / "prior_log.json"
    prior.write_text(json.dumps({"_meta": {}, **{a: {"name": a, "upa": r} for a, r in refs.items()}}))
    fake.prior_log, fake.old_refs = str(prior), refs
    monkeypatch.setattr(U, "WS", lambda token, url=None: fake)
    tok = tmp_path / "token"
    tok.write_text("tok")
    monkeypatch.setattr(U, "TOKEN_FILE", str(tok))
    return fake


def run(monkeypatch, gdir, *extra, prior=None, gs="new_set"):
    argv = ["upload_kbase_genomes.py", "--genome-dir", str(gdir), "--ws-id", str(WSID),
            "--new-version", "--gate-workers", "1", "--genomeset-name", gs, "--genomeset-description", "new"]
    if prior:
        argv += ["--prior-upload-log", prior]
    monkeypatch.setattr(sys, "argv", argv + list(extra))
    return U.main()


# --------------------------------------------------------------------------- the plan (pure)
def info(objid, name, typ, ver, meta):
    return [objid, name, typ, "t", ver, "me", WSID, "w", "c", 1, meta]


def test_plan_new_versions_decides_by_stamp_type_asv_and_prior_log():
    d = {a: f"sha_{a}" for a in ("a1", "a2", "a3", "a4", "a5", "a6", "a7")}
    g = U.GENOME_TYPE + "-17.2"
    infos = {"a1": info(11, "a1", g, 1, {"asv": "a1"}),                         # earlier build, no stamp
             "a2": info(12, "a2", g, 2, {"asv": "a2", U.STAMP_KEY: "sha_a2"}),  # already this build
             "a4": info(14, "a4", "KBaseFBA.FBAModel-1.0", 1, {"asv": "a4"}),   # not a Genome
             "a5": info(15, "a5", g, 1, {"asv": "other"}),                      # another ASV's object
             "a6": info(16, "a6", g, 3, {"asv": "a6", U.STAMP_KEY: "sha_old"}),  # not in prior log
             "a7": info(99, "a7", g, 1, {"asv": "a7", U.STAMP_KEY: "sha_a7"})}  # created by this mode
    prior = {f"{WSID}/11", f"{WSID}/12", f"{WSID}/14", f"{WSID}/15"}
    plan = U.plan_new_versions(list(d), d, infos, prior)
    assert plan["done"] == {"a2": f"{WSID}/12/2", "a7": f"{WSID}/99/1"}
    assert plan["todo"] == ["a1", "a3"]
    assert plan["prev"] == {"a1": f"{WSID}/11/1", "a3": None}
    assert plan["next_ver"] == {"a1": 2, "a3": 1}
    assert len(plan["refusals"]) == 3
    joined = "\n".join(plan["refusals"])
    assert "a4 -> 7/14/1: type KBaseFBA.FBAModel-1.0" in joined
    assert "a5 -> 7/15/1: metadata asv 'other'" in joined
    assert "a6 -> 7/16/3: not an object recorded by --prior-upload-log" in joined
    # without a prior log, a6 (stale stamp) is simply versioned again
    plan2 = U.plan_new_versions(["a6"], d, infos)
    assert plan2["todo"] == ["a6"] and plan2["next_ver"] == {"a6": 4}


def test_object_infos_pages_and_carries_metadata():
    fake = FakeWS()
    for i in range(5):
        fake.put(f"g{i}", U.GENOME_TYPE, {}, {"asv": f"g{i}"})
    fake.put("s", U.GENOMESET_TYPE, {"elements": {}}, {})
    got = fake.object_infos(WSID, None, page=2)
    assert sorted(got) == ["g0", "g1", "g2", "g3", "g4", "s"]
    assert got["g3"][10] == {"asv": "g3"}
    assert sorted(fake.object_infos(WSID, U.GENOME_TYPE, page=2)) == ["g0", "g1", "g2", "g3", "g4"]
    lists = [p[0] for m, p in fake.calls if m == "Workspace.list_objects"]
    assert all(p["includeMetadata"] == 1 for p in lists) and len(lists) >= 4


def test_gate_returns_file_digests_only_when_asked(tmp_path):
    gdir = make_build(tmp_path)
    man = json.load(open(gdir / "synthetic_genome_manifest.json"))
    assert U.require_own_patric_ids(str(gdir), man, workers=1, log=lambda m: None) == \
        {"n_genomes": 2, "n_checked": 3}
    out = U.require_own_patric_ids(str(gdir), man, workers=1, log=lambda m: None, return_digests=True)
    assert out["digests"] == {a: sha(gdir, a) for a in ("a1", "a2")}


# --------------------------------------------------------------------------- end to end
def test_new_version_saves_version_2_and_a_set_of_exact_versions(ws, monkeypatch, tmp_path):
    gdir = make_build(tmp_path)
    run(monkeypatch, gdir, prior=ws.prior_log)
    # each genome is version 2 of the SAME object; version 1 is kept, unchanged
    for a in ("a1", "a2"):
        assert ws.versions(a) == [1, 2]
        v1, v2 = ws.objects[a]
        assert v1["data"] == {"id": a, "old": True} and v1["info"][10] == {"asv": a}
        assert v2["info"][0] == v1["info"][0]
        assert v2["info"][10][U.STAMP_KEY] == sha(gdir, a)
        assert v2["info"][10]["previous_version"] == ws.old_refs[a]
        assert v2["info"][10]["function_annotation"] == "BV-BRC PATRIC (RASTtk) only"
        prov = v2["provenance"][0]
        assert prov["custom"][U.STAMP_KEY] == sha(gdir, a) and prov["custom"]["previous_version"] == ws.old_refs[a]
        assert prov["method_params"][0][U.STAMP_KEY] == sha(gdir, a)
    # nothing else received a version; the earlier set and the 20-hit-style object are untouched
    assert ws.versions("old_set") == [1] and ws.versions("Methanobacterium.1") == [1]
    assert sorted(ws.saved_names()) == ["a1", "a2", "new_set"]
    # the new set references the exact new versions
    new_set = ws.latest("new_set")
    refs = {k: e["ref"] for k, e in new_set["data"]["elements"].items()}
    assert refs == {a: f"{WSID}/{ws.objects[a][0]['info'][0]}/2" for a in ("a1", "a2")}
    assert all(VERSIONED.match(r) for r in refs.values())
    assert new_set["info"][10]["upload_mode"] == "new-version"
    assert new_set["info"][10]["genome_dir"] == str(gdir)
    assert sorted(new_set["provenance"][0]["input_ws_objects"]) == sorted(refs.values())
    old_refs = {k: e["ref"] for k, e in ws.latest("old_set")["data"]["elements"].items()}
    assert old_refs == ws.old_refs                                  # still version 1
    log = json.load(open(gdir / "kbase_upload_log.json"))
    assert {a: log[a]["upa"] for a in ("a1", "a2")} == refs
    assert all(log[a]["verified"] and log[a]["previous_version"] == ws.old_refs[a] for a in ("a1", "a2"))
    assert log["_genomeset"]["upa"] == U.upa_of(new_set["info"]) and log["_meta"]["mode"] == "new-version"


def test_complete_rerun_saves_nothing_and_reuses_the_set(ws, monkeypatch, tmp_path):
    gdir = make_build(tmp_path)
    run(monkeypatch, gdir, prior=ws.prior_log)
    n = len(ws.saved_names())
    run(monkeypatch, gdir, prior=ws.prior_log)
    assert len(ws.saved_names()) == n                               # no genome, no set saved
    assert ws.versions("a1") == [1, 2] and ws.versions("new_set") == [1]
    assert json.load(open(gdir / "kbase_upload_log.json"))["_genomeset"]["reused"] is True


def test_interrupted_run_resumes_by_stamp(ws, monkeypatch, tmp_path):
    """a2's save fails: exit non-zero, no GenomeSet; the re-run saves a2 only (a1's latest version
    carries this build's stamp, although a1 is also the name of the earlier object)."""
    gdir = make_build(tmp_path)
    ws.fail_names = {"a2"}
    with pytest.raises(SystemExit, match="INCOMPLETE"):
        run(monkeypatch, gdir, prior=ws.prior_log)
    assert ws.versions("a1") == [1, 2] and ws.versions("a2") == [1] and "new_set" not in ws.objects
    log = json.load(open(gdir / "kbase_upload_log.json"))
    assert "error" in log["a2"] and log["a1"]["upa"].endswith("/2")
    ws.fail_names = set()
    before = len(ws.saved_names())
    run(monkeypatch, gdir, prior=ws.prior_log)
    assert ws.saved_names()[before:] == ["a2", "new_set"]          # a1 not re-sent
    assert ws.versions("a1") == [1, 2] and ws.versions("a2") == [1, 2]
    refs = {k: e["ref"] for k, e in ws.latest("new_set")["data"]["elements"].items()}
    assert all(r.endswith("/2") for r in refs.values())


def test_lost_log_is_rebuilt_from_the_workspace(ws, monkeypatch, tmp_path):
    """The log is not trusted by name: deleting it (or a stale one) changes nothing."""
    gdir = make_build(tmp_path)
    run(monkeypatch, gdir, prior=ws.prior_log)
    (gdir / "kbase_upload_log.json").write_text(json.dumps({"a1": {"upa": "7/1/1", "adopted": True}}))
    n = len(ws.saved_names())
    run(monkeypatch, gdir, prior=ws.prior_log)
    assert len(ws.saved_names()) == n
    assert json.load(open(gdir / "kbase_upload_log.json"))["a1"]["upa"].endswith("/2")


def test_a_changed_build_is_saved_as_a_further_version(ws, monkeypatch, tmp_path):
    gdir = make_build(tmp_path)
    run(monkeypatch, gdir, prior=ws.prior_log)
    g = json.load(open(gdir / "a1.json"))
    g["features"][0]["functions"] = ["Prod changed"]
    (gdir / "a1.json").write_text(json.dumps(g))
    run(monkeypatch, gdir, prior=ws.prior_log)
    assert ws.versions("a1") == [1, 2, 3] and ws.versions("a2") == [1, 2]
    assert ws.latest("a1")["info"][10][U.STAMP_KEY] == sha(gdir, "a1")
    assert ws.versions("new_set") == [1, 2]                         # refs changed -> new set version
    assert ws.latest("new_set")["data"]["elements"]["a1"]["ref"].endswith("/3")


def test_a_new_asv_becomes_a_new_object(ws, monkeypatch, tmp_path):
    SOURCES["a3x"] = ["5.5"]
    try:
        gdir = make_build(tmp_path, name="build_b", asvs=("a1", "a2", "a3x"))
        run(monkeypatch, gdir, prior=ws.prior_log)
        assert ws.versions("a3x") == [1] and ws.latest("a3x")["info"][10]["previous_version"] == "none (new object)"
        n = len(ws.saved_names())
        run(monkeypatch, gdir, prior=ws.prior_log)                  # resume: not in the prior log, but stamped
        assert len(ws.saved_names()) == n
    finally:
        del SOURCES["a3x"]


@pytest.mark.parametrize("damage", ["other_asv_meta", "not_in_prior_log", "not_a_genome"])
def test_refuses_to_version_an_object_it_does_not_own(ws, monkeypatch, tmp_path, damage):
    gdir = make_build(tmp_path)
    if damage == "other_asv_meta":
        ws.latest("a2")["info"][10]["asv"] = "a3"
    elif damage == "not_in_prior_log":
        prior = json.load(open(ws.prior_log))
        del prior["a2"]
        open(ws.prior_log, "w").write(json.dumps(prior))
    else:
        ws.latest("a2")["info"][2] = "KBaseFBA.FBAModel-1.0"
    with pytest.raises(SystemExit, match="must not receive a new version") as ei:
        run(monkeypatch, gdir, prior=ws.prior_log)
    assert "a2 -> " in ei.value.code and "Nothing was uploaded" in ei.value.code
    assert ws.saved_names() == [] and not os.path.exists(gdir / "kbase_upload_log.json")


def test_refuses_an_existing_set_name_from_another_build(ws, monkeypatch, tmp_path):
    gdir = make_build(tmp_path)
    with pytest.raises(SystemExit, match="already holds 'old_set'"):
        run(monkeypatch, gdir, prior=ws.prior_log, gs="old_set")
    assert ws.saved_names() == []
    with pytest.raises(SystemExit, match="already holds 'Methanobacterium.1'"):
        run(monkeypatch, gdir, prior=ws.prior_log, gs="Methanobacterium.1")
    assert ws.saved_names() == []


def test_dry_run_reads_and_writes_nothing(ws, monkeypatch, tmp_path, capsys):
    gdir = make_build(tmp_path)
    before = sorted(os.listdir(gdir))
    run(monkeypatch, gdir, "--dry-run", prior=ws.prior_log)
    out = capsys.readouterr().out
    assert ws.saved_names() == [] and sorted(os.listdir(gdir)) == before
    assert {m for m, p in ws.calls} <= {"Workspace.get_workspace_info", "Workspace.list_objects",
                                        "Workspace.get_object_info3"}
    assert f"version 2 of {ws.old_refs['a1'].rsplit('/', 1)[0]} (supersedes {ws.old_refs['a1']})" in out
    assert "[dry] would save 2 genomes" in out and "versions to be created: {2: 2}" in out
    assert "[dry] GenomeSet 'new_set': new object, 2 elements" in out


def test_file_changed_after_the_gate_is_not_uploaded(ws, monkeypatch, tmp_path):
    gdir = make_build(tmp_path)
    real = U.require_own_patric_ids

    def gate(*a, **k):
        out = real(*a, **k)
        out["digests"]["a1"] = "0" * 64
        return out
    monkeypatch.setattr(U, "require_own_patric_ids", gate)
    with pytest.raises(SystemExit, match="INCOMPLETE"):
        run(monkeypatch, gdir, prior=ws.prior_log)
    assert ws.versions("a1") == [1] and "new_set" not in ws.objects
    assert "changed after the feature gate" in json.load(open(gdir / "kbase_upload_log.json"))["a1"]["error"]


def test_limit_saves_part_and_writes_no_set(ws, monkeypatch, tmp_path):
    gdir = make_build(tmp_path)
    run(monkeypatch, gdir, "--limit", "1", prior=ws.prior_log)     # exit 0, partial
    assert ws.versions("a1") == [1, 2] and ws.versions("a2") == [1] and "new_set" not in ws.objects
    run(monkeypatch, gdir, prior=ws.prior_log)
    assert ws.versions("a2") == [1, 2] and ws.versions("new_set") == [1]


def test_dataset_label_reaches_the_notes(ws, monkeypatch, tmp_path):
    gdir = make_build(tmp_path)
    run(monkeypatch, gdir, "--dataset-label", "EmilyKin", prior=ws.prior_log)
    assert "(EmilyKin)" in ws.latest("a1")["data"]["notes"]
    assert U.DEFAULT_DATASET in U.to_kbase_genome(genome("a1", ["7.1"]), asv="a1", sel={}, man={},
                                                  commit="c")["notes"]


def test_default_mode_still_adopts_by_name(ws, monkeypatch, tmp_path):
    """Without --new-version nothing changes: a same-named object is adopted (the hazard the new
    mode exists for), no genome is saved, and the set points at version 1."""
    gdir = make_build(tmp_path)
    monkeypatch.setattr(sys, "argv", ["upload_kbase_genomes.py", "--genome-dir", str(gdir), "--ws-id", str(WSID),
                                      "--gate-workers", "1", "--genomeset-name", "plain"])
    U.main()
    assert ws.saved_names() == ["plain"]
    refs = {k: e["ref"] for k, e in ws.latest("plain")["data"]["elements"].items()}
    assert refs == ws.old_refs
    assert "mode" not in json.load(open(gdir / "kbase_upload_log.json"))["_meta"]


def test_prior_log_needs_new_version(ws, monkeypatch, tmp_path):
    gdir = make_build(tmp_path)
    monkeypatch.setattr(sys, "argv", ["upload_kbase_genomes.py", "--genome-dir", str(gdir), "--ws-id", str(WSID),
                                      "--prior-upload-log", ws.prior_log])
    with pytest.raises(SystemExit) as ei:
        U.main()
    assert ei.value.code == 2 and ws.calls == []
