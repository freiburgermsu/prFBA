"""Offline tests for the E3 build integration: build_synthetic_genomes (build_for_asv hook,
prepare_pfunc, serial CLI), build_synthetic_genomes_parallel (--selection, --pfunc e3 through
a fork pool), reweight_build_e3.py (copy-on-write retro-apply) and the description strings.
Synthetic fixtures only (tiny genome cache, tiny PGFam cache, tiny clade catalogue); the real
codiffusion merge (bvbrc_to_kbase_genome.LocalGenomeConverter) runs, offline. Run:
    ~/Documents/py_venv/bin/python -m pytest -q tests/test_build_e3.py
"""
import copy
import hashlib
import json
import math
import os
import sys
import time

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import build_synthetic_genomes as BSG  # noqa: E402
import build_synthetic_genomes_parallel as BSGP  # noqa: E402
import pfunc_estimator as P  # noqa: E402
import reweight_build_e3 as R  # noqa: E402

BSG._patch_taxonomy_viz()

TAX = ("cellular organisms; Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; "
       "Enterobacteriaceae; Escherichia")


# --------------------------------------------------------------------------- fixtures
def patric_cds(gid, n, product, pg):
    """A fast-builder BV-BRC PATRIC CDS (patric_id 'fig|...'; pgfam None for the ~10% of
    PATRIC CDSs without a PGFam)."""
    al = [["PATRIC_id", f"fig|{gid}.peg.{n}"]]
    if pg:
        al.append(["PGFAM", pg])
    return {"type": "CDS", "functions": [product], "aliases": al, "pgfam_id": pg,
            "gene": None, "dna_sequence": "", "protein_translation": ""}


def refseq_cds(product):
    """A fast-builder RefSeq (NCBI PGAP) CDS as cached up to 3ce12d0: empty patric_id, no
    PGFam (what the unfiltered genome_feature query returned beside the PATRIC rows)."""
    return {"type": "CDS", "functions": [product], "aliases": [["PATRIC_id", ""]],
            "pgfam_id": None, "gene": None, "dna_sequence": "", "protein_translation": ""}


def src(gid, prods, refseq_first=(), refseq_last=()):
    """A cached source genome dict: prods = [(product, pgfam|None)] BV-BRC PATRIC CDSs, plus
    RefSeq CDSs (products) before and after them, as a mixed cache file holds them."""
    feats = [refseq_cds(p) for p in refseq_first]
    feats += [patric_cds(gid, i + 1, p, pg) for i, (p, pg) in enumerate(prods)]
    feats += [refseq_cds(p) for p in refseq_last]
    return {"id": gid, "scientific_name": f"Org {gid}", "taxonomy": TAX, "gc_content": 0.5,
            "features": feats}


# The genome cache is MIXED, as the real cache files written up to 3ce12d0 are: every build
# must use only the PATRIC CDSs (build_synthetic_genomes.bvbrc_only), so every expectation
# below is that of the PATRIC CDSs alone (PATRIC_SOURCES). The RefSeq CDSs would change the
# builds if they got in: 'Prod D' in 1.1 (RefSeq, listed first) would introduce asv_n2's
# 'Prod D' feature without its PGFAM alias and make its k/n 1.0; 'Prod B' in 2.1 and 3.1 and
# 'Prod F' in 2.1 repeat PATRIC products of other genomes (k/n 1.0 / 2/3); 'Prod R*' exist only
# as RefSeq CDSs (extra merged features).
SOURCES = {
    # 'Prod E' enters the merge from 1.1, whose PATRIC CDS has no PGFam for it; 2.1 carries
    # the same product string with PGF_E -> the product-match route.
    "1.1": src("1.1", [("Prod A", "PGF_A"), ("Prod B", None), ("Prod C", "PGF_C"),
                       ("Prod E", None)], refseq_first=["Prod D", "Prod R1"]),
    "2.1": src("2.1", [("Prod A", "PGF_A"), ("Prod C", "PGF_C2"), ("Prod D", "PGF_D"),
                       ("Prod E", "PGF_E")], refseq_last=["Prod B", "Prod F", "Prod R2"]),
    "3.1": src("3.1", [("Prod A", "PGF_A"), ("Prod D", "PGF_D"), ("Prod F", "PGF_F")],
               refseq_last=["Prod B", "Prod R3"]),
}
PATRIC_SOURCES = {g: BSG.bvbrc_only(s) for g, s in SOURCES.items()}
PGFAM_CACHE = {  # shared-cache format: genome_id -> sorted PGFam list
    "1.1": ["PGF_A", "PGF_C"], "2.1": ["PGF_A", "PGF_C2", "PGF_D", "PGF_E"],
    "3.1": ["PGF_A", "PGF_D", "PGF_F"], "4.1": [],          # empty list: treated as missing
    "11.1": ["PGF_A"], "12.1": ["PGF_D"], "13.1": ["PGF_X"],
}                                                         # 5.1: absent (fetch path)


def sel(gid, ident, species, rank):
    return {"rank": rank, "genome_id": gid, "taxon_id": int(gid.split(".")[0]),
            "organism": f"Org {gid}", "identity": ident, "align_score": ident * 500,
            "family": "Fam", "genus": None, "species": species, "role": "anchor",
            "confidence": "medium", "est_novel_fraction": 1.0, "est_genes": 100,
            "gene_set_empty": False, "equivalent_genomes": []}


def dec(rank, gid, ident, species):
    return {"rank": rank, "genome_id": gid, "organism": f"Org {gid}", "identity": ident,
            "align_score": ident * 500, "coverage": 1.0, "edlib_agreement": 0.0,
            "family": "Fam", "genus": None, "species": species, "disposition": "excluded",
            "stage": "identity_band", "reason": "test"}


def record(selected, decisions):
    return {"asv_len": 400, "best_identity": 0.99, "midas_taxonomy": None, "rel_ab": 0.01,
            "tier": "species", "consensus_family": "fam", "n_selected": len(selected),
            "flags": ["select_all:window=0.0"], "selected": selected, "provenance": {},
            "decisions": decisions}


DECS = [dec(1, "1.1", 0.99, "Sp one"), dec(2, "2.1", 0.99, "Sp two"),
        dec(3, "3.1", 0.97, "Sp three"), dec(4, "4.1", 0.96, "Sp four"),
        dec(5, "5.1", 0.95, "Sp five")]
AUDIT = {
    "asv_n2": record([sel("1.1", 0.99, "Sp one", 1), sel("2.1", 0.99, "Sp two", 2)], DECS),
    "asv_n1": record([sel("1.1", 0.99, "Sp one", 1)], DECS),
    "asv_n3": record([sel("1.1", 0.99, "Sp one", 1), sel("2.1", 0.99, "Sp two", 2),
                      sel("3.1", 0.97, "Sp three", 3)], DECS),
    "asv_abst": {**record([], []), "tier": "below", "flags": ["abstain_below_family"]},
}


@pytest.fixture
def world(tmp_path):
    cache = tmp_path / "genome_cache"
    cache.mkdir()
    for g, s in SOURCES.items():
        (cache / f"{g}.json").write_text(json.dumps(s))
    pg = tmp_path / "pgfam_cache.json"
    pg.write_text(json.dumps(PGFAM_CACHE))
    st = os.stat(pg)
    cat = tmp_path / "catalogue.json"
    cat.write_text(json.dumps({
        "_meta": {"source_cache_size": st.st_size,
                  "source_cache_mtime": time.strftime("%Y-%m-%dT%H:%M:%S%z",
                                                      time.localtime(st.st_mtime))},
        "genomes": {"1.1": [1, "Fam", "Ord"], "2.1": [2, "Fam", "Ord"], "3.1": [3, "Fam", "Ord"],
                    "11.1": [11, "Fam", "Ord"], "12.1": [12, "Fam", "Ord"],
                    "13.1": [13, "Fam", "Ord"]},
        "taxa": {}}))
    selp = tmp_path / "selection.json"
    selp.write_text(json.dumps({"_meta": {"gene_data": "bvbrc"}, **AUDIT}))
    return {"tmp": tmp_path, "cache": str(cache), "pgfam": str(pg), "catalogue": str(cat),
            "selection": str(selp)}


def e3_ctx(w, out_dir, records=None, **kw):
    ctx, meta = BSG.prepare_pfunc(records or AUDIT, str(out_dir), mode="e3",
                                  pgfam_cache=w["pgfam"], catalogue=w["catalogue"],
                                  fetch=False, taxdb=None, log=None, **kw)
    return ctx, meta


def build_one(out_dir, asv, pfunc, cache):
    os.makedirs(out_dir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(out_dir)
    try:
        for d in ("genome_objects", "taxonomies", "provenance"):
            os.makedirs(d, exist_ok=True)
        return BSG.build_for_asv(asv, copy.deepcopy(AUDIT[asv]), genome_cache_dir=cache,
                                 builder="fast", pfunc=pfunc)
    finally:
        os.chdir(cwd)


def rb(path):
    with open(path, "rb") as f:
        return f.read()


def jl(path):
    with open(path) as f:
        return json.load(f)


def files(d, asv):
    return {k: os.path.join(d, *([k] if k != "root" else []), f"{asv}.json")
            for k in ("root", "genome_objects", "provenance")}


def only_probabilities_differ(a_bytes, b_bytes):
    """True when b, with every feature/CDS probability set back to a's, dumps to a's bytes."""
    a, b = json.loads(a_bytes), json.loads(b_bytes)
    pa = {x["id"]: x.get("probability") for c in ("features", "cdss") for x in a[c]}
    for c in ("features", "cdss"):
        for x in b[c]:
            x["probability"] = pa[x["id"]]
    return json.dumps(b).encode() == a_bytes


def tree_digest(d):
    out = {}
    for r, _ds, fs in os.walk(d):
        for f in fs:
            p = os.path.join(r, f)
            out[os.path.relpath(p, d)] = (hashlib.sha256(rb(p)).hexdigest(), os.stat(p).st_mtime_ns)
    return out


# =========================================================================== build_for_asv
def test_kn_context_is_byte_identical_to_no_context(world):
    a, b = world["tmp"] / "a", world["tmp"] / "b"
    ea = build_one(a, "asv_n2", None, world["cache"])
    eb = build_one(b, "asv_n2", BSG.PfuncContext("kn"), world["cache"])
    assert ea == eb and "pfunc" not in eb
    for k in ("root", "genome_objects", "provenance"):
        assert rb(files(a, "asv_n2")[k]) == rb(files(b, "asv_n2")[k])
    assert not os.path.exists(b / BSG.PFUNC_DIR)


def test_e3_changes_only_probabilities_in_both_copies(world):
    kn_dir, e3_dir = world["tmp"] / "kn", world["tmp"] / "e3"
    build_one(kn_dir, "asv_n2", None, world["cache"])
    ctx, meta = e3_ctx(world, e3_dir)
    entry = build_one(e3_dir, "asv_n2", ctx, world["cache"])
    fk, fe = files(kn_dir, "asv_n2"), files(e3_dir, "asv_n2")
    for k in ("root", "genome_objects"):
        assert rb(fk[k]) != rb(fe[k])
        assert only_probabilities_differ(rb(fk[k]), rb(fe[k])), k
    root, go = jl(fe["root"]), jl(fe["genome_objects"])
    # the converter's copy keeps exactly its keys (annotation keys only in the root copy)
    assert all("selection_provenance" not in f for f in go["features"])
    assert all("selection_provenance" in f for f in root["features"])
    p_root = {f["id"]: f["probability"] for f in root["features"]}
    p_go = {f["id"]: f["probability"] for f in go["features"]}
    assert p_root == p_go
    for g in (root, go):
        cp = {c["id"]: c["probability"] for c in g["cdss"]}
        assert all(cp[f"{fid}.CDS"] == p for fid, p in p_root.items())
    # features: A(own), B(unresolved), C(own), E(product match via 2.1), D(own)
    fn = {f["functions"][0]: f["id"] for f in root["features"]}
    kn = {f["id"]: f["probability"] for f in jl(fk["root"])["features"]}
    assert p_root[fn["Prod B"]] == kn[fn["Prod B"]] == 0.5          # unresolved: keeps k/n
    # independent E3 for 'Prod A': reps 1.1, 2.1 (w 1), 3.1 (w e^-1); 4.1 empty, 5.1 absent;
    # family members = species 1, 2, 3, 11, 12, 13 -> PGF_A in 1.1, 2.1, 3.1, 11.1 -> 4/6
    w3 = math.exp(-(0.99 - 0.97) / P.TAU)
    den = 1 + 1 + w3 + P.ALPHA
    assert p_root[fn["Prod A"]] == pytest.approx((2 + w3 + P.ALPHA * 4 / 6) / den)
    assert p_root[fn["Prod E"]] == pytest.approx((1 + P.ALPHA * 1 / 6) / den)   # only 2.1 has PGF_E
    prov = jl(fe["provenance"])
    assert set(prov["probability_kn"]) == {fn[x] for x in ("Prod A", "Prod C", "Prod D", "Prod E")}
    assert prov["probability_kn"] == {fid: kn[fid] for fid in prov["probability_kn"]}
    assert prov["pfunc"]["applied"] and prov["pfunc"]["clade_rank"] == "family"
    assert prov["pfunc"]["feature_resolution"] == {"own_alias": 3, "product_match": 1,
                                                   "unresolved": 1, "features": 5}
    assert prov["pfunc"]["n_neighbours_missing"] == 1 and prov["pfunc"]["n_neighbours_empty"] == 1
    # the k/n provenance keys are unchanged
    pk = jl(fk["provenance"])
    assert {k: v for k, v in prov.items() if k not in ("probability_kn", "pfunc")} == pk
    assert entry["pfunc"]["applied"] and entry["pfunc"]["n_features_estimated"] == 4
    assert meta["pgfam_store"]["genomes_absent"] == 1 and meta["catalogue"]["current_vs_pgfam_cache"]


def test_e3_n1_and_failed_fetch_below_min_n_are_noops(world, monkeypatch):
    kn_dir, e3_dir = world["tmp"] / "kn", world["tmp"] / "e3"
    ctx, _ = e3_ctx(world, e3_dir)
    build_one(kn_dir, "asv_n1", None, world["cache"])
    entry = build_one(e3_dir, "asv_n1", ctx, world["cache"])
    for k in ("root", "genome_objects", "provenance"):
        assert rb(files(kn_dir, "asv_n1")[k]) == rb(files(e3_dir, "asv_n1")[k])
    assert entry["pfunc"]["applied"] is False and entry["pfunc"]["reason"].startswith("n_built 1")
    # n_selected 2, but 2.1 cannot be built -> n_built 1 -> k/n kept
    cache2 = world["tmp"] / "cache2"
    cache2.mkdir()
    (cache2 / "1.1.json").write_text(json.dumps(SOURCES["1.1"]))

    def boom(gid):
        raise RuntimeError("offline")
    monkeypatch.setattr(BSG, "fast_build_genome", boom)
    e = build_one(world["tmp"] / "e3b", "asv_n2", ctx, str(cache2))
    assert e["n_genomes"] == 1 and e["pfunc"]["applied"] is False
    assert "probability_kn" not in jl(files(world["tmp"] / "e3b", "asv_n2")["provenance"])


def test_prepare_pfunc_scope_store_fetch_and_shared_cache_untouched(world, monkeypatch):
    before = (rb(world["pgfam"]), os.stat(world["pgfam"]).st_mtime_ns)
    import prefetch_gene_families as PF
    calls = []

    def fake_fetch(gid):
        calls.append(gid)
        return ["PGF_A", "PGF_Z"]
    monkeypatch.setattr(PF, "fetch_pgfams", fake_fetch)
    out = world["tmp"] / "out"
    ctx, meta = BSG.prepare_pfunc(AUDIT, str(out), mode="e3", pgfam_cache=world["pgfam"],
                                  catalogue=world["catalogue"], fetch=True, taxdb=None, log=None)
    # in scope: asv_n2 and asv_n3 (n_selected >= 2); asv_n1 / asv_abst are not
    assert meta["n_asvs_in_scope"] == 2
    assert calls == ["5.1"]
    assert jl(out / BSG.PFUNC_DIR / BSG.PFUNC_FETCHED) == {"5.1": ["PGF_A", "PGF_Z"]}
    saved = jl(out / BSG.PFUNC_DIR / BSG.PFUNC_STORE)
    assert set(saved) == {"1.1", "2.1", "3.1", "4.1", "5.1", "11.1", "12.1", "13.1"}
    assert ctx.store.get("4.1") == frozenset()
    assert (rb(world["pgfam"]), os.stat(world["pgfam"]).st_mtime_ns) == before
    assert meta["pgfam_store"]["fetched"] == 1 and meta["estimator"] == "E3"
    # with n_by_asv (the retro path) the built n decides the scope
    _c, m2 = BSG.prepare_pfunc(AUDIT, str(world["tmp"] / "o2"), mode="e3",
                               pgfam_cache=world["pgfam"], catalogue=world["catalogue"],
                               fetch=False, taxdb=None, log=None, n_by_asv={"asv_n2": 2})
    assert m2["n_asvs_in_scope"] == 1
    # kn: nothing prepared, nothing written
    c3, m3 = BSG.prepare_pfunc(AUDIT, str(world["tmp"] / "o3"), mode="kn")
    assert not c3.e3 and m3["estimator"] == "kn" and not os.path.exists(world["tmp"] / "o3")


def test_pfunc_context_validation():
    with pytest.raises(ValueError):
        BSG.PfuncContext("e3")
    with pytest.raises(ValueError):
        BSG.PfuncContext("bogus")
    with pytest.raises(ValueError):
        BSG.PfuncContext("kn", dedup="strain")


# =========================================================================== CLIs
def run_serial(world, out, mode, monkeypatch):
    argv = ["build_synthetic_genomes.py", "--selection", world["selection"], "--out-dir", str(out),
            "--genome-cache-dir", world["cache"], "--pfunc", mode,
            "--pfunc-pgfam-cache", world["pgfam"], "--pfunc-catalogue", world["catalogue"],
            "--pfunc-no-fetch"]
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.chdir(world["tmp"])
    BSG.main()
    return jl(out / "synthetic_genome_manifest.json")


def run_parallel(world, out, mode, monkeypatch, procs=2):
    argv = ["build_synthetic_genomes_parallel.py", "--selection", world["selection"],
            "--out-dir", str(out), "--genome-cache-dir", world["cache"], "--procs", str(procs),
            "--io-workers", "2", "--pfunc", mode, "--pfunc-pgfam-cache", world["pgfam"],
            "--pfunc-catalogue", world["catalogue"], "--pfunc-no-fetch"]
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.chdir(world["tmp"])
    BSGP.main()
    return jl(out / "synthetic_genome_manifest.json")


def test_serial_and_parallel_builders_agree_in_both_modes(world, monkeypatch):
    built = ("asv_n2", "asv_n1", "asv_n3")
    for mode in ("kn", "e3"):
        ms = run_serial(world, world["tmp"] / f"s_{mode}", mode, monkeypatch)
        mp_ = run_parallel(world, world["tmp"] / f"p_{mode}", mode, monkeypatch)
        assert ms["_meta"]["pfunc"]["estimator"] == ("E3" if mode == "e3" else "kn")
        assert mp_["_meta"]["pfunc"]["estimator"] == ms["_meta"]["pfunc"]["estimator"]
        assert mp_["_meta"]["selection"] == world["selection"]
        for a in built:
            for k in ("root", "genome_objects", "provenance"):
                assert rb(files(world["tmp"] / f"s_{mode}", a)[k]) == \
                    rb(files(world["tmp"] / f"p_{mode}", a)[k]), (mode, a, k)
            assert ms[a] == mp_[a]
        assert ms["asv_abst"]["status"] == "abstain"
    # kn builds carry no pfunc entries / files; e3 builds do
    assert "pfunc" not in ms_kn_entry(world) and ms["asv_n2"]["pfunc"]["applied"]
    assert ms["_meta"]["pfunc"]["n_applied"] == 2          # asv_n2, asv_n3 (asv_n1: n = 1)
    assert ms["_meta"]["pfunc"]["not_applied_by_reason"] == {"n_built < min_n": 1}
    for a in built:   # e3 vs kn: only probabilities differ
        for k in ("root", "genome_objects"):
            assert only_probabilities_differ(rb(files(world["tmp"] / "s_kn", a)[k]),
                                             rb(files(world["tmp"] / "s_e3", a)[k]))


def ms_kn_entry(world):
    return jl(world["tmp"] / "s_kn" / "synthetic_genome_manifest.json")["asv_n2"]


def test_retro_cli_matches_build_time_e3_and_never_touches_source(world, monkeypatch):
    run_serial(world, world["tmp"] / "s_kn", "kn", monkeypatch)
    run_serial(world, world["tmp"] / "s_e3", "e3", monkeypatch)
    src_dir = world["tmp"] / "s_kn"
    before = tree_digest(src_dir)
    common = ["--build-dir", str(src_dir), "--genome-cache-dir", world["cache"],
              "--pfunc-pgfam-cache", world["pgfam"], "--pfunc-catalogue", world["catalogue"],
              "--pfunc-no-fetch"]
    for procs in (1, 2):
        out = world["tmp"] / f"r{procs}"
        assert R.main(common + ["--out-dir", str(out), "--procs", str(procs)]) == 0
        man = jl(out / "synthetic_genome_manifest.json")
        e3m = jl(world["tmp"] / "s_e3" / "synthetic_genome_manifest.json")
        for a in ("asv_n2", "asv_n3", "asv_n1"):
            for k in ("root", "genome_objects", "provenance"):
                assert rb(files(out, a)[k]) == rb(files(world["tmp"] / "s_e3", a)[k]), (procs, a, k)
            assert man[a]["pfunc"] == e3m[a]["pfunc"]
        for k in ("root", "genome_objects", "provenance"):      # n = 1: a byte copy
            assert rb(files(out, "asv_n1")[k]) == rb(files(src_dir, "asv_n1")[k])
        assert rb(out / "taxonomies" / "asv_n2.json") == rb(src_dir / "taxonomies" / "asv_n2.json")
        assert man["_meta"]["pfunc"]["n_applied"] == 2 and man["asv_abst"]["status"] == "abstain"
        assert man["_meta"]["reweighted_from"]["build_dir"] == str(src_dir)
        assert man["_meta"]["selection"] == world["selection"]
    assert tree_digest(src_dir) == before


def test_retro_cli_limit_asvs_and_refusals(world, monkeypatch):
    run_serial(world, world["tmp"] / "s_kn", "kn", monkeypatch)
    src_dir = world["tmp"] / "s_kn"
    common = ["--build-dir", str(src_dir), "--genome-cache-dir", world["cache"],
              "--pfunc-pgfam-cache", world["pgfam"], "--pfunc-catalogue", world["catalogue"],
              "--pfunc-no-fetch"]
    out = world["tmp"] / "lim"
    assert R.main(common + ["--out-dir", str(out), "--limit", "2"]) == 0
    man = jl(out / "synthetic_genome_manifest.json")
    assert [a for a in man if a != "_meta"] == ["asv_n2", "asv_n1"]
    assert man["_meta"]["built"] == 2 and man["_meta"]["abstain"] == 0
    lst = world["tmp"] / "list.txt"
    lst.write_text("asv_n3\nasv_abst\n")
    out2 = world["tmp"] / "lst"
    assert R.main(common + ["--out-dir", str(out2), "--asvs", str(lst)]) == 0
    assert [a for a in jl(out2 / "synthetic_genome_manifest.json") if a != "_meta"] == ["asv_n3"]
    with pytest.raises(SystemExit):                      # existing, non-empty out dir
        R.main(common + ["--out-dir", str(out)])
    with pytest.raises(SystemExit):                      # out dir inside the build dir
        R.main(common + ["--out-dir", str(src_dir / "e3")])
    run_serial(world, world["tmp"] / "s_e3", "e3", monkeypatch)
    with pytest.raises(SystemExit):                      # reweighting an E3 build is refused
        R.main(["--build-dir", str(world["tmp"] / "s_e3"), "--out-dir", str(world["tmp"] / "x"),
                "--pfunc-no-fetch", "--pfunc-pgfam-cache", world["pgfam"],
                "--pfunc-catalogue", world["catalogue"], "--genome-cache-dir", world["cache"]])
    # kn mode: a plain copy
    out3 = world["tmp"] / "kncopy"
    assert R.main(common + ["--out-dir", str(out3), "--pfunc", "kn"]) == 0
    for a in ("asv_n2", "asv_n1", "asv_n3"):
        for k in ("root", "genome_objects", "provenance"):
            assert rb(files(out3, a)[k]) == rb(files(src_dir, a)[k])
    assert "pfunc" not in jl(out3 / "synthetic_genome_manifest.json")["asv_n2"]


# =========================================================================== descriptions
OLD_NOTE = "Feature.probability = fraction of source genomes carrying the function (also in inference_data"
# the k/n method sentence; since the BV-BRC-only rule (2026-11-13) it also names the functions'
# annotation: FUNCTIONS_BVBRC for a stamped build, FUNCTIONS_MIXED for a build made before it
# (was: no functions clause at all)
def old_method(functions):
    return ("build_synthetic_genomes_parallel.py: per-ASV union of the selected genomes' genes by "
            "function (bvbrc_to_kbase_genome.LocalGenomeConverter.create_synthetic_genome), the "
            f"functions being {functions}, each "
            "feature tagged probability = fraction of source genomes carrying the function, "
            "selection_provenance core/accessory, source_genome_count")


OLD_METHOD = old_method(P.FUNCTIONS_MIXED)


def test_description_strings_follow_the_manifest():
    U = pytest.importorskip("upload_kbase_genomes")
    import write_mapping_provenance as W
    g = {"features": [], "cdss": []}
    s = {"selected": [{"genome_id": "1.1", "organism": "o", "identity": 0.99, "role": "anchor"}],
         "tier": "species", "best_identity": 0.99, "flags": []}
    kn_note = U.to_kbase_genome(g, asv="a", sel=s, man={}, commit="c")["notes"]
    assert OLD_NOTE in kn_note
    assert OLD_NOTE in U.to_kbase_genome(g, asv="a", sel=s, man={}, commit="c",
                                         pfunc_meta=P.estimator_metadata(mode="kn"))["notes"]
    meta = {**P.estimator_metadata(), "min_n": 2}
    on = U.to_kbase_genome(g, asv="a", sel=s, man={"pfunc": {"applied": True,
                                                             "n_features_estimated": 7}},
                           commit="c", pfunc_meta=meta)["notes"]
    assert "E3 estimate" in on and "the 7 features with a PGFam" in on and OLD_NOTE not in on
    off = U.to_kbase_genome(g, asv="a", sel=s, man={"pfunc": {"applied": False,
                                                              "reason": "n_built 1 < min_n 2"}},
                            commit="c", pfunc_meta=meta)["notes"]
    assert "E3 not applied (n_built 1 < min_n 2)" in off
    assert W.build_method_text({}) == OLD_METHOD
    assert W.build_method_text({"pfunc": P.estimator_metadata(mode="kn")}) == OLD_METHOD
    assert "E3 estimate" in W.build_method_text({"pfunc": meta})
    # BV-BRC PATRIC-only builds (manifest _meta.annotation) say so; unstamped builds say mixed
    st = {"annotation": BSG.ANNOTATION_STAMP}
    assert W.build_method_text(st) == old_method(P.FUNCTIONS_BVBRC)
    assert "BV-BRC (PATRIC) annotations" in W.build_method_text({**st, "pfunc": meta})
    assert "RefSeq" in OLD_METHOD and "3ce12d0" in OLD_METHOD
    g_new = U.to_kbase_genome(g, asv="a", sel=s, man={}, commit="c",
                              annotation=BSG.ANNOTATION_STAMP)
    g_old = U.to_kbase_genome(g, asv="a", sel=s, man={}, commit="c")
    assert f"Functions: {P.FUNCTIONS_BVBRC}." in g_new["notes"] and OLD_NOTE in g_new["notes"]
    assert f"Functions: {P.FUNCTIONS_MIXED}." in g_old["notes"]
    assert "BV-BRC (PATRIC) annotation functions" in g_new["ontology_events"][0]["method"]
    assert "RefSeq" in g_old["ontology_events"][0]["method"]
    assert "RefSeq" not in g_new["ontology_events"][0]["method"] + \
        g_new["ontology_events"][0]["description"]


# =========================================================================== fix round 2
import shutil  # noqa: E402
import socket  # noqa: E402
import subprocess  # noqa: E402

import prefetch_gene_families as PF  # noqa: E402

BUILT = ("asv_n2", "asv_n1", "asv_n3")
KINDS4 = ("root", "genome_objects", "provenance", "taxonomies")


def files4(d, asv):
    return {**files(d, asv), "taxonomies": os.path.join(d, "taxonomies", f"{asv}.json")}


def write_cache_and_catalogue(tmp, name, cache, genomes):
    """A PGFam cache + a catalogue that is current for it."""
    pg = tmp / f"{name}_pgfam.json"
    pg.write_text(json.dumps(cache))
    st = os.stat(pg)
    cat = tmp / f"{name}_catalogue.json"
    cat.write_text(json.dumps({
        "_meta": {"source_cache_size": st.st_size,
                  "source_cache_mtime": time.strftime("%Y-%m-%dT%H:%M:%S%z",
                                                      time.localtime(st.st_mtime))},
        "genomes": genomes, "taxa": {}}))
    return str(pg), str(cat)


@pytest.fixture
def offline(monkeypatch):
    """No network: every socket connect fails; BV-BRC PGFam fetches are recorded, then fail."""
    def no_net(*a, **k):
        raise OSError("network disabled in this test")
    monkeypatch.setattr(socket.socket, "connect", no_net)
    monkeypatch.setattr(socket, "create_connection", no_net)
    calls = []

    def fetch(gid, **kw):
        calls.append(gid)
        raise OSError("network disabled in this test")
    monkeypatch.setattr(PF, "fetch_pgfams", fetch)
    return calls


def cli(world, out, extra, monkeypatch, parallel=False, pgfam=None, catalogue=None,
        cache=None, selection=None):
    """Run a builder CLI on the fixture world; `extra` is appended verbatim (no --pfunc
    unless the test passes one, so the default is exercised)."""
    if parallel:
        argv = ["build_synthetic_genomes_parallel.py", "--procs", "2", "--io-workers", "2"]
    else:
        argv = ["build_synthetic_genomes.py"]
    argv += ["--selection", selection or world["selection"], "--out-dir", str(out),
             "--genome-cache-dir", cache or world["cache"],
             "--pfunc-pgfam-cache", pgfam or world["pgfam"],
             "--pfunc-catalogue", catalogue or world["catalogue"]] + list(extra)
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.chdir(world["tmp"])
    (BSGP.main if parallel else BSG.main)()
    return jl(out / "synthetic_genome_manifest.json")


def test_cli_default_is_e3():
    import argparse
    for default_call in (BSG.add_pfunc_args,):
        ap = argparse.ArgumentParser()
        default_call(ap)
        a = ap.parse_args([])
        assert a.pfunc == "e3" and a.pfunc_min_n == 2 and not a.pfunc_allow_stale_catalogue
        assert ap.parse_args(["--pfunc", "kn"]).pfunc == "kn"


def test_default_e3_build_runs_offline_when_every_genome_is_cached(world, offline, monkeypatch):
    full = {**PGFAM_CACHE, "5.1": ["PGF_A", "PGF_F"]}
    pg, cat = write_cache_and_catalogue(world["tmp"], "full", full, jl(world["catalogue"])["genomes"])
    ms = cli(world, world["tmp"] / "def_s", [], monkeypatch, pgfam=pg, catalogue=cat)
    mp_ = cli(world, world["tmp"] / "def_p", [], monkeypatch, parallel=True, pgfam=pg, catalogue=cat)
    me3 = cli(world, world["tmp"] / "exp_e3", ["--pfunc", "e3"], monkeypatch, pgfam=pg, catalogue=cat)
    assert offline == []                                   # nothing fetched, no network used
    for m in (ms, mp_):
        pm = m["_meta"]["pfunc"]
        assert pm["estimator"] == "E3" and pm["n_applied"] == 2
        assert pm["pgfam_store"]["genomes_absent"] == 0 and pm["pgfam_store"].get("fetched", 0) == 0
        assert m["asv_n2"]["pfunc"]["applied"] and m["asv_n3"]["pfunc"]["n_neighbours_missing"] == 0
    for a in BUILT:                                        # default == explicit e3, serial == parallel
        for k in ("root", "genome_objects", "provenance"):
            ref = rb(files(world["tmp"] / "exp_e3", a)[k])
            assert rb(files(world["tmp"] / "def_s", a)[k]) == ref
            assert rb(files(world["tmp"] / "def_p", a)[k]) == ref
    # a genome missing from the cache with the network down: the build still completes,
    # the failed fetch is counted (not silent) and the genome is left out of E3
    m2 = cli(world, world["tmp"] / "def_missing", [], monkeypatch)
    assert offline == ["5.1"]
    st = m2["_meta"]["pfunc"]["pgfam_store"]
    assert st["fetch_failed"] == 1 and st["genomes_absent"] == 1 and st["absent_examples"] == ["5.1"]
    assert m2["asv_n2"]["pfunc"]["applied"] and m2["asv_n2"]["pfunc"]["n_neighbours_missing"] == 1


def test_explicit_kn_reproduces_the_pre_e3_build(world, monkeypatch):
    m = cli(world, world["tmp"] / "kn", ["--pfunc", "kn"], monkeypatch)
    ref = world["tmp"] / "ref"
    for a in BUILT:                     # pfunc=None is the pre-E3 code path
        e = build_one(ref, a, None, world["cache"])
        assert m[a] == e and "pfunc" not in m[a]
        for k in KINDS4:
            assert rb(files4(world["tmp"] / "kn", a)[k]) == rb(files4(ref, a)[k]), (a, k)
    assert m["_meta"]["pfunc"] == {"estimator": "kn", "estimator_id": "k/n",
                                   "gating_rule": "k/n everywhere (--pfunc kn)", "min_n": 2}
    # the input stamps (every mode): only OUT/pfunc/source_digests.json and two _meta keys
    assert os.listdir(world["tmp"] / "kn" / BSG.PFUNC_DIR) == [BSG.SOURCE_DIGESTS]
    st = os.stat(world["selection"])
    assert m["_meta"]["selection_stamp"] == {"sha256": _sha(world["selection"]),
                                             "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    assert m["_meta"]["source_digests"]["file"] == "pfunc/source_digests.json"
    assert m["_meta"]["source_digests"]["n_genomes"] == 3
    # kn needs no catalogue at all
    m2 = cli(world, world["tmp"] / "kn2", ["--pfunc", "kn"], monkeypatch,
             catalogue=str(world["tmp"] / "absent_catalogue.json"))
    assert m2["_meta"]["pfunc"]["estimator"] == "kn"


def test_missing_or_stale_catalogue_fails_before_selection(world, monkeypatch):
    src = world["tmp"] / "r_src"                  # a k/n build for the retro tool
    run_serial(world, src, "kn", monkeypatch)

    def must_not_run(*a, **k):
        raise AssertionError("selection / warm / audit streaming ran before the catalogue check")
    absent = str(world["tmp"] / "absent_catalogue.json")
    monkeypatch.setattr(BSG, "get_selection", must_not_run)
    with pytest.raises(SystemExit) as ei:
        cli(world, world["tmp"] / "x1", [], monkeypatch, catalogue=absent)
    msg = str(ei.value.code)
    assert "not found" in msg and "build-catalogue" in msg and "systemd-run" in msg
    assert "--pfunc kn" in msg and "--pfunc-allow-stale-catalogue" in msg
    monkeypatch.setattr(BSGP, "warm_cache", must_not_run)
    monkeypatch.setattr(BSG, "_load", must_not_run)
    with pytest.raises(SystemExit) as ei:
        cli(world, world["tmp"] / "x2", [], monkeypatch, parallel=True, catalogue=absent)
    assert "build-catalogue" in str(ei.value.code)
    monkeypatch.setattr(R, "stream_records", must_not_run)
    with pytest.raises(SystemExit) as ei:
        R.main(["--build-dir", str(src), "--out-dir", str(world["tmp"] / "r_x"),
                "--genome-cache-dir", world["cache"], "--pfunc-catalogue", absent,
                "--pfunc-pgfam-cache", world["pgfam"]])
    assert "build-catalogue" in str(ei.value.code)
    # stale: the PGFam cache changed after the catalogue was built
    with open(world["pgfam"], "a") as f:
        f.write(" ")
    for parallel in (False, True):
        with pytest.raises(SystemExit) as ei:
            cli(world, world["tmp"] / f"x3{parallel}", [], monkeypatch, parallel=parallel)
        assert "STALE" in str(ei.value.code) and "build-catalogue" in str(ei.value.code)
    with pytest.raises(SystemExit) as ei:
        R.main(["--build-dir", str(src), "--out-dir", str(world["tmp"] / "r_y"),
                "--genome-cache-dir", world["cache"], "--pfunc-catalogue", world["catalogue"],
                "--pfunc-pgfam-cache", world["pgfam"]])
    assert "STALE" in str(ei.value.code)


def test_stale_catalogue_allowed_by_flag_and_recorded(world, monkeypatch):
    with open(world["pgfam"], "a") as f:
        f.write(" ")
    m = cli(world, world["tmp"] / "st", ["--pfunc-allow-stale-catalogue", "--pfunc-no-fetch"],
            monkeypatch)
    c = m["_meta"]["pfunc"]["catalogue"]
    assert c["current_vs_pgfam_cache"] is False and c["stale_allowed"] is True
    assert c["staled_by_this_run"] is False and m["asv_n2"]["pfunc"]["applied"]


def test_catalogue_staled_by_this_run_is_warned_not_fatal(world):
    pre = BSG.pfunc_preflight("e3", world["catalogue"], world["pgfam"], log=None)
    assert pre["current"] and BSG.pfunc_preflight("kn", "nope", "nope") is None
    # the run's own selection/prefetch writes a new genome into the PGFam cache
    time.sleep(1.1)
    with open(world["pgfam"], "w") as f:
        json.dump({**PGFAM_CACHE, "99.1": ["PGF_Q"]}, f)
    msgs = []
    ctx, meta = BSG.prepare_pfunc(AUDIT, str(world["tmp"] / "o"), mode="e3",
                                  pgfam_cache=world["pgfam"], catalogue=world["catalogue"],
                                  fetch=False, taxdb=None, log=msgs.append, preflight=pre)
    assert meta["catalogue"]["staled_by_this_run"] is True
    assert meta["catalogue"]["current_vs_pgfam_cache"] is False
    assert any("this run's selection/prefetch" in m and "build-catalogue" in m for m in msgs)
    with pytest.raises(SystemExit):              # the same staleness without the snapshot
        BSG.prepare_pfunc(AUDIT, str(world["tmp"] / "o2"), mode="e3", pgfam_cache=world["pgfam"],
                          catalogue=world["catalogue"], fetch=False, taxdb=None, log=None)


def _sha(p):
    return hashlib.sha256(rb(p)).hexdigest()


def test_e3_stamp_records_code_provenance(world, monkeypatch):
    ms = cli(world, world["tmp"] / "cs", ["--pfunc-no-fetch"], monkeypatch)
    mp_ = cli(world, world["tmp"] / "cp", ["--pfunc-no-fetch"], monkeypatch, parallel=True)
    for m, extra in ((ms, set()), (mp_, {"build_synthetic_genomes_parallel.py"})):
        cp = m["_meta"]["pfunc"]["code_provenance"]
        assert set(cp["files"]) == {"pfunc_estimator.py", "build_synthetic_genomes.py"} | extra
        for rel, rec in cp["files"].items():
            assert rec["sha256"] == _sha(os.path.join(ROOT, rel))
            assert rec["git_status"] is None or isinstance(rec["git_status"], str)
        assert cp["git_head"] is None or len(cp["git_head"]) == 40
        assert cp["dirty"] in (True, False, None)
        if cp["files"]["pfunc_estimator.py"]["git_status"]:
            assert cp["dirty"] is True
    # the kn stamp carries none of it
    mk = cli(world, world["tmp"] / "ck", ["--pfunc", "kn"], monkeypatch)
    assert "code_provenance" not in mk["_meta"]["pfunc"]


def test_build_time_e3_error_keeps_kn_and_the_genome(world, monkeypatch):
    kn_dir = world["tmp"] / "kn"
    build_one(kn_dir, "asv_n2", None, world["cache"])
    ctx, _ = e3_ctx(world, world["tmp"] / "e3")

    def boom(*a, **k):
        raise RuntimeError("boom")
    with monkeypatch.context() as mp:
        mp.setattr(P, "reweight_asv", boom)
        e = build_one(world["tmp"] / "e3a", "asv_n2", ctx, world["cache"])
    assert e["status"] == "built" and e["pfunc"]["applied"] is False
    assert e["pfunc"]["reason"] == "error: RuntimeError: boom" and e["pfunc"]["n_built"] == 2
    for k in ("root", "genome_objects", "provenance"):
        assert rb(files(world["tmp"] / "e3a", "asv_n2")[k]) == rb(files(kn_dir, "asv_n2")[k])
    # E3 applied in memory, then the genome_objects rewrite fails: every probability restored
    with monkeypatch.context() as mp:
        mp.setattr(BSG, "_dump_json", lambda *a, **k: (_ for _ in ()).throw(OSError("disk full")))
        e2 = build_one(world["tmp"] / "e3b", "asv_n2", ctx, world["cache"])
    assert e2["status"] == "built" and e2["pfunc"]["reason"] == "error: OSError: disk full"
    for k in ("root", "genome_objects", "provenance"):
        assert rb(files(world["tmp"] / "e3b", "asv_n2")[k]) == rb(files(kn_dir, "asv_n2")[k])
    assert not [f for f in os.listdir(world["tmp"] / "e3b" / "genome_objects") if ".tmp." in f]
    c = BSG.pfunc_counts([e, e2, {"pfunc": {"applied": True, "clade_rank": "family"}}])
    assert c["n_errors"] == 2 and c["not_applied_by_reason"] == {"error": 2} and c["n_applied"] == 1


def failed_rep_world(world):
    """Catalogue where 2.1 is in FamB (members 21.1-23.1), a genome cache without 1.1."""
    genomes = {**jl(world["catalogue"])["genomes"], "2.1": [2, "FamB", "Ord"],
               "21.1": [21, "FamB", "Ord"], "22.1": [22, "FamB", "Ord"], "23.1": [23, "FamB", "Ord"]}
    cache = {**PGFAM_CACHE, "21.1": ["PGF_A"], "22.1": ["PGF_D"], "23.1": ["PGF_X"]}
    pg, cat = write_cache_and_catalogue(world["tmp"], "frw", cache, genomes)
    gc2 = world["tmp"] / "cache_no_1_1"
    gc2.mkdir()
    for g in ("2.1", "3.1"):
        (gc2 / f"{g}.json").write_text(json.dumps(SOURCES[g]))
    return pg, cat, str(gc2)


def test_failed_first_rep_prior_anchors_on_first_built_rep(world, offline, monkeypatch):
    pg, cat, gc2 = failed_rep_world(world)
    recs = {"asv_n3": AUDIT["asv_n3"]}
    monkeypatch.setattr(BSG, "fast_build_genome", lambda gid: (_ for _ in ()).throw(
        RuntimeError("offline")))
    # before the merge every rep may anchor (a cached file can be unreadable, an uncached one
    # can fail to fetch); with the built ids known (retro tool) exactly the first built rep
    assert BSG.anchor_candidates(AUDIT["asv_n3"]) == ["1.1", "2.1", "3.1"]
    assert BSG.anchor_candidates(AUDIT["asv_n3"], built_ids=["3.1", "2.1"]) == ["2.1"]
    assert BSG.anchor_candidates({"selected": []}) == []
    ctx, meta = BSG.prepare_pfunc(recs, str(world["tmp"] / "o"), mode="e3", pgfam_cache=pg,
                                  catalogue=cat, fetch=False, taxdb=None, log=None)
    assert {"21.1", "22.1", "23.1"} <= ctx.store.genome_ids()
    e = build_one(world["tmp"] / "b", "asv_n3", ctx, gc2)
    assert e["n_genomes"] == 2 and e["failed"][0]["genome_id"] == "1.1" and e["pfunc"]["applied"]
    d = jl(files(world["tmp"] / "b", "asv_n3")["provenance"])["pfunc"]
    assert d["anchor"] == "2.1" and d["anchor_selected"] == "1.1" and d["anchor_built"] is False
    assert d["reps_not_built"] == ["1.1"] and d["anchor_family"] == "FamB"
    assert d["clade_rank"] == "family" and d["n_members_used"] == 4      # 2.1, 21.1, 22.1, 23.1
    assert d["n_neighbours"] == 5                  # 1.1 (not built) is still a neighbour
    assert d["source_digests"] == P.source_digests(["2.1", "3.1"], [PATRIC_SOURCES["2.1"],
                                                                     PATRIC_SOURCES["3.1"]])
    # a store scoped on selected[0] only (the pre-fix scope) would lose FamB's members: prior
    # dropped. A corrupt cached file of selected[0] (exists, does not parse) is the same case.
    (world["tmp"] / "cache_bad").mkdir()
    for g in ("2.1", "3.1"):
        shutil.copyfile(os.path.join(gc2, f"{g}.json"), world["tmp"] / "cache_bad" / f"{g}.json")
    (world["tmp"] / "cache_bad" / "1.1.json").write_text("{not json")
    e_bad = build_one(world["tmp"] / "bb", "asv_n3", ctx, str(world["tmp"] / "cache_bad"))
    d_bad = jl(files(world["tmp"] / "bb", "asv_n3")["provenance"])["pfunc"]
    assert e_bad["failed"][0]["genome_id"] == "1.1" and d_bad["anchor"] == "2.1"
    assert d_bad["clade_rank"] == "family" and d_bad["n_members_used"] == 4
    narrow = P.PGFamStore({g: ctx.store.get(g) for g in
                           P.required_genomes(AUDIT["asv_n3"], P.load_clade_catalogue(cat))
                           if ctx.store.get(g) is not None})
    ctx0 = BSG.PfuncContext("e3", store=narrow, catalogue=P.load_clade_catalogue(cat))
    build_one(world["tmp"] / "b0", "asv_n3", ctx0, gc2)
    d0 = jl(files(world["tmp"] / "b0", "asv_n3")["provenance"])["pfunc"]
    assert d0["clade_dropped_after_missing"] and d0["clade_rank"] is None


# --------------------------------------------------------------------------- retro guards
def retro(world, src, out, monkeypatch=None, cache=None, selection=None, extra=()):
    argv = ["--build-dir", str(src), "--out-dir", str(out), "--genome-cache-dir",
            cache or world["cache"], "--pfunc-pgfam-cache", world["pgfam"],
            "--pfunc-catalogue", world["catalogue"], "--pfunc-no-fetch"]
    if selection:
        argv += ["--selection", selection]
    rc = R.main(argv + list(extra))
    return rc, jl(out / "synthetic_genome_manifest.json")


def assert_byte_copy(out, src, asv):
    for k in KINDS4:
        assert rb(files4(out, asv)[k]) == rb(files4(src, asv)[k]), (asv, k)


MANIFEST = "synthetic_genome_manifest.json"


def make_legacy(src, audits=(), caches=()):
    """Turn a build into a pre-stamp (legacy) one: drop _meta.selection_stamp,
    _meta.source_digests and pfunc/source_digests.json; set the manifest's mtime to now and
    the given audits' / cache dirs' files 100 s older (os.utime, not write order)."""
    man = jl(src / MANIFEST)
    man["_meta"].pop("selection_stamp")
    man["_meta"].pop("source_digests")
    (src / MANIFEST).write_text(json.dumps(man, indent=1))
    os.remove(src / BSG.PFUNC_DIR / BSG.SOURCE_DIGESTS)
    t = time.time()
    os.utime(src / MANIFEST, (t, t))
    for p in list(audits) + [os.path.join(c, f) for c in caches for f in os.listdir(c)]:
        os.utime(p, (t - 100, t - 100))
    return t


def set_mtime(path, t):
    os.utime(path, (t, t))


def test_every_build_writes_the_audit_stamp_and_source_digests(world, monkeypatch):
    """Digests and stamp come from a real (fixture) build, in both modes and both builders."""
    want = P.source_digests(["1.1", "2.1", "3.1"], [PATRIC_SOURCES[g] for g in ("1.1", "2.1", "3.1")])
    st = os.stat(world["selection"])
    stamp = {"sha256": _sha(world["selection"]), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    ref = None
    for mode in ("kn", "e3"):
        for parallel in (False, True):
            out = world["tmp"] / f"st_{mode}_{parallel}"
            m = cli(world, out, ["--pfunc", mode, "--pfunc-no-fetch"], monkeypatch,
                    parallel=parallel)
            dfile = out / BSG.PFUNC_DIR / BSG.SOURCE_DIGESTS
            assert jl(dfile) == want                       # every built source genome
            ref = ref or rb(dfile)
            assert rb(dfile) == ref                        # same bytes in every mode / builder
            assert m["_meta"]["selection_stamp"] == stamp
            assert m["_meta"]["source_digests"]["n_genomes"] == 3
            assert m["_meta"]["source_digests"]["conflicting_genome_ids"] == []
            for a in BUILT:                                # the entries never carry them
                assert "source_digests" not in json.dumps(m[a])
            if mode == "e3":                               # == what E3 records per ASV
                d = jl(files(out, "asv_n2")["provenance"])["pfunc"]["source_digests"]
                assert d == {g: want[g] for g in ("1.1", "2.1")}
    run = {}
    conflicts = set()
    BSG.merge_source_digests(run, {"1.1": "a", "2.1": "b"}, conflicts)
    BSG.merge_source_digests(run, {"1.1": "a", "2.1": "c", "3.1": "d"}, conflicts)
    BSG.merge_source_digests(run, {"2.1": "b"}, conflicts)
    assert run == {"1.1": "a", "2.1": None, "3.1": "d"} and conflicts == {"2.1"}
    assert BSG.file_stamp(str(world["tmp"] / "nope.json")) is None


def test_retro_input_drift_guards(world, monkeypatch):
    src = world["tmp"] / "s_kn"
    run_serial(world, src, "kn", monkeypatch)
    # (a) a re-fetched product string in a source genome: function_prevalence differs
    drift = world["tmp"] / "cache_drift"
    shutil.copytree(world["cache"], drift)
    s2 = copy.deepcopy(SOURCES["2.1"])
    s2["features"][3]["functions"] = ["prod e"]
    (drift / "2.1.json").write_text(json.dumps(s2))
    rc, m = retro(world, src, world["tmp"] / "ra", cache=str(drift))
    assert rc == 0
    for a in ("asv_n2", "asv_n3"):
        assert m[a]["pfunc"]["reason"] == "inputs drifted since build: function_prevalence"
        assert_byte_copy(world["tmp"] / "ra", src, a)
    assert m["_meta"]["pfunc"]["n_inputs_drifted"] == 2 and m["_meta"]["pfunc"]["n_applied"] == 0
    # (b) the audit's selected genomes differ from what the build used; (c) record missing.
    # On a stamped build these are shadowed by the run-level sha256 check (test below), so
    # they are exercised on a legacy build, with --allow-audit-newer (the drifted audit is
    # newer than the build's manifest).
    leg = world["tmp"] / "s_leg"
    shutil.copytree(src, leg)
    aud = copy.deepcopy(AUDIT)
    # a regenerated audit that selects one more rep than the build used (the old subset
    # check, built genomes <= selected, would pass it)
    aud["asv_n2"]["selected"] = [sel("1.1", 0.99, "Sp one", 1), sel("2.1", 0.99, "Sp two", 2),
                                 sel("3.1", 0.99, "Sp three", 3)]
    del aud["asv_n3"]
    selp = world["tmp"] / "sel_drift.json"
    selp.write_text(json.dumps({"_meta": {"gene_data": "bvbrc"}, **aud}))
    t = make_legacy(leg, caches=[world["cache"]])
    set_mtime(selp, t + 50)
    with pytest.raises(SystemExit) as ei:                     # legacy + audit newer: refused
        retro(world, leg, world["tmp"] / "rb0", selection=str(selp))
    assert "--allow-audit-newer" in str(ei.value.code) and "Nothing was written" in str(ei.value.code)
    assert not os.path.exists(world["tmp"] / "rb0")
    rc, m = retro(world, leg, world["tmp"] / "rb", selection=str(selp),
                  extra=["--allow-audit-newer"])
    assert m["asv_n2"]["pfunc"]["reason"] == "inputs drifted since build: selected genomes"
    assert m["asv_n3"]["pfunc"]["reason"] == "record not in selection audit"
    for a in ("asv_n2", "asv_n3"):
        assert_byte_copy(world["tmp"] / "rb", leg, a)
    pm = m["_meta"]["pfunc"]
    assert pm["n_inputs_drifted"] == 1 and pm["n_record_not_in_audit"] == 1
    sa = pm["input_checks"]["selection_audit"]
    assert sa["allow_audit_newer"] is True and sa["audit_newer_than_build"] is True
    assert "allowed by --allow-audit-newer" in sa["result"]
    # legacy build, audit older than the manifest: runs without any flag
    rc, m = retro(world, leg, world["tmp"] / "rb_ok")
    assert rc == 0 and m["_meta"]["pfunc"]["n_applied"] == 2
    assert m["_meta"]["pfunc"]["input_checks"]["selection_audit"]["result"] == \
        "audit not newer than the build"


def test_retro_stamped_audit_mismatch_keeps_kn_everywhere(world, monkeypatch):
    """FINDING 1: a changed decisions[] (same selected reps) in a re-written audit."""
    src = world["tmp"] / "s_kn"
    run_serial(world, src, "kn", monkeypatch)
    aud = copy.deepcopy(AUDIT)
    aud["asv_n2"]["decisions"][2]["identity"] = 0.96          # one neighbour's identity only
    assert [s_["genome_id"] for s_ in aud["asv_n2"]["selected"]] == ["1.1", "2.1"]
    selp = world["tmp"] / "sel_dec.json"
    selp.write_text(json.dumps({"_meta": {"gene_data": "bvbrc"}, **aud}))
    calls = []
    monkeypatch.setattr(R, "stream_records", lambda *a, **k: calls.append(a) or ({}, {}))
    rc, m = retro(world, src, world["tmp"] / "rs", selection=str(selp), extra=["--procs", "2"])
    assert rc == 0 and calls == []                             # nothing streamed
    for a in BUILT:                                            # n = 1 included
        assert m[a]["status"] == "built"
        assert m[a]["pfunc"]["reason"] == "inputs drifted since build: selection audit"
        assert_byte_copy(world["tmp"] / "rs", src, a)
    pm = m["_meta"]["pfunc"]
    assert pm["n_applied"] == 0 and pm["n_inputs_drifted"] == 3
    assert pm["pgfam_store"]["genomes_needed"] == 0
    sa = pm["input_checks"]["selection_audit"]
    assert sa["result"].startswith("MISMATCH") and sa["stamp_in_build"]["sha256"] == \
        _sha(world["selection"]) and sa["audit_now"]["sha256"] == _sha(selp)
    # the same bytes at another path, newer than the build: the stamp matches -> applied
    monkeypatch.undo()
    cp = world["tmp"] / "sel_copy.json"
    shutil.copyfile(world["selection"], cp)
    set_mtime(cp, time.time() + 100)
    rc, m = retro(world, src, world["tmp"] / "rc", selection=str(cp))
    assert rc == 0 and m["_meta"]["pfunc"]["n_applied"] == 2
    assert m["_meta"]["pfunc"]["input_checks"]["selection_audit"]["result"] == "match"


def test_retro_source_digest_guard_from_a_real_build(world, monkeypatch):
    """FINDING 2: the digests come from the build itself (OUT/pfunc/source_digests.json), so
    a PGFAM-alias-only change in a cached source genome is caught on a k/n build."""
    src = world["tmp"] / "s_kn"
    run_serial(world, src, "kn", monkeypatch)
    recorded = jl(src / BSG.PFUNC_DIR / BSG.SOURCE_DIGESTS)
    assert recorded == {g: P.source_digest(PATRIC_SOURCES[g]) for g in ("1.1", "2.1", "3.1")}
    alias = world["tmp"] / "cache_alias"
    shutil.copytree(world["cache"], alias)
    s2 = copy.deepcopy(SOURCES["2.1"])
    s2["features"][2]["aliases"][1] = ["PGFAM", "PGF_D9"]          # 'Prod D': alias only
    s2["features"][2]["pgfam_id"] = "PGF_D9"
    (alias / "2.1.json").write_text(json.dumps(s2))
    prov = jl(files(src, "asv_n2")["provenance"])
    assert dict(BSG.function_prevalence([PATRIC_SOURCES["1.1"], BSG.bvbrc_only(s2)])) == \
        prov["function_prevalence"]
    rc, m = retro(world, src, world["tmp"] / "rd", cache=str(alias), extra=["--procs", "2"])
    assert rc == 0
    for a in ("asv_n2", "asv_n3"):                                 # both use 2.1
        assert m[a]["pfunc"]["reason"] == "inputs drifted since build: source digest"
        assert m[a]["pfunc"]["drift_genomes"] == ["2.1"]
        assert_byte_copy(world["tmp"] / "rd", src, a)
    assert m["_meta"]["pfunc"]["input_checks"]["source_genomes"]["mode"] == "digests"
    assert rb(world["tmp"] / "rd" / BSG.PFUNC_DIR / BSG.SOURCE_DIGESTS) == \
        rb(src / BSG.PFUNC_DIR / BSG.SOURCE_DIGESTS)              # copied along
    # an unchanged cache with fresh mtimes (a copy) is NOT drift on a digest build
    fresh = world["tmp"] / "cache_fresh"
    shutil.copytree(world["cache"], fresh)
    for f in os.listdir(fresh):
        set_mtime(fresh / f, time.time() + 100)
    rc, m = retro(world, src, world["tmp"] / "rf", cache=str(fresh))
    assert m["_meta"]["pfunc"]["n_applied"] == 2


def test_retro_legacy_build_cache_newer_than_build(world, monkeypatch):
    """FINDING 2, legacy builds (no source_digests.json): a cache file newer than the source
    manifest is drift for the ASVs that use it, unless --allow-cache-newer."""
    src = world["tmp"] / "s_kn"
    run_serial(world, src, "kn", monkeypatch)
    leg = world["tmp"] / "s_leg"
    shutil.copytree(src, leg)
    cache = world["tmp"] / "cache_leg"
    shutil.copytree(world["cache"], cache)
    t = make_legacy(leg, audits=[world["selection"]], caches=[str(cache)])
    set_mtime(cache / "3.1.json", t + 50)                         # used by asv_n3 only
    rc, m = retro(world, leg, world["tmp"] / "rl", cache=str(cache), extra=["--procs", "2"])
    assert rc == 0
    assert m["asv_n3"]["pfunc"]["reason"] == \
        "inputs drifted since build: source genome cache newer than build"
    assert m["asv_n3"]["pfunc"]["drift_genomes"] == ["3.1"]
    assert_byte_copy(world["tmp"] / "rl", leg, "asv_n3")
    assert m["asv_n2"]["pfunc"]["applied"]                         # 1.1, 2.1 not newer
    sg = m["_meta"]["pfunc"]["input_checks"]["source_genomes"]
    assert sg["mode"] == "mtime" and sg["allow_cache_newer"] is False
    assert not os.path.exists(world["tmp"] / "rl" / BSG.PFUNC_DIR / BSG.SOURCE_DIGESTS)
    rc, m = retro(world, leg, world["tmp"] / "rl2", cache=str(cache),
                  extra=["--allow-cache-newer"])
    assert m["asv_n3"]["pfunc"]["applied"]
    assert m["asv_n3"]["pfunc"]["cache_newer_allowed_genomes"] == ["3.1"]
    assert "cache_newer_allowed_genomes" not in m["asv_n2"]["pfunc"]
    sg = m["_meta"]["pfunc"]["input_checks"]["source_genomes"]
    assert sg["allow_cache_newer"] is True and sg["n_asvs_cache_newer_allowed"] == 1
    # the reweighted genome equals the build-time E3 one (the file was only touched)
    run_serial(world, world["tmp"] / "s_e3", "e3", monkeypatch)
    for k in ("root", "genome_objects"):
        assert rb(files(world["tmp"] / "rl2", "asv_n3")[k]) == \
            rb(files(world["tmp"] / "s_e3", "asv_n3")[k])


def test_retro_errors_keep_kn_and_built_always_has_files(world, monkeypatch):
    src = world["tmp"] / "s_kn"
    run_serial(world, src, "kn", monkeypatch)
    real = P.reweight_asv

    def flaky(genome_obj, rec, *a, **k):
        if rec and len(rec.get("selected") or []) == 2:      # asv_n2 only (n1's rec is None)
            raise RuntimeError("estimator blew up")
        return real(genome_obj, rec, *a, **k)
    with monkeypatch.context() as mp:
        mp.setattr(P, "reweight_asv", flaky)
        rc, m = retro(world, src, world["tmp"] / "re")
    assert rc == 1
    assert m["asv_n2"]["status"] == "built"
    assert m["asv_n2"]["pfunc"]["reason"] == "error: RuntimeError: estimator blew up"
    assert_byte_copy(world["tmp"] / "re", src, "asv_n2")
    assert m["asv_n3"]["pfunc"]["applied"]
    pm = m["_meta"]
    assert pm["pfunc"]["n_errors"] == 1 and pm["pfunc"]["n_errors_kept_kn"] == 1
    assert pm["pfunc"]["n_built_demoted_to_failed"] == 0 and pm["built"] == 3
    # an error after some E3 files were written: the fallback copy overwrites them
    real_dump = BSG._dump_json

    def dump_then_fail(obj, path):
        if os.sep + "provenance" + os.sep in path:
            raise OSError("disk full")
        return real_dump(obj, path)
    with monkeypatch.context() as mp:
        mp.setattr(BSG, "_dump_json", dump_then_fail)
        rc, m = retro(world, src, world["tmp"] / "rp")
    assert rc == 1 and m["_meta"]["pfunc"]["n_errors"] == 2
    for a in ("asv_n2", "asv_n3"):
        assert m[a]["status"] == "built" and m[a]["pfunc"]["reason"] == "error: OSError: disk full"
        assert_byte_copy(world["tmp"] / "rp", src, a)
    # a broken source build: missing provenance (n3), truncated root (n2), missing root (n1)
    br = world["tmp"] / "s_broken"
    shutil.copytree(src, br)
    os.remove(br / "provenance" / "asv_n3.json")
    (br / "asv_n2.json").write_bytes(rb(br / "asv_n2.json")[:100])
    os.remove(br / "asv_n1.json")
    rc, m = retro(world, br, world["tmp"] / "rbr")
    assert rc == 1
    assert {a: m[a]["status"] for a in BUILT} == {a: "failed" for a in BUILT}
    assert m["asv_n3"]["pfunc"]["reason"].startswith("error: FileNotFoundError")
    assert m["asv_n2"]["pfunc"]["reason"].startswith("error: JSONDecodeError")
    assert "root" in m["asv_n1"]["error"]
    assert m["_meta"]["built"] == 0 and m["_meta"]["failed"] == 3
    assert m["_meta"]["pfunc"]["n_built_demoted_to_failed"] == 3
    assert m["_meta"]["pfunc"]["n_built_demoted_unreadable_output"] == 1     # asv_n1
    assert m["_meta"]["pfunc"]["n_errors_kept_kn"] == 0
    # FINDING 3: a truncated n = 1 root raises nothing (it is only copied) but must not stay
    # 'built' with an unparseable file
    br1 = world["tmp"] / "s_broken1"
    shutil.copytree(src, br1)
    (br1 / "asv_n1.json").write_bytes(rb(br1 / "asv_n1.json")[: len(rb(br1 / "asv_n1.json")) // 3])
    rc, m = retro(world, br1, world["tmp"] / "rbr1", extra=["--procs", "2"])
    assert rc == 1
    assert m["asv_n1"]["status"] == "failed" and "root" in m["asv_n1"]["error"]
    assert "copied from the source build" in m["asv_n1"]["error"]
    assert m["asv_n1"]["pfunc"]["reason"] == "n_built 1 < min_n 2"
    assert m["asv_n2"]["status"] == "built" and m["asv_n2"]["pfunc"]["applied"]
    pm = m["_meta"]["pfunc"]
    assert pm["n_errors"] == 0 and pm["n_built_demoted_to_failed"] == 1
    assert pm["n_built_demoted_unreadable_output"] == 1
    assert m["_meta"]["built"] == 2 and m["_meta"]["failed"] == 1
    for out in ("re", "rp", "rbr", "rbr1"):        # every 'built' entry has readable files
        mm = jl(world["tmp"] / out / "synthetic_genome_manifest.json")
        for a, e in mm.items():
            if a != "_meta" and e.get("status") == "built":
                assert not R.unusable_files(a, str(world["tmp"] / out), parse=True)


def test_retro_code_provenance_and_kn_copy_meta(world, monkeypatch):
    src = world["tmp"] / "s_kn"
    run_serial(world, src, "kn", monkeypatch)
    rc, m = retro(world, src, world["tmp"] / "rr")
    files_ = m["_meta"]["pfunc"]["code_provenance"]["files"]
    assert "reweight_build_e3.py" in files_ and "pfunc_estimator.py" in files_
    assert files_["reweight_build_e3.py"]["sha256"] == _sha(os.path.join(ROOT, "reweight_build_e3.py"))
    d = jl(files(world["tmp"] / "rr", "asv_n2")["provenance"])["pfunc"]
    assert d["anchor_built"] is True and d["reps_not_built"] == []
    assert set(d["source_digests"]) == {"1.1", "2.1"}


# --------------------------------------------------------------------------- uploader / d100
def test_uploader_records_the_estimator_for_e3_builds_only(monkeypatch):
    U = pytest.importorskip("upload_kbase_genomes")
    monkeypatch.setattr(U, "now_kbase", lambda: "2026-01-01T00:00:00+0000")
    s = {"selected": [{"genome_id": "1.1", "organism": "o", "identity": 0.99, "role": "anchor",
                       "taxon_id": 1}], "tier": "species", "best_identity": 0.99, "flags": [],
         "asv_len": 400}
    man = {"n_features": 5, "n_equivalent_genomes": 0}
    prov = {"pipeline": {"git": {"commit": "c"}}}
    base_meta = U.ws_meta("a", s, man, "c", 0.1)
    base_prov = U.ws_provenance("a", s, prov, "c", "h", "s", "/g")
    kn = {"estimator": "kn", "estimator_id": "k/n", "gating_rule": "x", "min_n": 2}
    for pm in (None, kn):                                   # k/n uploads unchanged
        assert U.ws_meta("a", s, man, "c", 0.1, pfunc_meta=pm) == base_meta
        assert U.ws_provenance("a", s, prov, "c", "h", "s", "/g", pfunc_meta=pm) == base_prov
    assert "probability_estimator" not in base_meta
    assert "feature_probability" not in base_prov[0]["method_params"][0]
    e3 = {**P.estimator_metadata(), "min_n": 2, "n_applied": 7}
    on = {**man, "pfunc": {"applied": True}}
    off = {**man, "pfunc": {"applied": False, "reason": "n_built 1 < min_n 2"}}
    assert U.ws_meta("a", s, on, "c", 0.1, pfunc_meta=e3)["probability_estimator"] == "E3"
    assert U.ws_meta("a", s, off, "c", 0.1, pfunc_meta=e3)["probability_estimator"] == \
        "k/n (E3 not applied: n_built 1 < min_n 2)"
    mp0 = U.ws_provenance("a", s, prov, "c", "h", "s", "/g", pfunc_meta=e3)[0]["method_params"][0]
    assert mp0["feature_probability"] == e3
    assert {k: v for k, v in mp0.items() if k != "feature_probability"} == \
        base_prov[0]["method_params"][0]


def test_build_d100_selection_says_tieset_guard_not_applicable(world):
    tie = world["tmp"] / "tie.json"
    tie.write_text(json.dumps({"tieset": ["1.1"]}))
    out = world["tmp"] / "d100"
    cmd = [sys.executable, os.path.join(ROOT, "deployment_d100", "build_d100.py"),
           "--tieset", str(tie), "--selection", world["selection"], "--out-dir", str(out),
           "--genome-cache-dir", world["cache"], "--procs", "1", "--io-workers", "1",
           "--pfunc-pgfam-cache", world["pgfam"], "--pfunc-catalogue", world["catalogue"],
           "--pfunc-no-fetch"]
    r = subprocess.run(cmd, cwd=world["tmp"], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    assert "tie-set guard not applicable with --selection" in r.stdout
    assert "tie-set guard on" not in r.stdout
    m = jl(out / "synthetic_genome_manifest.json")
    assert m["_meta"]["pfunc"]["estimator"] == "E3"            # the default passes through
    assert "deployment_d100/build_d100.py" in m["_meta"]["pfunc"]["code_provenance"]["files"]
    assert m["_meta"]["annotation"] == BSG.ANNOTATION_STAMP      # BV-BRC PATRIC only
    for a in BUILT:                                              # from the mixed cache
        assert all(has_patric_id(f) for f in jl(out / f"{a}.json")["features"])


# --------------------------------------------------------------------------- d100_pipeline.sh
PIPE = os.path.join(ROOT, "deployment_d100", "d100_pipeline.sh")


def run_pipe(tmp, env_extra, script=PIPE):
    env = {k: v for k, v in os.environ.items() if k != "PFUNC"}
    env.update(NAME="t_d100", FASTA="/nonexistent.fasta", TAXCSV="/nonexistent.csv",
               STAGE="build", MIN_AVAIL_GB="0", **env_extra)
    return subprocess.run(["/bin/bash", script], env=env, capture_output=True, text=True,
                          timeout=120)


def test_d100_pipeline_default_e3_and_refuses_any_existing_syn(tmp_path):
    """Since the BV-BRC-only rule (2026-11-13) no PFUNC reproduces the shipped (mixed
    PATRIC + RefSeq) builds: the default is e3, every PFUNC refuses a non-empty SYN, and the
    audit and stats always go into the new SYN. (Was test_d100_pipeline_pfunc_default_kn_and_
    e3_guard: default kn writing the audit/stats into OUT, refusal for e3 only.)"""
    out, syn = tmp_path / "hits_d100", tmp_path / "syn_d100"
    out.mkdir()
    syn.mkdir()
    (syn / "synthetic_genome_manifest.json").write_text("{}")        # a shipped build
    before = tree_digest(syn)
    for pf in ("e3", "kn", None):               # any PFUNC into an existing, non-empty SYN
        env = {"OUT": str(out), "SYN": str(syn), **({"PFUNC": pf} if pf else {})}
        r = run_pipe(tmp_path, env)
        assert r.returncode == 2 and f"[abort] SYN={syn} exists and is not empty" in r.stdout, \
            r.stdout
        assert f"(PFUNC={pf or 'e3'})" in r.stdout and "_bvbrc_d100" in r.stdout
        assert "3ce12d0" in r.stdout
    r = run_pipe(tmp_path, {"OUT": str(out), "SYN": str(syn), "PFUNC": "bogus"})
    assert r.returncode == 2 and "PFUNC=bogus (must be kn or e3)" in r.stdout
    assert tree_digest(syn) == before and os.listdir(out) == []
    # which --pfunc / audit / stats paths the build stage passes: the script with PY=echo
    stub = tmp_path / "pipe_stub.sh"
    txt = open(PIPE).read()
    assert txt.count("\nPY=/home/freiburger/Documents/py_venv/bin/python\n") == 1
    stub.write_text(txt.replace("\nPY=/home/freiburger/Documents/py_venv/bin/python\n", "\nPY=echo\n"))
    (out / "asv_top20_alignment_hits.json").write_text("{}")
    for pf, name in ((None, "syn_new_d100"), ("kn", "syn_kn_d100")):   # a NEW SYN runs
        new = tmp_path / name
        env = {"OUT": str(out), "SYN": str(new), **({"PFUNC": pf} if pf else {})}
        r = run_pipe(tmp_path, env, script=str(stub))
        assert r.returncode == 0, r.stdout + r.stderr
        b = [ln for ln in r.stdout.splitlines() if "build_d100.py" in ln and "--tieset" in ln][-1]
        assert f"--pfunc {pf or 'e3'}" in b and f"--out-dir {new} " in b
        assert f"--selection-out {new}/asv_reference_selection_d100.json" in b
        assert f"--out {new}/deployment_stats_d100.json" in r.stdout
        assert f"{out}/asv_reference_selection_d100.json" not in r.stdout
    assert tree_digest(syn) == before


# =========================================================================== BV-BRC PATRIC only
# Since 2026-11-13 every source genome is BV-BRC's own PATRIC annotation only; the fixture
# genome cache above is mixed (PATRIC + RefSeq CDSs, as the cache files written up to 3ce12d0).
def has_patric_id(f):
    return any(a[0] == "PATRIC_id" and a[1] for a in f.get("aliases") or ())


def write_cache(d, sources):
    d.mkdir()
    for g, s_ in sources.items():
        (d / f"{g}.json").write_text(json.dumps(s_))
    return str(d)


def test_bvbrc_only_keeps_patric_features_only():
    g = SOURCES["1.1"]
    before = copy.deepcopy(g)
    f = BSG.bvbrc_only(g)
    assert g == before                                         # input not modified
    assert [x["functions"][0] for x in f["features"]] == ["Prod A", "Prod B", "Prod C", "Prod E"]
    assert all(has_patric_id(x) for x in f["features"])
    assert {k: v for k, v in f.items() if k != "features"} == \
        {k: v for k, v in g.items() if k != "features"}
    assert BSG.bvbrc_only(f) == f                              # idempotent
    # a RefSeq product identical to a PATRIC one is dropped too (2.1's RefSeq 'Prod B')
    assert "Prod B" not in {x["functions"][0] for x in PATRIC_SOURCES["2.1"]["features"]}
    # full-builder dicts: every feature collection; no PATRIC_id alias at all -> dropped
    full = {"id": "9.9",
            "features": [{"id": "a", "aliases": [["PATRIC_id", "fig|9.9.peg.1"], ["gene_name", "x"]]},
                         {"id": "b", "aliases": [["PATRIC_id", ""], ["RefSeq_locus_tag", "T_1"]]},
                         {"id": "c", "aliases": []}, {"id": "d"}],
            "cdss": [{"id": "a_CDS_1", "aliases": [["PATRIC_id", "fig|9.9.peg.1"]]},
                     {"id": "b_CDS_1", "aliases": [["PATRIC_id", ""]]}],
            "non_coding_features": [{"id": "r", "aliases": [["PATRIC_id", "fig|9.9.rna.1"]]},
                                    {"id": "s", "aliases": [["PATRIC_id", ""]]}],
            "mrnas": [], "feature_counts": {"CDS": 2}}
    ff = BSG.bvbrc_only(full)
    assert [x["id"] for x in ff["features"]] == ["a"] and [x["id"] for x in ff["cdss"]] == ["a_CDS_1"]
    assert [x["id"] for x in ff["non_coding_features"]] == ["r"] and ff["mrnas"] == []
    assert ff["feature_counts"] == {"CDS": 2}


def test_fast_build_genome_queries_patric_only_and_filters(world, monkeypatch):
    queries = []
    rows = [{"patric_id": "fig|7.1.peg.1", "product": "Prod A", "feature_type": "CDS",
             "pgfam_id": "PGF_A", "figfam_id": "FIG1"},
            {"patric_id": "fig|7.1.peg.2", "product": "Prod B", "feature_type": "CDS"},
            # a RefSeq row (were the API to ignore the filter): never gets in
            {"patric_id": "", "product": "Prod A", "feature_type": "CDS"},
            {"product": "Prod R", "feature_type": "CDS"}]

    def api(core, query, **kw):
        queries.append((core, query))
        if core == "genome":
            return [{"genome_name": "Org 7.1", "gc_content": 51.0,
                     "taxon_lineage_names": ["cellular organisms", "Bacteria"]}]
        return copy.deepcopy(rows)
    monkeypatch.setattr(BSG, "_api_json", api)
    g = BSG.fast_build_genome("7.1")
    q = dict(queries)["genome_feature"]
    assert "eq(annotation,PATRIC)" in q and "eq(feature_type,CDS)" in q
    assert [f["functions"] for f in g["features"]] == [["Prod A"], ["Prod B"]]
    assert all(has_patric_id(f) for f in g["features"])
    assert g["features"][0]["aliases"] == [["PATRIC_id", "fig|7.1.peg.1"], ["FIGFAM", "FIG1"],
                                           ["PGFAM", "PGF_A"]]
    # the fresh-build path writes the PATRIC-only dict to a NEW cache file
    cache = world["tmp"] / "fresh_cache"
    g2 = BSG.build_genome("7.1", str(cache), "fast")
    assert g2 == g and jl(cache / "7.1.json") == g


def test_cache_read_path_filters_on_read_and_never_rewrites(world, monkeypatch):
    p = os.path.join(world["cache"], "2.1.json")
    raw, mt = rb(p), os.stat(p).st_mtime_ns
    monkeypatch.setattr(BSG, "fast_build_genome", lambda gid: (_ for _ in ()).throw(
        AssertionError("a cached genome must not be fetched")))
    g = BSG.build_genome("2.1", world["cache"], "fast")
    assert g == PATRIC_SOURCES["2.1"] and len(g["features"]) == 4
    assert BSG.load_source_genome(p) == PATRIC_SOURCES["2.1"]
    assert len(jl(p)["features"]) == 7                          # the file still holds RefSeq
    assert rb(p) == raw and os.stat(p).st_mtime_ns == mt        # never rewritten
    # the d100 wrapper's atomic variant: the same filter on both paths
    sys.path.insert(0, os.path.join(ROOT, "deployment_d100"))
    try:
        import build_d100 as D
    finally:
        sys.path.remove(os.path.join(ROOT, "deployment_d100"))
    assert D._build_genome_atomic("2.1", world["cache"], "fast") == PATRIC_SOURCES["2.1"]
    monkeypatch.setattr(BSG, "fast_build_genome", lambda gid: copy.deepcopy(SOURCES["3.1"]))
    cache = world["tmp"] / "atomic_cache"
    g3 = D._build_genome_atomic("3.1", str(cache), "fast")
    assert g3 == PATRIC_SOURCES["3.1"] and jl(cache / "3.1.json") == PATRIC_SOURCES["3.1"]
    assert os.listdir(cache) == ["3.1.json"]


def test_full_builder_keeps_annotation_patric_rows_only(world, monkeypatch):
    KB = BSG.KB
    C = KB.BVBRCToKBaseConverter
    raw = [{"annotation": "PATRIC", "patric_id": "fig|8.1.peg.1", "feature_type": "CDS",
            "product": "Prod A", "pgfam_id": "PGF_A", "go": ["GO:1", "GO:2"]},
           {"annotation": "PATRIC", "patric_id": "fig|8.1.rna.1", "feature_type": "tRNA",
            "product": "tRNA-Ala"},
           {"annotation": "RefSeq", "patric_id": "", "feature_type": "CDS", "product": "Prod A",
            "refseq_locus_tag": "T_1"},
           {"annotation": "RefSeq", "feature_type": "gene", "product": "Prod R"}]
    monkeypatch.setattr(C, "fetch_genome_features", lambda self: copy.deepcopy(raw))
    monkeypatch.setattr(C, "fetch_genome_metadata", lambda self: {
        "genome_name": "Org 8.1", "taxon_lineage_names": ["cellular organisms", "Bacteria"],
        "gc_content": 0.5})
    monkeypatch.setattr(C, "fetch_genome_sequences", lambda self: [])
    monkeypatch.setattr(C, "fetch_feature_sequences", lambda self, md5s: {})
    g = BSG.fetch_source_genome("8.1", "full")        # applies the fetch patch itself
    assert C.fetch_genome_features._prfba_patched
    BSG._patch_go_field()                              # idempotent
    rows = C("8.1").fetch_genome_features()
    assert [r["patric_id"] for r in rows] == ["fig|8.1.peg.1", "fig|8.1.rna.1"]
    assert rows[0]["go"] == "GO:1,GO:2"
    feats = g["features"] + g["non_coding_features"]
    assert sorted(f["functions"][0] for f in feats) == ["Prod A", "tRNA-Ala"]
    assert all(has_patric_id(f) for f in feats + g["cdss"])
    cache = world["tmp"] / "full_cache"
    monkeypatch.setattr(BSG, "fetch_source_genome", lambda gid, b: (_ for _ in ()).throw(
        AssertionError("cached")))
    (cache).mkdir()
    (cache / "8.1.json").write_text(json.dumps({**g, "features": g["features"] + [
        {"id": "x", "functions": ["Prod R"], "aliases": [["PATRIC_id", ""]]}]}))
    assert BSG.build_genome("8.1", str(cache), "full")["features"] == g["features"]


@pytest.mark.parametrize("mode", ["kn", "e3"])
def test_merge_uses_patric_features_only(world, mode):
    """The merge, k/n (recomputed over PATRIC features only), core/accessory, function
    prevalence, E3's PGFam resolution and the digests: a build from the MIXED cache is
    byte-identical to one from a cache holding only the PATRIC CDSs."""
    pure = write_cache(world["tmp"] / "pure_cache", PATRIC_SOURCES)
    ctx = (e3_ctx(world, world["tmp"] / f"ctx_{mode}")[0] if mode == "e3" else None)
    for asv in BUILT:
        em = build_one(world["tmp"] / "mx", asv, ctx, world["cache"])
        ep = build_one(world["tmp"] / "pu", asv, ctx, pure)
        assert em == ep, asv
        for k in KINDS4:
            assert rb(files4(world["tmp"] / "mx", asv)[k]) == rb(files4(world["tmp"] / "pu", asv)[k])
    n3 = jl(files(world["tmp"] / "mx", "asv_n3")["root"])
    assert all(has_patric_id(f) for f in n3["features"])
    p = {f["functions"][0]: f["probability"] for f in n3["features"]}
    assert not any(fn.startswith("Prod R") for fn in p)                  # RefSeq-only products
    kn3 = jl(files(world["tmp"] / "mx", "asv_n3")["provenance"])
    # k/n over PATRIC CDSs: 'Prod B' (PATRIC in 1.1, RefSeq in 2.1 and 3.1) 1/3, not 1;
    # 'Prod F' (PATRIC in 3.1, RefSeq in 2.1) 1/3, not 2/3; 'Prod D' 2/3, not 1
    assert kn3["function_prevalence"] == {"Prod A": 3, "Prod B": 1, "Prod C": 2, "Prod E": 2,
                                          "Prod D": 2, "Prod F": 1}
    n2 = jl(files(world["tmp"] / "mx", "asv_n2")["root"])
    d = next(f for f in n2["features"] if f["functions"] == ["Prod D"])
    # introduced by 2.1's PATRIC CDS (with its PGFAM alias), not by 1.1's earlier RefSeq one
    assert d["aliases"] == [["PATRIC_id", "fig|2.1.peg.3"], ["PGFAM", "PGF_D"]]
    if mode == "kn":
        assert p["Prod B"] == p["Prod F"] == pytest.approx(1 / 3)
        assert p["Prod D"] == pytest.approx(2 / 3) and d["probability"] == 0.5
    else:
        pf = jl(files(world["tmp"] / "mx", "asv_n2")["provenance"])["pfunc"]
        assert pf["feature_resolution"] == {"own_alias": 3, "product_match": 1,
                                            "unresolved": 1, "features": 5}
        assert pf["source_digests"] == P.source_digests(
            ["1.1", "2.1"], [PATRIC_SOURCES["1.1"], PATRIC_SOURCES["2.1"]])
        fpg, _ = P.resolve_feature_pgfams(n2, [PATRIC_SOURCES["1.1"], PATRIC_SOURCES["2.1"]])
        assert fpg[d["id"]] == {"PGF_D"}


def test_every_manifest_records_the_annotation(world, monkeypatch):
    for mode in ("kn", "e3"):
        for parallel in (False, True):
            m = cli(world, world["tmp"] / f"an_{mode}_{parallel}", ["--pfunc", mode, "--pfunc-no-fetch"],
                    monkeypatch, parallel=parallel)
            assert m["_meta"]["annotation"] == BSG.ANNOTATION_STAMP == \
                "BV-BRC PATRIC only (eq(annotation,PATRIC); cached dicts filtered on read)"
    rc, m = retro(world, world["tmp"] / "an_kn_False", world["tmp"] / "an_retro")
    assert rc == 0 and m["_meta"]["annotation"] == BSG.ANNOTATION_STAMP


def test_retro_loads_patric_only_and_refuses_a_mixed_build(world, monkeypatch):
    srcs, missing = R.load_sources(["2.1", "9.9", "3.1"], world["cache"])
    assert srcs == [PATRIC_SOURCES["2.1"], PATRIC_SOURCES["3.1"]] and missing == ["9.9"]
    src = world["tmp"] / "s_kn"
    run_serial(world, src, "kn", monkeypatch)
    old = world["tmp"] / "s_old"                    # a build made at or before 3ce12d0
    shutil.copytree(src, old)
    man = jl(old / MANIFEST)
    man["_meta"].pop("annotation")
    (old / MANIFEST).write_text(json.dumps(man, indent=1))
    before = tree_digest(old)
    for mode in ("e3", "kn"):
        with pytest.raises(SystemExit) as ei:
            retro(world, old, world["tmp"] / f"ro_{mode}", extra=["--pfunc", mode])
        msg = str(ei.value.code)
        assert "not BV-BRC PATRIC-only" in msg and "3ce12d0" in msg and "Nothing was" in msg
        assert not os.path.exists(world["tmp"] / f"ro_{mode}")
    assert tree_digest(old) == before


def test_fetch_pgfams_queries_patric_only(monkeypatch):
    urls = []

    class Resp:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self, *a):
            return json.dumps([{"pgfam_id": "PGF_B"}, {"pgfam_id": "PGF_A"}, {}]).encode()

    def urlopen(req, timeout=None):
        urls.append(req.full_url)
        return Resp()
    monkeypatch.setattr(PF.urllib.request, "urlopen", urlopen)
    assert PF.fetch_pgfams("1.1") == ["PGF_A", "PGF_B"]
    assert "eq(annotation,PATRIC)" in urls[0] and "eq(feature_type,CDS)" in urls[0]
    PF.fetch_pgfams("1.1", family="plfam_id")
    assert "eq(annotation,PATRIC)" in urls[1] and "select(plfam_id)" in urls[1]
