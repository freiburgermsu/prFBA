"""Offline tests for pfunc_estimator (production E3). Synthetic fixtures only: no network,
no taxdump, no shared cache. Run:
    ~/Documents/py_venv/bin/python -m pytest -q tests/test_pfunc_estimator.py
"""
import copy
import json
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pfunc_estimator as P  # noqa: E402
import select_references as S  # noqa: E402

TAU = P.TAU


# --------------------------------------------------------------------------- fixtures
def dec(rank, gid, ident, score=None, cov=1.0, agree=0.0, species=None, stage=None):
    """A select_references decision (identity/coverage/agreement already 4 dp)."""
    return {"rank": rank, "genome_id": gid, "identity": ident,
            "align_score": score if score is not None else ident * 500,
            "coverage": cov, "edlib_agreement": agree, "species": species,
            "family": None, "genus": None, "disposition": "excluded", "stage": stage}


def sel(gid, ident, species=None):
    return {"genome_id": gid, "identity": ident, "species": species, "taxon_id": None}


def hit(rank, gid, ident, score=None, aligned=250, edlib=None, species=None):
    return {"rank": rank, "genome_id": gid, "identity": ident,
            "align_score": score if score is not None else ident * 500,
            "aligned_len": aligned, "edlib_identity": ident if edlib is None else edlib,
            "lineage": {"Species": species}}


def w(id1, x):
    return math.exp(-(id1 - x) / TAU)


# =========================================================================== constants
def test_constants_match_stage_e_and_select_references():
    assert (P.K, P.TAU, P.ALPHA, P.MIN_CLADE, P.MIN_N) == (100, 0.02, 2.0, 3, 2)
    assert P.COV_MIN == S.DEFAULT_KNOBS["cov_min"] == 0.90
    assert P.SWED_MAX == S.DEFAULT_KNOBS["swed_max"] == 0.02
    assert P.FAMILY_FLOOR == S.DEFAULT_KNOBS["family_floor"] == 0.865


# =========================================================================== stage 0
def test_stage0_decisions_filters():
    rec = {"selected": [], "decisions": [
        dec(1, "1.1", 0.99),                              # passes
        dec(2, "2.1", 0.98, cov=0.8999),                  # coverage < 0.90
        dec(3, "3.1", 0.97, agree=0.0201),                # |SW-edlib| > 0.02
        dec(4, "4.1", 0.96, agree=None),                  # None agreement passes
        dec(5, "5.1", 0.8649),                            # below family floor
        dec(6, "6.1", 0.865, stage="identity_floor"),     # 4-dp rounding edge: raw < floor
        dec(7, "7.1", 0.865),                             # exactly at floor passes
        dec(8, "8.1", 0.95, cov=0.90, agree=0.02),        # boundaries pass
    ]}
    info = {}
    nb = P.neighbours_from_decisions(rec, dedup="genome", reps_first=False, info=info)
    assert [g for g, _ in nb] == ["1.1", "4.1", "8.1", "7.1"]
    assert info["n_stage0"] == 4


def test_stage0_hits_matches_pe_rounding():
    L = 1000
    rec = {"asv_len": L, "top20": [
        hit(1, "1.1", 0.99, aligned=900),                       # cov 0.9 passes
        hit(2, "2.1", 0.99, aligned=899),                       # cov 0.899 fails
        hit(3, "3.1", 0.98, aligned=1000, edlib=0.96004),       # agree round(0.01996,4)=0.02 passes
        hit(4, "4.1", 0.98, aligned=1000, edlib=0.9594),        # agree 0.0206 fails
        hit(5, "5.1", 0.97, aligned=1200),                      # cov capped at 1.0
        hit(6, "6.1", 0.86499, aligned=1000),                   # raw identity below floor
    ]}
    rec["top20"][2]["edlib_identity"] = 0.96004
    rec["top20"][3]["edlib_identity"] = 0.9594
    nb = P.neighbours_from_hits(rec, dedup="genome", reps_first=False)
    assert [g for g, _ in nb] == ["1.1", "3.1", "5.1"]
    # coverage rounding: 899.96/1000 rounds to 0.9 and passes (PE: round(min(.,1),4))
    rec2 = {"asv_len": 100000, "top20": [hit(1, "9.1", 0.99, aligned=89996)]}
    assert [g for g, _ in P.neighbours_from_hits(rec2, reps_first=False)] == ["9.1"]
    # edlib_identity None passes
    rec3 = {"asv_len": 100, "top20": [hit(1, "9.2", 0.99, aligned=100)]}
    rec3["top20"][0]["edlib_identity"] = None
    assert [g for g, _ in P.neighbours_from_hits(rec3, reps_first=False)] == ["9.2"]


def test_hits_used_first_100_only():
    top = [hit(i + 1, f"{i + 1}.1", 0.99 - 0.0001 * i) for i in range(150)]
    nb = P.neighbours_from_hits({"asv_len": 250, "top20": top}, reps_first=False)
    assert len(nb) == 100 and nb[-1][0] == "100.1"


# =========================================================================== sort / dedup
def test_sort_by_identity_then_score_stable():
    rec = {"selected": [], "decisions": [
        dec(1, "a.1", 0.97, score=400),
        dec(2, "b.1", 0.98, score=300),
        dec(3, "c.1", 0.97, score=450),
        dec(4, "d.1", 0.97, score=400),   # exact tie with a.1 -> keeps rank order
        dec(5, "e.1", 0.98, score=300),   # exact tie with b.1 -> keeps rank order
    ]}
    nb = P.neighbours_from_decisions(rec, dedup="genome", reps_first=False)
    assert [g for g, _ in nb] == ["b.1", "e.1", "c.1", "a.1", "d.1"]


def test_genome_vs_species_dedup():
    decs = [
        dec(1, "10.1", 0.99, species="Escherichia coli"),
        dec(2, "10.1", 0.985, species="Escherichia coli"),   # multi-copy 16S, same genome
        dec(3, "10.2", 0.985, species="escherichia coli"),   # same species (case-folded)
        dec(4, "20.1", 0.98, species="Candidatus Foo bar"),
        dec(5, "20.2", 0.975, species="Foo bar"),            # same after prefix strip
        dec(6, "30.1", 0.97, species=""),                    # blank -> genome:30.1
        dec(7, "30.2", 0.97, species=None),                  # blank -> genome:30.2
    ]
    rec = {"selected": [], "decisions": decs}
    g = P.neighbours_from_decisions(rec, dedup="genome", reps_first=False)
    s = P.neighbours_from_decisions(rec, dedup="species", reps_first=False)
    assert [x for x, _ in g] == ["10.1", "10.2", "20.1", "20.2", "30.1", "30.2"]
    assert [x for x, _ in s] == ["10.1", "20.1", "30.1", "30.2"]
    # the first (highest-identity) copy of a genome supplies its identity
    assert dict(g)["10.1"] == pytest.approx(1.0)
    assert P.species_key("", "7.1") == "genome:7.1"
    assert P.species_key("Ca. Foo bar", "7.1") == "foo bar"


def test_reps_first_and_skip_rep_species():
    decs = [
        dec(1, "11.1", 0.99, score=500, species="Sp a"),
        dec(2, "12.1", 0.99, score=510, species="Sp a"),   # same species as rep 11.1
        dec(3, "13.1", 0.99, score=520, species="Sp b"),
        dec(4, "14.1", 0.97, species="Sp c"),
    ]
    rec = {"selected": [sel("13.1", 0.99, "Sp b"), sel("11.1", 0.99, "Sp a")], "decisions": decs}
    sp = P.neighbours_from_decisions(rec, dedup="species")
    assert [g for g, _ in sp] == ["13.1", "11.1", "14.1"]       # 12.1 skipped (rep species)
    gm = P.neighbours_from_decisions(rec, dedup="genome")
    assert [g for g, _ in gm] == ["13.1", "11.1", "12.1", "14.1"]
    pure = P.neighbours_from_decisions(rec, dedup="genome", reps_first=False)
    assert [g for g, _ in pure] == ["13.1", "12.1", "11.1", "14.1"]
    # q is order-free: same multiset of (gid, w)
    assert sorted(gm) == sorted(pure)


def test_K_cap_counts_reps():
    decs = [dec(i + 1, f"{i + 1}.1", 0.99 - 0.0001 * i) for i in range(120)]
    rec = {"selected": [sel("1.1", 0.99), sel("2.1", 0.99)], "decisions": decs}
    nb = P.neighbours_from_decisions(rec, dedup="genome", K=10)
    assert len(nb) == 10 and [g for g, _ in nb[:2]] == ["1.1", "2.1"]
    assert len(P.neighbours_from_decisions(rec, dedup="genome")) == 100


# =========================================================================== weights / fallback
def test_weights_exp_decay_and_reps_at_one():
    decs = [dec(1, "1.1", 0.99, species="A a"), dec(2, "2.1", 0.97, species="B b"),
            dec(3, "3.1", 0.95, species="C c")]
    rec = {"selected": [sel("1.1", 0.99, "A a")], "decisions": decs}
    nb = P.neighbours_from_decisions(rec)
    assert nb[0] == ("1.1", 1.0)
    assert nb[1][1] == pytest.approx(math.exp(-1.0))          # 0.02 below the best
    assert nb[2][1] == pytest.approx(math.exp(-2.0))
    nb5 = P.neighbours_from_decisions(rec, tau=0.005)
    assert nb5[1][1] == pytest.approx(math.exp(-4.0))


def test_fallback_to_reps_when_no_survivor():
    rec = {"asv_len": 250, "top20": [hit(1, "1.1", 0.99, aligned=100),     # coverage fails
                                     hit(2, "2.1", 0.80)]}                  # floor fails
    info = {}
    nb = P.neighbours_from_hits(rec, reps=["5.1", "6.1"], reps_first=False, info=info)
    assert nb == [("5.1", 1.0), ("6.1", 1.0)] and info["fallback_reps"] is True
    assert P.neighbours_from_hits(rec, reps=None, reps_first=False) == []
    # decisions: no survivor and no reps -> []
    assert P.neighbours_from_decisions({"selected": [], "decisions": [dec(1, "1.1", 0.80)]}) == []


def test_hits_and_decisions_agree_through_select_representatives():
    """End to end: a synthetic hit record -> select_references -> neighbours_from_decisions
    equals neighbours_from_hits on the raw record with the same reps (4-dp identities)."""
    sps = ["Alpha one", "Alpha one", "Beta two", "", "Gamma three", "Delta four"]
    top = []
    for i, (g, x) in enumerate([("1.1", 0.9921), ("1.2", 0.9921), ("2.1", 0.9921),
                                ("3.1", 0.9802), ("4.1", 0.9711), ("5.1", 0.9500)]):
        top.append({**hit(i + 1, g, x, score=round(x * 500)), "taxon_id": int(g.split(".")[0]),
                    "organism": sps[i] or "bacterium", "aligned_len": 250,
                    "lineage": {"Species": sps[i], "Family": "Famx"}})
    hits_rec = {"asv_len": 250, "best_identity": 0.9921, "top20": top}
    audit = S.select_representatives(hits_rec)
    reps = [s["genome_id"] for s in audit["selected"]]
    assert set(reps) <= {"1.1", "1.2", "2.1"} and len(reps) == 2
    for mode in ("species", "genome"):
        a = P.neighbours_from_decisions(audit, dedup=mode)
        b = P.neighbours_from_hits(hits_rec, reps=audit["selected"], dedup=mode)
        assert a == b


# =========================================================================== e3_probabilities
STORE = {
    "n1": {"P1", "P2"}, "n2": {"P1"}, "n3": {"P3"},
    "m1": {"P1", "P3"}, "m2": {"P1"}, "m3": {"P2"}, "m4": set(),
}


def test_e3_formula_with_clade():
    nb = [("n1", 1.0), ("n2", 0.5), ("n3", 0.25)]
    feats = {"f1": {"P1"}, "f2": {"P2"}, "f3": {"P2", "P3"}, "f4": {"P9"}, "f5": set()}
    q, d = P.e3_probabilities(feats, nb, ["m1", "m2", "m3"], STORE, alpha=2.0)
    den = 1.75 + 2.0
    assert q["f1"] == pytest.approx((1.5 + 2.0 * 2 / 3) / den)
    assert q["f2"] == pytest.approx((1.0 + 2.0 * 1 / 3) / den)
    assert q["f3"] == pytest.approx((1.25 + 2.0 * 2 / 3) / den)   # any of the set (m1, m3)
    assert q["f4"] == 0.0                                          # clamped (exact 0)
    assert "f5" not in q                                           # empty set -> keeps k/n
    assert d["clade_used"] and d["n_members_used"] == 3 and d["n_features_estimated"] == 4


def test_e3_no_clade_drops_alpha():
    nb = [("n1", 1.0), ("n2", 0.5), ("n3", 0.25)]
    q, d = P.e3_probabilities({"f1": {"P1"}}, nb, [], STORE, alpha=2.0)
    assert q["f1"] == pytest.approx(1.5 / 1.75)
    assert not d["clade_used"] and d["alpha_effective"] == 0.0


def test_e3_clamping_to_one_and_zero():
    nb = [("n1", 1.0), ("n2", 1.0)]
    q, _ = P.e3_probabilities({"a": {"P1"}}, nb, [], STORE)
    assert q["a"] == 1.0
    # numerator 1e-13 of 1 -> clamped to 0
    nb2 = [("n3", 1.0), ("n2", 1e-13)]
    q2, _ = P.e3_probabilities({"a": {"P1"}}, nb2, [], STORE)
    assert q2["a"] == 0.0
    # just under 1-1e-9 -> 1.0
    nb3 = [("n2", 1.0), ("n3", 5e-10)]
    q3, _ = P.e3_probabilities({"a": {"P1"}}, nb3, [], STORE)
    assert q3["a"] == 1.0


def test_e3_missing_and_empty_are_excluded_and_counted():
    nb = [("n1", 1.0), ("zz", 0.7), ("m4", 0.3)]           # zz absent, m4 empty set
    q, d = P.e3_probabilities({"f": {"P1"}}, nb, ["m1", "m2", "absent", "m4"], STORE)
    assert d["n_neighbours_missing"] == 1 and d["n_neighbours_empty"] == 1
    assert d["weight_missing"] == pytest.approx(1.0)
    assert d["n_members_missing"] == 1 and d["n_members_empty"] == 1
    # only 2 usable members < MIN_CLADE -> prior dropped
    assert not d["clade_used"] and d["clade_dropped_after_missing"]
    assert q["f"] == 1.0
    # no usable neighbour and no clade -> no estimate
    q0, d0 = P.e3_probabilities({"f": {"P1"}}, [("zz", 1.0)], [], STORE)
    assert q0 == {} and d0["no_estimate"]
    # members only (no usable neighbour): q = pi
    q1, _ = P.e3_probabilities({"f": {"P1"}}, [("zz", 1.0)], ["m1", "m2", "m3"], STORE)
    assert q1["f"] == pytest.approx(2 / 3)


def test_e3_accepts_pgfamstore():
    st = P.PGFamStore(STORE)
    q, _ = P.e3_probabilities({"f": {"P1"}}, [("n1", 1.0), ("n3", 1.0)], [], st)
    assert q["f"] == pytest.approx(0.5)


# =========================================================================== clade catalogue
def make_catalogue():
    # genome -> (species_taxid, family, order)
    g = {
        "100.1": (100, "FamA", "OrdX"), "100.2": (100, "FamA", "OrdX"),
        "100.10": (100, "FamA", "OrdX"),
        "200.1": (200, "FamA", "OrdX"), "300.5": (300, "FamA", "OrdX"),
        "400.1": (400, "FamB", "OrdX"), "500.1": (500, "FamB", "OrdX"),
        "600.1": (600, None, "OrdY"), "700.1": (700, "FamC", "OrdY"),
        "800.1": (None, "FamA", "OrdX"),                 # no species rank: never a member
        "900.1": (None, None, None),                     # unresolved taxid
    }
    taxa = {150: (150, "FamA", "OrdX")}
    return P.CladeCatalogue(g, taxa)


def test_clade_family_one_genome_per_species_smallest_id():
    rank, mem = P.clade_members("100.2", make_catalogue())
    assert rank == "family"
    assert sorted(mem) == ["100.1", "200.1", "300.5"]   # '100.1' < '100.10' < '100.2'
    assert "800.1" not in mem


def test_clade_backoff_to_order_and_none():
    cat = make_catalogue()
    rank, mem = P.clade_members("400.1", cat)        # FamB has 2 species -> order OrdX
    assert rank == "order" and len(mem) == 5         # 100, 200, 300, 400, 500
    rank, mem = P.clade_members("600.1", cat)        # family None, OrdY has 2 species
    assert (rank, mem) == (None, [])
    assert P.clade_members("900.1", cat) == (None, [])
    assert P.clade_members("12345.6", cat) == (None, [])     # unknown anchor, no taxdb


def test_clade_anchor_outside_cache_uses_taxa():
    rank, mem = P.clade_members("150.9", make_catalogue())
    assert rank == "family" and len(mem) == 3


def test_clade_exclusions_and_cap():
    cat = make_catalogue()
    rank, mem = P.clade_members("100.1", cat, exclude_genomes={"100.1"})
    assert rank == "family" and "100.10" in mem              # next smallest of species 100
    rank, mem = P.clade_members("100.1", cat, exclude_species={100})
    assert rank == "order"                                   # FamA left with 200, 300 only
    assert not any(g.startswith("100.") for g in mem)
    # cap: the species with the smallest sha1(str(taxid)) are kept
    import hashlib
    rank, mem = P.clade_members("100.1", cat, max_species=2, min_clade=2)
    want = sorted([100, 200, 300], key=lambda s: hashlib.sha1(str(s).encode()).hexdigest())[:2]
    assert sorted(cat.genomes[g][0] for g in mem) == sorted(want)


def test_build_clade_catalogue_offline(tmp_path, monkeypatch):
    cache = tmp_path / "cache.json"
    cache.write_text(json.dumps({"100.1": ["PGF_1"], "100.2": [], "abc": ["PGF_2"],
                                 "200.3": ["PGF_3", "PGF_4"]}))
    fake = {100: (100, "FamA", "OrdX"), 200: (None, None, None)}
    monkeypatch.setattr(P, "resolve_taxid", lambda t, db: fake.get(t, (None, None, None)))
    out = tmp_path / "cat.json"
    meta = P.build_clade_catalogue(str(cache), taxdb=object(), out_path=str(out), log=None)
    obj = json.loads(out.read_text())
    assert obj["genomes"] == {"100.1": [100, "FamA", "OrdX"], "abc": [None, None, None],
                              "200.3": [None, None, None]}
    assert meta["n_empty_pgfam_lists_excluded"] == 1 and meta["empty_pgfam_genomes"] == ["100.2"]
    assert meta["n_genome_id_without_numeric_taxid"] == 1
    cat = P.load_clade_catalogue(str(out))
    assert cat.ranks_for("100.7") == (100, "FamA", "OrdX")     # via taxa
    assert cache.read_text().startswith('{"100.1"')             # cache untouched
    assert P.catalogue_is_current(cat, str(cache))
    os.utime(cache, (1, 1))
    assert not P.catalogue_is_current(cat, str(cache))


# =========================================================================== PGFam store
def test_pgfamstore_ensure_side_fetch_save_load(tmp_path, monkeypatch):
    cache = tmp_path / "cache.json"
    cache.write_text(json.dumps({"1.1": ["PGF_a", "PGF_b"], "2.1": ["PGF_a"], "3.1": []}))
    side = tmp_path / "side.json"
    side.write_text(json.dumps({"4.1": ["PGF_c"]}))
    import prefetch_gene_families as PF

    def fake_fetch(g, **kw):
        if g == "5.1":
            return ["PGF_d"]
        raise RuntimeError("boom")
    monkeypatch.setattr(PF, "fetch_pgfams", fake_fetch)
    st = P.PGFamStore()
    stats = st.ensure(["1.1", "3.1", "4.1", "5.1", "6.1"], cache_path=str(cache),
                      side_path=str(side), fetch=True, workers=2, log=None)
    assert stats["from_side_file"] == 1 and stats["from_shared_cache"] == 2
    assert stats["fetched"] == 1 and stats["fetch_failed"] == 1 and stats["absent"] == 1
    assert st.get("1.1") == {"PGF_a", "PGF_b"} and st.get("3.1") == frozenset()
    assert st.get("6.1") is None and st.missing == {"6.1"}
    # the failure is never stored; the success is merged into the side file
    assert json.loads(side.read_text()) == {"4.1": ["PGF_c"], "5.1": ["PGF_d"]}
    # interned strings shared across genomes
    a1 = next(x for x in st.get("1.1") if x == "PGF_a")
    st.ensure(["2.1"], cache_path=str(cache), log=None)
    a2 = next(iter(st.get("2.1")))
    assert a1 is a2
    slim = tmp_path / "slim.json"
    st.save(str(slim))
    st2 = P.PGFamStore.load(str(slim))
    assert st2.get("5.1") == {"PGF_d"} and len(st2) == len(st)
    with pytest.raises(ValueError):
        P.PGFamStore().ensure(["9.9"], cache_path=None, fetch=True, log=None)


# =========================================================================== genome objects
def make_genome():
    feats, cdss = [], []
    spec = [("asv_1", "Prod A", [["PATRIC_id", "fig|1"], ["PGFAM", "PGF_A"]], 1.0),
            ("asv_2", "Prod B", [["PATRIC_id", ""]], 0.5),           # RefSeq-introduced
            ("asv_3", "Prod C", [["PATRIC_id", ""]], 0.5),           # no PGFam anywhere
            ("asv_4", "Prod D", [["PGFAM", "PGF_D1"], ["PGFAM", "PGF_D2"]], 1.0)]
    for fid, fn, al, p in spec:
        f = {"id": fid, "type": "gene", "functions": [fn], "aliases": al, "cdss": [f"{fid}.CDS"],
             "probability": p}
        feats.append(f)
        c = {k: copy.deepcopy(v) for k, v in f.items() if k != "cdss"}
        c.update(id=f"{fid}.CDS", type="CDS", parent_gene=fid)
        cdss.append(c)
    return {"id": "asv", "features": feats, "cdss": cdss}


def make_sources():
    return [
        {"id": "1.1", "features": [
            {"functions": ["Prod A"], "aliases": [["PGFAM", "PGF_A"]], "pgfam_id": "PGF_A"},
            {"functions": ["Prod B"], "aliases": [["PATRIC_id", ""]]}]},
        {"id": "2.1", "features": [
            {"functions": ["Prod B"], "aliases": [["PATRIC_id", "fig|2"]], "pgfam_id": "PGF_B1"},
            {"functions": ["Prod B"], "aliases": [["PGFAM", "PGF_B2"]]},    # alias only
            {"functions": ["Prod C"], "aliases": [["PATRIC_id", ""]]}]},
    ]


def test_resolve_feature_pgfams_routes():
    fpg, counts = P.resolve_feature_pgfams(make_genome(), make_sources())
    assert fpg == {"asv_1": {"PGF_A"}, "asv_2": {"PGF_B1", "PGF_B2"},
                   "asv_4": {"PGF_D1", "PGF_D2"}}
    assert counts == {"own_alias": 2, "product_match": 1, "unresolved": 1, "features": 4}


def test_apply_to_genome_updates_features_and_cds_only():
    g = make_genome()
    before = copy.deepcopy(g)
    old = P.apply_to_genome(g, {"asv_2": 0.25, "asv_4": 0.75})
    assert old == {"asv_2": 0.5, "asv_4": 1.0}
    fp = {f["id"]: f["probability"] for f in g["features"]}
    cp = {c["id"]: c["probability"] for c in g["cdss"]}
    assert fp == {"asv_1": 1.0, "asv_2": 0.25, "asv_3": 0.5, "asv_4": 0.75}
    assert cp == {"asv_1.CDS": 1.0, "asv_2.CDS": 0.25, "asv_3.CDS": 0.5, "asv_4.CDS": 0.75}
    assert len(g["features"]) == len(before["features"]) and len(g["cdss"]) == len(before["cdss"])
    # nothing but 'probability' changed
    for a, b in zip(g["features"] + g["cdss"], before["features"] + before["cdss"]):
        assert {k: v for k, v in a.items() if k != "probability"} == \
               {k: v for k, v in b.items() if k != "probability"}
    with pytest.raises(KeyError):
        P.apply_to_genome(g, {"nope": 0.1})


# =========================================================================== reweight_asv
def asv_fixture():
    rec = {"selected": [sel("1.1", 0.99, "Sp one"), sel("2.1", 0.99, "Sp two")],
           "decisions": [dec(1, "1.1", 0.99, species="Sp one"),
                         dec(2, "2.1", 0.99, species="Sp two"),
                         dec(3, "3.1", 0.97, species="Sp three")]}
    store = {"1.1": {"PGF_A", "PGF_B1"}, "2.1": {"PGF_A", "PGF_B1", "PGF_D1"},
             "3.1": {"PGF_A"}, "m1": {"PGF_A"}, "m2": {"PGF_D2"}, "m3": {"PGF_X"}}
    cat = P.CladeCatalogue({"1.1": (1, "Fam", "Ord"), "m1": (11, "Fam", "Ord"),
                            "m2": (12, "Fam", "Ord"), "m3": (13, "Fam", "Ord")})
    return rec, store, cat


def test_reweight_asv_kn_mode_and_n1_are_noops():
    rec, store, cat = asv_fixture()
    for kw in ({"mode": "kn", "n_built": 2}, {"mode": "e3", "n_built": 1}):
        g = make_genome()
        before = copy.deepcopy(g)
        diag, kn = P.reweight_asv(g, rec, make_sources(), store, cat, **kw)
        assert kn == {} and diag["applied"] is False and g == before
    g = make_genome()
    diag, kn = P.reweight_asv(g, {"selected": [], "decisions": []}, make_sources(), store, cat,
                              n_built=3)
    assert kn == {} and not diag["applied"]


def test_reweight_asv_e3_applied_at_n2():
    rec, store, cat = asv_fixture()
    g = make_genome()
    diag, kn = P.reweight_asv(g, rec, make_sources(), store, cat, n_built=2)
    assert diag["applied"] and diag["clade_rank"] == "family"
    assert diag["anchor_in_catalogue"] and diag["anchor_family"] == "Fam"
    assert kn == {"asv_1": 1.0, "asv_2": 0.5, "asv_4": 1.0}       # asv_3 unresolved: kept
    w3 = math.exp(-1.0)
    den = 2 + w3 + 2.0
    # members m1 (11), m2 (12), m3 (13): the clade also contains 1.1 (species 1) -> 4 species
    mem = P.clade_members("1.1", cat)[1]
    assert sorted(mem) == ["1.1", "m1", "m2", "m3"]
    pi = lambda S_: sum(1 for m in mem if store[m] & S_) / len(mem)   # noqa: E731
    exp_q = {"asv_1": (2 + w3 + 2 * pi({"PGF_A"})) / den,
             "asv_2": (2 + 2 * pi({"PGF_B1", "PGF_B2"})) / den,
             "asv_4": (1 + 2 * pi({"PGF_D1", "PGF_D2"})) / den}
    fp = {f["id"]: f["probability"] for f in g["features"]}
    cp = {c["id"]: c["probability"] for c in g["cdss"]}
    for fid, q in exp_q.items():
        assert fp[fid] == pytest.approx(q) and cp[f"{fid}.CDS"] == pytest.approx(q)
    assert fp["asv_3"] == 0.5 and cp["asv_3.CDS"] == 0.5
    assert diag["feature_resolution"]["unresolved"] == 1
    assert P.required_genomes(rec, cat) == {"1.1", "2.1", "3.1", "m1", "m2", "m3"}


def test_estimator_metadata():
    m = P.estimator_metadata()
    assert m["estimator"] == "E3" and m["neighbour_dedup"] == "species" and m["K"] == 100
    assert P.estimator_metadata(mode="kn")["estimator"] == "kn"


# =========================================================================== fix round 2
def test_decisions_capped_at_hits_used_by_rank():
    """Only the first HITS_USED (100) decisions by rank enter (PE top[:100])."""
    decs = [dec(i, f"{i}.1", round(0.95 - i * 1e-4, 4), species=f"Sp {i}") for i in range(1, 101)]
    decs.append(dec(101, "101.1", 0.999, species="Sp 101"))   # best identity, rank 101
    rec = {"selected": [], "decisions": list(reversed(decs))}   # stored order must not matter
    ids = [g for g, _ in P.neighbours_from_decisions(rec)]
    assert "101.1" not in ids and len(ids) == 100 and ids[0] == "1.1"
    ids101 = [g for g, _ in P.neighbours_from_decisions(rec, hits_used=101)]
    assert ids101[0] == "101.1" and len(ids101) == 100            # K = 100 still caps
    assert P.HITS_USED == 100


def test_estimator_metadata_kn_is_minimal():
    assert P.estimator_metadata(mode="kn") == {
        "estimator": "kn", "estimator_id": "k/n", "gating_rule": "k/n everywhere (--pfunc kn)"}
    m = P.estimator_metadata()
    assert m["hits_used"] == 100 and "source_digest" in m and "prior_anchor" in m


def test_source_digest_order_insensitive_and_sensitive_to_function_and_alias():
    s1, s2 = make_sources()[:2]
    d = P.source_digest(s2)
    rev = {**s2, "features": list(reversed(s2["features"]))}
    assert P.source_digest(rev) == d
    fn = copy.deepcopy(s2)
    fn["features"][2]["functions"] = ["prod c"]                     # product string re-fetched
    assert P.source_digest(fn) != d
    al = copy.deepcopy(s2)
    al["features"][1]["aliases"] = [["PGFAM", "PGF_B9"]]            # alias alone changed
    assert P.source_digest(al) != d
    top = copy.deepcopy(s2)
    top["features"][0]["pgfam_id"] = "PGF_B7"                       # top-level pgfam_id counts
    assert P.source_digest(top) != d
    # no PGFam -> the ('function', '') pair; still a stable hex digest
    assert len(P.source_digest({"features": [{"functions": ["x"]}]})) == 64
    assert P.source_digests(["1.1", "2.1"], [s1, s2]) == {"1.1": P.source_digest(s1), "2.1": d}


def anchor_fixture():
    """selected 1.1 (FamA), 2.1 (FamB), 3.1; FamA and FamB each have 3 other member species."""
    rec = {"selected": [sel("1.1", 0.99, "Sp one"), sel("2.1", 0.99, "Sp two"),
                        sel("3.1", 0.99, "Sp three")],
           "decisions": [dec(1, "1.1", 0.99, species="Sp one"),
                         dec(2, "2.1", 0.99, species="Sp two"),
                         dec(3, "3.1", 0.99, species="Sp three")]}
    store = {"1.1": {"PGF_A"}, "2.1": {"PGF_A", "PGF_D1"}, "3.1": {"PGF_A"},
             "a1": {"PGF_A"}, "a2": {"PGF_A"}, "a3": {"PGF_A"},
             "b1": {"PGF_D1"}, "b2": {"PGF_D1"}, "b3": {"PGF_X"}}
    cat = P.CladeCatalogue({"1.1": (1, "FamA", "Ord"), "a1": (11, "FamA", "Ord"),
                            "a2": (12, "FamA", "Ord"), "a3": (13, "FamA", "Ord"),
                            "2.1": (2, "FamB", "Ord"), "b1": (21, "FamB", "Ord"),
                            "b2": (22, "FamB", "Ord"), "b3": (23, "FamB", "Ord")})
    return rec, store, cat


def test_prior_anchored_on_first_built_rep_failed_rep_stays_neighbour():
    rec, store, cat = anchor_fixture()
    assert P.prior_anchor(rec) == "1.1" and P.prior_anchor(rec, ["2.1", "3.1"]) == "2.1"
    assert P.prior_anchor(rec, ["3.1"]) == "3.1" and P.prior_anchor({"selected": []}) is None
    srcs = [{"id": "2.1", "features": []}, {"id": "3.1", "features": []}]
    g = make_genome()
    diag, kn = P.reweight_asv(g, rec, srcs, store, cat, n_built=2, built_ids=["2.1", "3.1"])
    assert diag["applied"] and diag["anchor"] == "2.1" and diag["anchor_selected"] == "1.1"
    assert diag["anchor_built"] is False and diag["reps_not_built"] == ["1.1"]
    assert diag["anchor_family"] == "FamB" and diag["clade_rank"] == "family"
    assert diag["n_neighbours"] == 3 and diag["n_neighbours_used"] == 3   # 1.1 still a neighbour
    # PGF_A: neighbours 1.1, 2.1, 3.1 all carry it (w 1); FamB members 2.1, b1, b2, b3 -> 1/4
    fp = {f["id"]: f["probability"] for f in g["features"]}
    assert fp["asv_1"] == pytest.approx((3 + 2 * 1 / 4) / (3 + 2))
    assert diag["source_digests"] == P.source_digests(["2.1", "3.1"], srcs)
    # without built ids: the old rule (anchor = selected[0]), anchor_built unknown
    d0, _ = P.reweight_asv(make_genome(), rec, srcs, store, cat, n_built=2)
    assert d0["anchor"] == "1.1" and d0["anchor_built"] is None and d0["reps_not_built"] is None
    d1, _ = P.reweight_asv(make_genome(), rec, srcs, store, cat, n_built=3,
                           built_ids=["1.1", "2.1", "3.1"])
    assert d1["anchor"] == "1.1" and d1["anchor_built"] is True and d1["reps_not_built"] == []


def test_required_genomes_covers_every_anchor_candidate():
    rec, store, cat = anchor_fixture()
    base = {"1.1", "2.1", "3.1"}
    assert P.required_genomes(rec, cat) == base | {"a1", "a2", "a3"}
    assert P.required_genomes(rec, cat, anchors=["2.1"]) == base | {"b1", "b2", "b3"}
    assert P.required_genomes(rec, cat, anchors=["1.1", "2.1"]) == \
        base | {"a1", "a2", "a3", "b1", "b2", "b3"}


def test_no_feature_resolved_reason():
    rec, store, cat = asv_fixture()
    g = {"id": "asv", "features": [{"id": "x", "functions": ["Nothing"], "aliases": [],
                                    "probability": 0.5}], "cdss": []}
    before = copy.deepcopy(g)
    diag, kn = P.reweight_asv(g, rec, [{"id": "1.1", "features": []},
                                       {"id": "2.1", "features": []}], store, cat, n_built=2,
                              built_ids=["1.1", "2.1"])
    assert not diag["applied"] and diag["reason"] == "no feature resolved to a PGFam"
    assert diag["n_neighbours_used"] == 3 and not diag["no_estimate"] and kn == {} and g == before
    # no usable neighbour and no clade -> the no-estimate reason
    g2 = make_genome()
    d2, _ = P.reweight_asv(g2, rec, make_sources(), {}, P.CladeCatalogue({}), n_built=2)
    assert not d2["applied"] and d2["reason"].startswith("no estimate")


def test_catalogue_meta_currency_and_rebuild_cli(tmp_path, monkeypatch):
    cache = tmp_path / "cache.json"
    cache.write_text(json.dumps({"100.1": ["PGF_1"]}))
    monkeypatch.setattr(P, "resolve_taxid", lambda t, db: (100, "FamA", "OrdX"))
    monkeypatch.setattr(P, "load_taxdb", lambda *a, **k: object())
    out = tmp_path / "cat.json"
    assert P.main(["build-catalogue", "--cache", str(cache), "--out", str(out)]) == 0
    meta = P.catalogue_meta(str(out))
    assert meta["n_genomes"] == 1 and P.catalogue_meta_is_current(meta, str(cache))
    assert P.load_clade_catalogue(str(out)).path == str(out)
    os.utime(cache, (1, 1))
    assert not P.catalogue_meta_is_current(meta, str(cache))
    assert not P.catalogue_meta_is_current(meta, str(tmp_path / "absent.json"))
    cmd = P.catalogue_rebuild_command(str(cache), str(out))
    assert "systemd-run --user" in cmd and "MemoryMax" in cmd and "build-catalogue" in cmd
    assert f"--cache {cache}" in cmd and f"--out {out}" in cmd
    with pytest.raises(SystemExit):
        P.main(["build-catalogue", "--cache", str(cache), "--out", str(cache)])
