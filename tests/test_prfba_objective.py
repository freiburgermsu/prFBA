#!/usr/bin/env python
"""
test_prfba_objective.py — checks prfba_objective: the production reaction-probability rule,
the keyword arguments passed to MSProbability.prFBA, and the error on a ModelSEEDpy whose
prFBA lacks the fix branch's keywords.

Offline: the model is ModelSEEDpy's tests/test_data/e_coli_core.json, renamed as in the
fix branch's tests/core/test_msprobability_prfba.py (compartment suffix "0"; internal
reactions "rxnNNNNN_c0", so prFBA prices them; biomass and ATPM keep their ids).

    PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-prfba \\
        ~/Documents/py_venv/bin/python -m pytest -p no:cacheprovider tests/test_prfba_objective.py

Without the fix branch on the path, the tests that solve are skipped (the rule and the
keyword tests still run).

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import inspect
import json
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# modelseedpy.core.msprobability imports both at module level
pytest.importorskip("mscommunity")
pytest.importorskip("cobrakbase")

import cobra  # noqa: E402
import modelseedpy  # noqa: E402
from modelseedpy.core.msprobability import MSProbability  # noqa: E402

import prfba_objective as PO  # noqa: E402

BIOMASS = "BIOMASS_Ecoli_core_w_GAM_c0"
FIXED = all(k in inspect.signature(MSProbability.prFBA).parameters
            for k in PO.REQUIRED_KEYWORDS)
needs_fix = pytest.mark.skipif(not FIXED, reason="the imported ModelSEEDpy ({}) lacks the "
                               "fix branch; run with PYTHONPATH={}".format(
                                   modelseedpy.__file__, PO.FIX_WORKTREE))


def _fixture_path():
    for root in (os.path.dirname(os.path.dirname(os.path.abspath(modelseedpy.__file__))),
                 PO.FIX_WORKTREE):
        p = os.path.join(root, "tests", "test_data", "e_coli_core.json")
        if os.path.exists(p):
            return p
    pytest.skip("e_coli_core.json not found under modelseedpy or the fix worktree")


def _build():
    """e_coli_core renamed as the fix branch's test does; no probability notes (the tests
    write them)."""
    with open(_fixture_path()) as fh:
        model_json = json.load(fh)
    model_json["compartments"] = {k + "0": v for k, v in model_json["compartments"].items()}
    metabolites = {}
    for m in model_json["metabolites"]:
        m["id"] += "0"
        m["compartment"] += "0"
        metabolites[m["id"]] = m
    index = 0
    for r in model_json["reactions"]:
        r["metabolites"] = {i + "0": v for i, v in r["metabolites"].items()}
        comps = {metabolites[k]["compartment"] for k in r["metabolites"]}
        if r["id"].startswith("EX_"):
            r["id"] += "0"
        elif r["id"] in ("BIOMASS_Ecoli_core_w_GAM", "ATPM"):
            r["id"] += "_" + (list(comps)[0] if len(comps) == 1 else "c0")
        else:
            r["id"] = "rxn{:05d}_c0".format(index)
            index += 1
    return cobra.io.from_json(json.dumps(model_json))


@pytest.fixture(scope="module")
def base_model():
    return _build()


@pytest.fixture
def model(base_model):
    return base_model.copy()


def _feature_probs(model, drop_every=0):
    """A deterministic probability per gene (0, 0.2, ..., 1.0 by sorted position); with
    drop_every k, every k-th gene is left out of the map."""
    genes = sorted(g.id for g in model.genes)
    vals = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.5)
    return {g: vals[i % len(vals)] for i, g in enumerate(genes)
            if not (drop_every and i % drop_every == 0)}


def _production_rule(model, feature_probabilities):
    """codiffusion_bioreactor/load_genomes_to_modelseedpy.py l.331-355, verbatim except
    the counter and the comments."""
    for rxn in model.reactions:
        if rxn.genes:
            gene_probs = []
            for gene in rxn.genes:
                gene_prob = feature_probabilities.get(gene.id, 0)
                gene_probs.append(gene_prob)
            probability = max(gene_probs)
        else:
            probability = 0.0
        rxn.probability = probability
        if not hasattr(rxn, 'notes') or rxn.notes is None:
            rxn.notes = {}
        rxn.notes['probability'] = str(probability)


def _coefficients(model):
    """{variable name: objective coefficient} of the model's (linear) objective."""
    coefs = model.solver.objective.get_linear_coefficients(model.variables)
    return {v.name: c for v, c in coefs.items() if c != 0}


# ---------------------------------------------------------------- the reaction rule
def test_reaction_probability_is_max_over_genes_and_zero_without_genes(model):
    fp = _feature_probs(model)
    out = PO.set_reaction_probabilities(model, fp)
    assert set(out) == {r.id for r in model.reactions}
    n_gene, n_free = 0, 0
    for rxn in model.reactions:
        if rxn.genes:
            expect = max(fp[g.id] for g in rxn.genes)
            n_gene += 1
        else:
            expect = 0.0
            n_free += 1
        assert out[rxn.id] == expect
        assert rxn.probability == expect
        assert rxn.notes["probability"] == str(expect)
        assert float(rxn.notes["probability"]) == expect
    # e_coli_core has both kinds; exchanges, biomass and ATPM are gene-less
    assert n_gene > 50 and n_free > 20
    for rid in (BIOMASS, "ATPM_c0", "EX_glc__D_e0"):
        assert model.reactions.get_by_id(rid).notes["probability"] == "0.0"
    # a multi-gene reaction takes the maximum, not a product or a mean
    multi = [r for r in model.reactions if len(r.genes) >= 2
             and len({fp[g.id] for g in r.genes}) >= 2]
    assert multi
    for r in multi:
        assert r.probability == max(fp[g.id] for g in r.genes)


def test_gene_missing_from_the_map_takes_default(model):
    fp = _feature_probs(model, drop_every=3)
    missing = {g.id for g in model.genes} - set(fp)
    assert missing
    out = PO.set_reaction_probabilities(model, fp)
    out_half = PO.set_reaction_probabilities(model.copy(), fp, default=0.5)
    for rxn in model.reactions:
        if rxn.genes:
            assert out[rxn.id] == max(fp.get(g.id, 0.0) for g in rxn.genes)
            assert out_half[rxn.id] == max(fp.get(g.id, 0.5) for g in rxn.genes)
        else:
            assert out[rxn.id] == out_half[rxn.id] == 0.0      # default never applies here


def test_matches_the_production_loader(base_model):
    """The same probability on every reaction as the production loop, and the same note
    string wherever the string is defined. Where a reaction's maximum is a tie between a
    gene missing from the map (production: the int 0) and a gene at 0.0, production's note
    is "0" or "0.0" depending on the iteration order of rxn.genes (a frozenset; it differs
    between model copies), so only the number is compared there."""
    fp = _feature_probs(base_model, drop_every=4)
    ref, mine0, mine = base_model.copy(), base_model.copy(), base_model.copy()
    _production_rule(ref, fp)
    PO.set_reaction_probabilities(mine0, fp, default=0)
    PO.set_reaction_probabilities(mine, fp)
    n_strings = 0
    for r in ref.reactions:
        a, b, c = (m.reactions.get_by_id(r.id) for m in (ref, mine0, mine))
        for x in (b, c):
            assert x.probability == a.probability
            assert float(x.notes["probability"]) == float(a.notes["probability"])
        if all(g.id in fp for g in r.genes):          # no int-0 / 0.0 tie possible
            assert b.notes["probability"] == a.notes["probability"]
            assert c.notes["probability"] == a.notes["probability"]
            n_strings += 1
    assert n_strings > 50


def test_feature_probabilities_from_genome(tmp_path):
    genome = {"id": "g", "features": [
        {"id": "f1", "probability": 0.5}, {"id": "f2"}, {"id": "", "probability": 1.0},
        {"probability": 0.3}, {"id": "f3", "probability": 1.0}],
        "cdss": [{"id": "f1.CDS", "probability": 0.5}]}
    expect = {"f1": 0.5, "f2": 0, "f3": 1.0}
    assert PO.feature_probabilities_from_genome(genome) == expect
    p = tmp_path / "g.json"
    p.write_text(json.dumps(genome))
    assert PO.feature_probabilities_from_genome(str(p)) == expect
    assert PO.feature_probabilities_from_genome(p) == expect


# ---------------------------------------------------------------- keyword arguments
def _fixed_signature_fake(calls):
    def fake(model_s_, environment=None, abundances=None, min_prob=0.01, prob_exp=1,
             ex_weight=100, commkinetics=None, kinetics_coef=1000, printLP=False,
             expression=None, missing_probability=None, biomass_id=None, min_growth=1e-6):
        calls.append(dict(model=model_s_, environment=environment, min_prob=min_prob,
                          prob_exp=prob_exp, ex_weight=ex_weight,
                          missing_probability=missing_probability, biomass_id=biomass_id,
                          min_growth=min_growth))
        return "solution"
    return fake


def test_tuned_values():
    assert PO.TUNED_PRFBA_KWARGS == dict(ex_weight=1, min_prob=0.01, prob_exp=1,
                                         missing_probability=0.0, biomass_id="bio1")
    assert PO.RELEASED_PRFBA_KWARGS == dict(ex_weight=100, min_prob=0.01, prob_exp=1)


def test_kwargs_pass_through(monkeypatch):
    calls = []
    monkeypatch.setattr(MSProbability, "prFBA", staticmethod(_fixed_signature_fake(calls)))
    env = {"cpd00027": 100.0}
    sentinel = object()
    assert PO.run_prfba(sentinel, env) == "solution"
    c = calls.pop()
    assert c["model"] is sentinel and c["environment"] is env
    assert {k: c[k] for k in PO.TUNED_PRFBA_KWARGS} == PO.TUNED_PRFBA_KWARGS
    # tuned=False: the library's released terms; the fix-only keywords stay at its defaults
    PO.run_prfba(sentinel, env, tuned=False)
    c = calls.pop()
    assert {k: c[k] for k in PO.RELEASED_PRFBA_KWARGS} == PO.RELEASED_PRFBA_KWARGS
    assert c["missing_probability"] is None and c["biomass_id"] is None
    # overrides replace single keywords of either set
    PO.run_prfba(sentinel, tuned=True, ex_weight=100, biomass_id="bio2", min_growth=None)
    c = calls.pop()
    assert (c["ex_weight"], c["biomass_id"], c["min_growth"], c["environment"]) == \
        (100, "bio2", None, None)
    assert c["missing_probability"] == 0.0 and c["min_prob"] == 0.01
    PO.run_prfba(sentinel, tuned=False, biomass_id="bio1")
    c = calls.pop()
    assert (c["ex_weight"], c["biomass_id"]) == (100, "bio1")
    assert PO.prfba_kwargs(False, ex_weight=100) == PO.RELEASED_PRFBA_KWARGS


def test_unknown_override_is_refused(monkeypatch):
    calls = []
    monkeypatch.setattr(MSProbability, "prFBA", staticmethod(_fixed_signature_fake(calls)))
    with pytest.raises(TypeError, match="exchange_weight"):
        PO.run_prfba(object(), exchange_weight=1)
    assert not calls


def test_released_library_raises_a_clear_error(monkeypatch):
    def released(model_s_, environment=None, abundances=None, min_prob=0.01, prob_exp=1,
                 ex_weight=100, commkinetics=None, kinetics_coef=1000, printLP=False,
                 expression=None):
        raise AssertionError("must not be called")
    monkeypatch.setattr(MSProbability, "prFBA", staticmethod(released))
    with pytest.raises(ImportError) as ei:
        PO.run_prfba(object(), {})
    msg = str(ei.value)
    assert PO.FIX_BRANCH in msg and "PYTHONPATH=" + PO.FIX_WORKTREE in msg
    assert "missing_probability" in msg and "biomass_id" in msg
    assert modelseedpy.__file__ in msg


def test_missing_modelseedpy_raises_a_clear_error(monkeypatch):
    monkeypatch.setitem(sys.modules, "modelseedpy.core.msprobability", None)
    with pytest.raises(ImportError, match=PO.FIX_BRANCH):
        PO.require_fixed_prfba()


@needs_fix
def test_released_kwargs_are_the_library_defaults():
    params = inspect.signature(MSProbability.prFBA).parameters
    for k, v in PO.RELEASED_PRFBA_KWARGS.items():
        assert params[k].default == v, k
    for k in PO.TUNED_PRFBA_KWARGS:
        assert k in params
    assert params["missing_probability"].default is None
    assert params["biomass_id"].default is None


# ---------------------------------------------------------------- real solves
@needs_fix
def test_tuned_solve_prices_the_objective_as_declared(base_model):
    fp = _feature_probs(base_model)
    m = base_model.copy()
    probs = PO.set_reaction_probabilities(m, fp)
    sol = PO.run_prfba(m, tuned=True, biomass_id=BIOMASS)
    assert sol.status == "optimal"
    coef = _coefficients(m)
    for r in m.reactions:
        f, b = r.forward_variable.name, r.reverse_variable.name
        if r.id.startswith("rxn"):
            c = max(0.01, 1 - probs[r.id])
        elif r.id.startswith("EX_"):
            c = 1
        else:
            c = 0
        assert coef.get(f, 0) == c and coef.get(b, 0) == c, r.id
    # the floor is on the biomass reaction, at 95% of its maximum
    floor = m.constraints["min_value"]
    bio = m.reactions.get_by_id(BIOMASS)
    assert {v.name for v in floor.variables} == {bio.forward_variable.name,
                                                bio.reverse_variable.name}
    assert math.isclose(sol.fluxes[BIOMASS], floor.lb, rel_tol=1e-6)


@needs_fix
def test_wrapper_equals_the_direct_library_call(base_model):
    fp = _feature_probs(base_model, drop_every=5)
    m1, m2 = base_model.copy(), base_model.copy()
    PO.set_reaction_probabilities(m1, fp)
    PO.set_reaction_probabilities(m2, fp)
    s1 = PO.run_prfba(m1, tuned=True, biomass_id=BIOMASS)
    kw = dict(PO.TUNED_PRFBA_KWARGS, biomass_id=BIOMASS)
    s2 = MSProbability.prFBA(m2, environment=None, **kw)
    assert s1.objective_value == s2.objective_value
    assert _coefficients(m1) == _coefficients(m2)
    assert (s1.fluxes - s2.fluxes).abs().max() == 0.0


@needs_fix
def test_released_configuration_prices_exchanges_at_100(base_model):
    fp = _feature_probs(base_model)
    m = base_model.copy()
    PO.set_reaction_probabilities(m, fp)
    sol = PO.run_prfba(m, tuned=False)          # e_coli_core's objective is its biomass
    assert sol.status == "optimal"
    coef = _coefficients(m)
    ex = m.reactions.get_by_id("EX_glc__D_e0")
    assert coef[ex.forward_variable.name] == 100 and coef[ex.reverse_variable.name] == 100


@needs_fix
def test_unannotated_rxn_reaction_is_priced_by_missing_probability(base_model):
    """Tuned: a rxn reaction without a note (as a gapfilled reaction added after annotation)
    is priced at p = 0, cost 1. Released: KeyError before the model is changed."""
    fp = _feature_probs(base_model)
    m = base_model.copy()
    PO.set_reaction_probabilities(m, fp)
    rid = next(r.id for r in m.reactions if r.id.startswith("rxn") and r.genes
               and r.probability >= 0.5)
    del m.reactions.get_by_id(rid).notes["probability"]
    m_rel = m.copy()
    PO.run_prfba(m, tuned=True, biomass_id=BIOMASS)
    r = m.reactions.get_by_id(rid)
    assert _coefficients(m)[r.forward_variable.name] == 1.0
    with pytest.raises(KeyError, match="missing_probability"):
        PO.run_prfba(m_rel, tuned=False)
    assert "min_value" not in m_rel.constraints


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v", "-p", "no:cacheprovider"]))
