#!/usr/bin/env python
"""
prfba_objective.py — prFBA's recommended call of ModelSEEDpy's probability-weighted FBA
(MSProbability.prFBA), and the production rule that puts genome feature probabilities onto
model reactions.

WHAT prFBA SOLVES
  Hold biomass at 95% of its maximum (carbon uptake capped at 100), then minimize

      sum_r max(min_prob, 1 - p_r ** prob_exp) * |v_r|  +  sum_ex ex_weight * |v_ex|

  over the "rxn"-prefixed reactions r, with p_r = float(rxn.notes["probability"]), and the
  "EX_" exchanges. |v| is the forward plus the reverse variable.

THE TUNED CALL (TUNED_PRFBA_KWARGS) AND ITS EVIDENCE
  ex_weight=1, min_prob=0.01, prob_exp=1 (linear cost), missing_probability=0.0 (a reaction
  with no probability note, e.g. a gapfilled one, is priced like a reaction without genetic
  evidence), biomass_id="bio1" (the floor is on biomass, not on whatever objective a
  gapfilled model carries; models built through ATP correction carry ATPM_c0).

  It is the configuration Stage O selected, `est/lin/G0/X1`
  (region_validation/shared_reference/OBJECTIVE_TUNING_RESULTS.md; design
  PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md, section "Stage O"; harness
  scripts/objective_tuning.py). Selected on the tune half, confirmed on the held-out half:
    * Multi-genome selections (n >= 2, pooling A): truth-absent gene-backed flux falls by
      -0.0172 [-0.0221, -0.0118] of D against the shipped objective (k/n, ex_weight=100),
      past the declared SESOI of 0.01, and robust to leaving out any one holdout group.
      Exploratory split: the estimator (E3) alone at ex_weight 100 gives
      -0.0128 [-0.0186, -0.0088]; ex_weight 1 alone with k/n gives -0.0040.
    * The gain is ESTABLISHED ONLY AT n >= 2. At n = 1 it is not established:
      -0.0048 [-0.0113, +0.0063] on the confirm half, after -0.0124 on the tune half.
      The design therefore scopes the genome builder's E3 to multi-genome selections
      (n >= 2), with k/n (= 1.0) kept at n = 1. The production validation gate passed
      (10k panel, exclude_0.987, n >= 2: Brier -0.0530, AUROC +0.0538 against k/n), so E3 is
      the builder's default (--pfunc e3; --pfunc kn reproduces the earlier builds). This
      module does not decide that: it prices whatever feature probabilities the genome
      carries.
  The tuned LP differs from the released one only in ex_weight (100 -> 1); the estimator
  enters through the probabilities written onto the reactions.

  CAVEAT ON ex_weight=1 (Stage O's own). At ex_weight 100 the LP minimizes uptake first; at
  1 an import costs the same as a gene-less reaction, so the LP trades uptake against
  synthesis. Total exchange flux rises and never falls: +2.50 [+1.20, +3.39] units
  (pooling A, confirm half; OBJECTIVE_TUNING_RESULTS.md guard table) against a shipped
  median of about 165 units, up to +18% on single rows. WHICH nutrients are taken up more
  was not recorded (no per-exchange fluxes were stored). The results document recommends
  testing the lower exchange weight on measured exchange (uptake/secretion) data before
  relying on it. Pass tuned=False, or ex_weight=100, to keep the shipped exchange term.

  THE LIBRARY DEFAULT IS DELIBERATELY NOT CHANGED. ModelSEEDpy's MSProbability.prFBA keeps
  ex_weight=100 (RELEASED_PRFBA_KWARGS): other callers (community studies) use the library,
  and ex_weight=1 alters predicted uptake. prFBA's recommendation lives here, at prFBA's
  call site.

LIBRARY REQUIREMENT
  The released ModelSEEDpy (atp-safe-and-drain-fixes, c2e31a0, the tree the venv's editable
  install points at) cannot run prFBA on one model (TypeError on len(model), then
  PackageError from a second ElementUptakePkg) and has no missing_probability / biomass_id
  keywords. The fix is branch `fix-prfba-single-model` (pushed to origin, 09856f5), checked
  out as the worktree ~/Documents/ModelSEEDpy-prfba. Use it with
      PYTHONPATH=/home/freiburger/Documents/ModelSEEDpy-prfba
  run_prfba() raises ImportError on a library without those keywords.

EQUIVALENCE
  region_validation/shared_reference/scripts/objective_wrapper_equivalence.py checks that
  run_prfba() reproduces Stage O's harness LP (objective_tuning.prfba_tuned_copy) at
  est/lin/G0/X1, and the released configuration at kn/lin/G0/X100, on gapfilled Stage O
  models (output: data/objective_tuning/wrapper_equivalence.json).

USE
    from prfba_objective import (feature_probabilities_from_genome,
                                 set_reaction_probabilities, run_prfba)
    fp = feature_probabilities_from_genome("genome_objects/ASV.json")
    set_reaction_probabilities(model, fp)          # after gapfilling: gapfilled p = 0.0
    sol = run_prfba(model.copy(), environment={cpd: 100.0 for cpd in medium})

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import inspect
import json
import os

__all__ = ["TUNED_PRFBA_KWARGS", "RELEASED_PRFBA_KWARGS", "FIX_BRANCH", "FIX_WORKTREE",
           "REQUIRED_KEYWORDS", "feature_probabilities_from_genome",
           "set_reaction_probabilities", "require_fixed_prfba", "prfba_kwargs", "run_prfba"]

# Stage O's selected configuration est/lin/G0/X1 (ex_weight an int, as the harness's X_VALUE)
TUNED_PRFBA_KWARGS = dict(ex_weight=1, min_prob=0.01, prob_exp=1, missing_probability=0.0,
                          biomass_id="bio1")

# The library's own signature defaults for the objective terms (ModelSEEDpy c2e31a0 and the
# fix branch alike). missing_probability and biomass_id are left at the library's defaults
# (None: a reaction without a probability note raises KeyError; the floor applies to the
# objective the model carries on entry, so set model.objective = "bio1" first).
RELEASED_PRFBA_KWARGS = dict(ex_weight=100, min_prob=0.01, prob_exp=1)

FIX_BRANCH = "fix-prfba-single-model"
FIX_WORKTREE = "/home/freiburger/Documents/ModelSEEDpy-prfba"
REQUIRED_KEYWORDS = ("missing_probability", "biomass_id")


# ---------------------------------------------------------------- probabilities
def feature_probabilities_from_genome(genome_json) -> dict:
    """{feature id: probability} from a KBase-style genome (prFBA's synthetic genome).

    `genome_json` is the genome dict or a path to its JSON (`genome_objects/{asv}.json`, the
    copy the model builders read). The rule is the production loader's
    (codiffusion_bioreactor/load_genomes_to_modelseedpy.py l.276-282): every feature with a
    non-empty id, its `probability` or 0 when absent. Values are kept as stored."""
    if isinstance(genome_json, (str, os.PathLike)):
        with open(genome_json) as fh:
            genome_json = json.load(fh)
    feature_probabilities = {}
    for feature in genome_json.get("features", []):
        feature_id = feature.get("id", "")
        probability = feature.get("probability", 0)
        if feature_id:
            feature_probabilities[feature_id] = probability
    return feature_probabilities


def set_reaction_probabilities(model, feature_probabilities, default=0.0) -> dict:
    """Write each reaction's probability: the maximum over its genes' feature probabilities,
    0.0 for a reaction without genes. Returns {reaction id: probability}.

    The production rule (codiffusion_bioreactor/load_genomes_to_modelseedpy.py l.331-355;
    the same in build_prfba_community_models.py and EmilyKin's ek_pipeline.py): for EVERY
    reaction of the model, `rxn.probability = p` and `rxn.notes["probability"] = str(p)`.
    Gene-less reactions (exchanges, biomass, transport without a GPR, and gapfilled
    reactions when this runs after gapfilling) get 0.0, so prFBA prices them at 1.

    `default` is the probability of a gene that has no entry in `feature_probabilities`.
    The production loader uses the int 0 there; 0.0 is the same number, so every
    probability and every LP coefficient is identical, but a reaction whose maximum is such
    a missing gene has the note "0.0" here where production writes "0". Production's own
    string is not defined where a missing gene (int 0) ties a gene at 0.0: max() keeps the
    first, and rxn.genes is a frozenset whose order differs between model copies. Pass
    default=0 to reproduce production's strings wherever they are defined."""
    out = {}
    for rxn in model.reactions:
        if rxn.genes:
            probability = max(feature_probabilities.get(gene.id, default) for gene in rxn.genes)
        else:
            probability = 0.0
        rxn.probability = probability
        if not hasattr(rxn, "notes") or rxn.notes is None:
            rxn.notes = {}
        rxn.notes["probability"] = str(probability)
        out[rxn.id] = probability
    return out


# ---------------------------------------------------------------- the prFBA call
def require_fixed_prfba():
    """Return MSProbability from the imported modelseedpy, or raise ImportError if its
    prFBA lacks the fix branch's keywords (missing_probability, biomass_id)."""
    try:
        import modelseedpy
        from modelseedpy.core.msprobability import MSProbability
    except ImportError as exc:
        raise ImportError(
            "prfba_objective needs ModelSEEDpy's MSProbability ({}). Put the {} worktree "
            "first on the path: PYTHONPATH={}".format(exc, FIX_BRANCH, FIX_WORKTREE)) from exc
    params = inspect.signature(MSProbability.prFBA).parameters
    missing = [k for k in REQUIRED_KEYWORDS if k not in params]
    if missing:
        raise ImportError(
            "The imported ModelSEEDpy ({}) has an MSProbability.prFBA without {}: it is the "
            "released code, whose single-model path does not run. Use ModelSEEDpy branch "
            "'{}' (worktree {}): run with PYTHONPATH={}".format(
                getattr(modelseedpy, "__file__", "?"), ", ".join(missing), FIX_BRANCH,
                FIX_WORKTREE, FIX_WORKTREE))
    return MSProbability


def prfba_kwargs(tuned=True, **overrides) -> dict:
    """The keyword arguments run_prfba() passes to MSProbability.prFBA (beside the model and
    the environment): TUNED_PRFBA_KWARGS or RELEASED_PRFBA_KWARGS, then the overrides."""
    kwargs = dict(TUNED_PRFBA_KWARGS if tuned else RELEASED_PRFBA_KWARGS)
    kwargs.update(overrides)
    return kwargs


def run_prfba(model, environment=None, tuned=True, **overrides):
    """MSProbability.prFBA(model, environment, **kwargs) with prFBA's tuned objective
    (TUNED_PRFBA_KWARGS), or with the library's released defaults when tuned=False
    (RELEASED_PRFBA_KWARGS: ex_weight 100; the caller then sets model.objective to the
    biomass reaction and annotates every rxn reaction). `overrides` replace single keywords,
    e.g. ex_weight=100 or biomass_id="bio1".

    The model's reactions must carry notes["probability"] (set_reaction_probabilities).
    The model is MODIFIED IN PLACE (medium, carbon-uptake package, the 95% floor, the
    objective); pass model.copy() to keep the original. A second call on the same object
    raises PackageError. `environment` is a {compound or exchange id: uptake} dict or an
    MSMedia; None keeps the model's own exchange bounds. Returns the cobra Solution."""
    MSProbability = require_fixed_prfba()
    kwargs = prfba_kwargs(tuned, **overrides)
    params = inspect.signature(MSProbability.prFBA).parameters
    unknown = sorted(k for k in kwargs if k not in params)
    if unknown:
        raise TypeError("MSProbability.prFBA has no keyword(s) {}".format(unknown))
    return MSProbability.prFBA(model, environment=environment, **kwargs)
