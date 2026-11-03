"""Independent verification of Stage E (pfunc_estimators) -- shared context.

Written from the design (PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md, "Stage E") and the raw
inputs. Does NOT import pfunc_estimators, calibrate_pfunc or select_references.
"""
import gzip
import os
import re

import numpy as np
import orjson
import pandas as pd

SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference"
DATA = os.path.join(SR, "data")
OUT = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(OUT, "cache")
os.makedirs(CACHE, exist_ok=True)

COV_MIN, SWED_MAX, FLOOR = 0.90, 0.02, 0.865
UNIVERSE_EPS = 0.01


def jl(path):
    with open(path, "rb") as fh:
        return orjson.loads(fh.read())


def norm_ec(e):
    """score_folds' rule as documented: strip an 'EC:' prefix, drop partial ECs ('-')."""
    s = str(e).strip()
    if s[:3].upper() == "EC:":
        s = s[3:].strip()
    if not s or "-" in s:
        return None
    return s


# select_references._norm reimplemented (species key for the species-dedup check)
_PFX = re.compile(r"^(Candidatus\s+|Ca[._]\s*|midas_[a-z]_)\s*", re.I)


def sp_norm(s):
    if not s:
        return ""
    return _PFX.sub("", str(s).strip()).casefold()


class Ctx:
    def __init__(self, load_truth=True):
        pe = jl(os.path.join(DATA, "panel_eligible.json"))["genomes"]
        self.panel = pe
        self.gids = sorted(pe)
        self.gi = {g: i for i, g in enumerate(self.gids)}
        self.NG = len(self.gids)
        fp = jl(os.path.join(DATA, "fold_plan.json"))
        self.fold_plan = fp
        self.ggroup = fp["genome_group"]
        # per-genome taxonomy labels (raw strings, None kept as None)
        self.lab = {}
        for rk in ("genus", "family", "order", "species"):
            self.lab[rk] = [pe[g].get(rk) for g in self.gids]
        self.lab["taxon"] = [str(pe[g].get("taxon_id")) for g in self.gids]
        self.lab["group"] = [self.ggroup.get(g) for g in self.gids]
        # integer codes; a missing label is unique to its genome (never shared)
        self.code = {}
        for nm, vals in self.lab.items():
            v2 = [(v if v not in (None, "") else f"__missing__{nm}:{g}")
                  for v, g in zip(vals, self.gids)]
            u = {x: j for j, x in enumerate(sorted(set(v2)))}
            self.code[nm] = np.array([u[x] for x in v2], np.int64)
        # EC universe
        uni = jl(os.path.join(DATA, "panel_ec_universe.json"))["ec"]
        ul, seen = [], set()
        for e in uni:
            n = norm_ec(e)
            if n and n not in seen:
                seen.add(n)
                ul.append(n)
        self.ecs = ul
        self.ei = {e: j for j, e in enumerate(ul)}
        self.NU = len(ul)
        if load_truth:
            pec = jl(os.path.join(DATA, "_cache_panel_ec.json"))
            E = np.zeros((self.NG, self.NU), np.uint8)
            self.missing_truth = 0
            for g in self.gids:
                rec = pec.get(g)
                if rec is None:
                    self.missing_truth += 1
                    continue
                idx = [self.ei[n] for n in (norm_ec(x) for x in rec["ec"]) if n in self.ei]
                E[self.gi[g], idx] = 1
            self.E = E


def load_cells(density="1.0"):
    cols = ["amplicon", "md5", "region", "src_genome_id", "fold", "group_id", "density",
            "arm", "space", "outcome", "n_truth", "f1", "n_reps"]
    parts = []
    with gzip.open(os.path.join(DATA, "cells.csv.gz"), "rt") as fh:
        for ch in pd.read_csv(fh, usecols=cols, dtype=str, keep_default_na=False,
                              chunksize=400_000):
            m = ((ch.arm == "prfba_exact_tie") & (ch.space == "ec") &
                 (ch.outcome == "scored") & (ch.density == density))
            parts.append(ch[m])
    df = pd.concat(parts, ignore_index=True)
    df["fold"] = df.fold.astype(int)
    df = df.sort_values(["fold", "md5", "amplicon"]).reset_index(drop=True)
    return df
