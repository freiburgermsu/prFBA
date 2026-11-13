#!/usr/bin/env python
"""census.py -- part 1: every consumer of feature / reaction probability in the two deployments.

  sweep   grep codiffusion_bioreactor and EmilyKin (code: *.py, *.sh; notebooks: code-cell
          sources of *.ipynb) for probability / reaction_scores / MSProbability / prFBA, outside
          the data directories. Every hit is listed (results/census_sweep.json).
  table   the curated consumer table (CONSUMERS below). Each entry names file:line anchors with
          a regex that the line must match TODAY; the script re-reads every anchor and fails
          the entry if the line moved or changed. Entries also carry the evidence of which
          build the consumer reads (a file:line anchor of its genome / model path).
  scan    (--scan) the probability range of both E3 builds over every reweighted (n >= 2)
          genome: global min / max of the E3 values and of k/n, number of zeros, and how many
          replaced features cross p = 0.5 (the only fixed threshold a consumer applies besides
          p > 0). Reads SYN_e3/genome_objects/{asv}.json and SYN_e3/provenance/{asv}.json
          (`probability_kn`, the replaced k/n). Read-only; 8 processes.

Output: results/census.json (+ census_sweep.json).
  PYTHONDONTWRITEBYTECODE=1 ~/Documents/py_venv/bin/python census.py [--scan]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

PATTERN = re.compile(r"probability|reaction_scores|MSProbability|prFBA\(")
REPOS = {"codiffusion_bioreactor": C.CODIF, "EmilyKin": C.EK}
# data / generated directories (no code), plus vendored copies
SKIP_DIRS = re.compile(
    r"/(synthetic_genomes[^/]*|comm_models[^/]*|models[^/]*|media_gapfilled_models[^/]*|"
    r"gapfilled_models|genome_objects|genomes|runs|_pinned|\.ipynb_checkpoints|__pycache__|"
    r"__marimo__|bvbrc_alignment_hits[^/]*|bvbrc_embedding_hits|kbase_r232|\.git)(/|$)")

# ------------------------------------------------------------------ the curated consumer table
# anchors: (relative file, line, regex the line must match)
CODIF_BUILD = ("codiffusion_bioreactor/build_prfba_community_models.py", 41,
               r'^GDIR = "synthetic_genomes_54d1150"')
EK_CFG = "EmilyKin/community_modeling/config.json"
CONSUMERS = [
    # ------------------------------------------------------------- codiffusion_bioreactor
    {"id": "codif.build.notes",
     "anchors": [("codiffusion_bioreactor/build_prfba_community_models.py", 356,
                  r'probs = \{f\.get\("id"\): f\.get\("probability", 0\)'),
                 ("codiffusion_bioreactor/build_prfba_community_models.py", 361,
                  r"p = max\(\(probs\.get\(gg\.id, 0\) for gg in r\.genes\)"),
                 ("codiffusion_bioreactor/build_prfba_community_models.py", 362,
                  r'r\.notes\["probability"\] = str\(p\)')],
     "role": "WRITER: draft member model (stage build, _build_one): reaction p = max over GPR genes "
             "of feature p, 0.0 if gene-less, into reaction.notes['probability'] (SBML notes).",
     "path": "production (rebuild_pipeline.sh step 1, TAG v7e)",
     "reads_build": "synthetic_genomes_54d1150 (20-hit, k/n)", "build_anchor": CODIF_BUILD,
     "e3_changes_output": "only the notes values (no reaction, bound or GPR depends on p). The notes "
                          "are carried into mediagf and community SBML but no production step "
                          "reads them (rows below)."},
    {"id": "codif.mediagf",
     "anchors": [("codiffusion_bioreactor/build_prfba_community_models.py", 541,
                  r"gf = MSGapfill\(model, default_gapfill_templates=\[template\], test_conditions=tests,"),
                 ("codiffusion_bioreactor/build_prfba_community_models.py", 542, r"^\s+blacklist=blacklist\)"),
                 ("codiffusion_bioreactor/build_prfba_community_models.py", 558,
                  r"gf = MSGapfill\(model, default_gapfill_templates=\[_load_template"),
                 ("codiffusion_bioreactor/build_prfba_community_models.py", 559, r"^\s+test_conditions=tests\)")],
     "role": "member gapfill per (member, day medium): MSGapfill WITHOUT reaction_scores, so "
             "GapfillingPkg.compute_gapfilling_penalties multiplies every candidate by 2 - 0 = 2.",
     "path": "production (mediagf; the communities are assembled from its models)",
     "reads_build": "models_<TAG>/ drafts of synthetic_genomes_54d1150", "build_anchor": CODIF_BUILD,
     "e3_changes_output": "no: p is never passed to the gapfiller."},
    {"id": "codif.gapfill_models",
     "anchors": [("codiffusion_bioreactor/gapfill_models.py", 600, r"reaction_scores=\{\},")],
     "role": "legacy member gapfill (AuxoMedia + PyruvateMinimalAnaerobic): reaction_scores={} "
             "explicitly, factor 2 for every candidate.",
     "path": "legacy (stage gapfill; not part of --stage all, superseded by mediagf)",
     "reads_build": "models/ drafts (bundled-template run)", "build_anchor": CODIF_BUILD,
     "e3_changes_output": "no."},
    {"id": "codif.compartment_aware_gapfill",
     "anchors": [("codiffusion_bioreactor/compartment_aware_gapfill.py", 76,
                  r'p = float\(rxn\.notes\.get\("probability", 1\.0\)\)'),
                 ("codiffusion_bioreactor/compartment_aware_gapfill.py", 72, r"if not rxn\.genes:"),
                 ("codiffusion_bioreactor/compartment_aware_gapfill.py", 113,
                  r'score = table\[gene\]\["probability"\] if isinstance'),
                 ("codiffusion_bioreactor/compartment_aware_gapfill.py", 120, r"factor = 2 - highest_score")],
     "role": "per-member gapfill penalties: for each gene-bearing reaction of member i's draft, "
             "{gene: p} with p = notes['probability'] (a MISSING note defaults to 1.0, not 0); the "
             "candidate rxn_ci is priced base x (2 - max p).",
     "path": "PROTOTYPE (its own docstring); only caller demo_compartment_scoring.py",
     "reads_build": "models_<tag>/ drafts (default tag v7) of synthetic_genomes_54d1150",
     "build_anchor": ("codiffusion_bioreactor/compartment_aware_gapfill.py", 64,
                      r'models_dir = os\.path\.join\(repo, f"models_\{tag\}" if tag else "models"\)'),
     "e3_changes_output": "yes, for the reversal candidates of member-draft reactions whose p changes "
                          "(measured in part 2); the demo itself never solves a community gapfill."},
    {"id": "codif.demo_compartment_scoring",
     "anchors": [("codiffusion_bioreactor/demo_compartment_scoring.py", 109,
                  r'pkg\.parameters\["member_reaction_scores"\] = member_scores'),
                 ("codiffusion_bioreactor/demo_compartment_scoring.py", 140,
                  r"scores = CAG\.member_scores_from_drafts\(prov, tag=a\.tag\)")],
     "role": "demo: stock vs compartment-aware penalty vectors of a day-174 4-member subcommunity.",
     "path": "demo (output modeling_files/compartment_scoring_demo_day174_4members.json, 2026-09-24)",
     "reads_build": "models_v7/ drafts (54d1150)", "build_anchor": CODIF_BUILD,
     "e3_changes_output": "yes, penalties only (part 2)."},
    {"id": "codif.loader.build_model",
     "anchors": [("codiffusion_bioreactor/load_genomes_to_modelseedpy.py", 280,
                  r"probability = feature\.get\('probability', 0\)"),
                 ("codiffusion_bioreactor/load_genomes_to_modelseedpy.py", 355,
                  r"rxn\.notes\['probability'\] = str\(probability\)")],
     "role": "WRITER: GenomeLoader.build_model (bundled template, MSBuilder.build) sets "
             "rxn.probability / notes. The production builders use GenomeLoader only for "
             "load_genome / get_msgenome, not build_model.",
     "path": "old modeling path (modeling.py marimo notebook 2026-05, integrate_with_existing_code.py, tests/)",
     "reads_build": "GenomeLoader() default genome_objects/ (pre-prFBA genomes)",
     "build_anchor": ("codiffusion_bioreactor/modeling.py", 75, r"loader = GenomeLoader\(\)"),
     "e3_changes_output": "notes values only; nothing downstream reads them."},
    {"id": "codif.wl_audit",
     "anchors": [("codiffusion_bioreactor/methane_nongrowth/atp_wl/wl_audit/audit_merge.py", 187,
                  r'"probability": f\.get\("probability"\)'),
                 ("codiffusion_bioreactor/methane_nongrowth/atp_wl/wl_audit/audit_merge.py", 194,
                  r'strong = \[h for h in \(hits or \[\]\) if \(h\.get\("probability"\) or 0\) >= 0\.5\]'),
                 ("codiffusion_bioreactor/methane_nongrowth/atp_wl/wl_audit/audit_merge.py", 262,
                  r'c = "GAP" if any\(\(h\.get\("probability"\) or 0\) >= 0\.5'),
                 ("codiffusion_bioreactor/methane_nongrowth/atp_wl/wl_audit/audit.py", 254,
                  r'"probability": ft\.get\("probability"\)')],
     "role": "audit table of Wood-Ljungdahl / methanogenesis steps: an absent step is 'GAP' when a "
             "gene matching its roles has feature p >= 0.5, else '--(weak gene)'. The ONLY fixed "
             "p threshold other than p > 0 in either project's code.",
     "path": "analysis (one-off audit, methane_nongrowth/atp_wl, committed 2026-11)",
     "reads_build": "synthetic_genomes_54d1150",
     "build_anchor": ("codiffusion_bioreactor/methane_nongrowth/atp_wl/wl_audit/audit.py", 45,
                      r'GENOME_DIR = os\.path\.join\(ROOT, "synthetic_genomes_54d1150"\)'),
     "e3_changes_output": "only if re-pointed at an E3 build: a step can flip GAP <-> --(weak gene) "
                          "where E3 moves a matching gene across 0.5 (scan: crossing counts)."},
    {"id": "codif.util.MSProbability",
     "anchors": [("codiffusion_bioreactor/util.py", 19,
                  r"from modelseedpy\.core\.msprobability import MSProbability")],
     "role": "dead import: util.py never calls MSProbability; no prFBA( / MSProbability.prFBA call "
             "exists in either project (sweep).",
     "path": "none", "reads_build": "-", "build_anchor": None, "e3_changes_output": "no."},
    {"id": "codif.kbbasemodules",
     "anchors": [("codiffusion_bioreactor/kbbasemodules/basemodelingmodule.py", 263,
                  r'if "probability" in item\["scores"\]:')],
     "role": "KBase SDK base module: takes reaction probability from annotation-ontology event "
             "scores (a KBase workspace object), not from genome feature 'probability'.",
     "path": "vendored library, not imported by any pipeline script (sweep)",
     "reads_build": "-", "build_anchor": None, "e3_changes_output": "no."},
    {"id": "codif.bvbrc_converter",
     "anchors": [("codiffusion_bioreactor/bvbrc_to_kbase_genome.py", 1188,
                  r"functions\[function\]\['probability'\] /= num_genomes")],
     "role": "PRODUCER of k/n (the pre-prFBA BV-BRC -> KBase converter), not a consumer.",
     "path": "producer (old genome build)", "reads_build": "-", "build_anchor": None,
     "e3_changes_output": "n/a."},
    # ------------------------------------------------------------- EmilyKin
    {"id": "ek.genome_evidence.ppk",
     "anchors": [("EmilyKin/community_modeling/ek_pipeline.py", 769,
                  r'p = float\(ft\.get\("probability", 0\) or 0\)'),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 777, r"bucket\[fn\] = max\(p, bucket\[fn\]\)"),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 779,
                  r'"ppk_evidence": bool\(ppk1 or ppk2\)'),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 780, r'"ppk_max_probability": max\(')],
     "role": "stage templates: polyP evidence. ppk_evidence (the decision: polyP pool + sink added "
             "after gapfilling) is PRESENCE of a matching function (the defaultdict key exists for "
             "any p, even 0); ppk_max_probability is only recorded in template_assignment.json "
             "(no config or code thresholds it).",
     "path": "production (ek_pipeline stage templates)",
     "reads_build": "PATHS['genome_objects'] = EmilyKin/synthetic_genomes/genome_objects (20-hit)",
     "build_anchor": ("EmilyKin/community_modeling/ek_pipeline.py", 763,
                      r'g = json\.load\(open\(os\.path\.join\(PATHS\["genome_objects"\]'),
     "e3_changes_output": "recorded value only; the decision is p-independent."},
    {"id": "ek.build.notes",
     "anchors": [("EmilyKin/community_modeling/ek_pipeline.py", 879,
                  r'probs = \{f\.get\("id"\): f\.get\("probability", 0\)'),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 894,
                  r"p = max\(\(probs\.get\(gg\.id, 0\) for gg in r\.genes\)"),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 895, r'r\.notes\["probability"\] = str\(p\)')],
     "role": "WRITER: draft member models, the codiffusion rule (max over GPR genes, 0 gene-less).",
     "path": "production (ek_pipeline stage build)",
     "reads_build": "PATHS['genome_objects'] (config.json and config_r13_syn.json: synthetic_genomes/"
                    "genome_objects, 20-hit; config_r13.json: runs/r13/genome_objects_overlay, the "
                    "planned MAG-union overlay, not built)",
     "build_anchor": ("EmilyKin/community_modeling/ek_pipeline.py", 859,
                      r'_W\["loader"\] = GenomeLoader\(genome_dir=PATHS\["genome_objects"\]\)'),
     "e3_changes_output": "notes values only."},
    {"id": "ek.atpcorr.record",
     "anchors": [("EmilyKin/community_modeling/ek_pipeline.py", 2515,
                  r'prob = \{r\.id: float\(r\.notes\.get\("probability", 0\) or 0\)'),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 2555, r'"genome_probability": prob\.get\(r\.id, 0\.0\)'),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 2603,
                  r'"filtered_with_genome_evidence": \[f\["id"\] for f in filtered if f\["genome_probability"\] > 0\]')],
     "role": "stage atpcorr: records each ATP-filtered reaction's p and lists the filtered reactions "
             "with p > 0 (a REPORT field of <rid>.atp.json; nothing reads it back to decide).",
     "path": "production (record only)", "reads_build": "models/ drafts (20-hit)",
     "build_anchor": ("EmilyKin/community_modeling/ek_pipeline.py", 2511,
                      r'model = read_sbml_model\(P\("models", f"\{rid\}\.xml"\)\)'),
     "e3_changes_output": "recorded values; the p > 0 list is NOT strictly invariant: k/n is > 0 for every "
                          "feature, but E3 is exactly 0 for 67 features of the two builds (scan: "
                          "zero_features), so a filtered reaction whose genes are all among them would "
                          "drop out of the list. No decision reads the list."},
    {"id": "ek.assimilation.record",
     "anchors": [("EmilyKin/community_modeling/ek_pipeline.py", 2397,
                  r'"genome_probability": float\(r\.notes\.get\("probability", 0\) or 0\)')],
     "role": "assimilation_exemption: records p of each flux-carrying reaction (report).",
     "path": "production (record only)", "reads_build": "atp-corrected member models (20-hit)",
     "build_anchor": None, "e3_changes_output": "recorded values only."},
    {"id": "ek.gate.record",
     "anchors": [("EmilyKin/community_modeling/ek_pipeline.py", 4496,
                  r'p = float\(mr\.notes\.get\("probability", 0\) or 0\)'),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 5049,
                  r'key=lambda c: \(not c\.get\("gene_bearing"\), -rank\.get\(c\["tier"\], 0\), abs\(c\["flux"\]\), c\["comm_rxn"\]\)\)'),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 5857, r'"gpr", "genome_probability", "db_reversibility"')],
     "role": "community energy gate: every candidate/applied change records genome_probability; the "
             "candidate order (l.5049) uses gene_bearing, tier and flux, NOT p.",
     "path": "production (record only)", "reads_build": "community models (20-hit members)",
     "build_anchor": None, "e3_changes_output": "recorded values only."},
    {"id": "ek.mediagf",
     "anchors": [("EmilyKin/community_modeling/ek_pipeline.py", 2792,
                  r"gf = MSGapfill\(model0, default_gapfill_templates=\[template\], test_conditions=tests, blacklist=blacklist,"),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 3211,
                  r"gf = MSGapfill\(model, default_gapfill_templates=\[_load_template\(template_name\)\], test_conditions=tests,"),
                 ("EmilyKin/community_modeling/ek_pipeline.py", 3232,
                  r"gf2 = MSGapfill\(model, default_gapfill_templates=\[_load_template\(template_name\)\], test_conditions=tests,"),
                 ("EmilyKin/community_modeling/curate_templates.py", 433,
                  r'"reaction_scores": "MSGapfill\(reaction_scores=\.\.\.\) / GapfillingPkg parameter: not passed by the pipeline \(empty\)"')],
     "role": "member gapfill: no reaction_scores (the pipeline's own record says so, curate_templates "
             "l.433), factor 2 for every candidate.",
     "path": "production", "reads_build": "ATP-corrected drafts (20-hit)", "build_anchor": None,
     "e3_changes_output": "no."},
    {"id": "ek.report_fluxes",
     "anchors": [("EmilyKin/community_modeling/report_fluxes.py", 864, r"fmt\(c\.get\('genome_probability'\)\)")],
     "role": "prints the recorded genome_probability column of gate changes.",
     "path": "report", "reads_build": "gate.json of a run", "build_anchor": None,
     "e3_changes_output": "printed values only."},
    {"id": "ek.analysis.storage_census",
     "anchors": [("EmilyKin/community_modeling/plan_r13/maps/storage_work/gene_census.py", 91,
                  r'p = float\(ft\.get\("probability", 0\) or 0\)'),
                 ("EmilyKin/community_modeling/plan_r13/maps/storage_work/gene_census.py", 157,
                  r'"n_members_p05": sum\(1 for m in members if cls_of\[m\]\[k\] >= 0\.5\)')],
     "role": "r13 planning analysis: per member and storage class the max p of matching genes; "
             "reports members with any (p > 0) and with p >= 0.5.",
     "path": "analysis (plan_r13 maps)", "reads_build": "EmilyKin/synthetic_genomes/genome_objects (20-hit)",
     "build_anchor": ("EmilyKin/community_modeling/plan_r13/maps/storage_work/gene_census.py", 10,
                      r'GO = f"\{EK\}/synthetic_genomes/genome_objects"'),
     "e3_changes_output": "only if re-pointed: the p >= 0.5 counts can change (scan: 0.5 crossings); the "
                          "presence counts (best p > 0, best starts at 0.0) only through the 67 E3-zero "
                          "features."},
    {"id": "ek.analysis.other",
     "anchors": [("EmilyKin/community_modeling/plan_r13/maps/storage_work/vocab_scan.py", 28,
                  r'p = float\(ft\.get\("probability", 0\) or 0\)'),
                 ("EmilyKin/community_modeling/plan_r13/maps/evidence/members/scripts/genome_scan.py", 35,
                  r"pr=float\(ft\.get\('probability',0\) or 0\)"),
                 ("EmilyKin/community_modeling/plan_r13/maps/mags_probe/nif_probe.py", 38,
                  r'fn\[f\] = max\(fn\.get\(f, 0\), ft\.get\("probability", 0\)\)'),
                 ("EmilyKin/community_modeling/analysis/data_inconsistencies/verify/D14_data/d14_verify.py", 59,
                  r'"probability": round\(c\.get\("probability", 0\), 4\)')],
     "role": "planning / verification scripts that tabulate the max or mean p of matching genes "
             "(no thresholds besides presence).",
     "path": "analysis", "reads_build": "EmilyKin/synthetic_genomes (20-hit)",
     "build_anchor": ("EmilyKin/community_modeling/plan_r13/maps/storage_work/vocab_scan.py", 5,
                      r'GO = "/home/freiburger/Documents/EmilyKin/synthetic_genomes/genome_objects"'),
     "e3_changes_output": "tabulated values only, and only if re-pointed."},
]


# sweep hits outside the consumer table: why each file is not a consumer
NOT_CONSUMERS = {
    "codiffusion_bioreactor/BV_BRC.ipynb": "notebook statistics of the pre-prFBA build's probability distribution",
    "codiffusion_bioreactor/BV_BRC-Copy1.ipynb": "notebook statistics / plot of the pre-prFBA build's probabilities",
    "codiffusion_bioreactor/BV_BRC-Copy2.ipynb": "notebook PRODUCER of k/n (pre-prFBA build) + statistics",
    "codiffusion_bioreactor/bvbrc_to_kbase_genome-old.py": "old PRODUCER of k/n",
    "codiffusion_bioreactor/SpiecEasi/spiec_easi_core.py": "unrelated (edge-selection probability of SPIEC-EASI)",
    "EmilyKin/network_method_comparison/spiec_easi_core.py": "unrelated (SPIEC-EASI)",
    "EmilyKin/community_modeling/load_genomes_to_modelseedpy.py": "copy of the codiffusion loader (WRITER via build_model; ek_pipeline uses load_genome / get_msgenome only)",
    "EmilyKin/community_modeling/codiffusion_copies/build_prfba_community_models.py": "archived copy of the codiffusion builder (PROVENANCE.md), not imported",
    "EmilyKin/community_modeling/codiffusion_copies/load_genomes_to_modelseedpy.py": "archived copy, not imported",
    "EmilyKin/community_modeling/codiffusion_copies/util.py": "archived copy (dead MSProbability import), not imported",
    "EmilyKin/community_modeling/plan_evidence/round12/tamper/work/pins_code_edit/code/curate_templates.py": "tamper-test copy of curate_templates.py (documentation string)",
    "EmilyKin/community_modeling/plan_evidence/round12/tamper/work/thermo10_code_edit/code/curate_templates.py": "tamper-test copy of curate_templates.py (documentation string)",
    "EmilyKin/community_modeling/plan_r13/_revision/part2.py": "design text (MAG-union probability = bins / 3)",
    "EmilyKin/community_modeling/plan_r13/_revision/part3.py": "design text",
    "EmilyKin/community_modeling/plan_r13/maps/storage_work/gate_shares.py": "docstring only; reads gene_census output (presence, p > 0)",
    "EmilyKin/community_modeling/analysis/data_inconsistencies/verify/D14_prov/model_gpr.py": "prints p of CDS entries of the 20-hit build (display)",
    "EmilyKin/community_modeling/analysis/data_inconsistencies/dedup/parts/part2.py": "manuscript text strings",
}


def read_line(rel, line):
    p = os.path.join(C.DOC, rel)
    try:
        with open(p, errors="replace") as fh:
            for i, txt in enumerate(fh, 1):
                if i == line:
                    return txt.rstrip("\n")
    except FileNotFoundError:
        return None
    return None


def verify_table():
    out, n_bad = [], 0
    for e in CONSUMERS:
        e = dict(e)
        checks = []
        for rel, line, rx in e["anchors"] + ([e["build_anchor"]] if e.get("build_anchor") else []):
            txt = read_line(rel, line)
            ok = txt is not None and re.search(rx, txt) is not None
            n_bad += not ok
            checks.append({"file": rel, "line": line, "ok": ok, "text": (txt or "")[:200].strip()})
        e["anchors"] = [f"{rel}:{line}" for rel, line, _ in e["anchors"]]
        e["build_anchor"] = (f"{e['build_anchor'][0]}:{e['build_anchor'][1]}" if e.get("build_anchor") else None)
        e["verified"] = all(c["ok"] for c in checks)
        e["anchor_text"] = checks
        out.append(e)
    return out, n_bad


def ek_config_paths():
    rows = {}
    for f in ("config.json", "config_r13.json", "config_r13_syn.json"):
        p = f"{C.EK}/community_modeling/{f}"
        cfg = json.load(open(p))
        go = cfg["paths"]["genome_objects"]
        rows[f] = {"genome_objects": go, "exists": os.path.isdir(go),
                   "is_d100": "d100" in go}
    return rows


def sweep():
    hits = []
    for repo_name, root in REPOS.items():
        for dirpath, dirnames, filenames in os.walk(root):
            rel_dir = dirpath[len(C.DOC):]
            dirnames[:] = [d for d in dirnames if not SKIP_DIRS.search(f"{rel_dir}/{d}")]
            for fn in filenames:
                if fn.startswith("._"):
                    continue
                path = os.path.join(dirpath, fn)
                rel = os.path.relpath(path, C.DOC)
                if fn.endswith((".py", ".sh")):
                    try:
                        with open(path, errors="replace") as fh:
                            for i, line in enumerate(fh, 1):
                                if PATTERN.search(line):
                                    hits.append({"file": rel, "line": i, "kind": "code", "text": line.strip()[:220]})
                    except OSError:
                        continue
                elif fn.endswith(".ipynb"):
                    try:
                        nb = json.load(open(path))
                    except Exception:
                        continue
                    for ci, cell in enumerate(nb.get("cells", [])):
                        if cell.get("cell_type") != "code":
                            continue
                        for li, line in enumerate("".join(cell.get("source", "")).splitlines(), 1):
                            if PATTERN.search(line):
                                hits.append({"file": rel, "line": f"cell{ci}:l{li}", "kind": "notebook",
                                             "text": line.strip()[:220]})
    return hits


def classify_hits(hits, table):
    anchored = defaultdict(set)
    for e in table:
        for a in e["anchors"]:
            f, l = a.rsplit(":", 1)
            anchored[f].add(int(l))
    by_file = Counter(h["file"] for h in hits)
    covered_files = {a.rsplit(":", 1)[0] for e in table for a in e["anchors"]}
    return {"n_hits": len(hits), "n_files": len(by_file),
            "files": dict(sorted(by_file.items())),
            "files_with_a_table_entry": sorted(f for f in by_file if f in covered_files),
            "files_not_consumers": {f: NOT_CONSUMERS[f] for f in sorted(by_file)
                                    if f not in covered_files and f in NOT_CONSUMERS},
            "files_unclassified": sorted(f for f in by_file if f not in covered_files and f not in NOT_CONSUMERS)}


# ------------------------------------------------------------------ probability range scan
def _scan_one(args):
    import orjson
    dep, asv = args
    e3 = C.BUILDS[dep]["e3"]
    with open(f"{e3}/genome_objects/{asv}.json", "rb") as fh:
        g = orjson.loads(fh.read())
    with open(f"{e3}/provenance/{asv}.json", "rb") as fh:
        pk = orjson.loads(fh.read()).get("probability_kn") or {}
    ps = [f.get("probability") for f in g.get("features", [])]
    rep = [(pk[f["id"]], f["probability"]) for f in g.get("features", []) if f["id"] in pk]
    cds = {c["id"]: c.get("probability") for c in g.get("cdss", [])}
    cds_mismatch = sum(1 for f in g.get("features", [])
                       if f"{f['id']}.CDS" in cds and cds[f"{f['id']}.CDS"] != f.get("probability"))
    none = sum(1 for p in ps if p is None)
    ps = [p for p in ps if p is not None]
    up = sum(1 for k, e in rep if k < 0.5 <= e)
    down = sum(1 for k, e in rep if e < 0.5 <= k)
    return {"dep": dep, "asv": asv, "n_features": len(ps) + none, "n_none": none,
            "min_all": min(ps) if ps else None, "max_all": max(ps) if ps else None,
            "n_zero": sum(1 for p in ps if p == 0),
            "n_replaced": len(rep),
            "min_kn_replaced": min((k for k, _ in rep), default=None),
            "min_e3": min((e for _, e in rep), default=None),
            "max_e3": max((e for _, e in rep), default=None),
            "n_e3_zero": sum(1 for _, e in rep if e == 0),
            "n_changed": sum(1 for k, e in rep if k != e),
            "cross05_up": up, "cross05_down": down,
            "n_kn_ge05": sum(1 for k, _ in rep if k >= 0.5), "n_e3_ge05": sum(1 for _, e in rep if e >= 0.5),
            "cds_mismatch": cds_mismatch,
            "zero_features": [{"asv": asv, "feature": f["id"], "functions": f.get("functions"),
                               "probability_kn": pk.get(f["id"])}
                              for f in g.get("features", []) if f.get("probability") == 0]}


def scan(procs=8):
    jobs = []
    for dep, b in C.BUILDS.items():
        man = json.load(open(f"{b['e3']}/synthetic_genome_manifest.json"))
        for asv, v in man.items():
            if asv.startswith("_") or v.get("status") != "built":
                continue
            if (v.get("pfunc") or {}).get("applied"):
                jobs.append((dep, asv))
    C.log(f"scan: {len(jobs)} reweighted genomes")
    rows = []
    with Pool(procs, maxtasksperchild=50) as pool:
        for i, r in enumerate(pool.imap_unordered(_scan_one, sorted(jobs), chunksize=4), 1):
            rows.append(r)
            if i % 250 == 0:
                C.log(f"scan {i}/{len(jobs)}")
    out = {}
    for dep in C.BUILDS:
        R = [r for r in rows if r["dep"] == dep]
        out[dep] = {
            "n_genomes": len(R),
            "n_features": sum(r["n_features"] for r in R),
            "n_features_without_probability": sum(r["n_none"] for r in R),
            "n_replaced": sum(r["n_replaced"] for r in R),
            "n_replaced_value_changed": sum(r["n_changed"] for r in R),
            "min_probability_all_features": min(r["min_all"] for r in R),
            "max_probability_all_features": max(r["max_all"] for r in R),
            "n_probability_zero": sum(r["n_zero"] for r in R),
            "min_kn_of_replaced": min(r["min_kn_replaced"] for r in R if r["min_kn_replaced"] is not None),
            "min_e3": min(r["min_e3"] for r in R if r["min_e3"] is not None),
            "max_e3": max(r["max_e3"] for r in R if r["max_e3"] is not None),
            "n_e3_zero": sum(r["n_e3_zero"] for r in R),
            "cross_0.5_kn_below_to_e3_at_or_above": sum(r["cross05_up"] for r in R),
            "cross_0.5_kn_at_or_above_to_e3_below": sum(r["cross05_down"] for r in R),
            "replaced_kn_ge_0.5": sum(r["n_kn_ge05"] for r in R),
            "replaced_e3_ge_0.5": sum(r["n_e3_ge05"] for r in R),
            "genomes_with_any_0.5_crossing": sum(1 for r in R if r["cross05_up"] or r["cross05_down"]),
            "cds_probability_mismatch": sum(r["cds_mismatch"] for r in R),
            "zero_features": [z for r in sorted(R, key=lambda r: r["asv"]) for z in r["zero_features"]],
        }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--procs", type=int, default=8)
    a = ap.parse_args()
    table, n_bad = verify_table()
    hits = sweep()
    C.dump({"hits": hits}, f"{C.OUT}/results/census_sweep.json")
    res = {"_meta": {"script": "census.py", "code": C.code_record([os.path.abspath(__file__)]),
                     "pattern": PATTERN.pattern, "skipped_dirs_regex": SKIP_DIRS.pattern},
           "consumers": table, "n_anchor_failures": n_bad,
           "sweep_summary": classify_hits(hits, table),
           "emilykin_config_genome_objects": ek_config_paths(),
           "codiffusion_GDIR": read_line(CODIF_BUILD[0], CODIF_BUILD[1])}
    if a.scan:
        res["probability_scan"] = scan(a.procs)
    elif os.path.exists(f"{C.OUT}/results/census.json"):
        old = json.load(open(f"{C.OUT}/results/census.json"))
        if "probability_scan" in old:
            res["probability_scan"] = old["probability_scan"]
    C.dump(res, f"{C.OUT}/results/census.json")
    C.log(f"census: {len(table)} consumers, {n_bad} anchor failures, {len(hits)} sweep hits in "
          f"{res['sweep_summary']['n_files']} files")


if __name__ == "__main__":
    main()
