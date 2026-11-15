"""Re-score code items W4 and W5, offline on tiny fixtures.

W4: manuscript_figures/p2/render_comparison_d100.py (and the two stock render_comparison.py
    copies) compute the figure's scoreable n and median-dF1 range from the summary instead of
    hard-coding the first depth-20 run's values.
W5: function_maps.py --msd / --mst pin the ModelSEED inputs and record the resolved git commits.
"""
import importlib.util
import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
RENDERERS = ["manuscript_figures/p2/render_comparison_d100.py",
             "region_validation/picrust2_comparison/scripts/render_comparison.py",
             "region_validation_d100/picrust2_comparison/scripts/render_comparison.py"]
FMAPS = ["region_validation/picrust2_comparison/scripts/function_maps.py",
         "region_validation_d100/picrust2_comparison/scripts/function_maps.py"]


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _summary(meds, n=67713, p_holm=0.0):
    regions = list(meds)
    pt = {f"{r}|prfba_exclude_vs_picrust2|conditional": {"median_diff": m, "p_value": p_holm,
                                                         "p_holm": p_holm} for r, m in meds.items()}
    return {"by_region": {"ec": dict({r: {} for r in regions},
                                     all={"all": {"prfba_exclude": {"n": n}, "picrust2": {"n": n}}})},
            "paired_tests_ec_f1": pt}


@pytest.mark.parametrize("path", RENDERERS)
def test_renderer_labels_are_computed(path):
    mod = _load(path, "r_" + path.replace("/", "_")[:-3])
    lab = mod.labels(_summary({"V4": 0.0716, "V1-V3": 0.0951, "FullLength16S": 0.091}))
    assert lab["n_text"] == "67,713" and lab["n_regions"] == 3
    assert lab["dF1_text"] == "+0.07 to +0.10" and lab["p_text"] == "p ≈ 0"
    lab = mod.labels(_summary({"V4": 0.0803, "V1-V3": 0.096}, n=67038, p_holm=0.003))
    assert lab["dF1_text"] == "+0.08 to +0.10" and lab["n_text"] == "67,038"
    assert lab["p_text"].startswith("Holm p")
    src = open(os.path.join(ROOT, path)).read()
    assert "67,038" not in src.split("def labels")[1] and "+0.08 to +0.10" not in src


def test_renderer_refuses_mismatched_n():
    mod = _load(RENDERERS[0], "r_mismatch")
    s = _summary({"V4": 0.07})
    s["by_region"]["ec"]["all"]["all"]["picrust2"]["n"] = 1
    with pytest.raises(SystemExit):
        mod.labels(s)


def _tiny_modelseed(tmp_path):
    msd = tmp_path / "MSD"
    (msd / "Biochemistry" / "Aliases").mkdir(parents=True)
    (msd / "Biochemistry" / "reaction_00.tsv").write_text(
        "id\tname\tis_obsolete\nrxn00001\ta\t0\nrxn00002\tb\t1\n")
    (msd / "Biochemistry" / "Aliases" / "Unique_ModelSEED_Reaction_ECs.txt").write_text(
        "ModelSEED ID\tExternal ID\tSource\nrxn00001\t1.1.1.1\tEnzyme Class\n"
        "rxn00002\t1.1.1.1\tEnzyme Class\nrxn00001\t2.7.-.-\tEnzyme Class\n")
    mst = tmp_path / "MST" / "templates" / "v7.0"
    mst.mkdir(parents=True)
    tmpl = {"roles": [{"id": "ftr1", "name": "Alcohol dehydrogenase (EC 1.1.1.1)"}],
            "complexes": [{"id": "cpx1", "complexroles": [{"templaterole_ref": "~/roles/id/ftr1",
                                                           "triggering": 1}]}],
            "reactions": [{"id": "rxn00001_c", "templatecomplex_refs": ["~/complexes/id/cpx1"]}]}
    for name in ("GramNegModelTemplateV7.json", "GramPosModelTemplateV7.json", "ArchaeaTemplateV6.json"):
        (mst / name).write_text(json.dumps(tmpl))
    prod = tmp_path / "products.json"
    prod.write_text(json.dumps({"PGF_1": {"family_product": "Alcohol dehydrogenase (EC 1.1.1.1)"},
                                "PGF_2": {"family_product": "hypothetical protein"}}))
    return msd, mst, prod


@pytest.mark.parametrize("path", FMAPS)
def test_function_maps_pins_and_records_commits(tmp_path, path):
    msd, mst, prod = _tiny_modelseed(tmp_path)
    env = dict(os.environ, GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t", GIT_COMMITTER_NAME="t",
               GIT_COMMITTER_EMAIL="t@t")
    subprocess.run(["git", "init", "-q", str(msd)], check=True, env=env)
    subprocess.run(["git", "-C", str(msd), "add", "-A"], check=True, env=env)
    subprocess.run(["git", "-C", str(msd), "commit", "-qm", "x"], check=True, env=env)
    head = subprocess.run(["git", "-C", str(msd), "rev-parse", "HEAD"], capture_output=True,
                          text=True, check=True).stdout.strip()
    out = tmp_path / "maps.json"
    r = subprocess.run([PY, os.path.join(ROOT, path), "--products", str(prod), "--out", str(out),
                        "--msd", str(msd), "--mst", str(mst)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    d = json.load(open(out))
    assert d["pgfam_ec"] == {"PGF_1": ["1.1.1.1"]} and d["pgfam_rxn"] == {"PGF_1": ["rxn00001"]}
    assert d["ec_rxn"] == {"1.1.1.1": ["rxn00001"]}                 # obsolete + partial EC dropped
    m = d["_meta"]
    assert m["modelseed_db"] == str(msd) and m["modelseed_templates_dir"] == str(mst)
    assert m["modelseed_db_git"]["commit"] == head and m["modelseed_db_git"]["dirty"] == 0
    assert m["modelseed_templates_git"]["commit"] is None           # not a git checkout


@pytest.mark.parametrize("path", FMAPS)
def test_function_maps_rejects_a_wrong_pin(tmp_path, path):
    msd, mst, prod = _tiny_modelseed(tmp_path)
    r = subprocess.run([PY, os.path.join(ROOT, path), "--products", str(prod),
                        "--out", str(tmp_path / "m.json"), "--msd", str(mst), "--mst", str(mst)],
                       capture_output=True, text=True)
    assert r.returncode != 0 and "not a ModelSEEDDatabase checkout" in (r.stdout + r.stderr)


def test_driver_passes_fmaps_args_through():
    s = open(os.path.join(ROOT, "region_validation_d100", "run_picrust2_d100.sh")).read()
    line = [ln for ln in s.splitlines() if "function_maps.py" in ln and ln.startswith("run ")]
    assert len(line) == 1 and line[0].rstrip().endswith("${FMAPS_ARGS:-}")
