"""mag_validation/af_gate_counterfactual.py reads its full-length 16S levels from the paired
MAG-vs-16S comparisons instead of hard-coding them (they went stale when the comparisons were
re-scored in 2026-11). Offline, on tiny fixtures; main() is not run.
"""
import importlib.util
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load():
    spec = importlib.util.spec_from_file_location(
        "afc_levels", os.path.join(ROOT, "mag_validation", "af_gate_counterfactual.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # module level only: imports, constants, functions
    return mod


def test_levels_come_from_the_paired_comparisons(tmp_path):
    mod = _load()
    vals = {"paired_comparison_d100_legacy.json": 0.692, "paired_comparison_d100.json": 0.7116,
            "paired_comparison.json": 0.7097}
    for name, v in vals.items():
        other = {"conditional": {"amp": 0.1}, "population": {"amp": 0.2}}
        (tmp_path / name).write_text(json.dumps({
            mod.LEVEL_KEY: {"conditional": {"n": 5, "mag": 0.8, "amp": v}, "population": {"amp": 0.5}},
            "exclude|c0.5_x0|FullLength16S|f1": other, "include|c1_x0|FullLength16S|f1": other}))
    got = mod.fl16s_levels(str(tmp_path))
    assert got == [("legacy reducer, depth 100", 0.692, "paired_comparison_d100_legacy.json"),
                   ("exact-tie union, depth 100", 0.7116, "paired_comparison_d100.json"),
                   ("depth 20", 0.7097, "paired_comparison.json")]


def test_no_hard_coded_level_in_the_scoring_code():
    src = open(os.path.join(ROOT, "mag_validation", "af_gate_counterfactual.py")).read()
    code = src.split('"""', 2)[2]          # past the module docstring, which records the old values
    for old in ("0.6912", "0.7092", "0.70915", "0.709)", "> 0.71"):
        assert old not in code
