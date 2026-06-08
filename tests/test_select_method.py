#!/usr/bin/env python
"""
test_select_method.py — exhaustive check of align_hits.select_method / decide, the
pure hardware -> alignment-method decision tree. No GPU, references, or I/O required.

Runs standalone (`python tests/test_select_method.py`) or under pytest.

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from align_hits import Decision, Hardware, decide, select_method  # noqa: E402

MIN = 10.0  # default free-disk threshold (GB)


def _chk(d: Decision, method, k1, save):
    assert d.method == method, f"method: got {d.method!r}, want {method!r} :: {d.reason}"
    assert d.k1 == k1, f"k1: got {d.k1}, want {k1} :: {d.reason}"
    assert d.save_edlib == save, f"save: got {d.save_edlib}, want {save} :: {d.reason}"


def test_cuda_wins_regardless_of_cpu():
    # CUDA takes priority over everything — even with abundant cores/disk or none.
    for cores in (1, 7, 8, 16, 31, 32, 64, 256):
        for disk in (0.0, 5.0, 10.0, 1e6):
            _chk(select_method(True, cores, disk), "gpu", None, False)


def test_cpu_16_to_32_with_disk_caches_top5000():
    for cores in (16, 20, 31):
        _chk(select_method(False, cores, 10.0), "edlib_biopython", 5000, True)  # exactly 10 GB
        _chk(select_method(False, cores, 5000.0), "edlib_biopython", 5000, True)


def test_cpu_16_to_32_low_disk_falls_back_to_2000_no_cache():
    for cores in (16, 24, 31):
        _chk(select_method(False, cores, 9.99), "edlib_biopython", 2000, False)
        _chk(select_method(False, cores, 0.0), "edlib_biopython", 2000, False)


def test_cpu_8_to_16_passes_top2000():
    for cores in (8, 12, 15):
        # disk is irrelevant in this tier
        _chk(select_method(False, cores, 0.0), "edlib_biopython", 2000, False)
        _chk(select_method(False, cores, 1e6), "edlib_biopython", 2000, False)


def test_cpu_below_8_passes_top1000():
    for cores in (1, 2, 7):
        _chk(select_method(False, cores, 1e6), "edlib_biopython", 1000, False)


def test_cpu_32_plus_no_gpu_is_full_biopython():
    for cores in (32, 33, 48, 64, 128):
        _chk(select_method(False, cores, 1e6), "biopython_full", None, False)


def test_boundaries():
    # 32 is the >=32 full-biopython boundary, not the edlib tier
    _chk(select_method(False, 32, 1e6), "biopython_full", None, False)
    _chk(select_method(False, 31, 1e6), "edlib_biopython", 5000, True)
    # 16 is the bottom of the 5000/cache tier
    _chk(select_method(False, 16, 10.0), "edlib_biopython", 5000, True)
    _chk(select_method(False, 15, 10.0), "edlib_biopython", 2000, False)
    # 8 is the bottom of the 2000 tier
    _chk(select_method(False, 8, 0.0), "edlib_biopython", 2000, False)
    _chk(select_method(False, 7, 0.0), "edlib_biopython", 1000, False)


def test_custom_min_disk_threshold():
    # a stricter 50 GB threshold pushes a 30 GB-free, 20-core box down to the 2000 tier
    d = select_method(False, 20, 30.0, min_disk_gb=50.0)
    _chk(d, "edlib_biopython", 2000, False)
    d = select_method(False, 20, 60.0, min_disk_gb=50.0)
    _chk(d, "edlib_biopython", 5000, True)


def test_decide_force_overrides():
    hw = Hardware(cuda_available=False, n_cores=64, free_disk_gb=1e6)  # would auto -> full
    _chk(decide(hw, force_method="gpu"), "gpu", None, False)
    _chk(decide(hw, force_method="biopython_full"), "biopython_full", None, False)
    # forced edlib still derives k1/cache from the (cores, disk) tier
    _chk(decide(hw, force_method="edlib_biopython"), "edlib_biopython", 5000, True)
    hw2 = Hardware(False, 10, 1e6)
    _chk(decide(hw2, force_method="edlib_biopython"), "edlib_biopython", 2000, False)


def test_decide_auto_matches_select():
    hw = Hardware(cuda_available=True, n_cores=8, free_disk_gb=4.0)
    assert decide(hw) == select_method(True, 8, 4.0)


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"  PASS  {fn.__name__}")
    print(f"\nAll {len(fns)} select_method tests passed.")


if __name__ == "__main__":
    _run_all()
