#!/usr/bin/env python
"""
test_select_method.py — checks align_hits.choose_backend, the pure hardware -> SW-scoring
backend decision. The prefilter is ALWAYS edlib top-k (hardware-independent); only the
Stage-2 scoring engine depends on hardware: CUDA > Metal > CPU, with --force-backend override.

Runs standalone (`python tests/test_select_method.py`) or under pytest.

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from align_hits import Hardware, choose_backend  # noqa: E402


def _hw(cuda, metal, cores=16, disk=100.0):
    return Hardware(cuda, metal, cores, disk)


def test_cuda_wins_over_everything():
    for metal in (False, True):
        for cores in (1, 8, 31, 32, 64, 256):
            assert choose_backend(_hw(True, metal, cores)).backend == "cuda"


def test_metal_when_no_cuda():
    for cores in (1, 8, 64):
        assert choose_backend(_hw(False, True, cores)).backend == "metal"


def test_cpu_when_no_gpu_any_core_count():
    for cores in (1, 7, 8, 16, 32, 64, 256):
        assert choose_backend(_hw(False, False, cores)).backend == "cpu"


def test_force_backend_overrides_hardware():
    for b in ("cuda", "metal", "cpu"):
        assert choose_backend(_hw(False, False), force=b).backend == b   # force despite no hw
        assert choose_backend(_hw(True, True), force=b).backend == b      # force despite all hw


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print("ok", fn.__name__)
    print(f"all {len(fns)} tests passed")
