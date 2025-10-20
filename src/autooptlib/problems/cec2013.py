"""CEC2013 benchmark problems (subset) translated for the Python port."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Iterable, Sequence, Tuple

import numpy as np


_REPO_ROOT = Path(__file__).resolve().parents[3]
_DATA_DIR = (
    _REPO_ROOT
    / "AutoOptLib for Matlab"
    / "Problems"
    / "Numerial Benchmarks"
    / "CEC2013"
    / "data"
)
_SHIFT_DATA = None
_RNG = np.random.default_rng(42)


def _load_shift_data() -> np.ndarray:
    global _SHIFT_DATA
    if _SHIFT_DATA is None:
        path = _DATA_DIR / "shift_data.txt"
        _SHIFT_DATA = np.loadtxt(path, dtype=float)
    return _SHIFT_DATA


def _make_problem_entry(problem: Any, shift: np.ndarray) -> None:
    bound = np.vstack((np.full(shift.size, -100.0), np.full(shift.size, 100.0)))
    problem.type = ["continuous", "static", "certain"]
    problem.bound = bound
    problem.setting = getattr(problem, "setting", "")
    problem.dimension = shift.size
    problem.name = getattr(problem, "name", "cec2013_f1")

    def _evaluate(_data: Any, dec: np.ndarray) -> Tuple[float, float, None]:
        arr = np.asarray(dec, dtype=float)
        diff = arr - shift
        value = float(np.sum(diff**2) - 1400.0)
        return value, 0.0, None

    problem.evaluate = _evaluate


def cec2013_f1(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Python translation of the CEC2013 f1 benchmark."""
    mode = str(mode).lower()
    problems = list(problems)

    if mode == "construct":
        shift_data = _load_shift_data()
        data_entries = []
        for prob, inst in zip(problems, instances):
            dim = int(inst)
            if shift_data.shape[1] >= dim:
                shift = shift_data[0, :dim]
            else:
                shift = _RNG.uniform(-100.0, 100.0, size=dim)
            _make_problem_entry(prob, shift)
            data_entries.append(SimpleNamespace(o=shift))
        return problems, data_entries, None

    if mode == "repair":
        decs = np.asarray(instances, dtype=float)
        return decs, None, None

    if mode == "evaluate":
        data_obj = problems  # first argument is Data in MATLAB signature
        decs = np.asarray(instances, dtype=float)
        shift = np.asarray(getattr(data_obj, "o"), dtype=float)
        if decs.ndim == 1:
            diff = decs - shift
            value = float(np.sum(diff**2) - 1400.0)
            return value, 0.0, None
        diff = decs - shift
        values = np.sum(diff**2, axis=1) - 1400.0
        return values, np.zeros_like(values), None

    raise ValueError(f"Unsupported mode: {mode}")


__all__ = ["cec2013_f1"]

