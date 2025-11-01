"""CEC2013 benchmark problems translated for the Python port.

This module currently exposes a growing subset of the official
benchmark suite.  The implementations are direct ports of AutoOptLib's
Matlab definitions, sharing the same data files (shift vectors,
rotation matrices, etc.) so that fitness values remain comparable.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Iterable, List, Sequence, Tuple

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
_ROTATION_DATA: dict[int, np.ndarray] = {}
_RNG = np.random.default_rng(42)

# Rotation matrices are bundled for a fixed set of dimensions.  The text
# files store multiple stacked matrices (each block is D rows).
_ROTATION_FILES = {
    2: "M_D2.txt",
    5: "M_D5.txt",
    10: "M_D10.txt",
    20: "M_D20.txt",
    30: "M_D30.txt",
    40: "M_D40.txt",
    50: "M_D50.txt",
    60: "M_D60.txt",
    70: "M_D70.txt",
    80: "M_D80.txt",
    90: "M_D90.txt",
    100: "M_D100.txt",
}


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def _load_shift_data() -> np.ndarray:
    global _SHIFT_DATA
    if _SHIFT_DATA is None:
        path = _DATA_DIR / "shift_data.txt"
        _SHIFT_DATA = np.loadtxt(path, dtype=float)
        if _SHIFT_DATA.ndim == 1:
            _SHIFT_DATA = _SHIFT_DATA.reshape(1, -1)
    return _SHIFT_DATA


def _load_rotation_blocks(dim: int) -> np.ndarray | None:
    """Return stacked rotation matrices for the requested dimension."""
    if dim not in _ROTATION_FILES:
        return None

    if dim not in _ROTATION_DATA:
        path = _DATA_DIR / _ROTATION_FILES[dim]
        raw = np.loadtxt(path, dtype=float)
        if raw.ndim == 1:
            raw = raw.reshape(1, -1)
        rows, cols = raw.shape
        if cols < dim or rows % dim != 0:
            raise ValueError(
                f"Rotation data '{path.name}' has shape {raw.shape}, "
                f"cannot extract {dim}x{dim} blocks."
            )
        blocks = raw[:, :dim].reshape(rows // dim, dim, dim)
        _ROTATION_DATA[dim] = blocks
    return _ROTATION_DATA[dim]


def _random_orthonormal(dim: int) -> np.ndarray:
    """Fallback rotation when no precomputed matrix exists."""
    mat = _RNG.normal(size=(dim, dim))
    q, _ = np.linalg.qr(mat)
    # Ensure a proper rotation matrix (determinant +1).
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1.0
    return q


def _get_rotation(dim: int, index: int = 0) -> np.ndarray:
    blocks = _load_rotation_blocks(dim)
    if blocks is not None and index < blocks.shape[0]:
        return blocks[index].copy()
    return _random_orthonormal(dim)


def _construct_lambda(alpha: float, dim: int) -> np.ndarray:
    if dim <= 1:
        return np.eye(dim, dtype=float)
    indices = np.arange(dim, dtype=float)
    diag = alpha ** (indices / (2.0 * (dim - 1)))
    return np.diag(diag)


def _compute_t_osz(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    result = np.zeros_like(values)
    mask = values != 0.0
    if not np.any(mask):
        return result

    v = values[mask]
    x_hat = np.log(np.abs(v))
    sign = np.sign(v)
    sign[sign == 0.0] = 1.0

    c1 = np.where(v > 0.0, 10.0, 5.5)
    c2 = np.where(v > 0.0, 7.9, 3.1)
    transformed = np.exp(x_hat + 0.049 * (np.sin(c1 * x_hat) + np.sin(c2 * x_hat)))
    result[mask] = sign * transformed
    return result


def _compute_t_asym(values: np.ndarray, beta: float) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
        squeeze_back = True
    else:
        squeeze_back = False

    result = arr.copy()
    dim = arr.shape[1]
    if dim == 1:
        return result.squeeze(0) if squeeze_back else result

    exponents = 1.0 + beta * np.arange(dim) / (dim - 1)
    for j in range(dim):
        pos = arr[:, j] > 0.0
        if not np.any(pos):
            continue
        result[pos, j] = np.power(arr[pos, j], exponents[j]) * np.sqrt(arr[pos, j])

    return result.squeeze(0) if squeeze_back else result


def _prepare_shift(dim: int) -> np.ndarray:
    shift_data = _load_shift_data()
    if shift_data.shape[1] >= dim:
        return shift_data[0, :dim].copy()
    return _RNG.uniform(-100.0, 100.0, size=dim)


def _initialise_problem(problem: Any, shift: np.ndarray, name: str) -> None:
    bound = np.vstack((np.full(shift.size, -100.0), np.full(shift.size, 100.0)))
    problem.type = ["continuous", "static", "certain"]
    problem.bound = bound
    problem.setting = getattr(problem, "setting", "")
    problem.dimension = shift.size
    problem.name = getattr(problem, "name", name)


def _ensure_array(instances: Any) -> Tuple[np.ndarray, bool]:
    arr = np.asarray(instances, dtype=float)
    if arr.ndim == 1:
        return arr.reshape(1, -1), True
    return arr, False


def _finalise(values: np.ndarray, single: bool) -> Tuple[Any, Any, None]:
    values = np.asarray(values, dtype=float)
    if single:
        return float(values[0]), 0.0, None
    zeros = np.zeros_like(values)
    return values, zeros, None


def _construct_common(
    problems: Iterable[Any],
    instances: Sequence[int],
    name: str,
    make_data: Callable[[int, np.ndarray], SimpleNamespace],
) -> Tuple[List[Any], List[SimpleNamespace], None]:
    problem_list = list(problems)
    data_entries: List[SimpleNamespace] = []
    for prob, inst in zip(problem_list, instances):
        dim = int(inst)
        shift = _prepare_shift(dim)
        data = make_data(dim, shift)
        data_entries.append(data)
        _initialise_problem(prob, shift, name)

        def _make_eval(data_ref: SimpleNamespace) -> Any:
            def _evaluate(_data: Any, dec: np.ndarray) -> Tuple[Any, Any, None]:
                return globals()[name](data_ref, dec, "evaluate")

            return _evaluate

        prob.evaluate = _make_eval(data)
    return problem_list, data_entries, None


def _repair(instances: Any) -> Tuple[np.ndarray, None, None]:
    decs = np.asarray(instances, dtype=float)
    return decs, None, None


def _mode_guard(mode: str) -> str:
    mode = str(mode).lower()
    if mode not in {"construct", "repair", "evaluate"}:
        raise ValueError(f"Unsupported mode: {mode}")
    return mode


# ---------------------------------------------------------------------------
# Individual problems
# ---------------------------------------------------------------------------

def cec2013_f1(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted sphere."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f1",
            lambda dim, shift: SimpleNamespace(o=shift),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    diff = decs - shift
    values = np.sum(diff**2, axis=1) - 1400.0
    return _finalise(values, single)


def cec2013_f2(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted and rotated Schwefel 1.2 with osz transformation."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f2",
            lambda dim, shift: SimpleNamespace(o=shift, M=_get_rotation(dim)),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M = np.asarray(data_obj.M, dtype=float)

    z = (decs - shift) @ M
    z = _compute_t_osz(z)

    dim = z.shape[1]
    if dim == 1:
        weights = np.array([1.0], dtype=float)
    else:
        weights = (1e6) ** (np.arange(dim) / (dim - 1))
    values = np.sum(weights * (z**2), axis=1) - 1300.0
    return _finalise(values, single)


def cec2013_f3(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted and rotated Rosenbrock with asymmetry transformation."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f3",
            lambda dim, shift: SimpleNamespace(
                o=shift, M1=_get_rotation(dim, 0), M2=_get_rotation(dim, 1)
            ),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M1 = np.asarray(data_obj.M1, dtype=float)
    M2 = np.asarray(data_obj.M2, dtype=float)

    z = (decs - shift) @ M1
    z = _compute_t_asym(z, 0.5)
    z = z @ M2

    a = 1e6
    # f3 = z1^2 + a * sum_{i=2..D} z_i^2
    if z.shape[1] > 1:
        tail = np.sum(z[:, 1:] ** 2, axis=1)
    else:
        tail = np.zeros(z.shape[0])
    values = z[:, 0] ** 2 + a * tail
    values -= 1200.0
    return _finalise(values, single)


def cec2013_f4(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted Schwefel 1.2 with osz transformation."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f4",
            lambda dim, shift: SimpleNamespace(o=shift, M=_get_rotation(dim, 0)),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M = np.asarray(data_obj.M, dtype=float)

    z = (decs - shift) @ M
    z = _compute_t_osz(z)

    a = 1e6
    # f4 = a*z1^2 + sum_{i=2..D} z_i^2
    if z.shape[1] > 1:
        tail = np.sum(z[:, 1:] ** 2, axis=1)
    else:
        tail = np.zeros(z.shape[0])
    values = a * (z[:, 0] ** 2) + tail
    values -= 1100.0
    return _finalise(values, single)


def cec2013_f5(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted Ackley variant with non-uniform exponent scaling."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f5",
            lambda dim, shift: SimpleNamespace(o=shift),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)

    z = decs - shift
    dim = z.shape[1]
    if dim == 1:
        powers = np.array([2.0])
    else:
        powers = 2.0 + 4.0 * (np.arange(dim) / (dim - 1))
    values = np.sqrt(np.sum(np.abs(z) ** powers, axis=1)) - 1000.0
    return _finalise(values, single)


def cec2013_f6(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted and rotated Rosenbrock."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f6",
            lambda dim, shift: SimpleNamespace(o=shift, M=_get_rotation(dim, 0)),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M = np.asarray(data_obj.M, dtype=float)

    z = (2.048 * (decs - shift) / 100.0) @ M + 1.0
    left = z[:, :-1]
    right = z[:, 1:]
    values = np.sum(100.0 * (left**2 - right) ** 2 + (left - 1.0) ** 2, axis=1)
    values -= 900.0
    return _finalise(values, single)


def cec2013_f7(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted and rotated Expanded Schaffer's F7 function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f7",
            lambda dim, shift: SimpleNamespace(
                o=shift, M1=_get_rotation(dim, 0), M2=_get_rotation(dim, 1)
            ),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M1 = np.asarray(data_obj.M1, dtype=float)
    M2 = np.asarray(data_obj.M2, dtype=float)

    z = (decs - shift) @ M1
    z = _compute_t_asym(z, 0.5)
    z = z @ M2
    z = z @ _construct_lambda(10.0, z.shape[1])

    if z.shape[1] < 2:
        values = np.zeros(z.shape[0])
    else:
        left = z[:, :-1]
        right = z[:, 1:]
        r = np.sqrt(left**2 + right**2)
        sqrt_r = np.sqrt(r)
        inner = sqrt_r + sqrt_r * (np.sin(50.0 * r**0.2) ** 2)
        # f7 = ( mean(inner) )^2
        values = np.mean(inner, axis=1) ** 2
    values -= 800.0
    return _finalise(values, single)


def cec2013_f8(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted and rotated Ackley function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f8",
            lambda dim, shift: SimpleNamespace(
                o=shift, M1=_get_rotation(dim, 0), M2=_get_rotation(dim, 1)
            ),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M1 = np.asarray(data_obj.M1, dtype=float)
    M2 = np.asarray(data_obj.M2, dtype=float)

    z = (decs - shift) @ M1
    z = _compute_t_asym(z, 0.5)
    z = z @ M2
    z = z @ _construct_lambda(10.0, z.shape[1])

    dim = z.shape[1]
    term1 = np.sum(z**2, axis=1) / dim
    term2 = np.sum(np.cos(2.0 * np.pi * z), axis=1) / dim
    values = -20.0 * np.exp(-0.2 * np.sqrt(term1)) - np.exp(term2) + 20.0 + np.e
    values -= 700.0
    return _finalise(values, single)


def cec2013_f9(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted and rotated Weierstrass function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f9",
            lambda dim, shift: SimpleNamespace(
                o=shift, M1=_get_rotation(dim, 0), M2=_get_rotation(dim, 1)
            ),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M1 = np.asarray(data_obj.M1, dtype=float)
    M2 = np.asarray(data_obj.M2, dtype=float)

    z = 0.5 * (decs - shift) / 100.0
    z = z @ M1
    z = _compute_t_asym(z, 0.5)
    z = z @ M2
    z = z @ _construct_lambda(10.0, z.shape[1])

    a = 0.5
    b = 3.0
    k = np.arange(0, 21)
    a_k = a**k
    b_k = b**k
    inner = z[:, :, None] + 0.5
    term1 = np.sum(
        np.sum(a_k[np.newaxis, np.newaxis, :] * np.cos(2.0 * np.pi * b_k[np.newaxis, np.newaxis, :] * inner), axis=2),
        axis=1,
    )
    constant = np.sum(a_k * np.cos(2.0 * np.pi * b_k * 0.5))
    values = term1 - z.shape[1] * constant
    values -= 600.0
    return _finalise(values, single)


def cec2013_f10(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted and rotated Griewank function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f10",
            lambda dim, shift: SimpleNamespace(o=shift, M=_get_rotation(dim, 0)),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M = np.asarray(data_obj.M, dtype=float)

    z = 600.0 * (decs - shift) / 100.0
    z = z @ M
    z = z @ _construct_lambda(100.0, z.shape[1])

    sum_term = np.sum(z**2, axis=1) / 4000.0
    indices = np.sqrt(np.arange(1, z.shape[1] + 1, dtype=float))
    cos_terms = np.prod(np.cos(z / indices), axis=1)
    values = sum_term - cos_terms + 1.0
    values -= 500.0
    return _finalise(values, single)


def cec2013_f11(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted Rastrigin function with asymmetric transformation."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f11",
            lambda dim, shift: SimpleNamespace(o=shift),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)

    z = 5.12 * (decs - shift) / 100.0
    z = _compute_t_osz(z)
    z = _compute_t_asym(z, 0.2)
    z = z @ _construct_lambda(10.0, z.shape[1])

    values = np.sum(z**2 - 10.0 * np.cos(2.0 * np.pi * z) + 10.0, axis=1)
    values -= 400.0
    return _finalise(values, single)


def cec2013_f12(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted and rotated Rastrigin function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f12",
            lambda dim, shift: SimpleNamespace(
                o=shift, M1=_get_rotation(dim, 0), M2=_get_rotation(dim, 1)
            ),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M1 = np.asarray(data_obj.M1, dtype=float)
    M2 = np.asarray(data_obj.M2, dtype=float)

    z = 5.12 * (decs - shift) / 100.0
    z = z @ M1
    z = _compute_t_osz(z)
    z = _compute_t_asym(z, 0.2)
    z = z @ M2
    z = z @ _construct_lambda(10.0, z.shape[1])
    z = z @ M1

    values = np.sum(z**2 - 10.0 * np.cos(2.0 * np.pi * z) + 10.0, axis=1)
    values -= 300.0
    return _finalise(values, single)


def cec2013_f13(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Non-continuous shifted and rotated Rastrigin function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f13",
            lambda dim, shift: SimpleNamespace(
                o=shift, M1=_get_rotation(dim, 0), M2=_get_rotation(dim, 1)
            ),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M1 = np.asarray(data_obj.M1, dtype=float)
    M2 = np.asarray(data_obj.M2, dtype=float)

    z = 5.12 * (decs - shift) / 100.0
    z = z @ M1

    mask = np.abs(z) > 0.5
    rounded = np.round(2.0 * z[mask]) / 2.0
    z = z.copy()
    z[mask] = rounded

    z = _compute_t_osz(z)
    z = _compute_t_asym(z, 0.2)
    z = z @ M2
    z = z @ _construct_lambda(10.0, z.shape[1])
    z = z @ M1

    values = np.sum(z**2 - 10.0 * np.cos(2.0 * np.pi * z) + 10.0, axis=1)
    values -= 200.0
    return _finalise(values, single)


def _schwefel_core(z: np.ndarray, dim: int, penalty_sign: float) -> np.ndarray:
    """Helper for Schwefel-like functions with boundary handling."""
    abs_z = np.abs(z)
    result = np.zeros(z.shape[0], dtype=float)
    for i in range(dim):
        zi = z[:, i]
        term = np.zeros_like(result)
        cond1 = abs_z[:, i] <= 500.0
        term[cond1] = zi[cond1] * np.sin(np.sqrt(abs_z[cond1]))

        cond2 = zi > 500.0
        if np.any(cond2):
            temp = 500.0 - np.mod(zi[cond2], 500.0)
            penalty = penalty_sign * ((zi[cond2] - 500.0) ** 2) / (10000.0 * dim)
            term[cond2] = temp * np.sin(np.sqrt(np.abs(temp))) + penalty

        cond3 = zi < -500.0
        if np.any(cond3):
            temp = np.mod(abs_z[cond3], 500.0) - 500.0
            penalty = penalty_sign * ((zi[cond3] + 500.0) ** 2) / (10000.0 * dim)
            term[cond3] = temp * np.sin(np.sqrt(np.abs(temp))) + penalty

        result += term
    return result


def cec2013_f14(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted Schwefel function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f14",
            lambda dim, shift: SimpleNamespace(o=shift),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)

    z = 10.0 * (decs - shift)
    z = z @ _construct_lambda(10.0, z.shape[1])
    z = z + 4.209687462275036e2

    g_sum = _schwefel_core(z, z.shape[1], penalty_sign=-1.0)
    values = 418.9829 * z.shape[1] - g_sum
    values -= 100.0
    return _finalise(values, single)


def cec2013_f15(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted and rotated Schwefel function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f15",
            lambda dim, shift: SimpleNamespace(o=shift, M1=_get_rotation(dim, 0)),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M1 = np.asarray(data_obj.M1, dtype=float)

    z = 10.0 * (decs - shift)
    z = z @ M1
    z = z @ _construct_lambda(10.0, z.shape[1])
    z = z + 4.209687462275036e2

    g_sum = _schwefel_core(z, z.shape[1], penalty_sign=1.0)
    values = 418.9829 * z.shape[1] - g_sum
    values += 100.0
    return _finalise(values, single)


def cec2013_f16(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Shifted and rotated Katsuura function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f16",
            lambda dim, shift: SimpleNamespace(
                o=shift, M1=_get_rotation(dim, 0), M2=_get_rotation(dim, 1)
            ),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    shift = np.asarray(data_obj.o, dtype=float)
    M1 = np.asarray(data_obj.M1, dtype=float)
    M2 = np.asarray(data_obj.M2, dtype=float)

    z = (5.0 / 100.0) * (decs - shift)
    z = z @ M1
    z = z @ _construct_lambda(100.0, z.shape[1])
    z = z @ M2

    exponent = 10.0 / (z.shape[1] ** 1.2)
    product_val = np.ones(z.shape[0], dtype=float)
    j = np.arange(1, 33, dtype=float)
    two_pow_j = 2.0**j
    inv_two_pow_j = 1.0 / two_pow_j

    for i in range(z.shape[1]):
        zi = z[:, i][:, np.newaxis]
        frac = np.abs(two_pow_j * zi - np.round(two_pow_j * zi)) * inv_two_pow_j
        sum_j = np.sum(frac, axis=1)
        inner = 1.0 + (i + 1) * sum_j
        product_val *= inner**exponent

    values = 10.0 / (z.shape[1] ** 2) * product_val - 10.0 / (z.shape[1] ** 2)
    values += 200.0
    return _finalise(values, single)


def cec2013_f17(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Lunacek bi-Rastrigin function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f17",
            lambda dim, shift: SimpleNamespace(o=shift),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    o = np.asarray(data_obj.o, dtype=float)

    mu0 = 2.5
    d = 1.0
    s = 1.0 - 1.0 / (2.0 * np.sqrt(decs.shape[1] + 20.0) - 8.2)
    mu1 = -np.sqrt((mu0**2 - d) / s)

    y = 10.0 * (decs - o) / 100.0
    sign_vec = np.where(o > 0.0, 1.0, -1.0)
    xhat = 2.0 * sign_vec * y + mu0
    z = xhat - mu0
    z = z @ _construct_lambda(100.0, z.shape[1])

    term1 = np.sum((xhat - mu0) ** 2, axis=1)
    term2 = d * decs.shape[1] + s * np.sum((xhat - mu1) ** 2, axis=1)
    base = np.minimum(term1, term2)
    ras = 10.0 * (decs.shape[1] - np.sum(np.cos(2.0 * np.pi * z), axis=1))
    values = base + ras + 300.0
    return _finalise(values, single)


def cec2013_f18(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Rotated Lunacek bi-Rastrigin function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f18",
            lambda dim, shift: SimpleNamespace(o=shift, M1=_get_rotation(dim, 0), M2=_get_rotation(dim, 1)),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    o = np.asarray(data_obj.o, dtype=float)
    M1 = np.asarray(data_obj.M1, dtype=float)
    M2 = np.asarray(data_obj.M2, dtype=float)

    mu0 = 2.5
    d = 1.0
    s = 1.0 - 1.0 / (2.0 * np.sqrt(decs.shape[1] + 20.0) - 8.2)
    mu1 = -np.sqrt((mu0**2 - d) / s)

    y = 10.0 * (decs - o) / 100.0
    sign_vec = np.where(o > 0.0, 1.0, -1.0)
    xhat = 2.0 * sign_vec * y + mu0
    z = (xhat - mu0) @ M1
    z = z @ _construct_lambda(100.0, z.shape[1])
    z = z @ M2

    term1 = np.sum((xhat - mu0) ** 2, axis=1)
    term2 = d * decs.shape[1] + s * np.sum((xhat - mu1) ** 2, axis=1)
    base = np.minimum(term1, term2)
    ras = 10.0 * (decs.shape[1] - np.sum(np.cos(2.0 * np.pi * z), axis=1))
    values = base + ras + 400.0
    return _finalise(values, single)


def cec2013_f19(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Rotated Expanded Griewank's plus Rosenbrock's Function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f19",
            lambda dim, shift: SimpleNamespace(o=shift, M=_get_rotation(dim, 0)),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    o = np.asarray(data_obj.o, dtype=float)
    M = np.asarray(data_obj.M, dtype=float)

    z = 5.0 * (decs - o) / 100.0
    z = z @ M
    z = z + 1.0

    N, D = z.shape
    # Rosenbrock pair then Griewank applied to each scalar
    xi = z[:, :-1]
    xj = z[:, 1:]
    last_pair = (z[:, -1], z[:, 0])

    def rosen_pair(a, b):
        return 100.0 * (a**2 - b) ** 2 + (a - 1.0) ** 2

    parts = rosen_pair(xi, xj)
    parts_last = rosen_pair(last_pair[0], last_pair[1])[:, np.newaxis]
    g2_vals = np.hstack([parts, parts_last])

    # Griewank on each scalar u: u/4000 - cos(u) + 1
    g1_vals = g2_vals / 4000.0 - np.cos(g2_vals) + 1.0
    values = np.sum(g1_vals, axis=1) + 500.0
    return _finalise(values, single)


def cec2013_f20(problems: Iterable[Any], instances: Sequence[int], mode: str):
    """Rotated Expanded Scaffer's F6 function."""
    mode = _mode_guard(mode)
    if mode == "construct":
        return _construct_common(
            problems,
            instances,
            "cec2013_f20",
            lambda dim, shift: SimpleNamespace(o=shift, M1=_get_rotation(dim, 0), M2=_get_rotation(dim, 1)),
        )

    if mode == "repair":
        return _repair(instances)

    data_obj = problems if hasattr(problems, "o") else list(problems)[0]
    decs, single = _ensure_array(instances)
    o = np.asarray(data_obj.o, dtype=float)
    M1 = np.asarray(data_obj.M1, dtype=float)
    M2 = np.asarray(data_obj.M2, dtype=float)

    z = (decs - o) @ M1
    z = _compute_t_asym(z, 0.5)
    z = z @ M2

    def g(x, y):
        r2 = x**2 + y**2
        return 0.5 + (np.sin(np.sqrt(r2)) ** 2 - 0.5) / (1.0 + 0.001 * r2) ** 2

    xi = z[:, :-1]
    xj = z[:, 1:]
    term = g(xi, xj)
    last = g(z[:, -1], z[:, 0])
    values = np.sum(term, axis=1) + last + 600.0
    return _finalise(values, single)


# ---------------------------- Composition functions ------------------------

_OFFSETS = {
    1: -1400.0,
    2: -1300.0,
    3: -1200.0,
    4: -1100.0,
    5: -1000.0,
    6: -900.0,
    7: -800.0,
    8: -700.0,
    9: -600.0,
    10: -500.0,
    11: -400.0,
    12: -300.0,
    13: -200.0,
    14: -100.0,
    15: 100.0,
    16: 200.0,
    17: 300.0,
    18: 400.0,
    19: 500.0,
    20: 600.0,
}


def _eval_base(i: int, x: np.ndarray, o: np.ndarray, rot_idx: int) -> np.ndarray:
    # Evaluate base function f_i and return f' = f - f* (i.e. add offset magnitude)
    if i == 1:
        data = SimpleNamespace(o=o)
        val, _, _ = cec2013_f1(data, x, "evaluate")
        return val - _OFFSETS[1]
    if i == 2:
        data = SimpleNamespace(o=o, M=_get_rotation(x.shape[1], rot_idx))
        val, _, _ = cec2013_f2(data, x, "evaluate")
        return val - _OFFSETS[2]
    if i == 3:
        data = SimpleNamespace(
            o=o, M1=_get_rotation(x.shape[1], rot_idx), M2=_get_rotation(x.shape[1], rot_idx + 1)
        )
        val, _, _ = cec2013_f3(data, x, "evaluate")
        return val - _OFFSETS[3]
    if i == 4:
        data = SimpleNamespace(o=o, M=_get_rotation(x.shape[1], rot_idx))
        val, _, _ = cec2013_f4(data, x, "evaluate")
        return val - _OFFSETS[4]
    if i == 5:
        data = SimpleNamespace(o=o)
        val, _, _ = cec2013_f5(data, x, "evaluate")
        return val - _OFFSETS[5]
    if i == 6:
        data = SimpleNamespace(o=o, M=_get_rotation(x.shape[1], rot_idx))
        val, _, _ = cec2013_f6(data, x, "evaluate")
        return val - _OFFSETS[6]
    if i == 7:
        data = SimpleNamespace(
            o=o, M1=_get_rotation(x.shape[1], rot_idx), M2=_get_rotation(x.shape[1], rot_idx + 1)
        )
        val, _, _ = cec2013_f7(data, x, "evaluate")
        return val - _OFFSETS[7]
    if i == 8:
        data = SimpleNamespace(
            o=o, M1=_get_rotation(x.shape[1], rot_idx), M2=_get_rotation(x.shape[1], rot_idx + 1)
        )
        val, _, _ = cec2013_f8(data, x, "evaluate")
        return val - _OFFSETS[8]
    if i == 9:
        data = SimpleNamespace(
            o=o, M1=_get_rotation(x.shape[1], rot_idx), M2=_get_rotation(x.shape[1], rot_idx + 1)
        )
        val, _, _ = cec2013_f9(data, x, "evaluate")
        return val - _OFFSETS[9]
    if i == 10:
        data = SimpleNamespace(o=o, M=_get_rotation(x.shape[1], rot_idx))
        val, _, _ = cec2013_f10(data, x, "evaluate")
        return val - _OFFSETS[10]
    if i == 11:
        data = SimpleNamespace(o=o)
        val, _, _ = cec2013_f11(data, x, "evaluate")
        return val - _OFFSETS[11]
    if i == 12:
        data = SimpleNamespace(
            o=o, M1=_get_rotation(x.shape[1], rot_idx), M2=_get_rotation(x.shape[1], rot_idx + 1)
        )
        val, _, _ = cec2013_f12(data, x, "evaluate")
        return val - _OFFSETS[12]
    if i == 14:
        data = SimpleNamespace(o=o)
        val, _, _ = cec2013_f14(data, x, "evaluate")
        return val - _OFFSETS[14]
    if i == 15:
        data = SimpleNamespace(
            o=o, M1=_get_rotation(x.shape[1], rot_idx), M2=_get_rotation(x.shape[1], rot_idx + 1)
        )
        val, _, _ = cec2013_f15(data, x, "evaluate")
        return val - _OFFSETS[15]
    if i == 16:
        data = SimpleNamespace(
            o=o, M1=_get_rotation(x.shape[1], rot_idx), M2=_get_rotation(x.shape[1], rot_idx + 1)
        )
        val, _, _ = cec2013_f16(data, x, "evaluate")
        return val - _OFFSETS[16]
    if i == 19:
        data = SimpleNamespace(o=o, M=_get_rotation(x.shape[1], rot_idx))
        val, _, _ = cec2013_f19(data, x, "evaluate")
        return val - _OFFSETS[19]
    if i == 20:
        data = SimpleNamespace(
            o=o, M1=_get_rotation(x.shape[1], rot_idx), M2=_get_rotation(x.shape[1], rot_idx + 1)
        )
        val, _, _ = cec2013_f20(data, x, "evaluate")
        return val - _OFFSETS[20]
    raise ValueError(f"Unsupported base function id: {i}")


def _composition_value(x: np.ndarray, centers: np.ndarray, base_ids: Sequence[int], lambdas: Sequence[float], sigmas: Sequence[float], biases: Sequence[float], rot_start: int = 0) -> np.ndarray:
    N, D = x.shape
    n = len(base_ids)
    w = np.zeros((N, n), dtype=float)
    gvals = np.zeros((N, n), dtype=float)
    for i, fid in enumerate(base_ids):
        oi = centers[i, :D]
        diff = x - oi
        dist2 = np.sum(diff**2, axis=1)
        # raw g' value
        gvals[:, i] = _eval_base(fid, x, oi, rot_start + 2 * i)
        # weights
        denom = np.sqrt(dist2 + 1e-30)
        w[:, i] = (1.0 / denom) * np.exp(-dist2 / (2.0 * D * (sigmas[i] ** 2) + 1e-30))
    w_sum = np.sum(w, axis=1, keepdims=True) + 1e-30
    omega = w / w_sum
    lambdas = np.asarray(lambdas, dtype=float)
    biases = np.asarray(biases, dtype=float)
    values = np.sum(omega * (lambdas[np.newaxis, :] * gvals + biases[np.newaxis, :]), axis=1)
    return values


def _construct_composition(problems: Iterable[Any], instances: Sequence[int], name: str):
    problem_list = list(problems)
    data_entries: List[SimpleNamespace] = []
    for prob, inst in zip(problem_list, instances):
        dim = int(inst)
        shift = _prepare_shift(dim)  # for bounds only
        _initialise_problem(prob, shift, name)
        data_entries.append(SimpleNamespace(dim=dim))
        def _eval_closure(_):
            def inner(_d, dec):
                return globals()[name](SimpleNamespace(dim=dim), dec, 'evaluate')
            return inner
        prob.evaluate = _eval_closure(dim)
    return problem_list, data_entries, None


def cec2013_f21(problems: Iterable[Any], instances: Sequence[int], mode: str):
    mode = _mode_guard(mode)
    if mode == 'construct':
        return _construct_composition(problems, instances, 'cec2013_f21')
    if mode == 'repair':
        return _repair(instances)

    decs, single = _ensure_array(instances)
    shift_data = _load_shift_data()
    centers = shift_data[:5, :]
    base_ids = [6, 5, 3, 4, 1]
    sigmas = [10, 20, 30, 40, 50]
    lambdas = [1, 1e-6, 1e-26, 1e-6, 0.1]
    biases = [0, 100, 200, 300, 400]
    values = _composition_value(decs, centers, base_ids, lambdas, sigmas, biases) + 700.0
    return _finalise(values, single)


def cec2013_f22(problems: Iterable[Any], instances: Sequence[int], mode: str):
    mode = _mode_guard(mode)
    if mode == 'construct':
        return _construct_composition(problems, instances, 'cec2013_f22')
    if mode == 'repair':
        return _repair(instances)

    decs, single = _ensure_array(instances)
    centers = _load_shift_data()[:3, :]
    base_ids = [14, 14, 14]
    sigmas = [20, 20, 20]
    lambdas = [1, 1, 1]
    biases = [0, 100, 200]
    values = _composition_value(decs, centers, base_ids, lambdas, sigmas, biases) + 800.0
    return _finalise(values, single)


def cec2013_f23(problems: Iterable[Any], instances: Sequence[int], mode: str):
    mode = _mode_guard(mode)
    if mode == 'construct':
        return _construct_composition(problems, instances, 'cec2013_f23')
    if mode == 'repair':
        return _repair(instances)

    decs, single = _ensure_array(instances)
    centers = _load_shift_data()[:3, :]
    base_ids = [15, 15, 15]
    sigmas = [20, 20, 20]
    lambdas = [1, 1, 1]
    biases = [0, 100, 200]
    values = _composition_value(decs, centers, base_ids, lambdas, sigmas, biases) + 900.0
    return _finalise(values, single)


def cec2013_f24(problems: Iterable[Any], instances: Sequence[int], mode: str):
    mode = _mode_guard(mode)
    if mode == 'construct':
        return _construct_composition(problems, instances, 'cec2013_f24')
    if mode == 'repair':
        return _repair(instances)

    decs, single = _ensure_array(instances)
    centers = _load_shift_data()[:3, :]
    base_ids = [15, 12, 9]
    sigmas = [20, 20, 20]
    lambdas = [0.25, 1.0, 2.5]
    biases = [0, 100, 200]
    values = _composition_value(decs, centers, base_ids, lambdas, sigmas, biases) + 1000.0
    return _finalise(values, single)


def cec2013_f25(problems: Iterable[Any], instances: Sequence[int], mode: str):
    mode = _mode_guard(mode)
    if mode == 'construct':
        return _construct_composition(problems, instances, 'cec2013_f25')
    if mode == 'repair':
        return _repair(instances)

    decs, single = _ensure_array(instances)
    centers = _load_shift_data()[:3, :]
    base_ids = [15, 12, 9]
    sigmas = [10, 30, 50]
    lambdas = [0.25, 1.0, 2.5]
    biases = [0, 100, 200]
    values = _composition_value(decs, centers, base_ids, lambdas, sigmas, biases) + 1100.0
    return _finalise(values, single)


def cec2013_f26(problems: Iterable[Any], instances: Sequence[int], mode: str):
    mode = _mode_guard(mode)
    if mode == 'construct':
        return _construct_composition(problems, instances, 'cec2013_f26')
    if mode == 'repair':
        return _repair(instances)

    decs, single = _ensure_array(instances)
    centers = _load_shift_data()[:5, :]
    base_ids = [15, 12, 2, 9, 10]
    sigmas = [10, 10, 10, 10, 10]
    lambdas = [0.25, 1.0, 1e-7, 2.5, 10.0]
    biases = [0, 100, 200, 300, 400]
    values = _composition_value(decs, centers, base_ids, lambdas, sigmas, biases) + 1200.0
    return _finalise(values, single)


def cec2013_f27(problems: Iterable[Any], instances: Sequence[int], mode: str):
    mode = _mode_guard(mode)
    if mode == 'construct':
        return _construct_composition(problems, instances, 'cec2013_f27')
    if mode == 'repair':
        return _repair(instances)

    decs, single = _ensure_array(instances)
    centers = _load_shift_data()[:5, :]
    base_ids = [10, 12, 15, 9, 1]
    sigmas = [10, 10, 10, 20, 20]
    lambdas = [100, 10, 2.5, 25, 0.1]
    biases = [0, 100, 200, 300, 400]
    values = _composition_value(decs, centers, base_ids, lambdas, sigmas, biases) + 1300.0
    return _finalise(values, single)


def cec2013_f28(problems: Iterable[Any], instances: Sequence[int], mode: str):
    mode = _mode_guard(mode)
    if mode == 'construct':
        return _construct_composition(problems, instances, 'cec2013_f28')
    if mode == 'repair':
        return _repair(instances)

    decs, single = _ensure_array(instances)
    centers = _load_shift_data()[:5, :]
    base_ids = [19, 7, 15, 20, 1]
    sigmas = [10, 20, 30, 40, 50]
    lambdas = [2.5, 2.5e-3, 2.5, 5e-4, 0.1]
    biases = [0, 100, 200, 300, 400]
    values = _composition_value(decs, centers, base_ids, lambdas, sigmas, biases) + 1400.0
    return _finalise(values, single)

__all__ = [
    "cec2013_f1",
    "cec2013_f2",
    "cec2013_f3",
    "cec2013_f4",
    "cec2013_f5",
    "cec2013_f6",
    "cec2013_f7",
    "cec2013_f8",
    "cec2013_f9",
    "cec2013_f10",
    "cec2013_f11",
    "cec2013_f12",
    "cec2013_f13",
    "cec2013_f14",
    "cec2013_f15",
    "cec2013_f16",
    "cec2013_f17",
    "cec2013_f18",
    "cec2013_f19",
    "cec2013_f20",
    "cec2013_f21",
    "cec2013_f22",
    "cec2013_f23",
    "cec2013_f24",
    "cec2013_f25",
    "cec2013_f26",
    "cec2013_f27",
    "cec2013_f28",
]
