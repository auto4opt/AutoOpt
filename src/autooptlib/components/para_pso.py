"""Python translation of para_pso."""
from __future__ import annotations

from typing import Any, Iterable

import numpy as np

from ._utils import flex_get
from .update_pairwise import update_pairwise

try:  # noqa: SIM105
    from ..utils.solve import SolutionSet  # type: ignore
except Exception:  # pylint: disable=broad-except
    class SolutionSet(list):  # type: ignore
        def __init__(self, items: Iterable[Any]):
            super().__init__(items)

        def decs(self) -> np.ndarray:
            if not self:
                return np.zeros((0, 0))
            return np.vstack([np.asarray(flex_get(sol, "dec"), dtype=float) for sol in self])

        def objs(self) -> np.ndarray:
            return np.asarray([flex_get(sol, "obj") for sol in self], dtype=float).reshape(-1, 1)

        def cons(self) -> np.ndarray:
            return np.asarray([flex_get(sol, "con") for sol in self], dtype=float).reshape(-1, 1)

        def fits(self) -> np.ndarray:
            return np.asarray([flex_get(sol, "fit") for sol in self], dtype=float).reshape(-1, 1)


def _ensure_solution_set(solution: Any) -> SolutionSet:
    if isinstance(solution, SolutionSet):
        return solution
    return SolutionSet(solution)


def para_pso(*args):
    solution = _ensure_solution_set(args[0])
    problem = args[1]
    aux = args[2] if len(args) > 2 else None

    if aux is None:
        aux = {}
    if "Pbest" not in aux:
        return aux

    pbest_set = aux["Pbest"]
    combined = list(pbest_set) + list(solution)
    selected, _ = update_pairwise(combined, "execute")
    aux["Pbest"] = SolutionSet(selected)

    objs = np.asarray([flex_get(sol, "obj") for sol in aux["Pbest"]], dtype=float)
    best_idx = int(np.argmin(objs))
    aux["Gbest"] = aux["Pbest"][best_idx]
    return aux
