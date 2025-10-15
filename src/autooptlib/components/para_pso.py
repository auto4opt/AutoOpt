"""Python translation of para_pso."""
from __future__ import annotations

from typing import Any

import numpy as np

from ..utils.solve import SolutionSet
from ._utils import flex_get
from .update_pairwise import update_pairwise


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