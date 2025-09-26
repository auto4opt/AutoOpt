"""Python translation of update_round_robin."""

from __future__ import annotations

import numpy as np

from ._utils import ensure_rng, extract_fits, flex_get, solution_as_list


def update_round_robin(*args):
    mode = args[-1]
    if mode == "execute":
        solution = args[0]
        problem = args[1] if len(args) > 1 else None
        aux = args[3] if len(args) > 3 else None
        rng = ensure_rng(aux)

        sol_list = solution_as_list(solution)
        fitness = extract_fits(solution)
        pop_size = len(sol_list)
        if pop_size == 0:
            return [], None

        n = int(flex_get(problem, "N", pop_size))
        n = min(n, pop_size)
        k = min(10, n - 1)
        if k <= 0:
            selected = np.argsort(fitness)[:n]
            return [sol_list[i] for i in selected], None

        win = np.zeros(pop_size, dtype=int)
        for i in range(pop_size):
            opponents = rng.choice(pop_size, size=k, replace=False)
            win[i] = int(np.sum(fitness[i] <= fitness[opponents]))

        rank = np.argsort(win)[::-1]
        chosen = rank[:n]
        return [sol_list[i] for i in chosen], None

    if mode == "parameter":
        return None, None

    if mode == "behavior":
        return ["", ""], None

    raise ValueError(f"Unsupported mode: {mode}")
