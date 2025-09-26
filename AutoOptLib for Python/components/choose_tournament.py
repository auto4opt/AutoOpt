"""Python translation of choose_tournament."""

from __future__ import annotations

import numpy as np

from ._utils import flex_get, ensure_rng


def _extract_fitness(solution) -> np.ndarray:
    fitness = flex_get(solution, "fits")
    if fitness is None:
        raise ValueError("Solution object must provide fitness values")
    return np.asarray(fitness).reshape(-1)


def choose_tournament(*args):
    mode = args[-1]
    if mode == "execute":
        solution = args[0]
        problem = args[1] if len(args) > 1 else None
        aux = args[3] if len(args) > 3 else None
        rng = ensure_rng(aux)

        fitness = _extract_fitness(solution)
        pop_size = fitness.shape[0]
        k = 2
        n = int(flex_get(problem, "N", pop_size))

        match_index = rng.integers(0, pop_size, size=(n, k))
        match_fitness = fitness[match_index]
        row_min = match_fitness.min(axis=1)
        winners = []
        for i in range(n):
            candidates = np.where(match_fitness[i] == row_min[i])[0]
            winners.append(match_index[i, candidates[0]])
        return np.asarray(winners, dtype=int), None

    if mode == "parameter":
        return None, None

    if mode == "behavior":
        return ["", ""], None

    raise ValueError(f"Unsupported mode: {mode}")
