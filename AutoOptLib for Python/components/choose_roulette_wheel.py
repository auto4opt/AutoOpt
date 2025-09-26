"""Python translation of choose_roulette_wheel."""

from __future__ import annotations

import numpy as np

from ._utils import ensure_rng, flex_get


def _extract_fitness(solution) -> np.ndarray:
    fitness = flex_get(solution, "fits")
    if fitness is None:
        raise ValueError("Solution object must provide fitness values")
    arr = np.asarray(fitness, dtype=float)
    return arr.reshape(-1)


def choose_roulette_wheel(*args):
    mode = args[-1]
    if mode == "execute":
        solution = args[0]
        problem = args[1] if len(args) > 1 else None
        aux = args[3] if len(args) > 3 else None
        rng = ensure_rng(aux)

        fitness = _extract_fitness(solution)
        n = int(flex_get(problem, "N", fitness.shape[0]))

        fitness = fitness - min(fitness.min(), 0) + 1e-6
        inv_fit = 1.0 / fitness
        cdf = np.cumsum(inv_fit)
        cdf /= cdf[-1]

        random_values = rng.random(n)
        indices = np.searchsorted(cdf, random_values, side="left")
        return indices.astype(int), None

    if mode == "parameter":
        return None, None

    if mode == "behavior":
        return ["", ""], None

    raise ValueError(f"Unsupported mode: {mode}")
