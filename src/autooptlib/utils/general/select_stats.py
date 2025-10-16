from __future__ import annotations

import numpy as np
from typing import Tuple


def friedman_statistic(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("Input to friedman_statistic must be 2-D")
    n, k = matrix.shape
    ranks = np.argsort(np.argsort(matrix, axis=1), axis=1) + 1
    R = np.sum(ranks, axis=0)
    chi2 = (12 / (n * k * (k + 1))) * np.sum(R**2) - 3 * n * (k + 1)
    q_stat = chi2 * (k - 1) / (k * (k - 1))
    # Placeholder: no real multiple comparison run
    p_values = np.full((k, k), np.nan)
    return R / n, p_values
