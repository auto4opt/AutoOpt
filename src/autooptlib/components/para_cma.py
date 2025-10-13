"""Python translation of para_cma."""
from __future__ import annotations

from typing import Any

import numpy as np

from ._utils import flex_get


def _extract_fitness(obj: Any, mode: str) -> np.ndarray:
    if mode == "solution":
        fits = flex_get(obj, "fits")
        if callable(fits):
            return np.asarray(fits()).reshape(-1)
        if fits is not None:
            return np.asarray(fits).reshape(-1)
        raise ValueError("Solution must expose fits() for para_cma")
    if mode == "algorithm":
        vals = flex_get(obj, "avePerformAll")
        if callable(vals):
            return np.asarray(vals()).reshape(-1)
        if vals is not None:
            return np.asarray(vals).reshape(-1)
        raise ValueError("Algorithm must expose avePerformAll for para_cma")
    raise ValueError(f"Unknown type {mode}")


def para_cma(*args):
    solution = args[0]
    problem = args[1]
    aux = args[2] if len(args) > 2 else {}
    mode = args[3] if len(args) > 3 else "solution"

    if aux is None:
        aux = {}
    disturbance = np.asarray(aux.get("cma_Disturb"))
    if disturbance.size == 0:
        return aux

    half_n = int(aux.get("cma_halfN"))
    w = np.asarray(aux.get("cma_w"))
    better_n = float(aux.get("cma_betterN"))
    mean = np.asarray(aux.get("cma_mean"), dtype=float)
    sigma = np.asarray(aux.get("cma_sigma"), dtype=float)
    csigma = float(aux.get("cma_csigma"))
    dsigma = float(aux.get("cma_dsigma"))
    chi_n = float(aux.get("cma_chiN"))
    cc = float(aux.get("cma_cc"))
    ccov = float(aux.get("cma_ccov"))
    cmu = float(aux.get("cma_cmu"))
    hth = float(aux.get("cma_hth"))
    ps = np.asarray(aux.get("cma_ps"), dtype=float)
    pc = np.asarray(aux.get("cma_pc"), dtype=float)
    C = np.asarray(aux.get("cma_C"), dtype=float)

    fitness = _extract_fitness(solution, mode)
    rank = np.argsort(fitness)
    disturbance = disturbance[rank]
    disturbance_w = w @ disturbance[:half_n]
    mean = mean + sigma * disturbance_w

    try:
        chol_c = np.linalg.cholesky(C)
        inv_chol_T = np.linalg.solve(chol_c.T, disturbance_w)
    except np.linalg.LinAlgError:
        chol_c = np.linalg.cholesky(C + 1e-12 * np.eye(C.shape[0]))
        inv_chol_T = np.linalg.solve(chol_c.T, disturbance_w)

    ps = (1 - csigma) * ps + np.sqrt(csigma * (2 - csigma) * better_n) * inv_chol_T
    gmax = int(flex_get(problem[0] if isinstance(problem, (list, tuple)) else problem, "Gmax", 1))
    denom = np.sqrt(1 - (1 - csigma) ** (2 * (gmax + 1)))
    hs = float(np.linalg.norm(ps) / denom < hth)
    pc = (1 - cc) * pc + hs * np.sqrt(cc * (2 - cc) * better_n) * disturbance_w
    delta = (1 - hs) * cc * (2 - cc)
    C = (1 - ccov - cmu) * C + ccov * (np.outer(pc, pc) + delta * C)
    for i in range(half_n):
        C += cmu * w[i] * np.outer(disturbance[i], disturbance[i])
    C = 0.5 * (C + C.T)
    eigvals, eigvecs = np.linalg.eigh(C)
    eigvals = np.clip(eigvals, 0, None)
    C = eigvecs @ np.diag(eigvals) @ eigvecs.T

    sigma = sigma * np.exp(csigma / dsigma * (np.linalg.norm(ps) / chi_n - 1)) ** 0.3

    aux.update(
        {
            "cma_mean": mean,
            "cma_ps": ps,
            "cma_pc": pc,
            "cma_C": C,
            "cma_sigma": sigma,
        }
    )
    return aux