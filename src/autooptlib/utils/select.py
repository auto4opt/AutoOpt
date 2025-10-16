"""Selection utilities mirroring MATLAB Select.m (average mode)."""
from __future__ import annotations

from typing import Any, Iterable, List, Sequence, Tuple

import numpy as np

from .design import Design
from .design._helpers import get_flex
from .general.select_stats import friedman_statistic


def _ensure_design_list(algs: Iterable[Design]) -> List[Design]:
    if isinstance(algs, list):
        return algs
    return list(algs)


def _require_performance(algs: List[Design], problem: Any, data: Any, setting: Any, seed_instance: Sequence[int]) -> None:
    evaluate_mode = get_flex(setting, "evaluate", "exact")
    if evaluate_mode != "racing":
        return
    for alg in algs:
        performance = getattr(alg, "performance", None)
        if performance is None:
            continue
        for seed in seed_instance:
            row = np.asarray(performance)[seed]
            if np.sum(row) == 0:
                alg.evaluate(problem, data, setting, [seed])


def _collect_performance(algs: List[Design], setting: Any, seed_instance: Sequence[int]) -> np.ndarray:
    runs = int(get_flex(setting, "alg_runs", 1))
    matrix = np.zeros((len(seed_instance) * runs, len(algs)))
    for idx, alg in enumerate(algs):
        perf = alg.get_performance(setting, seed_instance)
        matrix[:, idx] = np.asarray(perf).reshape(-1)
    return matrix


def _compare_friedman(matrix: np.ndarray, alpha: float) -> np.ndarray:
    ranks, p_values = friedman_statistic(matrix)
    wins = np.zeros(matrix.shape[1], dtype=int)
    for (i, j), p_val in np.ndenumerate(p_values):
        if i >= j or np.isnan(p_val):
            continue
        diff = ranks[i] - ranks[j]
        if diff < 0 and p_val < alpha:
            wins[i] += 1
        elif diff > 0 and p_val < alpha:
            wins[j] += 1
    return wins


def select(algs: Iterable[Design], problem: Any, data: Any, setting: Any, seed_instance: Sequence[int]) -> List[Design]:
    alg_list = _ensure_design_list(algs)
    if not alg_list:
        return []

    _require_performance(alg_list, problem, data, setting, seed_instance)
    all_perf = _collect_performance(alg_list, setting, seed_instance)

    compare = get_flex(setting, "compare", "average")
    evaluate_mode = get_flex(setting, "evaluate", "exact")
    alg_n = int(get_flex(setting, "alg_n", len(alg_list)))
    alpha = float(get_flex(setting, "alpha", 0.05))

    if compare == "average":
        averages = np.mean(all_perf, axis=0)
        order = np.argsort(averages)
        if evaluate_mode in {"exact", "approximate", "racing"}:
            top = order[: min(alg_n, len(order))]
            return [alg_list[i] for i in top]
        if evaluate_mode == "intensification":
            old_avgs = averages[:alg_n]
            threshold = float(np.max(old_avgs)) if old_avgs.size else np.inf
            return [alg for alg, avg in zip(alg_list[alg_n:], averages[alg_n:]) if avg < threshold]
        raise NotImplementedError(f"Unsupported evaluate mode: {evaluate_mode}")

    if compare == "statistic":
        wins = _compare_friedman(all_perf, alpha=alpha)
        order = np.argsort(-wins)
        if evaluate_mode in {"exact", "approximate", "racing"}:
            top = order[: min(alg_n, len(order))]
            return [alg_list[i] for i in top]
        if evaluate_mode == "intensification":
            return [alg_list[i] for i in order if i >= alg_n and wins[i] > 0]
        raise NotImplementedError(f"Unsupported evaluate mode: {evaluate_mode}")

    raise NotImplementedError(f"Unsupported compare mode: {compare}")
