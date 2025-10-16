"""Selection utilities mirroring MATLAB Select.m (average mode)."""
from __future__ import annotations

from typing import Any, Iterable, List, Sequence

import numpy as np

from .design import Design
from .design._helpers import get_flex


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


def select(algs: Iterable[Design], problem: Any, data: Any, setting: Any, seed_instance: Sequence[int]) -> List[Design]:
    alg_list = _ensure_design_list(algs)
    if not alg_list:
        return []

    _require_performance(alg_list, problem, data, setting, seed_instance)
    all_perf = _collect_performance(alg_list, setting, seed_instance)
    averages = np.mean(all_perf, axis=0)

    compare = get_flex(setting, "compare", "average")
    evaluate_mode = get_flex(setting, "evaluate", "exact")
    alg_n = int(get_flex(setting, "alg_n", len(alg_list)))

    if compare != "average":
        raise NotImplementedError("Only compare='average' is supported in the current Python port")

    if evaluate_mode in {"exact", "approximate"}:
        order = np.argsort(averages)
        top = order[: min(alg_n, len(order))]
        return [alg_list[i] for i in top]

    if evaluate_mode == "intensification":
        old = alg_list[:alg_n]
        new = alg_list[alg_n:]
        if not new:
            return []
        old_avgs = averages[:alg_n]
        threshold = float(np.max(old_avgs)) if old_avgs.size else np.inf
        selected = [alg for alg, avg in zip(new, averages[alg_n:]) if avg < threshold]
        return selected

    if evaluate_mode == "racing":
        order = np.argsort(averages)
        top = order[: min(alg_n, len(order))]
        return [alg_list[i] for i in top]

    raise NotImplementedError(f"Unsupported evaluate mode: {evaluate_mode}")
