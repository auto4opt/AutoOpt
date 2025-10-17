"""Python translation of Utilities/General/Process.m."""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Iterable, List, Sequence, Tuple

import numpy as np

from ..design import Approximate, Design
from ..solve import input_algorithm, run_algorithm
from ..design._helpers import problem_list
from ..general.improve_rate import improve_rate
from ..space import space
from ...components import get_component
from .. import select as select_module


def _progress(app: Any, message: str) -> None:
    if app is None:
        return
    if hasattr(app, "TextArea"):
        try:
            app.TextArea.Value = message
        except Exception:  # pragma: no cover
            pass
    elif callable(app):
        app(message)


def _normalize_setting(setting: Any) -> SimpleNamespace:
    if isinstance(setting, SimpleNamespace):
        ns = setting
    elif isinstance(setting, dict):
        ns = SimpleNamespace(**setting)
    else:
        data = {
            name: getattr(setting, name)
            for name in dir(setting)
            if not name.startswith("__") and not callable(getattr(setting, name))
        }
        ns = SimpleNamespace(**data)

    for name in list(vars(ns)):
        if name[0].isupper():
            setattr(ns, name.lower(), getattr(ns, name))
    return ns


def _resolve_problem_callable(descriptor: Any) -> Callable[[Sequence[Any], Sequence[Any], str], Tuple[Any, Any, Any]]:
    if callable(descriptor):
        return descriptor
    if isinstance(descriptor, str):
        if ":" in descriptor:
            module_name, func_name = descriptor.split(":", 1)
            module = __import__(module_name, fromlist=[func_name])
            func = getattr(module, func_name)
            if callable(func):
                return func
        raise ValueError(
            "Problem descriptor must be a callable or 'module:function' string in Python translation."
        )
    raise TypeError("Unsupported problem descriptor format")


def _build_problem_struct(problem_descriptor: Any, instances: Sequence[Any], setting: SimpleNamespace):
    problems = []
    for _ in instances:
        problems.append(
            SimpleNamespace(
                name=problem_descriptor,
                setting="",
                N=int(getattr(setting, "ProbN", getattr(setting, "prob_n", 20))),
                Gmax=int(
                    math.ceil(
                        getattr(setting, "ProbFE", getattr(setting, "prob_fe", 5000))
                        / max(1, getattr(setting, "ProbN", getattr(setting, "prob_n", 20)))
                    )
                ),
            )
        )
    return problems


def _ensure_algorithm_list(obj: Iterable[Design]) -> List[Design]:
    return list(obj)


def _select(algs: Sequence[Design], problem: Any, data: Any, setting: Any, seeds: Sequence[int]) -> List[Design]:
    return select_module.select(algs, problem, data, setting, seeds)


def _mean_performance(algs: Sequence[Design]) -> np.ndarray:
    values = []
    for alg in algs:
        arr = alg.ave_perform_all()
        values.append(float(np.mean(arr)) if arr.size else float("inf"))
    return np.asarray(values, dtype=float)


def _update_cma_parameters(algs: Sequence[Design], problems: Sequence[Any], aux_list: List[Any]) -> List[Any]:
    for idx, aux in enumerate(aux_list):
        if isinstance(aux, dict) and "cma_Disturb" in aux:
            get_component("para_cma")(algs[idx], problems, aux, "algorithm")
    return aux_list


def _process_design(problem_descriptor: Any, instance_train: Sequence[Any], instance_test: Sequence[Any], setting: Any, app: Any) -> Tuple[List[Design], List[Design]]:
    setting_ns = _normalize_setting(setting)
    rng = np.random.default_rng(getattr(setting_ns, "random_seed", None))

    instance_train = list(instance_train)
    instance_test = list(instance_test)
    instances = instance_train + instance_test

    problems = _build_problem_struct(problem_descriptor, instances, setting_ns)
    construct_fn = _resolve_problem_callable(problem_descriptor)
    problems, data, _ = construct_fn(problems, instances, "construct")

    setting_ns = space(problems, setting_ns)
    alg_n = int(getattr(setting_ns, "AlgN", getattr(setting_ns, "alg_n", 1)))
    alg_fe = int(getattr(setting_ns, "AlgFE", getattr(setting_ns, "alg_fe", alg_n)))
    alg_runs = int(getattr(setting_ns, "AlgRuns", getattr(setting_ns, "alg_runs", 1)))
    alg_gmax = int(math.ceil(alg_fe / max(1, alg_n)))

    seed_train = rng.permutation(len(instance_train)).tolist()
    seed_test = (rng.permutation(len(instance_test)) + len(instance_train)).tolist()

    eval_mode = str(getattr(setting_ns, "Evaluate", getattr(setting_ns, "evaluate", "exact"))).lower()

    _progress(app, "Initializing...")

    algs: List[Design] = []
    surrogate = None
    if eval_mode in {"exact", "intensification"}:
        for _ in range(alg_n):
            alg = Design(problems, setting_ns)
            alg.evaluate(problems, data, setting_ns, seed_train)
            algs.append(alg)
    elif eval_mode == "racing":
        racing_k = int(getattr(setting_ns, "RacingK", getattr(setting_ns, "racingk", 1)))
        subset = seed_train[: racing_k or 1]
        for _ in range(alg_n):
            alg = Design(problems, setting_ns)
            alg.evaluate(problems, data, setting_ns, subset)
            algs.append(alg)
    elif eval_mode == "approximate":
        surrogate = Approximate(problems, data, setting_ns, seed_train)
        pool = surrogate.data
        indices = rng.choice(len(pool), size=min(alg_n, len(pool)), replace=False)
        algs = [pool[i] for i in indices]
    else:
        raise NotImplementedError(f"Evaluation mode '{eval_mode}' is not implemented in Python translation yet.")

    alg_trace: List[Design] = []
    G = 1

    while G <= alg_gmax:
        _progress(app, f"Designing... {100 * G / alg_gmax:.1f}%")
        improve = None
        inner_g = 1
        aux_list: List[Any] = [None] * len(algs)
        inner_gmax = int(math.ceil(alg_fe / max(alg_n, 1) / 10)) if alg_n == 1 else 1

        while (improve is None or improve[0] >= getattr(setting_ns, "IncRate", getattr(setting_ns, "inc_rate", 0.05))) and inner_g <= inner_gmax:
            new_algs: List[Design] = []
            for idx, alg in enumerate(algs):
                new_alg, aux_list[idx] = alg.get_new(problems, setting_ns, inner_g, aux_list[idx])
                new_algs.append(new_alg)

            if eval_mode == "exact":
                for new_alg in new_algs:
                    new_alg.evaluate(problems, data, setting_ns, seed_train)
                algs = _select(algs + new_algs, problems, data, setting_ns, seed_train)
            elif eval_mode == "approximate":
                for new_alg in new_algs:
                    new_alg.estimate(problems, setting_ns, seed_train, surrogate)
                algs = _select(algs + new_algs, problems, data, setting_ns, seed_train)
                if surrogate is not None and G in set(int(x) for x in np.asarray(surrogate.exact_g).reshape(-1)):
                    for new_alg in new_algs:
                        new_alg.evaluate(problems, data, setting_ns, seed_train)
                    surrogate.UpdateModel(new_algs, setting_ns)
            elif eval_mode == "racing":
                racing_k = int(getattr(setting_ns, "RacingK", getattr(setting_ns, "racingk", 1)))
                subset = seed_train[: racing_k or 1]
                for new_alg in new_algs:
                    new_alg.evaluate(problems, data, setting_ns, subset)
                algs = _select(algs + new_algs, problems, data, setting_ns, subset)
                remaining = seed_train[racing_k:]
                while len(algs) > alg_n and remaining:
                    algs = _select(algs, problems, data, setting_ns, [remaining[0]])
                    remaining = remaining[1:]
                seed_train = rng.permutation(len(instance_train)).tolist()
                if len(algs) > alg_n:
                    indices = rng.choice(len(algs), size=len(algs) - alg_n, replace=False)
                    algs = [alg for idx, alg in enumerate(algs) if idx not in indices]
            elif eval_mode == "intensification":
                while new_algs and seed_train:
                    subset = [seed_train[0]]
                    for alg in new_algs:
                        alg.evaluate(problems, data, setting_ns, subset)
                    new_algs = _select(algs + new_algs, problems, data, setting_ns, subset)
                    seed_train = seed_train[1:]
                seed_train = rng.permutation(len(instance_train)).tolist()
                for alg in new_algs:
                    for seed in seed_train:
                        if not np.any(alg.performance[seed, :]):
                            alg.evaluate(problems, data, setting_ns, [seed])
                if new_algs:
                    replace_idx = rng.choice(len(algs), size=min(len(new_algs), len(algs)), replace=False)
                    for dest, src in zip(replace_idx, new_algs):
                        algs[dest] = src
            else:
                raise NotImplementedError(f"Evaluation mode '{eval_mode}' is not implemented.")

            aux_list = _update_cma_parameters(algs, problems, aux_list)

            perf_values = _mean_performance(algs)
            wrapper = SimpleNamespace(avePerformAll=lambda: perf_values)
            improve = improve_rate(wrapper, improve, inner_g, "algorithm")
            inner_g += 1
            G += 1
            if G > alg_gmax:
                break

            best_idx = int(np.argmin(perf_values)) if len(perf_values) else 0
            alg_trace.append(algs[best_idx])

        if G > alg_gmax:
            break

    _progress(app, "Testing...")
    setting_ns.Evaluate = "exact"
    for alg in algs:
        alg.evaluate(problems, data, setting_ns, seed_test)
    algs = _select(algs, problems, data, setting_ns, seed_test)
    _progress(app, "Complete")
    return algs[:alg_n], alg_trace


def process(problem_descriptor: Any, *args, setting: Any, app: Any | None = None):
    """Main entry replicating MATLAB Process.m behaviour."""
    setting_ns = _normalize_setting(setting)
    mode = str(getattr(setting_ns, "Mode", getattr(setting_ns, "mode", "design"))).lower()

    if mode == "design":
        if len(args) < 2:
            raise ValueError("Design mode requires instance_train and instance_test arguments")
        instance_train, instance_test = args[:2]
        return _process_design(problem_descriptor, instance_train, instance_test, setting_ns, app)

    if mode == "solve":
        if len(args) < 1:
            raise ValueError("Solve mode requires instance list argument")
        instance_solve = args[0]
        problems = _build_problem_struct(problem_descriptor, instance_solve, setting_ns)
        construct_fn = _resolve_problem_callable(problem_descriptor)
        problems, data, _ = construct_fn(problems, instance_solve, "construct")
        alg, setting_ns = input_algorithm(setting_ns)
        _progress(app, "Solving...")
        best_solutions, all_solutions = run_algorithm(alg, problems, data, app, setting_ns)
        _progress(app, "Complete")
        return best_solutions, all_solutions

    raise ValueError("Mode must be 'design' or 'solve'.")
