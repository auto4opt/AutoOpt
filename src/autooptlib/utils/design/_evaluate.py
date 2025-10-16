from __future__ import annotations

from typing import Any, List, Sequence

import numpy as np

from ...components import get_component
from ..general.improve_rate import improve_rate
from ..solve import SolutionSet, make_solutions, repair_sol
from ._helpers import Pathway, PathwayParam, get_flex, get_problem_type, problem_list


def _init_population(problem: Any, data: Any, setting: Any) -> SolutionSet:
    ptype = get_problem_type(problem) or "continuous"
    pop_n = int(get_flex(problem, "N", get_flex(setting, "prob_n", 10)))
    if ptype == "continuous":
        bound = np.asarray(get_flex(problem, "bound"), dtype=float)
        lower = bound[0]
        upper = bound[1]
        decs = lower + (upper - lower) * np.random.rand(pop_n, lower.shape[-1])
    elif ptype == "discrete":
        bound = np.asarray(get_flex(problem, "bound"), dtype=int)
        lower = bound[0]
        upper = bound[1]
        decs = np.vstack([
            np.random.randint(lower[j], upper[j] + 1, size=pop_n)
            for j in range(lower.shape[-1])
        ]).T
    elif ptype == "permutation":
        bound = np.asarray(get_flex(problem, "bound"), dtype=int)
        if bound.ndim == 2 and bound.shape[1] > 0:
            d = bound.shape[1]
        else:
            d = int(get_flex(problem, "dimension", 0) or get_flex(problem, "D", 0))
            if d == 0:
                raise ValueError("Permutation problem requires explicit dimension or bound")
        domain = np.arange(1, d + 1)
        decs = np.vstack([np.random.permutation(domain) for _ in range(pop_n)])
    else:
        raise NotImplementedError(f"Unsupported problem type: {ptype}")
    return make_solutions(decs, problem, data)


def _split_indices(indices: Sequence[int], alg_p: int, total: int) -> List[List[int]]:
    array = np.asarray(indices, dtype=int).flatten()
    if array.size == 0:
        indices = list(range(total))
    else:
        indices = array.tolist()
    if alg_p <= 1:
        return [indices]
    splits: List[List[int]] = []
    chunk = max(1, len(indices) // alg_p)
    start = 0
    for i in range(alg_p):
        if i == alg_p - 1:
            splits.append(indices[start:])
        else:
            splits.append(indices[start:start + chunk])
            start += chunk
    if not splits[-1]:
        splits[-1] = list(range(total))
    return splits


def _execute_path(
    path: Pathway,
    params: PathwayParam,
    subset: SolutionSet,
    aux_states: List[Any],
    problem: Any,
    data: Any,
    setting: Any,
    generation: int,
) -> tuple[List[Any], List[Any], int]:
    if not aux_states or len(aux_states) != len(path.search):
        aux_states = [None] * len(path.search)
    current_parent = subset
    evals = 0
    for step_idx, step in enumerate(path.search):
        improve = None
        inner_g = 1
        limit = int(step.termination[1]) if step.termination.size > 1 else 1
        threshold = float(step.termination[0]) if step.termination.size > 0 else -np.inf
        primary_param = params.search[step_idx].primary if params.search else None
        secondary_param = params.search[step_idx].secondary if params.search else None
        aux_state = aux_states[step_idx] or {}
        while (improve is None or improve[0] >= threshold) and inner_g <= limit:
            primary_fn = get_component(step.primary)
            new_dec, aux_state = primary_fn(current_parent, problem, primary_param, aux_state, generation, inner_g, data, "execute")
            if step.secondary:
                new_dec = repair_sol(np.asarray(new_dec), problem)
                secondary_fn = get_component(step.secondary)
                new_dec, aux_state = secondary_fn(new_dec, problem, secondary_param, aux_state, generation, inner_g, data, "execute")
            new = make_solutions(np.asarray(new_dec), problem, data)
            if step.primary == "search_cma":
                aux_state = get_component("para_cma")(new, problem, aux_state, "solution")
            elif step.primary == "search_pso":
                aux_state = get_component("para_pso")(new, problem, aux_state)
            current_parent = new
            evals += len(new)
            improve = improve_rate(new, improve, inner_g, "solution")
            aux_states[step_idx] = aux_state
            inner_g += 1
    return list(current_parent), aux_states, evals


def evaluate(self, problem: Any, data: Any, setting: Any, seed_instance: Sequence[int]):
    if not self.operator_pheno:
        return self

    pathways: List[Pathway] = self.operator_pheno[0]
    params: List[PathwayParam] = self.parameter_pheno[0]
    alg_p = max(1, len(pathways))

    problems = problem_list(problem)
    data_list = problem_list(data) if data is not None else [None] * len(problems)

    evaluate_mode = get_flex(setting, "evaluate", "exact")
    metric = get_flex(setting, "metric", "quality")
    runs = int(get_flex(setting, "alg_runs", 1))

    max_seed = max(seed_instance) if seed_instance else 0
    rows_needed = max_seed + 1

    if evaluate_mode == "approximate":
        if self.performance_approx.shape != (rows_needed, runs):
            new_matrix = np.zeros((rows_needed, runs))
            if self.performance_approx.size:
                r, c = self.performance_approx.shape
                new_matrix[: min(r, rows_needed), : min(c, runs)] = self.performance_approx[: min(r, rows_needed), : min(c, runs)]
            self.performance_approx = new_matrix
        target_matrix = self.performance_approx
    else:
        if self.performance.shape != (rows_needed, runs):
            new_matrix = np.zeros((rows_needed, runs))
            if self.performance.size:
                r, c = self.performance.shape
                new_matrix[: min(r, rows_needed), : min(c, runs)] = self.performance[: min(r, rows_needed), : min(c, runs)]
            self.performance = new_matrix
        target_matrix = self.performance

    for seed in seed_instance:
        problem_obj = problems[seed] if seed < len(problems) else problems[-1]
        data_obj = data_list[seed] if seed < len(data_list) else None
        for run in range(runs):
            solutions = _init_population(problem_obj, data_obj, setting)
            aux_cache: List[List[Any]] = [[None] * len(path.search) for path in pathways]
            archive_names = list(get_flex(setting, "archive", []))
            archives: List[Any] = []
            for name in archive_names:
                if name == "archive_statistic":
                    archives.append(np.empty((0, 2)))
                else:
                    archives.append([])

            best_history: List[float] = []
            evaluations = 0
            elapsed = 0.0
            generation = 1

            if metric in ("quality", "auc"):
                Tmax = np.inf
                threshold = float(get_flex(setting, "thres", -np.inf))
            elif metric == "runtimeFE":
                Tmax = int(np.ceil(get_flex(setting, "tmax", np.inf) / max(get_flex(setting, "prob_n", 1.0), 1)))
                threshold = float(get_flex(setting, "thres", -np.inf))
            elif metric == "runtimeSec":
                Tmax = float(get_flex(setting, "tmax", np.inf))
                threshold = float(get_flex(setting, "thres", -np.inf))
            else:
                Tmax = np.inf
                threshold = float(get_flex(setting, "thres", -np.inf))

            while True:
                current_best = float(np.min(solutions.fits())) if len(solutions) else np.inf
                best_history.append(current_best)

                if generation > int(get_flex(problem_obj, "Gmax", 10)):
                    break
                if metric == "quality" and current_best <= threshold:
                    break
                if metric == "runtimeFE" and evaluations >= Tmax:
                    break
                if metric == "runtimeSec" and elapsed >= Tmax:
                    break

                choose_fn = get_component(pathways[0].choose)
                choose_param = getattr(params[0], "choose", None)
                selected, _ = choose_fn(solutions, problem_obj, choose_param, None, generation, 1, data_obj, "execute")
                selected = np.asarray(selected, dtype=int).reshape(-1)
                index_groups = _split_indices(selected, alg_p, len(solutions))

                new_items: List[Any] = []
                for p_idx, path in enumerate(pathways):
                    subset_idx = index_groups[p_idx] if p_idx < len(index_groups) else index_groups[-1]
                    subset = SolutionSet([solutions[i] for i in subset_idx]) if subset_idx else solutions
                    produced, aux_state, evals = _execute_path(path, params[p_idx], subset, aux_cache[p_idx], problem_obj, data_obj, setting, generation)
                    aux_cache[p_idx] = aux_state
                    new_items.extend(produced)
                    evaluations += evals

                update_fn = get_component(pathways[0].update)
                update_param = getattr(params[0], "update", None)
                updated, _ = update_fn(list(solutions) + new_items, problem_obj, update_param, None, generation, 1, data_obj, "execute")
                solutions = SolutionSet(updated)

                for j, name in enumerate(archive_names):
                    archive_fn = get_component(name)
                    if name == "archive_best":
                        archives[j], _ = archive_fn(solutions, archives[j], "execute")
                    elif name == "archive_statistic":
                        archives[j], _ = archive_fn(solutions, archives[j], "execute")
                    else:
                        archives[j], _ = archive_fn(solutions, archives[j], problem_obj, "execute")

                if metric == "runtimeFE":
                    elapsed = evaluations
                elif metric == "runtimeSec":
                    elapsed += 1.0

                generation += 1

            if metric == "quality":
                value = float(best_history[-1]) if best_history else np.inf
            elif metric == "auc":
                history = np.asarray(best_history, dtype=float)
                if history.size == 0:
                    value = np.inf
                else:
                    normalized = history - history.min()
                    value = float(np.trapz(normalized, dx=1.0))
            elif metric == "runtimeFE":
                value = evaluations
            elif metric == "runtimeSec":
                value = elapsed
            else:
                value = float(best_history[-1]) if best_history else np.inf

            target_matrix[seed, run] = value

    if evaluate_mode == "approximate":
        self.performance_approx = target_matrix
    else:
        self.performance = target_matrix

    return self
