from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, List, Optional, Sequence, Tuple, TYPE_CHECKING
import pickle

import numpy as np
import time

from ..design._helpers import ensure_rng, get_flex, get_problem_type
from ..design._helpers import Pathway, PathwayParam, SearchStep, SearchParam
from ..general.improve_rate import improve_rate
from ...components import get_component

if TYPE_CHECKING:
    from ..design import Design


@dataclass
class Solution:
    dec: np.ndarray
    obj: float
    con: float
    fit: float


class SolutionSet(Sequence[Solution]):
    def __init__(self, items: Iterable[Solution]):
        self._items: List[Solution] = list(items)

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, idx):
        return self._items[idx]

    # Aggregations compatible with components
    def decs(self) -> np.ndarray:
        if not self._items:
            return np.zeros((0, 0))
        return np.vstack([s.dec for s in self._items])

    def objs(self) -> np.ndarray:
        return np.asarray([s.obj for s in self._items]).reshape(-1, 1)

    def cons(self) -> np.ndarray:
        return np.asarray([s.con for s in self._items]).reshape(-1, 1)

    def fits(self) -> np.ndarray:
        return np.asarray([s.fit for s in self._items]).reshape(-1, 1)


def repair_sol(dec: np.ndarray, problem: Any) -> np.ndarray:
    ptype = get_problem_type(problem)
    bound = getattr(problem, "bound", None)
    if bound is None:
        return dec
    bound = np.asarray(bound)
    if ptype == "continuous":
        lower = bound[0]
        upper = bound[1]
        return np.minimum(np.maximum(dec, lower), upper)
    if ptype == "discrete":
        lower = bound[0]
        upper = bound[1]
        repaired = np.clip(np.rint(dec), lower, upper).astype(int)
        setting = get_flex(problem, "setting", "")
        if isinstance(setting, str) and "dec_diff" in setting:
            n, d = repaired.shape
            for i in range(n):
                seen = set()
                duplicates = []
                for j, val in enumerate(repaired[i]):
                    if val in seen:
                        duplicates.append(j)
                    else:
                        seen.add(val)
                if not duplicates:
                    continue
                for idx in duplicates:
                    choices = [v for v in range(int(lower[idx]), int(upper[idx]) + 1) if v not in repaired[i]]
                    if choices:
                        repaired[i, idx] = int(np.random.choice(choices))
            return repaired
        return repaired
    if ptype == "permutation":
        repaired = np.rint(dec).astype(int)
        if repaired.ndim == 1:
            rows = [repaired.copy()]
        else:
            rows = [row.copy() for row in repaired]
        d = len(rows[0])
        domain = list(range(1, d + 1))
        for row in rows:
            missing = [v for v in domain if v not in row]
            counts = {}
            duplicates = []
            for j, val in enumerate(row):
                counts[val] = counts.get(val, 0) + 1
                if counts[val] > 1:
                    duplicates.append(j)
            for j in duplicates:
                if missing:
                    row[j] = missing.pop(0)
        if repaired.ndim == 1:
            return rows[0]
        return np.vstack(rows)
    return dec


def _call_evaluator(problem: Any, data: Any, dec: np.ndarray):
    # Preferred: problem.evaluate(data, dec)
    eval_fn = getattr(problem, 'evaluate', None)
    if callable(eval_fn):
        out = eval_fn(data, dec)
    else:
        # Fallback: Problem.name referencing a callable in current scope
        name = getattr(problem, 'name', None)
        if callable(name):
            out = name(data, dec)
        else:
            raise NotImplementedError('Problem must provide an evaluate(data, dec) callable')
    # Normalize outputs: obj, con, acc(optional)
    if isinstance(out, tuple):
        if len(out) == 3:
            obj, con, acc = out
        elif len(out) == 2:
            obj, con = out
            acc = None
        else:
            obj, con, acc = out[0], 0.0, None
    else:
        obj, con, acc = out, 0.0, None
    return float(np.asarray(obj).reshape(1)[0]), float(np.sum(np.maximum(0.0, np.asarray(con)))), acc


def make_solutions(decs: np.ndarray, problem: Any, data: Any) -> SolutionSet:
    items: List[Solution] = []
    for i in range(decs.shape[0]):
        d = np.array(decs[i], dtype=float)
        d = repair_sol(d, problem)
        obj, con, _ = _call_evaluator(problem, data, d)
        feasible = con <= 0.0
        fit = obj if feasible else (con + 1e8)
        items.append(Solution(dec=d, obj=obj, con=con, fit=fit))
    return SolutionSet(items)


def _init_population(rng: np.random.Generator, problem: Any, data: Any, setting: Any) -> SolutionSet:
    ptype = get_problem_type(problem) or "continuous"
    pop_n = int(get_flex(problem, "N", get_flex(setting, "prob_n", 10)))
    if ptype == "continuous":
        bound = np.asarray(get_flex(problem, "bound"), dtype=float)
        lower = bound[0]
        upper = bound[1]
        decs = lower + (upper - lower) * rng.random((pop_n, lower.shape[-1]))
    elif ptype == "discrete":
        bound = np.asarray(get_flex(problem, "bound"), dtype=int)
        lower = bound[0]
        upper = bound[1]
        decs = np.vstack(
            [rng.integers(lower[j], upper[j] + 1, size=pop_n) for j in range(lower.shape[-1])]
        ).T
    elif ptype == "permutation":
        bound = np.asarray(get_flex(problem, "bound"), dtype=int)
        if bound.ndim == 2 and bound.shape[1] > 0:
            d = bound.shape[1]
        else:
            d = int(get_flex(problem, "dimension", 0) or get_flex(problem, "D", 0))
            if d == 0:
                raise ValueError("Permutation problem requires explicit dimension or bound")
        domain = np.arange(1, d + 1)
        decs = np.vstack([rng.permutation(domain) for _ in range(pop_n)])
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
    path: Any,
    params: Any,
    subset: SolutionSet,
    aux_states: List[Any],
    problem: Any,
    data: Any,
    setting: Any,
    generation: int,
) -> Tuple[List[Any], List[Any], int]:
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
            new_dec, aux_state = primary_fn(
                current_parent, problem, primary_param, aux_state, generation, inner_g, data, "execute"
            )
            if step.secondary:
                new_dec = repair_sol(np.asarray(new_dec), problem)
                secondary_fn = get_component(step.secondary)
                new_dec, aux_state = secondary_fn(
                    new_dec, problem, secondary_param, aux_state, generation, inner_g, data, "execute"
                )
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


def run_design(pathways: List[Any], params: List[Any], problem: Any, data: Any, setting: Any) -> dict:
    rng = ensure_rng(setting)
    solutions = _init_population(rng, problem, data, setting)
    aux_cache: List[List[Any]] = [[None] * len(path.search) for path in pathways]
    archive_names = list(get_flex(setting, "archive", []))
    archives: List[Any] = []
    for name in archive_names:
        if name == "archive_statistic":
            archives.append(np.empty((0, 2)))
        else:
            archives.append([])

    metric = get_flex(setting, "metric", "quality")
    prob_n = float(get_flex(setting, "prob_n", 1.0))
    tmax = get_flex(setting, "tmax", np.inf)
    if tmax is None:
        tmax = np.inf
    thres_val = get_flex(setting, "thres", -np.inf)
    thres = float(thres_val if thres_val is not None else -np.inf)
    gmax = int(get_flex(problem, "Gmax", 10))

    if metric == "quality":
        Tmax = np.inf
        threshold = thres
    elif metric == "auc":
        Tmax = np.inf
        threshold = -np.inf
    elif metric == "runtimeFE":
        Tmax = int(np.ceil(tmax / max(prob_n, 1)))
        threshold = thres
    elif metric == "runtimeSec":
        Tmax = float(tmax)
        threshold = thres
    else:
        Tmax = np.inf
        threshold = thres

    history: List[Optional[Solution]] = []
    fit_history: List[float] = []
    evaluations = 0
    elapsed_seconds = 0.0
    generation = 1

    while True:
        if len(solutions):
            best_solution = min(list(solutions), key=lambda s: s.fit)
            history.append(best_solution)
            fit_history.append(float(best_solution.fit))
        else:
            history.append(None)
            fit_history.append(np.inf)

        if generation > gmax:
            break
        if metric == "quality" and fit_history[-1] <= threshold:
            break
        if metric == "runtimeFE" and evaluations >= Tmax:
            break
        if metric == "runtimeSec" and elapsed_seconds >= Tmax:
            break

        iteration_start = time.perf_counter()

        choose_fn = get_component(pathways[0].choose)
        choose_param = getattr(params[0], "choose", None)
        selected, _ = choose_fn(solutions, problem, choose_param, None, generation, 1, data, "execute")
        selected = np.asarray(selected, dtype=int).reshape(-1)
        index_groups = _split_indices(selected, len(pathways), len(solutions))

        new_items: List[Any] = []
        for p_idx, path in enumerate(pathways):
            subset_idx = index_groups[p_idx] if p_idx < len(index_groups) else index_groups[-1]
            subset = SolutionSet([solutions[i] for i in subset_idx]) if subset_idx else solutions
            produced, aux_state, evals = _execute_path(
                path, params[p_idx], subset, aux_cache[p_idx], problem, data, setting, generation
            )
            aux_cache[p_idx] = aux_state
            new_items.extend(produced)
            evaluations += evals

        update_fn = get_component(pathways[0].update)
        update_param = getattr(params[0], "update", None)
        updated, _ = update_fn(list(solutions) + new_items, problem, update_param, None, generation, 1, data, "execute")
        solutions = SolutionSet(updated)

        for j, name in enumerate(archive_names):
            archive_fn = get_component(name)
            if name == "archive_best":
                archives[j], _ = archive_fn(solutions, archives[j], "execute")
            elif name == "archive_statistic":
                archives[j], _ = archive_fn(solutions, archives[j], "execute")
            else:
                archives[j], _ = archive_fn(solutions, archives[j], problem, "execute")

        if metric == "runtimeFE":
            pass
        elif metric == "runtimeSec":
            elapsed_seconds += time.perf_counter() - iteration_start

        generation += 1

    if history and history[-1] is None:
        final_solution = None
    else:
        final_solution = history[-1]

    return {
        "solutions": solutions,
        "history": history,
        "fit_history": fit_history,
        "evaluations": evaluations,
        "elapsed": elapsed_seconds,
        "archives": archives,
        "best_solution": final_solution,
    }


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _placeholder_solution() -> Solution:
    return Solution(dec=np.zeros((1, 0)), obj=float("inf"), con=float("inf"), fit=float("inf"))


def _ensure_solution_instance(obj: Any) -> Solution:
    if isinstance(obj, Solution):
        return obj
    if obj is None:
        return _placeholder_solution()
    if hasattr(obj, "dec") and hasattr(obj, "fit"):
        return Solution(
            dec=np.asarray(getattr(obj, "dec")),
            obj=float(getattr(obj, "obj", np.inf)),
            con=float(getattr(obj, "con", np.inf)),
            fit=float(getattr(obj, "fit", np.inf)),
        )
    return _placeholder_solution()


# ---------------------------------------------------------------------------
# Solve-mode helpers (InputAlg / RunAlg translation)
# ---------------------------------------------------------------------------

def _normalize_setting(setting: Any) -> Any:
    if isinstance(setting, dict):
        return type("Setting", (), setting)()
    for attr in dir(setting):
        if attr[0].isupper():
            setattr(setting, attr.lower(), getattr(setting, attr))
    return setting


def _to_array(value: Any) -> Optional[np.ndarray]:
    if value is None or value == []:
        return None
    arr = np.asarray(value, dtype=float)
    return arr.reshape(-1) if arr.ndim > 0 else arr


def _make_search_step(primary: str, secondary: Optional[str], termination: Sequence[float]) -> SearchStep:
    term = np.asarray(termination, dtype=float).reshape(-1)
    return SearchStep(primary=primary, secondary=secondary if secondary else None, termination=term)


def _make_search_param(primary: Any, secondary: Any) -> SearchParam:
    return SearchParam(primary=_to_array(primary), secondary=_to_array(secondary))


def _build_default_algorithm(name: str, setting: Any) -> "Design":
    name = name.strip().lower()
    tmpl = {
        "continuous genetic algorithm": [
            (
                "choose_tournament",
                [
                    ("cross_sim_binary", "search_mu_polynomial", (-np.inf, 1), [20.0], [0.2, 20.0]),
                ],
                "update_round_robin",
                [],
            ),
        ],
        "evolutionary programming": [
            (
                "choose_tournament",
                [
                    ("search_mu_gaussian", None, (-np.inf, 1), None, None),
                ],
                "update_round_robin",
                [],
            ),
        ],
        "fast evolutionary programming": [
            (
                "choose_tournament",
                [
                    ("search_mu_cauchy", None, (-np.inf, 1), None, None),
                ],
                "update_round_robin",
                [],
            ),
        ],
        "cma-es": [
            (
                "choose_traverse",
                [
                    ("search_cma", None, (-np.inf, 1), None, None),
                ],
                "update_greedy",
                [],
            ),
        ],
        "estimation of distribution": [
            (
                "choose_traverse",
                [
                    ("search_eda", None, (-np.inf, 1), None, None),
                ],
                "update_greedy",
                [],
            ),
        ],
        "particle swarm optimization": [
            (
                "choose_traverse",
                [
                    ("search_pso", None, (-np.inf, 1), [0.9], None),
                ],
                "update_pairwise",
                [],
            ),
        ],
        "differential evolution": [
            (
                "choose_traverse",
                [
                    ("search_de_current", None, (-np.inf, 1), [0.9, 0.1], None),
                ],
                "update_pairwise",
                [],
            ),
        ],
        "continuous random search": [
            (
                "choose_traverse",
                [
                    ("reinit_continuous", None, (-np.inf, 1), None, None),
                ],
                "update_greedy",
                [],
            ),
        ],
        "discrete genetic algorithm": [
            (
                "choose_tournament",
                [
                    ("cross_point_uniform", "search_reset_rand", (-np.inf, 1), [0.2], [0.2]),
                ],
                "update_round_robin",
                [],
            ),
        ],
        "discrete iterative local search": [
            (
                "choose_traverse",
                [
                    ("search_reset_one", None, (0.05, 10), None, None),
                    ("reinit_discrete", None, (-np.inf, 1), None, None),
                ],
                "update_greedy",
                [],
            ),
        ],
        "discrete simulated annealing": [
            (
                "choose_traverse",
                [
                    ("search_reset_one", None, (-np.inf, 1), None, None),
                ],
                "update_simulated_annealing",
                [],
                np.array([0.1]),
            ),
        ],
        "discrete random search": [
            (
                "choose_traverse",
                [
                    ("reinit_discrete", None, (-np.inf, 1), None, None),
                ],
                "update_greedy",
                [],
            ),
        ],
        "permutation genetic algorithm": [
            (
                "choose_tournament",
                [
                    ("cross_order_two", "search_swap", (-np.inf, 1), None, None),
                ],
                "update_round_robin",
                [],
            ),
        ],
        "permutation iterative local search": [
            (
                "choose_traverse",
                [
                    ("search_insert", None, (0.05, 10), None, None),
                    ("reinit_permutation", None, (-np.inf, 1), None, None),
                ],
                "update_greedy",
                [],
            ),
        ],
        "permutation simulated annealing": [
            (
                "choose_traverse",
                [
                    ("search_insert", None, (-np.inf, 1), None, None),
                ],
                "update_simulated_annealing",
                [],
                np.array([0.1]),
            ),
        ],
        "permutation variable neighborhood search": [
            (
                "choose_traverse",
                [
                    ("search_swap", None, (0.05, 10), None, None),
                    ("search_scramble", None, (0.05, 10), None, None),
                    ("search_insert", None, (0.05, 10), None, None),
                ],
                "update_greedy",
                [],
            ),
        ],
        "permutation random search": [
            (
                "choose_traverse",
                [
                    ("reinit_permutation", None, (-np.inf, 1), None, None),
                ],
                "update_greedy",
                [],
            ),
        ],
    }

    if name not in tmpl:
        raise NotImplementedError(f"Preset algorithm {name!r} is not available in the Python translation.")

    pathways_cfg = tmpl[name]
    pathways: List[Pathway] = []
    params: List[PathwayParam] = []
    for item in pathways_cfg:
        choose = item[0]
        search_entries = item[1]
        update = item[2]
        archive = item[3] if len(item) > 3 else []
        update_param = item[4] if len(item) > 4 else None

        search_steps = []
        search_params = []
        for entry in search_entries:
            primary, secondary, termination, primary_param, secondary_param = entry
            search_steps.append(_make_search_step(primary, secondary, termination))
            search_params.append(_make_search_param(primary_param, secondary_param))

        pathways.append(Pathway(choose=choose, search=search_steps, update=update, archive=archive))
        params.append(
            PathwayParam(
                choose=None,
                search=search_params,
                update=_to_array(update_param),
            )
        )

    from ..design import Design

    design = Design()
    design.construct([pathways], [params])
    return design


def input_algorithm(setting: Any) -> Tuple["Design", Any]:
    setting = _normalize_setting(setting)
    alg_file = getattr(setting, "AlgFile", getattr(setting, "alg_file", ""))
    if alg_file and str(alg_file).lower() != "none":
        path = Path(alg_file)
        if not path.exists():
            raise FileNotFoundError(f"Algorithm file {alg_file} not found.")
        with path.open("rb") as handle:
            data = pickle.load(handle)
        if isinstance(data, dict) and "algs" in data:
            alg = data["algs"][0]
        else:
            alg = data[0] if isinstance(data, (list, tuple)) else data
        from ..design import Design
        if not isinstance(alg, Design):
            raise TypeError("Loaded algorithm is not a Design instance.")
        setting.AlgP = len(getattr(alg, "operator_pheno", []) or [])
        return alg, setting

    alg_name = getattr(setting, "AlgName", getattr(setting, "alg_name", "")).strip()
    if not alg_name:
        raise ValueError("Please specify AlgFile or AlgName in solve mode.")
    design = _build_default_algorithm(alg_name.lower(), setting)
    setting.AlgP = len(design.operator_pheno or [])
    return design, setting


def _extract_mode(problem: Any) -> str:
    behavior = get_flex(problem, "type", ["continuous", "static"])
    if isinstance(behavior, (list, tuple)) and len(behavior) > 1:
        return str(behavior[1])
    return "static"


def run_algorithm(alg: 'Design', problems: Sequence[Any], data: Sequence[Any], app: Any, setting: Any):
    pathways = alg.operator_pheno[0]
    params = alg.parameter_pheno[0]
    alg_runs = int(get_flex(setting, "AlgRuns", 1))
    best_solutions: List[List[Solution]] = []
    all_solutions: List[List[List[Solution]]] = []

    for idx, problem in enumerate(problems):
        data_obj = data[idx] if idx < len(data) else None
        mode = _extract_mode(problem)
        run_histories: List[List[Solution]] = []
        run_best: List[Solution] = []
        fitness_sequence: List[float] = []

        for run in range(alg_runs):
            if mode == "static":
                result = run_design(pathways, params, problem, data_obj, setting)
                history = [sol for sol in result["history"] if isinstance(sol, Solution)]
                best = result["best_solution"]
                if best is None and history:
                    best = history[-1]
                run_histories.append(history)
                run_best.append(_ensure_solution_instance(best))
            elif mode == "sequential":
                current_problem = problem
                current_data = data_obj
                sequence_solutions: List[Solution] = []
                cumulative = 0.0
                while getattr(current_data, "continue", False):
                    result = run_design(pathways, params, current_problem, current_data, setting)
                    best = result["best_solution"]
                    if best is None and result["history"]:
                        history = [sol for sol in result["history"] if isinstance(sol, Solution)]
                        best = history[-1] if history else None
                    if best is None:
                        break
                    sequence_solutions.append(best)
                    cumulative += float(best.fit)
                    next_fn = getattr(current_problem, "advance_sequence", None)
                    if callable(next_fn):
                        current_problem, current_data = next_fn(best, current_data)
                    else:
                        name = getattr(current_problem, "name", None)
                        if callable(name):
                            current_problem, current_data, _ = name(current_problem, current_data, best, "sequence")
                        else:
                            break
                    if current_problem is problem and current_data is data_obj:
                        break
                run_histories.append(sequence_solutions)
                run_best.append(sequence_solutions[-1] if sequence_solutions else _placeholder_solution())
                fitness_sequence.append(cumulative)
            else:
                raise NotImplementedError(f"Unsupported problem mode: {mode}")

            if app is not None and hasattr(app, "TextArea"):
                pct = (idx * alg_runs + run + 1) / (len(problems) * alg_runs)
                app.TextArea.Value = f"Solving... {pct * 100:.1f}%"

        if mode == "static":
            best_idx = int(np.argmin([b.fit for b in run_best]))
            all_solutions.append(run_histories[best_idx])
        else:
            best_idx = int(np.argmin(fitness_sequence)) if fitness_sequence else 0
            all_solutions.append(run_histories[best_idx])
        best_solutions.append(run_best)

    return best_solutions, all_solutions


__all__ = [
    "Solution",
    "SolutionSet",
    "repair_sol",
    "make_solutions",
    "run_design",
    "input_algorithm",
    "run_algorithm",
]



