from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from autooptlib.components import __all__ as COMPONENT_NAMES, get_component


@dataclass
class SimpleSolution:
    dec: np.ndarray
    obj: float
    con: float
    fit: float


class SimpleSolutionSet(Sequence[SimpleSolution]):
    def __init__(self, items: Iterable[SimpleSolution]):
        self._items: List[SimpleSolution] = list(items)

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, idx):
        return self._items[idx]

    def decs(self) -> np.ndarray:
        if not self._items:
            return np.zeros((0, 0))
        return np.vstack([sol.dec for sol in self._items])

    def objs(self) -> np.ndarray:
        return np.asarray([sol.obj for sol in self._items]).reshape(-1, 1)

    def cons(self) -> np.ndarray:
        return np.asarray([sol.con for sol in self._items]).reshape(-1, 1)

    def fits(self) -> np.ndarray:
        return np.asarray([sol.fit for sol in self._items]).reshape(-1, 1)


class AuxState(dict):
    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError as exc:
            raise AttributeError(item) from exc

    def __setattr__(self, key, value) -> None:
        self[key] = value


@dataclass
class DummyProblem:
    kind: str
    dim: int = 4
    pop_size: int = 6

    def __post_init__(self) -> None:
        if self.kind == "continuous":
            lower = np.full(self.dim, -5.0)
            upper = np.full(self.dim, 5.0)
        elif self.kind == "discrete":
            lower = np.zeros(self.dim, dtype=int)
            upper = np.full(self.dim, 9, dtype=int)
        elif self.kind == "permutation":
            lower = np.ones(self.dim, dtype=int)
            upper = np.full(self.dim, self.dim, dtype=int)
        else:
            raise ValueError(f"Unsupported problem kind: {self.kind}")
        self.bound = np.vstack((lower, upper))
        self.N = self.pop_size
        self.Gmax = 10
        self.type = [self.kind, "static", "certain"]

    def evaluate(self, data, dec):  # noqa: ANN001
        arr = np.asarray(dec, dtype=float)
        obj = float(np.sum(arr ** 2))
        return obj, 0.0, None


def make_population(problem: DummyProblem, rng: np.random.Generator) -> SolutionSet:
    if problem.kind == "continuous":
        lower, upper = problem.bound
        decs = lower + (upper - lower) * rng.random((problem.N, problem.dim))
    elif problem.kind == "discrete":
        lower, upper = problem.bound.astype(int)
        decs = np.vstack(
            [rng.integers(lower[i], upper[i] + 1, size=problem.N) for i in range(problem.dim)]
        ).T
    elif problem.kind == "permutation":
        domain = np.arange(1, problem.dim + 1)
        decs = np.vstack([rng.permutation(domain) for _ in range(problem.N)])
    else:
        raise ValueError(f"Unknown problem type: {problem.kind}")

    items: List[SimpleSolution] = []
    for row in decs:
        dec = np.asarray(row, dtype=float if problem.kind == "continuous" else int)
        obj = float(np.sum(np.asarray(row, dtype=float) ** 2))
        items.append(SimpleSolution(dec=dec, obj=obj, con=0.0, fit=obj))
    return SimpleSolutionSet(items)


def solution_set_from_decs(decs: np.ndarray, problem: DummyProblem) -> SolutionSet:
    rng = np.random.default_rng(9876)
    items: List[SimpleSolution] = []
    for row in np.asarray(decs):
        dec = np.asarray(row, dtype=float if problem.kind == "continuous" else int)
        obj = float(np.sum(np.asarray(row, dtype=float) ** 2))
        items.append(SimpleSolution(dec=dec, obj=obj, con=0.0, fit=obj + rng.random() * 1e-6))
    return SimpleSolutionSet(items)


def unpack(result):
    if isinstance(result, tuple) and len(result) >= 1:
        return result[0]
    return result


def sample_parameter(spec):
    if spec is None:
        return None
    arr = np.asarray(spec, dtype=float)
    if arr.ndim == 0:
        return float(arr)
    if arr.ndim == 1:
        if arr.size == 0:
            return None
        if arr.size == 2:
            return float(arr.mean())
        return arr
    if arr.ndim == 2 and arr.shape[1] == 2:
        return arr.mean(axis=1)
    return arr


def make_aux(seed_counter: Iterable[int]) -> AuxState:
    seed = next(seed_counter)
    return AuxState(rng=np.random.default_rng(10_000 + seed))


def validate_components() -> None:
    seed_counter = itertools.count()
    rng = np.random.default_rng(12345)

    problem_cont = DummyProblem("continuous")
    problem_disc = DummyProblem("discrete")
    problem_perm = DummyProblem("permutation")

    pop_cont = make_population(problem_cont, rng)
    pop_disc = make_population(problem_disc, rng)
    pop_perm = make_population(problem_perm, rng)

    doubled_pop = list(pop_cont) + list(pop_cont)

    discrete_components = {
        "cross_point_one",
        "cross_point_two",
        "cross_point_n",
        "cross_point_uniform",
        "reinit_discrete",
        "search_reset_one",
        "search_reset_rand",
        "search_reset_creep",
    }
    permutation_components = {
        "cross_order_two",
        "cross_order_n",
        "reinit_permutation",
        "search_swap",
        "search_swap_multi",
        "search_scramble",
        "search_insert",
    }

    errors: List[Tuple[str, Exception]] = []

    def get_context(name: str) -> Tuple[SimpleSolutionSet, DummyProblem]:
        if name in permutation_components:
            return pop_perm, problem_perm
        if name in discrete_components:
            return pop_disc, problem_disc
        return pop_cont, problem_cont

    # Ensure search_cma and search_pso states are captured for para components
    extra_states: Dict[str, Tuple[np.ndarray, dict]] = {}

    for name in sorted(COMPONENT_NAMES):
        func = get_component(name)
        try:
            if name.startswith("para_"):
                continue  # handled later

            population, problem = get_context(name)

            param_spec = func(problem, "parameter") if name not in {"archive_tabu"} else func("parameter")
            param_spec = unpack(param_spec)
            param_value = sample_parameter(param_spec)

            behavior = func("behavior")
            unpack(behavior)  # ensure callable exists

            if name.startswith("archive_"):
                if name == "archive_best":
                    func(population, [], "execute")
                elif name == "archive_diversity":
                    func(population, [], problem, "execute")
                elif name == "archive_statistic":
                    func(population, [], "execute")
                elif name == "archive_tabu":
                    func(population, "execute")
            elif name.startswith("choose_"):
                aux = make_aux(seed_counter)
                func(population, problem, param_value, aux, 0, "execute")
            elif name.startswith("cross_"):
                aux = make_aux(seed_counter)
                func(population, problem, param_value, aux, "execute")
            elif name.startswith("search_"):
                aux = make_aux(seed_counter)
                offspring, aux_out = func(population, problem, param_value, aux, "execute")
                if name in {"search_cma", "search_pso"}:
                    extra_states[name] = (np.asarray(offspring), aux_out if isinstance(aux_out, dict) else {})
            elif name.startswith("reinit_"):
                aux = make_aux(seed_counter)
                func(population, problem, param_value, aux, "execute")
            elif name.startswith("update_"):
                aux = make_aux(seed_counter)
                g_val = 1
                func(doubled_pop, problem_cont, param_value, aux, g_val, "execute")
            else:
                raise ValueError(f"Unhandled component category for {name}")
        except Exception as exc:  # pylint: disable=broad-except
            errors.append((name, exc))

    # Validate para components separately
    try:
        search_cma_func = get_component("search_cma")
        if "search_cma" not in extra_states:
            aux = make_aux(seed_counter)
            offspring, aux_out = search_cma_func(pop_cont, problem_cont, None, aux, "execute")
            extra_states["search_cma"] = (np.asarray(offspring), aux_out if isinstance(aux_out, dict) else {})
        decs, aux_cma = extra_states["search_cma"]
        sol_set = solution_set_from_decs(decs, problem_cont)
        get_component("para_cma")(sol_set, problem_cont, aux_cma, "solution")
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(("para_cma", exc))

    try:
        search_pso_func = get_component("search_pso")
        if "search_pso" not in extra_states:
            aux = make_aux(seed_counter)
            offspring, aux_out = search_pso_func(pop_cont, problem_cont, None, aux, "execute")
            extra_states["search_pso"] = (np.asarray(offspring), aux_out if isinstance(aux_out, dict) else {})
        decs, aux_pso = extra_states["search_pso"]
        sol_set = solution_set_from_decs(decs, problem_cont)
        get_component("para_pso")(sol_set, problem_cont, aux_pso)
    except Exception as exc:  # pylint: disable=broad-except
        errors.append(("para_pso", exc))

    if errors:
        print("Component validation FAILED:")
        for name, err in errors:
            print(f" - {name}: {err}")
        raise SystemExit(1)

    print(f"All {len(COMPONENT_NAMES)} components executed successfully.")


if __name__ == "__main__":
    validate_components()
