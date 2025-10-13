from __future__ import annotations

from typing import Any, List, Sequence

import numpy as np

from ...components import get_component
from ..general.improve_rate import improve_rate
from ..solve import SolutionSet, make_solutions, repair_sol
from ._helpers import Pathway, PathwayParam, get_flex


def _init_population(problem: Any, data: Any, setting: Any) -> SolutionSet:
    ptype = get_flex(problem, "type", ["continuous"])
    pop_n = int(get_flex(problem, "N", 10))
    if isinstance(ptype, (list, tuple)):
        ptype = ptype[0]
    if ptype != "continuous":
        raise NotImplementedError("Only continuous problems supported in minimal evaluate()")
    bound = get_flex(problem, "bound")
    lower = np.asarray(bound[0], dtype=float)
    upper = np.asarray(bound[1], dtype=float)
    decs = lower + (upper - lower) * np.random.rand(pop_n, lower.shape[-1])
    return make_solutions(decs, problem, data)


def evaluate(self, problem: Any, data: Any, setting: Any, seed_instance: Sequence[int]):
    if not self.operator_pheno:
        return self

    pathway: Pathway = self.operator_pheno[0][0]
    path_params: PathwayParam = self.parameter_pheno[0][0]

    solutions = _init_population(problem, data, setting)
    aux_cache: List[Any] = [None] * len(pathway.search)
    archive_names = list(get_flex(setting, "archive", []))
    archives: List[Any] = [[] for _ in archive_names]

    metric = get_flex(setting, "metric", "quality")
    prob_n = float(get_flex(setting, "prob_n", 1.0))
    tmax = get_flex(setting, "tmax", np.inf)
    thres = float(get_flex(setting, "thres", -np.inf))
    gmax = int(get_flex(problem, "Gmax", 10))

    if metric in ("quality", "auc"):
        Tmax = np.inf
        Thres = -np.inf
    elif metric == "runtimeFE":
        Tmax = int(np.ceil(tmax / max(prob_n, 1)))
        Thres = thres
    else:
        Tmax = tmax
        Thres = thres

    G = 1
    t = 0.0

    while G <= gmax and t < Tmax and float(np.min(solutions.fits())) > Thres:
        for step_idx, step in enumerate(pathway.search):
            improve = None
            innerG = 1
            inner_limit = int(step.termination[1]) if step.termination.size > 1 else 1
            rate_threshold = float(step.termination[0]) if step.termination.size > 0 else -np.inf
            while (improve is None or improve[0] >= rate_threshold) and innerG <= inner_limit:
                choose_fn = get_component(pathway.choose)
                choose_param = getattr(path_params, "choose", None)
                ind, _ = choose_fn(solutions, problem, choose_param, None, G, innerG, data, "execute")
                ind = np.asarray(ind, dtype=int).reshape(-1)

                search_params = getattr(path_params, "search", [])
                primary_param = search_params[step_idx].primary if step_idx < len(search_params) else None
                secondary_param = search_params[step_idx].secondary if step_idx < len(search_params) else None

                parent_subset = SolutionSet([solutions[i] for i in ind]) if ind.size else solutions
                primary_fn = get_component(step.primary)
                aux_state = aux_cache[step_idx] or {}
                new_dec, aux_state = primary_fn(parent_subset, problem, primary_param, aux_state, G, innerG, data, "execute")

                if step.secondary:
                    new_dec = repair_sol(np.asarray(new_dec), problem)
                    secondary_fn = get_component(step.secondary)
                    new_dec, aux_state = secondary_fn(new_dec, problem, secondary_param, aux_state, G, innerG, data, "execute")

                new = make_solutions(np.asarray(new_dec), problem, data)

                update_fn = get_component(pathway.update)
                updated, _ = update_fn(list(solutions) + list(new), problem, getattr(path_params, "update", None), None, G, innerG, data, "execute")
                solutions = SolutionSet(updated)

                for j, name in enumerate(archive_names):
                    if name == "archive_best":
                        archive_fn = get_component(name)
                        archives[j], _ = archive_fn(solutions, archives[j], problem, "execute")

                if step.primary == "search_cma":
                    aux_state = get_component("para_cma")(new, problem, aux_state, "solution")

                improve = improve_rate(solutions, improve, innerG, "solution")
                aux_cache[step_idx] = aux_state
                innerG += 1
                G += 1
                if metric == "runtimeFE":
                    t = G
                if G > gmax or t >= Tmax or float(np.min(solutions.fits())) <= Thres:
                    break
            if G > gmax or t >= Tmax or float(np.min(solutions.fits())) <= Thres:
                break
        if metric != "runtimeFE":
            t += 1.0

    best = float(np.min(solutions.fits())) if len(solutions) else np.inf
    runs = int(get_flex(setting, "alg_runs", 1))
    for inst in seed_instance:
        if metric in ("quality", "auc"):
            self.performance[inst, :runs] = best
        else:
            self.performance[inst, :runs] = t
    return self
