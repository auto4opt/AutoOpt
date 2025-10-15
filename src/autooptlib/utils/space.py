from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np

from ..components import get_component
from .design._helpers import get_problem_type


def _to_namespace(setting: Any) -> SimpleNamespace:
    if isinstance(setting, SimpleNamespace):
        return setting
    if isinstance(setting, dict):
        return SimpleNamespace(**setting)
    # Fallback: copy attributes
    data = {k: getattr(setting, k) for k in dir(setting) if not k.startswith("__") and not callable(getattr(setting, k))}
    return SimpleNamespace(**data)


def space(problem: Any, setting: Any) -> SimpleNamespace:
    """Define operator/parameter space for the design problem (continuous only)."""
    set_obj = _to_namespace(setting)
    ptype = get_problem_type(problem)

    if ptype == "continuous":
        choose = ["choose_traverse", "choose_tournament", "choose_roulette_wheel", "choose_nich"]
        search = [
            "search_de_current",
            "search_de_current_best",
            "search_de_random",
            "search_cma",
            "search_pso",
            "cross_point_one",
            "cross_point_two",
            "cross_point_uniform",
            "cross_point_n",
            "search_mu_gaussian",
            "search_mu_cauchy",
            "search_mu_polynomial",
            "search_mu_uniform",
            "reinit_continuous",
        ]
        update = ["update_greedy", "update_round_robin", "update_pairwise", "update_always", "update_simulated_annealing"]
    elif ptype == "discrete":
        choose = ["choose_traverse", "choose_tournament", "choose_roulette_wheel", "choose_nich"]
        search = [
            "cross_point_one",
            "cross_point_two",
            "cross_point_uniform",
            "cross_point_n",
            "search_reset_one",
            "search_reset_rand",
            "search_reset_creep",
            "reinit_discrete",
        ]
        update = ["update_greedy", "update_round_robin", "update_pairwise", "update_always", "update_simulated_annealing"]
    elif ptype == "permutation":
        choose = ["choose_traverse", "choose_tournament", "choose_roulette_wheel", "choose_nich"]
        search = [
            "cross_order_two",
            "cross_order_n",
            "search_swap",
            "search_swap_multi",
            "search_scramble",
            "search_insert",
            "reinit_permutation",
        ]
        update = ["update_greedy", "update_round_robin", "update_pairwise", "update_always", "update_simulated_annealing"]
    else:
        raise NotImplementedError("space() currently supports continuous, discrete, and permutation problems")

    all_op = choose + search + update
    op_space = np.array(
        [
            [1, len(choose)],
            [len(choose) + 1, len(choose) + len(search)],
            [len(choose) + len(search) + 1, len(all_op)],
        ],
        dtype=int,
    )

    para_space = []
    behav_space = []
    para_local_space = []

    for name in all_op:
        try:
            comp = get_component(name)
        except KeyError:
            para_space.append(None)
            behav_space.append(None)
            para_local_space.append(None)
            continue

        try:
            params, _ = comp(problem, "parameter")
        except TypeError:
            params = None
        try:
            behavior, _ = comp("behavior")
        except TypeError:
            behavior = None

        if params is None:
            para_space.append(None)
        else:
            arr = np.asarray(params, dtype=float)
            if arr.ndim == 1:
                if arr.size % 2 != 0:
                    raise ValueError(f"Parameter bounds for {name} must contain pairs")
                arr = arr.reshape(-1, 2)
            elif arr.ndim > 1 and arr.shape[1] != 2:
                total = arr.size
                if total % 2 != 0:
                    raise ValueError(f"Parameter bounds for {name} must contain pairs")
                arr = arr.reshape(-1, 2)
            para_space.append(arr)
        behav_space.append(behavior)
        para_local_space.append(None)

    set_obj.AllOp = all_op
    set_obj.OpSpace = op_space
    set_obj.ParaSpace = para_space
    set_obj.BehavSpace = behav_space
    set_obj.ParaLocalSpace = para_local_space

    # Provide defaults for downstream code
    set_obj.TunePara = getattr(set_obj, "TunePara", False)
    set_obj.alg_p = getattr(set_obj, "alg_p", 1)
    set_obj.alg_q = getattr(set_obj, "alg_q", 3)
    set_obj.alg_n = getattr(set_obj, "alg_n", 1)
    return set_obj
