"""AutoOpt component library (Python port)."""

from __future__ import annotations

from importlib import import_module
from typing import Callable, Dict

_COMPONENT_MODULES = {
    "choose_traverse": "choose_traverse",
    "choose_tournament": "choose_tournament",
    "choose_roulette_wheel": "choose_roulette_wheel",
    "choose_nich": "choose_nich",
    "update_always": "update_always",
    "update_greedy": "update_greedy",
    "update_round_robin": "update_round_robin",
    "update_pairwise": "update_pairwise",
    "update_simulated_annealing": "update_simulated_annealing",
    "cross_point_one": "cross_point_one",
    "cross_point_two": "cross_point_two",
    "cross_point_uniform": "cross_point_uniform",
    "cross_point_n": "cross_point_n",
    "search_de_current": "search_de_current",
    "search_de_current_best": "search_de_current_best",
    "search_de_random": "search_de_random",
    "search_mu_gaussian": "search_mu_gaussian",
    "search_mu_cauchy": "search_mu_cauchy",
    "search_mu_uniform": "search_mu_uniform",
    "search_mu_polynomial": "search_mu_polynomial",
    "search_cma": "search_cma",
    "search_pso": "search_pso",
    "reinit_continuous": "reinit_continuous",
    "reinit_discrete": "reinit_discrete",
    "reinit_permutation": "reinit_permutation",
    "search_reset_one": "search_reset_one",
    "search_reset_rand": "search_reset_rand",
    "search_reset_creep": "search_reset_creep",
    "search_swap": "search_swap",
    "search_swap_multi": "search_swap_multi",
    "search_scramble": "search_scramble",
    "search_insert": "search_insert",
    "cross_order_two": "cross_order_two",
    "cross_order_n": "cross_order_n",
    "archive_best": "archive_best",
    "archive_diversity": "archive_diversity",
    "archive_statistic": "archive_statistic",
    "para_cma": "para_cma",
    "para_pso": "para_pso",
}

_cache: Dict[str, Callable] = {}

def get_component(name: str) -> Callable:
    if name not in _cache:
        module_name = _COMPONENT_MODULES.get(name)
        if module_name is None:
            raise KeyError(f"Component {name!r} not registered")
        module = import_module(f".{module_name}", package=__name__)
        _cache[name] = getattr(module, name)
    return _cache[name]

__all__ = sorted(_COMPONENT_MODULES)

