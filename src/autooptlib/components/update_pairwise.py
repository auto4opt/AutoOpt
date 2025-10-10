"""Python translation of update_pairwise."""

from __future__ import annotations

from ._utils import extract_fits, solution_as_list


def update_pairwise(*args):
    mode = args[-1]
    if mode == "execute":
        solution = args[0]
        sol_list = solution_as_list(solution)
        fitness = extract_fits(solution)
        if fitness.size % 2 != 0:
            raise ValueError("update_pairwise expects an even number of solutions")
        half = fitness.size // 2
        selected = []
        for i in range(half):
            old_fit = fitness[i]
            new_fit = fitness[half + i]
            if old_fit <= new_fit:
                selected.append(sol_list[i])
            else:
                selected.append(sol_list[half + i])
        return selected, None

    if mode == "parameter":
        return None, None

    if mode == "behavior":
        return ["", ""], None

    raise ValueError(f"Unsupported mode: {mode}")
