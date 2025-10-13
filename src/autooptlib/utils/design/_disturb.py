from __future__ import annotations

from typing import Any, List, Sequence, Tuple

import numpy as np

from ._helpers import ensure_rng, get_flex


def _as_int(value: Any) -> int:
    arr = np.asarray(value)
    if arr.size == 0:
        raise ValueError("Cannot convert empty array to int")
    return int(arr.reshape(-1)[0])


def _ensure_alg_list(algs: Any) -> List[Any]:
    if isinstance(algs, Sequence) and not isinstance(algs, (str, bytes)):
        return list(algs)
    return [algs]



def _copy_parameter_structure(parameter: Sequence[Any]) -> List[List[Any]]:
    copied: List[List[Any]] = []
    for entry in parameter:
        if entry is None:
            copied.append([None, None])
            continue
        values = entry[0] if len(entry) > 0 else None
        behavior = entry[1] if len(entry) > 1 else None
        values_copy = None if values is None else np.array(values, copy=True)
        copied.append([values_copy, behavior])
    return copied


def disturb(algs: Any, setting: Any, inner_g: int, aux: Any) -> Tuple[List[List[np.ndarray]], List[List[List[Any]]], Any]:
    alg_list = _ensure_alg_list(algs)
    rng = ensure_rng(setting)

    op_space = np.asarray(get_flex(setting, "op_space", required=True), dtype=int)
    para_space = list(get_flex(setting, "para_space", required=True))
    para_local_space = list(get_flex(setting, "para_local_space", default=[None] * len(para_space)))
    alg_p = int(get_flex(setting, "alg_p", required=True))
    alg_q = int(get_flex(setting, "alg_q", required=True))
    tune_para = bool(get_flex(setting, "tune_para", False))

    ind_non_empty_para = [idx for idx, space in enumerate(para_space) if space is not None and len(space) > 0]

    aux_list: List[Any]
    if aux is None:
        aux_list = [{} for _ in alg_list]
    elif isinstance(aux, Sequence):
        aux_list = list(aux)
        while len(aux_list) < len(alg_list):
            aux_list.append({})
    else:
        aux_list = [aux] + [{} for _ in range(len(alg_list) - 1)]

    new_ops: List[List[np.ndarray]] = []
    new_paras: List[List[List[Any]]] = []

    for idx, alg in enumerate(alg_list):
        this_op_raw = getattr(alg, "operator", None)
        this_para_raw = getattr(alg, "parameter", None)
        if not this_op_raw or this_para_raw is None:
            raise ValueError("Algorithm must provide operator and parameter fields for disturb().")

        if isinstance(this_op_raw[0], (list, tuple)):
            path_list = [np.array(path, copy=True) for path in this_op_raw[0]]
        else:
            path_list = [np.array(path, copy=True) for path in this_op_raw]
        this_op = path_list
        this_para = _copy_parameter_structure(this_para_raw)

        ind_search: List[int] = []
        for path in this_op:
            if path.shape[0] > 1:
                for val in path[1:, 0]:
                    ind_search.append(_as_int(val))
        ind_op: List[int] = []
        if this_op:
            ind_op.append(_as_int(this_op[0][0, 0]))
            ind_op.extend(ind_search)
            ind_op.append(_as_int(this_op[0][-1, -1]))

        ind_para = [p for p in ind_non_empty_para if (p + 1) in ind_op]

        aux_entry = aux_list[idx] if isinstance(aux_list[idx], dict) else {}

        if inner_g == 1:
            if not tune_para:
                total_ops = len(ind_op)
                total = total_ops + len(ind_para)
                if total == 0:
                    seed = 1
                else:
                    prob_op = 1.0 / total
                    prob_para = prob_op * len(ind_para)
                    probs = [prob_op] * total_ops + [prob_para]
                    # Disable impossible moves
                    if op_space[0, 1] - op_space[0, 0] + 1 <= 1 and probs:
                        probs[0] = 0.0
                    if op_space[1, 1] - op_space[1, 0] + 1 <= 1 and len(probs) >= 3:
                        for j in range(1, len(probs) - 1):
                            probs[j] = 0.0
                    if op_space[2, 1] - op_space[2, 0] + 1 <= 1 and probs:
                        probs[-2 if len(probs) > 1 else 0] = 0.0
                    if len(ind_para) == 0 and probs:
                        probs[-1] = 0.0
                    probs = np.asarray(probs, dtype=float)
                    if probs.sum() <= 0:
                        seed = int(rng.integers(max(1, total_ops))) + 1
                    else:
                        probs = probs / probs.sum()
                        seed = int(rng.choice(len(probs), p=probs)) + 1
            else:
                seed = len(ind_op) + 1
            aux_entry["seed"] = seed
        else:
            seed = int(aux_entry.get("seed", len(ind_op) + 1))

        if seed <= len(ind_op) and len(ind_op) > 0:
            # Operator disturbance
            if seed == 1:
                # Choose operator
                pool = np.arange(op_space[0, 0], op_space[0, 1] + 1)
                current = ind_op[0]
                pool = pool[pool != current]
                if pool.size > 0:
                    ind_new = _as_int(rng.choice(pool))
                    for path in this_op:
                        path[0, 0] = ind_new
            elif seed == len(ind_op):
                # Update operator
                pool = np.arange(op_space[2, 0], op_space[2, 1] + 1)
                current = ind_op[-1]
                pool = pool[pool != current]
                if pool.size > 0:
                    ind_new = _as_int(rng.choice(pool))
                    for path in this_op:
                        path[-1, -1] = ind_new
            else:
                if alg_p != 1:
                    raise NotImplementedError("disturb() currently supports single-path algorithms only")
                path = this_op[0]
                pool = np.arange(op_space[1, 0], op_space[1, 1] + 1)
                current = ind_op[seed - 1]
                pool = pool[pool != current]
                if pool.size == 0:
                    ind_new = current
                else:
                    ind_new = _as_int(rng.choice(pool))
                num_search = path.shape[0] - 1
                if num_search == 1 and num_search < alg_q:
                    sample_pool = list(pool) + [np.inf]
                elif num_search > 1 and num_search < alg_q:
                    sample_pool = list(pool) + [np.inf, -np.inf]
                elif num_search > 1 and num_search == alg_q:
                    sample_pool = list(pool) + [-np.inf]
                else:
                    sample_pool = list(pool)
                if not sample_pool:
                    sample_choice = ind_new
                else:
                    sample_choice = rng.choice(sample_pool)
                row_idx = seed - 1  # zero-based
                if sample_choice == np.inf:
                    insert_row = np.zeros(2, dtype=int)
                    insert_row[0] = ind_new
                    insert_row[1] = path[row_idx, 1]
                    path = np.insert(path, row_idx + 1, insert_row, axis=0)
                    path[row_idx, 1] = ind_new
                elif sample_choice == -np.inf and path.shape[0] > 2:
                    path[row_idx - 1, 1] = path[row_idx, 1]
                    path = np.delete(path, row_idx, axis=0)
                else:
                    path[row_idx, 0] = ind_new
                    path[row_idx - 1, 1] = ind_new
                this_op[0] = path
        else:
            # Parameter disturbance (resample within bounds)
            for para_idx in ind_para:
                entry = this_para[para_idx]
                values = entry[0]
                if values is None:
                    continue
                behavior = entry[1]
                space = para_local_space[para_idx]
                if behavior != "GS" and space is None:
                    space = para_space[para_idx]
                if space is None:
                    continue
                bounds = np.asarray(space, dtype=float)
                if bounds.ndim != 2 or bounds.shape[1] != 2:
                    continue
                lower = bounds[:, 0]
                upper = bounds[:, 1]
                new_vals = lower + (upper - lower) * rng.random(lower.shape)
                entry[0] = new_vals
                this_para[para_idx] = entry

        new_ops.append(this_op)
        new_paras.append(this_para)
        aux_list[idx] = aux_entry

    return new_ops, new_paras, aux_list
