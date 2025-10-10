"""Python translation of search_de_current."""

from __future__ import annotations

import numpy as np


def search_de_current(*args):
    mode = args[-1]
    if mode == "execute":
        parent_obj = args[0]
        para = args[2] if len(args) > 2 else None
        aux = args[3] if len(args) > 3 else None  # noqa: F841
        decs = getattr(parent_obj, "decs", None)
        if callable(decs):
            parent = decs()
        else:
            parent = decs if decs is not None else parent_obj
        parent = np.asarray(parent, dtype=float)
        n, d = parent.shape
        if para is None:
            f, cr = 0.5, 0.5
        else:
            arr = np.asarray(para).reshape(-1)
            f = float(arr[0])
            cr = float(arr[1]) if arr.size > 1 else 0.5
        p1 = parent.copy()
        p2 = parent[np.random.permutation(n)]
        p3 = parent[np.random.permutation(n)]
        mask = np.random.rand(n, d) < cr
        offspring = parent.copy()
        offspring[mask] = p1[mask] + f * (p2[mask] - p3[mask])
        return offspring, args[3] if len(args) > 3 else None
    if mode == "parameter":
        return [[0, 1], [0, 1]], None
    if mode == "behavior":
        return [["LS", "small", "small"], ["GS", "large", "large"]], None
    raise ValueError(f"Unsupported mode: {mode}")
