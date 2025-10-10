"""Python translation of cross_point_two."""

from __future__ import annotations

import numpy as np


def cross_point_two(*args):
    mode = args[-1]
    if mode == "execute":
        parent = args[0]
        aux = args[3] if len(args) > 3 else None  # noqa: F841
        decs = getattr(parent, "decs", None)
        if callable(decs):
            parent = decs()
        else:
            parent = decs if decs is not None else parent
        parent = np.asarray(parent)
        n, d = parent.shape
        if d < 2:
            raise ValueError("Two-point crossover requires D >= 2")
        p1 = parent[: (n + 1) // 2]
        p2 = parent[n // 2 :]
        nh = p1.shape[0]
        k = np.sort(np.random.randint(0, d, size=(nh, 2)), axis=1)
        off1 = p1.copy()
        off2 = p2.copy()
        for i in range(nh):
            a, b = k[i]
            off1[i, a : b + 1] = p2[i, a : b + 1]
            off2[i, a : b + 1] = p1[i, a : b + 1]
        offspring = np.vstack([off1, off2])
        return offspring[:n, :], args[3] if len(args) > 3 else None
    if mode == "parameter":
        return None, None
    if mode == "behavior":
        return ["", "GS"], None
    raise ValueError(f"Unsupported mode: {mode}")
