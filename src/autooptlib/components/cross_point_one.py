"""Python translation of cross_point_one."""

from __future__ import annotations

import numpy as np


def cross_point_one(*args):
    mode = args[-1]
    if mode == "execute":
        parent = args[0]
        aux = args[3] if len(args) > 3 else None  # noqa: F841 (kept for API parity)
        decs = getattr(parent, "decs", None)
        if callable(decs):
            parent = decs()
        else:
            parent = decs if decs is not None else parent
        parent = np.asarray(parent)
        n, d = parent.shape
        p1 = parent[: (n + 1) // 2]
        p2 = parent[n // 2 :]
        nh = p1.shape[0]
        k = np.random.randint(0, d, size=(nh,))
        off1 = p1.copy()
        off2 = p2.copy()
        for i in range(nh):
            off1[i, k[i] : ] = p2[i, k[i] : ]
            off2[i, k[i] : ] = p1[i, k[i] : ]
        offspring = np.vstack([off1, off2])
        return offspring[:n, :], args[3] if len(args) > 3 else None
    if mode == "parameter":
        return None, None
    if mode == "behavior":
        return ["", "GS"], None
    raise ValueError(f"Unsupported mode: {mode}")
