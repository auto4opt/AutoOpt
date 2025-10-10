"""Python translation of cross_point_uniform."""

from __future__ import annotations

import numpy as np


def cross_point_uniform(*args):
    mode = args[-1]
    if mode == "execute":
        parent = args[0]
        para = args[2] if len(args) > 2 else None
        aux = args[3] if len(args) > 3 else None  # noqa: F841
        prob = float(np.asarray(para).reshape(-1)[0]) if para is not None else 0.5
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
        mask = np.random.rand(nh, d) < prob
        off1 = p1.copy()
        off2 = p2.copy()
        off1[mask] = p2[mask]
        off2[mask] = p1[mask]
        offspring = np.vstack([off1, off2])
        return offspring[:n, :], args[3] if len(args) > 3 else None
    if mode == "parameter":
        return [0, 0.5], None
    if mode == "behavior":
        return [["LS", "small"], ["GS", "large"]], None
    raise ValueError(f"Unsupported mode: {mode}")
