"""Track mean and std of fitness."""
from __future__ import annotations

from typing import Any, Sequence

import numpy as np

from ._utils import flex_get


def archive_statistic(*args):
    mode = args[-1]
    if mode == "execute":
        solution = args[0]
        archive = args[1] if len(args) > 1 else []

        fits = flex_get(solution, "fits")
        if callable(fits):
            values = np.asarray(fits()).reshape(-1)
        else:
            values = np.asarray(fits).reshape(-1)
        stat = np.array([[np.mean(values), np.std(values)]])
        if archive is None or len(archive) == 0:
            return stat, None
        archive_arr = np.asarray(archive)
        return np.vstack([archive_arr, stat]), None

    if mode == "parameter":
        return None, None

    if mode == "behavior":
        return ["", ""], None

    raise ValueError(f"Unsupported mode: {mode}")