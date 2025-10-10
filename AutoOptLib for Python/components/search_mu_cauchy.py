"""Python translation of search_mu_cauchy."""

from __future__ import annotations

import numpy as np


def search_mu_cauchy(*args):
    mode = args[-1]
    if mode == "execute":
        parent = args[0]
        aux = args[3] if len(args) > 3 else None
        inner_g = int(args[5]) if len(args) > 5 else 1

        if not isinstance(parent, (np.ndarray, list, tuple)):
            decs = getattr(parent, "decs", None)
            parent = decs() if callable(decs) else (decs if decs is not None else parent)
        parent = np.asarray(parent, dtype=float)
        n, d = parent.shape

        if aux is None or not hasattr(aux, "cauchy_eta") or inner_g == 1:
            eta = np.random.rand(n, d)
            if aux is not None:
                setattr(aux, "cauchy_eta", eta)
        else:
            eta = getattr(aux, "cauchy_eta")
            if eta is None or np.shape(eta) != (n, d):
                eta = np.random.rand(n, d)
        disturb = eta * np.random.standard_cauchy(size=(n, d))
        offspring = parent + disturb

        tau1 = 1.0 / np.sqrt(2.0 * np.sqrt(d))
        tau2 = 1.0 / np.sqrt(2.0 * d)
        normal = np.random.randn(n, 1).repeat(d, axis=1)
        normal_j = np.random.randn(n, d)
        eta = eta * np.exp(tau2 * normal + tau1 * normal_j)
        if aux is not None:
            setattr(aux, "cauchy_eta", eta)
        return offspring, aux

    if mode == "parameter":
        return None, None

    if mode == "behavior":
        return ["", "GS"], None

    raise ValueError(f"Unsupported mode: {mode}")
