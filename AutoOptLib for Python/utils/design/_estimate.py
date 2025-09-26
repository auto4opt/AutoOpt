"""Python port of the MATLAB Estimate routine."""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np


def _call_use_embed(surrogate: Any, new_algs: Sequence[Any], setting: Any):
    if surrogate is None:
        raise ValueError("Surrogate must be provided")
    if hasattr(surrogate, "UseEmbed") and callable(getattr(surrogate, "UseEmbed")):
        return surrogate.UseEmbed(new_algs, setting)
    if hasattr(surrogate, "use_embed") and callable(getattr(surrogate, "use_embed")):
        return surrogate.use_embed(new_algs, setting)
    raise AttributeError("Surrogate does not provide a use_embed/UseEmbed method")


def _predict(model: Any, data: Any):
    if model is None:
        raise ValueError("Surrogate model is required for prediction")
    if hasattr(model, "predict"):
        return model.predict(data)
    raise AttributeError("Surrogate model lacks a predict method")


def estimate(new_algs: Sequence[Any], problem: Any, setting: Any, ind_instance: Sequence[int], surrogate: Any):
    """Estimate the designed algorithm's performance using a surrogate model."""
    embed_algs = _call_use_embed(surrogate, new_algs, setting)
    model = getattr(surrogate, "model", None)
    for alg_idx, alg in enumerate(new_algs):
        features = embed_algs[alg_idx]
        prediction = np.asarray(_predict(model, features))
        for instance_index in ind_instance:
            alg.performance_approx[instance_index, :] = prediction
    return new_algs
