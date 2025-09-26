"""Shared helpers for the AutoOpt component translations."""

from __future__ import annotations

import math
import re
from typing import Any, Iterable, Sequence

import numpy as np


def _name_variants(name: str) -> list[str]:
    variants = {name}
    variants.add(name.lower())
    variants.add(name.upper())
    snake = re.sub(r"(?<!^)(?=[A-Z])", "_", name).lower()
    variants.add(snake)
    if "_" in name:
        parts = name.split("_")
        camel = parts[0].lower() + "".join(part.title() for part in parts[1:])
        pascal = "".join(part.title() for part in parts)
        variants.add(camel)
        variants.add(pascal)
    else:
        parts = re.findall(r"[A-Za-z][^A-Z]*", name)
        if parts:
            variants.add("_".join(part.lower() for part in parts))
            camel = parts[0].lower() + "".join(part.title() for part in parts[1:])
            pascal = "".join(part.title() for part in parts)
            variants.add(camel)
            variants.add(pascal)
    return [v for v in variants if v]


def flex_get(obj: Any, name: str, default: Any = None) -> Any:
    if obj is None:
        return default
    for candidate in _name_variants(name):
        if isinstance(obj, dict) and candidate in obj:
            return obj[candidate]
        if hasattr(obj, candidate):
            attr = getattr(obj, candidate)
            return attr() if callable(attr) else attr
    return default


def ensure_rng(*candidates: Any) -> np.random.Generator:
    for candidate in candidates:
        if isinstance(candidate, np.random.Generator):
            return candidate
        if isinstance(candidate, dict):
            nested = candidate.get("rng")
            if isinstance(nested, np.random.Generator):
                return nested
        if hasattr(candidate, "rng"):
            nested = getattr(candidate, "rng")
            if isinstance(nested, np.random.Generator):
                return nested
    return np.random.default_rng()


def to_numpy(data: Any) -> np.ndarray:
    return np.asarray(data)


def pairwise_distances(a: np.ndarray) -> np.ndarray:
    diff = a[:, None, :] - a[None, :, :]
    return np.linalg.norm(diff, axis=-1)


def randperm(n: int, rng: np.random.Generator) -> np.ndarray:
    return rng.permutation(n)


def reshape_pairs(index: Sequence[int]) -> np.ndarray:
    array = np.asarray(index)
    if array.size % 2:
        raise ValueError("Index array must have even length for pairing")
    return array.reshape(-1, 2)


def as_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def ensure_column(vector: Any) -> np.ndarray:
    arr = np.asarray(vector)
    if arr.ndim == 1:
        return arr.reshape(-1, 1)
    return arr

