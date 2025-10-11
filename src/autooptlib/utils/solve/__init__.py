from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, List, Sequence

import numpy as np


@dataclass
class Solution:
    dec: np.ndarray
    obj: float
    con: float
    fit: float


class SolutionSet(Sequence[Solution]):
    def __init__(self, items: Iterable[Solution]):
        self._items: List[Solution] = list(items)

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, idx):
        return self._items[idx]

    # Aggregations compatible with components
    def decs(self) -> np.ndarray:
        if not self._items:
            return np.zeros((0, 0))
        return np.vstack([s.dec for s in self._items])

    def objs(self) -> np.ndarray:
        return np.asarray([s.obj for s in self._items]).reshape(-1, 1)

    def cons(self) -> np.ndarray:
        return np.asarray([s.con for s in self._items]).reshape(-1, 1)

    def fits(self) -> np.ndarray:
        return np.asarray([s.fit for s in self._items]).reshape(-1, 1)


def repair_sol(dec: np.ndarray, problem: Any) -> np.ndarray:
    bound = getattr(problem, 'bound', None)
    if bound is None:
        return dec
    lower = np.asarray(bound[0], dtype=float)
    upper = np.asarray(bound[1], dtype=float)
    return np.minimum(np.maximum(dec, lower), upper)


def _call_evaluator(problem: Any, data: Any, dec: np.ndarray):
    # Preferred: problem.evaluate(data, dec)
    eval_fn = getattr(problem, 'evaluate', None)
    if callable(eval_fn):
        out = eval_fn(data, dec)
    else:
        # Fallback: Problem.name referencing a callable in current scope
        name = getattr(problem, 'name', None)
        if callable(name):
            out = name(data, dec)
        else:
            raise NotImplementedError('Problem must provide an evaluate(data, dec) callable')
    # Normalize outputs: obj, con, acc(optional)
    if isinstance(out, tuple):
        if len(out) == 3:
            obj, con, acc = out
        elif len(out) == 2:
            obj, con = out
            acc = None
        else:
            obj, con, acc = out[0], 0.0, None
    else:
        obj, con, acc = out, 0.0, None
    return float(np.asarray(obj).reshape(1)[0]), float(np.sum(np.maximum(0.0, np.asarray(con)))), acc


def make_solutions(decs: np.ndarray, problem: Any, data: Any) -> SolutionSet:
    items: List[Solution] = []
    for i in range(decs.shape[0]):
        d = decs[i]
        d = repair_sol(d, problem)
        obj, con, _ = _call_evaluator(problem, data, d)
        feasible = con <= 0.0
        fit = obj if feasible else (con + 1e8)
        items.append(Solution(dec=d, obj=obj, con=con, fit=fit))
    return SolutionSet(items)