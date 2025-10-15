from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, List, Sequence

import numpy as np

from ..design._helpers import get_problem_type, get_flex


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
    ptype = get_problem_type(problem)
    bound = getattr(problem, "bound", None)
    if bound is None:
        return dec
    bound = np.asarray(bound)
    if ptype == "continuous":
        lower = bound[0]
        upper = bound[1]
        return np.minimum(np.maximum(dec, lower), upper)
    if ptype == "discrete":
        lower = bound[0]
        upper = bound[1]
        repaired = np.clip(np.rint(dec), lower, upper).astype(int)
        setting = get_flex(problem, "setting", "")
        if isinstance(setting, str) and "dec_diff" in setting:
            n, d = repaired.shape
            for i in range(n):
                seen = set()
                duplicates = []
                for j, val in enumerate(repaired[i]):
                    if val in seen:
                        duplicates.append(j)
                    else:
                        seen.add(val)
                if not duplicates:
                    continue
                for idx in duplicates:
                    choices = [v for v in range(int(lower[idx]), int(upper[idx]) + 1) if v not in repaired[i]]
                    if choices:
                        repaired[i, idx] = int(np.random.choice(choices))
            return repaired
        return repaired
    if ptype == "permutation":
        repaired = np.rint(dec).astype(int)
        if repaired.ndim == 1:
            rows = [repaired.copy()]
        else:
            rows = [row.copy() for row in repaired]
        d = len(rows[0])
        domain = list(range(1, d + 1))
        for row in rows:
            missing = [v for v in domain if v not in row]
            counts = {}
            duplicates = []
            for j, val in enumerate(row):
                counts[val] = counts.get(val, 0) + 1
                if counts[val] > 1:
                    duplicates.append(j)
            for j in duplicates:
                if missing:
                    row[j] = missing.pop(0)
        if repaired.ndim == 1:
            return rows[0]
        return np.vstack(rows)
    return dec


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
        d = np.array(decs[i], dtype=float)
        d = repair_sol(d, problem)
        obj, con, _ = _call_evaluator(problem, data, d)
        feasible = con <= 0.0
        fit = obj if feasible else (con + 1e8)
        items.append(Solution(dec=d, obj=obj, con=con, fit=fit))
    return SolutionSet(items)
