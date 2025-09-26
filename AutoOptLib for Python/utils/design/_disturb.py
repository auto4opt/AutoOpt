"""Placeholder for the MATLAB Disturb routine."""

from __future__ import annotations

from typing import Any, Sequence, Tuple


def disturb(algs: Sequence[Any], setting: Any, inner_g: int, aux: Any) -> Tuple[Any, Any, Any]:
    """Design new algorithms based on the current ones.

    The MATLAB implementation performs a sophisticated combination of
    structural mutations over operator graphs and parameter perturbations.
    This port will be completed once the dependent operator primitives are
    available on the Python side.
    """
    raise NotImplementedError("disturb() is pending a full translation from MATLAB")
