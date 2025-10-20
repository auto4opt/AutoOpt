from __future__ import annotations

try:  # noqa: SIM105
    from .utils.design import Design  # noqa: F401
except Exception:  # pylint: disable=broad-except
    Design = None  # type: ignore

try:  # noqa: SIM105
    from .components import get_component  # noqa: F401
except Exception:  # pylint: disable=broad-except
    get_component = None  # type: ignore

try:  # noqa: SIM105
    from .problems.cec2013 import cec2013_f1  # noqa: F401
except Exception:  # pylint: disable=broad-except
    cec2013_f1 = None  # type: ignore

__all__ = [
    name
    for name in ("Design", "get_component", "cec2013_f1")
    if name in globals() and globals()[name] is not None
]
