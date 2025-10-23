"""High-level entry point mirroring MATLAB's AutoOpt function."""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Iterable, List

from .utils.general import input_handler, output, process


def _to_arguments(kwargs: dict[str, Any]) -> List[Any]:
    args: List[Any] = []
    for key, value in kwargs.items():
        args.extend([key, value])
    return args


def _load_problem_callable(descriptor: Any):
    if callable(descriptor):
        return descriptor
    if isinstance(descriptor, str):
        module = import_module("autooptlib.problems")
        if hasattr(module, descriptor):
            return getattr(module, descriptor)
    raise ValueError(f"Unsupported problem descriptor: {descriptor!r}")


def autoopt(**kwargs: Any):
    """Python entry point for running AutoOpt-style workflows.

    Parameters mirror MATLAB usage, e.g.:

    autoopt(
        Mode="design",
        Problem="cec2013_f1",
        InstanceTrain=[30, 30, 30],
        InstanceTest=[30],
        AlgFE=400,
        ...
    )
    """

    if "Mode" not in kwargs and "mode" not in kwargs:
        raise ValueError("AutoOpt requires a Mode parameter ('design' or 'solve').")

    mode = kwargs.get("Mode", kwargs.get("mode"))
    setting = SimpleNamespace(Mode=str(mode))

    arguments = _to_arguments(kwargs)

    data = input_handler(arguments, setting, "data")
    if setting.Mode.lower() == "design":
        problem_descriptor, instance_train, instance_test = data
    else:
        problem_descriptor, instance_solve = data

    setting = input_handler(arguments, setting, "parameter")
    setting = input_handler(arguments, setting, "check")

    problem_callable = _load_problem_callable(problem_descriptor)

    if setting.Mode.lower() == "design":
        final_algs, alg_trace = process(
            problem_callable,
            instance_train,
            instance_test,
            setting=setting,
        )
        output(final_algs, alg_trace, instance_train, instance_test, setting=setting)
        return final_algs, alg_trace

    best_solutions, all_solutions = process(
        problem_callable,
        instance_solve,
        setting=setting,
    )
    output(best_solutions, all_solutions, instance_solve, setting=setting)
    return best_solutions, all_solutions


__all__ = ["autoopt"]
