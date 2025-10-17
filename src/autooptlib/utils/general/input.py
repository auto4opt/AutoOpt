"""Translation of MATLAB Utilities/General/Input.m."""

from __future__ import annotations

import warnings
from types import SimpleNamespace
from typing import Any, Iterable, Sequence


_PARAM_KEYS = {
    "AlgP",
    "AlgQ",
    "Archive",
    "LSRange",
    "IncRate",
    "ProbN",
    "ProbFE",
    "InnerFE",
    "AlgN",
    "AlgFE",
    "AlgRuns",
    "Metric",
    "Compare",
    "Evaluate",
    "Tmax",
    "Thres",
    "RacingK",
    "Surro",
    "AlgFile",
    "AlgName",
}


def _to_sequence(value: Any) -> Sequence[Any]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return value
    return (value,)


def _ensure_namespace(setting: Any) -> SimpleNamespace:
    if isinstance(setting, SimpleNamespace):
        return setting
    if isinstance(setting, dict):
        return SimpleNamespace(**setting)
    data = {
        name: getattr(setting, name)
        for name in dir(setting)
        if not name.startswith("__") and not callable(getattr(setting, name))
    }
    return SimpleNamespace(**data)


def _find_argument(arguments: Sequence[Any], name: str) -> tuple[bool, Any]:
    if not arguments:
        return False, None
    try:
        idx = arguments.index(name)
    except ValueError:
        return False, None
    if idx + 1 >= len(arguments):
        raise ValueError(f'Missing value for argument "{name}"')
    return True, arguments[idx + 1]


def _to_list(obj: Any) -> list[Any]:
    if isinstance(obj, list):
        return obj
    if isinstance(obj, tuple):
        return list(obj)
    return [obj]


def input_handler(arguments: Iterable[Any], setting: Any, mode: str):
    """Reimplementation of the MATLAB Input.m logic."""
    args = list(arguments)
    set_ns = _ensure_namespace(setting)

    if mode == "data":
        if "Problem" not in args:
            raise ValueError("Please set the targeted problem.")
        _, problem = _find_argument(args, "Problem")
        if set_ns.Mode == "design":
            if "InstanceTrain" not in args or "InstanceTest" not in args:
                raise ValueError("Please set the targeted problem instance indexes.")
            _, train = _find_argument(args, "InstanceTrain")
            _, test = _find_argument(args, "InstanceTest")
            return problem, _to_list(train), _to_list(test)

        if set_ns.Mode == "solve":
            if "InstanceSolve" not in args:
                raise ValueError("Please set the targeted problem instance indexes.")
            _, solve = _find_argument(args, "InstanceSolve")
            return problem, _to_list(solve)

        raise ValueError('Please set the mode to "design" or "solve".')

    if mode == "parameter":
        for key in _PARAM_KEYS:
            found, value = _find_argument(args, key)
            if found:
                setattr(set_ns, key, value)
        return set_ns

    if mode == "check":
        _check_setting(set_ns)
        return set_ns

    raise ValueError(f"Unsupported mode: {mode}")


def _check_setting(setting: SimpleNamespace) -> None:
    mode = getattr(setting, "Mode", None)
    if mode not in {"design", "solve"}:
        raise ValueError('Please set the mode to "design" or "solve".')

    if mode == "design":
        if getattr(setting, "AlgP", 0) > getattr(setting, "AlgQ", 0):
            raise ValueError("The number of pathways should not be larger than the number of operators.")
        if getattr(setting, "AlgN", 0) > getattr(setting, "AlgFE", 0):
            raise ValueError("The number of algorithms should not be larger than the evaluation budget.")
        if getattr(setting, "AlgRuns", 0) > getattr(setting, "ProbFE", 0):
            raise ValueError("The number of runs should not exceed problem evaluations.")

        evaluate = getattr(setting, "Evaluate", "exact")
        compare = getattr(setting, "Compare", "average")
        if evaluate == "racing" and compare != "statistic":
            raise ValueError(
                'The "racing" evaluation method should be used with the algorithm comparing method of "statistic".'
            )
        if evaluate == "racing" and not getattr(setting, "RacingK", None):
            raise ValueError(
                'Please set "Setting.K" as the number of instances evaluated before the first round of racing.'
            )
        if evaluate == "approximate" and not getattr(setting, "Surro", None):
            raise ValueError(
                'Please set "Setting.Surro" as the number of exact performance evaluations when using surrogate.'
            )
        if evaluate == "approximate" and compare == "statistic":
            raise ValueError(
                'It is not necessary to use the "statistic" algorithm comparing method when using the "approximate" evaluation method.'
            )
        if compare == "statistic" and getattr(setting, "AlgRuns", 1) == 1:
            raise ValueError(
                'Please run the design multiple times (Setting.AlgRuns>1) when using the "statistic" comparsion method.'
            )
        if getattr(setting, "ProbN", 0) < 5 and getattr(setting, "AlgP", 0) > 1:
            warnings.warn("It is better to have a large population size if involving the EDA operator", stacklevel=2)
        if getattr(setting, "AlgQ", 0) > 4:
            warnings.warn(
                "AlgQ is recommended to be larger than 4 for discrete and permutation problems due to the lack of so many search operators",
                stacklevel=2,
            )
        if compare == "statistic":
            warnings.warn(
                "It is better to have a large number of training intrances or have a large number of algorithm runs "
                "(set AlgRun to a large number), in order to make the statistical test discriminative.",
                stacklevel=2,
            )

    elif mode == "solve":
        if not getattr(setting, "AlgFile", None) and not getattr(setting, "AlgName", None):
            raise ValueError(
                "Please specify an algorithm file in Setting.AlgFile or specify an algorithm name in Setting.AlgName."
            )

        metric = getattr(setting, "Metric", "quality")
        if metric == "runtimeFE":
            if getattr(setting, "Tmax", None) in (None, []):
                setting.Tmax = getattr(setting, "ProbFE", None)
            if getattr(setting, "Thres", None) in (None, []):
                raise ValueError(
                    'Please set "Setting.Thres" as the lowest acceptable performance of the design algorithms, '
                    "the performance can be the solution quality."
                )
        if metric == "runtimeSec":
            if getattr(setting, "Tmax", None) in (None, []):
                raise ValueError('Please set "Setting.Tmax" as the maximum runtime (seconds).')
            if getattr(setting, "Thres", None) in (None, []):
                raise ValueError(
                    'Please set "Setting.Thres" as the lowest acceptable performance of the design algorithms, '
                    "the performance can be the solution quality."
                )
        if metric == "auc":
            tmax = getattr(setting, "Tmax", None)
            thres = getattr(setting, "Thres", None)
            if not isinstance(tmax, Sequence) or len(tmax) <= 1:
                raise ValueError(
                    '"Setting.Tmax" should contain multiple time points. The time points should the numbers of '
                    "function evaluations spent during the alorithm execution."
                )
            if not isinstance(thres, Sequence) or len(thres) != len(tmax):
                raise ValueError(
                    'The number of thresholds in "Setting.Thres" should be equal to the number of time points in '
                    '"Setting.Tmax". "Setting.Thres" refers to the lowest acceptable performance of the design '
                    "algorithms, the performance can be the solution quality."
                )
