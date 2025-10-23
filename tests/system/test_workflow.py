"""System-level tests for design and solve workflows."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from autooptlib import autoopt


@pytest.fixture()
def chdir_tmp(tmp_path):
    """Temporarily change working directory to an isolated tmp path."""
    prev = Path.cwd()
    os.chdir(tmp_path)
    try:
        yield tmp_path
    finally:
        os.chdir(prev)


def test_autoopt_design_creates_outputs(chdir_tmp):
    final_algs, alg_trace = autoopt(
        Mode="design",
        Problem="cec2013_f1",
        InstanceTrain=[10],
        InstanceTest=[10],
        AlgN=2,
        AlgFE=60,
        AlgRuns=2,
        ProbN=15,
        ProbFE=800,
        Evaluate="racing",
        Compare="statistic",
        RacingK=1,
        archive=["archive_best"],
        Seed=2024,
    )

    # basic structural checks
    assert final_algs, "No algorithms returned from design workflow."
    assert final_algs[0].performance.shape[0] == 2  # train + test instances
    assert final_algs[0].performance.shape[1] == 2
    assert alg_trace, "Algorithm trace should not be empty."

    # output files generated
    expected = [
        "Algs.pkl",
        "Algs_best_algs_iter.csv",
        "Algs_perf_final_algs.csv",
        "ConvergenceCurve.csv",
    ]
    for name in expected:
        assert Path(name).exists(), f"Missing expected design output {name}"


def test_autoopt_solve_generates_solution_logs(chdir_tmp):
    best_solutions, all_solutions = autoopt(
        Mode="solve",
        Problem="cec2013_f1",
        InstanceSolve=[10],
        AlgName="ICA",
        AlgRuns=1,
        ProbN=20,
        ProbFE=1000,
        Metric="quality",
        Seed=2025,
    )

    assert len(best_solutions) == 1
    assert len(best_solutions[0]) == 1
    assert len(all_solutions) == 1
    assert all(all_solutions[0]), "History for best run should not be empty."
    assert np.isfinite(best_solutions[0][0].fit)

    expected = [
        "Solutions.pkl",
        "Solutions.csv",
        "Fitness_all_runs.csv",
    ]
    for name in expected:
        assert Path(name).exists(), f"Missing expected solve output {name}"
