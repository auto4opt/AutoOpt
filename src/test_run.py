from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from autooptlib.utils.design import Design
from autooptlib.utils.design._helpers import Pathway, PathwayParam, SearchStep, SearchParam
from autooptlib.utils.general.process import process
from autooptlib.problems.cec2013 import cec2013_f1


def demo_design_operator() -> None:
    class DummyProblem:
        def __init__(self) -> None:
            self.bound = np.array([[-5.0, -5.0], [5.0, 5.0]])
            self.N = 6
            self.Gmax = 5
            self.type = ["continuous", "static", "certain"]

        def evaluate(self, data, dec):
            obj = float(np.sum(dec**2))
            return obj, 0.0, None

    class DummySetting:
        def __init__(self) -> None:
            self.metric = "quality"
            self.alg_runs = 1
            self.archive = ["archive_best"]
            self.prob_n = 1
            self.tmax = np.inf
            self.thres = -np.inf

    problem = DummyProblem()
    setting = DummySetting()

    step = SearchStep(
        primary="search_de_random",
        secondary="search_mu_uniform",
        termination=np.array([0.0, 3]),
    )
    params = SearchParam(primary=np.array([0.5, 0.8]), secondary=np.array([0.1]))
    pathway = Pathway(
        choose="choose_traverse",
        search=[step],
        update="update_greedy",
        archive=["archive_best"],
    )
    path_params = PathwayParam(choose=None, search=[params], update=None)

    design = Design()
    design.operator_pheno = [[pathway]]
    design.parameter_pheno = [[path_params]]
    design.performance = np.zeros((1, setting.alg_runs))
    design.performance_approx = np.zeros((1, setting.alg_runs))

    design.evaluate(problem, None, setting, [0])
    print("Design demo performance:", design.performance)


def demo_cec2013_solve() -> None:
    setting = SimpleNamespace(
        Mode="solve",
        AlgName="Continuous Random Search",
        AlgRuns=1,
        ProbN=30,
        ProbFE=5000,
        Metric="quality",
        Tmax=None,
        Thres=None,
    )
    best_solutions, _ = process(cec2013_f1, [50], setting=setting)
    best_fit = best_solutions[0][0].fit
    print("CEC2013 f1 best fitness (one run):", best_fit)


if __name__ == "__main__":
    demo_design_operator()
    demo_cec2013_solve()
