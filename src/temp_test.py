import numpy as np
from types import SimpleNamespace

from autooptlib.utils.design import Design
from autooptlib.utils.space import space

class DummyProblem:
    def __init__(self):
        self.bound = np.array([[-5.0, -5.0], [5.0, 5.0]])
        self.N = 6
        self.Gmax = 5
        self.type = ['continuous', 'static', 'certain']

    def evaluate(self, data, dec):
        obj = float(np.sum(dec ** 2))
        return obj, 0.0, None

problem = DummyProblem()
setting = space(
    problem,
    SimpleNamespace(metric='quality', alg_runs=1, archive=['archive_best'],
                    prob_n=1, tmax=np.inf, thres=-np.inf,
                    tune_para=False, alg_p=1, alg_q=3, alg_n=1),
)
setting.metric = 'quality'
setting.alg_runs = 1

design = Design(problem=[problem], setting=setting)
design.evaluate(problem, None, setting, [0])
print('Initial performance:', design.performance)

new_design, aux = design.get_new([problem], setting, inner_g=1, aux=None)
new_design.evaluate(problem, None, setting, [0])
print('New performance:', new_design.performance)
