import numpy as np
from autooptlib.utils.design import Design
from autooptlib.utils.design._helpers import Pathway, PathwayParam, SearchStep, SearchParam

class DummyProblem:
    def __init__(self):
        self.bound = np.array([[-5.0, -5.0], [5.0, 5.0]])
        self.N = 6
        self.Gmax = 5
        self.type = ['continuous', 'static', 'certain']

    def evaluate(self, data, dec):
        obj = float(np.sum(dec ** 2))
        return obj, 0.0, None  # (obj, con, acc)

class DummySetting:
    def __init__(self):
        self.metric = 'quality'
        self.alg_runs = 1
        self.archive = ['archive_best']
        self.prob_n = 1
        self.tmax = np.inf
        self.thres = -np.inf

problem = DummyProblem()
setting = DummySetting()

step = SearchStep(primary='search_de_random',
                  secondary='search_mu_uniform',
                  termination=np.array([0.0, 3]))
params = SearchParam(primary=np.array([0.5, 0.8]),
                     secondary=np.array([0.1]))
pathway = Pathway(choose='choose_traverse',
                  search=[step],
                  update='update_greedy',
                  archive=['archive_best'])
path_params = PathwayParam(choose=None,
                           search=[params],
                           update=None)

design = Design()
design.operator_pheno = [[pathway]]
design.parameter_pheno = [[path_params]]
design.performance = np.zeros((1, setting.alg_runs))
design.performance_approx = np.zeros((1, setting.alg_runs))

design.evaluate(problem, None, setting, [0])
print('Performance:', design.performance)
