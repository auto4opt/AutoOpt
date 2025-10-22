from autooptlib import autoopt

# 设计模式示例
final_algs, alg_trace = autoopt(
    Mode="design",
    Problem="cec2013_f1",
    InstanceTrain=[30, 30, 30],
    InstanceTest=[30],
    AlgN=5,
    AlgFE=400,
    ProbN=40,
    ProbFE=6000,
    Evaluate="racing",
    Compare="statistic",
    archive=["archive_best"],
)

# 求解模式示例
best_solutions, all_solutions = autoopt(
    Mode="solve",
    Problem="cec2013_f1",
    InstanceSolve=[30],
    AlgName="Continuous Random Search",
    AlgRuns=2,
    ProbN=40,
    ProbFE=6000,
)
