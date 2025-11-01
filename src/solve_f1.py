from types import SimpleNamespace
from autooptlib.utils.general.process import process
from autooptlib.problems.cec2013 import cec2013_f1

setting = SimpleNamespace(
    Mode="design",
    AlgN=6,             # 同时维护的候选算法数，适当增大保持多样性
    AlgFE=1000,          # 设计阶段的总评估次数（原来 50 太低）
    AlgRuns=3,          # 每个候选算法重复运行 2 次，结果更稳定
    ProbN=40,           # 问题侧的种群规模
    ProbFE=5000,        # 每次 run_design 的函数评估额度
    Metric="quality",
    Evaluate="racing",  # racing 会先用少量实例筛选，省掉大量全量评估
    Compare="statistic",
    RacingK=1,          # 每次先用 1 个实例（可视需要调大）
    IncRate=0.03,
    archive=["archive_best"],
)

final_algs, alg_trace = process(
    cec2013_f1,
    [30, 30, 30],   # 训练实例（可以保留 3 个 30 维）
    [30],           # 测试实例
    setting=setting,
)

print(f"设计完成的算法数: {len(final_algs)}")
print("测试集性能:\n", final_algs[0].performance)
