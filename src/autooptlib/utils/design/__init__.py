# --------------------------------------------------------------------------
# Copyright (C) <2025>  <Swarm Intelligence Lab>
#
# This file is part of AutoOptLib_py.
# AutoOptLib_py is a free software. You can use, redistribute, and/or modify
# it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or any later
# version.
# --------------------------------------------------------------------------

import numpy as np

# 从同级目录下的其他模块（文件）中导入我们即将实现的函数
# 这些函数对应于您提供的 .m 文件
from ._initialize import initialize
from ._repair import repair
from ._decode import decode
from ._disturb import disturb
from ._evaluate import evaluate
from ._estimate import estimate


class Design:
    """
    该类负责设计算法。
    它对应于 MATLAB 版本中的 @DESIGN 类文件夹。
    """

    # 1. 构造函数 (对应 DESIGN.m 中的 DESIGN function)
    def __init__(self, problem=None, setting=None):
        """
        构造函数，用于初始化一个算法设计实例。
        如果提供了 problem 和 setting，则会创建一个全新的、随机的算法实例。
        如果不提供参数，则会创建一个空的 Design 对象。
        """
        # 初始化所有属性 (对应 MATLAB 的 properties)
        self.operator = None
        self.parameter = None
        self.operator_pheno = None
        self.parameter_pheno = None
        self.performance = np.array([])
        self.performance_approx = np.array([])

        # 如果传入了参数，则执行初始化逻辑 (对应 nargin > 0)
        if problem is not None and setting is not None:
            # 调用静态方法进行初始化
            # 注意：在Python中，静态方法可以通过类名直接调用
            operators, paras = Design.initialize(setting, 1)
            operators, paras = Design.repair(operators, paras, problem, setting)
            op_pheno, para_pheno = Design.decode(operators, paras, problem, setting)

            self.operator = operators
            # Paras 在MATLAB中是一个1x1的cell，里面是真正的参数cell array
            # 这里我们直接取其内容
            self.parameter = paras[0]
            self.operator_pheno = op_pheno
            self.parameter_pheno = para_pheno

            # 假设 problem 是一个列表或类似结构，可以len()
            # 假设 setting 是一个对象，可以访问属性
            self.performance = np.zeros((len(problem), setting.alg_runs))
            self.performance_approx = np.zeros((len(problem), setting.alg_runs))

    # 2. 将外部函数“挂载”为类的方法
    # ------------------------------------
    # 静态方法 (对应 MATLAB 的 methods(Static))
    # 这些方法不依赖于某个具体的实例(self)，因此声明为静态方法
    initialize = staticmethod(initialize)
    repair = staticmethod(repair)
    decode = staticmethod(decode)

    # 实例方法 (对应 MATLAB 的 methods)
    # 这些方法需要访问实例的属性(self)，例如 self.operator
    disturb = disturb
    evaluate = evaluate
    estimate = estimate

    # 3. 直接在类中实现的其他方法 (来自 DESIGN.m)
    # ---------------------------------------------
    def get_new(self, problem, setting, inner_g, aux):
        """
        根据当前算法设计新算法。对应 GetNew 方法。
        """
        # 创建一个新的、空的 Design 实例
        new_obj = Design()

        # 调用实例方法 disturb
        new_op, new_para, aux = self.disturb(setting, inner_g, aux)

        # 调用静态方法 repair 和 decode
        operators, paras = Design.repair(new_op, new_para, problem, setting)
        curr_op, curr_para = Design.decode(operators, paras, problem, setting)

        # 填充新实例的属性
        new_obj.operator = operators
        new_obj.parameter = paras[0]
        new_obj.operator_pheno = curr_op
        new_obj.parameter_pheno = curr_para
        new_obj.performance = np.zeros((len(problem), setting.alg_runs))
        new_obj.performance_approx = np.zeros((len(problem), setting.alg_runs))

        return new_obj, aux

    def get_performance(self, setting, seed_instance):
        """
        获取算法的所有性能结果。对应 GetPerformance 方法。
        """
        # 注意：Python使用0-based索引，并且切片不包含末尾元素
        # seed_instance 在Python中应为一个列表或numpy数组
        if setting.evaluate == 'approximate' and np.sum(self.performance_approx) != 0 and np.sum(self.performance) == 0:
            # reshape 的 'F' 参数表示按列优先（Fortran-style），以匹配MATLAB的行为
            return self.performance_approx[seed_instance, :].reshape(-1, 1, order='F')
        else:
            return self.performance[seed_instance, :].reshape(-1, 1, order='F')

    def construct(self, operator, parameter):
        """对应 Construct 方法"""
        self.operator_pheno = operator
        self.parameter_pheno = parameter

    # ... 此处可以继续添加 avePerformAll 等其他简短的方法 ...