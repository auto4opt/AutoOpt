function T_asym = computeTAsym(Decs, beta)
    % 计算 T_asym 的函数
    % 输入:
    %   Decs - n x D 矩阵，每行是一个样本，每列是一个维度
    %   beta - 参数 beta
    % 输出:
    %   T_asym - n x D 矩阵，对应的 T_asym 值

    [n, D] = size(Decs); % 获取样本数和维度
    T_asym = Decs;       % 初始化 T_asym

    % 遍历每个样本和维度
    for i = 1:n
        for j = 1:D
            if Decs(i, j) > 0
                % 应用公式
                T_asym(i, j) = Decs(i, j)^(1 + beta * (j - 1) / (D - 1)) * sqrt(Decs(i, j));
            end
        end
    end
end
