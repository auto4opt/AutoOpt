function Tosz = computeTosz(Decs)
    % 初始化矩阵
    [n, D] = size(Decs); % 获取样本数和维度
    Tosz = zeros(size(Decs)); % 用于存储结果

    % 遍历每个样本和维度
    for i = 1:n
        for j = 1:D
            x_i = Decs(i, j); % 当前元素
            if x_i == 0
                Tosz(i, j) = 0;
                continue;
            end

            x_hat = log(abs(x_i));

            % 计算 sign(x_i)
            if x_i < 0
                sign_x = -1;
            else
                sign_x = 1;
            end

            % 计算 c1 和 c2
            if x_i > 0
                c1 = 10;
                c2 = 7.9;
            else
                c1 = 5.5;
                c2 = 3.1;
            end

            % 计算 T_osz
            Tosz(i, j) = sign_x * exp(x_hat + 0.049 * (sin(c1 * x_hat) + sin(c2 * x_hat)));
        end
    end
end
