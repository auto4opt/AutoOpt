function Lambda_alpha = constructLambda(alpha, D)
% constructLambda 构造给定alpha和维数D的对角矩阵Lambda_alpha
%
%   输入参数:
%       alpha: 标量，表示构造矩阵的基数
%       D    : 标量，表示目标对角矩阵的维度
%
%   输出参数:
%       Lambda_alpha: D x D的对角矩阵，其对角元素满足
%                     Lambda_alpha(i,i) = alpha^{(i-1)/(2*(D-1))}， i=1,...,D
%
%   示例:
%       Lambda_alpha = constructLambda(1e6, 10);

    if D < 2
        error('D must be greater than or equal to 2.');
    end

    indices = 0:(D-1);
    lambda = alpha.^(indices/(2*(D-1)));
    Lambda_alpha = diag(lambda);
end
