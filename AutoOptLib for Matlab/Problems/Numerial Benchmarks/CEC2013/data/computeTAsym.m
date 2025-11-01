function T_asym = computeTAsym(Decs, beta)
% T_asym transformation as defined in CEC2013 definitions:
%   if x_i > 0, y_i = x_i^(1 + beta * (i-1)/(D-1) * sqrt(x_i))
%   else       y_i = x_i

    [n, D] = size(Decs);
    T_asym = Decs;

    for i = 1:n
        for j = 1:D
            x = Decs(i, j);
            if x > 0
                expo = 1 + beta * (j - 1) / (D - 1) * sqrt(x);
                T_asym(i, j) = x ^ expo;
            end
        end
    end
end

