function T_asym = computeTAsym(Decs, beta)

    [n, D] = size(Decs); 
    T_asym = Decs;      

    for i = 1:n
        for j = 1:D
            if Decs(i, j) > 0
                T_asym(i, j) = Decs(i, j)^(1 + beta * (j - 1) / (D - 1)) * sqrt(Decs(i, j));
            end
        end
    end
end
