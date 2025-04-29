function Tosz = computeTosz(Decs)
    [n, D] = size(Decs); 
    Tosz = zeros(size(Decs)); 

    for i = 1:n
        for j = 1:D
            x_i = Decs(i, j); 
            if x_i == 0
                Tosz(i, j) = 0;
                continue;
            end

            x_hat = log(abs(x_i));

            if x_i < 0
                sign_x = -1;
            else
                sign_x = 1;
            end

            if x_i > 0
                c1 = 10;
                c2 = 7.9;
            else
                c1 = 5.5;
                c2 = 3.1;
            end

            Tosz(i, j) = sign_x * exp(x_hat + 0.049 * (sin(c1 * x_hat) + sin(c2 * x_hat)));
        end
    end
end
