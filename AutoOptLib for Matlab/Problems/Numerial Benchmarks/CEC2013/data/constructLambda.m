function Lambda_alpha = constructLambda(alpha, D)
    if D < 2
        error('D must be greater than or equal to 2.');
    end

    indices = 0:(D-1);
    lambda = alpha.^(indices/(2*(D-1)));
    Lambda_alpha = diag(lambda);
end
