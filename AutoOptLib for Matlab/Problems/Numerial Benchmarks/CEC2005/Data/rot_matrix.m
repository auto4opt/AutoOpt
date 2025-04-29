function M = rot_matrix(D,c)
A = normrnd(0,1,D,D);
P = cGram_Schmidt(A);
A = normrnd(0,1,D,D);
Q = cGram_Schmidt(A);
u = rand(1,D);
D = c.^((u-min(u))./(max(u)-min(u)));
D = diag(D);
M = P*D*Q;
end