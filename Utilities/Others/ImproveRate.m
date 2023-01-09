function improve = ImproveRate(Solution,improve,innerG,type)
% Calculate designed algorithms' performance improvement rate during k
% consecutive iterations.

k = 3;

% initialize
if innerG == 1
    improve = [1,zeros(1,k)]; % [fitness improve rate, fitness 1, 2, ..., k]
end

% keep the best fitness found at the latest k iteration in the improve vector
switch type
    case 'solution'
        Con      = sum(max(0,Solution.cons),2);
        Feasible = Con <= 0;
        Fitness  = Feasible.*Solution.objs + ~Feasible.*(Con+1e8);
    case 'algorithm'
        Fitness    = Solution.avePerformAll; 
end
improve    = [improve,min(Fitness)];
improve(2) = [];

% calculate the fitness improve rate
if innerG >= k
    tempRate = zeros(1,k-1);
    for i = 2:k
        tempRate(i-1) = (improve(i)-improve(i+1))./improve(i);
    end
    improve(1) = max(tempRate);
end
end