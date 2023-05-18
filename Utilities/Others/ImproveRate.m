function improve = ImproveRate(Solution,improve,innerG,type)
% Calculate the designed algorithm's performance improvement rate during k
% consecutive iterations.

%----------------------------Copyright-------------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 

% Please reference the paper below if using AutoOptLib in your publication:
% @article{zhao2023autooptlib,
%  title={AutoOptLib: A Library of Automatically Designing Metaheuristic 
%         Optimization Algorithms in Matlab},
%  author={Zhao, Qi and Yan, Bai and Hu, Taiwei and Chen, Xianglong and 
%          Yang, Jian and Shi, Yuhui},
%  journal={arXiv preprint 	arXiv:2303.06536},
%  year={2023}
% }
%--------------------------------------------------------------------------

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