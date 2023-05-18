function decs = RepairSol(decs,Problem)
% Repair infeasible solutions.

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

switch Problem.type{1}
    case 'continuous'
        Lower = Problem.bound(1,:);
        Upper = Problem.bound(2,:);
        decs  = max(min(decs,Upper),Lower); % limit solutions within decision space
    case 'discrete'
        if contains(Problem.setting,'dec_diff') % if elements of a solution should be different with respect to each other
            [N,D] = size(decs);
            for i = 1:N
                [~,ind] = unique(decs(i,:),'stable'); 
                while numel(ind) < D
                    DupInd = setdiff(1:D,ind);
                    for j = 1:numel(DupInd)
                        decs(i,DupInd(j)) = randperm(Problem.bound(2,DupInd(j)),1);
                    end
                    [~,ind] = unique(decs(i,:),'stable'); 
                end
            end
        end
    case 'permutation'
        % don't need to do anything       
end
end