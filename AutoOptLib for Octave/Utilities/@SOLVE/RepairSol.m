function decs = RepairSol(decs,Problem)
% Repair infeasible solutions.

%----------------------------Copyright-------------------------------------
% Copyright (C) <2025>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.
%--------------------------------------------------------------------------

switch Problem.type{1}
    case 'continuous'
        Lower = Problem.bound(1,:);
        Upper = Problem.bound(2,:);
        decs  = max(min(decs,Upper),Lower); % limit solutions within decision space
    case 'discrete'
        % if contains(Problem.setting,'dec_diff') % if elements of a solution should be different with respect to each other
        if strcmp(Problem.setting, 'dec_diff')
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
