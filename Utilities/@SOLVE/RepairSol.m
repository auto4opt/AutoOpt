function decs = RepairSol(decs,Problem)
% Repair infeasible solutions.

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