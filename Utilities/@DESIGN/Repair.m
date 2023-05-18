function [Operators,Paras] = Repair(Operators,Paras,Problem,Setting)
% Ensure the designed algorithm(s) to be reasonable.

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

AllOp = Setting.AllOp;
ParaLocalSpace = Setting.ParaLocalSpace;
BehavSpace = Setting.BehavSpace;


switch Problem(1).type{1}
    case 'continuous'
        indMu = find(contains(AllOp,'search_mu'));
    case {'discrete','permutation'}
        indMu = find(contains(AllOp,'search'));
end
indCross = find(contains(AllOp,'cross'));

for i = 1:size(Operators,1)
    for j = 1:Setting.AlgP
        % crossover should be followed by mutation
        for k = 2:size(Operators{i,j},1)-1
            if ismember(Operators{i,j}(k,1),indCross) && ~ismember(Operators{i,j}(k,2),indMu)
                indThisMu = datasample(indMu,1);
                Operators{i,j}(k,2) = indThisMu;
                Operators{i,j}(k+1,1) = indThisMu;
            end
        end
        % crossover should not be the last search operator
        if ismember(Operators{i,j}(end,1),indCross)
            thisIndMu = datasample(indMu,1);
            Operators{i,j}(end,1)   = thisIndMu;
            Operators{i,j}(end-1,2) = thisIndMu;
        end

        % delete edges that the starting and end points are the same
        rowDelete = Operators{i,j}(:,1)==Operators{i,j}(:,2);
        Operators{i,j}(rowDelete,:) = [];

        % if search_pso is involved, it should be incorporated with choose_traverse and update_always
        indPSO = find(strcmp(AllOp,'search_pso'));
        if ismember(indPSO,Operators{i,j})
            indChoose = find(strcmp(AllOp,'choose_traverse'));
            indUpdate = find(strcmp(AllOp,'update_always'));
            % only use the first pathway's choose and update operators in algorithm performance evaluation
            Operators{i,1}(1,1) = indChoose; 
            Operators{i,1}(end,end) = indUpdate;
            Operators{i,j} = [indChoose,indPSO;indPSO,indUpdate];
        end

        % there should be at most one global search operator in a pathway
        indSearch = Operators{i,j}(2:end,1); % indices of search operators
        for k = 1:numel(indSearch)
            Paras{i}{indSearch(k),2} = 'LS'; % set all search operators' behaviors as local search        
        end
        for k = 1:numel(indSearch)
            if isempty(BehavSpace{indSearch(k)}{2,1})
                indSearch(k) = 0; 
            end
        end
        indSearch(indSearch==0) = []; % delete operators that only behave local search
        indGS = []; % indices of global search operators
        for k = 1:numel(indSearch)
            if isempty(BehavSpace{indSearch(k)}{1,1}) % if operator only behaves global search
                indGS = [indGS,indSearch(k)];
            elseif any(Paras{i}{indSearch(k),1} < ParaLocalSpace{indSearch(k)}(:,1)) || any(Paras{i}{indSearch(k),1} > ParaLocalSpace{indSearch(k)}(:,2))
                % if parameters are smaller than the lower bound or larger
                % than the upper bound of the parameter space for local
                % search
                indGS = [indGS,indSearch(k)];
            end
        end
        if numel(indGS) == 1
            Paras{i}{indGS,2} = 'GS'; % revise the only global search operator's behavior to global
        elseif numel(indGS) > 1 % repair the global search operator to behave local search
            indRetain = indGS(randperm(numel(indGS),1)); % index of the retained global search operator
            rowRetain = Operators{i,j}(:,1)==indRetain;
            if ismember(indRetain,indCross) && ismember(Operators{i,j}(rowRetain,2),indMu) && ismember(Operators{i,j}(rowRetain,2),indGS)
                % if the global search operator is a crossover and the
                % crossover is followed by a global mutation, then retain 
                % the mutation.
                indRetain = [indRetain,Operators{i,j}(rowRetain,2)];
            end
            for k = 1:numel(indRetain)
                Paras{i}{indRetain(k),2} = 'GS'; % change the retained global search operator's behavior to global
                indGS(indGS==indRetain(k)) = []; % delete the retained global search operator
            end

            for k = 1:numel(indGS)
                if ~isempty(ParaLocalSpace{indGS(k)})
                    % if can perform local search, initialize the parameter
                    % values from the parameter space for local search
                    lower = ParaLocalSpace{indGS(k)}(:,1);
                    upper = ParaLocalSpace{indGS(k)}(:,2);
                    Paras{i}{indGS(k),1} = lower+(upper-lower).*rand(numel(lower),1);
                else
                    % if cannot perform local search, delete the operator
                    indDel = find(Operators{i,j}(:,2) == indGS(k));
                    Operators{i,j}(indDel,2) = Operators{i,j}(indDel+1,2);
                    Operators{i,j}(indDel+1,:) = [];
                end
            end
        end
    end
end