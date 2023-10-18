function [Operators,Paras] = Initialize(Setting,N)
% Initialize the designed algorithm(s).

%----------------------------Copyright-------------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 
%--------------------------------------------------------------------------

OpSpace   = Setting.OpSpace;
ParaSpace = Setting.ParaSpace;
Operators = cell(N,Setting.AlgP); % N algorithms, each with AlgP search 
% pathways. Each cell contains a number of AlgP matrix within edges of the graph 
% representation of a designed ealgorithm
Paras     = cell(N,1); % parameters of N algorithms
for i = 1:N
    % initialize operators
    indChoose = randi([OpSpace(1,1),OpSpace(1,2)]); % different search pathways use the same choose operator
    indUpdate = randi([OpSpace(3,1),OpSpace(3,2)]); % different search pathways use the same update operator
    for j = 1:Setting.AlgP
        currQ = randi(Setting.AlgQ); % number of search operators in the current search pathway
        Operators{i,j} = zeros(currQ+1,2);

        indSearch = randi([OpSpace(2,1),OpSpace(2,2)]);
        Operators{i,j}(1,:) = [indChoose,indSearch];

        indSearchStart = indSearch;
        for k = 2:currQ
            indSearchEnd = randi([OpSpace(2,1),OpSpace(2,2)]);
            Operators{i,j}(k,:) = [indSearchStart,indSearchEnd];
            indSearchStart = indSearchEnd;
        end

        Operators{i,j}(end,:) = [indSearchStart,indUpdate];
    end
    % initialize parameters
    tempPara = cell(length(ParaSpace),2); % first column of each row
    % contains a cloumn vector of an operator's parameter values,
    % second column will contain a string of whether the operator
    % performs local or global search behavior.
    indNonEmptyPara = find(~cellfun(@isempty,ParaSpace)==1);
    for j = 1:numel(indNonEmptyPara)
        k = indNonEmptyPara(j);
        tempPara{k,1} = ParaSpace{k}(:,1)+(ParaSpace{k}(:,2)-ParaSpace{k}(:,1)).*rand(size(ParaSpace{k},1),1);
    end
    Paras{i} = tempPara;
end
end