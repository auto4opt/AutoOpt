function [opPheno,paraPheno] = Decode(Operators,Paras,Problem,Setting)
% Decode the designed algorithm from the graph representation.
% Operator: 1*P, P search pathways
% Para    : 1*P, P search pathways

%----------------------------Copyright-------------------------------------
% Copyright (C) <2025>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.
%--------------------------------------------------------------------------

AllOp     = Setting.AllOp;
rate      = Setting.IncRate;
innerGmax = ceil(Setting.InnerFE/Setting.ProbN);

switch Problem(1).type{1}
    case 'continuous'
        % indMu = find(contains(AllOp,'search_mu'));
        indMu = find(~cellfun('isempty', strfind(AllOp, 'search_mu')));
    case {'discrete','permutation'}
        % indMu = find(contains(AllOp,'search'));
        indMu = find(~cellfun('isempty', strfind(AllOp, 'search')));
end
% indCross = find(contains(AllOp,'cross'));
indCross = find(~cellfun('isempty', strfind(AllOp, 'cross')));

N         = size(Operators,1);
opPheno   = cell(N,Setting.AlgP);
paraPheno = cell(N,Setting.AlgP);
for i = 1:N % for each algorithm
    for j = 1:Setting.AlgP  % for each search pathway
        opPheno{i,j}.Choose   = AllOp{Operators{i,j}(1,1)};
        opPheno{i,j}.Search   = cell(size(Operators{i,j},1)-1,3); % number of search operators * 3
        opPheno{i,j}.Update   = AllOp{Operators{i,j}(end,end)};
        opPheno{i,j}.Archive  = Setting.Archive;

        paraPheno{i,j}.Choose = Paras{i}{strcmp(opPheno{i,j}.Choose,AllOp),1};
        paraPheno{i,j}.Search = cell(size(Operators{i,j},1)-1,2); % number of search operators * 2
        paraPheno{i,j}.Update = Paras{i}{strcmp(opPheno{i,j}.Update,AllOp),1};

        k = 2;
        while k <= size(Operators{i,j},1)
            % set search operators and their parameters
            thisSearchInd = Operators{i,j}(k,1);
            thisSearchIndNext = Operators{i,j}(k,2);
            opPheno{i,j}.Search{k-1,1}   = AllOp{thisSearchInd};
            paraPheno{i,j}.Search{k-1,1} = Paras{i}{thisSearchInd,1};
            if ismember(thisSearchInd,indCross) && ismember(thisSearchIndNext,indMu)
                % put mutation operator to the second column
                opPheno{i,j}.Search{k-1,2}   = AllOp{thisSearchIndNext};
                paraPheno{i,j}.Search{k-1,2} = Paras{i}{thisSearchIndNext,1};
                % set search operators' termination conditions
                if strcmp(Paras{i}{thisSearchInd,2},'GS') || strcmp(Paras{i}{thisSearchIndNext,2},'GS')
                    opPheno{i,j}.Search{k-1,3} = [-inf,1]; % global search operator terminates after 1 iteration
                else
                    opPheno{i,j}.Search{k-1,3} = [rate,innerGmax];
                end
                % jump to the row after the next row of the Operator matrix
                k = k+2;
            else
                % set search operators' termination conditions
                if strcmp(Paras{i}{thisSearchInd,2},'GS')
                    opPheno{i,j}.Search{k-1,3} = [-inf,1];
                else
                    opPheno{i,j}.Search{k-1,3} = [rate,innerGmax];
                end
                k = k+1;
            end
        end

        % delete empty rows
        rowDelete = cellfun(@isempty,opPheno{i,j}.Search(:,1));
        opPheno{i,j}.Search(rowDelete,:) = [];
        paraPheno{i,j}.Search(rowDelete,:) = [];
    end
end
end