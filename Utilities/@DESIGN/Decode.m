function [currOp,currPara] = Decode(Operator,Para,Problem,Setting)
% Decode the designed algorithm from the graph representation.
% Operator: 1*P, P search pathways
% Para    : 1*P, P search pathways

[AllOp,~,~,~,~] = Space(Problem,Setting);
rate      = Setting.IncRate;
innerGmax = ceil(Setting.InnerFE/Setting.ProbN);

switch Problem(1).type{1}
    case 'continuous'
        indMu = find(contains(AllOp,'search_mu'));
    case {'discrete','permutation'}
        indMu = find(contains(AllOp,'search'));
end
indCross = find(contains(AllOp,'cross'));
currOp   = cell(1,Setting.AlgP);
currPara = cell(1,Setting.AlgP);

for i = 1:Setting.AlgP % for each search pathway
    currOp{i}.Choose   = AllOp{Operator{i}(1,1)};
    currOp{i}.Search   = cell(size(Operator{i},1)-1,3); % number of search operators * 3
    currOp{i}.Update   = AllOp{Operator{i}(end,end)};
    currOp{i}.Archive  = Setting.Archive;

    currPara{i}.Choose = Para{strcmp(currOp{i}.Choose,AllOp),1};
    currPara{i}.Search = cell(size(Operator{i},1)-1,2); % number of search operators * 2
    currPara{i}.Update = Para{strcmp(currOp{i}.Update,AllOp),1};
    
    j = 2;
    while j <= size(Operator{i},1)
        % set search operators and their parameters
        thisSearchInd = Operator{i}(j,1);
        thisSearchIndNext = Operator{i}(j,2);
        currOp{i}.Search{j-1,1}   = AllOp{thisSearchInd};
        currPara{i}.Search{j-1,1} = Para{thisSearchInd,1};
        if ismember(thisSearchInd,indCross) && ismember(thisSearchIndNext,indMu)
            % put mutation operator to the second column
            currOp{i}.Search{j-1,2}   = AllOp{thisSearchIndNext}; 
            currPara{i}.Search{j-1,2} = Para{thisSearchIndNext,1};
            % set search operators' termination conditions
            if strcmp(Para{thisSearchInd,2},'GS') || strcmp(Para{thisSearchIndNext,2},'GS')
                currOp{i}.Search{j-1,3} = [-inf,1]; % global search operator terminates after 1 iteration
            else
                currOp{i}.Search{j-1,3} = [rate,innerGmax];
            end
            % jump to the row after the next row of the Operator matrix
            j = j+2;
        else
            % set search operators' termination conditions
            if strcmp(Para{thisSearchInd,2},'GS') 
                currOp{i}.Search{j-1,3} = [-inf,1]; 
            else
                currOp{i}.Search{j-1,3} = [rate,innerGmax];
            end
            j = j+1;
        end
    end

    % delete empty rows
    rowDelete = cellfun(@isempty,currOp{i}.Search(:,1));
    currOp{i}.Search(rowDelete,:) = [];
    currPara{i}.Search(rowDelete,:) = [];
end
end