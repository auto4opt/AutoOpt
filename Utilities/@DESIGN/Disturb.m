function [NewOp,NewPara,Aux]= Disturb(Algs,Setting,innerG,Aux)
% Design new algorithm(s) based on the current ones.

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

OpSpace   = Setting.OpSpace;
ParaSpace = Setting.ParaSpace;
ParaLocalSpace = Setting.ParaLocalSpace;

% get indices of non-empty parameters
indNonEmptyPara = find(~cellfun(@isempty,ParaSpace)==1);

NewOp   = cell(Setting.AlgN,Setting.AlgP);
NewPara = cell(Setting.AlgN,1);
for i = 1:Setting.AlgN
    thisOp   = Algs(i).operator;  % 1*P cells
    thisPara = Algs(i).parameter; % 1*1 cells

    % extract indices of operators and non-empty parameters
    indSearch = [];
    for j = 1:Setting.AlgP
        indSearch = [indSearch;thisOp{j}(2:end,1)];     
    end    
    indOp   = [thisOp{1}(1,1);indSearch;thisOp{1}(end,end)]; % indices of operators
    indPara = indNonEmptyPara(ismember(indNonEmptyPara,indOp)); % indices of parameters

    % determine where to disturb
    if innerG == 1
        if Setting.TunePara == false
            probOp   = 1/(numel(indOp)+numel(indPara)); % probability of disturbing an operator
            probPara = probOp*numel(indPara); % probability of disturbing all parameters
            prob     = [repmat(probOp,1,numel(indOp)),probPara];
            if numel(OpSpace(1,1):OpSpace(1,2)) == 1
                prob(1) = 0;
            end
            if numel(OpSpace(2,1):OpSpace(2,2)) == 1
                prob(2:numel(indOp)-1) = 0;
            end
            if numel(OpSpace(3,1):OpSpace(3,2)) == 1
                prob(numel(indOp)) = 0;
            end
            Aux{i}.seed = randsrc(1,1,[1:numel(indOp)+1;prob./sum(prob)]);
        else % always disturb parameters
            Aux{i}.seed = numel(indOp)+1;
        end
    end

    % disturb
    seed = Aux{i}.seed;
    if seed <= numel(indOp)
        % disturb an operator
        if seed == 1
            % disturb the choose operator
            indPool = OpSpace(1,1):OpSpace(1,2);
            indPool(indPool == indOp(seed)) = [];
            indNew = datasample(indPool,1);
            for j = 1:Setting.AlgP
                thisOp{j}(1,1) = indNew;
            end
        elseif seed == numel(indOp)
            % disturb the update operator
            indPool = OpSpace(3,1):OpSpace(3,2);
            indPool(indPool == indOp(seed)) = [];
            indNew = datasample(indPool,1);
            for j = 1:Setting.AlgP
                thisOp{j}(end,end) = indNew;
            end
        else
            % disturb the search operator
            if Setting.AlgP == 1
                % a single search pathway
                indPool = OpSpace(2,1):OpSpace(2,2);
                indPool(indPool==indOp(seed)) = [];
                indNew  = datasample(indPool,1);
                if numel(indSearch) == 1 && numel(indSearch) < Setting.AlgQ
                    samplePool = [indPool,+inf]; % +inf: add an operator after the current one
                elseif numel(indSearch) > 1 && numel(indSearch) < Setting.AlgQ
                    samplePool = [indPool,+inf,-inf]; % -inf: delete the current operator
                elseif numel(indSearch) > 1 && numel(indSearch) == Setting.AlgQ
                    samplePool = [indPool,-inf];
                else
                    samplePool = indPool;
                end
                sampleInd = datasample(samplePool,1);
                if sampleInd == +inf
                    thisOp{1} = [thisOp{1}(1:seed,:);zeros(1,2);thisOp{1}(seed+1:end,:)];
                    thisOp{1}(seed+1,:) = [indNew,thisOp{1}(seed,2)];
                    thisOp{1}(seed,2)   = indNew;
                elseif sampleInd == -inf
                    thisOp{1}(seed-1,2) = thisOp{1}(seed,2);
                    thisOp{1}(seed,:)   = [];
                else
                    thisOp{1}(seed,1)   = indNew;
                    thisOp{1}(seed-1,2) = indNew;
                end
            else
                % multiple search pathways             
                indPath = [];
                for j = 1:Setting.AlgP
                    numSearch = size(thisOp{j},1)-1; % number of search operators in pathway j
                    indPath = [indPath,repmat(j,1,numSearch)]; % indices of search pahways that search operators belong to
                end
                thisPath = indPath(seed-1);
                indChoose = thisOp{thisPath}(1,1);
                indUpdate = thisOp{thisPath}(end,end);
                currQ = randi(Setting.AlgQ); % number of search operators in the (seed-1)-th search pathway
                thisOp{thisPath} = zeros(currQ+1,2);
                indSearch = randi([OpSpace(2,1),OpSpace(2,2)]);
                thisOp{thisPath}(1,:) = [indChoose,indSearch];
                indSearchStart = indSearch;
                for k = 2:currQ
                    indSearchEnd = randi([OpSpace(2,1),OpSpace(2,2)]);
                    thisOp{thisPath}(k,:) = [indSearchStart,indSearchEnd];
                    indSearchStart = indSearchEnd;
                end      
                thisOp{thisPath}(end,:) = [indSearchStart,indUpdate];
            end
        end
    else
        % disturb all parameters
        currPara = [];
        currParaSpace = [];
        for j = 1:numel(indPara) % extract parameter values and parameter space
            for k = 1: numel(thisPara{indPara(j),1}) % for each parameter of operator indPara(j)
                currPara = [currPara,thisPara{indPara(j),1}(k)];
            end
            if strcmp(thisPara{indPara(j),2},'LS')
                currParaSpace = [currParaSpace;ParaLocalSpace{indPara(j)}];
            else
                currParaSpace = [currParaSpace;ParaSpace{indPara(j)}];
            end
        end
        if Setting.AlgN == 1
            [currPara,Aux{i}] = search_cma(currPara,currParaSpace',Aux{i},'algorithm');
        else
            [currPara,~] = search_mu_polynomial(currPara,currParaSpace','algorithm');
        end
        % insert new parameters to NewPara
        for j = 1:numel(indPara)
            for k = 1:numel(thisPara{indPara(j),1})
                thisPara{indPara(j),1}(k) = currPara(1);
                currPara(1) = [];
            end
        end
    end

    for j = 1:Setting.AlgP
        NewOp{i,j} = thisOp{j};
    end
    NewPara{i} = thisPara;
end