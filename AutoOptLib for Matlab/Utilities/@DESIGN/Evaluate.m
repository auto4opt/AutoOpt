function NewAlgs = Evaluate(NewAlgs,Problem,Data,Setting,seedInstance)
% Evaluate the designed algorithm's performance.

%----------------------------Copyright-------------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 
%--------------------------------------------------------------------------

for i = 1:length(NewAlgs)
    Operator    = NewAlgs(i).operatorPheno;
    Parameter   = NewAlgs(i).parameterPheno;
    Performance = NewAlgs(i).performance;
    for j = 1:length(Problem(seedInstance))
        switch Problem(seedInstance(j)).type{2}
            case 'static'
                for k = 1:Setting.AlgRuns
                    [ArchSolution,t] = RunDesign(Operator,Parameter,Problem(seedInstance(j)),Data(seedInstance(j)),Setting);
                    switch Setting.Metric
                        case 'quality' % solution quality as performance metric
                            Performance(seedInstance(j),k) = min(ArchSolution.fits); % best performance 
                        case {'runtimeFE','runtimeSec'} % running time as performance metric
                            Performance(seedInstance(j),k) = t;
                        case 'auc' % AUC as performance metric
                            TimePoints = ceil(Setting.Tmax./Setting.ProbN); % time (Tmax) should be function evaluations when using auc
                            Performance(seedInstance(j),k) = 1/(sum(ArchSolution(TimePoints).fits <= Setting.Thres')/numel(TimePoints)+eps); % 1/AUC, the smaller the better performance
                    end
                end

            case 'sequential'
                for k = 1:Setting.AlgRuns
                    while Data(seedInstance(j)).continue == true
                        [ArchSolution,t] = RunDesign(Operator,Parameter,Problem(seedInstance(j)),Data(seedInstance(j)),Setting);
                        switch Setting.Metric
                            case 'quality'
                                Performance(seedInstance(j),k) = Performance(seedInstance(j),k) + min(ArchSolution.fits);
                            case {'runtimeFE','runtimeSec'}
                                Performance(seedInstance(j),k) = Performance(seedInstance(j),k) + t;
                            case 'auc'
                                TimePoints = ceil(Setting.Tmax./Setting.ProbN);
                                Performance(seedInstance(j),k) = Performance(seedInstance(j),k) + 1/sum(ArchSolution(TimePoints).fits <= Setting.QuaThres)/numel(TimePoints);
                        end
                        % update the current instance of the problem sequence
                        [~,best] = min(ArchSolution.fits);
                        [thisProblem,thisData,~] = feval(str2func(Problem(seedInstance(j)).name),Problem(seedInstance(j)),Data(seedInstance(j)),ArchSolution(best),'sequence');
                        Problem(seedInstance(j)) = thisProblem;
                        Data(seedInstance(j)) = thisData;
                    end
                end
        end
    end
    NewAlgs(i).performance = Performance;
end
end

function [ArchSolution,t] = RunDesign(Operator,Para,Problem,Data,Setting)
% Problem.bound: 2*D for continuous and discrete problems, 2*1 for permutation problem
% Problm.N     : number of solutions to the targeted problem instance
% Operator     : 1*P, P search pathways
% Para         : 1*P, P search pathways

% initialize solution(s)
switch Problem.type{1}
    case 'continuous'
        Lower = Problem.bound(1,:);
        Upper = Problem.bound(2,:);
        dec   = unifrnd(repmat(Lower,Problem.N,1),repmat(Upper,Problem.N,1));
    case 'discrete'
        D     = size(Problem.bound,2); % size of the solution space
        dec   = zeros(Problem.N,D);
        for j = 1:D
            dec(:,j) = randi([Problem.bound(1,j),Problem.bound(2,j)],Problem.N,1);
        end
    case 'permutation'
        [~,dec] = sort(rand(Problem.N,Problem.bound(end)),2);
end
Solution = SOLVE;
for i = 1:Problem.N
    Solution(i) = SOLVE(dec(i,:),Problem,Data);
end

% initialize archive(s)
Archive = cell(length(Operator{1}.Archive),1);
for i = 1:length(Operator{1}.Archive)
    [Archive{i},~] = feval(str2func(Operator{1}.Archive{i}),Solution,Problem,'execute');
end
ArchSolution = Solution(randi(Problem.N)); % the best solution found at each iteration, Gmax*1
Aux = cell(1,Setting.AlgP); 
for i = 1:Setting.AlgP
    Aux{i} = struct; % auxiliary structure array
end

% iterate
G = 1;
t = 0;
switch Setting.Metric
    case {'quality','auc'}
        Tmax  = +inf;
        Thres = -inf;
    case 'runtimeFE'
        Tmax  = ceil(Setting.Tmax./Setting.ProbN); % time changes to interation G
        Thres = Setting.Thres;
    case 'runtimeSec'
        Tmax  = Setting.Tmax;
        Thres = Setting.Thres;
end
if Setting.AlgP == 1 % if the designed algorithm has a single search pathway
    while G <= Problem.Gmax && t < Tmax && min(ArchSolution.fits) > Thres
        tic;
        for i = 1:size(Operator{1}.Search,1) % for each search operation
            improve = 1;
            innerG  = 1;
            while improve(1) >= Operator{1}.Search{i,end}(1) && innerG <= Operator{1}.Search{i,end}(2) % termination condition of search operator i
                % choose where to search from
                [ind,~] = feval(str2func(Operator{1}.Choose),Solution,Problem,Para{1}.Choose,Aux{1},G,innerG,Data,'execute');

                % search from the chosen solution(s)
                [new_dec,Aux{1}] = feval(str2func(Operator{1}.Search{i,1}),Solution(ind),Problem,Para{1}.Search{i,1},Aux{1},G,innerG,Data,'execute');
                if ~isempty(Operator{1}.Search{i,2}) % if designed a sexual evolutonary algorithm with crossover and mutation
                    new_dec = SOLVE.RepairSol(new_dec,Problem);
                    [new_dec,Aux{1}] = feval(str2func(Operator{1}.Search{i,2}),new_dec,Problem,Para{1}.Search{i,2},Aux{1},G,innerG,Data,'execute');
                end
                New = SOLVE;
                for j = 1:size(new_dec,1)
                    New(j) = SOLVE(new_dec(j,:),Problem,Data);
                end

                if strcmp(Operator{1}.Search{i,1},'search_pso') % update pbest and gbest for PSO' particle fly operator
                    Aux{1} = para_pso(New,Problem,Aux{1});
                elseif strcmp(Operator{1}.Search{i,1},'search_cma') % update parameters for CMA-ES
                    Aux{1} = para_cma(New,Problem,Aux{1},'solution');
                end
                
                % update solution(s)
                [Solution,~] = feval(str2func(Operator{1}.Update),[Solution,New],Problem,Para{1}.Update,Aux{1},G,innerG,Data,'execute');

                % update archive(s)
                for j = 1:length(Operator{1}.Archive)
                    [Archive{i},~] = feval(str2func(Operator{1}.Archive{j}),Solution,Archive{i},Problem,'execute');
                end
                [ArchSolution,~] = archive_best(Solution,ArchSolution,'execute');
                
                improve = ImproveRate(Solution,improve,innerG,'solution');
                innerG  = innerG+1;
                G       = G+1;
                if strcmp(Setting.Metric,'runtimeFE')
                    t = G;
                elseif strcmp(Setting.Metric,'runtimeSec')
                    t = t +toc;
                end
                if G > Problem.Gmax || t >= Tmax || min(ArchSolution.fits) <= Thres
                    break
                end
            end
            if G > Problem.Gmax || t >= Tmax || min(ArchSolution.fits) <= Thres
                break
            end
        end
    end

elseif Setting.AlgP > 1 % if the designed algorithm has multiple search pathways
    eachN  = round(Setting.ProbN*(1/Setting.AlgP));
    innerG = 1;
    while G <= Problem.Gmax && t < Tmax && min(ArchSolution.fits) > Thres
        tic;
        % choose where to search from
        [ind,~] = feval(str2func(Operator{1}.Choose),Solution,Problem,Para{1}.Choose,Aux{i},G,innerG,Data,'execute');

        % search from the chosen solution(s)
        allNew = [];
        for i = 1:Setting.AlgP
            % determine indices of solutions to be searched from
            if i == Setting.AlgP
                currInd = ind;
            else
                currInd = ind(1:eachN);
                ind(1:eachN) = []; % delete used indices
            end
            % search
            [new_dec,Aux{i}] = feval(str2func(Operator{i}.Search{1}),Solution(currInd),Problem,Para{i}.Search{1},Aux{i},G,innerG,Data,'execute');
            New = SOLVE;
            for j = 1:size(new_dec,1)
                New(j) = SOLVE(new_dec(i,:),Problem,Data);
            end
            if strcmp(Operator{i}.Search{1},'search_pso') % update pbest and gbest for PSO' particle fly operator
                Aux{i} = para_pso(New,Problem,Aux{i});
            elseif strcmp(Operator{i}.Search{1},'search_cma') % update parameters for CMA-ES
                Aux{i} = para_cma(New,Problem,Aux{i},'solution');
            end
            allNew = [allNew,New];
        end

        % update solution(s)
        [Solution,~] = feval(str2func(Operator{1}.Update),[Solution,allNew],Problem,Para{1}.Update,Aux{i},G,innerG,Data,'execute');

        % update archive(s)
        for i = 1:length(Operator{1}.Archive)
            [Archive{i},~] = feval(str2func(Operator{1}.Archive{i}),Solution,Archive{i},Problem,'execute');
        end
        [ArchSolution,~] = archive_best(Solution,ArchSolution,'execute');

        innerG = innerG+1;
        G      = G+1;
        if strcmp(Setting.Metric,'runtimeFE')
            t = G;
        elseif strcmp(Setting.Metric,'runtimeSec')
            t = t +toc;
        end
    end
end
end
