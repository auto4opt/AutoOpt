function [bestSolutions,allSolutions] = RunAlg(Alg,Problem,Data,app,Setting)
% Solve the targeted problem instances by a designed algorithm.

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

Operator  = Alg.operatorPheno;
Parameter = Alg.parameterPheno;
bestSolutions = SOLVE; % the best solution at the final iteration of each algorithm run
allSolutions  = SOLVE; % the best solution at each iteration of the best algorithm run (obtains the best solutions among all runs)
tic;
str = 'Initializing...';
if ~isempty(app)
    app.TextArea.Value = str;
    drawnow;
else
    bar = waitbar(0,str);
end
counter = 1;
for i = 1:length(Problem)
    currSolutions   = cell(Setting.AlgRuns,1); 
    fitnessSequence = zeros(Setting.AlgRuns,1);
    for j = 1:Setting.AlgRuns
        switch Problem(i).type{2}
            case 'static'
                [ArchSolution,t] = RunDesign(Operator,Parameter,Problem(i),Data(i),Setting);
                currSolutions{j} = ArchSolution(2:end); % best solutions at each iteration of the current algorithm run
                [~,best] = min(ArchSolution.fits);
                bestSolutions(i,j) = ArchSolution(best);
                
            case 'sequential'
                currProblem  = Problem(i);
                currData     = Data(i);
                currSolution = SOLVE;
                k = 1;
                while currData.continue == true
                    [ArchSolution,t] = RunDesign(Operator,Parameter,currProblem,currData,Setting);
                    [~,best] = min(ArchSolution.fits);
                    currSolution(k) = ArchSolution(best);
                    k = k+1;

                    % update the current instance of the problem sequence
                    [~,best] = min(ArchSolution.fits);
                    [currProblem,currData,~] = feval(str2func(Problem(i).name),currProblem,currData,ArchSolution(best),'sequence');
                end
                currSolutions{j} = currSolution; % final solutions of all subproblems obtained in the current algorithm run
                fitnessSequence(j) = sum(currSolution.fits);
        end        
        str = ['Solving... ',num2str(100*counter/(length(Problem)*Setting.AlgRuns)),'%'];
        if ~isempty(app)
            app.TextArea.Value = str;
            drawnow;
        else
            waitbar(counter/(length(Problem)*Setting.AlgRuns),bar,str);
        end
        counter = counter+1;
    end
    switch Problem(i).type{2}
        case 'static'
            [~,best] = min(bestSolutions(i,:).fits);
            allSolutions(i,1:length(currSolutions{1})) = currSolutions{best}; % solutions at each iteration of the best algorithm run
        case 'sequential'
            [~,best] = min(fitnessSequence);
            allSolutions(i,1:length(currSolutions{1})) = currSolutions{best}; % solutions at the final iteration of each subproblem of the best algorithm run         
    end
end
str = 'Complete';
if ~isempty(app)
    app.TextArea.Value = str;
    drawnow;
else
    waitbar(100,bar,str);
end
toc;
end

function [ArchSolution,t] = RunDesign(Operator,Para,Problem,Data,Setting)
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
Solution = SOLVE(dec,Problem,Data);

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
        for i = 1:size(Operator{1}.Search,1)
            improve = 1;
            innerG  = 1;
            while improve(1) >= Operator{1}.Search{i,end}(1) && innerG <= Operator{1}.Search{i,end}(2) % termination condition of search operator i
                % choose where to search from
                [ind,~] = feval(str2func(Operator{1}.Choose),Solution,Problem,Para{1}.Choose,Aux{1},G,innerG,Data,'execute');


                % search from the chosen solution(s)
                [New,Aux{1}] = feval(str2func(Operator{1}.Search{i,1}),Solution(ind),Problem,Para{1}.Search{i,1},Aux{1},G,innerG,Data,'execute');
                
                if ~isempty(Operator{1}.Search{i,2}) % if designed a sexual evolutonary algorithm with crossover and mutation
                    New = SOLVE.RepairSol(New,Problem);
                    [New,Aux{1}] = feval(str2func(Operator{1}.Search{i,2}),New,Problem,Para{1}.Search{i,2},Aux{1},G,innerG,Data,'execute');
                end
                New = SOLVE(New,Problem,Data);

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
            [New,Aux{i}] = feval(str2func(Operator{i}.Search{1}),Solution(currInd),Problem,Para{i}.Search{1},Aux{i},G,innerG,Data,'execute');
            New = SOLVE(New,Problem,Data);
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