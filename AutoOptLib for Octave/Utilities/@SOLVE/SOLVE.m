% Class for solving the targeted problem.

%----------------------------Copyright-------------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 
%--------------------------------------------------------------------------

classdef SOLVE < handle
    properties(SetAccess = private)
        dec; % decision variables
        obj; % objective values
        con; % constraint violations
        fit; % solution's fitness (considering constraint violations)
        acc; % accessory data 
    end
    
    methods
        % construct solutions
        function solution = SOLVE(dec,Problem,Data)
            if nargin > 0
                solution = SOLVE;
                dec = solution.RepairSol(dec,Problem);
                switch Problem.type{3}
                    case 'certain' 
                        dec = feval(str2func(Problem.name),Data,dec,'repair'); % repair solutions       
                        [obj,con,acc] = feval(str2func(Problem.name),Data,dec,'evaluate'); % evaluate fitness                    
                    case 'uncertain' % for uncertain problems
                        obj = zeros(Problem.sampleN,1);
                        con = zeros(Problem.sampleN,1);
                        for i = 1:ProblemsampleN
                            [obj(i),con(i),acc] = feval(str2func(Problem.name),Data,dec,'evaluate');
                        end
                        if contains(Problem.setting,'uncertain_average')
                            obj = mean(obj,1);
                            con = mean(con,1);
                        elseif contains(Problem.setting,'uncertain_worst')
                            obj = max(obj,[],1);
                            con = max(con,[],1);
                        end 
                end

                solution.dec = dec;
                solution.obj = obj;

                if ~isempty(con)
                    solution.con = con;
                else
                    solution.con = 0;
                end

                Con = sum(max(0,solution.con));
                feasible = Con<= 0;
                solution.fit = feasible.*solution.obj+~feasible.*(Con+1e8);

                if ~isempty(acc)
                    solution.acc = cell(1,length(acc));
                    for j = 1:length(acc)
                        if iscell(acc{j})
                            solution.acc{j} = acc{j}{i};
                        elseif ismatrix(acc{j})
                            solution.acc{j} = acc{j}(i,:);
                        else
                            error('Accessory data should be saved in cell or matrix format.');
                        end
                    end
                end
                
            end     
        end

    end

    methods(Static)
        % ensure solutions to be feasible 
        dec = RepairSol(dec,Problem)

        % input algorithm
        [Alg,Setting] = InputAlg(Setting)

        % solve the targeted problem instances by the designed algorithm(s).
        [bestSolutions,allSolutions] = RunAlg(Alg,ProblemSolve,DataSolve,app,Setting)

        function value = decs(Solution)
            % value = cat(1,solution.dec);
            [~, N] = size(Solution);
            [~, M] = size(Solution(1).dec);
            value = zeros(N, M);
            for i = 1:N
                value(i,:) = Solution(i).dec;
            end
        end

        function value = fits(Solution)
            % value = cat(1,solution.fit);
            [~, N] = size(Solution);
            [~, M] = size(Solution(1).fit);
            value = zeros(N, M);
            for i = 1:N
                value(i,:) = Solution(i).fit;
            end
        end

        function value = objs(Solution)
            % value = cat(1,solution.obj);
            [~, N] = size(Solution);
            [~, M] = size(Solution(1).obj);
            value = zeros(N, M);
            for i = 1:N
                value(i,:) = Solution(i).obj;
            end
        end

        function value = cats(Solution1, Solution2)
            [~, N1] = size(Solution1);
            [~, N2] = size(Solution2);
            value = SOLVE;
            for i = 1:N1
                value(i) = Solution1(i);
            end
            for j = 1:N2
                value(j+N1) = Solution2(j);
            end
        end

        function value = cons(Solution)
            % value = cat(1,solution.con);
            [~, N] = size(Solution);
            [~, M] = size(Solution(1).con);
            value = zeros(N, M);
            for i = 1:N
                value(i,:) = Solution(i).con;
            end
        end

    end
end