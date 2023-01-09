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
                solution(1,size(dec,1)) = SOLVE;
                dec = solution.RepairSol(dec,Problem);
                switch Problem.type{3}
                    case 'certain' 
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
                for i = 1:length(solution)
                    solution(i).dec = dec(i,:);
                    solution(i).obj = obj(i,:);
                    if ~isempty(con)
                        solution(i).con = con(i,:);
                    else
                        solution(i).con = 0;
                    end
                    Con = sum(max(0,solution(i).con));
                    feasible = Con<= 0;
                    solution(i).fit = feasible.*solution(i).obj+~feasible.*(Con+1e8);
                    if ~isempty(acc)
                        solution(i).acc = cell(1,length(acc));
                        for j = 1:length(acc)
                            if iscell(acc{j})
                                solution(i).acc{j} = acc{j}{i};
                            elseif ismatrix(acc{j})
                                solution(i).acc{j} = acc{j}(i,:);
                            else
                                error('Accessory data should be saved in cell or matrix format.');
                            end
                        end
                    end
                end
            end     
        end
        
        function value = decs(solution)
            value = cat(1,solution.dec);
        end
        
        function value = objs(solution)
            value = cat(1,solution.obj);
        end

        function value = cons(solution)
            value = cat(1,solution.con);
        end

        function value = fits(solution)
            value = cat(1,solution.fit);
        end
    end

    methods(Static)
        % ensure solutions to be feasible 
        dec = RepairSol(dec,Problem)

        % input algorithm
        [Alg,Setting] = InputAlg(Setting)

        % solve the targeted problem instances by the designed algorithm(s).
        [bestSolutions,allSolutions] = RunAlg(Alg,ProblemSolve,DataSolve,app,Setting)
    end
end