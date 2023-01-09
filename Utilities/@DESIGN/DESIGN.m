classdef DESIGN < handle
    properties (SetAccess = private)
        operator;
        parameter;
        operatorPheno;
        parameterPheno;
        performance;
        performanceApprox;
    end

    methods
        % initialize the designed algorithms
        function obj = DESIGN(varargin)
            if nargin > 0
                Problem = varargin{1};
                Setting = varargin{2};
                if nargin == 3
                    N = varargin{3};
                else
                    N = Setting.AlgN;
                end

                obj(1,N) = DESIGN;
                [Operators,Paras] = obj.Initialize(Problem,Setting,N);
                [Operators,Paras] = obj.Repair(Operators,Paras,Problem,Setting);
                for i = 1:N
                    obj(i).operator  = Operators(i,:);
                    obj(i).parameter = Paras{i};
                    [currOp,currPara]     = obj.Decode(Operators(i,:),Paras{i},Problem,Setting);
                    obj(i).operatorPheno  = currOp;
                    obj(i).parameterPheno = currPara;
                    obj(i).performance       = zeros(length(Problem),Setting.AlgRuns);
                    obj(i).performanceApprox = zeros(length(Problem),Setting.AlgRuns);
                end
            end
        end

        % design new algorithms based on the current ones
        function [objNew,Aux] = GetNew(obj,Problem,Setting,innerG,Aux)
            [NewOp,NewPara,Aux] = obj.Disturb(Problem,Setting,innerG,Aux);
            [Operators,Paras]   = obj.Repair(NewOp,NewPara,Problem,Setting);
            objNew(1,Setting.AlgN) = DESIGN;
            for i = 1:Setting.AlgN
                objNew(i).operator  = Operators(i,:);
                objNew(i).parameter = Paras{i};
                [currOp,currPara]     = objNew.Decode(Operators(i,:),Paras{i},Problem,Setting);
                objNew(i).operatorPheno  = currOp;
                objNew(i).parameterPheno = currPara;
                objNew(i).performance       = zeros(length(Problem),Setting.AlgRuns);
                objNew(i).performanceApprox = zeros(length(Problem),Setting.AlgRuns);
            end
        end

        % get algorithms' average performance or statistically comparing results
        function value = GetPerformance(Algs,Setting,indInstance)
            allPerform = zeros(numel(indInstance)*Setting.AlgRuns,length(Algs));
            for i = 1:length(Algs)
                % reshape algorithm i's all performance values (each run on each instance) to a column vector
                if strcmp(Setting.Evaluate,'approximate') && sum(Algs(i).performanceApprox,'all') ~= 0 && sum(Algs(i).performance,'all') == 0
                    allPerform(:,i) = reshape(Algs(i).performanceApprox(indInstance,:)',size(allPerform,1),1);
                else
                    allPerform(:,i) = reshape(Algs(i).performance(indInstance,:)',size(allPerform,1),1);
                end
            end
            switch Setting.Compare
                case 'average' 
                    value = mean(allPerform,1);
                case 'statistic' 
                    [~,~,stats] = friedman(allPerform,1,'off');
                    value = multcompare(stats,'Display','off');
            end
        end

        % design new algorithms based on the current ones
        [NewOp,NewPara,Aux] = Disturb(Algs,Problem,Setting,innerG,Aux)

        % exactly evaluate algorithms' performance
        Algs = Evaluate(Algs,Problem,Data,Setting,indInstance)

        % approximatly estimate algorithms' performance
        Algs = Estimate(Algs,Problem,Setting,indInstance,Surrogate)

        function Alg = Construct(Alg,operator,parameter)
            Alg.operatorPheno  = operator;
            Alg.parameterPheno = parameter;
        end
        
        function value = avePerformAll(Algs)
            value = zeros(length(Algs),1);
            for i = 1:length(Algs)
                value(i) = mean(Algs(i).performance,'all');
            end
        end

        function value = avePerformApproxAll(Algs)
            value = zeros(length(Algs),1);
            for i = 1:length(Algs)
                value(i) = mean(Algs(i).performanceApprox,'all');
            end
        end

        function value = avePerformPer(Algs,ind)
            value = zeros(length(Algs),1);
            for i = 1:length(Algs)
                value(i) = mean(Algs(i).performance(ind,:));
            end
        end

        function value = avePerformApproxPer(Algs,ind)
            value = zeros(length(Algs),1);
            for i = 1:length(Algs)
                value(i) = mean(Algs(i).performanceApprox(ind,:));
            end
        end
    end

    methods(Static)
        % main process
        

        % initialize graph representations of the designed algorithms
        [Operators,Paras] = Initialize(Problem,Setting,N)

        % repair the designed algorithms to ensure the algorithms' reasonability
        [Operators,Paras] = Repair(Operators,Paras,Problem,Setting)
        
        % decode the designed algorithms from their graph representations
        [currOp,currPara] = Decode(Operators,Paras,Problem,Setting) 

        % select algorithms
        Algs = Select(Algs,NewAlgs,Problem,Data,Setting,indInstance)
    end
end