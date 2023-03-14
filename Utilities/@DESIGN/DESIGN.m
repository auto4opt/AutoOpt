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
        %% initialize the designed algorithms
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
                [Operators,Paras] = obj.Initialize(Setting,N);
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

        %% design new algorithms based on the current ones
        function [objNew,Aux] = GetNew(obj,Problem,Setting,innerG,Aux)
            [NewOp,NewPara,Aux] = obj.Disturb(Setting,innerG,Aux);
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

        %% get algorithms' average performance or statistically comparing results
        function value = GetPerformance(obj,Setting,seedInstance)
            allPerform = zeros(numel(seedInstance)*Setting.AlgRuns,length(obj));
            for i = 1:length(obj)
                % reshape algorithm i's all performance values (each run on each instance) to a column vector
                if strcmp(Setting.Evaluate,'approximate') && sum(obj(i).performanceApprox,'all') ~= 0 && sum(obj(i).performance,'all') == 0
                    allPerform(:,i) = reshape(obj(i).performanceApprox(seedInstance,:)',size(allPerform,1),1);
                else
                    allPerform(:,i) = reshape(obj(i).performance(seedInstance,:)',size(allPerform,1),1);
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

        %% design new algorithms based on the current ones
        [NewOp,NewPara,Aux] = Disturb(obj,Problem,Setting,innerG,Aux)

        %% exactly evaluate algorithms' performance
        obj = Evaluate(obj,Problem,Data,Setting,indInstance)

        %% approximatly estimate algorithms' performance
        obj = Estimate(obj,Problem,Setting,indInstance,Surrogate)

        %% select algorithms
        Algs = Select(obj,Problem,Data,Setting,indInstance)

        function obj = Construct(obj,operator,parameter)
            obj.operatorPheno  = operator;
            obj.parameterPheno = parameter;
        end
        
        function value = avePerformAll(obj)
            value = zeros(length(obj),1);
            for i = 1:length(obj)
                value(i) = mean(obj(i).performance,'all');
            end
        end

        function value = avePerformApproxAll(obj)
            value = zeros(length(obj),1);
            for i = 1:length(obj)
                value(i) = mean(obj(i).performanceApprox,'all');
            end
        end

        function value = avePerformPer(obj,ind)
            value = zeros(length(obj),1);
            for i = 1:length(obj)
                value(i) = mean(obj(i).performance(ind,:));
            end
        end

        function value = avePerformApproxPer(obj,ind)
            value = zeros(length(obj),1);
            for i = 1:length(obj)
                value(i) = mean(obj(i).performanceApprox(ind,:));
            end
        end
    end

    methods(Static)        
        % initialize graph representations of the designed algorithms
        [Operators,Paras] = Initialize(Problem,Setting,N)

        % repair the designed algorithms to ensure the algorithms' reasonability
        [Operators,Paras] = Repair(Operators,Paras,Problem,Setting)
        
        % decode the designed algorithms from their graph representations
        [currOp,currPara] = Decode(Operators,Paras,Problem,Setting) 
    end
end