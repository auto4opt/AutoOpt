% Class for designing algorithms

%----------------------------Copyright-------------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 
%--------------------------------------------------------------------------

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
        %% initialize an algorithm
        function obj = DESIGN(varargin)
            if nargin > 0
                Problem = varargin{1};
                Setting = varargin{2};

                obj = DESIGN;
                [Operators,Paras]   = obj.Initialize(Setting,1);
                [Operators,Paras]   = obj.Repair(Operators,Paras,Problem,Setting);
                [OpPheno,ParaPheno] = obj.Decode(Operators,Paras,Problem,Setting);
                obj.operator          = Operators;
                obj.parameter         = Paras{1};
                obj.operatorPheno     = OpPheno;
                obj.parameterPheno    = ParaPheno;
                obj.performance       = zeros(length(Problem),Setting.AlgRuns);
                obj.performanceApprox = zeros(length(Problem),Setting.AlgRuns);
            end
        end

        %% design new algorithms based on the current ones
        function [objNew,Aux] = GetNew(obj,Problem,Setting,innerG,Aux)
            objNew = DESIGN;
            [NewOp,NewPara,Aux] = obj.Disturb(Setting,innerG,Aux);
            [Operators,Paras]   = objNew.Repair(NewOp,NewPara,Problem,Setting);
            [currOp,currPara]   = objNew.Decode(Operators,Paras,Problem,Setting);
            objNew.operator          = Operators;
            objNew.parameter         = Paras{1};
            objNew.operatorPheno     = currOp;
            objNew.parameterPheno    = currPara;
            objNew.performance       = zeros(length(Problem),Setting.AlgRuns);
            objNew.performanceApprox = zeros(length(Problem),Setting.AlgRuns);
        end

        %% get algorithm's all performance results
        function value = GetPerformance(obj,Setting,seedInstance)
            % reshape algorithm i's all performance values (each run on each instance) to a column vector
            if strcmp(Setting.Evaluate,'approximate') && sum(obj.performanceApprox,'all') ~= 0 && sum(obj.performance,'all') == 0
                value = reshape(obj.performanceApprox(seedInstance,:)',numel(seedInstance)*Setting.AlgRuns,1);
            else
                value = reshape(obj.performance(seedInstance,:)',numel(seedInstance)*Setting.AlgRuns,1);
            end
        end

        %% design new algorithms based on the current ones
        [NewOp,NewPara,Aux] = Disturb(obj,Problem,Setting,innerG,Aux)

        %% exactly evaluate algorithms' performance
        obj = Evaluate(obj,Problem,Data,Setting,indInstance)

        %% approximatly estimate algorithms' performance
        obj = Estimate(obj,Problem,Setting,indInstance,Surrogate)

        %% select algorithms
        % Algs = Select(obj,Problem,Data,Setting,indInstance)  % Select
        % should be a public method, in order to be compatible with Octave
        
        %%
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

        function value = cats(Algs1, Algs2)
            [~, N1] = size(Algs1);
            [~, N2] = size(Algs2);
            value = DESIGN;
            for i = 1:N1
                value(i) = Algs1(i);
            end
            for j = 1:N2
                value(j+N1) = Algs2(j);
            end
        end
    end
end