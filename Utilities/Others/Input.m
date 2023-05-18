function varargout = Input(varargin)
% Process the input of algorithm design.

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

switch varargin{end}
    case 'data'
        value   = varargin{1};
        Setting = varargin{2};
        % input problem
        if any(strcmp(value,'Problem'))
            ind = find(strcmp(value,'Problem'));
            prob  = value{ind+1};
            varargout{1} = prob;
        else
            error('Please set the targeted problem.');
        end
        if strcmp(Setting.Mode,'design') && any(strcmp(value,'InstanceTrain')) && any(strcmp(value,'InstanceTest'))
            ind = find(strcmp(value,'InstanceTrain'));
            instanceTrain = value{ind+1}; % indexes of problem instances for designing algorithm
            ind = find(strcmp(value,'InstanceTest'));
            instanceTest  = value{ind+1}; % indexes of problem instances for testing the designed algorithm
            varargout{2} = instanceTrain;
            varargout{3} = instanceTest;
        elseif strcmp(Setting.Mode,'solve') && any(strcmp(value,'InstanceSolve'))
            ind = find(strcmp(value,'InstanceSolve'));
            instanceSolve = value{ind+1}; % indexes of problem instances to be solved
            varargout{2} = instanceSolve;
        else
            error('Please set the targeted problem instance indexes.');
        end

    case 'parameter'
        value   = varargin{1};
        Setting = varargin{2};

        if any(strcmp(value,'AlgP'))
            ind = find(strcmp(value,'AlgP'));
            Setting.AlgP = value{ind+1};
        end
        if any(strcmp(value,'AlgQ'))
            ind = find(strcmp(value,'AlgQ'));
            Setting.AlgQ = value{ind+1};
        end
        if any(strcmp(value,'Archive'))
            ind = find(strcmp(value,'Archive'));
            Setting.Archive = value{ind+1};
        end
        if any(strcmp(value,'LSRange'))
            ind = find(strcmp(value,'LSRange'));
            Setting.LSRange = value{ind+1};
        end
        if any(strcmp(value,'IncRate'))
            ind = find(strcmp(value,'IncRate'));
            Setting.IncRate = value{ind+1};
        end
        if any(strcmp(value,'ProbN'))
            ind = find(strcmp(value,'ProbN'));
            Setting.ProbN = value{ind+1};
        end
        if any(strcmp(value,'ProbFE'))
            ind = find(strcmp(value,'ProbFE'));
            Setting.ProbFE = value{ind+1};
        end
        if any(strcmp(value,'InnerFE'))
            ind = find(strcmp(value,'InnerFE'));
            Setting.InnerFE = value{ind+1};
        end
        if any(strcmp(value,'AlgN'))
            ind = find(strcmp(value,'AlgN'));
            Setting.AlgN = value{ind+1};
        end
        if any(strcmp(value,'AlgFE'))
            ind = find(strcmp(value,'AlgFE'));
            Setting.AlgFE = value{ind+1};
        end
        if any(strcmp(value,'AlgRuns'))
            ind = find(strcmp(value,'AlgRuns'));
            Setting.AlgRuns = value{ind+1};
        end
        if any(strcmp(value,'Metric'))
            ind = find(strcmp(value,'Metric'));
            Setting.Metric = value{ind+1};
        end
        if any(strcmp(value,'Compare'))
            ind = find(strcmp(value,'Compare'));
            Setting.Compare = value{ind+1};
        end
        if any(strcmp(value,'Evaluate'))
            ind = find(strcmp(value,'Evaluate'));
            Setting.Evaluate = value{ind+1};
        end
        if any(strcmp(value,'Tmax'))
            ind = find(strcmp(value,'Tmax'));
            Setting.Tmax = value{ind+1};
        end
        if any(strcmp(value,'Thres'))
            ind = find(strcmp(value,'Thres'));
            Setting.Thres = value{ind+1};
        end
        if any(strcmp(value,'RacingK'))
            ind = find(strcmp(value,'RacingK'));
            Setting.RacingK = value{ind+1};
        end
        if any(strcmp(value,'Surro'))
            ind = find(strcmp(value,'Surro'));
            Setting.Surro = value{ind+1};
        end
        if any(strcmp(value,'AlgFile'))
            ind = find(strcmp(value,'AlgFile'));
            Setting.AlgFile = value{ind+1};
        end
        if any(strcmp(value,'AlgName'))
            ind = find(strcmp(value,'AlgName'));
            Setting.AlgName = value{ind+1};
        end
        varargout{1} = Setting;

    case 'check'
        Setting = varargin{1};
        switch Setting.Mode
            case 'design'
                % AlgP or AlgQ should be equal to 1
                if Setting.AlgP > 1 && Setting.AlgQ > 1
                    error(['Setting.AlgP or Setting.AlgQ should be equal to 1 - ' ...
                        'For algorithms with multiple search pathways (AlgP>1), ' ...
                        'each search pathway should have only one search operator (AlgQ=1). ' ...
                        'For algorithms with a single search pathway (AlgP=1), ' ...
                        'the search pathway can have multiple search operators (AlgQ>=1).'])
                end

                % set Tmax and Thres when using runtimeFE or runtimeSec as the performance metric
                if strcmp(Setting.Metric,'runtimeFE') && isempty(Setting.Tmax)
                    Setting.Tmax = Setting.ProbFE;
                end
                if strcmp(Setting.Metric,'runtimeFE') && isempty(Setting.Thres)
                    error(['Please set "Setting.Thres" as the lowest acceptable performance ' ...
                        'of the design algorithms, the performance can be the solution quality.'])
                end
                if strcmp(Setting.Metric,'runtimeSec') && isempty(Setting.Tmax)
                    error('Please set "Setting.Tmax" as the maximum runtime (seconds).')
                end
                if strcmp(Setting.Metric,'runtimeSec') && isempty(Setting.Thres)
                    error(['Please set "Setting.Thres" as the lowest acceptable performance ' ...
                        'of the design algorithms, the performance can be the solution quality.'])
                end

                %  time (Tmax) should be function evaluations when using auc
                if strcmp(Setting.Metric,'auc') && numel(Setting.Tmax) <= 1
                    error(['"Setting.Tmax" should contain multiple time points. The time ' ...
                        'points should the numbers of function evaluations spent during ' ...
                        'the alorithm execution.'])
                end

                % thresholds should be corresponding to time points when using auc
                if strcmp(Setting.Metric,'auc') && numel(Setting.Thres) ~= numel(Setting.Tmax)
                    error(['The number of thresholds in "Setting.Thres" should be equal to ' ...
                        'the number of time points in "Setting.Tmax". "Setting.Thres" ' ...
                        'refers to the lowest acceptable performance of the design' ...
                        ' algorithms, the performance can be the solution quality.'])
                end

                % the "racing" evaluation method should be used with the "statictic" algorithm comparing method
                if strcmp(Setting.Evaluate,'racing') && ~strcmp(Setting.Compare,'statistic')
                    error(['The "racing" evaluation method should be used with the ' ...
                        'algorithm comparing method of "statistic". '])
                end

                % should set RacingK when using the "racing" evaluation method
                if strcmp(Setting.Evaluate,'racing') && isempty(Setting.RacingK)
                    error(['Please set "Setting.K" as the number of instances evaluated ' ...
                        'before the first round of racing.'])
                end

                % should set Surro when using the "approximate" evaluation method
                if strcmp(Setting.Evaluate,'approximate') && isempty(Setting.Surro)
                    error(['Please set "Setting.Surro" as the number of exact performance evaluations' ...
                        ' when using surrogate.'])
                end

                % it is not necessary to use "statistic" algorithm comparing method when using the "approximate" evaluation method
                if strcmp(Setting.Evaluate,'approximate') && strcmp(Setting.Compare,'statistic')
                    error(['It is not necessary to use the "statistic" algorithm ' ...
                        'comparing method when using the "approximate" evaluation method.'])
                end

                % should run the design multiple times when using the "statistic" comparsion method
                if strcmp(Setting.Compare,'statistic') && Setting.AlgRuns == 1
                    error(['Please run the design multiple times (Setting.AlgRuns>1)' ...
                        ' when using the "statistic" comparsion method.'])
                end

                % better to have a large population size when involving the EDA operator
                if Setting.ProbN < 5 && Setting.AlgP > 1
                    warning(['It is better to have a large population size if ' ...
                        'involving the EDA operator'])
                end

                % AlgP cannot be very large
                if Setting.AlgQ > 4
                    warning(['AlgQ is recommended to be larger than 4 for ' ...
                        'discrete and permutation problems due to the lack' ...
                        ' of so many search operators'])
                end

                % better to have a large number of training intrances or have a large number of algorithm runs, in order to make the statistical test discriminative
                if strcmp(Setting.Compare,'statistic')
                    warning(['It is better to have a large number of training intrances ' ...
                        'or have a large number of algorithm runs (set AlgRun to a large number), ' ...
                        'in order to make the statistical test discriminative.']);
                end

            case 'solve'
                % should specify the file or name of the algorithm being used
                if isempty(Setting.AlgFile) && isempty(Setting.AlgName)
                    error(['Please specify an algorithm file in Setting.AlgFile or ' ...
                        'specify an algorithm name in Setting.AlgName.'])
                end

                % set Tmax and Thres when using runtimeFE or runtimeSec as the performance metric
                if strcmp(Setting.Metric,'runtimeFE') && isempty(Setting.Tmax)
                    Setting.Tmax = Setting.ProbFE;
                end
                if strcmp(Setting.Metric,'runtimeFE') && isempty(Setting.Thres)
                    error(['Please set "Setting.Thres" as the lowest acceptable performance ' ...
                        'of the design algorithms, the performance can be the solution quality.'])
                end
                if strcmp(Setting.Metric,'runtimeSec') && isempty(Setting.Tmax)
                    error('Please set "Setting.Tmax" as the maximum runtime (seconds).')
                end
                if strcmp(Setting.Metric,'runtimeSec') && isempty(Setting.Thres)
                    error(['Please set "Setting.Thres" as the lowest acceptable performance ' ...
                        'of the design algorithms, the performance can be the solution quality.'])
                end

                %  time (Tmax) should be function evaluations when using auc
                if strcmp(Setting.Metric,'auc') && numel(Setting.Tmax) <= 1
                    error(['"Setting.Tmax" should contain multiple time points. The time ' ...
                        'points should the numbers of function evaluations spent during ' ...
                        'the alorithm execution.'])
                end

                % thresholds should be corresponding to time points when using auc
                if strcmp(Setting.Metric,'auc') && numel(Setting.Thres) ~= numel(Setting.Tmax)
                    error(['The number of thresholds in "Setting.Thres" should be equal to ' ...
                        'the number of time points in "Setting.Tmax". "Setting.Thres" ' ...
                        'refers to the lowest acceptable performance of the design' ...
                        ' algorithms, the performance can be the solution quality.'])
                end
        end
        varargout{1} = Setting;
end
end