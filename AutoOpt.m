function AutoOpt(varargin)
% The main function of AutoOpt.
% -------------------------------------------------------------------------
% Problem settings:
% P       : input problem name
% prob    : input problem properity file name (Stack.m/Place.m)
% instance: indexes of the targeting problem instances
% type    : {continuous\discrete\permutation, static\sequential, certain\uncertain}
% -------------------------------------------------------------------------
% Settings of the designed algorithm(s):
% Setting.Mode    : design/solve, i.e., the aim is to designing algorithm or solving problem
% Setting.AlgP    : number of searth pathways in a designed algorithm
% Setting.AlgQ    : maximum number of search operators in a search pathway
% Setting.Archive : name of the archive(s) that will be used in the designed algorithm(s)
% Setting.LSRange : range of parameter values that make the algorithm perform local search
% Setting.IncRate : minimum rate of solutions' fitness increase during 3 consecutive iterations
% Setting.InnerFE : maximum number of function evaluations for each call of local search
% -------------------------------------------------------------------------
% Settings of the design process:
% Setting.AlgN    : number of algorithms to be designed
% Setting.AlgRuns : number of algorithm runs on each problem instance
% Setting.ProbN   : population size of the designed algorithms for solving the targeted problem instances
% Setting.ProbFE  : number of fitness evaluations of the designed algorithms for solving the targeted problem instances
% Setting.Estimate: {quality/runtimeFE/runtimeSec/auc,
%                    average/statistic,
%                    default/approximate/intensification/racing}
%                   {algorithm's performance metric,
%                    method for evaluating designed algoritm's performance,
%                    method for comparing designed algoritm's performance}
% Setting.AlgFE   : number of function evaluations for designing algorithms
% Setting.Tmax    : maximum running time, time can be function evaluation or wall clock time (in second)
% Setting.Thres   : the lowest acceptable performance of the design algorithms, the performance can be the solution quality
% Setting.RacingK : number of instances evaluated before the first round of racing
% -------------------------------------------------------------------------
% Settings of solving the targeted problem:
% Setting.Alg     : algorithm file name, e.g., Algs
% -------------------------------------------------------------------------
% Example of running AutoOpt: 
%    AutoOpt()
%    AutoOpt('mode','design','problem','CEC2005_f1','instanceTrain',[1,2],'instanceTest',3)
%    AutoOpt('mode','solve','problem','CEC2005_f1','instanceSolve',[1,2],'AlgFile','Algs')

% cd(fileparts(mfilename('fullpath')));
% addpath(genpath(cd));

if nargin == 0 % call the GUI
    APP;
else 
    % get mode
    if any(strcmp(varargin,'mode'))
        Setting = struct;
        ind = find(strcmp(varargin,'mode'));
        Setting.Mode  = varargin{ind+1};
    else
        error('Please set the mode to "design" or "solve".');
    end

    % get problem
    if strcmp(Setting.Mode,'design')
        [prob,instanceTrain,instanceTest] = Input(varargin,Setting,'data');
    elseif strcmp(Setting.Mode,'solve')
        [prob,instanceSolve] = Input(varargin,Setting,'data');
    end

    % default parameters
    switch Setting.Mode
        case 'design'
            Setting.AlgP     = 1;
            Setting.AlgQ     = 2;
            Setting.Archive  = '';
            Setting.LSRange  = 0.2;
            Setting.IncRate  = 0.05;
            Setting.ProbN    = 50;
            Setting.ProbFE   = 1000;
            Setting.InnerFE  = 100;
            Setting.AlgN     = 10;
            Setting.AlgFE    = 100;
            Setting.AlgRuns  = 1;
            Setting.Metric   = 'quality'; % quality/runtimeFE/runtimeSec/auc
            Setting.Compare  = 'average'; % average/statistic
            Setting.Evaluate = 'default'; % default/approximate/intensification/racing
            Setting.Tmax     = [];
            Setting.Thres    = [];
            Setting.RacingK  = max(1,round(length(instanceTrain)*0.2));
            Setting.Surro    = Setting.ProbFE*0.2; % number of exact performance evaluations
            Setting.TunePara = false; % true/false
            Setting          = Input(varargin,Setting,'parameter'); % replace default parameters with user-defined ones
            Setting          = Input(Setting,'check'); % avoid conflicting parameter settings
            [algs,performance,performTrend] = Process(prob,instanceTrain,instanceTest,Setting);
            Output(algs,performance,performTrend,instanceTest,Setting);

        case 'solve'
            Setting.Mode = 'solve';
            Setting.AlgFile  = '';
            Setting.AlgName  = 'Continuous Genetic Algorithm';
            Setting.Metric   = 'quality';
            Setting.Tmax     = [];
            Setting.Thres    = [];
            Setting.ProbN    = 10;
            Setting.ProbFE   = 500;
            Setting.AlgRuns  = 2;
            Setting          = Input(varargin,Setting,'parameter');
            Setting          = Input(Setting,'check');
            [bestSolutions,allSolutions,~] = Process(prob,instanceSolve,Setting);
            Output(bestSolutions,allSolutions,instanceSolve,Setting);
    end
end
end