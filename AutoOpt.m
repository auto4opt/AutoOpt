function AutoOpt(varargin)
% --------------------------Introduction-----------------------------------
% AutoOptLib is a MATLAB library for automatically designing metaheuristic
% optimization algorithms. 
% AutoOptLib is licensed under the GNU General Public License v3.0 and is
% free for academic use.
% Please read the documentation at https://github.com/qz89/AutoOpt for a 
% step-by-step user guidance.
% Please reference the following paper if using AutoOptLib in your publication:
% @article{zhao2023autooptlib,
%  title={AutoOptLib: A Library of Automatically Designing Metaheuristic Optimization Algorithms in Matlab},
%  author={Zhao, Qi and Yan, Bai and Hu, Taiwei and Chen, Xianglong and Shi, Yuhui},
%  journal={arXiv preprint arXiv:2303.06536},
%  year={2023}
% }
% ----------------------------Settings-------------------------------------
% Settings of the targeted problem:
% Problem      : problem name
% InstanceTrain: indexes of training instances
% InstanceTest : indexes of test instances
% -------------------------------------------------------------------------
% Settings of the designed algorithm(s):
% Setting.Mode    : design/solve, i.e., the aim is to designing algorithm or solving problem
% Setting.AlgP    : number of search pathways in a designed algorithm
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
% Setting.Metric  : quality/runtimeFE/runtimeSec/auc, i.e., metric for evaluating algorithms' performance 
% Setting.Evaluate: exact/approximate/intensification/racing, i.e., method for evaluating algoritm's performance
% Setting.Compare : average/statistic, i.e., method for comparing the performance of algorithms
% Setting.AlgFE   : number of function evaluations for designing algorithms
% Setting.Tmax    : maximum running time, time can be function evaluations or wall clock time (in second)
% Setting.Thres   : the lowest acceptable performance of the design algorithms, the performance can be the solution quality
% Setting.RacingK : number of instances evaluated before the first round of racing
% Setting.Surro   : number of exact performance evaluations when using surrogate 
% -------------------------------------------------------------------------
% Settings of solving the targeted problem:
% Setting.Alg     : algorithm file name, e.g., Algs
% -------------------------------------------------------------------------
% Example of running AutoOpt: 
%    AutoOpt()
%    AutoOpt('Mode','design','Problem','CEC2005_f1','InstanceTrain',[1,2],'InstanceTest',3)
%    AutoOpt('Mode','design','Problem','CEC2005_f1','InstanceTrain',[1,3],'InstanceTest',2,'AlgN',2,'AlgFE',4,'AlgRuns',1,'Compare','average')
%    AutoOpt('Mode','solve','Problem','CEC2005_f1','InstanceSolve',[1,2],'AlgFile','Algs','ProbN',10,'ProbFE',100,'AlgRuns',2)
% -------------------------------------------------------------------------
cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));

if nargin == 0 % call the GUI
    APP;
else 
    % get mode
    if any(strcmp(varargin,'Mode'))
        Setting = struct;
        ind = find(strcmp(varargin,'Mode'));
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
            Setting.AlgQ     = 3;
            Setting.Archive  = '';
            Setting.IncRate  = 0.05;
            Setting.ProbN    = 20;
            Setting.ProbFE   = 2000;
            Setting.InnerFE  = 200;
            Setting.AlgN     = 5;
            Setting.AlgFE    = 2000;
            Setting.AlgRuns  = 5;
            Setting.Metric   = 'quality'; % quality/runtimeFE/runtimeSec/auc
            Setting.Evaluate = 'exact';   % exact/approximate/intensification/racing
            Setting.Compare  = 'statistic'; % average/statistic
            Setting.Tmax     = [];
            Setting.Thres    = [];
            Setting.LSRange  = 0.3;
            Setting.RacingK  = max(1,round(length(instanceTrain)*0.2)); 
            Setting.Surro    = Setting.ProbFE*0.3;
            Setting          = Input(varargin,Setting,'parameter'); % replace default parameters with user-defined ones
            Setting          = Input(Setting,'check'); % avoid conflicting parameter settings
            [algs,algTrace]  = Process(prob,instanceTrain,instanceTest,Setting);
            Output(algs,algTrace,instanceTrain,instanceTest,Setting);

        case 'solve'
            Setting.Mode = 'solve';
            Setting.AlgFile  = '';
            Setting.AlgName  = 'Continuous Genetic Algorithm';
            Setting.Metric   = 'quality';
            Setting.Tmax     = [];
            Setting.Thres    = [];
            Setting.ProbN    = 100;
            Setting.ProbFE   = 50000;
            Setting.AlgRuns  = 31;
            Setting          = Input(varargin,Setting,'parameter');
            Setting          = Input(Setting,'check');
            [bestSolutions,allSolutions] = Process(prob,instanceSolve,Setting);
            Output(bestSolutions,allSolutions,instanceSolve,Setting);
    end
end
end
