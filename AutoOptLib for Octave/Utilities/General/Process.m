function [output1,output2] = Process(varargin)
% The main process of algorithm design and problem solving.

%----------------------------Copyright-------------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.
%--------------------------------------------------------------------------

if strcmp(varargin{1}(end:-1:end-1),'m.')
    varargin{1}(end:-1:end-1)=[]; % delete '.m'
end
Setting = varargin{end};
switch Setting.Mode
    case 'design'
        tic;
        str = 'Initializing...';
        if nargin > 4
            app = varargin{end-1};
            app.TextArea.Value = str;
            drawnow;
        else
            bar = waitbar(0,str);
        end
        %% construct training problem properties
        Problem  = struct('name',[],'type',[],'bound',[],'setting',{''},'N',[],'Gmax',[]);
        instanceTrain = varargin{2};
        instanceTest  = varargin{3};
        seedTrain     = randperm(numel(instanceTrain));
        seedTest      = randperm(numel(instanceTest))+length(seedTrain);
        instance      = [instanceTrain,instanceTest];
        for i = 1:numel(instance)
            Problem(i).name    = varargin{1};
            Problem(i).setting = '';
            Problem(i).N       = Setting.ProbN;
            Problem(i).Gmax    = ceil(Setting.ProbFE/Setting.ProbN);
        end
        [Problem,Data,~] = feval(str2func(Problem(1).name),Problem,instance,'construct'); % infill problems' constraints and search boundary, construct data properties

        %% initialize algorithms and evaluate their performance
        Setting = Space(Problem,Setting); % get design space
        AlgGmax = ceil(Setting.AlgFE/Setting.AlgN);
        Algs = DESIGN;
        switch Setting.Evaluate
            case {'exact','intensification'}
                for i = 1:Setting.AlgN  % traverse each algorithm for being compatible with Octave
                    Algs(i) = DESIGN(Problem,Setting);
                    Algs(i) = Algs(i).Evaluate(Problem,Data,Setting,seedTrain);
                end
            case 'racing'
                for i = 1:Setting.AlgN
                    Algs(i) = DESIGN(Problem,Setting);
                    Algs(i) = Algs(i).Evaluate(Problem,Data,Setting,seedTrain(1:Setting.RacingK));
                end
            case 'approximate'
                Surrogate = Approximate(Problem,Data,Setting,seedTrain); % initialize surrogate model
                Algs = Surrogate.data(randperm(length(Surrogate.data),Setting.AlgN)); % get initial algorithms
        end

        %% iterate
        G = 1;
        AlgTrace = DESIGN; % for save best algorithms found at each iteration
        while G <= AlgGmax
            str = ['Designing... ',num2str(100*G/AlgGmax),'%'];
            if nargin > 4
                app.TextArea.Value = str;
                drawnow;
            else
                waitbar(G/AlgGmax,bar,str);
            end
            improve = 1;
            innerG  = 1;
            Aux     = cell(Setting.AlgN,1); % auxiliary data
            if Setting.AlgN == 1
                innerGmax = ceil(Setting.AlgFE/Setting.AlgN/10);
            else % global search for designing multiple algorithms
                innerGmax = 1;
            end
            while improve(1) >= Setting.IncRate && innerG <= innerGmax
                %% design new algorithms
                NewAlgs = DESIGN;
                for i = 1:length(Algs) % traverse each algorithm for being compatible with Octave
                    [NewAlgs(i),Aux{i}] = Algs(i).GetNew(Problem,Setting,innerG,Aux{i});
                end

                %% performance evaluation and algorithm selection
                switch Setting.Evaluate
                    case 'exact'
                        for i = 1:length(NewAlgs) % traverse each algorithm for being compatible with Octave
                            NewAlgs(i) = NewAlgs(i).Evaluate(Problem,Data,Setting,seedTrain);
                        end
                       % Algs = Select([Algs,NewAlgs],Problem,Data,Setting,seedTrain);
                       Algs = Select(DESIGN.cats(Algs,NewAlgs),Problem,Data,Setting,seedTrain);

                    case 'approximate'
                        % get surrogate
                        for i = 1:length(NewAlgs)
                            NewAlgs(i) = NewAlgs(i).Estimate(Problem,Setting,seedTrain,Surrogate);
                        end
                        Algs = Select([Algs,NewAlgs],Problem,Data,Setting,seedTrain);
                        % update surrogate
                        if ismember(G,Surrogate.exactG)
                            for i = 1:length(NewAlgs)
                                NewAlgs(i) = NewAlgs(i).Evaluate(Problem,Data,Setting,seedTrain);
                            end
                            Surrogate = Surrogate.UpdateModel(NewAlgs,Setting);
                        end

                    case 'intensification'
                        % screen survivals from new algorithms
                        while ~isempty(NewAlgs) && ~isempty(seedTrain)
                            for i = 1:length(NewAlgs)
                                NewAlgs(i) = NewAlgs(i).Evaluate(Problem,Data,Setting,seedTrain(1));
                            end
                            NewAlgs = Select([Algs,NewAlgs],Problem,Data,Setting,seedTrain(1));
                            seedTrain(1) = [];
                        end
                        % restore instance indices
                        seedTrain = randperm(numel(instanceTrain));
                        % evaluate new incumbents (NewAlgs)' performance on all instances
                        for i = 1:length(NewAlgs)
                            for j = 1:numel(seedTrain)
                                if sum(NewAlgs(i).performance(seedTrain(j),:)) == 0 % if haven't evaluated on instance j
                                    NewAlgs(i) = NewAlgs(i).Evaluate(Problem,Data,Setting,seedTrain(j));
                                end
                            end
                        end
                        % update incumbent algorithms
                        Algs(randperm(Setting.AlgN,length(NewAlgs))) = NewAlgs;

                    case 'racing'
                        % screen survivals from all algorithms (racing)
                        for i = 1:length(NewAlgs)
                            NewAlgs(i) = NewAlgs(i).Evaluate(Problem,Data,Setting,seedTrain(1:Setting.RacingK));
                        end
                        Algs = Select([Algs,NewAlgs],Problem,Data,Setting,seedTrain(1:Setting.RacingK));
                        seedTrain(1:Setting.RacingK) = [];
                        while length(Algs) > Setting.AlgN && ~isempty(seedTrain)
                            Algs = Select(Algs,Problem,Data,Setting,seedTrain(1));
                            seedTrain(1) = [];
                        end
                        % restore instance indices
                        seedTrain = randperm(numel(instanceTrain));
                        % delete redundant algorithms after racing
                        if length(Algs) > Setting.AlgN
                            ind = randperm(length(Algs),length(Algs)-Setting.AlgN);
                            Algs(ind) = [];
                        end
                end

                %% update auxiliary data
                for i = 1:Setting.AlgN
                    if isfield(Aux{i},'cma_Disturb') % if use CMA-ES
                        Aux{i} = para_cmaes(Algs(i),Problem,Aux{i},'algorithm'); % update CMA-ES's parameters
                    end
                end

                %% record best algorithms at each iteration
                allPerform = zeros(numel(seedTrain)*Setting.AlgRuns,length(Algs));
                for i = 1:length(Algs)
                    allPerform(:,i) = Algs(i).GetPerformance(Setting,seedTrain);
                end
                [~,best] = min(mean(allPerform,1));
                AlgTrace(G) = Algs(best);

                improve = ImproveRate(Algs,improve,innerG,'algorithm');
                innerG  = innerG+1;
                G       = G+1;
                if G > AlgGmax
                    break
                end
            end
        end

        %% test the designed algorithm(s)
        str = 'Testing... ';
        if nargin > 4
            app.TextArea.Value = str;
            drawnow;
        else
            % waitbar(100,bar,str);
            waitbar(1,bar,str);
        end
        Setting.Evaluate = 'exact';
        for i = 1:length(Algs)
            Algs(i) = Algs(i).Evaluate(Problem,Data,Setting,seedTest);
        end
        Algs = Select(Algs,Problem,Data,Setting,seedTest); % sort algorithms in descending order in terms of their performance
        output1 = Algs; % final algorithms
        output2 = AlgTrace; % best algorithms found at each iteration of design

        str = 'Complete';
        if nargin > 4
            app.TextArea.Value = str;
            drawnow;
        else
            % waitbar(100,bar,str);
            waitbar(1,bar,str);
        end
        toc;

    case 'solve'
        %% construct problem properties
        ProblemSolve  = struct('name',[],'type',[],'bound',[],'setting',{''},'N',[],'Gmax',[]);
        instanceSolve = varargin{2};
        for i = 1:numel(instanceSolve)
            ProblemSolve(i).name = varargin{1};
            ProblemSolve(i).setting = '';
            ProblemSolve(i).N    = Setting.ProbN;
            ProblemSolve(i).Gmax = ceil(Setting.ProbFE/Setting.ProbN);

        end
        [ProblemSolve,DataSolve,~] = feval(str2func(ProblemSolve(1).name),ProblemSolve,instanceSolve,'construct'); % infill problems' constraints and search boundary, construct data properties

        %% solve the problem
        Solution = SOLVE;
        [Alg,Setting] = Solution.InputAlg(Setting); % algorithm profile
        if nargin > 3
            app = varargin{end-1};
            [bestSolutions,allSolutions] = Solution.RunAlg(Alg,ProblemSolve,DataSolve,app,Setting);
        else
            [bestSolutions,allSolutions] = Solution.RunAlg(Alg,ProblemSolve,DataSolve,[],Setting);
        end
        output1 = bestSolutions; % the best solution at the final iteration of each algorithm run
        output2 = allSolutions;  % the best solution at each iteration of the best algorithm run
end
end
