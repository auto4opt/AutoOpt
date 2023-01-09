function [output1,output2,output3] = Process(varargin)
% The main process of algorithm design and problem solving.

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
        % construct training problem properties
        ProblemTrain  = struct('name',[],'type',[],'bound',[],'setting',{''},'N',[],'Gmax',[]);
        instanceTrain = varargin{2};
        for i = 1:numel(instanceTrain)
            ProblemTrain(i).name = varargin{1};
            ProblemTrain(i).N    = Setting.ProbN;
            ProblemTrain(i).Gmax = ceil(Setting.ProbFE/Setting.ProbN);
        end
        [ProblemTrain,DataTrain,~] = feval(str2func(ProblemTrain(1).name),ProblemTrain,instanceTrain,'construct'); % infill problems' constraints and search boundary, construct data properties
        
        % construct test problem properties
        ProblemTest  = struct('name',[],'type',[],'bound',[],'setting',{''},'N',[],'Gmax',[]);
        instanceTest = varargin{3};
        for i = 1:numel(instanceTest)
            ProblemTest(i).name = varargin{1};
            ProblemTest(i).N    = Setting.ProbN;
            ProblemTest(i).Gmax = ceil(Setting.ProbFE/Setting.ProbN);
        end
        [ProblemTest,DataTest,~] = feval(str2func(ProblemTest(1).name),ProblemTest,instanceTest,'construct');
        
        % initialize algorithms and evaluate their performance
        AlgGmax = ceil(Setting.AlgFE/Setting.AlgN);
        indInstance = randperm(numel(instanceTrain));
        switch Setting.Evaluate
            case {'default','intensification'}
                Algs = DESIGN(ProblemTrain,Setting);
                Algs = Algs.Evaluate(ProblemTrain,DataTrain,Setting,indInstance);
            case 'racing'
                Algs = DESIGN(ProblemTrain,Setting);
                Algs = Algs.Evaluate(ProblemTrain,DataTrain,Setting,indInstance(1:Setting.RacingK));
            case 'approximate'
                Surrogate = Approximate(ProblemTrain,DataTrain,Setting,indInstance); % initialize surrogate model
                Algs = Surrogate.data(randperm(length(Surrogate.data),Setting.AlgN)); % get initial algorithms
        end
        
        % iterate
        G = 1;
        performTrend = zeros(AlgGmax,1);
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
                % design new algorithms
                [NewAlgs,Aux] = Algs.GetNew(ProblemTrain,Setting,innerG,Aux);
                
                % performance evaluation and algorithm selection
                switch Setting.Evaluate
                    case 'default'
                        NewAlgs = NewAlgs.Evaluate(ProblemTrain,DataTrain,Setting,indInstance);
                        Algs = Algs.Select(Algs,NewAlgs,ProblemTrain,DataTrain,Setting,indInstance);
                        
                    case 'approximate'
                        % get surrogate
                        NewAlgs = NewAlgs.Estimate(ProblemTrain,Setting,indInstance,Surrogate);
                        Algs = Algs.Select(Algs,NewAlgs,ProblemTrain,DataTrain,Setting,indInstance);
                        % update surrogate
                        if ismember(G,Surrogate.exactG)
                            NewAlgs = NewAlgs.Evaluate(ProblemTrain,DataTrain,Setting,indInstance);
                            Surrogate = Surrogate.UpdateModel(NewAlgs,Setting);
                        end
                        
                    case 'intensification'
                        % screen survivals from new algorithms
                        while ~isempty(NewAlgs) && ~isempty(indInstance)
                            NewAlgs = NewAlgs.Evaluate(ProblemTrain,DataTrain,Setting,indInstance(1));
                            NewAlgs = Algs.Select(Algs,NewAlgs,ProblemTrain,DataTrain,Setting,indInstance(1));
                            indInstance(1) = [];
                        end
                        % restore instance indices
                        indInstance = randperm(numel(instance));
                        % evaluate new incumbents (NewAlgs)' performance on all instances
                        for i = 1:length(NewAlgs)
                            for j = 1:numel(indInstance)
                                if sum(NewAlgs(i).performance(indInstance(j),:)) == 0 % if haven't evaluated on instance j
                                    NewAlgs(i) = NewAlgs(i).Evaluate(ProblemTrain,DataTrain,Setting,indInstance(j));
                                end
                            end
                        end
                        % update incumbent algorithms
                        Algs(randperm(Setting.AlgN,length(NewAlgs))) = NewAlgs;
                        
                    case 'racing'
                        % screen survivals from all algorithms (racing)
                        NewAlgs = NewAlgs.Evaluate(ProblemTrain,DataTrain,Setting,indInstance(1:Setting.RacingK));
                        Algs = Algs.Select(Algs,NewAlgs,ProblemTrain,DataTrain,Setting,indInstance(1:Setting.RacingK));
                        indInstance(1:Setting.RacingK) = [];
                        while length(Algs) > Setting.AlgN && ~isempty(indInstance)
                            Algs = Algs.Select(Algs,[],ProblemTrain,DataTrain,Setting,indInstance(1));
                            indInstance(1) = [];
                        end
                        % restore instance indices
                        indInstance = randperm(numel(instance));
                        % delete redundant algorithms after racing
                        if length(Algs) > Setting.AlgN
                            ind = randperm(length(Algs),length(Algs)-Setting.AlgN);
                            Algs(ind) = [];
                        end
                end
                
                % update auxiliary data
                for i = 1:Setting.AlgN
                    if isfield(Aux{i},'cma_Disturb') % if use CMA-ES
                        Aux{i} = para_cmaes(Algs(i),ProblemTrain,Aux{i},'algorithm'); % update CMA-ES's parameters
                    end
                end
                
                % record best algorithms at each iteration
                currCompare = Setting.Compare;
                Setting.Compare = 'average';
                currPerform = Algs.GetPerformance(Setting,indInstance);
                Setting.Compare = currCompare;
                [~,best] = min(currPerform);
                currPerform = currPerform(best);
                performTrend(G) = currPerform;
                          
                improve = ImproveRate(Algs,improve,innerG,'algorithm');
                innerG  = innerG+1;
                G       = G+1;
                if G > AlgGmax
                    break
                end
            end
        end
        
        % test the designed algorithm(s)
        str = 'Testing... ';
        if nargin > 4
            app.TextArea.Value = str;
            drawnow;
        else
            waitbar(100,bar,str);
        end
        Setting.Evaluate = 'default';
        indInstance = 1:numel(instanceTest);
        Algs = Algs.Evaluate(ProblemTest,DataTest,Setting,indInstance);
        Algs = Algs.Select(Algs,[],ProblemTest,DataTest,Setting,indInstance); % sort algorithms in descending order in terms of their performance
        allPerform = zeros(numel(indInstance)*Setting.AlgRuns,length(Algs));
        for i = 1:length(Algs)
            % reshape algorithm i's all performance values (each run on each instance) to a column vector
            allPerform(:,i) = reshape(Algs(i).performance(indInstance,:)',size(allPerform,1),1);
        end
        output1 = Algs;
        output2 = allPerform;
        output3 = performTrend;

        str = 'Complete';
        if nargin > 4
            app.TextArea.Value = str;
            drawnow;
        else
            waitbar(100,bar,str);
        end
        toc;
        
    case 'solve'
        % construct problem properties
        ProblemSolve  = struct('name',[],'type',[],'bound',[],'setting',{''},'N',[],'Gmax',[]);
        instanceSolve = varargin{2};
        for i = 1:numel(instanceSolve)
            ProblemSolve(i).name = varargin{1};
            ProblemSolve(i).N    = Setting.ProbN;
            ProblemSolve(i).Gmax = ceil(Setting.ProbFE/Setting.ProbN);
        end
        [ProblemSolve,DataSolve,~] = feval(str2func(ProblemSolve(1).name),ProblemSolve,instanceSolve,'construct'); % infill problems' constraints and search boundary, construct data properties
        
        % solve the problem
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
        output3 = [];
end
end