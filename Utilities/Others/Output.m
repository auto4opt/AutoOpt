function Output(varargin)
% Save and output the results.

Setting = varargin{end};
switch Setting.Mode
    case 'design'
        algs         = varargin{1};
        performance  = varargin{2};
        performTrend = varargin{3};
        instance     = varargin{4};
        
        % delete previous files
        delete('Algs.xlsx')
        
        % save algorithm representations
        save('Algs.mat','algs');
        
        % save algorithm performance
        firstRow = cell(1,Setting.AlgN+1);
        firstRow{1} = 'Instance Index';
        for i = 1:Setting.AlgN
            firstRow{i+1} = ['Alg',num2str(i)];
        end
        writecell(firstRow,'Algs.xlsx','Sheet','Performance','Range','A1');
        firstColumn = repmat(instance',1,Setting.AlgRuns);
        firstColumn = reshape(firstColumn',length(instance)*Setting.AlgRuns,1);
        writematrix(firstColumn,'Algs.xlsx','Sheet','Performance','Range','A2');
        writematrix(performance,'Algs.xlsx','Sheet','Performance','Range','B2');
        
        % save algorithm profiles
        for i = 1:Setting.AlgN
            Alg = cell(Setting.AlgQ+5,Setting.AlgP*2);
            Alg{1,1}     = 'Components';
            Alg{2,1}     = algs(i).operatorPheno{1}.Choose;
            Alg{end-1,1} = algs(i).operatorPheno{1}.Update;
            Alg{end,1}   = algs(i).operatorPheno{1}.Archive;
            for j = 1:Setting.AlgP
                for k = 1:length(algs(i).operatorPheno{j}.Search)-1
                    Alg{k+2,j} = algs(i).operatorPheno{j}.Search{k};
                end
                Alg{Setting.AlgQ+3,j} = num2str(algs(i).operatorPheno{j}.Search{end});
            end
            Alg{1,Setting.AlgP+1}     = 'Parameters';
            Alg{2,Setting.AlgP+1}     = num2str(algs(i).parameterPheno{1}.Choose);
            Alg{end-1,Setting.AlgP+1} = num2str(algs(i).parameterPheno{1}.Update);
            for j = 1:Setting.AlgP
                for k = 1:length(algs(i).operatorPheno{j}.Search)-1
                    Alg{k+2,j+Setting.AlgP} = num2str(algs(i).parameterPheno{j}.Search{k});
                end
            end
            sheetName = ['Algorithm ',num2str(i)];
            writecell(Alg,'Algs.xlsx','Sheet',sheetName,'Range','A1');
        end
        
        % depict the convergence curve of the design process
        curve = plot(1:length(performTrend),performTrend);
        xlabel('Iterations');
        ylabel('Performance');
        title('Convergence Curve of the Design Process');
        saveas(curve,'Convergence Curve of the Design Process');
        close;

    case 'solve'
        bestSolutions = varargin{1};
        allSolutions  = varargin{2};
        instance      = varargin{3};
        
        % delete previous files
        delete('Solutions.xlsx')
        
        % save original format of solutions
        save('Solutions.mat','bestSolutions','allSolutions');
        
        % save solutions and fitness obtained in the best run on each problem instance
        numIterations  = size(allSolutions,2);
        solutions      = cell(1+length(instance),2+numIterations);
        solutions{1,1} = 'Instance Index';
        solutions{1,2} = 'Best';
        fitness        = cell(1+length(instance),2+numIterations);
        fitness{1,1}   = 'Instance Index';
        fitness{1,2}   = 'Best Solution';
        for i = 1:numIterations
            solutions{1,2+i} = ['Iteration ',num2str(i)];
            fitness{1,2+i}   = ['Iteration ',num2str(i)];
        end
        for i = 1:length(instance)
            solutions{1+i,1} = instance(i);
            fitness{1+i,1}   = instance(i);
            [~,best]         = min(allSolutions(i,:).fits);
            solutions{1+i,2} = num2str(allSolutions(i,best).dec);
            fitness{1+i,2}   = num2str(allSolutions(i,best).fit);
            for j = 1:numIterations
                solutions{1+i,2+j} = num2str(allSolutions(i,j).dec);
                fitness{1+i,2+j}   = num2str(allSolutions(i,j).fit);
            end
            % depict convergence curves
            if i <= 20
                curve = plot(1:numIterations,allSolutions(i,:).fits);
                xlabel('Iterations');
                ylabel('Fitness');
                title(['Convergence Curve of Instance ',num2str(instance(i))]);
                saveas(curve,['Convergence Curve of Instance ',num2str(instance(i))]);
                close;
            else
                warning('Only display the first 20 instances for saving computational resources.')
            end
        end
        writecell(solutions,'Solutions.xlsx','Sheet','Solutions','Range','A1');
        writecell(fitness,'Solutions.xlsx','Sheet','Fitness','Range','A1');
        
        % save fitness of all runs on each problem instance
        fitnessAll = cell(1+length(instance),3+Setting.AlgRuns);
        fitnessAll{1,1} = 'Instance Index';
        fitnessAll{1,2} = 'Mean';
        fitnessAll{1,3} = 'Std';
        for i = 1:Setting.AlgRuns
            fitnessAll{1,3+i} = ['Run ',num2str(i)];
        end
        for i = 1:length(instance)
            fitnessAll{1+i,1} = instance(i);
            fitnessAll{1+i,2} = num2str(mean(bestSolutions(i,:).fits));
            fitnessAll{1+i,3} = num2str(std(bestSolutions(i,:).fits));
            for j = 1:Setting.AlgRuns
                fitnessAll{1+i,3+j} = num2str(bestSolutions(i,j).fit);
            end
        end
        writecell(fitnessAll,'Solutions.xlsx','Sheet','Fitness of All Runs','Range','A1');
end