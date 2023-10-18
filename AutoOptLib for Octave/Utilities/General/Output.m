function Output(varargin)
% Save and output the results.

%----------------------------Copyright-------------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.
%--------------------------------------------------------------------------

Setting = varargin{end};
switch Setting.Mode
    case 'design'
        algs          = varargin{1};
        algTrace      = varargin{2};
        instanceTrain = varargin{3};
        instanceTest  = varargin{4};
        if nargin > 5
            app = varargin{end-1};
        end

        % delete previous result files
        delete('Algs.xlsx','Algs_final_algs.csv','Algs_perf_final_algs.csv','Algs_best_algs_iter.csv','Algs_perf_best_algs_iter.csv')

        % save algorithm representations
        save('Algs.mat','algs', '-mat');

        %% save final algorithms' pseudocode
        allCode = cell(50,Setting.AlgN);
        for i = 1:Setting.AlgN
            code = pseudocode(algs(i),Setting,i,'final');
            if length(code) > size(allCode,1)
                error('Please set the number of rows of "allCode" larger than the length of "code".');
            end
            allCode(1:length(code),i) = code;
            % display on GUI
            if i == 1 && exist('app','var')
                app.TextArea2.Value = code;
            end
        end
        % writecell(allCode,'Algs.xlsx','Sheet','Final Algs');
        % writecell(allCode,'Algs_final_algs.csv');
        xlswrite('Algs.xlsx', allCode, 'Final Algs');
        cell2csv('Algs_final_algs.csv', allCode, ',');

        %% save final algorithms' performance
        firstRow = cell(1,Setting.AlgN+1);
        firstRow{1} = 'Instance Index';
        firstRow{2} = 'Best Algorithm';
        for i = 2:Setting.AlgN
            firstRow{i+1} = ['Algorithm ',num2str(i)];
        end
        sheetName = 'Perf of Final Algs';
        % writecell(firstRow,'Algs.xlsx','Sheet',sheetName,'Range','A1');
        xlswrite('Algs.xlsx', firstRow, sheetName, 'A1');

        firstColumn = repmat(instanceTest',1,Setting.AlgRuns);
        firstColumn = reshape(firstColumn',length(instanceTest)*Setting.AlgRuns,1);
        % writematrix(firstColumn,'Algs.xlsx','Sheet',sheetName,'Range','A2');
        xlswrite('Algs.xlsx', firstColumn, sheetName, 'A2');

        performance = zeros(numel(instanceTest)*Setting.AlgRuns,length(algs));
        for i = 1:length(algs)
            % reshape algorithm i's all performance values (each run on each instance) to a column vector
            performance(:,i) = reshape(algs(i).performance(numel(instanceTrain)+1:numel(instanceTrain)+numel(instanceTest),:)',size(performance,1),1);
        end

        % writematrix(performance,'Algs.xlsx','Sheet',sheetName,'Range','B2');
        xlswrite('Algs.xlsx', performance, sheetName, 'B2');
        % data = readcell('Algs.xlsx','Sheet',sheetName);
        [~, ~, data] = xlsread('Algs.xlsx',sheetName);
        % writecell(data,'Algs_perf_final_algs.csv');
        cell2csv('Algs_perf_final_algs.csv', data, ',');


        %% save pseudocode of the best algorithms found at each iteration of design
        allCode = cell(50,length(algTrace));
        for i = 1:length(algTrace)
            code = pseudocode(algTrace(i),Setting,i,'iterate');
            if length(code) > size(allCode,1)
                error('Please set the number of rows of "allCode" larger than the length of "code".');
            end
            allCode(1:length(code),i) = code;
        end
        % writecell(allCode,'Algs.xlsx','Sheet','Best Algs at Iteration');
        % writecell(allCode,'Algs_best_algs_iter.csv');
        xlswrite('Algs.xlsx', allCode,'Best Algs at Iteration');
        cell2csv('Algs_best_algs_iter.csv', allCode, ',');

        %% save performance of the best algorithms found at each iteration of design
        firstRow = cell(1,length(algTrace)+1);
        firstRow{1} = 'Instance Index';
        for i = 1:length(algTrace)
            firstRow{i+1} = ['Algorithm at Iteration ',num2str(i)];
        end
        sheetName = 'Perf of Best Algs at Iteration';
        % writecell(firstRow,'Algs.xlsx','Sheet',sheetName,'Range','A1');
        xlswrite('Algs.xlsx', firstRow,sheetName,'A1');

        firstColumn = repmat(instanceTrain',1,Setting.AlgRuns);
        firstColumn = reshape(firstColumn',length(instanceTrain)*Setting.AlgRuns,1);
        % writematrix(firstColumn,'Algs.xlsx','Sheet',sheetName,'Range','A2');
        xlswrite('Algs.xlsx', firstColumn, sheetName, 'A2');

        performTrace = zeros(numel(instanceTrain)*Setting.AlgRuns,length(algTrace));
        for i = 1:length(algTrace)
            % reshape algorithm i's all performance values (each run on each instance) to a column vector
            performTrace(:,i) = reshape(algTrace(i).performance(1:numel(instanceTrain),:)',size(performTrace,1),1);
        end
        % writematrix(performTrace,'Algs.xlsx','Sheet',sheetName,'Range','B2');
        xlswrite('Algs.xlsx', performTrace, sheetName, 'B2');
        % data = readcell('Algs.xlsx','Sheet',sheetName);
        [~, ~, data] = xlsread('Algs.xlsx',sheetName);
        % writecell(data,'Algs_perf_best_algs_iter.csv');
        cell2csv('Algs_perf_best_algs_iter.csv', data, ',');

        %% depict the convergence curve of the design process
        performTrace = mean(performTrace,1);
        curve = plot(1:length(performTrace),performTrace);
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
        delete('Solutions.xlsx','Solutions.csv','Fitness.csv','Constraint_iolation.csv','Fitness_all_runs.csv')

        % save original format of solutions
        save('Solutions.mat','bestSolutions','allSolutions', '-mat');

        %% save solutions and fitness obtained in the best run on each problem instance
        numIterations   = size(allSolutions,2);
        solutions       = cell(1+length(instance),2+numIterations);
        solutions{1,1}  = 'Instance Index';
        solutions{1,2}  = 'Best Solution';
        fitness         = cell(1+length(instance),2+numIterations);
        fitness{1,1}    = 'Instance Index';
        fitness{1,2}    = 'Best Solution';
        constraint      = cell(1+length(instance),2+numIterations);
        constraint{1,1} = 'Instance Index';
        constraint{1,2} = 'Best Solution';
        for i = 1:numIterations
            solutions{1,2+i}  = ['Iteration ',num2str(i)];
            fitness{1,2+i}    = ['Iteration ',num2str(i)];
            constraint{1,2+i} = ['Iteration ',num2str(i)];
        end
        for i = 1:length(instance)
            solutions{1+i,1}  = instance(i);
            fitness{1+i,1}    = instance(i);
            constraint{1+i,1} = instance(i);
            [~,best]          = min(SOLVE.fits(allSolutions(i,:)));
            solutions{1+i,2}  = num2str(allSolutions(i,best).dec);
            fitness{1+i,2}    = num2str(allSolutions(i,best).fit);
            constraint{1+i,2} = num2str(allSolutions(i,best).con);
            for j = 1:numIterations
                solutions{1+i,2+j}  = num2str(allSolutions(i,j).dec);
                fitness{1+i,2+j}    = num2str(allSolutions(i,j).fit);
                constraint{1+i,2+j} = num2str(allSolutions(i,j).con);
            end
            % depict convergence curves
            if i <= 20
                curve = plot(1:numIterations,SOLVE.fits(allSolutions(i,:)));
                xlabel('Iterations');
                ylabel('Fitness');
                title(['Convergence Curve of Instance ',num2str(instance(i))]);
                saveas(curve,['Convergence Curve of Instance ',num2str(instance(i))]);
                close;
            else
                warning('Only display the first 20 instances for saving computational resources.')
            end
        end
        % writecell(solutions,'Solutions.xlsx','Sheet','Solutions','Range','A1');
        xlswrite('Solutions.xlsx', solutions,'Solutions','A1');
        % writecell(solutions,'Solutions.csv');
        cell2csv('Solutions.csv', solutions, ',');
        % writecell(fitness,'Solutions.xlsx','Sheet','Fitness','Range','A1');
        xlswrite('Solutions.xlsx', fitness,'Fitness', 'A1');
        % writecell(fitness,'Fitness.csv');
        cell2csv('Fitness.csv', fitness,',');
        % writecell(constraint,'Solutions.xlsx','Sheet','Constraint Violation','Range','A1');
        xlswrite('Solutions.xlsx', constraint,'Constraint Violation','A1');
        % writecell(constraint,'Constraint_iolation.csv');
        cell2csv('Constraint_iolation.csv', constraint, ',');


        %% save fitness of all runs on each problem instance
        fitnessAll = cell(1+length(instance),3+Setting.AlgRuns);
        fitnessAll{1,1} = 'Instance Index';
        fitnessAll{1,2} = 'Mean';
        fitnessAll{1,3} = 'Std';
        for i = 1:Setting.AlgRuns
            fitnessAll{1,3+i} = ['Run ',num2str(i)];
        end
        for i = 1:length(instance)
            fitnessAll{1+i,1} = instance(i);
            fitnessAll{1+i,2} = num2str(mean(SOLVE.fits(bestSolutions(i,:))));
            fitnessAll{1+i,3} = num2str(std(SOLVE.fits(bestSolutions(i,:))));
            for j = 1:Setting.AlgRuns
                fitnessAll{1+i,3+j} = num2str(bestSolutions(i,j).fit);
            end
        end
        % writecell(fitnessAll,'Solutions.xlsx','Sheet','Fitness of All Runs','Range','A1');
        xlswrite('Solutions.xlsx', fitnessAll,'Fitness of All Runs','A1');
        % writecell(fitnessAll,'Fitness_all_runs.csv');
        cell2csv('Fitness_all_runs.csv', fitnessAll, ',');
end

end

function cell2csv(filename, cellArray, delimiter)
    % 打开要写入的CSV文件
    fileID = fopen(filename, 'w');

    % 遍历元胞数组，按行优先的顺序写入CSV文件
    [rows, cols] = size(cellArray);
    for row = 1:rows
        for col = 1:cols
            % 将元胞内容转换为字符串，并写入CSV文件
            cellValue = cellArray{row, col};
            if ischar(cellValue)
                fprintf(fileID, '"%s"', cellValue);
            else
                fprintf(fileID, '%s', num2str(cellValue));
            end
            if col < cols
                fprintf(fileID, delimiter);
            end
        end
        fprintf(fileID, '\n'); % 换行以分隔行
    end

    % 关闭文件
    fclose(fileID);
end


