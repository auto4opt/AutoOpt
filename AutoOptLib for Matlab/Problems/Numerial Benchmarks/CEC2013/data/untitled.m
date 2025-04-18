% 初始化
clc;
clear;

% 参数设置
numIterations = 20; % 总调用次数
problems = arrayfun(@(i) ['CEC2013_f', num2str(i)], 1:numIterations, 'UniformOutput', false); % 问题名称数组
AlgFile = 'Algs'; % 固定算法名称
currentFile = mfilename('fullpath'); % 获取当前文件的完整路径
baseDir = fileparts(currentFile); % 获取当前文件的目录; % 当前工作目录
disp(baseDir)
% 主循环
for i = 1:numIterations
    try
        % 创建文件夹名称
        folderName = fullfile(baseDir, ['f', num2str(i)]);
        
        % 创建文件夹
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
        
        % 切换到新文件夹
        cd(folderName);
        disp(['进入文件夹: ', folderName]);
        
        % 动态设置参数
        mode = 'solve'; % 固定模式
        problem = problems{i}; % 当前问题名称
        instanceSolve = 50; % 当前 'InstanceSolve' 值
        
        % 调用 AutoOpt 函数
        AutoOpt('Mode', mode, 'Problem', problem, 'InstanceSolve', instanceSolve, 'AlgFile', AlgFile);
        
        % 返回到初始工作目录
        cd(baseDir);
    catch ME
        disp(['第 ', num2str(i), ' 次调用 AutoOpt 出现错误: ', ME.message]);
        
        % 确保返回到初始工作目录
        cd(baseDir);
    end
end

% 提示完成
disp('所有调用完成！');
