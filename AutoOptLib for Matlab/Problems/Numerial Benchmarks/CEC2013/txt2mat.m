% 设置数据文件夹路径
folderPath = 'data';  % 请根据实际路径调整

% 获取所有 .txt 文件
txtFiles = dir(fullfile(folderPath, '*.txt'));

% 遍历所有 .txt 文件并转换为 .mat 文件
for i = 1:length(txtFiles)
    % 获取文件的完整路径
    txtFilePath = fullfile(txtFiles(i).folder, txtFiles(i).name);
    
    % 读取 .txt 文件数据
    data = load(txtFilePath);
    
    % 创建保存的 .mat 文件路径
    [~, fileName, ~] = fileparts(txtFiles(i).name);  % 获取文件名（去除扩展名）
    matFilePath = fullfile(txtFiles(i).folder, [fileName, '.mat']);
    
    % 保存为 .mat 文件
    save(matFilePath, 'data');
    
    % 输出转换的信息
    disp(['已将 ', txtFiles(i).name, ' 转换为 ', [fileName, '.mat']]);
end
