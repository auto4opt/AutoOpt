folderPath = 'data'; 

txtFiles = dir(fullfile(folderPath, '*.txt'));

for i = 1:length(txtFiles)
    txtFilePath = fullfile(txtFiles(i).folder, txtFiles(i).name);
    
    data = load(txtFilePath);

    [~, fileName, ~] = fileparts(txtFiles(i).name); 
    matFilePath = fullfile(txtFiles(i).folder, [fileName, '.mat']);
    
    save(matFilePath, 'data');
end
