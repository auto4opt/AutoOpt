% prepare
days = 7;
filename = 'DataOrg.xlsx';

% elminate null
[~,~,raw] = xlsread(filename);
raw(strcmp(raw(:,3),'ActiveX VT_ERROR: '),:) = []; % null in size
raw(strcmp(raw(:,4),'ActiveX VT_ERROR: '),:) = []; % null in size
raw(strcmp(raw(:,5),'ActiveX VT_ERROR: '),:) = []; % null in size
raw(strcmp(raw(:,6),'ActiveX VT_ERROR: '),:) = []; % null in frequency
raw(strcmp(raw(:,7),'ActiveX VT_ERROR: '),:) = []; % null in SKU

% elminate zeros
for i = 3:7 
    tempLoc = false(size(raw,1),1);
    for j = 1:size(raw,1)
        if raw{j,i}==0
            tempLoc(j) = true;
        end
    end
    raw(tempLoc,:) = [];
end

xlswrite('DataClearned.xlsx',raw,'Sheet1'); % for restore to original names of products and materials

raw(1,:) = [];
raw(:,8) = [];

% keep all cells being strings
for i = 1:size(raw,1)
    for j = 1:size(raw,2)
        raw{i,j} = num2str(raw{i,j}); 
    end
end

% give products numerical indexes
IndProd = raw(:,1);
temp = unique(IndProd);
NumIndProd = zeros(length(IndProd),1);
for i = 1:length(temp)
    NumIndProd(strcmp(IndProd,temp(i))) = i; % numerical indexes of products
end

% give materials numerical indexes
IndMater = raw(:,2);
temp = unique(IndMater);
NumIndMater = zeros(length(IndMater),1);
for i = 1:length(temp)
    NumIndMater(strcmp(IndMater,temp(i))) = i; % numerical indexes of materials
end

% load size and frequency data
NumSizeFr = str2double(raw(:,3:7)); % strings in multiple cells to doubles in a single matrix

data = [NumIndProd,NumIndMater,NumSizeFr]; % data with indexes of products and materials

% calculate SKU (boxes)
data(:,6) = data(:,6)./30.*days; % picking up frequencies
filename = 'Size.xlsx';
[~,~,rawSize] = xlsread(filename);
for i = 1:length(temp)
    NumMaterBox = rawSize{strcmp(rawSize(:,1),temp(i)),2}; % number of materials in one box
    data(strcmp(IndMater,temp(i)),7) = ceil(data(strcmp(IndMater,temp(i)),7)./NumMaterBox); % update SKU (boxes), common material has the same SKU value in each line of the material
end
writecell
xlswrite('DataProcessed.xlsx',data,'Sheet1');