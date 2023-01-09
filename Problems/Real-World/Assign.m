%% load data
filename = 'DataProcessed.xlsx';
[Data,~,~] = xlsread(filename);

%% rank products according to their future picking up frequencies, ranks are products' locations
ProdInd = unique(Data(:,1)); % all products
NumProd = length(ProdInd); % number of products
MaterInd = unique(Data(:,2)); % all materials
fr = zeros(NumProd,1); % future picking up frequency of p
for i = 1: NumProd
   fr(i) = sum(Data(Data(:,1)==ProdInd(i),6)); 
end
[~,loc] = sort(fr,'descend'); % locations of products, corresponding to products' indexes
temp = zeros(NumProd,1);
for i = 1:NumProd
    temp(i) = find(loc==i);
end
loc = temp;

%% assign each common material to a product, such that this common material's overall picking up distance is minimum
for i = 1:length(MaterInd)
    frMater = Data(Data(:,2)==MaterInd(i),6); % future picking up frequencies of material i 
    if length(frMater) > 1 % if common material
        ProdInd_ComMater = Data(Data(:,2)==MaterInd(i),1); % indexes of products that use the common material
        dis = sum(repmat(frMater,1,NumProd) .* abs(repmat(loc(ProdInd_ComMater),1,NumProd) - repmat(loc',length(frMater),1)),1); % 1:NumProd
        [~,best] = min(dis);
        Data(Data(:,2)==MaterInd(i),1) = ProdInd(best); % assign a new product for the common material
    end
end

%% delete repetitive material records
for i = 1:length(MaterInd)
    repMaterLine = find(Data(:,2)==MaterInd(i));
    Data(repMaterLine(2:end),:) = [];
end

%% save data
xlswrite('DataAssigned.xlsx',Data,'Sheet1');
xlswrite('DataAssigned.xlsx',loc,'Sheet2');