function [Objs,WasteL,NumRacks] = Stack(Data,Decs)
N = size(Decs,1); % population size
D = size(Decs,2); % deminsion of solution, i.e., number of materials
Objs = zeros(N,1);
WasteL = zeros(N,1); % totally waste lengthes of all racks used for the given product
NumRacks = cell(N,1); % number of eack kind of used racks

materSize = Data.orgData(:,3:5); % |materials|*3, all materials' sizes
materVol = prod(materSize,2).*Data.orgData(:,7); % |materials|*1, volumn of all boxes of each material (SKU boxes)
numRack = size(Data.rackSize,1);
materL = zeros(numRack,D); % 28*|materials|, length side of each material in each kind of racks
NL = zeros(numRack,D); % 28*|materials|, number of materials stacked along the length side of each kind of racks
for d = 1:D % for each material
    materL(Data.stackSt{d}(:,4)==1,d) = materSize(d,2);
    materL(Data.stackSt{d}(:,4)==2,d) = materSize(d,1);
    NL(:,d) = Data.stackSt{d}(:,1);
end
rackL = Data.rackSize(:,1); % 28*1, lengthes of all kinds of racks
rackVol = prod(Data.rackSize,2); % 28*1, volumns of all kinds of racks

for i = 1: N  
    solution = Decs(i,:); % index of rack being used by all materials of a given product
    rackInd = unique(solution); % indexes of kinds of racks being used
    C = zeros(length(rackInd),1); 
    Clength = zeros(length(rackInd),1); 
    NumRack =  zeros(length(rackInd),1); 
    for j = 1:length(rackInd)
        j_materL   = materL(rackInd(j),solution==rackInd(j)); % length side of a single box of each material (belong to the given product) that uses rack j
        j_NL       = NL(rackInd(j),solution==rackInd(j)); % number of materials stacked along the length side of  rack j
        j_materVol = materVol(solution==rackInd(j));
        j_rackL    = rackL(rackInd(j)); % length of rack j
        j_rackVol  = rackVol(rackInd(j)); % volumn of rack j
        
        if sum(j_materL.*j_NL) <= j_rackL
            C(j) = j_rackVol-sum(j_materVol);
            Clength(j) = j_rackL-sum(j_materL.*j_NL);
            NumRack(j) = NumRack(j)+1;
        else
            j_materSize = materSize(solution==rackInd(j),:); % size of single box of the materials that use rack j
            [~,ind] = sort(prod(j_materSize,2),'descend');
            k = 1;
            stackedMaterInd = [];
            currL = 0;
            currVol = 0;
            C(j) = 0;
            while ~isempty(ind) % while each material
                while currL < j_rackL && k <= numel(ind) % while each unstacked materials
                    currL = currL+j_materL(ind(k))*j_NL(ind(k));
                    currVol = currVol+j_materVol(ind(k));
                    stackedMaterInd = [stackedMaterInd,ind(k)]; 
                    if currL > j_rackL
                        currL = currL-sum(j_materL(ind(k))*j_NL(ind(k)));
                        currVol = currVol-j_materVol(ind(k));
                        stackedMaterInd(end) = [];
                    end
                    k = k+1;
                end

                C(j) = C(j)+(j_rackVol-currVol);
                Clength(j) = Clength(j)+(j_rackL-currL);
                NumRack(j) = NumRack(j)+1;
                for p = 1:length(stackedMaterInd)
                    ind(ind==stackedMaterInd(p)) = [];
                end
                k = 1;
                currL = 0;
                currVol = 0;
            end
        end
    end
    Objs(i) = sum(C);
    WasteL(i) = sum(Clength);
    NumRacks{i} = NumRack;
end
end