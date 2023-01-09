function [output1,output2,output3] = Stack(varargin)
switch varargin{end}
    case 'construct' % % construct problem and data
        type     = {'discrete','static','certain'};
        Problem  = varargin{1};
        instance = varargin{2};

        orgData = readmatrix('DataAssigned.xlsx','Sheet',1);
        prodInd = unique(orgData(:,1)); % all products
        numProd = length(prodInd);
        Data = struct('orgData',[],'rackSize',[],'stackSt',[],'inFesMater',[]);

        lowerH = 20;
        upperH = 150;
        intervalH = 10;
        H = (lowerH:intervalH:upperH)';
        numH = length(H);
        rackSize = [repmat(180,numH,1),repmat(60,numH,1),H;repmat(240,numH,1),repmat(120,numH,1),H;]; % sizes of all kinds of racks, [long,width,height]
        
        for i = 1:numProd
            Problem(i) = Problem(1); % a number of numProd problem instances
            Problem(i).type = type;
            Data(i).orgData = orgData(orgData(:,1)==prodInd(i),:); % original data about the materials of product i
            
            % materials' stacking strategies when using each kind of racks
            boxSize = Data(i).orgData(:,3:5);  % size of boxes of all materials of product i, [long,width,height]
            numMater = size(Data(i).orgData,1); % number of materials
            stackSt = cell(numMater,1);
            cons = cell(numMater,1);
            inFesMater = false(numMater,1);
            for j = 1:numMater
                NH = floor(rackSize(:,3)./boxSize(j,3));
                NW = zeros(size(rackSize,1),1);
                NL = zeros(size(rackSize,1),1);
                control = zeros(size(rackSize,1),1); % control materials' stacking directions, 1: length side of material stacks in width side of rack, 2: width side of material stacks in width side of rack
                for k = 1:size(rackSize,1)
                    if rackSize(k,2) < min(boxSize(j,1),boxSize(j,2))
                        NW(k) = 0;
                    elseif rackSize(k,2) >= max(boxSize(j,1),boxSize(j,2))
                        if rem(rackSize(k,2),boxSize(j,1)) < rem(rackSize(k,2),boxSize(j,2)) % length side of material stacks in width side of rack
                            NW(k) = floor(rackSize(k,2)./boxSize(j,1));
                            control(k) = 1;
                        else
                            NW(k) = floor(rackSize(k,2)./boxSize(j,2));
                            control(k) = 2;
                        end
                    elseif rackSize(k,2) < boxSize(j,1) && rackSize(k,2) >= boxSize(j,2) % width side of material stacks in width side of rack
                        NW(k) = floor(rackSize(k,2)./boxSize(j,2));
                        control(k) = 2;
                    elseif rackSize(k,2) >= boxSize(j,1) && rackSize(k,2) < boxSize(j,2) % length side of material stacks in width side of rack
                        NW(k) = floor(rackSize(k,2)./boxSize(j,1));
                        control(k) = 1;
                    end
                    if NW(k) == 0 || NH(k) == 0
                        NL(k) = 0;
                    else
                        NL(k) = ceil(Data(i).orgData(j,7)./NW(k)./NH(k));
                        if control(k) == 1
                            total_L =  boxSize(j,2)*NL(k);
                        elseif control(k) == 2
                            total_L =  boxSize(j,1)*NL(k);
                        end
                        if total_L > rackSize(k,1)
                            NL(k) = 0;
                        end
                    end
                end
                stackSt{j} = [NL,NW,NH,control]; % 28*4
                cons{j} = NL.*NW.*NH==0; % 28*1 logical
                if prod(sum(stackSt{j}(:,1:3))) == 0
                    inFesMater(j) = true; % indexes of materials that cannot be stacked on all kinds of racks
                end
            end
            Data(i).rackSize = rackSize;
            Data(i).stackSt = stackSt; % product i's materials' stacking strategies
            Data(i).inFesMater = inFesMater; % indexes of materials that cannot be stacked on all kinds of racks
            Problem(i).cons = cons;

            Data(i).orgData(inFesMater,:) = []; % delete unstacked materials
            Data(i).stackSt(inFesMater) = [];
            Problem(i).cons(inFesMater) = [];

            D = size(Data(i).orgData,1); % dimension of decision space
            Problem(i).bound = [ones(1,D);repmat(size(rackSize,1),1,D)]; % product i's decision space's boundry, i.e., (min, max) indexes of kinds of racks being used
        end

        output1 = Problem(instance);
        output2 = Data(instance);

    case 'evaluate' % % evaluate fitness
        Data = varargin{1};
        Decs = varargin{2};

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
        for i = 1:D % for each material
            materL(Data.stackSt{i}(:,4)==1,i) = materSize(i,2);
            materL(Data.stackSt{i}(:,4)==2,i) = materSize(i,1);
            NL(:,i) = Data.stackSt{i}(:,1);
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
                j_materL   = materL(rackInd(j),solution==rackInd(j)); % length side of the materials (belong to the given product) that use rack j
                j_NL       = NL(rackInd(j),solution==rackInd(j)); % number of materials stacked along the length side of rack j
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
                    k = 1; % counter of materials that use rack j
                    currL = 0;
                    currVol = 0;
                    C(j) = 0;
                    while k <= length(j_materL) % while each material
                        while currL < j_rackL && k <= length(j_materL)
                            currL = currL+j_materL(ind(k))*j_NL(ind(k));
                            currVol = currVol+j_materVol(ind(k));
                            k = k+1;
                        end
                        if currL > j_rackL
                            currL = currL-sum(j_materL(ind(k-1))*j_NL(ind(k-1)));
                            currVol = currVol-j_materVol(ind(k-1));
                            k = k-1;
                        end
                        C(j) = C(j)+(j_rackVol-currVol);
                        Clength(j) = Clength(j)+(j_rackL-currL);
                        NumRack(j) = NumRack(j)+1;
                        currL = 0;
                        currVol = 0;
                    end
                end
            end
            Objs(i) = sum(C);
            WasteL(i) = sum(Clength);
            NumRacks{i} = NumRack;
        end
        output1 = Objs;
        output3 = {WasteL,NumRacks};
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end