function [output1,output2,output3] = Place(varargin)
maxH = 200;
switch varargin{end}
    case 'construct' % construct problem and data
        type     = {'permutation','sequential','certain'};
        Problem  = varargin{1};
        instance = varargin{2};

        orgData = load('DataStacked.mat');
        Data = struct('orgData',[],'rackH',[],'maxH',[],'rackInd',[],'prodInd',[],'preSolRack',[],'preSolProd',[],'continue',[]);
        
        lowerH = 20;
        upperH = 150;
        intervalH = 10;
        H = (lowerH:intervalH:upperH)';
        numH = length(H);
        rackSize = [repmat(180,numH,1),repmat(60,numH,1),H;repmat(240,numH,1),repmat(120,numH,1),H;]; % sizes of all kinds of racks, [long,width,height]
        rackH = rackSize(:,3);
        
        numInstance = size(orgData.usedRacks,2);
        for i = 1:numInstance
            Problem(i) = Problem(1); % a number of numProd problem instances
            Problem(i).type = type;
            Data(i).orgData = orgData.usedRacks(:,i); % cells, each collects indexes (can be repetitive) of all racks used for one product         
            Data(i).rackH   = rackH;
            Data(i).maxH    = maxH;
            Data(i).prodInd = 1; % index of the current considered product
            Data(i).rackInd = Data(i).orgData{Data(i).prodInd}; % indexes of racks that the current product used

            totalH = sum(Data(i).rackH(Data(i).rackInd));
            Data(i).preSolRack = [];
            Data(i).preSolProd = [];
            while totalH < maxH && Data(i).prodInd < length(Data(i).orgData)
                Data(i).preSolRack = [Data(i).preSolRack,Data(i).rackInd]; % a part of the solution (rack indexes) to the current problem
                Data(i).preSolProd = [Data(i).preSolProd,repmat(Data(i).prodInd,1,numel(Data(i).rackInd))]; % a part of the solution (product indexes) to the current problem
                Data(i).maxH = maxH-totalH; % vertical height being considered in the current problem
                Data(i).prodInd = Data(i).prodInd+1; % index of the current considered product
                Data(i).rackInd = Data(i).orgData{Data(i).prodInd};
                totalH = totalH+sum(Data(i).rackH(Data(i).rackInd));
            end

            Problem(i).bound = [1;numel(Data(i).rackInd)];

            if Data(i).prodInd == length(Data(i).orgData) && totalH <= maxH
                Data(i).continue = false;
            else
                Data(i).continue = true;
            end
        end

        output1 = Problem(instance);
        output2 = Data(instance);

    case 'evaluate' % evaluate fitness
        Data = varargin{1};
        Decs = varargin{2};

        N = size(Decs,1); % population size
        D = size(Decs,2); % deminsion of solution, i.e., number of racks
        Objs = zeros(N,1);
        acc1 = zeros(N,1);
        acc2 = cell(N,1); % explicit solutions, i.e., racks' placement
        for i = 1: N
            currH = 0;
            j = 1;
            k = [];
            while currH < Data.maxH && j <= D
                currH = currH+Data.rackH(Data.rackInd(Decs(i,j)));
                k = [k,Data.rackInd(Decs(i,j))]; % indexes of racks that have been placed
                if currH > Data.maxH
                    currH = currH-Data.rackH(Data.rackInd(Decs(i,j)));
                    k(end) = []; % delete the rack that cannot be placed
                end
                j = j+1;
            end

            Objs(i) = Data.maxH - currH;
            acc2{i}(1,:) = [Data.preSolRack,k]; % all racks placed in the current column
            acc2{i}(2,:) = [Data.preSolProd,repmat(Data.prodInd,1,numel(k))]; % each rack's product index
        end
        output1 = Objs;
        output3 = {acc1,acc2};

    case 'sequence' % change time step in the problem sequence
        Problem = varargin{1};
        Data = varargin{2};
        solution = varargin{3};
        
        Data.maxH = maxH;
        usedRackInd = solution.acc{2}(1,:);
        usedRackInd = usedRackInd(solution.acc{2}(2,:)==Data.prodInd); % current product's used racks
        for i = 1:numel(usedRackInd)
            sameInd = find(Data.rackInd==usedRackInd(i));
            Data.rackInd(sameInd(1)) = [];
        end

        totalH = sum(Data.rackH(Data.rackInd));
        Data.preSolRack = [];
        Data.preSolProd = [];
        while totalH < maxH && Data.prodInd < length(Data.orgData)
            Data.preSolRack = [Data.preSolRack,Data.rackInd]; % a part of the solution (rack indexes) to the current problem
            Data.preSolProd = [Data.preSolProd,repmat(Data.prodInd,1,numel(Data.rackInd))]; % a part of the solution (product indexes) to the current problem
            Data.maxH = maxH-totalH; % vertical height being considered in the current problem
            Data.prodInd = Data.prodInd+1; % index of the current considered product
            Data.rackInd = Data.orgData{Data.prodInd};
            totalH = totalH+sum(Data.rackH(Data.rackInd));
        end

        while numel(Data.rackInd) == 1 && Data.prodInd < length(Data.orgData) % if only one rack for placement
            Data.prodInd = Data.prodInd+1;
            Data.rackInd = [Data.rackInd,Data.orgData{Data.prodInd}];
            totalH = sum(Data.rackH(Data.rackInd));
        end

        Problem.bound = [1;numel(Data.rackInd)];

        if Data.prodInd == length(Data.orgData) && totalH <= maxH
            Data.continue = false;
        else
            Data.continue = true;
        end

        output1 = Problem;
        output2 = Data;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end