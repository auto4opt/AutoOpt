function [output1,output2,output3] = CEC2013_f18(varargin)
% Rotated Lunacek bi-Rastrigin Function (CEC2013)

switch varargin{end}
    case 'construct'
        type     = {'continuous','static','certain'};
        Problem  = varargin{1};
        instance = varargin{2};
        orgData  = load('shift_data.mat');
        o        = orgData.data(1,:);
        Data     = struct('o',[],'M1',[],'M2',[]);
        for i = 1:length(instance)
            Problem(i).type = type;
            D = instance(i);
            Problem(i).bound = [ones(1,D)*-100; ones(1,D)*100];
            if length(o) >= D
                Data(i).o = o(1:D);
            else
                Data(i).o = -100+200*rand(1,D);
            end
            fileMap = containers.Map([2,5,10,20,30,40,50,70,80,90,100], ...
                {'M_D2.mat','M_D5.mat','M_D10.mat','M_D20.mat','M_D30.mat','M_D40.mat','M_D50.mat','M_D70.mat','M_D80.mat','M_D90.mat','M_D100.mat'});
            if isKey(fileMap,D)
                M = load(fileMap(D));
                Data(i).M1 = M.data(1:D,1:D);
                Data(i).M2 = M.data(D+1:2*D,1:D);
            else
                A = normrnd(0,1,D,D);
                [Data(i).M1,~] = cGram_Schmidt(A);
                [Data(i).M2,~] = cGram_Schmidt(A);
            end
        end
        output1 = Problem; output2 = Data;

    case 'repair'
        output1 = varargin{2};

    case 'evaluate'
        Data = varargin{1};
        o    = Data.o; M1 = Data.M1; M2 = Data.M2;
        X    = varargin{2};
        [N,D] = size(X);

        mu0 = 2.5; d = 1; s = 1 - 1/(2*sqrt(D+20)-8.2);
        mu1 = -sqrt((mu0^2 - d)/s);

        Y = 10*(X - repmat(o,N,1))/100;
        sgn = ones(1,D); sgn(o<=0) = -1;
        Xhat = 2*repmat(sgn,N,1).*Y + mu0;
        Z = (Xhat - mu0)*M1; Z = Z*constructLambda(100,D); Z = Z*M2;

        term1 = sum((Xhat-mu0).^2,2);
        term2 = d*D + s*sum((Xhat-mu1).^2,2);
        base  = min(term1, term2);
        ras   = 10*(D - sum(cos(2*pi*Z),2));
        output1 = base + ras + 400;
end

if ~exist('output2','var'), output2 = []; end
if ~exist('output3','var'), output3 = []; end
end

