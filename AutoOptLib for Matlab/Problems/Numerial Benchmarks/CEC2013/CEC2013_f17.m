function [output1,output2,output3] = CEC2013_f17(varargin)
% Lunacek bi-Rastrigin Function (CEC2013)

switch varargin{end}
    case 'construct'
        type     = {'continuous','static','certain'};
        Problem  = varargin{1};
        instance = varargin{2};
        orgData  = load('shift_data.mat');
        o        = orgData.data(1,:);
        Data     = struct('o',[]);
        for i = 1:length(instance)
            Problem(i).type = type;
            D = instance(i);
            Problem(i).bound = [ones(1,D)*-100; ones(1,D)*100];
            if length(o) >= D
                Data(i).o = o(1:D);
            else
                Data(i).o = -100+200*rand(1,D);
            end
        end
        output1 = Problem; output2 = Data;

    case 'repair'
        output1 = varargin{2};

    case 'evaluate'
        Data = varargin{1};
        o    = Data.o;
        X    = varargin{2};
        [N,D] = size(X);

        mu0 = 2.5; d = 1;
        s = 1 - 1/(2*sqrt(D+20)-8.2);
        mu1 = -sqrt((mu0^2 - d)/s);

        Y = 10*(X - repmat(o,N,1))/100;
        sgn = ones(1,D); sgn(o<=0) = -1;
        Xhat = 2*repmat(sgn,N,1).*Y + mu0;
        Z = Xhat - mu0;
        Z = Z*constructLambda(100,D);

        term1 = sum((Xhat-mu0).^2,2);
        term2 = d*D + s*sum((Xhat-mu1).^2,2);
        base  = min(term1, term2);
        ras   = 10*(D - sum(cos(2*pi*Z),2));
        output1 = base + ras + 300;
end

if ~exist('output2','var'), output2 = []; end
if ~exist('output3','var'), output3 = []; end
end

