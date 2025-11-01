function [output1,output2,output3] = CEC2013_f19(varargin)
% Rotated Expanded Griewank's plus Rosenbrock's Function (CEC2013)

switch varargin{end}
    case 'construct'
        type     = {'continuous','static','certain'};
        Problem  = varargin{1};
        instance = varargin{2};
        orgData  = load('shift_data.mat');
        o        = orgData.data(1,:);
        Data     = struct('o',[],'M',[]);
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
                Data(i).M = M.data(1:D,1:D);
            else
                A = normrnd(0,1,D,D);
                [Data(i).M,~] = cGram_Schmidt(A);
            end
        end
        output1 = Problem; output2 = Data;

    case 'repair'
        output1 = varargin{2};

    case 'evaluate'
        Data = varargin{1};
        o    = Data.o; M = Data.M;
        X    = varargin{2};
        [N,D] = size(X);

        Z = 5*(X - repmat(o,N,1))/100; Z = Z*M; Z = Z + 1;

        % pairwise Rosenbrock, then Griewank per scalar, ring connection
        fit = zeros(N,1);
        for i = 1:D
            j = i+1; if j>D, j=1; end
            u = 100*(Z(:,i).^2 - Z(:,j)).^2 + (Z(:,i)-1).^2; % Rosenbrock
            g = u/4000 - cos(u) + 1;                          % Griewank
            fit = fit + g;
        end
        output1 = fit + 500;
end

if ~exist('output2','var'), output2 = []; end
if ~exist('output3','var'), output3 = []; end
end

