function [output1,output2,output3] = CEC2013_f9(varargin)

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
            lower = zeros(1,D)-100;
            upper = zeros(1,D)+100;
            Problem(i).bound = [lower;upper];
            
            if length(o) >= D
                curr_o = o(1:D);
            else 
                curr_o = -100+200*rand(1,D);
            end
            Data(i).o = curr_o;
            
            fileMap = containers.Map([2, 5, 10, 20, 30, 40, 50, 70, 80, 90, 100], ...
                         {'M_D2.mat', 'M_D5.mat', 'M_D10.mat', 'M_D20.mat', 'M_D30.mat', 'M_D40.mat', 'M_D50.mat', 'M_D70.mat', 'M_D80.mat', 'M_D90.mat', 'M_D100.mat'});
            if isKey(fileMap, D)
                M = load(fileMap(D));
                M1 = M.data(1:D,1:D);
                M2 = M.data(D+1:2*D,1:D);
            else
                A = normrnd(0, 1, D, D);
                [M1, ~] = cGram_Schmidt(A);
                [M2, ~] = cGram_Schmidt(A);
            end
            Data(i).M1 = M1;
            Data(i).M2 = M2;
        end      
        output1 = Problem;
        output2 = Data;
        
   case 'repair'
        Decs = varargin{2};
        output1 = Decs;
    
    case 'evaluate'
        Data = varargin{1};
        o    = Data.o;
        M1   = Data.M1;
        M2   = Data.M2;
        Decs = varargin{2};
        
        [N,D] = size(Decs);
        Decs  = Decs-repmat(o,N,1);
        % disp(size(Decs))
        % disp(size(M))
        Decs = 0.5*Decs/100;
        Decs  = Decs*M1;
        Decs  = computeTAsym(Decs,0.5);
        Decs  = Decs*M2;
        Decs  = Decs*constructLambda(10,D);

        a = 0.5;
        b = 3;
        kmax = 20;

        term1 = 0;
        term2 = 0;
        for i = 1:D
            inner_sum = 0;  % 对第 i 维度的内部 k 循环累加
            for k = 0:kmax
                z = Decs(:,i);
                inner_sum = inner_sum + a^k * cos(2 * pi * b^k * (z + 0.5));
            end
            term1 = term1 + inner_sum;
        end
        for k = 0:kmax
            term2 = term2 + a^k * cos(2 * pi * b^k * 0.5);
        end
        fit = term1-D*term2;
        output1= fit-600;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end