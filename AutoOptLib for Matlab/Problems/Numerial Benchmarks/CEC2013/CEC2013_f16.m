function [output1,output2,output3] = CEC2013_f16(varargin)

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
        Decs = Decs - repmat(o, N, 1);
        Decs = (5/100)*Decs;
        Decs = Decs * M1;
        Decs = Decs * constructLambda(100, D);
        z = Decs * M2;

        exponent = 10 / (D^(1.2));
        
        productVal = 1;
        
        for i = 1:D
            zi = z(:,i);

            sum_j = 0;
            for j = 1:32
                factor = 2^j * zi;
                sum_j = sum_j + abs(factor - round(factor))/ (2^j);
            end

            innerTerm = 1 + i * sum_j;

            factor_i = innerTerm.^exponent;
            
            productVal = productVal .* factor_i;
        end
        
        fit = (10/(D^2))*productVal - (10/(D^2));
        
        output1 = fit + 200;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end