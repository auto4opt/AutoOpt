function [output1,output2,output3] = CEC2013_f17(varargin)

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
        Decs = varargin{2};
        
        [N,D] = size(Decs);
        mu0   = 2.5;
        d     = 1;
        x     = Decs;
        
        s = 1 - 1/(2*sqrt(D+20)-8.2);
        mu1 = -sqrt((mu0^2 - d)/s);

        Decs = Decs - repmat(o, N, 1);
        y = 0.1 * Decs;
        y = 2 * y;
        
        x_prime = zeros(D);
        z = (x_prime - mu0) * constructLambda(100, D);
        z_prime = zeros(D);


        for i=1:D
            x_prime(i) = sign(x(i)) * y(i) + mu0;
            z_prime(i) = sign(z(i)) * y(i) + mu0;
        end
        

        term1 = 0;
        term2 = 0;
        term3 = 0;
        for i=1:D
            term1 = term1 + (x_prime(i) - mu0).^2;
            term2 = term2 + (x_prime(i) - mu1).^2;
            term3 = term3 + cos(2*pi*z_prime(i));
        end

        fit = min(term1, D + s*term2) + 10*(D-term3);
        output1 = fit + 300;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end