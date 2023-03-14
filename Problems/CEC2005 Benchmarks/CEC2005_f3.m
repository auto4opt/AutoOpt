function [output1,output2,output3] = CEC2005_f3(varargin)
% The Shifted Rotated High Conditioned Elliptic Function 1.2 from the
% benchmark for the CEC 2005 Special Session on Real-Parameter Optimization.

switch varargin{end}
    case 'construct'
        type     = {'continuous','static','certain'};
        Problem  = varargin{1};
        instance = varargin{2};
        orgData  = load('high_cond_elliptic_rot_data');
        o        = orgData.o;
        Data     = struct('o',[],'M',[]);
        for i = 1:length(instance)
             Problem(i).type = type;

            if instance(i) == 1
                D = 10;
            elseif instance(i) == 2
                D = 30;
            elseif instance(i) == 3
                D = 50;
            else
                error('Only instances 1, 2, and 3 are available.')
            end
            lower = zeros(1,D)-100;
            upper = zeros(1,D)+100;
            Problem(i).bound = [lower;upper];
            
            if length(o) >= D
                curr_o = o(1:D);
            else
                curr_o = -100+200*rand(1,D);
            end
            Data(i).o = curr_o;
            
            if D == 2
                M = load('elliptic_M_D2');
                M = M.M;
            elseif D == 10
                M = load('elliptic_M_D10');
                M = M.M;
            elseif D == 30
                M = load('elliptic_M_D30');
                M = M.M;
            elseif D == 50
                M = load('elliptic_M_D50');
                M = M.M;
            else
                A = normrnd(0,1,D,D);
                [M,~] = cGram_Schmidt(A);
            end
            Data(i).M = M;
        end      
        output1 = Problem;
        output2 = Data;
        
   case 'repair'
        Decs = varargin{2};
        output1 = Decs;
    
    case 'evaluate'
        Data = varargin{1};
        o    = Data.o;
        M    = Data.M;
        Decs = varargin{2};
        
        [N,D] = size(Decs);
        Decs  = Decs-repmat(o,N,1);
        Decs  = Decs*M;
        
        a = 1e+6;
        fit = 0;
        for i = 1:D
            fit = fit+a.^((i-1)/(D-1)).*Decs(:,i).^2;
        end
        output1= fit-450;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end