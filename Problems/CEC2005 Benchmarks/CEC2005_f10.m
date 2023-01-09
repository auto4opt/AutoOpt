function [output1,output2,output3] = CEC2005_f10(varargin)
% The Shifted Rotated Rastrigin's Function from the benchmark for the CEC
% 2005 Special Session on Real-Parameter Optimization.

switch varargin{end}
    case 'construct'
        type     = {'continuous','static','certain'};
        Problem  = varargin{1};
        instance = varargin{2};
        orgData  = load('rastrigin_func_data');
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
            lower = zeros(1,D)-5;
            upper = zeros(1,D)+5;
            Problem(i).bound = [lower;upper];
            
            if length(o) >= D
                curr_o = o(1:D);
            else
                curr_o = -5+10*rand(1,D);
            end
            Data(i).o = curr_o;
            
            if D == 2
                M = load('rastrigin_M_D2');
                M = M.M;
            elseif D == 10
                M = load('rastrigin_M_D10');
                M = M.M;
            elseif D == 30
                M = load('rastrigin_M_D30');
                M = M.M;
            elseif D == 50
                M = load('rastrigin_M_D50');
                M = M.M;
            else
                M = rot_matrix(D,2);
            end
            Data(i).M = M;
        end
        output1 = Problem;
        output2 = Data;
        
    case 'evaluate'
        Data = varargin{1};
        o    = Data.o;
        M    = Data.M;
        Decs = varargin{2};
        
        [N,~] = size(Decs);
        Decs  = Decs-repmat(o,N,1);
        Decs  = Decs*M;
        
        fit = sum(Decs.^2-10.*cos(2.*pi.*Decs)+10,2);
        output1 = fit-330;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end