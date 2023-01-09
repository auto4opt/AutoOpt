function [output1, output2, output3] = CEC2005_f5(varargin)
% The Schwefel's Problem 2.6 with Global Optimum on Bounds from the
% benchmark for the CEC 2005 Special Session on Real-Parameter Optimization.

switch varargin{end}
    case 'construct'
        type     = {'continuous','static','certain'};
        Problem  = varargin{1};
        instance = varargin{2};
        orgData  = load('schwefel_206_data');
        o        = orgData.o;
        A        = orgData.A;
        Data     = struct('A',[],'B',[]);
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
                curr_A = A(1:D,1:D);
            else
                curr_o = -100+200*rand(1,D);
                curr_A = round(-100+2*100.*rand(D,D));
                while det(curr_A) == 0
                    curr_A = round(-100+2*100.*rand(D,D));
                end
            end
            curr_o(1:ceil(D/4)) = -100;
            curr_o(max(floor(0.75*D),1):D) = 100;
            curr_B = curr_A*curr_o';
            Data(i).A = curr_A;
            Data(i).B = curr_B;
        end
        output1 = Problem;
        output2 = Data;
        
    case 'evaluate'
        Data  = varargin{1};
        A     = Data.A;
        B     = Data.B;
        Decs  = varargin{2};
        [N,~] = size(Decs);
        
        fit = zeros(N,1);
        for i = 1:N
            fit(i) = max(abs(A*(Decs(i,:)')-B));
        end
        output1 = fit-310;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end