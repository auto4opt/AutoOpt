function [output1,output2,output3] = CEC2013_f4(varargin)
% The f4 Function from the benchmark for the CEC 2013 Special
% Session on Real-Parameter Optimization.

%------------------------------Reference-----------------------------------
% Liang J J, Qu B Y, Suganthan P N, et al. Problem definitions and 
% evaluation criteria for the CEC 2013 special session on real-parameter 
% optimization[R]. Computational Intelligence Laboratory, Zhengzhou 
% University, Zhengzhou, China and Nanyang Technological University, 
% Singapore, Technical Report, 2013, 201212(34): 281-295.
%------------------------------Copyright-----------------------------------
% Copyright (C) <2025>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 
%--------------------------------------------------------------------------

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
            else
                A = normrnd(0, 1, D, D);
                [M1, ~] = cGram_Schmidt(A);
            end
            Data(i).M1 = M1;
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
        Decs  = Decs-repmat(o,N,1);
        % disp(size(Decs))
        % disp(size(M))
        Decs  = Decs*M1;
        Decs  = computeTosz(Decs);

        a = 1e+6;
        fit = a.*Decs(:,1).^2;
        for i = 1:D
            fit = fit+Decs(:,i).^2;
        end
        output1= fit-1100;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end