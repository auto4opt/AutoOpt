function [output1,output2,output3] = CEC2013_f11(varargin)
% The f11 Function from the benchmark for the CEC 2013 Special
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
        end      
        output1 = Problem;
        output2 = Data;
        
   case 'repair'
        Decs = varargin{2};
        output1 = Decs;
    
    case 'evaluate'
        Data = varargin{1};
        o    = Data.o;
        Decs = varargin{2};
        
        [N,D] = size(Decs);
        Decs  = Decs-repmat(o,N,1);
        % disp(size(Decs))
        % disp(size(M))
        Decs = 5.12*Decs/100;
        Decs = computeTosz(Decs);
        Decs = computeTAsym(Decs, 0.2);
        z = Decs*constructLambda(10, D);
        fit = 0;
        for i = 1:D
           fit = fit + (z(:,i).^2 - 10.*cos(2.*pi.*z(:,i)) + 10);
        end
        
        output1= fit-400;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end