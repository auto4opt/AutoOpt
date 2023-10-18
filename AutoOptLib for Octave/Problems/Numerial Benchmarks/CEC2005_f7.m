function [output1,output2,output3] = CEC2005_f7(varargin)
% The Shifted Rotated Griewank's Function without Bounds from the benchmark
% for the CEC 2005 Special Session on Real-Parameter Optimization.

%------------------------------Reference-----------------------------------
% Suganthan P N, Hansen N, Liang J J, et al. Problem definitions and 
% evaluation criteria for the CEC 2005 special session on real-parameter 
% optimization[R]. KanGAL report, 2005.
%------------------------------Copyright-----------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 

% Please reference the paper below if using AutoOptLib in your publication:
% @article{zhao2023autooptlib,
%  title={AutoOptLib: A Library of Automatically Designing Metaheuristic 
%         Optimization Algorithms in Matlab},
%  author={Zhao, Qi and Yan, Bai and Hu, Taiwei and Chen, Xianglong and 
%          Yang, Jian and Shi, Yuhui},
%  journal={arXiv preprint 	arXiv:2303.06536},
%  year={2023}
% }
%--------------------------------------------------------------------------

switch varargin{end}
    case 'construct'
        type     = {'continuous','static','certain'};
        Problem  = varargin{1};
        instance = varargin{2};
        orgData  = load('griewank_func_data');
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
            lower = zeros(1,D) - 600;
            upper = zeros(1,D) + 600;
            Problem(i).bound = [lower;upper];
            
            if length(o) >= D
                curr_o = o(1:D);
            else
                curr_o = -600+0*rand(1,D);
            end
            Data(i).o = curr_o;
            
            if D == 2
                M = load('griewank_M_D2');
                M = M.M;
            elseif D == 10
                M = load('griewank_M_D10');
                M = M.M;
            elseif D == 30
                M = load('griewank_M_D30');
                M = M.M;
            elseif D == 50
                M = load('griewank_M_D50');
                M = M.M;
            else
                M = rot_matrix(D,3);
                M = M.*(1+0.3.*normrnd(0,1,D,D));
            end
            Data(i).M = M;
        end      
        output1 = Problem;
        output2 = Data;
        
    case 'repair'
        Decs = varargin{2};
        output1 = Decs;
    
    case 'evaluate'
        Data  = varargin{1};
        o     = Data.o;
        M     = Data.M;
        Decs  = varargin{2};
        
        [N,D] = size(Decs);
        Decs  = Decs-repmat(o,N,1);
        Decs  = Decs*M;
        
        fit = 1;
        for i = 1:D
            fit = fit.*cos(Decs(:,i)./sqrt(i));
        end
        fit = sum(Decs.^2,2)./4000-fit+1;
        output1= fit-180;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end