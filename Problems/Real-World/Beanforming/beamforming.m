function [output1,output2,output3] = beamforming(varargin)
% The beanforming problem in RIS-aided communications.

%------------------------------Reference-----------------------------------
% Yan B, Zhao Q, Li M, et al. Fitness landscape analysis and niching 
% genetic approach for hybrid beamforming in RIS-aided communications[J]. 
% Applied Soft Computing, 2022, 131: 109725.
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
    case 'construct' % define problem properties
        Problem  = varargin{1};
        instance = varargin{2};       

        orgData  = load('Beanforming.mat','Data');
        Data = orgData.Data((instance));
        for i = 1:length(instance)
            D = size(Data(i).G,1);
            phases_cnt = 2^Data(i).b-1;
            lower = zeros(1,D); % 1*D, lower bound of the D-dimension decision space
            upper = repmat(phases_cnt,1,D); % 1*D, upper bound of the D-dimension decision space
            Problem(i).type = {'discrete','static','certain'};
            Problem(i).bound = [lower;upper];
        end
        
        output1 = Problem;
        output2 = Data;

    case 'repair' % repair solutions
        Decs = varargin{2};
        output1 = Decs;
    
    case 'evaluate' % evaluate solution's fitness
        Data = varargin{1}; % load problem data
        m    = varargin{2}; % load the current solution(s)
        
        b     = Data.b;
        PT    = Data.PT;
        G     = Data.G;
        Hd    = Data.Hd;
        Hr    = Data.Hr;
        omega = Data.omega;
      
        sR = get_sum_rate(m,b,Hd, Hr,G,PT,omega); % calculate objective value
        sR = sR'; % N*1

        output1 = sR; % matrix for saving objective function values
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end