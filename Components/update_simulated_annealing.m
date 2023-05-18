function [output1,output2] = update_simulated_annealing(varargin)
% Simulated annealing's update mechanism, i.e., accept worse soluion with a
% probability.

%------------------------------Reference-----------------------------------
% Kirkpatrick S, Gelatt Jr C D, Vecchi M P. Optimization by simulated 
% annealing[J]. Science, 1983, 220(4598): 671-680.
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
mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};
        Para     = varargin{3};
        G        = varargin{5};

        T_initial = Para; % initial temperture
        T_final   = 0.01; % final temperture
        Rate      = nthroot(T_final/T_initial,Problem.Gmax); % temperture decrease rate
        T         = T_initial*Rate^G; % current temperture

        Old = Solution(1:Problem.N);
        New = Solution(Problem.N+1:end);
      
        accept = rand(Problem.N,1) < exp((Old.fits-New.fits)./abs(Old.fits+1e-6)./T); % accept worse solutions by annealing
        accept(Old.fits>New.fits) = true; % always accept better solutions      
        Old(accept) = New(accept);

        output1 = Old;

    case 'parameter'
        output1 = [0.1,1]; % initial temperture 

    case 'behavior'
        output1 = {'';''}; 
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end