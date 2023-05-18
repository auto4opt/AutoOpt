function [output1,output2] = search_pso(varargin)
% Particle swarm optimization's particle fly operator.

%------------------------------Reference-----------------------------------
% Shi Y, Eberhart R. A modified particle swarm optimizer[C]//1998 IEEE 
% International Conference on Evolutionary Computation Proceedings. IEEE 
% World Congress on Computational Intelligence (Cat. No. 98TH8360). IEEE, 
% 1998: 69-73.
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
        Para     = varargin{3};
        Aux      = varargin{4};

        [N,D] = size(Solution.decs);
        W = Para; % inertia weight

        % initialize pbest, gbest, and velocity
        if ~isfield(Aux,'Pbest')
            Aux.Pbest = Solution; % personal best solutions
            Aux.Gbest = Solution(randi(N)); % global best solution
            Aux.V = zeros(N,1);
        end

        % particle fly
        Dec   = Solution.decs;
        Pbest = Aux.Pbest.decs;
        Gbest = Aux.Gbest.dec;

        r1 = repmat(rand(N,1),1,D);
        r2 = repmat(rand(N,1),1,D);
        V  = W.*Aux.V+2*r1.*(Pbest-Dec)+2*r2.*(Gbest-Dec); % N*D
        output1 = Dec+V; % N*D

        % update volecity
        Aux.V   = V;
        output2 = Aux;

    case 'parameter'
        output1 = [0,0.5]; % inertia weight

    case 'behavior'
        output1 = {'';'GS'}; % always perform global search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end