function [output1,output2] = cross_sim_binary(varargin)
% Simulated binary crossover.

%------------------------------Reference-----------------------------------
% Deb K, Agrawal R B. Simulated binary crossover for continuous search 
% space[J]. Complex systems, 1995, 9(2): 115-148.
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
        Parent = varargin{1};
        Para   = varargin{3};
        Aux    = varargin{4};

        Parent = Parent.decs;
        DisC = Para;
            
        [N,D] = size(Parent);
        Parent1 = Parent(1:ceil(N/2),:);
        Parent2 = Parent(floor(N/2)+1:end,:);
        Nhalf = size(Parent1,1);

        beta = zeros(Nhalf,D);
        mu = rand(Nhalf,D);

        beta(mu<=0.5) = (mu(mu<=0.5)*2).^(1/(1+DisC));
        beta(mu>0.5)  = (2-mu(mu>0.5)*2).^(-1/(1+DisC));
        Offspring = [0.5*((1+beta).*Parent1+(1-beta).*Parent2);
            0.5*((1-beta).*Parent1+(1+beta).*Parent2)];        
        output1   = Offspring(1:N,:);
        output2   = Aux;

    case 'parameter'
        output1 = [20,40]; % crossover distribution

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