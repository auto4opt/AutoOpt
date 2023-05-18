function [output1,output2] = search_de_current(varargin)
% The "current/1" differential mutation.
%------------------------------Reference-----------------------------------
% Storn R, Price K. Differential evolution-a simple and efficient heuristic
% for global optimization over continuous spaces[J]. Journal of Global 
% Optimization, 1997, 11(4): 341-359.
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
        F  = Para(1);
        CR = Para(2);      

        [N,D]   = size(Parent);
        Parent1 = Parent;
        Parent2 = Parent(randperm(N),:);
        Parent3 = Parent(randperm(N),:);

        ind = rand(N,D) < CR;
        Offspring = Parent;
        Offspring(ind) = Parent1(ind) + F*(Parent2(ind)-Parent3(ind));
        output1 = Offspring;
        output2 = Aux;

    case 'parameter'
        output1 = [0,1;0,1]; % F and CR 

    case 'behavior'
        output1 = {'LS','small','small';'GS','large','large'}; % small F and CR values perform local search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end