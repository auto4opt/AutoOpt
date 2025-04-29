function [output1,output2] = search_de_random(varargin)
% The "random/1" differential mutation.

%------------------------------Reference-----------------------------------
% Storn R, Price K. Differential evolution-a simple and efficient heuristic
% for global optimization over continuous spaces[J]. Journal of Global 
% Optimization, 1997, 11(4): 341-359.
%------------------------------Copyright-----------------------------------
% Copyright (C) <2025>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 
%--------------------------------------------------------------------------

mode = varargin{end};
switch mode
    case 'execute'
        Parent = varargin{1};
        Para   = varargin{3};
        Aux    = varargin{4};

        F       = Para(1);
        CR      = Para(2);
        Parent  = SOLVE.decs(Parent);

        [N,D]   = size(Parent);
        Parent1 = Parent(randperm(N),:);
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