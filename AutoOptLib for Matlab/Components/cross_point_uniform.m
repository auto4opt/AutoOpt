function [output1,output2] = cross_point_uniform(varargin)
% Uniform crossover.

%------------------------------Copyright-----------------------------------
% Copyright (C) <2023>  <Qi Zhao>

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

        Prob = Para;
        Parent = Parent.decs;
        [N,D] = size(Parent);
        Parent1 = Parent(1:ceil(N/2),:);
        Parent2 = Parent(floor(N/2)+1:end,:);
        Nhalf = size(Parent1,1);

        ind = rand(Nhalf,D) < Prob;
        Offspring1 = Parent1;
        Offspring2 = Parent2;
        Offspring1(ind) = Parent2(ind);
        Offspring2(ind) = Parent1(ind);
        Offspring = [Offspring1;Offspring2];
        output1   = Offspring(1:N,:);
        output2   = Aux;

    case 'parameter'
        output1 = [0,0.5]; % crossover probability

    case 'behavior'
        output1 = {'LS','small';'GS','large'}; % small probabilities perform local search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end