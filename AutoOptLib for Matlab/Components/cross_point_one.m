function [output1,output2] = cross_point_one(varargin)
% One-point crossover.

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
        Aux    = varargin{4};
        
        Parent = Parent.decs;
        [N,D] = size(Parent);
        Parent1 = Parent(1:ceil(N/2),:);
        Parent2 = Parent(floor(N/2)+1:end,:);
        Nhalf = size(Parent1,1);

        k = randi(D,Nhalf,1);
        Offspring1 = Parent1;
        Offspring2 = Parent2;
        Offspring1(:,k:end) = Parent2(:,k:end);
        Offspring2(:,k:end) = Parent1(:,k:end);
        Offspring = [Offspring1;Offspring2];
        output1   = Offspring(1:N,:);
        output2   = Aux;

    case 'parameter'
        % n/a

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