function [output1,output2] = cross_point_n(varargin)
% Randomly select n indices and exchange the n elements in these n indices
% between a pair of solutions. Similar with the cross_uniform operator.

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
        
        Parent = Parent.decs;
        n = round(Para);
        [N,D] = size(Parent);
        Parent1 = Parent(1:ceil(N/2),:);
        Parent2 = Parent(floor(N/2)+1:end,:);
        Nhalf = size(Parent1,1);

        k = zeros(Nhalf,n);
        for i = 1:Nhalf
            k(i,:) = randperm(D,n); % the n points to be exchanged
        end

        Offspring1 = Parent1;
        Offspring2 = Parent2;
        Offspring1(:,k(:,n)) = Parent2(:,k(:,n));
        Offspring2(:,k(:,n)) = Parent1(:,k(:,n));
        Offspring = [Offspring1;Offspring2];
        output1   = Offspring(1:N,:);
        output2   = Aux;

    case 'parameter'
        % number of problem's decision variable
        Problem  = varargin{1};
        D = zeros(length(Problem),1);
        for i = 1:length(Problem)
            D(i) = size(Problem(i).bound,2);
        end
        D = min(D);
        n_max = max(1,round(D*0.5)); % the maximum n in n-point crossover;
        output1 = [1,n_max]; % the n point

    case 'behavior'
        output1 = {'LS','small';'GS','large'}; % small n values perform local search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end