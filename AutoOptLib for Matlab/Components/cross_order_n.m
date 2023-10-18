function [output1,output2] = cross_order_n(varargin)
% N-order crossover.

%------------------------------Copyright-----------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

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

        n = round(Para);
        Parent = Parent.decs;
        [N,D] = size(Parent);
        Parent1 = Parent(1:ceil(N/2),:);
        Parent2 = Parent(floor(N/2)+1:end,:);
        Nhalf = size(Parent1,1);

        Offspring1 = Parent1;
        Offspring2 = Parent2;

        seed = zeros(Nhalf,n);
        for i = 1:Nhalf
            seed(i,:) = randperm(D,n); % the n points to be exchanged
        end

        for i = 1:Nhalf
            k = seed(i,:);
            Temp = setdiff(Parent2(i,:),Parent1(i,k),'stable'); % variables of parent 2 that are not appeared in the selected segement (k) of parent 1
            ind  = setdiff(1:D,k,'stable'); % indexes of variables of parent 2 that are not appeared in the selected segement (k) of parent 1
            Offspring1(i,ind) = Temp;

            Temp = setdiff(Parent1(i,:),Parent2(i,k),'stable'); % variables of parent 1 that are not appeared in the selected segement (k) of parent 2
            ind  = setdiff(1:D,k,'stable'); % indexes of variables of parent 1 that are not appeared in the selected segement (k) of parent 2
            Offspring2(i,ind) = Temp;
        end
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
        n_max = max(1,round(D*0.5)); % the maximum n in n-point crossover
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