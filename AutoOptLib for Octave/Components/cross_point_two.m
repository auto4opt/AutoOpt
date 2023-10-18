function [output1,output2] = cross_point_two(varargin)
% The two point crossover operation in genetic algorithm.

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
        Aux    = varargin{4};

        Parent  = SOLVE.decs(Parent);
        [N,D]   = size(Parent);
        Parent1 = Parent(1:ceil(N/2),:);
        Parent2 = Parent(floor(N/2)+1:end,:);
        Nhalf   = size(Parent1,1);

        if D < 2
            error(['Two-point crossover is unavailable for 1-dimensional ' ...
                'problems. Please remove the two-point crossover from the ' ...
                'Design_Space.'])
        end
        
        k = zeros(Nhalf,2);
        for i = 1:Nhalf
            k(i,:) = randperm(D,2);
        end
        k = sort(k,2);

        Offspring1 = Parent1;
        Offspring2 = Parent2;
        Offspring1(:,k(:,1):k(:,2)) = Parent2(:,k(:,1):k(:,2));
        Offspring2(:,k(:,1):k(:,2)) = Parent1(:,k(:,1):k(:,2));
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