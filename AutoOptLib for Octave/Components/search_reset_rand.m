function [output1,output2] = search_reset_rand(varargin)
% Reset each element of solutions to a random value with a probability.

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
        Solution = varargin{1};
        Problem  = varargin{2};
        Para     = varargin{3};
        Aux      = varargin{4};

        Prob = Para;
        if ~isnumeric(Solution)
            New = SOLVE.decs(Solution);
        else
            New = Solution;
        end
        [N,D] = size(New);

        Temp = zeros(N,D);
        for j = 1:D
            Temp(:,j) = randi([Problem.bound(1,j),Problem.bound(2,j)],N,1);
        end

        ind = rand(N,D) < Prob; % N*D
        New(ind) = Temp(ind);
        output1 = New;
        output2 = Aux;

    case 'parameter'
        output1 = [0,0.5]; % reset probability 

    case 'behavior'
        output1 = {'LS','small';'GS','large'}; % samll probabilities perform local search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end