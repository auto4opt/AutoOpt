function [output1,output2] = update_round_robin(varargin)
% Conduct update via round-robin tournament.

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

        Fitness  = SOLVE.fits(Solution);
        K   = min(10,Problem.N-1);
        ind = 1:length(Solution);
        win = zeros(length(Solution),1);
        for i = 1:length(Solution)
            win(i) = sum(Fitness(i)<=Fitness(randperm(end,K)));
        end
        [~,rank] = sort(win,'descend');
        ind = ind(rank(1:Problem.N));
        output1 = Solution(ind);

    case 'parameter'
        % n/a

    case 'behavior'
        output1 = {'';''}; 
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end