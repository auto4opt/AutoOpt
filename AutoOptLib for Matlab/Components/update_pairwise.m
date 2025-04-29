function [output1,output2] = update_pairwise(varargin)
% Select the better solution from each pair of old and new solutions.

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

        Fitness  = Solution.fits;
        tempN    = size(Fitness,1)/2;
        Old      = Fitness(1:tempN);
        New      = Fitness(tempN+1:end);

        compare1 = Old<=New;
        compare2 = Old>New;

        ind      = 1:length(Solution);
        ind      = ind([compare1;compare2]);
        output1  = Solution(ind);

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