function [output1,output2] = archive_best(varargin)
% Collect the best (in terms of fitness) solution of each iteration.

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
        Solution    = varargin{1};
        CurrArchive = varargin{2};

        Fitness  = SOLVE.fits(Solution);
        [~,best] = min(Fitness);
        output1  = SOLVE.cats(CurrArchive,Solution(best));
    
    case 'parameter'
        % no parameter

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