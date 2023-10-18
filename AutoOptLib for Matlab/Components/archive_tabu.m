function [output1,output2] = archive_tabu(varargin)
% Collect N solutions to the tabu list.

%------------------------------Reference-----------------------------------
% Glover F. Tabu search—part I[J]. ORSA Journal on computing, 1989, 1(3): 
% 190-206.
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
        Solution = varargin{1};
        output1  = Solution; % prevent the search visiting the last N solutions.

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