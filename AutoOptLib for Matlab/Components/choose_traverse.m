function [output1,output2] = choose_traverse(varargin)
% Select each of the cureent soutions.

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
        Problem  = varargin{2};
        index = 1:Problem.N;
        output1 = index;

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