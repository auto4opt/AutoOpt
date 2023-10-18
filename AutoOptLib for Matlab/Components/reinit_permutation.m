function [output1,output2] = reinit_permutation(varargin)
% Reinitialization for permutation problems.

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
        Aux      = varargin{4};

        [N,D] = size(Solution.decs);
        [~,output1] = sort(rand(N,D),2);
        output2 = Aux;
    
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