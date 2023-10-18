function [output1,output2] = archive_statistic(varargin)
% Collect the average and standard deviation of solution's fitness of each
% iteration.

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
        Solution    = varargin{1};
        CurrArchive = varargin{2};

        Fitness = SOLVE.fits(Solution);
        output1 = [mean(Fitness),std(Fitness)];
        output1 = [CurrArchive;output1];
    
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