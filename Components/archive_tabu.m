function [output1,output2] = archive_tabu(varargin)
% Collect N solutions to the tabu list.

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