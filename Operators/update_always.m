function [output1,output2] = update_always(varargin)
% Always select newly generated solutions.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};

        output1 = Solution(end-Problem.N+1:end);

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