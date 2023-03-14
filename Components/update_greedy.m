function [output1,output2] = update_greedy(varargin)
% Select the best (in terms of fitness) N solutions.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};

        [~,ind]  = sort(Solution.fits,'ascend');
        output1 = Solution(ind(1:Problem.N));

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