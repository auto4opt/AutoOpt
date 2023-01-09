function [output1,output2] = archive_best(varargin)
% Collect the best (in terms of fitness) solution of each iteration.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        ArchBest = varargin{2};

        Fitness  = Solution.fits;
        [~,best]  = min(Fitness);
        output1  = [ArchBest,Solution(best)];
    
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