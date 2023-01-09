function [output1,output2] = search_swap(varargin)
% Swap two randomly selected elements.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};

        if ~isnumeric(Solution)
            Solution = Solution.decs;
        end      
        [N,D] = size(Solution);
        New = Solution;

        for i = 1:N
            k = randperm(D,2);
            New(i,k(1)) = Solution(i,k(2));
            New(i,k(2)) = Solution(i,k(1));
        end
        output1 = New;
        output2 = varargin{5};

    case 'parameter'
        % n/a

    case 'behavior'
        output1 = {'LS';''}; % always perform local search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end