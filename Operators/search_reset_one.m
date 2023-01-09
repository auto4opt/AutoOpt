function [output1,output2] = search_reset_one(varargin)
% Randomly select an element of each solution and reset it to a random
% value.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};

        if ~isnumeric(Solution)
            New = Solution.decs;
        else
            New = Solution;
        end

        [N,D] = size(New);
        ind = randi(D,N,1); % N*1
        
        for i = 1:N
            curr = New(i,ind(i));
            while New(i,ind(i)) == curr
                New(i,ind(i)) = randi([Problem.bound(1,ind(i)),Problem.bound(2,ind(i))]);
            end
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