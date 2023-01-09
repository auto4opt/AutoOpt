function [output1,output2] = search_reset_rand(varargin)
% Reset each element of solutions to a random value with a probability.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};
        Para     = varargin{3};

        Prob = Para;
        if ~isnumeric(Solution)
            New = Solution.decs;
        else
            New = Solution;
        end
        [N,D] = size(New);

        Temp = zeros(N,D);
        for j = 1:D
            Temp(:,j) = randi([Problem.bound(1,j),Problem.bound(2,j)],N,1);
        end

        ind = rand(N,D) < Prob; % N*D
        New(ind) = Temp(ind);
        output1 = New;
        output2 = varargin{5};

    case 'parameter'
        output1 = [0,0.5]; % reset probability 

    case 'behavior'
        output1 = {'LS','small';'GS','large'}; % samll probabilities perform local search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end