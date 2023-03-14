function [output1,output2] = reinit_discrete(varargin)
% Reinitialization for discrete problems.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};
        Aux      = varargin{4};

        [N,D] = size(Solution.decs);
        output1 = zeros(N,D);
        for j = 1:D
            output1(:,j) = randi([Problem.bound(1,j),Problem.bound(2,j)],N,1);
        end
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