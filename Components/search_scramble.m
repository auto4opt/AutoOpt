function [output1,output2] = search_scramble(varargin)
% Scramble all the elements between two randomly selected indices of each
% solution.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Aux      = varargin{4};

        if ~isnumeric(Solution)
            New = Solution.decs;
        else
            New = Solution;
        end
        [N,D] = size(New);

        for i = 1:N
            k = randperm(D,2);
            k = sort(k);
            n = k(2)-k(1)+1; % number of elements to be scarmbled

            Temp = New(i,k(1):k(2));
            ind = randperm(n);
            Temp = Temp(ind);
            New(i,k(1):k(2)) = Temp;
        end
        output1 = New;
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