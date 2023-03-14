function [output1,output2] = search_insert(varargin)
% Randomly select two elements, then insert the second element to the
% position next (right side) to the first element.

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
            CurrentNew = New(i,:);
            k = randperm(D,2); % the two points particulated in insertion should be different
            k = sort(k,'ascend');
            Temp = CurrentNew(k(2));
            CurrentNew(k(2)) = [];
            New(i,:) = [CurrentNew(1:k(1)),Temp,CurrentNew(k(1)+1:end)];
        end
        output1 = New;
        output2 = Aux;
    
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