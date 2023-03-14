function [output1,output2] = update_round_robin(varargin)
% Conduct update via round-robin tournament.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};

        Fitness  = Solution.fits;
        K   = min(10,Problem.N-1);
        ind = 1:length(Solution);
        win = zeros(length(Solution),1);
        for i = 1:length(Solution)
            win(i) = sum(Fitness(i)<=Fitness(randperm(end,K)));
        end
        [~,rank] = sort(win,'descend');
        ind = ind(rank(1:Problem.N));
        output1 = Solution(ind);

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