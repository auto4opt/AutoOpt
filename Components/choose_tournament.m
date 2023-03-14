function [output1,output2] = choose_tournament(varargin)
% K-tournament selection.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};
        Fitness = Solution.fits;
        K = 2;

        MatchIndex = randi(size(Fitness,1),Problem.N,K);
        Fitness = repmat(Fitness,1,K);
        MatchFitness = Fitness(MatchIndex);
        index = MatchIndex(MatchFitness == min(MatchFitness,[],2));
        index = index(1:Problem.N);
        output1 = index;

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