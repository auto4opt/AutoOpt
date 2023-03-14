function [output1,output2] = choose_roulette_wheel(varargin)
% Roulette wheel selection.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};
        Fitness = Solution.fits;

        % according to fitness rank
        % S = 2;
        % [~,rank] = sort(Fitness); % for single-objective only
        % Prob = (2-S)/N + 2.*rank.*(S-1)./(N*(N-1));

        % according to fitness proportation
        Fitness = reshape(Fitness,1,[]);
        Fitness = Fitness-min(min(Fitness),0)+1e-6;
        Fitness = cumsum(1./Fitness);
        prob    = Fitness./max(Fitness);
        index   = arrayfun(@(S)find(rand<=prob,1),1:Problem.N);
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