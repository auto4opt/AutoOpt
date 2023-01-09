function [output1,output2]  = search_de_current_best(varargin)
% The "current-to-best/1" differental mutation.

mode = varargin{end};
switch mode
    case 'execute'
        Parent = varargin{1};
        Para   = varargin{3};

        F        = Para(1);
        CR       = Para(2);
        [~,best] = min(Parent.fits);
        GbestDec = Parent(best).dec; % the best solution
        Parent   = Parent.decs;

        [N,D]   = size(Parent);
        Parent1 = Parent;
        Parent2 = repmat(GbestDec,N,1);
        Parent3 = Parent1(randperm(N),:);

        ind = rand(N,D) < CR;
        Offspring = Parent;
        Offspring(ind) = Parent1(ind) + F*(Parent2(ind)-Parent3(ind));
        output1 = Offspring;
        output2 = varargin{5};

    case 'parameter'
        output1 = [0,1;0,1]; % F and CR

    case 'behavior'
        output1 = {'LS','small','small';'GS','large','large'}; % small F and CR values perform local search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end