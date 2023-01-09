function [output1,output2] = reinit_permutation(varargin)
% Peinitialization for permutation problems

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};

        [N,D] = size(Solution.decs);
        [~,output1] = sort(rand(N,D),2);
        output2 = varargin{5};
    
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