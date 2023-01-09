function [output1,output2] = reinit_continuous(varargin)
% Reinitialization for continuous problems.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};

        N = length(Solution);
        Lower = Problem.bound(1,:);
        Upper = Problem.bound(2,:);
        output1 = unifrnd(repmat(Lower,N,1),repmat(Upper,N,1));
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