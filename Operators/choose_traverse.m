function [output1,output2] = choose_traverse(varargin)
% Select each of the cureent soutions.

mode = varargin{end};
switch mode
    case 'execute'
        Problem  = varargin{2};
        index = 1:Problem.N;
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