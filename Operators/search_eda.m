function [output1,output2] = search_eda(varargin)
% The estimation of distribution.

mode = varargin{end};
switch mode
    case 'execute'
        Parent = varargin{1};

        Parent = Parent.decs;
        [N,D]  = size(Parent);
        
        output1 = zeros(N,D);
        for i = 1:D
            pd = fitdist(Parent(:,i),'normal'); % fit a normal distribution
            output1(:,i) = random(pd,[N,1]); % sample for the fitted distribution
        end
        output2 = varargin{5};
  
    case 'parameter'
        % n/a

    case 'behavior'
        output1 = {'';'GS'}; % always performs global search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end