function [output1,output2] = search_eda(varargin)
% The estimation of distribution.

%------------------------------Reference-----------------------------------
% Baluja S, Caruana R. Removing the genetics from the standard genetic 
% algorithm[M]//Machine Learning Proceedings 1995. Morgan Kaufmann, 1995: 
% 38-46.
%------------------------------Copyright-----------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 
%--------------------------------------------------------------------------

mode = varargin{end};
switch mode
    case 'execute'
        Parent = varargin{1};
        Aux    = varargin{4};

        Parent = SOLVE.decs(Parent);
        [N,D]  = size(Parent);
        
        output1 = zeros(N,D);
        for i = 1:D
            % pd = fitdist(Parent(:,i),'normal'); % fit a normal distribution
            [meanEst, stdEst] = normfit(Parent(:,i));
            % output1(:,i) = random(pd,[N,1]); % sample for the fitted distribution
            output1(:,i) = normrnd(meanEst, stdEst, N, 1);
        end
        output2 = Aux;
  
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