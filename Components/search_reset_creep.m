function [output1,output2] = search_reset_creep(varargin)
% Add a small value (positive or negative) to each element of solutions
% with a probability, for discrete problems with ordinal attributes.

%------------------------------Copyright-----------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 

% Please reference the paper below if using AutoOptLib in your publication:
% @article{zhao2023autooptlib,
%  title={AutoOptLib: A Library of Automatically Designing Metaheuristic 
%         Optimization Algorithms in Matlab},
%  author={Zhao, Qi and Yan, Bai and Hu, Taiwei and Chen, Xianglong and 
%          Yang, Jian and Shi, Yuhui},
%  journal={arXiv preprint 	arXiv:2303.06536},
%  year={2023}
% }
%--------------------------------------------------------------------------

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};
        Para     = varargin{3};
        Aux      = varargin{4};

        Prob = Para(1);
        Amp  = Para(2);
        Amp  = max(1,round((Problem.bound(2,:)-Problem.bound(1,:)+1).*Amp)); % 1*D
        if ~isnumeric(Solution)
            New = Solution.decs;
        else
            New = Solution;
        end
        [N,D] = size(New);

        StepSize = zeros(N,D);
        for i = 1:D
            StepSize(:,i) = randi(Amp(i),N,1); % sample step size for each dimension, N*D
        end
        ind = rand(N,D) < 0.5; % N*D logical indices
        StepSize(ind) = -StepSize(ind);
        Temp = New+StepSize;
        Lower = repmat(Problem.bound(1,:),N,1);
        Upper = repmat(Problem.bound(2,:),N,1);
        Temp = max(min(Temp,Upper),Lower);

        ind = rand(N,D) < Prob; % N*D
        New(ind) = Temp(ind);
        output1 = New;
        output2 = Aux;

    case 'parameter'
        output1 = [0,0.5;0,0.5]; % probability and amplitude 

    case 'behavior'
        output1 = {'LS','small','small';'GS','large','large'}; % small probabilities and amplitudes perform local search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end