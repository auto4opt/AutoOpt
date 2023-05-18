function [output1,output2] = search_swap_multi(varargin)
% Swap each pair of elements between two randomly selected indices.

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
        Aux      = varargin{4};

        if ~isnumeric(Solution)
            Solution = Solution.decs;
        end
        [N,D] = size(Solution);
        New = Solution;

        for i = 1:N
            k = randperm(D,2);
            k = sort(k);
            n = floor((k(2)-k(1)+1)/2); % number of pairs of elements to be swapped

            for j = 1:n
                New(i,k(1)+j-1) = Solution(i,k(2)-j+1);
                New(i,k(2)-j+1) = Solution(i,k(1)+j-1);
            end
        end
        output1 = New;
        output2 = Aux;

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