function [output1,output2] = archive_diversity(varargin)
% Collect the N most diversified solutions found so far to an external
% archive.

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
        Solution    = varargin{1};
        CurrArchive = varargin{2};
        Problem     = varargin{3};

        Solution = [Solution,CurrArchive];

        % distance between each pair of solutions
        if contains(Problem.type{1},'continuous')
            dist = pdist2(Solution.decs,Solution.decs,'euclidean'); % n*n, n is the number of solutions
        else
            dist = pdist2(Solution.decs,Solution.decs,'hamming');
        end

        ind = [randi(length(Solution));zeros(Problem.N-1,1)];
        for i = 2:Problem.N
            % sort the total distance between each solution and the archive solutions in descending order
            [~,rank] = sort(sum(dist(ind(1:i-1),:),1),'descend');
            j = 1;
            ind(i) = rank(j); % select the solution with the largest distance to the archive solutions
            while ismember(ind(i),ind(1:i-1))
                j = j+1;
                ind(i) = rank(j);
            end
        end
        output1 = Solution(ind);
        
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