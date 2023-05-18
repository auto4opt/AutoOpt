function [output1,output2] = choose_brainstorm(varargin)
% Brain storm optimization's idea picking up for selecting solutions.

%------------------------------Reference-----------------------------------
% Shi Y. An optimization algorithm based on brainstorming process[M]//
% Emerging Research on Swarm Intelligence and Algorithm Optimization. IGI 
% Global, 2015: 1-35.
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

        K     = max(2,round(Para(1))); % number of clusters
        gamma = Para(2); % probability of selecting cluster centers

        % clustering
        if contains(Problem(1).type{1},'continuous')
            ClusterInd = kmeans(Solution.decs,K,'Distance','sqeuclidean'); % clustering in continuous search space
        else
            error('choose_cluster is only available for continuous problems.'); 
        end

        Cluster = cell(K,1);
        ClusterCenter = zeros(K,1);
        SelRateCenter = zeros(K,1);
        Objs = Solution.objs;
        for i = 1:K
            Cluster{i} = find(ClusterInd == i); % indexes of solutions belonging to cluster i
            [~,ind] = min(Objs(Cluster{i},:)); % fittest solution
            ClusterCenter(i) = Cluster{i}(ind); % record the fittest solution as cluster center
            SelRateCenter(i) = length(Cluster{i})/length(Solution); % select centers according to their scales
            Cluster{i}(ind) = []; % reserve solutions except for the center
        end
        temp = zeros(K,1);
        for i = 1:K
            temp(i) = SelRateCenter(i)./sum(SelRateCenter);
        end
        SelRateCenter = temp;

        % choose
        index = randsrc(Problem(1).N,1,[ClusterCenter';SelRateCenter']); % select N cluster centers from 1:K, some of them may be the same
        for i = 1:Problem(1).N % select cluster center with probability gamma, then replace the center with a random solution from the cluster
            if rand <= gamma && ~isempty(Cluster{ClusterCenter==index(i)})
                ind = randi(numel(Cluster{ClusterCenter==index(i)}));
                index(i) = Cluster{ClusterCenter==index(i)}(ind);
            end
        end
        output1 = index;

    case 'parameter'
        Problem = varargin{1};
        k_max   = max(1,round(Problem(1).N/5)); % maximum number of clusters in brain storm optimization's idea picking up
        output1 = [1,k_max;0,1]; % number of clusters and probability of selecting cluster center in BSO's exploration

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