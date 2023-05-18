function [output1,output2] = choose_nich(varargin)
% Adaptive niching based on the nearest-better clustering. Crossover 
% (if available) will be conducted between solutions from the same specie.

%------------------------------Reference-----------------------------------
% Preuss M. Niching the CMA-ES via nearest-better clustering[C]//
% Proceedings of the 12th annual conference companion on Genetic and 
% evolutionary computation. 2010: 1711-1718.

% Yan B, Zhao Q, Li M, et al. Fitness landscape analysis and niching 
% genetic approach for hybrid beamforming in RIS-aided communications[J].
% Applied Soft Computing, 2022, 131: 109725.
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
        G        = varargin{5};
        
        disance = pdist2(Solution.decs,Solution.decs); % compute distance
        species = mNBC(disance,1,-1,G,Problem.Gmax); % divide species
        
        index = [];
        for i = 1:length(species)
            currInd = species(i).idx;
            currInd = currInd(randperm(numel(currInd)));
            index = [index;currInd]; % gather all solutions' indexes 
        end
        
        if mod(length(index),2) == 1 % if the number of solutions is odd 
            odd = true;
            index = [index;index(end)];
        else
            odd = false;
        end

        % crossover between solutions from the same specie
        index = reshape(index,2,length(index)/2);
        index = index';
        index = [index(:,1);index(:,2)];
        if odd == true
            index(end) = [];
        end

        output1 = index;

    case 'parameter'
        % n/a

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

function species = mNBC(disance,fai,min_size,g,gmax)
% matdis: the distance between each pair of the individuals
% fai: the weight in NBC
% min_size: the min size of the species
% D: the dimension ofthe problem
% g: the current generation

    % reset 'min_size' if 'min_size' is equal to -1
    if min_size == -1
        min_size = round(3+g/gmax*7);
    end 
    
    % find the nearest better neighour and the distance of each individual
    NP = size(disance, 1);
    nbc = zeros(NP, 3);
    nbc(1, :) = [1 -1 0]; % the best individual do not have the nearest better neighbour
    for i = 2:NP
        nbc(i, 1) = i;
        [nbc(i, 3), nbc(i, 2)] = min(disance(i, 1:i-1));
    end
    
    % set the follow value
    follow = ones(NP, 1);
    for i = NP:-1:2 
        follow(nbc(i, 2)) = follow(nbc(i, 2)) + follow(i); % the subtree rooted at nearest better neighbour individual must contain all current individual's nodes
    end
    
    % cut the edge from the longest to the shortest
    meandis = fai * mean(nbc(2:NP, 3));
    seeds = 1;
    
    [~, sort_index] = sort(nbc(:, 3), 'descend');
    for i = 1:NP
        if nbc(sort_index(i), 3) > meandis % one of the cut conditions
            inf_index = sort_index(i); % the inferior individual
            sup_index = nbc(sort_index(i), 2); % the superior individual
            top_index = sup_index; % the root of the subtree which contains the inferior and superior individuals
            while nbc(top_index, 2) ~= -1
                top_index = nbc(top_index, 2);
                sup_index = [sup_index, top_index];
            end
            
            if follow(inf_index) >= min_size && follow(top_index) - follow(inf_index) >= min_size % the other of the cut conditions
                % cut operator
                nbc(inf_index, 2) = -1;
                nbc(inf_index, 3) = 0;
                % put the current seed into the set
                seeds = [seeds; sort_index(i)];
                follow(sup_index) = follow(sup_index) - follow(inf_index);
            end
        end
    end
    
    % set the root of subtree which contain the current individual
    m = zeros(NP, 2);
    m(1:NP, 1) = 1:NP;
    for i = 1:NP
       j = nbc(i, 2);
       k = j;
       while j ~= -1
           k =j;
           j = nbc(j, 2);
       end
       if k == -1
           m(i, 2) = i;
       else
           m(i, 2) = k;
       end
    end
    
    % construct the result
    species = struct();
    for i=1:length(seeds)
       species(i).seed = seeds(i);
       species(i).idx = m(m(:, 2) == seeds(i), 1);
       species(i).len = length(species(i).idx);
    end  
end