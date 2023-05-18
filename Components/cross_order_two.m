function [output1,output2] = cross_order_two(varargin)
% Two-order crossover.

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
        Parent = varargin{1};
        Aux    = varargin{4};
        Parent = Parent.decs;
        [N,D] = size(Parent);
        Parent1 = Parent(1:ceil(N/2),:);
        Parent2 = Parent(floor(N/2)+1:end,:);
        Nhalf = size(Parent1,1);

        Offspring1 = Parent1;
        Offspring2 = Parent2;
        for i = 1:Nhalf
            k = randperm(D,2);
            k = sort(k,'ascend');
            Temp = setdiff(Parent2(i,:),Parent1(i,k(1):k(2)),'stable'); % variables of parent 2 that are not appeared in the selected segement (k(1):k(2)) of parent 1
            ind  = setdiff(1:D,k(1):k(2),'stable'); % indexes of variables of parent 2 that are not appeared in the selected segement (k(1):k(2)) of parent 1
            Offspring1(i,ind) = Temp;


            k = randperm(D,2);
            k = sort(k,'ascend');
            Temp = setdiff(Parent1(i,:),Parent2(i,k(1):k(2)),'stable'); % variables of parent 1 that are not appeared in the selected segement (k(1):k(2)) of parent 2
            ind  = setdiff(1:D,k(1):k(2),'stable'); % indexes of variables of parent 1 that are not appeared in the selected segement (k(1):k(2)) of parent 2
            Offspring2(i,ind) = Temp;
        end
        Offspring = [Offspring1;Offspring2];
        output1   = Offspring(1:N,:);
        output2   = Aux;

    case 'parameter'
        % no parameter

    case 'behavior'
        output1 = {'';'GS'};  % always perform global search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end