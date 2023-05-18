function output = Select(AllAlgs,Problem,Data,Setting,seedInstance)
% Select promising algorithms.

%----------------------------Copyright-------------------------------------
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

% get algorithms' performance
if strcmp(Setting.Evaluate,'racing')
    for i = 1:length(AllAlgs)
        for j = 1:numel(seedInstance)
            if sum(AllAlgs(i).performance(seedInstance(j),:)) == 0 % if haven't evaluate on instance j
                AllAlgs(i) = AllAlgs(i).Evaluate(Problem,Data,Setting,seedInstance(j));
            end
        end
    end
end
c = AllAlgs.GetPerformance(Setting,seedInstance);

% select algorithms
switch Setting.Evaluate
    case {'exact','approximate'}
        if strcmp(Setting.Compare,'average')
            [~,ind] = sort(c,'ascend'); % small s values refer to better performance
        else
            win = zeros(length(AllAlgs),1);
            for i = 1:length(AllAlgs)
                ind1 = find(c(:,1)==i); % indices of pairwise comparisions that Alg i stands on the first position
                ind2 = find(c(:,2)==i);
                for j = 1:numel(ind1)
                    if c(ind1(j),4) < 0 && c(ind1(j),6) < 0.05 % Alg i's result is better than the opponent && Alg i and the opponent are different
                        win(i) = win(i)+1;
                    end
                end
                for j = 1:numel(ind2)
                    if c(ind2(j),4) > 0 && c(ind2(j),6) < 0.05
                        win(i) = win(i)+1;
                    end
                end
            end
            [~,ind] = sort(win,'descend');
        end
        AllAlgs = AllAlgs(ind(1:Setting.AlgN));
        output = AllAlgs;

    case 'intensification'
        if strcmp(Setting.Compare,'average')
            [~,ind] = sort(c,'ascend');
            deleteInd = [];
            for i = 1:length(ind)
                if ind(i) > length(OldAlgs) && ~any(ind(i+1:end)<=length(OldAlgs))
                    deleteInd = [deleteInd,ind(i)];
                end
            end 
            NewAlgs(deleteInd-length(OldAlgs)) = []; % delete new algorithms that are not better than all incumbents
        else
            rowInd = find(c(:,1)==length(OldAlgs));
            rowInd = rowInd(end);
            c = c(1:rowInd,:); % reserve comparisons between old and new algorithms
            win = zeros(length(AllAlgs),1);
            for i = length(OldAlgs)+1:length(AllAlgs) % for each new algorithm
                ind1 = find(c(:,1)==i); 
                ind2 = find(c(:,2)==i);
                for j = 1:numel(ind1)
                    if c(ind1(j),4) < 0 && c(ind1(j),6) < 0.05 
                        win(i) = win(i)+1;
                    end
                end
                for j = 1:numel(ind2)
                    if c(ind2(j),4) > 0 && c(ind2(j),6) < 0.05
                        win(i) = win(i)+1;
                    end
                end
            end
            win(1:length(OldAlgs),:) = []; % delete information of old algorithms
            NewAlgs(win==0) = []; % delete new algorithms that are not better than all incumbents
        end
        output = NewAlgs;

    case 'racing'
        win = zeros(length(AllAlgs),1);
        for i = 1:length(AllAlgs)
            ind1 = find(c(:,1)==i);
            ind2 = find(c(:,2)==i);
            for j = 1:numel(ind1)
                if c(ind1(j),4) <= 0 && c(ind1(j),6) < 0.05 % <=: not worse than
                    win(i) = win(i)+1;
                end
            end
            for j = 1:numel(ind2)
                if c(ind2(j),4) >= 0 && c(ind2(j),6) < 0.05
                    win(i) = win(i)+1;
                end
            end
        end
        notWorseInd = find(win>0);
        if numel(notWorseInd) < Setting.AlgN
            worseInd    = find(win==0);
            notWorseInd = [notWorseInd;datasample(worseInd,Setting.AlgN-numel(notWorseInd),'Replace',false)];
        end
        AllAlgs = AllAlgs(notWorseInd);
        output = AllAlgs;
end
end