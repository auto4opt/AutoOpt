function output = Select(OldAlgs,NewAlgs,Problem,Data,Setting,indInstance)
% Select promising algorithms.

% get algorithms' performance
Algs = [OldAlgs,NewAlgs];
if strcmp(Setting.Evaluate,'racing')
    for i = 1:length(Algs)
        for j = 1:numel(indInstance)
            if sum(Algs(i).performance(indInstance(j),:)) == 0 % if haven't evaluate on instance j
                Algs(i) = Design_Evaluate(Algs(i),Problem,Data,Setting,indInstance(j));
            end
        end
    end
end
c = Algs.GetPerformance(Setting,indInstance);

% select algorithms
switch Setting.Evaluate
    case {'default','approximate'}
        if strcmp(Setting.Compare,'average')
            [~,ind] = sort(c,'ascend'); % small s values refer to better performance
        else
            win = zeros(length(Algs),1);
            for i = 1:length(Algs)
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
        Algs = Algs(ind(1:Setting.AlgN));
        output = Algs;

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
            win = zeros(length(Algs),1);
            for i = length(OldAlgs)+1:length(Algs) % for each new algorithm
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
        win = zeros(length(Algs),1);
        for i = 1:length(Algs)
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
        Algs = Algs(notWorseInd);
        output = Algs;
end
end