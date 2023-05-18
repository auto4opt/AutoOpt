function Code = pseudocode(algorithm,Setting,algInd,state)
% Output the pseudocode of the designed algorithm.

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

% prepare choose
choose    = algorithm.operatorPheno{1}.Choose;
choose    = [choose,'('];
parameter = algorithm.parameterPheno{1}.Choose;
for j = 1:numel(parameter)
    temp   = num2str(parameter(j));
    choose = [choose,temp,','];
end
choose = ['  S = ',choose,'S)'];

% prepare update
update    = algorithm.operatorPheno{1}.Update;
update    = [update,'('];
parameter = algorithm.parameterPheno{1}.Update;
for j = 1:numel(parameter)
    temp   = num2str(parameter(j));
    update = [update,temp,','];
end
update = ['  S = ',update,'S,S_new)'];

% prepare archive
if ~isempty(algorithm.operatorPheno{1}.Archive)
    archive = algorithm.operatorPheno{1}.Archive;
    archive = ['  A = ',archive,'(A,S)'];
end

if Setting.AlgP == 1
    Code = cell(4+7*size(algorithm.operatorPheno{1}.Search,1),1);
    k = 1;
    
    % first row
    switch state
        case 'final'
            if algInd == 1
                Code{k} = 'Best algorithm:'; k = k+1;
            else
                Code{k} =  ['Algorithm ',num2str(algInd),':']; k = k+1;
            end
        case 'iterate'
            Code{k} = ['Algorithm at Iteration ',num2str(algInd),':']; k = k+1;
    end

    % initialize
    Code{k} = 'S = initialize()'; k = k+1;

    % archive
    if ~isempty(algorithm.operatorPheno{1}.Archive)
        Code{k} = archive; k = k+1;
    end

    % while
    Code{k} = 'while algorithm termination condition not met'; k = k+1;

    % for each search operation
    for j = 1:size(algorithm.operatorPheno{1}.Search,1)
        % prepare search
        termninate1 = num2str(algorithm.operatorPheno{1}.Search{j,end}(1));
        termninate2 = num2str(algorithm.operatorPheno{1}.Search{j,end}(2));
        search      = algorithm.operatorPheno{1}.Search{j,1};
        parameter   = algorithm.parameterPheno{1}.Search{j,1};
        search      = [search,'('];
        for l = 1:numel(parameter)
            temp = num2str(parameter(l));
            search = [search,temp,','];
        end
        search = ['  S_new = ',search,'S)'];
        % if is a sexual evolutonary algorithm with crossover and mutation
        if ~isempty(algorithm.operatorPheno{1}.Search{j,2})
            search2    = algorithm.operatorPheno{1}.Search{j,2};
            parameter2 = algorithm.parameterPheno{1}.Search{j,2};
            search2    = [search2,'('];
            for l = 1:numel(parameter2)
                temp = num2str(parameter2(l));
                search2 = [search2,temp,','];
            end
            search2 = ['  S_new = ',search2,'S_new)'];
        end

        % search
        if algorithm.operatorPheno{1}.Search{j,end}(2) == 1 % global search
            Code{k} = choose; k = k+1;
            Code{k} = search; k = k+1;
            if ~isempty(algorithm.operatorPheno{1}.Search{j,2})
                Code{k} = search2; k = k+1;
            end
            Code{k} = update; k = k+1;
            if ~isempty(algorithm.operatorPheno{1}.Archive)
                Code{k} = archive; k = k+1;
            end
        else % ierative local search
            Code{k} = ['  while solution improvement>',termninate1,' or inner_iteration<=',termninate2]; k = k+1;
            Code{k} = ['  ',choose]; k = k+1;
            Code{k} = ['  ',search]; k = k+1;
            if ~isempty(algorithm.operatorPheno{1}.Search{j,2})
                Code{k} = ['  ',search2]; k = k+1;
            end
            Code{k} = ['  ',update]; k = k+1;
            if ~isempty(algorithm.operatorPheno{1}.Archive)
                Code{k} = ['  ',archive]; k = k+1;
            end
            Code{k} = '  end while'; k = k+1;
        end
    end
    Code{k} = 'end while';
else
    Code = cell(7+2*Setting.AlgP,1);
    k = 1;

    % first row
    switch state
        case 'final'
            if algInd == 1
                Code{k} = 'Best algorithm:'; k = k+1;
            else
                Code{k} =  ['Algorithm ',num2str(algInd),':']; k = k+1;
            end
        case 'iterate'
            Code{k} = ['Algorithm at Iteration ',num2str(algInd),':']; k = k+1;
    end

    % initialize
    Code{k} = 'S = initialize()'; k = k+1;

    % archive
    if ~isempty(algorithm.operatorPheno{1}.Archive)
        Code{k} = archive; k = k+1;
    end

    % while
    Code{k} = 'while algorithm termination condition not met'; k = k+1;

    % choose
    Code{k} = choose; k = k+1;

    % prepare subpopulation
    dividePop = '  {S1,';
    for j = 2:Setting.AlgP
        dividePop = [dividePop,'S',num2str(j),','];
    end
    dividePop = [dividePop(1:end-1),'} = S'];
    mergePop = '  S_new = {S1_new,';
    for j = 2:Setting.AlgP
        mergePop = [mergePop,'S',num2str(j),'_new,'];
    end
    mergePop = [mergePop(1:end-1),'}'];

    % divide population
    Code{k} = dividePop; k = k+1;

    % for each search pathway
    for j = 1:Setting.AlgP
        % prepare search
        search      = algorithm.operatorPheno{j}.Search{1};
        parameter   = algorithm.parameterPheno{j}.Search{1};
        search      = [search,'('];
        for l = 1:numel(parameter)
            temp = num2str(parameter(l));
            search = [search,temp,','];
        end
        search = ['  S',num2str(j),'_new = ',search,'S',num2str(j),')'];
        % if is a sexual evolutonary algorithm with crossover and mutation
        if ~isempty(algorithm.operatorPheno{j}.Search{2})
            search2    = algorithm.operatorPheno{j}.Search{2};
            parameter2 = algorithm.parameterPheno{j}.Search{2};
            search2    = [search2,'('];
            for l = 1:numel(parameter2)
                temp = num2str(parameter2(l));
                search2 = [search2,temp,','];
            end
            search2 = ['  S',num2str(j),'_new = ',search2,'S',num2str(j),')'];
        end

        % search
        Code{k} = search; k = k+1;
        if ~isempty(algorithm.operatorPheno{j}.Search{2})
            Code{k} = search2; k = k+1;
        end
    end

    % merge subpopulation
    Code{k} = mergePop; k = k+1;

    % update
    Code{k} = update; k = k+1;

    % archive
    if ~isempty(algorithm.operatorPheno{1}.Archive)
        Code{k} = archive; k = k+1;
    end

    Code{k} = 'end while';
end

% delete empty rows
rowDelete = cellfun(@isempty,Code(:,1));
Code(rowDelete,:) = [];
end