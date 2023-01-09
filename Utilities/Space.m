function [AllOp,OpSpace,ParaSpace,ParaLocalSpace,BehavSpace] = Space(Problem,Setting)
% Define the design space.

% operater space
switch Problem(1).type{1}
    case 'continuous'
        Choose  = {'choose_traverse';'choose_tournament';'choose_roulette_wheel';'choose_cluster'}; 
        Search  = {'search_pso';'search_de_current';'search_de_current_best';'search_de_random';
            'cross_arithmetic';'cross_sim_binary';'cross_point_one';'cross_point_two';
            'cross_point_uniform';'search_mu_gaussian';'search_mu_cauchy';'search_mu_polynomial';
            'search_mu_uniform';'search_eda';'search_cma';'reinit_continuous'};
        Update  = {'update_greedy';'update_round_robin';'update_pairwise';'update_always';'update_simulated_annealing'};
       
    case 'discrete'
        Choose  = {'choose_traverse';'choose_tournament';'choose_roulette_wheel'};
        Search  = {'cross_point_one';'cross_point_two';'cross_point_uniform';'search_reset_one';
            'search_reset_rand';'reinit_discrete'};
        Update  = {'update_greedy';'update_round_robin';'update_pairwise';'update_always';'update_simulated_annealing'};

    case 'permutation'
        Choose  = {'choose_traverse';'choose_tournament';'choose_roulette_wheel'};
        Search  = {'cross_order_two';'cross_order_n';'search_swap';'search_swap_multi';
            'search_scramble';'search_insert';'reinit_permutation'};
        Update  = {'update_greedy';'update_round_robin';'update_pairwise';'update_always';'update_simulated_annealing'};
end
OpSpace = [1,length(Choose);
    length(Choose)+1,length(Choose)+length(Search);
    length(Choose)+length(Search)+1,length(Choose)+length(Search)+length(Update)]; % each row reports the search space of one kind of operaters, e.g., the search space of Choose is intergers from 1 to 4.
AllOp = [Choose;Search;Update];

% parameter and behavior space
ParaSpace      = cell(length(AllOp),1); % parameter space
ParaLocalSpace = cell(length(AllOp),1); % parameter space for performing local search
BehavSpace     = cell(length(AllOp),1); % search behaviors
for i = 1:length(AllOp)
    [thisParaSpace,~]  = feval(str2func(AllOp{i}),Problem,'parameter');
    [thisBehavSpace,~] = feval(str2func(AllOp{i}),'behavior');
    ParaSpace{i} = thisParaSpace;
    BehavSpace{i} = thisBehavSpace;
    if ~isempty(thisBehavSpace{1,1}) && ~isempty(thisParaSpace) % if the operator can perform local search and has parameter(s)
        if isempty(thisBehavSpace{2,1}) % if only performs local search
            thisParaLocalSpace = thisParaSpace;
        else
            thisParaLocalSpace = zeros(size(thisParaSpace));
            for j = 1:size(thisParaLocalSpace,1) % for each parameter of operator i
                lower = thisParaSpace(j,1);
                upper = thisParaSpace(j,2);
                trend = thisBehavSpace{1,j+1}; % parameter setting for local search
                if strcmp(trend,'small')
                    thisParaLocalSpace(j,:) = [lower,lower+(upper-lower)*Setting.LSRange];
                elseif strcmp(trend,'large')
                    thisParaLocalSpace(j,:) = [upper-upper*Setting.LSRange,upper];
                end
            end
        end
        ParaLocalSpace{i} = thisParaLocalSpace;
    end
end