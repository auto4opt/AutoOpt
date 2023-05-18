function [output1,output2] = search_mu_polynomial(varargin)
% The polynomial mutation.

%------------------------------Reference-----------------------------------
% Deb K, Goyal M. A combined genetic adaptive search (GeneAS) for 
% engineering design[J]. Computer Science and Informatics, 1996, 26: 30-45.
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
        Parent  = varargin{1};
        Problem = varargin{2};
        Para    = varargin{3};
        Aux     = varargin{4};
        ProbM   = Para(1);
        DisM    = Para(2);
        
        if ~isnumeric(Parent)
            Offspring = Parent.decs;           
        else
            Offspring = Parent;
        end    
        [N,D] = size(Offspring);

        Lower = repmat(Problem.bound(1,:),N,1);
        Upper = repmat(Problem.bound(2,:),N,1);

        Site  = rand(N,D) < ProbM;
        mu    = rand(N,D);
        temp  = Site & mu<=0.5;
        Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
            (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(DisM+1)).^(1/(DisM+1))-1);
        temp = Site & mu>0.5;
        Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
            (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(DisM+1)).^(1/(DisM+1)));
        output1 = Offspring;
        output2 = Aux;
    
    case 'parameter'
        output1 = [0,0.3;20,40]; % mutation probability and distribution

    case 'behavior'
        output1 = {'LS','small','large';'GS','large','small'}; % small probabilities and large distributions perform local search

    case 'algorithm'
        Parent    = varargin{1};
        bound     = varargin{2};
        ProbM     = 1;
        DisM      = 30;
        Offspring = Parent;
        [N,D] = size(Offspring);

        Lower = repmat(bound(1,:),N,1);
        Upper = repmat(bound(2,:),N,1);
        Site  = rand(N,D) < ProbM;
        mu    = rand(N,D);
        temp  = Site & mu<=0.5;
        Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
            (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(DisM+1)).^(1/(DisM+1))-1);
        temp = Site & mu>0.5;
        Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
            (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(DisM+1)).^(1/(DisM+1)));
        output1 = Offspring;
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end