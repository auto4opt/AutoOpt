function [output1,output2] = search_cma(varargin)
% The evolution strategy with convariance matrix adaption.

%------------------------------Reference-----------------------------------
% N. Hansen and A. Ostermeier, Completely derandomized selfadaptation in
% evolution strategies, Evolutionary Computation, 2001, 9(2): 159-195.
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
        lower   = Problem.bound(1,:);
        upper   = Problem.bound(2,:);
        Aux     = varargin{4};
        if ~isnumeric(Parent)
            Parent  = Parent.decs;
        end
        [N,D] = size(Parent);

        % initialize parameters
        if ~isfield(Aux,'cma_halfN')
            Aux.cma_halfN   = round(N/2);
            w               = log(Aux.cma_halfN+0.5)-log(1:Aux.cma_halfN);
            Aux.cma_w       = w./sum(w);
            Aux.cma_betterN = 1/sum(Aux.cma_w.^2);

            Aux.cma_csigma  = (Aux.cma_betterN+2)/(D+Aux.cma_betterN+5);
            Aux.cma_dsigma  = Aux.cma_csigma+2*max(sqrt((Aux.cma_betterN-1)/(D+1))-1,0)+1;
            Aux.cma_chiN     = sqrt(D)*(1-1/(4*D)+1/(21*D^2));

            Aux.cma_cc    = (4+Aux.cma_betterN/D)/(4+D+2*Aux.cma_betterN/D);
            Aux.cma_ccov  = 2/((D+1.3)^2+Aux.cma_betterN);
            Aux.cma_cmu   = min(1-Aux.cma_ccov,2*(Aux.cma_betterN-2+1/Aux.cma_betterN)/((D+2)^2+2*Aux.cma_betterN/2));
            Aux.cma_hth   = (1.4+2/(D+1))*Aux.cma_chiN;

            Aux.cma_mean  = unifrnd(lower,upper);
            Aux.cma_ps    = zeros(1,D);
            Aux.cma_pc    = zeros(1,D);
            Aux.cma_C     = eye(D);
            Aux.cma_sigma = 0.1*(upper-lower);
        end

        % load parameters
        C     = Aux.cma_C;
        mean  = Aux.cma_mean;
        sigma = Aux.cma_sigma;

        % search
        Disturb = zeros(N,D);
        for i = 1:N
            Disturb(i,:) = mvnrnd(zeros(1,D),C);
        end
        output1 = mean+sigma.*Disturb;
        Aux.cma_Disturb = Disturb;
        output2 = Aux;

    case 'parameter'
        % n/a

    case 'behavior'
        output1 = {'';'GS'}; % always performs global search

    case 'algorithm'
        Parent = varargin{1};
        bound  = varargin{2};
        Aux    = varargin{3};
        lower  = bound(1,:);
        upper  = bound(2,:);
        [N,D] = size(Parent);

        % initialize parameters
        if ~isfield(Aux,'cma_halfN')
            Aux.cma_halfN   = round(N/2);
            w               = log(Aux.cma_halfN+0.5)-log(1:Aux.cma_halfN);
            Aux.cma_w       = w./sum(w);
            Aux.cma_betterN = 1/sum(Aux.cma_w.^2);

            Aux.cma_csigma  = (Aux.cma_betterN+2)/(D+Aux.cma_betterN+5);
            Aux.cma_dsigma  = Aux.cma_csigma+2*max(sqrt((Aux.cma_betterN-1)/(D+1))-1,0)+1;
            Aux.cma_chiN     = sqrt(D)*(1-1/(4*D)+1/(21*D^2));

            Aux.cma_cc    = (4+Aux.cma_betterN/D)/(4+D+2*Aux.cma_betterN/D);
            Aux.cma_ccov  = 2/((D+1.3)^2+Aux.cma_betterN);
            Aux.cma_cmu   = min(1-Aux.cma_ccov,2*(Aux.cma_betterN-2+1/Aux.cma_betterN)/((D+2)^2+2*Aux.cma_betterN/2));
            Aux.cma_hth   = (1.4+2/(D+1))*Aux.cma_chiN;

            Aux.cma_mean  = unifrnd(lower,upper);
            Aux.cma_ps    = zeros(1,D);
            Aux.cma_pc    = zeros(1,D);
            Aux.cma_C     = eye(D);
            Aux.cma_sigma = 0.1*(upper-lower);
        end

        % load parameters
        C     = Aux.cma_C;
        mean  = Aux.cma_mean;
        sigma = Aux.cma_sigma;

        % search
        Disturb = zeros(N,D);
        for i = 1:N
            Disturb(i,:) = mvnrnd(zeros(1,D),C);
        end
        output1 = mean+sigma.*Disturb;
        Aux.cma_Disturb = Disturb;
        output2 = Aux;
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end