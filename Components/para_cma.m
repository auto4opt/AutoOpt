function Aux = para_cma(varargin)
% Update CMA-ES's parameters.

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

Solution = varargin{1};
Problem  = varargin{2};
Aux      = varargin{3};
type     = varargin{4};

% prepare parameters
Disturb = Aux.cma_Disturb;
halfN   = Aux.cma_halfN;
w       = Aux.cma_w;
betterN = Aux.cma_betterN;
mean    = Aux.cma_mean;
sigma   = Aux.cma_sigma;
csigma  = Aux.cma_csigma;
dsigma  = Aux.cma_dsigma;
chiN     = Aux.cma_chiN;
cc      = Aux.cma_cc;
ccov    = Aux.cma_ccov;
cmu     = Aux.cma_cmu;
hth     = Aux.cma_hth;
ps      = Aux.cma_ps;
pc      = Aux.cma_pc;
C       = Aux.cma_C;

% update parameters
switch type
    case 'solution'
        Fitness  = Solution.fits;
    case 'algorithm'
        Fitness  = Solution.avePerformAll;
end
[~,rank] = sort(Fitness,'ascend');
Disturb  = Disturb(rank,:);
DisturbW = w*Disturb(1:halfN,:);
mean     = mean+sigma.*DisturbW;

ps    = (1-csigma)*ps+sqrt(csigma*(2-csigma)*betterN)*DisturbW/chol(C)';
hs    = norm(ps)/sqrt(1-(1-csigma)^(2*(Problem(1).Gmax+1))) < hth;
pc    = (1-cc)*pc+hs*sqrt(cc*(2-cc)*betterN)*DisturbW;
delta = (1-hs)*cc*(2-cc);
C     = (1-ccov-cmu)*C+ccov*(pc'*pc+delta*C);
for i = 1:halfN
    C = C+cmu*w(i)*Disturb(i,:)'*Disturb(i,:);
end
[V,E] = eig(C);
if any(diag(E)<0)
    C = V*max(E,0)/V;
end
sigma = sigma*exp(csigma/dsigma*(norm(ps)/chiN-1))^0.3;

Aux.cma_mean  = mean;
Aux.cma_ps    = ps;
Aux.cma_pc    = pc;
Aux.cma_C     = C;
Aux.cma_sigma = sigma;
end