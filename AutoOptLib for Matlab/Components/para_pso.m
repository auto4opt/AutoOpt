function Aux = para_pso(varargin)
% Update pbest and gbest for PSO' particle fly operator.

%------------------------------Reference-----------------------------------
% Shi Y, Eberhart R. A modified particle swarm optimizer[C]//1998 IEEE 
% International Conference on Evolutionary Computation Proceedings. IEEE 
% World Congress on Computational Intelligence (Cat. No. 98TH8360). IEEE, 
% 1998: 69-73.
%------------------------------Copyright-----------------------------------
% Copyright (C) <2023>  <Swarm Intelligence Lab>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 
%--------------------------------------------------------------------------

Solution = varargin{1};
Problem  = varargin{2};
Aux      = varargin{3};
Aux.Pbest = update_pairwise([Aux.Pbest,Solution],Problem,'execute');
[~,best]  = min(Aux.Pbest.objs);
Aux.Gbest = Aux.Pbest(best);
end