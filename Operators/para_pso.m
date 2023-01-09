function Aux = para_pso(varargin)
% Update pbest and gbest for PSO' particle fly operator.

Solution = varargin{1};
Problem  = varargin{2};
Aux      = varargin{3};
Aux.Pbest = update_pairwise([Aux.Pbest,Solution],Problem,'execute');
[~,best]  = min(Aux.Pbest.objs);
Aux.Gbest = Aux.Pbest(best);
end