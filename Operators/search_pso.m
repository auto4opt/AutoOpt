function [output1,output2] = search_pso(varargin)
% Particle swarm optimization's particle fly operator

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Para     = varargin{3};
        Aux      = varargin{5};

        [N,D] = size(Solution.decs);
        W = Para; % inertia weight

        % initialize pbest, gbest, and velocity
        if ~isfield(Aux,'Pbest')
            Aux.Pbest = Solution; % personal best solutions
            Aux.Gbest = Solution(randi(N)); % global best solution
            Aux.V = zeros(N,1);
        end

        % particle fly
        ParticleDec = Solution.decs;
        PbestDec = Aux.Pbest.decs;
        GbestDec = Aux.Gbest.dec;

        r1 = repmat(rand(N,1),1,D);
        r2 = repmat(rand(N,1),1,D);
        V = W.*Aux.V + 2*r1.*(PbestDec-ParticleDec) + 2*r2.*(GbestDec-ParticleDec); % N*D
        output1 = ParticleDec + V; % N*D

        % update volecity
        Aux.V   = V;
        output2 = Aux;

    case 'parameter'
        output1 = [0,0.5]; % inertia weight

    case 'behavior'
        output1 = {'';'GS'}; % always perform global search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end