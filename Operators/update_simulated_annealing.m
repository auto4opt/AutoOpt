function [output1,output2] = update_simulated_annealing(varargin)
% Simulated annealing's update mechanism, i.e., accept worse soluion with a
% probability.

mode = varargin{end};
switch mode
    case 'execute'
        Solution = varargin{1};
        Problem  = varargin{2};
        Para     = varargin{3};
        G        = varargin{4};

        T_initial = Para; % initial temperture
        T_final   = 0.01; % final temperture
        Rate      = nthroot(T_final/T_initial,Problem.Gmax); % temperture decrease rate
        T         = T_initial*Rate^G; % current temperture

        Old = Solution(1:Problem.N);
        New = Solution(Problem.N+1:end);
        accept = rand(Problem.N,1) < exp((Old.fits-New.fits)./(abs(Old.fits+1e-6)/T));
        Old(accept) = New(accept);
        output1 = Old;

    case 'parameter'
        output1 = [0.1,1]; % initial temperture 

    case 'behavior'
        output1 = {'';''}; 
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end