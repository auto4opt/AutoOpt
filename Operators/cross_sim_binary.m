function [output1,output2] = cross_sim_binary(varargin)
% Simulated binary crossover.

mode = varargin{end};
switch mode
    case 'execute'
        Parent = varargin{1};
        Para   = varargin{3};

        Parent = Parent.decs;
        DisC = Para;
            
        [N,D] = size(Parent);
        Parent1 = Parent(1:ceil(N/2),:);
        Parent2 = Parent(floor(N/2)+1:end,:);
        Nhalf = size(Parent1,1);

        beta = zeros(Nhalf,D);
        mu = rand(Nhalf,D);

        beta(mu<=0.5) = (mu(mu<=0.5)*2).^(1/(1+DisC));
        beta(mu>0.5)  = (2-mu(mu>0.5)*2).^(-1/(1+DisC));
        Offspring = [0.5*((1+beta).*Parent1+(1-beta).*Parent2);
            0.5*((1-beta).*Parent1+(1+beta).*Parent2)];        
        output1   = Offspring(1:N,:);
        output2   = varargin{5};

    case 'parameter'
        output1 = [20,40]; % crossover distribution

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