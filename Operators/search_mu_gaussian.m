function [output1,output2] = search_mu_gaussian(varargin)
% The Gaussian mutation.

mode = varargin{end};
switch mode
    case 'execute'
        Parent  = varargin{1};
        Aux     = varargin{5};
        innerG  = varargin{6};

        if ~isnumeric(Parent)
            Parent = Parent.decs;
        end
        [N,D] = size(Parent);
        
        % initialize eta
        if innerG == 1
            Aux.gaussian_eta = rand(N,D); 
        end
     
        % search
        Disturb = Aux.gaussian_eta.*randn(N,D);
        output1 = Parent+Disturb;

        % update eta
        tau1      = 1/sqrt(2*sqrt(D));
        tau2      = 1/sqrt(2*D);
        normal    = repmat(randn(N,1),1,D);
        normal_j  = randn(N,D);
        Aux.gaussian_eta = Aux.gaussian_eta.*exp(tau2*normal+tau1*normal_j);
        output2   = Aux;


    case 'parameter'
        % n/a

    case 'behavior'
        output1 = {'';'GS'}; % always performs global search
end

if ~exist('output1','var')
    output1 = [];
end
if ~exist('output2','var')
    output2 = [];
end
end