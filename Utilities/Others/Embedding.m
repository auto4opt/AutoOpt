function output = Embedding(Algs,Setting,Surrogate,mode)
% Embed the designed algorithm into a compact vector representation by 
% auto-encoder.

%----------------------------Copyright-------------------------------------
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

randSeed = Surrogate.randSeed;
layer    = 2;   % number of auto-encoder layers
p        = 0.3; % corrupt probability
lower    = 1;
upper    = size(Algs(1).parameter,1);

% get vector representations of the designed algorithms
if Setting.AlgP == 1 % one search pathway
    for i = 1:length(Algs)
        operator = zeros(1,Setting.AlgQ+2); % vector representation of operators
        operator(1) = Algs(i).operator{1}(1,1);
        operator(2:size(Algs(i).operator{1},1)) = Algs(i).operator{1}(2:end,1);
        operator(end) = Algs(i).operator{1}(end,2);

        parameter = zeros(2,Setting.AlgQ+2); % vector representation of parameters
        for j = 1:Setting.AlgQ+2
            if operator(j) ~= 0
                for k = 1:numel(Algs(i).parameter{operator(j),1})
                    parameter(k,j) = Algs(i).parameter{operator(j),1}(k);
                end
            end
        end
        parameter = reshape(parameter,1,(Setting.AlgQ+2)*2);

        operator = (operator-min(operator))./(max(operator)-min(operator)); % normilize to [0,1]
        tempAlg  = [operator,parameter];
        tempAlg  = tempAlg(randSeed); % disrupt the order of elements in the vector representation
        tempAlg  = reshape(tempAlg,numel(randSeed)/3,3); % for reducting dimensionality of the vector representation
        VectorAlgs(:,3*i-2:3*i) = tempAlg;
    end
else % multiple search pathways
    for i = 1:length(Algs)
        operator = zeros(1,Setting.AlgP+2);
        operator(1) = Algs(i).operator{1}(1,1);
        for j = 1:Setting.AlgP
            operator(j+1) = Algs(i).operator{j}(2,1);
        end
        operator(end) = Algs(i).operator{1}(end,2);
        
        parameter = zeros(2,Setting.AlgP+2); % each operator has at most 2 parameters
        for j = 1:numel(Algs(i).parameter{operator(1),1}) % choose operator's parameter(s)
            parameter(j,1) = Algs(i).parameter{operator(1),1}(j);
        end
        for j = 1:Setting.AlgP % search operators' parameter(s)
            for k = 1:numel(Algs(i).parameter{operator(j+1),1})
                parameter(k,j+1) = Algs(i).parameter{operator(j+1),1}(k);
            end
        end
        for j = 1:numel(Algs(i).parameter{operator(end),1}) % update operator's parameter(s)
            parameter(j,end) = Algs(i).parameter{operator(end),1}(j);
        end
        parameter = reshape(parameter,1,(Setting.AlgP+2)*2);

        operator = (operator-lower)./(upper-lower); % normilize to [0,1]
        tempAlg  = [operator,parameter];
        tempAlg  = tempAlg(randSeed); % disrupt the order of elements in the vector representation
        tempAlg  = reshape(tempAlg,numel(randSeed)/3,3); % for reducting dimensionality of the vector representation
        VectorAlgs(:,3*i-2:3*i) = tempAlg;
    end
end

switch mode
    case 'get' % get the embedding mapping by mSDA
        [EmbedMap,~] = mSDA(VectorAlgs,p,layer);
        output = EmbedMap; 
    case 'use' % use the mapping to embed algorithms
        EmbedMap = Surrogate.embedding;
        temp     = [VectorAlgs;ones(1,length(Algs)*3)];
        wx       = EmbedMap(:,:,1)*temp;
        argslist = zeros(numel(randSeed)/3,length(Algs));
        for i = 1:length(Algs)
            argslist(:,i) = wx(:,3*i-2)+wx(:,3*i-1)+wx(:,3*i);
        end
        wx = argslist/3;
        wx = tanh(wx);
        layer = size(EmbedMap,3);
        for i = 2:layer
            wx = [wx;ones(1,length(Algs))];
            wx = EmbedMap(:,:,i)*wx;
            wx = tanh(wx);
        end
        EmbedAlgs = wx';
        output = EmbedAlgs;
end
end

function [Ws,hs] = mSDA(X,p,l)
[d,n] = size(X);
n  = n/3;
Ws = zeros(d,d+1,l);
hs = zeros(d,n,l+1);
%hs(:,:,1) = X;
[Ws(:,:,1),hs(:,:,1+1)] = mDA_de(X,p);
for t = 2:l
    [Ws(:,:,t),hs(:,:,t+1)] = mDA(hs(:,:,t),p);
end
end

function [W,h] = mDA_de(X,p)
X = [X;ones(1,size(X,2))];
d = size(X,1);
q = [ones(d-1,1).*(1-p); 1];
S = X*X';
Q = S.*(q*q');
Q(1:d+1:end) = q.*diag(S);
P = S.*repmat(q',d,1);
W = P(1:end-1,:)/(Q+1e-5*eye(d));

length = size(X,2)/3;
h = zeros(d-1,length);
wx = W*X;
for i = 1:length
    temp = wx(:,3*i-2:3*i);
    temp = temp(:,1)+temp(:,2)+temp(:,3);
    h(:,i) = temp/3;
end
h = tanh(h);
%h = tanh(W*X);
end

function [W,h] = mDA(X,p)
X = [X;ones(1,size(X,2))];
d = size(X,1);
q = [ones(d-1,1).*(1-p); 1];
S = X*X';
Q = S.*(q*q');
Q(1:d+1:end) = q.*diag(S);
P = S.*repmat(q',d,1);
W = P(1:end-1,:)/(Q+1e-5*eye(d));
h = tanh(W*X);
end
