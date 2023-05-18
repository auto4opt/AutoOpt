% Class for estimating the designed algorithm's performance.

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
classdef Approximate < handle
    properties
        data;
        embedding;
        model;
        randSeed;
        exactG;
    end

    methods
        function obj = Approximate(Problem,Data,Setting,indInstance)
            if nargin > 0
                obj = Approximate;
                % seed for disrupting the order of elements in the vector representation of algorithms
                if Setting.AlgP == 1 % one search pathway
                    obj.randSeed = randperm((Setting.AlgQ+2)*3);
                else % multiple search pathways
                    obj.randSeed = randperm((Setting.AlgP+2)*3);
                end

                % determine iterations with exact performance estimations
                ExactGmax  = Setting.Surro/Setting.AlgN;
                AlgGmax    = ceil(Setting.AlgFE/Setting.AlgN);
                obj.exactG = 1:AlgGmax/ExactGmax:AlgGmax;

                % get training data of surrogate
                TrainAlgs1 = DESIGN(Problem,Setting,500);
                obj.data = TrainAlgs1.Evaluate(Problem,Data,Setting,indInstance); % exactly evaluate performance

                % train embedding
                TrainAlgs2 = DESIGN(Problem,Setting,1000);
                obj.embedding = obj.GetEmbed(TrainAlgs2,Setting);
                
                % train surrogate
                obj.model = obj.GetModel(obj.data,Setting);
            end
        end

        function EmbedMap = GetEmbed(obj,data,Setting)
            EmbedMap  = Embedding(data,Setting,obj,'get');
        end

        function EmbedAlgs = UseEmbed(obj,data,Setting)
            EmbedAlgs = Embedding(data,Setting,obj,'use');
        end

        function Model = GetModel(obj,data,Setting)
            % train a random forest surrogate to estimate algorithms' performance
            EmbedAlgs = Embedding(data,Setting,obj,'use'); % get the embeded algorithms
            Labels    = data.avePerformAll; % all algorithms' average performance on all instances
            Model     = TreeBagger(1000,EmbedAlgs,Labels,'Method','regression','MinLeafSize',5); % train a random forest surrogate
        end

        function obj = UpdateModel(obj,Algs,Setting)
            % Update the surrogate model during the design process.
            N = Setting.AlgN;
            a = Algs.avePerformAll;
            b = Algs.avePerformApproxAll;
            total  = 0;
            error  = 0;
            indNew = [];
            for i = 1:N
                for j = i+1:N
                    res = xor(a(i)<a(j),b(i)<b(j)); % res=xor(A,B)=1 if A is different with B; res=0 otherwise
                    if res == 1
                        error = error+1;
                        indNew = [indNew,i,j];
                    end
                    total = total+1;
                end
            end
            newData = Algs(unique(indNew));
            obj.data  = [obj.data,newData];
            obj.model = obj.GetModel(obj.data,Setting);
        end
    end
end