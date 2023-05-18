function [output1,output2,output3] = prob_name(varargin)
% Template for writing problem file. 

%----------------------------Copyright-------------------------------------
% Copyright (C) <2023>  <Qi Zhao>

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

switch varargin{end}
    case 'construct' % define problem properties
        Problem = varargin{1};
        % define problem type in the following three cells.
        % first cell : 'continuous'\'discrete'\'permutation'
        % second cell: 'static'\'sequential'
        % third cell : 'certain'\'uncertain'
        Problem.type = {'','',''}; 
        
        % define the bound of solution space
        lower = []; % 1*D, lower bound of the D-dimension decision space
        upper = []; % 1*D, upper bound of the D-dimension decision space
        Problem.bound = [lower;upper];        
        
        % define specific settings (optional), choices: 
        % 'dec_diff': elements of the solution should be different w.r.t 
        %             each other for discrete problems
        % 'uncertain_average': averaging the fitness over multiple fitness 
        %             evaluations for uncertain problems
        % 'uncertain_worst'  : use the worse fitness among multiple fitness
        %             evaluations as the fitness for uncertain problems
        Problem.setting = {''}; % put choice(s) into the cell        
        
        % set the number of samples for uncertain problems (optional)
        Problem.sampleN = [];
               
        output1 = Problem;
        
        % load/construct data file
        Data = load(''); % for .mat format
        Data = readmatrix('','Sheet',1); % for .xlsx format
        output2 = Data;

    case 'repair' % repair solutions
        Data = varargin{1};
        Decs = varargin{2};

        % define methods for repairing solutions

        
        output1 = Decs;
    
    case 'evaluate' % evaluate solution's fitness
        Data = varargin{1}; % load problem data
        Decs = varargin{2}; % load the current solution(s)
        
        % define the objective function in the following

        
        % define the inequal constraint(s) in the following, equal 
        % constraints should be transformed to inequal ones

       
        % calculate the constraint violation in the following


        % collect accessory data for understanding the solutions in the 
        % following (optional)

         output1 = ; % matrix for saving objective function values
         output2 = ; % matrix for saving constraint violation values (optional)
         output3 = ; % matrix or cells for saving accessory data (optional), a solution's accessory data should be saved in a row
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end