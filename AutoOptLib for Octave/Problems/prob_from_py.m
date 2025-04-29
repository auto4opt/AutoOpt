function [output1,output2,output3] = prob_from_py(varargin)
% Interface with Python problem file.
%
% if the modification of .py module can not be synchronized, try:
% >>
% clear classes
% obj = py.importlib.import_module('myfun');
% py.importlib.reload(obj);
%
% or restart the matlab software.

%----------------------------Copyright-------------------------------------
% Copyright (C) <2025>  <Qi Zhao>

% AutoOptLib is a free software. You can use, redistribute, and/or modify
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or any later 
% version. 
%--------------------------------------------------------------------------

switch varargin{end}
    case 'construct' % define problem properties
        Problem  = varargin{1};
        instance = varargin{2};
        Data     = struct('instanceInd',[]);
        
        % get problem type
        type = string(py.prob_interface.get_type());
        type = {type{1},type{2},type{3}};

        for i = 1:length(instance)
            Problem(i).type = type;
            
            % get solution space boundary
            [lower,upper] = double(py.prob_interface.get_bound());
            Problem(i).bound = [lower;upper];

            Data(i).instanceInd = i; % index of problem instance
        end
        output1 = Problem;
        output2 = Data;
        
    case 'repair' % repair solutions
        Decs = varargin{2};
        output1 = Decs;
        
    case 'evaluate' % evaluate solution's fitness
        Data = varargin{1};
        Decs = varargin{2};

        instanceInd = Data.instanceInd;
        [obj,con,acc] = double(py.prob_interface.evaluate(Decs,instanceInd)); % evaluate from Python
        output1 = obj;
        output2 = con; % onstraint violation values (optional)
        output3 = acc; % accessory data (optional)
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end