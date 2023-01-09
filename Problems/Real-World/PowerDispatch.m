function [output1,output2,output3] = PowerDispatch(varargin)
switch varargin{end}
    case 'construct' % construct problem and data
        Problem = varargin{1};
        instance = varargin{2};

        orgData = readmatrix('DataPowerDispatch.xlsx','Sheet',1);
        Data = struct('orgData',[]);
        
        for i = 1:numel(instance)
            Problem(i) = Problem(1); 
            Data(i).orgData = orgData;
            Problem(i).bound = [orgData(:,4)';orgData(:,5)']; % bound constraint
        end

        output1 = Problem(instance);
        output2 = Data(instance);

    case 'evaluate' % evaluate fitness
        Data = varargin{1}.orgData;
        Decs = varargin{2}; % solution set

        N = size(Decs,1); % population size
        Objs = zeros(N,1); % fitness of solutions
        for i = 1:N
            solution = Decs(i,:); % 1*|generator|
            F = zeros(1,size(Data,1)); % fuel cost of generators, 1*|generator|
            for j = 1:size(Data,1) % for each generator
                F(j) = Data(j,1) + Data(j,2)*solution(j) + Data(j,3)*solution(j)^2; % objective function
            end 
            Objs(i) = sum(F);
            if sum(solution) < Data(1,6) % power balance constraint 
               Objs(i) = 10^6; % assign an extreme large fitness value to infeasible solution  
            end
        end
        output1 = Objs;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end