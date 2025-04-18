function [output1,output2,output3] = CEC2013_f13(varargin)

switch varargin{end}
    case 'construct'
        type     = {'continuous','static','certain'};
        Problem  = varargin{1};
        instance = varargin{2};
        orgData  = load('shift_data.mat');
        o        = orgData.data(1,:);
        Data     = struct('o',[],'M',[]);
        for i = 1:length(instance)
             Problem(i).type = type;

            D = instance(i);
            lower = zeros(1,D)-100;
            upper = zeros(1,D)+100;
            Problem(i).bound = [lower;upper];
            
            if length(o) >= D
                curr_o = o(1:D);
            else 
                curr_o = -100+200*rand(1,D);
            end
            Data(i).o = curr_o;
            
            fileMap = containers.Map([2, 5, 10, 20, 30, 40, 50, 70, 80, 90, 100], ...
                         {'M_D2.mat', 'M_D5.mat', 'M_D10.mat', 'M_D20.mat', 'M_D30.mat', 'M_D40.mat', 'M_D50.mat', 'M_D70.mat', 'M_D80.mat', 'M_D90.mat', 'M_D100.mat'});
            if isKey(fileMap, D)
                M = load(fileMap(D));
                M1 = M.data(1:D,1:D);
                M2 = M.data(D+1:2*D,1:D);
            else
                A = normrnd(0, 1, D, D);
                [M1, ~] = cGram_Schmidt(A);
                [M2, ~] = cGram_Schmidt(A);
            end
            Data(i).M1 = M1;
            Data(i).M2 = M2;
        end      
        output1 = Problem;
        output2 = Data;
        
   case 'repair'
        Decs = varargin{2};
        output1 = Decs;
    
    case 'evaluate'
        Data = varargin{1};
        o    = Data.o;
        M1   = Data.M1;
        M2   = Data.M2;
        Decs = varargin{2};
        
        [N,D] = size(Decs);
        % 1. 移位并缩放到[-5.12,5.12]
        Decs = Decs - repmat(o, N, 1);     % x - o
        Decs = (5.12/100)*Decs;            % 5.12*(x-o)/100
        Decs = Decs * M1;                  % 得到 x̄
        
        % 2. 非连续处理 y_i
        % 若 |x̄_i| > 0.5 则 y_i = round(2*x̄_i)/2，否则 y_i = x̄_i
        for i = 1:D
            mask = abs(Decs(:,i)) > 0.5;
            Decs(mask,i) = round(2*Decs(mask,i))/2; 
        end
        % 此时 Decs即为 y
        
        % 3. T_osz 变换
        Decs = computeTosz(Decs);
        
        % 4. T_asy^{0.2} 变换
        Decs = computeTAsym(Decs, 0.2);
        
        % 5. 旋转：M2、Lambda^{10}、M1
        Decs = Decs * M2;
        Decs = Decs * constructLambda(10, D);
        z    = Decs * M1;
        
        % 6. 计算 Rastrigin 基本形式的值并加上偏移
        fit = 0;
        for i = 1:D
            fit = fit + (z(:,i).^2 - 10.*cos(2.*pi.*z(:,i)) + 10);
        end
        
        % 根据题中设定对 fit 进行整体偏移
        % 这里示例使用 fit - 200，实际请根据题中给定的 f_{13}^* 调整
        output1 = fit - 200;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end