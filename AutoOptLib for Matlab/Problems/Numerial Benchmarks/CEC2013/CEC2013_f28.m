function [output1,output2,output3] = CEC2013_f28(varargin)

switch varargin{end}
    case 'construct'
    type     = {'continuous','static','certain'};
    Problem  = varargin{1};
    instance = varargin{2};
    orgData  = load('shift_data.mat');
    o       = orgData.data;
    Data     = struct('Os',[],'M',[]);
    for i = 1:length(instance)
        Problem(i).type = type;

        D = instance(i); % 维度
        lower = zeros(1,D)-100;
        upper = zeros(1,D)+100;
        Problem(i).bound = [lower;upper];
        
        Data(i).o = o;
        
        fileMap = containers.Map([2, 5, 10, 20, 30, 40, 50, 70, 80, 90, 100], ...
                     {'M_D2.mat', 'M_D5.mat', 'M_D10.mat', 'M_D20.mat', 'M_D30.mat', 'M_D40.mat', 'M_D50.mat', 'M_D70.mat', 'M_D80.mat', 'M_D90.mat', 'M_D100.mat'});
        if isKey(fileMap, D)
            M = load(fileMap(D)).data;
        else
            error("Unsuportted dimensions!")
        end
        Data(i).M = M;
    end      
    output1 = Problem;
    output2 = Data;
        
    case 'repair'
        Decs = varargin{2};
        % 这里可以对 Decs 进行边界裁剪或其它修复
        % 也可以什么都不做, 直接返回
        output1 = Decs;
        output2 = [];
        output3 = [];
    
    case 'evaluate'
        %% ============ 0) 读取问题数据与决策变量 ============
        Data = varargin{1};       % 构造阶段保存的数据
        Decs = varargin{2};       % 当前要评估的一组解 (N×D)
        
        % 你可能需要存储多个子函数的平移/旋转信息（o_i, M1_i, M2_i等）。
        % 或者有时只存储一个大矩阵，再在内部拆分给各子函数使用。
        % 这里仅示例，把它们用 cell{} 分开存储：
        Os  = Data.o;
        Ms = Data.M;
        
        [N, D] = size(Decs);
        
        %% ============ 1) 定义组合参数 ============
        sigma   = [10, 20, 30, 40, 50];
        lambda  = [2.5, 2.5e-3, 2.5, 5e-4, 0.1];
        biasVal = [0, 100, 200, 300, 400];
        
        % 一般可设定最终常数偏移 f^* = 0 (或其它值)
        f_star  = 1400;
        
        % 初始化用于存储每个子函数值 g_i(x) 的向量
        g_vals  = zeros(N,5);
        
        %% ============ 2) 逐个子函数计算 g_i(x) ============
        i1 = 1;
        x_shifted = Decs - repmat(Os(i1,1:D), N, 1); 
        z         = (5/100)*x_shifted*Ms(1:D,1:D)+1;
        g_sum = zeros(N,1);
        for i=1:D
            zi = z(:,i);
            j = i + 1;
            if j>D
                j=1;
            end
            zj = z(:,j);
            s = 100*(zi.^2 - zj).^2 + (zi-1).^2;
            g_sum = g_sum + (s.^2/4000 - cos(s)/sqrt(i) + 1);
        end

        g_vals(:,i1) = g_sum;
        
        
        i2 = 2;
        x_shifted = Decs - repmat(Os(i2,1:D), N, 1);
        x         = x_shifted*Ms(1:D,1:D);
        x         = computeTAsym(x,0.5);
        x         = x*Ms(D+1:i2*D,1:D)*constructLambda(10,D);

        g_sum = zeros(N,1);
        for i = 1:D-1
            y_i = x(:,i);
            y_j = x(:, i+1);
            z = sqrt(y_i.^2 + y_j.^2);
            g_sum = g_sum+(sqrt(z)+sqrt(z)*sin(50*(z.^0.2)).^2).^2;
        end
        g_sum = g_sum/(D-1);
        
        g_vals(:,i2) = g_sum;


        i3 = 3;
        x_shifted = Decs - repmat(Os(i3,1:D), N, 1); 
        z         = 10*( (x_shifted * Ms(2*D+1:i3*D,1:D)) * constructLambda(10,D) ) + 4.209687462275036e+002;  % (N×D)
        g_sum = zeros(N,1);
        for i = 1:D
            zi = z(:,i);
            abs_zi = abs(zi);
            
            % 条件分支计算g(zi)
            cond1 = (abs_zi <= 500);
            cond2 = (zi > 500);
            cond3 = (zi < -500);
            
            g_val = zeros(N,1);
            
            % 条件1: |z_i| <= 500
            g_val(cond1) = zi(cond1).*sin(abs_zi(cond1).^(1/2));
            
            % 条件2: z_i > 500
            if any(cond2)
                z_mod = mod(zi(cond2), 500);  % mod(z_i,500)
                temp = 500 - z_mod;
                g_val(cond2) = temp .* sin(sqrt(abs(temp))) + ((zi(cond2)-500).^2)/(10000*D);
            end
            
            % 条件3: z_i < -500
            if any(cond3)
                z_abs_mod = mod(abs_zi(cond3), 500); 
                temp = z_abs_mod - 500;
                g_val(cond3) = temp .* sin(sqrt(abs(temp))) + ((zi(cond3)+500).^2)/(10000*D);
            end
            
            g_sum = g_sum + g_val;
        end

        g_vals(:,i3) = 418.9829 * D - sum(g_sum);



        i4 = 4;
        x_shifted = Decs - repmat(Os(i4,1:D), N, 1);
        z         = x_shifted*Ms(i3*D+1:i4*D,1:D);
        z         = computeTAsym(z,0.5)*Ms(i4*D+1:5*D, 1:D);
        g_sum     = zeros(N,1);

        for i=1:D
            zi = z(:,i);
            j = i + 1;
            if j>D
                j=1;
            end
            zj = z(:,j);
            s = 0.5 + (sin(sqrt(zi.^2+zj.^2)).^2-0.5) / (1+0.001*(zi.^2+zj.^2)).^2;
            g_sum = g_sum + s;
        end
        
        g_vals(:,i4) = g_sum;

        i5 = 5;
        x_shifted = Decs - repmat(Os(i3,1:D), N, 1); 
        z         = x_shifted;

        g_vals(:,i5) = sum(z.^2,2);
        
        %% ============ 3) 计算组合权重 w_i ============
        %   w_i = 1/sqrt( sum( (x - o_i)^2 ) ) * exp( - sum( (x - o_i)^2 ) / (2 * D * sigma_i^2 ) )
        %   然后归一化: omega_i = w_i / sum(w_k)
        
        w = zeros(N,5);
        for i = 1:5
            diff_i = Decs - repmat(Os(i,1:D), N, 1);
            dist2  = sum(diff_i.^2, 2);  % (N×1)
            
            w_i = 1./sqrt(dist2+1e-30) .* ...
                  exp(-dist2/(2*D*(sigma(i)^2)+1e-30));
            
            w(:,i) = w_i;
        end
        w_sum = sum(w,2) + 1e-30;  % 防止除0
        omega = w ./ w_sum;        % (N×5)
        
        %% ============ 4) 计算组合函数 f(x) ============
        % f(x) = sum_{i=1..5} [ omega_i * (lambda_i*g_i(x) + bias_i ) ] + f_star
        compVal = zeros(N,1);
        for i = 1:5
            compVal = compVal + ...
                omega(:,i) .* ( lambda(i)*g_vals(:,i) + biasVal(i) );
        end
        
        fit = compVal + f_star;
        
        %% ============ 5) 返回结果 ============
        output1 = fit;    % 最终适应值
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end