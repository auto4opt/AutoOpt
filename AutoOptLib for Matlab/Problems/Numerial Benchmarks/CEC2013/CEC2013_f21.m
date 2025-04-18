function [output1,output2,output3] = CEC2013_f21(varargin)

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
        lambda  = [1, 1e-6, 1e-26, 1e-6, 0.1];
        biasVal = [0, 100, 200, 300, 400];
        
        % 一般可设定最终常数偏移 f^* = 0 (或其它值)
        f_star  = 700;
        
        % 初始化用于存储每个子函数值 g_i(x) 的向量
        g_vals  = zeros(N,5);
        
        %% ============ 2) 逐个子函数计算 g_i(x) ============
        
        % -- 2.1) g1: Rotated Rosenbrock’s Function f6’
        %    公式: f6’(x) = sum_{i=1..(D-1)} [ 100 (z_i^2 - z_{i+1})^2 + (z_i - 1)^2 ] + f6*
        %    其中: z = M1( 2.048*(x - o)/100 ) + 1
        %           (假设 f6* 已内置成使得子函数在平移中心处为 0)
        i1 = 1;
        x_shifted = Decs - repmat(Os(i1,1:D), N, 1); 
        z         = (2.048/100)*( x_shifted * Ms(1:i1*D,1:D) ) + 1;  % (N×D)
        g_tmp = zeros(N,1);
        for r = 1:(D-1)
            g_tmp = g_tmp + 100*(z(:,r).^2 - z(:,r+1)).^2 + (z(:,r)-1).^2;
        end
        % g_tmp此时就是 f6’(x)
        g_vals(:,i1) = g_tmp;
        
        % -- 2.2) g2: Rotated Different Powers Function f5’
        %    公式: f5’(x) = sqrt( sum_{i=1..D} |z_i|^{ 2+4*(i-1)/(D-1) } ) + f5*
        %    “Rotated”版本若有需要则： z = M( x-o )。示例只做 M1:
        i2 = 2;
        x_shifted = Decs - repmat(Os(i2,1:D), N, 1);
        z         = x_shifted * Ms(i1*D+1:i2*D,1:D);
        
        exponents = linspace(2, 6, D);
        %   其中 6 = 2 + 4*(D-1)/(D-1).
        tmpSum = zeros(N,1);
        for dd = 1:D
            tmpSum = tmpSum + abs(z(:,dd)).^exponents(dd);
        end
        g_vals(:,i2) = sqrt(tmpSum);  % + f5* (假设内部已减去常数,同上)
        
        % -- 2.3) g3: Rotated Bent Cigar Function f3’
        %    公式: f3’(x) = z1^2 + 10^6 * sum_{i=2..D} z_i^2
        %           其中 z = M2( T^0.5_asy( M1( x-o ) ) )
        i3 = 3;
        x_shifted = Decs - repmat(Os(i3,1:D), N, 1);
        z1 = x_shifted * Ms(i2*D+1:i3*D,1:D);           % 先乘 M1
        z2 = computeTAsym(z1, 0.5);         % 例如 T^0.5_asym
        z  = z2 * Ms(i3*D+1:4*D,1:D);                  % 再乘 M2
    
        g_tmp = z(:,1).^2 + 1e6*sum(z(:,2:end).^2, 2);
        g_vals(:,i3) = g_tmp;
        
        % -- 2.4) g4: Rotated Discus Function f4’
        %    公式: f4’(x) = 10^6 * z1^2 + sum_{i=2..D} z_i^2
        %           z = T_0.5Z( M1( x - o ) )  (文献中略有不同写法)
        i4 = 4;
        x_shifted = Decs - repmat(Os(i4,1:D), N, 1);
        z1 = x_shifted * Ms(i3*D+1:i4*D,1:D);           
        z = computeTosz(z1);  % 假设有个 T_0.5Z 函数
        g_tmp = 1e6*z(:,1).^2 + sum(z(:,2:end).^2,2);
        g_vals(:,i4) = g_tmp;
        
        % -- 2.5) g5: Sphere Function f1’
        %    公式: f1’(x) = sum_{i=1..D} z_i^2
        %           z = x - o
        i5 = 5;
        x_shifted = Decs - repmat(Os(i5,1:D), N, 1);
        z = x_shifted;  % Sphere 无旋转
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