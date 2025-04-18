function [output1,output2,output3] = CEC2013_f14(varargin)

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
        end      
        output1 = Problem;
        output2 = Data;
        
    case 'repair'
        Decs = varargin{2};
        output1 = Decs;
    
    case 'evaluate'
        Data = varargin{1};
        o    = Data.o;
        Decs = varargin{2};
        
        [N,D] = size(Decs);
        % z = Λ^{10}(10(x-o)) + 420.9687462275036
        Decs = Decs - repmat(o, N, 1);  % (x - o)
        Decs = 10 * Decs;               % 10(x-o)
        Decs = Decs*constructLambda(10, D);
        z = Decs + 4.209687462275036e+002; % 加上常数偏移
        
        % 计算g(z_i)
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
                g_val(cond2) = temp .* sin(sqrt(abs(temp))) - ((zi(cond2)-500).^2)/(10000*D);
            end
            
            % 条件3: z_i < -500
            if any(cond3)
                z_abs_mod = mod(abs_zi(cond3), 500); 
                temp = z_abs_mod - 500;
                g_val(cond3) = temp .* sin(sqrt(abs(temp))) - ((zi(cond3)+500).^2)/(10000*D);
            end
            
            g_sum = g_sum + g_val;
        end
        
        fit = 418.9829 * D - sum(g_sum);
        
        output1 = fit - 100;
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end