function [output1,output2,output3] = CEC2013_f26(varargin)

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

        D = instance(i);
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
        output1 = Decs;
        output2 = [];
        output3 = [];
    
    case 'evaluate'
        Data = varargin{1};       
        Decs = varargin{2};      
        Os  = Data.o;
        Ms = Data.M;
        
        [N, D] = size(Decs);
        
        sigma   = [10, 10, 10, 10, 10];
        lambda  = [0.25, 1, 1e-7, 2.5, 10];
        biasVal = [0, 100, 200, 300, 400];
        
        f_star  = 1200;
        
        g_vals  = zeros(N,5);
        
        i1 = 1;
        x_shifted = Decs - repmat(Os(i1,1:D), N, 1); 
        z         = 10*( (x_shifted * Ms(1:D,1:D)) * constructLambda(10,D) ) + 4.209687462275036e+002;  % (N×D)
        g_sum = zeros(N,1);
        for i = 1:D
            zi = z(:,i);
            abs_zi = abs(zi);

            cond1 = (abs_zi <= 500);
            cond2 = (zi > 500);
            cond3 = (zi < -500);
            
            g_val = zeros(N,1);
            
            g_val(cond1) = zi(cond1).*sin(abs_zi(cond1).^(1/2));
            
            if any(cond2)
                z_mod = mod(zi(cond2), 500); 
                temp = 500 - z_mod;
                g_val(cond2) = temp .* sin(sqrt(abs(temp))) + ((zi(cond2)-500).^2)/(10000*D);
            end
            
            if any(cond3)
                z_abs_mod = mod(abs_zi(cond3), 500); 
                temp = z_abs_mod - 500;
                g_val(cond3) = temp .* sin(sqrt(abs(temp))) + ((zi(cond3)+500).^2)/(10000*D);
            end
            
            g_sum = g_sum + g_val;
        end

        g_vals(:,i1) = 418.9829 * D - sum(g_sum);
        
        i2 = 2;
        x_shifted = Decs - repmat(Os(i2,1:D), N, 1); 
        z         = 10* ((5.12/100)*x_shifted * Ms(i1*D+1:i2*D,1:D)); 
        z         = computeTosz(z);
        z         = computeTAsym(z,0.2);
        z         = z*Ms(i2*D+1:3*D,1:D)*constructLambda(10,D)*Ms(i1*D+1:i2*D,1:D);

        g_sum = zeros(N,1);
        for i = 1:D
           g_sum = g_sum + (z(:,i).^2 - 10.*cos(2.*pi.*z(:,i)) + 10);
        end
        
        g_vals(:,i2) = g_sum;

        i3 = 3;
        x_shifted = Decs - repmat(Os(i3,1:D), N, 1); 
        z         = x_shifted * Ms(i2*D+1:i3*D,1:D); 
        z         = computeTosz(z);

        g_sum = zeros(N,1);
        for i = 1:D
           g_sum = g_sum + ((10^6)^((i-1)/(D-1)) * z(:,i).^2);
        end
        
        g_vals(:,i3) = g_sum;

        i4 = 4;
        x_shifted = Decs - repmat(Os(i4,1:D), N, 1);
        z         = (0.5/100)*( (x_shifted * Ms(i3*D+1:i4*D,1:D))); 
        z         = computeTAsym(z,0.5);
        z         = z*Ms(i4*D+1:5*D,1:D)*constructLambda(10,D);
        g_sum = zeros(N,1);
        a = 0.5;
        b = 3;
        kmax = 20;

        term2 = 0;
        for i = 1:D
            inner_sum = 0;
            for k = 0:kmax
                zi = z(:,i);
                inner_sum = inner_sum + a^k * cos(2 * pi * b^k * (zi + 0.5));
            end
            g_sum = g_sum + inner_sum;
        end
        for k = 0:kmax
            term2 = term2 + a^k * cos(2 * pi * b^k * 0.5);
        end
        g_sum = g_sum - D*term2;
        
        g_vals(:,i4) = g_sum;

        i5 = 5;
        x_shifted = Decs - repmat(Os(i5,1:D), N, 1); 
        z         = 6*x_shifted*Ms(i4*D+1:i5*D,1:D)*constructLambda(100,D);

        sum_term = sum(z.^2, 2) / 4000;

        cos_terms = ones(N,1);
        for i = 1:D
            cos_terms = cos_terms .* cos(z(:,i)./sqrt(i));
        end

        g_vals(:,i5) = sum_term - cos_terms + 1;
        
        w = zeros(N,5);
        for i = 1:5
            diff_i = Decs - repmat(Os(i,1:D), N, 1);
            dist2  = sum(diff_i.^2, 2);  
            
            w_i = 1./sqrt(dist2+1e-30) .* ...
                  exp(-dist2/(2*D*(sigma(i)^2)+1e-30));
            
            w(:,i) = w_i;
        end
        w_sum = sum(w,2) + 1e-30;  
        omega = w ./ w_sum;        

        compVal = zeros(N,1);
        for i = 1:5
            compVal = compVal + ...
                omega(:,i) .* ( lambda(i)*g_vals(:,i) + biasVal(i) );
        end
        
        fit = compVal + f_star;
        
        output1 = fit;   
end

if ~exist('output2','var')
    output2 = [];
end
if ~exist('output3','var')
    output3 = [];
end
end