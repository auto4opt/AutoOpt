function sumRs = get_sum_rate(ms,b,Hd,Hr,G,PT,omega)
[~, Nt]=size(G);
m_num=size(ms,1);
sumRs=zeros(1, m_num);
Qs = exp(1i*ms/(2^b)*2*pi);

for mi=1:m_num
    %% get hybrid CSI: F        
    F= Hd+Hr*diag(Qs(mi,:)')*G;
    
    %% get BS beamforming: VD
    W=F'/(F*F');
    vk=diag(W'*W);
    [vk_fill] = water_filling(PT, vk);
    PKs=vk_fill./vk;
    VD=W.*repmat(sqrt(PKs'),Nt,1);
    
    %% get sum rate
    He = abs(F * VD).^2;
    FV_K_others = sum(He,2) - diag(He);
    R = diag(He)./( FV_K_others + 1);
    sumRs(mi) = 1./sum(log(R+1).*omega');
end
end