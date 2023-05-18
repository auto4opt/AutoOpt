function generateInstance
K    = 4;
Nt   = 4;
NR2  = 40:40:400;
b    = 2;
SNR  = 2;
PT   = 10.^(SNR/10);
Data = struct;
for i = 1:numel(NR2)
    [G,Hd,Hr,omega] = get_loc_pathloss_csi(K,Nt,NR2(i)); % produce the targeted problem, for b=2
    Data(i).b     = b;
    Data(i).PT    = PT;
    Data(i).G     = G;
    Data(i).Hd    = Hd;
    Data(i).Hr    = Hr;
    Data(i).omega = omega;
end
save('Beanforming.mat','Data');
end