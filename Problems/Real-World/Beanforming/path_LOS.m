function [ loss ] = path_LOS( d )
%     d=d/1000;
%     loss=89.5 + 16.9*log10(d); 
%     loss=38.46 + 20*log10(d);
%     loss=35.6 + 22*log10(d); %默认。数值越大，代表path loss越大，衰减越严重，sum rate越小。
% loss=20 + 20*log10(d); %很难收敛，对于WSR
loss=20 + 20*log10(d);%20收敛很慢
end 

