function [ loss ] = path_LOS( d )
%     d=d/1000;
%     loss=89.5 + 16.9*log10(d); 
%     loss=38.46 + 20*log10(d);
%     loss=35.6 + 22*log10(d); %Ĭ�ϡ���ֵԽ�󣬴���path lossԽ��˥��Խ���أ�sum rateԽС��
% loss=20 + 20*log10(d); %��������������WSR
loss=20 + 20*log10(d);%20��������
end 

