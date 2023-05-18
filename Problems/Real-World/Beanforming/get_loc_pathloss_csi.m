%% 
function  [G, Hd, Hr, omega] = get_loc_pathloss_csi(K, Nt, NR2)
%% generate new location for the K users
Pt=zeros(K,2);
%%
Lroom=100; %默认200
Wroom=30; %默认30
k1=[1,0];
k2=[0,1];
R=10;
%%
for k0=1:K
    r=rand(1,1)*R;
    theta=rand(1,1)*2*pi;
    px=r*cos(theta);
    py=r*sin(theta);
    pt=[Lroom,Wroom]+px*k1+py*k2;
    Pt(k0,:)=pt;
end
%%
%save('user_location.mat','K','Pt','Lroom','Wroom','R');

plot_ind=0;
if plot_ind==1
    figure
    plot(Pt(:,1),Pt(:,2),'ro');
    xlim([Lroom-R,Lroom+R]);ylim([Wroom-R,Wroom+R]);
    hold on
    theta=linspace(0,1,100).*2.*pi; %半径100
    hold on
    plot(Lroom+R*cos(theta),Wroom+R*sin(theta),'r.')
end

%-------------------------------------------------------------------------
%% compute the pathloss according to the users' locations
%% 
AP=[0,0];
%%
IRS=[100,0];
d_g=sqrt(sum(abs(AP-IRS).^2));
L_g=path_LOS( d_g );
%% IRS-assist link
Lu=zeros(1,K);
for k0=1:K
    pt=Pt(k0,:)-IRS;
    du=sqrt(sum(abs(pt).^2));
    Lu(k0)=path_LOS( du );
end
Lu=Lu+L_g;
%% direct link
Ld=zeros(1,K);
for k0=1:K
    pt=Pt(k0,:)-AP;
    dk=sqrt(sum(abs(pt).^2));
    Ld(k0)=path_NLOS( dk );
%     Ld_test(k0)=path_LOS( dk );
end
%Lu-Ld
%%
%save('user_pathloss.mat','K','Lu','Ld');

%-------------------------------------------------------------------------
%% Generate the channel coefficients
noise=-170+10*log10(180*1e3); %传输带宽越大，sum rate越低
path_d=10.^((-noise-Ld)/10);
path_i=10.^((-noise-Lu)/10);
%%
ite=1; 
%%
%%
pd=sqrt(path_d);
pd=repmat(pd.',1,Nt);
ps=sqrt(path_i);
ps=repmat(ps.',1,NR2);
%% theta_init, channel Hd
Hd_w=zeros(K,Nt,ite);

for j0=1
    Hd=sqrt(1/2).*(randn(K,Nt)+1j.*randn(K,Nt));
    %%
    Hd_w(:,:,j0)=Hd;
end
% load('CSI8(20,20)Hd_NR2_10.mat','K','NR2','Nt','Pt',...
%   'pd','ps','Hd','AP_angle','IRS_angle',...
% 'User_angle','Hd');

%%
eb=10;
eb2=1/(1+eb);
eb1=1-eb2;
eb1=sqrt(eb1);
eb2=sqrt(eb2);
%% channel G
AP_angle=rand(1,1);
IRS_angle=rand(1,1);
G_sig=zeros(NR2,Nt,ite);
for i0=1:ite
    G_sig(:,:,i0)=sqrt(1/2).*(randn(NR2,Nt)+1j.*randn(NR2,Nt));
end
%% channel Hr_w
User_angle=rand(1,K);
Hr_sig=zeros(K,NR2,ite);
for i0=1:ite
    Hr_sig(:,:,i0)=sqrt(1/2).*(randn(K,NR2)+1j.*randn(K,NR2));
end

Hd=pd.*Hd_w(:,:,1);
G=channel_G(AP_angle,IRS_angle,G_sig(:,:,1),eb1,eb2,NR2,Nt);
Hr=ps.*channel_Hr(User_angle,Hr_sig(:,:,1),eb1,eb2,K,NR2);

% Generate weights of users ——omega
% weight=1./((path_d));
% omega=weight./sum(weight);
omega=ones(1, K);
%%
% save('CSI8(20,20)Hd_NR2_10.mat','K','NR2','Nt','Pt',...
%   'pd','ps','Hd','AP_angle','IRS_angle',...
% 'User_angle','Hd','Hr','G');
end