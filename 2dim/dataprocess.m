function [S_iftxyz,Nf2,S_XFT,k1,Ky,kx]=dataprocess(S,deltaX,f,Nx,Nf,R0)
c=3e8;
kx_max=pi/deltaX;
kx=linspace(-kx_max,kx_max,Nx);
k=2*pi*f/c;
k1=2*k;
Kx=kx'*ones(1,Nf);

for f_index=1:Nf
    Ky1=sqrt((k1(f_index)).^2-(Kx(:,1)).^2);
    Ky(:,f_index)=Ky1;
end
% figure,plot(Kx,Ky,'.');
% xlabel('Kx');ylabel('Ky');
% title('波数域');
% grid on;
S_FT=[];
S_FT=fftshift(fft(fftshift(S,1),[],1),1);    %傅里叶变换


%% 相位补偿
W=exp(1j.*Ky.*R0);                    %补偿项
S_XFT=S_FT.*W;                       %相位补偿


%% 距离维插值
[ax,ay]=find(Ky==max(max(Ky)));
ky_mid=Ky(ax(1),1);
ky_min=min(min(Ky));
ky_max=max(max(Ky));

Nf2=fix((ky_max-ky_min)/(ky_max-ky_mid)*Nf)+1;
ky1=linspace(ky_min,ky_max,Nf2);

ky2=repmat(ky1,Nx,1);                %% 均匀的ky
kx2=repmat(kx',1,Nf2);
k2=sqrt(kx2.^2+ky2.^2);              %% 非均匀的k
dk=(max(k1)-min(k1))/(Nf-1);         %% k的间隔
S1=zeros(Nx,Nf2);                    %% 插值后的回波矩阵

% 归一化
k_un=((k1-min(k1))/dk)+1;            %% 之前的k
k_un=repmat(k_un,Nx,1);


k2_un=(k2-min(k1))/dk+1;             %% 插值后的k


for i3=1:Nx
    [x1,z1]=find(k2_un(i3,:)>k_un(1,Nf));
    if isempty(z1)==1
        z1=Nf2+1;
    end
    zz2=min(min(z1));

    [x3,z3]=find(k2_un(i3,:)<k_un(1,1));
    if isempty(z3)==1
        z3=0;
    end
    zz3=max(max(z3));

    x11=k_un(i3,:);
    y11=S_XFT(i3,:);
    xs=k2_un(i3,zz3+1:zz2-1);
    ys=interp1(x11,y11,xs,'spline');
    S1(i3,zz3+1:zz2-1)=ys;
end



S_iftxyz=zeros(Nx,Nf2);
S_iftxyz=fftshift(ifft2(fftshift((S1))));
end