clc;
clear all;
close all;
tic;
%收发机参数
c=3e8;              %光速
fs=75e9;            %发射信号起始频率，单位：Hz
B=35e9;             %带宽
fc=fs+B/2;          %发射信号中心频率，单位：Hz
lambda_c=c/fc;      %中心频率波长
Nx=255;             %收发机在高度维采样点数
Nf=201;                          %频率采样点数
delta_f=B/(Nf-1);                %频率采样间隔
f=fs+(0:Nf-1)*delta_f;           %频率序列

theta_1ant1=60;                       %天线波束角，单位：°
theta_ant=theta_1ant1*pi/180;            %转化为弧度制
Kxmax=(4*pi*(fs+B/2)/c)*sin(theta_ant/2);
deltaX=pi/Kxmax;                         % x方向上的采样间隔
Lx=(Nx-1)*deltaX;
R0=Lx/2/tan(theta_ant/2);                %天线阵列到目标区域中心的距离，单位：m 
x_tr = ((-(Nx-1)/2:(Nx-1)/2)*deltaX).';   % 收发机高度维的采样坐标
y_tr=-R0.*ones(Nx,1);                      %收发机在距离维采样坐标
             

%构建回波信号 
j=sqrt(-1);
 object=[0,0,1
     0.05,0.1,1
     0.1,0.2,1
  
    ];      %点目标

S=zeros(Nx,Nf); 
num=size(object,1);
for  i1=1:Nx
     s=zeros(1,Nf);
         for i2=1:num
             x=object(i2,1);                                 %目标的横坐标
             y=object(i2,2);                                 %目标的纵坐标
             A=object(i2,3);                                 %目标的幅度
             R=sqrt((x_tr(i1)-x).^2+(y_tr(i1)-y).^2);        %天线到目标的距离
             s=s+A*exp(-j*2*pi*f*2*R/c);                    
         end
         S(i1,:)=s;                                           %回波
end


 
kx_max=pi/deltaX;
kx=linspace(-kx_max,kx_max,Nx);
k=2*pi*f/c;
k1=2*k;
Kx=kx'*ones(1,Nf); 
 
for f_index=1:Nf
    Ky1=sqrt((k1(f_index)).^2-(Kx(:,1)).^2);
    Ky(:,f_index)=Ky1;
end
figure,plot(Kx,Ky,'.');
xlabel('Kx');ylabel('Ky');
title('波数域');
grid on;
S_FT=[];
S_FT=fftshift(fft(fftshift(S,1),[],1),1); 


%% 相位补偿
W=exp(j.*Ky.*R0);                    %补偿项
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



Dy=(c/2/B)*(Nf-1);
ObjectX_pos=((-(Nx-1)/2:(Nx-1)/2)*deltaX);  % 收发机方位维X维的采样坐标
ObjectY_pos=linspace(-Dy/2,Dy/2,Nf2);  % 收发机方位维Y维的采样坐标


figure,
imagesc(ObjectX_pos',ObjectY_pos,abs(S_iftxyz'));
xlabel('高度维/m'),ylabel('距离维/m');
set(gca, 'YDir', 'normal');

toc;