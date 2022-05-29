%% 不使用fft 进行3d rma 成像，未补零  通过变换角度进行调整波数域通过变换角度进行调整波数域 @ty
%% 二维仿真回波成像
clc;
clear;
close all;
tic;
%% 场景参数设置
% 收发机参数
c=3e8;                           % 波的传播速度，单位：m/s
f0=75e9;                         % 雷达发射信号中心频率，单位：Hz
B=35e9;                         % 带宽，单位：Hz
fs=110e9;
% R0=2;
% 探测器位置参数
lambda=c/(f0+0.5*B);
% %% 111
% Nx=128;
% Ny = 128;                   % 收发机在Y维方位维上的采样点个数
% Nf=64;

% 222
Nx=256;
Ny = 256;                   % 收发机在Y维方位维上的采样点个数
Nf=128;

% %% 333
% Nx=512;
% Ny = 512;                   % 收发机在Y维方位维上的采样点个数
% Nf=256;

%% 水平方向
th=20*pi/180;
Kxmax=4*pi*fs/c*sin(th/2);
x_step=pi/Kxmax;          % x方向上的采样间隔
% if x_step>(lambda/2)     
%     x_step=lambda/2;
% %     Kxmax=pi/x_step;
% end
deltaX=x_step;                  % 收发机在X维方位维上的采样间距  默认半波长，单位：m
lambda2=lambda/x_step;

Lx=Nx*x_step;
R0=Lx/2/tan(th/2);         % 方位维X维的孔径长度


% Lx=R0*tan(th/2)*2;         % 方位维X维的孔径长度
% Lx2=R0*tan(th/2);  

% Nx=2*fix(Lx2/x_step)+1;                   % 天线阵数
TRX_pos = ((-(Nx-1)/2:(Nx-1)/2)*deltaX).';  % 收发机方位维X维的采样坐标

%% 竖直方向

% Ny=2^nextpow2(Nx);

% deltaY=lambda/2;                 % 收发机在Y维方位维上的采样间距，单位：m
deltaY=deltaX;
TRY_pos = ((-(Ny-1)/2:(Ny-1)/2)*deltaY).';  % 收发机方位维Y维的采样坐标
Ly=(Ny-1)*deltaY;                              % 方位维Y维的孔径长度

% 扫频参数

delta_f = B/(Nf-1);               % 收发机频率采样间隔，单位：Hz
rhoz = c/2/B;                     % 距离维理论分辨率
TRZ =-R0;                         % 收发机距离维位置，单位：m
%% 理论分辨率计算
rhoy=deltaX;
rhox=deltaY;
%% 波数域计算
fvec=zeros(Nf,1);
kvec=zeros(Nf,1);
fvec = (f0:delta_f:f0+B).';
kvec = 2*pi*fvec/c;                             %波数域K的序列                                 

Kz=zeros(Nx,Ny,Nf);
Kx=zeros(Nx,Ny);
Ky=zeros(Nx,Ny);

kx_max=Kxmax;   % 求kx的波域范围
kx=(linspace(-kx_max,kx_max,Nx)).';   % 设置kx的波域坐标
ky_max=pi/deltaY;   % 求ky的波域范围
ky = (linspace(-ky_max,ky_max,Ny)).';   % 设置ky的波域坐标
Kx = kx * ones(1,Ny);  
Ky = ones(Nx,1) * (ky.');

for index_f=1:Nf                                %用Kx,Ky和kvec得到Kz
    k = kvec(index_f);
    kz = sqrt((2*k).^2 - Kx.^2 - Ky.^2); 
    Kz(:,:,index_f) = kz;
end
k1=2*kvec;
K=zeros(Nx,Ny,Nf);
for i=1:Nx
    for j=1:Ny
        K(i,j,:)=k1.';
    end
end  
kx2=zeros(Nx,Nf);
kz2=zeros(Nx,Nf);
for i=1:Nf
    kx2(:,i)=kx;
end
kkk=permute(Kz,[1 3 2]);
kz2=(kkk(:,:,1));
%% 投影c.v
num=numel(Kz);
z_map1=reshape(Kz,1,num) ;
z_map2=unique(z_map1);
tic;
z_map=sort(z_map2);
toc;
% Mzmap=std(z_map,1);
Mzmap=std(z_map);
uzmap=abs(mean(z_map));
cvkz=Mzmap./uzmap;
%% 读取实验数据并构建回波信号
j = sqrt(-1);
%% 目标设置
% Ptar=[0,0,0,1]*diag([rhox rhoy rhoz 1]);     % [x轴坐标，y轴坐标，z轴坐标，反射能量]
%   Ptar=[15,35,15,1;16,35,15,1;17,35,15,1;18,35,15,1;19,35,15,1;20,35,15,1;21,35,15,1;22,35,15,1;23,35,15,1;24,35,15,1;25,35,15,1;26,35,15,1;
%       27,35,15,1;28,35,15,1;29,35,15,1;30,35,15,1;31,35,15,1;32,35,15,1;33,35,15,1;34,35,15,1;35,35,15,1;36,35,15,1;37,35,15,1;38,35,15,1;39,35,15,1;40,35,15,1;
%       41,35,15,1;42,35,15,1;43,35,15,1;44,35,15,1;45,35,15,1;46,35,15,1;47,35,15,1;48,35,15,1;49,35,15,1;50,35,15,1;
%       30,15,15,1;30,16,15,1;30,17,15,1;30,17,15,1;30,18,15,1;30,19,15,1;30,20,15,1;30,21,15,1;30,22,15,1;30,23,15,1;30,24,15,1;30,25,15,1;30,25,15,1;30,26,15,1;30,27,15,1;30,28,15,1;30,29,15,1;30,30,15,1;
%       30,31,15,1;30,32,15,1;30,33,15,1;30,34,15,1;30,35,15,1;30,15,15,1;30,14,15,1;30,13,15,1;30,12,15,1;30,11,15,1;30,10,15,1;30,9,15,1;30,8,15,1;30,7,15,1;30,6,15,1;30,5,15,1;30,4,15,1;30,3,15,1;30,2,15,1;30,1,15,1;30,0,15,1;
%       ]*diag([rhox/100*Nx/2 rhoy/75*Ny/2 rhoz/30*Nf/2 1]);     % [x轴坐标，y轴坐标，z轴坐标，反射能量]
  
    Ptar=[0,0,0.5,1;0.05,0,0.5,1;0.1,0,0.5,1;0.15,0,0.5,1;0.2,0,0.5,1;0.25,0,0.5,1;0.3,0,0.5,1;0.35,0,0.5,1;0.4,0,0.5,1;
            -0.05,0,0.5,1;-0.1,0,0.5,1;-0.15,0,0.5,1;-0.2,0,0.5,1;-0.25,0,0.5,1;-0.3,0,0.5,1;-0.35,0,0.5,1;-0.4,0,0.5,1;
          0.4,-0.4,0.5,1;0.4,-0.35,0.5,1;0.4,-0.3,0.5,1;0.4,-0.25,0.5,1;0.4,-0.2,0.5,1;0.4,-0.15,0.5,1;0.4,-0.1,0.5,1;0.4,-0.05,0.5,1;0.4,0,0.5,1;
          0.4,0.4,0.5,1;0.4,0.35,0.5,1;0.4,0.3,0.5,1;0.4,0.25,0.5,1;0.4,0.2,0.5,1;0.4,0.15,0.5,1;0.4,0.1,0.5,1;0.4,0.05,0.5,1;
          0,0,0.45,1;0,0,0.40,1;0,0,0.35,1;0,0,0.30,1;0,0,0.25,1;0,0,0.20,1;0,0,0.15,1;0,0,0.10,1;0,0,0.05,1;0,0,0,1;
          ]*diag([rhox*Nx/2 rhoy*Ny/2 rhoz*Nf/2 1]);     % [x轴坐标，y轴坐标，z轴坐标，反射能量]
Object_num = length(Ptar(:,1));                                  % 计算目标点个数（像素点数）

%% 仿真回波信号
j = sqrt(-1);
S = zeros(Nx,Ny,Nf);                   %构建每个探测点的空回波

for index_y = 1:Ny
    for index_x = 1:Nx
        s = 0;
        for i=1:1:Object_num
            xn = Ptar(i,1);
            yn = Ptar(i,2);
            zn = Ptar(i,3);
            R = sqrt((xn - TRX_pos(index_x)).^2 + (yn - TRY_pos(index_y)).^2 + (zn-TRZ).^2);  %计算每个收发器到每个目标点的距离
            sigma = Ptar(i,4);                                                     %假设所有的点产生的幅值一样
%             s = s + sigma*exp(-1j*2*kvec*R);   %加入时延产生回波，并将所有目标点的回波叠加在一起，构建收发器接收信号的真实情况
%            s = s + sigma*exp(-1j*2*kvec*R).*exp(j*2*kvec*R0);
                s = s + sigma*exp(-1j*2*kvec*(R-R0));
        end
        S(index_x,index_y,:) = s;
    end
end

S = cell2mat(struct2cell(load('./echo_data/S2.mat')));
% S = S(1:2:end,1:2:end,1:2:end);

% load('F:\实测数据\data20181110剪刀剪线刀夹子块1z0.26x+0.01\matlab_3D.mat');

%% 补零计算
% S_lin=permute(S_lin,[2 3 1]);
tic;
%% 解回波
%% fft
S_ftxy=zeros(Nx,Ny,Nf);                                    % 对回波的x,y两个方位维进行fft
S_ftxy=fftshift(fft2(fftshift(S)));
%% 相位补偿
S_ftxy=S_ftxy.*exp(-1j*K*R0).*exp(1j*Kz*R0);
% S_ftxy=S_ftxy.*exp(1j*Kz*R0);
%% 插值部分（spaine）



a=zeros(Nx,Ny);
a=Kz(:,:,Nf);
[ax,ay]=find(a==max(max(a)));
KzM=max(max(max(Kz)));
Kzm=min(min(min(Kz)));                    %找到kz最大值,最小值，以及其坐标
Kzm=abs(Kzm);
Kzmid=Kz(ax(1),ay(1),1);
Nf2=fix(Nf*(KzM-Kzm)/(KzM-Kzmid))+1;      %建立的新矩阵的z轴大小，新矩阵的点数
% Nf2=Nf;
% Nf2=2^nextpow2(Nf2);
dnf=Nf2-Nf;
kz2=linspace(Kzm,KzM,Nf2);
Kz2=zeros(Nx,Ny,Nf2);
Kx2=zeros(Nx,Ny,Nf2);
Ky2=zeros(Nx,Ny,Nf2);
K2=zeros(Nx,Ny,Nf2);
for i=1:Nx
    for j=1:Ny
          Kz2(i,j,:)=kz2;                  %新的Kz矩阵（插值后的均匀矩阵）
    end
end
Ky2=repmat(Ky,[1 1 Nf2]);
Kx2=repmat(Kx,[1 1 Nf2]);
K2=sqrt(Kz2.^2+Kx2.^2+Ky2.^2); 



S2=zeros(Nx,Ny,Nf2);                      %
K3=zeros(Nx,Nf);                         % 原来的Kz矩阵
K4=zeros(Nx,Nf2);                        % 新的Kz矩阵       
S4=zeros(Nx,Nf);        
S6=zeros(Nx,Nf2);                                 
K3_uni=zeros(Nx,Nf);                    %插值过程中用量的预留空间
K4_uni=zeros(Nx,Nf2);                    %插值过程中用量的预留空间

dk1=(K(1,1,Nf)-K(1,1,1))/(Nf-1);                   %新矩阵的平均间隔
tic;
h=waitbar(0,'please wait');               %进度条
for i=1:Ny                                %将三维分为Ny个二维矩阵（因为kx，ky两个维度均匀）
K3(:,:)=K(:,i,:);
K4(:,:)=K2(:,i,:);
S4(:,:)=S_ftxy(:,i,:);        %fft后的回波的一个二维矩阵
     
K3_uni=((K3-K(1,1,1)*ones(Nx,Nf))/dk1)+1;   %原来矩阵的归一化
K4_uni=((K4-K(1,1,1)*ones(Nx,Nf2))/dk1)+1;  %新矩阵的归一化


    for j=1:Nx
    [x3,z3]=find(K4_uni(j,:)>K3_uni(j,Nf));
    if isempty(z3)==1
        z3=Nf2+1;
    end
    zz3=min(min(z3));
    [x4,z4]=find(K4_uni(j,:)<K3_uni(j,1));
    if isempty(z4)==1
        z4=0;
    end
    zz4=max(max(z4));
    x=K3_uni(j,:);
    y=S4(j,:);
    xs=K4_uni(j,zz4+1:zz3-1);
    ys=interp1(x,y,xs,'spline');    
    S6(j,zz4+1:zz3-1)=ys;
    end 
%computation here%
    waitbar(i/Ny,h);
    S2(:,i,:)=S6;                         %插值后的回波
end
toc;
delete(h);



%% ifft
S_iftxyz=zeros(Nx,Ny,Nf);
S_iftxyz=fftshift(ifftn(fftshift((S2))));
% S_iftxyz=fftshift(ifftn(fftshift((S_ftxy))));
toc;
%% 反演网格计算
ObjectX_pos=((-(Nx-1)/2:(Nx-1)/2)*deltaX).';  % 收发机方位维Y维的采样坐标
ObjectY_pos=((-(Ny-1)/2:(Ny-1)/2)*deltaY).';  % 收发机方位维Y维的采样坐标

%% 2
Dy=(Nf-1)*rhoz;               % 距离深度
Nf2=Nf;
%% 4

ObjectZ_pos=linspace(-Dy/2,Dy/2,Nf2)+R0;
ObjectZ_pos1=ObjectZ_pos-R0;
ObjectZ2=0.3;
Nf3=round(ObjectZ2/Dy*Nf2);
if mod(Nf3,2)==1
    Nf3=Nf3+1;
else
    Nf3=Nf3;
end
ObjectZ_pos=linspace(-ObjectZ2/2,ObjectZ2/2,Nf2)+R0;
%% 构建成像图像
toc;% 计时结束

%% 采用高频带图像融合的方式进行图像处理；之后将图像融合的方式变为function，从不同利用不同方向上的切片来实现目标的三维成像
%% 正视图
% xy平面图像
% N1=(Nf2-Nf3)/2+1;
% N1=fix(N1);
% N2=N1+Nf3-1;
% N2=fix(N2);
% S_iftxyz(:,:,1:N1)=0;
% S_iftxyz(:,:,N2:Nf2)=0;

%% 数据调试接口
% S_iftxyz = load('S_iftxyz1.mat');

max_xy = max(permute(S_iftxyz,[1 2 3]),[],3);
aaa=abs(max_xy.');
% bbb=zeros(Ny,Nx);
% bbb(1:Nymove,:)=aaa((125-Nymove+1):125,:);
% bbb((Nymove+1):125,:)=aaa(1:(125-Nymove),:);
figure
% imagesc(ObjectX_pos,ObjectY_pos,fliplr(aaa));             
imagesc(ObjectX_pos,ObjectY_pos,abs(max_xy.'));
% colormap(gray);
xlabel('x(m)'),ylabel('y(m)');
set(gca,'YDir','normal')
title('主视图');
% frame=getframe(gca);
% imwrite(frame.cdata,'正视图.png');


%% 

% xz平面图像
max_xz = max(permute(S_iftxyz,[1 3 2]),[],3);
% max_xz2=max_xz(:,N1:N2).';
figure
% imagesc(ObjectX_pos,ObjectZ_pos,abs(fliplr(max_xz2)));          % 距离维指定区域
imagesc(ObjectX_pos,ObjectZ_pos1,abs(max_xz.'));                         % 距离维全部区域
% colormap(gray);
xlabel('x(m)'),ylabel('z(m)');
set(gca,'YDir','normal')
title('侧视图');
% set(gca,'XTick',[],'YTick',[])%去掉横纵坐标的刻度线
% frame=getframe(gca);
% imwrite(frame.cdata,'侧视图.png');


% yz平面图像
max_yz = max(permute(S_iftxyz,[2 3 1]),[],3);
aa=abs(max_yz.');
% bb=zeros(Nf2,Ny);
% bb(:,1:Nymove)=aa(:,(125-Nymove+1):125);
% bb(:,(Nymove+1):125)=aa(:,1:(125-Nymove));
% aa=aa(N1:N2,:);
figure
% imagesc(ObjectY_pos,ObjectZ_pos1,aa);                 % 距离维指定区域
imagesc(ObjectX_pos,ObjectZ_pos1,abs(max_yz.'));                        % 距离维全部区域
% % colormap(gray);
xlabel('y(m)'),ylabel('z(m)');
set(gca,'YDir','normal')
title('俯视图');
% frame=getframe(gca);
% imwrite(frame.cdata,'俯视图.png');

toc;


% figure;scatter3(Ptar(:,1),Ptar(:,2),Ptar(:,3),'.');
% xlabel('x(m)'),ylabel('y(m)'),zlabel('z(m)');
% save .\echo_data\S_iftxyz1.mat S_iftxyz
% save .\echo_data\S1.mat S