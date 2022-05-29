clc
clear
close all
%% main
img_dir = 'D:\E\浏览器下载\第一期\20220110_1_Code\LR\LRBI\pseudo1\x8\';
img_name = 'image1_LRBI_x8.png';
resolx = 0.01;
resoly = 0.01;
target = create_2decho(img_dir,img_name,resolx,resoly);
num=size(target,1);
%% 成像参数设置
c = 3e8;
R0 = 0.57;   %场景中心距离阵列平面的距离
fmin = 75e9;
% fmax = 20e9;
B = 35e9;   %信号带宽
fc = fmin+0.5*B;
lambda = c/fc;   %中心频率处波长
Nx = 200;
thetaX_span = 60;    %天线波束角，单位：°
rho_x = lambda/2;
Lx = (Nx-1)*rho_x;   %x方向维孔径大小
Nf = 200;
rho_f = c/2/B;
Lf = (Nf-1)*rho_f;
f_step = B/(Nf-1);            % 采样间隔
f = fmin+(0:Nf-1)*f_step;     % 采样频率

x_tr = ((-(Nx-1)/2:(Nx-1)/2)*rho_x).';     
y_tr=-R0.*ones(Nx,1);                        

%% 构建回波
j=sqrt(-1);
for i=1:Nx
    s=zeros(1,Nf);
            for j1=1:num                                       %目标的数量
               x=target(j1,2);                                 %目标的横坐标
               y=target(j1,1);                                 %目标的纵坐标
               A=target(j1,3);                                 %目标的幅度
               R=sqrt((x_tr(i)-x).^2+(y_tr(i)-y).^2);        %天线到目标的距离
               s=s+A*exp(-j*2*pi*f*2*R/c);                    
            end
     S(i,:)=s;
end

[S_iftxyz,Nf1]=dataprocess(S,rho_x,f,Nx,Nf,R0);

% deltaX=(c/(fmin+0.5*B))/2;
% Dy=(c/2/B)*(Nf-1);
% deltaX = 0.0019;
% Dy = 0.0416;
deltaX = rho_x;
Dy = Lf;
ObjectX_pos=((-(Nx-1)/2:(Nx-1)/2)*deltaX);  % 收发机方位维X维的采样坐标
ObjectY_pos=linspace(-Dy/2,Dy/2,Nf1);  % 收发机方位维Y维的采样坐标


figure,
imagesc(ObjectX_pos',ObjectY_pos,flipud(abs(S_iftxyz')));

% imagesc(ObjectX_pos1',ObjectY_pos1,abs(S_iftxyz1'));
title('反演回波成像');
xlabel('方位维/m'),ylabel('距离维/m');
set(gca, 'YDir', 'normal');
axis off;   %不显示坐标轴
% frame=getframe(gca);
% imwrite(frame.cdata,'1.png');
% colormap(gray);


function target = create_2decho(img_dir,img_name,resolx,resoly)
%% 图像反演回波点位置
im_RGB = imread(fullfile(img_dir,img_name));
im_gray = rgb2gray(im_RGB);
im = double(im_gray);
im_norm = im/max(max(im));
[Px,Py] = size(im_norm);
pixelx_pos = ((-(Px-1)/2:(Px-1)/2)*resolx).'; 
pixely_pos = ((-(Py-1)/2:(Py-1)/2)*resoly).'; 
target = [];
for i= 1:Px
    for j = 1:Py
        if im_norm(i,j)>0 && im_norm(i,j)~=im_norm(1,1) && im_norm(i,j)~=im_norm(end,end)
            target1 = [pixelx_pos(i),pixely_pos(j),im_norm(i,j)];
            target = [target;target1];
        end
    end
end

end


function [S_iftxyz,Nf2]=dataprocess(S,deltaX,f,Nx,Nf,R0)
%%   RMA成像关键部分
j=sqrt(-1);
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
S_FT=[];
S_FT=fftshift(fft(fftshift(S,1),[],1),1);    %傅里叶变换   方位维傅里叶变换

%% 调试图
%         figure,
%         plot(Kx,Ky,'-');
%         xlabel('Kx');ylabel('Ky');
%         title('波数域');
%          hold on

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

%% 调试图
%         figure,
%         plot(kx2,ky2,'.');
%         xlabel('Kx');ylabel('Ky');
%         title('波数域');


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
    ys=interp1(abs(x11),y11,abs(xs),'spline');    %%加了abs 转换为实数
    S1(i3,zz3+1:zz2-1)=ys;
end

S_iftxyz=zeros(Nx,Nf2);
S_iftxyz=fftshift(ifft2(fftshift((S1))));
end