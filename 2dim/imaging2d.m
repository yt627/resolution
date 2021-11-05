%% 20211025编写，二维回波频带融合仿真 步进频
clc
clear
close all

%% 收发机参数设置

%收发机1参数设置
c=3e8;                                         %光速
fs1=10e9;                                      %频带1起始频率
B1=0.5e9;                                      %频带1带宽
Tr1=3.3e-7;                                    %频带1脉宽
K1=B1/Tr1;                                     %频带1调频斜率
Nx1=100;                                       %方位向采样点数
Nf1=128;                                       %频带1距离采样点数（频点）
%delta_x1=
%delta_f1=c/(2*B1);                             %频带1距离为分辨率
delta_f1=B1/(Nf1-1);                            %频带1频率采样间隔
f1=fs1+(0:Nf1-1)*delta_f1;                      %频带1频率序列

%% 1 1
theta_1ant1=60;                                 %天线波束角，单位：°
theta_ant1=theta_1ant1*pi/180;                  %转化为弧度制
Kxmax1=(4*pi*(fs1+B1/2)/c)*sin(theta_ant1/2);
deltaX1=pi/Kxmax1;                              % x方向上的采样间隔
Lx1=(Nx1-1)*deltaX1;
R1=Lx1/2/tan(theta_ant1/2);                     %天线阵列到目标区域中心的距离，单位：m
x_tr1 = ((-(Nx1-1)/2:(Nx1-1)/2)*deltaX1).';     %收发机高度维的采样坐标
y_tr1=-R1.*ones(Nx1,1);                         %收发机在距离维采样坐标


%收发机2参数设置
c=3e8;                                         %光速
fs2=10.35e9;                                   %频带2起始频率
B2=0.3e9;                                      %频带2带宽
Tr2=3.3e-7;                                    %频带2脉宽
K2=B2/Tr2;                                     %频带2调频斜率
Nx2=100;                                       %方位向采样点数
Nf2=128;                                       %频带2距离采样点数（频点）
%delta_x2=
%delta_f2=c/(2*B2);                            %频带2距离为分辨率
delta_f2=B2/(Nf2-1);                           %频带2频率采样间隔
f2=fs2+(0:Nf2-1)*delta_f2;                     %频带2频率序列


%% 2 2
theta_1ant2=60;                                %天线波束角，单位：°
theta_ant2=theta_1ant2*pi/180;                 %转化为弧度制
Kxmax2=(4*pi*(fs2+B2/2)/c)*sin(theta_ant2/2);
deltaX2=pi/Kxmax2;                             %x方向上的采样间隔
Lx2=(Nx2-1)*deltaX2;
R2=Lx2/2/tan(theta_ant2/2);                    %天线阵列到目标区域中心的距离，单位：m
x_tr2 = ((-(Nx2-1)/2:(Nx2-1)/2)*deltaX2).';    %收发机高度维的采样坐标
y_tr2=-R2.*ones(Nx2,1);                        %收发机在距离维采样坐标

%全频带参数设置
c=3e8;                                         %光速
fs=10e9;                                       %全频带起始频率
B=0.65e9;                                      %全频带带宽
Tr=6.6e-7;                                     %全频带脉宽
K=B/Tr;                                        %全频带调频斜率
Nx=100;                                        %方位向采样点数
Nf=300;                                       %全频带距离采样点数（频点）
%delta_x=
%delta_f=c/(2*B);                              %全频带距离为分辨率
delta_f=B/(Nf-1);                              %全频带频率采样间隔
f=fs+(0:Nf-1)*delta_f;                         %全频带频率序列

%% 3 3
theta_1ant=60;                       %天线波束角，单位：°
theta_ant=theta_1ant*pi/180;            %转化为弧度制
Kxmax=(4*pi*(fs+B/2)/c)*sin(theta_ant/2);
deltaX=pi/Kxmax;                         % x方向上的采样间隔
Lx=(Nx-1)*deltaX;
R=Lx/2/tan(theta_ant/2);                %天线阵列到目标区域中心的距离，单位：m
x_tr = ((-(Nx-1)/2:(Nx-1)/2)*deltaX).';   % 收发机高度维的采样坐标
y_tr=-R.*ones(Nx,1);                      %收发机在距离维采样坐标



%% 目标参数设置
j=sqrt(-1);
 object=[
         % 0.3,   0,   1
          0,     0,   1
          0,     -0.2, 1
          0,     -0.5, 1
          0,     -0.9, 1
        % -0.3,  0.7, 1
   ];      %点目标
num=size(object,1);
%% 回波仿真
% 频带1回波
S1=zeros(Nx1,Nf1);                            %空回波
for i1=1:Nx1
    s1=zeros(1,Nf1);
            for j1=1:num                                       %目标的数量
               x=object(j1,1);                                 %目标的横坐标
               y=object(j1,2);                                 %目标的纵坐标
               A=object(j1,3);                                 %目标的幅度
               R=sqrt((x_tr1(i1)-x).^2+(y_tr1(i1)-y).^2);        %天线到目标的距离
               s1=s1+A*exp(-j*2*pi*f1*2*R/c);
            end
     S1(i1,:)=s1;
end

% 频带2回波
S2=zeros(Nx2,Nf2);                            %空回波
for i1=1:Nx2
    s2=zeros(1,Nf2);
            for j1=1:num
               x=object(j1,1);                                 %目标的横坐标
               y=object(j1,2);                                 %目标的纵坐标
               A=object(j1,3);                                 %目标的幅度
               R=sqrt((x_tr2(i1)-x).^2+(y_tr2(i1)-y).^2);        %天线到目标的距离
               s2=s2+A*exp(-j*2*pi*f2*2*R/c);
            end
     S2(i1,:)=s2;
end

% 全频带回波
S=zeros(Nx,Nf);                            %空回波
for i1=1:Nx
    s=zeros(1,Nf);
            for j1=1:num
               x=object(j1,1);                                 %目标的横坐标
               y=object(j1,2);                                 %目标的纵坐标
               A=object(j1,3);                                 %目标的幅度
               R=sqrt((x_tr(i1)-x).^2+(y_tr(i1)-y).^2);        %天线到目标的距离
               s=s+A*exp(-j*2*pi*f*2*R/c);
            end
     S(i1,:)=s;
end




%% 回波数据融合
ff=[f1 f2];
[ff_sort,ff_index]=sort(ff);

SS=[S1 S2];
[m,n]=size(SS);
SSS=zeros(m,n);
for i=1:n
    SSS(:,i)=SS(:,ff_index(i));
end
%距离维插值为均匀序列
fff=linspace(ff_sort(1),ff_sort(end),Nf);
SSSS=zeros(m,Nf);
for i=1:m
    SSSS(i,:)=interp1(ff_sort,SSS(i,:),fff);
end

%% 成像
% 1
[S_iftxyz1,Nf11,~,~,~,~]=dataprocess(S1,deltaX1,f1,Nx1,Nf1,R1);   %调用dataprocess函数， 成像预处理，傅里叶变换+距离维插值

Dy1=(c/2/B1)*(Nf1-1);
ObjectX_pos1=((-(Nx1-1)/2:(Nx1-1)/2)*deltaX1);      % 收发机方位维X维的采样坐标
ObjectY_pos1=linspace(-Dy1/2,Dy1/2,Nf11);           % 收发机方位维Y维的采样坐标


figure,
imagesc(ObjectX_pos1',ObjectY_pos1,abs(S_iftxyz1'));
title('频带1');
xlabel('方位维/m'),ylabel('距离维/m');
set(gca, 'YDir', 'normal');
colormap(gray);

% 2
[S_iftxyz2,Nf22,~,~,~,~]=dataprocess(S2,deltaX2,f2,Nx2,Nf2,R2);   % 调用dataprocess函数，成像预处理，傅里叶变换+距离维插值

Dy2=(c/2/B2)*(Nf2-1);
ObjectX_pos2=((-(Nx2-1)/2:(Nx2-1)/2)*deltaX2);        % 收发机方位维X维的采样坐标
ObjectY_pos2=linspace(-Dy2/2,Dy2/2,Nf22);             % 收发机方位维Y维的采样坐标


figure,
imagesc(ObjectX_pos2',ObjectY_pos2,abs(S_iftxyz2'));
title('频带2');
xlabel('方位维/m'),ylabel('距离维/m');
set(gca, 'YDir', 'normal');
colormap(gray);

% 全
[S_iftxyz3,Nfff2,~,~,~,~]=dataprocess(S,deltaX,f,Nx,Nf,R);   %调用dataprocess函数， 成像预处理，傅里叶变换+距离维插值

Dy=(c/2/B)*(Nf-1);
ObjectX_pos=((-(Nx-1)/2:(Nx-1)/2)*deltaX);        % 收发机方位维X维的采样坐标
ObjectY_pos=linspace(-Dy/2,Dy/2,Nfff2);           % 收发机方位维Y维的采样坐标


figure,
imagesc(ObjectX_pos',ObjectY_pos,abs(S_iftxyz3'));
title('全频带');
xlabel('方位维/m'),ylabel('距离维/m');
set(gca, 'YDir', 'normal');
colormap(gray);


% 融合
[S_iftxyz33,Nfff,~,~,~,~]=dataprocess(SSSS,deltaX,f,Nx,Nf,R); % 调用dataprocess函数，成像预处理，傅里叶变换+距离维插值

Dy=(c/2/B)*(Nf-1);
ObjectX_pos=((-(Nx-1)/2:(Nx-1)/2)*deltaX);     % 收发机方位维X维的采样坐标
ObjectY_pos=linspace(-Dy/2,Dy/2,Nfff);         % 收发机方位维Y维的采样坐标


figure,
imagesc(ObjectX_pos',ObjectY_pos,abs(S_iftxyz33'));
title('融合频带');
xlabel('方位维/m'),ylabel('距离维/m');
set(gca, 'YDir', 'normal');
colormap(gray);


%% 强度图
figure
mesh(ObjectX_pos1',ObjectY_pos1,abs(S_iftxyz1'));
xlabel('方位维/m'),ylabel('距离维/m'),zlabel('强度');title('频带1');

figure
mesh(ObjectX_pos2',ObjectY_pos2,abs(S_iftxyz2'));
xlabel('方位维/m'),ylabel('距离维/m'),zlabel('强度');title('频带2');

figure
mesh(ObjectX_pos',ObjectY_pos,abs(S_iftxyz3'));
xlabel('方位维/m'),ylabel('距离维/m'),zlabel('强度');title('全频带');
figure
mesh(ObjectX_pos',ObjectY_pos,abs(S_iftxyz33'));
xlabel('方位维/m'),ylabel('距离维/m'),zlabel('强度');title('融合频带');