%% 基于前处理方式的SAR成像仿真(参数不同无重叠)
clc
clear all
close all
%% 参数配置
c=3e8;
%% 脉冲1
fc1=10e9;
B1=500e6;
Tr1=3.3e-7;
K1=B1/Tr1;  %调频斜率
Ny1=1024;
delta_y1=c/(2*B1);%分辨率

%分辨率
y_grid1=(0:Ny1-1)*delta_y1; 

%% 脉冲2
fc2=10.5e9;
B2=300e6;
Tr2=3.3e-7;
K2=B2/Tr2;    %调频斜率
Ny2=512;
delta_y2=c/(2*B2);%分辨率

%分辨率
y_grid2=(0:Ny2-1)*delta_y2;

%% 合成脉冲
fc=10.15e9;
B=0.8e9;
% Ny=Ny1+Ny2;
Ny=2000;
delta_y=c/(2*B);

%分辨率、
y_grid=(0:Ny-1)*delta_y;


%% 目标点位置设置
%   0  2   7.8   8
%% 设定收发天线的位置
x_TR=0;
y_TR=0;
z_TR=0;

R_ref=0;%设参考目标位置

% 时间序列
t1=linspace(0,Tr1,Ny1);
t11=linspace(0,Tr1,1250);
f1=fc1+K1*t1;
% t2=linspace(Tr1,Tr1+Tr2,Ny2);
t2=linspace(0,Tr2,Ny2);
t22=linspace(0,Tr2,750);
f2=fc2+K2*t2;

x_target=0;
y_target=0;

%% 脉冲1
%% 解调频参考信号
s_if1=0;%初始化中频信号 
phase_ref1=2*pi*(fc1*(t1-(2*R_ref/c))+K1*((t1-(2*R_ref/c)).^2)/2);
s_ref1=exp(-1j*phase_ref1);

%% 目标回波构造

z_target=[10;11;17.65;18];  % 目标位置 
sigma=ones(4,1);       % 目标强度
for i=1:4
   
R_i=sqrt((x_TR-x_target)^2+(y_TR-y_target)^2+(z_TR-z_target(i)).^2);%目标到收发天线之间的距离
phase_r1=2*pi*(fc1*(t1-(2*R_i/c))+K1*((t1-(2*R_i/c)).^2)/2);
s_r1=sigma(i)*exp(-1j*phase_r1);

s_if1=s_if1+s_r1;
end  

s_dcp1=s_if1.*conj(s_ref1);%进行差频处理

% f1=fc1+K1*t1;
s_dcp1=fft(s_dcp1);
% s_compa1=s_dcp1.*exp(-1j*pi*f1.^2/K1);%%%
s_ift1=ifft(s_dcp1);

%% 子脉冲1距离像
G_pc1=fft(s_ift1)./max(fft(s_ift1));
figure
plot(y_grid1,(abs(G_pc1)));
title('子带1一维距离维成像仿真图');
xlabel('距离/m');
ylabel('归一化幅度');


%% 脉冲2
%% 解调频参考信号
s_if2=0;%初始化中频信号 
phase_ref2=2*pi*(fc2*(t2-(2*R_ref/c))+K2*((t2-(2*R_ref/c)).^2)/2);
s_ref2=exp(-1j*phase_ref2);

%% 目标回波构造

% z_target=[10;11;17.8;18];  % 目标位置 
sigma=ones(4,1);       % 目标强度
for i=1:4
   
R2_i=sqrt((x_TR-x_target)^2+(y_TR-y_target)^2+(z_TR-z_target(i)).^2);%目标到收发天线之间的距离
phase_r2=2*pi*(fc2*(t2-(2*R2_i/c))+K2*((t2-(2*R2_i/c)).^2)/2);
s_r2=sigma(i)*exp(-1j*phase_r2);

s_if2=s_if2+s_r2;
end  

s_dcp2=s_if2.*conj(s_ref2);%进行差频处理

% f2=fc2+K2*t2;
s_dcp2=fft(s_dcp2);
% s_compa2=s_dcp2.*exp(-1j*pi*f2.^2/K2);%%%
s_ift2=ifft(s_dcp2);

%% 子脉冲2距离像
G_pc2=fft(s_ift2)./max(fft(s_ift2));
figure
plot(y_grid2,(abs(G_pc2)));
title('子带2一维距离维成像仿真图');
xlabel('距离/m');
ylabel('归一化幅度');



%% 合成宽带及插值处理 距离维采样点数128
ff1=linspace(10e9,10.5e9,1250);
ff2=linspace(10.5e9,10.8e9,750);
s_ift1=interp1(f1,s_ift1,ff1,'nearest');
s_ift2=interp1(f2,s_ift2,ff2,'nearest');

S_ift=[s_ift1 s_ift2];


s_dcp3=fft(S_ift);
% s_compa2=s_dcp2.*exp(-1j*pi*f2.^2/K2);%%%
s_ift3=ifft(s_dcp3);

%% 合成脉冲距离像
G_pc3=fft(s_ift3)./max(fft(s_ift3));
figure
plot(y_grid,(abs(G_pc3)));
title('合成宽带一维距离维成像仿真图');
xlabel('距离/m');
ylabel('归一化幅度');

