%% 基于前处理方式的SAR成像仿真   （参数(采样率)不同不重叠）调频斜率相同只需要插值即可
clc
clear all
close all
%% 参数配置
c=3e8;
%% 脉冲1
fc1=10e9;
B1=500e6;
Tr1=3.3e-7;
K1=B1/Tr1;  %调频斜率  ：  1.5152e15
Ny1=1024;
delta_y1=c/(2*B1);%分辨率

%分辨率
y_grid1=(0:Ny1-1)*delta_y1; 

%% 脉冲2
fc2=10.5e9;
B2=500e6;
Tr2=3.3e-7;
K2=B2/Tr2;    %调频斜率  ：  1.1364e15
Ny2=512;
delta_y2=c/(2*B2);%分辨率

%分辨率
y_grid2=(0:Ny2-1)*delta_y2;

%% 合成脉冲
fc=10.225e9;
B=1e9;
Tr=2*3.3e-7;
K=B/Tr;
% Ny=Ny1+Ny2-1;
% Ny=Ny1+Ny2;
Ny=2048;
delta_y=c/(2*B);

%分辨率、
y_grid=(0:Ny-1)*delta_y;


%% 目标点位置设置
%   10  10.5   17.7   18
%% 设定收发天线的位置
x_TR=0;
y_TR=0;
z_TR=0;

R_ref=0;%设参考目标位置

% 时间序列
t1=linspace(0,Tr1,Ny1);
t2=linspace(0,Tr2,Ny2);
t22=linspace(0,Tr2,1024);
t=linspace(0,Tr,2048);
% t2=linspace(0,Tr2,Ny2);

x_target=0;
y_target=0;

%% 脉冲1
%% 解调频参考信号
s_if1=0;%初始化中频信号 
phase_ref1=2*pi*(fc1*(t1-(2*R_ref/c))+K1*((t1-(2*R_ref/c)).^2)/2);
s_ref1=exp(-1j*phase_ref1);

%% 目标回波构造

z_target=[10;11;17.7;18];  % 目标位置 
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





% %% 真实宽带
% %% 解调频参考信号
% s_if=0;%初始化中频信号 
% phase_ref=2*pi*(fc*(t-(2*R_ref/c))+K*((t-(2*R_ref/c)).^2)/2);
% s_ref=exp(-1j*phase_ref);
% 
% %% 目标回波构造
% 
% % z_target=[10;11;17.8;18];  % 目标位置 
% sigma=ones(4,1);       % 目标强度
% for i=1:4
%    
% R_i=sqrt((x_TR-x_target)^2+(y_TR-y_target)^2+(z_TR-z_target(i)).^2);%目标到收发天线之间的距离
% phase_r=2*pi*(fc*(t-(2*R_i/c))+K*((t-(2*R_i/c)).^2)/2);
% s_r=sigma(i)*exp(-1j*phase_r);
% 
% s_if=s_if+s_r;
% end  
% 
% s_dcp=s_if.*conj(s_ref);%进行差频处理
% 
% % f2=fc2+K2*t2;
% s_dcp=fft(s_dcp);
% % s_compa2=s_dcp2.*exp(-1j*pi*f2.^2/K2);%%%
% s_ift=ifft(s_dcp);
% 
% %% 真实宽带距离像
% G_pc=fft(s_ift)./max(fft(s_ift));
% figure
% plot(y_grid,(abs(G_pc)));
% title('真实宽带一维距离维成像仿真图');
% xlabel('距离/m');
% ylabel('归一化幅度');



%% 合成宽带及插值处理 距离维采样点数1024

% [x1,y1]=size(s_ift1);
% [x2,y2]=size(s_ift2);
% S_ift=zeros(1,y1+y2-1);
% S_ift(:,1:1024)=s_ift1;
% S_ift(:,1025:end)=s_ift2(:,2:end);

%将脉冲2的回波采样间隔插值与脉冲1一致，这样合成之后的就是均匀采样了

s_ift2=interp1(t2,s_ift2,t22,'spline');

S_ift=[s_ift1 s_ift2];

% 插值变成均匀序列



% s_dcp=[s_dcp1 s_dcp2];

s_dcp=fft(S_ift);
% s_compa3=s_dcp3.*exp(-1j*pi*fc.^2/K);%%%
s_ift=ifft(s_dcp);

%% 合成脉冲距离像
G_pc=fft(s_ift)./max(fft(s_ift));
figure
plot(y_grid,(abs(G_pc)));
title('合成脉冲一维距离维成像仿真图');
xlabel('距离/m');
ylabel('归一化幅度');

