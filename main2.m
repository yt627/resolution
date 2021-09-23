%% 基于前处理方式的SAR成像仿真    (参数*（调频斜率和采样率）不同有重叠)需要对重叠部分进行排序再插值
clc
clear all
close all
%% 参数配置
c=3e8;
%% 脉冲1
fc1=10e9;
B1=500e6;
Tr1=5.3e-7;
K1=B1/Tr1;  %调频斜率
Ny1=1024;
delta_y1=c/(2*B1);%分辨率

%分辨率
y_grid1=(0:Ny1-1)*delta_y1; 

%% 脉冲2
fc2=10.35e9;
B2=300e6;
Tr2=3.3e-7;
K2=B2/Tr2;    %调频斜率
Ny2=512;
delta_y2=c/(2*B2);%分辨率

%分辨率
y_grid2=(0:Ny2-1)*delta_y2;

%% 合成脉冲
fc=9.75e9;
B=0.8e9;
Tr3=6.6e-7;
K=B/Tr3;
Ny=2500;
Ny3=2500;
delta_y=c/(2*B);

%分辨率、
y_grid=(0:Ny-1)*delta_y;
y_grid3=(0:Ny3-1)*delta_y;


%% 目标点位置设置
%   10  11   17.75   18
%% 设定收发天线的位置
x_TR=0;
y_TR=0;
z_TR=0;

R_ref=0;%设参考目标位置

% 时间序列
t1=linspace(-Tr1/2,Tr1/2,Ny1);
% t11=linspace(0,Tr1,1250);
f1=fc1+K1*t1;

t2=linspace(-Tr2/2,Tr2/2,Ny2);
% t22=linspace(0,Tr2,750);
f2=fc2+K2*t2;

t=linspace(-(Tr1+Tr2)/2,(Tr1+Tr2)/2,Ny1+Ny2);
tt=linspace(-(Tr1+Tr2)/2,(Tr1+Tr2)/2,Ny);
t3=linspace(-Tr3/2,Tr3/2,Ny3);
f3=fc+K*t3;

x_target=0;
y_target=0;

%% 脉冲1
%% 解调频参考信号
s_if1=0;%初始化中频信号 
phase_ref1=2*pi*(fc1*(t1-(2*R_ref/c))+K1*((t1-(2*R_ref/c)).^2)/2);
s_ref1=exp(-1j*phase_ref1);

%% 目标回波构造
z_target=[17.7;18];  % 目标位置 
sigma=ones(2,1);       % 目标强度
for i=1:2
   
R_i=sqrt((x_TR-x_target)^2+(y_TR-y_target)^2+(z_TR-z_target(i)).^2);%目标到收发天线之间的距离
phase_r1=2*pi*(fc1*(t1-(2*R_i/c))+K1*((t1-(2*R_i/c)).^2)/2);
s_r1=sigma(i)*exp(-1j*phase_r1);

s_if1=s_if1+s_r1;
end  

s_dcp1=s_if1.*conj(s_ref1);%进行差频处理

s_dcp1=fft(s_dcp1);
s_compa1=s_dcp1.*exp(-1j*pi*fc1.^2/K1);%%%
s_ift1=ifft(s_compa1);

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
sigma=ones(2,1);       % 目标强度
for i=1:2
   
R2_i=sqrt((x_TR-x_target)^2+(y_TR-y_target)^2+(z_TR-z_target(i)).^2);%目标到收发天线之间的距离
phase_r2=2*pi*(fc2*(t2-(2*R2_i/c))+K2*((t2-(2*R2_i/c)).^2)/2);
s_r2=sigma(i)*exp(-1j*phase_r2);

s_if2=s_if2+s_r2;
end  

s_dcp2=s_if2.*conj(s_ref2);%进行差频处理

s_dcp2=fft(s_dcp2);
s_compa2=s_dcp2.*exp(-1j*pi*fc2.^2/K2);%%%
s_ift2=ifft(s_compa2);

%% 子脉冲2距离像
G_pc2=fft(s_ift2)./max(fft(s_ift2));
figure
plot(y_grid2,(abs(G_pc2)));
title('子带2一维距离维成像仿真图');
xlabel('距离/m');
ylabel('归一化幅度');



%% 合成宽带及插值处理 距离维采样点数128
[x1,y1]=size(f1);
[x2,y2]=size(f2);

f=[f1 f2];% 索引号大于y1为脉冲2的排序
[f_sort,index]=sort(f);%对频率序列进行排序
S_ift=[s_ift1 s_ift2];
[x,y]=size(S_ift);
S_ift_sort=zeros(1,y);
for i=1:y
    S_ift_sort(i)=S_ift(index(i));
    
end

ff=linspace(f_sort(1),f_sort(end),Ny);


% real=real(S_ift_sort);
% imag=imag(S_ift_sort);
% real1=interp1(f_sort,real,ff);
% imag1=interp1(f_sort,imag,ff);
% S_ift_sort1=real1+1j*imag1;

S_ift_sort1=interp1(f_sort,S_ift_sort,ff,'spline');



% s_ift1=interp1(t1,s_ift1,t11,'cubic');
% s_ift2=interp1(t2,s_ift2,t22,'cubic');

% S_ift=[s_ift1 s_ift2];


s_dcp3=fft(S_ift_sort1);
% s_compa2=s_dcp2.*exp(-1j*pi*f2.^2/K2);%%%
s_ift3=ifft(s_dcp3);

%% 合成脉冲距离像
G_pc3=fft(s_ift3)./max(fft(s_ift3));
figure
plot(y_grid,(abs(G_pc3)));
title('合成宽带一维距离维成像仿真图');
xlabel('距离/m');
ylabel('归一化幅度');





%% 真实宽带距离像
%% 解调频参考信号
s_if3=0;%初始化中频信号 
phase_ref3=2*pi*(fc*(t3-(2*R_ref/c))+K*((t3-(2*R_ref/c)).^2)/2);
s_ref3=exp(-1j*phase_ref3);

%% 目标回波构造
sigma=ones(2,1);       % 目标强度
for i=1:2
   
R3_i=sqrt((x_TR-x_target)^2+(y_TR-y_target)^2+(z_TR-z_target(i)).^2);%目标到收发天线之间的距离
phase_r3=2*pi*(fc*(t3-(2*R3_i/c))+K*((t3-(2*R3_i/c)).^2)/2);
s_r3=sigma(i)*exp(-1j*phase_r3);

s_if3=s_if3+s_r3;
end  

s_dcp3=s_if3.*conj(s_ref3);%进行差频处理

s_dcp3=fft(s_dcp3);
s_compa3=s_dcp3.*exp(-1j*pi*fc.^2/K);%%%
s_ift3=ifft(s_dcp3);

%% 真实带宽距离像
G_pc3=fft(s_ift3)./max(fft(s_ift3));
figure
plot(y_grid3,(abs(G_pc3)));
title('真实带宽一维距离维成像仿真图');
xlabel('距离/m');
ylabel('归一化幅度');


%%
s_ift3=s_ift3(1:10:end);
S_ift_sort1=S_ift_sort1(1:10:end);

phase1=imag(s_ift3);
real1=real(s_ift3);
phase2=imag(S_ift_sort1);
real2=real(S_ift_sort1);
atan1=atan(phase1./real1);
% atan1=atan1(1:10:end);
atan2=atan(phase2./real2);
% atan2=atan2(1:10:end);
figure
f33=f3(1:10:end);
stem(f33,atan1);
hold on
fff=ff(1:10:end);
stem(fff,atan2);
xlabel('频率（Hz）');
ylabel('相位');
h=legend('真实宽带','合成宽带');

%% 误差分析
[~,ii]=size(atan1);
erro1=zeros(1,ii);
for i=1:ii
    erro1(i)=(atan1(i)-atan2(i))/atan1(i);
end
figure
stem(fff,erro1);
xlabel('频率（Hz）');
ylabel('相位误差');

% 幅度
erro2=zeros(1,ii);
for i=1:ii
    erro2(i)=(abs(s_ift3(i))-abs(S_ift_sort1(i)))/abs(s_ift3(i));
end
figure
stem(fff,erro2);
xlabel('频率（Hz）');
ylabel('幅度误差');



% figure
% stem(f3,abs(s_ift3));
% hold on
% stem(ff,abs(S_ift_sort1));
% h=legend('真实宽带','合成宽带');
