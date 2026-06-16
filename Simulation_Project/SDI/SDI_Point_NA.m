%% 仿真带NA的SDI信号
clear;
close all;
clc;
%% 定义光源
lam_peri = 3e-4; %波长采样周期,um制
lam_lim = [0.3,1.1]; %波长采样范围
lam = lam_lim(1):lam_peri:lam_lim(2); %波长采样
lam = lam(:); %转为列向量
vk0 = 1./lam; %转为波数
vsk_ini = gen_lightsource(vk0,1);%获取光源信号
vsk = vsk_ini./(vk0.^2);%波长波数域转换
vsk=vsk./(max(abs(vsk)));

figure();
tiledlayout(3,1);

nexttile;
plot(lam,vsk_ini,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('$\lambda$/$\mu$m','Interpreter','latex');

nexttile;
plot(vk0,vsk_ini,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu m^{-1}$','Interpreter','latex');

nexttile;
plot(vk0,vsk,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu m^{-1}$','Interpreter','latex');
%% 定义角度
NA = 0.3; %系统NA
theta_max=asin(NA); %最大NA对应的空气中光线角度theta
theta_peri=0.01; %角度theta的采样周期
theta_array = 0:theta_peri:theta_max; %theta采样数组
%% 定义样品结构
% sample_stru = {'SiO2',5;...
%                'Si',inf}; %样品结构
sample_stru = {'Si',inf}; %样品结构
[r_Se,r_Sm] = CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru); %计算样品TE场、TM场反射率
%% 参考镜反射率
[r_Me,r_Mm] = CalcMirrorAmplitudeReflectivity(vk0,theta_array); %计算参考镜TE场、TM场反射率号
%% 选择偏振模式，生成白光干涉信号
system_pol = 'unpolar';%非偏振模式

signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,sample_dis,system_pol);

figure();
plot(vk0,signal,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu m^{-1}$','Interpreter','latex');

figure();
plot(lam,signal.*vk0.^2,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('$\lambda$/$\mu m$','Interpreter','latex');
%%
interD=(signal+eps)./(vsk+eps);

figure();
plot(vk0,interD,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu m^{-1}$','Interpreter','latex');
