%% 这个脚本用来仿真NA=0的白光干涉过程，但是采用当前代码库的函数进行编写
clear;
close all;
clc;
%% 定义光源
z_samp=0.0719;%PZT单次移动间隔，um制
% z_samp=0.05;%PZT单次移动间隔，um制

N_half = round(10 / z_samp); % 计算单侧点数
z_scan = (-N_half : N_half) * z_samp; % 这样生成的数组中心绝对是 0

fre=make_axis_freq(length(z_scan),z_samp,'101');
fre=fre(fre>2/1.1);%只选取大于1.1的部分；
vk0=(fre/2)';
vk0_samp=1/(z_samp*length(z_scan));

vsk_ini=gen_lightsource(vk0,1);%获取光源信号
vsk=vsk_ini./(vk0.^2);%波长波数域转换

figure(1);
tiledlayout(3,1);

nexttile;
plot(1./vk0,vsk_ini,'LineWidth',1.5,'Color',GetColor(1,1));
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