%% %这个脚本用来仿真NA=0的白光干涉过程
clear;
close all;
clc;
%% 1. 定义光谱曲线
z_samp=0.0719;%PZT单次移动间隔，um制
N_half = round(10 / z_samp); % 计算单侧点数
z_scan = (-N_half : N_half) * z_samp; % 这样生成的数组中心绝对是 0
fre=make_axis_freq(length(z_scan),z_samp,'101');
fre=fre(fre>2/1.1);%只选取大于1.1的部分；
vk0_ini=(fre/2)';
vk0_samp=1/(z_samp*length(z_scan));
[vk0,vSk_ini]=genLightSource(vk0_ini,10);%获取光源信号
vSk=vSk_ini./(vk0.^2);%波长波数域转换

% 转换为 single 以节省后续计算内存
vk0 = single(vk0);
vSk = single(vSk);


figure(1);
plot(vk0,vSk,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu m^{-1}$','Interpreter','latex');     