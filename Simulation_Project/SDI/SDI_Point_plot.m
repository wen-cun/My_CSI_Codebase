%% 这个脚本用来绘制有无NA、SDI的信号变化。
clear;
close all;
clc;
%% %光源部分

lam_samp=3e-4;%光谱仪分辨率0.3nm;
lam_lim=[0.3,0.8];%光谱仪采样范围
lam=lam_lim(1):lam_samp:lam_lim(2);
vk0=(1./lam)';
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
%% 系统参数
A1=1;%反射光振幅
A2=0.25;%参考光振幅
NA=0.3;%带NA的系统
d=0.5;%样品与参考镜的距离
%% NA=0的SDI信号

Modu_0=A1^2+A2^2+2*A1*A2.*cos(4*pi*vk0*d);
I_0=vsk.*Modu_0;
I_0=I_0./max(I_0);

Color=GetColor(2,1);
figure(2);
plot(1./vk0,I_0,'--','LineWidth',1.5,'Color',Color(1,:));
defaultAxes(2);
xlabel('$\lambda$/$\mu$m','Interpreter','latex');
%% 带NA的SDI信号
theta_samp=0.01;
theta=0:theta_samp:asin(NA);
theta_mat=repmat(theta,length(vk0),1);%行向量复制为矩阵
vk0_mat=repmat(vk0,1,length(theta));%列向量复制为矩阵
Modu_NA=(A1^2+A2^2+2*A1*A2.*cos(4*pi.*vk0_mat*d.*cos(theta_mat))).*sin(theta_mat).*cos(theta_mat);
Modu_NA=sum(Modu_NA,2);%每行求和
I_NA=vsk.*Modu_NA;
I_NA=I_NA./max(I_NA);

figure(2);
hold on;
plot(1./vk0,I_NA,'LineWidth',1.5,'Color',Color(2,:));
hold off;
defaultAxes(2);
xlabel('$\lambda$/$\mu$m','Interpreter','latex');