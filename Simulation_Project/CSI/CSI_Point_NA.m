%% 这个脚本用来仿真带NA的白光干涉过程
clear;
close all;
clc;
%% 定义光源
z_peri=0.0719;%PZT单次移动间隔，um制
% z_samp=0.05;%PZT单次移动间隔，um制

N_half = round(10 / z_peri); % 计算单侧点数
z_scan = (-N_half : N_half) * z_peri; % 这样生成的数组中心绝对是 0

fre = make_axis_freq(length(z_scan),z_peri,'101');
fre = fre(fre> 2/1.1);%只选取大于1.1的部分；
vk0 = (fre /2)';

vsk_ini = gen_lightsource(vk0,1);%获取光源信号
vsk = vsk_ini./(vk0.^2); %波长波数域转换
vsk = vsk./(max(abs(vsk))); %归一化光谱
k0 = sum(vk0 .* vsk) / sum(vsk); %0NA对应的光谱中心

figure();
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
NA = 0; %系统NA
theta_max=asin(NA); %最大NA对应的空气中光线角度theta
theta_peri=0.01; %角度theta的采样周期
theta_array = 0:theta_peri:theta_max; %theta采样数组
%% 样品反射率
% sample_stru = {'SiO2',5;...
%                'Si',inf}; %样品结构
sample_stru = {'Si',inf}; %样品结构
[r_Se,r_Sm] = CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru); %计算样品TE场、TM场反射率

% figure();
% surf(repmat(vk0,1,length(theta_array)),repmat(theta_array,length(vk0),1),abs(r_Se).^2);
% shading interp;
% defaultAxes(3);
% defaultColor(1);
% xlabel('k/$\mu m^{-1}$','Interpreter','latex');
% ylabel('$\theta$/rad','Interpreter','latex');

%% 参考镜反射率
[r_Me,r_Mm] = CalcMirrorAmplitudeReflectivity(vk0,theta_array); %计算参考镜TE场、TM场反射率

% figure();
% surf(repmat(vk0,1,length(theta_array)),repmat(theta_array,length(vk0),1),abs(r_Me).^2);
% shading interp;
% defaultAxes(3);
% defaultColor(1);
% xlabel('k/$\mu m^{-1}$','Interpreter','latex');
% ylabel('$\theta$/rad','Interpreter','latex');
%% 选择偏振模式，生成白光干涉信号
% system_pol = 'unpolar';%非偏振模式
system_pol = 'ideal';%理想偏振模式

sample_dis = 1;%样品与参考镜之间的距离;
tic;
signal = CSIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,theta_array,sample_dis,system_pol);%生成CSI信号
toc;

figure();
plot(z_scan,signal,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');

%信号加高斯白噪声
SNR = 40; %40dB的噪声
signal = awgn(signal,SNR,'measured');
signal = signal/max(abs(signal)); %归一化

figure();
plot(z_scan,signal,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');
%% 点式FDA算法
[FDA,LinearFit,Hilb]=CSIPointFDA(signal,z_scan,z_peri,k0,NA);
%%
disp(['预设高度h:',num2str(sample_dis)]);
disp(['希尔伯特变换+重心法:',num2str(Hilb)]);
disp(['线性拟合高度h:',num2str(LinearFit)]);
disp(['FDA处理高度h:',num2str(FDA)]);