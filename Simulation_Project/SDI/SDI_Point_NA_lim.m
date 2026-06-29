%% 这个脚本用判断0.05~10um内是否能正确提取测量结果
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
sample_dis=[0:0.05:1.1,1.5:0.5:10,10.1:0.1:11];
z_pre = nan*ones(size(sample_dis)); %预先分配内存
tic;
for ii=1:length(sample_dis)
    disp(['正在计算:',num2str(ii),'/',num2str(length(sample_dis))]);
    signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,sample_dis(ii),system_pol);
    SNR = 40; %40dB的噪声
    signal = awgn(signal,SNR,'measured');
    signal = signal/max(abs(signal)); %归一化
    [z_coa,rsquare] = SDIPointModulFit(signal,lam,vsk_ini,NA);
    z_pre(ii) = SDIPointModelFit(signal,z_coa,rsquare,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol);
end
toc;
%% 展示结果
Color=GetColor(2,1);

figure();
plot(sample_dis,sample_dis,'LineWidth',1.5,'Color',Color(1,:));
hold on;
plot(sample_dis,z_pre,'*','LineWidth',1.5,'Color',Color(2,:));
hold off;
defaultAxes(2);
legend('Real','Measurement','EdgeColor','none','Location','NorthWest');
xlabel('z/$\mu$m','Interpreter','latex');

figure();
plot(sample_dis,z_pre-sample_dis,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu$m','Interpreter','latex');
