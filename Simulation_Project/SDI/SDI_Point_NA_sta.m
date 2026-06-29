%% 这个脚本用来仿真SDI处理算法的准确度、重复测量稳定性
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
sample_dis=6;
signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,sample_dis,system_pol);
signal_pur=signal;

num_meas = 20; %重复测量的次数
z_pre = nan*ones(num_meas,1); %预先分配内存

tic;
for ii=1:num_meas
disp(['正在计算:',num2str(ii),'/',num2str(num_meas)]);
SNR = 40; %40dB的噪声
signal = awgn(signal_pur,SNR,'measured');
signal = signal/max(abs(signal)); %归一化
[z_coa,rsquare] = SDIPointModulFit(signal,lam,vsk_ini,NA);
z_pre(ii) = SDIPointModelFit(signal,z_coa,rsquare,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol);
end
toc;
%% 展示结果
Color=GetColor(2,1);

figure();
plot(1:num_meas,sample_dis*ones(1,num_meas),'LineWidth',1.5,'Color',Color(1,:));
hold on;
scatter(1:num_meas,z_pre,25,Color(2,:),'filled');
hold off;
defaultAxes(2);
legend('Real','Measurement','EdgeColor','none');
xlabel('num');

z_pre_mean = mean(z_pre); %计算平均值
z_pre_std = std(z_pre); %计算标准差
z_error = abs(z_pre_mean-sample_dis) / sample_dis; %计算平均偏差

disp(['预设值:',num2str(sample_dis)]);
disp(['平均测量值:',num2str(z_pre_mean)]);
disp(['平均偏差:',num2str(z_error)]);
disp(['标准差:',num2str(z_pre_std)]);