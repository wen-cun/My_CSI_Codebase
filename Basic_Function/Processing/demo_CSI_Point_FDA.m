%% 这个脚本用来确定点式的CSI FDA算法
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
vsk = vsk_ini./(vk0.^2);%波长波数域转换
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

sample_dis = 3;%样品与参考镜之间的距离;
signal = CSIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,theta_array,sample_dis,system_pol);%生成CSI信号

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
% signal = (signal-mean(signal(1:10)))'; %减去直流分量，并转为列向量
% signal=signal(:);
% signal=signal-mean(signal);
% 
% clip_num = 15; %两侧各取10各点，注意和信号包络宽度有关
% env_signal = abs(hilbert(signal)); % 希尔伯特变换提取信号包络
% [~, I_maxindex] = max(env_signal); % 信号包络最大值对应的索引
% 
% env_clip = env_signal(I_maxindex-clip_num:I_maxindex+clip_num); %选取部分包络，使用重心法计算信号中心（样品与参考镜之间的距离）
% z_clip = (z_scan(I_maxindex-clip_num:I_maxindex+clip_num))'; %选取用于计算的信号
% z_center = sum(z_clip.*env_clip.^2)/sum(env_clip.^2); % 希尔伯特变换+重心拟合算法计算的信号中心
% 
% Color = GetColor(2,1); %用于绘制图形的颜色
% figure();
% plot(z_scan,signal,'LineWidth',1.5,'Color',Color(1,:)); %信号本身
% hold on;
% plot(z_scan,env_signal,'--','LineWidth',1.5,'Color',Color(2,:)); %信号包络
% scatter(z_scan(I_maxindex),env_signal(I_maxindex),25,'filled');  %信号包络最大值对应的索引
% plot(z_clip,env_clip,'--','LineWidth',1.5); %用于重心法计算的信号
% hold off;
% defaultAxes(2);
% xlabel('z/$\mu m$','Interpreter','latex');
% 
% tau = 1;%超高斯窗标准差
% order = 3;%超高斯窗阶数
% window = exp(-((z_scan'-z_center).^2/(2*tau^2)).^order);
% 
% figure();
% tiledlayout(2,1);
% nexttile;
% plot(z_scan,signal,'LineWidth',1.5,'Color',Color(1,:));
% hold on;
% plot(z_scan,window,'--','LineWidth',1.5,'Color',Color(2,:));
% hold off;
% defaultAxes(2);
% xlabel('z/$\mu m$','Interpreter','latex');
% 
% signal = signal.*window;   %信号添加超高斯窗
% 
% nexttile;
% plot(z_scan,signal,'LineWidth',1.5,'Color',GetColor(1,1));
% defaultAxes(2);
% xlabel('z/$\mu m$','Interpreter','latex');
% 
% 
% 
% spectrm = fftshift(fft(ifftshift(signal))); %傅里叶变换
% spectrm = spectrm/max(abs(spectrm));
% fre = make_axis_freq(length(signal),z_peri,'101'); %生成频域坐标
% 
% figure();
% tiledlayout(2,1);
% 
% nexttile;
% plot(fre,abs(spectrm),'LineWidth',1.5,'Color',GetColor(1,1));
% defaultAxes(2);
% xlabel('fre/Hz','Interpreter','latex');
% 
% nexttile;
% plot(fre,unwrap(angle(spectrm)),'LineWidth',1.5,'Color',GetColor(1,1));
% defaultAxes(2);
% xlabel('fre/Hz','Interpreter','latex');
% 
% vk0_choose_num = 9;%波数域中，选择中心波数左右个9个点，共计19个点进行拟合
% vk0_fre = (fre/2)';
% c = sqrt(1-NA^2);
% k0_NA = k0*2*(1+c+c^2)/(3*(1+c));
% [~,index] = min(abs(vk0_fre-k0_NA));%找到最接近中心的1.8的数据对应的索引
% k1 = vk0_fre(index);%用最接近1.8数据的波数作为中心
% vk0_fre_fit = vk0_fre(index-vk0_choose_num:index+vk0_choose_num);
% spectrm_fit = spectrm(index-vk0_choose_num:index+vk0_choose_num); %对频谱进行选取
% % w = abs(spectrm_fit).^2;
% w = abs(spectrm_fit);
% w = w / max(w);
% %理想中，模的平方作为权重能让估计的方差最小，但在实际中，受到非均匀采样、样品反射率等误差的影响，以模作为权重更稳健
% 
% angle_fit = unwrap(angle(spectrm_fit)); %提取相位
% 
% figure();
% tiledlayout(2,1);
% 
% nexttile;
% plot(vk0_fre_fit,abs(spectrm_fit),'LineWidth',1.5,'Color',GetColor(1,1));
% defaultAxes(2);
% xlabel('k/$\mu^{-1}$m','Interpreter','latex');
% 
% nexttile;
% plot(vk0_fre_fit,angle_fit,'LineWidth',1.5,'Color',GetColor(1,1));
% defaultAxes(2);
% xlabel('k/$\mu^{-1}$m','Interpreter','latex');
% %% 最小二乘拟合求系数a,b (使用高效矩阵左除 \ )
% x = 4*pi*(vk0_fre_fit-k1); % 拟合的x项
% 
% % 1. 构造设计矩阵 X (第一列是常数项1对应截距，第二列是x对应斜率)
% X_mat = [ones(length(x), 1), x];
% 
% % 2. 提取权重的平方根
% W_sqrt = sqrt(w);
% 
% % 3. 将权重应用到 X_mat 和 angle_fit (利用MATLAB隐式扩展)
% X_w = X_mat .* W_sqrt; 
% Y_w = angle_fit .* W_sqrt;
% 
% % 4. 使用左除求解 (底层采用极其稳健的 QR 分解)
% coeffs = X_w \ Y_w;
% 
% % 5. 提取系数
% phi = coeffs(1); % 截距 (载波相位)
% h1  = -coeffs(2); % 负斜率 (粗略的高度)
% 
% % 还原拟合曲线用于绘图
% angle_fitlm = phi - h1 * x; % 等价于 X_mat * coeffs
% 
% figure();
% plot(vk0_fre_fit,angle_fit,'--*','LineWidth',1.5,'Color',Color(1,:));
% hold on;
% plot(vk0_fre_fit,angle_fitlm,'LineWidth',1.5,'Color',Color(2,:));
% hold off;
% defaultAxes(2);
% legend('Original Data','Linear fit','EdgeColor','none');
% xlabel('k/$\mu m^{-1}$','Interpreter','latex');
% 
% 
% M = round(-2*k1*h1-phi/(2*pi));%整数M
% h2 = -phi/(4*pi*k1)-M/(2*k1);%精确的高度

[h2,h1,z_center]=CSIPointFDA(signal,z_scan,z_peri,k0,NA);
%%
disp(['预设高度h:',num2str(sample_dis)]);
disp(['希尔伯特变换+重心法:',num2str(z_center)]);
disp(['线性拟合高度h:',num2str(h1)]);
disp(['FDA处理高度h:',num2str(h2)]);
%% 封装FDA算法
function [FDA,LinearFit,Hilb]=CSIPointFDA(signal,z_scan,z_peri,k0,NA)
% 这个函数用FDA算法求解点式CSI信号

% 输入 signal 为CSI单像素轴向序列，double型1×N维向量，N为帧数
% 输入 z_scan 为位置扫描序列，double型1×N维向量，N为帧数
% 输入 z_peri 为z_scan的间隔，等于PZT扫描步长，double型 1×1维向量
% 输入 k0 为光源光谱的质心, double型 1×1维向量
% 输入 NA 为光学系统的NA， double型1×1维向量

% 输出 FDA 为使用FDA算法计算的样品与参考镜之间的距离
% 输出 LinearFit 为使用相位拟合方法计算的样品与参考镜之间的距离
% 输出 Hilb为使用希尔伯特变换+重心法计算的样品与参考镜之间的距离

%     例：
%     z_peri=0.0719;%PZT单次移动间隔，um制
%     N_half = round(10 / z_peri); % 计算单侧点数
%     z_scan = (-N_half : N_half) * z_peri; % 这样生成的数组中心绝对是 0
%     NA=0.3;
%     k0=1.8;
%     [FDA,LinearFit,Hilb]=CSIPointFDA(signal,z_scan,z_peri,k0,NA)

signal=signal(:);
signal=signal-mean(signal);

clip_num = 15; %两侧各取10各点，注意和信号包络宽度有关
env_signal = abs(hilbert(signal)); % 希尔伯特变换提取信号包络
[~, I_maxindex] = max(env_signal); % 信号包络最大值对应的索引

env_clip = env_signal(I_maxindex-clip_num:I_maxindex+clip_num); %选取部分包络，使用重心法计算信号中心（样品与参考镜之间的距离）
z_clip = (z_scan(I_maxindex-clip_num:I_maxindex+clip_num))'; %选取用于计算的信号
z_center = sum(z_clip.*env_clip.^2)/sum(env_clip.^2); % 希尔伯特变换+重心拟合算法计算的信号中心

tau = 1;%超高斯窗标准差
order = 3;%超高斯窗阶数
window = exp(-((z_scan'-z_center).^2/(2*tau^2)).^order);
signal = signal.*window;   %信号添加超高斯窗

spectrm = fftshift(fft(ifftshift(signal))); %傅里叶变换
spectrm = spectrm/max(abs(spectrm));
fre = make_axis_freq(length(signal),z_peri,'101'); %生成频域坐标

vk0_choose_num = 9;%波数域中，选择中心波数左右个9个点，共计19个点进行拟合
vk0_fre = (fre/2)';
c = sqrt(1-NA^2);
k0_NA = k0*2*(1+c+c^2)/(3*(1+c));
[~,index] = min(abs(vk0_fre-k0_NA));%找到最接近中心的1.8的数据对应的索引
k1 = vk0_fre(index);%用最接近1.8数据的波数作为中心
vk0_fre_fit = vk0_fre(index-vk0_choose_num:index+vk0_choose_num);
spectrm_fit = spectrm(index-vk0_choose_num:index+vk0_choose_num); %对频谱进行选取
% w = abs(spectrm_fit).^2;
w = abs(spectrm_fit);
w = w / max(w);
%理想中，模的平方作为权重能让估计的方差最小，但在实际中，受到非均匀采样、样品反射率等误差的影响，以模作为权重更稳健
angle_fit = unwrap(angle(spectrm_fit)); %提取相位

x = 4*pi*(vk0_fre_fit-k1); % 拟合的x项

% 1. 构造设计矩阵 X (第一列是常数项1对应截距，第二列是x对应斜率)
X_mat = [ones(length(x), 1), x];
% 2. 提取权重的平方根
W_sqrt = sqrt(w);
% 3. 将权重应用到 X_mat 和 angle_fit (利用MATLAB隐式扩展)
X_w = X_mat .* W_sqrt; 
Y_w = angle_fit .* W_sqrt;
% 4. 使用左除求解 (底层采用极其稳健的 QR 分解)
coeffs = X_w \ Y_w;
% 5. 提取系数
phi = coeffs(1); % 截距 (载波相位)
h1  = -coeffs(2); % 负斜率 (粗略的高度)
M = round(-2*k1*h1-phi/(2*pi));%整数M
h2 = -phi/(4*pi*k1)-M/(2*k1);%精确的高度

FDA=h2;
LinearFit=h1;
Hilb=z_center;
end