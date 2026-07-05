%% 这个脚本用来恢复指定序列的SDI信号
clear;
close all;
clc;
%% 生成扫描序列
z_peri=0.0719;%PZT单次移动间隔，um制
% z_samp=0.05;%PZT单次移动间隔，um制

N_half = round(10 / z_peri); % 计算单侧点数
z_scan_pur = (-N_half : N_half) * z_peri; %不带振动的纯净扫描序列

t_peri = 1e-3; %时间采样
t = (0:length(z_scan_pur)-1)*t_peri;

vib_freq_hz = [100,200];
vib_amp = [0.2,0.1];

vib = zeros(size(t));
for ii = 1:numel(vib_freq_hz)
    vib = vib + vib_amp(ii).*cos(2*pi*vib_freq_hz(ii).*t);
end

z_scan = z_scan_pur + vib; %扫描序列加上振动量

figure();
plot(t,vib,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('x/$\mu$m','Interpreter','latex');

Color=GetColor(2,1);
figure();
plot(z_scan_pur,z_scan_pur,'LineWidth',1.5,'Color',Color(1,:));
hold on;
plot(z_scan_pur,z_scan,'-','LineWidth',1.5,'Color',Color(2,:));
hold off;
legend('Pur','Vib','EdgeColor','none','Location','NorthWest');
defaultAxes(2);
xlabel('x/$\mu$m','Interpreter','latex');

vib_spe = fftshift(fft(vib)); %振动的频谱
vib_spe = vib_spe./max(abs(vib_spe));
vib_fre=make_axis_freq(length(t),t_peri,'101');
vib_pow = abs(vib_spe).^2; %振动的功率谱
vib_val = vib_fre >= 0;

figure();
plot(vib_fre(vib_val),vib_pow(vib_val),'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('x/Hz');
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
NA = 0.4; %系统NA 
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
%% 选择偏振模式，生成光谱干涉信号
system_pol = 'unpolar';%非偏振模式

tic;
signal=nan*ones(length(vk0),length(z_scan));%信号序列，每一列代表指定离焦量位移下的SDI信号
for ii=1:length(z_scan)
    signal(:,ii) = SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,z_scan(ii),system_pol); %生成信号
    %信号加高斯白噪声
    SNR = 40; %40dB的噪声
    signal(:,ii) = awgn(signal(:,ii),SNR,'measured');
    signal(:,ii) = signal(:,ii)/max(abs(signal(:,ii))); %归一化
end
toc;

figure();
for ii=1:15:length(z_scan)
plot3(z_scan(ii)*ones(size(vk0)),lam,signal(:,ii),'LineWidth',1.5);
hold on;
end
hold off;
defaultAxes(3);
ylabel('$\lambda$/$\mu$m','Interpreter','latex');
%% 粗细拟合恢复扫描序列
z_scan_rec = nan*ones(size(z_scan)); %恢复的扫描序列，预先分配内存
source_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行粗略定位，以降低噪声
valid = vsk_ini > source_thr*max(vsk_ini);
tic;
for ii=1:length(z_scan_rec)
    if mod(ii,10)==0
        disp(['正在计算',num2str(ii),'/',num2str(length(z_scan_rec))]);
    end
    z_coa = SDIPointModulFit(signal(:,ii),lam,vsk_ini,valid,NA);
    z_scan_rec(ii) = SDIPointModelFit(signal(:,ii),z_coa,valid,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol);
    
end
toc;

vib_rec = z_scan_rec-z_scan_pur; %振动噪声
%% 展示结果

figure();
scatter(z_scan,z_scan_rec,25,GetColor(1,1),'filled');
defaultAxes(2);
xlabel('z/$\mu$m','Interpreter','latex');
ylabel('z$_{rec}$/$\mu$m','Interpreter','latex');

figure();
plot(z_scan,z_scan-z_scan_rec,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu$m','Interpreter','latex');

figure();
plot(t,vib_rec,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu$m','Interpreter','latex');
