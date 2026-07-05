%% 这个脚本用来仿真SDI补偿CSI振动问题
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
%% 选择偏振模式，生成CSI信号
system_pol = 'unpolar';%非偏振模式
% system_pol = 'ideal';%理想偏振模式

sample_dis = 3;%样品与参考镜之间的距离;

signal_CSI = CSIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,theta_array,sample_dis,system_pol);%生成CSI信号

figure();
plot(z_scan_pur,signal_CSI,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');

%信号加高斯白噪声
SNR = 40; %40dB的噪声
signal_CSI = awgn(signal_CSI,SNR,'measured');
signal_CSI = signal_CSI/max(abs(signal_CSI)); %归一化

figure();
plot(z_scan_pur,signal_CSI,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');

%% 生成SDI信号
tic;
signal_SDI=nan*ones(length(vk0),length(z_scan));%信号序列，每一列代表指定离焦量位移下的SDI信号
for ii=1:length(z_scan)
    signal_SDI(:,ii) = SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,z_scan(ii),system_pol); %生成信号
    %信号加高斯白噪声
    SNR = 40; %40dB的噪声
    signal_SDI(:,ii) = awgn(signal_SDI(:,ii),SNR,'measured');
    signal_SDI(:,ii) = signal_SDI(:,ii)/max(abs(signal_SDI(:,ii))); %归一化
end
toc;

figure();
for ii=1:15:length(z_scan)
plot3(z_scan(ii)*ones(size(vk0)),lam,signal_SDI(:,ii),'LineWidth',1.5);
hold on;
end
hold off;
defaultAxes(3);
ylabel('$\lambda$/$\mu$m','Interpreter','latex');
%% 从SDI信号恢复扫描序列
z_scan_rec = nan*ones(size(z_scan)); %恢复的扫描序列，预先分配内存
source_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行粗略定位，以降低噪声
valid = vsk_ini > source_thr*max(vsk_ini);
tic;
for ii=1:length(z_scan_rec)
    if mod(ii,10)==0
        disp(['正在计算',num2str(ii),'/',num2str(length(z_scan_rec))]);
    end
    z_coa = SDIPointModulFit(signal_SDI(:,ii),lam,vsk_ini,valid,NA);
    z_scan_rec(ii) = SDIPointModelFit(signal_SDI(:,ii),z_coa,valid,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol);
    
end
toc;

vib_rec = z_scan_rec-z_scan_pur; %振动噪声

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
%% 从恢复的非均匀采样序列，处理CSI信号
%粗定位 能量重心法


signal_CSI = signal_CSI(:);
signal_CSI = signal_CSI-mean(signal_CSI); %去掉直流

z_scan = z_scan(:);
valid_z = isfinite(z_scan);

energy = signal_CSI.^2;

d_coa = sum(energy(valid_z).*z_scan(valid_z)) ...
      / sum(energy(valid_z));
  
% d_coa = sum((signal.^2.*z_scan_un'),'omitnan')./sum(signal.^2,'omitnan');

z_scan=z_scan';

figure();
scatter(z_scan,signal_CSI,25,Color(1,:),'filled');
hold on;
scatter(d_coa,0,25,Color(2,:),'filled');
hold off;
defaultAxes(2);
legend('Signal','Center','EdgeColor','none');
xlabel('z/$\mu m$','Interpreter','latex');



disp(['预设高度h:',num2str(sample_dis)]);
disp(['粗定位法:',num2str(d_coa)]);

% 精确定位,粗网格修正粗定位结果
tic;
d_add = 0.5; %扫描的上范围
d_min = max(d_coa-d_add,0); %搜索的下限
d_max = d_coa+d_add; %搜索的上限
d_peri = 1e-2; %10nm采样

d_minus_min = -d_coa-d_add; %搜索负半部分
z_minus_max = min(-d_coa+d_add,0); %搜索负半部分

d_gra = [d_minus_min:d_peri:z_minus_max,d_min:d_peri:d_max]; %搜索的区间
cost = nan*ones(size(d_gra)); %预先分配内存
for ii=1:length(d_gra)
    signal_gra = CSIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,theta_array,d_gra(ii),system_pol); %计算信号
    signal_gra = signal_gra(:);
    signal_gra = signal_gra-mean(signal_gra,'omitnan'); %去掉直流
    cost(ii) = sum(abs(signal_CSI-signal_gra).^2,'all','omitnan');%计算误差
end
toc;

local_mask = islocalmin(cost,'FlatSelection','center','SamplePoints',d_gra); %判断是否为局部极小值
idx_local = find(local_mask);

% 防止严格单调或边缘最小导致没有检测到局部极小值
[~,idx_global] = min(cost);

if isempty(idx_local)
    idx_local = idx_global;
elseif ~ismember(idx_global,idx_local)
    idx_local = [idx_local;idx_global];
end

figure();
plot(d_gra,cost,'LineWidth',1.5,'Color',Color(1,:));
hold on;
scatter(d_gra(local_mask),cost(local_mask),25,Color(2,:),'filled');
hold off;
defaultAxes(2);
legend('Cost','LocalMin','EdgeColor','none');
xlabel('z/$\mu$ m','Interpreter','latex');

% 按损失值排序
[~,order] = sort(cost(idx_local),'ascend');
idx_sorted = idx_local(order);




% 精确定位，细网格搜索正负半轴分别进行细网格搜索

d_num = 6; %仅选取前六个极小值

disp('粗网格候选位置：');
for ii=1:d_num
    disp(num2str(d_gra(idx_sorted(ii))));
end

d_add_fine = 1.5*d_peri; %精确定位范围采用粗网格定义的1个半以内。
d_step_fine = 1e-3;
num_fine_half = ceil(d_add_fine/d_step_fine);
fine_offsets = (-num_fine_half:num_fine_half)*d_step_fine;


cost_cand = nan*ones(1,numel(fine_offsets)) ; %给每个候选位置损失函数分配内存

d_para = nan*ones(1,d_num) ; %给每个候选位置抛物线位置分配内存
cost_para = nan*ones(1,d_num) ; %给每个候选位置抛物线拟合误差分配内存

tic;
figure();
for ii = 1:d_num
    d_cand = fine_offsets + d_gra(idx_sorted(ii));
    for jj = 1:numel(d_cand)
     signal_gra = CSIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,theta_array,d_cand(jj),system_pol); %计算信号
    signal_gra = signal_gra(:);
    signal_gra = signal_gra-mean(signal_gra,'omitnan'); %去掉直流
    cost_cand(jj) = sum(abs(signal_CSI-signal_gra).^2,'all','omitnan');%计算误差
        
    end
    plot(d_cand,cost_cand,'LineWidth',1.5,'Color',Color(1,:));
    hold on;
    [d_para(ii),cost_para(ii),~] = RefineMinimumByParabola(d_cand,cost_cand);
end

toc;

[cost_pre,min_index] = min(cost_para); %确定最小抛物线候选中心
d_pre = d_para(min_index); %精确搜索的结果



hold on;
scatter(d_pre,cost_pre,25,Color(2,:),'filled');
hold off;
defaultAxes(2);
xlabel('z/$\mu$ m','Interpreter','latex');

disp(['细网格拟合法:',num2str(d_pre)]);
%% 展示恢复的CSI信号
signal_rec = CSIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan_pur,theta_array,d_pre,system_pol); %计算信号
signal_rec = signal_rec-mean(signal_rec,'omitnan'); %去掉直流

figure();
plot(z_scan,signal_CSI+1,'LineWidth',1.5,'Color',Color(1,:));
hold on;
plot(z_scan_pur,signal_rec,'LineWidth',1.5,'Color',Color(2,:));
hold off;
legend('Original','Recover','EdgeColor','none');
defaultAxes(2);
xlabel('z/$\mu$ m','Interpreter','latex');