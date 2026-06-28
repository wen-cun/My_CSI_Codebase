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
sample_dis=5;
signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,sample_dis,system_pol);

figure();
plot(lam,signal,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('$\lambda$/$\mu m$','Interpreter','latex');

%信号加高斯白噪声
SNR = 40; %40dB的噪声
signal = awgn(signal,SNR,'measured');
signal = signal/max(abs(signal)); %归一化

% figure();
% plot(lam,signal,'LineWidth',1.5,'Color',GetColor(1,1));
% defaultAxes(2);
% xlabel('$\lambda$/$\mu m$','Interpreter','latex'); 
Color=GetColor(2,1);
figure();
plot(lam,signal,'LineWidth',1.5,'Color',Color(1,:));
hold on;
plot(lam,vsk_ini,'--','LineWidth',1.5,'Color',Color(2,:));
hold off;
defaultAxes(2);
legend('Singal','Source','EdgeColor','none');
xlabel('$\lambda$/$\mu m$','Interpreter','latex');

%% 求解调制项
% interD = (signal+eps)./(vsk_ini+eps);
source_thr = 0.01;
valid_div = vsk_ini > source_thr*max(vsk_ini);

interD = nan(size(signal));
interD(valid_div) = signal(valid_div)./vsk_ini(valid_div);
figure();
plot(lam,interD,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('$\lambda$/$\mu m$','Interpreter','latex'); 

%% 粗略定位
source_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行粗略定位，以降低噪声
valid = vsk_ini > source_thr*max(vsk_ini);

figure();
plot(lam(valid),interD(valid),'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('$\lambda$/$\mu m$','Interpreter','latex'); 

x_cho = 4*pi*1./lam(valid);
y_cho = interD(valid);
weight = vsk_ini(valid).^2; %权重
weight = weight./ max(weight); %权重归一化

tic;
[fitResult, gof] = fit(x_cho,y_cho,'fourier1','Weights',weight);
toc;
coeff = coeffvalues(fitResult);
y_fit = coeff(1)+coeff(2)*cos(coeff(4)*x_cho)+coeff(3)*sin(coeff(4)*x_cho);

c = sqrt(1-NA^2);
mu = 2*(1+c+c^2)/(3*(1+c));

z_coa = abs(coeff(4))/mu;

Color = GetColor(2,1);
figure();
scatter(x_cho,y_cho,20,Color(1,:),'filled');
defaultAxes(2);
hold on;
plot(x_cho,y_fit,'LineWidth',1.5,'Color',Color(2,:));
hold off;
legend('Original Data','Fourier fit');
title(['Rsquare=',num2str(gof.rsquare)]);
xlabel('4$\pi/\lambda$','Interpreter','latex');

disp(['预设高度h:',num2str(sample_dis)]);
disp(['粗定位法:',num2str(z_coa)]);
%% 精确定位
tic;
if gof.rsquare<=0.75
    warning('粗定位拟合效果较差，已改为粗细网格搜索法');
    z_add = 2.5; %扫描的上范围
    z_min = max(z_coa-z_add,0); %搜索的下限
    z_max = z_coa+z_add; %搜索的上限
    z_peri = 5e-2; %50nm采样
    z_gra = z_min:z_peri:z_max; %搜索的区间
    cost = nan*ones(size(z_gra)); %预先分配内存
    for ii=1:length(z_gra)
        signal_gra = SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,z_gra(ii),system_pol); %计算信号
        cost(ii) = sum(abs(signal(valid)-signal_gra(valid)).^2,'all');%计算误差
    end
    [~,index]=min(cost);
    figure();
    plot(z_gra,cost,'LineWidth',1.5,'Color',Color(1,:));
    hold on;
    scatter(z_gra(index),cost(index),25,Color(2,:),'filled');
    hold off;
    defaultAxes(2);
    xlabel('z/$\mu$ m','Interpreter','latex');
    z_coa=z_gra(index); %粗网格修正粗定位结果
    disp(['粗网格拟合法:',num2str(z_coa)]);
end
z_add = 0.5; %扫描的上下范围
z_min = max(z_coa-z_add,0); %搜索的下限
z_max = z_coa+z_add; %搜索的上限
z_peri = 1e-3; %1nm采样
z_gra = z_min:z_peri:z_max; %搜索的区间
cost = nan*ones(size(z_gra)); %预先分配内存
for ii=1:length(z_gra)
    signal_gra = SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,z_gra(ii),system_pol); %计算信号
    cost(ii) = sum(abs(signal(valid)-signal_gra(valid)).^2,'all');%计算误差
end
toc;
figure();
plot(z_gra,cost,'LineWidth',1.5,'Color',Color(1,:));
defaultAxes(2);
xlabel('z/$\mu$ m','Interpreter','latex');

[~,min_index] = min(cost); %误差最小值对应的索引
z_pre = z_gra(min_index);
cost_pre = cost(min_index);

if min_index == 1 || min_index == length(z_gra)
    warning('误差最小值位于搜索边界，建议扩大 搜索范围');
else
    %抛物线三点拟合
    % 最小值左右相邻的三个代价函数值
    cost_left   = cost(min_index-1);
    cost_center = cost(min_index);
    cost_right  = cost(min_index+1);
    
    denominator = cost_left - 2*cost_center + cost_right;
    
    % denominator>0 表示局部抛物线开口向上
    if isfinite(denominator) && denominator > eps(max(abs(cost)))
        delta_index = 0.5 * ...
            (cost_left - cost_right) / denominator;
        
        % 顶点原则上应该落在相邻两个采样点之间
        if abs(delta_index) <= 1
            z_pre = z_gra(min_index) + delta_index*z_peri;
            
            % 抛物线顶点对应的代价值
            cost_pre = cost_center ...
                - (cost_left-cost_right)^2/(8*denominator);
        else
            warning('抛物线顶点超出相邻采样区间，使用离散最小值。');
        end
    else
        warning('局部代价函数不满足开口向上的抛物线条件，使用离散最小值。');
    end
end

hold on;
scatter(z_pre,cost_pre,25,Color(2,:),'filled');
hold off;
%% 输出结果
disp(['模型拟合法:',num2str(z_pre)]);