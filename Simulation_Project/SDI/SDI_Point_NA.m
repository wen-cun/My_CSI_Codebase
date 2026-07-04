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

C=angle((r_Se.*conj(r_Me)+r_Sm.*conj(r_Mm))); %判断干涉项的相位是否随波长、角度变化明显
[theta_array_surf,vk0_surf]=meshgrid(theta_array,vk0);
figure();

surf(theta_array_surf,vk0_surf,C);
shading interp;
defaultAxes(3);
defaultColor(1);
ylabel('k/$\mu m^{-1}$','Interpreter','latex');
xlabel('$\theta$/rad','Interpreter','latex');
%% 选择偏振模式，生成光谱干涉信号
system_pol = 'unpolar';%非偏振模式
sample_dis = 0.432495;
signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,sample_dis,system_pol);

figure();
plot(lam,signal,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('$\lambda$/$\mu m$','Interpreter','latex');


%信号加高斯白噪声
SNR = 40; %40dB的噪声
signal = awgn(signal,SNR,'measured');
signal = signal/max(abs(signal)); %归一化

figure();
plot(lam,signal,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
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
%% 精确定位,粗网格修正粗定位结果
tic;
z_add = 2.5; %扫描的上范围
z_min = max(z_coa-z_add,0); %搜索的下限
z_max = z_coa+z_add; %搜索的上限
z_peri = 2.5e-2; %50nm采样
z_gra = z_min:z_peri:z_max; %搜索的区间
cost = nan*ones(size(z_gra)); %预先分配内存
for ii=1:length(z_gra)
    signal_gra = SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,z_gra(ii),system_pol); %计算信号
    cost(ii) = sum(abs(signal(valid)-signal_gra(valid)).^2,'all');%计算误差
end
[~,index]=min(cost);

z_coa_pos=z_gra(index); %粗网格修正粗定位结果
figure();
plot(z_gra,cost,'LineWidth',1.5,'Color',Color(1,:));
hold on;
scatter(z_gra(index),cost(index),25,Color(2,:),'filled');
hold off;
defaultAxes(2);
xlabel('z/$\mu$ m','Interpreter','latex');

z_minus_min = -z_coa-z_add; %搜索负半部分
z_minus_max = min(-z_coa+z_add,0); %搜索负半部分
z_gra = z_minus_min:z_peri:z_minus_max; %搜索的区间
cost_minus = nan*ones(size(z_gra)); %预先分配内存
for ii=1:length(z_gra)
    signal_gra = SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,z_gra(ii),system_pol); %计算信号
    cost_minus(ii) = sum(abs(signal(valid)-signal_gra(valid)).^2,'all');%计算误差
end
[~,index]=min(cost_minus);

z_coa_minus=z_gra(index); %粗网格修正粗定位结果

hold on;
plot(z_gra,cost_minus,'LineWidth',1.5,'Color',Color(1,:));
scatter(z_gra(index),cost_minus(index),25,Color(2,:),'filled');
hold off;
toc;
disp(['粗网格拟合法/正半轴:',num2str(z_coa_pos)]);
disp(['粗网格拟合法/负半轴:',num2str(z_coa_minus)]);

%% 精确定位，细网格搜索正负半轴分别进行细网格搜索
z_add_fine = 0.3;
z_step_fine = 1e-3;

% 正半轴细网格
z_pos_min = max(z_coa_pos-z_add_fine,0);
z_pos_max = z_coa_pos+z_add_fine;
z_fine_pos = z_pos_min:z_step_fine:z_pos_max;

% 负半轴细网格
z_neg_min = z_coa_minus-z_add_fine;
z_neg_max = min(z_coa_minus+z_add_fine,0);
z_fine_neg = z_neg_min:z_step_fine:z_neg_max;

cost_fine_pos = nan(size(z_fine_pos));
cost_fine_neg = nan(size(z_fine_neg));

tic;

% 正半轴细搜索
for ii = 1:numel(z_fine_pos)

    signal_model = SDIPointSignalGenerate( ...
        NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm, ...
        theta_array,z_fine_pos(ii),system_pol);

    cost_fine_pos(ii) = sum( ...
        abs(signal(valid)-signal_model(valid)).^2);
end

% 负半轴细搜索
for ii = 1:numel(z_fine_neg)

    signal_model = SDIPointSignalGenerate( ...
        NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm, ...
        theta_array,z_fine_neg(ii),system_pol);

    cost_fine_neg(ii) = sum( ...
        abs(signal(valid)-signal_model(valid)).^2);
end

toc;
figure();
plot(z_fine_pos,cost_fine_pos,'LineWidth',1.5,'Color',Color(1,:));
hold on;
plot(z_fine_neg,cost_fine_neg,'LineWidth',1.5,'Color',Color(1,:));
hold off;
defaultAxes(2);
xlabel('z/$\mu$ m','Interpreter','latex');

[z_pos_refined,cost_pos_refined,pos_valid] = ...
    RefineMinimumByParabola(z_fine_pos,cost_fine_pos);

[z_neg_refined,cost_neg_refined,neg_valid] = ...
    RefineMinimumByParabola(z_fine_neg,cost_fine_neg);

if cost_pos_refined <= cost_neg_refined
    z_pre = z_pos_refined;
    cost_pre = cost_pos_refined;
    selected_branch = "positive";
else
    z_pre = z_neg_refined;
    cost_pre = cost_neg_refined;
    selected_branch = "negative";
end

disp(['选择分支: ',char(selected_branch)]);
hold on;
scatter(z_pre,cost_pre,25,Color(2,:),'filled');
hold off;
%% 输出结果
disp(['模型拟合法:',num2str(z_pre)]);
%% 抛物线拟合的函数
function [z_refined,cost_refined,success] = ...
    RefineMinimumByParabola(z_grid,cost_grid)

z_grid = z_grid(:);
cost_grid = cost_grid(:);

[cost_min,index_min] = min(cost_grid);

z_refined = z_grid(index_min);
cost_refined = cost_min;
success = false;

% 最小值位于边缘时无法做三点拟合
if index_min == 1 || index_min == numel(z_grid)
    warning('局部最小值位于细搜索边界，使用离散最小值。');
    return;
end

idx = index_min-1:index_min+1;

% 以中心点为原点，提高数值稳定性
z0 = z_grid(index_min);
z_local = z_grid(idx)-z0;
cost_local = cost_grid(idx);

% cost = a*z^2+b*z+c
p = polyfit(z_local,cost_local,2);

% 必须开口向上
if ~all(isfinite(p)) || p(1) <= 0
    warning('局部代价函数不满足开口向上的抛物线条件。');
    return;
end

z_offset = -p(2)/(2*p(1));

% 顶点必须位于左右两个相邻采样点之间
if z_offset < z_local(1) || z_offset > z_local(3)
    warning('抛物线顶点超出局部三点范围。');
    return;
end

z_refined = z0+z_offset;
cost_refined = polyval(p,z_offset);

% SSE理论上不应小于0
cost_refined = max(cost_refined,0);

success = true;

end