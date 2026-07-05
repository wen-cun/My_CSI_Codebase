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
sample_dis = 9.0004;
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
z_peri = 1e-2; %10nm采样

z_minus_min = -z_coa-z_add; %搜索负半部分
z_minus_max = min(-z_coa+z_add,0); %搜索负半部分

z_gra = [z_minus_min:z_peri:z_minus_max,z_min:z_peri:z_max]; %搜索的区间
cost = nan*ones(size(z_gra)); %预先分配内存
for ii=1:length(z_gra)
    signal_gra = SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,z_gra(ii),system_pol); %计算信号
    cost(ii) = sum(abs(signal(valid)-signal_gra(valid)).^2,'all');%计算误差
end
toc;

local_mask = islocalmin(cost,'FlatSelection','center','SamplePoints',z_gra); %判断是否为局部极小值
idx_local = find(local_mask);

% 防止严格单调或边缘最小导致没有检测到局部极小值
[~,idx_global] = min(cost);

if isempty(idx_local)
    idx_local = idx_global;
elseif ~ismember(idx_global,idx_local)
    idx_local = [idx_local;idx_global];
end

figure();
plot(z_gra,cost,'LineWidth',1.5,'Color',Color(1,:));
hold on;
scatter(z_gra(local_mask),cost(local_mask),25,Color(2,:),'filled');
hold off;
defaultAxes(2);
legend('Cost','LocalMin','EdgeColor','none');
xlabel('z/$\mu$ m','Interpreter','latex');

% 按损失值排序
[~,order] = sort(cost(idx_local),'ascend');
idx_sorted = idx_local(order);




%% 精确定位，细网格搜索正负半轴分别进行细网格搜索

z_num = 6; %仅选取前六个极小值

disp('粗网格候选位置：');
for ii=1:z_num
    disp(num2str(z_gra(idx_sorted(ii))));
end

z_add_fine = 1.5*z_peri; %精确定位范围采用粗网格定义的1个半以内。
z_step_fine = 1e-3;
num_fine_half = ceil(z_add_fine/z_step_fine);
fine_offsets = (-num_fine_half:num_fine_half)*z_step_fine;


cost_cand = nan*ones(1,numel(fine_offsets)) ; %给每个候选位置损失函数分配内存

z_para = nan*ones(1,z_num) ; %给每个候选位置抛物线位置分配内存
cost_para = nan*ones(1,z_num) ; %给每个候选位置抛物线拟合误差分配内存

tic;
figure();
for ii = 1:z_num
    z_cand = fine_offsets + z_gra(idx_sorted(ii));
    for jj = 1:numel(z_cand)
        
        signal_model = SDIPointSignalGenerate( ...
            NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm, ...
            theta_array,z_cand(jj),system_pol);
        cost_cand(jj) = sum( ...
            abs(signal(valid)-signal_model(valid)).^2);
    end
    plot(z_cand,cost_cand,'LineWidth',1.5,'Color',Color(1,:));
    hold on;
    [z_para(ii),cost_para(ii),~] = RefineMinimumByParabola(z_cand,cost_cand);
end

toc;

[cost_pre,min_index] = min(cost_para); %确定最小抛物线候选中心
z_pre = z_para(min_index); %精确搜索的结果



hold on;
scatter(z_pre,cost_pre,25,Color(2,:),'filled');
hold off;
defaultAxes(2);
xlabel('z/$\mu$ m','Interpreter','latex');


%% 调用封装函数
z_coa_fun = SDIPointModulFit(signal,lam,vsk_ini,valid,NA);
z_pre_fun = SDIPointModelFit(signal,z_coa_fun,valid,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol);
%% 输出结果
disp(['模型拟合法:',num2str(z_pre)]);

disp(['封装函数粗定位:',num2str(z_coa_fun)]);
disp(['封装函数模型拟合:',num2str(z_pre_fun)]);
%%
function z_coa=SDIPointModulFit(signal,lam,vsk_ini,valid,NA)
% 这个函数用调制拟合的方法给出位移值的粗定位结果

% 输入 signal 为等波长空间的信号强度，double型N×1维向量
% 输入 lam 为等波长空间的采样点波长，double型N×1维向量
% 输入 vsk_ini 为等波长空间的光谱强度，double型N×1维向量
% 输入 valid 为计算信号内部的有效位置逻辑数组，Logical型N×1维向量
% 输入 NA 为光学系统的NA，double型1×1维向量

% 输出 z_coa 为位移值的粗定位结果，double型1×1维向量

%     例：
%     lam_peri = 3e-4; %波长采样周期,um制
%     lam_lim = [0.3,1.1]; %波长采样范围
%     lam = lam_lim(1):lam_peri:lam_lim(2); %波长采样
%     lam = lam(:); %转为列向量
%     vk0 = 1./lam; %转为波数
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     NA=0;
%     sample_dis=5;
%     system_pol='ideal';
%     signal=SDIPointSignalGenerate(NA,vk0,nan,nan,nan,nan,0,sample_dis,system_pol);
%     source_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行粗略定位，以降低噪声
%     valid = vsk_ini > source_thr*max(vsk_ini);
%     z_coa=SDIPointModulFit(signal,lam,vsk_ini,valid,NA);

source_thr = 0.01;  %结算调制相位的光谱强度阈值
valid_div = vsk_ini > source_thr*max(vsk_ini); %光谱阈值上的数据对应的索引

interD = nan(size(signal)); %调制项分配内存
interD(valid_div) = signal(valid_div)./vsk_ini(valid_div); %计算调制项


x_cho = 4*pi*1./lam(valid);
y_cho = interD(valid);

weight = vsk_ini(valid).^2; %将光谱数据强度当做权重
weight = weight./ max(weight); %权重归一化


fitResult = fit(x_cho,y_cho,'fourier1','Weights',weight);
coeff = coeffvalues(fitResult);

c = sqrt(1-NA^2);
mu = 2*(1+c+c^2)/(3*(1+c)); %小NA修正系数

z_coa = abs(coeff(4))/mu; %粗定位结果

end
%% 
function cost=CalcSDISeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal,valid)
% 这个函数用来计算搜索范围内模型信号与目标信号的损失值

% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入system_pol为系统的偏振状态
% 输入z_gra为样品与参考镜的距离, double型1×K维向量，K为搜索的点数
% 输入 signal 为要匹配的目标信号，double型 N×1维向量
% 输入 valid 为计算信号内部的有效位置逻辑数组，Logical型N×1维向量

% 输出 cost为每个搜索位置，与目标信号的损失误差，double型1×K维向量

%     例：
%     lam_peri = 3e-4; %波长采样周期,um制
%     lam_lim = [0.3,1.1]; %波长采样范围
%     lam = lam_lim(1):lam_peri:lam_lim(2); %波长采样
%     lam = lam(:); %转为列向量
%     vk0 = 1./lam; %转为波数
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vsk=vsk./(max(abs(vsk)));
%     NA = 0.3; %系统NA
%     theta_max=asin(NA); %最大NA对应的空气中光线角度theta
%     theta_peri=0.01; %角度theta的采样周期
%     theta_array = 0:theta_peri:theta_max; %theta采样数组
%     sample_stru={'SiO2',5;...
%                   'Si',inf};
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru);
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array);
%     system_pol='ideal';
%     signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,3,system_pol));
%     z_gra=0:1e-3:5;
%     source_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行粗略定位，以降低噪声
%     valid = vsk_ini > source_thr*max(vsk_ini);
%     cost=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal,valid);


cost = nan*ones(size(z_gra)); %预先分配内存
for ii=1:length(z_gra)
    signal_gra = SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,z_gra(ii),system_pol); %计算信号
    cost(ii) = sum(abs(signal(valid)-signal_gra(valid)).^2,'all');%计算误差
end
end
%%
function [z_pre,cost_pre]=SDIPointModelFit(signal,z_coa,valid,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol)
% 这个函数用精细网格给出SDI信号位移值的精确定位结果

% 输入 signal 为要匹配的目标信号，double型 N×1维向量
% 输入 z_coa 为粗定位的结果，double型1×1维向量
% 输入 valid 为计算信号内部的有效位置逻辑数组，Logical型N×1维向量
% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入system_pol为系统的偏振状态

% 输出 z_pre 为精确定位的结果，double型1×1维向量
% 输出 cost_pre 为精确定位结果对应的损失值,double 型1×1维向量

%     例：
%     lam_peri = 3e-4; %波长采样周期,um制
%     lam_lim = [0.3,1.1]; %波长采样范围
%     lam = lam_lim(1):lam_peri:lam_lim(2); %波长采样
%     lam = lam(:); %转为列向量
%     vk0 = 1./lam; %转为波数
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vsk=vsk./(max(abs(vsk)));
%     NA = 0.3; %系统NA
%     theta_max=asin(NA); %最大NA对应的空气中光线角度theta
%     theta_peri=0.01; %角度theta的采样周期
%     theta_array = 0:theta_peri:theta_max; %theta采样数组
%     sample_stru={'SiO2',5;...
%                   'Si',inf};
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru);
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array);
%     system_pol='ideal';
%     signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,3,system_pol));
%     source_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行粗略定位，以降低噪声
%     valid = vsk_ini > source_thr*max(vsk_ini);
%     z_coa=SDIPointModulFit(signal,lam,vsk_ini,valid,NA);
%     [z_pre,cost_pre]=SDIPointModelFit(signal,z_coa,valid,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol)

  
z_add = 2.5; %扫描的上范围
z_min = max(z_coa-z_add,0); %搜索的下限
z_max = z_coa+z_add; %搜索的上限
z_peri = 1e-2; %10nm采样

z_minus_min = -z_coa-z_add; %搜索负半部分
z_minus_max = min(-z_coa+z_add,0); %搜索负半部分

z_gra = [z_minus_min:z_peri:z_minus_max,z_min:z_peri:z_max]; %搜索的区间

cost = CalcSDISeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal,valid);

local_mask = islocalmin(cost,'FlatSelection','center','SamplePoints',z_gra); %判断是否为局部极小值
idx_local = find(local_mask);

% 防止严格单调或边缘最小导致没有检测到局部极小值
[~,idx_global] = min(cost);

if isempty(idx_local)
    idx_local = idx_global;
elseif ~ismember(idx_global,idx_local)
    idx_local = [idx_local;idx_global];
end
% 按损失值排序
[~,order] = sort(cost(idx_local),'ascend');
idx_sorted = idx_local(order);
    
z_num = 6; %仅选取前六个极小值

z_add_fine = 1.5*z_peri; %精确定位范围采用粗网格定义的1个半以内。
z_step_fine = 1e-3;
num_fine_half = ceil(z_add_fine/z_step_fine);
fine_offsets = (-num_fine_half:num_fine_half)*z_step_fine;



z_para = nan*ones(1,z_num) ; %给每个候选位置抛物线位置分配内存
cost_para = nan*ones(1,z_num) ; %给每个候选位置抛物线拟合误差分配内存

for ii = 1:z_num
    z_cand = fine_offsets + z_gra(idx_sorted(ii));
    cost_cand = CalcSDISeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_cand,signal,valid);
    [z_para(ii),cost_para(ii),~] = RefineMinimumByParabola(z_cand,cost_cand);
end


[cost_pre,min_index] = min(cost_para); %确定最小抛物线候选中心
z_pre = z_para(min_index); %精确搜索的结果    
    
end
%% 抛物线拟合的函数
function [z_refined,cost_refined,success] =RefineMinimumByParabola(z_grid,cost_grid)
%这个函数用来给损失函数最小值处三点抛物线拟合

% 输入 z_grid 是要进行抛物线拟合的位移搜索范围，double型1×N维向量
% 输入 cost_grid 是要进行抛物线拟合的损失值，double型1×N维向量

% 输出 z_refined 是抛物线拟合后的中心，double型1×1维向量
% 输出 cost_refined 是抛物线中心对应的损失函数，double型1×1维向量
% 输出 success 是成功的标志，double型1×1维Logical 向量，返回1时表明成功进行抛物线拟合，返回0时表明失效

%     例：
%     z_grid=[0:1e-3:2.5];
%     lam_peri = 3e-4; %波长采样周期,um制
%     lam_lim = [0.3,1.1]; %波长采样范围
%     lam = lam_lim(1):lam_peri:lam_lim(2); %波长采样
%     lam = lam(:); %转为列向量
%     vk0 = 1./lam; %转为波数
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vsk=vsk./(max(abs(vsk)));
%     NA = 0.3; %系统NA
%     theta_max=asin(NA); %最大NA对应的空气中光线角度theta
%     theta_peri=0.01; %角度theta的采样周期
%     theta_array = 0:theta_peri:theta_max; %theta采样数组
%     sample_stru={'SiO2',5;...
%                   'Si',inf};
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru);
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array);
%     system_pol='ideal';
%     signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,3,system_pol));
%     z_grid=0:1e-3:5;
%     source_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行粗略定位，以降低噪声
%     valid = vsk_ini > source_thr*max(vsk_ini);
%     cost_grid=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_grid,signal,valid);
%     [z_refined,cost_refined,success]=RefineMinimumByParabola(z_grid,cost_grid);

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