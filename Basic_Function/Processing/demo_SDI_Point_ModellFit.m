%% 采用调试拟合法给出SDI的粗定位结果的函数
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

figure();
plot(lam,signal,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('$\lambda$/$\mu m$','Interpreter','latex'); 

%% 粗定位结果
tic;
[z_coa,rsquare]=SDIPointModulFit(signal,lam,vsk_ini,NA);
toc;

disp(['预设高度h:',num2str(sample_dis)]);
disp(['粗定位法:',num2str(z_coa),'置信程度：',num2str(rsquare)]);
%% 精确定位结果

tic;
[z_pre,cost_pre]=SDIPointModelFit(signal,z_coa,rsquare,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol);
toc;

disp(['模型拟合法:',num2str(z_pre)]);
%%
function [z_coa,rsquare]=SDIPointModulFit(signal,lam,vsk_ini,NA)
% 这个函数用调制拟合的方法给出位移值的粗定位结果

% 输入 signal 为等波长空间的信号强度，double型N×1维向量
% 输入 lam 为等波长空间的采样点波长，double型N×1维向量
% 输入 vsk_ini 为等波长空间的光谱强度，double型N×1维向量
% 输入 NA 为光学系统的NA，double型1×1维向量

% 输出 z_coa 为位移值的粗定位结果，double型1×1维向量
% 输出 rsquare 为粗定位结果的置信程度，double型1×1维向量

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
%     [z_coa,rsquare]=SDIPointModulFit(signal,lam,vsk_ini,NA);

source_thr = 0.01;  %结算调制相位的光谱强度阈值
valid_div = vsk_ini > source_thr*max(vsk_ini); %光谱阈值上的数据对应的索引

interD = nan(size(signal)); %调制项分配内存
interD(valid_div) = signal(valid_div)./vsk_ini(valid_div); %计算调制项

fitting_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行拟合
valid = vsk_ini > fitting_thr*max(vsk_ini); %拟合阈值上的数据对应的索引

x_cho = 4*pi*1./lam(valid);
y_cho = interD(valid);

weight = vsk_ini(valid).^2; %将光谱数据强度当做权重
weight = weight./ max(weight); %权重归一化


[fitResult, gof] = fit(x_cho,y_cho,'fourier1','Weights',weight);
coeff = coeffvalues(fitResult);

c = sqrt(1-NA^2);
mu = 2*(1+c+c^2)/(3*(1+c)); %小NA修正系数

z_coa = abs(coeff(4))/mu; %粗定位结果
rsquare=gof.rsquare; %将拟合的有效系数当做粗定位的置信程度

end
%% 
function cost=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal)
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
%     cost=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal);


vsk_ini = vsk.*vk0.^2; %等波长域的光谱强度
fitting_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行拟合
valid = vsk_ini > fitting_thr*max(vsk_ini); %拟合阈值上的数据对应的索引
cost = nan*ones(size(z_gra)); %预先分配内存
for ii=1:length(z_gra)
    signal_gra = SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,z_gra(ii),system_pol); %计算信号
    cost(ii) = sum(abs(signal(valid)-signal_gra(valid)).^2,'all');%计算误差
end
end
%%
function [z_pre,cost_pre]=SDIPointModelFit(signal,z_coa,rsquare,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol)
% 这个函数用精细网格给出SDI信号位移值的精确定位结果

% 输入 signal 为要匹配的目标信号，double型 N×1维向量
% 输入 z_coa 为粗定位的结果，double型1×1维向量
% 输入 rsquare为粗定位结果的置信程度，double型1×1维向量，当rsquare<0.75时将先采用粗网格修正粗定位结果
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
%     [z_coa,rsquare]=SDIPointModulFit(signal,lam,vsk_ini,NA);
%     [z_pre,cost_pre]=SDIPointModelFit(signal,z_coa,rsquare,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol)

if rsquare<=0.75
    warning('粗定位效果较差，改用粗细网格搜索法精确定位');
    
    z_add = 2.5; %扫描的上范围
    z_min = max(z_coa-z_add,0); %搜索的下限
    z_max = z_coa+z_add; %搜索的上限
    z_peri = 5e-2; %50nm采样
    z_gra = z_min:z_peri:z_max; %搜索的区间
    
    cost=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal);
    
    [~,index]=min(cost);
    z_coa=z_gra(index); %粗网格修正粗定位结果
end

z_add = 0.5; %扫描的上下范围
z_min = max(z_coa-z_add,0); %搜索的下限
z_max = z_coa+z_add; %搜索的上限
z_peri = 1e-3; %1nm采样
z_gra = z_min:z_peri:z_max; %搜索的区间

cost=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal);

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
end