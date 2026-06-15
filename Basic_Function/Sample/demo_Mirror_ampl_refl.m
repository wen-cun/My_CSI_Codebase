%% 这个脚本用来计算反射镜振幅反射率，采用传递矩阵来计算
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
vk0_samp = 1/(z_peri *length(z_scan));

vsk_ini = gen_lightsource(vk0,1);%获取光源信号
vsk = vsk_ini./(vk0.^2);%波长波数域转换

figure(1);
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
NA = 0.15; %系统NA
theta_max=asin(NA); %最大NA对应的空气中光线角度theta
theta_peri=0.01; %角度theta的采样周期
theta_array = (0:theta_peri:theta_max)'; %theta采样数组

% sample_stru = {'SiO2',0.05;...
%                'Ag',0.2;...
%                'Boro33',inf}; %样品结构
% %% 读取样品每层介质的折射率
% layernum=size(sample_stru,1);%介质的数量
% Index_init_name = {'RefractiveIndex_'};
% Index_last_name = {'.csv'};
% Index_name = strcat(Index_init_name,sample_stru(:,1),Index_last_name);
% Index = cell(layernum,2);%折射率元胞数组
% for ii=1:layernum
%     Index_table = readtable(Index_name{ii});%读取介质折射率
%     Index{ii,1} = interp1(Index_table.wl,Index_table.n,1./vk0,'linear');%线性插值查询光源波长处的折射率
%     Index{ii,2} = interp1(Index_table.wl,Index_table.k,1./vk0,'linear');%线性插值查询光源波长处的消光系数
% end
% %% 计算每个波长，每个角度下的样品振幅反射率
% 
% r_te = nan*ones(length(vk0),length(theta_array)); %TE场振幅反射率矩阵
% r_tm = nan*ones(length(vk0),length(theta_array)); %Tm场振幅反射率矩阵
% 
% for ii=1:length(vk0)
%     for jj=1:length(theta_array)
%         p0_te = cos(theta_array(jj)); %空气中的TE场导纳
%         p0_tm = 1/cos(theta_array(jj)); %空气中的TM场导纳
%         Index_now = cellfun(@(v) v(ii), Index, 'UniformOutput', false); %该波长下，所有介质的折射率
%         [M_te,M_tm] = CalcTransferMatrix(vk0(ii),theta_array(jj),Index_now,sample_stru);
%         
%         N2 = (Index{layernum,1}(ii)+1i*Index{layernum,2}(ii));%基底介质的复折射率
%         q1=Calcq(N2,theta_array(jj)); %基底介质的q值
%         p2_te = q1; %基底的导纳，这里要用到折射率的消光系数
%         r_te(ii,jj) = CalcAmplitudeReflectivity(M_te,p0_te,p2_te);%振幅反射率
%         
%         p2_tm = N2^2/q1; %介质的导纳，这里要用到折射率的消光系数
%         r_tm(ii,jj) = CalcAmplitudeReflectivity(M_tm,p0_tm,p2_tm);%振幅反射率
%     end
% end
[r_te,r_tm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array);
figure(2);
plot(1./vk0,abs(r_tm(:,2)).^2,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('$\lambda$/$\mu$m','Interpreter','latex');
%% 计算反射镜振幅反射率的函数
function [r_te,r_tm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array)
% 这个脚本用来计算所有波长，所有角度下的参考镜振幅反射率

% 输入 vk0 为光源的波数序列，double型M*1维变量
% 输入 theta_array 为角度序列，double型N*1维变量

% 输出 r_te 为TE偏振下的振幅反射率，complex double型M*N维变量
% 输出 r_tm 为TM偏振下的振幅反射率，complex double型M*N维变量

%     例：
%     vk0=1./1.1:0.01:1./0.3;
%     theta=0:0.01:asin(0.15);
%     [r_te,r_tm]=CalcSampleAmplitudeReflectivity(vk0,theta_array)

sample_stru = {'SiO2',0.05;...
               'Ag',0.2;...
               'Boro33',inf}; %反射镜结构

layernum=size(sample_stru,1);%介质的数量
Index_init_name = {'RefractiveIndex_'};
Index_last_name = {'.csv'};
Index_name = strcat(Index_init_name,sample_stru(:,1),Index_last_name);
Index = cell(layernum,2);%折射率元胞数组
for ii=1:layernum
    Index_table = readtable(Index_name{ii});%读取介质折射率
    Index{ii,1} = interp1(Index_table.wl,Index_table.n,1./vk0,'linear');%线性插值查询光源波长处的折射率
    Index{ii,2} = interp1(Index_table.wl,Index_table.k,1./vk0,'linear');%线性插值查询光源波长处的消光系数
end

r_te = nan*ones(length(vk0),length(theta_array)); %TE场振幅反射率矩阵
r_tm = nan*ones(length(vk0),length(theta_array)); %Tm场振幅反射率矩阵

for ii=1:length(vk0)
    for jj=1:length(theta_array)
        p0_te = cos(theta_array(jj)); %空气中的TE场导纳
        p0_tm = 1/cos(theta_array(jj)); %空气中的TM场导纳
        Index_now = cellfun(@(v) v(ii), Index, 'UniformOutput', false); %该波长下，所有介质的折射率
        [M_te,M_tm] = CalcTransferMatrix(vk0(ii),theta_array(jj),Index_now,sample_stru);
        
        N2 = (Index{layernum,1}(ii)+1i*Index{layernum,2}(ii));%基底介质的复折射率
        q1=Calcq(N2,theta_array(jj)); %基底介质的q值
        p2_te = q1; %基底的导纳，这里要用到折射率的消光系数
        r_te(ii,jj) = CalcAmplitudeReflectivity(M_te,p0_te,p2_te);%振幅反射率
        
        p2_tm = N2^2/q1; %介质的导纳，这里要用到折射率的消光系数
        r_tm(ii,jj) = CalcAmplitudeReflectivity(M_tm,p0_tm,p2_tm);%振幅反射率
    end
end
end