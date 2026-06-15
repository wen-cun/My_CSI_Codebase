%% 这个脚本用来计算样品振幅反射率，采用传递矩阵来计算
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

sample_stru = {'SiO2',5;...
               'Si',inf}; %样品结构
% sample_stru = {'Si',inf}; %样品结构
[r_te,r_tm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru);
%% 读取样品每层介质的折射率
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

% %双层循环计算每一个波长，每个角度下的te场振幅反射率
% for ii=1:length(vk0)
%     for jj=1:length(theta_array)
%         p0_te = cos(theta_array(jj)); %空气中的TE场导纳
%         p0_tm = 1/cos(theta_array(jj)); %空气中的TM场导纳
%         M_te = eye(2);%初始化TE场传递矩阵
%         M_tm = eye(2);%初始化TM场传递矩阵
%         if layernum==1  %仅有一层介质的情况
%             N=Index{1,1}(ii)+1i*Index{1,2}(ii);%介质复折射率
%             q=sqrt(N^2-sin(theta_array(jj))^2);
%             if imag(q)<0
%                     q=-q;
%             end
%             p2_te =q; %介质的导纳，这里要用到折射率的消光系数
%             r_te(ii,jj) = (p0_te-p2_te)/(p0_te+p2_te);
%             p2_tm = N^2/q; %介质的导纳，这里要用到折射率的消光系数
%             r_tm(ii,jj) = (p0_tm-p2_tm)/(p0_tm+p2_tm);
%             continue;
%         else
%             for kk=1:layernum-1 %内层循环计算每层介质传递矩阵
%                 %计算第kk层薄膜传递矩阵
%                 N = Index{kk,1}(ii)+1i*Index{kk,2}(ii);%介质复折射率
%                 q = sqrt(N^2-sin(theta_array(jj))^2);
%                 if imag(q)<0
%                     q=-q;
%                 end
%                 p1_te = q; %介质的TE场导纳，这里要用到折射率的消光系数
%                 delta = 2*pi*vk0(ii)*sample_stru{kk,2}*q; %介质内部的TE场delta
%                 M1_te = [cos(delta),      -1i/p1_te*sin(delta);...
%                     -1i*p1_te*sin(delta), cos(delta)]; %该介质的传递矩阵
%                 M_te=M_te*M1_te;%更新TE场传递矩阵
%                 
%                 p1_tm = N^2/q; %介质的TM场导纳，这里要用到折射率的消光系数
%                 M1_tm = [cos(delta),      -1i/p1_tm*sin(delta);...
%                     -1i*p1_tm*sin(delta), cos(delta)]; %该介质的传递矩阵
%                 M_tm=M_tm*M1_tm;%更新TM场传递矩阵
%             end
%             N2 = (Index{layernum,1}(ii)+1i*Index{layernum,2}(ii));%基底介质的负折射率
%             q1 = sqrt(N2^2-sin(theta_array(jj))^2);
%             if imag(q1)<0
%                 q1 = -q1;
%             end
%             p2_te = q1; %基底的导纳，这里要用到折射率的消光系数
%             gamma_te = (M_te(2,1)+M_te(2,2)*p2_te)/(M_te(1,1)+M_te(1,2)*p2_te); %结构等效介质导纳
%             r_te(ii,jj)=(p0_te-gamma_te)/(p0_te+gamma_te);%振幅反射率
%             
%             p2_tm = N2^2/q1; %介质的导纳，这里要用到折射率的消光系数
%             gamma_tm = (M_tm(2,1)+M_tm(2,2)*p2_tm)/(M_tm(1,1)+M_tm(1,2)*p2_tm); %结构等效介质导纳
%             r_tm(ii,jj)=(p0_tm-gamma_tm)/(p0_tm+gamma_tm);%振幅反射率
%         end
%     end
% end
%%
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
%% 计算膜层结构传递矩阵
function [M_te,M_tm] = CalcTransferMatrix(vk0,theta,index,sample_stru)
% 这个函数用来计算膜层结构的传递矩阵

% 输入 vk0 为波数，double型1×1变量
% 输入 theta 为角度，double型1×1变量
% 输入 index 为每层结构的单波长复折射率，指定为layernum×2型cell元胞数组，每个子胞内都为double型1×1变量
% 输入 sample_stru 为该膜层结构，指定为layernum×2型cell元胞数字，第一列子胞为每层膜的介质，第二列子胞为每层膜的厚度

% 输出 M_te 为TE偏振下，该膜层结构的总传递矩阵，complex double型2×2型变量
% 输出 M_tm 为TM偏振下，该膜层结构的总传递矩阵，complex double型2×2型变量

%     例：
%        vk0=1/0.5;
%        theta=asin(0.15);
%        index{1,1}=1.5；
%        index{1,2}=0.02；
%        index{2,1}=1.75;
%        index{2,2}=0.1;
%        sample_stru={'SiO2',5;'Si',inf};
%        [M_te,M_tm]=CalcTransferMatrix(vk0,theta,index,samp_stru)

layernum=size(sample_stru,1); %介质的总个数
M_te=eye(2);
M_tm=eye(2);
if layernum==1
    %直接返回单位矩阵
else
    for kk=1:layernum-1 %内层循环计算每层介质传递矩阵
                %计算第kk层薄膜传递矩阵
                N = index{kk,1}+1i*index{kk,2};%介质复折射率
                q = Calcq(N,theta);
                p1_te = q; %介质的TE场导纳，这里要用到折射率的消光系数
                delta = 2*pi*vk0*sample_stru{kk,2}*q; %介质内部的TE场delta
                M1_te = CalcSingleTransferMatrix(p1_te,delta); %该介质的传递矩阵 
                M_te=M_te*M1_te;%更新TE场传递矩阵
                
                p1_tm = N^2/q; %介质的TM场导纳，这里要用到折射率的消光系数
                M1_tm = CalcSingleTransferMatrix(p1_tm,delta); %该介质的传递矩阵 
                M_tm=M_tm*M1_tm;%更新TM场传递矩阵
   end
end
end

%% 计算q值的函数
function q=Calcq(N,theta)
% 这个函数用来计算q值

% 输入 N 为该介质的复折射率，complex double型1×1变量
% 输入 theta 为空气中，光线的角度，double型1×1变量

% 输出 q 为该介质内部的q值，complex double型1×1变量

%      例：
%      N=1.5+0.02*1i;
%      theta=asin(0.15);
%      q=Calcq(N,theta)

q = sqrt(N^2-sin(theta)^2);
if imag(q)<0 || (abs(imag(q)) < 1e-12 && real(q) < 0)
    q=-q; %如果q值的虚部小于0，说明平方根多添加了一个负号，为了确保吸收介质内，光线的振幅随传播距离增加而减小，做此判断
end
end

%% 计算单层介质传递矩阵的函数
function M=CalcSingleTransferMatrix(p,delta)
% 这个函数用来计算单层介质传递矩阵，注：无论TE偏振，还是TM偏振，均满足此形式

% 输入 p 为介质的导纳,complex double型1×1型变量
% 输入 delta 为介质内的delta值，complex double型1×1型变量

% 输出 M为该介质的传递矩阵，complex double 2×2型变量

%     例：
%      N = 1.5+0.02i;%介质复折射率
%      q = Calcq(N,theta);
%      p1 = q; 
%      delta = 2*pi*vk0*5*q; 
%      M = CalcSingleTransferMatrix(p1,delta); %该介质的传递矩阵

M = [cos(delta),      -1i/p*sin(delta);...
         -1i*p*sin(delta), cos(delta)]; %该介质的传递矩阵
end
%% 计算振幅反射率的函数
function r=CalcAmplitudeReflectivity(M,p0,p2)
% 这个函数用来计算整个膜层的振幅反射率，注：无论TE偏振，还是TM偏振，均满足此形式

% 输入 M 为薄膜结构的总传输矩阵，complex double 2×2型变量
% 输入 p0 为空气中的导纳,complex double 1×1型变量
% 输入 p1 为基底中的导纳,complex double 1×1型变量

% 输出 r 为薄膜结构总的振幅反射率，complex double 1×1型变量

%     例：
%     p0=cos(theta);
%     N=1.5+0.02i;
%     theta=asin(0.15);
%     q=Calcq(N,theta);
%     M=[1,0;0,1];
%     p2=2*pi*1/0.2*5*q;
%     r=CalcAmplitudeReflectivity(M,p0,p2)

Ye = (M(2,1)+M(2,2)*p2)/(M(1,1)+M(1,2)*p2); %等效介质导纳
r=(p0-Ye)/(p0+Ye);%电场振幅反射率
end

%% 计算样品所有波长，波数下振幅反射率的函数
function [r_te,r_tm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)
% 这个脚本用来计算所有波长，所有角度下的样品振幅反射率

% 输入 vk0 为光源的波数序列，double型M*1维变量
% 输入 theta_array 为角度序列，double型N*1维变量
% 输入 sample_stru 为膜层结构，K*2型元胞数组变量，K为总介质个数，其中第一列为子胞为膜层结构介质的名称，第二列为该膜层介质的厚度

% 输出 r_te 为TE偏振下的振幅反射率，complex double型M*N维变量
% 输出 r_tm 为TM偏振下的振幅反射率，complex double型M*N维变量

%     例：
%     vk0=1./1.1:0.01:1./0.3;
%     theta=0:0.01:asin(0.15);
%     sample_stru={'SiO2',5;...
%                   'Si',inf};
%     [r_te,r_tm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)

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