%% %这个脚本用来仿真NA=0的白光干涉过程
clear;
close all;
clc;
%% %定义光源
z_samp=0.0719;%PZT单次移动间隔，um制

N_half = round(10 / z_samp); % 计算单侧点数
z_scan = (-N_half : N_half) * z_samp; % 这样生成的数组中心绝对是 0

fre=make_axis_freq(length(z_scan),z_samp,'101');
fre=fre(fre>2/1.1);%只选取大于1.1的部分；
vk0_ini=(fre/2)';
vk0_samp=1/(z_samp*length(z_scan));
[vk0,vSk_ini]=genLightSource(vk0_ini,10);%获取光源信号
vSk=vSk_ini./(vk0.^2);%波长波数域转换

figure(1);
tiledlayout(3,1);

nexttile;
plot(1./vk0,vSk_ini,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('$\lambda$/$\mu$m','Interpreter','latex');

nexttile;
plot(vk0,vSk_ini,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu m^{-1}$','Interpreter','latex');

nexttile;
plot(vk0,vSk,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu m^{-1}$','Interpreter','latex');
%% %生成白光信号
nktable=readtable("C:\Users\lenovo\Desktop\MATLAB\codebase_cct\project\w260201\RefractiveIndexINFO_HK9L.csv");%读取分光棱镜折射率数据
nk.lambda=nktable{:,1};
nk.n=nktable{:,2};
nk.k=1./nk.lambda;
n=interp1(nk.k,nk.n,vk0,'spline');%利用插值获得所需波数的BS折射率

figure(2);
plot(vk0,n,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu m^{-1}$','Interpreter','latex');

vk0_lim=[1.6,1.95];%选取部分光谱区域
k0 = 1.8; %光谱中心波数
%拟合折射率系数a,b;
vk0_fit=double(vk0(vk0>vk0_lim(1)&vk0<vk0_lim(2))-k0);%选取拟合的波数
n_fit=double(n(vk0>vk0_lim(1)&vk0<vk0_lim(2)))-1;%选取拟合的折射率
[fitresult,gof]=fit([vk0_fit,vk0_fit.^2],n_fit,'poly11');
coeff=coeffvalues(fitresult);%计算系数a,b,c
c=coeff(3);
b=coeff(2);
a=coeff(1);
n_result=coeff(3)*vk0_fit.^2+coeff(2)*vk0_fit+coeff(1);%计算拟合结果
% n=c*(vk0-k0).^2+b*(vk0-k0)+a+1;

Color=GetColor(2,1);%双绘图的颜色
figure(3);
plot(vk0_fit,n_fit,'--*','LineWidth',1.5,'Color',Color(1,:));
hold on;
plot(vk0_fit,n_result,'LineWidth',1.5,'Color',Color(2,:));
hold off;
defaultAxes(2);
legend('Original data','Linearfit','EdgeColor','none','Location','NorthWest');
xlabel('k-k0/$\mu m^{-1}$','Interpreter','latex');   

flaw=0;%BS缺陷水平
% h=a*flaw;%样品与参考镜的距离，确保白光信号包络中心基本在扫描中心位置
h=3;

A1=1;
A2=1;

%多重反射项
env_location=7;%多重反射包络的位置
flaw_multi_rig=0;
flaw_multi_lef=-flaw_multi_rig;
h_multi_rig=b*flaw_multi_rig+env_location;%主包络右侧的峰
h_multi_lef=b*flaw_multi_lef-env_location;%主包络左侧的峰
A3=0;
A4=0;%无寄生反射信号


z_mat=repmat(z_scan,length(vk0),1);%行向量复制多次成矩阵
vk0_mat=repmat(vk0,1,length(z_scan));%列向量复制多次成矩阵
vSk_mat=repmat(vSk,1,length(z_scan));%列向量复制多次成矩阵
n_mat=repmat(n,1,length(z_scan));%行向量复制多次成矩阵

%加入GDD、TOD、FOD
GDD=0;%二阶色散
TOD=0;%三阶色散 %无色散
phase1=1/2*4*pi^2*GDD*(vk0_mat-k0).^2;
phase2=1/6*8*pi^3*TOD*(vk0_mat-k0).^3;
phase=4*pi*vk0_mat.*(z_mat-h+(n_mat-1)*flaw)+phase1+phase2;%每一列代表一个Z项，每一行代表一个k
phase_multi_rig=4*pi*vk0_mat.*(z_mat-h_multi_rig+(n_mat-1)*flaw_multi_rig);%主包络右侧多重反射对应的相位
phase_multi_lef=4*pi*vk0_mat.*(z_mat-h_multi_lef+(n_mat-1)*flaw_multi_lef);%主包络左侧多重反射对应的相位

interf=A1.^2+A2.^2+2*A1*A2*cos(phase)+2*A3*A4*cos(phase_multi_rig)+2*A3*A4*cos(phase_multi_lef);%干涉项
I=sum(vSk_mat.*interf*vk0_samp,1);%离散积分求和
I=I/max(I);%归一化

%信号加高斯白噪声
SNR=40;%40dB的噪声
I=awgn(I,SNR,'measured');
I=I/max(I);%归一化
%40dB的高斯白噪声和实验中观察到的信号非常接近

figure(4);
plot(z_scan,I,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');
%% 信号加窗
I=(I-mean(I(1:10)))';%减去直流分量，并转为列向量

clip_num=15;%两侧各取10各点，注意和信号包络宽度有关
env_I = abs(hilbert(I));%信号包络
[~, I_maxindex] = max(env_I);
env_clip=env_I(I_maxindex-clip_num:I_maxindex+clip_num); %选取信号
z_clip=(z_scan(I_maxindex-clip_num:I_maxindex+clip_num))'; %选取坐标
z_center=sum(z_clip.*env_clip.^2)/sum(env_clip.^2);%质心所在位置

figure(5);
plot(z_scan,I,'LineWidth',1.5,'Color',Color(1,:));
hold on;
plot(z_scan,env_I,'--','LineWidth',1.5,'Color',Color(2,:));
scatter(z_scan(I_maxindex),env_I(I_maxindex),25,'filled');




% [~,I_maxindex]=max(I);
% env_clip=I(I_maxindex-clip_num:I_maxindex+clip_num); %选取信号
% z_clip=(z_scan(I_maxindex-clip_num:I_maxindex+clip_num))'; %选取坐标
% z_center=sum(z_clip.*env_clip.^2)/sum(env_clip.^2);%质心所在位置

plot(z_clip,env_clip,'--*','LineWidth',1.5);
hold off;
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');

tau=1;%超高斯窗标准差
order=3;%超高斯窗阶数
window=exp(-((z_scan'-z_center).^2/(2*tau^2)).^order);


figure(6);
tiledlayout(2,1);

nexttile;
plot(z_scan,I,'LineWidth',1.5,'Color',Color(1,:));
hold on;
plot(z_scan,window,'--','LineWidth',1.5,'Color',Color(2,:));
hold off;
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');

I=I.*window;

nexttile;
plot(z_scan,I,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');
%% 傅里叶变换，选取拟合的波段

spectrm=fftshift(fft(ifftshift(I)));%傅里叶变换
spectrm=spectrm/max(abs(spectrm));
fre=make_axis_freq(length(I),z_samp,'101');%生成频域坐标

figure(7);
tiledlayout(2,1);

nexttile;
plot(fre,abs(spectrm),'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('fre/Hz','Interpreter','latex');

nexttile;
plot(fre,unwrap(angle(spectrm)),'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('fre/Hz','Interpreter','latex');

vk0_choose_num=9;%波数域中，选择中心波数左右个9个点，共计19个点进行拟合
vk0_fre=(fre/2)';
[~,index]=min(abs(vk0_fre-k0));%找到最接近中心的1.8的数据对应的索引
k1=vk0_fre(index);%用最接近1.8数据的波数作为中心
vk0_fre_fit=vk0_fre(index-vk0_choose_num:index+vk0_choose_num);
spectrm_fit=spectrm(index-vk0_choose_num:index+vk0_choose_num);%对频谱进行选取
w = abs(spectrm_fit).^2;
w = w / max(w);

angle_fit=unwrap(angle(spectrm_fit));%提取相位

figure(8);
tiledlayout(2,1);

nexttile;
plot(vk0_fre_fit,abs(spectrm_fit),'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu^{-1}$m','Interpreter','latex');

nexttile;
plot(vk0_fre_fit,angle_fit,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu^{-1}$m','Interpreter','latex');
%% 最小二乘拟合求系数a,b (使用高效矩阵左除 \ )
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

% 还原拟合曲线用于绘图
angle_fitlm = phi - h1 * x; % 等价于 X_mat * coeffs

figure(9);
plot(vk0_fre_fit,angle_fit,'--*','LineWidth',1.5,'Color',Color(1,:));
hold on;
plot(vk0_fre_fit,angle_fitlm,'LineWidth',1.5,'Color',Color(2,:));
hold off;
defaultAxes(2);
legend('Original Data','Linear fit','EdgeColor','none');
xlabel('k/$\mu m^{-1}$','Interpreter','latex');


M=round(-2*k1*h1-phi/(2*pi));%整数M
h2=-phi/(4*pi*k1)-M/(2*k1);%精确的高度
%%
disp(['预设高度h:',num2str(h)]);
disp(['线性拟合高度h:',num2str(h1)]);
disp(['FDA处理高度h:',num2str(h2)]);