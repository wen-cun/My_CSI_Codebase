%% 这个脚本用来绘制课题规划0525PPT中的图
clear;
close all;
clc;
%% %绘制带振动的CSI信号
z_samp=0.0719;%PZT单次移动间隔，um制
% z_samp=0.05;%PZT单次移动间隔，um制

N_half = round(10 / z_samp); % 计算单侧点数
z_scan = (-N_half : N_half) * z_samp; % 这样生成的数组中心绝对是 0

fre=make_axis_freq(length(z_scan),z_samp,'101');
fre=fre(fre>2/1.1);%只选取大于1.1的部分；
vk0=(fre/2)';
vk0_samp=1/(z_samp*length(z_scan));

vsk_ini=gen_lightsource(vk0,1);%获取光源信号
vsk=vsk_ini./(vk0.^2);%波长波数域转换

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
%%
A1=1;
A2=1;
h_stand=3;
% h=h_stand*ones(1,length(z_scan))+0.5*rand(1,length(z_scan));%改动这里，这里一直在抖
f1=200;%频段在200Hz
f2=400;
h=h_stand*ones(1,length(z_scan))+0.25*sin(2*pi*f1*z_scan);%改动这里，这里一直在抖
h_mat=repmat(h,length(vk0),1);%行向量复制多次成矩阵
z_mat=repmat(z_scan,length(vk0),1);%行向量复制多次成矩阵
vk0_mat=repmat(vk0,1,length(z_scan));%列向量复制多次成矩阵
vsk_mat=repmat(vsk,1,length(z_scan));%列向量复制多次成矩阵
phase=4*pi*vk0_mat.*(z_mat-h_mat);%每一列代表一个Z项，每一行代表一个k
interf=A1.^2+A2.^2+2*A1*A2*cos(phase);%干涉项

%% %绘制SDI图
SDI=vsk_mat.*interf*vk0_samp;
figure(100);
for ii=1:40:length(h)
plot3((z_scan(ii)-h(ii))*ones(size(vk0)),vk0,SDI(:,ii),'LineWidth',1.5);
hold on;
end
hold off;
defaultAxes(3);
ylabel('k/$\mu m^{-1}$','Interpreter','latex');
%% 绘制CSI图
I=sum(vsk_mat.*interf*vk0_samp,1);%离散积分求和
I=I/max(I);%归一化
I=I-(mean(I(1:10)))';

%信号加高斯白噪声
SNR=100;%40dB的噪声
I=awgn(I,SNR,'measured');
I=I/max(I);%归一化
%40dB的高斯白噪声和实验中观察到的信号非常接近

figure(4);%带噪声CSI信号
plot(z_scan,I,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');

%% 给信号重排
z_scan_corr=z_scan-h;%修正后的高度
figure(5);
plot(z_scan_corr,'LineWidth',1.5);
xlabel('z/$\mu m$','Interpreter','latex');
I_corr=I;
for ii=1:length(z_scan)       %冒泡排序
    for jj=ii:length(z_scan)-1
        if z_scan_corr(jj)>z_scan_corr(jj+1)
           var_a=z_scan_corr(jj);
           var_b=I_corr(jj);
           z_scan_corr(jj)=z_scan_corr(jj+1);
           I_corr(jj)=I_corr(jj+1);
           z_scan_corr(jj+1)=var_a;
           I_corr(jj+1)=var_b;
        else
            continue;
        end
    end
end

figure(6);%修正后的CSI信号
plot(z_scan_corr+h_stand,I_corr,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');

figure(7);
plot(z_scan_corr,'LineWidth',1.5);
xlabel('z/$\mu m$','Interpreter','latex');

%均匀插值
z_scan_uni=z_scan-h_stand;%目标均匀的采样位置
I_uni=interp1(z_scan_corr,I_corr,z_scan_uni,'spline');%查询该点的强度
I_uni=(I_uni-mean(I_uni(1:10)))';%减去直流分量，并转为列向量

z_uni_mat=repmat(z_scan_uni,length(vk0),1);%行向量复制多次成矩阵
phase=4*pi*vk0_mat.*z_uni_mat;%每一列代表一个Z项，每一行代表一个k
interf=A1.^2+A2.^2+2*A1*A2*cos(phase);%干涉项
I_real=sum(vsk_mat.*interf*vk0_samp,1);%离散积分求和
I_real=(I_real-mean(I_real(1:10)))';%减去直流分量，并转为列向量

Color=GetColor(2,1);
figure(8);%修正后的CSI信号
plot(z_scan_uni,I_uni,'LineWidth',1.5,'Color',Color(1,:));
hold on;
plot(z_scan_uni,I_real,'--','LineWidth',1.5,'Color',Color(2,:));
hold off;
defaultAxes(2);
legend('spline','real');
xlabel('z/$\mu m$','Interpreter','latex');

figure(9);%修正后的CSI信号
plot(z_scan_uni+h_stand,I_real,'LineWidth',1.5,'Color',Color(1,:));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');

figure(10);%修正后的CSI信号
plot(z_scan_uni+h_stand,I_uni,'LineWidth',1.5,'Color',Color(1,:));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');
