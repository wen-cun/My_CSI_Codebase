%% 这个脚本是生成点CSI干涉信号的demo
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
NA = 0; %系统NA
theta_max=asin(NA); %最大NA对应的空气中光线角度theta
theta_peri=0.01; %角度theta的采样周期
theta_array =0:theta_peri:theta_max; %theta采样数组
%% 样品反射率
% sample_stru = {'SiO2',5;...
%                'Si',inf}; %样品结构
sample_stru = {'Si',inf}; %样品结构
[r_Se,r_Sm] = CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru); %计算样品TE场、TM场反射率
%% 参考镜反射率
[r_Me,r_Mm] = CalcMirrorAmplitudeReflectivity(vk0,theta_array); %计算参考镜TE场、TM场反射率
% %% 选择偏振模式，生成白光干涉信号
% % system.por = 'unpolar';%非偏振模式
% system.por = 'ideal';%非偏振模式
% 
% sample_dis = 3;%样品与参考镜之间的距离;
% switch system.por
%     case 'ideal'  %理想情况
%         A1 = 1; %样品振幅反射率
%         A2 = 1; %参考镜振幅反射率
%         if NA<eps
%             %0NA
%             vk0_mat = repmat(vk0,1,length(z_scan));
%             vsk_mat = repmat(vsk,1,length(z_scan)); %将列向量的光谱数据复制为多列的矩阵
%             z_scan_mat =repmat(z_scan,length(vk0),1);  %将行向量的扫描序列复制为多行的矩阵
%             inter=A1^2+A2^2+2*A1*A2*real(exp(1i*4*pi.*vk0_mat.*(z_scan_mat-sample_dis)));
%             signal = sum(vsk_mat.*inter,1); %每一列求和，生成信号
%         else
%             vk0_mat = repmat(vk0,1,length(theta_array));
%             vsk_mat = repmat(vsk,1,length(theta_array)); %将列向量的光谱数据复制为多列的矩阵
%             theta_mat = repmat(theta_array,length(vk0),1);
%             inter_D = vsk_mat.*(A1^2+A2^2).*sin(theta_mat).*cos(theta_mat);%干涉的直流项
%             I_D = sum(inter_D,'all'); %光强直流项
%             signal = nan*ones(size(z_scan)); %预先分配内存
%             inter_A = vsk_mat.*(r_Se.*conj(r_Me)+r_Sm.*conj(r_Mm)).*sin(theta_mat).*cos(theta_mat); %干涉的交流项
%             for ii = 1:length(z_scan)
%                 phase = exp(1i*4*pi.*vk0_mat.*cos(theta_mat).*(z_scan(ii)-sample_dis));
%                 signal(ii) = I_D + 2*real(sum(inter_A.*phase,'all'));
%             end
%             signal = signal+I_D; %信号加上直流项
%         end
%     case 'unpolar' %自然光/非偏振模式
%         if NA<eps
%             %0NA
%             vk0_mat = repmat(vk0,1,length(z_scan));
%             vsk_mat = repmat(vsk,1,length(z_scan)); %将列向量的光谱数据复制为多列的矩阵
%             r_Se_mat = repmat(r_Se,1,length(z_scan));
%             r_Sm_mat = repmat(r_Sm,1,length(z_scan));
%             r_Me_mat = repmat(r_Me,1,length(z_scan));
%             r_Mm_mat = repmat(r_Mm,1,length(z_scan)); %将列向量的反射率复制为多行的矩阵
%             z_scan_mat =repmat(z_scan,length(vk0),1);  %将行向量的扫描序列复制为多行的矩阵
%             inter = abs(r_Se_mat).^2+abs(r_Sm_mat).^2+abs(r_Me_mat).^2+abs(r_Mm_mat).^2 ...
%                 +2*real(r_Se_mat.*conj(r_Me_mat).*exp(1i*4*pi.*vk0_mat.*(z_scan_mat-sample_dis))...
%                 +r_Sm_mat.*conj(r_Mm_mat).*exp(1i*4*pi.*vk0_mat.*(z_scan_mat-sample_dis)));
%             signal = sum(vsk_mat.*inter,1); %每一列求和，生成信号
%         else
%             vk0_mat = repmat(vk0,1,length(theta_array));
%             vsk_mat = repmat(vsk,1,length(theta_array)); %将列向量的光谱数据复制为多列的矩阵
%             theta_mat = repmat(theta_array,length(vk0),1);
%             inter_D = vsk_mat.*(abs(r_Se).^2+abs(r_Sm).^2+abs(r_Me).^2+abs(r_Mm).^2).*sin(theta_mat).*cos(theta_mat);%干涉的直流项
%             
%             I_D = sum(inter_D,'all'); %光强直流项
%             signal = nan*ones(size(z_scan)); %预先分配内存
%             inter_A = vsk_mat.*A1*conj(A2).*sin(theta_mat).*cos(theta_mat); %干涉的交流项 存
%             for ii = 1:length(z_scan)
%                 phase = exp(1i*4*pi.*vk0_mat.*cos(theta_mat).*(z_scan(ii)-sample_dis));
%                 signal(ii) = I_D + 2*real(sum(inter_A.*phase,'all'));
%             end
%         end
% end
% signal = signal./max(abs(signal),[],'all'); %信号归一化
%% 选择偏振模式，生成白光干涉信号
% system_pol = 'unpolar';%非偏振模式
system_pol = 'unpolar';%理想偏振模式

sample_dis = 3;%样品与参考镜之间的距离;

signal=CSIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,theta_array,sample_dis,system_pol);
%% 查看信号

figure(2);
plot(z_scan,signal,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('z/$\mu m$','Interpreter','latex');
%%  CSI点式信号的生成
function signal=CSIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,theta_array,sample_dis,system_pol)
% 这个函数用来计算点式CSI信号

% 输入NA为系统的NA
% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入z_scan为扫描序列, double型1×K维向量，K为图像帧数
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量
% 输入system_pol为系统的偏振状态

%输出signal为CSI信号， double 1×K维向量

%     例：
%     z_peri=0.0719;%PZT单次移动间隔，um制
%     N_half = round(10 / z_peri); % 计算单侧点数
%     z_scan = (-N_half : N_half) * z_peri; % 这样生成的数组中心绝对是 0
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vk0=1./1.1:0.01:1./0.3;
%     NA=0.15；
%     theta_array=0:0.01:asin(NA);
%     sample_stru={'SiO2',5;...
%                   'Si',inf};
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array)
%     sample_dis=3;
%     system_pol='unpolar';
%     signal=CSIPointSingalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,theta_array,sample_dis,system_pol)

switch system_pol
     case 'ideal'  %理想情况
        A1 = 1; %样品振幅反射率
        A2 = 1; %参考镜振幅反射率
        if NA<eps
            %0NA
            signal=CSIPointIdeal0NA(vk0,vsk,A1,A2,z_scan,sample_dis);
        else
            signal=CSIPointIdealNA(vk0,vsk,theta_array,A1,A2,z_scan,sample_dis);
        end
    case 'unpolar' %自然光/非偏振模式
        if NA<eps
            %0NA
            signal=CSIPointUnpolar0NA(vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,sample_dis);
        else
            signal=CSIPointUnpolarNA(vk0,vsk,theta_array,r_Se,r_Sm,r_Me,r_Mm,z_scan,sample_dis);
        end
    otherwise
        error('未包含该偏振状态,仅支持: ideal,unpolar');
end
end
%% 自然光偏振0NA函数
function signal=CSIPointUnpolar0NA(vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,sample_dis)
% 这个函数用来计算自然光偏振NA=0下的点式CSI信号

% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入z_scan为扫描序列, double型1×K维向量，K为图像帧数
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量

%输出signal为CSI信号， double 1×K维向量

%     例：
%     z_peri=0.0719;%PZT单次移动间隔，um制
%     N_half = round(10 / z_peri); % 计算单侧点数
%     z_scan = (-N_half : N_half) * z_peri; % 这样生成的数组中心绝对是 0
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vk0=1./1.1:0.01:1./0.3;
%     theta_array=0:0.01:asin(0.15);
%     sample_stru={'SiO2',5;...
%                   'Si',inf};
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array)
%     sample_dis=3;
%     signal=CSIPointUnpolar0NA(vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,sample_dis)

vk0_mat = repmat(vk0,1,length(z_scan));
vsk_mat = repmat(vsk,1,length(z_scan)); %将列向量的光谱数据复制为多列的矩阵

r_Se_mat = repmat(r_Se,1,length(z_scan));
r_Sm_mat = repmat(r_Sm,1,length(z_scan));
r_Me_mat = repmat(r_Me,1,length(z_scan));
r_Mm_mat = repmat(r_Mm,1,length(z_scan)); %将列向量的反射率复制为多行的矩阵

z_scan_mat =repmat(z_scan,length(vk0),1);  %将行向量的扫描序列复制为多行的矩阵

inter = abs(r_Se_mat).^2+abs(r_Sm_mat).^2+abs(r_Me_mat).^2+abs(r_Mm_mat).^2 ...
    +2*real(r_Se_mat.*conj(r_Me_mat).*exp(1i*4*pi.*vk0_mat.*(z_scan_mat-sample_dis))...
    +r_Sm_mat.*conj(r_Mm_mat).*exp(1i*4*pi.*vk0_mat.*(z_scan_mat-sample_dis)));

signal = sum(vsk_mat.*inter,1); %每一列求和，生成信号

signal = signal./max(abs(signal),[],'all'); %信号归一化
end
%% 自然光偏振有限NA函数
function signal=CSIPointUnpolarNA(vk0,vsk,theta_array,r_Se,r_Sm,r_Me,r_Mm,z_scan,sample_dis)
% 这个函数用来计算自然光偏振有限NA下的点式CSI信号

% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入z_scan为扫描序列, double型1×K维向量，K为图像帧数
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量

%输出signal为CSI信号， double 1×K维向量

%     例：
%     z_peri=0.0719;%PZT单次移动间隔，um制
%     N_half = round(10 / z_peri); % 计算单侧点数
%     z_scan = (-N_half : N_half) * z_peri; % 这样生成的数组中心绝对是 0
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vk0=1./1.1:0.01:1./0.3;
%     theta_array=0:0.01:asin(0.15);
%     sample_stru={'SiO2',5;...
%                   'Si',inf};
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array)
%     sample_dis=3;
%     signal=CSIPointUnpolarNA(vk0,vsk,theta_array,r_Se,r_Sm,r_Me,r_Mm,z_scan,sample_dis)

vk0_mat = repmat(vk0,1,length(theta_array));
vsk_mat = repmat(vsk,1,length(theta_array)); %将列向量的光谱数据复制为多列的矩阵
theta_mat = repmat(theta_array,length(vk0),1);
inter_D = vsk_mat.*(abs(r_Se).^2+abs(r_Sm).^2+abs(r_Me).^2+abs(r_Mm).^2).*sin(theta_mat).*cos(theta_mat);%干涉的直流项
I_D = sum(inter_D,'all'); %光强直流项

signal = nan*ones(size(z_scan)); %预先分配内存
inter_A = vsk_mat.*(r_Se.*conj(r_Me)+r_Sm.*conj(r_Mm)).*sin(theta_mat).*cos(theta_mat); %干涉的交流项
for ii = 1:length(z_scan)
    phase = exp(1i*4*pi.*vk0_mat.*cos(theta_mat).*(z_scan(ii)-sample_dis));
    signal(ii) = I_D + 2*real(sum(inter_A.*phase,'all'));
end

signal = signal./max(abs(signal),[],'all'); %信号归一化
end
%% 理想偏振0NA函数
function signal=CSIPointIdeal0NA(vk0,vsk,A1,A2,z_scan,sample_dis)
% 这个函数用来计算理想偏振偏振NA=0下的点式CSI信号

% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入A1为样品振幅反射率, complex double型1×1维向量
% 输入A2为参考镜振幅反射率, complex double型1×1维向量
% 输入z_scan为扫描序列, double型1×K维向量，K为图像帧数
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量

%输出signal为CSI信号， double 1×K维向量

%     例：
%     z_peri=0.0719;%PZT单次移动间隔，um制
%     N_half = round(10 / z_peri); % 计算单侧点数
%     z_scan = (-N_half : N_half) * z_peri; % 这样生成的数组中心绝对是 0
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vk0=1./1.1:0.01:1./0.3;
%     theta_array=0:0.01:asin(0.15);
%     A1=1;
%     A2=2;
%     sample_dis=3;
%     signal=CSIPointIdeal0NA(vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,z_scan,sample_dis)

vk0_mat = repmat(vk0,1,length(z_scan));
vsk_mat = repmat(vsk,1,length(z_scan)); %将列向量的光谱数据复制为多列的矩阵
z_scan_mat =repmat(z_scan,length(vk0),1);  %将行向量的扫描序列复制为多行的矩阵
inter=abs(A1)^2+abs(A2)^2+2*real(A1*conj(A2)*exp(1i*4*pi.*vk0_mat.*(z_scan_mat-sample_dis)));
signal = sum(vsk_mat.*inter,1); %每一列求和，生成信号

signal = signal./max(abs(signal),[],'all'); %信号归一化
end
%% 理想偏振有限NA函数
function signal=CSIPointIdealNA(vk0,vsk,theta_array,A1,A2,z_scan,sample_dis)
% 这个函数用来计算理想偏振有限NA下的点式CSI信号

% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入A1为样品振幅反射率, complex double型1×1维向量
% 输入A2为参考镜振幅反射率, complex double型1×1维向量
% 输入z_scan为扫描序列, double型1×K维向量，K为图像帧数
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量

%输出signal为CSI信号， double 1×K维向量

%     例：
%     z_peri=0.0719;%PZT单次移动间隔，um制
%     N_half = round(10 / z_peri); % 计算单侧点数
%     z_scan = (-N_half : N_half) * z_peri; % 这样生成的数组中心绝对是 0
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vk0=1./1.1:0.01:1./0.3;
%     theta_array=0:0.01:asin(0.15);
%     A1=1;
%     A2=1;
%     sample_dis=3;
%     signal=CSIPointIdealNA(vk0,vsk,theta_array,A1,A2,z_scan,sample_dis)

vk0_mat = repmat(vk0,1,length(theta_array));
vsk_mat = repmat(vsk,1,length(theta_array)); %将列向量的光谱数据复制为多列的矩阵
theta_mat = repmat(theta_array,length(vk0),1);
inter_D = vsk_mat.*(abs(A1)^2+abs(A2)^2).*sin(theta_mat).*cos(theta_mat);%干涉的直流项
I_D = sum(inter_D,'all'); %光强直流项
signal = nan*ones(size(z_scan)); %预先分配内存
inter_A = vsk_mat.*A1*conj(A2).*sin(theta_mat).*cos(theta_mat); %干涉的交流项 存
for ii = 1:length(z_scan)
    phase = exp(1i*4*pi.*vk0_mat.*cos(theta_mat).*(z_scan(ii)-sample_dis));
    signal(ii) = I_D + 2*real(sum(inter_A.*phase,'all'));
end

signal = signal./max(abs(signal),[],'all'); %信号归一化
end