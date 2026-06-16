%% SDI信号计算函数的demo
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
% system_pol = 'ideal';%非偏振模式

sample_dis = 0.5;%样品与参考镜之间的距离;
% switch system_pol
%     case 'ideal'  %理想情况
%         A1 = 1; %样品振幅反射率
%         A2 = 1; %参考镜振幅反射率
%         if NA<eps
%             %0NA
%             vk0_mat = vk0; %光谱数据保持不变
%             inter=abs(A1)^2+abs(A2)^2+2*real(A1*conj(A2)*exp(1i*4*pi.*vk0_mat.*sample_dis));
%             signal = vsk.*inter; %光源直接乘以调制项
%         else
%             vk0_mat = repmat(vk0,1,length(theta_array)); %将列向量的波数数据复制为多列的矩阵
%             theta_mat = repmat(theta_array,length(vk0),1);
%             inter = (abs(A1)^2+abs(A2)^2+2*real(A1*conj(A2)*...
%                 exp(1i*4*pi.*vk0_mat.*sample_dis.*cos(theta_mat)))).*sin(theta_mat).*cos(theta_mat);
%             signal=vsk.*sum(inter,2);%每一行求和为调制项，光源直接乘以调制项
%         end
%     case 'unpolar' %自然光/非偏振模式
%         if NA<eps
%             %0NA
%             vk0_mat = vk0; %光谱数据保持不变
%             inter=abs(r_Se).^2+abs(r_Sm).^2+abs(r_Me).^2+abs(r_Mm).^2+...
%                 2*real((r_Se.*conj(r_Me)+r_Sm.*conj(r_Mm)).*exp(1i*4*pi.*vk0_mat.*sample_dis));
%             signal = vsk.*inter; %光源直接乘以调制项
%         else
%             vk0_mat = repmat(vk0,1,length(theta_array)); %将列向量的光谱数据复制为多列的矩阵
%             theta_mat = repmat(theta_array,length(vk0),1);
%             inter= (abs(r_Se).^2+abs(r_Sm).^2+abs(r_Me).^2+abs(r_Mm).^2+...
%                 2*real((r_Se.*conj(r_Me)+r_Sm.*conj(r_Mm)).*exp(1i*4*pi.*vk0_mat.*sample_dis.*cos(theta_mat))))...
%                 .*sin(theta_mat).*cos(theta_mat);
%             signal = vsk.*sum(inter,2); % 每一行求和为调制项，光源直接乘以调制项
%         end
% end
% signal = signal./max(abs(signal),[],'all'); %信号归一化


signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,sample_dis,system_pol);
%% 查看信号
figure();
plot(vk0,signal,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu m^{-1}$','Interpreter','latex');

figure();
plot(lam,signal.*vk0.^2,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('$\lambda$/$\mu m$','Interpreter','latex');
%%
interD=(signal+eps)./(vsk+eps);

figure();
plot(vk0,interD,'LineWidth',1.5,'Color',GetColor(1,1));
defaultAxes(2);
xlabel('k/$\mu m^{-1}$','Interpreter','latex');
%%  SDI点式信号的生成
function signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,sample_dis,system_pol)
% 这个函数用来计算点式SDI信号

% 输入NA为系统的NA
% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量
% 输入system_pol为系统的偏振状态

%输出signal为CSI信号， double 型N×1维向量


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
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array)
%     sample_dis=3;
%     system_pol='ideal';
%     signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,sample_dis,system_pol))

switch system_pol
    case 'ideal'  %理想情况
        A1 = 1; %样品振幅反射率
        A2 = 1; %参考镜振幅反射率
        if NA<eps
            %0NA
            signal = SDIPointIdeal0NA(vk0,vsk,A1,A2,sample_dis);
        else
            signal = SDIPointIdealNA(vk0,vsk,theta_array,A1,A2,sample_dis);
        end
    case 'unpolar' %自然光/非偏振模式
        if NA<eps
            %0NA
            signal=SDIPointUnpolar0NA(vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,sample_dis);
        else
            signal=SDIPointUnpolarNA(vk0,vsk,theta_array,r_Se,r_Sm,r_Me,r_Mm,sample_dis);
        end
    otherwise
        error('未包含该偏振状态,仅支持: ideal,unpolar');
end
end
%%
function signal=SDIPointIdeal0NA(vk0,vsk,A1,A2,sample_dis)
% 这个函数用来计算理想偏振偏振NA=0下的点式SDI信号

% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入A1为样品振幅反射率, complex double型1×1维向量
% 输入A2为参考镜振幅反射率, complex double型1×1维向量
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量

%输出signal为SDI信号， double 型N×1维向量

%     例：
%     lam_peri = 3e-4; %波长采样周期,um制
%     lam_lim = [0.3,1.1]; %波长采样范围
%     lam = lam_lim(1):lam_peri:lam_lim(2); %波长采样
%     lam = lam(:); %转为列向量
%     vk0 = 1./lam; %转为波数
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vsk=vsk./(max(abs(vsk)));
%     A1=1;
%     A2=2;
%     sample_dis=3;
%     signal=SDIPointIdeal0NA(vk0,vsk,A1,A2,sample_dis)

vk0_mat = vk0; %光谱数据保持不变
inter=abs(A1)^2+abs(A2)^2+2*real(A1*conj(A2)*exp(1i*4*pi.*vk0_mat.*sample_dis));
signal = vsk.*inter; %光源直接乘以调制项
signal = signal./max(abs(signal),[],'all'); %信号归一化
end
%%
function signal=SDIPointIdealNA(vk0,vsk,theta_array,A1,A2,sample_dis)
% 这个函数用来计算理想偏振偏振有限NA下的点式SDI信号

% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入A1为样品振幅反射率, complex double型1×1维向量
% 输入A2为参考镜振幅反射率, complex double型1×1维向量
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量

%输出signal为SDI信号， double 型N×1维向量

%     例：
%     lam_peri = 3e-4; %波长采样周期,um制
%     lam_lim = [0.3,1.1]; %波长采样范围
%     lam = lam_lim(1):lam_peri:lam_lim(2); %波长采样
%     lam = lam(:); %转为列向量
%     vk0 = 1./lam; %转为波数
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vsk=vsk./(max(abs(vsk)));
%     A1=1;
%     A2=2;
%     sample_dis=3;
%     NA = 0.3; %系统NA
%     theta_max=asin(NA); %最大NA对应的空气中光线角度theta
%     theta_peri=0.01; %角度theta的采样周期
%     theta_array = 0:theta_peri:theta_max; %theta采样数组
%     signal=SDIPointIdealNA(vk0,vsk,A1,A2,sample_dis)

vk0_mat = repmat(vk0,1,length(theta_array)); %将列向量的波数数据复制为多列的矩阵
theta_mat = repmat(theta_array,length(vk0),1);
inter = (abs(A1)^2+abs(A2)^2+2*real(A1*conj(A2)*...
    exp(1i*4*pi.*vk0_mat.*sample_dis.*cos(theta_mat)))).*sin(theta_mat).*cos(theta_mat);
signal=vsk.*sum(inter,2);%每一行求和为调制项，光源直接乘以调制项
signal = signal./max(abs(signal),[],'all'); %信号归一化
end
%%
function signal=SDIPointUnpolar0NA(vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,sample_dis)
% 这个函数用来计算自然偏振NA=0下的点式SDI信号

% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量

%输出signal为SDI信号， double 型N×1维向量

%     例：
%     lam_peri = 3e-4; %波长采样周期,um制
%     lam_lim = [0.3,1.1]; %波长采样范围
%     lam = lam_lim(1):lam_peri:lam_lim(2); %波长采样
%     lam = lam(:); %转为列向量
%     vk0 = 1./lam; %转为波数
%     vsk_ini = gen_lightsource(vk0,1);%获取光源信号
%     vsk = vsk_ini./(vk0.^2);%波长波数域转换
%     vsk=vsk./(max(abs(vsk)));
%     sample_stru={'SiO2',5;...
%                   'Si',inf};
%     NA = 0; %系统NA
%     theta_max=asin(NA); %最大NA对应的空气中光线角度theta
%     theta_peri=0.01; %角度theta的采样周期
%     theta_array = 0:theta_peri:theta_max; %theta采样数组
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array)
%     sample_dis=3;
%     signal=SDIPointUnpolar0NA(vk0,vsk,A1,A2,sample_dis)

vk0_mat = vk0; %光谱数据保持不变
inter=abs(r_Se).^2+abs(r_Sm).^2+abs(r_Me).^2+abs(r_Mm).^2+...
    2*real((r_Se.*conj(r_Me)+r_Sm.*conj(r_Mm)).*exp(1i*4*pi.*vk0_mat.*sample_dis));
signal = vsk.*inter; %光源直接乘以调制项
signal = signal./max(abs(signal),[],'all'); %信号归一化
end
%%
function signal=SDIPointUnpolarNA(vk0,vsk,theta_array,r_Se,r_Sm,r_Me,r_Mm,sample_dis)
% 这个函数用来计算自然偏振NA=0下的点式SDI信号

% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入sample_dis为样品与参考镜的距离, double型1×1维向量

%输出signal为SDI信号， double 型N×1维向量

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
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array)
%     sample_dis=3;
%     signal=SDIPointUnpolarNA(vk0,vsk,A1,A2,sample_dis)

vk0_mat = repmat(vk0,1,length(theta_array)); %将列向量的光谱数据复制为多列的矩阵
theta_mat = repmat(theta_array,length(vk0),1);
inter= (abs(r_Se).^2+abs(r_Sm).^2+abs(r_Me).^2+abs(r_Mm).^2+...
    2*real((r_Se.*conj(r_Me)+r_Sm.*conj(r_Mm)).*exp(1i*4*pi.*vk0_mat.*sample_dis.*cos(theta_mat))))...
    .*sin(theta_mat).*cos(theta_mat);
signal = vsk.*sum(inter,2); % 每一行求和为调制项，光源直接乘以调制项
signal = signal./max(abs(signal),[],'all'); %信号归一化
end