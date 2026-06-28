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

tic;
[fitResult, gof] = fit(x_cho,y_cho,'fourier1','Weights',weight);
toc;
coeff = coeffvalues(fitResult);

c = sqrt(1-NA^2);
mu = 2*(1+c+c^2)/(3*(1+c)); %小NA修正系数

z_coa = abs(coeff(4))/mu; %粗定位结果
rsquare=gof.rsquare; %将拟合的有效系数当做粗定位的置信程度

end