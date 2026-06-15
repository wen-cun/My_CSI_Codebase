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

