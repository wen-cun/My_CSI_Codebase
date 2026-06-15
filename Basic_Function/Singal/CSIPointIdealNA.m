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
% signal=signal*2/sin(theta_array(1,end))^2;

signal = signal./max(abs(signal),[],'all'); %信号归一化
end
