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

vk0 = vk0(:);
vsk = vsk(:);
theta_array = theta_array(:).';

z_shape = size(z_scan);
z_scan = z_scan(:).';

N = numel(vk0);
M = numel(theta_array);
K = numel(z_scan);

if K == 0
    signal = reshape([],z_shape);
    return;
end

%% 基本输入检查

if numel(vsk) ~= N
    error('vk0和vsk的长度不一致。');
end

expected_size = [N,M];

if ~isequal(size(r_Se),expected_size) || ...
   ~isequal(size(r_Sm),expected_size) || ...
   ~isequal(size(r_Me),expected_size) || ...
   ~isequal(size(r_Mm),expected_size)

    error('反射率矩阵尺寸均应为N×M。');
end

if ~isscalar(sample_dis) || ~isfinite(sample_dis)
    error('sample_dis必须为有限标量。');
end

%% 预计算与扫描位置无关的量

cos_theta = cos(theta_array);

% 角度权重
angular_weight = sin(theta_array).*cos_theta;

% N×M光谱角度权重矩阵
weight_mat = vsk*angular_weight;

% 直流项
direct_term = ...
    abs(r_Se).^2 + abs(r_Sm).^2 + ...
    abs(r_Me).^2 + abs(r_Mm).^2;

I_D = real(sum(weight_mat.*direct_term,'all'));

% 复干涉系数
inter_A = weight_mat.*( ...
    r_Se.*conj(r_Me) + ...
    r_Sm.*conj(r_Mm));

% N×M相位系数
phase_factor = 4*pi*(vk0*cos_theta);

% 展开为列向量，便于矩阵乘法
phase_factor_vec = phase_factor(:);
inter_A_vec = inter_A(:);

num_phase = numel(phase_factor_vec);

%% 分块计算非均匀扫描位置对应的CSI信号

signal = zeros(1,K);

% 控制临时相位矩阵大小
% 可根据计算机内存调整，例如5e5、1e6、2e6
target_phase_elements = 1e6;

block_size = floor(target_phase_elements/max(num_phase,1));
block_size = max(block_size,1);
block_size = min(block_size,K);

for first_index = 1:block_size:K

    last_index = min(first_index+block_size-1,K);
    index_block = first_index:last_index;

    % 当前分块的非均匀扫描位置，1×B
    delta_z = z_scan(index_block)-sample_dis;

    % num_phase×B相位矩阵
    phase_block = exp( ...
        1i*(phase_factor_vec*delta_z));

    % 注意使用.'，避免对inter_A进行复共轭
    interference_block = ...
        inter_A_vec.'*phase_block;

    signal(index_block) = ...
        I_D + 2*real(interference_block);
end

%% 归一化

finite_signal = isfinite(signal);

if any(finite_signal)

    normalizer = max(abs(signal(finite_signal)));

    if isfinite(normalizer) && normalizer > 0
        signal(finite_signal) = ...
            signal(finite_signal)./normalizer;
    else
        warning('CSI信号归一化因子无效。');
    end

else
    warning('CSI信号中不存在有限值。');
end

%% 恢复输入扫描序列的形状

signal = reshape(signal,z_shape);

end


