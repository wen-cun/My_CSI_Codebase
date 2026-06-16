function [r_te,r_tm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru)
% 这个脚本用来计算所有波长，所有角度下的样品振幅反射率

% 输入 vk0 为光源的波数序列，double型M*1维变量
% 输入 theta_array 为角度序列，double型1*N维变量
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
