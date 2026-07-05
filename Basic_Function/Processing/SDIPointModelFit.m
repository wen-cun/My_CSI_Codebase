function [z_pre,cost_pre]=SDIPointModelFit(signal,z_coa,valid,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol)
% 这个函数用精细网格给出SDI信号位移值的精确定位结果

% 输入 signal 为要匹配的目标信号，double型 N×1维向量
% 输入 z_coa 为粗定位的结果，double型1×1维向量
% 输入 valid 为计算信号内部的有效位置逻辑数组，Logical型N×1维向量
% 输入vk0为光谱采样点波数, double型N×1维向量，N为查询点个数
% 输入vsk为采样点强度（等波数空间）, double型N×1维向量，N为查询点个数
% 输入r_Se为样品TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为样品TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Me为参考镜TE偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入r_Sm为参考镜TM偏振反射率, complex double型N×M维向量，N为光谱采样点个数,M为角度序列的长度
% 输入theta_array为角度采样序列, double型1×M维向量
% 输入system_pol为系统的偏振状态

% 输出 z_pre 为精确定位的结果，double型1×1维向量
% 输出 cost_pre 为精确定位结果对应的损失值,double 型1×1维向量

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
%     [r_Se,r_Sm]=CalcSampleAmplitudeReflectivity(vk0,theta_array,sample_stru);
%     [r_Me,r_Mm]=CalcMirrorAmplitudeReflectivity(vk0,theta_array);
%     system_pol='ideal';
%     signal=SDIPointSignalGenerate(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,3,system_pol));
%     source_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行粗略定位，以降低噪声
%     valid = vsk_ini > source_thr*max(vsk_ini);
%     z_coa=SDIPointModulFit(signal,lam,vsk_ini,valid,NA);
%     [z_pre,cost_pre]=SDIPointModelFit(signal,z_coa,valid,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol)

  
z_add = 2.5; %扫描的上范围
z_min = max(z_coa-z_add,0); %搜索的下限
z_max = z_coa+z_add; %搜索的上限
z_peri = 1e-2; %10nm采样

z_minus_min = -z_coa-z_add; %搜索负半部分
z_minus_max = min(-z_coa+z_add,0); %搜索负半部分

z_gra = [z_minus_min:z_peri:z_minus_max,z_min:z_peri:z_max]; %搜索的区间

cost = CalcSDISeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal,valid);

local_mask = islocalmin(cost,'FlatSelection','center','SamplePoints',z_gra); %判断是否为局部极小值
idx_local = find(local_mask);

% 防止严格单调或边缘最小导致没有检测到局部极小值
[~,idx_global] = min(cost);

if isempty(idx_local)
    idx_local = idx_global;
elseif ~ismember(idx_global,idx_local)
    idx_local = [idx_local;idx_global];
end
% 按损失值排序
[~,order] = sort(cost(idx_local),'ascend');
idx_sorted = idx_local(order);
    
z_num = 6; %仅选取前六个极小值

z_add_fine = 1.5*z_peri; %精确定位范围采用粗网格定义的1个半以内。
z_step_fine = 1e-3;
num_fine_half = ceil(z_add_fine/z_step_fine);
fine_offsets = (-num_fine_half:num_fine_half)*z_step_fine;



z_para = nan*ones(1,z_num) ; %给每个候选位置抛物线位置分配内存
cost_para = nan*ones(1,z_num) ; %给每个候选位置抛物线拟合误差分配内存

for ii = 1:z_num
    z_cand = fine_offsets + z_gra(idx_sorted(ii));
    cost_cand = CalcSDISeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_cand,signal,valid);
    [z_para(ii),cost_para(ii),~] = RefineMinimumByParabola(z_cand,cost_cand);
end


[cost_pre,min_index] = min(cost_para); %确定最小抛物线候选中心
z_pre = z_para(min_index); %精确搜索的结果    
    
end