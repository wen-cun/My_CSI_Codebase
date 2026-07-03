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
    z_peri = 5e-2; %50nm采样
    z_gra = z_min:z_peri:z_max; %搜索的区间
    
    cost_pos=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal,valid);
    
    [~,index]=min(cost_pos);
    z_coa_pos=z_gra(index); %粗网格正半轴修正结果
    
    z_minus_min = -z_coa-z_add; %搜索负半部分
    z_minus_max = min(-z_coa+z_add,0); %搜索负半部分
    z_gra = z_minus_min:z_peri:z_minus_max; %搜索的区间
    
    cost_neg=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal,valid);
    
    [~,index]=min(cost_neg);
    z_coa_neg=z_gra(index); %粗网格负半轴修正结果
    
    z_add_fine = 0.25;
    z_step_fine = 1e-3;
    
    % 正半轴细网格
    z_pos_min = max(z_coa_pos-z_add_fine,0);
    z_pos_max = z_coa_pos+z_add_fine;
    z_fine_pos = z_pos_min:z_step_fine:z_pos_max;
    cost_fine_pos=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_fine_pos,signal,valid);
    
    
    % 负半轴细网格
    z_neg_min = z_coa_neg-z_add_fine;
    z_neg_max = min(z_coa_neg+z_add_fine,0);
    z_fine_neg = z_neg_min:z_step_fine:z_neg_max;
    cost_fine_neg=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_fine_neg,signal,valid);
    
    [z_pos_refined,cost_pos_refined,~] = ...
        RefineMinimumByParabola(z_fine_pos,cost_fine_pos);
    
    [z_neg_refined,cost_neg_refined,~] = ...
        RefineMinimumByParabola(z_fine_neg,cost_fine_neg);
    
    if cost_pos_refined <= cost_neg_refined
        z_pre = z_pos_refined;
        cost_pre = cost_pos_refined;
    else
        z_pre = z_neg_refined;
        cost_pre = cost_neg_refined;
    end
    
end