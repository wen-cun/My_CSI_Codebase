function [z_pre,cost_pre]=SDIPointModelFit(signal,z_coa,rsquare,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol)
% 这个函数用精细网格给出SDI信号位移值的精确定位结果

% 输入 signal 为要匹配的目标信号，double型 N×1维向量
% 输入 z_coa 为粗定位的结果，double型1×1维向量
% 输入 rsquare为粗定位结果的置信程度，double型1×1维向量，当rsquare<0.75时将先采用粗网格修正粗定位结果
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
%     [z_coa,rsquare]=SDIPointModulFit(signal,lam,vsk_ini,NA);
%     [z_pre,cost_pre]=SDIPointModelFit(signal,z_coa,rsquare,NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol)

if rsquare<=0.75
    warning('粗定位效果较差，改用粗细网格搜索法精确定位');
    
    z_add = 2.5; %扫描的上范围
    z_min = max(z_coa-z_add,0); %搜索的下限
    z_max = z_coa+z_add; %搜索的上限
    z_peri = 5e-2; %50nm采样
    z_gra = z_min:z_peri:z_max; %搜索的区间
    
    cost=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal);
    
    [~,index]=min(cost);
    z_coa=z_gra(index); %粗网格修正粗定位结果
end

z_add = 0.5; %扫描的上下范围
z_min = max(z_coa-z_add,0); %搜索的下限
z_max = z_coa+z_add; %搜索的上限
z_peri = 1e-3; %1nm采样
z_gra = z_min:z_peri:z_max; %搜索的区间

cost=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_gra,signal);

[~,min_index] = min(cost); %误差最小值对应的索引
z_pre = z_gra(min_index);
cost_pre = cost(min_index);

if min_index == 1 || min_index == length(z_gra)
    warning('误差最小值位于搜索边界，建议扩大 搜索范围');
else
    %抛物线三点拟合
    % 最小值左右相邻的三个代价函数值
    cost_left   = cost(min_index-1);
    cost_center = cost(min_index);
    cost_right  = cost(min_index+1);
    
    denominator = cost_left - 2*cost_center + cost_right;
    
    % denominator>0 表示局部抛物线开口向上
    if isfinite(denominator) && denominator > eps(max(abs(cost)))
        delta_index = 0.5 * ...
            (cost_left - cost_right) / denominator;
        
        % 顶点原则上应该落在相邻两个采样点之间
        if abs(delta_index) <= 1
            z_pre = z_gra(min_index) + delta_index*z_peri;
            
            % 抛物线顶点对应的代价值
            cost_pre = cost_center ...
                - (cost_left-cost_right)^2/(8*denominator);
        else
            warning('抛物线顶点超出相邻采样区间，使用离散最小值。');
        end
    else
        warning('局部代价函数不满足开口向上的抛物线条件，使用离散最小值。');
    end
end
end
