function [z_refined,cost_refined,success] =RefineMinimumByParabola(z_grid,cost_grid)
%这个函数用来给损失函数最小值处三点抛物线拟合

% 输入 z_grid 是要进行抛物线拟合的位移搜索范围，double型1×N维向量
% 输入 cost_grid 是要进行抛物线拟合的损失值，double型1×N维向量

% 输出 z_refined 是抛物线拟合后的中心，double型1×1维向量
% 输出 cost_refined 是抛物线中心对应的损失函数，double型1×1维向量
% 输出 success 是成功的标志，double型1×1维Logical 向量，返回1时表明成功进行抛物线拟合，返回0时表明失效

%     例：
%     z_grid=[0:1e-3:2.5];
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
%     z_grid=0:1e-3:5;
%     source_thr = 0.05; %仅选取光源强度在最大值0.05以上的值进行粗略定位，以降低噪声
%     valid = vsk_ini > source_thr*max(vsk_ini);
%     cost_grid=CalcSeekCost(NA,vk0,vsk,r_Se,r_Sm,r_Me,r_Mm,theta_array,system_pol,z_grid,signal,valid);
%     [z_refined,cost_refined,success]=RefineMinimumByParabola(z_grid,cost_grid);

z_grid = z_grid(:);
cost_grid = cost_grid(:);

[cost_min,index_min] = min(cost_grid);

z_refined = z_grid(index_min);
cost_refined = cost_min;
success = false;

% 最小值位于边缘时无法做三点拟合
if index_min == 1 || index_min == numel(z_grid)
    warning('局部最小值位于细搜索边界，使用离散最小值。');
    return;
end

idx = index_min-1:index_min+1;

% 以中心点为原点，提高数值稳定性
z0 = z_grid(index_min);
z_local = z_grid(idx)-z0;
cost_local = cost_grid(idx);

% cost = a*z^2+b*z+c
p = polyfit(z_local,cost_local,2);

% 必须开口向上
if ~all(isfinite(p)) || p(1) <= 0
    warning('局部代价函数不满足开口向上的抛物线条件。');
    return;
end

z_offset = -p(2)/(2*p(1));

% 顶点必须位于左右两个相邻采样点之间
if z_offset < z_local(1) || z_offset > z_local(3)
    warning('抛物线顶点超出局部三点范围。');
    return;
end

z_refined = z0+z_offset;
cost_refined = polyval(p,z_offset);

% SSE理论上不应小于0
cost_refined = max(cost_refined,0);

success = true;

end

