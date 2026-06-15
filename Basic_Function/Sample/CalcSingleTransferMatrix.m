function M=CalcSingleTransferMatrix(p,delta)
% 这个函数用来计算单层介质传递矩阵，注：无论TE偏振，还是TM偏振，均满足此形式

% 输入 p 为介质的导纳,complex double型1×1型变量
% 输入 delta 为介质内的delta值，complex double型1×1型变量

% 输出 M为该介质的传递矩阵，complex double 2×2型变量

%     例：
%      N = 1.5+0.02i;%介质复折射率
%      q = Calcq(N,theta);
%      p1 = q; 
%      delta = 2*pi*vk0*5*q; 
%      M = CalcSingleTransferMatrix(p1,delta); %该介质的传递矩阵

M = [cos(delta),      -1i/p*sin(delta);...
         -1i*p*sin(delta), cos(delta)]; %该介质的传递矩阵
end

