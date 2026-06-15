function q=Calcq(N,theta)
% 这个函数用来计算q值

% 输入 N 为该介质的复折射率，complex double型1×1变量
% 输入 theta 为空气中，光线的角度，double型1×1变量

% 输出 q 为该介质内部的q值，complex double型1×1变量

%      例：
%      N=1.5+0.02*1i;
%      theta=asin(0.15);
%      q=Calcq(N,theta)

q = sqrt(N^2-sin(theta)^2);
if imag(q)<0 || (abs(imag(q)) < 1e-12 && real(q) < 0)
    q=-q; %如果q值的虚部小于0，说明平方根多添加了一个负号，为了确保吸收介质内，光线的振幅随传播距离增加而减小，做此判断
end
end
