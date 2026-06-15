function r=CalcAmplitudeReflectivity(M,p0,p2)
% 这个函数用来计算整个膜层的振幅反射率，注：无论TE偏振，还是TM偏振，均满足此形式

% 输入 M 为薄膜结构的总传输矩阵，complex double 2×2型变量
% 输入 p0 为空气中的导纳,complex double 1×1型变量
% 输入 p1 为基底中的导纳,complex double 1×1型变量

% 输出 r 为薄膜结构总的振幅反射率，complex double 1×1型变量

%     例：
%     p0=cos(theta);
%     N=1.5+0.02i;
%     theta=asin(0.15);
%     q=Calcq(N,theta);
%     M=[1,0;0,1];
%     p2=2*pi*1/0.2*5*q;
%     r=CalcAmplitudeReflectivity(M,p0,p2)

Ye = (M(2,1)+M(2,2)*p2)/(M(1,1)+M(1,2)*p2); %等效介质导纳
r=(p0-Ye)/(p0+Ye);%电场振幅反射率
end

