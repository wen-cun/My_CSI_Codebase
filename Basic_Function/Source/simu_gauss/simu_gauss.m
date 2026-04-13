function vsk=simu_gauss(vk0,mu,sigma,SNR)
%这个函数用来返回模拟的高斯光源

%输入vk0为查询点的波数,  变量类型：double型向量
%输入mu为高斯光源的中心，  变量类型：double
%输入sigma为高斯光源的标准差，  变量类型：double
%输入SNR为光源的信噪比，当指定为inf时，不在光源中添加噪声，  变量类型：double or inf

%输出vsk为查询点的强度，与输入的vk0同维度，  变量类型：double 型变量

%  例：
%  vk0=1./(0.3:0.001:0.8)';
%  vsk=simu_gauss(vk0,1/0.5,2,40);
%  返回中心在500nm，标准差为2，信噪比为40的模拟高斯光源

vsk_ini=1/sqrt(sigma*sqrt(2*pi)).*exp(-(vk0-mu).^2/(2*sigma.^2));
if SNR==inf
    vsk=vsk_ini;
else
    vsk=awgn(vsk_ini,SNR,'measured');
end
vsk=vsk./max(abs(vsk));
end

